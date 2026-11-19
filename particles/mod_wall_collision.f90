!> Module for checking collisions between test particles and triangular
!> mesh representing the first wall. Can be used to estimate heat loads
!> or FILD signals.
!>
!> To use this module in a particle simulation, begin by initializing
!> the wall data.
!>
!> type(octree_triangle) :: octree ! The struct containing the wall data.
!> call mod_wall_collision_init('wall.h5', octree, max_depth)
!>
!> The HDF5 file contains a real array where wall elements are stored:
!> 'nodes' ! Triangle vertices as [x1_1, y1_1, z1_1, x2_1, y2_1, z2_1, ...
!>                                 x3_1, y3_1, z3_1, ..., x3_n, y3_n, z3_n]
!> 'n_triangle' ! Number of elements.
!>
!> For optimization, the computational volume will be divided recursively
!> into 8**max_depth subvolumes and collisions are only checked against
!> the elements that are in the same volume as particle ini/end points.
!> max_depth = 5 is a good value to be used.
!>
!> To check wall collisions at each time step, use (inside the simulation loop)
!> 
!> call mod_wall_collision_check(p, q, octree, wall_id, wall_pos)
!>
!> where p is the marker position at the beginning of the time step
!> (in cylindrical coordinates) and q at the end of the time step. The subroutine
!> returns id of the wall element that was hit (zero otherwise) and
!> the impact point as p + t * (q - p) where parameter t is determined by the collision
!> check. The wall element id corresponds to the position of the element in the input file.
!> 
!> Free resources with
!>
!> call mod_wall_collision_free(octree)
!>
!> You can export the data with
!>
!> call mod_wall_collision_export(sim, file)
!>
!> Note that this assumes that (for lost markers) you have stored the data as
!> i_elm = (-1) * wall_id. The minus sign is there just to notify that the marker is
!> lost.
!>
!> The exported data contains IDs of the wetted triangles and associated particle
!> (sum of weights) and energy loads (sum of weight * energy).
!>
module mod_wall_collision
use mpi
use hdf5_io_module
use constants, only: ATOMIC_MASS_UNIT,SPEED_OF_LIGHT
use mod_particle_types
use mod_particle_sim
use mod_coordinate_transforms, only: cartesian_to_cylindrical, cylindrical_to_cartesian
use mod_math_operators, only: cross_product

implicit none


private
public :: octree_node, mod_wall_collision_init, mod_wall_collision_check, mod_wall_collision_free, mod_wall_collision_export
public :: mod_wall_collision_check_restart, wall_impact_energy
public :: octree_triangle, octree_init, octree_free, octree_find, &
          mod_wall_collision_intersect ! For unit tests

!> Triangle representing one wall element.
type octree_triangle
   integer triangle_id   !< Unique identifier
   real(kind=8) :: v0(3) !< x,y,z coordinates of vertex A
   real(kind=8) :: v1(3) !< x,y,z coordinates of vertex B
   real(kind=8) :: v2(3) !< x,y,z coordinates of vertex C
end type octree_triangle

!> Node / volume belonging to the octree
type octree_node
   integer      :: depth         !< The depth of this node, for root depth = 1
   real(kind=8) :: boundary(2,3) !< [x_min, x_max; y_min, y_max, z_min, z_max] of the volume belonging this node
   type(octree_triangle), allocatable :: contained(:) !< Triangles contained within this node.
   type(octree_node), pointer         :: children(:) => null() !< Child nodes this element has. Leaf nodes never
                                                               !< assign this, so it must start out nullified.
end type octree_node

contains


function wall_impact_energy(sim, particle, mass, time) result(energy)
  type(particle_sim), intent(in) :: sim
  class(particle_base), intent(in) :: particle
  real*8, intent(in) :: mass, time
  real*8 :: energy, E(3), B(3), psi, U, gamma, st(2), r_found, z_found
  integer :: i_elm, ierr

  energy = 0.d0
  select type (particle)
  type is (particle_kinetic_relativistic)
    gamma = sqrt(1.d0 + dot_product(particle%p, particle%p)/(mass*SPEED_OF_LIGHT)**2)
  type is (particle_gc_relativistic)
    i_elm = 0
    call find_RZ(sim%fields%node_list, sim%fields%element_list, particle%x(1), particle%x(2), &
                r_found, z_found, i_elm, st(1), st(2), ierr)
    if (ierr /= 0) return
    call sim%fields%calc_EBpsiU(time, i_elm, st, particle%x(3), E, B, psi, U)
    gamma = sqrt(1.d0 + particle%p(1)**2/(mass*SPEED_OF_LIGHT)**2 &
                 + 2*particle%p(2)*norm2(B)/(mass*SPEED_OF_LIGHT**2))
  class default
    return
  end select
  energy = (gamma - 1.d0)*mass*ATOMIC_MASS_UNIT*SPEED_OF_LIGHT**2
end function wall_impact_energy

subroutine mod_wall_collision_check_restart(sim, file)
  type(particle_sim), intent(in) :: sim
  character(len=*), intent(in) :: file
  integer(HID_T) :: file_id
  integer :: ierr
  real*8 :: time
  logical :: has_time

  if (sim%my_id /= 0) return
  call HDF5_open(file, file_id, ierr)
  if (ierr /= 0) then
    write(*,*) 'Cannot open wall-load restart file ', trim(file)
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if
  has_time = .false.
  call h5lexists_f(file_id, '/time', has_time, ierr)
  if (ierr /= 0 .or. .not. has_time) then
    write(*,*) 'Missing checkpoint time in ', trim(file)
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if
  call HDF5_real_reading(file_id, time, '/time')
  call HDF5_close(file_id)
  if (abs(time - sim%time) > 1.d-14) then
    write(*,*) 'Particle and wall-load checkpoint times differ: ', sim%time, time
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if
end subroutine mod_wall_collision_check_restart

subroutine mod_wall_collision_export(sim, file, iangle_groups, impact_wall, impact_energy, previous_file)
  type(particle_sim), intent(in) :: sim    !< The particle sim struct
  character(len=*),   intent(in) :: file   !< Filename for output
  real*8, intent(in) :: iangle_groups(:,:) !< Angle of incidence (ngroup,nprt) for all markers in this process (dummy for confined markers)

  ! Explicit impacts select only the collisions accumulated since the previous export.
  integer, optional, intent(in) :: impact_wall(:,:)
  real*8, optional, intent(in) :: impact_energy(:,:)
  character(len=*), optional, intent(in) :: previous_file

  type :: saved_wall_load
     integer, allocatable :: ids(:)
     real*8, allocatable :: particles(:), energy(:), angle(:)
  end type saved_wall_load
  type(saved_wall_load), allocatable :: previous_load(:)
  integer :: n_previous
  integer, allocatable :: previous_ids(:)
  real*8, allocatable :: previous_particles(:), previous_energy(:), previous_angle(:)
  real*8, allocatable :: total_particles(:), total_energy(:), angle_energy(:)
  logical, allocatable :: wetted(:)
  logical :: have_previous

  integer              :: n_wetted               ! Number of wetted triangles
  integer, allocatable :: wetted_id(:)           ! IDs of wetted triangles (i.e. their position index in wall input)
  real*8, allocatable  :: particle_deposition(:) ! Deposited particles (assuming weight = number of real particles)
  real*8, allocatable  :: energy_deposition(:)   ! Deposited energy (J)
  real*8, allocatable  :: iangle_mean(:)         ! Mean angle of incidence weight by energy deposition of each marker

  integer, allocatable :: particles_per_proc(:), wall_id_all(:), wall_id(:)
  real*8, allocatable  :: weight_all(:), energy_all(:), weight(:), energy(:), iangle(:), iangle_all(:)

  character(len=5) :: group_name
  integer(HID_T) :: file_id, group_id, previous_id
  integer :: ierr, my_id, n_mpi, n_total, n_here, maxid, i, j, k

  call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, n_mpi, ierr)

  have_previous = .false.
  if (present(previous_file)) have_previous = len_trim(previous_file) > 0

  if ( my_id .eq. 0 ) then
     if (have_previous) then
        call HDF5_open(previous_file, previous_id, ierr)
        if (ierr /= 0) then
           write(*,*) 'Cannot open previous wall-load checkpoint ', trim(previous_file)
           call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
        end if
        if (allocated(sim%groups)) then
           allocate(previous_load(size(sim%groups)))
           do i=1,size(sim%groups)
              write(group_name, '(A,i0.3,A)') '/', i, '/'
              call HDF5_integer_reading(previous_id, n_previous, group_name//'nwetted')
              allocate(previous_load(i)%ids(n_previous), previous_load(i)%particles(n_previous), &
                       previous_load(i)%energy(n_previous), previous_load(i)%angle(n_previous))
              if (n_previous == 0) cycle
              call HDF5_array1D_reading_int(previous_id, previous_load(i)%ids, group_name//'wallid')
              call HDF5_array1D_reading(previous_id, previous_load(i)%particles, group_name//'particledepot')
              call HDF5_array1D_reading(previous_id, previous_load(i)%energy, group_name//'energydepot')
              call HDF5_array1D_reading(previous_id, previous_load(i)%angle, group_name//'angleofincidence')
           end do
        end if
        ! Release the input before creating output, including when the paths alias.
        call HDF5_close(previous_id)
     end if
     write(*,*) 'Creating wall-load checkpoint ', trim(file)
     call HDF5_create(file, file_id, ierr)
     if( ierr .ne. 0 ) then
        write(*,*) "Could not create wall-load checkpoint ", trim(file)
        call report_open_wallload_files()
        call MPI_Abort(MPI_COMM_WORLD, -1, ierr)
        return
     end if
     call HDF5_real_saving(file_id, sim%time, '/time')
  end if

  if (allocated(sim%groups)) then
     do i=1,size(sim%groups,1)
        if (.not. allocated(sim%groups(i)%particles)) then
           if( my_id .eq. 0 ) write(*,*) "WARNING: group ", i, " not allocated, exiting"
           call MPI_ABORT(MPI_COMM_WORLD, 1, ierr)
           return
        end if

        if( my_id .eq. 0 ) then
           write(group_name,"(A,i0.3,A)") "/", i, "/"
           call h5gcreate_f(file_id, group_name, group_id, ierr)
           call h5gclose_f(group_id, ierr)
        end if

        ! Find total number of markers and number of markers in this process
        allocate(particles_per_proc(0:n_mpi-1))
        particles_per_proc = 0
        n_here = size(sim%groups(i)%particles,1)
        call MPI_Gather(n_here,1,MPI_INTEGER,&
             particles_per_proc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        n_total = sum(particles_per_proc,1)

        ! Collect wall IDs, weights, and evaluate energy. TODO add other particle types
        allocate( wall_id(n_here), weight(n_here), energy(n_here), iangle(n_here) )
        wall_id = 0
       weight = 0.d0
       energy = 0.d0
       iangle = 0.d0
       do k=1,n_here
          if (present(impact_wall)) then
             wall_id(k) = impact_wall(i,k)
          else
             wall_id(k) = -sim%groups(i)%particles(k)%i_elm
          end if
          if (wall_id(k) <= 0) cycle
          weight(k) = sim%groups(i)%particles(k)%weight
          iangle(k) = iangle_groups(i,k)
          if (present(impact_energy)) then
             energy(k) = impact_energy(i,k)
          else
             energy(k) = wall_impact_energy(sim, sim%groups(i)%particles(k), sim%groups(i)%mass, sim%time)
          end if
       end do

       ! Gather data from other processes
        allocate( wall_id_all(n_total), weight_all(n_total), energy_all(n_total), iangle_all(n_total) )
        call MPI_Gatherv(wall_id(:), n_here, MPI_INTEGER, &
             wall_id_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_mpi-1)], &
             MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

        call MPI_Gatherv(weight(:), n_here, MPI_REAL8, &
             weight_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_mpi-1)], &
             MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        call MPI_Gatherv(energy(:), n_here, MPI_REAL8, &
             energy_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_mpi-1)], &
             MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        call MPI_Gatherv(iangle(:), n_here, MPI_REAL8, &
             iangle_all(:), particles_per_proc, [(sum(particles_per_proc(0:i-1),1), i=0,n_mpi-1)], &
             MPI_REAL8, 0, MPI_COMM_WORLD, ierr)

        deallocate(particles_per_proc, wall_id, energy, weight, iangle)

        if( my_id .eq. 0 ) then
           n_previous = 0
          if (have_previous) n_previous = size(previous_load(i)%ids)

          maxid = max(0, maxval(wall_id_all))
          if (n_previous > 0) then
             call move_alloc(previous_load(i)%ids, previous_ids)
             call move_alloc(previous_load(i)%particles, previous_particles)
             call move_alloc(previous_load(i)%energy, previous_energy)
             call move_alloc(previous_load(i)%angle, previous_angle)
             maxid = max(maxid, maxval(previous_ids))
          end if

          allocate(wetted(maxid), total_particles(maxid), total_energy(maxid), angle_energy(maxid))
          wetted = .false.
          total_particles = 0.d0
          total_energy = 0.d0
          angle_energy = 0.d0
          if (n_previous > 0) then
             wetted(previous_ids) = .true.
             total_particles(previous_ids) = previous_particles
             total_energy(previous_ids) = previous_energy
             angle_energy(previous_ids) = previous_energy*previous_angle
             deallocate(previous_ids, previous_particles, previous_energy, previous_angle)
          end if

          do k=1,n_total
             j = wall_id_all(k)
             if (j <= 0) cycle
             wetted(j) = .true.
             total_particles(j) = total_particles(j) + weight_all(k)
             total_energy(j) = total_energy(j) + weight_all(k)*energy_all(k)
             angle_energy(j) = angle_energy(j) + weight_all(k)*energy_all(k)*iangle_all(k)
          end do

          n_wetted = count(wetted)
          call HDF5_integer_saving(file_id, n_wetted, group_name//'nwetted')
          if (n_wetted > 0) then
             wetted_id = pack([(j, j=1,maxid)], wetted)
             particle_deposition = pack(total_particles, wetted)
             energy_deposition = pack(total_energy, wetted)
             iangle_mean = pack(angle_energy, wetted)
             where (energy_deposition /= 0.d0) iangle_mean = iangle_mean/energy_deposition
             call HDF5_array1D_saving_int(file_id, wetted_id, n_wetted, group_name//'wallid')
             call HDF5_array1D_saving(file_id, particle_deposition, n_wetted, group_name//'particledepot')
             call HDF5_array1D_saving(file_id, energy_deposition, n_wetted, group_name//'energydepot')
             call HDF5_array1D_saving(file_id, iangle_mean, n_wetted, group_name//'angleofincidence')
             deallocate(wetted_id, particle_deposition, energy_deposition, iangle_mean)
          end if
          deallocate(wetted, total_particles, total_energy, angle_energy)
       end if
       deallocate(wall_id_all, weight_all, energy_all, iangle_all)
     end do
  end if

  if (my_id .eq. 0) then
     call HDF5_close(file_id)
     write(*,*) 'Wall loads stored at ', file
  end if
  
end subroutine mod_wall_collision_export

!> Identify handles retaining files when checkpoint creation fails.
subroutine report_open_wallload_files()
  integer(HID_T), allocatable :: ids(:)
  integer(SIZE_T) :: n_open, name_length
  integer :: ierr, i, object_types
  character(len=1024) :: filename

  object_types = ior(H5F_OBJ_FILE_F, ior(H5F_OBJ_DATASET_F, H5F_OBJ_GROUP_F))
  call h5fget_obj_count_f(int(H5F_OBJ_ALL_F, HID_T), object_types, n_open, ierr)
  if (ierr /= 0 .or. n_open == 0) return
  allocate(ids(n_open))
  call h5fget_obj_ids_f(int(H5F_OBJ_ALL_F, HID_T), object_types, n_open, ids, ierr)
  if (ierr /= 0) return
  do i=1,int(n_open)
     call h5fget_name_f(ids(i), filename, name_length, ierr)
     if (ierr == 0) write(*,*) 'Open HDF5 handle ', ids(i), ' file ', trim(filename)
  end do
end subroutine report_open_wallload_files
!> Initialize the wall data.
!>
!> Wall data is read from HDF5 file and octree is initialized.
subroutine mod_wall_collision_init(file, max_depth, octree)
  implicit none

  character(len=*), intent(in)     :: file      !< Filename of the wall data
  integer, intent(in)              :: max_depth !< Maximum number of subdivisions to be made
  type(octree_node), intent(inout) :: octree    !< Initialized wall data / root node of the octree

  integer(HID_T)                     :: file_id
  real*8, allocatable                :: wall(:)
  type(octree_triangle), allocatable :: triangles(:)
  real*8  :: boundary(2,3), offset
  integer :: n, i, err

  call HDF5_open(file,file_id)
  call HDF5_integer_reading(file_id,n,'ntriangle')

  allocate(wall(9*n))
  call HDF5_array1D_reading(file_id, wall, 'nodes')
  call HDF5_close(file_id)

  ! Convert wall array into triangles and initialize the octree
  allocate(triangles(n))
  do i=1,n
     triangles(i)%triangle_id = i
     triangles(i)%v0 = (/ wall( (i-1)*9 + 1 ), wall( (i-1)*9 + 2 ), wall( (i-1)*9 + 3 ) /)
     triangles(i)%v1 = (/ wall( (i-1)*9 + 4 ), wall( (i-1)*9 + 5 ), wall( (i-1)*9 + 6 ) /)
     triangles(i)%v2 = (/ wall( (i-1)*9 + 7 ), wall( (i-1)*9 + 8 ), wall( (i-1)*9 + 9 ) /)
  end do

  ! Bounding box that contains all wall elements (with a little bit of padding)
  offset = 1.0d-3
  boundary(1,1) = minval( wall(1:9*n-2:3) ) - offset
  boundary(2,1) = maxval( wall(1:9*n-2:3) ) + offset
  boundary(1,2) = minval( wall(2:9*n-1:3) ) - offset
  boundary(2,2) = maxval( wall(2:9*n-1:3) ) + offset
  boundary(1,3) = minval( wall(3:9*n-0:3) ) - offset
  boundary(2,3) = maxval( wall(3:9*n-0:3) ) + offset
  deallocate(wall)

  call octree_init(triangles, max_depth, boundary, octree, err)
  deallocate(triangles)

end subroutine mod_wall_collision_init

!< Release resources used by the wall data / octree
subroutine mod_wall_collision_free(octree)
  implicit none

  type(octree_node), intent(inout) :: octree
  call octree_free(octree)
  
end subroutine mod_wall_collision_free

!< Check for wall collisions between particle ini and end points.
subroutine mod_wall_collision_check(p, q, octree, wall_id, wall_pos, iangle)
  implicit none

  real(kind=8), intent(in)      :: p(3), q(3) !< Particle initial and final positions in (r,z,phi)
  type(octree_node), intent(in) :: octree     !< Initialized octree data.
  integer, intent(out) :: wall_id     !< Zero when there is no collision or ID (position in the input array) of the 
                                      !< wall element the particle collided. Negative if particle outside octree volume.
  real*8, intent(out)  :: wall_pos(3) !< The intersection point (r,z,phi) if applicable
  real*8, intent(out)  :: iangle      !< Angle of incidence [rad] if applicable

  type(octree_node), pointer :: nodep, nodeq
  real*8                     :: pxyz(3), qxyz(3), t
  integer                    :: i, err

  wall_id = 0

  ! Find nodes for both initial and end points.
  pxyz = cylindrical_to_cartesian(p)
  qxyz = cylindrical_to_cartesian(q)
  call octree_find(octree, pxyz, nodep, err)
  if(err .gt. 0) then
     wall_id = -1
     return
  end if

  call octree_find(octree, qxyz, nodeq, err)
  if(err .gt. 0) then
     wall_id = -1
     return
  end if

  ! Check collisions for the initial point's node
  do i = 1,size(nodep%contained)
     call mod_wall_collision_intersect(pxyz, qxyz, nodep%contained(i)%v0, &
          nodep%contained(i)%v1, nodep%contained(i)%v2, t, iangle)
     if (t .ge. 0.d0) then
        wall_id = nodep%contained(i)%triangle_id
        wall_pos = cartesian_to_cylindrical( pxyz + t * ( qxyz - pxyz ) )
        exit
     end if
  end do
  

  ! If no collision and ini and end points are in separate nodes,
  ! check collisions also for the other node.
  if( wall_id .eq. 0 .and. .not. associated(nodep, nodeq) ) then
  
     do i = 1,size(nodeq%contained)
        call mod_wall_collision_intersect(pxyz, qxyz, nodeq%contained(i)%v0, &
             nodeq%contained(i)%v1, nodeq%contained(i)%v2, t, iangle)
        if (t .ge. 0.d0) then
           wall_id = nodeq%contained(i)%triangle_id
           wall_pos = cartesian_to_cylindrical( pxyz + t * ( qxyz - pxyz ) )
           exit
        end if
     end do
     
  end if
  
end subroutine mod_wall_collision_check

!< Check for collision between line segment and triangle in 3D.
!< 
!< This subroutine implements Moeller-Trumbore algorithm.
subroutine mod_wall_collision_intersect(p, q, v0, v1, v2, t, iangle)
  implicit none

  real(kind=8), dimension(3), intent(in) :: p, q       !< Ini and end (x,y,z) positions defining the segment.
  real(kind=8), dimension(3), intent(in) :: v0, v1, v2 !< Triangle vertices (x,y,z).
  real(kind=8), intent(out) :: t !< Parameter defining the intersection point as p + t * (q -p) or negative if
                                 !< there is no intersection.
  real(kind=8), intent(out) :: iangle !< Angle of incidence [rad] if applicable

  real(kind=8) :: e1(3), e2(3), e3(3), pq(3)
  real(kind=8) :: ao(3), dao(3)
  real(kind=8) :: det, invdet, u, v, length, epsilon

  epsilon = 1.d-6

  e1 = v1 - v0
  e2 = v2 - v0
  e3 = cross_product(e1,e2)

  pq     = q - p
  length = norm2(pq)
  det    = -dot_product( pq / length, e3 )
  invdet = 1.d0 / det

  ao  = p - v0
  dao = cross_product(ao, pq / length )
  u =  dot_product(e2,dao) * invdet
  v = -dot_product(e1,dao) * invdet
  t =  dot_product(ao,e3)  * invdet

  if (abs(det) >= epsilon .and. t >= 0.d0 .and. u >= 0.d0 .and. v >= 0.d0 .and. &
         (u+v) <= 1.d0 .and. t .lt. length ) then
     t = t / length

     ! Make sure triangle normal and incident vector point to same direction 
     if( dot_product(e3, pq) .lt. 0.d0 ) e3 = -e3
     ! Incident angle is then just the angle between these two
     iangle = acos( dot_product(pq, e3) / ( length * norm2(e3) ) )
  else
     t = -1.d0
  end if

end subroutine mod_wall_collision_intersect

!< Initialize octree
subroutine octree_init(triangles, max_depth, boundary, node, err)
  implicit none

  type(octree_triangle), intent(in) :: triangles(:)  !< Array of wall triangles
  integer, intent(in)               :: max_depth     !< Maximum number of divisions 
  real(kind = 8), intent(in)        :: boundary(2,3) !< Bounding box
  type(octree_node), intent(inout)  :: node          !< Root node
  integer, intent(out)              :: err           !< Non zero if some triangles are outside bounding box

  integer i, j

  node%depth = 1
  node%boundary = boundary

  ! Check that all triangles are within the bounding box before constructing the octree
  err = 0
  do i = 1, size(triangles)
     if ( .not. contains_triangle(node%boundary, triangles(i)) ) then
        err = 1
        return
     end if
  end do

  call octree_construction(triangles, max_depth, node)

end subroutine octree_init

!< Routine for initializing nodes recursively
recursive subroutine octree_construction(triangles, max_depth, node)
  implicit none

  type(octree_triangle), intent(in) :: triangles(:)
  integer, intent(in)               :: max_depth
  type(octree_node), intent(inout)  :: node
  integer :: i, j
  integer :: contained_tri

  ! Find how many triangles are contained within this node
  contained_tri = 0

  do i = 1, size(triangles)
     if ( .not. contains_triangle(node%boundary, triangles(i)) ) cycle
     contained_tri = contained_tri + 1
  end do

  allocate(node%contained(contained_tri))

  ! No triangles in this node -> no need to make children
  if( contained_tri .eq. 0 ) return

  ! Repeat the loop, this time assigning all the triangles that are within this node
  j=1
  do i = 1, size(triangles)
     if ( .not. contains_triangle(node%boundary, triangles(i)) ) cycle
     node%contained(j)%triangle_id = triangles(i)%triangle_id
     node%contained(j)%v0 = triangles(i)%v0
     node%contained(j)%v1 = triangles(i)%v1
     node%contained(j)%v2 = triangles(i)%v2
     j = j + 1
  end do

  ! Return if we have reached the max depth
  if(node%depth==max_depth) return

  ! Make children and free the memory allocated for triangles as those will be contained within child nodes
  call node_division(node)
  do i = 1, 8
     call octree_construction(node%contained, max_depth, node%children(i))
  end do
  deallocate(node%contained)

end subroutine octree_construction

!< Function for checking if a given axis aligned box contains the given triangle
!<
!< This function is conservative. In reality it doesn't check if the triangle is
!< inside the box but simply checks if the bounding box of the triangle and the given
!< box overlap. For now this results in some extra collision checks so feel free to improve
!< this algorithm if you want.
logical function contains_triangle(nbox, triangle)
  real*8, intent(in)                :: nbox(2,3)
  type(octree_triangle), intent(in) :: triangle

  real*8 :: bbox(2,3) ! Bounding box for the triangle
  bbox(1,1) = minval((/triangle%v0(1), triangle%v1(1), triangle%v2(1)/))
  bbox(2,1) = maxval((/triangle%v0(1), triangle%v1(1), triangle%v2(1)/))
  bbox(1,2) = minval((/triangle%v0(2), triangle%v1(2), triangle%v2(2)/))
  bbox(2,2) = maxval((/triangle%v0(2), triangle%v1(2), triangle%v2(2)/))
  bbox(1,3) = minval((/triangle%v0(3), triangle%v1(3), triangle%v2(3)/))
  bbox(2,3) = maxval((/triangle%v0(3), triangle%v1(3), triangle%v2(3)/))

  contains_triangle = ( &
       ( ( bbox(1,1) .le. nbox(1,1) .and. nbox(1,1) .le. bbox(2,1) ) .or.    &
       ( bbox(1,1) .le. nbox(2,1) .and. nbox(2,1) .le. bbox(2,1) ) .or.    &
       ( nbox(1,1) .le. bbox(1,1) .and. bbox(1,1) .le. nbox(2,1) ) .or.    &
       ( nbox(1,1) .le. bbox(2,1) .and. bbox(2,1) .le. nbox(2,1) ) ) .and. &
       ( ( bbox(1,2) .le. nbox(1,2) .and. nbox(1,2) .le. bbox(2,2) ) .or.    &
       ( bbox(1,2) .le. nbox(2,2) .and. nbox(2,2) .le. bbox(2,2) ) .or.    &
       ( nbox(1,2) .le. bbox(1,2) .and. bbox(1,2) .le. nbox(2,2) ) .or.    &
       ( nbox(1,2) .le. bbox(2,2) .and. bbox(2,2) .le. nbox(2,2) ) ) .and. &
       ( ( bbox(1,3) .le. nbox(1,3) .and. nbox(1,3) .le. bbox(2,3) ) .or.    &
       ( bbox(1,3) .le. nbox(2,3) .and. nbox(2,3) .le. bbox(2,3) ) .or.    &
       ( nbox(1,3) .le. bbox(1,3) .and. bbox(1,3) .le. nbox(2,3) ) .or.    &
       ( nbox(1,3) .le. bbox(2,3) .and. bbox(2,3) .le. nbox(2,3) ) )       &
       )

end function contains_triangle

!< Divide node volume and allocated and initialize the children
subroutine node_division(node)
  implicit none

  type(octree_node), intent(inout), target :: node
  integer :: i, j, k, l
  real(kind=8) :: offset

  ! Add just a little bit of offset so that the boxes overlap
  ! and we don't end up with (numerical) gaps
  offset = 1.0d-3

  allocate(node%children(8))
  l = 1
  do k = 1,2
     do j = 1,2
        do i = 1,2
           node%children(l)%depth = node%depth + 1
           
           node%children(l)%boundary(1,1) = node%boundary(1,1) &
                + (i-1) * (node%boundary(2,1) - node%boundary(1,1)) * 0.5d0 - offset
           node%children(l)%boundary(2,1) = node%boundary(2,1) &
                - (2-i) * (node%boundary(2,1) - node%boundary(1,1)) * 0.5d0 + offset
           node%children(l)%boundary(1,2) = node%boundary(1,2) &
                + (j-1) * (node%boundary(2,2) - node%boundary(1,2)) * 0.5d0 - offset
           node%children(l)%boundary(2,2) = node%boundary(2,2) &
                - (2-j) * (node%boundary(2,2) - node%boundary(1,2)) * 0.5d0 + offset
           node%children(l)%boundary(1,3) = node%boundary(1,3) &
                + (k-1) * (node%boundary(2,3) - node%boundary(1,3)) * 0.5d0 - offset
           node%children(l)%boundary(2,3) = node%boundary(2,3) &
                - (2-k) * (node%boundary(2,3) - node%boundary(1,3)) * 0.5d0 + offset

           l = l + 1
        end do
     end do
  end do
end subroutine node_division

!< Find the node containing the point
subroutine octree_find(node, x, nodeout, err)
  implicit none

  type(octree_node), target, intent(in)     :: node    !< The octree root node
  real(kind=8), intent(in)                  :: x(3)    !< x,y,z coordinates of the query point
  type(octree_node), pointer, intent(inout) :: nodeout !< Node containg the point
  integer, intent(out)                      :: err     !< Non-zero if the point is outside the octree

  if(  (x(1)) > node%boundary(1,1) .and. &
       (x(1)) < node%boundary(2,1) .and. &
       (x(2)) > node%boundary(1,2) .and. &
       (x(2)) < node%boundary(2,2) .and. &
       (x(3)) > node%boundary(1,3) .and. &
       (x(3)) < node%boundary(2,3)) then

     err = 0
     call octree_traverse(node, x, nodeout)
  else
     err = 1
  end if

end subroutine octree_find

!< Find recursively the node where given point belongs to
recursive subroutine octree_traverse(node, x, nodeout)
  implicit none

  type(octree_node), target, intent(in)     :: node
  real(kind=8), intent(in)                  :: x(3)
  type(octree_node), pointer, intent(inout) :: nodeout

  integer :: i

  ! This node has no children, return the node
  if( .not. associated(node%children) ) then
     nodeout => node
  else
     ! Find the child to which the given point belongs to
     do i = 1,8
        if(  (x(1)) > node%children(i)%boundary(1,1) .and. &
             (x(1)) < node%children(i)%boundary(2,1) .and. &
             (x(2)) > node%children(i)%boundary(1,2) .and. &
             (x(2)) < node%children(i)%boundary(2,2) .and. &
             (x(3)) > node%children(i)%boundary(1,3) .and. &
             (x(3)) < node%children(i)%boundary(2,3)) then

           ! Child found (no need to loop further), repeat the process for the child node
           call octree_traverse(node%children(i), x, nodeout)
           exit
        end if
     end do
  end if

end subroutine octree_traverse

!< Free resources used by the octree
!<
!< This routine should be called only for the root node.
subroutine octree_free(node)
  implicit none
  type(octree_node), intent(inout) :: node

  call clean_node(node)
  node%depth    = 0
  node%boundary = 0

end subroutine octree_free

!< Free node resources (children and triangles) recursively
recursive subroutine clean_node(node)
  implicit none

  type(octree_node), intent(inout) :: node
  integer i

  if(associated(node%children)) then
     if(size(node%children)==8) then
       do i = 1,8
         call clean_node(node%children(i))

         if( allocated(node%children(i)%contained) ) then
           deallocate(node%children(i)%contained)
         end if
       end do
       deallocate(node%children)
     end if
  end if
end subroutine clean_node

end module mod_wall_collision
