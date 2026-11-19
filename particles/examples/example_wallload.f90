!> Assess wall loads using relativistic guiding centers.
!>
!> This example reads input markers from part_restart.h5.
!> Checkpoints are written at numbered fluid snapshots divisible by nout_particles.
!> With restart=.true., copy a matching pair to part_restart.h5 and wallload_restart.h5.
!> nstep_particles sets the duration of each invocation from the loaded particle time.
!> You will also need the wall input. See the wiki page "How to assess wall loads with the particle tracker?"
!> for details.
!>
Program  example_wallload

use particle_tracer
use mod_particle_io
use mod_particle_diagnostics
use mod_fields_linear   
use mod_fields_hermite_birkhoff
use mod_gc_relativistic
use mod_kinetic_relativistic
use mod_wall_collision
use constants, only: PI
use phys_module, only: tstep_particles, nstep_particles, nout_particles, restart, restart_particles, &
                       wall_collision_model, max_depth_wall
use mpi, only: wallload_mpi_integer => MPI_INTEGER, wallload_comm => MPI_COMM_WORLD, &
               wallload_abort => MPI_Abort
                  
implicit none

interface
  subroutine find_RZ(node_list, element_list, R_find, Z_find, R_out, Z_out, ielm_out, s_out, t_out, ifail)
    use data_structure, only: type_node_list, type_element_list
    implicit none
    type(type_node_list), intent(in) :: node_list
    type(type_element_list), intent(in) :: element_list
    real*8, intent(in) :: R_find, Z_find
    real*8, intent(out) :: R_out, Z_out, s_out, t_out
    integer, intent(inout) :: ielm_out
    integer, intent(out) :: ifail
  end subroutine find_RZ
end interface

external :: MPI_Bcast

! Set up the simulation variables
real(kind=8)                      :: timesteps(1)
real(kind=8)                      :: target_time, t_particle, sim_time, push_dt
integer(kind=4)                   :: n_part, i, j, k, n_steps, ifail, wall_id, fluid_step, ierr
type(write_particle_diagnostics)  :: diag
real(kind=8),dimension(3)         :: pos_prev, wall_pos
real*8, allocatable :: iangle(:,:), impact_energy(:,:)
integer, allocatable :: impact_wall(:,:)
character(len=64) :: previous_wallload, particle_file, wallload_file
character(len=16) :: step_name
logical :: at_fluid_step, fields_exhausted

type(particle_kinetic_relativistic), allocatable :: prtkin(:)
real*8 :: rnd(1), psi, U, B(3), E(3)

type(octree_node) :: wall
class(action), pointer :: field_action

! Read and broadcast input parameters before loading the wall model.
call sim%initialize(num_groups=1)
call mod_wall_collision_init(wall_collision_model,max_depth_wall,wall)

! Read particle file (atomic rates are not needed for orbit tracking).
call read_simulation_hdf5(sim, 'part_restart.h5', load_atomic_data=.false.)
n_part = size(sim%groups(1)%particles)
allocate(prtkin(n_part))
allocate(iangle(1,n_part), impact_energy(1,n_part), impact_wall(1,n_part))
impact_wall = 0
impact_energy = 0.d0
iangle = 0.d0
previous_wallload = ''
if (restart) then
  previous_wallload = 'wallload_restart.h5'
  call mod_wall_collision_check_restart(sim, previous_wallload)
end if

! Schedule field updates and the end of particle tracking.
! Start with the numbered JOREK restart at or before the particle time.
timesteps(1) = tstep_particles
sim_time = real(nstep_particles,8) * timesteps(1)
events = [event(read_jorek_fields_interp_linear(i=last_file_before_time(sim%time), stop_at_end=.false.)), &
          event(stop_action(),start=sim%time+sim_time)]

! Run first event to read the JOREK fields
call with(sim, events(1))
events(2)%start = sim%time + sim_time

if (nout_particles <= 0 .or. timesteps(1) <= 0.d0) then
  if (sim%my_id == 0) write(*,*) 'nout_particles and tstep_particles must be positive.'
  call wallload_abort(wallload_comm, 1, ierr)
end if

if (.not. (restart .or. restart_particles)) then
  do i=1,size(sim%groups(1)%particles)
     select type (p=>sim%groups(1)%particles(i))
     type is (particle_gc_relativistic)
        call find_RZ(sim%fields%node_list, sim%fields%element_list, &
             p%x(1), p%x(2), &
             p%x(1), p%x(2), p%i_elm, p%st(1), p%st(2), ifail)
     end select
  end do


! Uncomment this block to turn guiding centers to gyro-orbit particles
  !allocate(prtkin(n_part))
  !do i=1,size(sim%groups(1)%particles)
  !   select type (p=>sim%groups(1)%particles(i))
  !   type is (particle_gc_relativistic)
  !      call sim%fields%calc_EBpsiU(sim%time, p%i_elm, p%st, p%x(3), E, B, psi, U)
  !      prtkin(i) = relativistic_gc_to_relativistic_kinetic(sim%fields%node_list,sim%fields%element_list, &
  !           p,sim%groups(1)%mass, B, 2*PI*rnd(1))
  !   end select
  !end do
  !deallocate(sim%groups(1)%particles)
  !allocate(particle_kinetic_relativistic::sim%groups(1)%particles(n_part))
  !do i=1,size(sim%groups(1)%particles)
  !   select type (p=>sim%groups(1)%particles(i))
  !   type is (particle_kinetic_relativistic)
  !      p%x      = prtkin(i)%x
  !      p%p      = prtkin(i)%p
  !      p%st     = prtkin(i)%st
  !      p%i_elm  = prtkin(i)%i_elm
  !      p%weight = prtkin(i)%weight
  !      p%q      = prtkin(i)%q
  !   end select
  !end do
  !deallocate(prtkin)

end if


call update_field_status()
do while (.not. sim%stop_now .and. .not. fields_exhausted .and. sim%time < events(2)%start)
  target_time = next_event_at(sim, events)
  at_fluid_step = events(1)%run_at(target_time)

  do i=1,1
    n_steps = max(1, ceiling((target_time - sim%time)/timesteps(i) - 1.d-12))
    push_dt = (target_time - sim%time)/n_steps

    select type (particles => sim%groups(i)%particles)
    type is (particle_gc_relativistic)
      !$omp parallel do default(private) &
      !$omp shared (i, n_steps, push_dt, sim, wall, iangle, impact_wall, impact_energy)
       do j=1,size(particles,1)
          do k=1,n_steps
             if(particles(j)%i_elm .le. 0) exit

             pos_prev = particles(j)%x

             t_particle = sim%time+real(k,8)*push_dt
             call runge_kutta_fixed_dt_gc_push_jorek(sim%fields,t_particle,push_dt, &
                  sim%groups(i)%mass,particles(j))
             if (particles(j)%i_elm .le. 0) exit
             call mod_wall_collision_check(pos_prev, particles(j)%x, wall, wall_id, wall_pos, iangle(i,j))
             if(wall_id .gt. 0) then
                particles(j)%x      = wall_pos
                particles(j)%i_elm  = -wall_id
                particles(j)%t_loss = t_particle
                impact_wall(i,j) = wall_id
                impact_energy(i,j) = wall_impact_energy(sim, particles(j), sim%groups(i)%mass, t_particle)
             end if
          end do
       end do
       !$omp end parallel do

    type is (particle_kinetic_relativistic)
       !$omp parallel do default(private) &
       !$omp shared (i, n_steps, push_dt, sim, wall, iangle, impact_wall, impact_energy)
       do j=1,size(particles,1)
          do k=1,n_steps
             if (particles(j)%i_elm .le. 0) exit
             pos_prev = particles(j)%x

             t_particle = sim%time+real(k,8)*push_dt
             call volume_preserving_push_jorek(particles(j),sim%fields,sim%groups(i)%mass,t_particle,push_dt,ifail)
             if (particles(j)%i_elm .le. 0) exit
             call mod_wall_collision_check(pos_prev, particles(j)%x, wall, wall_id, wall_pos, iangle(i,j))
             if(wall_id .gt. 0) then
                particles(j)%x      = wall_pos
                particles(j)%i_elm  = -wall_id
                particles(j)%t_loss = t_particle
                impact_wall(i,j) = wall_id
                impact_energy(i,j) = wall_impact_energy(sim, particles(j), sim%groups(i)%mass, t_particle)
             end if
          end do
       end do
       !$omp end parallel do
    end select

  enddo

  sim%time = target_time
  if (at_fluid_step) then
    ! Match event_run: the action pointer takes precedence over stored_action.
    if (associated(events(1)%action)) then
      field_action => events(1)%action
    else
      field_action => events(1)%stored_action
    end if
    select type (reader => field_action)
    type is (read_jorek_fields_interp_linear)
      fluid_step = reader%i
    end select
    call MPI_Bcast(fluid_step, 1, wallload_mpi_integer, 0, wallload_comm, ierr)
    if (mod(fluid_step, nout_particles) == 0) then
      write(step_name, '(i0.5)') fluid_step
      particle_file = 'part_out_fluid'//trim(step_name)//'.h5'
      wallload_file = 'wallload_fluid'//trim(step_name)//'.h5'
      call write_checkpoint(particle_file, wallload_file)
    end if
  end if
  call with(sim, events, at=sim%time)
  call update_field_status()
enddo

call mod_wall_collision_free(wall)

call write_checkpoint('part_out.h5', 'wallload.h5')

deallocate(iangle)

! Finalize the simulation
call sim%finalize

contains

subroutine update_field_status()
  fields_exhausted = .true.
  select type (fields => sim%fields)
  type is (jorek_fields_interp_linear)
    fields_exhausted = fields%static
  end select
end subroutine update_field_status

subroutine write_checkpoint(particle_path, wall_path)
  character(len=*), intent(in) :: particle_path, wall_path

  call mod_wall_collision_export(sim, wall_path, iangle, impact_wall, impact_energy, previous_wallload)
  call write_simulation_hdf5(sim, particle_path)
  previous_wallload = wall_path
  impact_wall = 0
end subroutine write_checkpoint

end program example_wallload
