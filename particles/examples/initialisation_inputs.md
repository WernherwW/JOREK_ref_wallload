# Runtime particle initializer inputs

Rebuild once after updating the source:

```bash
make test_initialisation_gc_phase_space test_initialisation_phase_space
```

Copy the appropriate example namelist into the run directory, edit it, and pass
its path as the first argument. Keep the normal JOREK namelist on standard input:

```bash
mpirun -np "$SLURM_NTASKS" ./test_initialisation_gc_phase_space gc_initialization.nml \
  < input_CQ_particles > output_initREs

mpirun -np "$SLURM_NTASKS" ./test_initialisation_phase_space full_orbit_initialization.nml \
  < input_CQ_particles > output_initREs
```

Without a command-line argument, both binaries use their built-in defaults.
A supplied file overrides only the fields it contains. MPI rank 0 reads it and
broadcasts the resulting configuration, so all ranks receive the same values.
The effective configuration is printed once. Unknown fields and invalid input
stop the MPI job with an error.

The group name is &particle_initialization. The example .nml files list all
supported fields. A job script can generate a small override file directly:

```bash
MARKERS_PER_RANK=1000000
cat > initializer.nml <<EOF
&particle_initialization
  n_particles = $MARKERS_PER_RANK
  start_time = 0.044d0
  Ekinbound = 1.99d6, 2.01d6
/
EOF
mpirun -np "$SLURM_NTASKS" ./test_initialisation_gc_phase_space initializer.nml \
  < input_CQ_particles > output_initREs
```

- n_particles is the marker count **per MPI rank**.
- Energy bounds are in eV, time in seconds, mass in AMU, positions in metres,
  and angles in radians. Namelists require numeric angles, not PI expressions.
- jorek_filename is the restart basename, using the existing reader's naming
  convention, rather than a complete .h5 path.
- test_case supports jorek_recurrent_density_re,
  jorek_current_density_re, uniform_weight, and uniform.
- The RE-current case requires WITH_Refluid, a full toroidal interval, and,
  for full-orbit particles, a full uniform gyro interval. This matches the
  domain of Ipre_tot and the sampling assumption.
- n_variables must remain 6 for GC or 7 for full orbit. nchi and Chibound
  are full-orbit options only. All mesh sizes must be at least 2.
- charge is a nonzero integer in elementary-charge units. Chargebound is
  derived from this fixed charge; internal PDF parameter counts are not inputs.
- do_verification enables histogram/error diagnostics and their output files.
  write_txt controls text output; the particle restart is always written.
- Relative filenames are resolved from the working directory.
- Filename buffers support up to 512 characters. The namelist path supports
  up to 1024 characters.

The existing executable-stack workaround in Makefile.inc remains in effect.
