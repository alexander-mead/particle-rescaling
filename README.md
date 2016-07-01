# Particle rescaling

This code carries out the Angulo & White (2010; http://arxiv.org/abs/0912.4277) rescaling algorithm on n-body dark matter particles. It uses the FFTW libraries, which need to be linked with the compile command. I use

gfortran modpart.f90 -L/usr/local/lib -lfftw3

and I think the code should compile with any Fortran compiler, but you may need to edit the FFTW library location in the first few lines of modpart.f90 to be appropriate for your machine.

The code needs to be run with an 'ini' parameters file. Two examples are given: 'example_zs.ini', which only carries out the redshift and size scaling parts of the method, and 'example_zsd.ini', which also does the displacement field correction.

The stuff in the ini files should be fairly obvious, but it needs to be pointed towards files that contain tabulated growth functions and growth rates, as well as tabulated matter power spectra.

growth_initial - tablulated growth function in the initial cosmology (a vs. g(a) in ascending a order)

growth_target  - tablulated growth function in the target cosmology (a vs. g(a) in ascending a order)

rate_initial   - tablulated logarithmic growth rate in the initial cosmology (a vs. f(a) in ascending a order)

rate_target    - tablulated logarithmic growth rate in the target cosmology (a vs. f(a) in ascending a order)

pk_initial     - tablulate k vs. P(k) from CAMB

pk_target      - tablulate k vs. P(k) from CAMB

sig8_initial   - initial cosmology sigma_8

sig8_target    - target cosmology sigma_8

L_target       - target box size

Om_m target    - target Omega_m

Om_v target    - target Omega_v

izsd           - 0/1 : Decide to do displacement field correction

ifudge         - 0/1 : Decide to fudge the displacement field to have the 'correct' variance (not sure if this is a good idea)

idisp          - 0/1 : If 1 then the displacement field can be read in directly, otherwise it is 'reconstructed' from the particle data

n_files        - number of files to rescale

file1          - First gadget-2 file location, output location, z for output, k_nl

file2          - etc.
