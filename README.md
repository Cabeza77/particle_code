# particle_code
Particle code for FARGO-ADSG

This code is for calculating particle orbits in FARGO-ADSG simulations. It might also work for other FARGO versions, but it is not yet tested for them.


COMPILING
--------------------

To compile the code use
    make

If you want to compile the parallel version of the code use
    make BUILD=parallel


EXECUTION
--------------------

To execute use
    ./particle_code input_file.nml
where you have to give the input file as first argument.

If you compiled in parallel mode use
    ./particle_code input_file.nml setthreads N
to run the code with N threads. If you don't give the number of threads, N=1 will be used.
