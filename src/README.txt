This is the CLOUDSC mini-application (version 1.0)
--------------------------------------------------

A test case has NGPTOT = 160000

To run a test case with varying NPROMAs = 1000 64 12
use the following command :

./run COMP 160000 1000 64 12 2>&1 | tee COMP.$(hostname -s).out

where COMP is one of INTEL, GNU, PGI or CCE

This command will rebuild from scratch using the 'makefile'
and run test cases with varying number of OpenMP-threads

Once your executable is around, you could
also run the executable (caveat: stack size not set, nor binding/affinity)
manually by

./main.x.COMP OMP NGPTOT NPROMA-list

where OMP is the number of threads for OMP-parallel regions,
NGPTOT the number of grid point columns, NPROMA-list
is a list of NPROMAs to use.  

Some example outputs can be seen under example_outputs/ directory.

Please observe also the directory bin/ for some tools to discover 
f.ex. your system properties (core count, number of sockets etc.) 
These tools are invoked by the run-script

The input data cloudsc.bin is a Fortran unformatted stream binary
(no record delimiters). In contains data for just 100 grid point columns
and will be inflated to full spectre of NGPTOT.

This application does not need MPI nor BLAS libraries for performance.
Just a compiler that understands OpenMP directives.
Fortran must be at least level F2003. 


Sami Saarinen, ECMWF
20-Oct-2015
email: sami.saarinen@ecmwf.int
