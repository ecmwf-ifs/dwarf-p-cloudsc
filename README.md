Dwarf-P-cloudMicrophysics-IFSScheme
-----------------------------------
Contact: Gianmarco Mengaldo (gianmarco.mengaldo@ecmwf.int), 
Sami Saarinen (sami.saarinen@ecmwf.int), 
Willem Deconinck (willem.deconinck@ecmwf.int)

Dwarf-P-cloudMicrophysics-IFSScheme is intended to test the cloud micro physics.

The code is written in Fortran 2003 and it has been tested using the various compilers, including:

    GCC 4.8.
    Cray 8.4.
    PGI.
    INTEL. 

This application does not need MPI nor BLAS libraries for performance. Just a compiler that understands 
OpenMP directives. Fortran must be at least level F2003.

Inside the dwarf directory you can find some example of outputs inside the example-outputs/ directory.

In addition, to run the dwarf it is necessary to use an input file that can be found inside the config-files/ 
directory winthin the dwarf folder.


Prototypes available
--------------------
- prototype1: includes the cloud scheme from IFS that is naturally 
suited to host-type machines and optimized on the Cray system at 
ECMWF.


Third Party Requirements
------------------------
Requirements to compile this dwarf:

    CMake: for use and installation see http://www.cmake.org/
    ecbuild

Recommended:

    OpenMP: shared memory parallelisation


Download and Installation
-------------------------
Please refer to: https://software.ecmwf.int/stash/projects/ESCAPE/repos/escape/browse


Running and testing
-------------------
Please refer to the documentation: https://software.ecmwf.int/wiki/display/ESCAPE/Dwarf+3+-+cloud+scheme 
