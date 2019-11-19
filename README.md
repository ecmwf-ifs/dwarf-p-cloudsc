Dwarf-P-cloudMicrophysics-IFSScheme
-----------------------------------
Contact: Michael Lange (michael.lange@ecmwf.int),
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
- **Prototype 1**: The original cloud scheme from IFS that is naturally
  suited to host-type machines and optimized on the Cray system at ECMWF.
- **Prototype 2**: A cleaned up version of the CLOUDSC prototype that is
  validates runs against platform and language agnostic off-line
  reference data via serialbox. The kernel code also is slightly
  cleaner than the one in prototype 1. Note: To run this prototype,
  please install with serialbox enabled (default in bundle) and
  run protoype1 once in the build directory.
- **Prototype 3**: Standalone C version of the kernel that has been generated
  by ECMWF tools. This requires the serialbox validation mechanism as above.

Download and Installation
-------------------------
The preferred method to install the CLOUDSC dwarf uses the bundle
definition shipped in the main repository. For this please
```
git clone ssh://git@git.ecmwf.int/escape/dwarf-p-cloudmicrophysics-ifsscheme.git cloudsc_bundle
cd cloudsc_bundle
<git checkout develop>  # For the latest version please use the `develop` branch
```

The individual protoype variants of the dwarf can be enable or disabled in the `bundle.yml`
file by simply (un)commenting the individual project definitions.

Then simply install the bundle via:
```
./cloudsc-bundle create  # Checks out dependecy packages
./cloudsc-bundle build <--build-type=debug|bit|release> <--toolchain=toolchain_file>
```

Running and testing
-------------------

The different prototype variants of the dwarf create different binaries that all behave similar.
The basic three arguments define (in this order):
* Number of OpenMP threads
* Size of overall working set in columns
* Block size (NPROMA) in columns

An example:
```
cd build
./bin/dwarf-P-cloudMicrophysics-IFSScheme 4 16384 32  # The original
./bin/dwarf-cloudsc-v2 4 16384 32   # The cleaned-up Fortran
./bin/dwarf-cloudsc-v3 4 16384 32   # The standalone C version
```
