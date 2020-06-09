Since CMake 3.16 there is support for OpenACC::OpenACC_Fortran target which we like to use.

The contributed file cmake/compat/3.16/FindOpenACC.cmake file is copied from a CMake 3.16 distribution and will be used
instead upon find_package(OpenACC), only when CMake version < 3.16
