# Lifting Line Theory (C++)
By Benjamin Moore

## Aim of the code:
- provide the base lifting line theory code but written in and for c++
- produce the data required to design a wing for UoN MMME3080 Coursework 1
- display the results in an understandable manner

## Dependencies
This code relies on:
- matplot++ - https://alandefreitas.github.io/matplotplusplus/
- gnuplot (required by matplot++) - https://www.gnuplot.info/
- eigen3 - https://eigen.tuxfamily.org/index.php?title=Main_Page

Vcpkg has been used to set these up, however they can be installed independently 
- vcpkg - https://vcpkg.io/en/

## Compiling
The code has been designed to be compiled with CMake. \
Once the dependencies have been installed the entire program should compile

The program has been set up using vcpkg and it is recommended that this method is used for both installing the dependencies and helping CMake link to the dependencies (Especially for windows machines).

If a build other than *Visual Studio 2022* is required, you can change the cmake generator in `CMakePresets.json`. This can be set to any of the cmake generators depending on the build desired.

To compile using vcpkg (This uses the cmake generator that has been chosen):\
&nbsp;&nbsp;&nbsp;&nbsp;`cmake --preset=default`\
&nbsp;&nbsp;&nbsp;&nbsp;`mkdir build`\
&nbsp;&nbsp;&nbsp;&nbsp;`cmake --build build`
