# Concave Generalized Bézier Patch #

This is a minimal library for evaluating CGB patches.

## Installation ##

The program was tested under Linux and Windows.

You will need the following libraries:
- libgeom, which is part of the [transfinite](https://bitbucket.org/salvipeter/transfinite) library
- [libharmonic](https://github.com/salvipeter/harmonic)
- Shewchuk's [Triangle](http://www.cs.cmu.edu/%7Equake/triangle.html)

### Building under Linux ###

Download and compile the dependencies.
Note that you have to build Triangle as a library with the `-fpic` flag.

```
> cd concave-gb
> mkdir build
> cd build
> cmake -D LIBGEOM_ROOT=/your/path/to/transfinite \
        -D LIBHARMONIC_ROOT=/your/path/to/harmonic \
		-D LIBTRIANGLE_ROOT=/your/path/to/triangle \
		-D CMAKE_BUILD_TYPE=Debug # or Release \
		../src
> make
```

### Building under Windows ###

Download the dependencies into `dependencies/{transfinite,harmonic,triangle}`.
A Visual Studio 2015 solution file is supplied for compilation.

## Usage ##

The library is documented in the header file `cgb.h`.

There is also a stand-alone program in the `example` directory that can be used from the command line to generate meshes from .cgb files.
Note that a high-density resolution also needs a higher level of detail for the parameterization.
