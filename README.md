# A vertex-centric representation for adaptive diamond-kite meshes

This is the code accompanying the paper "A vertex-centric representation for adaptive diamond-kite meshes" submitted for publication.

Adaptive diamond-kite meshes were introduced in the paper [Diamond-kite adaptive quadrilateral meshing](https://doi.org/10.1007/s00366-013-0327-9) by Eppstein (2014).

## Running the code

Do `make` and it outputs a sample adaptive diamond-kite mesh in four formats: eps, csv, off, obj. The eps file is converted to pdf and displayed.

## Dependencies

My main platform is macOS.
The code is written in Python 2 but the only change for Python 3 is fixing `print` statements. 
The Makefile uses [epstopdf](https://tug.org/epstopdf/) from [MacTeX](https://tug.org/mactex/) but the builtin `/usr/bin/pstopdf ` works as well.


