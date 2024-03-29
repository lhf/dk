# A vertex-centric representation for adaptive diamond-kite meshes

This is the code accompanying the paper [A vertex-centric representation for adaptive diamond-kite meshes](https://doi.org/10.1016/j.cag.2024.103910) (2024).

Adaptive diamond-kite meshes were introduced in the paper [Diamond-kite adaptive quadrilateral meshing](https://doi.org/10.1007/s00366-013-0327-9) by Eppstein (2014).

Running `make` outputs a sample adaptive diamond-kite mesh in four formats: eps, csv, off, obj. The eps file is converted to pdf using pstopdf.

The file `in.csv` is a sample adaptive diamond-kite mesh generated by random refinement. It can be loaded with `make I=in.csv`.

The program `kf.py` is a variant of `dk.py` that computes kite fractals described in the paper [Fractal tilings based on kite- and dart-shaped prototiles](https://doi.org/10.1016/S0097-8493%2800%2900134-5) by Fathauer (2001). Running `make kf` outputs a sample kite fractal.

The code should work in both Python 2 and Python 3.

