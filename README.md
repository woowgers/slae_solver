![Tests](https://github.com/woowgers/slae_solver/actions/workflows/c-cpp.yml/badge.svg)

# Simple SLAE solver
## Solver uses Gauss' method to solve the system.
## SLAE format
```
<number of equations> <number of variables>
<coefficients of row #1>
<coefficients of row #2>
...
<coefficients of row #N>
```
You can either specify a system via `stdin`, or pass a file with a system as a program parameter.


## TODO
* Add rational fractions support
* Add more tests
* Add CI
