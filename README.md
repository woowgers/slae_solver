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
* Figure out why tests are failing on stepik.org
* Add rational fractions support
