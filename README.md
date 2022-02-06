# System of linear equations solver

A Python module that can solve systems of linear equations using Gauss-Jordan elimination.

*Developed for fun and as a coding exercise; better alternatives certainly exist.*

**Tested with Python 3.10.1, untested with other versions.**

## Usage:

```python
# (Assuming the file is named "solver.py")
import solver

# Solve the system {x-2y=4, -2x+3y=0}
solution = solver.solve([ [1, -2, 4], [-2, 3, 0] ])
# solution = ('-12.0', '-8.0') ; x = -12, y = -8

# Solve the system {y+4z=2}
solution = solver.solve([ [0, 1, 4, 2] ])
# solution = ('X0', '-4.0*X2+2.0', 'X2') ; x = x, y = -4z+2, z = z
# Variables in a solution are named as 'Xn', where n is an integer
# starting at 0

# Example of an impossible system {x=1, x=2}
solution = solver.solve([ [1, 1], [1, 2] ])
# solution = ()
```


