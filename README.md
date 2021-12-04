# simplex

Implementation of Simplex Method in Python

## Notes
This implementation can be used for [canonical or standard form](https://en.wikipedia.org/wiki/Linear_programming).

For now, this implementation avoids degeneracy (by using [Bland's rule](https://en.wikipedia.org/wiki/Bland%27s_rule)). I might later modify the implementation to choose the entering variable to be the one with the largest coefficient, then resort to Bland's rule after multiple iterations of degeneracy.

The revised simplex method avoids tableaus but is a bit less readable.

## Requirements
This implementation uses Python and requires [numpy](https://numpy.org/install/).

## Usage
Clone the repository and change into it.
```
$ git clone https://github.com/jli0108/simplex.git
$ cd simplex
```
Modify either `simplex.py` or `revised_simplex.py`.

Modify `maximize` depending on whether you want to solve a maximization or minimiation problem.

Modify the arrays `c`, `A_B`, and `b` with your favorite editor, i.e.
```
$ code simplex.py
```
Run `simplex.py` or `revised_simplex.py` with either `python` or `python3`.
