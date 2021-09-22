# simplex

Simple Implementation of Simplex Method using Tableaus

## Notes
This implementation requires the problem to be in [canonical form](https://en.wikipedia.org/wiki/Linear_programming).

That is, you want to maximize/minimize c^T x subject to the constraints Ax <= b and x >= 0.

For now, this implementation requires an initial feasible solution at the origin.

For now, this implementation avoids degeneracy (by using [Bland's rule](https://en.wikipedia.org/wiki/Bland%27s_rule)). I might later modify the implementation to choose the entering variable to be the one with the largest coefficient, then resort to Bland's rule after multiple iterations of degeneracy.

## Requirements
This implementation uses Python and requires [numpy](https://numpy.org/install/).

## Usage
Clone the repository and change into it.
```
$ git clone https://github.com/jli0108/simplex.git
$ cd simplex
```
Modify `maximize` depending on whether you want to solve a maximization or minimiation problem.

Modify the arrays `c`, `A_B`, and `b` in the `simplex.py` file with your favorite editor, i.e.
```
$ code simplex.py
```
Run `simplex.py` with either `python simplex.py` or `python3 simplex.py`.
