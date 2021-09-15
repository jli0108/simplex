# simplex

Simple Implementation of Simplex Method using Tableaus

## Notes
This implementation assumes you want to solve a maximization problem in canonical form.

That is, you want to maximize c^T x subject to the constraints Ax <= b and x >= 0.

I think the implementation avoids degeneracy (because of [Bland's rule](https://en.wikipedia.org/wiki/Bland%27s_rule)) , but I have not thoroughly tested it.

## Requirements
This implementation uses Python and requires [numpy](https://numpy.org/install/).

## Usage
Clone the repository and change into it
```
$ git clone https://github.com/jli0108/simplex.git
$ cd simplex
```
Modify the arrays `c`, `A_B`, and `b` in the `simplex.py` file with your favorite editor, i.e.
```
$ code simplex.py
```
Run `simplex.py` with either `python simplex.py` or `python3 simplex.py`.
