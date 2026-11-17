# Exact bounded sparse exclusion

This supplement enumerates all supports of degree at most four with at most four
monomials in Gaussian Witt coordinates. It uses exact Wick coefficients and
rational Gröbner bases in Singular. The three all-order family arguments are in
Section 4 of the article; a finite replay does not replace those arguments.

Run from this directory with Python 3.12, SymPy 1.14.0 and Singular 4.3.2:

~~~sh
python3 verify_three_real_gmc_sparse_pareto.py
~~~

The two Python files are copied without mathematical changes from the retained
research checker. The legacy script and output use “Pareto” as a name; their
certified exclusion is exactly the bounded rectangle stated above, not global
Pareto minimality of degree/support (4,5).

The program enumerates supports rather than trusting a stored survivor list,
checks coefficient-symmetry normalization on every support, and verifies the
three surviving coefficient conditions and product factorization. A nonzero
exit code or incomplete run is a failed or incomplete replay. Independent
implementation and external mathematical review are separate obligations.
