# Archived two-transfer example

This is the one exact local-algebra example worth extracting from the former
transfer programme.

For

\[
 S=Z^2+pZ+q,qquad V=S^2+XZ+Y,
\]

triangular elimination from `U^2=V^3` gives

\[
 X^3=0,qquad 2XY=pX^2,qquad Y^2=qX^2.                    \tag{1}
\]

Hence the resulting algebra is finite flat of rank four over `k[p,q]`, with
basis

\[
 1,quad X,quad Y,quad X^2.                              \tag{2}
\]

At `p=q=0` its fiber is

\[
 k[X,Y]/(X^3,XY,Y^2),                                     \tag{3}
\]

of Hilbert function `(1,2,1)` and socle dimension two.  This corrects the
naive tensor product of two dual-number factors, but it does not extend to an
all-`k` theorem.

The original derivation, together with the isolated `k=3,4` computations, is
preserved in
[three-four-transfer-computations.md](three-four-transfer-computations.md).
The failure of both proposed uniform models is documented in
[failed-all-k-filtrations.md](failed-all-k-filtrations.md).
