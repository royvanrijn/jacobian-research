# Clean-room strengthening of the explicit quartic model

This audit independently reconstructs the central exact algebra of
[QUARTIC_WEIGHTED_GEOMETRY.md](QUARTIC_WEIGHTED_GEOMETRY.md).  It imports
neither SymPy nor the project weighted-map implementation.

Starting from

\[
u=1+3xy,\qquad \gamma=1-4xy-x^2z,
\]

the checker performs the exact divisions defining the three map coordinates,
recomputes the full Jacobian, and obtains `det DG=-6`.  It also verifies both
rational collision points and the marked inverse equation

\[
E(W)=W^2-W^4-2BCW+AC^2.
\]

Using the quartic discriminant formula in an independently implemented sparse
polynomial ring gives

\[
\operatorname{disc}_W(E)=-16C^2Q_4,
\]

with the same `Q_4` as the main proof.  The repeated-root parametrization is
checked through

\[
E_{r-2r^3,\,r^2-3r^4}(W)
=-(W-r)^2(W^2+2rW+3r^2-1).                         \tag{1}
\]

Finally, the audit independently checks the `C=0` relation

\[
(B^2-A)x^2=1
\]

on the `gamma=0` chart and the omitted-node factorization

\[
W^2-W^4-\tfrac14=-(W^2-\tfrac12)^2.
\]

Together with the reconstruction and properness arguments in the main note,
these identities supply a second implementation of the algebra controlling
the image and boundary theorem.

Run

```bash
python3 scripts/audit_c14_independent.py
```

The checker uses only `fractions` and custom sparse-polynomial arithmetic.
An unrelated CAS reproduction of the full singular-locus radical and every
formal escape path would be the remaining step before calling the whole C14
package independently verified.
