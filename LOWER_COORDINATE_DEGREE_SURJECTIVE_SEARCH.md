# Lower-coordinate-degree surjective search

## Result of the first exact pass

No degree below \(17\) is produced in the tested classes.  The search does
find a sharp false lead: the cancellation type \((m,r)=(1,2)\) has coordinate
degrees

\[
(\deg F_1,\deg F_2,\deg F_3)=(9,8,13),
\]

but it is not surjective.  Its inverse polynomial at
\((P,Q,R)=(1,2,1/5)\) is

\[
\Psi(T)-R=\frac{(T-1)^5}{5},\qquad
D(T)=1-T(2-T)=(T-1)^2.
\]

The only inverse root is therefore entirely on the reconstruction boundary.
This rules out the apparent degree-\(13\) record by the exact full-contact
criterion, rather than by a numerical image test.

The strongest surviving candidate from another family is the cancellation
type \((m,r)=(3,1)\).  Its coordinate degrees are

\[
(13,10,24).
\]

It is surjective: comparison of its inverse antiderivative with
\(\lambda(T-a)^3(T-b)^2\) gives the unit ideal.  Thus no degree-five inverse
polynomial in this family can lose all of its simple roots.  The normalized
formula is too large to improve the current record.

## Weighted quintics

For the canonical seed

\[
H(W)=W^4(1-W),
\]

the weighted map has degree profile \((17,16,4)\) and determinant one.  The
full-contact test is particularly small.  A monic quintic without a simple
root has the form

\[
(W-a)^3(W-b)^2,
\]

including the quintuple-root case \(a=b\).  Matching the unaffected
\(W^4,W^3,W^2\) coefficients of

\[
W^5-W^4+sW-t
\]

gives an ideal equal to \((1)\).  Hence the map is surjective.

The target \((-16,-16,1)\) has the two distinct rational preimages

\[
\left(-\frac1{15},16,-3880\right),\qquad
\left(\frac1{32},-30,30624\right).
\]

Thus this simple canonical member already supplies determinant,
noninjectivity, and surjectivity at the standing bound.

The leading terms of a general genuine weighted quintic are forced by its
nonzero quintic coefficient, so coefficient sparsity does not lower
\((17,16,4)\).  The checker additionally exhausts every single elementary
source shear

\[
v\longmapsto v+\lambda m(v_1,v_2),
\qquad 1\le\deg m\le2,
\]

on the canonical member and finds no specialization with maximum coordinate
degree below \(17\).

## Cancellation quintics

The degree equation

\[
r(m+1)+1=5
\]

has only the two types

\[
(m,r)=(1,2),\qquad(3,1).
\]

The first expands to \((9,8,13)\) but has the quintuple-contact omission
displayed above.  The second expands to \((13,10,24)\) and passes the exact
full-contact test.  Every single elementary source shear of degree at most
two was also checked on the surjective \((3,1)\) map; none lowers its maximum
degree below \(24\), hence none challenges \(17\).

Low-degree target shears do not create a missing cancellation in the
normalized \((3,1)\) profile: the degrees of \(P,Q\) are \(13,10\), and the
degrees of their monomials of target degree at most three are all distinct.
In particular no such monomial has degree \(24\), the leading degree of
\(R\).

## Sparse quadratic gauges

The sparsest genuine quintic seed is

\[
G(S)=S+\mu S^3+\nu S^5,\qquad \mu\nu\ne0.
\]

Its normalized quadratic-gauge map still has coordinate degrees

\[
(7,32,30).
\]

More decisively, every member of this sparse two-torus is nonsurjective.
Choose \(p,a,b\) over the algebraic closure satisfying

\[
15\nu p^3=4\mu^2,\qquad
\mu p a^2=-1,\qquad
2b+3a=0.
\]

After choosing the two free target coefficients \(B,C\), direct coefficient
comparison gives

\[
E_{p,B,C}(S)=\nu p^5(S-a)^3(S-b)^2.
\]

This is an exact full-contact omitted target.  Sparse coefficients therefore
do not create a competitive surjective quadratic gauge.

## Alternative primitive generators

The existing bounded primitive-root search was rerun.  In the degree-two
mixed-root ansatz, primitiveness reduces the candidate to a Möbius function
of the old root.  The remaining strata fail as follows:

- a root-only non-affine Möbius change has an unavoidable order-two pole;
- a \(1/P\) translation makes the controlled-divisor coefficient
  nonpolynomial;
- the reciprocal scaling \(T=S/P\) leaves the forced term \(QS^2/P\);
- a unimodular nonmonomial denominator requires a nonpolynomial
  third-order shear.

Thus the bounded one-primitive search does not yield a fourth polynomial
suspension.  Determinant-supported multi-prime ledgers and a second primitive
reconstruction variable remain genuinely open; this pass does not claim to
exclude them.

## Reproduction

Run

```bash
.venv/bin/python scripts/search_lower_coordinate_degree_surjective.py
.venv/bin/python scripts/explore_root_changing_incidence_degree_four.py
.venv/bin/python scripts/verify_fourth_suspension_valuation_fan.py
.venv/bin/python scripts/verify_fourth_suspension_primitive_root.py
```

The first checker verifies the three degree-five family comparisons,
determinants where competitive, exact full-contact decisions, the rational
weighted collision, the sparse quadratic-gauge omission, and the bounded
elementary source-shear search.
