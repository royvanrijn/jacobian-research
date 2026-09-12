# The one-pair Schur obstruction for the Dvorsky cubic

## 1. Statement

Work over a characteristic-zero field.  Let

\[
 P=(t+c)(ad+bt)
\tag{1.1}
\]

be the Dvorsky--Long cubic, and add one variable \(s\).  The most direct
quadraticization pairs the previously cubic operator factor
\(\partial_t\) with \(\partial_s\):

\[
 \widetilde\Delta
 =\partial_a\partial_d-\partial_b\partial_c+\partial_t\partial_s.
\tag{1.2}
\]

The symbol of (1.2) is a nondegenerate quadratic form in six variables.
Thus it is linearly equivalent over an algebraic closure to an ordinary
Laplacian.  Nevertheless every polynomial or formal hyperplane lift for
this canonical one-pair completion is impossible.

> **Theorem 1.1 — complete hyperplane-lift obstruction.** Let
> \(F(a,b,c,d,t,s)\) be any polynomial, or formal power series in \(s\)
> with polynomial coefficients, such that
> \[
>  F(a,b,c,d,t,0)=P.
> \tag{1.3}
> \]
> If \(\widetilde\Delta F=0\), write
> \[
>  \left.\partial_sF\right|_{s=0}=-ct+R_0(a,b,c,d)
> \tag{1.4}
> \]
> and put \(\rho=(\partial_cR_0)(0)\).  Then
> \[
>  \left.
>  \widetilde\Delta^2(F^2)
>  \right|_{a=b=c=d=s=0}
>  =12t^2-8\rho t.
> \tag{1.5}
> \]
> In particular, the second pure moment is nonzero, so no such lift satisfies
> \(\widetilde\Delta^m(F^m)=0\) even for \(m=1,2\).

This is an obstruction to one precise but natural quadraticization class.
It is not an obstruction to a different nondegenerate quadratic operator,
a lift with more auxiliary variables, or a nonlinear specialization that
does not recover \(P\) on the hyperplane \(s=0\).

## 2. The unrestricted transverse jet

Expand an arbitrary lift along the new direction:

\[
 F=P+sR_1+s^2R_2+O(s^3).
\tag{2.1}
\]

\[
 D=\partial_a\partial_d-\partial_b\partial_c.
\]

Since

\[
 DP=c,
\tag{2.2}
\]

The coefficients of \(s^0\) and \(s^1\) in
\(\widetilde\Delta F=0\) give

\[
\partial_tR_1=-c,\qquad
DR_1+2\partial_tR_2=0.
\]

Consequently

\[
 R_1=-ct+R_0(a,b,c,d)
\tag{2.3}
\]

and
\[
 R_2=-\frac{t}{2}DR_0+S_0(a,b,c,d)
\tag{2.4}
\]

for arbitrary polynomials \(R_0,S_0\).  Higher coefficients of \(F\) are
constrained recursively, but they cannot enter the restriction of
\widetilde\Delta^2(F^2)\) to \(s=0\), since that operator differentiates
at most twice with respect to \(s\).

## 3. The second-moment obstruction

Only the two-jet of \(R_0,S_0\) at the origin can affect the required
axis restriction.  Substituting (2.3)--(2.4) and applying the product rule
gives

\[
 \left.
 \widetilde\Delta^2(F^2)
 \right|_{a=b=c=d=s=0}
 =4t(3t-2\rho),
\tag{3.1}
\]

where \(\rho=(\partial_cR_0)(0)\).  Every other transverse jet cancels.
The coefficient of \(t^2\) is therefore always \(12\), proving (1.5).

This identifies the first structural gate in the proposed Schur route:
merely completing the determinant block by pairing \(t\) with one new
direction cannot work while preserving the cubic and the restriction
(1.3), regardless of degree mixing.  The next viable ansatz must change
at least one of those requirements:

1. add a second auxiliary block and use its Schur rank to cancel the
   invariant \(12t^2\); or
2. replace hyperplane restriction by a nonlinear polarization that
   changes the transverse jet identity; or
3. use a different nondegenerate quadratic completion whose cross terms
   are not equivalent to (1.2) under a specialization-preserving change.

Any replacement should still be designed to preserve the all-order pure
and mixed contractions.

Any candidate must still preserve the fixed-multiplier defect; solving only
the pure moment equations does not produce an ordinary-Laplacian GVC
counterexample.

## 4. Reproduction

Run

```bash
.venv/bin/python scripts/verify_dvorsky_one_pair_schur_obstruction.py
```

The checker retains the homogeneous cubic normal-form regression, then
parametrizes the complete two-jets of the arbitrary polynomials \(R_0\)
and \(S_0\).  It verifies the first two transverse harmonic equations and
extracts the exact axis formula (3.1).
