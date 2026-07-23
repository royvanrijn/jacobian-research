# Puncture-rank frontier for Keller boundaries

This note audits the proposed assertion that a uniquely selected rational
ramified boundary of a nonproper Keller map has unit rank at most one.  It
separates a formulation issue from a new low-degree obstruction.

The outcome is:

1. the literal statement with `E=A^1` or `G_m` is not a statement about the
   canonical boundary prime used in this repository: that prime is a
   surface;
2. the natural corrected examples are
   `A^2` and `G_m x A^1`, whose unit ranks are zero and one;
3. finite normalization, tame generic ramification, rationality, and absence
   of a second boundary prime over the same target divisor do not by
   themselves turn Jelonek's uniruledness theorem into a puncture bound;
4. nevertheless, the most direct two-unit generalization of the reciprocal
   chart is impossible in every geometric degree from four through seven.

The last item is an exact bounded theorem, not heuristic search.

Work over an algebraically closed field `k` of characteristic zero.

## 1. The dimensional correction

For a dominant quasi-finite map

\[
 F:\mathbb A^3\longrightarrow\mathbb A^3
\]

the canonical finite normalization

\[
 \bar X_F=\operatorname{Norm}_{\mathbb A^3}k(\mathbb A^3)
\]

is three-dimensional.  Every prime component of
`\bar X_F - A^3` has codimension one, hence dimension two.  Thus a canonical
boundary prime `E` cannot literally be isomorphic to `A^1` or `G_m`.

In the two established suspension families, the relevant normalizations
have the form

\[
 E_{\rm wt}\simeq\mathbb A^1\times\mathbb A^1,\qquad
 E_{\rm rec}\simeq\mathbb G_m\times\mathbb A^1,
\]

or have the same unit lattices after the intrinsic affine-line parameter is
retained.  Consequently

\[
 \operatorname{rank}\mathcal O(E_{\rm wt})^*/k^*=0,
 \qquad
 \operatorname{rank}\mathcal O(E_{\rm rec})^*/k^*=1.
\]

A corrected puncture statement must therefore do one of the following:

- state only the unit-rank inequality for the boundary surface;
- canonically extract an `A^1`-fibration `E -> C` and state the conclusion
  for the relative curve `C`; or
- assume from the start that `E=C x A^1` and conclude
  `C=A^1` or `G_m`.

The last conclusion also needs normality (equivalently smoothness for a
curve here).  A geometrically integral rational affine curve of unit rank
zero need not itself be `A^1`; singular rational curves give immediate
counterexamples to that inference.

## 2. Why the suggested general inputs do not yet prove the bound

Two of the proposed hypotheses carry no extra force in the present setting.
For a generically finite polynomial self-map of affine space, normalization
in the finite function-field extension is finite because finitely generated
algebras over a field are excellent (in particular Nagata).  Generic DVR
ramification is tame automatically in characteristic zero.  The substantive
hypotheses are therefore the intrinsic selection, rationality, and the
no-additional-boundary condition; the last must specify whether it excludes
only boundary primes over the same target divisor or also transverse affine
divisors meeting the selected boundary.

Jelonek--Lasoń prove that the nonproperness hypersurface of a generically
finite polynomial map from affine space is covered by polynomial curves of
bounded degree.  For a curve, being covered by polynomial curves would force
a one-puncture normalization.  Here, however, the nonproperness set is a
surface and the selected canonical boundary is another surface finite over
it.  Uniruledness of the target surface does not bound the unit rank of a
curve factor in the source surface.

The generic DVR identity

\[
 \sum_{E_i\mid Z}e(E_i/Z)f(E_i/Z)=\operatorname{gdeg}(F)
\]

does not repair this gap.  The residue degree `f(E/Z)` is the generic degree
of one surface over another.  It is not the degree of a chosen punctured
curve over `A^1`.  In particular, a three-punctured relative curve need not
consume three residue sheets.

Likewise, tame log crepancy

\[
 K_{\bar X_F}+D_F=\pi^*(K_{\mathbb A^3}+B_F)
\]

is a codimension-one identity.  After adjunction it gives the expected
Riemann--Hurwitz identity on any genuinely extracted curve quotient, but it
does not construct that quotient or bound its punctures.

The missing general-theorem input is therefore precise:

> **Relative-curve extraction problem.**  Canonically extract from the
> selected boundary surface a curve quotient `E -> C` such that independent
> units of `E` descend to puncture units of `C`, and prove that the Keller
> reconstruction gives at most one independent puncture character.

Without this input, the proposed theorem is stronger than the listed tools.

## 3. Universal two-center reciprocal chart

The counterexample strategy in the question has a clean universal form.
Let `f in k[y]` have degree `d>=1`, let `r>=1`, and put

\[
 A=1+xf(y),\qquad
 B=A^{r+1}z+g(y,A),
\]

\[
 s=\frac{x}{A},\qquad P=AB,\qquad Q=y+xB,\qquad
 Y=Q-sP.
\]

Then identically

\[
 Y=y,\qquad
 D=1-sf(Y)=A^{-1}.
\]

The normal critical boundary of the coordinate-preserving core has ring

\[
 \mathcal O(E)
 =k[Y,P,f(Y)^{-1}].
\]

If `f` has `ell` distinct roots, then

\[
 \operatorname{rank}\mathcal O(E)^*/k^*=\ell.
\]

Thus two distinct roots give exactly the desired rank-two surface

\[
 \left(\mathbb P^1\setminus\{\alpha,\beta,\infty\}\right)
 \times\mathbb A^1.
\]

The chart determinant is independent of `f` and `g`.  First,

\[
 \det\frac{\partial(s,P,Q)}{\partial(x,y,B)}=-A^{-1}.
\]

Since `partial_z B=A^{r+1}`,

\[
 \boxed{
 \det\frac{\partial(s,P,Q)}{\partial(x,y,z)}=-A^r=-D^{-r}.
 }
\tag{3.1}
\]

Define the core output

\[
 R=C\int_0^s\{1-tf(Q-Pt)\}^r\,dt.
\tag{3.2}
\]

Then `R_s=CD^r`, so (3.1) gives

\[
 \det\frac{\partial(P,Q,R)}{\partial(x,y,z)}=-C
\]

whenever (3.2) is polynomial.  Its generic inverse degree is

\[
 N=r(d+1)+1.
\tag{3.3}
\]

This is the smallest direct realization of
`D=1-sU_1^aU_2^b`: after affine normalization of the two roots,

\[
 f(y)=y^a(y-1)^b,\qquad a,b\ge1.
\]

The determinant ledger therefore poses no obstruction.  Polynomiality is
the entire issue.

## 4. The boundary-moment obstruction

Write

\[
 g_0(y)=g(y,0).
\]

In (3.2), substitute `t=sv`.  The identities above give

\[
 R=
 \frac{Cx}{A^{r+1}}
 \int_0^1
 \left[
 A-xv f\bigl(y+xB(1-v)\bigr)
 \right]^r\,dv.
\tag{4.1}
\]

If `R` is polynomial, the numerator in (4.1) is divisible by `A^{r+1}`.
Reducing its constant term modulo `A`, where

\[
 x=-f(y)^{-1},\qquad B=g_0(y),
\]

gives the necessary identity

\[
 \boxed{
 \mathcal M_{f,r}(y,g_0):=
 \int_0^1v^r
 f\left(
 y-\frac{g_0(y)}{f(y)}(1-v)
 \right)^r\,dv=0.
 }
\tag{4.2}
\]

After clearing denominators, (4.2) is a polynomial equation in
`k[y,g_0]`.  It is only the first cancellation equation, so failure of
(4.2) rules out every higher jet `g(y,A)` at once.

### Theorem 4.1 -- no two-center reciprocal lift in degrees four through seven

Let

\[
 f(y)=y^a(y-1)^b,\qquad a,b\ge1,
\]

and form (3.1)--(3.2).  If

\[
 4\le N=r(a+b+1)+1\le7,
\]

then there is no polynomial `g(y,A)` for which `R` is polynomial.
Consequently this entire two-center reciprocal family contains no
rank-two Keller map of geometric degree four, five, six, or seven.

#### Proof

The possible triples `(r,a,b)` are:

\[
\begin{array}{c|c}
r& (a,b)\\ \hline
1&(1,1),\\
1&(1,2),(2,1),\\
1&(1,3),(2,2),(3,1),\\
1&(1,4),(2,3),(3,2),(4,1),\\
2&(1,1).
\end{array}
\tag{4.3}
\]

For each row, clear the denominator in (4.2), replace `g_0(y)` by an
indeterminate `G`, and obtain

\[
 H_{r,a,b}(y,G)\in\mathbb Q[y,G].
\]

Exact absolute factorization shows that every polynomial in (4.3) is
absolutely irreducible and has `G`-degree greater than one.  If a polynomial solution
`G=g_0(y)` existed, then `G-g_0(y)` would divide `H_{r,a,b}` in
`k[y,G]`, a contradiction.  The checker listed in Section 7 performs these
factorizations with Singular's characteristic-zero absolute-factorization
algorithm and verifies the complete list (4.3).  It separately factors over
`Q` as a regression against the generated inputs.

For the first and smallest case no computer algebra is needed.  With
`r=a=b=1`, equation (4.2) is equivalent to

\[
 G^2-2y(y-1)(2y-1)G+6y^3(y-1)^3=0.
\tag{4.4}
\]

Its discriminant is

\[
 4y^2(y-1)^2(1+2y-2y^2),
\]

which is not a square in `k[y]`.  Hence (4.4) has no polynomial root.
QED

### Scope

The theorem allows:

- arbitrary positive multiplicities at the two boundary centers;
- every geometric degree from four through seven in this chart;
- an arbitrary polynomial center `g(y,A)`, not merely a monomial or a
  bounded coefficient ansatz; and
- arbitrary higher `A`-jets and arbitrary `z`-terms, because the obstruction
  occurs already modulo `A`.

It does not exclude:

- a nonmultiplicative rational source chart;
- two independent reconstruction variables;
- a nontrivial target Jacobian ledger;
- a critical boundary not presented by `1-sf(Q-Ps)`; or
- a construction in which the two puncture units do not lift from two
  affine-linear boundary centers.

Those are exactly the escape routes already suggested by the
controlled-boundary programme.

## 5. What this says about the conjecture

The proposed rank bound should be split into two assertions.

### Corrected structural conjecture

Let `E` be the canonically selected normal ramified boundary **surface**.
Assume the canonical package also extracts an `A^1`-fibration

\[
 E\longrightarrow C
\]

onto a smooth rational affine curve, identifies
`\mathcal O(E)^*/k^*` with `\mathcal O(C)^*/k^*`, and supplies a single
primitive reconstruction character.  Under the saturated-link,
boundary-monotonicity, no-extra-boundary, and closed-point conductor
conditions, conjecture that

\[
 \operatorname{rank}\mathcal O(E)^*/k^*\le1.
\]

Then `C=A^1` or `G_m`; it is `C`, not the surface `E`, that has this
isomorphism type.

### Proven low-degree theorem

For the direct reciprocal realization with two affine-linear centers,
the conjectural conclusion holds in degrees at most seven for the stronger
reason that polynomiality already fails at its first boundary moment.

This changes the smallest computational target.  Searching more coefficients
inside the same degree-four or degree-five reciprocal chart cannot succeed.
A meaningful counterexample search must leave that chart, most economically
by introducing a target ledger or a second reconstruction variable.

## 6. Next proof target

Equation (4.2) exists in every degree.  For

\[
 f=y^a(y-1)^b
\]

the bounded factorizations suggest:

> **Two-center moment conjecture.**  For all `a,b,r>=1`, the cleared
> polynomial `H_{r,a,b}(y,G)` has no factor linear in `G`.

This statement would exclude the entire direct two-unit reciprocal family in
all degrees.  It is much narrower than the general puncture-rank conjecture
and is now a concrete problem about one family of beta-type moments.

A proof can start with two elementary valuation reductions.  Any rational
solution `G=g_0(y)` forces `f|g_0`; otherwise (4.2) has a nonzero leading pole
at one of the two roots.  Writing `g_0=fc`, comparison at infinity forces
`deg c<=1`.  The remaining problem is therefore finite-dimensional:
exclude an affine `c(y)=lambda y+mu` from the moment identity.  The leading
coefficient is a Jacobi/beta moment; the lower coefficients encode the
incompatibility of the two distinct centers.  This is the most promising
route to an all-degree theorem.

## 7. Exact verification

Run

```bash
.venv/bin/python scripts/verify_puncture_rank_frontier.py
```

The checker verifies:

1. the universal chart determinant;
2. the rank-two boundary ring for two distinct roots;
3. the degree list (4.3);
4. the hand quadratic (4.4) and its nonsquare discriminant; and
5. absolute irreducibility and absence of a linear `G`-factor for every
   degree-four through degree-seven boundary moment.
