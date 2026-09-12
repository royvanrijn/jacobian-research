# Coprime carrier jets: saturation, discriminant, and the first missing invariant

> **Status.**  Sections 1, 2, and 4 are general exact identities or parameter
> counts.  Section 3 is an exact symbolic audit in bidegrees `(2,3)`, `(2,5)`,
> `(3,4)`, and `(3,5)`; the displayed all-bidegree determinant formula is a
> conjecture.  Sections 5 and 6 are exact for the F2 `k=1` bidegree `(3,5)`.
> They show that the seven-center carrier packet ends exactly at raw-parameter
> saturation, while its nonimmersion divisor generically replaces one node by
> one cusp and changes the conductor scheme from two reduced points to one
> double point.  No identification of this jet ramification with a localized
> `ch_2`/`Fitt_1` class is claimed.

The exact calculations are replayed by
[`verify_coprime_carrier_jet_discriminant_pattern.py`](../scripts/verify_coprime_carrier_jet_discriminant_pattern.py).

## 1. Universal primitive carrier coordinates

Let `1<m<n` be coprime.  Choose the unique positive integers `r,s` with

\[
 sn-rm=1,\qquad 1\le s<m.                       \tag{1.1}
\]

For a polynomial parametrization with leading terms

\[
 p(t)=t^m+\cdots,\qquad q(t)=t^n+\cdots,
\]

put

\[
 x=\frac{p^r}{q^s},\qquad y=\frac{p^n}{q^m}.   \tag{1.2}
\]

At `t=infinity`, `x` is a uniformizer and `y` is a unit.  Indeed, with
`u=t^{-1}`, `\bar p=u^mp(u^{-1})`, and
`\bar q=u^nq(u^{-1})`,

\[
 x=u\frac{\bar p^r}{\bar q^s},\qquad
 y=\frac{\bar p^n}{\bar q^m}.                  \tag{1.3}
\]

The Bezout determinant in (1.1) gives the exact inverse monomials

\[
 \boxed{p=x^{-m}y^s,\qquad q=x^{-n}y^r.}        \tag{1.4}
\]

Here and below `p,q` in (1.4) denote the corresponding functions in the
punctured formal carrier chart; leading target scales can be restored
separately.

## 2. Universal three-parameter affine transport

After removing the two leading scales, every lower-triangular affine target
change has the form

\[
 P=A(p+\mu),\qquad
 Q=B(q+\eta p+\nu).                             \tag{2.1}
\]

Let `X,Y` be the dimensionless carrier coordinates formed from `P,Q` by the
same exponents as in (1.2).  Equations (1.4) give

\[
\boxed{\begin{aligned}
X={}&x\frac{(1+\mu x^m/y^s)^r}
 {(1+\eta x^{n-m}y^{s-r}+\nu x^n/y^r)^s},\\
Y={}&y\frac{(1+\mu x^m/y^s)^n}
 {(1+\eta x^{n-m}y^{s-r}+\nu x^n/y^r)^m}.
\end{aligned}}                                  \tag{2.2}
\]

Thus the dimensionless carrier transport always has exactly three
parameters, independently of `(m,n)`.  For the carrier grading
`wt(x)=1, wt(y)=0`, their weights are

\[
 \boxed{\operatorname{wt}(\eta,\mu,\nu)=(n-m,m,n).} \tag{2.3}
\]

For `(m,n)=(3,5)`, one has `(r,s)=(3,2)`, and (2.2) is exactly the fixed
transport in
[`F2_AFFINE_K1_CARRIER_JET_FACTORIZATION.md`](F2_AFFINE_K1_CARRIER_JET_FACTORIZATION.md).

## 3. The raw-jet discriminant pattern

The standard monic, centered, lower-triangular normal-form slice is

\[
\begin{aligned}
p(t)&=t^m+\sum_{i=1}^{m-2}a_i t^i,\\
q(t)&=t^n+\sum_{\substack{1\le j\le n-1\\j\ne m}}b_jt^j.
\end{aligned}                                    \tag{3.1}
\]

The missing `t^(m-1)` coefficient comes from translating `t`; the constants
and leading coefficients come from target translations and scales; and the
missing `t^m` coefficient of `q` comes from the target shear.  This slice has

\[
 M=(m-2)+(n-2)=m+n-4                         \tag{3.2}
\]

displayed coefficients.  It deliberately retains the residual weighted
rescaling of `t`, just as the F2 four-parameter slice retains the weights of
`a,b,c,d`.

Write the transported carrier graph as

\[
 Y=1+J_1X+J_2X^2+\cdots .                       \tag{3.3}
\]

Adding the three parameters in (2.2) gives the square candidate map

\[
 \Phi_{m,n}:\mathbb A^{m+n-1}longrightarrow
 \mathbb A^{m+n-1},\qquad
 (a_i,b_j,\mu,\eta,\nu)\longmapsto
 (J_1,\ldots,J_{m+n-1}).                        \tag{3.4}
\]

Exact formal-series calculation gives

\[
\begin{array}{c|c}
(m,n)&\det(d\Phi_{m,n})/\operatorname{Res}(p',q')\\ \hline
(2,3)& 2\\
(2,5)&-2\\
(3,4)&-3\\
(3,5)& 3.
\end{array}                                      \tag{3.5}
\]

This supports the following general pattern.

> **Primitive carrier-discriminant conjecture.**  For every coprime
> `1<m<n`, the determinant in (3.4) is
> `epsilon_(m,n)*m*Res(p',q')`, where
> `epsilon_(m,n)` is a nonzero sign depending on conventions.

Only the four rows in (3.5) are currently proved.  In particular, the
general equality of divisors, its multiplicity one, and its sign remain
open.  The invariant formulation should identify the ramification divisor
of the raw carrier-jet chart with the nonimmersion divisor of the polynomial
parametrization.

## 4. Saturation and the one-order gap

The count (3.2) plus the universal three parameters gives the candidate
saturation order

\[
 \boxed{N_*(m,n)=m+n-1.}                        \tag{4.1}
\]

Whenever the determinant in (3.4) is nonzero, the first `N_*` raw carrier
coefficients are local coordinates.  Consequently:

- no compatibility relation among raw jets of order at most `N_*` can hold
  on a dense open set;
- the closure of the `(N_*+1)`-jet image is an irreducible hypersurface; and
- the first target-normalization-invariant raw compatibility can appear only
  at order `N_*+1`.

For `(3,5)`, `N_*=7`.  The F2 carrier refinement has seven prescribed
centers, and contact `8` asks only for equality of `J_1,...,J_7`.  Therefore

\[
 \boxed{\text{the fan stops exactly one coefficient before the first raw
 invariant relation.}}                          \tag{4.2}
\]

This explains the failure of further coefficient factoring more precisely
than parameter counting alone: the existing packet reaches saturation but
never becomes overdetermined.  Computing `J_8` could exhibit the first raw
hypersurface, but that equation is not tested by the present seven-center
fan and cannot by itself create a boundary obstruction there.

## 5. Exact `(3,5)` nonimmersion stratum

For

\[
p=t^3+at,\qquad q=t^5+bt^4+ct^2+dt,
\]

the raw seven-jet discriminant is

\[
\Delta_{\rm imm}=\operatorname{Res}(p',q')
=25a^4+48a^3b^2-144a^2bc+90a^2d+108ac^2+81d^2. \tag{5.1}
\]

A dense parametrization of this irreducible hypersurface marks its unique
common critical point `r`:

\[
 \boxed{a=-3r^2,\qquad
 d=-5r^4-4br^3-2cr.}                            \tag{5.2}
\]

Substitution in the collision quartic gives the exact factorization

\[
 R(u)=(u-2r)T(u),                               \tag{5.3}
\]

where

\[
\boxed{
T(u)=u^3+(b+2r)u^2+(2br+r^2)u
      -2br^2-c+2r^3.}                           \tag{5.4}
\]

The factor `u-2r` is the diagonal pair `(r,r)`.  The ordinary-cusp open
condition is

\[
 \chi=6br^2-c+20r^3\ne0,                       \tag{5.5}
\]

because

\[
 \det\!\begin{pmatrix}p''(r)&p'''(r)\\q''(r)&q'''(r)\end{pmatrix}
 =12\chi.                                       \tag{5.6}
\]

The second possible common critical point is excluded by
`2br^2+c != 0`.  The residual cubic is squarefree when

\[
 (2r^3-3br^2-c)
 (4b^3-12b^2r-69br^2-27c+50r^3)\ne0.           \tag{5.7}
\]

These conditions are simultaneously nonempty.  At

\[
 (r,b,c)=(1,0,1),\qquad(a,b,c,d)=(-3,0,1,-7),  \tag{5.8}
\]

one has

\[
 T=u^3+2u^2+u+1,quad
 \operatorname{Disc}(T)=-23,quad
 \operatorname{Res}(T,3u^2-12)=-513.           \tag{5.9}
\]

The exact tangent and distinct-image tests show that its three residual
roots give three distinct ordinary nodes, all disjoint from the cusp.
Hence a dense open subset of `Delta_imm=0` has precisely

\[
 \boxed{A_2+3A_1}                               \tag{5.10}
\]

as its affine singularity packet.  Its delta ledger is still `1+3=4`.
At the witness (5.8), the raw seven-jet Jacobian has rank six, so its generic
corank along the irreducible discriminant is exactly one.

## 6. The conductor factor records the degeneration, not its `ch_2` length

Let

\[
 C(t)=\operatorname{Res}_u
 \left(R(u),t^2-ut+(u^2+a)\right)                \tag{6.1}
\]

be the degree-eight affine conductor polynomial.  Equations (5.2)--(5.4)
give

\[
 \boxed{C(t)=(t-r)^2C_6(t),\qquad \deg C_6=6.} \tag{6.2}
\]

Thus the two simple conductor points of one generic node collide to the
double conductor point of the cusp, while the other six points remain the
three nodal pairs.  At (5.8),

\[
C_6=t^6+2t^5-6t^4-12t^3+11t^2+22t+1,           \tag{6.3}
\]

with nonzero discriminant and `C_6(1)=19`.  The conductor scheme therefore
has exactly the claimed `2+6` split on a nonempty open subset of the
nonimmersion divisor.

This is the precise local pattern relevant to the logarithmic `ch_2`
program:

\[
\begin{array}{c|c|c|c}
\text{delta unit}&\text{normalization branches}&\text{conductor preimage}
&\text{raw jet rank loss}\\ \hline
A_1&2&1+1&0\\
A_2&1&2&1.
\end{array}                                      \tag{6.4}
\]

The last column is proved here only for the generic `(3,5)` raw seven-jet
chart.  It belongs to the parameter map of target curves, not to the
logarithmic cotangent complex of the Keller map.  In fact the subsequent
[`unibranch attachment theorem`](LOG_UNIBRANCH_ATTACHMENT_FITTING.md) proves
that a minimal transverse boundary attachment above a branch of multiplicity
`m_C` and local residue index `q_p` has point correction `q_p*m_C`.  For the
ordinary cusp this is `2q_p`, not the jet corank one.  Affine source preimages contribute
zero because the logarithmic complex is exact on `A^2`; only source-boundary
attachments can carry this class.

The
[`affine strict-log-étale resolution theorem`](AFFINE_KELLER_STRICT_LOG_ETALE_RESOLUTION.md)
makes that last statement resolution-independent: the entire embedded
resolution of the `A_2+3A_1` packet is an étale base change, so its relative
logarithmic cokernel is zero.  The coefficient-space corank and affine
conductor split cannot become a source-surface Fitting class without a
separately located compactification-boundary attachment.

<!-- status-consumer: PAER1 60eb24b2232d159e -->

## 7. Generalized gaps and the useful next computations

The factor pattern separates six genuinely different gaps.

1. **All-bidegree determinant theorem.**  Prove the conjecture after
   quotienting the residual weighted reparametrization correctly.  A proof
   should explain multiplicity one, not just reproduce more symbolic rows.
2. **First invariant versus available fan.**  For `(3,5)`, compute the raw
   eighth-jet hypersurface only if a global source refinement actually tests
   it.  The current seven-center fan does not.
3. **Global normalization fixing.**  Relate `mu,eta,nu` to the four affine
   singular fibers or to the completed fixed-coordinate source pair.  Until
   then the first seven raw centers are coordinates, not obstructions.
4. **Locate and test the boundary attachment.**  The local algebra is now
   exact: a minimal transverse SNC attachment contributes `q_p*m_C`, hence
   `2q_p` over the ordinary cusp and total `2f` over a complete
   residue-degree-`f` cusp fiber.  The missing source compiler must locate the
   boundary points, determine their residue indices, and verify transverse
   unimodularity and minimal residual-Jacobian order.  Jet corank and
   conductor degree cannot substitute for those checks.
5. **Weighted source fiber counts.**  Determine the numbers `N_i` above the
   target singular values.  The conserved total delta four only supplies the
   bound `sum N_i delta_i<=4(d-1)`.
6. **Complete global filtration.**  Put the boundary divisorial Smith
   packets and boundary node/cusp gluing quotients on one common SNC model
   before applying positivity or parity to the global `ch_2` remainder.
   Retain affine conductor lengths only in the separate fiber/escape ledger.

The immediate productive branch is Gap 4, not a ninth blind Laurent
coefficient: use the cusp value to locate its boundary preimages, then extract
the residue index and the two transverse first derivatives required by the
unibranch attachment theorem.

<!-- status-consumer: LUAF1 b0279670ffbd3fa5 -->

## Reproduction

```bash
.venv/bin/python scripts/verify_coprime_carrier_jet_discriminant_pattern.py
```

The checker verifies the universal Bezout/transport exponents, recomputes
the first three new determinant rows in (3.5) (the `(3,5)` row is supplied
by `PF2K1JF1`), factors the complete `(3,5)` nonimmersion collision and
conductor packets, proves the cusp-plus-three-node witness, and checks raw
seven-jet corank one there.
