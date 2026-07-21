# The weighted family as a normalized marked-root space

This note generalizes the marked-root viewpoint from the exceptional cubic to
the admissible weighted inverse-pencil family.  It packages formulas already
used in [WEIGHTED_SEED_THEOREM.md](WEIGHTED_SEED_THEOREM.md) and the divisor
classification in
[DICRITICAL_COMPACTIFICATION.md](../experimental/geometry/DICRITICAL_COMPACTIFICATION.md).  The main
new point is the correct global formulation: the raw simple-root locus works
over `C!=0`, but over `C=0` one must normalize the root incidence and retain
exactly the branches on which reconstruction is regular.

Universally, this is the pullback of the marked linear-factor map
`mu_(1,n-1)` along

\[
 (A,B,C)\longmapsto H(W)-BCW+cAC^2.
\]

The parameter space is three-dimensional even though the polynomial depends
on it through two pencil coordinates: `C` retains the reconstruction scale.
This factorization-slice interpretation is developed in
[the archived factorization-slice note](../archive/transfer-program/cubic-factorization-obstruction.md).

Work over a characteristic-zero field `k`.  Let `H in k[W]` be an admissible
primitive of degree `n`, with

\[
H(0)=H(1)=0,\qquad H'(0)=0,\qquad H'(1)=-c\ne0.
\]

Put

\[
\kappa={H''(1)\over c}\ne-2,\qquad
a_0=-{1+\kappa\over2+\kappa},\qquad b_0\ne0.
\]

For source coordinates `(x,y,z)`, the weighted construction uses

\[
v=xy,\quad S=x^2z,\quad u=1+v,
\quad\gamma=1+a_0v+b_0S,\quad W=u\gamma.            \tag{1}
\]

Write the associated polynomial map as `G_H(x,y,z)=(A,B,C)`.

## Polynomiality and constant Jacobian

Put `p=H'` and

\[
q(W)={Wp(W)-H(W)\over c}.
\]

Because `p(0)=0`, the polynomial `p(W)` is divisible by `W`; hence
`p(W)/gamma` is polynomial after `W=u gamma`.  Because
`H(0)=H'(0)=0`, `q(W)` is divisible by `W^2`, so
`q(W)/gamma^2` is polynomial as well.  It remains to justify the displayed
divisions by `x` and `x^2` in (4).

At `x=0` one has `(u,gamma,W)=(1,1,1)`.  Thus

\[
c+{p(W)\over\gamma}=c+p(1)=0,
\]

so the numerator defining `B` is divisible by `x`.  Write
`kappa=p'(1)/c`.  Since `q(1)=-1` and `q'(1)=kappa`, differentiation at
`x=0` gives

\[
{\partial\over\partial x}
\left(u+{q(W)\over\gamma^2}\right)_{x=0}
=\bigl(1+\kappa+a_0(2+\kappa)\bigr)y=0.
\]

The numerator itself also vanishes at `x=0`, so it is divisible by `x^2`.
Consequently all three expressions in (4) are polynomials.

There is also a denominator-free structural determinant proof.  On the dense
open `C=x gamma ne 0`, introduce

\[
s=BC=c\gamma+p(W),\qquad
t=cAC^2=cW\gamma+Wp(W)-H(W).
\]

The three coordinate changes have determinants

\[
\det{\partial(x,v,S)\over\partial(x,y,z)}=x^3,
\]

\[
\det{\partial(W,\gamma,C)\over\partial(x,v,S)}
=b_0\gamma^2,
\qquad
\det{\partial(s,t,C)\over\partial(W,\gamma,C)}
=-c^2\gamma,
\]

and

\[
\det{\partial(s,t,C)\over\partial(A,B,C)}=-cC^3.
\]

Since `C=x gamma`, the chain rule gives

\[
\boxed{\det DG_H=b_0c.}                            \tag{2a}
\]

Both sides are polynomials, so equality on `C ne 0` proves the identity on
all of `A^3`.  In particular `G_H` is etale and quasi-finite without any
appeal to later boundary theory.

## The finite marked-root incidence

For a target `(A,B,C)`, define

\[
E_{A,B,C}(W)=H(W)-BCW+cAC^2.                        \tag{2}
\]

Its projective homogenization has value `h_n W^n` on the line at infinity,
where `h_n` is the nonzero leading coefficient of `H`.  It therefore has no
projective root at infinity.  The incidence

\[
\mathcal I_H=
\{((A,B,C),W)\in\mathbb A^3\mathbin\times\mathbb A^1:
E_{A,B,C}(W)=0\}                                   \tag{3}
\]

is the entire projective root incidence.  Its projection to the target is
finite flat of degree `n`: after division by `h_n`, equation (2) is monic in
`W`, so its coordinate ring is free with basis `1,W,...,W^(n-1)`.

The source marks the root (1).  Indeed, with `q` as above, the weighted map is
defined in invariant coordinates by

\[
C=x\gamma,\qquad
B={c+H'(W)/\gamma\over x},\qquad
A={u+q(W)/\gamma^2\over x^2}.                      \tag{4}
\]

Equations (1) and (4) give the polynomial identity

\[
E_{G_H(x,y,z)}(W(x,y,z))=0.                         \tag{5}
\]

Thus there is a marked-root morphism

\[
\iota_0:\mathbb A^3_{x,y,z}\longrightarrow\mathcal I_H.
                                                                    \tag{6}
\]

## Reconstruction on `C!=0`

On the incidence, put

\[
E_W=H'(W)-BC,
\qquad
\gamma=-{E_W\over c}={BC-H'(W)\over c}.             \tag{7}
\]

Whenever `C!=0` and the marked root is simple, so `E_W!=0`, reconstruct

\[
x={C\over\gamma},\qquad
u={W\over\gamma},\qquad
v=u-1,
\qquad y={v\over x}={W-\gamma\over C},              \tag{8}
\]

and

\[
S={\gamma-1-a_0v\over b_0},
\qquad z={S\over x^2}.                              \tag{9}
\]

These expressions have no denominator beyond the units `C`, `E_W`, `b_0`,
and `c`.  Conversely, on a source point with `C=x gamma!=0`, both `x` and
`gamma` are nonzero and

\[
E_W(W(x,y,z))=-c\gamma\ne0.                         \tag{10}
\]

Substitution of (8)--(9) in (4), using `E=0`, returns `(A,B,C)` and recovers
`W=u gamma`.  Hence

\[
G_H^{-1}(C\ne0)\ \cong\
\mathcal I_H^{\rm simp}\cap\{C\ne0\}.              \tag{11}
\]

Under (11), `G_H` forgets the marked simple root.  The repeated-root divisor
`E=E_W=0` is the ramification divisor of the finite root cover, while the
affine weighted Keller map is etale.  It is therefore not a branch divisor of
`G_H`.

## Why the raw simple-root locus is not the global answer

At `C=0`, equation (2) specializes to `H(W)=0`.  Two effects prevent a global
replacement of the source by `I_H^simp`.

First, the root `W=0` has multiplicity `m_0>=2`, but after normalization some
of its branches reconstruct to finite points on the source chart `gamma=0`.
Thus finite source points can lie over a repeated root of the specialized
pencil.

Second, an additional root `rho!=0,1` of `H` can be simple while (8) still has
a pole: generically `W-gamma` does not vanish with `C`, so `y` diverges.
Thus a simple specialized root can represent a boundary point rather than an
affine source point.

Both phenomena are resolved by normalizing the marked-root incidence and
testing regularity of all reconstruction coordinates.

## The global normalized marked-root theorem

The polynomial (2) is irreducible in `k[A,B,C,W]`: as a polynomial in `A`
over the UFD `k[B,C,W]`, it is primitive and linear over the fraction field.
Indeed, `cC^2` and `H(W)-BCW` are relatively prime.  Let

\[
\nu:\widetilde{\mathcal I}_H\longrightarrow\mathcal I_H              \tag{12}
\]

be its normalization.  This morphism is finite because the incidence is of
finite type over a field, so \(\widetilde{\mathcal I}_H\) remains finite over
the target.  The rational functions (7)--(9) belong to the common function
field of `A^3` and `I_H`.  Let `R_H` be their simultaneous regularity locus on
\(\widetilde{\mathcal I}_H\):

\[
\mathcal R_H=
\{p\in\widetilde{\mathcal I}_H:x,y,z
\text{ are regular at }p\}.                         \tag{13}
\]

This is an open subset.

**Theorem.**  The marked-root morphism (6) lifts uniquely to an isomorphism

\[
\iota:\mathbb A^3_{x,y,z}\xrightarrow{\sim}\mathcal R_H.
                                                                    \tag{14}
\]

Under this identification, `G_H` is the finite marked-root projection
restricted to the regular-reconstruction open.  Over `C!=0`, equation (14)
is exactly the simple-root isomorphism (11).

**Proof.**  The source is normal, and (6) is dominant and birational by the
explicit inverse on the dense open `C!=0`.  It therefore factors uniquely
through the normalization (12).  The lift is quasi-finite because the
construction lemma gives `det(DG_H)=b_0c!=0`, so `G_H` is etale and
quasi-finite, and every fiber of the lift is contained in a fiber of `G_H`.
A quasi-finite birational morphism to a normal variety is an open immersion by
Zariski's Main Theorem.

On its image, the rational functions (7)--(9) pull back to the source
coordinates, so the image lies in `R_H`.  Conversely, wherever these three
functions are regular they define a morphism to `A^3`.  The identities
`G_H(x,y,z)=(A,B,C)` and `W=(1+xy)gamma`, verified on the dense open `C!=0`,
hold in the common function field and hence on that regularity neighborhood.
This morphism is inverse to the open immersion.  Therefore its domain is
exactly `R_H`, proving (14).  \(\square\)

This theorem makes the global regularity criterion explicit.  It does not
assert that every simple root at `C=0` reconstructs, nor that every repeated
root there is lost.

## The finite charts over `C=0`

The distinguished root `W=1` is simple because `H'(1)=-c`.  In the completed
etale local incidence chart, equation (2) gives

\[
W=1-{B\over c}C+
\left(A+{2+\kappa\over2c^2}B^2\right)C^2+O(C^3).   \tag{15}
\]

Equations (7)--(9) then show

\[
x\longrightarrow0,qquad
y\longrightarrow-{2+\kappa\over c}B.              \tag{16}
\]

The definition of `a_0` cancels the coefficient of `C` in
`gamma-1-a_0v`; consequently `S=O(C^2)` and `z=S/x^2` is regular.  This is the
general finite `x=0` chart.

Write

\[
H(W)=W^{m_0}L(W),\qquad L(0)=h_0\ne0.
\]

The finite zero-cluster chart uses `W=CR`.  After division by `C^2`, its
strict equation is

\[
C^{m_0-2}R^{m_0}L(CR)-BR+cA=0.                    \tag{17}
\]

For `m_0=2`, its generic special fiber is

\[
h_0R^2-BR+cA=0,
\]

and gives the two generic finite `gamma=0` branches.  For `m_0>=3`, its
generic special fiber is `R=cA/B` and gives the single finite zero-cluster
branch.  Repeated points of (17) can still meet the discriminant boundary;
the theorem retains only its regular-reconstruction open.

## Pole divisors and C16

Because \(\widetilde{\mathcal I}_H\) is normal, the failure locus of the rational
reconstruction functions is the union of their polar prime divisors.  The
valuation calculation of Claim C16 identifies exactly those divisors:

1. the discriminant divisor `D_Delta`, where `E_W=0` over `C!=0`;
2. for every additional nonzero primitive root `rho_j!=1`, the divisor
   `D_(rho_j)` over `C=0`, including the case in which `rho_j` is simple;
3. when `m_0>=3`, the escaping zero-cluster divisor `D_0`.

The distinguished root-one divisor and the finite component(s) of (17) lie
in `R_H` generically.  The remaining zero-cluster normalization has scale

\[
C=\delta^{m_0-1},\qquad W=\eta\delta,
\qquad \eta^{m_0-1}=B/h_0,                          \tag{18}
\]

and is `D_0`; additional roots have the Kummer scales recorded in C16.  Thus
C16 is precisely the divisorial classification of the complement
\(\widetilde{\mathcal I}_H\setminus\mathcal R_H\).  On the normalized
inverse-graph compactification, these
same valuations are the dicritical boundary divisors.

The Kummer exponents in (18) are ramification indices.  Each corresponding
prime divisor has residue degree one over `C=0`; the explicit blow-ups and
the intervening nondicritical `A`-chains are computed in
[C16_BLOWUP_GEOMETRY.md](../experimental/geometry/C16_BLOWUP_GEOMETRY.md).

The exact image and omission statements remain seed-dependent.  The theorem
organizes their inverse geometry but does not replace the boundary arguments
in the canonical, deformed, and repeated-root notes.

## Relation with the exceptional cubic

For the exceptional cubic, the alternative binary root `[U:V]` used in
[MARKED_ROOT_MODEL.md](MARKED_ROOT_MODEL.md) already packages the affine
source as a raw simple-root locus, including its root at infinity.  For a
general weighted seed the natural root `W` has no projective root at infinity,
but its incidence is singular over the multiple root `W=0`; normalization and
the regular-reconstruction condition in (13) are the required replacements.

Run the exact algebraic audit with

```bash
.venv/bin/python scripts/verify_weighted_marked_root_model.py
```

It checks the universal reconstruction identities, the root-one cancellation,
the normalized zero-cluster equations, and the complete construction for the
canonical inverse quartic.
