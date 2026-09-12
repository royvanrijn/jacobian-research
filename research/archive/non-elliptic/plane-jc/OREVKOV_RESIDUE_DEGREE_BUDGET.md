# Orevkov's residue-degree budget and the clean F2 cusp atlas

> **Status.**  The residue-normalized global identity in Sections 1--3 is an
> unconditional consequence of Orevkov's Euler-multiplicity formula and
> Riemann--Hurwitz.  For every dicritical row of transverse/residue degrees
> `(e,f)`, the generic term `e` and the forced finite residue ramification
> together cost at least `e*f`.  More precisely,
>
> \[
> d-1=\sum_i e_i f_i+\sum_{i,x}\epsilon_{i,x},
> \qquad \epsilon_{i,x}\ge0.
> \]
>
> For F2, the rightmost sum is exactly the numerical residual
> `d-1-sum_i(e_i*f_i)` in the logarithmic complete-chain identity.  This
> removes the former sign hypothesis from that **total numerical
> inequality**, but it does not identify the individual logarithmic modules.
> The branch-multiplicity refinement in Section 5 requires the displayed
> regular/Cartier local hypothesis.  Section 6 is an exact specialization of
> Orevkov's 2026 classification of smooth, bijective, irreducibly ramified
> germs; it shows that the only nontrivial clean `(3,5)` cusp row is
> `(local degree, ramification order)=(15,7)`.  Singular normalization points,
> non-bijective residue maps, and reducible ramification remain outside that
> classification.  The same arithmetic enumeration classifies all five
> irreducible cusp types on the `k=1` F2 atlas: `A2`, `A4`, `A6`, `E6`, and
> `E8`.  In particular, a single clean ramification row cannot cover both
> the `A4` and `A2` points of the `A4+A2+A1` target, because their possible
> nontrivial ramification orders are disjoint.

The arithmetic is replayed by
[`verify_orevkov_residue_degree_budget.py`](../scripts/verify_orevkov_residue_degree_budget.py).

## 1. Orevkov's global identity

Let

\[
 F:\mathbb C^2\longrightarrow\mathbb C^2
\]

be a polynomial local biholomorphism of geometric degree `d`.  Use the
regularization and topological contractions in Orevkov's construction.
Orevkov writes `L_infinity` for the inverse image of the target line at
infinity, `L_FC=L-L_infinity`, `L_C` for the components constant on the
affine target, and

\[
 L_F=L_{FC}-L_C.                                                \tag{1.0}
\]

The curves denoted below by `l` are **only** the components
`l` contained in `L_F`:
the affine exceptional/dicritical curves mapping nonconstantly to the
nonproperness set.  The terminal F2 row over target infinity is in
`L_infinity` and is not a summand in this formula.  This scope is explicit in
Orevkov's Lemma 4.2 and prevents the terminal `(e,f)=(1,6)` packet from being
incorrectly added to the affine budget.

Each such `l` is rational, has one point over infinity, and
`l^o=l-{infinity}` is isomorphic to `A1`.  Write `mu_l` for the generic local
surface degree of the contracted map along `l`, and `mu_x` for its local
degree at a finite point of `l`.  Orevkov's Lemma 4.2 is

\[
 \boxed{
 d-1=\sum_l\left(
  \mu_l+\sum_{x\in l^o}(\mu_x-\mu_l)
 \right).}                                      \tag{1.1}
\]

Only finitely many inner summands are nonzero, and upper semicontinuity gives
`mu_x>=mu_l`.

At a generic smooth point of `l`, the surface map has transverse local form
`(s,r)->(s,r^e)`.  Hence

\[
 \mu_l=e_l,                                     \tag{1.2}
\]

where `e_l` is the transverse index used in the F2 boundary ledger.

## 2. The residue cover supplies the missing factor `f`

Let `C_l` be the image component and let `C_l^nu=A1` be its affine
normalization.  The restriction to `l^o` factors as a polynomial map

\[
 g_l:\mathbb A^1\longrightarrow\mathbb A^1
\]

of degree `f_l`, the residue degree of the dicritical row.  Its extension
`P1->P1` has a unique point over infinity and is totally ramified there.
Riemann--Hurwitz therefore gives

\[
 \sum_{x\in l^o}(q_x-1)=f_l-1,                 \tag{2.1}
\]

where `q_x` is the local degree of `g_l` and equals one away from finitely
many points.

Local-degree conservation gives

\[
 \mu_x\ge e_lq_x.                              \tag{2.2}
\]

Indeed, after moving to a nearby smooth value of the image branch, the
residue point splits into `q_x` points on `l`, each carrying transverse
degree `e_l`; any other local sheets only increase `mu_x`.

Define the residue-normalized local excess

\[
 \epsilon_{l,x}:=\mu_x-e_lq_x\ge0.             \tag{2.3}
\]

Substituting (1.2), (2.1), and (2.3) into the contribution of `l` in (1.1)
gives the exact decomposition

\[
\begin{aligned}
 e_l+\sum_x(\mu_x-e_l)
 &=e_l+e_l\sum_x(q_x-1)+\sum_x\epsilon_{l,x}\\
 &=e_lf_l+\sum_x\epsilon_{l,x}.                \tag{2.4}
\end{aligned}
\]

Thus residue degree is not an optional correction to Orevkov's formula: it
is exactly the finite Riemann--Hurwitz part of that formula.

## 3. The global residue-degree theorem

Summing (2.4) gives

\[
 \boxed{
 d-1=\sum_l e_lf_l+\sum_{l,x}\epsilon_{l,x},
 \qquad \epsilon_{l,x}\ge0.}                  \tag{3.1}
\]

In particular,

\[
 \boxed{\sum_l e_lf_l\le d-1.}                 \tag{3.2}
\]

This is stronger than retaining only Orevkov's generic inequality
`sum_l e_l<=d-1`.  Equality in (3.2) means precisely that every finite
special point has no excess beyond the ramification of the residue cover.

Group the rows over affine nonproperness components `C_j` and write

\[
 A_j=\sum_{l\mapsto C_j}e_lf_l.
\]

Then (3.1) becomes

\[
 \boxed{d-1-\sum_jA_j=\sum_{l,x}\epsilon_{l,x}\ge0.} \tag{3.3}
\]

No cyclicity, Fitting-ideal sign, or logarithmic determinant normalization
is used in this deduction.

## 4. Identification with the F2 complete-chain remainder

Every F2 nonproperness component has rational one-puncture normalization,
so its logarithmic canonical degree is `beta_j=-1`.  The complete-chain
calculation in
[`F2_AFFINE_GLOBAL_RAMIFICATION_BUDGET.md`](F2_AFFINE_GLOBAL_RAMIFICATION_BUDGET.md)
gives the numerical remainder

\[
 B_f-\left(\deg K_{\rm aff}+\frac12D_{\log}^2\right)
 =d-1-\sum_jA_j.                                \tag{4.1}
\]

Combining (3.3) and (4.1) identifies it with

\[
 \boxed{\sum_{l,x}\epsilon_{l,x}.}             \tag{4.2}
\]

This is the sought conceptual bridge between Orevkov's topology and the
localized `ch_2` budget.  It proves that the **total number** is
nonnegative even when a contracted logarithmic cokernel is generically
noncyclic.  It does not say that a scalar determinant module splits into the
correct local summands, nor does it identify which contracted curve carries
which `epsilon`.  The glued-versus-split node distinction still belongs to
the full perfect complex and its `Fitt_1` filtration.

For one target component, if a meridian fixes `u` sheets, then
`A=d-u`, and (3.3) reads

\[
 \sum_{l,x}\epsilon_{l,x}=u-1.                 \tag{4.3}
\]

The degree-six cubic-inertia E8 row has `(e,f)=(3,1)`, `u=3`, and hence
exact excess budget two.  This is the same equality row isolated by the
logarithmic calculation, but now its nonnegativity is unconditional.

## 5. What ordinary branch multiplicity can force

There is a useful local refinement, with an important hypothesis.  Suppose
`x` is represented by a finite map of regular surface germs, the relevant
source branch is the Cartier divisor `E=(r=0)`, and the image branch has
multiplicity `m_C`.  If its normalization parameter pulls back with local
degree `q_x`, then in `R/(r)=k[[t]]` the target maximal ideal has order
`q_x*m_C`.  The quotient map

\[
 R/(f^*\mathfrak m)\twoheadrightarrow
 R/(r,f^*\mathfrak m)
\]

gives

\[
 \mu_x\ge q_xm_C.                              \tag{5.1}
\]

Together with (2.2),

\[
 \epsilon_{l,x}\ge q_x\max(0,m_C-e_l).        \tag{5.2}
\]

For the `(3,5)` cusp, `m_C=3`.  A simple-inertia row `e=2` therefore has
at least `q_x` excess at every regular Cartier point above the cusp.  A
complete cusp fiber on a residue-degree-`f` row costs at least

\[
 ef+f=3f.                                       \tag{5.3}
\]

Applied to the already enumerated one-component simple-inertia E8 actions
through the F2 degree ceiling, the necessary inequality is

\[
 R:=\sum f_l\le u-1.                            \tag{5.4}
\]

Every row fails it: `(d,u,R)` is one of

\[
 (6,2,2),(10,2,4),(12,4,4),(15,3,6),
 (20,4,8),(24,4,10).
\]

This recovers their exclusion whenever the cusp points are regular Cartier
points.  It must not be applied silently at a singular point of Orevkov's
contracted source.  Such a singular/contracted attachment is precisely where
the positive-dimensional logarithmic Smith packet can live.

For cubic inertia `e=3`, (5.2) is neutral.  The second cusp exponent five,
not just branch multiplicity three, must then enter.  That is the sharp local
gap exposed by the degree-six equality row.

## 6. The smooth clean cusp germs form a finite atlas

Orevkov's 2026 theorem classifies finite analytic germs between smooth
surface germs which are ramified along one irreducible curve, whose image is
`u^d1=v^d2`, and whose restriction to every component over that curve is
bijective.  Solving both arithmetic families in his Theorem 2, including the
swap `d1<->d2`, gives the following complete list for the irreducible cusp
types occurring in the `k=1` target atlas.  The identity row `(N,n)=(1,1)`
occurs in every line and is suppressed.

\[
\begin{array}{c|c|c|c}
\text{type}&(d_1,d_2)&\text{nontrivial }(N,n)&
 \text{passport }\alpha\mid\beta\mid(n,1^{N-n})\\ \hline
A_2&(2,3)&(3,2)&(3)\mid(2,1)\mid(2,1)\\
   &&(6,4)&(3,3)\mid(2,2,2)\mid(4,1,1)\\
A_4&(2,5)&(5,3)&(5)\mid(2,2,1)\mid(3,1,1)\\
   &&(10,6)&(5,5)\mid(2^5)\mid(6,1^4)\\
A_6&(2,7)&(7,4)&(7)\mid(2,2,2,1)\mid(4,1^3)\\
   &&(14,8)&(7,7)\mid(2^7)\mid(8,1^6)\\
E_6&(3,4)&(4,2)&(4)\mid(3,1)\mid(2,1,1)\\
   &&(6,3)&(3,3)\mid(4,2)\mid(3,1^3)\\
   &&(12,6)&(4,4,4)\mid(3^4)\mid(6,1^6)\\
E_8&(3,5)&(15,7)&(5,5,5)\mid(3^5)\mid(7,1^8).
\end{array}                                                    \tag{6.1}
\]

Here `N` is the local covering degree and `n` is the ramification order.
Swapping target coordinates interchanges the first two passport columns.
For one clean germ, its Orevkov contribution is exactly `N`: the generic
term is `n` and the cusp jump is `N-n`.  Thus a row of local degree `N`
forces the global inequality `d>=N+1`.

The `E8` line is especially rigid.  A clean smooth-source bijective `(3,5)`
germ cannot have ramification order `2,3,4,5`, or `6`; its only nonidentity
row is `(N,n)=(15,7)` and forces `d>=16`.  Any lower-order F2 `E8` packet
must leave at least one hypothesis of the classification: it has a singular
finite-normalization point, a non-bijective residue branch, reducible
ramification, or another local branch component.  The explicit cubic F2
model does exactly this after resolution: it contains a contracted divisor
with generic logarithmic Smith form `diag(t,t)`.

There is also a useful compatibility obstruction inside the severe atlas.
A single irreducible ramification divisor has one generic transverse order
along its target component.  The possible nontrivial clean orders at an
`A4` point are `{3,6}`, whereas those at an `A2` point are `{2,4}`.  Hence

\[
 \boxed{\text{the `A4+A2+A1` component cannot be clean at both cusps on
 one bijective irreducible ramification row}.}             \tag{6.2}
\]

At least one cusp must instead be singular/non-Cartier upstairs,
non-bijective, or supported by a different ramification component.  This is
a classification of the clean escape routes, not a closure of all cusp
attachments.

## 7. Current classification and next finite target

The cusp problem now has three disjoint regimes.

| regime | disposition |
| --- | --- |
| residue ramification only | booked exactly by `e*f` in (3.1) |
| regular Cartier cusp point with `e<m_C` | positive excess from (5.2) |
| smooth, bijective, irreducibly ramified `k=1` cusp germ | finite atlas (6.1) |
| singular/non-Cartier or reducible ramification | remaining local frontier |

For `(75,125)`, the next useful computation is therefore not another scalar
determinant expansion.  It is a finite classification of the normal surface
singularities, non-bijective residue maps, and reducible ramification graphs
over the five rows in (6.1), with local degree at most `28`, followed by
their compatibility with the fixed terminal A6 section and the compiled
source intersection graph.

## Sources

- S. Yu. Orevkov,
  [*On three-sheeted polynomial mappings of C2*](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  Lemma 4.2 and Corollary 4.3.
- Nguyen Van Chau,
  [*Non-zero constant Jacobian polynomial maps of C2*](https://matwbn.icm.edu.pl/ksiazki/apm/apm71/apm7135.pdf),
  Remark 4.9.
- S. Yu. Orevkov,
  [*On germs of mappings C2 to C2*](https://www.math.univ-toulouse.fr/~orevkov/k-en.pdf),
  Theorem 2.

## Reproduction

```bash
.venv/bin/python scripts/verify_orevkov_residue_degree_budget.py
```
