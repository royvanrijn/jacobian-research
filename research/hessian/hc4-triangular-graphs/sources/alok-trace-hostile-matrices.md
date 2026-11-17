# Hostile matrices for the surviving trace lattices

## Claim boundary

The determinant, symmetry, normalization-module cokernel, corank, theta, and
trace-unit constraints derived in
[`finite-flat-trace-lattices.md`](finite-flat-trace-lattices.md) do **not** by
themselves eliminate the surviving one-dicritical passports.  This note gives
exact polynomial matrices realizing all of those quadratic constraints for:

- an `A6` branch with one finite `T(2,5)` cusp and smooth normalization
  collisions; and
- an `S6` branch with the full `1+1+1` jump partition, realized by three
  `T(2,3)` cusps and three smooth normalization collisions.

These are hostile trace-lattice models.  They are not rank-six commutative
algebras, finite covers, or Keller maps.  Their purpose is to prove that the
listed quadratic data alone are insufficient.  A later complement-group audit
now excludes both particular curves globally.  The `S6` curve's four-meridian
projection cannot generate a transitive six-sheet group using transpositions.
The `A6` curve has affine complement group `Z`, so its 3-cycle meridians have
only cyclic image.  These facts do not invalidate either matrix calculation,
but they show that connected monodromy is another essential datum alongside
multiplication, associativity, the regular representation, and
compactification data.

## 1. A square-discriminant `A6` model

Let

\[
R=\mathbb C[P,Q]
\]

and parametrize a rational plane curve by

\[
\nu(t)=(P(t),Q(t))=(t^2+t^3,t^5).
\tag{1.1}
\]

Its implicit equation is

\[
b=-P^5+5P^3Q-5PQ^2+Q^3+Q^2.
\tag{1.2}
\]

At `t=0`, (1.1) has orders `(2,5)`.  Its only other finite singularities are
two smooth two-preimage collisions at

\[
Q=1,
\qquad
P^2+P-1=0.
\]

Indeed, the paired parameters are nontrivial fifth roots `t` and `t^{-1}`.

Consider the symmetric matrix

\[
K=
\begin{pmatrix}
-P^2 & P^2-Q & P-Q\\
P^2-Q & -Q & -P\\
P-Q & -P & P-1
\end{pmatrix}.
\tag{1.3}
\]

Exact calculation gives

\[
\det K=b,
\qquad
K(P(t),Q(t))
\begin{pmatrix}1\\t\\t^2\end{pmatrix}=0.
\tag{1.4}
\]

The generators `1,t,t^2` span `C[t]` over `R`: the identities

\[
t^3=P-t^2,
\qquad
t^4=Pt-P+t^2,
\qquad
t^5=Q
\]

reduce all higher powers.  The standard Bezoutian presentation then gives

\[
\operatorname{coker}K\cong\nu_*\mathcal O_{\mathbb A^1}.
\tag{1.5}
\]

Now set

\[
\Phi=\operatorname{diag}(K,I_3),
\qquad
H=
\begin{pmatrix}0&I_3\\I_3&0\end{pmatrix},
\qquad
T=\Phi^tH\Phi=
\begin{pmatrix}0&K\\K&0\end{pmatrix}.
\tag{1.6}
\]

Then

\[
\det\Phi=b,
\qquad
\det H=-1,
\qquad
\det T=-b^2,
\tag{1.7}
\]

which is the `A6` factorization up to a nonzero scalar unit.  The exact ranks
are:

| target position | `rank K` | `rank Phi` | `rank T` |
|---|---:|---:|---:|
| off `b=0` | 3 | 6 | 6 |
| generic point of `b=0` | 2 | 5 | 4 |
| `T(2,5)` cusp | 1 | 4 | 2 |
| either normalization collision | 1 | 4 | 2 |

Thus `Phi` has the required coranks one and two, while the trace form has
coranks two and four.  Even the distinguished unit line can be modeled.  In
the two three-dimensional blocks, take

\[
u=(e_3,-3(e_2+e_3)).
\]

It is primitive and

\[
u^tTu=6.
\tag{1.8}
\]

Consequently the unit line and a free rank-five orthogonal complement cannot
provide the missing contradiction either.

## 2. A full `S6` `1+1+1` model

Take

\[
P(t)=t^4-6t^2,
\qquad
Q(t)=t^5-5t^3.
\tag{2.1}
\]

The common derivative factor is

\[
\gcd(P',Q')=t(t^2-3).
\]

The three critical parameters `0,+sqrt(3),-sqrt(3)` give `T(2,3)` cusps at

\[
(0,0),
\qquad
(-9,-6\sqrt3),
\qquad
(-9,6\sqrt3).
\tag{2.2}
\]

The curve also has three ordinary normalization collisions at

\[
(-5,0),
\qquad
(-4,2\sqrt2),
\qquad
(-4,-2\sqrt2).
\tag{2.3}
\]

Its implicit equation is

\[
\begin{aligned}
b={}&P^5+10P^4+25P^3+10P^2Q^2\\
   &+90PQ^2-Q^4+216Q^2.
\end{aligned}
\tag{2.4}
\]

Define

\[
B=
\begin{pmatrix}
0&6Q&-5P&-Q&P\\
6Q&-5P&-Q&P&0\\
-5P&-Q&P-30&0&6\\
-Q&P&0&1&0\\
P&0&6&0&-1
\end{pmatrix}.
\tag{2.5}
\]

It satisfies

\[
\det B=b,
\qquad
B(P(t),Q(t))
\begin{pmatrix}1\\t\\t^2\\t^3\\t^4\end{pmatrix}=0,
\qquad
\operatorname{coker}B\cong\nu_*\mathcal O_{\mathbb A^1}.
\tag{2.6}
\]

The rank-six symmetric matrix

\[
T_{S_6}=\operatorname{diag}(B,1)
\tag{2.7}
\]

has determinant `b`, rank five at a generic branch point, and rank four at
all six points in (2.2)--(2.3).  Thus its coranks are exactly one generically
and two at every cusp or normalization collision.  The vector
`sqrt(6)*e_6` has norm six, so the trace-unit constraint also survives.

## 3. Compactified theta symmetry also survives

Every entry of (2.5) is affine linear in `P,Q`.  Homogenizing the constants
produces a symmetric `5 x 5` matrix of projective linear forms whose
determinant is the projective closure of (2.4).  Its cokernel is the pushforward
of

\[
\mathcal O_{\mathbb P^1}(-1),
\]

the unique theta characteristic on the rational normalization.  Its
restriction to `A1` is trivial, exactly as required by the affine
normalization module.

This compatibility is structural.  Beauville's Theorem B turns self-dual ACM
sheaves into symmetric resolutions, and his Remark 4.4 notes that an integral
plane curve obtains such a representation by pushing forward an ineffective
theta characteristic from its normalization.  Therefore symmetry or theta
parity cannot exclude these rational branch curves without an additional
Keller-specific pole filtration or projective extension.

## 4. What the matrices omit

For an actual finite flat algebra `A/R`, the trace matrix is the quadratic
shadow of multiplication.  Choose a basis

\[
e_0=1,e_1,\ldots,e_5
\]

and write `X_i` for multiplication by `e_i`.  Beyond the matrices above, one
needs all of the following polynomial identities:

\[
X_0=I,
\qquad
X_ie_0=e_i,
\qquad
X_ie_j=X_je_i,
\tag{4.1}
\]

\[
X_iX_j=\sum_k(X_ie_j)_kX_k,
\tag{4.2}
\]

\[
X_i^tT=TX_i,
\qquad
\operatorname{tr}(X_iX_j)=T_{ij}.
\tag{4.3}
\]

Equivalently, the symmetric cubic trace tensor

\[
C_{ijk}=\operatorname{Tr}(e_ie_je_k)
\]

must satisfy `C_{0ij}=T_{ij}`.  Recovering `e_i e_j` from
`T(e_ie_j,e_k)=C_{ijk}` gives the concrete integrality condition

\[
\operatorname{adj}(T)C_{ij,*}
\in(\det T)R^6
\tag{4.4}
\]

for every `i,j`, followed by unitality, associativity, and regular-trace
normalization.  None of (4.1)--(4.4) follows from a symmetric determinant or
normalization-module cokernel.  This multiplicative enhancement is the
narrowest current algebraic target.

The exact checker
[`trace_hostile_matrices.py`](../scripts/trace_hostile_matrices.py) verifies
both implicit equations, determinants, parametrized kernel vectors, all
listed ranks, the trace factorizations, and the norm-six vectors.  A hostile
entry perturbation breaks the `A6` certificate.  The independent
[`S6` topology checker](../scripts/s6_trace_curve_topology.py) exhausts all
transposition-valued representations of the explicit curve's affine
complement and proves that none is transitive; see the
[topology note](s6-trace-curve-topology.md).  The
[`A6` one-pair checker](../scripts/a6_one_pair_infinity.py) proves that the
explicit `A6` curve's complement is cyclic and, under additional infinity
hypotheses, excludes its entire collision-delta-at-most-two family; see the
[one-pair note](a6-one-pair-infinity.md).

## Primary sources

- Arnaud Beauville,
  [“Determinantal Hypersurfaces”](https://math.univ-cotedazur.fr/u/beauvill/pubs/det.pdf),
  *Michigan Math. J.* 48 (2000), 39--64.  Theorem B gives symmetric
  resolutions for self-dual ACM sheaves; Proposition 4.2 and Remark 4.4 treat
  theta characteristics and integral plane curves.
- Branko Curgus and Aad Dijksma,
  [“A proof of the main theorem on Bezoutians”](https://ems.press/journals/em/articles/12441),
  *Elemente der Mathematik* 58 (2003), 153--165.  The Bezoutian nullity theorem
  identifies the matrix nullity with the degree of the polynomial gcd,
  matching the normalization-fiber multiplicity formula used here.

The displayed matrices and their specialization to the degree-six trace
tables are exact constructions derived and checked in this repository.
