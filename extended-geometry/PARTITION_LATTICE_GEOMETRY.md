# Weighted partition-lattice geometry of the exceptional locus

Work over a field of characteristic zero.  This note upgrades the
set-theoretic common-coarsening statement to a generic scheme-theoretic
calculation on every **minimal** intersection boundary of two atomic
components.  It also records the complete census through degree twenty and
separates component intersections from self-intersections of normalization
sheets.

The result is a weighted partition complex rather than a normal-crossings
complex.  Its primitive weight is the length-two square/cube transfer block.
The component counts, interval structure, and binomial branch profiles used
here are recorded separately in the
[enumerative skeleton](NONSURJECTIVE_ENUMERATIVE_GEOMETRY.md); that note
deliberately stops before the scheme-theoretic weights computed below.

## 1. Atomic vertices and the collision poset

The vertices in degree `N` are

\[
 \mathcal C_{a,b},\qquad N=2a+3b,
\]

with atomic partition `2^a3^b` and dimension `a+b-1`.  The boundary poset of
one vertex is the coarsening poset of its atomic partition.  By the exact
closure theorem,

\[
 \mathcal E_\nu\subset\mathcal C_{a,b}
 \quad\Longleftrightarrow\quad
 2^a3^b\preceq\nu.
\]

Set-theoretic intersections are therefore common upper ideals in the
coarsening order.  Scheme structure requires one additional decoration: an
allocation of every collision part `m_rho` by

\[
 2i_\rho+3j_\rho=m_\rho.
\]

Distinct allocations are distinct points of the smooth component
normalization, even when they lie above the same collision partition.

## 2. Minimal common coarsening

Let `C_(a,b)` and `C_(c,d)` be distinct degree-`N` components.  Equal degree
gives

\[
 2(a-c)+3(b-d)=0,
\]

so there is a nonzero integer `k` with

\[
 a-c=3k,\qquad b-d=-2k.
\]

Put `t=|k|`.  After common double and triple atoms are cancelled, the only
primitive equality is

\[
 3\cdot2=2\cdot3=6.
\]

Consequently the unique minimal common integer-partition coarsening is

\[
 \boxed{\nu_{\min}
 =6^t3^{\min(b,d)}2^{\min(a,c)}.}                 \tag{1}
\]

This uniqueness is stronger than merely finding one collision: any common
block cancels its common twos and threes, and every remaining relation splits
into primitive sixes.  Keeping the sixes separate is precisely minimality in
the coarsening order.

## 3. The universal collision weight

At one primitive sixfold transfer, cancel the common factors and translate
the collision root.  The two allocations compare a monic cubic `Q` and a
monic quadratic `T` through

\[
 Q^2=T^3.
\]

Write

\[
 Q=Z^3+uZ+v,\qquad T=Z^2+\alpha Z+\beta.
\]

Coefficient comparison has Groebner basis

\[
 \alpha,\qquad v,\qquad2u-3\beta,\qquad\beta^2.
\]

Thus one block is a smooth root-position coordinate times

\[
 \boxed{Z_1=k[\epsilon]/(\epsilon^2),\qquad
        \operatorname{length}Z_1=2.}             \tag{2}
\]

The number `2` is the primitive weight of the partition complex.

## 4. Affine-difference rigidity

The leading seed coefficient recovers `D(M)`.  Equality of normalized seeds
therefore makes their monic omitted polynomials differ by at most a constant
and a linear term.  It remains to check that this does not enlarge an
intersection.  Let

\[
 M=Q^2R^3,\qquad P=S^2T^3,qquad M-P=\lambda W+\mu.
\]

The Wronskian

\[
 PM'-MP'=P(M-P)'-(M-P)P'
\]

is divisible by

\[
 QR^2ST^2.
\]

Its divisor has degree

\[
 a+2b+c+2d=N+{b+d\over2}>N,                     \tag{3}
\]

because two distinct equal-degree atomic types cannot both have `b=d=0`.
The right side has degree at most `N`, so the Wronskian vanishes.  Its two
leading coefficients give `(1-N)lambda=0` and `-N mu=0`.  Hence

\[
 \boxed{M-P\in\langle1,W\rangle\Longrightarrow M=P.}            \tag{4}
\]

This argument is valid over rings with nilpotents: it uses monicity and the
invertibility of nonzero integers, not cancellation by a possibly
zero-divisorial coefficient.

## 5. Generic minimal-intersection theorem

Complete on the ordered-root chart at the geometric generic point of
`E_(nu_min)`.  Formal Hensel factorization separates its distinct roots.
The common double and triple atoms give smooth diagonal factors, while the
`t` distinct sixfold roots give `t` independent copies of (2).  Equation
(4) identifies the actual normalized-seed correspondence with this strong
factorization correspondence.  Finally, the splitting-and-rescaling
construction from the contact-strata theorem gives a point where `Phi=0` has
nonzero derivative in the common scaling/root-position direction (the
derivative is controlled by the unit `D`).  Thus on a dense open, `Phi=0`
eliminates one smooth support parameter and does not alter the transverse
Artin algebra.

Therefore

\[
 \boxed{
 D_{(a,b),(c,d)}^{\mathrm{tr}}
 \simeq
 k[\epsilon_1,\ldots,\epsilon_t]/
 (\epsilon_1^2,\ldots,\epsilon_t^2).}            \tag{5}
\]

In particular,

\[
\boxed{I((a,b),(c,d))=2^t,
\qquad t={|a-c|\over3}={|b-d|\over2}.}          \tag{6}
\]

Here `I` means the transverse length of the normalization-branch
correspondence at the geometric generic point of the minimal stratum.  The
component intersections are generally excess-dimensional, so (6) is not a
classical proper-intersection number in the ambient seed space.

The transverse embedding dimension is `t`, its tangent-cone Hilbert vector
is

\[
 \left({t\choose0},{t\choose1},\ldots,{t\choose t}\right),      \tag{7}
\]

and the full intersection tangent dimension is

\[
 \dim E_{\nu_{\min}}+t=\ell(\nu_{\min})-1+t.                   \tag{8}
\]

Thus multiplicity factors over independent primitive collision clusters:

\[
 I=\prod_{\rho=1}^t\operatorname{length}Z_1=2^t.                \tag{9}
\]

### Weighted interval series

As a corollary, apply the transverse length `2^t` to the unweighted interval
series (12) of the [enumerative skeleton](NONSURJECTIVE_ENUMERATIVE_GEOMETRY.md).
In the notation of that series, `y` records the component-index gap `k`
(equal to the transfer-block count called `t` above), while the exponent of
the series variable `t` records the dimension of the top contact stratum
plus one.  Retaining the labels that mark these as refinements of (12), the
scheme-weighted refinement is

\[
 \boxed{
 \sum_{N,k,i}2^k x^N t^{\dim\mathcal E_{\nu_i}+1}y^k
 ={1\over(1-tx^2)(1-tx^3)}
   {2ytx^6\over1-2ytx^6}.}                                    \tag{12a}
\]

More finely, (7) replaces the scalar weight `2^k` by the tangent-cone
Hilbert polynomial `(1+z)^k`, giving

\[
 \boxed{
 \sum_{N,k,i}(1+z)^k x^N t^{\dim\mathcal E_{\nu_i}+1}y^k
 ={1\over(1-tx^2)(1-tx^3)}
   {y(1+z)tx^6\over1-y(1+z)tx^6}.}                            \tag{12b}
\]

This is a generic statement on the exact minimal stratum.  If primitive
sixfold clusters collide with one another, the relevant higher transfer
block is not literally a tensor product of dual numbers; its multiplication
can change even when its length is preserved in computed low orders.

## 6. Normalization self-pairs

A component has an off-diagonal normalization pair as soon as `a>=3` and
`b>=2`.  Its first such boundary is

\[
 \boxed{6^2 3^{b-2}2^{a-3}.}                    \tag{10}
\]

The two sheets exchange the allocations `(3,0)` and `(0,2)` at the two
sixfold roots.  Hence (5) applies with two independent blocks: the transverse
length is four, transverse embedding dimension is two, and tangent-cone
Hilbert vector is `(1,2,1)`.  This recovers the degree-twelve `(6,6)` theorem
and explains why the same local singularity repeats after adjoining untouched
atomic roots.

The word *self* below means this off-diagonal normalization correspondence;
it is not the formal intersection of a scheme with itself.

## 7. Complete census through degree twenty

The table records every unordered pair of distinct components and every
first normalization self-pair.  `dim E` is the reduced support dimension,
`t` is the number of primitive transfer clusters, `I` is transverse length,
and `edim` is transverse embedding dimension.  Degrees with no nontrivial
row are omitted.

| N | pair | minimal coarsening | dim E | t | I | edim | Hilbert | factors? |
|---:|---|---|---:|---:|---:|---:|---|---|
| 6 | `C_(3,0) ∩ C_(0,2)` | `6` | 0 | 1 | 2 | 1 | `(1,1)` | yes |
| 8 | `C_(4,0) ∩ C_(1,2)` | `6 2` | 1 | 1 | 2 | 1 | `(1,1)` | yes |
| 9 | `C_(3,1) ∩ C_(0,3)` | `6 3` | 1 | 1 | 2 | 1 | `(1,1)` | yes |
| 10 | `C_(5,0) ∩ C_(2,2)` | `6 2^2` | 2 | 1 | 2 | 1 | `(1,1)` | yes |
| 11 | `C_(4,1) ∩ C_(1,3)` | `6 3 2` | 2 | 1 | 2 | 1 | `(1,1)` | yes |
| 12 | `C_(6,0) ∩ C_(3,2)` | `6 2^3` | 3 | 1 | 2 | 1 | `(1,1)` | yes |
| 12 | `C_(6,0) ∩ C_(0,4)` | `6^2` | 1 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 12 | `C_(3,2) ∩ C_(0,4)` | `6 3^2` | 2 | 1 | 2 | 1 | `(1,1)` | yes |
| 12 | `C_(3,2) self` | `6^2` | 1 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 13 | `C_(5,1) ∩ C_(2,3)` | `6 3 2^2` | 3 | 1 | 2 | 1 | `(1,1)` | yes |
| 14 | `C_(7,0) ∩ C_(4,2)` | `6 2^4` | 4 | 1 | 2 | 1 | `(1,1)` | yes |
| 14 | `C_(7,0) ∩ C_(1,4)` | `6^2 2` | 2 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 14 | `C_(4,2) ∩ C_(1,4)` | `6 3^2 2` | 3 | 1 | 2 | 1 | `(1,1)` | yes |
| 14 | `C_(4,2) self` | `6^2 2` | 2 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 15 | `C_(6,1) ∩ C_(3,3)` | `6 3 2^3` | 4 | 1 | 2 | 1 | `(1,1)` | yes |
| 15 | `C_(6,1) ∩ C_(0,5)` | `6^2 3` | 2 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 15 | `C_(3,3) ∩ C_(0,5)` | `6 3^3` | 3 | 1 | 2 | 1 | `(1,1)` | yes |
| 15 | `C_(3,3) self` | `6^2 3` | 2 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 16 | `C_(8,0) ∩ C_(5,2)` | `6 2^5` | 5 | 1 | 2 | 1 | `(1,1)` | yes |
| 16 | `C_(8,0) ∩ C_(2,4)` | `6^2 2^2` | 3 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 16 | `C_(5,2) ∩ C_(2,4)` | `6 3^2 2^2` | 4 | 1 | 2 | 1 | `(1,1)` | yes |
| 16 | `C_(5,2) self` | `6^2 2^2` | 3 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 17 | `C_(7,1) ∩ C_(4,3)` | `6 3 2^4` | 5 | 1 | 2 | 1 | `(1,1)` | yes |
| 17 | `C_(7,1) ∩ C_(1,5)` | `6^2 3 2` | 3 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 17 | `C_(4,3) ∩ C_(1,5)` | `6 3^3 2` | 4 | 1 | 2 | 1 | `(1,1)` | yes |
| 17 | `C_(4,3) self` | `6^2 3 2` | 3 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 18 | `C_(9,0) ∩ C_(6,2)` | `6 2^6` | 6 | 1 | 2 | 1 | `(1,1)` | yes |
| 18 | `C_(9,0) ∩ C_(3,4)` | `6^2 2^3` | 4 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 18 | `C_(9,0) ∩ C_(0,6)` | `6^3` | 2 | 3 | 8 | 3 | `(1,3,3,1)` | yes |
| 18 | `C_(6,2) ∩ C_(3,4)` | `6 3^2 2^3` | 5 | 1 | 2 | 1 | `(1,1)` | yes |
| 18 | `C_(6,2) ∩ C_(0,6)` | `6^2 3^2` | 3 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 18 | `C_(3,4) ∩ C_(0,6)` | `6 3^4` | 4 | 1 | 2 | 1 | `(1,1)` | yes |
| 18 | `C_(6,2) self` | `6^2 2^3` | 4 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 18 | `C_(3,4) self` | `6^2 3^2` | 3 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 19 | `C_(8,1) ∩ C_(5,3)` | `6 3 2^5` | 6 | 1 | 2 | 1 | `(1,1)` | yes |
| 19 | `C_(8,1) ∩ C_(2,5)` | `6^2 3 2^2` | 4 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 19 | `C_(5,3) ∩ C_(2,5)` | `6 3^3 2^2` | 5 | 1 | 2 | 1 | `(1,1)` | yes |
| 19 | `C_(5,3) self` | `6^2 3 2^2` | 4 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 20 | `C_(10,0) ∩ C_(7,2)` | `6 2^7` | 7 | 1 | 2 | 1 | `(1,1)` | yes |
| 20 | `C_(10,0) ∩ C_(4,4)` | `6^2 2^4` | 5 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 20 | `C_(10,0) ∩ C_(1,6)` | `6^3 2` | 3 | 3 | 8 | 3 | `(1,3,3,1)` | yes |
| 20 | `C_(7,2) ∩ C_(4,4)` | `6 3^2 2^4` | 6 | 1 | 2 | 1 | `(1,1)` | yes |
| 20 | `C_(7,2) ∩ C_(1,6)` | `6^2 3^2 2` | 4 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 20 | `C_(4,4) ∩ C_(1,6)` | `6 3^4 2` | 5 | 1 | 2 | 1 | `(1,1)` | yes |
| 20 | `C_(7,2) self` | `6^2 2^4` | 5 | 2 | 4 | 2 | `(1,2,1)` | yes |
| 20 | `C_(4,4) self` | `6^2 3^2 2` | 4 | 2 | 4 | 2 | `(1,2,1)` | yes |

There are `36` distinct-component rows and `10` first self-pair rows.  The
only transverse lengths through degree twenty are `2,4,8`.

## 8. Relation to the degree-eighteen triple point

The row `C_(3,4) self` has generic reduced support `E_(6,6,3,3)` and length
four.  Specializing the two common triple roots to another sixfold root gives
`E_(6,6,6)`, where the component has three normalization sheets.  The known
ordered triple intersection has transverse algebra

\[
 k[\epsilon_1,\epsilon_2,\epsilon_3]/
 (\epsilon_1^2,\epsilon_2^2,\epsilon_3^2)
\]

of length eight.  That is a genuinely higher face of the weighted complex:
its cubic socle class is invisible in every pairwise slice.  The
excess-conormal theorem in the degree-eighteen note proves that the pairwise
data recover the entire three-sheet completed component ring on a nonempty
open.

## 9. Proof boundary and next theorem

Equations (1)--(9) settle every minimal intersection of two distinct atomic
components in all degrees, not only through degree twenty.  The bounded
census additionally records the first self-pair of every component in the
requested range.

This does **not** yet assign a universal algebra to every nonminimal
coarsening or every multi-sheet face.  The remaining uniform singularity
theorem needs two extra layers:

1. classify the higher local transfer block `Z_k` when several primitive
   sixfold transfers collide at one root; and
2. prove effective Cech descent for three or more normalization sheets, so
   pairwise conductor congruences recover the completed coefficient ring.

The partition lattice supplies the correct indexing.  The transfer-block
algebras supply its weights, and the conductor equalizers supply its
higher-face compatibility.

## Executable census

Run

```bash
python scripts/verify_exceptional_partition_lattice.py
python scripts/verify_exceptional_partition_lattice.py --json
```

The checker exhausts all full-contact partitions through degree twenty,
verifies the unique minimal common coarsening for every pair, checks the
universal dual-number Groebner basis and Wronskian degree gap, audits every
self-pair allocation, and emits the complete machine-readable census.
