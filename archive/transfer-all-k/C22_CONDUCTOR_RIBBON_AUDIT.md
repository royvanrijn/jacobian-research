# C22 conductor-ribbon audit: counterexample to the canonical norm theorem

## Verdict

The proposed statement

> the completed degree-`k` square/cube factorization fiber is canonically
> isomorphic to the divided-power symmetric product of the split conductor
> ribbon

is false with the natural norm morphism.  It already fails for `k=2` over
every characteristic-zero field.  The divided-power space and the
factorization fiber both have length four at the double collision, but the
norm-generated subalgebra has length three.  Thus equality of ranks and
tangent dimensions does not provide an inverse.

This invalidates the independent proof in
[C22_DEFORMATION_AUDIT.md](C22_DEFORMATION_AUDIT.md) and the Boolean-norm
identification used in
[ALL_K_TRANSFER_BLOCK_THEOREM.md](ALL_K_TRANSFER_BLOCK_THEOREM.md).  It does
not disprove the numerical conjecture that the factorization block is finite
flat of rank `2^k`; it disproves the asserted ribbon model and removes the
current all-`k` proof of that conjecture.

Throughout the counterexample, `K` is a field in which `2` and `3` are
invertible.  In particular it applies in characteristic zero.

## 1. Precise rings and completion ideals

Let

\[
 U=Z^{3k}+u_1Z^{3k-1}+\cdots+u_{3k},\qquad
 V=Z^{2k}+v_1Z^{2k-1}+\cdots+v_{2k}
\]

and put

\[
 D_k=K[u_1,\ldots,u_{3k},v_1,\ldots,v_{2k}]
     /(\operatorname{coeff}_Z(U^2-V^3)).                       \tag{1}
\]

Because `2` is invertible, there is a unique monic

\[
 S=Z^k+s_1Z^{k-1}+\cdots+s_k                                  \tag{2}
\]

such that

\[
 V=S^2+A,\qquad \deg_Z A<k.                                  \tag{3}
\]

The coefficients `s_i` are obtained recursively from the coefficients of
degrees `2k-1,...,k` of `V`; the coefficients `a_0,...,a_(k-1)` of `A`
generate an ideal `J` in `D_k`.  The completed factorization block is

\[
 \mathfrak Z_k=\operatorname{Spf}\widehat{D_k}^{,J}.          \tag{4}
\]

Its base ring is

\[
 A_k=K[s_1,\ldots,s_k],                                      \tag{5}
\]

with the map to (4) defined by the extracted square root (2), not by a
separately chosen root cycle.

For the ribbon, put

\[
 E=K[z,\epsilon]/(\epsilon^2),\qquad E_0=K[z].                \tag{6}
\]

The affine divided-power symmetric product is

\[
 \mathfrak B_k=\operatorname{Spec}\Gamma_K^k(E).              \tag{7}
\]

Since `E` is `K`-flat, in characteristic zero its coordinate ring identifies
with

\[
 \Gamma_K^k(E)\simeq
 \left(K[r_1,\ldots,r_k,\epsilon_1,\ldots,\epsilon_k]
 / (\epsilon_1^2,\ldots,\epsilon_k^2)\right)^{S_k}.            \tag{8}
\]

The natural support map is induced by `E_0 -> E` and has base
`Gamma_K^k(E_0)=K[e_1(r),...,e_k(r)]`.  Its vertical ideal is nilpotent of
order at most `k+1`, so completing along the zero section changes no
Artinian calculation below.

## 2. Representing functors on Artinian local algebras

Let `R` be an Artinian local `K`-algebra.

The functor represented by (4) consists of pairs `(U_R,V_R)` of monic
polynomials of degrees `(3k,2k)` such that

\[
 U_R^2=V_R^3,                                                \tag{9}
\]

with the coefficients of `A_R=V_R-S_R^2` nilpotent, where `S_R` is the
unique polynomial extracted by (2)--(3).  Equivalently, these are continuous
maps from the `J`-adic ring in (4) to `R`.

The functor represented by (7) consists of multiplicative homogeneous
polynomial laws

\[
 n:E\longrightarrow R                                      \tag{10}
\]

of degree `k`.  Its natural support cycle is the restriction of `n` to
`E_0`.  This formulation remains meaningful when support points collide and
is the precise divided-power replacement for an unordered root list.

## 3. The norm-polynomial morphism

For a law (10), define

\[
 V_n(Z)=n((Z-z)^2+\epsilon),                                \tag{11}
\]

\[
 U_n(Z)=n((Z-z)^3+\tfrac32\epsilon(Z-z)).                   \tag{12}
\]

Multiplicativity and

\[
 ((Z-z)^3+\tfrac32\epsilon(Z-z))^2
 =((Z-z)^2+\epsilon)^3
\]

in `E[Z]` give `U_n^2=V_n^3`.  Hence (11)--(12) define a
coordinate-free norm morphism

\[
 N_k:\mathfrak B_k\longrightarrow\mathfrak Z_k.             \tag{13}
\]

This proves that the norm is well-defined at collisions.  It does **not**
prove that it is an isomorphism.

There is a first compatibility failure.  The support polynomial obtained by
restricting (10) to `E_0` need not equal the factorization base (2) extracted
from `V_n`.  For `k=2`, write

\[
 S_0=(Z-r_1)(Z-r_2)=Z^2+pZ+q,qquad
 \theta_0=\epsilon_1+\epsilon_2.
\]

Then direct coefficient comparison in (11) gives

\[
 S_{\mathrm{fact}}=Z^2+pZ+q+\frac12\theta_0.                \tag{14}
\]

Thus (13) is not a morphism over the natural base
`Gamma^2(A^1)`.  Any valid comparison has to use the sheared base (14).

## 4. Failure of an inverse at `k=2`

Impose the factorization collision base `S_fact=Z^2`.  In the ordered ribbon
cover this gives

\[
 r_2=-r_1,qquad r_1^2={\epsilon_1+\epsilon_2\over2},qquad
 \epsilon_1^2=\epsilon_2^2=0.                               \tag{15}
\]

Write

\[
 V_n=Z^4+XZ+Y.                                               \tag{16}
\]

Reduction using (15) gives

\[
 X=2r_1(\epsilon_2-\epsilon_1),\qquad
 Y={5\over2}\epsilon_1\epsilon_2,                            \tag{17}
\]

and consequently

\[
 X^2=XY=Y^2=0.                                               \tag{18}
\]

The norm-generated collision algebra therefore has basis `1,X,Y` and length
three.

The actual square/cube factorization fiber is

\[
 K[X,Y]/(X^3,XY,Y^2),                                      \tag{19}
\]

with basis `1,X,Y,X^2` and length four.  The class `X^2` is the missing
collision deformation.  Hence the pullback on coordinate rings induced by
(13) kills a nonzero element, and no inverse construction exists.

An explicit adversarial Artin point makes the failure visible.  Over

\[
 R=K[t]/(t^3)
\]

put

\[
 V=Z^4+tZ,qquad
 U=Z^6+\frac32tZ^3+\frac38t^2.                              \tag{20}
\]

Then `U^2=V^3`, the extracted base is exactly `S=Z^2`, and `X^2=t^2` is
nonzero.  Equations (17)--(18) show that (20) cannot factor through the
ribbon norm functor.

This example has the same reduced point and first-order tangent direction as
a ribbon deformation.  It is detected only at second order, which is why
length and tangent comparisons did not expose the error.

## 5. Exact checks for `k=1,2,3`

The script `scripts/verify_c22_ribbon_functor.py` performs the following
calculation over `QQ`.

1. Form the ordered Boolean ribbon and the exact norm polynomials.
2. Extract the factorization base from `V_n`, rather than imposing the
   forgetful support base.
3. Specialize that base to `S=Z^k`.
4. Compute the full invariant subspace of the ordered collision algebra.
5. Compute the subalgebra generated by the norm coefficients.
6. Independently eliminate `U` from `U^2=V^3` and compute the factorization
   fiber length.

The exact results are:

| `k` | ordered ribbon length | invariant ribbon length | norm-generated length | factorization length |
|---:|---:|---:|---:|---:|
| 1 | 2 | 2 | 2 | 2 |
| 2 | 8 | 4 | 3 | 4 |
| 3 | 48 | 8 | 4 | 8 |

Thus `k=1` is the only checked case in which the norm is an isomorphism.
At `k=3` it misses four of the eight invariant/factorization directions.
These are bounded exact computations, not an all-`k` theorem.

## 6. What remains of the rank-`2^k` argument

The divided-power ribbon itself is finite locally free of rank `2^k` over
its **natural** support base when `k!` is invertible.  Indeed, over ordered
roots it is the free Boolean algebra on `k` square-zero directions; Reynolds
splitting and descent give a finite projective module of generic rank `2^k`.
Over a field, Quillen--Suslin makes it free.

This does not prove the same statement for (4):

- the norm is not over the natural base, by (14);
- after the base shear, the norm is not surjective, by (18)--(20); and
- equal generic rank would not exclude torsion supported on the
  discriminant.

The existing straightening argument does not repair this.  Its distinct-root
Vandermonde determinant proves generic independence only.  The asserted
monic confluent straightening relations over the full base were not written
down, and the counterexample shows that the proposed norm generators do not
span even the `k=2` collision fiber.

For clarity, the abstract split-surjection/Nakayama step would be valid under
the following hypotheses.  Let `(A,m)` be local, let `M` be a finite
`A`-module, let `B` be finite free of rank `n`, and suppose `M -> B` is
surjective and `M` can be generated by at most `n` elements.  Projectivity of
`B` splits the surjection, so `M=B direct-sum K`; additivity of minimal
generator numbers and Nakayama force `K=0`.  In C22, the missing hypotheses
are precisely the all-base `n=2^k` generator bound and a genuine surjection.
The norm construction supplies neither.

The numerical rank-`2^k` statement for the factorization block remains:

- proved directly for `k=1,2`;
- supported by exact relative Groebner presentations for `k=3,4`; and
- supported at the maximally collided fiber for `k=5,6`.

It is not presently proved uniformly in `k`.

## 7. Characteristic ledger

Every characteristic-sensitive step is as follows.

1. Extracting `S` from `V=S^2+A` and defining the cusp jet uses `1/2`.
2. Solving the one-root coefficient equations gives
   `b=3d/2`; the next equation is `-3d^2/4=0`.  Inferring `d^2=0`
   therefore requires both `2` and `3` invertible.
3. Identifying divided powers with symmetric tensors and using a Reynolds
   splitting requires `k!` invertible.  The divided-power functor itself does
   not require this identification.
4. The triangular coefficient recovery in the rejected norm proof divides
   by `1,2,...,k`.
5. Its Laurent straightening also used nonvanishing of
   `binom(3/2,d)`; this is proof-specific and cannot now support the theorem.
6. The separate affine-difference Wronskian proof divides by `6k` and
   `1-6k`; over a field that particular proof excludes characteristics
   dividing `6k(6k-1)`.
7. The explicit counterexample (20) requires only that `2` and `3` be
   invertible.

Consequently characteristic zero is a safe common hypothesis, but statements
in positive characteristic must name which of the seven operations they use.

## 8. Corrected open problem

The canonical ribbon theorem should not be repaired by changing only the
word “canonical.”  The second-order class (20) proves that the underlying
norm functor is too small.  A viable replacement must enlarge the local cusp
object so that its collision moduli retain higher principal-part classes such
as `X^2`.

The immediate all-`k` problem is:

> Prove directly, in the completed factorization ring (4), a monic
> straightening basis of size `2^k`, or find a value of `k` at which the
> rank or flatness fails.

Until that is done, C21 and the allocation-intersection statements that use
the all-`k` Boolean relation have a material unresolved dependency.

The direct elimination identity, exact bases through `k=4`, and bounded
counterexample search are now recorded in
[DIRECT_TRANSFER_BASIS.md](DIRECT_TRANSFER_BASIS.md).
