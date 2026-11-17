# Three-torsion independence from every field unramified at three

The three-torsion representation of the fixed elliptic curve is
surjective. More particularly, its inertia at three normally generates
the whole Galois group. The three-torsion field is therefore disjoint
from every finite Galois extension unramified at three, including the
prime-to-three torsion and rational division fields used in the Kummer
arguments. This gives exact compatible Frobenius densities; it does not
control reduction modulo the square of a varying residue prime.

## 1. The full three-torsion group

Use the fixed short model
\[
 E:y^2=X^3+AX+B,\quad A=-10758382,\quad B=6600759744,
 \qquad K_3=\mathbb Q(E[3]).
\]
Its discriminant is nonzero and prime to three. The four x-coordinates
of the nonzero three-torsion pairs are the roots of
\[
 \psi_3(X)=3X^4+6AX^2+12BX-A^2.
\]
Indeed the duplication formula gives
\[
 X(2P)=X(P)\quad\Longleftrightarrow\quad
 (3X^2+A)^2-12X(X^3+AX+B)=-\psi_3(X)=0.
\]
For these points `y!=0`; duplication cannot give `2P=P` for a
nonidentity point, so it gives `2P=-P`. Conversely every nonzero
three-torsion point satisfies this equation.

For an integral monic polynomial with the same splitting field, put
\[
 g(Z)=27\psi_3(Z/3)
      =Z^4-193650876Z^2+712882052352Z-3125055147963948.
                                                               \tag{1.1}
\]
Its reductions have the exact factorizations
\[
 g(Z)\bmod11=Z^4-Z^2-3Z-1\quad\text{irreducible},
\]
\[
 g(Z)\bmod13=(Z-6)(Z^3+6Z^2-6Z+5),             \tag{1.2}
\]
where the cubic is irreducible. Both reductions are squarefree. The
checker proves the quartic assertion by excluding every monic linear
and quadratic divisor and the cubic assertion by excluding all thirteen
roots; these are complete finite irreducibility criteria.

**Theorem ETD1 (full mod-three image).**
\[
             \operatorname{Gal}(K_3/\mathbb Q)
                        =\operatorname{GL}_2(\mathbb F_3).
                                                               \tag{1.3}
\]

**Proof.** Dedekind's Frobenius factorization theorem applied to the
monic integral polynomial (1.1) supplies a four-cycle and a three-cycle
in its root group. The group is transitive and has order divisible by
twelve. As a subgroup of `S_4`, its order is twelve or twenty-four.
The only index-two subgroup is `A_4`: any homomorphism from `S_4` to
`C_2` sends its conjugate generating transpositions to the same value.
A four-cycle is odd, so the group is `S_4`.

The four roots are the four one-dimensional subspaces of `E[3]`:
each consists of the two points `P,-P` with the same x-coordinate.
Thus the projective image is `PGL_2(F_3)=S_4`. Its full linear image
has order twenty-four or forty-eight. In the smaller case it would
be a subgroup isomorphic to `S_4` missing the scalar `-I`. But in
characteristic three the only nonidentity involution of determinant
one is `-I`: an involution is diagonalizable with eigenvalues in
`{1,-1}`.
The three involutions in a Klein four subgroup of such an embedded
`S_4` would consequently all have determinant minus one. The product
of two would have determinant one and be the third nonidentity
involution, a contradiction. The linear image is full. \(\square\)

The exact factorization import is Milne, *Algebraic Number Theory*,
version 3.00, Corollary 8.22 and Theorem 8.23, printed p.138
(PDF p.146), in the [retained primary](sources/milne-ant300.pdf).
Using (1.1) rather than a nonmonic polynomial verifies its integral
monic hypothesis directly.

## 2. No nontrivial unramified-at-three Galois subextension

**Theorem ETD2 (universal unramified disjointness).** For every
finite Galois extension `C/Q` unramified at three,
\[
                  \boxed{K_3\cap C=\mathbb Q.}     \tag{2.1}
\]
In particular `K_3` and `C` are linearly disjoint over `Q`.

**Proof.** The determinant on `E[3]` is the cyclotomic character,
by the Galois-module identity
\[
                       \bigwedge^2E[3]\simeq\mu_3.
\]
In particular `Q(zeta_3) subset K_3`. This quadratic extension is
ramified at three: `zeta_3-1` satisfies `T^2+3T+3`, with valuation
one half. Consequently the inertia subgroup `I_3` at three in
`Gal(K_3/Q)` contains an element of determinant minus one.

On `GL_2(F_3)` the determinant equals the sign of the action on the
four projective lines. For instance `SL_2(F_3)` maps to an index-two
subgroup of `S_4`, hence `A_4`. Let `N` be the normal closure of
`I_3`. Its projective image is normal in `S_4` and contains an odd
permutation, so it is all `S_4`. To see the last elementary fact,
a normal subgroup containing a transposition contains all of them.
If it contains a four-cycle `c=(1234)`, its commutator with `(12)`
is a three-cycle, so it contains `A_4` and also the odd element.
The smaller-lift obstruction in the proof of ETD1 now gives
`N=GL_2(F_3)`.

The intersection `K_3 cap C` is Galois over `Q` and unramified at
three. Its defining normal subgroup therefore contains `I_3`, hence
its normal closure, which is the whole group. This proves (2.1).
The product description of the compositum follows because both
extensions are Galois. \(\square\)

The determinant identity is stated in Milne, *Elliptic Curves*
(2006), Chapter V, printed p.224 (PDF p.232), the paragraph after
Conjecture 8.5: [primary book](https://www.jmilne.org/math/Books/ectext6.pdf),
[retained copy](sources/milne-elliptic-curves-2006.pdf).
It is the stated Weil-pairing identity that is used, not the
conjecture discussed on that page.

## 3. Compatible exact prime densities

Let
\[
 \Sigma_3=\{g\in\operatorname{GL}_2(\mathbb F_3):
                     \det g=-1,\ \det(I-g)\ne0\}.
\]
There are twelve such matrices. Indeed twenty-four matrices have
determinant minus one. Those with an eigenvalue one have the distinct
eigenvalues `1,-1` and form one conjugacy class, represented by
`diag(1,-1)`. Its centralizer has four elements, so that class has
twelve elements. The remaining twelve constitute `Sigma_3`.

**Corollary ETD3 (exact joint density).** Let `C/Q` be finite Galois
and unramified at three, let `G_C=Gal(C/Q)`, and let
`Omega subset G_C` be conjugacy-stable. The good rational primes
unramified in the relevant fields satisfying
\[
 \operatorname{Frob}_p(C/\mathbb Q)\in\Omega,\qquad
             p\equiv2\pmod3,\qquad3\nmid\#E(\mathbb F_p)
\]
have exact natural prime density
\[
                  \boxed{\frac{|\Omega|}{4|G_C|}.}            \tag{3.1}
\]
With `C=Q`, the density is `1/4`. With complete splitting in a
fixed `C`, it is `1/(4[C:Q])`.

**Proof.** At a good prime different from three, reduction identifies
the prime-to-residue-characteristic three-torsion, equivariantly for
Frobenius. Thus its determinant is `p mod3`, and `E(F_p)` has order
divisible by three exactly when that Frobenius fixes a nonzero
three-torsion point. The latter equivalence is also the elementary
order-three criterion for a finite abelian group.

ETD1--ETD2 give the product Galois group
`GL_2(F_3) times G_C`. The required set is `Sigma_3 times Omega`.
It is conjugacy-stable, and Chebotarev gives
`12|Omega|/(48|G_C|)`, proving (3.1). Omitting finitely many primes
does not change it. \(\square\)

The good-reduction torsion interface is Milne, *Abelian Varieties*
(January 2, 2022), Proposition 20.7 and Corollary 20.8, printed/PDF
p.44, in [the retained source](sources/milne-abelian-varieties-2022.pdf).
The density theorem is Milne ANT300, Definition 8.30 and Theorem
8.31, printed p.140 (PDF p.148), in the source cited above.
For example `C` can be any fixed finite Galois compositum of torsion
and division fields of rational points with multipliers prime to
three: good reduction at three and the finite-etale division theorem
give the required unramifiedness.

The disjointness is stronger than just disjointness from `Q(zeta_3)`.
It permits arbitrary conjugacy-compatible data in `C` while imposing
the displayed three-torsion conditions independently. It does not
assert independence of the reductions at different residue primes,
nor does it control the varying formal depths `e_p,k_p`.

## 4. Exact finite controls and assurance

Run from the repository root:
```sh
python3 research/arithmetic/verify_elliptic_exception_three_torsion.py
```
The checker verifies the duplication identity, the monic transform,
complete irreducibility controls at eleven and thirteen, and all
forty-eight invertible matrices and their projective actions. It
exhausts the finite group facts used in the lift, sign and density
arguments. The
[receipt](elliptic_exception_three_torsion_checks.json) pins the
source and checker. Frobenius factorization, the Weil pairing,
local inertia and Chebotarev are written/imported arguments with
the exact interfaces above; their deep proofs are not certified by
the finite controls. No formal verification, external review or
independent whole-proof replay is claimed.
