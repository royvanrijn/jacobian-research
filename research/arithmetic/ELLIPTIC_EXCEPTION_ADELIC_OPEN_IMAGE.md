# Finite-exception torsion images and the meaning of the class-group density gap

The explicit elliptic curve is non-CM. Serre's open-image theorem
therefore supplies the full torsion-image hypotheses used in the
strict Selmer converse outside one fixed finite set of primes.
Together with the retained Kummer theorem, it also gives the full
rank-two affine torsion/division image at every sufficiently large
auxiliary prime, after any fixed number-field extension.

The resulting class-group criterion is unconditional outside finite
exceptions. It shows that the proposed second-copy density condition
asks simultaneously for control of deep formal images and of primes
with nonzero `Sha[p]`; it is not a way to bypass either question.

## 1. The primary torsion-image input

Retain the fixed curve
\[
 E:y^2=X^3-10758382X+6600759744.
\]
[EKX1](ELLIPTIC_EXCEPTION_KUMMER_BRANCH_EXCLUSION.md) proves
`End(E_Qbar)=Z`: its `j` invariant is nonintegral, with valuation
minus one at `1889`, whereas a CM `j` invariant is an algebraic
integer. This is an actual non-CM proof, not numerical evidence.

The precise imported statement is Serre, *Proprietes galoisiennes
des points d'ordre fini des courbes elliptiques*, Inventiones
mathematicae **15** (1972), 259--331, section 4.4, Theorem 3,
printed p.299, and Corollary 1, printed p.300. For a non-CM
elliptic curve over a number field `L`, its adelic torsion image
is open in `GL_2(Zhat)` and contains the entire `ell`-adic factor
outside finitely many `ell`. Section 4.2, Theorem 2, printed
p.294, also states full mod-`ell` image outside finitely many primes.

[Author-hosted primary PDF](https://www.college-de-france.fr/media/jean-pierre-serre/UPL5874918517843398173_Serre_proprie_te_s_galoisiennes_des_courbes_elliptiques.pdf),
[retained PDF](sources/serre-1972-galois-torsion.pdf),
[exact-version metadata](sources/serre-1972-galois-torsion.json).
The retained file has one added cover page: printed pages
259,260,294,299,300 are PDF pages 2,3,37,42,43. The scanned
theorem and corollary pages were rendered and read directly.
Serre's proof is imported, not reproved here.

**Theorem EOI1 (full torsion factors outside fixed exceptions).**
For every fixed number field `L`, there is an integer `N_L>=1`
such that the image on all torsion of `E/L` contains
\[
 \ker\!\left(\operatorname{GL}_2(\widehat{\mathbb Z})
              \longrightarrow\operatorname{GL}_2(\mathbb Z/N_L)
       \right).
                                                               \tag{1.1}
\]
In particular it contains the whole factor `GL_2(Z_ell)` for
every `ell` not dividing `N_L`, with identity in all other
factors. Its image on `E[ell^n]` is `GL_2(Z/ell^n)` for every
such `ell` and every `n>=1`.

**Proof.** Non-CM is unchanged by base extension. Apply the stated
open-image theorem to `E/L`. The kernels of reduction modulo
positive integers form a neighbourhood basis of the identity in
`GL_2(Zhat)`, giving (1.1). For `ell` not dividing `N_L`, the
whole `ell`-adic factor lies in that kernel. Reduction proves the
last assertion. \(\square\)

No value of `N_L` or explicit exceptional-prime list is computed.
This theorem does not assert that the complete adelic image is
all of `GL_2(Zhat)`; the rational two-torsion already prevents that.

## 2. The strict Selmer converse applies at every sufficiently large prime

Use the notation of
[ECS1--ECS3](ELLIPTIC_EXCEPTION_CLASS_GROUP_SELMER_CONVERSE.md):
\[
 K_p=\mathbb Q(E[p]),\qquad V_p=E[p],\qquad
 \mathcal C_p=\operatorname{Hom}_{G_p}(\operatorname{Cl}(K_p)/p,V_p),
\]
and let `e_p` be the depth of the actual full rational formal image.

**Corollary EOI2 (an unconditional finite-exception converse).**
There is a fixed finite set `Sigma` such that, for every prime
`p` outside `Sigma`, ECS1--ECS3 apply and
\[
 \boxed{\operatorname{Cl}(K_p)/p\twoheadrightarrow_{G_p}V_p^2
       \quad\Longleftrightarrow\quad
       e_p>1\ \text{or}\ \Sha(E/\mathbb Q)[p]\ne0.}              \tag{2.1}
\]
The class-group map space is exactly the strict Selmer group
`ker(Sel_p(E/Q)->H^1(Q_p,V_p))`, with its exact rational/Sha
sequence from ECS2.

Consequently the following density-zero assertions are equivalent:

1. Primes with the actual equivariant two-copy quotient in (2.1)
   have natural density zero.
2. Primes with `e_p>1` have natural density zero, and primes with
   `Sha(E/Q)[p]!=0` have natural density zero.

**Proof.** Let `Sigma` include primes dividing `N_Q`, primes below
seven, bad reduction primes, and prime divisors of the finite
bad-local quotient orders `m_ell` in ECS1. EOI1 makes
`G_p=GL_2(F_p)` outside this set; every other ECS hypothesis is
included by construction. Apply ECS1--ECS3.

Outside `Sigma`, the prime set in assertion 1 is exactly the union
of the two sets in assertion 2. For any two sets of primes, their
union has density zero if and only if each has density zero: bound
each counting function by that of the union, and the union by their
sum. Finite changes have no effect. \(\square\)

There is no assumption here that `Sha` is finite. If one separately
assumes its finiteness, its nonzero `p`-torsion occurs at only finitely
many primes, so the two-copy and deep-image conditions agree outside
a finite set. Without that additional input, the second-copy density
condition is stronger than the deep-image density condition alone.

This sharpens the interpretation of
[EWC3](ELLIPTIC_EXCEPTION_WIEFERICH_CLASS_GROUP_BRIDGE.md): a
density-zero theorem for the extra class-group copy would indeed give
the stated weighted first moments, but the exact converse shows the
two arithmetic contributions that such a theorem must control.
Neither density is proved here, and first moments still leave the
conditional correlation and recurrence obligations open.

## 3. The full affine image at a large auxiliary prime

Fix an abstract free basis `P_1,P_2` of `E(Q)/torsion`, with rational
lifts, and let `L` be any fixed number field. For `ell` and `n>=1`,
write
\[
 L_n=L(E[\ell^n]),\qquad
 M_n=L(E[\ell^n],[\ell^n]^{-1}P_1,[\ell^n]^{-1}P_2).
\]

**Theorem EOI3 (full rank-two affine division image).** For every
fixed `L`, all sufficiently large `ell` have, simultaneously for
every `n>=1`, the full affine image
\[
 \boxed{\operatorname{Gal}(M_n/L)
     = E[\ell^n]^2\rtimes\operatorname{GL}_2(\mathbb Z/\ell^n)}.
                                                               \tag{3.1}
\]
Here the equality means the image in the natural affine action
after choosing compatible division points. Passing to the inverse
limit gives `(T_ell E)^2 rtimes GL_2(Z_ell)`.

More generally there is a fixed finite set of primes `Sigma_L`
such that, for every positive integer `m` supported outside it,
the corresponding image is
`E[m]^2 rtimes GL_2(Z/m)`. Thus finitely many such auxiliary
primary factors can be used simultaneously with their full product
image; the same fixed exceptional set works for every `m`.

**Proof.** Since `End(E_Qbar)=Z`, independence of the two basis
points makes `(P_1,P_2)` nondegenerate on `E^2`, by the argument of
EKX1; this remains true over `L`. Apply the retained uniform Kummer
index theorem to this fixed abelian variety and point over `L`:
Corvaja--Demeio--Javanpeykar--Lombardo--Zannier, *Distribution of
rational points on ramified covers of abelian varieties*, Compositio
Mathematica **158** (2022), Theorem 6.9 and Remark 6.10, printed
pp.2144--2145/PDF pp.37--38:
[primary DOI](https://doi.org/10.1112/S0010437X22007746),
[retained PDF](sources/cdjlz-2022-ramified-covers.pdf).

It bounds the index of `Gal(M_n/L_n)` in `E[ell^n]^2` by a constant
`c_L` independent of `ell,n`. That index is a power of `ell`.
For `ell>c_L` it must therefore be one, at every `n`. EOI1 makes
the torsion projection full outside another fixed finite set.
An affine subgroup with full translation kernel and full linear
projection is the entire affine group: over each linear element,
one lift and all translations give every affine lift. This proves
(3.1). Compatible choices and compactness give the inverse limit.
For the general `m`, the same Kummer index divides `m^4` and is
at most `c_L`. If every prime factor of `m` exceeds `c_L`, the
index must be one. Include the primes dividing `N_L` and the
primes at most `c_L` in `Sigma_L`, then use (1.1) for the full
linear projection modulo `m`. The same affine-group argument
proves the simultaneous statement.
\(\square\)

The theorem also applies after adjoining any fixed finite Galois
data to `L`; the sufficiently-large threshold can change with that
extension. It supplies full images for fixed auxiliary primes,
not independence of reductions at varying residue primes. It
does not bound a growing family of coefficient characters or
the large common-quotient tail in the recurrence problem.

These are written deductions using the stated primary theorems.
No numerical exceptional threshold, artificial finite checker,
formal verification or external whole-proof review is claimed.

## 4. A single open affine image, including the exceptional factors

**Corollary EOI4 (adelic rank-two affine openness).** Choose
compatible division points of the two basis points at every level.
For every fixed number field `L`, the resulting adelic affine image
is an open subgroup of
\[
 \mathcal A=(\widehat{\mathbb Z}^{\,2})^2
                         \rtimes\operatorname{GL}_2(\widehat{\mathbb Z})
          =\prod_\ell\bigl((T_\ell E)^2
                                  \rtimes\operatorname{GL}_2(\mathbb Z_\ell)\bigr).
                                                               \tag{4.1}
\]
Consequently it contains the full product of the auxiliary affine
factors outside a fixed finite set, with identity on every factor
in that set. Its Haar distribution is the product of the full
outside factors and its possibly proper image on the exceptional
factors.

**Proof.** Let `I_L` bound the index of the linear image at every
level, by Serre's theorem, and retain the uniform translation-kernel
index bound `c_L`. At every positive integer `m`, the affine image
has index at most `I_L c_L` in
`E[m]^2 rtimes GL_2(Z/m)`, by the subgroup index formula for its
linear projection and translation kernel. These are projections of
one compact adelic image `H`.

The index `[A:H]` is the supremum of these finite indices. Indeed
any finitely many distinct cosets of the closed subgroup `H` are
separated in a common finite quotient, using compactness and the
congruence neighbourhood basis. Thus `[A:H]<=I_L c_L` is finite.
A closed subgroup of finite index in a compact group is open:
each coset is closed, and its complement is a finite union of closed
cosets. Hence `H` contains one congruence kernel in (4.1), which
contains all factors outside the prime divisors of its modulus.

Since every outside-factor element with identity elsewhere is in
`H`, it can cancel the outside coordinates of any element of `H`.
This proves the asserted direct-product decomposition of `H` and
therefore of normalized Haar measure. The argument uses the full
finite-level affine indices; it does not identify a finite-level
translation kernel with a kernel over the infinite torsion field.
\(\square\)

Openness is not a product law for reductions at different residue
primes. It describes their joint auxiliary division fields and the
associated Frobenius distribution. Passing from fixed auxiliary
information to recurrence still requires a uniform estimate.
