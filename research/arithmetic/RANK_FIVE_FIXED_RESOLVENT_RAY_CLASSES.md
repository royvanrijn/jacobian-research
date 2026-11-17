# Ray-class characters and Hasse failures with a fixed quadratic resolvent

For every quadratic number field \(K\), infinitely many non-Galois
cubic fields \(L\) make \(K\times L\) Hasse-failing. This follows by
restricting the ray-class characters to conductors supported on primes
that split in \(K\) and are one modulo three. The argument concerns
abstract field algebras. It does not place infinitely many of them on
the existing elliptic slice of the fixed map.

## 1. The imported class-field interface

Fix a quadratic number field \(K\) of fundamental discriminant \(d_K\),
with nontrivial automorphism \(\sigma\). For a finite
\(\sigma\)-invariant modulus \(\mathfrak m\), let
\(\operatorname{Cl}_{\mathfrak m}(K)\) denote the ordinary finite ray
class group. A cubic character is a homomorphism to \(\mu_3\). Call it
**anti-invariant** if
\[
 \chi(\sigma a)=\chi(a)^{-1}.                         \tag{1.1}
\]
Its primitive conductor is the least modulus through which it factors.

Global class field theory associates to a nontrivial character (1.1)
a cyclic cubic extension \(N/K\). Its kernel is stable under \(\sigma\),
and the induced action on its order-three quotient is inversion.
Consequently \(N/\mathbb Q\) is Galois with group \(S_3\).
Its degree-three subfields form one \(\mathbb Q\)-isomorphism class
\(L\), with quadratic resolvent \(K\). The characters
\(\chi,\chi^{-1}\) give the same field, and these are the only
characters with that kernel.

The exact imported counting interface is Cohen--Thorne,
*Dirichlet series associated to cubic fields with given quadratic
resolvent*, arXiv:1301.3563v2 (14 August 2013), equation (3.3), printed
p.6, and Proposition 3.6, equation (3.8), printed p.8:
[primary PDF](https://arxiv.org/pdf/1301.3563v2),
[retained PDF](sources/cohen-thorne-2013-fixed-quadratic-resolvent-v2.pdf),
[metadata](sources/cohen-thorne-2013-fixed-quadratic-resolvent-v2.json).
For a positive integer \(m\), it identifies the number of cubic-field
isomorphism classes with discriminant \(d_K f^2,\ f\mid m\), with
\[
 \frac{3^{r_m}-1}{2},\qquad
 r_m=\dim_{\mathbb F_3}
   \left(\operatorname{Cl}_{(m)}(K)/
                     \operatorname{Cl}_{(m)}(K)^3\right)^-.       \tag{1.2}
\]
Here minus denotes inversion under \(\sigma\). These are fields,
not cubic orders or automorphism-weighted counts. The general
character correspondence and its equivariance under \(\sigma\) are
the global class-field input recalled there; they are not reproved here.
The brief printed proof of Proposition 3.6 writes \(D^*\) where the
definition of \(H'\) and the stated formula use \(D\). Formula (1.2)
follows the stated formula and its definition, with \(K=\mathbb Q(\sqrt D)\).

For this correspondence the conductor ideal of \(N/K\) is
\[
 \mathfrak f(N/K)=f\mathcal O_K,\qquad D_L=d_K f^2.     \tag{1.3}
\]
The precise primary interface is Cohen--Morra, *Counting Cubic
Extensions with given Quadratic Resolvent*, arXiv:1003.1869v1,
Introduction, printed pp.1--2:
[primary PDF](https://arxiv.org/pdf/1003.1869v1),
[retained PDF](sources/cohen-morra-2010-fixed-quadratic-resolvent-v1.pdf),
[metadata](sources/cohen-morra-2010-fixed-quadratic-resolvent-v1.json).
The retained version identifier is used rather than its inconsistent
internal title-page date. Its conductor statement is imported; the
unfiltered counting constant in its Theorem 1.1 is not an asymptotic
constant for the Hasse subfamily.

## 2. Which anti-invariant characters give Hasse failures?

**Theorem RRC1 (character criterion).** Let \(\chi\) be a nontrivial
anti-invariant cubic ray-class character of \(K\), and let \(L\) be
its associated cubic field. Then
\[
 \boxed{\operatorname{Spec}(K\times L)
     \text{ is Hasse-failing}
 \quad\Longleftrightarrow\quad
 \chi\text{ is unramified at every finite prime }
       \mathfrak p\text{ with }\sigma\mathfrak p=\mathfrak p.}    \tag{2.1}
\]
Equivalently, the primitive conductor of \(\chi\) is supported only
on pairs of primes exchanged by \(\sigma\).

**Proof.** Neither field factor is \(\mathbb Q\), so there is no
rational point. The real cubic always has a real root.

At a rational prime \(p\) split in \(K\), the quadratic algebra is
\(\mathbb Q_p^2\), which already supplies a point, independently of
the cubic character.

Otherwise there is a unique prime \(\mathfrak p\) of \(K\) above \(p\),
and it is fixed by \(\sigma\). Suppose first that \(\chi\) is
unramified there. Its Artin value satisfies
\[
 \chi(\mathfrak p)=\chi(\sigma\mathfrak p)
                         =\chi(\mathfrak p)^{-1}.
\]
Since the value has order dividing three, it is one. Equivalently,
the Frobenius in \(N/K\) is trivial. The extension \(N/K\) therefore
splits completely at \(\mathfrak p\). This argument evaluates the
primitive unramified character, even if a larger chosen ray modulus
contains \(\mathfrak p\). The decomposition subgroup in
\(S_3=\operatorname{Gal}(N/\mathbb Q)\) has order two. In its action on
the three roots it fixes one root, so \(L\) has a
\(\mathbb Q_p\)-factor.

If \(\chi\) is ramified at this fixed prime, the local cyclic cubic
extension \(N/K\) has ramification degree three. Its local degree
over \(\mathbb Q_p\) is six, so the decomposition subgroup is all
of \(S_3\). It fixes none of the three roots. Thus \(L\) has no
\(\mathbb Q_p\)-factor, and the nonsplit quadratic has none either.
This proves failure of local solubility at \(p\).

These alternatives prove both directions at every finite place,
including two and three. They are equivalent to the
[field-discriminant criterion RDF1](RANK_FIVE_DISCRIMINANT_CLASSIFICATION.md).
In particular ramification at three is permitted only when \(K\)
splits over \(\mathbb Q_3\); in that case the possible factor in
(1.3) is \(3^2\). \(\square\)

The fixed-prime criterion concerns ramification, not all local
character values: at a fixed prime, anti-invariance already forces an
unramified order-three character to have trivial Frobenius.

## 3. A quantitative finite-modulus lower bound and infinitude

Set
\[
 u_K=
 \begin{cases}
  1,&K\text{ real, or }K=\mathbb Q(\sqrt{-3}),\\
  0,&K\text{ imaginary and }K\ne\mathbb Q(\sqrt{-3}).
 \end{cases}                                                    \tag{3.1}
\]

**Theorem RRC2 (fixed-resolvent ray-class lower bound).** Let
\(p_1,\ldots,p_t\) be distinct rational primes which split in \(K\)
and satisfy \(p_i\equiv1\pmod3\). Suppose \(t\ge u_K\), and put
\(m=\prod_{i=1}^t p_i\). There are at least
\[
 \boxed{\frac{3^{\,t-u_K}-1}{2}}                     \tag{3.2}
\]
distinct cubic-field isomorphism classes \(L\) such that
\(K\times L\) is Hasse-failing and
\[
 f\mid m,\qquad |D_L|\le |d_K|m^2.                   \tag{3.3}
\]
In particular every quadratic number field is the quadratic
resolvent of infinitely many Hasse-failing rank-five algebras.

**Proof.** Consider the ray-class exact sequence
\[
 \mathcal O_K^*\longrightarrow
 R=(\mathcal O_K/m\mathcal O_K)^*
 \longrightarrow\operatorname{Cl}_{(m)}(K)
 \longrightarrow\operatorname{Cl}(K)\longrightarrow1.           \tag{3.4}
\]
It is the exact sequence (3.7) on printed p.7 of the retained
Cohen--Thorne text. Taking three-primary parts and then minus parts
preserves exactness. For the second operation this follows from the
idempotents \((1\pm\sigma)/2\), since two is invertible on a finite
abelian three-group.

Each split pair above \(p_i\) gives
\[
 (\mathcal O_K/p_i\mathcal O_K)^*
                   \simeq\mathbb F_{p_i}^*\times\mathbb F_{p_i}^*,
\]
with \(\sigma\) exchanging the factors. Its three-primary minus
part consists of \((a,a^{-1})\), and is a nontrivial cyclic
three-group because \(3\mid p_i-1\). Thus the three-primary minus
part of \(R\) is a product of \(t\) nontrivial cyclic groups and
has three-rank \(t\).

The unit image can remove at most \(u_K\) generators in that part.
Indeed the global unit group has free rank one in the real
quadratic case and zero in the imaginary case; among quadratic
fields its torsion has a nontrivial three-part only for
\(\mathbb Q(\sqrt{-3})\), where that part is cyclic of order three.
Prime-to-three torsion contributes nothing. Projecting the unit
image to the minus part cannot increase this generator bound.

Consequently the quotient of the three-primary minus residue group
by the unit image has three-rank at least \(t-u_K\).
It injects into the three-primary minus ray class group by (3.4).
For finite abelian three-groups \(A\subset B\), the inequality
\(\dim A[3]\le\dim B[3]\) shows that the three-rank cannot decrease
on passing from the subgroup to the ambient group. Hence \(r_m\)
in (1.2) satisfies
\[
 r_m\ge t-u_K.                                      \tag{3.5}
\]
This uses the rank inequality for subgroups, not an assertion that
\(A/3A\) injects into \(B/3B\).

There are \(3^{r_m}-1\) nontrivial anti-invariant cubic characters.
Each inverse pair has one kernel and gives one cubic-field
isomorphism class, producing (3.2). Different kernels give different
normal closures and hence different cubic-field isomorphism
classes. Their conductors divide \(m\mathcal O_K\); (1.3) gives
\(f\mid m\) and (3.3). Every ramified prime is in a split pair, so
RRC1 makes all these products Hasse-failing.

Finally apply Chebotarev to the finite Galois extension
\(K\mathbb Q(\zeta_3)/\mathbb Q\). Its completely split primes have
positive density. Outside the finitely many ramified primes they
split in both \(K\) and \(\mathbb Q(\zeta_3)\), the latter being
exactly the condition \(p\equiv1\pmod3\). Thus there are infinitely
many eligible primes. Letting \(t\) grow makes (3.2) unbounded.
No asymptotic in \(m\) or in \(D_L\) is asserted. \(\square\)

The remaining classical inputs are grounded in Milne, *Algebraic
Number Theory*, version 3.00 (11 February 2008),
[author-hosted primary text](https://www.jmilne.org/math/CourseNotes/ANT300.pdf),
[retained copy](sources/milne-ant300.pdf). Theorem 5.1, printed p.78
(PDF p.86), states that the unit rank is \(r_1+r_2-1\), giving the
unit cost used above. Theorem 8.31 and Corollary 8.32, printed p.140
(PDF p.148), give Chebotarev and positive density
\(1/[K\mathbb Q(\zeta_3):\mathbb Q]\) for the completely split primes.
These retained statements were inspected directly. The equivalence
between splitting in \(\mathbb Q(\zeta_3)\) and \(p\equiv1\pmod3\),
for \(p\ne3\), follows from a nontrivial cube root of unity in the
residue field. No effective prime bound is required.

The bound is deliberately only a finite-modulus lower bound. The
unit image and the class-group contribution can make the actual
number larger.

## 4. Consequences and the elliptic limitation

For fixed \(K\), distinct cubic fields yield distinct complete
algebras \(K\times L\): their field factors of degrees two and three
are intrinsic. Common Arithmetic Fibers realizes each such
rank-five finite étale algebra as a full fibre of some
Jacobian-one map on \(\mathbb A^3_{\mathbb Q}\); see
[FMA7](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md).
The realizing map is allowed to vary.

Each of the four quadratic fields on the existing elliptic slice
therefore has infinitely many abstract Hasse-failing algebras.
This does not prove that infinitely many associated ray-class
characters occur in that slice's cubic factors. The
[explicit elliptic criterion](ELLIPTIC_EXCEPTION_EXPLICIT_LOCAL_CRITERION.md)
still requires simultaneous valuation conditions at all varying
primes. Neither (3.2) nor an unfiltered fixed-resolvent counting
asymptotic controls this particular image.

The deductions here use the cited global class-field correspondence,
conductor identity, exact finite-group argument, and the existing
arithmetic and realization results. The retained PDFs and metadata
pin the imported interfaces. Global class field theory, the unit
theorem and Chebotarev are not replayed or formally
verified here. No claim of novelty, external review, an elliptic
Hasse-infinitude theorem, or a filtered Cohen--Morra asymptotic is made.
