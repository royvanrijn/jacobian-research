# Wild affine fields from a single Pell prime path

The [Pell prime family PI3](PELL_PERIOD_RING_CLASS_INERTIA.md#a-simultaneous-family-with-genuine-prime-classes)
gives an actual prime of high ordinary order in a real quadratic
field, with all designated ramified prime directions buried in
its ring class groups. The same element also gives scalar
affine fields whose wild discriminants can be determined exactly.
This uses its local radical depth at the ramified primes and the
orientation of the moving prime. At three the local conductor
calculation extends the mechanism of
[WD2–WD4](GENUS_THREE_PELL_WILD_DICKSON_FIELDS.md)
to this Pell family. The scalar extension over the cyclotomic
mirror is a separate extension from PI3's ring class extension
over the real quadratic field.

## A composite-degree different formula

**Theorem AW1.** Fix an odd integer \(m>1\) and a real quadratic
field \(K\) of discriminant \(\Delta\), disjoint from
\(\mathbb Q(\zeta_m)\). Suppose every \(p\mid m\) ramifies
in \(K\). Let \(\gamma\in K^\times\) be positive of norm one,
with ideal valuations divisible by \(m\). At each \(p\mid m\)
assume
\[
 v_{\mathfrak p}(\gamma-1)=1
 \quad\text{for the normalized valuation of }K_{\mathfrak p}.
 \tag{1}
\]
Put \(t=\gamma+\gamma^{-1}\).
The Dickson polynomial \(D_m(X,1)-t\) defines a degree-\(m\)
number field \(L\) of signature \((1,(m-1)/2)\).
Its normal closure \(W\) has group
\((\mathbb Z/m)\rtimes(\mathbb Z/m)^\times\), acting by
scalar affine maps. If \(B\) is the fixed field of translations,
then \(W/B\) is unramified at every rational prime outside
\(m\). The absolute discriminant of \(L\) is
\[
       |\operatorname{Disc}(L)|
       =\Delta^{(m-1)/2}\prod_{p^a\parallel m} p^{am}.
 \tag{2}
\]

Write \(C=\mathbb Q(\zeta_m)\), \(E=KC\) and \(v^m=\gamma\),
using the positive real root. For \(p\mid m\), (1) makes
\(\gamma\) a non-\(p\)th power in \(K\). Indeed a local
\(p\)th root of \(\gamma\) would have residue one, and the
\(p\)th power of a principal unit has difference from one
of valuation at least three in the ramified quadratic
completion. The prime-degree polynomial \(X^p-\gamma\)
is therefore irreducible over \(K\). One elementary check
is that a proper factor of degree \(e\) would give, by the
product of its roots, an element with \(p\)th power \(\gamma^e\)
up to a sign which is a \(p\)th power. Since \(p\nmid e\)
this would force \(\gamma\) to be a \(p\)th power.
The real-root extension of degree \(p\) is nonnormal over
the real field \(K\). It cannot lie in the abelian
extension \(E/K\). As \(E\) contains the \(p\)th roots of
unity, \(\gamma\) stays a non-\(p\)th power in \(E\).

Thus \(E(v)/E\) is cyclic of degree \(m\).
An automorphism of \(E\) with cyclotomic exponent \(a\)
and action \(\eta=\pm1\) on \(K\) lifts by \(v\mapsto v^\eta\),
and acts on translations by \(\eta a\). This gives a
split extension over \(\mathbb Q\). Its central involution
\(\tau\) combining both conjugations fixes the roots
\[
 \zeta_m^j v+\zeta_m^{-j}v^{-1},\qquad j\in\mathbb Z/m.
\]
The roots are distinct since \(\gamma^2\ne1\).
Their permutation action has kernel exactly \(\langle\tau\rangle\)
and the full scalar affine group as image. In particular
\(B=E^{\langle\tau\rangle}\) and the polynomial is irreducible
over \(\mathbb Q\). Complex conjugation fixes \(v\) and
negates the index \(j\), giving the signature.

Outside \(m\), Kummer theory for virtual units gives
unramifiedness over \(E\), and the quotient \(W/B\) inherits
it. A rational prime outside \(m\) that ramifies in \(K\)
has quadratic inertia acting by \(-1\) on the \(m\) roots.
The permutation representation has \((m-1)/2\) copies of
the local sign direction, giving its quadratic discriminant
contribution even if the quadratic inertia is wild at two.

Fix \(p^a\parallel m\), put \(n=p^a\), \(m_0=m/n\), and take
a finite unramified extension \(F/\mathbb Q_p\) containing
\(\mu_{m_0}\). Unramified base change preserves the exponent
of the discriminant of a separable local algebra. In the
ramified quadratic extension \(K_{\mathfrak p}F/F\) there is
a unique \(\gamma_0=\gamma^{1/m_0}\) close to one, with
\(v_{\mathfrak p}(\gamma_0-1)=1\) and norm one over \(F\).
Let \(v_0^n=\gamma_0\). The polynomial
\((1+X)^n-\gamma_0\) is Eisenstein over \(K_{\mathfrak p}F\).
Consequently \(v_0-1\) generates that extension as a
uniformizer, and its relative different exponent is \(2an\),
from the derivative \(n v_0^{n-1}\).

The roots with index zero on the \(\mu_{m_0}\)-factor
give the polynomial \(D_n(X,1)-(\gamma_0+\gamma_0^{-1})\).
Its shift by \(X=2\) is Eisenstein over \(F\):
\(D_n(X,1)\equiv X^n\bmod p\), and
\[
  \gamma_0+\gamma_0^{-1}-2=(\gamma_0-1)^2/\gamma_0
\]
has \(F\)-valuation one. At a root \(x=v_0+v_0^{-1}\),
\[
 D_n'(x,1)=n\frac{\gamma_0-\gamma_0^{-1}}{v_0-v_0^{-1}}.
\]
The numerator other than \(n\) has \(F\)-valuation \(1/2\);
the denominator has valuation \(1/(2n)\) because
\((v_0-v_0^{-1})^2=x^2-4\).
Thus this degree-\(n\) factor has different exponent
\(an+(n-1)/2\).

For a pair \(\{\eta,\eta^{-1}\}\subset\mu_{m_0}\) with
\(\eta\ne1\), the root \(x_\eta=\eta v_0+\eta^{-1}v_0^{-1}\)
generates \(F(v_0)\). The quadratic equation for \(v_0\)
has a unique simple solution near one over \(F(x_\eta)\),
since \(\eta-\eta^{-1}\) is a unit. That field contains
\(K_{\mathfrak p}F\) because \(v_0^n=\gamma_0\notin F\);
it has degree \(2n\) over \(F\).
The reductions \(x_\eta\bmod p\) distinguish the
\((m_0-1)/2\) opposite pairs and the zero index. The degrees
\(n\) and \(2n\) account for all \(m\) roots. The quadratic
different exponent one and the relative binomial different
give \(n+2an\) for each degree-\(2n\) factor.
The full exponent at \(p\) is therefore
\[
 an+(n-1)/2+\frac{m_0-1}{2}(2a+1)n=am+(m-1)/2 .
\]
This gives (2). The argument concerns field discriminants:
the local generators used for the different calculations
are uniformizers, and square factors in a polynomial
discriminant do not enter the formula.

## Exact conductor at every prime-power layer

**Theorem AW2.** In AW1 take \(m=p^a\). If \(p=3\) assume
\(K_{\mathfrak p}\ne\mathbb Q_3(\zeta_3)\).
For \(j\le a\), let \(W_j/B_j\) denote the normal extension
and translation-fixed field at degree \(p^j\), using the
same \(\gamma\). The cyclic \(W_a/B_a\) is totally ramified
at every prime above \(p\).
Let \(\delta=1\) if \(p\equiv1\bmod4\) and \(\delta=2\)
if \(p\equiv3\bmod4\). Every character of order \(p^j\)
has conductor exponent at each such prime
\[
             c_j=\frac{(p+3)p^{j-1}+2(\delta-1)}{2\delta}.
 \tag{3}
\]
The relative discriminant exponent there is
\[
 \sum_{j=1}^a\varphi(p^j)c_j
 =\frac{p+3}{2\delta(p+1)}(p^{2a}-1)
             +\frac{\delta-1}{\delta}(p^a-1).
 \tag{4}
\]
Outside \(p\) all these characters are unramified.

First, \(\gamma\) stays non-\(p\)th in the local \(E_a\).
Over \(E_1=K_{\mathfrak p}(\zeta_p)\) it is still non-\(p\)th
by the prime-to-\(p\) degree. If it became a \(p\)th power
over \(E_a\), its Kummer line over \(E_1\) would be that
of \(\zeta_p\), since \(E_1(\zeta_{p^2})/E_1\) is the
unique degree-\(p\) step. There is a local automorphism
fixing \(K_{\mathfrak p}\) and acting nontrivially on
\(\zeta_p\). For \(p\ge5\) this follows since
\([E_1:K_{\mathfrak p}]\ge(p-1)/2\), and for \(p=3\)
it is the stated local hypothesis. That automorphism
fixes the class of \(\gamma\) and acts nontrivially
on the entire nonzero line of \(\zeta_p\), a contradiction.
The extension over local \(E_a\) therefore has the full
group \(A=\mathbb Z/p^a\).
The same automorphism lifts and acts by a scalar
nontrivial modulo \(p\) on \(A\). Any abelian unramified
quotient would have trivial action under this conjugation
because it is a quotient of the cyclic residue field group.
It must be zero. Thus translations are totally ramified,
including after taking the central quotient over local \(B_a\).

For the local base, write \(C_j=\mathbb Q_p(\zeta_{p^j})\)
and \(K_{\rm cyc}\) for its unique quadratic subfield.
It is the ramified quadratic field of squareclass
\((-1)^{(p-1)/2}p\).
If \(K_{\mathfrak p}\ne K_{\rm cyc}\), the compositum
with \(C_j\) adds an unramified quadratic factor.
The restriction of \(\tau\) combines inversion on \(C_j\)
and the nontrivial action on \(K_{\mathfrak p}\).
Inversion on \(C_j\) acts on \(K_{\rm cyc}\) by
\((-1)^{(p-1)/2}\). It follows that a completion of \(B_j\)
has ramification index \(\varphi(p^j)/\delta\) and residue
degree \(f\), where \(f=2\) precisely if
\(p\equiv3\bmod4\) and \(K_{\mathfrak p}\ne K_{\rm cyc}\).
There are \(\delta/f\) such completions.
For \(\delta=1\), the local base is \(C_j\) or has the
same discriminant exponent after unramified quadratic
base change. For \(\delta=2\) it is the real cyclotomic
layer \(C_j^+\), with the unramified factor when \(f=2\).
Putting \(d_j=(j(p-1)-1)p^{j-1}\), its discriminant
exponent over \(\mathbb Q_p\) is
\[
              d(B_j)=f(d_j-\delta+1)/\delta .
\]
Here \(d_j\) comes from the derivative of
\(\Phi_{p^j}\) at a primitive root. The quadratic
step \(C_j/C_j^+\) is tame of different exponent one.

The permutation representation on the degree-\(p^j\)
roots contains each translation character once.
The primitive characters form \(\delta/f\) orbits under
the local scalar group, since they are faithful and
the local group has order \(\varphi(p^j)f/\delta\)
modulo translations. Each primitive orbit representation
is induced from a character over the local \(B_j\).
All these characters have the same conductor because
they are generators of the same character group.
By [Artin conductor induction](sources/schwein-formal-degree-2024.pdf),
the conductor contribution of one primitive orbit is
\(d(B_j)+f c_j\).
The imprimitive part is the representation of the
degree-\(p^{j-1}\) field. The difference of the exponents
in (2) at levels \(j\) and \(j-1\), with the level-zero
exponent zero, yields
\[
 \frac{\delta}{f}(d(B_j)+f c_j)
 = j p^j+\frac{p^j-1}{2}
       -(j-1)p^{j-1}-\frac{p^{j-1}-1}{2}.
\]
Solving gives (3) first over \(B_j\).

For \(a>j\), the first step \(B_{j+1}/B_j\) has conductor
exponent \(\lambda_j=(p^j+\delta-1)/\delta\). This follows by
transitivity of the different and the just-computed
quadratic-cyclotomic discriminants, then the
conductor-discriminant formula for that cyclic step.
Every nontrivial character of \(B_a/B_j\) has conductor
at least \(\lambda_j\ge c_j\).
Over the local \(B_j\), \(W_j\) and \(B_a\) are linearly
disjoint. A local scalar automorphism acting nontrivially
modulo \(p\) on translations exists: the local scalar
group is full, or the index-two subgroup of units in
the case \(p\equiv3\bmod4,\ K_{\mathfrak p}=K_{\rm cyc}\)
with \(p\ge7\). Such an automorphism acts trivially on
the abelian \(B_a/\mathbb Q_p\). Any common translation
quotient would have trivial action and hence be zero.

Let \(\chi\) be a primitive translation character at level
\(j\), and \(\lambda\ne1\) a character from \(B_a/B_j\).
If \(a(\lambda)>c_j\), the twist \(\chi\lambda\) has
conductor \(a(\lambda)\). If \(a(\lambda)=c_j\), their
last principal-unit characters lie in different
eigenspaces for that scalar automorphism. One is
multiplied by a scalar different from one modulo \(p\),
and the other is fixed. They cannot cancel, so the
twist has conductor \(a(\lambda)\) in this case too.
Induction of the base-changed \(\chi\) is the sum of
all the twists. Subtract the discriminant contribution
\(\sum_\lambda a(\lambda)\) in the induction formula.
The residue degree of the base change is one, so the
conductor remaining over \(B_a\) is exactly \(c_j\).
This is the mechanism of WD3 at \(p=3\).
Finally summing the characters' conductor exponents
gives (4), by the [conductor-discriminant formula](sources/milne-cft403.pdf).

## One actual split prime in all the fields

**Corollary AW3.** Use the fixed set \(S\), exponents \(a_p\),
\(B=\prod_{p\in S}p\), and \(m=\prod_{p\in S}p^{a_p}\)
of PI3. For every divisor \(e\mid m\), there are
\(\gg_{S,(a_p)} X^{1/(2m)}/(\log X)^2\) distinct real
quadratic fields \(K\) of discriminant at most \(X\) on
which all the following hold.

- The genuine moving prime \(q\) has exact order \(m\)
  in the ordinary class group and exact order \(Bm\)
  after lifting to the ring class group of conductor \(B\),
  with the inertia depths of PI3.
- The same norm element gives AW1's degree-\(m\) field
  with (2). The prime \(q\) splits in its translation-fixed
  field and has Frobenius of order \(e\) in the cyclic
  relative extension.
- For each \(p\in S\), the degree-\(p^{a_p}\) field from
  the same norm element satisfies AW2, with a unique
  prime over \(p\) in its cyclotomic mirror. The \(q\)
  Frobenius there has order \(p^{v_p(e)}\).

In PI3 put \(\gamma=\alpha/\bar\alpha\). Both factors are
positive for large \(q\), and
\((\gamma)=(\mathfrak q/\bar{\mathfrak q})^m\).
At every designated \(p\) we have (1) because the
coefficient of \(\sqrt D\) in \(\alpha\) is a unit and
the scalar coefficient is a unit. Also
\(n/p\equiv2\bmod p\) for PI3's integer parameter \(n\):
\(q^m-c^2\equiv-4p\bmod p^2\) and
\(A=2(Bc-1)\equiv-2\bmod p\).
Hence \(D/p\) is a square unit at \(p\), and
\(K_{\mathfrak p}=\mathbb Q_p(\sqrt p)\). This satisfies
the AW2 local hypothesis at three, and gives \(f=\delta\)
throughout this family. The quadratic field ramifies
at two, so it is disjoint from the odd cyclotomic field.

Add \(q\equiv1\bmod m\) to the progression in PI3.
This is compatible with the other prescribed residues,
which are prime to \(B\) except for \(q\equiv1\bmod B\).
At the oriented prime above \(q\), the normalization
\(\gamma/q^m=\bar\alpha^{-2}\) is a local unit. Its residue
is \(C^{-2}\), where
\[
                   C=\frac{c(Bc-2)}{Bc-1}.
\]
It is independent of \(q\). Because \(Bc-1=2^{e_{\rm aux}}s\)
in the auxiliary choice of PI3, \(C\) has valuation
\(-1\) at the fixed prime \(s\), which is prime to \(m\).
Thus \(C\) gives full Kummer degree \(m\) over
\(\mathbb Q(\zeta_m)\), and its normal closure over
\(\mathbb Q\) has the usual scalar affine group.
The commutator contains every translation, since
commuting with the scalar \(-1\) multiplies a translation
by \(-2\). The maximal abelian subfield is therefore
exactly \(\mathbb Q(\zeta_m)\). The cyclotomic field
imposing the progression intersects this Kummer field
in that subfield. [Chebotarev](sources/milne-cft403.pdf) gives a positive density of
primes satisfying the progression and any prescribed
translation order \(e\). This conjugacy condition fixes
the order; the inverse square in \(C^{-2}\) preserves
orders because \(m\) is odd. At the prime \(q\) the
quadratic and cyclotomic fields both split, so the
same residue computation gives the mirror Frobenius
and all its primary projections.

For the count, all moduli and the fixed auxiliary Kummer
field are chosen for the given \(S,(a_p)\). The Chebotarev
prime count remains \(\gg T/\log T\) up to \(T\).
The uniform positive-unit multiplicity bound from PI3
gives the claimed number of quadratic fields. Also \(K\)
is determined by the degree-\(m\) affine field for fixed
\(m\) in this family. Its normal closure recovers \(B_m\)
as the fixed field of the commutator. Adjoining
\(\mathbb Q(\zeta_m)\) recovers \(K(\zeta_m)\).
Any other quadratic subfield outside \(\mathbb Q(\zeta_m)\)
differs from \(K\) by a nontrivial cyclotomic quadratic
character. Such a character ramifies at some \(p\in S\),
where its product with \(K\) is unramified, since both
quadratic factors there are ramified. This cannot satisfy
PI3. In particular this gives
\(\gg_{S,(a_p)}Y^{1/(m(m-1))}/(\log Y)^2\) distinct
degree-\(m\) fields with discriminant at most \(Y\) on
each prescribed Frobenius slice.

## Verification

The [exact checker](verify_pell_prime_path_wild_affine_fields.py)
verifies finite Dickson field discriminants including a
mixed degree-fifteen case, the conductor arithmetic, and a
joint true-prime fixture from PI3 with \(m=45\) and Frobenius
order \(45\). Run it with Sage's Python.
The local different and conductor proof, Chebotarev
selection, and uniform field count are written arguments.
The bounds fix the ramified primes and their powers.
They do not estimate the unrestricted sum over Pell kernels.

## All wild conductors in the composite mirror

**Corollary AW4.** In AW3 suppose \(m\) has at least two
distinct prime divisors. For \(p^a\parallel m\) put
\(m_0=m/p^a\) and \(d=\operatorname{ord}_{m_0}(p)\).
Use \(\delta\) as in AW2. Every prime of the composite
translation-fixed field \(B\) above \(p\) has ramification
index \(\varphi(p^a)\) over \(\mathbb Q\). Its residue
degree \(f_B\) is \(\operatorname{lcm}(\delta,d)\), except
that it is half this value when \(p\equiv3\bmod4\),
\(4\mid d\), and \(p^{d/2}\equiv-1\bmod m_0\).
The local relative group of \(W/B\) at each such
prime is exactly the \(p\)-primary part of the global
translation group and is totally ramified.

For a character of the composite translation group,
its conductor exponent there is zero if its
\(p\)-primary part is trivial. For a \(p\)-primary
part of order \(p^j>1\) the exponent is
\[
                           \frac{p+3}{2}p^{j-1}.
 \tag{7}
\]
Consequently the relative discriminant exponent of
\(W/B\) at each prime over \(p\) is
\[
                    \frac{p+3}{2(p+1)}(p^{2a}-1).
 \tag{8}
\]
This gives the conductor of each character simultaneously
at every designated prime. The extension is unramified
at the remaining finite primes.

With \(C=\mathbb Q(\zeta_m)\) and \(E=KC\),
consider the local tower from the prime-power cyclotomic
mirror \(B_a\) of AW2 to a completion of the composite
mirror \(B\). Adjoining \(\mu_{m_0}\) to \(E_a\) is
unramified of cyclotomic residue degree \(d\) over
\(\mathbb Q_p\). The completion of \(E\) thus has
residue degree \(\operatorname{lcm}(\delta,d)\).
The central involution acts nontrivially on
\(\mu_{m_0}\), so it does not belong to inertia there.
Thus \(e(B/\mathbb Q_p)=\varphi(p^a)\), and
\(e(B/B_a)=\delta\). If \(p\equiv1\bmod4\),
\(K_{\mathfrak p}=K_{\rm cyc}\) and the central
involution is already absent from the local
decomposition group of \(E_a\). If \(p\equiv3\bmod4\),
its restriction to \(E_a\) exists and fixes that
field's unramified quadratic factor. It extends to
the local decomposition group in \(E\) precisely
if inversion of \(\mu_{m_0}\) is an even power of
Frobenius. This is the stated condition on \(d\)
and \(p^{d/2}\). When it extends, quotienting removes
an unramified quadratic factor from \(E\), proving
the formula for \(f_B\).

The \(m_0\)th radical of \(\gamma\) is already locally
in \(K_{\mathfrak p}\) as used in AW1. Over \(E\) all
its conjugates are present at \(p\), so the only
local translations are the full \(p\)-part of AW2,
which stays totally ramified under this unramified
extension from \(E_a\). Their images in \(W\) have
the same degrees. Equivalently \(W/B\) at this prime
is obtained by adjoining \(W_a\) locally to \(B\).
The extension \(B/B_a\) is tame, whereas \(W_a/B_a\)
is totally wildly ramified, so the local layers are
disjoint. Every nontrivial character of \(B/B_a\)
has conductor at most one. For a wild character of
conductor \(c_j\) from AW2, all of its twists by
these characters have conductor \(c_j\). By Artin
induction and the tame discriminant \(f(\delta-1)\)
of the base change of residue degree \(f\),
the new conductor is \(\delta c_j-(\delta-1)\).
This equals (7), and summing over the local
translations gives (8).
