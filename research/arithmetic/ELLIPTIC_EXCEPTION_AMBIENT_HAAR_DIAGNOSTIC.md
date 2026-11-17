# Ambient and rational-closure Haar measures on the elliptic slice

The independent ambient product gives measure zero to the all-primes
soluble set. This does not determine the measure in the actual rational
closure: that closure itself has ambient measure zero. The difference
can be quantified at every finite set of good primes. A local formula
also identifies the formal depths that affect the actual pole
obstructions.

## 1. The curve and the two measures

Retain
\[
 E_0:\ y^2=x^3-10656x^2+27091730x,\quad
 F(x)=128x-25865,\quad q_0=2709173,
\]
\[
 A=\{1,2,5q_0,10q_0\},\qquad
 \Gamma=E_0(\mathbb Q),\qquad
 \alpha:\Gamma\longrightarrow A.
\]
The [arithmetic theorem](ELLIPTIC_EXCEPTION_ARITHMETIC.md) gives
rank two, rational torsion \(C_2\), and the indicated descent image.
The elliptic discriminant is
\[
 \Delta_E=16(27091730)^2(10656^2-4\cdot27091730)
         =2^9 5^2 7^3 1889\,q_0^2.                  \tag{1.1}
\]
For a finite set \(S\) of good primes, put
\[
 \mathcal G_S=A\times\prod_{p\in S}E_0(\mathbb Q_p),\qquad
 H_S=\overline{\{(\alpha(P),(P)_p):P\in\Gamma\}}
                                      \subseteq\mathcal G_S.
\]
Let \(\lambda_S\) be normalized Haar measure on the **ambient**
group \(\mathcal G_S\); the intrinsic measure on \(H_S\) is the
different normalized Haar measure of EFA1--EFA3.

Likewise let \(\mathcal G_\infty=A\times\prod_p E_0(\mathbb Q_p)\)
carry product Haar measure \(\lambda_\infty\), and let
\(H_\infty\) be the rational closure of
[EFA3](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md), with its
intrinsic normalized Haar measure \(\nu_\infty\).
No product independence is assumed for \(\nu_\infty\).

Write
\[
 \mathcal P=\{p\ge11:p\equiv2\pmod3,\ p\nmid\Delta_E\},\quad
 d=51730,\quad x_0=\frac{25865}{128},\quad c_0=\frac{639735}{128}.
\]
The exact identity
\[
 f_E(x_0)=x_0c_0^2,\qquad x_0=d/16^2                 \tag{1.2}
\]
shows that the zero divisor of \(F\) consists of
\[
 R_\pm=(x_0,\ \pm639735\sqrt d/2048),
 \qquad \operatorname{div}(F)=R_++R_--2O.            \tag{1.3}
\]
The two zeros form a closed point of degree two, with field
\(\mathbb Q(\sqrt d)\). For \(p\in\mathcal P\), all the constants
in (1.2) are units: their odd prime factors other than three and
five are \(7,739,42649\), each one modulo three.

## 2. The ambient all-primes measure is zero

**Theorem EHD1 (ambient local fraction and product measure).**
Let \(N_p=\#E_0(\mathbb F_p)\), and let \(V_p\) be the admissible
local soluble set of EFA1. For \(p\in\mathcal P\), its normalized
Haar measure in the full local group \(E_0(\mathbb Q_p)\) is
\[
 \boxed{\lambda_p(V_p)
 =1-\frac{b_p}{N_p}
          \left(1-\frac1{p^2+p+1}\right),\qquad
 b_p=1+2\,\mathbf1_{(d/p)=1}.}                       \tag{2.1}
\]
Consequently the set of tuples in \(\mathcal G_\infty\) soluble
at every prime of \(\mathcal P\), and hence the everywhere locally
soluble set, has \(\lambda_\infty\)-measure zero.

**Proof of (2.1).** Good reduction gives \(N_p\) residue discs,
each of Haar measure \(1/N_p\). Outside the pole disc \(O\) and
the two possible zero discs in (1.3), \(F\) is a unit. The zero
discs exist exactly when \((d/p)=1\); their \(y\)-coordinates are
units. Hensel lifting gives their two local centers \(R_\pm\).
In either zero disc, the valuation of \(F\) is the formal depth
from the corresponding center, since \(F\) has a simple zero
and \(dx/(2y)\) is a unit differential there. In the pole disc,
a formal point of depth \(n\ge1\) has \(v_p(F)=-2n\).
By ELV1, local solubility in each of these \(b_p\) discs is exactly
\(3\mid n\). Its conditional mass is
\[
 \sum_{j\ge1}(1-p^{-1})p^{1-3j}
                         =\frac1{p^2+p+1}.
\]
The finitely many excluded local points have measure zero.
This proves (2.1).

**Proof of the ambient zero measure.** The elementary bound
\(N_p\le2p+1\) follows by counting at most two \(y\)-values for
each affine \(x\), and adding \(O\). Thus the insoluble mass in
(2.1) is at least \(1/(6p)\).

The set \(\mathcal P\) has natural density \(1/2\) among primes:
omit finitely many primes from the nonsplit primes in
\(\mathbb Q(\zeta_3)\), and apply Chebotarev. In particular
\(\sum_{p\in\mathcal P}1/p=\infty\). For completeness, positive
relative prime density and partial summation reduce this to
Euler's divergence \(\sum_p1/p=\infty\). The latter follows
by taking logarithms of the Euler product for \(\zeta(s)\) as
\(s\downarrow1\); the higher prime-power terms stay bounded,
while \(\zeta(s)\) diverges.

The ambient measure is a product by definition. Its finite
soluble-cylinder measures are at most
\(\prod_{p\in S\cap\mathcal P}(1-1/(6p))\), which tend to zero.
Continuity of probability measure proves the assertion.
This is an ambient product calculation, not a calculation on
\(H_\infty\). \(\square\)

The exact density input is Milne, *Algebraic Number Theory*,
version 3.00, Theorem 8.31, printed p.140 (PDF p.148):
[primary text](https://www.jmilne.org/math/CourseNotes/ANT300.pdf),
[retained copy](sources/milne-ant300.pdf).
The local root and branch arguments are those of ELV1 and
[the eligibility theorem](ELLIPTIC_EXCEPTION_PRIME_ELIGIBILITY.md).

## 3. The rational closure is also ambient-null

**Theorem EHD2 (support index and singularity).** For any finite set
\(S\) of \(s\) good odd primes, put
\(r_p=|E_0(\mathbb F_p)/2E_0(\mathbb F_p)|\). Then
\[
 \boxed{[\mathcal G_S:H_S]\ge\frac12\prod_{p\in S}r_p
                    \ge2^{s-1},\qquad
 \lambda_S(H_S)\le2^{1-s}.}                         \tag{3.1}
\]
An index may be regarded as infinite in this inequality; in fact
\(H_S\) is open for every finite \(S\), as explained below.
The full closure satisfies
\[
 \boxed{\lambda_\infty(H_\infty)=0,\qquad
                         \nu_\infty(H_\infty)=1.}    \tag{3.2}
\]
Thus its intrinsic Haar measure is singular with respect to
ambient product Haar measure.

**Proof.** The finite quotient
\[
 \mathcal G_S\longrightarrow
       A\times\prod_{p\in S}
              E_0(\mathbb F_p)/2E_0(\mathbb F_p)     \tag{3.3}
\]
is surjective by good reduction. The rational two-torsion point
\(T_0=(0,0)\) remains nonzero at every good odd prime. Therefore
each \(r_p\ge2\).
The rational image in (3.3) factors through
\[
 \Gamma/2\Gamma,\qquad
 |\Gamma/2\Gamma|=2^{2+1}=8,
\]
using rank two and rational torsion \(C_2\). Taking closure cannot
enlarge an image in a finite quotient. Hence the image of \(H_S\)
has at most eight elements, while the ambient quotient has
\(4\prod r_p\) elements. Both the index and Haar inequalities
in (3.1) follow. The measure inequality is valid even without
using openness. For small \(s\), the trivial bound by one also
applies.

Choose arbitrarily many distinct good primes. The inverse image
of \(H_S\) under projection of the full ambient product has measure
\(\lambda_S(H_S)\), and contains \(H_\infty\). Letting \(s\)
grow proves the first identity of (3.2). The second is its defining
intrinsic normalization. \(\square\)

For clarity, finite-place openness is compatible with (3.1); see also
[EQA1](ELLIPTIC_EXCEPTION_QUANTITATIVE_FINITE_LOCAL_ASYMPTOTICS.md).
Choose deep formal subgroups at the finitely many primes so each
is additively isomorphic to \(\mathbb Z_p\). Multiply one rational
nontorsion point by an integer killing \(A\) and the corresponding
finite local quotients. Its formal coordinate is nonzero at each
prime. Integer Chinese remainder approximation shows that its
cyclic closure contains a product of open formal subgroups.
This product lies in \(H_S\), so \(H_S\) is open. Its index must
nevertheless grow at least as in (3.1).
The index measures the support of the rational image. It is not
the local soluble fraction \(\mu_{S,a}\).

Ambient nullness alone therefore places no upper bound on the
intrinsic measure of the all-soluble subset of \(H_\infty\). EHD1 therefore proves neither vanishing nor positivity
of the all-primes coefficient in EFA3.

## 4. The actual local image and formal depth

**Proposition EHD3 (no anomalous good primes and the pole fraction).**
At every good prime \(p\ge7\),
\[
 \boxed{p\nmid N_p.}                                \tag{4.1}
\]
Let \(J_p\) be the actual finite image of
\(\Gamma\to A\times E_0(\mathbb F_p)\).
Under the prime-to-\(p\) torsion splitting of the good-reduction
sequence, the closure of the joint rational image in
\(A\times E_0(\mathbb Q_p)\) is exactly
\[
 \boxed{J_p\times E_{e_p}(p)}                        \tag{4.2}
\]
for some integer \(e_p\ge1\).
It is characterized by
\[
 e_p=\min_{Q\in\Gamma,\ Q\text{ nontorsion}}
                   v_p\!\left(t([N_p]Q)\right),
 \qquad t=-X/y.                                    \tag{4.3}
\]
For \(p\in\mathcal P\), put
\(\delta_p=3\lceil e_p/3\rceil-e_p\in\{0,1,2\}\).
For each pair \((a,O)\in J_p\), the conditional soluble fraction
in its formal component is
\[
 \boxed{s_{p,e_p}=\frac{p^{2-\delta_p}}{p^2+p+1}.}    \tag{4.4}
\]
Thus its absolute insoluble pole mass in the whole joint local
closure is
\[
 \frac1{|J_p|}\left(1-s_{p,e_p}\right),             \tag{4.5}
\]
and is zero for class \(a\) if \((a,O)\notin J_p\).

**Proof of (4.1).** Hasse's bound gives
\(0<N_p\le p+1+2\sqrt p<2p\) for \(p\ge7\).
If \(p\mid N_p\), it would follow that \(N_p=p\).
But \(N_p\) is even because of \(T_0\), while \(p\) is odd.
This proves (4.1).
The precise imported bound is Milne, *Elliptic Curves* (2006),
Chapter IV, Theorem 9.4, printed pp.149--150 (PDF pp.157--158):
[primary book](https://www.jmilne.org/math/Books/ectext6.pdf),
[retained copy](sources/milne-elliptic-curves-2006.pdf).
The retained statement and proof were inspected.

For (4.2), the good-reduction sequence has formal kernel
\(E_1(p)\simeq\mathbb Z_p\), and its quotient has order \(N_p\)
prime to \(p\). It splits uniquely into prime-to-\(p\) torsion
and the formal group. Explicitly, subtracting
\([N_p]^{-1}[N_p]Q\) in the formal group from a lift \(Q\) gives
the unique torsion lift of its residue. Multiplication by \(N_p\)
is a unit automorphism on the formal group.

The closed formal image of the finitely generated group is
\(E_{e_p}(p)\). It is nonzero because a rational nontorsion point
cannot become torsion in \(E_0(\mathbb Q_p)\). Its minimum depth
is (4.3), since multiplication by \(N_p\) preserves formal depth.
The finite group \(J_p\) has order prime to \(p\), so integer
Chinese remainder approximation separates its finite coordinate
from the pro-\(p\) formal coordinate. This proves the product
description (4.2).

Within \(E_{e_p}(p)\), depth \(n\ge e_p\) has conditional mass
\((p-1)p^{e_p-n-1}\). At the pole, ELV1 requires \(3\mid n\).
Summing from \(n=3\lceil e_p/3\rceil\) in steps of three gives
(4.4), and the uniform finite-coordinate mass gives (4.5).
\(\square\)

Without fixing \(a\), let \(R_p\) be the reduction image of
\(\Gamma\) in \(E_0(\mathbb F_p)\). Summing (4.5) over the pairs
above \(O\) gives total pole obstruction mass
\[
 \frac1{|R_p|}(1-s_{p,e_p}).                        \tag{4.6}
\]
If \(e_p=1\), this is
\(|R_p|^{-1}(1-1/(p^2+p+1))\), of order at least \(1/p\).
If \(3\mid e_p\), its conditional factor is instead
\((p+1)/(p^2+p+1)\), of order \(1/p\) before division by
\(|R_p|\). No distribution theorem for \(e_p\) or \(R_p\) is
asserted.

The minimum in (4.3) can be taken over any finite set of
Mordell--Weil generators. A known finite-index subgroup gives
the same depth outside the finitely many primes dividing its
index. Existence of one rational witness of depth one already
proves \(e_p=1\), without a generator basis. The retained witnesses
at 11 and 17 supply such examples.

## 5. What this does and does not settle

These statements give an exact ambient zero-measure theorem,
a quantitative support-index obstruction to transferring that
theorem to the rational closure, and a local formula involving
the actual reduction image and formal depth. They do not prove
that the intrinsic all-primes measure vanishes or is positive.

Even divergent marginal obstruction masses would not by themselves
give a product law for the actual diagonal closure. The relevant
losses after imposing earlier prime conditions are conditional
losses in that closure. The exact product criterion in
[EFA4](ELLIPTIC_EXCEPTION_FINITE_LOCAL_ASYMPTOTICS.md#6-conditional-prime-losses-and-the-limiting-measure)
uses these actual conditional losses. The ambient lower bound of
EHD1 does not bound them. EFA3's arithmetic tail condition and the
prime eligibility restrictions of EPE2 remain separate obligations.

Run

```sh
python3 research/arithmetic/verify_elliptic_exception_ambient_haar.py --write-receipt
```

The [checker](verify_elliptic_exception_ambient_haar.py) verifies the
discriminant and branch identities, finite good-reduction groups,
the finite two-quotient image bound, exact local mass factors and
bounded formal-depth sums, and two depth-one rational witnesses.
The finite-field controls cover every good prime from 7 through 251;
branch and ambient-mass controls use those primes at least 11 that
are two modulo three. Joint two-quotient images use the six primes
11, 13, 17, 19, 23, 29. Formal-depth controls use every active prime
at most 101 and depths 1 through 12. These are stated finite domains,
not evidence for an unrestricted depth distribution.
The [receipt](elliptic_exception_ambient_haar_checks.json) states
the finite coverage and source hashes. The unrestricted measure,
group and density statements are proved above using the stated
classical inputs. No finite computation is promoted to an
all-primes distribution assertion, and no formal verification
or external review is claimed.
