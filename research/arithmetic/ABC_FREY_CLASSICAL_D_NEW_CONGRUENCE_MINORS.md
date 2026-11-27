# Classical congruence minors for the fixed Frey norm

Use the [fixed-level family](ABC_FREY_FIXED_LEVEL_UNWEIGHTED_MONODROMY.md)
\(a=3x^2,\ b=y^2,\ c=z^3\) with \(32\cdot3\cdot5\cdot7\cdot11\mid a\).
Let \(R=\operatorname{rad}(abc)\), \(M\in\{110,770\}\) as in (FL6),
\(D=R/M\), \(n_p=v_p(\Delta_E^{\min})\), and \(H=H_{D,M,5}\).
Write \(\delta_{D'}=\delta_{D',R/D'}\) for the optimal Shimura degree for
each \(D'\mid D\) with an even number of primes. Use \(f\) for the
normalized rational newform of conductor \(R\) attached to \(E\).

The prime-to-\(\prod_{\ell\le163}\ell\) part of \(H\) is an exact gcd of
congruence indices from classical cusp forms after division by known
Tate factors. The comparison uses the original \(D\)-new space away
from an explicit set of primes dividing \(R\). Including the classical
full space in the gcd gives a purely classical integer whose ratio to
\(H\) has rational height \(O(1+\omega(D))\).

## The congruence index from Fourier minors

For \(D'\mid D\) as above, let \(V_{D'}\subset S_2(\Gamma_0(R),\mathbb Q)\)
be the subspace new at every \(p\mid D'\), allowing oldforms at other
primes. Over \(\mathbb Q\), decompose it as \(\mathbb Qf\oplus W_{D'}\)
using the idempotent \(e_f\) which is one on \(f\) and zero on all
other systems of good Hecke eigenvalues. Strong multiplicity one and
the exact conductor of \(f\) make its summand a line.

Let \(T_{D'}\) be the image of the full integral level \(R\) Hecke algebra
on this space, including the operators at primes dividing \(R\). Put
\(r_{D'}=\min\{r>0:r e_f\in T_{D'}\}\). This is Hamidi's \(r^1_{D'}(R)\).
Define the Sturm bound and the coefficient map by
\[
 B=\left\lfloor\frac R6\prod_{p\mid R}(1+p^{-1})\right\rfloor,
 \qquad A_B(g)=(a_1(g),\ldots,a_B(g))^{\mathsf T}.
\]
Set \(F=A_B(f)\) and
\(K_{D'}= A_B(W_{D'})\cap\mathbb Z^B\). Use any integer basis matrix
\(\mathbf K_{D'}\) of this saturated lattice and put
\(d=1+\operatorname{rank}K_{D'}\). Then
\[
 \boxed{\displaystyle
   r_{D'}=\gcd\{\text{absolute values of the }d\times d
                  \text{ minors of }[\mathbf K_{D'}\mid F]\}.}
 \tag{CD1}
\]
Equivalently, \(v_\ell(r_{D'})\ge k\) if and only if
\(F\bmod\ell^k\) belongs to the span of the columns of
\(\mathbf K_{D'}\bmod\ell^k\). This tests congruence to an integral
form in the rational complement, which can be a sum of several
packets.

Indeed the integral \(a_1\)-pairing is perfect between the full Hecke
order and the integral cusp forms. It identifies \(T_{D'}^\vee\) with
the integral forms in \(V_{D'}\), since the latter are the annihilator
of the kernel of the Hecke quotient. The operators \(T_n\) for
\(1\le n\le B\) generate the full order additively, and their images
generate the quotient. Hence \(A_B\) identifies \(T_{D'}^\vee\)
with \(A_B(V_{D'})\cap\mathbb Z^B\).
[Agashe–Stein, Sturm generation and the integral pairing](sources/agashe-stein-hecke-generation-sturm.pdf).

The primitive intersection \(T_{D'}\cap\mathbb Q e_f\) has generator
\(r_{D'}e_f\). Every functional on this primitive rank-one
sublattice extends to \(T_{D'}\). Thus projection of \(T_{D'}^\vee\)
to the \(f\)-line has image \((1/r_{D'})\mathbb Zf\), with kernel the
integral forms in \(W_{D'}\). The cyclic quotient by \(K_{D'}+\mathbb ZF\)
in its saturation has order \(r_{D'}\), proving (CD1) by Smith normal
form. The congruence test follows because \(A_B(V_{D'})\cap\mathbb Z^B\)
is saturated and \(K_{D'}\) is the kernel of the projection. These
are finite matrices at \(R\). One can obtain \(W_{D'}\) using the old
and new decomposition or the sum of the images of finitely many
\(T_q-a_q(f)\) with \(q\nmid R\), and obtain \(\mathbf K_{D'}\) by
saturating its truncated rational span.

## Exact comparison at large coefficient primes

For odd \(\ell\), Hamidi's Theorems A and B imply
\(v_\ell(r_{D'})\ge v_\ell(\delta_{D'})\), with equality for
\(\ell>3\) if the residual representation is absolutely irreducible
and there is at most one noncontrollable prime in \(D'\). The
criterion includes \(\ell\mid R\), and \(p=\ell\) is controllable.
Corollary 4.1 gives the divisibilities
\(\delta_{D'}\mid r_{D'}\mid r_1=\delta_1\) at every prime. The
[Frey residual image calculation](ABC_FREY_JOINT_DIVISION_FIELDS.md#the-image-over-the-hecke-order)
gives absolute irreducibility for \(\ell\ge5\), in particular for
\(\ell>163\).
[Hamidi, Theorems A–B, Corollary 4.1 and Definition 5.2](sources/hamidi-d-new-modular-shimura-degree-2607.23244v1.pdf).

At a multiplicative \(p\ne\ell\) the mod-\(\ell\) Tate extension,
up to unramified quadratic twist, has diagonal characters the
cyclotomic character and \(1\). It is ramified if \(\ell\nmid n_p\).
Indeed the Kummer class of the Tate parameter is ramified in that
case by its valuation. An unramified residual Frobenius with
\(p\not\equiv1\bmod\ell\) has distinct eigenvalues. All possibly
noncontrollable primes in \(D\) therefore belong to
\[
 S_\ell=\{p\mid D:\ \ell\mid \gcd(n_p,p-1)\}.
\]
This set never contains \(p=\ell\). A residue test identifies the
actual noncontrollable subset \(S_\ell^*\). For each \(p\mid D\),
let \(t\) be the letter in \(\{a,b,c\}\) divisible by \(p\). Take
\(s=b\) if \(t=a\) and \(s=a\) otherwise, and put
\[
 t_0=t/p^{v_p(t)},\qquad u_p=t_0/(16s)\pmod p\quad\text{in }\mathbb F_p^\times.
\]
For \(p\ne\ell\) and \(k\ge1\) the local action on \(E[\ell^k]\) is
scalar precisely when
\[
 \ell^k\mid\gcd(n_p,p-1)
           \quad\text{and}\quad
 u_p^{(p-1)/\ell^k}=1\pmod p .
 \tag{CD2a}
\]
Taking \(k=1\) gives \(S_\ell^*\). To prove the test, remove the
unramified quadratic twist of the Tate representation. It is then
scalar precisely when it is trivial, since the quotient character
is trivial. The cyclotomic character is trivial if and only if
\(\ell^k\mid p-1\). With this condition the extension is trivial if
and only if the Tate parameter \(q_p\) is an \(\ell^k\)-th power
in \(\mathbb Q_p^\times\). This holds exactly when its valuation
\(n_p\) is divisible by \(\ell^k\) and its unit residue is an
\(\ell^k\)-th power. The principal units are uniquely
\(\ell^k\)-divisible for \(p\ne\ell\). Finally,
\(q_p=j(E)^{-1}+O(j(E)^{-2})\) and
\[
 j(E)^{-1}=\frac{16a^2b^2c^2}{(16(a^2+ab+b^2))^3}
\]
show \(p^{-n_p}q_p\equiv u_p^2\bmod p\).
For example if \(p\mid c\) then \(b\equiv-a\bmod p\), giving the
same unit square. Squaring is invertible on the
\(\ell\)-primary quotient of \(\mathbb F_p^\times\) for odd
\(\ell\). The scalar criterion is used only for odd \(\ell\).
[Tate, the local curve and its multiplicative uniformization](sources/tate-nonarchimedean-uniformization.pdf).

For \(\ell>163\), remove an even number of members of \(S_\ell^*\)
leaving at most one member, choosing the least member to retain if
\(|S_\ell^*|\) is odd. Denote the resulting divisor by \(D_\ell\).

In fact at a good coefficient prime \(\ell>163\), meaning \(\ell\nmid R\),
one can always use \(D\) itself. In the notation of Böckle–Khare–Manning,
Theorem 7.4 gives equal congruence lengths for the local anemic order
\(T^{\mathrm{st},Q}\) and its cohomology multiplicity module
\(M^{\mathrm{st},Q}(N)\), from their equal defects over the same order.
Here \(Q\) is the set of primes in \(D\), \(N=R/D\), and their
coefficient prime is \(\ell\nmid 2NQ\). The congruence length on
cohomology is the valuation of the optimal Shimura degree, with
valuation normalized in the coefficient ring. This identification
is given in the proof of Theorem 7.4 by the pullback and pushforward
of the optimal quotient. The full localized order contains the
anemic image on this cohomology, so the valuation of its idempotent
denominator is at most that same valuation. The reverse inequality
is \(\delta_D\mid r_D\). Extending the coefficient ring scales both
lengths by its ramification index. For \(D=1\) use
\(r_1=\delta_1\) directly.
[Böckle–Khare–Manning, §7, Theorem 7.4 and its proof](sources/bockle-khare-manning-local-global-wiles-defect-2108.09729v2.pdf).

Define
\[
 \mathcal E=\{\ell>163:\ \ell\mid R,\ |S_\ell^*|\ge2\}.
\]
Take \(D_\ell=D\) for \(\ell\notin\mathcal E\). For \(\ell>163\),
the preceding comparisons give
\[
 v_\ell(H)=v_\ell(r_{D_\ell})-v_\ell(n_5)
                     -\sum_{p\mid D/D_\ell}v_\ell(n_p).
 \tag{CD2}
\]
To see the normalization and the direction of the inequality, the
[uniform degree transfer](ABC_CUBIC_ADJOINT_SHIMURA_PROJECTORS.md#uniform-shimura-transfer)
writes \(\delta_1=\gamma_{D'}\delta_{D'}\prod_{p\mid D'}n_p\), with
\(\gamma_{D'}\) supported on coefficient primes at most \(163\) and
of rational height \(O(\omega(D'))\). Its spectator hypotheses hold
since every \(R/D'\) contains \(M\). Also
\(H=(\delta_D/n_5)\kappa\), where \(\kappa\) has bounded rational
height and the same support by
[(MN1)–(MN2)](ABC_FREY_MONODROMY_NORMS.md#exact-integral-projector-on-the-torus).
Consequently at \(\ell>163\) the right side of (CD2) with any
\(D'\mid D\) of even prime cardinality is at least \(v_\ell(H)\).
For \(D_\ell\) it equals \(v_\ell(H)\).

The support of the extra computations satisfies
\[
  \left(\prod_{\ell\in\mathcal E}\ell\right)^2
          \ \Bigm|\ \prod_{p\mid D}\gcd(n_p,p-1),\qquad
  2|\mathcal E|\le
       \sum_{p\mid D}\omega\bigl(\gcd(n_p,p-1)\bigr).
 \tag{CD3}
\]
Thus at most \(1+|\mathcal E|\) spaces suffice for the exact large
part. One can also use only Hamidi for (CD2) by allowing every
\(\ell>163\) with \(|S_\ell^*|\ge2\) in \(\mathcal E\), whether or not
\(\ell\mid R\).

## Scalar primes occur at arbitrary depth in the fixed image

Fix any finite set of primes \(\ell>163\), and for each \(\ell\)
prescribe a depth \(k\ge1\) and a finite number of distinct primes
\(p\equiv1\bmod\ell^k\) at which to require scalar action.
The primes \(p\) may be chosen disjoint from the coefficient primes
and from all other lists. Such choices can be made repeatedly by
taking a prime factor of \(\Phi_{\ell^k}(t)\) with \(t\) a positive
multiple of \(\ell\) and of all excluded primes. Every prime factor
then has an element of order \(\ell^k\).

There are triples in the image of the fixed \(n=60\) map with each
coefficient prime \(\ell\mid R\) and with \(E[\ell^k]\) scalar at every
prescribed \(p\mid D\). To construct them take the input
\((m,1,m+1)\) with \(m\) even, so \(U=(m-1)^n\) and \(W=(m+1)^n\)
in (FL3) for \(m>1\). Write \(P=U(W-U)W\). At any prescribed
\(p\), let \(d=\ell^k\), choose \(m\equiv p^d m_0\bmod p^{d+1}\)
with \(m_0\ne0\bmod p\), and expand \(W-U\) at \(m=0\). Since
\(p>163\) and \(p\nmid 2n\),
\[
 v_p(P)=d,\qquad P/p^d\equiv 2n m_0\bmod p.
\]
The output satisfies \(v_p(a)=2d\), \(n_p=4d\), and \(b\equiv1\bmod p\).
Thus with \(C_p=27n^2/16\) in \(\mathbb F_p^\times\), the unit in
(CD2a) is \(u_p=C_pm_0^2\). Choose \(m_0=C_p^{(d-1)/2}\), so
\(u_p=C_p^d\). The CRT combines the requirements for \(m\) at
the prescribed primes with \(m\equiv0\bmod2\prod_\ell\ell\).
Each coefficient prime then divides the output letter \(a\).
Also every prescribed \(p\mid a\), and \(p>163\) puts it in \(D\)
for either choice of \(M\). Equation (CD2a) gives scalar action
at the requested depth.

The number of exceptional coefficient primes and the number of
noncontrollable primes over each can both be arbitrarily large.
This construction controls the local scalar conditions. The
difference between the congruence index and the Shimura degree
at such a prime remains a separate comparison. For a finite
check of the construction, \(\ell=167\) and \(p=2339,5011\)
give depth \(1\). Direct integer evaluation of (FL3) verifies the
valuation and unit-power tests at both primes.

## A purely classical integer and the remaining bound

For \(\mathcal D=\{1,D\}\cup\{D_\ell:\ell\in\mathcal E\}\), define
\[
 C_{D'}=\frac{r_{D'}}{n_5\prod_{p\mid D/D'}n_p},\qquad
 G=\gcd_{D'\in\mathcal D} C_{D'},\qquad J=\operatorname{num}(G).
 \tag{CD4}
\]
The gcd of positive rationals takes the minimum of the valuations
at each prime, and the numerator is taken in lowest terms.
By (CD2) and its inequality, \(v_\ell(J)=v_\ell(H)\) for every
\(\ell>163\). The smaller primes can be compared without any
multiplicity-one assertion. At each of them, by (CD4) and the
same degree transfer,
\[
  \min_{D'\in\mathcal D}v_\ell\!\left(\frac{\gamma_D}
                                                {\gamma_{D'}\kappa}\right)
      \le v_\ell(G/H)
      \le v_\ell(\gamma_D/\kappa).
\]
The lower bound uses \(\delta_{D'}\mid r_{D'}\). The upper bound uses
\(1\in\mathcal D\) and \(r_1=\delta_1\). Clipping negative
valuations of \(G\) to form \(J\) cannot increase the absolute
difference from \(v_\ell(H)\). The prime support of the corrections
is fixed, so
\[
 h_{\mathbb Q}(J/H)=O(1+\omega(D)),\qquad
 \log J=\log c+\log R+
           O\!\left(\frac{\log c}{\log(\log c+2)}\right).
 \tag{CD5}
\]
This uses at most one additional space for the small primes.
For an exact reconstruction of \(H\), its small-prime valuations
can instead be supplied by the fixed-level Brandt norm.

Full rational abc is equivalent to \(J\ll_\varepsilon R^{2+\varepsilon}\)
for every \(\varepsilon>0\), either on the whole fixed-level family
or on the images of the single \(n=60\) map in (FL3). This follows
from (CD5) with exponent slack, then (FL5). The gcd in (CD1) gives
integer linear congruences for testing the bound. The remaining
obligation is a uniform estimate for the size of the normalized gcd
in (CD4), including the powers at small primes. The support bound
(CD3) gives no estimate for the depth of the congruences it lists.
A sequence with \(J/R^{2+\varepsilon}\) unbounded for one fixed
\(\varepsilon>0\) would disprove abc.


## A single matrix sufficient for the sharp estimate

The unconditional divisibilities \(\delta_D\mid r_D\mid\delta_1\)
and the uniform degree transfer also give
\[
 \log r_D=\log c+\log R+
           O\!\left(\frac{\log c}{\log(\log c+2)}\right).
 \tag{CD6}
\]
Here the difference between the two bounding degrees is
\(O(\omega(D)+\sum_{p\mid D}\log n_p)\) on the logarithmic scale.
Thus the sharp bound for just \(r_D\) in (CD1), using the original
\(D\)-new space for every triple, is already equivalent to abc.
The adjusted gcd is useful when an exact valuation of \(H\) is
needed at a large prime.

Let \(L=\{u\in\mathbb Z^B:\mathbf K_D^{\mathsf T}u=0\}\).
There is a dual version of (CD1),
\[
                    F^{\mathsf T} L=r_D\mathbb Z .
 \tag{CD7}
\]
Indeed \(\sum_{n\le B}u_nT_n\) is supported on the \(f\)-line
precisely for \(u\in L\), and all elements of the order have such
coefficient representations. A single relation \(u\in L\) with
\(F^{\mathsf T}u\ne0\) gives \(r_D\le |F^{\mathsf T}u|\).
In particular a relation of \(\ell^1\)-norm
\(O_\eta(R^{3/2+\eta})\) for each \(\eta>0\), uniformly over the
family, would prove abc. The Hecke coefficient bound gives
\(\max_{n\le B}|a_n(f)|\le\max_{n\le B}\tau(n)\sqrt n=R^{1/2+o(1)}\),
since \(B=R^{1+o(1)}\). The outstanding step is to find a relation
with small evaluation on \(F\) and prove that evaluation is nonzero.


### Finite arithmetic checks

PARI calculations at level \(37\) gave \(r_1=2\) for both rational
newforms using \(B=6\). For the curve \(y^2+y=x^3-x\), PARI's
modular-degree calculation also returned \(2\). At level \(102\),
using \(B=36\), the \(6\)-new congruence indices for the three
rational newforms of exact conductor \(102\) are \(2,2,2\),
while the full indices are \(2,12,4\). The order of their first
three Fourier coefficients is \((1,1,1),(1,-1,1),(1,-1,-1)\).
For these examples the complementary spans came from the images
of \(T_q-a_q(f)\) for prime \(q\le23\) away from the level.
Saturation by two integer kernels, followed by Smith normal form,
gave the same indices as (CD7). These examples check the lattice
formulas on their stated inputs.
