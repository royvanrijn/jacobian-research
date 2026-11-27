# Cubic ray coordinates and Cassels–Tate pairings in inert-prime twists

This note uses the family, isogeny conventions and exact Selmer calculation
of [MS1](MORDELL_INERT_PRIME_SELMER_MOMENTS.md). Put
\(E_R:y^2=x^3-2R^2\), \(E'_R:y^2=x^3+54R^2\), and
\(S_R=\operatorname{Sel}_\psi(E_R/\mathbf Q)\).
The pairing considered here is the Cassels–Tate pairing on
\(\Sha(E_R)[\psi]\), pulled back to \(S_R\).
It is also the restriction of the full three-Selmer pairing under MS1.

**Result and scope.** There is a complete Kummer basis for every ray space,
a general pairwise variation formula, and an implemented exact full-matrix
evaluator from norm certificates, including the coordinate at 3.
Seven full local-symbol certificates are retained. They give new
nondegenerate matrices at \(R=71\) and \(355\), proving rank zero
and \(\Sha(E_R)[3]\cong(\mathbf Z/3\mathbf Z)^2\).
The explicit 5–23 entry includes its correction at 3 and
holds whenever \(115\mid R\), without a restriction on \(R/115\bmod9\).
There is also an explicit entry between the wild character and the 5-ray
character. The 5–23 entry is nonzero for two thirds
of the allowed moving primes \(q\equiv-1\pmod9\) in \(R=115q\).
Thus infinitely many admissible curves have rational rank at most one and
\(\dim_{\mathbf F_3}\Sha(E_R)[3]\ge2\).
Complete normalized matrices are obtained below for selected products of
two, three and four allowed primes, including spaces with the wild
coordinate. The cubic local descent has the explicit scalar formula
(24), including in split cubic algebras; direct local certificates
also verify complete tame and wild controls. The full evaluator computes
every off-diagonal entry without a rational-rank or rational-point
oracle. A uniform runtime guarantee and a Frobenius-only rule with all
pairwise reference constants already evaluated are not claimed.
These are internal written proofs with an exact algebra checker; no
independent mathematical review or formal verification is claimed.

## Pairing, kernel and literature match

For an isogeny \(\phi:A\to B\), the pairing on
\(\operatorname{Sel}_{\hat\phi}(B)\) has radical equal to the image of
\(\operatorname{Sel}_3(A)\).
[van Beek–Fisher, introduction, diagram preceding (2), and Definition 1.2](https://arxiv.org/pdf/1711.02432v1).
Apply this with \(A=E'_R\), \(B=E_R\), \(\phi=\hat\psi\).
In particular, \(\operatorname{Sel}_{\hat\psi}(E'_R)=0\) does not imply
that the pairing on \(S_R\) vanishes. Its radical is controlled by the
**full** three-Selmer group of \(E'_R\). The kernel theorem concerns this
isogeny-Selmer pairing; it must not be replaced by a same-curve
three-Selmer surjectivity assertion.

The closest directly applicable comparison is
[Koymans–Smith, Proposition 4.1, equations (4.1)–(4.8) and its proof](https://arxiv.org/pdf/2405.09311v1).
Take their \(d=54\). Then \(E_{d,n}=E'_n\), while
\(E_{-27d,n}\cong E_n\) through \((x,y)\mapsto(x/9,y/27)\).
Their coefficient line \(M_{-3d}\) is \(M_{-2}\), the sign line of
\(k=\mathbf Q(\sqrt{-2})\). Their pairing \(\mathrm{CTP}_{54,n}\)
is precisely the pairing needed here, up to fixed isogeny/Weil-pairing
sign conventions. The proof below specializes their cochain comparison
and gives its governing extension explicitly.

[Wu–Zhu, Lemma 6.1, Lemma 6.2 and Corollary 6.3](https://arxiv.org/pdf/2611.29541v1)
give the relevant governing-field construction and explain why inert
moving primes can also contribute central Frobenius coordinates. Their
Theorems A and B concern the Fermat cubic, whose coefficient line is
\(\mu_3\), rather than \(M_{-2}\). Their six-prime construction uses two
moving primes split in \(\mathbf Q(\zeta_3)\), excluded here. Those
theorems are not applied to this family. The argument here is a
specialization of the Koymans–Smith comparison method, not a claim that
this method or governing-field variation is new.

## Coordinates and an exact norm solution

Let \(L=\mathbf Q(s)\), \(s^2=6\). Fix the usual quadratic-descent
identification
\[
 H^1(\mathbf Q,M_{-2})
   \cong \ker\bigl(N:L^\times/L^{\times3}
                         \longrightarrow\mathbf Q^\times/\mathbf Q^{\times3}\bigr).
 \tag{1}
\]
It is obtained by restricting to \(L\), where \(M_{-2}\cong\mu_3\),
and taking the anti-invariant Kummer subspace. Evaluation of a rational
point \((x,y)\in E'_R(\mathbf Q)\) is \([y-3Rs]\).
This is the same forward Kummer map used in
[MH1](MORDELL_SELMER_LARGE_SUPPORT_HEIGHT.md).

Use the following oriented ray characters, specified by their Kummer
representatives in (1):

| Character | Representative \(a_p\) | \(a_p/p\) |
|---|---|---|
| \(\chi_5\) | \(-35+15s\) | \(-7+3s\) |
| \(\chi_{23}\) | \(667+276s\) | \(29+12s\) |
| \(\chi_{53}\) | \(53+159s\) | \(1+3s\) |
| \(\chi_{359}\) | \(-25489+10770s\) | \(-71+30s\) |
| \(\chi_{503}\) | \(9557+6036s\) | \(19+12s\) |

Each last-column norm is \(-p\); hence \(N(a_p)=(-p)^3\), and the
valuations at the two primes over \(p\) are 1 and 2. There is no other
valuation support modulo three. If \(a_p=A+Bs\), then
\(A\equiv\pm1\pmod9\), \(B\equiv0\pmod3\). At 3 the field
\(k\) splits and \(L_3=\mathbf Q_3(\zeta_3)\). Consequently
\(a_p\in\pm(1+\pi^3\mathcal O_{L_3})\), where
\(\pi=\sqrt{-3}\). These are unramified cubic Kummer classes:
substituting \(1+\pi T\) into \(X^3=1+\pi^3 c\) gives
\(T^3-\pi T^2-T=c\), whose reduction has derivative \(-1\) and is
either split or irreducible of degree three over \(\mathbf F_3\).
Thus the corresponding characters of \(k\) are unramified at 3.
The class-number-one ray calculation in MS1 identifies these with the
individual tame ray coordinates. In particular they are independent
whenever their primes divide \(R\), and they are the entire basis if
\(\epsilon_R=0\) and the displayed primes exhaust its support.

Set
\[
 a=-35+15s,\qquad b=667+276s,\qquad r^3=a,
\]
and write
\[
 X(s,r)=\left(-\frac75s-\frac{64}{15}\right)r^2
         +\left(\frac13s+\frac{13}{3}\right)r
         +\frac{14}{3}s-\frac83.
\]
There is an exact norm identity
\[
             N_{L(r)/L}\bigl(X(s,r)\bigr)=b.                \tag{2}
\]
For \(A+Br+Cr^2\) this norm is
\(A^3+aB^3+a^2C^3-3aABC\); substitution verifies (2) directly.
The identity, rather than trust in the norm solver that found it, is
the retained computational input.

## Explicit central symbol

For an allowed prime \(q\notin\{5,23\}\), choose
\(s_q\in\mathbf F_q\) with \(s_q^2=6\), and work in
\(\mathbf F_{q^2}=\mathbf F_q[\zeta]\), \(\zeta^2+\zeta+1=0\).
Both choices exist: \((6/q)=1\) and \(q\equiv2\pmod3\).
Let \(r_q\) be the unique cube root of \(-35+15s_q\) in
\(\mathbf F_q\). Define
\[
 \begin{split}
 T(s,r)&=X(s,\zeta r)X(s,\zeta^2r)^2,\\
 D_q&=T(s_q,r_q)T(-s_q,-5/r_q),\\
 D_q^{(q^2-1)/3}&=\zeta^{\kappa(q)},\qquad
                   \kappa(q)\in\mathbf F_3.              \tag{3}
 \end{split}
\]
All denominators and factors in (3) are nonzero: the only rational
primes in denominators of \(X\) are 3 and 5, and its norm is supported
at 23. The symbol is unchanged on replacing \(s_q\) by \(-s_q\) or
\(\zeta\) by \(\zeta^2\), as also follows from the governing-field
description below. Exact values needed here are
\[
 \kappa(53)=2,\quad \kappa(71)=1,\quad
 \kappa(359)=2,\quad \kappa(503)=0.                       \tag{4}
\]

**Proposition MCT1 (fixed-block formula).** Let \(t\) be a squarefree
product of allowed primes different from 5 and 23, with
\(t\equiv\pm1\pmod9\), including \(t=1\). Put \(R=115t\).
There is one fixed nonzero scalar \(u\in\mathbf F_3^\times\),
independent of \(t\), such that
\[
 3\langle\chi_5,\chi_{23}\rangle_{R}
       =u\left(\sum_{q\mid t}\kappa(q)-\kappa(503)\right).
                                                               \tag{5}
\]
Here \(3\langle\ ,\ \rangle\) means the identification
\(j/3\mapsto j\), not multiplication inside \(\mathbf Q/\mathbf Z\).
Define the normalized pairing \(B_R=u^{-1}3\langle\ ,\ \rangle_R\).
All matrices in this note use this **single common normalization**.
The usual convention differs by at most an overall sign; the remaining
Weil-sign calibration is not asserted to have been computed.

### The governing extension and its descent

Take \(F=L(\zeta)\), \(v^3=b\),
\(\xi=X(s,r)\), and let \(\sigma r=\zeta r\) fix \(F,v\).
Put
\[
 \theta=\sigma\xi(\sigma^2\xi)^2,\qquad
 \tau(s)=-s,\quad\tau(r)=-5/r,\quad\tau(v)=-23/v,
 \quad\tau(\zeta)=\zeta,
 \quad \eta=\theta\tau\theta.
\]
The norm identities \(a\bar a=(-5)^3\), \(b\bar b=(-23)^3\)
make \(\tau\) an involution of \(F(r,v)\).
Let \(\mathcal H=F(r,v,w)\), \(w^3=\eta\).

The radicands \(a,b\) are independent over \(F\) modulo cubes by
their separate valuations at 5 and 23. Thus \([F(r,v):F]=9\).
Equation (2) gives
\[
 \frac{\sigma\theta}{\theta}=\frac{\xi^3}{b},\qquad
 \frac{\sigma\eta}{\eta}
   =\left(\frac{-23\xi}{v^2\tau\sigma^2\xi}\right)^3.
 \tag{6}
\]
Consequently \(\sigma\) extends by
\(\sigma(w)=Aw\), \(A=-23\xi/(v^2\tau\sigma^2\xi)\).
The norm of \(A\) under \(\sigma\) is one. The automorphism
\(\delta(v)=\zeta v\) fixing \(F,r\) extends by \(\delta(w)=w\).
Since \(\delta(A)=\zeta A\), the two extensions have a nontrivial
central commutator. In particular \(\eta\) cannot be a cube in
\(F(r,v)\): otherwise the two commuting automorphisms of that field,
applied to a hypothetical cube root, would give this contradictory
commutator even after multiplying either lift by a cube root of unity.
Thus \(\operatorname{Gal}(\mathcal H/F)\) is the Heisenberg group
\(P\) of order 27 and exponent three, with center generated by
\(z(w)=\zeta w\).

Extend \(\tau\) by \(\tau(w)=w\). Let \(c\) fix \(s,r,v\) and
invert \(\zeta\), and put
\[
 C=(\sigma\xi\,\sigma^2\xi)
       (\tau\sigma\xi\,\tau\sigma^2\xi),\qquad c(w)=C/w.
\]
Indeed \(\eta c(\eta)=C^3\). Both \(\tau\) and \(c\) are
involutions and commute, since \(\tau(C)=c(C)=C\).
Each acts by negation on \(P/Z(P)\) and fixes \(Z(P)\).
Their product acts trivially on this quotient and its center. The
automorphisms of \(P\) acting trivially on both form a group of order
nine (a lift of either quotient generator can only be multiplied by a
central element). This group has odd order, so the involution
\(\tau c\) centralizes \(P\). Therefore
\[
 H=\mathcal H^{\langle\tau c\rangle},\qquad
 \operatorname{Gal}(H/\mathbf Q)=P\rtimes\langle c\rangle,
 \qquad H^{P}=k.                                            \tag{7}
\]
The extension is unramified outside \(2,3,5,23\): this follows from
the norms and denominators of \(a,b,\xi\), and the elementary tame
Kummer ramification criterion outside 3.

The two quotient coordinates in \(P/Z(P)\) are precisely
\(\chi_5,\chi_{23}\). In Heisenberg coordinates there is a cochain
\(h:G_{\mathbf Q}\to\mathbf F_3\), unramified outside those four
primes, with \(dh\) a fixed nonzero multiple of
\(\chi_5\smile\chi_{23}\). One obtains this directly by taking the
central coordinate in the multiplication law
\((x,y,z)(x',y',z')=(x+x',y+y',z+z'+xy')\).
The involution acts by \((x,y,z)\mapsto(-x,-y,z)\).

At an allowed prime \(q\), Frobenius can be conjugated to \(z^j c\):
conjugating by a quotient vector removes its two noncentral
coordinates, since 2 is invertible modulo three. Its square is the
central element \(z^{2j}\). Choosing the roots \(r,v\) in
\(\mathbf F_q\), its action on \(w\) after squaring is exactly the
power in (3). Thus \(\kappa(q)=2j\). This also proves the choice
independence asserted after (3).

### Comparison of twists

Use the reference parameter \(R_0=115\cdot503\). Its entire Selmer
space is generated by rational Kummer classes, as the exact witnesses
below prove, so its pairing is zero.

Put \(a_t=\chi_{3,t/503}\in H^1(\mathbf Q,\mu_3)\). At each of
\(2,3,5,23,\infty\) it is zero: units are cubes at 2, 5 and 23,
and \(t/503\equiv\pm1\pmod9\) at 3. Hence local isomorphisms of
the two cubic twists identify their full local Kummer conditions there.
At every prime dividing \(503t\) outside that fixed set, both
\(\chi_5\) and \(\chi_{23}\) vanish locally. Indeed they are
unramified, and the unramified cohomology of the inert sign line is
\(M_{-2}/(\mathrm{Frob}-1)M_{-2}=0\). Choose embeddings making
their cocycles zero there and use zero local full lifts.

These facts check all local conditions in Koymans–Smith Definition 2.2.
At fixed places the first cocycle is zero; at varying places the two
ray cocycles are zero. Every permuted pairwise cup product is locally
trivial. Outside these places all classes are unramified, and at 3 the
two ray classes are unramified, so their cup product is also zero.

Now apply the proof of their Proposition 4.1, with \(d=54\) and
twists \(115t\), \(115\cdot503\). For a single moving prime this
is exactly its stated specialization with \(m=115\), \(k=1\).
For several moving primes the same proof applies: its only uses of
\(p,q\) are to identify local conditions on the fixed set and use
zero local lifts on the varying set, which were just checked at every
prime of \(503t\). For completeness, the comparison of the two
three-torsion extensions, their equations (4.3)–(4.4), changes the
differential of the global two-cochain by a fixed nonzero scalar times
\(a_t\smile\chi_5\smile\chi_{23}\), since
\((503/t)^2\equiv t/503\) modulo cubes. A primitive for this change
is the corresponding scalar times \(-a_t\smile h\), using
\(d(a_t\smile h)=-a_t\smile dh\).
The ambiguity is a global two-cocycle and has reciprocity sum zero.

The local-lift terms cancel at fixed places. At varying places the
cochain \(h\) is an unramified character, so the local invariant is
its Frobenius value times \(v_q(t/503)\), up to the one common
coefficient/Weil sign. By (7) this value is a fixed nonzero multiple
of \(\kappa(q)\). Summing proves (5). No equality of local symbols
at different primes or unproved finiteness of \(\Sha\) is used.

## Infinite admissible family

**Theorem MCT2.** Among primes satisfying
\[
 q\equiv-1\pmod9,\qquad (-2/q)=-1,
\]
each value of \(\kappa(q)\) occurs with relative natural density
\(1/3\). For the two thirds with \(\kappa(q)\ne0\), the family
\(R=115q\) satisfies
\[
 \operatorname{rank}E_R(\mathbf Q)\le1,\qquad
 \dim_{\mathbf F_3}\Sha(E_R/\mathbf Q)[3]\ge2.              \tag{8}
\]
The four nonzero lines in \(\langle\chi_5,\chi_{23}\rangle\)
have no rational representative at any height.

To prove the density assertion, the commutator subgroup of
\(P\rtimes\langle c\rangle\) contains \(P/Z(P)\), since \(c\)
acts by negation, and contains \(Z(P)\) by commutators within \(P\).
Its maximal abelian subextension is therefore exactly \(k\).
It follows that \(H\cap\mathbf Q(\zeta_9)=\mathbf Q\), since
the only quadratic subfield of \(\mathbf Q(\zeta_9)\) is
\(\mathbf Q(\sqrt{-3})\), different from \(k\).
In the compositum, of degree \(54\cdot6=324\), prescribe
Frobenius in the conjugacy class of \((z^jc,-1\bmod9)\).
For each \(j\) the class has size nine.
[Chebotarev, in Lenstra's Theorem 3.1, pp. 3–4](sources/lenstra-chebotarev-density-notes.pdf), gives rational
prime density \(9/324=1/36\) for each central value, hence total
density \(1/18\) for the nonzero values. These primes are all
allowed. The assertion is infinitude with an exact prime density,
not an inference from the finite symbol controls.

For these curves \(R\equiv2\pmod9\), so MS1 gives
\(\dim S_R=3\). Equation (5), with (4), gives a nondegenerate
two-by-two block, so the alternating pairing has rank two. Rational
Kummer classes lie in its radical. Since \(E_R(\mathbf Q)[3]=0\),
the Kummer exact sequence then gives (8). Every nonzero vector in
that two-dimensional block pairs nontrivially with another vector,
so none is a rational Kummer class. These are genuine global
obstructions, beyond the height cutoff of MH1.

One of these obstructions has the particularly simple plane-cubic model
\[
 C_q:\quad U^3-7U^2V+18UV^2-14V^3+23qW^3=0.             \tag{8a}
\]
For every prime in the nonzero-symbol family, \(C_q\) has points over
every completion of \(\mathbf Q\) and has no rational point. Indeed
the coefficient of \(s\) in
\(a_5(U+Vs)^3\) is
\(15(U^3-7U^2V+18UV^2-14V^3)\). Equation (8a) is therefore the
\(\chi_5\)-covering defined by
\(a_5(U+Vs)^3=(y-3Rs)W^3\), with \(R=115q\).
Its map to \(E'_R\), on \(W\ne0\), has
\(x=-5(U^2-6V^2)/W^2\) and \(y\) the rational coefficient of
\(a_5(U+Vs)^3/W^3\). The norm identity proves the point equation.
There is no rational point at \(W=0\): that would make the nonzero
Kummer class \(a_5\) rational modulo cubes, and its norm condition
would then make it a cube. The binary cubic has discriminant \(-200\),
so the plane cubic is smooth. Local solubility follows from its Selmer
membership; global insolubility follows from the nonzero pairing.
For \(q=53\), this gives the explicit Hasse failure
\[
       U^3-7U^2V+18UV^2-14V^3+1219W^3=0.
\]

More generally a pairing of rank \(2r\) on the MS1 ray space gives
\(\dim\Sha(E_R)[3]\ge2r\) and
\(\operatorname{rank}E_R(\mathbf Q)\le\omega(R)+\epsilon_R-2r\),
without assuming \(\Sha\) finite.

## Complete exact matrices and rational controls

In the ordered bases indicated by the prime factors, the normalized
matrices are
\[
 \begin{array}{c|c|c|c}
 R & \text{basis primes} & B_R & \operatorname{rank}E_R(\mathbf Q)\\ \hline
 115 &(5,23)&0_{2\times2}&2\\
 57845 &(5,23,503)&0_{3\times3}&3\\
 6095 &(5,23,53)&
 \begin{pmatrix}0&2&0\\1&0&1\\0&2&0\end{pmatrix}&1\\
 2188105 &(5,23,53,359)&
 \begin{pmatrix}0&1&1&0\\2&0&1&2\\2&2&0&2\\0&1&1&0\end{pmatrix}&2
 \end{array}                                                \tag{9}
\]
The last two curves have
\(\Sha(E_R)[3]\cong(\mathbf Z/3\mathbf Z)^2\).
These rank and Sha conclusions use (5) plus the exact rational
certificates, not the exploratory PARI rank outputs.

The certificates are points on \(E'_R\). The vector \(v\) in each
row means the exact identity
\[
       y-3Rs=d^3\prod_{p\mid R}a_p^{v_p},                 \tag{10}
\]
with \(d\in L\). All coordinates and cube roots \(d\) are retained
in [the exact input](mordell_ray_pairing_exact.json).

| \(R\) | Point \((x,y)\) on \(E'_R\) | Ray vector |
|---|---|---|
| 115 | \((-5,845)\) | \((2,0)\) |
| 115 | \((-69,621)\) | \((0,1)\) |
| 57845 | \((-1509,421011)\) | \((0,0,1)\) |
| 57845 | \((7545/4,3463155/8)\) | \((2,0,2)\) |
| 57845 | \((11569/4,3621097/8)\) | \((0,2,1)\) |
| 6095 | \((795,50085)\) | \((2,0,1)\) |
| 2188105 | \((-26925,15460335)\) | \((1,0,0,2)\) |
| 2188105 | \((530265/16,1099129635/64)\) | \((2,1,2,0)\) |

For the two zero controls these vectors span the entire ray space.
For 6095, the block entry is \(\kappa(53)=2\); its single radical
vector determines the other entries of an alternating three-by-three
matrix. For 2188105, the block entry is
\(\kappa(53)+\kappa(359)=1\); the two independent radical vectors
uniquely determine the alternating four-by-four matrix. Their rational
dimensions equal the resulting radical dimensions. All these curves
have \(\epsilon_R=0\). There is no assertion that a zero pairing in
general implies rational solubility.

## General coordinates and pairwise variation

Write \(\varepsilon=5+2s\), and denote the wild character by
\(\chi_0=[\varepsilon]\) in (1). The index 0 denotes its conductor
at 3, not a rational prime factor of \(R\).

**Proposition MCT3 (basis for every admissible parameter).** For each
allowed prime \(p\), choose \(\pi_p=A_p+B_ps\) of norm \(-p\),
and choose \(e_p\in\{0,1,2\}\) so that the coefficient of \(s\) in
\[
             a_p=p\pi_p\varepsilon^{e_p}
\]
is divisible by three. The classes \(\chi_p=[a_p]\), for \(p\mid R\),
together with \(\chi_0\) when \(\epsilon_R=1\), are a basis of
\(S_R\). Changing the generator of the same prime ideal changes the
normalized representative only by a cube. Interchanging the two primes
over \(p\) reverses that basis vector.

Here is an existence and uniqueness check. The quadratic field of
discriminant 24 has class number one: Minkowski's bound is
\(\sqrt{24}/2<3\), and its prime ideal of norm 2 is principal,
generated by \(2+s\). Thus a prime over \(p\) has a generator of
norm \(\pm p\). Norm \(+p\) is impossible modulo 3, so the norm is
\(-p\). The units are \(\pm\varepsilon^{\mathbf Z}\): a unit of
norm \(-1\) is impossible modulo 3, and the least positive Pell unit
greater than one is \(5+2s\). If \(p\pi_p=A+Bs\), then
\(A\ne0\pmod3\) and multiplication by \(\varepsilon^e\) changes
\(B\) modulo 3 to \(2^e(B+eA)\). There is exactly one suitable
\(e\pmod3\). Its norm condition then forces \(A\equiv\pm1\pmod9\).
The local argument following (1) proves that this class is unramified
at 3. It has just the one tame ray coordinate at \(p\).
The unit \(\varepsilon\) has no tame support and is ramified at 3:
its linear coefficient is nonzero modulo 3, whereas a unit cube has
linear coefficient divisible by three. MS1 now supplies the stated
dimensions and proves that these independent directions exhaust the
Selmer space. This also specifies the extra wild basis direction.

### Normalizing the pairwise governing field

For any two distinct basis classes \(f,g\), their cup product with
values in the trivial module \(\mathbf F_3\) vanishes globally.
At each allowed tame prime, at least one class is unramified and hence
locally zero in the inert sign line. Away from the support the two
classes are unramified. At 2 they are locally zero. At 3,
\(H^2(\mathbf Q_3,\mathbf F_3)=0\), by local duality and
\(\mu_3(\mathbf Q_3)=0\), even when one class is wild.
Restriction to \(\mathbf Q(\zeta_3)\), the Brauer local–global
injection and corestriction by two therefore prove global vanishing.

In particular there is an exact norm solution
\(N_{L(\sqrt[3]{a_f})/L}\xi=a_g\).
Over \(F=L(\zeta_3)\) this follows from the split cubic symbol algebra.
It descends to \(L\) without a Hasse-norm assumption: if a solution
over \(F\) is \(x\), its norm to \(L(\sqrt[3]{a_f})\) has norm
\(a_g^2\); dividing \(a_g\) by that norm gives a solution of norm
\(a_g\).

If \(N_{L/\mathbf Q}(a_f)=c_f^3\) and
\(N_{L/\mathbf Q}(a_g)=c_g^3\), repeat (6) with
\(\tau r=c_f/r\), \(\tau v=c_g/v\), and
\(A=c_g\xi/(v^2\tau\sigma^2\xi)\). The commutator is the same
for every ordered pair: \(\delta(A)=\zeta A\).
This fixes the relative normalization of all pairwise central symbols;
there is only the one overall Weil sign already present in (5).

A particular norm solution may introduce auxiliary central ramification.
It must be removed before using an unrestricted prime-sum formula.
Outside \(2,3\) and the supports of \(f,g\), inertia is central and
tame. Frobenius raises tame inertia to the rational-prime power, so
nonzero central inertia can occur only at \(\ell\equiv1\pmod3\).
Adding a suitable cubic Dirichlet character of conductor \(\ell\)
to the primitive cochain removes it. Only finitely many such primes
occur. In a residue formula this operation adds the corresponding
cubic character of \(q\); it is not legitimate simply to discard
the auxiliary primes of a norm solution.

For two tame classes one can also make the governing extension
unramified at 3. Their quotient coordinates there factor through the
same unramified cyclic line. An unramified lift to the exponent-three
Heisenberg group exists, and the difference from any given lift on
inertia is an order-three character of \(\mathbf Q_3^\times\).
Its ramified part is one-dimensional. The global cubic character
\(\chi_9\) of conductor 9 spans it, so adding exactly one multiple
of \(\chi_9\) removes this ramification. Explicitly,
\[
 \chi_9(2)=1,\qquad
 \chi_9(n)=\begin{cases}0&n\equiv\pm1\pmod9,\\
                         1&n\equiv\pm2\pmod9,\\
                         2&n\equiv\pm4\pmod9.
           \end{cases}                                      \tag{11}
\]
There is no remaining global-character ambiguity after this normalization:
an order-three rational character unramified outside 2 and the two
allowed primes, and unramified at 3, is unramified everywhere and zero.

Let \(\widetilde\kappa_{p,r}(q)\) be the squared-Frobenius central
coordinate of this normalized field at an allowed \(q\ne p,r\).
At a regular prime of the norm solution it is given by (3), with the
appropriate radicands, norm solution and the Dirichlet corrections just
described. At an exceptional allowed prime, use its local Kummer
representative: its valuation is divisible by three, so divide by a
cube of a uniformizer before applying the residue power. This is a
local computation, not substitution into a denominator that vanishes.

**Proposition MCT4 (general pairwise variation).** In the common
normalization of (5), for two distinct tame coordinates \(p,r\mid R\),
\[
 B_R(\chi_p,\chi_r)
   =B_{pr}(\chi_p,\chi_r)
       +\sum_{q\mid R/(pr)}\widetilde\kappa_{p,r}(q).       \tag{12}
\]
There is no condition on \(R/(pr)\bmod9\), and the formula also
holds when \(\epsilon_R=1\).
For a wild–tame pair \((\chi_0,\chi_p)\), choose any admissible
reference \(R_0\) containing \(p\) with \(\epsilon_{R_0}=1\).
Normalize away auxiliary ramification outside \(2,3,p\), retaining
the necessary ramification at 3. Then, whenever \(p\mid R\) and
\(\epsilon_R=1\),
\[
 B_R(\chi_0,\chi_p)=B_{R_0}(\chi_0,\chi_p)
    +\sum_{q\mid R/p}\kappa_{0,p}(q)
    -\sum_{q\mid R_0/p}\kappa_{0,p}(q).                  \tag{13}
\]
The remaining freedom to add \(\chi_9\) does not change (13), since
\(\chi_9(R/R_0)=0\). The base values in (12) and (13) are genuine
inputs; this proposition does not assume they vanish.

To prove (12), adapt the cochain comparison of MCT1 after normalizing
at 3. The only new point is the lack of a local cube ratio there.
The restrictions of two tame characters to \(G_3\) are proportional.
By alternation one may replace one argument by a linear combination
that is zero on \(G_3\), or interchange arguments if necessary.
The corresponding local full-lift term vanishes for both twists, even
if their local curves are not isomorphic. The primitive changes by a
multiple of the diagonal primitive \(-f^2/2\), whose coboundary is
\(f\smile f\). That primitive is unramified at 3 and zero at the
moving primes where \(f\) is zero. The normalized cochain at 3 is
therefore an unramified character, and it pairs to zero with the
unit twist ratio. All other fixed-place comparisons are by local
cubes, as in MCT1; the remaining terms are exactly (12).
For (13), the twist ratio is already a cube at 3 because both parameters
are \(\pm1\pmod9\), so the original comparison proof applies.
The universal commutator calculation above gives the same scalar in
both formulas. No per-pair sign is fitted from the examples.

### The explicit correction for 5 and 23

**Theorem MCT5 (5–23 entry for every admissible multiple of 115).**
The original norm solution (2) needs exactly the correction
\[
 \widetilde\kappa_{5,23}(q)=\kappa(q)+\chi_9(q),\qquad
 B_R(\chi_5,\chi_{23})=
       \sum_{q\mid R/115}\bigl(\kappa(q)+\chi_9(q)\bigr)
                                                        \tag{14}
\]
for every admissible \(R\) divisible by 115.

There are no auxiliary primes outside \(2,3,5,23\) in (2). Thus
the normalization argument leaves just one unknown coefficient
\(m\in\mathbf F_3\) of \(\chi_9\). A central twist of the
Kummer radicand from \(\eta\) to \(\eta\zeta^m\) changes the
regular-prime residue symbol by
\[
             m(q^2-1)/3=m\chi_9(q)\pmod3.                \tag{15}
\]
The rational-control curve \(R=115\cdot173=19895\) has
\(\epsilon_R=0\). Its three retained rational Kummer vectors,
in the basis \((5,23,173)\), are
\((0,1,2),(1,2,1),(1,2,2)\), which span the whole space.
Its pairing is zero. Since \(\kappa(173)=2\) and
\(\chi_9(173)=1\), formula (12) and the zero base pairing at 115
force \(m=1\). This determines the correction from an exact rational
certificate, independently of all nonzero matrices. It is not a fit
to an assumed rank or a bounded point search. Equations (12) and
(15) prove (14), including parameters with a wild coordinate.

### An explicit wild–tame entry

There is a second small norm identity. Put \(r^3=\varepsilon\) and
\[
 X_0(s,r)=(25s-53)r^2+(10s-28)r-s-28.
\]
Then
\[
           N_{L(r)/L}X_0(s,r)=-35+15s=a_5.               \tag{16}
\]
For allowed \(q\ne5\), define
\[
 T_0(s,r)=X_0(s,\zeta r)X_0(s,\zeta^2r)^2,\qquad
 D_{0,q}=T_0(s_q,r_q)T_0(-s_q,1/r_q),
 \quad r_q^3=5+2s_q,
\]
and let \(D_{0,q}^{(q^2-1)/3}=\zeta^{\kappa_{0,5}(q)}\).
The coefficients are integral and the norm is supported at 5, so
there are no auxiliary or exceptional allowed primes in this formula.
Its governing field is unramified outside \(2,3,5\).

**Theorem MCT6 (wild–5 entry).** For every admissible \(R\) divisible
by 5 with \(\epsilon_R=1\),
\[
 B_R(\chi_0,\chi_5)
       =\sum_{q\mid R/5}\kappa_{0,5}(q)-2.              \tag{17}
\]
Indeed \(\kappa_{0,5}(47)=2\), and the reference curve
\(R_0=235=5\cdot47\) has full rational Kummer span, with vectors
\((2,0,1),(1,2,1),(0,1,2)\) in the basis \((0,5,47)\).
Its pairing is zero. Apply (13). In particular this formula includes
the extra coordinate at 3 without treating it as a tame prime.

## A local-symbol reduction for every matrix entry

Equations (12)–(13) reduce variation to fixed pairwise fields, but
their absolute reference constants still have to be supplied. The
following specialization of van Beek–Fisher gives an alternative
finite local reduction, including every wild entry. It specifies the
remaining lifting data rather than claiming those data have been
computed uniformly.

**Proposition MCT7 (full local reduction).** Let \(R>1\), let
\(a_i\) be the basis of MCT3, put \(\rho^3=R\), and work over
\(L=\mathbf Q(\sqrt6)\). Set \(d=\sqrt{-3}\), and on \(E'_R\)
take
\[
 S=(0,3Rs),\qquad T=(-6\rho^2,3Rsd),\qquad
 \sigma(\rho)=\zeta^2\rho.
\]
Then \(\sigma T=S+T\). For each \(i\), solve
\[
 N_{L(\rho)/L}\xi_i=a_i,\qquad
 b_i=\sigma(\xi_i)^2\sigma^2(\xi_i).                    \tag{18}
\]
Such solutions exist because the classes are Selmer classes, by
van Beek–Fisher Lemma 1.1 and Theorem 2.11.
For each place \(w\) of \(L\), choose
\(P_{i,w}=(X_{i,w},Y_{i,w})\in E'_R(L_w)\) with
\([Y_{i,w}-3Rs]=[a_i]\), avoiding the zeros of the tangent functions.
Put
\[
 t_T(P)=Y+sd\rho X+3Rsd.
\]
The local difference class has a representative
\(c_{i,w}\in L_w(\zeta)^\times/L_w(\zeta)^{\times3}\) whose
restriction to \(L_w(\zeta,\rho)\) is
\[
                    [b_i/t_T(P_{i,w})].                  \tag{19}
\]
If \(\zeta_W\) is the Weil-compatible root in van Beek–Fisher
Lemma 1.3, then the usual \(j/3\mapsto j\) pairing is
\[
 3\langle\chi_i,\chi_j\rangle_{\mathbf Q}
   =2\sum_w [L_w(\zeta):L_w]^{-1}
       \operatorname{Ind}_{\zeta_W}(c_{i,w},a_j)_w
                  \quad\text{in }\mathbf F_3.            \tag{20}
\]
The root \(\zeta_W\) and its relation to the chosen residue-symbol
normalization must be retained; the standard Weil-sign calibration
is still separate. In (20), \(\operatorname{Ind}_{\zeta_W}\) denotes
the exponent in \(\mathbf F_3\), rather than its value divided by
three. The case \(R=1\) has a one-dimensional Selmer
space and zero alternating matrix.

Here (18) is Theorem 2.11 for the reverse isogeny over \(L\), whose
kernel point \(S\) is rational. Its quotient kernel is \(\mu_3\),
and its three-torsion field is \(L(\zeta,\sqrt[3]R)\).
Direct addition verifies the displayed \(\sigma T\); the tangent
slope at \(T\) is \(-sd\rho\), giving \(t_T\).
Definition 1.2 and Lemmas 1.3/2.10 then give the local formula over
\(L\). Restriction to \(L\) multiplies Cassels–Tate values by two,
whose inverse in \(\mathbf F_3\) is also two, giving (20).
There is no division by a degree divisible by three.

The descent in (19) is essential: taking a Hilbert symbol directly
in the cubic extension would multiply its invariant by three and
erase it. Its existence follows from exactness of the local
cohomology sequence. Any ambiguity from restriction to the cubic
extension is a multiple of \([R]\); it pairs trivially with
\(a_j\) because the same local obstruction vanishes for the Selmer
class \(a_j\). Thus (20) is independent of this choice.

Only finitely many finite places are needed: those over \(6R\)
and the rational denominators of the coefficients of the \(\xi_i\)
in the basis \(1,\rho,\rho^2\). Outside them the norm identities
make all conjugates of \(\xi_i\) units, the local full Kummer
conditions are unramified, and the local invariant is zero.
At 3, \(L_3=\mathbf Q_3(\zeta)\), and the local symbol is the
explicit four-dimensional calculation from van Beek–Fisher (9).
With basis \((\lambda,1-\lambda,1-\lambda^2,1-\lambda^3)\),
\(\lambda=1-\zeta_W\), its exponent matrix is
\[
 \begin{pmatrix}
 0&0&0&2\\0&0&1&0\\0&2&0&0\\1&0&0&0
 \end{pmatrix}.                                          \tag{21}
\]
This includes the wild coordinate. Producing the local points and
their coordinates in (21) is still a required computational step.
MCT8 below gives an explicit scalar representative for (19). The
pairwise residue checker does not perform these local computations;
the direct checker below verifies its own retained certificates.

## An explicit scalar for the cubic local descent

The local class in (19) can be computed by arithmetic in its base field.
This is a Hilbert–90 realization of the difference of the two full
Kummer lifts in van Beek–Fisher Definition 1.2 and Theorem 2.11. It
also applies when the local cubic algebra splits. It does not require
restriction of the local invariant to an extension of degree three.

**Proposition MCT8 (explicit local scalar descent).** In MCT7, fix a
place, put \(K=L_w(\zeta)\), and use the cyclic étale algebra
\(D=K[\rho]/(\rho^3-R)\), with \(\sigma\rho=\zeta^2\rho\).
Let \(P=(X,Y)\) be the chosen local point, and choose \(c\in K^\times\)
such that
\[
 \alpha=Y-3Rs=a_i c^3,\qquad u=3Rsd.
\]
Avoid the zeros of \(\alpha\) and \(Y^2+162R^2\). Define
\[
 \begin{split}
 g&=\frac{\alpha+s\zeta\rho X}{X+6\zeta\rho^2},&
 \xi_P&=\sigma^2g=
          \frac{\alpha+s\zeta^2\rho X}{X+6\rho^2},\\
 \theta&=c\xi_i/\xi_P,&
 \gamma_h&=h+\theta\sigma h+
                    \theta\sigma\theta\,\sigma^2h .
 \end{split}                                                   \tag{23}
\]
For at least one \(h\in\{1,\rho,\rho^2\}\), the element
\(\gamma_h\) is a unit of \(D\). For any such choice,
\[
 C_{i,w}=\frac{a_i}{(Y-u)N_{D/K}(\gamma_h)}\in K^\times,
 \qquad
 \frac{b_i}{t_T(P)}=C_{i,w}(\sigma\gamma_h)^3.             \tag{24}
\]
Consequently \([C_{i,w}]\) supplies the class in (19), for tame
and wild coordinates alike.

Here is an algebraic proof, including the normalizing factor in (24).
The secant from \(S\) to \(T\), divided by the vertical line at
\(S+T\), gives \(g\). The curve equation and \(s^2=6,d^2=-3\)
give the exact identities
\[
 N(g)=\alpha,\qquad
 \frac{\sigma t_T(P)}{t_T(P)}=\frac{\alpha}{g^3},\qquad
 g^2\sigma g=\frac{\alpha}{Y-u}\,t_T(P).                 \tag{25}
\]
For example the denominator norm in \(g\) is
\(X^3+216R^2=Y^2+162R^2\). These identities hold in the
coordinate ring after clearing the displayed denominators, so they
hold on every admitted local point. In particular \(N(\xi_P)=\alpha\)
and \(N(\theta)=1\).

The explicit Hilbert–90 sum in (23) satisfies
\(\gamma_h=\theta\sigma\gamma_h\). After extending scalars to a
splitting field, \(\sigma\) cyclically permutes the three components
of \(D\). Since \(\theta\) is a unit, this identity implies that a
nonzero \(\gamma_h\) is a unit. The linear operator
\(h\mapsto\gamma_h\) is nonzero: in the split algebra its matrix has
diagonal entries one. Thus one of the three basis choices is nonzero.
Writing \(\gamma=\gamma_h\), one obtains
\[
 \xi_i=\frac{\xi_P}{c}\frac{\gamma}{\sigma\gamma},
 \qquad
 b_i=\frac{\sigma(\xi_P)^2\sigma^2(\xi_P)}{c^3}
               \frac{(\sigma\gamma)^3}{N\gamma}.
\]
Since \(\sigma\xi_P=g\), the last identity in (25) and
\(\alpha=a_i c^3\) prove (24). In particular the cube factor in
(24) is \((\sigma\gamma)^3\); the matching root \(c\) has already
cancelled. Changing the point, root or Hilbert–90 seed changes no
pairing value, by the same local exactness used in MCT7.

### Computing the wild local coordinates

**Lemma MCT9 (finite cube test at three).** In
\(K=\mathbf Q_3(\zeta)\), let \(\lambda=1-\zeta\) and normalize
\(v(\lambda)=1\). A unit is a cube if and only if it is congruent
to \(1\) or \(-1\) modulo \(\lambda^4\). A nonzero element is a
cube if and only if its valuation is divisible by three and its unit
part has this property, using any fixed uniformizer.

For the necessity, use \(\pi=\sqrt{-3}\), which generates the same
prime ideal as \(\lambda\). Every integral unit is \(a+b\pi\),
with \(a\in\mathbf Z_3^\times,b\in\mathbf Z_3\). Its cube is
\(a^3+(b^3-a^2b)\pi^3\) modulo \(\pi^4\). The second coefficient
is divisible by three, while \(a^3\equiv\pm1\pmod9\).
Conversely, the logarithm and exponential identify
\(1+\lambda^2\mathcal O_K\) with \(\lambda^2\mathcal O_K\),
and cubing with multiplication by 3. Thus every element of
\(1+\lambda^4\mathcal O_K\) is a cube; both signs are cubes.

This gives the coordinates in (21) without a numerical root decision:
remove the valuation modulo three, enumerate the 27 products of the
three unit basis elements, and apply the unit test. Exactly one
product works, by the basis statement preceding (21). Relative
precision strictly greater than four suffices to distinguish the
classes. The checker below requires a larger precision margin and
rejects unresolved valuations. Local matching cube roots are obtained
by the same logarithm argument at 3 and by unit Hensel lifting away
from 3.

### Direct local certificate replay

The separate direct checker uses (23)–(24), without reconstructing
an alternating matrix from rational radical vectors:

```sh
sage -python research/arithmetic/verify_mordell_ray_pairing_local.py
```

Its input is [the direct local certificate](mordell_ray_pairing_local_exact.json).
The original retained cases are \(R=145\), in basis \((0,5,29)\), and
\(R=115\), in basis \((5,23)\). They give, respectively,
\[
 \begin{pmatrix}0&2&1\\1&0&1\\2&2&0\end{pmatrix},
 \qquad 0_{2\times2}.                                      \tag{26}
\]
They test both a split local cubic algebra and a cubic field at 3.
The input also retains \(R=71\), in basis \((0,71)\), with matrix
\(\left(\begin{smallmatrix}0&1\\2&0\end{smallmatrix}\right)\).
It cross-checks the new full evaluator below by the explicit scalar
Hilbert–90 descent, using the same exact norm witnesses but different
local arithmetic. Its eight records cover \(3,11,71\).
The input contains exact global norm coefficients and local point
seeds; the checker derives the matching cube roots and scalar classes.
It verifies the norm identities over \(\mathbf Q(\sqrt6)\) before
performing local arithmetic. The three identities in (25) are also
checked by exact polynomial reduction for symbolic \(X,Y,\rho\),
subject to the curve and quadratic relations. This is not interpolation
from a finite set of curves.

The global norm solutions for 145 were found using PARI `rnfisnorm`
with an auxiliary prime search and then checked exactly; the solver's
completeness is not used. The norm solutions for 115 can be recovered
from its retained rational point identities using \(\xi_P\) and the
quadratic norm down to \(L(\rho)\). Those points supply inputs to the
global lift, rather than any entries or ranks of the resulting matrix.

For each row the checker covers both places at every split prime of
\(L\), and the unique place at a nonsplit prime, over all primes
dividing \(3R\) or a norm-coefficient denominator. At 2 all columns
are locally zero by MS1, so this place contributes zero. The real
places have no odd-order Brauer invariant. For 145 the retained
rational primes are \(3,5,17,19,29,43\), with twenty local records.
For 115 they are \(3,5,23,59,114259\), with thirteen local records.
The auxiliary primes are part of the proof support. For example, the
row for \(\chi_0\) at 145 has local symbol vector \((1,2,2)\) at 3
and \((2,2,0)\) at 17, each with weight 2 in (20); together they give
the first row \((0,2,1)\). Omitting 17 would even give a nonzero
diagonal entry.

Internal replay uses SageMath 10.9 with 240 digits of rational-prime
precision, exact global coefficients and capped-relative local
arithmetic. Scalar and matching-root checks require more than twelve
uniformizer digits of relative accuracy, and the unit-class decision
requires strictly more than four. Unresolved valuations or a precision
failure abort the replay. Corrupting an exact norm or removing an
auxiliary place is rejected. The local comparison, reciprocity and
the interpretation as a Cassels–Tate pairing remain the written inputs
in MCT7–MCT9; this is internal computational replay, not formal
verification or independent mathematical review.

The direct checker accepts other admissible parameters with exact norm
coefficients and local point seeds in the same format. Its `--generate`
mode fills local records from supplied global norms by a bounded
point search. Exhausting that search proves no insolubility. A
uniformly validated certificate generator from \(R\) alone, including
the global norm search, is not supplied by the two retained examples.

## Exact controls beyond the original congruence class

The retained input now also gives the following complete normalized
matrices. A basis index 0 denotes \(\chi_0\).
\[
 \begin{array}{c|c|c}
 R&\text{basis}&B_R\\\hline
 19895&(5,23,173)&0_{3\times3}\\
 235&(0,5,47)&0_{3\times3}\\
 3335&(5,23,29)&\begin{pmatrix}0&1&1\\2&0&2\\2&1&0\end{pmatrix}\\
 145&(0,5,29)&\begin{pmatrix}0&2&1\\1&0&1\\2&2&0\end{pmatrix}\\
 27485&(5,23,239,0)&
 \begin{pmatrix}0&1&2&0\\2&0&1&2\\1&2&0&1\\0&1&2&0\end{pmatrix}\\
 33695&(5,23,293,0)&
 \begin{pmatrix}0&2&1&1\\1&0&2&1\\2&1&0&1\\2&2&2&0\end{pmatrix}
 \end{array}                                             \tag{22}
\]
The first two are rational controls of full Selmer dimension three.
The other four have pairing rank two; their rational Kummer dimensions
are respectively 1,1,2,2, giving those exact rational ranks and
\(\dim\Sha(E_R)[3]=2\). In the last two matrices both (14) and
(17) are checked. For \(R=33695\) both corresponding entries are
nonzero, so this checks their relative normalization rather than
only agreement on a zero entry.

For example \(a_{29}=145+87s\), \(a_{47}=-611+282s\),
\(a_{173}=24047+9861s\), \(a_{239}=-5975+2868s\), and
\(a_{293}=-23147+9669s\). All twenty rational Kummer identities,
including those that determine the two reference zero matrices,
are given as exact point/vector/cube-root certificates in the input.
The calibration uses the full rational span at 19895. The other new
matrices follow from the formulas and their own rational radical
vectors; no analytic rank assumption is used.

## Exact full-matrix evaluation from norm certificates

The scalar formula (24) supplies a general local descent. For this
particular family one can also evaluate the full matrix using exact
number-field local cube tests, with only nine candidate classes at 3.
This second route supplies an automatic, bounded norm-search front end
and avoids fitting any entry from a rational radical vector.

**Theorem MCT10 (full-matrix certificate rule).** Let \(R>1\) be
admissible and let \(a_i\) be the complete basis of MCT3. Given exact
norm witnesses \(\xi_i\) in (18), the finite rules (27)–(30) below
compute all off-diagonal entries, including every wild entry, up to
one common global Weil normalization. Diagonals are zero by
alternation. All candidate-class searches are exhaustive, and all
required prime contributions are included. The implementation checks
skew symmetry after computing both directions separately.

### A global lift that descends to the rational curve

Let \(F=L(\zeta)\), \(M=F(\rho)\), and let \(\tau\) fix
\(\rho,\zeta\) and send \(s\) to \(-s\). Starting with (18), put
\[
                 \beta_i=\left(\frac{b_i}{\tau b_i}\right)^2.
                                                               \tag{27}
\]
This averaging step permits using just one place of \(F\) over each
rational prime. It cannot be omitted for an arbitrary \(L\)-lift.
Indeed \(\tau S=-S\), \(\tau T=-T\), so the conjugate cohomology
class has coordinates \(((\tau a_i)^{-1},(\tau b_i)^{-1})\).
Twice the sum of the original class and its conjugate has second
coordinate (27) and first coordinate \((a_i/\tau a_i)^2\).
Since \(a_i\tau a_i=c_i^3\) for rational \(c_i\), this first
coordinate is \(a_i\) modulo cubes. The average is invariant under
\(\operatorname{Gal}(L/\mathbf Q)\). Inflation–restriction, whose
order-two cohomology vanishes on these three-torsion modules, therefore
descends it to a global full-three-torsion lift over \(\mathbf Q\).
Local invariants can now be computed at any chosen place of \(F\)
and divided by its local degree over \(\mathbf Q_p\), as in
van Beek–Fisher Lemma 1.3. This explains why the two conjugate places
of \(L\) no longer have to be summed separately.

### The coordinate at 3 reduces to a three- or nine-class test

Choose the place of \(F\) above 3 where
\(sd/3=\sqrt{-2}\equiv4\pmod9\). Its completion is
\(K_3=\mathbf Q_3(\zeta)\), of degree two over \(\mathbf Q_3\).
With \(\lambda=1-\zeta\) and \(\eta_3=1-\lambda^3\), there are
unique \((\alpha_j,\upsilon_j)\in\mathbf F_3^2\) with
\[
                     [a_j]=[\zeta]^{\alpha_j}[\eta_3]^{\upsilon_j}.
\]
The norm-one subspace is the anti-invariant Kummer subspace. Its two
basis classes are \(\zeta,\eta_3\); the implementation also checks
existence and uniqueness by all nine local cube tests. For a tame
class \(\alpha_j=0\).

Choose a local rational point \(P_i\) with forward class \(a_i\)
and form \(D_i=\beta_i/t_T(P_i)\). If the class is locally zero,
take \(P_i=O\) and the tangent class is 1. Since the rational kernel
of the reverse isogeny over \(\mathbf Q_3\) is of cyclotomic type,
the descended difference class comes from
\(\mathbf Q_3^\times/\mathbf Q_3^{\times3}\), with basis 3 and 4.

If \(\epsilon_R=1\), there is a cube root of \(R\) in
\(\mathbf Q_3\). Hence \(M\) splits locally over \(K_3\), and
there are unique \(e_i,f_i\in\mathbf F_3\) such that
\(D_i/(3^{e_i}4^{f_i})\) is a cube in \(K_3\). All nine tests
are made. Thus the cubic descent that appears in (19) requires no
extension-field root search in precisely the case with a wild
coordinate.

If \(\epsilon_R=0\), the extension \(K_3(\rho)/K_3\) is cubic.
The restriction kernel on \(K_3^\times/K_3^{\times3}\) is generated
by \([R]\), which spans the rational unit direction \([4]\).
There is therefore a unique \(e_i\) such that \(D_i/3^{e_i}\) is
a cube in that cubic completion. Only three tests are needed. Put
\(f_i=0\); all columns are tame, so the omitted unit direction
annihilates them. In both cases the contribution over \(\mathbf Q_3\) is
\[
       \ell_{ij}(3)=2\bigl(e_i\upsilon_j+2f_i\alpha_j\bigr)
                      \quad(i\ne j).                      \tag{28}
\]
To see this directly from (21), use
\([3]=[\lambda]^2[\zeta]^2\) and
\([4]=[1-\lambda^2][\eta_3]\). The Hilbert exponent is
\(e_i\upsilon_j+2f_i\alpha_j\), and the local restriction degree
is two. The checker verifies these two basis identities by local
cube tests as well.

Local point seeds do not require a global rational point search.
One unramified generator is
\(P_u=(1,\sqrt{1+54R^2})\), choosing the root congruent to 1
modulo 3. If \(\epsilon_R=1\), a complementary generator is
\(P_w=(3,9\sqrt{(1+2R^2)/3})\). The radicand in the second root
is \(1\pmod3\). After removing a cube, its Kummer value has a
nonzero linear coefficient modulo 3, so it is ramified; the first
point gives a nonzero unramified class. Thus
\(mP_u+nP_w\), \(m,n\in\{0,1,-1\}\), supplies all nine local
classes. Without the wild coordinate use \(O,\pm P_u\).

The implementation retains rational approximations \((X,Y_0)\).
For \(h=v_3(Y_0^2-X^3-54R^2)-v_3(Y_0)\), it requires the strict
Hensel hypothesis \(v_3(Y_0^2-X^3-54R^2)>2v_3(Y_0)\).
It also verifies that the matching Kummer value and every tangent or
norm-factor evaluation is stable modulo a cube. A relative error of
valuation strictly greater than \(2v(3)\) is sufficient, by Hensel
applied to \(Z^3-u\) at 1. This is an exact integer-valuation bound,
not agreement between two numerical precision settings. A failed
bound in verification is rejected. Generation doubles precision and
uses the finite point enlargement proved in MCT13 below.

### Tame and auxiliary contributions

At an allowed prime \(q\mid R\), an off-diagonal second argument
is locally zero unless it is \(\chi_q\). If \(j\) indexes
\(\chi_q\) and \(i\ne j\), the first argument is locally zero,
so use \(P_i=O\). Let \(\mathfrak P\) be a prime of \(M\) above
a chosen prime \(\mathfrak q\) of \(F\) over \(q\), and normalize
\(v_{\mathfrak P}(\rho)=1\). The value
\(v=v_{\mathfrak P}(\beta_i)\) is divisible by three. The removed
factor \(\rho^v\) is therefore a cube. Its residue-field exponent gives
\[
 \ell_{ij}(q)=2v_{\mathfrak q}(a_j)
   \operatorname{Ind}_{\zeta}
       \overline{(\beta_i/\rho^v)}^{(q^2-1)/3}.           \tag{29}
\]
The residue field is \(\mathbf F_{q^2}\). The descended local
kernel class is unramified, since an order-three character of
\(\mathbf Q_q^\times\) cannot ramify for \(q\equiv2\pmod3\).
This proves that its unit residue gives the needed symbol. A wild
second argument is locally zero at every such \(q\).

For an auxiliary rational prime \(q\nmid6R\), the curve has good
reduction and all local point lifts are unramified. Their contribution
to the ramified coordinate of the discrepancy is zero. The extension
\(M/F\) is also unramified there. Let
\(\nu_i=v_{\mathfrak P}(\beta_i)\bmod3\),
\(f_q=[F_{\mathfrak q}:\mathbf Q_q]\in\{1,2\}\), and put
\[
 \gamma_j=\operatorname{Ind}_{\zeta}
            \overline{a_j}^{(q^{f_q}-1)/3},\qquad
       \ell_{ij}(q)=-f_q^{-1}\nu_i\gamma_j.               \tag{30}
\]
The residue power in (30) must be taken in \(F\), not in a degree-three
unramified extension, which would erase the unramified cubic class.
Formula (30) is the tame Hilbert-symbol formula with the two factors
interchanged. Its sign is essential.

The auxiliary set consists of all prime divisors of all denominators
of the \(\xi_i\) coefficients, with \(2,3\) and the primes of
\(R\) handled separately. The certificate retains and verifies a
complete prime factorization of each denominator. Outside this set,
\(N\xi_i=a_i\) makes the norm factors units and the invariant is
zero. The place 2 contributes zero by MS1, and the real place has
no odd-primary invariant. Summing (28)–(30) proves MCT10.

### Implementation and retained complete matrices

[The full evaluator](mordell_ray_full_pairing.py) uses exact algebraic
numbers and PARI's exact local power test, documented in
[the retained official function catalogue](sources/pari-general-number-fields-local-powers.html)
under `nfislocalpower` and `nfmodpr`.
[The certificates](mordell_ray_full_pairing_certificates.json) retain
the norm coefficients, denominator factorizations, rational Hensel
seeds, the cube-root approximation when needed, all local classes at
3 and every local contribution. Verification never invokes a norm
solver. The software pins are SageMath 10.9 and PARI 2.17.3.

```sh
sage -python research/arithmetic/mordell_ray_full_pairing.py \
  --verify research/arithmetic/mordell_ray_full_pairing_certificates.json \
  --compare-controls --self-test
```

The seven retained cases are
\[
 \begin{array}{c|c|c}
 R&\text{basis}&\text{pairing rank}\\\hline
 53&(0,53)&0\\71&(0,71)&2\\115&(5,23)&0\\
 145&(0,5,29)&2\\235&(0,5,47)&0\\
 6095&(5,23,53)&2\\355&(5,71)&2
 \end{array}
\]
The four overlapping complete matrices agree with the earlier controls.
At 71, the two off-diagonal entries receive contributions \((2,2)\)
from 3, \((0,0)\) from 71, and \((2,0)\) from the auxiliary prime
11. Their sum is \((1,2)\). Thus even this smallest new wild example
requires the auxiliary contribution; the local matrix at 3 alone is
not the answer.
The scalar-descent checker also gives the same new matrix at 71,
without the averaging or number-field local-cube route. This is a
cross-check of implementations, not independent mathematical review.
The full evaluator rejects a corrupted norm, an invalid Hensel point,
an omitted auxiliary-prime record, a corrupted output matrix and a
corrupted prime-block witness.

Automatic input generation is available, for example:

```sh
sage -python research/arithmetic/mordell_ray_full_pairing.py \
  --generate 71 --output .scratch/mordell-pairing/full-71.json
```

Each norm-search attempt runs in an isolated process with a timeout.
The generator accepts only exact identities, so no class-group
completeness or GRH assumption used by a search routine enters the
retained theorem. A failed search, precision bound, factorization or
resource limit produces no matrix claim. In exploration, a norm search
at 27485 exhausted its PARI stack and a 45-second attempt at 7685
did not complete; neither is part of the retained full-certificate
coverage. MCT12–MCT13 below supply a prime-block reduction and a terminating
unbounded fallback. They supply no useful uniform runtime bound.

**Corollary MCT11 (two rank-zero curves and explicit Hasse failures).**
The full certificates give
\[
 B_{71}=\begin{pmatrix}0&1\\2&0\end{pmatrix}_{(0,71)},\qquad
 B_{355}=\begin{pmatrix}0&2\\1&0\end{pmatrix}_{(5,71)}.
\]
Consequently, for both parameters,
\[
       \operatorname{rank}E_R(\mathbf Q)=0,\qquad
       \Sha(E_R)[3]\cong(\mathbf Z/3\mathbf Z)^2.
\]
Both ray spaces have dimension two by MS1. Their nondegenerate pairing
leaves no rational Kummer direction, proving the rank assertion; the
Kummer exact sequence gives the stated group. No global point search
or finiteness conjecture for \(\Sha\) is used.
In particular the following smooth plane cubics have points over every
completion of \(\mathbf Q\) and have no rational point:
\[
 \begin{split}
 2U^3+15U^2V+36UV^2+30V^3+213W^3&=0,\\
 U^3-7U^2V+18UV^2-14V^3+71W^3&=0.
 \end{split}                                             \tag{31}
\]
The first is the class \(a_0=5+2s\) at 71; the second is the class
\(a_5=-35+15s\) at 355. They come from the coefficient of \(s\)
in \(a_i(U+Vs)^3=(Y-3Rs)W^3\). Their binary discriminants are
\(-648\) and \(-200\), respectively. Local solubility is Selmer
membership; the nonzero pairing proves global insolubility. As in
(8a), a rational point at infinity would make the nonzero Kummer
class a cube and is also excluded.

## Prime-block norm assembly and a terminating general algorithm

The closest norm reduction in the primary literature is van Beek–Fisher
Theorem 3.1, which applies over every number field. Their Algorithm 3.6
and its iteration bound, Theorem 3.7, are stated over the rationals.
The paragraph following Example 3.10 explicitly leaves termination of
analogues over other number fields unproved. We use the norm-swapping
identity over \(L=\mathbf Q(\sqrt6)\), and supply a different, elementary
termination argument. We claim no extension of their iteration bound.

**Theorem MCT12 (prime-block assembly of every norm row).** Fix an
allowed prime \(p\), the tame representative \(a_p\) of MCT3, and
\(F_p=L[r_p]/(r_p^3-a_p)\). Every allowed prime \(q\), including
\(q=p\), belongs to \(N_{F_p/L}(F_p^\times)\). Consequently the
reverse witness for any admissible \(R\) can be assembled as
\[
 \eta_{p,R}=\prod_{q\mid R}\eta_{p,q},\qquad
 N_{F_p/L}(\eta_{p,q})=q.                                  \tag{32}
\]
The field and the individual prime witnesses depend on \(p,q\),
independently of the other factors of \(R\).

For the wild row put \(F_0=L[r_0]/(r_0^3-a_0)\). For a squarefree
product \(b>1\) of allowed primes,
\[
 b\in N_{F_0/L}(F_0^\times)
       \quad\Longleftrightarrow\quad b\equiv\pm1\pmod9.    \tag{33}
\]
Whenever the wild coordinate occurs, the factors of \(R\) can be
partitioned into blocks \(b\) of at most three primes satisfying
(33). Multiplying their witnesses constructs \(\eta_{0,R}\).

For any row, write \(\eta_{i,R}=A+B r_i+C r_i^2\). The exact
forward witness is
\[
 \xi_{i,R}=
 \frac{AB-C^2a_i-B\rho}{-(B^2-AC)-C\rho},\qquad \rho^3=R,
 \qquad N_{L(\rho)/L}(\xi_{i,R})=a_i.                     \tag{34}
\]
Thus the expensive reverse number field depends on the row rather
than the full parameter. This is an algebraic reduction; it does not
guarantee smaller coefficients or easily factored denominators.

Proof. MS1 and MCT7 give a forward norm witness for the tame class
\(a_p\) at parameter \(p\). The equivalence (i)–(ii) in van
Beek–Fisher Theorem 3.1 therefore makes \(p\) a norm from \(F_p\).
For \(q\ne p\), apply the same argument at the squarefree admissible
parameter \(pq\); the same representative \(a_p\) is a Selmer
class there. Hence \(pq\) is a norm. Dividing the two witnesses makes
\(q\) a norm. Multiplicativity proves (32). No assertion about the
rational rank of either curve is used.

If \(b\equiv\pm1\pmod9\), MS1 includes \(a_0\) at parameter
\(b\); MCT7 and norm swapping give sufficiency in (33). Conversely,
at 3 the completion of \(L\) is \(\mathbf Q_3(\zeta)\).
The class \(a_0\) has a nonzero \([\zeta]\) coordinate, while
\([\eta_3]\) is unramified. A rational 3-adic unit is a power of
\([4]\) modulo cubes. The table (21) gives a nonzero symbol between
\([\zeta]\) and \([4]\), and zero between \([\eta_3]\) and
\([4]\). The cyclic local norm criterion therefore forces \(b\)
to be a rational cube at 3, equivalently \(b\equiv\pm1\pmod9\).
This proves necessity without a global integral-norm assumption.

Let \(\chi_9\) be the cubic character used in (14). Keep each prime
with \(\chi_9(q)=0\) as a singleton, pair as many primes with
characters 1 and 2 as possible, and group the remaining primes into
triples of the same character. Their count is divisible by three
because \(\chi_9(R)=0\). This proves the block bound. Finally,
(34) is the explicit norm swap (22) in van Beek–Fisher, with its
variables interchanged. Substitution, or the cubic norm identity
below, verifies the signs. Its denominator cannot vanish: \(R\) is
not a cube in \(L\), and \(B=C=0\) would imply that \(R=A^3\).
This completes the proof.

**Theorem MCT13 (unrestricted effective full-matrix rule).** For every
admissible \(R>1\), there is a terminating exact-arithmetic algorithm
for its entire ray-basis pairing matrix, including the coordinate at
3 when present, up to the one common Weil normalization in MCT10.
No rank oracle, rational-point oracle, GRH assumption or conjecture
about finiteness of \(\Sha\) is required. No useful uniform runtime
or coefficient-height bound is asserted. The case \(R=1\) is the zero
one-dimensional alternating matrix.

The algorithm first tries bounded norm searches, using MCT12 if
convenient. If these fail, it enumerates positive denominators \(d\)
and six integer coefficients of \(A,B,C\in\mathbf Z[s]\), by
increasing maximum absolute height, until
\[
 A^3+R B^3+R^2C^3-3RABC=a_i d^3.                          \tag{35}
\]
It then takes \(\xi_i=(A+B\rho+C\rho^2)/d\). The norm identity
(35) is checked exactly. Every witness assured by MCT7 has such a
representation after clearing rational denominators; dividing their
common integer gcd gives a primitive tuple in the enumeration.
Each height has finitely many tuples, so the search terminates.
Factoring the finitely many resulting denominators is an effective
finite step; no bound on its cost is needed for this assertion.

For local computation, double the rational 3-adic precision until the
strict Hensel and cube-stability inequalities in MCT10 hold. There is
a small point-selection issue: a fixed representative can lie at a
zero of a tangent function. The following uniform finite enlargement
of the point list removes that issue. Put
\[
 Q=\left(\frac19,\frac1{27}\sqrt{1+54R^2 3^6}\right),
 \qquad \sqrt{1+54R^2 3^6}\equiv1\pmod3.                 \tag{36}
\]
For each of the three or nine seed points in MCT10, try
\(P+3hQ\), \(0\le h\le5\), as precision increases.

Indeed \(Q\) lies in the formal group, with parameter
\(-x/y\) of valuation one. The formal logarithm has that same
valuation (its terms of degree \(n\ge2\) have valuation at least
\(n-v_3(n)>1\)), so \(Q\) is not torsion. All six translates
are distinct, and \(3Q\) has zero isogeny Kummer class because
\([3]\) factors through the isogeny. They therefore preserve the
required local class. Among affine points, Hensel lifting of \(Y\)
only excludes the at most three points with \(Y=0\); the tangent
function excludes only \(T\), since its divisor is \(3T-3O\).
Including \(O\) leaves at most five forbidden points. Also
\(Y-3Rs\ne0\) for a rational affine point, since
\(s\notin\mathbf Q_3\). At least one of the six translates for
each class is consequently usable.

The rational approximations obtained by Hensel roots and addition
converge to those exact points; removable intermediate failures are
skipped. At a usable point the Kummer and tangent values are nonzero,
so their relative errors tend to infinite valuation. Each of the four
norm factors at a local cube root of \(R\) is nonzero, by its
nonzero exact norm. Their error margins also tend to infinity.
Eventually all the strict inequalities exceed \(2v(3)\).
The finite local cube tests then certify the selected classes, rather
than infer them from convergence. For a supplied certificate, a bad
bound is rejected instead of silently increasing its recorded precision.
Finally MCT10 supplies all remaining finite local contributions and
proves that their sum is the required pairing. This proves termination
and correctness for arbitrary finite prime support.

The implementation offers the exhaustive fallback with `--complete`;
its ordinary mode keeps the global searches bounded and reports
incompletion when none succeeds. Local precision is adaptive in both
generation modes. The exhaustive fallback is deliberately simple and
can be prohibitively slow. It closes the existence/termination issue,
while practical scalability and a formula using only Frobenius symbols
remain different obligations. MCT12 also does not identify any of the
individual pairing entries with an additive function of \(R\): it
assembles norm witnesses, and all subsequent local terms still apply.

### Optimization evidence and practical limits

The retained 145 certificate now uses (32)–(34) for all three rows:
the wild row has the block \(5\cdot29\), and each tame row has the
two singleton prime targets. Its five block witnesses are verified
individually and multiplied before checking (34). Starting generation
at precision 6 triggers an increase to 12, where all strict bounds
pass. The resulting matrix is unchanged from (26) and the separate
scalar-descent calculation. A deliberately corrupted block witness
is rejected by full replay.

For 355, using an auxiliary-prime cutoff 11 in the direct PARI search
instead of the previous cutoff 30 reduces the complete retained case
from about 139 KB to 31 KB in compact JSON. In this internal run,
generation and full local evaluation took about two seconds. The
matrix remains the nondegenerate matrix of MCT11. This is a measured
improvement on that parameter, not a monotonicity theorem about the
cutoff: different cutoffs can give much larger witnesses.

The bounded searches also put denominator factorization inside each
isolated worker, so a difficult denominator is covered by its timeout.
Prime-block generation at 27485 still exceeded a 25-second worker
budget in row 239; the reverse-only alternative exceeded 30 seconds
there. At 7685 the prime-block, direct and reverse alternatives each
exceeded a 15-second worker budget in row 29. Neither run produced a
complete matrix, and neither is claimed as a completed certificate.
They do not contradict MCT13: the exhaustive fallback was not invoked.
The retained seven cases cover prime-block assembly, adaptive precision,
split and nonsplit cubic local algebras at 3, nonzero auxiliary terms,
and both zero and nonzero full matrices. They do not constitute a
performance bound for arbitrary parameters.

## Evidence, replay and remaining scope

Run from the repository root:

```sh
python3 research/arithmetic/verify_mordell_ray_pairing.py
```

The checker uses Python's standard library and exact fractions. It
checks (2) and (16), all four embedding choices in each retained symbol,
the central twist (15), all twenty rational Kummer identities,
independence of the rational vectors, the correction calibrated by
the zero reference matrix, and uniqueness and ranks of all ten complete
matrices. It also rejects a
corrupted norm value. It does not certify the global comparison,
the descent of the governing field, Chebotarev or MS1. Those are
written proof inputs, with their exact dependencies above.

Internal execution: Python 3.12.3, all checks passed. Exploration used
the repository's SageMath 10.9 / PARI 2.17.3 environment to find a norm
solution and rational points; neither solver is needed to replay the
retained certificate. The source/input digests are pinned in the
mathematical catalogue. Retained primary versions are
[Koymans–Smith v1](sources/koymans-smith-cubic-twist-pairing-2405.09311v1.pdf),
[van Beek–Fisher v1](sources/van-beek-fisher-cubic-norm-pairing-1711.02432v1.pdf),
and [Wu–Zhu v1](sources/wu-zhu-fermat-pairing-2611.29541v1.pdf).

The full-matrix evaluator now supplies exact local matching and descent
from verified norm inputs, with a bounded automatic search front end.
MCT13 gives an unrestricted terminating algorithm with an exhaustive
fallback. Remaining work concerns practical norm generation and scaling
to larger supports, plus a Frobenius-only formula with all pairwise
reference constants evaluated. The usual
\(j/3\mapsto j\) Weil-sign calibration, independent replay and
external mathematical review remain separate tasks.
