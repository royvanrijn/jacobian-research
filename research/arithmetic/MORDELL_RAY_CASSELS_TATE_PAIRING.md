# A nonzero cubic ray pairing block on infinitely many inert-prime twists

This note uses the family, isogeny conventions and exact Selmer calculation
of [MS1](MORDELL_INERT_PRIME_SELMER_MOMENTS.md). Put
\(E_R:y^2=x^3-2R^2\), \(E'_R:y^2=x^3+54R^2\), and
\(S_R=\operatorname{Sel}_\psi(E_R/\mathbf Q)\).
The pairing considered here is the Cassels–Tate pairing on
\(\Sha(E_R)[\psi]\), pulled back to \(S_R\).
It is also the restriction of the full three-Selmer pairing under MS1.

**Result and scope.** There is an explicit residue-field formula for the
entry between the ray characters at 5 and 23. It is nonzero for two thirds
of the allowed moving primes \(q\equiv-1\pmod9\) in \(R=115q\).
Thus infinitely many admissible curves have rational rank at most one and
\(\dim_{\mathbf F_3}\Sha(E_R)[3]\ge2\).
Complete normalized matrices are obtained below for selected products of
two, three and four allowed primes. A formula for **every** entry for
arbitrary \(R\), including the extra coordinate at 3, remains open here.
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

## Evidence, replay and remaining scope

Run from the repository root:

```sh
python3 research/arithmetic/verify_mordell_ray_pairing.py
```

The checker uses Python's standard library and exact fractions. It
checks (2), all four embedding choices in each retained symbol (4),
the eight identities (10), independence of the rational vectors, and
uniqueness and ranks of all four matrices (9). It also rejects a
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

The next extension is the arbitrary ray-basis matrix: construct the
corresponding exact norm/cochain data for every pair, retain the local
term at 3 when the twist ratio is not a local cube, and include the
extra wild coordinate when \(\epsilon_R=1\). The fixed 5–23 formula
and infinite-family theorem do not discharge that larger obligation.
The usual \(j/3\mapsto j\) Weil-sign calibration, independent replay
and external mathematical review also remain separate tasks.
