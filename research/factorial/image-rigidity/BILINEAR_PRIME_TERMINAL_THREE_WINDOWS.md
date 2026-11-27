# Prime and doubled-prime pure-cross three-windows

Sections 1--7 exclude every noncentral resonance of odd prime terminal order. Combining this with the canonical opposite-axis theorem closes the pure-cross three-window at every odd prime terminal order. Sections 8--12 prove doubled-prime windows away from the two fourth-root axis ratios. Section 13 gives a uniform good-prime lift from a smaller Chebyshev array, including further explicit families at terminals \(3p,4p,5p\). The exceptional axis families and arbitrary terminal orders remain open.

## 1. Statement and normalization

Define monic integer polynomials
\[
D_N(u)=2T_N(u/2),\qquad E_n(u)=U_n(u/2),
\]
and put
\[
F_n(u,A)=\sum_{j=0}^n E_n^{(j)}(u)A^j.
\tag{1}
\]
Thus
\[
F_n(u,A)=G_n(u/2,A/2),\qquad
G_n(c,b)=\sum_{j=0}^n U_n^{(j)}(c)b^j.
\tag{2}
\]
The leading coefficient of \(F_n\) in \(A\) is \(n!\), independent of \(u\).

**Theorem.** Let \(p\) be an odd prime and \(n=p-1\). If \(c\ne0\) and \(T_p(c)=0\), then \(G_n(c,b)\) is squarefree as a polynomial in \(b\) over \(\mathbb C\). Equivalently, \(F_n(u,A)\) is squarefree at every nonzero zero \(u=2c\) of \(D_p\).

The theorem is algebraic and therefore remains true after extension to any characteristic-zero field. It does not assert irreducibility of \(G_n\).

## 2. The noncentral resonance polynomial is Eisenstein

The integer recurrences
\[
D_0=2,\ D_1=u,\ D_{N+1}=uD_N-D_{N-1},
\qquad
E_0=1,\ E_1=u,\ E_{n+1}=uE_n-E_{n-1}
\]
give monicity, parity, and \(D_N'=NE_{N-1}\). The Laurent identity
\[
D_N(z+z^{-1})=z^N+z^{-N}
\]
follows from the same recurrence. In characteristic \(p\),
\[
D_p(z+z^{-1})=(z+z^{-1})^p.
\]
Substitution \(u\mapsto z+z^{-1}\) is injective on \(\mathbb F_p[u]\), as seen from the highest Laurent exponent. Hence
\[
D_p(u)\equiv u^p\pmod p.
\tag{3}
\]
Since \(p\) is odd, \(D_p\) is odd, and
\[
M_p(u):=D_p(u)/u
\]
is a monic integer polynomial of degree \(p-1=n\). All nonleading coefficients are divisible by \(p\), while
\[
M_p(0)=D_p'(0)=pE_{p-1}(0)
=p(-1)^{(p-1)/2}
\tag{4}
\]
is not divisible by \(p^2\). Eisenstein's criterion proves that \(M_p\) is irreducible over \(\mathbb Q\).

Consequently every nonzero zero of \(D_p\) has minimal polynomial \(M_p\), and
\[
\overline{M_p}(u)=u^n.
\tag{5}
\]

## 3. A formal discriminant-order lemma

Let \(k\) be an algebraically closed field of characteristic different from two. Suppose \(P(u,A)\in k[u,A]\) has degree \(n\) in \(A\) and a leading coefficient which is a nonzero constant. Assume that at \(u=0\):

* every multiple root \(a\) of \(P(0,A)\) is nonzero;
* \(P_{AA}(0,a)\ne0\);
* \(P_u(0,a)\ne0\).

Then every multiple root is exactly double. If there are \(r\) such distinct roots, then
\[
\operatorname{ord}_{u=0}\operatorname{disc}_A P(u,A)=r\le\lfloor n/2\rfloor.
\tag{6}
\]

Here is the complete local argument. After dividing by the constant leading coefficient, factor \(P(0,A)\) into its monic factors for distinct roots. Those factors are pairwise coprime. The elementary coefficient-lifting argument using their Bézout identities uniquely factors \(P\) over \(k[[u]][A]\) into monic factors \(Q_a(u,A)\) of the same degrees: at each power of \(u\), the coefficient error is uniquely distributed among the factors by reduction modulo their coprime constant factors.

A simple-root factor has unit discriminant. At a double root, translate \(A=a+y\) and write
\[
Q_a=y^2+\alpha(u)y+\beta(u),\qquad
\alpha(0)=\beta(0)=0.
\]
The other factors and the leading coefficient give a unit \(V(u,y)\), so \(P=VQ_a\) locally. At \((0,0)\),
\[
P_u=V(0,0)\beta'(0).
\]
The hypothesis \(P_u\ne0\) implies \(\beta'(0)\ne0\), and therefore
\[
\operatorname{disc}_y Q_a
=\alpha(u)^2-4\beta(u)
=-4\beta'(0)u+O(u^2)
\]
has order exactly one. Resultants between factors belonging to distinct constant roots are units. The product formula for discriminants gives (6), including its assertion that the discriminant is not identically zero. This proof uses only formal coefficient lifting and \(2\ne0\); it does not assume analytic or characteristic-zero Puiseux theory.

## 4. Every multiple root at \(u=0\) is a simple fold modulo \(p\)

Take \(n=p-1=2h\), and reduce \(F_n\) modulo \(p\). Its degree remains exactly \(n\), since \(n!\ne0\) in \(\mathbb F_p\). Write this reduction as \(\overline F\).

At \(u=0\),
\[
E_n(0)=(-1)^h=:\epsilon,\qquad
\overline F(0,0)=\epsilon\ne0.
\tag{7}
\]
Thus none of its roots in \(A\) is zero.

Telescoping (1) gives the universal integer identity
\[
A\,\partial_u F_n=F_n-E_n.
\tag{8}
\]
At a multiple root \((0,a)\), (7)--(8) yield
\[
\partial_u\overline F(0,a)=-\epsilon/a\ne0.
\tag{9}
\]

For clarity, the second-derivative assertion can be obtained entirely from the native polynomial at \(u=0\). Its exact integer expansion is
\[
F_n(0,A)=R_h(A^2),\qquad
R_h(X)=\sum_{i=0}^{h}
(-1)^{h-i}\frac{(h+i)!}{(h-i)!}X^i.
\tag{10}
\]
The coefficient ratio in (10) is
\[
c_i=-(h+i)(h-i+1)c_{i-1}.
\]
If \(\theta=X\,d/dX\), this is equivalent to
\[
R_h-\epsilon=-X\{h(h+1)-\theta-\theta^2\}R_h,
\]
or
\[
X^3R_h''+2X^2R_h'
-\{h(h+1)X+1\}R_h=-\epsilon.
\tag{11}
\]
At a multiple nonzero root \(a\) of \(F_n(0,A)\), put \(x=a^2\). Since \(2a\ne0\), one has \(R_h(x)=R_h'(x)=0\). Equation (11), reduced modulo \(p\), gives
\[
R_h''(x)=-\epsilon/x^3,\qquad
\partial_A^2\overline F(0,a)
=4xR_h''(x)=-4\epsilon/a^4\ne0.
\tag{12}
\]
All operations are legitimate in characteristic \(p\): the factorials in (10) have indices at most \(2h=p-1\), and no division by a multiple of \(p\) occurs.

Applying Section 3 over \(\overline{\mathbb F}_p\) proves
\[
\operatorname{ord}_{u=0}
\operatorname{disc}_A\overline F(u,A)
\le n/2<n.
\tag{13}
\]

## 5. Characteristic-zero contradiction

Let
\[
\Delta_n(u)=\operatorname{disc}_A F_n(u,A)\in\mathbb Z[u].
\]
Discriminants commute with reduction here because their defining universal polynomial has integer coefficients and the leading coefficient \(n!\) remains nonzero modulo \(p\). Thus (13) says
\[
0\ne\overline{\Delta_n}(u),\qquad
\operatorname{ord}_{u=0}\overline{\Delta_n}<n.
\tag{14}
\]

Suppose \(u_0\ne0\), \(D_p(u_0)=0\), and \(F_n(u_0,A)\) has a multiple root. Then \(\Delta_n(u_0)=0\). By irreducibility of \(M_p\), it divides \(\Delta_n\) over \(\mathbb Q[u]\). Because \(M_p\) is monic and both polynomials have integer coefficients, monic polynomial division shows that the quotient is in \(\mathbb Z[u]\). Reducing modulo \(p\) and using (5) forces
\[
u^n\mid\overline{\Delta_n}(u),
\]
contrary to (14). This proves the theorem.

## 6. All phases at every odd prime terminal

Let
\[
f=BXY+dX+eY,\qquad Bde\ne0,
\]
over any characteristic-zero field, with the factorial functional \(\Lambda(X^iY^j)=i!j!\). The native pure-cross reduction in the canonical
[BILINEAR_CHEBYSHEV_SQUAREFREE_TRANSFER.md](BILINEAR_CHEBYSHEV_SQUAREFREE_TRANSFER.md), says that a zero triple ending at moment \(N\) has, after adjoining \(s^2=de\),
\[
c=(d+e)/(2s),\quad b=B/(2s),\quad
T_N(c)=0,\quad G_{N-1}(c,b)=G_{N-1,b}(c,b)=0.
\tag{15}
\]

**Pure-cross prime-terminal theorem.** For every odd prime \(p\), the three moments
\[
\Lambda(f^{p-2}),\quad\Lambda(f^{p-1}),\quad\Lambda(f^p)
\]
cannot all vanish.

If \(d+e\ne0\), this follows from (15) and Sections 1--5. If \(d+e=0\), use exactly the opposite-axis theorem in the cited canonical source. Its proof identifies, for \(n=2m\),
\[
G_{2m}(0,b)=(2b)^{2m}
 B_{m,m+1}\bigl(-(2b)^{-2}\bigr),\qquad
B_{m,\alpha}(z)=\sum_{j=0}^{m}\binom mj(\alpha)_j j!\,z^{m-j}.
\tag{16}
\]
For \(m=1\) the polynomial \(B_{1,2}\) is linear. For \(m\ge2\), that canonical proof imports the squarefreeness theorem for \(B_{m,\alpha}\) at \(\alpha>m-1\), from
[MULTIVARIATE_ONE_PRIVATE_ROBIN_RESOLVENT.md](MULTIVARIATE_ONE_PRIVATE_ROBIN_RESOLVENT.md). The parameter here is \(\alpha=m+1\). The change \(b\mapsto-(2b)^{-2}\) has nonzero derivative at \(b\ne0\), so it preserves simplicity. This is an existing theorem being used for the central case, not an assertion that the reduction-modulo-\(p\) argument alone handles it.

For clarity on fields, the noncentral result is a coprimality statement over the algebraic coefficient field of \(c\), and hence over every characteristic-zero extension. The central assertion has the same characteristic-zero scope in its canonical source. Alternatively, any putative finite collection of coefficients can be placed in a finitely generated characteristic-zero field, embedded in \(\mathbb C\), and (15) preserved; the theorem then excludes that image.

In the necessary equation \((e/d)^p=-1\), the noncentral cases are the primitive \(2p\)-th roots, while \(e/d=-1\) is the central case. Thus no possible phase ratio has been omitted. For \(p=3\) the starting order \(p-2=1\) is positive and covered.

## 7. Scope

No finite computation is a premise of Sections 1--5. The noncentral proof is self-contained. The final all-phase prime-terminal theorem additionally uses the exact existing canonical opposite-axis squarefreeness theorem specified in Section 6.

This is not an all-start pure-cross theorem and does not close general binary quadratic SFC, SFC(2), or SFC(3). It is an infinite class of starting indices \(p-2\), with no coefficient restriction beyond the actual three-term hypothesis \(Bde\ne0\). The following doubled-prime result has an additional, explicit axis restriction.

## 8. Doubled primes away from the fourth-root axis ratios

**Theorem (doubled-prime nonexceptional windows).** Let \(p\) be any odd prime, \(n=2p-1\), and let \(c\) satisfy
\[
                 T_{2p}(c)=0,\qquad T_2(c)\ne0.
\]
Then \(G_n(c,b)\) is squarefree over \(\mathbb C\). Consequently, over every characteristic-zero field, if
\[
 f=BXY+dX+eY,\qquad Bde\ne0,\qquad d^2+e^2\ne0,
\]
the three factorial moments at powers \(2p-2,2p-1,2p\) cannot all vanish.

This is uniform in the prime and in all three coefficients. The excluded axis quotients are exactly \(e/d=\pm i\) after scalar extension. The theorem does not restrict the normalized cross coefficient to be real or to have a real conjugate.

Set \(u=2c\), \(A=2b\) and use \(F_n\) from (1). Composition of the Laurent identities for \(D_N\) gives
\[
 D_{2p}(u)=D_p(u^2-2),\qquad
 D_{2p}(u)=(u^2-2)M_p(u^2-2),\quad M_p(w)=D_p(w)/w.       \tag{17}
\]
Thus the nonexceptional resonances are exactly
\[
                         M_p(u^2-2)=0.                     \tag{18}
\]
It suffices to prove that \(F_{2p-1}(u,A)\) has no repeated root at any such algebraic \(u\).

## 9. A tame coefficient field and two possible root valuations

Embed the algebraic parameters in \(\overline{\mathbb Q}_p\), and normalize \(v(p)=1\). Put \(w=u^2-2\). Section 2 proves that \(M_p\) is Eisenstein, of degree \(p-1\), and hence
\[
                           v(w)=1/(p-1).                   \tag{19}
\]
For completeness this valuation follows directly from the constant coefficient of valuation one and the other nonleading coefficients of valuation at least one: the unique possible tie of minimal terms is between the leading term and the constant term.

The reduction \(\bar u=a\) satisfies \(a^2=2\). Choose a finite unramified extension \(L/\mathbb Q_p\) containing a square root \(a_0\) of 2 with that residue. In the complete field \(L(w)\), the simple residue root of \(X^2-(2+w)\) lifts uniquely, since \(2a\ne0\). Therefore the chosen \(u\) belongs to \(L(w)\). The polynomial \(M_p\) is still Eisenstein over \(L\), and
\[
                   v(L(w)^\times)=\frac1{p-1}\mathbb Z.    \tag{20}
\]
Indeed the powers \(1,w,\ldots,w^{p-2}\) form a basis; their valuations are distinct modulo \(\mathbb Z\), so a linear combination has the valuation of its unique least term. The same observation proves their independence and the Eisenstein degree assertion. Also \(u\) cannot belong to an unramified extension: \(u^2-2=w\) has the nonintegral valuation (19).

Write
\[
                    F_n(u,A)=\sum_{j=0}^n a_j(u)A^j,
              \qquad a_j(u)=E_n^{(j)}(u),\quad n=2p-1.
\]
Every coefficient of \(a_j(u)\) as an integer polynomial is divisible by \(j!\). Thus \(v(a_j(u))\ge0\) for \(j<p\) and \(v(a_j(u))\ge1\) for \(p\le j\le2p-1\).
The Laurent formula for \(E_n\), reduced modulo \(p\), gives
\[
 E_{2p-1}(u)\equiv u^p(u^2-4)^{(p-1)/2}\pmod p.           \tag{21}
\]
For example substitute \(u=z+z^{-1}\):
\((z^{2p}-z^{-2p})/(z-z^{-1})=u^p(z-z^{-1})^{p-1}\) in characteristic \(p\). Injectivity of that substitution proves (21).
Since \(a^2=2\), both \(a\) and \(a^2-4=-2\) are nonzero. Consequently
\[
 v(a_0(u))=0,\quad v(a_{p-1}(u))=0,\quad
                         v(a_{2p-1}(u))=1.                \tag{22}
\]
For the middle assertion differentiate (21) \(p-1\) times: the result is \((p-1)!u^p\) modulo \(p\). The last assertion uses \(a_n=n!\) and \(v((2p-1)!)=1\).

These bounds show directly that every root \(A\) has valuation
\[
                            v(A)=0\quad\hbox{or}\quad -1/p. \tag{23}
\]
If \(v(A)>0\), the constant term is uniquely least. If \(t=v(A)<0\), every term other than those with indices \(p-1\) and \(2p-1\) has valuation strictly above at least one of those two terms. Their valuations can tie only when
\((p-1)t=1+(2p-1)t\), or \(t=-1/p\). This argument requires no degree-preserving reduction of the whole polynomial modulo \(p\).

If \(v(A)=-1/p\), the powers \(1,A,\ldots,A^{p-1}\) are linearly independent over \(L(w)\): their valuations are distinct modulo the group (20). Thus \([L(w,A):L(w)]\ge p\). In characteristic zero an irreducible factor containing a repeated root must divide \(F_n\) at least twice. That would require degree at least \(2p\), whereas \(\deg_A F_n=2p-1\). No root of negative valuation can repeat.

## 10. Integral double roots would force an unramified parameter

It remains to exclude a repeated root of valuation zero. Let its residue be \(\alpha\ne0\). The telescoping identity (8) remains valid for this \(n\). The differentiated Chebyshev equation also gives the integer identity
\[
 -A^4F_{AA}-(3A^3+2uA^2)F_A
 +[n(n+2)A^2-uA+4-u^2]F
       =(4-u^2)E_n(u)-(n+1)A D_{n+1}(u).                  \tag{24}
\]
It follows by collecting the coefficient of every power of \(A\); in particular its two initial terms are included. This is also (A.5) of the real-conjugate-ratio note under \(u=2c,A=2b\).

Reduce at \((u,A)=(a,\alpha)\). Since \(n+1=2p\), \(\bar F=\bar F_A=0\) gives
\[
 \bar F_u=-\frac{E_n(a)}{\alpha}\ne0,\qquad
 \bar F_{AA}=-\frac{(4-a^2)E_n(a)}{\alpha^4}\ne0.           \tag{25}
\]
Here \(E_n(a)\ne0\) follows from (21), and \(4-a^2=2\ne0\). Therefore the two equations \(F=F_A=0\), in the two unknowns \(u,A\), have invertible Jacobian at the residue point: its determinant is \(\bar F_u\bar F_{AA}\).

Choose a further finite unramified extension \(L'\) containing the residues \(a,\alpha\). The simple-system version of Hensel lifting gives a solution \((\widehat u,\widehat A)\) in its valuation ring with these residues. It is the only such solution even after passing to any algebraic extension. The elementary argument is useful here. Newton correction with the inverse Jacobian doubles the positive valuation of the error at each step, giving a convergent solution in the complete unramified field. For uniqueness, subtract two solutions and let \(s>0\) be the minimum valuation of their coordinate differences. The invertible integral linear term has minimum valuation exactly \(s\), while every Taylor remainder has valuation at least \(2s\); they cannot sum to zero.

The hypothesized pair \((u,A)\) must consequently equal this unramified pair. But \(u\) cannot be unramified by (19). This excludes every remaining repeated root and proves squarefreeness in (18).

The ramification comparison is indispensable. The two-equation Jacobian argument forces a *critical parameter* to be unramified; it does not assert squarefreeness for every integral parameter with the same residue.

## 11. Phase conclusion and the exact remaining axis family

For a hypothetical zero triple ending at \(2p\), the native reduction (15) applies. With \(s^2=de\) and \(u=(d+e)/s\),
\[
                    u^2-2=\frac{d^2+e^2}{de}.              \tag{26}
\]
If \(d^2+e^2\ne0\), equations (17)--(18) put it in the proved squarefree sector. This contradicts the repeated-root requirement. The argument is algebraic: the resonant \(u\) and any repeated \(A\) are algebraic over \(\mathbb Q\), so their field embeds into \(\overline{\mathbb Q}_p\); the contradiction rules them out over every characteristic-zero field.

If \(d^2+e^2=0\), then \(u^2=2\), and \(u\) is unramified at the odd prime. Section 10 gives no contradiction there. The remaining doubled-prime obligation is precisely squarefreeness of
\[
                  F_{2p-1}(u,A),\qquad u^2=2,
                         \quad p\text{ any odd prime}.     \tag{27}
\]
This is a fixed axis-ratio family with unrestricted cross coefficient, not a finite set of coefficient profiles. Its \(p=3\) case was already proved by the exact quintic certificate in the first-five-starts note. No all-prime conclusion for (27) is inferred from that case.

The degree and tame-value-group hypotheses in Section 9 cannot be dropped. Over \(\mathbb Q_p\), \((1+pA^p)^2\) has repeated roots of valuation \(-1/p\), but has degree \(2p\). Over \(\mathbb Q_3(\pi)\), \(\pi^3=3\), the degree-three polynomial \((A-\pi^{-1})^2(A+\pi^{-1})\) has a repeated root of valuation \(-1/3\); its coefficient value group already contains thirds. In particular the proof above does not extend to prime-power terminals by replacing \(p\) with \(p^r\).

## 12. Reproducible supplemental checks

At \(p=3\), the resonance polynomial is \(M_3(u^2-2)=(u^2-2)^2-3\). For \(a^2=2\) in characteristic three, direct differentiation of \(E_5=u^5-4u^3+3u\) gives \(\bar F_5(a,A)=-a-A+aA^2\). Its unique double root is \(A=a\), with \(\bar F_u=1\) and \(\bar F_{AA}=2a\). This illustrates the unit fold explicitly; its lift is unramified and cannot equal a root of the displayed Eisenstein-in-\(u^2-2\) resonance. The negative-valuation branch has root degree at least three and cannot repeat in degree five.

The [doubled-prime checker](verify_bilinear_doubled_prime_ramification.py) independently constructs the integer Chebyshev recurrences. It checks (21), both differential identities, coefficient divisibility, and the three endpoint valuations for the explicitly recorded primes. It also checks exact discriminant coprimality with \(M_p(u^2-2)\) at three small primes, and rejects removal of the degree bound using \((1+3A^3)^2\). The [receipt](bilinear_doubled_prime_ramification_verification.json) records the environment, source hashes, precise finite domains and actual outcome. These checks are internal transcription and small-case controls; the proof for all odd primes is Sections 8--11. No independent replay, external review, or formal verification is claimed.

## 13. A uniform good-prime lift of the smaller factorial array

For \(m\ge2\), set
\[
 \Delta_{m-1}(u)=\operatorname{disc}_A F_{m-1}(u,A),\qquad
 R_m=\operatorname{Res}_u(D_m(u),\Delta_{m-1}(u))\in\mathbb Z,
                                                               \tag{28}
\]
and put \(R_1=1\). The discriminant of a degree-one polynomial is one, so \(R_2=1\). A prime is called good here only when \(p\nmid R_m\); no assertion that \(R_m\ne0\) for every \(m\) is assumed.

**Theorem (good-prime lifting).** Let \(m\ge1\), let \(p>m\) be an odd prime, and suppose \(p\nmid R_m\). If
\[
                    D_{mp}(u)=0,\qquad D_m(u)\ne0,
\]
then \(F_{mp-1}(u,A)\) is squarefree in characteristic zero. Consequently, for \(Bde\ne0\) and \(d^m+e^m\ne0\), the three factorial moments of \(BXY+dX+eY\) ending at \(mp\) cannot all vanish, over every characteristic-zero field.

Here \(m,p\) both vary. The explicitly testable smaller-array condition in (28) is a hypothesis. It is not replaced by an unproved all-index squarefreeness assertion.

To prove the theorem put \(n=mp-1\) and \(w=D_m(u)\). Composition gives \(D_p(w)=0\), and \(w\ne0\) gives \(M_p(w)=0\). Thus \(v(w)=1/(p-1)\). The residue \(a=\bar u\) is a root of \(D_m\). The identities
\[
 D_m'=mE_{m-1},\qquad
 D_m^2+(4-u^2)E_{m-1}^2=4
\]
show that \(a^2\ne4\) and \(E_{m-1}(a)\ne0\); indeed \(D_m(\pm2)=2(\pm1)^m\ne0\). As \(p>m\), \(D_m'(a)\ne0\). After an unramified extension \(L\) containing the residue \(a\), Hensel lifting of \(D_m(X)=w\) puts \(u\) in \(L(w)\), with value group \((p-1)^{-1}\mathbb Z\). The parameter \(u\) is not unramified because \(D_m(u)=w\) has nonintegral valuation.

The Laurent formula now gives the uniform congruence
\[
 E_{mp-1}(u)\equiv E_{m-1}(u)^p(u^2-4)^{(p-1)/2}\pmod p.  \tag{29}
\]
For \(a_j(u)=E_n^{(j)}(u)\), the factorial divisibility and \(n<p^2\) give \(v(a_j(u))\ge\lfloor j/p\rfloor\). Also \(a_0\) and \(a_{p-1}\) are units, while \(v(a_n)=m-1\). When \(m\ge2\), the lower coefficient polygon consists of the horizontal interval from 0 to \(p-1\), followed by the line of slope \(1/p\) ending at \(mp-1\). One can check this without a polygon theorem: for a negative candidate root valuation \(t\), the least lower bound in each block \(kp\le j<(k+1)p\) occurs at its right endpoint, and the successive endpoint heights differ by \(1+pt\). Hence a negative root must have \(t=-1/p\). Positive root valuations are excluded by the unit constant term. For \(m=1\), the leading coefficient is a unit and all roots have valuation zero.

For every unit-valued repeated root, equations (8) and (24), with \(n+1=mp\), again give the two unit derivatives in (25): \(E_n(a)\ne0\) and \(4-a^2\ne0\) by (29). The Hensel argument of Section 10 would make \(u\) unramified. Thus every possible repeated root has valuation \(-1/p\).

The negative-valuation roots retain the smaller factorial array exactly. For \(1\le k\le m\), coefficient extraction in (29) gives
\[
 \frac{a_{kp-1}(u)}{p^{k-1}}
       \equiv(-1)^k\bigl(E_{m-1}^{(k-1)}(a)\bigr)^p\pmod{\mathfrak m}.
                                                               \tag{30}
\]
Here division is legitimate coefficientwise over \(\mathbb Z_p\). To check the factor, the unit part of \((kp-1)!\) is
\[
 (kp-1)!/p^{k-1}\equiv(-1)^k(k-1)!\pmod p.
\]
In the Hasse derivative of order \(kp-1\), the only contributing term of \((u+t)^2-4\) raised to \((p-1)/2\) has degree \(p-1\), and the other factor contributes degree \((k-1)p\). Its coefficient is \((E_{m-1}^{[k-1]}(a))^p\). Since \(k-1<p\), multiplying by \((k-1)!\) gives (30).

Let \(Y=pA^p\). After division by \(A^{p-1}\), the terms lying on the negative-valuation edge therefore reduce to
\[
 r_a(Y)=-\sum_{j=0}^{m-1}(-1)^j
                         \bigl(E_{m-1}^{(j)}(a)\bigr)^pY^j. \tag{31}
\]
This is the coefficientwise Frobenius twist of \(-F_{m-1}(a,-Y)\), of degree \(m-1\), with nonzero constant term. The hypothesis \(p\nmid R_m\), and the fact that \((m-1)!\) is a unit, imply that \(F_{m-1}(a,A)\) is squarefree at every root \(a\) of \(D_m\). Thus \(r_a\) has \(m-1\) distinct nonzero roots.

Here is the full root-separation step, including the wild denominator \(p\). Enlarge \(L\) unramifiedly to contain the residue roots of \(r_a\), and set \(K=L(w)\); its value group is still \((p-1)^{-1}\mathbb Z\). Adjoin \(\pi\) with \(\pi^p=p\). The polynomial
\[
                    \pi^{p-1}F_n(u,Z/\pi)
\]
has integral coefficients and a unit leading coefficient. Its reduction is
\(Z^{p-1}r_a(Z^p)\). Over the residue field this has one factor \(Z^{p-1}\), and \(m-1\) pairwise coprime factors \((Z-\zeta)^p\), each associated to one root \(\eta=\zeta^p\) of \(r_a\). Coprime-factor Hensel lifting, by the coefficient-lifting argument of Section 3, shows that each nonzero residue cluster contains exactly \(p\) roots counted with multiplicity.

The cluster labeled by \(\eta\) is intrinsically the condition \(\overline{pA^p}=\eta\). It is preserved by every automorphism of \(\overline{\mathbb Q}_p/K\), because \(\eta\) lies in the residue field of \(K\). Its monic degree-\(p\) factor is therefore in \(K[A]\). Each root has valuation \(-1/p\); the distinct valuations modulo \(v(K^\times)\) of \(1,A,\ldots,A^{p-1}\) force its minimal polynomial over \(K\) to have degree at least \(p\). The cluster factor is consequently irreducible of degree \(p\), and is separable in characteristic zero. Its roots are simple in that cluster, and the different clusters are disjoint. This excludes all negative-valuation repeated roots. Separability of the residual polynomial alone, without this value-group and degree argument, would not suffice in a wildly ramified coefficient field.

Finally \(u=(d+e)/s\), \(s^2=de\), gives
\[
                       D_m(u)=\frac{d^m+e^m}{s^m}.
\]
The native zero-triple condition and the algebraic embedding argument of Section 11 prove the phase conclusion.

Direct exact evaluation of the small base resultants gives
\[
 \begin{array}{c|c}
 m&R_m\\ \hline
 2&1\\
 3&2^7\\
 4&2^{18}3^4\,23^2\\
 5&2^{51}3^{15}5^4\,7^4.
 \end{array}                                                \tag{32}
\]
Thus the theorem includes every terminal \(3p\) with odd prime \(p>3\), every \(4p\) with prime \(p>4\), \(p\ne23\), and every \(5p\) with prime \(p>7\), provided respectively \(d^3+e^3\), \(d^4+e^4\), or \(d^5+e^5\) is nonzero. The excluded lower-axis families are not closed by this lift. A bad prime in (32) is a failure of this sufficient modular test, not a counterexample to a factorial window.

The retained checker also verifies (29)--(31) on its listed \((m,p)\) pairs, computes (32) exactly, and rejects the bad-prime applications \((m,p)=(4,23),(5,7)\). These are finite checks of the uniform mechanism and the four literal constants in (32); they do not establish an unrestricted theorem by an initial-range extrapolation.

A concrete control shows why residual separability needs the tame coefficient field. Let \(\pi^3=3\), \(b=\pi^{-1}\), and work over \(\mathbb Q_3(\pi)\). The degree-five polynomial
\[
 (3A^3-3\pi A+2)(A^2+1)
                  =3(A-b)^2(A+2b)(A^2+1)
\]
has the same coefficient lower bounds and the same two root valuations as the doubled-prime case \(p=3,m=2\). Its negative-edge polynomial is the separable linear polynomial \(2+Y\), \(Y=3A^3\), but \(b\) is a double root. Its coefficient value group already contains thirds. This refutes the broader residual-separability shortcut, not a factorial-window conjecture. The checker verifies the double root and the simple residual polynomial exactly.

## 14. The exact bad-prime square-cluster obstruction

**Proposition (perfect-square cluster reduction).** Let \(m\ge1\), let \(p>m\) be an odd prime, and let \(u\) satisfy \(D_{mp}(u)=0\ne D_m(u)\). Do not assume \(p\nmid R_m\). After a finite unramified base extension, put \(K=L(D_m(u))\) as in Section 13, so \(v(K^\times)=(p-1)^{-1}\mathbb Z\), and include all roots of the residual polynomial \(r_a\) in the residue field. Then every repeated root of \(r_a\) is exactly double. A repeated root of \(F_{mp-1}(u,A)\) exists if and only if one of the monic degree-\(2p\) Hensel factors attached to those double residual roots is a perfect square in \(K[A]\). In that event the factor is \(Q(A)^2\), with \(Q\) irreducible of degree \(p\). All its \(p\) distinct roots are double roots of \(F_{mp-1}\).

Here the cluster factor is the factor containing exactly the roots with a prescribed nonzero residue \(\overline{pA^p}=\eta\). It is an actual factor of the characteristic-zero polynomial, not a freely chosen polynomial with the same leading reduction.

To prove the residual assertion, set \(F=F_{m-1}\) and work in characteristic \(p\) at a residue \(a\) of \(u\). As in Section 13, \(D_m(a)=0\), \(E_{m-1}(a)\ne0\), and \(4-a^2\ne0\). Also \(F(a,0)\ne0\). At a nonzero repeated root \(\beta\), identities (8) and (24), now with \(n=m-1\), give
\[
 F_u(a,\beta)=-E_{m-1}(a)/\beta\ne0,\qquad
 F_{AA}(a,\beta)=-(4-a^2)E_{m-1}(a)/\beta^4\ne0.          \tag{34}
\]
Since \(p\ne2\), the multiplicity is exactly two. Equation (31) is a coefficientwise Frobenius twist and a sign change of variable, so it preserves root multiplicities. Thus each root of \(r_a\) has multiplicity one or two. The first derivative in (34) concerns the untwisted smaller array; differentiating its Frobenius twist with respect to \(a\) would instead give zero in characteristic \(p\).

The scaled polynomial \(\pi^{p-1}F_{mp-1}(u,Z/\pi)\), \(\pi^p=p\), reduces to \(Z^{p-1}r_a(Z^p)\). A root \(\eta\) of \(r_a\) of multiplicity \(e\in\{1,2\}\) consequently has a coprime residue cluster of degree \(ep\). The Hensel and Galois argument of Section 13 puts its monic cluster factor \(C_\eta\) in \(K[A]\), of that exact degree. Every root in this cluster has valuation \(-1/p\), hence minimal-polynomial degree at least \(p\) over \(K\).

Integral repeated roots have already been excluded without a good-prime hypothesis. Simple residual clusters have degree \(p\) and are irreducible and separable by the same degree argument. Suppose now a degree-\(2p\) cluster contains a repeated root with monic minimal polynomial \(Q\) of degree \(d\ge p\). Separability in characteristic zero implies \(Q^2\mid C_\eta\). Hence \(2d\le2p\), forcing \(d=p\) and \(C_\eta=Q^2\). Conversely, if a monic \(C_\eta\) is a square in \(K[A]\), its monic square root has degree \(p\). The valuation degree bound makes that square root irreducible; its characteristic-zero roots are simple, and each becomes exactly double in the cluster. This proves both directions. It is a reduction, not an exclusion of the square factors.

### A low-order discriminant shortcut fails

The order argument of Sections 3--5 does not extend to every lifted terminal merely by removing the integer content of the discriminant. There is an exact obstruction at the bad pair \((m,p)=(5,7)\) in (32).

**Proposition (bad-lift discriminant order, computational).** Put
\[
 \Delta(u)=\operatorname{disc}_A F_{34}(u,A)\in\mathbb Z[u],
 \qquad \Delta=c\Delta_{\rm prim},
\]
where \(c>0\) is the gcd of its integer coefficients. Then \(\deg\Delta=1122\), \(v_7(c)=136\), and
\[
 \operatorname{ord}_{q}\bigl(\overline{\Delta_{\rm prim}}\bigr)=49,
 \qquad q(u)=u^4+2u^2-2\in\mathbb F_7[u].                  \tag{33}
\]
Here \(q\) is the noncentral factor of \(\overline{D_5}=u q\); the order means the largest exponent for divisibility by the entire polynomial \(q\). In particular it is not less than \(p-1=6\). The reduction at the central factor \(u\) is a unit.

The [exact checker](verify_bilinear_bad_prime_discriminant_order.py) constructs \(E_{34}\) from the integer Chebyshev recurrence, forms \(F_{34}\), computes its discriminant over \(\mathbb Z[u]\), removes its full integer content, and performs exact repeated division modulo seven. The [receipt](bilinear_bad_prime_discriminant_order_verification.json) retains the input coefficients, primitive discriminant modulo seven, quotient by \(q^{49}\), and nonzero next remainder. Multiplying the quotient back by \(q^{49}\) checks the stated order independently of the repeated-division loop. This is one exact finite computational assertion, with the stated environment and no independent-replay or formal-verification claim.

At a noncentral \(D_5\) residue, the Eisenstein lift \(M_7(D_5(u))=0\) has local ramification degree six. Its mod-seven multiplicity therefore cannot contradict the order 49 in (33). A higher-precision argument retaining the actual lifted factorial array is needed if the good-prime assumption is to be removed this way. Neither (33) nor \(7\mid R_5\) asserts a repeated characteristic-zero root or a factorial counterexample.

