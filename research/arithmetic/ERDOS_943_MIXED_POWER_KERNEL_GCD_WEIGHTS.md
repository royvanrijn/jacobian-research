# Recovering the kernel gcd for fixed mixed powers

Fix integers \(2\leq k<j\), \(3\mid j\), \(k\nmid j\).
Let \(n\geq2\). A simple bipartite graph has squarefree positive
left kernels \(b\leq K_1\), right kernels \(d\leq K_2\), and one
selected positive representation on every edge:
\[
 \kappa_1a^kb^j+\kappa_2c^kd^j=n,\qquad
 \gcd(\kappa_1a^kb^j,\kappa_2c^kd^j)=1,\quad
 a\leq R_1,\quad c\leq R_2.                              \tag{1}
\]
The positive integer coefficient pair is fixed during the graph
count. All four real caps are in \([1,n]\).
Constants may depend on \(k,j,\epsilon\), but not on the target,
caps or coefficient pair. No aggregate over coefficient pairs
or exponents is counted.

For a shared-right path \(b_1-d-b_2\), with distinct \(b_1,b_2\), put
\[
 a_i=\alpha r_i,\quad\alpha=\gcd(a_1,a_2),\qquad
 b_i=g u_i,\quad g=\gcd(b_1,b_2),
\]
\[
 L=\frac{c_1^k-c_2^k}{\alpha^kg^j},\qquad
 e=\operatorname{pf}_k(L),\qquad
 \mathcal V_{k,j}=\frac{\alpha^3g^{\,j/k}}{e^{\,1-1/k}}.     \tag{2}
\]
Here \(\operatorname{pf}_k(L)\) is the positive \(k\)-power-free part:
\(|L|=e v^k\), \(v>0\) integral, and \(0\leq v_p(e)<k\).
The positive real powers in weights are only counting weights.
The quantities in (2) are defined on every selected path.
For the fixed left coefficient define
\[
 \chi_k(\kappa_1)=\prod_p p^{\lceil v_p(\kappa_1)/k\rceil},
 \qquad \chi_k(1)=1.
\]

**Theorem 1 (No kernel-gcd index).** For \(A,G,E\geq1\),
\[
 \#\{\text{selected paths}: A\leq\alpha<2A,\ G\leq g<2G,\ e\leq E\}
 \ll_{k,j,\epsilon}
 \frac{n^\epsilon R_1^2R_2 E^{\,1-1/k}}{\chi_k(\kappa_1)A^3G^{\,j/k}}.       \tag{3}
\]
In particular
\[
 \sum_{\text{selected shared-right paths}}\mathcal V_{k,j}
                 \ll_{k,j,\epsilon}n^\epsilon R_1^2R_2/\chi_k(\kappa_1). \tag{4}
\]
These path estimates need no small-height gate.

**Theorem 2 (Actual cycle weight).** If in addition
\[
                         2(R_1R_2)^2<n,                  \tag{5}
\]
then
\[
 \sum_{C_6}\sum_{\text{shared-right paths in }C_6}\mathcal V_{k,j}
                 \ll_{k,j,\epsilon}n^\epsilon K_1R_1^2R_2/\chi_k(\kappa_1), \tag{6}
\]
and, for every real \(Y>0\),
\[
 \#\{C_6:\text{some shared-right path has }\mathcal V_{k,j}\geq Y\}
       \ll_{k,j,\epsilon} n^\epsilon K_1R_1^2R_2/(\chi_k(\kappa_1)Y). \tag{7}
\]
Cycles are simple and use the selected actual labels.
Side-interchanged versions hold. These are sector counts, not a
full mixed-power representation theorem.

## 1. The integral reduced cubic

Actual coprimality in (1) implies all variables and coefficients
on an edge are units modulo \(n\). It also gives
\(\gcd(\alpha g,\kappa_2d)=\gcd(\kappa_1,\kappa_2d)=1\).
Subtracting the edge equations first shows
\(\alpha^k\mid c_1^k-c_2^k\), and then
\(g^j\mid(c_1^k-c_2^k)/\alpha^k\). Hence \(L\) is integral and
\[
 \kappa_1r_1^ku_1^j-\kappa_1r_2^ku_2^j+\kappa_2Ld^j=0.   \tag{8}
\]
Equation(8) and \(\gcd(\kappa_1,\kappa_2d)=1\) also give
\(\kappa_1\mid L\), an indispensable coefficient restriction below.
It is nonzero. Indeed, for a positive squarefree \(b\), the prime
valuations of \(a^kb^j\) are \(k m+j\delta\),
\(\delta\in\{0,1\}\). The two residues \(0,j\pmod k\) differ
precisely under \(k\nmid j\). They uniquely determine \(b\) and
then \(a\), even if \(\gcd(k,j)>1\) or \(a,b\) are not coprime.
Thus \(c_1=c_2\) would force \(b_1=b_2\), a contradiction.

The integers \(u_1,u_2\) are distinct, coprime, positive and
squarefree; each is coprime to \(d\). Therefore
\[
 (u_1^{j/3},u_2^{j/3},d^{j/3})
\]
is a primitive integer triple. The powers are integral because
\(3\mid j\), and their maximum height is at most \(n^{j/3}\).

The first classical input is the exact Heath-Brown
arXiv2601.07817v1 Theorem1.3: primitive zeros of
\(AX^3+BY^3+CZ^3=0\), for arbitrary nonzero signed integer
coefficients, of maximum height at most \(H\) number
\(O_\eta(H^\eta)\), uniformly in the coefficients.
Thus, after fixing ONLY \((r_1,r_2,L)\), equation(8) gives
\(O_{k,j,\epsilon}(n^\epsilon)\) reduced kernel triples
\((u_1,u_2,d)\), by resizing the exponent.
Every coefficient is nonzero. This supplies completions of the
coefficient slice; the frequency of those slices is still to be counted.

## 2. An exact fixed-level lifting problem

Fix \((r_1,r_2,L)\) and one reduced kernel triple capable of
producing an actual path. Write
\[
 A_i=\kappa_1r_i^ku_i^j>0,\qquad S=\alpha^kg^j.
\]
The actual equations are
\[
 A_iS+\kappa_2c_i^kd^j=n.
\]
Using (8) to eliminate \(S\) gives
\[
                         A_2c_1^k-A_1c_2^k=nL.           \tag{9}
\]
The right side is nonzero, of absolute value at most \(n^{k+1}\),
because \(|L|<R_2^k\leq n^k\).

We prove that (9) has at most \(O_{k,j,\epsilon}(n^\epsilon)\)
actual-capable positive root pairs, uniformly in all the fixed
coefficient and reduced-kernel data.

### 2.1. Degrees at least three: exact uniform Thue--Mahler

Suppose \(k\geq3\).
Let \(c=\gcd(A_1,A_2)\). If \(c\nmid nL\), no lift exists.
Otherwise put
\[
 F(X,Y)=(A_2X^k-A_1Y^k)/c,\qquad T=nL/c.
\]
This is a primitive integral binary form of degree exactly \(k\),
with \(k\) DISTINCT complex projective roots: the coefficients
are nonzero and characteristic zero makes the roots of a nonzero
binomial separable. Reducibility over \(\mathbb Q\) causes no problem.

For a root pair, let \(\gamma=\gcd(c_1,c_2)>0\) and
\(c_i=\gamma s_i\). Homogeneity gives
\(\gamma^k\mid T\) and
\[
              \gcd(s_1,s_2)=1,\qquad F(s_1,s_2)=T/\gamma^k.
\]
There are at most \(\tau(|T|)\) possible \(\gamma\).

The precise second classical input is Evertse's1985 Theorem6:
for an integral binary form of degree \(d\geq3\) with at least
three distinct complex linear factors, the number of primitive
integer pairs whose nonzero absolute form value is supported on
a fixed set of \(t\) primes is at most
\[
                              7^{d^3(2t+3)},              \tag{10}
\]
independently of the coefficients. This includes reducible and
repeated-root forms satisfying the distinct-factor condition.
Apply it with \(d=k\) and the prime divisors of \(|T|\).
Every exact value \(T/\gamma^k\), of either sign, is included.
Thus each \(\gamma\) contributes at most
\(7^{k^3(2\omega(|T|)+3)}\).

For fixed \(k,\eta\), both this bound and \(\tau(|T|)\) are
\(O_{k,\eta}(|T|^\eta)\). For example a fixed constant to the
power \(\omega(m)\) is \(O_\eta(m^\eta)\), by separating the
finitely many small primes and using \(C\leq p^\eta\) for the
remaining primes. Since \(|T|\leq n^{k+1}\), choosing smaller
initial exponents proves the desired \(n^\epsilon\) bound.
No bound on the heights of the coefficients of \(F\) is needed,
and no primitive coordinate, polynomial content, target-zero,
or leading-degree-drop case is discarded.

### 2.2. Degree two: the real-quadratic lemma

If \(k=2\), the hypothesis \(k\nmid j\) makes \(j\) odd.
Actual coprimality gives \(\gcd(\kappa_1,\kappa_2d)=1\), so (8)
implies \(\kappa_1\mid L\) on any actual-capable reduced triple.
Set
\[
 D=u_1u_2,\quad
 x=r_2u_2^{(j+1)/2}c_1,\quad
 y=r_1u_1^{(j-1)/2}c_2,\quad
 M=u_2nL/\kappa_1.
\]
Multiplying (9), after division by \(\kappa_1\), by \(u_2\)
gives exactly
\[
                              x^2-Dy^2=M.                \tag{11}
\]
Here \(D>1\) is squarefree and nonsquare, since \(u_1,u_2\) are
distinct, coprime squarefree integers. The target \(M\) is a
nonzero integer. Bounds adequate for uniformity are
\[
 D\leq n^2,\qquad |M|\leq n^4,\qquad
 |x|,|y|\leq n^{(j+5)/2}.                                \tag{12}
\]

The needed elementary norm lemma, proved in the companion
*Recovering the kernel gcd by an elementary real-quadratic norm count*,
Lemma1, says uniformly for squarefree \(D>1,M\ne0\),
\[
 \#\{x^2-Dy^2=M:\ |x|,|y|\leq H\}
 \ll \tau(|M|)2^{\omega(|M|)}
             [1+\log(H(1+\sqrt D))].                    \tag{13}
\]
Its argument removes \(\gcd(x,y)\), uses at most
\(2^{\omega(|M'|)+1}\) primitive residue classes
\(x/y\pmod{|M'|}\), and shows two solutions in one class have
an integral norm-one ratio. Its positive absolute embedding
ratios are separated by \((3+\sqrt5)/2\), from their integral
trace. This proves the logarithmic factor without ideal
factorization, a class-number estimate or Dirichlet's unit theorem.
Ramified primes, the dyadic case, either sign of \(M\), and
\(|M'|=1\) are included in that lemma.

By (12)--(13), there are \(O_{j,\epsilon}(n^\epsilon)\) pairs
\(x,y\), and hence at most that many positive \(c_1,c_2\), since
the displayed fixed nonzero factors recover them injectively.
This proves the same root-lift bound as in the \(k\geq3\) case.

## 3. Recovering the omitted kernel gcd

For either case, each root pair fixes
\[
 S=\frac{n-\kappa_2c_1^kd^j}{A_1}.                        \tag{14}
\]
Discard nonintegral or nonpositive values.
A positive integer has at most one expression
\(S=\alpha^kg^j\) with \(\alpha\) positive integral and \(g\)
positive squarefree: the two distinct residues \(0,j\pmod k\)
determine \(g\) prime by prime, and then determine \(\alpha\).
This uses \(k\nmid j\), not the stronger \(\gcd(k,j)=1\).

Consequently no \(g\), \(\alpha\), root sign or representation
multiplicity is added after the reduced cubic and the root lifts.
Recover \(b_i=gu_i\) and retain only candidates satisfying every
original squarefree, gcd, cap, positivity and selected-label condition.
Discarding invalid candidates only decreases the count.
In particular all completions of a fixed \((r_1,r_2,L)\)
number \(O_{k,j,\epsilon}(n^\epsilon)\).

This lifting lemma counts actual completions AFTER the coefficient
data are fixed. It is the reason the following frequency estimate
does not sum separately over the integers \(g\).

## 4. Coefficient frequency and weighted paths

On the dyadic conditions in Theorem1,
\[
             r_i\leq R_1/A,\qquad
             0<|L|<L_0=R_2^k/(A^kG^j).
\]
The number of possible nonzero signed \(L\) with
\(\operatorname{pf}_k(L)\leq E\) is at most
\[
 2\sum_{1\leq t\leq E}\lfloor(L_0/t)^{1/k}\rfloor
       \leq\frac{2k}{k-1}L_0^{1/k}E^{\,1-1/k}.            \tag{15}
\]
Dropping \(k\)-power-freeness of \(t\) only enlarges this bound.
Terms with \(t>L_0\) vanish; \(L_0<1\) and \(A>R_1\) give
empty sets. All inequalities hold for real \(E\geq1\).

To obtain the coefficient gain, retain the necessary divisibility
\(\kappa_1\mid L\). Write
\[
 v_p(\kappa_1)=k r_p+s_p,\quad0\leq s_p<k,\qquad
 h=\prod_p p^{r_p},\quad K=\prod_{s_p>0}p,\quad
 \chi=\chi_k(\kappa_1)=hK.
\]
For \(e\) \(k\)-power-free define
\[
 Q(e)=\prod_{\substack{p\mid K\\v_p(e)\geq s_p}}p.
\]
The prime-valuation condition
\(\kappa_1\mid e v^k\) is equivalent to
\[
                         \chi/Q(e)\mid v.                \tag{15a}
\]
For \(s_p=0\), it forces \(v_p(v)\geq r_p\).
For \(s_p>0\), it forces \(v_p(v)\geq r_p+1\) when
\(v_p(e)<s_p\), and at least \(r_p\) otherwise.
Thus every local case in (15a) is accounted for.

For \(a\mid K\), put \(B_a=\prod_{p\mid a}p^{s_p}\).
The exact finite divisor expansion is
\[
 Q(e)=\sum_{a\mid K}\varphi(a)\,\mathbf1_{B_a\mid e},
\]
because \(K\) is squarefree. Dropping \(k\)-power-freeness only
enlarges the following positive sum. For every real \(E\geq1\),
\[
 \begin{split}
 \sum_{\substack{e\leq E\\e\ k\text{-power-free}}}
                 e^{-1/k}Q(e)
 &\leq\sum_{a\mid K}\varphi(a)
                  \sum_{\substack{e\leq E\\B_a\mid e}}e^{-1/k}\\
 &\leq\frac{k}{k-1}E^{\,1-1/k}
                   \sum_{a\mid K}\frac{\varphi(a)}{B_a}\\
 &\leq\frac{k}{k-1}E^{\,1-1/k}\tau(K).                   \tag{15b}
 \end{split}
\]
Indeed the inner sum is
\(B_a^{-1/k}\sum_{m\leq E/B_a}m^{-1/k}
 \leq[k/(k-1)]E^{1-1/k}/B_a\), with an empty sum if \(B_a>E\).
Also
\(\sum_{a\mid K}\varphi(a)/B_a
=\prod_{p\mid K}(1+(p-1)/p^{s_p})\leq2^{\omega(K)}=\tau(K)\).

For each such \(e\), (15a) bounds possible positive \(v\) by
\(L_0^{1/k}e^{-1/k}Q(e)/\chi\).
Therefore the number of signed \(L\) relevant to actual paths is
\[
 \ll_k\frac{\tau(K)}{\chi}\,
                 L_0^{1/k}E^{\,1-1/k}.                  \tag{15c}
\]
This is valid for either sign, for \(\kappa_1=1\), and for all
empty/floor ranges. A nonempty graph has \(\kappa_1<n\), so
\(\tau(K)\ll_\eta n^\eta\), uniformly, for each \(\eta>0\).

There are at most \((R_1/A)^2\) ordered pairs \(r_1,r_2\).
Multiply their count by (15c) and the completion bound of
Sections1--3, resizing the initial epsilons to absorb \(\tau(K)\).
This proves (3). There is no \(g\)-index in this coefficient set.

For (4), split \(\alpha,g,e\) into dyadic ranges.
They are bounded respectively by \(n,n,n^k\), giving
\(O_k((1+\log n)^3)\) range triples.
On a block with \(e\in[E_0,2E_0)\), the weight is
\(O_{k,j}(A^3G^{j/k}/E_0^{1-1/k})\).
This cancels the dyadic variables in (3), used with \(E=2E_0\),
while retaining the fixed factor \(1/\chi_k(\kappa_1)\).
A smaller initial epsilon absorbs the logarithms.
This proves the path theorem without (5).

## 5. Actual codegrees and cycles

For completeness the general actual-codegree argument under (5)
has no new deep input. With
\(c_i=\gamma s_i\), actual within-edge coprimality gives
\[
 \gcd(r_1,r_2)=\gcd(s_1,s_2)=\gcd(r_1,s_1)=\gcd(r_2,s_2)=1.
\]
The reduced positive fraction
\[
             u/v=(a_1c_2)/(a_2c_1)
                 =(r_1s_2)/(r_2s_1)
\]
has the coprime numerator and denominator exactly \(r_1s_2,r_2s_1\).
They are units modulo \(n\) of height at most \(R_1R_2\), and
\[
                         b_1^ju^k\equiv b_2^jv^k\pmod n.
\]
For fixed \(k\), unit \(k\)-th-root classes number at most
\(C_k^{\omega(n)}\ll_{k,\epsilon}n^\epsilon\).
Each class contains at most one reduced fraction by (5).
Divisor-many allocations recover \(r_i,s_i\), and
\[
 \kappa_1\alpha^k(b_1^ju^k-b_2^jv^k)=n(s_2^k-s_1^k)
\]
recovers \(\alpha\) at most once.
The bracket cannot vanish: squarefree kernel valuation differences
are \(0,\pm1\), and \(k\nmid j\) makes the zero congruence force
\(b_1=b_2\).
The first edge then determines \(\gamma^kd^j\), whose
squarefree-\(d\) decomposition is unique by the same residue
argument. Thus every same-side pair has \(n^\epsilon\) selected
actual common neighbors, uniformly; interchange the sides for the
other orientation.

Each fixed selected path has \(O_{k,j,\epsilon}(K_1n^\epsilon)\)
simple-cycle continuations by choosing the third left kernel and
two such common neighbors. This and (4) prove (6).
Positivity of the weight proves (7) for every \(Y>0\).

## 6. Exact gain and boundaries

Let the preceding general weight be
\[
 \mathcal W_{k,j}
       =\alpha^3g^{j/k-1}/e^{1-1/k}.
\]
The new one is exactly
\[
                         \mathcal V_{k,j}=g\mathcal W_{k,j}. \tag{16}
\]
Put \(\ell=(c_1^k-c_2^k)/\alpha^k=g^jL\).
For \(\gamma=\gcd(c_1,c_2)\), actual coprimality gives
\(\gcd(\gamma,\alpha g)=1\), so \(\gamma^k\mid L\).
Writing \(|L|=e v^k\) then implies \(\gamma\mid v\).
The natural unrefined quotient \(e_0=|\ell|/\gamma^k\) and its
weight \(W_0=\alpha^3/e_0^{1-1/k}\) satisfy
\[
       \frac{\mathcal V_{k,j}}{W_0}
                    =g^j(v/\gamma)^{k-1}\geq g^j.        \tag{17}
\]
This ratio is a positive integer.

The positive root cap gives
\[
 \alpha g^{j/k}e^{1/k}v<R_2,\qquad
 \alpha^{k+2}g^j<R_2^{k-1}\mathcal V_{k,j},\qquad
 \mathcal V_{k,j}<\alpha^2R_2.                           \tag{18}
\]
Thus every actual weight is \(>R_2^{-(k-1)}\), and membership in
\(\mathcal V_{k,j}\geq Y\) requires
\(\alpha>\sqrt{Y/R_2}\).
These are genuine necessary bounds and counted sectors; they do
not force a large weight on every residual cycle.

For \(k=2,j=3\), the weight is precisely
\(\alpha^3g^{3/2}/\sqrt{e}\), the companion Pell theorem.
Its actual two-edge example
\(67193=38^3+111^2=34^3+167^2\) has ratio8 to \(W_0\),
compared with ratio4 for \(\mathcal W_{2,3}\).
It is an actual path, not a six-cycle or a sharpness construction.

## 7. The exact square-cube specialization

**Corollary 3.** For \(k=2,j=3\), write uniquely
\(\kappa_1=h^2q\) with \(q\) positive squarefree.
Then \(\chi_2(\kappa_1)=hq\), and Theorems1--2 give
\[
 \sum_{\text{paths}}\frac{\alpha^3g^{3/2}}{\sqrt e}
       \ll_\epsilon \frac{n^\epsilon R_1^2R_2}{hq},
 \qquad
 \sum_{C_6}\sum_{\text{paths in }C_6}
          \frac{\alpha^3g^{3/2}}{\sqrt e}
       \ll_\epsilon \frac{n^\epsilon K_1R_1^2R_2}{hq}.
\]
The tail has the same numerator divided by \(hqY\), for every
\(Y>0\). These square-cube bounds use the local uniform
real-quadratic norm lemma and the primitive cubic theorem.
The casewise proof in Sections1--5 does not require the
degree-at-least3 Thue theorem when \(k=2\), nor any separate
square-cube path or cycle theorem.

For the existing shaped square-cube caps the resource is
\(K_1R_1^2R_2
=16n^{131/153}/(\kappa_1^{19/9}\kappa_2^{14/9})\).
Under the existing height gate, the \(n^{26/51+\epsilon}\) cycle
budget follows on
\[
 \frac{\alpha^3g^{3/2}}{\sqrt e}\geq
 \frac{n^{53/153}}
 {\kappa_1^{19/9}\kappa_2^{14/9}hq}.
\]
For sufficiently large \(n\), those shaped caps satisfy (5)
uniformly, since \(2(R_1R_2)^2\leq32n^{44/51}<n\).
The bounded initial nonempty range is uniform because all labels,
kernels and coefficients are bounded by its target bound.
The balanced cap specialization has tail
\(O_\epsilon(n^\epsilon T^4/(hqY))\).
No statement that every residual cycle meets these sectors is made.

For a squarefree coefficient \(\kappa_1\), the general theorem
has \(\chi_k(\kappa_1)=\kappa_1\), for every allowed \(k,j\).
This is a coefficient gain in a fixed graph, not permission to
sum arbitrary coefficient pairs.

## 8. Exact inputs and scope

1. Heath-Brown, arXiv2601.07817v1, Theorem1.3, printed4:
primitive diagonal-cubic zeros, arbitrary nonzero signed integer
coefficients, uniform height-subpower count. It is applied only
after fixing \((r_1,r_2,L)\). No global triple theorem is used.
2. J.-H. Evertse, *Upper bounds for the numbers of solutions of certain
Diophantine equations* (1985), Theorem6, printed84, with its supplied
derivation printed85: degree at least3, at least three DISTINCT
complex projective roots, primitive integer coordinates, uniform
bound(10). This is the retained1985 three-root theorem, not the
different1997 irreducible-form theorem. Its supplied proof imports
the stated S-unit Theorem5; that underlying proof is not re-audited here.
3. The companion elementary real-quadratic norm lemma, Lemma1,
is used only at \(k=2\); its complete proof includes every local
root class and unit-spacing case. The actual general codegrees
are re-derived in Section5.

The fixed exact-level lifting problem recovers the kernel gcd
instead of summing it as a free index, and the necessary divisor
\(\kappa_1\mid L\) supplies the explicit coefficient gain.
There is no uniform-in-unbounded-\(k,j\) assertion, unrestricted
coefficient-pair sum, unselected-edge multiplicity estimate,
full six-cycle theorem or full \(k\)-full prescribed-sum conclusion.
