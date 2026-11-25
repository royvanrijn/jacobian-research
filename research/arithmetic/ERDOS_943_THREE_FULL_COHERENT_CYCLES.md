# Three-full sums from coherent six-cycles

Let \(r_3(N)\) be the number of ordered pairs of positive 3-full integers with sum \(N\). An integer is 3-full if every prime in its factorization has exponent at least three; the integer one is included.

**Theorem 1.** For every \(\epsilon>0\) and every positive integer \(N\),
\[
                         r_3(N)\ll_\epsilon N^{33/134+\epsilon}.       \tag{1}
\]
More precisely, for fixed positive integers \(\kappa_1,\kappa_2,n\), the number of canonical tuples in
\[
 \kappa_1a_0^3a_1^4a_2^5+\kappa_2c_0^3c_1^4c_2^5=n,\qquad
 \gcd(\kappa_1a_0^3a_1^4a_2^5,\kappa_2c_0^3c_1^4c_2^5)=1             \tag{2}
\]
is \(O_\epsilon(n^{33/134+\epsilon})\), uniformly in the positive coefficients. Here \(a_1,a_2\) are squarefree and coprime, and \(c_1,c_2\) are squarefree and coprime. The cube coordinates may share primes with their own fourth- or fifth-power coordinates.

The cycle argument retains the actual positive representations. A coherent cycle has primitive integer vertex potentials. Complementing these potentials turns a missing-root recovery into a unit cubic congruence. Counting both the original and complementary gcd parameters gives a cycle bound in every dyadic root rectangle. Combining it with two omitted-coordinate graphs and a mixed quartic curve gives (1).

## 1. Canonical factors and the counting inputs

Every positive 3-full integer has a unique presentation
\[
                              m=a_0^3a_1^4a_2^5.                    \tag{3}
\]
For a prime exponent \(v\ge3\), its residue modulo three determines the presentation: if \(v\equiv0\), put the whole exponent in \(a_0^3\); if \(v\equiv1\), put one prime in \(a_1\) and the remaining exponent in \(a_0^3\); if \(v\equiv2\), put one prime in \(a_2\) and the remainder in \(a_0^3\). The respective cube exponents are \(v/3,(v-4)/3,(v-5)/3\). They are nonnegative in their applicable cases. Replacing a coordinate by a divisor preserves the canonical conditions.

We use the following three counting theorems.

* Evertse's 1985 Theorem 6 bounds primitive solutions to a fixed nonzero value of an integer binary form with at least three distinct complex projective linear factors by a degree-dependent constant to the power \(1+\omega(n)\). The bound is uniform in the coefficients. In particular, for fixed \(j=3\) or \(4\), the primitive positive solutions of \(P x^j+Qy^j=n\), with \(P,Q>0\), number \(O_\epsilon(n^\epsilon)\), uniformly in \(P,Q\). No irreducibility of this binary form is assumed. The form has \(j\) distinct projective roots because \(P,Q\ne0\). This is the [1985 theorem](sources/evertse-uniform-thue-mahler-1985.pdf), rather than the distinct 1997 irreducible-form theorem.
* Bombieri--Pila's 1989 Theorem 5 bounds integer points on an absolutely irreducible affine curve of fixed degree \(d\) in a square of side \(H\) by \(O_{d,\eta}(H^{1/d+\eta})\), uniformly in its coefficients. The degree-dependent bounded range of \(H\) is absorbed in the constant. We use its [degree-four case](sources/bombieri-pila-arcs-ovals-1989.pdf).
* Morris--Saxton's Theorem 1.5(a), with cycle parameter three, implies that a simple graph on \(v\) vertices with \(E\) edges satisfies
\[
                       E\ll v^{4/3}+v\,C_6(G)^{1/6}.               \tag{4}
\]
Above a constant multiple of \(v^{4/3}\), their theorem supplies a constant multiple of \((E/v)^6\) simple six-cycles; below that threshold the first term suffices. We use the [simple-cycle theorem](sources/morris-saxton-1309.2927v3.pdf).

Throughout, \(n^\epsilon\) can absorb any fixed number of divisor bounds, fixed-order factorization counts, unit-root counts and powers of \(\log(2n)\). All varying integers to which these bounds are applied have size at most a fixed power of \(n\). The finitely many smallest targets are absorbed in the constants. Nonempty positive equations have \(\kappa_i<n\).

We also record two elementary facts. A unit modulo \(m\) has \(O_j(C_j^{\omega(m)})=O_{j,\epsilon}(m^\epsilon)\) \(j\)-th roots for fixed \(j\). For odd prime powers this follows from the cyclic unit group; at powers of two its kernel has size at most \(2j\). More specifically the cube-map kernel has at most \(3^{\omega(m)}\) elements, also when the modulus is \(m^3\). At powers of two cubing is bijective. At powers of three, writing a cube root of one as \(z=1+3w\) gives
\[
                v_3(z^3-1)=v_3(z-1)+1,
\]
and at most three classes at every exponent. CRT completes the assertion, with the modulus-one case interpreted as one class.

Second, one unit residue class modulo \(n\) contains \(O(1+H^2/n)\) positive reduced fractions \(u/v\) with \(u,v\le H\). Distinct primitive positive vectors in that class have nonzero determinant divisible by \(n\). Their lengths are at most \(\sqrt2 H\), so their angular separation is at least \(n/(2H^2)\), unless two such vectors cannot exist. Ordering their angles in the first quadrant proves the bound. These elementary facts are also used in the [mixed-power graph argument](ERDOS_943_DUAL_GRAPHS_MIXED_POWER_PRODUCTS.md).

## 2. Actual boxes and five elementary bounds

Fix (2) and dyadic intervals for its six coordinates. Write their lower endpoints as
\[
 (A_0,A_1,A_2;C_0,C_1,C_2)
                  =(n^a,n^b,n^c;n^d,n^e,n^f).
\]
All exponents are nonnegative, and a nonempty box satisfies
\[
                 3a+4b+5c\le1,\qquad 3d+4e+5f\le1.                \tag{5}
\]
The positive coefficients only shrink this region. There are \(O((\log(2n))^6)\) boxes. Let \(\mathcal R\) count the primitive representations in one box, and put
\[
 p=b+c,\quad q=e+f,\quad r=a+c,\quad s=d+f,\quad A=a+d.
\]
Counting a complete left or right canonical tuple determines the other full summand, so
\[
                       \mathcal R\ll \min(n^{a+p},n^{d+q}).      \tag{6}
\]
Fixing the four kernel coordinates leaves a primitive positive binary cubic. Thus also
\[
                              \mathcal R\ll_\epsilon n^{p+q+\epsilon}.
                                                                        \tag{7}
\]
If \(A>1/2\), (5) gives
\[
                       p+q\le(2-3A)/4<1/8.
\]
This already proves a stronger bound than required for that box. Henceforth assume \(A\le1/2\).

For completeness, the omitted-coordinate graph argument is valid for the actual canonical factors in (3), without replacing the fourth/fifth kernel by an arbitrary powerfree variable. Omit exponent \(j=3\) or \(4\); a vertex consists of the other two coordinates. A kernel pair is joined when a representation in the box exists. Keep one actual omitted-root pair on each edge. The edge multiplicity in the full representation count is \(n^\epsilon\), by the corresponding primitive binary-form bound of Section 1.

Fix two distinct left vertices with products \(P_1,P_2\), including their common coefficient, and a common right vertex of product \(Q\). The selected equations have the form
\[
                         P_i x_i^j+Qy_i^j=n\quad(i=1,2).
\]
Write \(x_i=g\alpha_i,\ y_i=h\gamma_i\) using the respective pair gcds. Full-summand coprimality gives all required within-edge cross-coprimalities. Therefore
\[
             u/v=\alpha_1\gamma_2/(\alpha_2\gamma_1)
\]
is reduced, has unit denominator modulo \(n\), and has height at most four times the product of the two omitted-coordinate lower scales. Elimination gives
\[
 (u/v)^j\equiv P_2/P_1\pmod n,\qquad
       g^j(P_1u^j-P_2v^j)=n(\gamma_2^j-\gamma_1^j).             \tag{8}
\]
At omitted-coordinate product at most \(\sqrt n\), the elementary unit-root and fraction bounds give only \(n^\epsilon\) fractions, including the constant-width dyadic boundary. Each fraction has divisor-many allocations to its four reduced factors.

The bracket in (8) cannot vanish in an actual collision. Its vanishing forces \(\gamma_1=\gamma_2=1\), so the two omitted right coordinates coincide. The complete right summands then coincide, hence so do the complete left summands. After dividing by the fixed coefficient, uniqueness of (3) contradicts distinct left vertices. A nonzero bracket fixes \(g\) at most once. The first edge then determines \(Qh^j\). Divide by the known right coefficient and use (3): replacing the omitted coordinate by its divisor \(h\) preserves the canonical conditions, so the complete common right vertex and \(h\) are recovered uniquely. Positivity, integrality and stored-edge restrictions only discard candidates. This proves codegree \(O_\epsilon(n^\epsilon)\) on both sides.

For \(j=3\), the omitted-coordinate product is \(A_0C_0\le\sqrt n\). For \(j=4\), it is \(A_1C_1\le\sqrt n\), since \(b,e\le1/4\). The two graph vertex exponents are respectively \(p,q\) and \(r,s\). Cauchy--Schwarz and pair counting give, for side sizes \(m_1,m_2\) and maximum codegree \(\Lambda\),
\[
 E\ll m_1+m_2+\sqrt\Lambda\min(m_1\sqrt{m_2},m_2\sqrt{m_1}).
\]
Indeed \(E^2/m_2-E\le\Lambda m_1^2\) and its symmetric inequality suffice. Restoring the edge multiplicities proves
\[
 \mathcal R\ll_\epsilon n^{G_3+\epsilon},\quad
 G_3=\max\{p,q,\min(p+q/2,q+p/2)\},
\]
\[
 \mathcal R\ll_\epsilon n^{G_4+\epsilon},\quad
 G_4=\max\{r,s,\min(r+s/2,s+r/2)\}.                            \tag{9}
\]
This is the actual omitted-coordinate construction of [Powerful sum graphs](POWERFUL_SUM_GRAPHS.md); the details above specify the needed coefficient-uniform canonical scope.

Finally fix \(a_2,c_2\), and then \(a_0,c_1\). The remaining equation is \(PX^4+QY^3=n\), with positive \(P,Q\). It is absolutely irreducible: a root \(\xi\) of \(n-PX^4\) is nonzero and simple, so Eisenstein at \(X-\xi\), over \(\mathbb C[X]\), applies to \(Y^3-(n-PX^4)/Q\). Bombieri--Pila in a containing square, and the other orientation obtained by fixing \(a_1,c_0\), give
\[
 \mathcal R\ll_\epsilon n^{J+\epsilon},\quad
 J=c+f+\min\{a+e+\max(b,d)/4,\ b+d+\max(a,e)/4\}.              \tag{10}
\]
Thus, when \(A\le1/2\), the minimum of \(L=a+p,L'=d+q,G_3,G_4,J\) is a valid box exponent.

## 3. Noncoherent cycles: a uniform product count

Use the graph omitting the cube coordinates in one fixed box. Its left and right kernel products are
\[
             K_i=b_i^4d_i^5,\qquad L_j=s_j^4t_j^5,
\]
with squarefree coprime coordinate pairs. Each edge stores one actual solution
\[
                  \kappa_1K_i a_{ij}^3+\kappa_2L_j c_{ij}^3=n. \tag{11}
\]
A deterministic choice, such as the lexicographically first root pair in the box, suffices. All arguments below hold uniformly for any such selection. Let \(E\) be the number of edges. The full representation count is at most \(n^\epsilon E\). Write \(A_0=B,\ C_0=D\) in this section and the next two, so \(B\le a_{ij}<2B,\ D\le c_{ij}<2D\).

Order a simple six-cycle by its edge positions \(11,22,33,12,23,31\). Put
\[
 U=a_{11}a_{22}a_{33}c_{12}c_{23}c_{31},\qquad
 V=a_{12}a_{23}a_{31}c_{11}c_{22}c_{33}.                    \tag{12}
\]
They are units modulo \(n\), and multiplying (11) modulo \(n\) gives
\(U^3\equiv V^3\pmod n\). If \(R_1,R_2\) are the actual clipped upper root caps, then \(U,V\le R=(R_1R_2)^3\ll(BD)^3\). Call the cycle coherent if \(U=V\).

We first count the noncoherent sector. For all \(R\ge1\),
\[
 \#\{U,V\le R:(UV,n)=1,\ U^3\equiv V^3\pmod n,\ U\ne V\}
       \ll_\epsilon n^\epsilon(Rn^{-1/3}+R^2/n).              \tag{13}
\]
Reduce \(U=gu,V=gv\), with \((u,v)=1\). The multiplier \(g\) is a unit modulo \(n\), so \(u^3\equiv v^3\pmod n\). Since \(u\ne v\), the nonzero integer \(u^3-v^3\) is divisible by \(n\), giving \(\max(u,v)\ge n^{1/3}\). In a dyadic height range \(h<\max(u,v)\le2h\), the cube-root residue count and angular packing give
\(O_\epsilon(n^\epsilon(1+h^2/n))\) reduced pairs. Restoring \(g\le R/h\) and summing the geometric height ranges \(n^{1/3}\ll h\le R\) proves (13). A shift of the first dyadic endpoint affects only an absolute constant. If the range is empty the count is zero.

For fixed \(U,V\), there are at most \(d_6(U)d_6(V)\ll_\epsilon n^\epsilon\) allocations to the twelve ordered roots. These roots determine the six kernel products at most once when \(U\ne V\): the six equations (11), as a linear system in \(K_1,K_2,K_3,L_1,L_2,L_3\), have determinant
\[
                    \pm\kappa_1^3\kappa_2^3(U^3-V^3)\ne0.
\]
A valid kernel has only prime exponents four or five, and its squarefree coprime fourth/fifth factors are unique. Reject nonpositive, nonintegral or noncanonical answers and repeated vertices. Consequently
\[
              C_{6,\mathrm{noncoh}}\ll_\epsilon
                 n^\epsilon\{(BD)^3n^{-1/3}+(BD)^6/n\}.        \tag{14}
\]
Ordered choices of the cycle orientation change only an absolute factor.

## 4. Coherent potentials and the original divisor walk

Assume \(U=V\). There are unique positive integer labels
\(P_1,P_2,P_3,Q_1,Q_2,Q_3\) with global gcd one such that
\[
 \frac{a_{ij}}{c_{ij}}=\frac{P_i}{Q_j},\quad
 g_{ij}=\gcd(P_i,Q_j),\quad
             a_{ij}=P_i/g_{ij},\quad c_{ij}=Q_j/g_{ij}.        \tag{15}
\]
Propagating the ratios along five tree edges gives rational labels, and (12) supplies the closing edge. The graph is connected, so all such labels differ by one rational scale. Primitive integer normalization gives uniqueness. Full-summand coprimality implies \((a_{ij},c_{ij})=1\), which justifies the reduced identities in (15).

Their sizes have the sharper bounds
\[
                       P_i\le R_1^3R_2^2,\qquad
                       Q_j\le R_1^2R_2^3.                    \tag{16}
\]
For an explicit clearing of denominators on the tree \(11,31,33,23,22\), take
\(P_1=a_{11}c_{31}a_{33}c_{23}a_{22}\).
The successive right/left labels are
\[
\begin{split}
 Q_1&=c_{11}c_{31}a_{33}c_{23}a_{22},&
 P_3&=c_{11}a_{33}c_{23}a_{22}a_{31},\\
 Q_3&=c_{11}c_{23}a_{22}a_{31}c_{33},&
 P_2&=c_{11}a_{22}a_{31}c_{33}a_{23},\\
 Q_2&=c_{11}a_{31}c_{33}a_{23}c_{22}.
\end{split}
\]
Every left label has three \(a\)- and two \(c\)-factors, every right label two \(a\)- and three \(c\)-factors. Dividing by the global gcd only decreases them. All labels are units modulo \(n\), and all are at most \(n^5\).

Define integer weighted potentials
\[
                X_i=\kappa_1K_iP_i^3,\qquad
                Y_j=\kappa_2L_jQ_j^3.
\]
On every edge,
\[
                  X_i+Y_j=ng_{ij}^3,\qquad
                  \gcd(X_i,Y_j)=g_{ij}^3.                    \tag{17}
\]
This follows by multiplying the two actual coprime summands in (11) by \(g_{ij}^3\). The potentials are at most \(n^{16}\).

A selected anchor edge and a prescribed value of its intrinsic gcd determine only \(n^\epsilon\) coherent cycles. Index the anchor as \(11\), with prescribed gcd \(s\). It fixes
\(P_1=s a_{11},Q_1=s c_{11}\), and hence \(X_1,Y_1\). Traverse \(31,33,23,22\). First choose \(t=g_{31}\mid Q_1\). Equation (17) fixes
\[
                             X_3=nt^3-Y_1.
\]
Reject nonpositive answers or failures of divisibility by \(\kappa_1\). Decompose \(X_3/\kappa_1\) by (3); it uniquely recovers \(P_3\) as the cube coordinate and \(K_3\) as its fourth/fifth kernel. Next choose \(g_{33}\mid P_3\), which fixes \(Y_3=ng_{33}^3-X_3\); divide by \(\kappa_2\) and canonically recover \(Q_3,L_3\). Repeat using \(g_{23}\mid Q_3\) and \(g_{22}\mid P_2\) to recover \(P_2,K_2,Q_2,L_2\).

After each step impose (16) on the new label before using it as a divisor dividend. Check that the chosen divisor equals the actual endpoint gcd, that the recovered roots and kernels meet all caps and canonical restrictions, and that each edge agrees with the stored actual root pair and full-summand primitivity. Check edge \(12\) at the end, together with distinct vertices and global label gcd one. All actual cycles survive; all other choices may be discarded. Four divisor choices from integers bounded by \(n^5\) cost \(n^\epsilon\). The initial prescribed \(s\) is not silently counted: its anchor labels are checked at once, and the parameter \(s\) is summed explicitly below.

Writing \(g_*=\min_{\text{cycle edges}}g_{ij}\), orient at an edge attaining the minimum. There are only an absolute number of orientations. For every real \(Y>0\),
\[
                         \#\{C:g_*(C)\le Y\}
                            \ll_\epsilon E Y n^\epsilon.    \tag{18}
\]
When \(Y<1\) the sector is empty. This is an actual count of varying anchors and their gcds, not merely a completion bound for fixed uncounted data.

## 5. Complementary labels and cubic-root recovery

Let \(M\) be the lcm of the six primitive labels. The complementary labels
\[
               \bar P_i=M/P_i,\qquad \bar Q_j=M/Q_j
\]
are again primitive, with lcm \(M\), and complementation is an involution. At a prime, let \(t_v\) be its six label valuations, whose minimum is zero. The product identity in (12) and (15) gives
\[
 v_p(U)=\tfrac12\sum_{vw\ \mathrm{on\ the\ cycle}}|t_v-t_w|
                      \ge\max_v t_v .
\]
Indeed the two arcs from a minimum to a maximum each have variation at least the maximum. Thus
\[
                              M\mid U,\qquad M\le(R_1R_2)^3. \tag{19}
\]
Put \(\bar g_{ij}=\gcd(\bar P_i,\bar Q_j)\). Direct calculation gives
\[
 a_{ij}=\bar Q_j/\bar g_{ij},\quad
 c_{ij}=\bar P_i/\bar g_{ij},\quad
 g_{ij}\bar g_{ij}=\frac{M}{a_{ij}c_{ij}}
                      \le D_0:=\frac{(R_1R_2)^3}{BD}
                      \le64(BD)^2.                          \tag{20}
\]
The tree argument with exchanged root roles also gives
\[
               \bar P_i\le R_2^3R_1^2,\qquad
               \bar Q_j\le R_2^2R_1^3,                       \tag{21}
\]
so these divisor dividends too are at most \(n^5\).

Fix an anchor \(11\) and its complementary gcd \(s\). Now
\(\bar P_1=s c_{11},\bar Q_1=s a_{11}\), and both kernels are known. In the same four-step order, choose \(t=\bar g_{31}\mid\bar Q_1\); this fixes the original root \(a_{31}=\bar Q_1/t\). The original full-summand primitivity supplies the crucial units
\[
 (a_{31},\kappa_2L_1n)=1,\qquad
       c_{31}^3\equiv n(\kappa_2L_1)^{-1}\pmod{a_{31}^3}.     \tag{22}
\]
These units follow from the full two summands, not merely from coprimality of the two roots. If the fixed \(a_{31}\) is in its box, (22) gives at most
\[
                  O_\epsilon(n^\epsilon(1+D/B^3))
\]
candidates for \(c_{31}\) in its interval. There are at most \(3^{\omega(a_{31})}\) residue classes, and each has at most \(1+2D/a_{31}^3\) representatives below \(2D\). This also covers root one: modulo one there is one class and its complete interval length is retained.

For each candidate compute exactly
\[
              K_3=\frac{n-\kappa_2L_1c_{31}^3}
                              {\kappa_1a_{31}^3},\qquad
              \bar P_3=t c_{31}.                            \tag{23}
\]
The kernel must be a positive integer with prime exponents exactly four or five (or be one); its squarefree fourth/fifth factors are unique. In particular this step does not allow a free cube multiplier in the kernel. It recovers the missing original root before testing the integer kernel.

At edge \(33\), choose \(\bar g_{33}\mid\bar P_3\). This fixes \(c_{33}\); the symmetric unit congruence modulo \(c_{33}^3\) gives at most
\(n^\epsilon(1+B/D^3)\) candidates for \(a_{33}\). Compute \(L_3\) exactly and set \(\bar Q_3=\bar g_{33}a_{33}\). Edges \(23,22\) repeat the two types of step, recovering \(K_2,\bar P_2,L_2,\bar Q_2\). Apply the label bounds (21) before every subsequent divisor choice. Impose positivity, coefficient divisibility, the exact kernel prime-exponent test, actual gcd identities, all caps, stored-edge agreement, the closing edge and global primitive normalization. Thus a fixed anchor and prescribed complementary gcd have at most
\[
 n^\epsilon\mathcal L,\qquad
             \mathcal L=(1+D/B^3)^2(1+B/D^3)^2               \tag{24}
\]
completions. No uncounted root, kernel or scaling loop remains.

It follows, for every \(Y>0\), that
\[
                 \#\{C:\bar g_*(C)\le Y\}
                      \ll_\epsilon E Y\mathcal L n^\epsilon. \tag{25}
\]
As before, a threshold below one gives an empty integer-gcd sector.

If \(g_*>Y\), (20) puts every complementary gcd at most \(D_0/Y\). Equations (18) and (25) therefore cover all coherent cycles. Choose \(Y=\sqrt{D_0\mathcal L}\) to obtain
\[
 C_{6,\mathrm{coh}}\ll_\epsilon E n^\epsilon
           BD(1+D/B^3)(1+B/D^3).                            \tag{26}
\]
If \(D_0/Y<1\), that residual sector is empty and the same upper bound holds. Formula (26) is uniform in every root rectangle and in the coefficients. It needs no short-interval condition such as \(2D<B^3\), and no prior bound for \(E\).

## 6. The cycle self-bound in exponent form

Return to the six exponent coordinates and put
\[
        V_0=\max(p,q),\quad A=a+d,\quad
        w=A+(d-3a)^++(a-3d)^+
                        =\max\{A,2a-2d,2d-2a\}.
\]
The graph has \(O(n^{V_0})\) vertices. Formula (26) bounds its coherent cycles by \(E n^{w+\epsilon}\), whereas (14) bounds its noncoherent cycles by
\[
                  n^\epsilon\{n^{3A-1/3}+n^{6A-1}\}.
\]
Taking sixth roots in (4) gives
\[
 E\ll n^{4V_0/3}
       +n^{V_0+w/6+\epsilon}E^{1/6}
       +n^{V_0+A/2-1/18+\epsilon}
       +n^{V_0+A-1/6+\epsilon}.                             \tag{27}
\]
If \(x\le C(B_0+D_1x^{1/6})\), either \(x\le2CB_0\) or
\(x^{5/6}\le2CD_1\). Therefore \(x\ll_C B_0+D_1^{6/5}\). Applying this elementary absorption to (27), and restoring the \(n^\epsilon\) multiplicity on each selected edge, proves
\[
 \mathcal R\ll_\epsilon n^{F_{\mathrm{cyc}}+\epsilon},\qquad
 F_{\mathrm{cyc}}=
 \max\left\{\frac{4V_0}{3},\frac{6V_0+w}{5},
                  V_0+\frac A2-\frac1{18},V_0+A-\frac16\right\}. \tag{28}
\]
This is a solved self-bound, rather than a use of the desired representation estimate as an input.

On the closed cone \(a\le3d,d\le3a\), one has \(w=A\). There we use (28) together with (6), (9), (10). Outside this cone, the old five bounds alone will suffice. This weaker use of the available bounds keeps the exact optimization below short and explicit.

## 7. A rational proof of the exponent \(33/134\)

Assume (5) and \(A\le1/2\). Set
\[
 F_*=\min\{L,L',G_3,G_4,J\}.
\]
On the cone \(a\le3d,d\le3a\), put
\[
 C_*=\max\{4V_0/3,(6V_0+A)/5,
                     V_0+A/2-1/18,V_0+A-1/6\}
\]
and \(F=\min(F_*,C_*)\); outside the cone put \(F=F_*\).
We prove \(F\le33/134\) on the full stated real domain, without a symmetry assumption or a finite grid. The actual box exponent from all the proved bounds is no larger than \(F\).

Fix any \(0\le t\le F\). Define
\[
\begin{split}
 R&=1-3a-4b-5c,&S&=1-3d-4e-5f,\\
 H&=a+p-t,&H'&=d+q-t,\\
 M_1&=2p+q-2t,&M_2&=p+2q-2t,\\
 N_1&=2r+s-2t,&N_2&=r+2s-2t,\\
 I_p&=p-t,&I_r&=r-t,\qquad I_s=s-t .
\end{split}                                                        \tag{29}
\]
The always nonnegative slacks are \(R,S,H,H'\). The other slacks are used only in their specified branches. Also put
\[
\begin{split}
 J_b&=c+f+a+e+b/4-t,&J_d&=c+f+a+e+d/4-t,\\
 J_a&=c+f+b+d+a/4-t,&J_e&=c+f+b+d+e/4-t .
\end{split}
\]
Since \(J\ge t\), the first two are nonnegative respectively when \(b\ge d\) or \(d\ge b\), and the last two respectively when \(a\ge e\) or \(e\ge a\). An unselected slack need not be nonnegative.

Every identity below has a mirror under exchanging
\((a,b,c)\) and \((d,e,f)\). This swaps \(p,q\), \(r,s\), \(R,S\), \(H,H'\), \(M_1,M_2\), \(N_1,N_2\), \(J_b,J_e\), and \(J_d,J_a\). All resources and bounds are preserved.

### Individual terms of the cubic graph

Suppose \(p\ge t\). The following four exhaustive branches of \(G_4\ge t\) give nonnegative right-hand sides:
\[
\begin{array}{ll}
 r\ge t:&
 25-104t=25R+71H+21c+29I_p+4I_r,\\[1mm]
 s\ge t,\ b\ge d:&
 100-416t=58R+42S+90H+84H'+5b+116I_p+42I_s+84J_b,\\[1mm]
 s\ge t,\ d\ge b:&
 200-832t=137R+63S+243H+84H'+31b+274I_p+63I_s+168J_d,\\[1mm]
 N_1,N_2\ge0:&
 75-312t=58R+17S+116H+22H'+46e+34f+116I_p+29N_1 .
\end{array}                                                        \tag{30}
\]
If neither individual term \(r,s\) is at least \(t\), both mixed inequalities are required by \(G_4\ge t\). Thus these cases exhaust the domain, and each gives \(t\le25/104<33/134\). The mirror handles \(q\ge t\). Otherwise \(p,q<t\), so \(G_3\ge t\) forces \(M_1,M_2\ge0\).

### Individual terms of the quartic graph

If \(r\ge t\), choose the appropriate maximum in the second part of \(J\):
\[
\begin{array}{ll}
 a\ge e:&
 99-405t=67R+32S+126H+64H'+52M_1+6M_2+67I_r+32J_a,\\[1mm]
 e\ge a:&
 11-45t=8R+3S+16H+5H'+6M_1+8I_r+4J_e .
\end{array}                                                        \tag{31}
\]
Both give \(t\le11/45<33/134\). The mirror handles \(s\ge t\). In the remaining case, \(G_4\ge t\) forces \(N_1,N_2\ge0\). Both graphs are now in their mixed branches.

### The region outside the cubic root-height cone

If \(a\ge3d\), the identity
\[
 61-256t=15R+46S+158H'+26M_1+8b+23N_2+22(a-3d)              \tag{32}
\]
gives \(t\le61/256<33/134\). The mirror handles \(d\ge3a\), including the boundaries. Thus the only remaining region is inside the cone, where \(w=A\).

### The four cycle terms

By exchanging the two sides assume \(q\ge p\), so \(V_0=q\). Since \(C_*\ge t\), at least one of
\[
 T_1=4q/3-t,\quad T_2=A+6q-5t,\quad
 T_3=q+A/2-1/18-t,\quad T_4=q+A-1/6-t
\]
is nonnegative. Direct expansion gives
\[
\begin{split}
 12-49t&=4R+8S+8H+16H'+4M_1+4N_2+9T_1,\\
 33-134t&=11R+22S+16H+38H'+14M_1+11N_2+6T_2,\\
 25-108t&=9R+18S+18H'+18M_1+9N_2+36T_3,\\
 5-22t&=3R+3S+6M_1+N_1+N_2+6T_4 .
\end{split}                                                        \tag{33}
\]
Use the identity corresponding to a nonnegative \(T_i\). Every other term in that identity is already nonnegative. The resulting upper bounds are \(12/49,33/134,25/108,5/22\), each at most \(33/134\). The mirrors cover \(p\ge q\).

Equations (30)--(33) exhaust all cases and prove \(F\le33/134\). These are polynomial identities on all seven real variables; their use as inequalities depends only on the explicitly stated nonnegative slacks.

For this specified sub-menu, the constant is attained at
\[
                  (a,b,c;d,e,f)=(12,7,14;9,13,11)/134.
\]
Both resources are one, and \(p=21/134,q=24/134,r=26/134,s=20/134,A=21/134\). Then \(L=L'=G_3=G_4=33/134\), while \(J=177/536>33/134\). The four cycle terms are \(32/134,33/134,487/2412,34/201\), so their maximum is \(33/134\). The point is strictly inside the cone. This establishes the exact maximum of the particular bounds used here; it does not assert optimality after adding every other omitted-coordinate or mixed-curve bound.

## 8. Summing boxes and removing common divisors

The high-\(A\) boxes were covered by (7). All other boxes have exponent at most \(33/134\) by Section 7. Summing the logarithmically many boxes and allocating the subpower losses proves the uniform primitive coefficient assertion of Theorem 1. The constants in the graph reconstruction, unit congruences, binary forms, curves, cycle theorem and divisor estimates are independent of \(\kappa_1,\kappa_2\). No coefficient pairs were summed in proving that assertion.

For arbitrary positive 3-full summands \(m_1+m_2=N\), let \(h=(m_1,m_2)\). The integer \(h\) is 3-full, since every positive minimum of two exponents at least three is at least three. Put \(n=N/h\). The quotients \(m_1/h,m_2/h\) are coprime. From each quotient extract the product \(\kappa_i\) of its prime powers of exponent one or two. Such a prime must divide \(h\); otherwise its exponent in the original 3-full summand would still be one or two. Since \(v_p(h)\ge3\), every extracted power divides \(h\). Thus
\[
        \kappa_i\mid h,\qquad m_i/h=\kappa_i u_i,
\]
where \(u_i\) is 3-full. Canonically present \(u_i\) by (3). This gives exactly an equation (2), with coprime full summands. For fixed \(h,\kappa_1,\kappa_2\), the original \(m_i=h\kappa_i u_i\) are uniquely recovered.

There are at most \(\tau(N)^3\) choices for \(h,\kappa_1,\kappa_2\). Apply the coefficient-uniform primitive estimate at \(n=N/h\le N\), sum these choices and absorb the divisor factor in \(N^\epsilon\). The case \(n<2\) has no positive two-summand representation. This proves (1) for every positive \(N\), including noncoprime original summands.

## 9. Relation to the existing graph bounds

The canonical factorization and omitted-coordinate graphs are those used in [the general \(k\)-full sum bounds](POWERFUL_SUM_GRAPHS.md), whose \(k=3\) specialization gives exponent \(1/4\). The mixed-curve step is the degree-four instance of the coefficient-uniform curve argument in [the mixed-power product theorem](ERDOS_943_DUAL_GRAPHS_MIXED_POWER_PRODUCTS.md). Their canonical sets and actual edge conditions are kept distinct here: a general 3-full integer may contain prime exponent five and is not necessarily a product of a cube and a fourth power.

The present proof adds the actual coherent-cycle reconstruction and complementary unit-cubic recovery, then the explicit exponent calculation. It does not use a fixed-shift powerful-pair theorem, an average over the target, the former quarter-bound principal, a conjectured cycle frequency, or an unproved simultaneous-values estimate as a premise. The result is a pointwise theorem for the fixed positive sum, with its nonprimitive normalization supplied in Section 8.
