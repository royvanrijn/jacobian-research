# Prime-terminal Chebyshev squarefreeness and pure-cross three-windows

This proof excludes every noncentral resonance of odd prime terminal order. Combining it with the already proved canonical opposite-axis theorem also closes the pure-cross three-window at every odd prime terminal order. General composite terminal orders are not covered by the new argument.

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

This is not an all-start pure-cross theorem and does not close general binary quadratic SFC, SFC(2), or SFC(3). It is an infinite class of starting indices \(p-2\), with no coefficient restriction beyond the actual three-term hypothesis \(Bde\ne0\). Composite terminal orders remain separate.

