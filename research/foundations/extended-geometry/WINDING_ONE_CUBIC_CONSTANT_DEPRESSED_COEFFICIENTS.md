# Constant depressed coefficients in the quartic-leading septic family

Work over a field of characteristic zero on \(t^2+2xy=1\), with the
normalized sphere functional \(\Lambda\). Write \(U=1+x\) and \(v=Ut\).
A regular winding-one phase in the normalized quadruple-root family is
\[
 f=\frac Ux P(x,v),\qquad P(x,v)=v^3+d(x)v^2+b(x)v+c(x).
\]
Put \(d_0=d(0)\). Regularity gives
\(P(0,v)=P_0(v)=(v^2-1)(v+d_0)\). The depressed coefficients are
\[
 p=b-d^2/3,\qquad q=c-bd/3+2d^3/27.
\]

**Theorem ConstantDepressedCoefficientExclusion.** If \(f\) has least
sphere degree seven and one of \(p,q\) is constant in \(x\), then the
positive pure moments of \(f\) do not vanish eventually.

The same statement holds before normalization when
\(E=e(1+\lambda x)^4\), \(e\lambda\ne0\), and the three
[quadruple-root boundary divisibilities](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#a-quadruple-root)
hold. Scaling \(x\) by \(\lambda\), scaling \(y\) inversely, and
multiplying the phase by a nonzero scalar give the displayed form.

## Reduction to one polynomial of the base variable

The [winding-one degree formula](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md)
gives \(\deg d\le3\), \(\deg b\le5\), and \(\deg c\le7\).
If \(\deg d=3\), then \(p\) has degree six and \(q\) has degree nine.
Their highest powers come from \(-d^2/3\) and \(2d^3/27\), respectively.
Thus either constant-coefficient hypothesis implies \(\deg d\le2\),
and consequently \(\deg p\le5\) and \(\deg q\le7\).

**Lemma DepressedTranslationMoments.** For any normalized regular
cubic above, replacing \(P(x,v)\) by \(P(x,v+(d_0-d(x))/3)\)
preserves every positive sphere moment and makes the quadratic
coefficient constant \(d_0\).

Let \(\delta=(d_0-d(x))/3\) and replace \(P(x,v)\) by
\(P(x,v+\delta)\). The fixed-endpoint identity for the sphere moments is
\[
 \Lambda(f^m)=[x^m]U^{m-1}\frac12\int_{-1}^1P(x,v)^m\,dv.
 \tag{1}
\]
In the change from \(t\) to \(v\), the initial endpoints are \(\pm U\).
The discarded strips from \(\pm1\) to \(\pm U\) contribute a
multiple of \(x^{m+1}\), since \(P(0,\pm1)=0\).
The translation by \(\delta\in xK[x]\) introduces strips with
the same divisibility. This proves that it preserves every positive
moment. Both strip arguments are polynomial identities, using a
parametrization of each strip by \(v=\pm1+x z\) or
\(v=\pm1+\delta z\) with \(0\le z\le1\).

The translated cubic is \((v+d_0/3)^3+p(x)(v+d_0/3)+q(x)\).
Thus if \(p\) is constant, its difference from \(P_0\) is
\(H(x)=q(x)-q(0)\), of degree at most seven. If \(q\) is constant,
the difference is \(H(x)(v+d_0/3)\), where \(H=p-p(0)\) has degree
at most five.

Set \(X=x/(1+x)\) and \(h(X)=H(X/(1-X))\). The coefficient
substitution in (1) uses \(U^{m-1}dx/x^{m+1}=dX/X^{m+1}\).
Write \(\langle A\rangle=\tfrac12\int_{-1}^1A(v)\,dv\). For
\(R=1\) or \(R=v+d_0/3\), the resulting moments are
\[
 M_m=\sum_{j\ge1}\binom mj\langle P_0^{m-j}R^j\rangle
                         [X^m]h(X)^j.                       \tag{2}
\]
The sum ends at \(m\). If \(h\) has first term \(a X^k\), then for
\(k\le m<2k\) only the linear term contributes. If \(N+1<2k\)
and the linear coefficient is nonzero for every \(k<m\le N+1\),
vanishing of these moments is
impossible for \(\deg H\le N\). Indeed it would give
\(h(X)=aX^k\bmod X^{N+2}\). The coefficient of \(x^{N+1}\) in
\(a(x/(1+x))^k\) is
\(a(-1)^{N+1-k}\binom{N}{k-1}\ne0\).

## The direction \(R=1\)

Here vanishing of the first twelve moments suffices. Suppose \(H\ne0\)
and write \(k=\operatorname{ord}_x H\), \(a=[x^k]H\), and
\(z=d_0^2\). The first possible moment forces
\(\langle P_0^{k-1}\rangle=0\). In particular \(k\ge2\).

If \(d_0=0\), the odd powers of \(P_0\) integrate to zero and
its even powers have nonzero rational integrals. This forces \(k\)
even. At \(m=2k\) the linear term still integrates to zero and
the quadratic term is \(\binom{2k}{2}a^2\langle P_0^{2k-2}\rangle\),
which is nonzero. Here \(2k\le12\).

For \(d_0\ne0\), up to nonzero rational factors the linear
coefficients \(\langle P_0^n\rangle\) are as follows.
\[
\begin{array}{c|l}
n&\langle P_0^n\rangle\text{ up to scale}\\\hline
1&d_0\\
2&z+1/7\\
3&d_0(z+1/3)\\
4&z^2+6z/11+3/143\\
5&d_0(z^2+10z/13+1/13)\\
6&z^3+z^2+3z/17+1/323\\
7&d_0(z^3+21z^2/17+105z/323+5/323)
\end{array}                                                   \tag{3}
\]
Thus \(k=2\) is impossible. For \(k=5,6,7\), the four monic
polynomials in the last four rows are pairwise coprime. Their
resultants modulo \(29\), taken in lexicographic order of pairs of
rows, are \(11,3,11,17,2,17\). The denominators are units and the
polynomials are monic. Hence every coefficient required for
\(m=k+1,\ldots,8\) is nonzero, and (2) gives the truncation
contradiction with \(N=7\).

For the other two values of \(k\), the initial equations and
the final tests reduce to the following table. It lists the
coefficients \(([x^k]H,\ldots,[x^7]H)\) and gives the last moments
up to nonzero rational factors.
\[
\begin{array}{c|c|l|l}
k&z&([x^k]H,\ldots,[x^7]H)&\text{tests}\\\hline
3&-1/7&(a,-3a,6a,-10a-56d_0a^2,15a+336d_0a^2)&
 M_8\sim a(a-d_0/8),\quad M_{10}\sim a(a-5d_0/28)\\
4&-1/3&(a,-4a,10a,-20a)& M_9\sim a
\end{array}                                                    \tag{4}
\]
The equations used to obtain the coefficient tuples are
\(M_{k+1}=\cdots=M_7=0\), successively linear in the newest
coefficient. The tests in either row contradict \(a\ne0\), since
\(d_0\ne0\). This proves \(H=0\) in the first direction.

## The direction \(R=v+d_0/3\)

For \(d_0\ne0\), vanishing of the first seven moments suffices.
The linear coefficient \(\langle R\rangle=d_0/3\) shows that
the first index of a nonzero \(H\) must be \(k\ge2\). Again put
\(z=d_0^2\). The needed linear coefficients are
\[
\begin{array}{c|l}
n&\langle P_0^nR\rangle\text{ up to scale}\\\hline
1&z+3/5\\
2&d_0(z+1)\\
3&z^2+4z/3+1/11\\
4&d_0(z^2+18z/11+3/11)\\
5&z^3+25z^2/13+7z/13+3/221 .
\end{array}                                                   \tag{5}
\]
The three monic polynomials in the last three rows have pairwise
resultants \(11,10,10\) modulo \(29\). For \(k=4,5\), the
truncation argument applies with \(N=5\). For \(k=2,3\), the
corresponding substitutions from \(M_{k+1}=\cdots=M_5=0\) are
\[
\begin{array}{c|c|l|l}
k&z&([x^k]H,\ldots,[x^5]H)&\text{tests}\\\hline
2&-3/5&(a,-2a,3a+11a^2/6,-4a-22a^2/3)&
 M_6\sim a(221a^2+880a+240),\quad M_7\sim d_0 a(55a+18)\\
3&-1&(a,-3a,6a)& M_7\sim d_0a .
\end{array}                                                    \tag{6}
\]
In the first row \(a=-18/55\) is not a root of the quadratic.
The last row rules out any nonzero \(a\).

For direct arithmetic verification of (3)–(6), put
\[
 C_{n,i}=\sum_{r=0}^n\sum_{s=0}^n
  \frac{(-1)^{n-s}\binom nr\binom ns\,d_0^{n-r}
       \,\mathbf 1_{r+i\ {\rm even}}}{2s+r+i+1}.
\]
Then \(C_{n,i}=\langle P_0^n v^i\rangle\). In (2) use
\(C_{m-j,0}\) or
\(\sum_{i=0}^j\binom ji(d_0/3)^{j-i}C_{m-j,i}\),
respectively. If \(H=\sum a_i x^i\), use
\([X^r]h=\sum_{i\le r}\binom{r-1}{i-1}a_i\).
All substitutions in the two tables take place in the quadratic
field specified by \(d_0^2=z\). The resultant checks can be made
by reducing the displayed polynomials and their Sylvester
determinants in \(\mathbb F_{29}\).

## The septic conclusion

By the [logarithmic boundary corollary](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#consequences-in-degree-seven),
eventual purity would give purity at every positive order and
a fixed surviving marker \(x\). If \(q\) is constant and
\(d_0=0\), then \(q=0\) by regularity. The projective cubic
invariant is constant, including the zero pair in its definition.
The [constant-invariant safety theorem](CONSTANT_CUBIC_INVARIANT_ELLIPTIC_MODULI_BRIDGE.md)
would make \(f\) safe, contrary to that marker.

In the other cases (2)–(6) force both \(p,q\) to be constant.
Since \(\deg d\le2\), expansion of
\((v+d/3)^3+p(v+d/3)+q\) and the degree formula give
least sphere degree at most six. This contradicts the septic
hypothesis.
