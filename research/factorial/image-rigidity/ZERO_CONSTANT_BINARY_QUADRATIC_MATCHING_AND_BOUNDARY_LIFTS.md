# Matching and boundary lifts for zero-constant binary quadratics

Write
\[
 f=\frac{\alpha}{2}x^2+\beta xy+\frac{\gamma}{2}y^2+dx+ey,\qquad
 P_j=\frac{\Lambda_2(f^j)}{j!},\qquad E(T)=\sum_{j\ge0}P_jT^j .
\]
The natural coefficients are \(\alpha,\beta,\gamma,d,e\). Over a finite
extension of \(\mathbb Q_2\), let \(v(2)=1\), with \(v(0)=+\infty\),
and divide the phase by a natural coefficient of minimum valuation.
Then all natural coefficients are integral and at least one is a unit.

**Theorem.** A block of four zero positive-index moments can occur
only in one of these residue types, up to exchanging the axes.

1. Only \(\beta\) is a unit. The block then starts at an index at
   least two.
2. The units are \(\alpha,\beta\) and their residues agree. In this
   case, let \(r\) be the unique unit root of
   \(1+\alpha T+d^2T^2\). The four-zero block also requires
   \[
                   v\bigl(\beta-\alpha+d+e+d(e-d)r\bigr)\ge1 .
                   \tag{1}
   \]

Every binary quadratic in the zero-constant family with
\(\gamma=0,d=e\), or with \(\gamma=0,\beta=\alpha\), satisfies
its actual-support Strong Factorial window in characteristic zero.
The genuine three-term bilinear family \(\beta xy+dx+ey\) has
its three-window if \(v(\beta)\ge\min(v(d),v(e))\) at some
place over two.

The full five-parameter family necessarily has nonzero phases with
four consecutive zero moments at every start. Its actual support-five
target is a five-window.

## The matching reduction and the unit pole

For integral natural coefficients the [weighted matching analysis](BINARY_ZERO_CONSTANT_QUADRATIC_WEIGHTED_MATCHING_LOCAL.md)
proves the two residue types in arbitrary characteristic-two residue
fields. It gives the universal coefficientwise congruence
\[
 E(T)\equiv N(T)/D(T)\pmod{2\mathcal O[[T]]},\qquad
 N=1+(d+e+\beta)T+deT^2,\quad
 D=(1+\alpha T+d^2T^2)(1+\gamma T+e^2T^2).
 \tag{2}
\]
This is a congruence of universal integral coefficient polynomials
and holds at precision one even over ramified fields. The unit-factor
divisibility in the weighted analysis gives \(N(r)\equiv0\bmod
2\mathcal O\) in the second type, which is (1).

## The exact interior and boundary recurrence

Use the univariate boundary moments
\[
 A_j=\Lambda_1((\alpha x^2/2+dx)^j)/j!,\qquad
 B_j=\Lambda_1((\gamma y^2/2+ey)^j)/j!.
\]
They start at \(A_0=B_0=1\), \(A_1=\alpha+d,\ B_1=\gamma+e\).
For \(j\ge2\),
\[
 A_j=\alpha(2j-1)A_{j-1}+d^2A_{j-2},\qquad
 B_j=\gamma(2j-1)B_{j-1}+e^2B_{j-2}.
 \tag{3}
\]
This is the native single-variable recurrence in the matching
proof. Define
\[
\begin{split}
 \Delta&=\alpha\gamma-\beta^2,\qquad D_{\rm lin}=\alpha+\gamma-2\beta,\qquad
 C_{\rm lin}=\gamma d^2-2\beta de+\alpha e^2,\\
 T_j&=(\gamma-\beta) B_j+(\alpha-\beta) A_j
       +(\gamma d-\beta e)B_{j-1}+(\alpha e-\beta d)A_{j-1}.
\end{split}
\]
For every \(j\ge2\) the exact relation is
\[
        D_{\rm lin} P_j=2j\Delta P_{j-1}+C_{\rm lin} P_{j-2}+T_j. \tag{4}
\]
To check it put \(H=\begin{pmatrix}\alpha&\beta\\\beta&\gamma\end{pmatrix}\),
\(a=(d,e)^{\mathsf T}\), \({\bf1}=(1,1)^{\mathsf T}\), and
\(V_j=\Lambda_2((x,y)^{\mathsf T} f^j)/j!\). Monomial factorial
evaluation gives
\[
 H V_{j-1}={\bf1}P_j-aP_{j-1}-(B_j,A_j)^{\mathsf T},
 \qquad {\bf1}^{\mathsf T}V_{j-1}=2jP_{j-1}-a^{\mathsf T} V_{j-2}.
\]
The first identity uses the two boundary terms in integration by
parts. The second uses \(x f_x+y f_y=2f-dx-ey\) in the weighted
Euler identity. Multiply the second by \(\det H\) and use both
Stein identities after multiplying by \(\operatorname{adj}H\).
The terms in \({\bf1}^{\mathsf T}\operatorname{adj}(H)a P_{j-1}\)
cancel, yielding (4) even when \(\Delta=0\). Thus a block of
\(L\ge3\) zero moments forces \(L-2\) consecutive zero \(T_j\),
starting two positions after the moment block starts.

## The support with one diagonal absent

Take \(\gamma=0\). Since \(B_j=e^j\), the forcing becomes
\[
              T_j=(\alpha-\beta)A_j+(\alpha e-\beta d)A_{j-1}
                        -2\beta e^j.                         \tag{5}
\]

First suppose \(d=e\ne0\), and take \(d=e=1\) by scaling the
phase. If \(\alpha\beta\ne0\), the desired window has length
four. If \(\alpha=\beta\), (5) is \(-2\beta\ne0\).
Otherwise two consecutive zero \(T_j,T_{j+1}\) would give
\[
 0=T_{j+1}-T_j=(\alpha-\beta)\alpha(2j+1)A_j
\]
by (3). Thus \(A_j=0\). For \(v(\alpha)>0\), (3) gives
\(A_j\equiv1\) in the residue field. For \(v(\alpha)<0\),
the leading term of \(A_j\) as a polynomial in \(\alpha\)
has odd coefficient \((2j-1)!!\), and the remaining coefficients
are integral, so the leading term is uniquely least in valuation.
It follows that \(v(\alpha)=0\). Equation (3) run backwards modulo
the maximal ideal, with unit last coefficient, shows
\(A_{j-1}\) is a unit. Then
\[
                (\alpha-\beta)A_{j-1}=2\beta
\]
forces \(v(\beta)=0\): for \(v(\beta)>0\) the left side is a unit
and for \(v(\beta)<0\) it has valuation \(v(\beta)\).
With \(d=e=1,\ \alpha,\beta\) units, (2) already excludes four
zeros. If \(\alpha=0,\beta\ne0\), formula (5) is \(-4\beta\)
and gives the three-window. If \(\beta=0\), the phase is diagonal
and satisfies the [matching theorem](QUADRATIC_STRONG_FACTORIAL_WINDOWS_FROM_MATCHINGS.md).
When \(d=e=0\) use the
[homogeneous quadratic theorem](BINARY_HOMOGENEOUS_QUADRATIC_STRONG_FACTORIAL.md).

Next suppose \(\alpha=\beta\ne0\). If both linears are nonzero
then there are four terms. A four-zero block forces the second
matching residue type, so normalize \(\alpha=\beta=1\) and
\(v(d),v(e)>0\). Now every \(A_j\) is a unit by (3), while
\[
               T_j=(e-d)A_{j-1}-2 e^j .
\]
If both \(T_j,T_{j+1}\) are zero, then \(e-d\ne0\) and their
ratio gives \(A_j/A_{j-1}=e\), contrary to \(v(e)>0\).

If exactly one linear is present, the same plane has a
three-window by a direct calculation. For \(e=0\), (5)
would require \(A_{j-1}=0\) at the end of a three-zero run.
After normalizing \(\alpha=\beta=1\), the case \(v(d)<0\) gives
\(P_l/d^l\) with unit residue for every \(l\), since the only
unit natural coefficient is the linear coefficient in \(x\).
If \(v(d)>0\), then \(A_{j-1}\) is a unit. If \(v(d)=0\),
(2) has proper denominator of degree two, excluding even
two consecutive zeros. For \(d=0\) and \(e\ne0\),
\(A_{j-1}=(2j-3)!!\) is a unit after the same normalization
and (5) becomes \(e A_{j-1}-2e^j\). If \(v(e)\ge0\), this
has two different term valuations. If \(v(e)<0\), normalize
by \(e\) instead and (2) reduces to \(1/(1+T)\), whose
coefficients are all nonzero. When both linears vanish use
the homogeneous theorem again. The case \(\alpha=\beta=0\)
is diagonal.

## A valuation-one lift for the bilinear three-term phase

Assume \(\alpha=\gamma=0\) and \(\beta de\ne0\). Formula (4) gives
\[
                 P_j=j\beta P_{j-1}+deP_{j-2}+d^j+e^j .      \tag{6}
\]
Three consecutive zeros ending at \(k\ge3\) therefore require
\(d^k+e^k=0\), so \(v(d)=v(e)\). Under the theorem's local
hypothesis take \(d=1\), write \(q=e\) and \(b=\beta\), and put
\(s=1+q,\ H_0(T)=(1-T)(1-qT)\). Here \(b\) is integral and
\(q^k=-1\). Expanding the cross count \(i\) gives the exact
identity
\[
                    E(T)=\sum_{i\ge0} i!\,b^i T^i H_0(T)^{-i-1}.
                    \tag{7}
\]
For example the coefficient of \(b^i\) in \(P_l\) is the sum of
\(i!\binom{i+u}{i}\binom{i+w}{i}q^w\) over \(u+w=l-i\), which
is the claimed series coefficient.

Write \(h_l(1,q)=1+q+\cdots+q^l\). In characteristic two the
even moments in (7) reduce to \(h_l(1,\bar q)\): the \(i=1\)
term has only odd indices and every \(i\ge2\) vanishes.
If \(\bar q=1\), every even moment is a unit. Otherwise
\(\bar q^k=1\), and \(k\) must be odd since for even \(k\) we
would have \(\bar P_k=(1+\bar q)/(1+\bar q)=1\). At odd \(k\)
the reduction of \(P_k\) is
\[
               (\bar s+\bar b)\,h_{(k-1)/2}(1,\bar q)^2 .
\]
The second factor is nonzero because \(\bar q^{k+1}=\bar q\ne1\).
Thus \(b+s\) has positive valuation. Put \(m=k-1=2l\).
The even coefficient polynomial \(P_m(b)\) reduces independently
of \(b\) as a universal polynomial, so
\(P_m(b)-P_m(-s)\in 2(b+s)\mathcal O\).

It remains to evaluate \(P_m(-s)\) to valuation one. The terms
with \(i\ge4\) in (7) are zero modulo \(4\mathcal O\), and the
\(i=3\) term contributes zero modulo \(4\mathcal O\) at the
even index, since \(H_0^{-4}\) has only exponents divisible by
four modulo two. The first two terms at \(b=-s\) use
\[
 \frac{H_0-sT}{H_0^2}
     = (T/H_0)' + \frac{2T(-s+qT)}{H_0^2}.
\]
They and the \(i=2\) term give
\[
 \frac{P_m(-s)}2\equiv
 \frac1{\bar s}
    +\bar q\,h_{l-1}(1,\bar q)^2
    +\bar s^2[T^{2l-2}]\frac{1+\bar qT^2}{\overline H_0^4}
        \pmod{\mathfrak m_{\mathcal O}}.                     \tag{8}
\]
Indeed \([T^m](T/H_0)'=k h_m(1,q)=2k/(1-q)\), and the
even-index part of \(\overline H_0^{-3}\) is the displayed
numerator over \(\overline H_0^4\). The middle term in (8)
is \(1/\bar s\), using \(\bar q^{2l}=\bar q^{-1}\).
When \(l\) is even the last term is
\(\bar s^2\bar q\,h_{(l-2)/2}(1,\bar q)^4=1/\bar s\).
When \(l\) is odd it is
\(\bar s^2h_{(l-1)/2}(1,\bar q)^4=1/\bar s\).
These use \(\bar q^{2l}=\bar q^{-1}\) and
\(\bar q^{2l+2}=\bar q\), respectively.
The total residue is \(1/\bar s\ne0\), hence \(v(P_m(b))=1\),
a contradiction. The sector \(v(\beta)<v(d)=v(e)\) is open in
this three-window argument. Two adjacent zero forcing terms
in (6) are impossible for any nonzero \(d,e\), so the bilinear
family has an unrestricted four-window.

## Field scope and the four-window obstruction

For the unconditional assertions, a characteristic-zero
counterexample at a fixed start would give an algebraic
counterexample. Adjoin inverses for required nonzero coefficients
and take an algebraic point of the finitely generated system by
the weak Nullstellensatz. An embedding into
\(\overline{\mathbb Q}_2\) puts the resulting coefficients in a
finite local extension, so the reductions above apply.

For each fixed positive \(n\), the equations
\(P_n=\cdots=P_{n+3}=0\) are four homogeneous equations in
the five independent natural coefficients. The
[homogeneous common-zero lemma](https://stacks.math.columbia.edu/tag/03RB)
gives a nonzero solution over \(\overline{\mathbb Q}\).
Indeed their affine zero set has dimension
at least one by the height bound for an ideal generated by four
elements of a polynomial ring in five variables. Thus it
contains a point other than the origin. Equivalently, four
homogeneous hypersurfaces in \(\mathbb P^4\) always intersect.
This argument uses no computation and applies to every start.
It does not prescribe the support of a witness or resolve
the five-window for the general support.
