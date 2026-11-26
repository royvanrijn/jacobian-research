# Dyadic residue symmetry for the bilinear three-term phase

Let \(f=\beta XY+dX+eY\), where \(\beta,d,e\ne0\), and
\(P_n=\Lambda_2(f^n)/n!\), with
\(\Lambda_2(X^iY^j)=i!j!\). Work over a finite extension of
\(\mathbb Q_2\), with \(v(2)=1\).
Any three consecutive zero moments ending at \(k\) require
\(v(d)=v(e)\). Write \(s=v(d)-v(\beta)\).
The three-window holds at every start if \(1/2<s<3/4\).
It also holds if \(3/4<s<7/8\) and \(e/d\) has residue
different from \(1\). The proof below gives further explicit
necessary terminal congruences in higher valuation bands.

## Reduction to a root of unity

Divide by \(d\) and write \(b=\beta/d,\ q=e/d\),
so the phase becomes \(bXY+X+qY\) and \(v(b)=-s\).
Continue to use \(P_n\) after normalization. The boundary
identity in [matching and boundary lifts](ZERO_CONSTANT_BINARY_QUADRATIC_MATCHING_AND_BOUNDARY_LIFTS.md#a-valuation-one-lift-for-the-bilinear-three-term-phase)
gives
\[
   P_n=nbP_{n-1}+qP_{n-2}+1+q^n,\qquad n\ge2 .
\]
Hence a zero triple requires \(q^k=-1\). In particular \(q\)
is a unit. The same source gives
\[
  P_n=\sum_{j=0}^n j!b^j p_{n,j},\qquad
  p_{n,j}=[T^n] F_j(T),\qquad
  F_j(T)=\frac{T^j}{H(T)^{j+1}},\quad H=(1-T)(1-qT).
                                                        \tag{1}
\]
Let \(K=2^{v_2(k)}\), and write \(\bar q\) for the residue
of \(q\). Its odd order divides \(k\), since \(\bar q^k=1\).
If \(\bar q=1\), the order of \(q\) is \(2K\).

## The two parity layers

Suppose \(1-2^{-(m-1)}<s<1-2^{-m}\), with \(m\ge2\),
and put \(Q=2^m\). The multiplier \(j!b^j\) has valuation
\[
  \phi(j)=v_2(j!)-js
          =\sum_{r\ge0} j_r\big((1-s)2^r-1\big),
                                                        \tag{2}
\]
where \(j_r\) are the binary digits. The unique unrestricted
minimum is \(j=Q-1\), with odd-index residue series
\(T^{Q-1}/\bar H^Q\), up to normalization by a unit.
At even indices the smallest possible contribution comes
uniquely from \(j=Q-2\). This is the least cost among
even \(j\), and every odd \(j\) has a valuation gain of at
least \(1\) in \(p_{n,j}\) at an even \(n\), by the expansion
of \(H^{-j-1}\) modulo \(2\). The gain is greater than
\(\phi(Q-2)-\phi(Q-1)=s\). Its even residue series is
\(T^{Q-2}(1+\bar qT^2)/\bar H^Q\). All these expressions come
from congruences of universal integer polynomials.

For \(\bar q\ne1\), on the three residue classes of \(n\)
below, the corresponding residue coefficient is nonzero
unless the indicated power of \(\bar q\) equals \(1\):
\[
\begin{array}{c|ccc}
 n\bmod Q& Q-2&Q-1&0\\ \hline
 \text{power of }\bar q& n+2&n+1&n .
\end{array}                                               \tag{3}
\]
Here \(n>0\). The coefficients are geometric sums in
\(\bar q^Q\), with an extra factor \(\bar q\) in the last
class. The positive indices in each class are at or after
its first possible coefficient. It follows that a zero
triple cannot end at \(k\equiv\pm1,\pm2\bmod Q\), since a
relevant exponent in (3) is \(k\pm1\) or \(k\pm2\).

For \(\bar q=1\), the same geometric sums are nonzero
exactly when their lengths are odd. They exclude
\(k\bmod 2Q\in\{Q-2,Q-1,Q,Q+1,Q+2\}\).

## Exact residue at the middle index

If \(q\ne1\), the substitution \(T=1/(qS)\) interchanges
the poles at \(1\) and \(q^{-1}\) and takes \(F_j(T)\,dT\)
to \(-F_j(S)\,dS\). There is no residue at infinity of
\(T^{-n-1}F_j(T)\,dT\) for \(n\ge0\). The sum of residues gives
\[
 p_{n,j}
   =-\operatorname{Res}_{T=1}
        \left(T^{-n-1}-q^{n+1}T^{n+1}\right)F_j(T)\,dT .
                                                        \tag{4}
\]
This is a formal partial fraction identity. At \(n=k-1\)
the multiplier in parentheses becomes \(T^{-k}+T^k\).
With \(u=T-1\),
\[
     T^{-k}+T^k=2+(T^k-1)^2/T^k
                  \equiv 2 \pmod{(4,u^{2K})}.             \tag{5}
\]
Indeed \(T^k-1\in(2,u^K)\) by Frobenius, so its square is
zero modulo \((4,u^{2K})\).

Write \(\epsilon=1-q\). Up to sign,
\[
 \operatorname{Res}_{T=1}(F_j\,dT)
     =\epsilon^{-j-1}
       \sum_{r=0}^{j}\binom jr\binom{j+r}{r}(q/\epsilon)^r .
                                                        \tag{6}
\]
For \(r>0\) the binomial product is even: oddness of the
first factor requires all bits of \(r\) to occur in \(j\),
whereas oddness of the second requires those bits to be
disjoint. Both follow by expanding powers in characteristic
two.

If \(\bar q\ne1\), then \(\epsilon\) is a unit. Equations
(4)--(6) give \(v(p_{k-1,j})=1\) for every \(j<2K\).
When \(K\ge Q\), the minimizer \(Q-1\) in (2) is below
\(2K\), and every index \(j\ge2K\) has cost at least
\(\phi(Q-1)+2K(1-s)-1>\phi(Q-1)+1\). Thus \(P_{k-1}\)
has a unique term of least valuation. This excludes
\(Q\mid k\) in every stated band for \(\bar q\ne1\).
In particular, at \(Q=4\) it completes the exclusion
together with (3).

If \(\bar q=1\) and \(K\ge2\), then \(v(\epsilon)=1/K\):
\(\epsilon\) is a root of \((1-Z)^K+1\), which is Eisenstein
at \(2\). For \(j<K\), the \(r=0\) summand in (6) dominates,
by the evenness just proved. In (5), a correction coefficient
at \(u^h\) for \(0\le h\le j\) has valuation at least \(2\).
The contribution to the residue therefore has valuation at
least \(2-(2j+1-h)/K > 1-(j+1)/K\). It follows that
\[
                   v(p_{k-1,j})=1-(j+1)/K,\qquad j<K.  \tag{7}
\]
This shifts the least-cost comparison below \(K\) from
\(s\) to \(s+1/K\).

## The third open band when \(\bar q\ne1\)

For \(3/4<s<7/8\) it remains to address the classes
\(k\equiv3,4,5\bmod 8\) after (3) and the middle-index
residue argument.

Consider an even index \(n\equiv2,4\bmod8\). At these
indices the least possible cross part is \(j=2\). Here
are all the comparisons. Among even \(j\), those with
\(j\equiv6\bmod8\) have residue zero, since
\(\bar F_j=T^j\bar H/\bar H^{8a}\) for an integer \(a\ge1\);
at an even index its support lies in \(6,0\bmod8\).
The resulting gain of \(1\) exceeds \(\phi(2)-\phi(6)=4s-3\).
All other even \(j\) have least cost at \(j=2\), uniquely.
For odd \(j\equiv1\bmod4\), cost is at least \(\phi(5)\)
and the gain is at least \(1\). For odd \(j\equiv3\bmod4\),
cost is at least \(\phi(7)\) and the gain is at least \(2\).
The last gain follows because an odd coefficient of
\((1-T)^{-a}\) is divisible by \(2^{v_2(a)}\), by
\(\binom{a+r-1}{r}=(a/r)\binom{a+r-1}{r-1}\) for odd
\(r\), and the same holds in \(H^{-a}\). Both totals exceed
\(\phi(2)\). In the \(j=2\) summand,
\(\bar F_2=T^2\bar H/\bar H^4\). At \(n=k-1\) with
\(k\equiv3,5\bmod8\), its coefficient is a nonzero
geometric sum, since the relevant exponents are \(k+1,k-1\)
respectively and \(\bar q^k=1\). This excludes these two
classes. For \(\bar q=1\), the same comparison shows that
every even index \(2,4\bmod8\) has a nonzero exact moment,
because the geometric sum length is odd. Together with
the earlier parity layer, the remaining terminal classes
for \(\bar q=1\) in this band are \(0,1,15\bmod16\).

For \(k\equiv4\bmod8\) and \(\bar q\ne1\), one has \(K=4\).
At the middle index, all \(j<8\) have their coefficient
valuation \(1\); all \(j\ge16\) have cost greater than
\(\phi(7)+1\). Among \(8\le j<16\), the only term
of smaller cost than \(j=11\) is \(j=15\). Its residue at
\(n=k-1\) is zero because \(\bar F_{15}=T^{15}/\bar H^{16}\)
is supported on \(15\bmod16\), so it gains at least \(1\)
and its resulting valuation exceeds \(\phi(11)\).
Also \(\phi(11)-\phi(7)=4(1-s)<1\).
Finally \(\bar F_{11}=T^{11}\bar H^4/\bar H^{16}\)
gives a nonzero geometric sum at \(n=k-1\), with exponent
\(k+4\) if \(k\equiv12\bmod16\), and \(k-4\) if
\(k\equiv4\bmod16\). The latter case starts at \(k=20\):
\(k=4\) is impossible because the odd order of \(\bar q\)
exceeds one and divides \(k\). These exponents are
not multiples of that odd order. Thus \(j=11\) is the unique
least term and excludes the final class.

## The second open band when \(\bar q=1\)

If \(k\) is odd, then \(q=-1\). Put \(k=2r+1\) and \(X=b^2\).
Since \(H=1-T^2\), (1) gives
\[
 E_r(X):=P_{2r}=\sum_{t=0}^r c_t X^t,\qquad
 c_t=\frac{(r+t)!}{(r-t)!},\qquad
 P_{2r+1}/b=(r+1+X\partial_X)E_r .
\]
Thus a zero triple requires a repeated zero of \(E_r\).
The valuation of \(c_t/c_{t-1}=(r+t)(r+1-t)\) is at least
one. In each pair \(t=2i-1,2i\), the two halved even
factors are consecutive in parity, since they sum to an
odd integer. At least one ratio gains another factor \(2\),
so \(v(c_t)\ge t+\lfloor t/2\rfloor\).
For \(1<2s<3/2\), a nonconstant term can have valuation
at most zero at \(X=b^2\) only if \(t\) is odd and
\(v(c_t)=(3t-1)/2\). Among these candidates the term
valuations strictly increase with \(t\). Thus a zero of
\(E_r\) requires a unique nonconstant term of valuation
zero to cancel the constant. Its coefficient in
\(X E_r'\) has the same valuation since \(t\) is odd;
all other terms have higher valuation. There is no
repeated zero.

For even \(k\), the elementary parity layer leaves only
\(k\equiv0\bmod8\), so \(K\ge8\). At the middle index
equation (7) holds, and all \(j\ge K\) have cost at least
\(\phi(3)+K(1-s)-1>\phi(3)+1\). Below \(K\), the minimum
is unique unless \(s+1/K=3/4\): the least index is \(3\)
or \(7\) as follows by applying (2) at \(s+1/K\).
These two indices are both below \(K\), and already
\(j=3\) has total valuation less than \(\phi(3)+1\).
This excludes the strict alternatives.

### The tie at \(s=3/4-1/K\)

At \(n=k\), (4) has multiplier \(M=T^{-k-1}+qT^{k+1}\).
Split \(M_0=T^{-1}+qT\) from \(M-M_0\). The residue of
\(M_0 F_j\,dT\), up to sign, is \(N_j(\epsilon)/\epsilon^{2j+1}\),
where for \(j\ge1\)
\[
 N_j(Z)=\sum_{r=0}^j
  \binom{j+r}{r}
  \left(\binom{j-1}{j-r}+(1-Z)\binom{j+1}{j-r}\right)
  (1-Z)^r Z^{j-r}.                                      \tag{8}
\]
Every \(N_j\) is even coefficientwise. Indeed, over
\(\mathbb F_2(q)\) the residue identity (4) at \(n=0\)
identifies this residue with \(-[T^0]F_j=0\).
For the following low indices (8) gives the indicated
lower bounds on its valuation, including the denominator:
\[
\begin{array}{c|rrrrrr}
j&1&2&3&5&6&7\\ \hline
v(N_j(\epsilon)/\epsilon^{2j+1})\ge
&2-3/K&1-4/K&3-7/K&2-7/K&1-8/K&4-15/K .
\end{array}
\]
For direct verification of the table, \(N_1=4-4Z\),
\(N_2=6Z^2-18Z+12\), and \(N_3,N_7\) are coefficientwise
divisible by \(8,16\), respectively. In \(N_5\) all
coefficients are divisible by \(4\), and those of degree
below four are divisible by \(8\). In \(N_6\) all are
divisible by \(2\), and those of degree below five
are divisible by \(4\). The table entries plus \(\phi(j)\)
are strictly positive at the tie for \(K\ge8\).
The \(j=4\) cost is positive and \(p_{k,4}\) is integral.

Write \(A=T^k-1\). The coefficients of \(u^h\) in
\(T^{\pm k}-1\), for \(1\le h<K\), have valuation at least
\(\log_2 K-\lfloor\log_2 h\rfloor\), by the binomial identity
used above. In particular the contribution of \(M-M_0\)
at \(j<K\) has valuation at least
\(\log_2 K-\lfloor\log_2 j\rfloor-2j/K\).
This follows by expanding \(F_j\,dT\) locally as in (6),
and noting that the smallest exponent \(h\) supplied by
\(M-M_0\) is one. For \(K\ge16\) and \(1\le j\le7\)
this is enough after adding \(\phi(j)\), since even its
minimum at \(j=7\) plus \(\phi(3)\) is positive.
Every \(j\ge16\) has positive cost at this tie. For
\(8\le j\le15\) the bound
\(\phi(j)\ge-1/4+11/K\) is enough at \(K\le32\).
For \(K\ge64\) use the numerator evenness bound
\(1-31/K\) and the correction bound
\(\log_2 K-3-30/K\); both still give strictly positive
total valuation. Indices outside \(0\le j\le k\) need
no comparison.

For \(K=8\) the cost of \(j=6\) is already positive, and
the correction estimate suffices at \(j=1,2,3\).
At \(j=5\), terms from a multiplier exponent \(h\le3\)
gain at least \(2\) before denominator loss at most
\(10/8\). For \(h\ge4\), the local denominator exponent
is \(6+r\) with \(r\le1\). The correction gains at least
\(1\), with another gain of \(1\) when \(r=1\) from
\(\binom{5+r}{r}\). The resulting bound is at least \(1/4\),
which exceeds \(-\phi(5)=1/8\).
At \(j=7\), for \(r>0\) in the local expansion,
\(\binom{7+r}{r}\) is divisible by \(2\). For \(h\le3\)
the multiplier gains at least \(2\) and the denominator
loss is at most \(14/8\), so the combined bound exceeds
\(3/8\). For \(h\ge4\) one has \(r\le3\) and, at \(r>0\),
\(\binom{7+r}{r}\) is divisible by \(4\), which gives
the same strict inequality. For \(r=0\), use
\[
     M-M_0=A(qT-T^{-1})+A^2 T^{-1}/(1+A).
\]
The \(A^2\) part gains at least \(2\) in every relevant
degree. In the first part the \(r=0\) numerator is
\(A(qT^8-T^6)\). For the \(A\) coefficients of \(u^i\)
with \(i\le3\) the gain is at least \(2\). For
\(i=4,5,6,7\), the required coefficients of
\(qT^8-T^6\) have degrees \(3,2,1,0\), respectively.
Their valuations are at least \(2,0,1,1/8\), while the
\(A\) coefficients gain \(1,3,2,3\), respectively.
Thus the correction in the \(r=0\) part has valuation
at least \(2-8/8=1\). This exceeds \(-\phi(7)=3/8\).

All \(j\ge1\) contributions to \(P_k\) in (1) therefore
have strictly positive valuation at this tie.
The \(j=0\) coefficient is \((1+q)/(1-q)\), a unit since
\(v(1-q)=1/K<1\). This excludes the final tie and proves
the three-window for \(1/2<s<3/4\).
