# Exclusion of a fourth-order tail at a quadruple leading root

Work over a characteristic-zero field on \(2xy+t^2=1\).
Let \(\Lambda\) take the angular constant term and integrate in
\(t\) over \([-1,1]\) with measure \(dt/2\). Put \(U=1+x\) and
\(v=Ut\).

## Proposition (Fourth-order tail exclusion)

An eventually pure regular winding-one cubic of least sphere degree
seven cannot have the form
\[
 f=\frac Ux\left(P_0(v)+
             x^4(a_4v+b_4)+x^5(a_5v+b_5)+b_6 x^6+b_7 x^7\right),
 \qquad P_0(v)=(v^2-1)(v+d).
\]

Suppose it has this form and is eventually pure. The
[quadruple-root boundary corollary](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#consequences-in-degree-seven)
makes every positive moment zero. The [high-tail change of variables](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_HIGH_TAIL.md#a-linear-moment-block)
preserves the moments and writes the phase as
\[
 \frac{P_0(T)}X+\sum_{j=4}^7
          \frac{X^{j-1}}{(1-X)^j}(a_jT+b_j),\qquad a_6=a_7=0.
 \tag{1}
\]
Use the interval moments \(I_n=\int P_0(T)^n\,dT/2\),
\(J_n=\int T P_0(T)^n\,dT/2\), and \(K_n=\int T^2 P_0(T)^n\,dT/2\).
Write \(B_j(T)=a_jT+b_j\), and use the convention that
\(\binom{n}{k}=0\) for \(k>n\).

### The linear rows

The moments at \(m=4,5,6,7\) have only linear contributions from
the tail: a term with at least two factors from the tail has
order in \(X\) at least \(8-m>0\). The row at \(m=n+1\)
divided by \(m\) is therefore
\[
 R_n=\left(\binom n3J_n,\binom n4J_n,\binom n3I_n,
                 \binom n4I_n,\binom n5I_n,\binom n6I_n\right)
 \tag{2}
\]
on the column \((a_4,a_5,b_4,b_5,b_6,b_7)^t\).

Two other linear equations follow from the
[inverse-period criterion](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_INVERSE_PERIOD.md).
With \(s=(x/U)^{1/3}\) and \(p,q\) the depressed coefficients of
the cubic in \(v\), the coefficients of \(w^{-1},w^{-2}\) in
the inverse branch are \(-p/(3Us)\) and \(-q/(3Us^2)\).
Their constant terms at \(x=\infty\) vanish. Equivalently,
\[
 [x^1](1+x^{-1})^{-2/3}p=0,\qquad
 [x^1](1+x^{-1})^{-1/3}q=0.
\]
Since \(d\) is constant, the coefficients of \(x^j\) in \(p,q\)
for \(j\ge4\) are \(a_j,b_j-da_j/3\), respectively. The
first row becomes \(a_4-(11/12)a_5=0\). The second is
\(\sum_{j=4}^7\beta_j(b_j-da_j/3)=0\), where
\(\beta_j=\binom{-1/3}{j-1}\).

Let \(M\) have first row \((1,-11/12,0,0,0,0)\), followed by
the four rows \(R_3,R_4,R_5,R_6\). For column \(j\) set
\[
 k_j=\frac{(-1)^j}{d}\det M_{\widehat j}.
 \tag{3}
\]
These are polynomials in \(d\) over \(\mathbb Q\). At \(d=0\) the
rows \(R_3,R_5\) have entries only in the two \(a\)-columns,
and \(R_4,R_6\) have entries only in the four \(b\)-columns.
With the first row, the rank is at most four, so all the maximal
minors are divisible by \(d\).

### Two higher moments and arithmetic verification

For \(d\ne0\), the column \(k\) in (3) is nonzero and spans
the kernel of \(M\), as verified below. On this column the second
inverse-period row is
\[
 D=\sum_{j=4}^7\beta_j(k_{j-1}-(d/3)a_j(k))
   =\frac{268435456}{470507587106131575}(d^2-1)^2S(d^2),
 \tag{4}
\]
where \(a_4(k)=k_1,\ a_5(k)=k_2,\ a_6(k)=a_7(k)=0\), and
\[
 S(z)=196911z^7+446947z^6+223883z^5-94449z^4
                    -7507z^3+17161z^2+3417z+69.
\]
For \(m=8,9\), let \(L_m=R_{m-1}k\), and let \(Q_m\) be the
quadratic contribution on \(k\) divided by \(\binom m2\).
Three tail factors have order in \(X\) at least \(12-m>0\).
From (1), the only contributing index pairs for \(m=8\) are
\((4,4)\), and for \(m=9\) they are \((4,4),(4,5),(5,4)\).
Their respective angular coefficients are \(1\) and \(8,1,1\).
Thus
\[
 \begin{aligned}
 Q_8&=\frac12\int_{-1}^1 B_4(T)^2 P_0(T)^6\,dT,\\
 Q_9&=\frac12\int_{-1}^1
     (8B_4(T)^2+2B_4(T)B_5(T))P_0(T)^7\,dT,
 \end{aligned}
\]
with all tail coefficients evaluated on \(k\).
If a nonzero multiple \(\lambda k\) is pure, its linear
condition is \(D=0\), and its next two moment equations after
division by \(\lambda\) imply
\[
 W=8\cdot36\,L_8 Q_9-9\cdot28\,L_9 Q_8=0.              \tag{5}
\]

For exact replay of (3)–(5) one can compute all interval moments
from
\[
 \frac12\int_{-1}^1T^\ell P_0(T)^n\,dT
 =\sum_{\substack{0\le j\le n\\j+\ell\ {\rm even}}}\ \sum_{h=0}^n
 \frac{(-1)^{n-h}\binom nj\binom nh}{2h+j+\ell+1}\,d^{n-j}
 \qquad (\ell=0,1,2).
\]
Expanding (2) and (3) and applying the row with entries
\(\beta_j\) gives (4) over \(\mathbb Q\). The following
smaller arithmetic verifies the required coprimality tests.
All denominators in these calculations are prime to \(29\).
Modulo \(29\), the coefficient lists for \(k_3,k_6\) in
increasing powers of \(z=d^2\) are
\[
 k_3:(7,3,8,-7,-12,-3,7,7,-11),\qquad
 k_6:(-9,10,-4,1,-4,-8,7,-11,-11).
\]
Starting from their monic versions in that order and making each
nonzero remainder monic after division, the successive
unnormalized leading coefficients are
\(1,13,11,7,24,17,8,10\), at degrees \(7,6,\ldots,0\).
In particular the two reductions are coprime. Each cofactor
omitting a \(b\)-column has degree in \(d\) at most \(17\),
because its determinant has one \(J_n\) entry with degree
at most \(n-1\) and the other moment entries have degree
at most \(n\) for \(n=3,4,5,6\). After division by \(d\),
both reduced \(k_3,k_6\) attain the bound \(16\) in \(d\).
This proves their coprimality over \(\mathbb Q\), and hence
the claim about \(M\) at \(d\ne0\).

Modulo \(29\), the remainder of \(W\) upon division by \(S(d^2)\)
is \(-5d^{12}+10d^{10}+8d^8-10d^6-10d^4+11d^2+9\).
On the monic versions of \(S(z)\) and this remainder in \(z\),
with the same normalization rule, the successive unnormalized
leading coefficients are \(24,3,9,19,14,17\), at degrees
\(5,4,\ldots,0\). Substituting \(z=d^2\) preserves the two
coprimality relations. The leading coefficient of \(S(d^2)\)
is a \(29\)-adic unit, so a common factor over \(\mathbb Q[d]\)
of \(S(d^2)\) and \(W\) would give a common factor of positive
degree after reduction. This argument does not require a degree
comparison for \(W\). In particular (4) and (5) leave only
\(d=\pm1\) when \(d\ne0\).

### The exceptional values

At \(d=\pm1\) the interval moments satisfy
\(J_n/I_n=dn/(3n+2)\). For \(d=1\), substituting
\(T=2u-1\) reduces this to the mean of \(2u-1\) with weight
\(u^{2n}(1-u)^n\) on \([0,1]\). Reflection gives the other
sign. The column \((11,12,-3d,-4d,0,0)^t\) satisfies the
first row of \(M\) and all rows (2), since
\(11\binom n3+12\binom n4=(3n+2)\binom n3\) and
\(-3d\binom n3-4d\binom n4=-dn\binom n3\).
It spans the kernel. In particular its linear eighth moment is
zero. Its quadratic eighth moment is a nonzero rational multiple
of \(\int_{-1}^1(11T-3d)^2P_0(T)^6\,dT>0\). Thus any pure
tail at \(d=\pm1\) is zero, also over an arbitrary
characteristic-zero field.

At \(d=0\), the fourth moment forces \(a_4=0\) because
\(I_3=0\) and \(J_3\ne0\). The first inverse-period row then
gives \(a_5=0\). The linear eighth moment vanishes by parity,
and its quadratic contribution is \(\binom82 b_4^2 I_6\).
Since \(I_6\ne0\), we have \(b_4=0\). The fifth moment yields
\(b_5=0\), and the seventh gives \(6b_6+b_7=0\).
The second inverse-period row now gives
\(b_6-(8/9)b_7=0\), hence \(b_6=b_7=0\).

In every case eventual purity forces all six tail coefficients
to vanish. The resulting phase has sphere degree at most six,
as in the [high-tail proposition](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_HIGH_TAIL.md#proposition-high-tail-exclusion).
This contradicts the degree-seven hypothesis.
