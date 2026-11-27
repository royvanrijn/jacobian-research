# Constant shifts of binary phases with prime-power axis monomials

Let \(E\) have characteristic zero and a real-valued nonarchimedean
valuation extending \(v_p\), normalized by \(v(p)=1\). Write
\(R_n=\Lambda(q^n)/n!\), where \(\Lambda(X^iY^j)=i!j!\).
Assume \(R_0=1\), the \(R_n\) are integral, and their residue series is
\[
 {\cal R}(T)=\frac{1+cT}{(1-aT)(1-bT)}.                                      \tag{1}
\]
Here \(a,b,c\) lie in the residue field and either factor may be \(1\).
All pole tests refer to the reduced fraction. We use
\(e_j(z)=\sum_{i=0}^j z^i/i!\) for \(0\le j<p\) over the residue field.

**Window theorem.** Suppose \({\cal R}\) is nonpolynomial. Every
two-window of \(q\) at a positive start contains a nonzero moment.
For any \(k\ne0\), the moments of \(q+k\) have a three-window, with the
following allowed enlargement when \(v(k)=0\): use four if \(p>3\) and
\(c\ne0\). The constant may belong to any real-valued extension of
\(E\). For this assertion assume either of the following.

- The coefficient field has \(v(E^\times)\subseteq (1/e)\mathbb Z\)
  for an integer \(e>0\) prime to \(p\).
- There are two distinct geometric poles, with frequencies \(a,b\).
  For the specified constant, if
  \(\delta=v(k)=(1-1/L)/(p-1)\), \(L=p^h,\ h\ge1\), set
  \(\mu=\overline{k^L/(-p)^{(L-1)/(p-1)}}\).
  There is no \(j<p\) for which both
  \(e_j(\mu/a^L)=e_j(\mu/b^L)=0\).
  At every other valuation there is no extra condition on the constant.

The two-pole alternative is uniform in the constant at \(p=2\).
At \(p=3\) it is uniform if \((a/b)^2\ne-1\). If that equality holds,
a failure of the two-pole test requires the simultaneous zeros in
the \(j=2\) row. The coefficient-field alternative also covers that row.

**Binary support consequence.** Write the nonzero terms of \(q\) as
\(u_{ij}X^iY^j/(i!j!)\) with integral natural coefficients. Suppose unit
axis monomials occur only at \(X^{Q_x}\), \(Y^{Q_y}\), where \(Q_x,Q_y\)
are powers of \(p\) including \(1\). Every other unit monomial satisfies
\[
             i<Q_x,\qquad j<Q_y,\qquad i/Q_x+j/Q_y\ge1.                     \tag{2}
\]
Either of the designated axis coefficients may be missing or a nonunit.
Then (1) holds with \(a=\bar u_{Q_x,0}\), \(b=\bar u_{0,Q_y}\) and
\(c\) the sum of the other unit coefficient residues. If the series is
nonpolynomial, the window theorem under its stated hypotheses implies
the actual-support Strong Factorial window for \(q+k\). In particular
the criterion applies to homogeneous binary forms of any prime-power
degree plus a constant and to inhomogeneous forms satisfying (2).
Nonunit monomials have no exponent restriction.

At \(p=3\) the two-pole case at \(v(k)=1/3\) also has its
actual-support window without a grid bound if (2) holds for every
other coefficient with valuation less than one, and the designated
axis monomials are the only axis terms of valuation less than one.
This uses the stronger reduction in Section 6.

## 1. Getting the residue fraction

The multinomial coefficient of \(\prod u_\nu^{m_\nu}\) in \(R_n\) is
the integer
\[
 \frac{\prod_i(\sum_\nu m_\nu\nu_i)!}
 {\prod_\nu m_\nu!\prod_\nu(\nu!)^{m_\nu}},\qquad \sum_\nu m_\nu=n.
\]
Use the [native digit valuation and unit-factorial formula](FACTORIAL_MOMENT_SERIES_PRIME_POWER_RATIONALITY.md#degree-bounds-in-characteristic-p).
Only unit coefficients can contribute. Every counted nonsingleton
has count less than \(p\), and the combined base-\(p\) loads must be
carry-free. The contributions of all nondesignated unit monomials are
below the slots of the corresponding axis monomials. Thus they
contribute total loads \(I<Q_x,\ J<Q_y\) when carry-free. By (2) their
total count is at most \(I/Q_x+J/Q_y<2\). Zero insertions give numerator
\(1\), one gives its own natural coefficient residue times \(T\).
The designated axis monomials do not interact with those lower slots
or with one another. Their separate normalized pure-power coefficients
are one modulo the maximal ideal by the same unit-factorial formula.
Summing their counts gives the denominator of (1).

## 2. The normalized constant-shift rows

Put \(E_m(k)=\Lambda((q+k)^m)/m!\). The binomial theorem gives
\(E_m(k)=\sum_{j=0}^m R_{m-j}k^j/j!\).
With \(\delta=v(k)\), the cost of a term before including \(R_{m-j}\)
is
\[
                      \frac{s_p(j)-[1-(p-1)\delta]j}{p-1}.                 \tag{3}
\]
If \(\delta<0\), this cost strictly decreases at every step, since
\(s_p(j+1)-s_p(j)\le1\). The monic term using \(R_0\) is uniquely least.
For \(\delta\ge1/(p-1)\) every \(j>0\) has positive cost, and the
residue is \(\bar R_m\).

If \(0<\delta<1/(p-1)\) is noncritical, choose the unique power \(L\)
of \(p\) with \(1<L[1-(p-1)\delta]<p\). In \(j=Lu+t\) with \(0\le t<L\),
the cost in (3) decomposes into a strictly decreasing function of
\(u\) and a sum of strictly positive digit weights for nonzero \(t\).
Thus the unique least index is \(L\lfloor m/L\rfloor\). A zero moment
forces \(\bar R_{m\bmod L}=0\).

At a critical value, put \(j_0=pL\lfloor m/(pL)\rfloor\) and
\(m=j_0+dL+t\) with \(0\le d<p,\ 0\le t<L\).
In (3) each place below \(L\) has positive weight, the place \(L\)
has zero weight, and the high block at multiples of \(pL\) strictly
decreases with its block index. The least indices are exactly
\(j_0+iL\) for \(0\le i\le d\).
Legendre's formula gives
\[
 v\big((j_0+iL)!/j_0!\big)=i(L-1)/(p-1).
\]
After dividing by \((-p)^{i(L-1)/(p-1)}\) the quotient reduces to \(i!\).
Indeed the higher digits of the two factorial indices agree, and the
only changed digit is \(i\) at \(L\), so the unit-factorial formula
cancels their other digit factorials. Dividing \(E_m(k)\) by
\(k^{j_0}/j_0!\) now gives the following residue, evaluated at \(W=1\):
\[
 F_{d,t}(W)=\sum_{i=0}^d \frac{(\mu W)^i}{i!}\bar R_{t+(d-i)L}.             \tag{4}
\]
The same calculation permits \(\delta=0\) by taking \(L=1\). In that
case \(t=0\), \(\mu=\bar k\), and the test is the prefix through degree
\(p-1\) of \(e_{p-1}(\mu T){\cal R}(T)\), repeated modulo \(p\).
The test always has entry \(1\) at index zero.

A nonzero proper rational series with denominator degree at most two
has no two consecutive zero coefficients, since the recurrence
propagates them backwards and forwards. In (1) a polynomial part
can occur only for denominator degree one, and has degree zero.
Its nonzero proper part then supplies nonzero coefficients at all
positive indices. This proves the unshifted two-window, the
noncritical two-window by the circular test and \(\bar R_0=1\),
and the assertion for \(\delta\ge1/(p-1)\).

At \(\delta=0\) a zero block cannot cross index zero modulo \(p\).
Multiplying the test series by \((1-aT)(1-bT)\) shows that a triple
ending at index \(2\le j<p\) requires
\[
              [T^j]e_{p-1}(\mu T)(1+cT)
                    =\frac{\mu^{j-1}}{(j-1)!}(c+\mu/j)=0.
\]
If \(c=0\) this is impossible. Two consecutive such vanishing
conditions are always impossible, so four suffices for \(c\ne0\).
For \(p\le3\) the prefix already rules out any three-window of zeros.

## 3. Two distinct poles at a critical valuation

Both partial-fraction coefficients \(A,B\) are nonzero. With
\(x=\mu/a^L,\ y=\mu/b^L\) we have \(x\ne y\) and
\[
       F_{d,t}(W)=A a^{t+dL}e_d(xW)+B b^{t+dL}e_d(yW).                    \tag{5}
\]
If \(e_d(x),e_d(y)\) are not both zero, this row has no adjacent zero
test entries by the two distinct frequencies. A three-window
meets at most two rows because \(L\ge2\), and so it cannot vanish
if all rows it meets obey this test. At \(p=2\) no \(e_d\) has two
roots. At \(p=3\) only \(e_2\) can have two roots, and they are
\(-1\pm i\), \(i^2=-1\). Their ratios are the two primitive fourth
roots of unity. Since \(L\) is a power of three, \(x/y=(b/a)^L\)
can equal such a ratio only if \((a/b)^2=-1\).

## 4. Multiplicity on a coefficient grid prime to p

Use the monic polynomial \(P_m(Z)=m!E_m(Z)\). If the moments
vanish from \(n\) through \(m\), the derivative identity
\(P_m^{(r)}=(m)_rP_{m-r}\) gives root multiplicity at least \(m-n+1\).
The following consequence of the coefficient grid retains the
multiplicity of the actual root in the lattice variable \(W\).

Assume \(\delta=(L-1)/(L(p-1))\), \(L=p^h\), \(h\ge1\), and that
\(F_{d,t}(W)\) is nonzero. For every polynomial over \(E\), the
exponents of terms of least \(\delta\)-cost are congruent modulo
\(L\). Indeed equality of two costs puts their difference times
\(\delta\) in \((1/e)\mathbb Z\), and the reduced denominator of
\(\delta\) is \(L\). After substituting \(Z=kV\) and dividing by a
least term's scalar coefficient, the reduction is a monomial in
\(V\) times a polynomial in \(W=V^L\). Such reductions multiply
up to residue units, by the integral-domain property of the residue
polynomial ring. A root \(k\) is algebraic over \(E\) if it is a root
of \(P_m\). Its separable minimal polynomial to the power of the
root multiplicity divides \(P_m\). The reduction of that minimal
polynomial vanishes at \(V=1\), so its polynomial in \(W\) vanishes
at \(1\). Formula (4) is exactly the polynomial in \(W\) for \(P_m\),
up to an overall unit and a monomial factor, because it is nonzero
before evaluation. Consequently the multiplicity of \(k\) in
\(P_m\) is at most \(\operatorname{ord}_{W=1} F_{d,t}(W)\).

For two distinct poles, suppose a three-zero block exists. In at
least one row it contains two entries, so that row has both
\(e_d(x)=e_d(y)=0\). In particular \(d\ge2\), hence \(L\ge3\).
If the last two entries of the block belong to this row, both
must give \(F'_{d,t}(1)=0\). At roots of \(e_d\) we have
\(e'_{d}(z)=e_{d-1}(z)=-z^d/d!\ne0\); thus the derivative in (5)
is a sum of two nonzero multiples of \(a^t,b^t\). It cannot
vanish at two adjacent \(t\). These polynomials in \(W\) are
nonzero because \(e_d(xW),e_d(yW)\) are linearly independent for
\(d\ge1\), by their constant and linear coefficients.

The remaining split has two entries in a row with \(e_d(x)=e_d(y)=0\)
and the last entry at the start of the next row. A wrap to row zero
already gives a unit. Otherwise put \(j=d+1<p\). The last entry
requires a root of multiplicity at least three in \(F_{j,0}\).
At \(z=x,y\) we have
\[
 e_j(z)=z^j/j!,\qquad
 \left.\frac{d^2}{dW^2} e_j(zW)\right|_{W=1}=-z^{j+1}/(j-1)! .
\]
The first identity uses \(e_{j-1}(z)=0\), and the second uses
\(e_{j-2}(z)=-z^{j-1}/(j-1)!\).
The equations from value and second derivative in (5) have
nonzero determinant because \(x,y\) are nonzero and distinct.
This contradicts the two nonzero partial-fraction coefficients.

For only one geometric pole, if the pole is double the coefficients
are \((A+Bn)a^n\). Among the last two orders of a three-window
there is a unit anchor \(\bar R_t\) because \(L\) is divisible by
\(p\) and \(A+Bt\) cannot vanish at successive residues modulo \(p\).
At such an anchor \(\bar R_{t+hL}=\bar R_t a^{hL}\) for every \(h\).
For a simple pole use any \(t>0\), which exists among the last two
orders, and the same formula holds even with a constant polynomial
part. In both cases (4) is a nonzero scalar times \(e_d(\mu W/a^L)\).
Every root of \(e_d\) is simple by the relation used above, so this
anchor cannot have multiplicity two. Its order is at least \(n+1\),
giving the required contradiction.

For the binary support consequence a nonpolynomial residue needs an
axis term. If \(q\) is a monomial, its moments are nonzero. If
\(q+k\) has only two terms with \(k\ne0\), use the already proved
[constant-plus-monomial two-window](MULTIVARIATE_CONSTANT_PLUS_MONOMIAL_DEGREE_FOUR.md).
The supported three-window covers all larger supports. When four
is needed, two distinct poles require two unit axis columns and
\(c\ne0\) requires a third unit column. Under the tame alternative
one axis may be absent. In that case \(c\ne0\) and a three-term
phase has only a single linear denominator factor in (1). At
\(\delta=0\), a triple forces two adjacent zero coefficients of
\(e_{p-1}(\mu T)(1+cT)\) by the order-one recurrence. The explicit
coefficient in Section 2 forbids this. Thus four is needed only
when the actual support is at least four.

## 5. The odd-prime cancellation is realized

Take \(p=3\), \(i^2=-1\), and \(q=(X^3+iY^3)/3!\) over
\(\mathbb Q_3(i)\). Here \({\cal R}(T)=1/((1-T)(1-\bar iT))\)
has two poles. For every \(L=3^h\), \(h\ge1\), put
\(\rho=\bar i^L\), \(\mu=-1-\rho\). Then
\(e_2(\mu)=e_2(\mu/\rho)=0\), since
\(\mu/\rho=-1+\rho\). Realize the constant by choosing
\(\pi^L=3\) and \(k=\pi^{(L-1)/2}u\), with a unit \(u\) whose
residue satisfies \(\bar u^L=(-1)^{(L-1)/2}\mu\).
Such a residue exists by Frobenius in \(\mathbb F_9\).
The entire \(d=2\) row of \(L\) consecutive leading tests vanishes.
The coefficient field \(\mathbb Q_3(i)\) itself has value group
\(\mathbb Z\), since \(X^2+1\) has separable irreducible reduction
modulo three. The actual three-window follows from Section 4.
Thus length-\(L\) first-residue blocks alone cannot be excluded by
having two poles at odd primes.

**Exact ternary-prime boundary.** For any two distinct poles in the
\(p=3\) canceled-row case, the maximum leading zero-block length is
\(L+1\) if \(c=0\), and \(L\) otherwise. In particular there is an
actual \((L+2)\)-window when \(c=0\), and an actual \((L+1)\)-window
otherwise, without any coefficient-grid assumption.

To see this, the \(d=2\) row vanishes throughout. Row zero has
entry \(1\) at its start, so the block does not extend right.
In the \(d=1\) row the two multipliers \(e_1(x),e_1(y)\) are
nonzero and opposite because \(x,y=-1\pm i\) in either order.
Hence no two entries in row one vanish consecutively. At its
last entry, \(t=L-1\), cancellation occurs exactly when
\(A a^{2L-1}=B b^{2L-1}\). Put \(r=b/a\), so \(r^2=-1\) and
\(r^{2L-1}=r\) because \(L\) is odd. The condition is \(A=Br\),
equivalent to \(c=-(Ab+Ba)=0\). Every block disjoint from the
vanishing row has length at most two. This proves the exact
leading lengths and their stated actual-moment consequence.
For example \(L=3\) already excludes a five-zero block in the
two-pole case over arbitrary real-valued \(3\)-adic fields.

For the same \(q\), the separate choice \(k=-(1+i)\) illustrates
sharpness of the three-window. Direct expansion gives first and
second moments zero, and the third moment \(1622(1-i)\).
Indeed for one axis \(X^3/3!\) the first three moments are
\(1,20,1680\). Centering gives third moment \(1680-3\cdot20+2=1622\),
and independence of the two axes multiplies this by \(1+i^3\).

## 6. Lifting the first canceled row over the full modulo-three ring

**First critical row.** At \(p=3,\ v(k)=1/3\), suppose the two
distinct-pole support consequence has the enhanced exponent
hypothesis stated above. Then the actual-support window holds
over arbitrary real-valued extensions, including wild
coefficient grids. In this sector four consecutive zeros are
impossible. When the support has size three, three consecutive
zeros are impossible.

First, the finite reduction preceding (4) has a version modulo
\(p\mathcal O\), for every \(p,L\) and each \(m=pLH+i,\ 0\le i<pL\)
at the critical valuation. In the normalized sum for \(E_m(k)\),
each \(j=pLH+\ell\), \(0\le\ell\le i\), has exactly the contribution
\[
              \frac{k^\ell}{\ell!}\,R_{i-\ell}
                       \binom{pLH+\ell}{\ell}^{-1}.
\]
The binomial is \(1\bmod p\) by the no-carry digit formula. The
factor \(k^\ell/\ell!\) is integral by (3) and the nonnegative
low-place weights. For every \(j<pLH\) the cost relative to the
normalization is at least one. Indeed the adjacent high-block
cost drop from \(H-1\) to \(H\) is \(1+v_p(H)\), and earlier
drops are positive. The remaining low-place cost is nonnegative.
Consequently an actual zero at \(m\) gives
\[
                         E_i(k)=0\pmod{p\mathcal O},                     \tag{6}
\]
where \(\mathcal O\) is the valuation ring of the field containing
\(k\) and the coefficients.

Now put \(p=L=3\) and work in \(\mathcal A=\mathcal O/3\mathcal O\),
which can have nilpotents. The proof in Section 1 gives (1) in
\(\mathcal A[[T]]\), with \(a,b,c\) the reductions in this ring of
the designated coefficients and the sum of the other coefficients.
In fact divisibility by \(3\) of a native integer count coefficient
annihilates its entire term in \(\mathcal A\); the remaining
count vectors are still carry-free. All coefficients of valuation
at least one disappear, and the stronger exponent condition puts
the remaining ones in the same disjoint slots. The one-insertion
numerator and the two geometric sums therefore hold in the full
ring.

With \(U=k^3/(-3)\), which is a unit, the binomial expansion for
\(E_i(k)\) through \(i=8\) gives
\[
 \sum_{i=0}^8 E_i(k)T^i
 \equiv \frac{(1+cT)e_2(kT)e_2(UT^3)}{(1-aT)(1-bT)}
             \pmod{T^9}\quad\text{in }\mathcal A[[T]].                 \tag{7}
\]
For \(0\le d,t\le2\), the constant coefficient at \(T^{3d+t}\)
in the truncated exponential \( \sum_{j=0}^8 k^j T^j/j!\) agrees
with \(U^d k^t/(d!t!)\) modulo \(3\mathcal O\). This is the
unit-factorial formula with valuation \(v_3((3d+t)!)=d\),
and the quotient after extracting \(U^d k^t\) is an ordinary
rational unit. In particular the coefficients at degrees seven
and eight of the numerator in (7) are
\[
                  (U^2/2)(k+c),\qquad (U^2/2)(k^2/2+ck).              \tag{8}
\]

A four-zero block in the leading test can start only at index
five modulo nine, by the boundary in Section 5. By (6) this
forces \(E_5=E_6=E_7=E_8=0\) in \(\mathcal A\). Then multiplying
the actual prefix of degree at most four by the denominator
in (7) gives a polynomial of degree at most six, so both
coefficients in (8) must vanish. Since \(U\) is a unit, \(c=-k\)
and \(k^2/2+ck=0\) in \(\mathcal A\). This contradicts \(v(k^2)=2/3\).

With actual support size three there are exactly the two unit
axis terms and the nonzero constant, so \(c=0\) already in the
coefficient field. The candidate triple starts at five or six
modulo nine. At five the actual prefix has degree at most four
through index seven, and the numerator coefficient at seven
must vanish. At six the prefix has degree at most five through
index eight, and the numerator coefficient at eight must vanish.
Equation (8) contradicts each requirement when \(c=0\), since
\(v(k)=1/3\).

Without the enhanced exponent hypothesis, the exact leading
length in Section 5 already covers supports of size at least five.
Support size three has just been covered because there is no
foreign column. Thus only size four can remain at \(L=3\).
If its one additional nonconstant column is outside (2), or is
an axis term at another position, its natural coefficient \(u\)
has \(v(u)>0\) by the original residue hypothesis. In fact a
four-zero block requires \(v(u)\le1/3\). To prove this, let
\(q_0\) consist of the two designated axis monomials. By native
integrality, \(R_i(q)-R_i(q_0)\) has valuation at least \(v(u)\).
For \(i\le8\) all factors \(k^\ell/\ell!\) in the normalized shift
are integral, so this estimate survives the shift and
multiplication by \((1-aT)(1-bT)\). For \(q_0+k\) the numerator
coefficient at seven is \(U^2 k/2\bmod 3\mathcal O\) by (8).
It has valuation \(1/3\), and therefore cannot vanish
modulo \(3\mathcal O\) after a correction of higher valuation.
As the leading four-zero block must start at five, (6) and
the degree bound force that coefficient to vanish for \(q+k\).
This proves the necessary bound on the foreign coefficient.

## 7. A digit restriction in the larger canceled rows

**Larger critical rows.** Keep the \(p=3\) two-pole canceled-row
case, and the full modulo-three exponent hypothesis of Section 6.
Write the possible support-size block as \(pLH+s,\ldots,pLH+s+w-1\),
where \(w\ge3\) and \(s\) is its initial residue. Every integer
\[
              j\in[s+2-2L,\ s+w-1-2L]
\]
is in \([1,L-1]\). If any such \(j\) has \(s_3(j)\le2\), then
\(3\nmid j\) and \(v(c)=v(k)\) for the actual coefficient sum \(c\).
Moreover all such \(j\) must lie in the same nonzero class modulo
three. In particular when \(c=0\) the indicated range can contain
only integers with \(s_3(j)\ge3\).

For the proof, use the full ring \(\mathcal A\) from Section 6.
For \(j=dL+t<pL\), extracting \(U^d k^t\), with \(U=k^L/(-3)^{(L-1)/2}\),
the factorial unit formula gives \(k^j/j!=U^d k^t/(d!t!)\)
in \(\mathcal A\). Thus (7) holds through degree \(3L-1\) with
\(e_2(UT^L)\) and with \(e_2(kT)\) replaced by
\(K(T)=\sum_{t=0}^{L-1} k^t T^t/t!\). By (6) the actual prefix
through the block's last index has degree at most \(s-1\).
Multiplication by the two denominator factors forces its
numerator coefficients at \(s+2\) through \(s+w-1\) to vanish.
The leading boundary puts \(s\ge2L-1\), and the block ends before
\(3L\), as claimed. At degree \(2L+j\) in this range the coefficient
is \(U^2/2\) times the coefficient of \((1+cT)K(T)\) at \(j\).
Only the highest \(e_2(UT^L)\) term contributes because
\(\deg((1+cT)K)\le L\). Therefore the actual \(c,k\) satisfy
\[
         v(c+k/j)\ \ge\ 1-v(k^{j-1}/(j-1)!).                            \tag{9}
\]
By (3), \(v(k^j/j!)<1\) is equivalent to \(s_3(j)\le2\) in
\(1\le j<L\). At such \(j\), (9) demands a value strictly higher
than \(v(k/j)\), and hence \(v(c)=v(k/j)\). Integrality of \(c\)
rules out \(3\mid j\), since \(v(k)<1/2\). If two such unit indices
were distinct modulo three, \(v(k/j-k/j')=v(k)\), contradicting
the strict valuation bounds from (9) at the two indices.

For pure axes plus a constant, \(c=0\). At \(L=9\), the three-window
can therefore only start at \(21,23,24\bmod 27\) in the canceled-row
case. These are the three indices \(s=16+j\) for which
\(1\le j\le8\) and \(s_3(j)\ge3\).

The three possible starts impose more precise pole conditions. If
\(r=b/a\), an actual triple at \(s=21\) requires \(v(r^2+1)=2/27\).
At \(s=23\) or \(24\) it requires \(v(r^2+1)\ge1/9\).
In particular \(s=21\) is excluded whenever the coefficient-grid
denominator has 3-part less than \(27\).

Here is the calculation in the pure-axis case. Still work modulo
\(3\mathcal O\), with \(H(T)=e_2(UT^9)\). In \(K(T)\) the last
surviving coefficient is at degree six and has valuation \(2/3\).
The product \(HK\) has degree \(24\) with leading coefficient
\(w=U^2 k^6/(2\cdot 6!)\) of the same valuation. At \(s=23,24\),
the denominator \(D=(1-aT)(1-bT)\) divides \(HK\): the actual
prefix times \(D\) agrees through degree at least \(25\) and both
sides have degree at most \(25\). At \(s=21\) the same argument
through degree \(23\) gives divisibility of \(HK-wT^{24}\).
Since \(K(a^{-1}),K(b^{-1})\) both have residue one, this proves
\(v(e_2(U/a^9)),v(e_2(U/b^9))\ge1\) at the later starts. At
\(s=21\) both values are \(2/3\), and the ratio of the two
evaluations of \(e_2\) is one in the residue field. For the last
assertion use \(\bar r^{24}=1\) in the two expressions
\(w a^{-24}/K(a^{-1})\) and \(w b^{-24}/K(b^{-1})\).

In a valued extension write the two roots of \(e_2\) as
\(\xi_a,\xi_b=-1\pm i\), in the order of their residues at the
actual arguments \(z_a=U/a^9,\ z_b=U/b^9\). The different roots
are separated by a unit, so \(v(z_i-\xi_i)=v(e_2(z_i))\).
At \(s=21\), the errors \(\Delta_i=z_i-\xi_i\) have valuation
\(2/3\) and \(\Delta_a/\Delta_b\) has residue \(-1\), by the
opposite unit derivatives \(e_2'(\xi_i)=1+\xi_i\).
For the fourth root \(\zeta=\xi_a/\xi_b\) it follows that
\[
 v(r^9-\zeta)=v((\Delta_a-\zeta\Delta_b)/z_b)=2/3,
\]
because \(-1-\bar\zeta\ne0\). At the other starts this difference
has valuation at least one. In all cases \(r\) reduces to the
fourth root \(\zeta\) itself, since ninth power is Frobenius
squared and fixes fourth roots. If \(0<v(r-\zeta)<1/9\), the
ninth power of \(r-\zeta\) is the unique least term in
\(r^9-\zeta\) by the binomial coefficients, which are divisible
by three except at the ends. If \(v(r-\zeta)\ge1/9\), the
difference has valuation at least one. Also \(r+\zeta\) is a
unit, so \(v(r-\zeta)=v(r^2+1)\). This gives both asserted
constraints in the original coefficient field.

## 8. The foreign monomial at the first critical radius

**Four-term foreign-column reduction.** In the remaining \(p=L=3\)
case in Section 6, write the third nonconstant monomial as
\(uX^\alpha Y^\beta/(\alpha!\beta!)\). A four-zero block requires a
nonsingleton foreign monomial. Every digit of \(\alpha,\beta\) must
be \(0\) or \(1\), the total number of \(1\)'s must be even, and
\[
      \lfloor \alpha/Q_x\rfloor\equiv
      \lfloor \beta/Q_y\rfloor\equiv0\pmod9.
                                                                  \tag{10}
\]
Here a singleton means an axis monomial with prime-power exponent.
After dividing the whole phase by the unit axis coefficient \(a\),
continue to write \(b,k,u\) for the normalized coefficients, and
put \(U=k^3/(-3)\). An actual block starting at \(9H+5\) further
requires
\[
 \begin{split}
 v(u)&=1/3,\qquad v(u+k)\ge4/3,\\
 v\big((b^2+1)^3-6H b^3\big)&\ge4/3.
 \end{split}                                                       \tag{11}
\]
In particular \(v(b^2+1)=1/3\) if \(3\nmid H\), and
\(v(b^2+1)\ge4/9\) otherwise.

### The finite modulo-three calculation

Normalize \(a=1\). At the canceled row let \(i\) be the fourth root
reducing to \(\bar b\) in a valued extension. Then
\(\bar U=\bar i-1\). Indeed the roots of \(e_2\) are
\(-1\pm\bar i\), and the actual residues of the arguments at the
two poles are \(\bar U,\bar U/\bar b^3\).
The ring \(\mathcal A=\mathcal O/3\mathcal O\) therefore receives
the localization of \(\mathbb F_9[B,V,k,u]\) at its origin by
\(B=b-i,\ V=U-(i-1)\). Any polynomial with nonzero constant
term maps to a unit, since the four variables have positive
valuation. Also \(k^3=0\) in \(\mathcal A\) and \(k^2\ne0\).

Let \(x,y,m\) be the counts of the two axis and the foreign terms.
For a nonsingleton foreign term all native integer coefficients
with \(m\ge3\) are divisible by three, by the digit valuation
formula in Section 1. A carry on doubling \((\alpha,\beta)\)
also kills all \(m=2\) terms. Otherwise let \(w\) be the number of
\(1\) digits of \((\alpha,\beta)\). The foreign-only factorial
coefficient for \(m=2\) is \((-1)^{w-1}\) modulo three,
by the product of digit factorials and the division by \(2!\).
When the doubling term is killed assign it value zero instead.
Thus, with \(r=\lfloor\alpha/Q_x\rfloor\bmod9\),
\(s=\lfloor\beta/Q_y\rfloor\bmod9\) and that value denoted by
\(\eta\), the remaining native coefficient for \(x+y+m\le8\)
is
\[
   \binom{x+mr}{x}\binom{y+ms}{y}
       \begin{cases}1&m=0,1,\\ \eta&m=2.\end{cases}                 \tag{12}
\]
To check the binomial factors, first remove the three
within-column factorial coefficients. In the first coordinate
the remaining multinomial is the binomial for adding
\(Q_x x\) to \(m\alpha\). The digits below \(Q_x\) are disjoint
from \(x\), and the digit binomial is
\(\binom{x+\lfloor m\alpha/Q_x\rfloor}{x}\bmod3\).
When the \(m=2\) coefficient survives, the lower \(\alpha\)
digits can be doubled without carry across \(Q_x\), so the
quotient equals \(m\lfloor\alpha/Q_x\rfloor\).
The same argument applies to the other coordinate.
Counts less than nine read only two low digits in these
quotients. If either \(r\) or \(s\) has a \(2\) digit, only
\(\eta=0\) is used. The other residue pairs allow both signs and
zero, which gives a finite cover of all nonsingleton exponents.

For a foreign singleton, the two occupied positions in its
coordinate either differ by one or by at least two. In the
second case counts through eight cannot overlap and the native
coefficient is one. In the first case let \(\ell,h\) be the
counts at the lower and higher position. The only shared
load is that of \(\lfloor\ell/3\rfloor\) and \(h\), so the
coefficient is \(\binom{h+\lfloor\ell/3\rfloor}{h}\bmod3\).
There are four adjacent-position profiles, one for each choice
of coordinate and orientation, and the one separated profile.

These rules give 118 profiles. In each, \(R_n\) is obtained
by summing the native coefficients times \(b^y u^m\) at total
count \(n\). The four equations are the coefficients at
\(5,\ldots,8\) of
\[
       \left(\sum_{n=0}^8 R_nT^n\right)
                   e_2(kT)e_2(UT^3)
       \pmod{T^9}.                                                \tag{13}
\]
Sixteen profiles have a direct polynomial argument. With
\(D=(1-T)(1-bT)\), their reduction at \(k=0,\ U+1+b^3=0,\
-U^2-b^3=0\) has series
\[
  D^2\left(1+
     \sum_{m=1,2} \eta_m u^m T^m
                      (1-T)^{-mr}(1-bT)^{-ms}\right),
  \qquad \eta_1=1,\quad \eta_2=\eta.                                 \tag{14}
\]
For these profiles, replacing each negative exponent modulo
nine in (14) gives degree at most four in all present terms.
This replacement is valid modulo \(T^9\) by Frobenius. In
particular the four coefficients in (13) belong to the ideal
generated by \(k,\ A=U+1+b^3,\ C=-U^2-b^3\). The coefficient
at seven of \(D\) times (13) has nonzero \(k\)-coefficient
\(-(\bar i-1)^2\) at the origin and zero \(A,C\) coefficients.
The coefficients at five and six of (13) have an invertible
two-by-two \(A,C\) matrix at the origin. Indeed their \(A,C\)
contribution is from \((A T^3+C T^6)/D\), and
\([T^2]D^{-1}=\bar i\) there. The other \(A,C,k\) terms have
coefficients in the maximal ideal or are part of the two
unneeded \(k\)-coefficients at orders five and six. These three
linear combinations of \(k,A,C\) have a unit determinant
locally, so the four equations force \(k=0\) in \(\mathcal A\).

For the other profiles the [exact checker](verify_binary_first_wild_three_adic_row.py)
computes local standard bases for (13) and \(k^3\) in
\(\mathbb F_9[B,V,k,u]_{(B,V,k,u)}\) with Singular's
negative-degree local order. Its zero normal form tests local ideal
membership. [Singular 4.3.2, §5.1.149](https://www.singular.uni-kl.de/Manual/4-3-2/sing_357.htm)
and [§5.1.129](https://www.singular.uni-kl.de/Manual/4-3-2/sing_337.htm).
The reduction of \(k^2\) is zero
in 101 of those 102 profiles, establishing ideal membership
and contradicting \(v(k^2)=2/3\). The single surviving profile
is \((r,s,\eta)=(0,0,-1)\), giving (10) and the digit condition.
As a separate check of (12) and the singleton formulas,
the checker compares 19,470 individual coefficients with
integer factorial quotients. Run from the repository root
with Python 3.12 and Singular 4.3.2 using
\[
 \texttt{python3 research/factorial/image-rigidity/verify_binary_first_wild_three_adic_row.py}.
\]
This finite computation covers exponent profiles of unbounded
size by the digit reduction above.

In the surviving profile, \(DR=1+uT-u^2T^2\bmod T^9\) in
\(\mathcal A\). The degree-seven condition for \(D\) times
(13) gives \(-U^2(k+u)+U k^2u^2=0\). If \(v(u)<1/3\)
the term \(-U^2u\) has strictly least valuation below one.
Thus \(v(u)=1/3\) by Section 6, and \(u+k=0\) in \(\mathcal A\).
Now \(e_2(kT)(1+uT-u^2T^2)=1\) in \(\mathcal A[T]\).
The vanishing at five and six in \(e_2(UT^3)/D\) forces
\(e_2(U),e_2(U/b^3)=0\) in \(\mathcal A\). This follows from
the proper part of the fraction after its degree-four
polynomial part, whose two coefficients at successive
orders determine both partial-fraction coefficients
because the poles are distinct modulo the maximal ideal.
In particular \(v(U-(i-1))\ge1\) and \(v(b-i)\ge1/3\).

There is an explicit solution of the four congruences
in the full ring modulo three. Take \(Q_x,Q_y\ge3\),
the foreign term \(XY\), \(b=i\), \(k^3=-3(i-1)\), and
\(u=-k\). Then \(e_2(UT^3)=D^3\) in \(\mathcal A[T]\),
and (13) is \(D^2\). Its coefficients at five through eight
vanish, with \(v(k)=1/3\) as required. This example establishes
the need to lift the actual equations beyond this modulus.

### Lifting the survivor to valuation \(4/3\)

Write \(I_\gamma=\{z:v(z)\ge\gamma\}\). Let
\(h_n=\sum_{y=0}^n b^y\), the pure geometric coefficient with
\(a=1\), and set \(h_n=0\) at negative indices. For a pure axis
of exponent \(Q=3^h\), its normalized factorial coefficient
\((Qn)!/(n!(Q!)^n)\) is \(1\bmod9\). To see this, splitting
the multiples of three from \((3n)!\) gives the unit factor
\(\prod_{j=0}^{n-1}(3j+1)(3j+2)\equiv2^n\bmod9\). Iterating
the splitting for \(Q\) and normalizing by the \(n\)-th power
of the corresponding factors for \(Q!\) cancels them all.
Consequently the pure part of \((DR)_n\) is
\(\mathbf 1_{n=0}\bmod9\). For \(n\le8\), the coefficients of
the \(u,u^2\) parts may be reduced modulo three at precision
\(I_{4/3}\), because their valuations are at least \(1/3,2/3\).
Terms of degree at least three in \(u\) have coefficients
divisible by three by the nonsingleton estimate. Therefore
\[
           D R(T)\equiv 1+uT-u^2T^2
                \pmod{(T^9,I_{4/3})}.                              \tag{15}
\]

Let \(\widetilde E_i=E_{9H+i}(k)/(k^{9H}/(9H)!)\).
Its denominator is the same for all \(i\) in the block.
In \(\widetilde E_7-(1+b)\widetilde E_6+b\widetilde E_5\),
all terms with constant-power index below \(9H\) drop out
at precision \(I_{4/3}\): their scalar coefficients have
valuation at least one and the accompanying coefficient
of \(DR\) has valuation at least \(1/3\) at positive
indices, by the pure-axis congruence and \(v(u)=1/3\).
For the other terms, (15) and the binomial coefficient in
Section 6 give exactly
\[
                  k^7/7!+u k^6/6!-u^2 k^5/5!
                         \pmod{I_{4/3}} .
\]
The binomial correction costs at least one additional
valuation unit, and each retained term here has valuation
at least \(1/3\). At \(u=-k\) this expression is
\(-48k^7/7!\) of valuation \(4/3\). Its difference from the
expression at the actual \(u\) is \((u+k)\) times a unit:
the linear factor has \(k^6/6!\) of valuation zero and
the other contribution has positive valuation.
The three zero equations therefore give \(v(u+k)\ge4/3\).

We next compute \(\widetilde E_i\) itself for \(5\le i\le8\)
modulo \(I_{4/3}\). In the main block the inverse binomial
is \(1-3H H_d\bmod9\) for \(\ell=3d,\ d=0,1,2\), by separating
the factors at \(3,6\) in its product. Here \(H_0=0,\ H_1=1,\
H_2=1+1/2=0\bmod3\) are harmonic sums. The other \(\ell\)
have positive constant cost at least \(1/3\), so their
binomial correction vanishes at this precision.
The main block thus contributes \(E_i(k)-3HU h_{i-3}\),
using \(k^3/3!\equiv U\bmod3\).
If \(3\mid H\), the earlier blocks have cost at least two
and contribute nothing. Otherwise only the three indices
\(j=9H-3d,\ d=1,2,3\), can contribute from the previous block.
The ratio \((9H)!/j!\) reduced after dividing by \(3^{d+1}\)
is \(-H(d-1)!\): the nonmultiples of three give \((-1)^d\),
and the multiples of three after scaling give
\(H(-1)^{d-1}(d-1)!\). Dividing also by \(k^{3d}\)
gives the tail contribution
\[
             3H\big(U^{-1}h_{i+3}
                    -U^{-2}h_{i+6}-U^{-3}h_{i+9}\big).
\]
It vanishes at precision \(I_{4/3}\).
Indeed modulo the maximal ideal \(h_n\) is a linear
combination of the two frequencies \(z=1,\bar b\).
For each frequency \(x=\bar U/z^3\) satisfies \(x^2=1+x\),
and so \(x^{-1}-x^{-2}-x^{-3}=0\). This checks both terms in
the spectral decomposition for every \(H\).

By (15) we may substitute \(u=-k\). The numerator from the
constant expansion of \(E_i\) now has coefficient
\((1-j^2)k^j/j!\) at degree \(j\). Through degree eight,
only degrees zero, three and six remain modulo \(I_{4/3}\),
at which the coefficients are \(1,4U,-7U^2/16\).
Combining the main and tail calculations gives
\[
   \widetilde E_i\equiv [T^i]\frac{g_H(UT^3)}{D(T)}
                   \pmod{I_{4/3}},\qquad
      g_H(z)=1+(4-3H)z-\frac7{16}z^2 .                             \tag{16}
\]
If these normalized moments vanish, the two successive
coefficients at five and six again force
\(g_H(U),g_H(U/b^3)=0\bmod I_{4/3}\). The polynomial
part has degree at most four and both poles are unit-separated
in this quotient. The arguments \(z_1=U,\ z_2=U/b^3\)
are units with a unit difference, so Vieta's identities
hold in the quotient by \(I_{4/3}\).
They give
\[
 \begin{split}
 b^3+b^{-3}
 &=\frac{(z_1+z_2)^2}{z_1z_2}-2\\
 &\equiv -\frac{16}{7}(4-3H)^2-2\equiv 6H
                                    \pmod{I_{4/3}} .
 \end{split}
\]
The last congruence is ordinary reduction modulo nine of
a rational \(3\)-adic integer. Since \(v(b^2+1)\ge1/3\),
the terms in \(b^6+1-(b^2+1)^3=-3b^2(b^2+1)\)
have valuation at least \(4/3\). This proves (11).

The congruence at precision \(I_{4/3}\) is realizable for every
\(H\ge0\). In \(\mathbb Q_3(i)\), lift the two simple roots of
\(e_2\bmod3\) to roots \(z_1,z_2\) of \(g_H\), choosing
\(z_1\equiv i-1\bmod3\). Set \(U=z_1\) and in a valued extension
set \(b^3=z_1/z_2,\ k^3=-3U,\ u=-k\).
Then \(v(b-i)\ge1/3\), since \(z_1/z_2\equiv i^3\bmod3\).
Take the foreign monomial \(XY\) and axis powers at least three.
Equation (16) applies and \(g_H(UT^3)/D\) is an actual polynomial
of degree four because the two linear factors divide its
numerator. Thus \(\widetilde E_5,\ldots,\widetilde E_8\) all
vanish modulo \(I_{4/3}\). This family realizes the lifted
necessary test, with the four exact moment equations still open.

### The next necessary digit congruence

**Digit-run restriction.** In the surviving profile (10), let
\(d_x,d_y\) be the ternary digits of \(\alpha,\beta\) immediately
below \(Q_x,Q_y\), using zero when a slot is the units digit.
Let \(N_{\rm run}\) be the total number of maximal nonempty
runs of \(1\)'s in the two ternary expansions. An actual
four-zero block also requires
\[
                   N_{\rm run}+d_x+d_y\equiv1\pmod3.                 \tag{17}
\]
In particular the \(XY\) family with both axis powers at least
nine is excluded at the next precision. The family with both
axis powers equal to three still passes (17).

We prove the congruence by retaining all contributions of
valuation below two to
\(\widetilde E_8-(1+b)\widetilde E_7+b\widetilde E_6\).
Put
\[
 C_2=\binom{2\alpha}{\alpha}\binom{2\beta}{\beta}/2,\qquad
 e=(C_2+1)/3,\qquad
 r_x=\lfloor\alpha/Q_x\rfloor/9,\quad r_y=\lfloor\beta/Q_y\rfloor/9 .
\]
Here \(C_2\) is an integer, \(C_2\equiv-1\bmod3\), and \(r_x,r_y\)
are integers by (10). Write \(J_D,J_R\) for the series
\[
 \begin{split}
 J_D&=\frac{d_xT}{1-T}+\frac{d_y bT}{1-bT},\\
 J_R&=r_x(T^3-T^6)+r_y(b^3T^3-b^6T^6).
 \end{split}
\]
Only their coefficients modulo three will be needed, but \(b\)
can be retained as the actual unit. Through order eight the
native expansion modulo \(I_2\) is
\[
 D R(T)=1+uT-u^2T^2
           +3uT(J_D+J_R)+3u^2T^2(e-2J_R).                          \tag{18}
\]

Here is the factorial computation for (18). For one insertion
count \(m=1,2\), put \(Q=3^h,\ n=\alpha\) in the \(x\)
coordinate or \(n=\beta\) in the \(y\) coordinate, and let
\(x'\le8\) be its axis count. The shifted binomial factor
\(\binom{Qx'+mn}{Qx'}\) can first be divided by
\(\binom{x'+m\lfloor n/Q\rfloor}{x'}\), a unit.
To compute the ratio modulo nine use the unit product
\(\prod_{1\le t\le 3q+d,\,3\nmid t}t\equiv
2^q\prod_{t=1}^d(3q+t)\bmod9\).
At every level below \(Q\) the factors \(2^q\) cancel.
All residual linear factors cancel except those at the level
immediately below \(Q\), since at the lower levels the
binomial's axis input has high quotient divisible by three.
At that last level the ratio is \(1+3x' H_{md}\), where
\(d\) is the digit immediately below \(Q\) and \(H_j\) is
the harmonic sum modulo three. It is \(1\) if \(Q=1\).
The upper binomial equals
\(1+3m r H_{\lfloor x'/3\rfloor}\bmod9\) by expanding its
factors \(1+9mr/t\) for \(1\le t\le x'\), where
\(r=\lfloor n/Q\rfloor/9\).
Since \(H_1=1,\ H_2=0\bmod3\), the lower digit produces
the term in \(J_D\) only when \(m=1\).
The upper binomial correction for \(m=1\) produces the
term in \(J_R\), since \(H_{\lfloor x'/3\rfloor}\) is one
exactly at \(x'=3,4,5\) within \(0\le x'\le8\), while for
\(m=2\) that correction doubles. Multiply the two coordinate
factors, use the pure-axis quotient modulo nine and the
foreign-only coefficients \(1,C_2\), and sum through order
eight. This gives (18). Counts \(m\ge3\) have valuation at
least two after multiplication by \(u^m\).

In the main constant-power block, the inverse-binomial
normalization corrections for indices \(\ell=3,4,5\)
multiply \((DR)_{8-\ell}\) of valuation at least \(4/3\),
so are discarded. At all other indices its normalization
factor is \(1\bmod9\). To use (18), the constant expansion
is \(H_0(T)(1+kT-k^2T^2)\) modulo \(3\mathcal O\), where
\(H_0=e_2(UT^3)\equiv D^3\bmod 3\mathcal O\).
The coefficient at eight from \(1+uT-u^2T^2\) is discarded:
at \(u=-k\) its valuation is \(8/3\) because it equals
\(-63k^8/8!\), and the perturbation \(u+k\) of valuation
at least \(4/3\) gives no term below two. In fact its linear
multiplier at \(u=-k\) is \(15k^7/7!\) of valuation \(4/3\),
and the quadratic multiplier is integral.
For the terms involving \(J_D,J_R\), set
\(A_D=[T^6]D^3 J_D,\ A_R=[T^6]D^3 J_R\) modulo
the maximal ideal. The \(T^7\) coefficients of \(D^3J_D\)
and \(D^3J_R\) are zero modulo \(3\mathcal O\):
the first series is a polynomial of degree at most six,
and the second is supported at multiples of three in
characteristic three. The remaining contribution at order
eight from (18) is therefore
\[
             3uk(A_D+A_R)+3u^2(-eU^2-2A_R)
         \equiv -3k^2(eU^2+A_D)\pmod{I_2}.
\]
Using \(\bar b=\bar i,\ \bar U^2=\bar i\) gives
\(A_D=\bar i(d_x+d_y)\) by the top degree coefficients of
\(D^3J_D\).

For the earlier constant-power blocks, only the immediately
preceding one can contribute below valuation two, and only
when \(3\nmid H\). At precision \(I_1\) the coefficients of
\(DR\) through degree seventeen are
\[
             1+uT-u^2T^2+u C_*T^{10}-2u^2 C_* T^{11},
           \qquad C_* =r_x+b^9r_y.
\]
Indeed at these indices \(x'\le17\) and the shifted binomial
at \(m=1\) or \(2\) modulo three is one for \(x'<9\) and
\(1+mr\) for \(x'\ge9\) by the first additional digit.
Terms using that extra digit in both coordinates require
at least eighteen axis choices. Thus only \(DR_{10}\)
and \(DR_{11}\) contribute to the earlier block at order
eight. They use constant-power indices \(9H-2\) and \(9H-3\).
The respective normalization scalars are \(3Hk/U\) and
\(3H/U\) to the required precision, by the preceding
calculation of the tail and the ratio \(k/(9H-2)\).
Their total is \(3H C_*(ku-2u^2)/U\), which vanishes modulo
\(I_2\) since \(u=-k\bmod I_{4/3}\). Hence the exact zeros
force \(e+d_x+d_y=0\bmod3\).

Finally let \(w\) be the total number of \(1\) digits and
\(J\) count adjacent \(11\) pairs across the two separate
expansions (with no pair between them).
The same unit product, applied to each doubling binomial,
gives \(C_2\equiv2^{w-1}(1+3J)\bmod9\).
Indeed a digit \(1\) with next digit \(q\bmod3\) contributes
\(2(1+3q)\), while a digit \(0\) contributes one.
With \(w\) even, this reads
\(e\equiv w-1-J=N_{\rm run}-1\bmod3\), giving (17).

### The adjacent digits must vanish

**Adjacent-digit restriction.** An actual block in the first
wild four-term case satisfies \(d_x=d_y=0\) and
\(N_{\rm run}\equiv1\bmod3\). In particular the foreign monomial
\(XY\) is excluded at \(v(k)=1/3\) for every pair of prime-power
axis exponents under the two-pole hypothesis.

More precisely, suppose (10), (11), (16) and (17) hold and keep
the variables of Section 8 with \(a=1\). Then
\[
 \widetilde E_8-(1+b)\widetilde E_7+b\widetilde E_6
       \equiv -9U(i d_x+d_y)\pmod{I_{7/3}}.                         \tag{19}
\]
The right-hand side has valuation two unless both digits
vanish, since \(1,\bar i\) are independent over \(\mathbb F_3\).
We compute (19) directly from the normalized moments.

First the pure-axis factorial quotient has the next
integer congruence. For \(Q=3^h\), \(h\ge1\), one has
\[
  \frac{(Qn)!}{n!(Q!)^n}
       \equiv 1+9\left\lfloor\frac{n+1}{3}\right\rfloor
                                                  \pmod{27}.      \tag{20}
\]
Indeed the \(3\)-free product through \(3n\) is congruent to
\(2^n(1+9\binom{n+1}{3})\bmod27\), by multiplying the factors
\(2+9j(j+1)\). Iterating the splitting and normalizing by
the fixed factorial gives a sum of differences
\(\binom{3^r n+1}{3}-n\binom{3^r+1}{3}\), \(0\le r<h\), in
the coefficient of \(9\). For \(r\ge1\) the difference
vanishes modulo three by the second ternary digit. The
remaining summand is \(\lfloor(n+1)/3\rfloor\bmod3\).
For \(Q=1\) the pure coefficient is one exactly.
The pure correction from (20) in \(DR(T)\) is therefore
\[
  9T^2\left(\frac{\mathbf 1_{Q_x>1}}{1-T^3}
      +\frac{b^2\mathbf 1_{Q_y>1}}{1-b^3T^3}\right)
                                                    \pmod{I_{7/3}}.
\]
Its degree-eight constant-shift contribution vanishes:
at this precision the factor \(9\) lets us multiply by
\(H_0\equiv D^3\bmod3\mathcal O\), leaving a polynomial of
degree at most five in characteristic three.

In the main block the \(m=1,2\) foreign insertions from
(18) already have full precision \(I_{7/3}\) on each of
their terms. In fact their integer coefficient errors in
the modulo-nine computation have valuation at least two,
while their scalar factors have valuations \(1/3,2/3\).
The contribution from \(m=3\) is irrelevant unless the
foreign exponents have altogether exactly two digits equal
to \(1\). Its own normalized factorial quotient has
valuation \(w-1\) by Legendre's formula, with \(w\) the
total ternary digit count. In the remaining case \(w=2\),
the coefficient after multiplication by \(u^3\) has value
at least two and the shifted binomials modulo three add
axis loads \(d_x,d_y\) from \(3\alpha,3\beta\). The
degree-eight coefficient after convolution is thus a
scalar times
\[
        [T^5](1-T)^{3-d_x}(1-bT)^{3-d_y}.
\]
For \(w=2\), (17) forces \(d_x+d_y=0\) when the two
digits are adjacent in one coordinate and \(d_x+d_y=2\)
otherwise. The displayed coefficient is zero in both
cases, by the exponents being both three or both two.
All \(m\ge4\) have value at least \(7/3\), using the
own factorial valuation bound \(v_3(m!)\) and \(v(u)=1/3\).

The inverse binomials for \(\ell=3,4,5\) in the main
block differ from one by a multiple of three, and multiply
\((DR)_{8-\ell}\) of valuation at least \(4/3\) by (15).
For \(\ell=6,7,8\) the differences have valuation at least
two, and the three retained scalar products have
valuation at least \(2/3\), even before these differences.
For \(\ell=0,1,2\) the differences also have valuation at
least two and multiply \((DR)_{8-\ell}\) of valuation at
least \(4/3\). Thus the main-block ratios may be replaced
by one in (19). For the earlier blocks, only the terms
from the immediately previous block can affect (19).
Even the next block back has cost at least two and
multiplies positive-order coefficients of \(DR\) of
valuation at least \(1/3\). For the previous block, the
native expansion through order seventeen used after
(18) is valid at precision \(I_{4/3}\): the errors in
the foreign coefficients modulo three gain at least
\(1/3\) from \(u\), and the pure error has value at least
two. Again only \((DR)_{10},(DR)_{11}\) contribute and
the sum \(3H C_*(ku-2u^2)/U\) has valuation at least \(7/3\),
since \(u=-k\bmod I_{4/3}\).

It remains to compute the next constant convolution in
(18). The base term \(1+uT-u^2T^2\) still contributes
zero at precision \(I_{7/3}\), since its coefficient at
\(u=-k\) and its change under \(u+k\) have valuations at
least \(8/3\). In the terms with factor \(3u\) or
\(3u^2\), replace the constant expansion by
\(D^3(1+kT-k^2T^2)\); the replacement error has
valuation at least one before these factors. Keep
\(A_D,A_R\) as in Section 8, with the unit \(b\)
retained. The degree-seven coefficients are divisible
by three as before. The degree-six contribution modulo
\(I_{7/3}\) is
\[
                       -3k^2(eU^2+A_D),\qquad
                 A_D=-b^3(d_x+d_y).
\]
The degree-six \(J_R\) terms cancel by the same
coefficient \(-3\) as previously. Moreover
\(U^2\equiv i,\ b^3\equiv -i\bmod3\mathcal O\) by the
full modulo-three equations, so \(eU^2+A_D\equiv0\bmod3\mathcal O\)
by (17). The degree-five contribution comes only from
the \(3u\) term with the factor \(-k^2T^2\), since
\(D^3(e-2J_R)\) has no degree-five coefficient modulo
three. In \(D^3(J_D+J_R)\) that coefficient is
\(-b^3d_x-b^2d_y\equiv i d_x+d_y\bmod3\mathcal O\).
Therefore the last contribution is
\(-3u k^2(i d_x+d_y)\equiv-9U(i d_x+d_y)\bmod I_{7/3}\),
as asserted.

For \(XY\), if one axis exponent is \(1\), its associated
quotient in (10) is \(1\), already excluded. Otherwise
\(N_{\rm run}=2\), so (17) would require both digits
to be \(1\), contrary to (19). This includes the two
cubic axes and uses no tameness or discreteness of
the coefficient value group.

### A complete exclusion for the foreign quartic on one axis

**First critical quartic exclusion.** Under the two-pole
hypothesis, the four-term support with foreign monomial
\(X^4\) or \(Y^4\) has its four-window at \(v(k)=1/3\) for
arbitrary prime-power axis exponents. The proof applies to
every starting index and every real-valued extension of
\(v_3\).

By symmetry take \(X^4\). Condition (10) rules out
\(Q_x=1,3\), and the adjacent-digit condition rules out
\(Q_x=9\). The only starting indices to consider are
\(9H+5\), \(H\ge0\). For the lift adjoin \(\pi\) with
\(\pi^3=3\), and use integral parameters \(A,B,C\) given by
\[
        k=\pi(1+i)+\pi^2 A,\qquad b=i+\pi B,\qquad
        u=-k+\pi^4 C .                                             \tag{21}
\]
Indeed \(v(b-i)\ge1/3,\ v(u+k)\ge4/3\) by Section 8.
Also \(U\equiv i-1\bmod3\mathcal O\) gives
\((k/\pi)^3\equiv(1+i)^3\bmod3\mathcal O\). If the two
cube bases differed by valuation strictly between zero
and \(1/3\), the cube of their difference would be
the unique least term below valuation one. Their
residues agree, so \(A\) is integral.

We use precision \(I_4\). Write \(E_i=\widetilde E_i\)
for \(5\le i\le8\) in the following calculation and
\(Z_i=E_i-(1+b)E_{i-1}+bE_{i-2}\) at \(i=7,8\).
Two uniform reductions make the arithmetic finite.

1. In each normalized moment only offsets
   \(-27\le\ell\le i\) from \(j_0=9H\) can contribute at
   precision \(I_4\). A block four or more nine-index
   steps back loses at least four valuation units, by
   (3). For nonnegative \(\ell\), its coefficient is
   \(k^\ell/\prod_{t=1}^{\ell}(9H+t)\). The ratio after
   changing \(H\) modulo \(27\) is \(1\bmod I_4\)
   because every \(1\le t\le8\) has valuation at most
   one. For a negative offset \(-d\) use
   \(k^{-d}\prod_{t=0}^{d-1}(9H-t)\), a product that is
   zero if the index falls below zero. Telescoping the
   difference when \(H\) changes modulo \(27\) gives
   valuation at least
   \[
                   3+\lceil d/3\rceil+\lceil d/9\rceil-d/3
                           \ \ge\ 4 .
   \]
   Indeed each original factor has its lower bound
   \(\mathbf1_{3\mid t}+\mathbf1_{9\mid t}\le2\), and a
   changed factor has at least five. The subtraction of
   \(d/3\) is from \(k^{-d}\). Thus \(H\bmod27\) suffices.

2. Use \(Q_x\in\{27,81,243\}\) and \(Q_y\in\{1,3\}\)
   as representatives. The axis quotient in (20)
   satisfies, for every power \(Q>1\), the stronger
   congruence \(1+9\binom{n+1}{3}\bmod81\). In its proof
   all cross terms of two factors divisible by nine
   still disappear. Each term for \(r\ge1\) in the
   iterated sum is \(3^{3r-1}(n^3-n)/2\), divisible
   by nine. This explains why \(Q_y=3\) represents
   every nontrivial power. For \(Q_x\ge243\), factor
   the extra \(x\)-binomial
   \(\binom{Q_x x+4m}{4m}=\prod_{j=1}^{4m}(1+Q_x x/j)\)
   away from the axis and foreign-only coefficients.
   For \(m=1,2\) it is \(1\bmod81\). For \(3\le m\le5\)
   it is \(1\bmod27\), and the terms already have
   valuation at least two from \(u^m\) and the native
   factorial quotient. For \(m\ge6\), the bound
   \(m/3+v_3(m!)\ge4\) discards the whole term. Thus
   \(Q_x=243\) represents all \(Q_x\ge243\) at this
   precision.

The [exact Python checker](verify_binary_wild_x4_first_radius.py)
expands the 162 representative systems in
\(\mathbb Z[i,\pi,A,B,C]/(i^2+1,\pi^3-3,\pi^{12})\).
It uses the \(3\)-free factorial parts and integer valuations
for each coefficient. For a negative offset it uses the
falling product above, with the unit \((1+i+\pi A)^{-d}\)
expanded by its finite binomial truncation. Thus the
formal substitutions are valid at integral parameters in
arbitrary valued extensions, including fractional
valuations within the coefficient ring. In the following
initial coefficients, \(\bar i\) and the other residues are
written without bars.

The first four initial polynomials are as follows, where
\(\varepsilon=0\) if \(Q_x=27\), and \(1\) otherwise:
\[
\begin{array}{c|c|l}
\text{form}&\pi\text{ order}&\text{initial coefficient}\\ \hline
 E_5&3&i B^3-i A^3+(1+i)H\\
 E_6&3&-B^3+(1-i)A^3\\
 Z_7&4&2+i-iC\\
 Z_8&8&iC^2+(2+i)C+\varepsilon-H .
\end{array}                                                        \tag{22}
\]
The first two equations have an invertible matrix in the
cubes, and the third is affine with unit slope. Thus a
block forces the residues \(A=(i-1)H,\ B=H,\ C=1+i\), and
then \(H\equiv\varepsilon+2\bmod3\). Furthermore the
degree-four row of \(Z_7\) gives \(C=1+i+\pi W\) with
\(W\) integral, by the actual \(\pi^5\) error bound.
Make this substitution for all remaining rows.

If \(Q_x=27\), the order-nine initial coefficient of \(Z_8\)
is independent of \(W\) and evaluates to \(-1-i\) at the
forced \(A,B\) for \(H\equiv2\). This is impossible. For
\(Q_x\ge81\) we have \(H\equiv0\), so \(A,B\) have positive
valuation. If \(Q_y=1\), the order-five row of \(Z_7\) is
\(i(A-W)\), hence \(W\) also has positive valuation. The
initial order-ten coefficient of
\(Z_8+i\pi^6E_5+2\pi^6E_6\) at zero residues is \(-i\),
again impossible.

If \(Q_x\ge81,\ Q_y>1\), put
\[
 G=Z_8+(2+i)\pi^6 E_5+iZ_7^2
        +\pi^7B\big((1+2i)E_5+iE_6\big)+(2+2i)\pi^6 Z_7 .
\]
The checker verifies its first order is eleven and its
coefficient there is
\(c_H+B+(2+2i)B^5+A^3B^2\), with \(c_H\) a constant.
If \(c_H\ne0\) the positive valuations of \(A,B\)
exclude the block. Otherwise the coefficient is \(B\)
times a unit, so \(v(B)\ge1/3\) by the order-twelve
error. At order four in \(E_5,E_6\) neither coefficient
depends on \(W\). Their constants are respectively zero
and \(i\). From \(E_5/\pi^3\) and \(v(B)\ge1/3\) we
obtain \(v(A^3)>1/3\), since the degree-four row has
zero residue at \(A=B=0\). Now \(E_6/\pi^3\) has its
order-four constant \(i\pi\) as the unique term of
valuation \(1/3\), a contradiction.

Run the finite verification from the repository root with
Python 3.12 via
\(\texttt{python3 research/factorial/image-rigidity/verify_binary_wild_x4_first_radius.py}\).
The representative reductions above are independent of
the computation and supply its all-start and all-axis scope.

### Bounded digit windows at the first radius

**Finite upper profiles at precision four.** For the
remaining first-radius four-term families satisfying
(10), (11), (17) and (19), the four normalized
moment rows modulo \(I_4\) depend on \(H\bmod27\)
and on finitely many realizable exponent profiles.
This includes unbounded upper quotients
\(\lfloor\alpha/Q_x\rfloor,\lfloor\beta/Q_y\rfloor\).
For each axis it suffices to retain five digits
immediately below it and seven immediately above,
together with the axis height truncated at five.
In addition the foreign-only factorial residues
\(F_m\bmod81\) are needed for \(m\le5\), with only
\(m=1,2\) needed at total digit count at least four.
The profile reduction leaves the exact exclusion of
the remaining first-radius rows open.

Indeed the scalar bounds in the quartic calculation
apply to all integral native coefficients with
\(v(k)=1/3\). At precision \(I_4\) they give
\(-27\le\ell\le8\) and \(H\bmod27\), so all axis
counts are at most \(35\). The foreign-only factorial
ratio for \(m\) insertions is
\[
 F_m=\frac{(m\alpha)!(m\beta)!}
           {m!(\alpha!\beta!)^m}.
\]
With total foreign digit count \(\sigma\),
Legendre's formula gives
\(v_3(F_m)\ge(\sigma-1)v_3(m!)\) as in Section 9.
Since \(v(u)=1/3\), counts \(m\ge6\) do not enter
\(I_4\). At \(\sigma\ge4\) the same holds for
\(3\le m\le5\). The pure axis factor uses only the
height-zero flag, by the modulus \(81\) identity in
the quartic calculation.

For a cross-binomial factor on either axis write
\(Q=3^q\), \(N=m\alpha\) for that axis and
\(d=\lfloor N/Q\rfloor\). The other coordinate uses
\(\beta\) in place of \(\alpha\). Separate its upper
binomial by repeatedly splitting each factorial
into its \(3\)-free part and its multiples of three.
For \(g(n)=\prod_{1\le k\le n,\,3\nmid k}k\) this gives
the exact factorization
\[
 \binom{Qx+N}{Qx}
   =\binom{x+d}{x}
      \prod_{j=0}^{q-1}
         \frac{g(x3^{q-j}+N_j)}{g(x3^{q-j})g(N_j)},
       \qquad N_j=\lfloor N/3^j\rfloor .
\]
The powers of three cancel in the \(q\) splitting
steps since \(Qx\) is divisible by \(3^q\).
Modulo \(81\) the relation \(g(n+81)=-g(n)\)
makes each product factor one when \(q-j\ge4\).
The remaining factors use \(N_j\bmod81\) for
\(q-3\le j<q\), omitting negative \(j\).
The upper binomial modulo \(81\) depends only on
\(d\bmod3^7\). To prove this, in Vandermonde's
identity for \(\binom{x+d+\delta}{x}\) take
\(3^7\mid\delta\). For \(1\le a\le x\le35\),
the factor \(\binom{\delta}{a}\) has valuation at
least \(7-v_3(a)\ge4\).

For \(m\le5\), multiplication by \(m\) in a word
of zeros and ones has carry at most two, and a
zero clears it. Along successive ones the carry
stabilizes after at most two digits. Thus two
preceding digits determine the incoming carry at
any position. The five digits below the axis give
the three positions \(q-3,q-2,q-1\) and the
preceding digits needed by their incoming carries.
The seven above it suffice for \(d\bmod3^7\),
and include the upper digits needed by the unit
product windows. Absent digits are zero.

For the \(F_m\) data, the finite transfer with
bounded carries described in the low-exponent
profile lemma of Section 9 works across the axis
boundary as well. Mark the boundary in each input
digit word, keep the last five lower digits and
the first seven upper digits, and retain the unit
residue product and valuation sum for each required
\(m\). The valuation sums can be saturated at
\(4+v_3(m!)\), since reaching that bound on either
axis already gives \(F_m=0\bmod81\).
One tracks the total digit count up to four and its
parity, and counts the runs modulo three, keeping
the previous digit across the boundary. The last
lower digit and first two upper digits enforce
(10) and (19). These checks also enforce (17)
using the run count. All states are finite, every
allowed word has a path, and every accepting path
gives a pair of words. The four rows at \(I_4\)
can thus be tested through an exact finite
reachability problem followed by the linked
moment calculation.

## 9. Low foreign exponents at larger critical radii

Suppose \(p=3,\ L=3^h\ge9\), and the foreign term of a
four-term phase is \(u X^\alpha Y^\beta/(\alpha!\beta!)\)
with \(v(u)>0,\ \alpha<Q_x,\ \beta<Q_y\). Assume the
foreign exponents have only digits \(0,1\) and an even
total digit count. This includes \(XY\) when both axis
powers are at least three. Keep the two unit axes and
the canceled leading row. In particular the ratio
\(r=b/a\) reduces to a primitive fourth root.
Normalize \(a=1\) and put \(\delta=v(k)=(1-1/L)/2\).

**Quadratic numerator filter.** A four-zero block may
only start at \(3LH+s\) with \(2L-1\le s\le3L-4\).
Let \(j=s+2-2L\), and apply the following filter to both
indices \(j,j+1\). An index \(t\) with \(s_3(t)\le2\)
cannot be divisible by three. At a unit index \(t\)
of digit sum at most two, one has \(v(u)=\delta\) and
\[
 v(u+k)\ \ge\ \delta+\begin{cases}
          1-v(k^t/t!) &t\equiv1\bmod3,\\
          (1-v(k^t/t!))/2 &t\equiv2\bmod3.
          \end{cases}                                              \tag{23}
\]

Indeed the native series modulo the full \(3\mathcal O\)
ideal is \((1+uT-u^2T^2)/D\). For one or two foreign
insertions all digits fit below both axes without
carries, since each digit is at most one. The
foreign-only second coefficient is \(-1\bmod3\)
by its even digit count, and counts of at least three
have integer coefficients divisible by three.
Use (6) and the low constant polynomial
\(K(T)=\sum_{t=0}^{L-1}k^tT^t/t!\). Its effective degree
in \(\mathcal O/3\mathcal O\) is \(2L/3\): positive
coefficients survive exactly at digit sum at most two.
The degree of \((1+uT-u^2T^2)K(T)\) is therefore
less than \(L\). Multiplying by \(e_2(UT^L)\) shows
that the only contributions at \(2L+j,2L+j+1\) come
from its highest term. These two numerator coefficients
must vanish modulo \(3\mathcal O\), by the degree at most
\(s+1\) of the actual prefix times \(D\).
After removing \(U^2/2\), the condition for \(t=j,j+1\)
is
\[
   \frac{k^t}{t!}\big(1+t z-t(t-1)z^2\big)\equiv0
                  \pmod{3\mathcal O},\qquad z=u/k,                   \tag{24}
\]
where fair interpretation of (24) at \(t=1\) simply
omits the quadratic term.

If \(3\mid t\), the linear and quadratic terms in
brackets have positive valuation since \(v(z)>-\delta>-1/2\).
Thus a small digit sum would contradict (24). At a unit
index of small digit sum, \(v(z)>0\) leaves the constant
uniquely least. For \(v(z)<0\), the linear term is
uniquely least when \(t\equiv1\), and the quadratic
term uniquely least when \(t\equiv2\). For the first
assertion use \(1+2v(z)>v(z)\), again by \(v(z)>-1/2\).
Hence \(z\) has value zero and its residue is \(-1\),
the unique root in either case. With \(z=-1+y\), the
constant term \(1-t^2\) has valuation at least one.
The linear coefficient \(t(2t-1)\) is a unit at
\(t\equiv1\), while at \(t\equiv2\) it is divisible
by three and the quadratic coefficient is a unit.
The least-term comparison in this expansion gives
(23).

**Second critical radius.** For \(L=9\), any such
exact four-zero block must start in one of the
following three classes, with the indicated
necessary conditions:
\[
\begin{array}{c|c|c}
s\bmod27 & \text{constraint on }u,k & v(r^2+1)\\ \hline
17 & v(u+k)\ge1 & 1/27\\
20 & v(u+k)\ge2/3 & 2/27\\
23 & v(u)\ge1/3 & \ge1/9 .
\end{array}                                                        \tag{25}
\]
In the first two rows (23) also gives \(v(u)=4/9\).
To see there are only these starts, the candidate
values of \(j\) are \(1,\ldots,7\), and \(3,6\) are the
positive multiples of three of digit sum at most two
among the test indices. Their exclusion leaves
\(j=1,4,7\). At the first two, (23) yields the
displayed bounds on \(u+k\). At \(j=7\) the only
potential term below valuation one in the direct
coefficient
\(K_7+uK_6-u^2K_5\) is \(uK_6\), where
\(K_t=k^t/t!\) and \(v(K_6)=2/3\). This proves the
third-row bound.

For the pole constraints, put \(H_0=e_2(UT^9)\) and
\(S=(1+uT-u^2T^2)K(T)\), and work in
\(\mathcal O/3\mathcal O\). At \(s=17\) the bound
\(v(u+k)\ge1\) gives \(S=1+P T^3+Q T^6\), where
\(P=-8k^3/6,\ Q=-35 k^6/6!\) have valuations
\(1/3,2/3\). Indeed at \(u=-k\) the degree-\(t\)
coefficient is \((1-t^2)k^t/t!\). Agreement of
\(D\) times the actual prefix with \(H_0 S\) through
degree twenty then gives the exact divisibility
\[
                D\mid S(1+UT^9)+(U^2/2)T^{18}
                       \quad\text{in }\mathcal A[T].
\]
Let \(z_a=U,\ z_b=U/b^9\), and let \(\xi_a,\xi_b\)
be the corresponding distinct roots \(-1\pm i\)
of \(e_2\) in a valued extension. Evaluation at
the two poles gives
\(e_2(z_h)=-(1+z_h)(S(h^{-1})-1)\)
modulo \(3\mathcal O\) for \(h=1,b\). Consequently
the errors \(\Delta_h=z_h-\xi_h\) have valuation
\(1/3\) and the residues of \(\Delta_h/P\) are
\(-\bar h^{-3}\), using \(e_2'(\xi_h)=1+\xi_h\).
Since \(\zeta=\xi_a/\xi_b\) is the fourth root
congruent to \(b\), the residue of
\((\Delta_a-\zeta\Delta_b)/P\) equals
\(-1+\bar b^{-2}=-2\ne0\). Thus \(v(b^9-\zeta)=1/3\).

At \(s=20\), write \(w=u+k\) so \(v(w)\ge2/3\).
Expanding \(S\) gives \(1+wT+P T^3+Q T^6\) modulo
\(3\mathcal O\). In fact the \(w\)-correction beyond
degree one is \(wK_{t-1}+(2kw-w^2)K_{t-2}\), of
valuation at least one for \(3\le t\le8\); in
degree two combine it as \(3kw-w^2\). The last
omitted term through degree twenty-four of \(H_0S\)
is \((U^2/2)Q T^{24}\). The actual zero block forces
\[
                     D\mid H_0S-(U^2/2)Q T^{24}.
\]
Since \(S\) reduces to one at each of the two poles
and \(\bar b^{24}=1\), the two evaluations of
\(e_2\) at the \(z_h\) have value \(2/3\) with ratio
of residue one. The opposite unit derivatives of
\(e_2\) at the two roots imply
\(\overline{\Delta_a/\Delta_b}=-1\).
It follows that \(v(b^9-\zeta)=2/3\).

At \(s=23\), the product \(H_0S\) has degree at most
twenty-six in \(\mathcal A[T]\), and the same
comparison through degree twenty-six gives
\(D\mid H_0S\), now with no term removed. Both evaluations
of \(S\) are units, so \(v(b^9-\zeta)\ge1\).
In each of the three cases, if \(0<v(b-\zeta)<1/9\)
the ninth power of \(b-\zeta\) is the unique least
term in \(b^9-\zeta\), and otherwise
\(v(b^9-\zeta)\ge1\). Since \(b+\zeta\) is a unit,
this proves the three claims on \(v(r^2+1)\) in
(25). They hold in the original coefficient field,
even though the roots of \(e_2\) were taken in a
valued extension.

**Late second-radius coefficient alignment.** In the
\(s=23\) row of (25) one has the stronger bounds
\[
                        v(u)=4/9,\qquad v(u+k)\ge10/9.                \tag{26}
\]
This lift applies to all the low even-digit foreign
monomials of Section 9.

To prove it, normalize the four moments by
\(k^{27H}/(27H)!\) and use the degree-twenty-five
combination \(D\widetilde E\). Write \(K_\ell=k^\ell/\ell!\)
for constant-power indices in its main block. Its native
series, to precision \(I_{16/9}\) once \(v(u)=4/9\), is
\[
       DR(T)=1+uT-u^2T^2+
          3uT\left(\frac{d_xT}{1-T}+\frac{d_y bT}{1-bT}\right).
\]
This follows from the modulo-nine single-insertion
calculation for (18), now with both upper quotients zero.
The pure error has value at least two, the double
insertion error has value at least \(1+2v(u)=17/9\),
and counts of at least three have value at least
\(1+3v(u)\). Initially use only \(v(u)\ge1/3\) from
(25). At precision \(I_{4/3}\) the same formula without
the \(3u\) term is valid.

Every constant-power block before the main one loses
at least one valuation unit, and at degree twenty-five
it multiplies a native coefficient of \(DR\) at an
index greater than two. Such a coefficient has value
at least \(1+v(u)\). Hence all earlier blocks may be
discarded even at precision \(I_{16/9}\). The main
inverse-binomial factors are \(1\bmod3\). In the three
retained products from \(1+uT-u^2T^2\) their errors
have value at least two, since \(v(K_{24})=2/3\),
\(v(K_{25})=10/9,\ v(K_{23})=11/9\) and
\(v(u)\ge1/3\). Thus at precision \(I_{4/3}\) the
zero condition reads
\((u+k/25)K_{24}=0\bmod I_{4/3}\). This gives
\(v(u)=4/9\).

We can now retain the \(3u\) term, with value \(13/9\).
It sees only the constant-power indices \(0,9,18\)
at precision \(I_{16/9}\), because every other
index costs at least \(1/3\). Their residue polynomial
is \(e_2(UT^9)\equiv D^9\bmod3\) at the canceled
row. Its product with the two displayed simple
fractions has degree at most eighteen over the
residue field, so its degree-twenty-four coefficient
vanishes. Thus the same zero condition for
\((u+k/25)K_{24}\) holds at precision \(I_{16/9}\).
Subtracting \(k/25-k\) of valuation \(13/9\) proves
(26).

### Realization of the three second-radius congruences

**Second-radius residue realizations.** Each of the three
classes in (25) occurs as a simultaneous zero of the four
normalized equations modulo \(3\mathcal O\), for any
\(H\ge0\), even with cubic axes and foreign \(XY\).
Thus the full modulo-three test for this slice has an
explicit lift in all three classes.

Take \(u=-k,\ L=9\) and let \(\psi=k^3/3\), so
\(\psi^3=3U\). The exact polynomial
\[
 S_0(T)=1-4\psi T^3-\frac7{16}\psi^2T^6
\]
represents \(S\bmod3\mathcal O\). Use the truncated
polynomials
\[
 \begin{split}
 F_{17}(T)&=S_0(T)(1+UT^9)+(U^2/2)T^{18},\\
 F_{20}(T)&=S_0(T)(1+UT^9)
               +(U^2/2)T^{18}(1-4\psi T^3).
 \end{split}
\]
For either \(s=17,20\), choose \(\psi\) of value \(1/3\)
so that \(F_s(1)=0\) with \(U=\psi^3/3\).
After substituting \(\psi=\pi Z,\ \pi^3=3\), the resulting
equation for \(Z\) has integral coefficients and reduces
to \(e_2(Z^3)=0\). Choose a unit root lifting either
residue satisfying this equation. To justify the lifts
even for multiple residue roots, factor an integral
polynomial of nonzero reduction over a valued algebraic
closure. Normalize each factor with a negative-valuation
root by its constant term. Such factors have constant
nonzero reduction, and each other factor reduces to
\(T\) minus the residue of its root. The product can
be scaled to the original polynomial by a unit, since
the product of the reduced factors is nonzero. Thus
every root of the reduction is the residue of a root.

The reduction of \(F_s(T)\) is \(e_2(\bar U T^9)\).
Besides \(T=1\) it has a distinct finite residue root.
Choose a root \(1/b\) lifting that residue by the
same argument, and choose \(k^3=3\psi\). Then
\(D=(1-T)(1-bT)\) divides \(F_s(T)\) exactly,
and \(F_s/D\) has degree at most \(s-1\).
Through index \(s+3\), \(F_s\) agrees with \(e_2(UT^9)S\)
modulo three, so the four low-prefix coefficients
of indices \(s,\ldots,s+3\) vanish. Equation (6) in its
congruence form transfers this test to every \(3LH+s\).

For \(s=23\), take \(U\) to be a root of \(e_2\) and
\(b^9\) the ratio of that root to the other root.
Set \(k^9=81U,\ u=-k\). Then \(D\) divides \(e_2(UT^9)\)
exactly and \(e_2(UT^9)S_0/D\) has degree twenty-two.
This realizes the third sequence of congruences.
All three families concern the normalized residue
tests; the exact moment equations retain higher
valuation terms.

### The foreign product at the second radius

**Second critical \(XY\) exclusion.** With two unit
prime-power axes \(Q_x,Q_y\ge3\), a foreign term \(uXY\)
of positive valuation, and distinct residue poles,
there are no four consecutive zero moments at
\(v(k)=4/9\). This holds for all starting indices and
arbitrary real-valued extensions of \(v_3\).

The uncanceled rows are already excluded in Section 3.
For the canceled row take \(a=1\), order the axes as
\(Q_x\le Q_y\), and use (25) and (26). Let \(i^2=-1\)
reduce to \(\bar b\). In the proof of (25) the root
\(\xi_a\) is then \(-1-i\), since its ratio to the
other root is \(i\). The same proof gives
\(v(U-\xi_a)\ge1/3,2/3,1\) for \(s=17,20,23\),
respectively. Adjoin \(\pi\) with \(\pi^e=3\) and
write
\[
 k=\pi^{4e/9}(-1-i+\pi^j A),\qquad
 b=i+\pi^j B,\qquad u=-k+\pi^\rho C,
 \qquad
 \begin{array}{c|rrr}
 s&e&j&\rho\\\hline
 17&27&1&27\\
 20&27&2&18\\
 23&9&1&10
 \end{array}                                                        \tag{27}
\]
with \(A,B,C\) integral. For \(A\), note that
\((-1-i)^9=16(-1-i)\equiv\xi_a\bmod3\mathcal O\).
The residues of \(k/\pi^{4e/9}\) and \(-1-i\)
agree by injectivity of the ninth power in the
residue field. If the bases of the ninth powers
differ by a value strictly between zero and \(1/9\),
the ninth power of their difference is the unique
least term and has valuation less than one. This
gives the bound \(j/e\). The parameters \(B,C\)
are integral by (25) and (26).

Set \(E_j=\widetilde E_j\) after dividing by
\(k^{27H}/(27H)!\), and
\(Z_j=E_j-(1+b)E_{j-1}+bE_{j-2}\). The exact block
forces the four forms \(E_s,E_{s+1},Z_{s+2},Z_{s+3}\)
to vanish. The following reductions justify the
finite calculation of their coefficients.

1. Let the working precision be \(I_{N/e}\), with
   \(N/e\le4\), and write \(h=\lceil N/e\rceil\).
   Formula (3) shows that each block of twenty-seven
   indices before \(27H\) costs at least one further
   valuation unit. Thus the only offsets \(\ell\)
   from \(27H\) needed here satisfy
   \(-27(h-1)\le\ell\le s+3\).
   The scalars are invariant modulo \(I_h\) when
   \(H\) changes modulo \(3^{h-1}\). For \(\ell\ge0\),
   all denominators \(27H+t\), \(1\le t\le\ell\le26\),
   have value at most two, and the scalar is integral.
   For \(\ell=-d<0\), use
   \(k^{-d}\prod_{t=0}^{d-1}(27H-t)\), which also gives
   zero when the previous index is negative. A changed
   factor has difference of value at least \(h+2\).
   Each original factor has lower bound
   \(\mathbf1_{3\mid t}+\mathbf1_{9\mid t}+\mathbf1_{27\mid t}\).
   Removing one of these lower bounds and dividing by
   \(k^d\) gives for the telescoped difference
   \[
   h-1+\lceil d/3\rceil+\lceil d/9\rceil
       +\lceil d/27\rceil-4d/9\ \ge\ h .
   \]

2. At precision \(I_h\), use representatives
   \(Q\in\{3,9,27\}\) if \(h\le3\) and
   \(Q\in\{3,9,27,81\}\) if \(h=4\). Indeed both pure
   axis quotients have the common reduction modulo
   \(81\) proved after (21). With \(m\) insertions
   of \(XY\), the additional factors are \(m!\) and
   the two binomials \(\binom{Qx+m}{m}\) for the
   corresponding axis counts \(x\). Since \(v(u)=4/9\),
   counts \(m\ge6\) cost more than four. For \(m=1,2\)
   the product \(\prod_{j=1}^m(1+Qx/j)\) is
   \(1\bmod3^h\) as soon as \(Q\ge3^h\).
   For \(3\le m\le5\) it is \(1\bmod 3^{h-1}\)
   and \(m!\,u^m\) already has valuation at least one.
   This proves stabilization at the largest representative.

The [checker](verify_binary_xy_second_radius.py) works
in \(\mathbb Z[i,\pi,A,B,C]/(i^2+1,\pi^e-3,\pi^N)\).
It computes factorial ratios by valuations and unit
parts. For \(\ell<0\) it expands the inverse of the
unit \(-1-i+\pi^jA\) by the finite binomial
truncation. The first two rows have an invertible
initial matrix in \(A^9,B^9\). At \(s=17,20\) their
unique residue centers \((a_0,b_0)\) are
\((-1+i,-1)\), \((i,1)\), respectively. At \(s=23\)
with \(Q_x=3\) and \(3\mid H\), both centers are zero.
For these cases take \(t^3=\pi\), with \(v_t(t)=1\).
The checker subtracts the centers and verifies that
the first two rows can be combined to start with
unit multiples of \((A-a_0)^9,(B-b_0)^9\), both at
the same coefficient order \(p\). It checks, at
all remaining coefficient orders \(\ell<p+9r\),
that each monomial of total degree \(d\) in the two
centered variables has \(\ell+rd\ge p+9r\).
Here \(r=1,2,1\) for \(s=17,20,23\), respectively.
If the least of the two centered valuations were
\(0<z<r\), the row for the variable attaining that
minimum would have its ninth power uniquely least:
for \(d<9\) the checked inequality gives
\(\ell+dz>p+9z\), and for \(d\ge9\) this follows from
\(\ell>p\). The minimum is therefore at least \(r\).
Make the substitutions \(A=a_0+t^r A'\),
\(B=b_0+t^r B'\), with integral \(A',B'\).
The other cases at \(s=23\) need no substitution.

To reduce \(Z_{s+3}\), the checker combines the
first two rows to give leading monomials \(A'^9,B'^9\),
each with a unit coefficient and possibly a constant
term. For the unrecentered cases these are \(A^9,B^9\).
The initial coefficient of \(Z_{s+2}\) is a unit
multiple of \(C\), again allowing a constant.
In increasing uniformizer order the checker cancels
monomials divisible by one of these three leading
monomials, by subtracting the corresponding monomial
multiple of its row. At a fixed order this decreases
degree, since only the constant term of the pivot
row can contribute another term there. All replacements
are congruences modulo the ideal generated by the
first three forms and the working precision. The
first nonzero coefficient left is a nonzero constant
in every system:
\[
\begin{array}{c|l|r|c}
s & \text{axis condition} & \text{valuation} & \text{residue}\\\hline
17 & Q_x=Q_y=3 & 17/9 & 1\\
17 & Q_x=3<Q_y & 7/3 & -1-i\\
17 & Q_x\ge9 & 17/9 & -1\\
20 & Q_x=Q_y=3 & 20/9 & -1+i\\
20 & Q_x=3<Q_y & 20/9 & 1-i\\
20 & Q_x=Q_y=9 & 29/9 & 1-i\\
20 & Q_x=9<Q_y & 11/3 & 1\\
20 & Q_x\ge27 & 29/9 & -1+i\\
23 & Q_x=Q_y=3 & 3 & i\\
23 & Q_x=3<Q_y & 23/9 & i\\
23 & Q_x\ge9 & 23/9 & -i
\end{array}                                                        \tag{28}
\]
Each constant has value below its working precision,
so the four forms cannot vanish. The arithmetic is
over formal integral parameters throughout, including
when the coefficient field has fractional valuations
between successive powers of the uniformizer. The
291 systems use only the Python 3.12 standard library.
Run the checker from the repository root with
\(\texttt{python3 research/factorial/image-rigidity/verify_binary_xy_second_radius.py}\).

### Finite native data for the remaining low foreign terms

**Low-digit native profiles at precision four.** Suppose
\(Q_x=3^{q_x}, Q_y=3^{q_y}\), with foreign exponents
\(\alpha<Q_x,\beta<Q_y\) consisting of digits \(0,1\)
and positive even total digit count. If \(v(u)\ge4/9\),
the whole native moment series modulo \(I_4\) is
determined by the unit axis coefficients, \(u\), and
a finite set of exponent profiles. The realizable
profiles can be enumerated by finite-state transfers
on the two digit words. Thus at \(L=9\) the all-axis
reduction for \(XY\) can be replaced by a finite
profile reduction for this low even-digit sector.
The exact four-row lift for the other profiles remains
open.

For axis counts \(x,y\) and foreign count \(m\), the
native coefficient of \(b^y u^m\) factors as the two
pure axis quotients times
\[
  F_m\binom{Q_x x+m\alpha}{m\alpha}
       \binom{Q_y y+m\beta}{m\beta},\qquad
  F_m=\frac{(m\alpha)!\,(m\beta)!}{m!\,(\alpha!\beta!)^m}.
\]
Here \(a=1\) as before. Write
\(\sigma=s_3(\alpha)+s_3(\beta)\). The digit bound
\(s_3(mn)\le s_3(m)s_3(n)\) follows by expanding the
product in base three and carrying. Hence
\[
             v_3(F_m)\ge(\sigma-1)v_3(m!).
\]
In particular \(m\ge6\) gives value above four
after multiplication by \(u^m\). If \(\sigma\ge4\),
then \(3\le m\le5\) also gives value above four,
since \(v_3(m!)=1\) in that range. Counts \(m=0,1\)
have \(F_m=1\); the other counts can use \(F_m\bmod81\).

For either exponent \(a<3^q\) of this kind, let
\(q_*=\min(q,3)\), \(d=q-q_*\), and
\(b_m=\lfloor ma/3^d\rfloor\). For \(1\le m\le5\),
\(ma<3\cdot3^q\), so each factor \(1+3^q x/j\),
\(1\le j\le ma\), is integral. Factors with \(3^d\nmid j\)
are \(1\bmod81\): if \(d>0\), then
\(q-v_3(j)\ge q-d+1=4\). The remaining factors
give
\[
             \binom{3^q x+ma}{ma}
          \equiv \binom{3^{q_*}x+b_m}{b_m}\pmod {81}.
\]
The pure quotient depends only on whether \(q=0\),
by the identity following (21); the remaining axis
data are \(q_*\) and the five \(b_m\), each bounded
independently of \(q\).

For exact reachability of the \(F_m\) data, read
each axis word from the least significant end.
At position \(j\), put \(n_j=\lfloor a/3^j\rfloor\),
\(c_{m,j}=\lfloor m(a\bmod3^j)/3^j\rfloor\); then
\(\lfloor ma/3^j\rfloor=mn_j+c_{m,j}\).
If the next digit is \(d_j\), the carry updates by
\(c_{m,j+1}=\lfloor(md_j+c_{m,j})/3\rfloor\).
Legendre's formula gives the valuation of
\((ma)!/(a!)^m\) as \(\sum_{j\ge1}c_{m,j}\).
For a fixed \(m\) this sum need only be stored up
to \(4+v_3(m!)\), since reaching that threshold
already makes \(F_m=0\bmod81\).

Let \(g(n)=\prod_{1\le k\le n,\ 3\nmid k} k\bmod81\).
The product of the units modulo \(81\) is \(-1\),
so \(g(n+81)=-g(n)\). The unit of the same factorial
ratio is the product over \(j\ge0\) of
\(g(mn_j+c_{m,j})/g(n_j)^m\). Each factor depends
only on the carry and \(n_j\bmod81\), since changing
\(n_j\) by \(81\) cancels the two signs. Four digits
at a time give \(n_j\bmod81\), so the product is
updated with a bounded delay. At the end, append
zeros to flush the last digits and carries. The
\(b_m\) can be read from the top \(q_*\) digits and
the carry at their lower edge, obtained by retaining
the last three digits and their incoming carries.
Length truncated at three, digit count truncated at
four, and digit-count parity also have finite state.
The two resulting axis states give the \(F_m\) by
combining the valuation sums and dividing the
unit products by the unit of \(m!\).
Every finite path that ends at a word boundary
constructs an admissible axis word, and all admissible
words occur. Keeping the positive even digit counts
therefore gives exactly the attainable profiles.

### The second-radius lift for digit counts at least four

**High digit-count second-radius exclusion.** In the
canceled \(L=9\) row of (25), if the total number of
ones in the two foreign digit words is even and at
least four, the four-window holds at every starting
index. Let \(d_x,d_y\) be the leading digits of the two
words, including leading zeros to give lengths \(q_x,q_y\).
Let \(e_x,e_y\) be their next digits. Use zero for an
absent digit. With
\[
                  F_2=\frac12\binom{2\alpha}{\alpha}
                                  \binom{2\beta}{\beta},
\]
the precision-four test leaves the following
conditions on an exact four-zero block:
\[
\begin{array}{c|l}
s & \text{necessary conditions at precision } I_4\\ \hline
17 & \text{excluded}\\
20 & d_x=d_y=e_x=e_y=0,\quad F_2\equiv11\pmod {27}\\
23 & d_x=d_y=0,\quad F_2\equiv17-9(e_x+e_y)\pmod {27}.
\end{array}                                                        \tag{29}
\]
For each exponent profile meeting the conditions in the
last two rows, and for every \(H\), the four normalized
moment congruences at that \(s\) admit a simultaneous
solution modulo \(I_4\) in a finite valued extension
of \(\mathbb Q_3\). The first three moments can be
made exactly zero in such a solution. The start \(s=20\)
is excluded at the higher precision \(I_{38/9}\).
At \(s=23\), higher precision first forces
\(e_x=e_y=0\) and \(Q_x,Q_y>1\), then gives
incompatible residues in the linked moment rows.

The possible exponent data at either precision
\(I_h\), \(h=3,4\), simplify because only one and
two foreign insertions survive. Their binomial factors
use \(b_m=m b_1\) when the leading \(h-1\) digits are
retained, since \(m\le2\) and there is no carry.
An axis \(Q=1\) has a separate type with no foreign
digits. For other axes, use length \(h-1\), padding
shorter words with zeros at the low end. At precision
\(I_4\) keep \(Q=3\) with foreign digit \(1\) as a
separate type. To check the padding rule, the
extra factors after scaling a binomial by three
are \(1+3^{q+1}x/j\) for \(3\nmid j\) through a
multiple of three. At the padding steps from \(q=1\)
to \(2\) at precision three and from \(q=2\) to \(3\)
at precision four, their linear terms sum to zero
at the required precision because each full block
has reciprocal sum zero modulo three. The quadratic
terms are already beyond that precision. The zero
foreign word on \(Q=3\) requires no separate type.

For any pair of the resulting axis types with at least
one ordinary prefix type and for any \(F_2\equiv-1\bmod3\),
the data occur modulo \(3^h\). Indeed take an arbitrarily
long axis with that prefix and append isolated ones
well below the prefix and at distance at least four
from one another. The product for the unit factorial
ratio of doubling a word depends on windows of four
digits, by the \(g\)-formula above. An isolated one
in a position at least one with these separations
contributes
\(\binom{2\cdot3^j}{3^j}\equiv20\bmod81\)
for \(j\ge1\). This follows by multiplying the
four \(g(2\cdot3^a)/g(3^a)^2\) values for
\(a=0,1,2,3\), which are \(2,10,1,1\) modulo
\(81\), using just the available values when \(j<3\).
The residue \(20\) generates the units modulo
\(81\) since \(20^{18}=28\) and \(20^{27}=-1\)
there. Choose the number of ones to obtain the
target residue of \(F_2\), adding a full period to
ensure at least four ones. The combined digit count
is even since the doubled factorial ratios each
reduce to \((-1)\) to the power of their digit count.
Pairs with both axes of the exceptional types cannot
have enough ones. Interchanging the axes and
renormalizing the first coefficient to one does not
change \(F_2\) or the start. There are therefore
\(126\) unordered axis-pair and \(F_2\) profiles at
precision \(I_3\), and \(1404\) at precision \(I_4\).

For the linked rows use (27) and the normalized
scalar reduction of the \(XY\) calculation. In
fact \(H\bmod9\) suffices even at precision \(I_4\)
for this calculation. If \(H\) changes by a
multiple of nine, every scalar changes by an element
of \(I_3\). For the nonnegative offsets \(\ell<9\)
its change belongs to \(I_4\) owing to the valuations
of the denominators. For \(18\le\ell\le26\) the
only denominators of valuation two are at \(9,18\).
The two reciprocal residues after multiplication by
nine sum to zero modulo three, so their first-order
changes cancel at precision \(I_4\). For the middle
offsets \(9\le\ell\le17\), the second-order
combination \(Z_j\) uses coefficients of \(DR\) of
index at least two, as \(j\ge19\). At index above
two those coefficients belong to \(I_1\). The only
index two case is \(\ell=17,j=19\), where the
scalar itself has valuation \(14/9\), so its change
has valuation at least \(3+14/9\).
The negative offsets also multiply indices above
two in \(DR\). Thus \(Z_{s+2},Z_{s+3}\) are invariant
modulo \(I_4\). For \(E_s,E_{s+1}\) the parameter
bounds established by their first orders use
precision below \(I_3\). Each multiplier used to
cancel a term with one of these two rows in the
final reduction has value at least one, as checked
in the calculation. Their changes in \(I_3\) therefore
have no effect on the reducing identity modulo \(I_4\).

The [checker](verify_binary_low_digit_second_radius.py)
uses the axis factors just described and the pure
axis quotient modulo \(81\). It first performs the
center and valuation checks from the \(XY\) calculation
and then the same reduction of the final row. At
precision \(I_3\) the only zero remainders satisfy
\[
\begin{array}{c|l}
s & \text{conditions}\\ \hline
17 & d_x=d_y=0,\quad F_2\equiv23-9(e_x+e_y)\pmod {27}\\
20 & d_x=d_y=0,\quad F_2\equiv2\pmod9\\
23 & F_2\equiv8-3(d_x+d_y)\pmod9.
\end{array}                                                        \tag{30}
\]
The computation at \(I_4\) needs only the refinements
of (30). It gives zero remainders exactly in (29);
all other remainders have a unit constant initial
coefficient. This proves necessity for exact
four-zero blocks. To replay both stages, run
\(\texttt{python3 research/factorial/image-rigidity/verify_binary_low_digit_second_radius.py}\)
from the repository root.

To prove the realization statement, use the exact
first three forms at the chosen \(H\), with the
parameters and refinements checked in the reduction.
Over the complete discrete valuation ring generated
by the uniformizer and \(i\), clear denominators
that are units as formal series, using powers of
the unit part of \(k\). Divide each form by its
verified common uniformizer factor. Row combinations
and multiplication by unit constants give three
polynomials with reductions
\[
        A'^9-a_0,\qquad B'^9-b_0,\qquad C-c_0 .
\]
Here the primes on \(A',B'\) also denote the unchanged
parameters in the cases with no refined scaling.
These monic polynomials in distinct variables form
a regular sequence. The quotient by the three
lifted polynomials is torsion-free over the discrete
valuation ring. If a uniformizer multiple is in
their ideal, its representation
reduces to a syzygy of the three initial polynomials.
These syzygies are generated by the pairwise
relations, by induction using that each successive
monic polynomial is a nonzerodivisor modulo the
earlier ones. Lift and remove the pairwise relations
and divide the remaining representation by the
uniformizer. The quotient's uniformizer completion
is still torsion-free. It is finite over the base
ring, since the \(81\) residue basis monomials
generate by successive uniformizer approximation.
Inverting the uniformizer therefore gives a nonzero
finite algebra over the base field. A maximal
ideal gives a point over a finite field extension
with all three parameters integral. It solves the
first three forms exactly. The fourth form is in
their ideal modulo \(I_4\) by the zero remainder
identity. This proves the claimed simultaneous
congruence for the actual exponent word, with no
limit in its exponents.

One further precision still gives a finite test here.
At precision \(I_5\), the counts \(m=3,4\) may enter
for total digit count four. If the digit count is
at least six, only \(m=1,2\) contribute as before.
For \(I_5\) the pure axis quotient is independent of
the positive axis height. To see this, decompose it
into successive units at heights \(j\). The factor
for a single block at height \(j\) is
\[
        \prod_{\substack{1\le r\le 3^j\\3\nmid r}}
                        (1+3^j x/r).
\]
For \(j\ge2\), the sum of \(r^{-2}\) is divisible by
three by counting the units. Pairing \(r\) with
\(3^j-r\) shows that the sum of \(r^{-1}\) is
divisible by \(3^{j+1}\), since the number of
pairs is \(3^{j-1}\) and each product
\(r(3^j-r)\) reduces to \(-1\).
The linear and quadratic terms therefore vanish
modulo \(3^{2j+1}\), and the terms of higher
degree already have valuation at least \(3j\).
The block is \(1\bmod243\).
The axis interaction reduction above now retains
four leading digits instead of three, with the
bounded carries for \(m=3,4\) retained when needed.
Together with the same unit product modulo \(243\),
this gives the inputs for a lift above precision four.

For \(s=20\) it suffices to use \(I_{114/27}\).
At this precision the third and later foreign
insertions are still invisible. For the first two
insertions the new coefficient errors at modulus
\(81\) have value at least \(4+4/9>114/27\).
Only the pure axis factors need the higher modulus
\(243\). By (29) there are five unordered axis
pairs represented by \((0,0),(3,0),(3,1)\) as
height-prefix types, excluding \((0,0)\) paired
with itself. The three possible residues of \(F_2\)
modulo \(81\) are \(11,38,65\).
When \(H\) changes by a multiple of \(27\), the
normalized scalar differences gain one valuation
unit over the preceding \(H\bmod9\) comparison.
The first two forms are invariant modulo \(I_4\),
and the two \(Z\) forms are invariant modulo \(I_5\).
The final reduction again uses first-row multipliers
of value at least one. Thus \(H\bmod27\) suffices.
In all \(405\) representatives the initial
remainder is constant and nonzero at a valuation
at most \(113/27\). This excludes the middle start.

At \(s=23\), take \(I_{37/9}\) if either second
digit is one. The same modulus and omission bounds
apply as for \(s=20\). There are nine eligible
axis pairs among the five types \((0,0),(3,0),(3,1),
(3,3),(3,4)\), with at least one of the last two
types present. Equation (29) fixes \(F_2\bmod27\).
Every one of the \(729\) representatives has its
constant obstruction at valuation four.

If both first and second digits are zero and one
axis has degree one, take \(I_{41/9}\). Now the
first insertion uses four leading digits for its
binomial, and the second insertion can still use
\(F_2\bmod81\). The four prefixes for the nontrivial
axis are \(0,1,3,4\), with axis height four.
Padding a height-three word with a low zero preserves
the binomial factors modulo \(243\) by the same
reciprocal-sum argument. With the two leading zeros,
every shorter word is zero.

The third insertion term has value at least \(13/3\)
and can be discarded from the reducing identity
at this precision. For the first two forms its value
plus the multiplier bound of one exceeds \(41/9\).
For the \(Z\) forms, \(3\alpha<Q_x\) and \(3\beta<Q_y\)
because the leading digits are zero. The axis
binomials for the third insertion are then one
modulo three, and its contribution to \(DR\)
is \(F_3u^3T^3\) modulo \(I_{16/3}\). The constant-shift
coefficients at offsets \(22,23\) have values \(7/9\)
and \(11/9\), respectively. Their products with
the third insertion also have values above \(41/9\).
All higher insertions are beyond the working
precision. The reduction at the \(324\) resulting
representatives gives a nonzero constant at valuation
\(40/9\).

For two axes of degree above one with both leading
digits zero on each, the same input reduction and
third-insertion omission hold through \(I_{43/9}\).
In fact the third insertion errors in the first two
rows have value at least \(13/3\) before multiplication
by the row coefficients of value at least one.
Its errors in the \(Z\) forms remain above
\(I_{43/9}\), and the fourth insertion has value
at least \(43/9\). Thus there are ten unordered
pairs of the ordinary prefixes \(0,1,3,4\).
Write \(r_x,r_y\) for the third foreign digits.
The \(810\) residue tests at \(I_{42/9}\) consider
the three \(F_2\bmod81\) residues allowed by (29).
For \(3\mid H\) all final initial coefficients
evaluate to a nonzero residue at the forced
parameters. For \(3\nmid H\) the only remaining
cases satisfy \(F_2\equiv44-27(r_x+r_y)\pmod {81}\).
In these \(180\) cases the last reduction at \(I_{43/9}\)
has initial order \(41\) in \(\pi\), with coefficient
\(i(B-h)\), where \(h=1\) for \(H\equiv1\bmod3\)
and \(h=-1\) for \(H\equiv2\bmod3\). Thus
\(v(B-h)\ge1/9\) in an exact block.

The first two moment rows give
\(A\equiv h(1-i),\ B\equiv h\) in the residue field.
The ninth-power comparison used for (28) permits
\(A=h(1-i)+tA'\), \(B=h+tB'\), \(t^3=\pi\), with
\(A',B'\) integral. In the \(h=-1\) case their
order-\(36\) rows in \(t\) give \(\bar B'\ne0\),
contradicting \(v(B-h)\ge1/9\).

In the \(h=1\) case those rows give
\(\bar A'=-1-i\), \(\bar B'=0\). Use
\(A'=-1-i+A_1\), \(B'=t^2B_1\). In suitable
combinations \(X,Y\) of the first two rows the
checker obtains, through the orders indicated,
\[
\begin{split}
 X &= t^{36}A_1^9
       +t^{39}\bigl(-1-i+\text{terms divisible by }A_1\bigr)
       \pmod{t^{40}},\\
 Y+t^3X &= t^{42}(2+i-B_1)A_1^9
       +t^{45}\bigl((1+i)B_1+(i-1)B_1^3
                        +\text{terms divisible by }A_1\bigr)
       \pmod{t^{46}} .
\end{split}
\]
The first equation with \(X=0\) gives
\(t^{-3}A_1^9\equiv1+i\) in the residue field.
Substitute this into the second and take its
order-\(45\) residue to obtain
\(1+(i-1)\bar B_1^3=0\). Thus \(\bar B_1=-1+i\).
In contrast, the final reduced row is
\[
       i\pi^{41}(B-1) +(1+i)\pi^{42}A(B+1)
                   \pmod{\pi^{43}} .
\]
With \(B=1+\pi B_1\) and \(\bar A=1-i\),
its order-\(42\) coefficient is \(i\bar B_1+1\).
This would require \(\bar B_1=i\), a contradiction.
The argument uses exact congruences only through
\(I_{43/9}\), so higher foreign insertions cannot
restore a zero block.

### The second-radius lift for two ternary ones

**Two-digit second-radius exclusion.** In the hypotheses of
Section 9, the four-window at \(L=9\) holds when the
two foreign exponents have altogether exactly two
ternary ones. This includes shifted products
\(X^{3^r}Y^{3^s}\) and axis terms \(X^{3^r+3^s}\)
with \(r<s\), including \(X^4\).
**Second-radius low even-digit exclusion.** Together
with the high digit-count exclusion, this gives the
four-window for the entire sector of Section 9 at
\(L=9\) with \(v(u)>0\).

Use the three starts (25) and the parameters (27).
For a target precision \(I_W\), let
\(w_m=4m/9+v_3(m!)\). The term with \(m\) foreign
insertions has value at least \(w_m\). Here the
largest \(W\) used is \(38/9\), so \(m\ge6\)
is omitted. For \(w_m<W\) the required integer
accuracy is \(h_m=\lceil W-w_m\rceil\).
At that accuracy the insertion data consist of
\(F_m/m!\) modulo \(3^{h_m}\) and the two
cross-binomial polynomials from the finite native
profile calculation, taken modulo \(3^{h_m}\).
Legendre's bound in that calculation gives
\(v_3(F_m/m!)\ge0\). Multiplying the stored residue
by \(m!\) and \(u^m\) gives the required precision
for the moment. The pure axis quotient uses at
most modulus \(243\), at which it depends only on
whether the axis height is zero or positive.

For each foreign word with \(j\) ones, enumerate
the \(j+1\) lengths of the zero runs before,
between and after the ones, including zeros up to
the declared axis height. Each length can be
capped at \(\lceil W\rceil+2\) without changing any
stored data. For the factorial units, use the
product \(g(mn_i+c_{m,i})/g(n_i)^m\) at digit
position \(i\) from the preceding finite-profile
argument. Its window has at most \(h_m\) digits.
Every carry for \(m\le5\) is at most two and
therefore clears at a zero digit. In a long
zero run only the first position can have an
incoming carry, and only the last \(h_m-1\)
positions can see a subsequent one. At each
other position the valuation contribution and
the carry are zero and the unit contribution is
one. The cap retains both edges and keeps their
windows disjoint. For the axis cross binomial,
the retained polynomial uses the leading
\(h_m-1\) digits and the multiplication carry at
their lower edge. With at most two ones this
carry is determined by the two preceding digits
if there are that many. The stated cap preserves
this upper window and whether the height is zero.
Thus it preserves all relevant binomial
polynomials and all residues of \(F_m/m!\).
Every bounded zero-run word is itself allowed.

The [checker](verify_binary_two_digit_second_radius.py)
enumerates the splits \((1,1)\) and \((0,2)\) of
the two ones and interchanges the axes if needed.
At precision three there are eighteen profiles
in the first split and thirty in the second.
The program uses the separated axis factors for
\(R\) and computes \(DR\) by taking the first
difference of each of its axis factors, extended
by zero at axis count \(-1\). Convolution with the
normalized shift coefficients gives
\(\widetilde E_s,\widetilde E_{s+1},Z_{s+2},Z_{s+3}\)
directly. The first two moment rows give the
residue centers and, in the ramified charts, the
lower bounds on parameter deviations by the same
ninth-power comparison as for (28). The reduction
of the last row uses the first two ninth-power
pivots and the linear third pivot. A nonzero value
of its initial coefficient at the forced residues
excludes the block.

For \(W=3\), the scalar estimates following (27)
give the period \(H\bmod9\). The comparison in the
high digit-count lift applies unchanged at \(W=4\)
with this period, since \((DR)_n\in I_1\) for \(n>2\)
also with two foreign digits. It gives period
\(H\bmod27\) for \(4<W\le5\). In the latter two
comparisons the two uncombined rows are transported
at precision one unit below the two \(Z\) rows.
The checker verifies that every multiplier using
either of the first two rows has valuation at
least one, which suffices for the reducing identity
at the target precision. The residue-center and
initial parameter-bound checks use orders below
the transfer cutoffs.

At each later stage the checker tests only
refinements of the preceding profiles and
starting classes that survived. The refinement map
uses the same exponent word with the earlier cutoffs.
Smaller \(h_m\) remove factors already one at that
precision, and simply reduce the foreign residues,
so this map is independent of the representative.
The resulting counts across both splits are
\[
\begin{array}{c|c|r|r}
s & W & \text{tested} & \text{retained}\\ \hline
17 & 3 & 432 & 0\\
20 & 3 & 432 & 189\\
20 & 4 & 486 & 135\\
20 & 38/9 & 405 & 0\\
23 & 3 & 432 & 99\\
23 & 4 & 135 & 18\\
23 & 37/9 & 54 & 0
\end{array}
\]
The precision-four profiles at \(s=23\) are
represented by the axis quartic with axis height
three, and the other axis either linear or
nonlinear. Their final constant obstruction has
valuation four. The precision-four profiles
retained at \(s=20\) have terminal obstructions at
valuation at most \(113/27\). All three starts are
therefore excluded. Replay the reductions from the
repository root with
\(\texttt{python3 research/factorial/image-rigidity/verify_binary_two_digit_second_radius.py}\).

### A unit low foreign coefficient

**Unit foreign coefficient.** Keep the low even-digit
foreign monomial of Section 9 and the two unit axes, but
now take \(u\) a unit. Assume that the reduced native
residue retains both distinct poles. Then the four-window
holds at every noncritical constant valuation and at
the first two critical radii \(1/3,4/9\), with arbitrary
coefficient ramification. If the nonconstant coefficient
grid is prime to \(3\), it holds at all constant valuations.

Here the polynomial part of
\(\overline R(T)=(1+\bar u T-\bar u^2T^2)/
((1-\bar a T)(1-\bar b T))\) is the nonzero constant
\(J=-\bar u^2/(\bar a\bar b)\). Both exponential
coefficients of its proper part are nonzero, so
\(\overline R\) has a two-window at positive indices.
The zero, negative, and noncritical positive constant
arguments in Section 2 apply unchanged to this
two-window. At \(v(k)=0\), index zero of the
three-entry test is a unit, giving a three-window.

At a critical radius use (4). In any row \(d\),
the polynomial part affects only \(t=0\), where its
additional contribution at \(W=1\) is \(J\mu^d/d!\).
At positive \(t\), the two-exponential test of (5)
applies. Any row whose two exponential multipliers
do not both vanish has at most one consecutive zero
at positive indices, and at most two starting at
its first entry. A block that crosses such a row
boundary can have at most one zero from the end
of the preceding row and two from the beginning
of the next. The only exception can be the canceled
row \(d=2\), where entry zero is a unit because of
\(J\), and all positive entries vanish. This zero
run ends before the next period because its initial
entry is \(R_0=1\). For \(L=3\) even this run has
only two zeros. For \(L=9\) a candidate four-window
inside it has \(s=19,\ldots,23\). The coefficient
test (24) remains valid there: with \(u\) a unit,
\(v(z)=-\delta\), so a small digit-sum index is
impossible in every residue class, by the same
least-term comparison. This leaves \(s=23\).
At its index \(t=7\), the numerator coefficient
\(K_7+uK_6-u^2K_5\) has value \(2/3\) from
the unique term \(uK_6\). Thus it cannot vanish
modulo \(3\mathcal O\).

Finally, in a prime-to-three coefficient grid the
multiplicity bound in Section 4 covers any canceled
row at every radius. Any four-zero block inside it
has positive \(t\) at its last two entries, so the
proper part gives (5) at both. The derivative in
\(W\) at \(1\) has nonzero weights at both frequencies
because the roots of \(e_2\) are simple, and cannot
vanish at two consecutive \(t\). At either index the
actual block would force multiplicity at least two.
This contradicts the lattice-variable multiplicity
bound. The mixed term in the polynomial part never
enters these two tests, since only \(t=0\) sees it.

## 10. Scope and remaining cases

[Edo–van den Essen, Definition 2.7 and Conjecture 2.8](https://arxiv.org/pdf/1304.3956)
give the actual-support and every-start target.
[Rocks, Lemma 3.8, Theorem 3.10 and Proposition 3.11](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf)
give the pure prime-power factorial quotient, the unshifted
prime-power family (including powers with the paper's specified
window), and the constant shift when \(p>m\) and \(n<p-m\),
respectively. The argument here uses the local digit formula,
allows mixed and inhomogeneous supports, and treats unbounded
starting orders under the stated residue and valuation hypotheses.
[Liu–Sun, Theorem 2.6](https://doi.org/10.1017/S0004972719000546)
and [Long, Theorem 1.1 and Corollary 1.3](https://hexagonmath.org/pdf/2611.00813v1)
concern all-powers vanishing and nonvanishing at infinitely many
orders. They supply no finite support-size estimate used here.

The next local obstruction is an actual moment block above a
canceled critical row on a coefficient grid whose denominator
is divisible by \(p\). At \(p=3\) with two geometric poles the
exact leading data are a primitive fourth-root pole ratio, the
specified critical value and a canceled row \(d=2\).
At \(L=3\), only support size four with a nonsingleton foreign
column satisfying (10) remains. The explicit congruence solution
shows that the full modulo-three equations alone leave this
profile. Its actual solutions must satisfy (11), (16) and
the vanishing adjacent digits in (19). Within those
conditions the fourth power on an axis is also excluded by
(21)–(22) and the subsequent linked lifts. The bounded
upper profiles give a finite first-radius reduction
at precision four for the other foreign-column cases.
At \(L\ge9\), a candidate block of length at least three is
confined to the canceled row, with its one allowed left
extension when \(c=0\). Only support lengths at most
\(L+\mathbf 1_{c=0}\) need further work there. Equation (6)
retains the linked small native prefixes modulo \(3\mathcal O\),
and Section 7 discards further starts on the full modulo-three
support class. For the low even-digit foreign term,
Section 9 gives the quadratic numerator test and the
three possible starts with precise pole constraints
at \(L=9\). The linked higher-precision equations
exclude all three starts for \(XY\) at this radius.
Other low even-digit foreign terms reduce to finitely
many realizable native profiles at precision four.
For digit counts at least four, (29) gives exactly
the surviving congruence profiles at precision four.
The subsequent linked lifts exclude all three starts
for exact blocks. The two-digit exclusion completes
the low even-digit, positive-valuation sector at
\(L=9\). For a unit low foreign coefficient the analysis
also excludes the first two critical radii and all
radii on prime-to-three grids, when both distinct
poles survive. One must use the linked
higher-precision native coefficients or another exact constraint. For example, a bound on the least term at \(v(k)\)
alone loses a factor \(p\) as soon as the coefficient grid permits
denominators divisible by \(p\). The polynomial native residue
sector requires a separate argument, as do supports whose
residue series has higher degree. Full SFC(2) remains open.
