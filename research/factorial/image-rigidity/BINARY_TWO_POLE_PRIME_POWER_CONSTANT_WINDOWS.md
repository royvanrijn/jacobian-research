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

## 8. Scope and remaining cases

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
At \(L=3\), the only remaining support length is four, with one
foreign column of natural valuation \(0<v(u)\le1/3\).
At \(L\ge9\), a candidate block of length at least three is
confined to the canceled row, with its one allowed left
extension when \(c=0\). Only support lengths at most
\(L+\mathbf 1_{c=0}\) need further work there. Equation (6)
retains the linked small native prefixes modulo \(3\mathcal O\),
and Section 7 discards further starts on the full modulo-three
support class. One must use the linked
higher-precision native coefficients or another exact constraint. For example, a bound on the least term at \(v(k)\)
alone loses a factor \(p\) as soon as the coefficient grid permits
denominators divisible by \(p\). The polynomial native residue
sector requires a separate argument, as do supports whose
residue series has higher degree. Full SFC(2) remains open.
