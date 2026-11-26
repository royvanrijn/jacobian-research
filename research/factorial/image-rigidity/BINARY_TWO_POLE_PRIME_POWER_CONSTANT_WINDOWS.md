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

## 9. Scope and remaining cases

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
the vanishing adjacent digits in (19).
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
