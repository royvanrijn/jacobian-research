# Deep prime certificates for opposed rows with a common upper exponent

Consider an original ordinary factorial binomial. Suppose a
group of its unequal rows all have larger exponent \(M\), and
both orientations occur in the group. The gaps may differ
arbitrarily, with arbitrary multiplicities. The retained rows
outside the group are unrestricted. Write \(d_+,d_-\) for the
least gaps in the increasing and decreasing orientations within
the group, respectively. Thus the group contains rows
\((M-d_+,M)\) and \((M,M-d_-)\). All comparisons use the
physical argument of the original binomial.

For adjacent degrees let \(m\ge4\) be either marked order and let
\(m'\) be the other. For each retained unequal row \(i\) put
\(\ell_i=\min(a_i,b_i)\) and \(d_i=|b_i-a_i|\), and set
\[
 B_j=j!\prod_{i\ {\rm retained}}\prod_{h=1}^{j d_i}(j\ell_i+h)
                      \qquad(j=m,m').
\]
Let \(\mathcal U_m\) consist of the nonmultiples \(u\) of \(m\)
in \(0\le u<m\min(d_+,d_-)\) for which
\[
     \left\lfloor\frac u{d_+}\right\rfloor
           +\left\lfloor\frac u{d_-}\right\rfloor
                            \in\{m-3,m-2\}.
\]
For \(u\in\mathcal U_m\) write
\[
     D(u)=1+\max\!\left(
                   \left\lfloor u/d_+\right\rfloor,
                   \left\lfloor u/d_-\right\rfloor\right).
\]
Use the gaps \(d_i\) also to index members of the selected group, and
define the integer constants
\[
 \begin{split}
 C_m(u)&=B_m
       \prod_{i\ {\rm group}}\ \prod_{\substack{0\le h<m d_i\\h\ne u}}
                                      |h-u|,\\
 C_{m'}(u)&=B_{m'}
       \prod_{i\ {\rm group}}\ \prod_{0\le h<m'd_i}|m'u-mh|.
 \end{split}
\]
These constants are positive and independent of \(M\) for fixed
retained rows and group gaps.

**Theorem.** The set \(\mathcal U_m\) is nonempty. A physical common
root of the two original moments forces, simultaneously for
every \(u\in\mathcal U_m\),
\[
                mM-u\ \mid\ C_m(u)C_{m'}(u)^{D(u)}.           \tag{1}
\]
A failure of any one divisibility certifies physical coprimality at
the marked pair. The test includes primes dividing \(m\).
Applying it at the pairs below the
[effective absolute ordinary cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
gives all-degree SFC(2) for the original profile whenever \(M\) is
larger than the resulting computable threshold. The gaps of the
two orientations need no comparability assumption.

The simultaneous tests improve the dependence on the retained
profile. With fixed group gaps, let \(t=|\mathcal U_m|\),
\(D_*=\max_{u\in\mathcal U_m}D(u)\), and suppose the retained
rows have total absolute gap \(G\) and smaller exponents at
most \(B\). There is a computable positive constant \(K\) depending
only on the group gaps, \(m,m'\), and \(G\) such that a common
root forces
\[
    mM-\max\mathcal U_m
         \le K(1+B)^{G(m+m'D_*)/t}.                           \tag{2}
\]
The exponent uses the product of all eligible central marks.

## A common middle interval

Put \(w(u)=\lfloor u/d_+\rfloor+\lfloor u/d_-\rfloor\).
On the indicated offset range it starts at zero and ends at
least \(m-1\). Its jumps are at most two. Hence the preimage of
\(\{m-3,m-2\}\) is a nonempty interval of integers in the
interior of the offset range. If it contains two consecutive
integers, one is a nonmultiple of \(m\). If it consists of one
integer \(u\), both the entry and exit are jumps, and at least
one of them is a jump of two. Thus both \(d_+,d_-\) divide
one of \(u,u+1\), and at least one divides the other. That
gap is one, and \(0<u\le w(u)\le m-2\). This also gives a
nonmultiple of \(m\).

Strip the factorial of the smaller exponent from each row
at the marked order. Its remaining factors can be written
\(mM-h\), \(0\le h<m d_i\), in the group. The factor
\(mM-u\) first occurs at index \(m-\lfloor u/d_i\rfloor\)
of an increasing row. For a decreasing row it occurs up to
index \(\lfloor u/d_i\rfloor\). Thus after marking the
common factor in every row of the group, no marked factor
occurs exactly on the interval
\[
 I=\left[\left\lfloor u/d_-\right\rfloor+1,\
                   m-\left\lfloor u/d_+\right\rfloor-1\right].
\]
Its width is \(m-w(u)-2\in\{0,1\}\). Outside \(I\) at
least one marked factor occurs. The greatest distance
of an index from \(I\) is \(D(u)\). The least gap in
each orientation determines these endpoints even when
there are many more rows of that orientation.

## Comparing the two actual Newton polygons

Take a prime \(p\mid mM-u\) and put
\(E=v_p(mM-u)\), \(U=v_p(C_m(u))\), \(L=v_p(C_{m'}(u))\).
Suppose \(E>U+D(u)L\). At order \(m\), every other factor of
a group row has valuation \(v_p(h-u)\), by writing
\(mM-h=(mM-u)+(u-h)\). After deleting the marked factors,
the coefficient valuations all lie between zero and \(U\).
At order \(m'\), the identity
\[
        m(m'M-h)=m'(mM-u)+(m'u-mh)
\]
has a nonzero second term independent of \(M\).
Indeed \(m\nmid u\) and \(\gcd(m,m')=1\).
Its valuation is at most \(L<E\), so the valuation of
\(m'M-h\) is \(v_p(m'u-mh)-v_p(m)\). It is nonnegative
and at most \(v_p(m'u-mh)\). This bounds every stripped
coefficient valuation at order \(m'\) between zero and \(L\).
In particular all slopes in that polygon lie in \([-L,L]\).

For any real slope in \([-L,L]\), a marked ordinate outside
\(I\), after subtraction of the slope times its index, is
strictly greater than the corresponding value at the nearest
index in \(I\). The difference is at least \(E-U-D(u)L\).
Thus all edges of the marked polygon with a slope also
available in the other polygon lie in \(I\), with total
width at most one. The
[deep-block degree bound](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md#marks-in-arbitrary-blocks)
then gives common factor degree at most one. A rational
common factor must have even degree, because one of the two
moments is an even-order expectation strictly positive
on the real axis. This proves that \(E>U+D(u)L\) is impossible
at a physical common root. Taking the inequalities at every
prime proves (1).

## Multiplying the offset restrictions

The least common multiple of the integers \(mM-u\) divides
the least common multiple of the right sides of (1).
Put \(W=\max\mathcal U_m-\min\mathcal U_m\).
At any prime \(p\), choose a mark of maximal \(p\)-valuation.
For \(p^r\) up to that valuation, all the marks divisible by
\(p^r\) lie in a single progression of spacing \(p^r\) in an
interval of width \(W\). Excluding the maximal mark leaves
at most \(\lfloor W/p^r\rfloor\) such marks. Hence the quotient
of the product of all marks by their least common multiple
divides \(W!\), by Legendre's formula. This gives
\[
 \prod_{u\in\mathcal U_m}(mM-u)
 \ \mid\
       \operatorname{lcm}_{u\in\mathcal U_m}
                    \bigl(C_m(u)C_{m'}(u)^{D(u)}\bigr)\,W!.    \tag{3}
\]
The least common multiple on the right divides
\[
          B_m B_{m'}^{D_*}\bigl((m m' d_{\max})!\bigr)^
                                      {2R(1+D_*)},             \tag{4}
\]
where \(R\) is the number of rows in the group and \(d_{\max}\)
is their largest gap. Indeed each offset product for a single
row has values at most \(mm'd_{\max}\), with each absolute value
occurring at most twice. Thus its product divides
\(\bigl((mm'd_{\max})!\bigr)^2\) independently of \(u\). From the product definition of
\(B_j\), one has \(B_j\le j![j(1+G)]^{jG}(1+B)^{jG}\).
Taking \(t\)-th roots in (3) proves (2). In particular the
retained height occurs only once in the least common multiple.
The offset set can grow with the gaps. For example, if
\(s=1/d_++1/d_-\), every nonmultiple of \(m\) among the integers
in \([(m-2)/s,(m-1)/s)\) is eligible.

Suppose in addition that \(M\) is the largest exponent anywhere
in the profile. For each low pair let
\(\alpha=G(m+m'D_*)/t\) in (2). If \(\alpha<1\), every possible
collision at that pair has an effective bound on \(M\) depending
only on the group gaps, \(G\) and the pair. Indeed \(B\le M\)
and \(mM-\max\mathcal U_m\ge mM/2\) as soon as
\(M\ge2\max\mathcal U_m/m\). Since \(1+M\le2M\), (2) then gives
\(M^{1-\alpha}\le 2^{1+\alpha}K/m\).
This applies at all low pairs, for example, under the simple
sufficient bound \(\min(d_+,d_-)\ge4J^2(G+1)\).
To check it, the interval \([(m-2)/s,(m-1)/s)\) of length \(1/s\)
contains at least \(\lfloor 1/s\rfloor\) integers, with at most
\(\lfloor 1/(ms)\rfloor+1\) multiples of \(m\). Thus
\(t\ge 3\min(d_+,d_-)/8-2\) for \(m\ge4\), and
\(m+m'D_*<J(J-1)\) since \(D_*\le m-1\).
This gives the claimed \(\alpha<1\).

An exact example with unequal gaps is the original pair
\((2,3),(3,1)\), at orders \(5,6\) and \(p=2\). After the same
scalar factorial divisions, the two arrays of valuations are
\((8,7,8,7,4,3)\) and \((12,12,8,9,9,9,6)\).
The first polygon has only slope \(-1\), and the second has
slopes \(-2,-1/2\). Thus the direct polygon test separates the
pair at a prime dividing the marked order. Taking the actual
ordinates can improve on the sufficient depth bound in (1).

For the all-degree assertion use (1) or (2) at every
\(4\le N<J\), marking the even order in the pair \(N-1,N\).
The first two pairs follow from
[real positivity and quadratic nondivisibility](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients).
The absolute cutoff \(J\) handles larger orders, and its
numerical value is unevaluated. A rational Bézout identity
transfers the result to every characteristic-zero field.
