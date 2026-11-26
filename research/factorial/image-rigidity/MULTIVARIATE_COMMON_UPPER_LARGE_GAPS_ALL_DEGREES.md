# All-degree coprimality for large gaps at a common global upper exponent

For distinct nonnegative integer exponent columns \(a,b\), let
\[
 Q_j(z)=\sum_{k=0}^j\binom jk
                \prod_i((j-k)a_i+kb_i)!\,z^k .
\]
Discard equal rows, which contribute only scalar factors, and let
\(M\) be the largest exponent in the remaining rows. Select a
nonempty group of \(R\) rows within which every larger exponent
equals \(M\). Write \(D\) for the least gap in the group,
\(G\) for the sum of the absolute gaps outside it, and suppose
the group gaps belong to \([D,AD]\), for some integer \(A\ge1\).

**Theorem.** For each \(A\) there are computable constants
\(\delta>0\) and \(D_0\) such that
\[
                      D\ge D_0,\qquad RG\le\delta D
\]
imply coprimality of \(Q_n(z)\) and \(Q_{n+1}(z)\) for every \(n\ge1\)
over every characteristic-zero field. The constants are independent
of \(R\) and of the exponents outside the group. In particular,
with bounded group size the result applies to \(G=o(D)\), even if
some of the retained exponents stay small.

Fix an integer \(J\ge4\) at least the
[effective absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
The constants below depend on \(A,J\). The numerical value of \(J\)
is unevaluated. The argument uses one progression of integers at any
possible low-order contact. Deep prime marks constrain its large
prime factors. Removing terms of largest valuation for the small
primes then gives an incompatible upper bound on the remaining
product.

## Reduction to a progression of eligible marks

By the
[curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md),
only
\[
                 S=\sum_i\frac{(b_i-a_i)^2}{\max(a_i,b_i)}<12
\]
needs consideration. In particular \(M>D^2/12\). If the whole
selected group has the same orientation, the signed total degree
difference has magnitude at least \(RD-G\). The
[signed degree difference theorem](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
handles that case as soon as \(D\ge292J\) and \(\delta\le1/2\).
Assume the group has both orientations.

The pairs of terminal indices two and three are covered by
[real positivity and quadratic nondivisibility](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients).
By the absolute cutoff it remains to exclude a common root for
\(4\le m<J\), with \(m'=m-1\). Let \(d_+,d_-\) be the least
group gap in each orientation and set \(s=1/d_++1/d_-\).
Select the integers
\[
 u\in[(m-2)/s,(m-1)/s),\qquad u\equiv1\pmod m .
\]
They form an arithmetic progression of length
\(k\ge D/(2m)-1\), because \(1/s\ge D/2\). They are all members
of the eligible set in the
[deep mark theorem](MULTIVARIATE_COMMON_UPPER_EXPONENT_DEEP_MARKS.md#comparing-the-two-actual-newton-polygons).
Indeed \(su\) belongs to \([m-2,m-1)\), so the sum of its two
floors belongs to \(\{m-3,m-2\}\), and \((m-1)/s<mD\).
In particular the marks are nonmultiples of \(m\) and are
strictly below \(mD\).

Write \(x_u=mM-u>0\) for the corresponding descending progression,
and \(W\) for its span. Then \(\gcd(x_u,m)=1\), \(W\le JD\),
and \(k\le JD\). Take \(D\ge48J^2\), so \(k\ge D/(4J)\)
and \(x_u>m(M-D)\ge M>D^2/12\).

## Bounding the possible height

For retained rows let \(\ell_i=\min(a_i,b_i)\), \(d_i=|b_i-a_i|\),
and write
\[
                    B_j=j!\prod_{i\ {\rm retained}}
                                  \prod_{h=1}^{j d_i}(j\ell_i+h).
\]
Every retained factor is at most \(jM\), so
\(B_j\le j!(JM)^{jG}\).
Set
\[
 H=J^2A,\qquad C=8J(2J+4JH+2J^2),\qquad
 \delta=\frac{1}{128J^2(C+1)}.
\]
The deep mark theorem says that a common root forces
\(x_u\mid C_m(u)C_{m'}(u)^{D(u)}\) for all selected \(u\).
Here \(D(u)\le m-1\), and the factors of \(C_m(u)/B_m\)
and \(C_{m'}(u)/B_{m'}\) are the nonzero absolute offset
differences \(|h-u|\) and \(|m'u-mh|\), respectively.
All are at most \(HD\). In each row a given positive offset
magnitude occurs at most twice, since it is the absolute
value of a nonconstant affine function of \(h\). Therefore
\[
       \operatorname{lcm}_u x_u
          \mid B_m B_{m'}^{m-1}((HD)!)^{2Rm}.
\]
For every prime \(p\), remove an index of largest \(p\)-valuation
from the valuation sum for \(\prod_u x_u\). Among the other
indices, at most \(\lfloor W/p^e\rfloor\) are divisible by
\(p^e\) for each \(e\ge1\). Distinct multiples are at distance
at least \(p^e\) and the removed term is divisible by every
\(p^e\) that occurs. Summing over \(e\) shows
\[
 \frac{\prod_u x_u}{\operatorname{lcm}_u x_u}\mid W!.
\]
Combining the two divisibilities gives
\[
 M^k\le W!m!(m'!)^{m-1}(JM)^{T}((HD)!)^{2Rm},
 \qquad T=G(m+m'(m-1))\le J^2G\le k/2 .
\]
For \(D\ge H\), both \(W\) and \(HD\) are at most \(D^2\).
The logarithm of the product on the right after removing \(M^T\)
is at most
\[
 [2J+4RJH+2J^2]D\log D .
\]
Indeed \(W!\) costs at most \(2JD\log D\), the group factorial
costs at most \(4RJHD\log D\), and the two order factorials
and \(J^T\) cost at most \(2J^2D\log D\) together.
Since \(k/2\ge D/(8J)\), the claimed contact requires
\[
                              M\le D^{CR}.                 \tag{1}
\]

## A prime that excludes contact

A prime \(p>HD\) dividing one of the marks \(x_u\) and
not dividing \(B_mB_{m'}\) also avoids every group offset
difference. The divisibility required by the deep mark theorem
would then fail at that \(u\). Thus under a hypothetical
contact, every prime \(p>HD\) dividing a mark divides
\(B_mB_{m'}\).

Such primes divide at most one mark each because \(W\le HD\).
None divides the order factorials \(m!,m'!\). Each of the at
most \(2JG\) retained factors is at most \(JM\), and by (1)
has at most \(R(C+1)\) distinct prime factors above \(HD\).
Indeed \(\log(JM)\le (CR+1)\log D\) and \(\log(HD)\ge\log D\).
Remove every mark divisible by any of these large primes.
The number removed is at most
\[
             2JGR(C+1)\le D/(64J)\le k/16 .
\]

Next remove, for each prime \(p\le HD\) dividing a mark,
one mark of largest \(p\)-valuation in the original progression.
These choices may coincide with previous deletions. The elementary
bound \(\pi(x)\le9x/\log x\) for \(x>1\) suffices. To see it,
write \(\vartheta(x)=\sum_{p\le x}\log p\). Primes in
\((2^{i-1},2^i]\) divide \(\binom{2^i}{2^{i-1}}\).
Summing the log bounds from the binomial theorem yields
\(\vartheta(x)\le4x\log 2\). Splitting the prime count at
\(\sqrt{x}\), and using \(\log x\le\sqrt{x}\), gives the
claim. If \(\log D\ge288JH\), then
\(\pi(HD)\le D/(32J)\le k/8\).
At least \(3k/4\) marks remain after both deletions.

For \(p\mid m\) the marks are already units. For each other
\(p\le HD\), the progression step \(m\) is invertible modulo
\(p^e\). Thus the indices of \(p^e\)-divisible terms form
one congruence class, and after removing a maximal-valuation
term at most \(\lfloor(k-1)/p^e\rfloor\) such indices remain.
There are no \(p>HD\) in the remaining product. Its valuation
at each prime is therefore at most that of \((k-1)!\).
But each remaining mark exceeds \(D^2/12\), which gives
\[
 \prod_{\text{remaining } u}x_u
      > (D^2/12)^{3k/4}>(JD)^k\ge(k-1)! .
\]
The middle inequality follows from \(D\ge48J^2\) and
\(48>12^{3/2}\). This contradiction excludes a common root.

One can take \(D_0=\max\{H,48J^2,292J,2^{1024JH}\}\), which
ensures all the above bounds, including \(\log D\ge288JH\)
by \(\log 2>1/2\). A rational Bézout identity transfers
coprimality to every characteristic-zero field.

For example, the original three-row profile
\((D^2-D,D^2),(D^2,D^2-D),(1,2)\) has \(A=1\), \(R=2\),
\(G=1\), curvature \(5/2\), and signed total gap \(1\).
The theorem applies for all sufficiently large \(D\), even though
the retained row has fixed positive curvature and fixed relative gap.

## Upper exponents close to the global maximum

Exact agreement of the larger exponents in the group can be
relaxed. Suppose its larger exponents are \(M_i=M-e_i\), with
\[
                  0\le e_i\le D/4,\qquad D\le d_i\le AD .
\]
The same condition \(RG\le\delta D\) gives the all-degree conclusion
with the same constants. The retained rows need only have
exponents at most \(M\).

Only the construction of the progression and verification of the
offsets change. If the group has both orientations, put
\(E_0=\max_i e_i\), \(L_0=\min_i(e_i+d_i)\), and define
\[
 L_\pm(u)=\max_{i\ {\rm of\ orientation}\ \pm}\frac{u-me_i}{d_i}
\]
for real \(u\) on \([mE_0,mL_0]\).
Here \(L_0\ge D>E_0\). Both functions are strictly increasing
and each is \(1/D\)-Lipschitz. At the left endpoint their sum
is at most \(m/2\le m-2\). At the right endpoint a row attaining
\(L_0\) contributes \(m\), and the other orientation is
nonnegative. Hence
\[
                         m-2\le L_+(u)+L_-(u)<m-1
\]
defines a half-open interval of length at least \(D/2\) contained
in \([mE_0,mL_0)\). Choose its integers \(u\equiv1\pmod m\).
Since a gap-\(D\) row gives \(L_0\le E_0+D\le5D/4\), the
progression has \(k\ge D/(2m)-1\), \(W\le JD\), \(k\le JD\) and
\(x_u=mM-u>M\) under the size assumptions in the proof.

In every group row the factor \(x_u\) has local offset
\(u-me_i\) at order \(m\), between zero and \(md_i-1\).
It is present at increasing indices \(k\ge m-\lfloor L_+(u)\rfloor\)
in at least one increasing row, and at decreasing indices
\(k\le\lfloor L_-(u)\rfloor\) in at least one decreasing row.
These exactly describe the union of the marked indices.
The sum of the two floors is \(m-3\) or \(m-2\), so the
unmarked interval has width one or zero.

At the marked order use the absolute offsets \(|u-me_i-h|\)
for the other factors. At the other order use
\[
 m(m'M_i-h)=m'(mM-u)+(m'u-mm'e_i-mh).
\]
The last term cannot be zero since \(m\nmid u\).
It is the difference of \(m'(u-me_i)\) and \(mh\), both
nonnegative and less than \(mm'd_i\). Thus all offset
magnitudes are at most \(HD\) with the original \(H=J^2 A\),
and each occurs at most twice per row. After deleting every copy of \(x_u\) where it occurs, the
cleaned ordinates lie in \([0,U]\), where \(U\) bounds the
valuation of the retained factors and the other marked-order
offsets. If \(E=v_p(x_u)\) exceeds both this bound and the
other-order offset bound \(L\), the other array has slopes
in \([-L,L]\), by the displayed identity. Each marked ordinate
outside the unmarked interval exceeds the nearest inside
ordinate, after subtraction of any slope in this range, by
at least \(E-U-D(u)L\), where \(D(u)\le m-1\) is the largest
distance to that interval. Thus \(E>U+D(u)L\) again confines
every shared edge to width at most one. The
[deep-block comparison](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md#marks-in-arbitrary-blocks)
then excludes a common factor by even-degree parity. This
gives the same necessary divisors and the same factorial bound
on their least common multiple. The height estimate and deletion
of prime multiples apply with the same constants, proving the
corollary.
