# Fixed six-row obstructions to bounded high-prime root screens

Fix any terminal \(N=n+1\ge3\) and any integer gap \(D>N\).
For every finite set of primes \(p>ND\) and every choice of positive
precisions \(L_p\), there is a fully shared ordinary profile of
exactly six unequal rows, consisting of mirror pairs with gaps
\(D,1,1\), with the following properties at every selected prime.

- Both normalized factorial polynomials of orders \(n,N\) have
  endpoint coefficient valuations one and interior coefficient
  valuations zero. In particular they have one root of valuation
  \(1\) and one root of valuation \(-1\).
- Write the roots of valuation \(1\) as \(p w_n,p w_N\).
  They are in \(\mathbb Q_p\) and \(w_n,w_N\) are units with
  \(w_n\equiv w_N\pmod{p^{L_p}}\). The roots of valuation \(-1\)
  are their reciprocals, and their normalized units also agree
  modulo \(p^{L_p}\).

The bases can be arbitrarily large, hence the profiles can be chosen
in the [uniform weak-curvature range](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md)
where ordinary coprimality holds at every terminal. Thus the matched
local roots do not give a characteristic-zero common root. The
[fixed-prime residual obstruction](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md#a-simultaneous-obstruction-at-finitely-many-primes)
already treated horizontal Newton polygons with unrestricted numbers
of rows. The present obstruction has fixed dimension and matches
both nonhorizontal faces to arbitrarily prescribed precision.

## The marked and neutral rows

Write \(a=D-1\). Choose the smaller exponent \(A\) on the gap-\(D\)
mirror pair by
\[
                        A+a\equiv p\pmod {p^2}
\]
at every selected \(p\), using the Chinese remainder theorem.
Let \(B,C\) be the smaller exponents on the two unit-gap mirror
pairs. At orders \(j=n,N\), divide the original factorial cells
by the row-minimum factorials, as in the
[translated resultant proof](MULTIVARIATE_TRANSLATED_RESULTANT_DEGREE.md).
The resulting ordinary polynomial \(F_j\) is reciprocal.

Since \(p>ND\), the only sampled factor \(jA+h\),
\(1\le h\le jD\), divisible by \(p\) is the one with \(h=ja\).
Its valuation is one. As \(D>j\), this index is strictly
above \((j-1)D\). In either orientation it occurs only at the
endpoint using the larger entry to order \(j\).
For a neutral unit-gap pair, require its base \(x\) to belong
modulo \(p\) to the set \(S_p\) for which
\[
               jx+h\ne0\pmod p
           \quad(j=n,N,\ 1\le h\le j).
\]
There are at most \(n+N\) excluded residues. Since
\(p>ND>2(n+N)\), one has \(|S_p|>(p-1)/2\).
With \(B,C\in S_p\), this proves the stated coefficient valuations.

## Matching both exterior roots

Fix such a prime \(p\). After substituting \(z=pw\) and
dividing by \(p\), the two reduced polynomials are affine
in \(w\), each with nonzero endpoints. Each has a unique
nonzero residue root. If \(k\) is its index in the original
array, a neutral unit-gap pair contributes
\((jx+1)_k(jx+1)_{j-k}\). Its ratio at \(k=0\) versus \(k=1\)
therefore equals \((jx+j)/(jx+1)\). It multiplies the ratio of
the two reduced root values, for the lower order versus the
upper, by
\[
                       T(x)=\frac{n}{N}\frac{Nx+1}{nx+1}.
\]
This map is injective on \(S_p\) and has nonzero values there.
Any prescribed member of \(\mathbb F_p^\times\) is a product
of two members of \(T(S_p)\): the two subsets \(T(S_p)\) and
the prescribed element times \(T(S_p)^{-1}\) must intersect.
Thus choose residues \(B,C\in S_p\) so that the two reduced
roots have the same value \(w_0\). Do this independently at
every selected prime. Then fix one ordinary positive \(C\)
in the selected classes, leaving \(B\) to be lifted.

For fixed \(A,C,p\), write
\[
                   G_j(B,w)=p^{-1}F_j(A,B,C;pw).
\]
The marked pair's endpoint factors show that \(G_j\) is a
polynomial over \(\mathbb Z_p\) in \(B,w\).
Near the chosen residue \(B\), reduce it to \(A_j(B)+B_j(B)w\),
where \(B_j\) is a unit at the selected residue. The residue
root is the rational function \(r_j(B)=-A_j(B)/B_j(B)\).
The quotient \(r_n(B)/r_N(B)\) is a nonzero constant times
\(T(B)\), as above. Its derivative at the matched point is
nonzero, since
\[
                        T'(B)=\frac{n}{N(nB+1)^2}.
\]
It follows that the derivative of \(r_n-r_N\) is nonzero at
the matched point. The Jacobian determinant of the two
equations \(G_n=G_N=0\) in \((B,w)\) reduces there to
\(-B_nB_N(r_n'-r_N')\), a unit. Successively lifting the two
digits of \(B,w\) uses this invertible linear system at each
power of \(p\). Thus there are \(p\)-adic values \(B_p,w_p\)
solving both equations, with \(w_p\) a unit.

Choose one ordinary positive \(B\) congruent to \(B_p\) modulo
\(p^{L_p}\) at every selected prime. For this \(B\) both
polynomials \(G_j(B,w)\) have a simple residue root, so their
unique unit roots satisfy \(w_j\equiv w_p\pmod{p^{L_p}}\).
Their reciprocal units agree at the same precision.
This proves the claims.

## Size and resultant valuation

To have arbitrarily weak curvature, first choose \(A,C\)
arbitrarily large in their CRT classes, and after lifting do
the same for \(B\). Their gaps remain fixed. If all primes
\(ND<p\le Y\) are selected with a common precision \(L\),
the three CRT moduli divide \(\prod p^2,\prod p,\prod p^L\),
respectively. Taking positive representatives above any fixed
lower bound gives maximal exponent at most
\(C\,2^{2Y\max(2,L)}\), by the elementary estimate
\(\prod_{p\le Y}p\le2^{2Y}\) used in the
[quantitative prime-screen obstruction](MULTIVARIATE_QUANTITATIVE_PRIME_SCREEN_OBSTRUCTION.md#unit-precision-and-curvature-within-the-size-bound).
Here \(C\) depends on \(N,D\) and the chosen lower bound.
The number of rows stays six.

The product formula also yields
\(v_p(\operatorname{Res}(F_n,F_N))\ge2+2L_p\) when the
resultant is nonzero. Indeed the two leading coefficient
valuations contribute \(n+N\). Summing the minima of each
pair of root valuations contributes \(-(n+N-2)\):
there are \(n+N-1\) pairs involving a valuation-\(-1\) root
and one pair of valuation-\(1\) roots. The two matched
exterior pairs improve their respective minima by at least
\(L_p\) each. This valuation bound concerns the scalar-free
polynomials \(F_n,F_N\) defined above.
