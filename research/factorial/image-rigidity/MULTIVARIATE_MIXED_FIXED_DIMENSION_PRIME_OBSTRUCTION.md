# Mixed horizontal and exterior root obstructions at bounded prime ranges

Fix \(N=n+1\ge3\). Choose a gap \(D>N\) large enough that the
successive distinct fractions from
\[
            \{k/j:j=n,N,\ 1\le k\le\lfloor j/2\rfloor\}
\]
are separated by at least one strictly intermediate fraction of
the form \(u/D\), including between zero and the first fraction.
Fix positive precisions \(L_p\) for the finitely many primes
\(p\le ND\). There is a row count depending only on these fixed
data such that, for any further finite set of primes \(p>ND\)
and any positive precisions \(L_p\) there as well, a fully shared
ordinary mirrored profile of that row count satisfies both of
the following properties.

- For \(p\le ND\), the two moments divided by their constant
  coefficients are congruent to \((1+z)^j\) modulo \(p^{L_p}\),
  for \(j=n,N\).
- For each selected \(p>ND\), the scalar-free polynomials have
  endpoint coefficient valuations one and interior valuations
  zero. Their roots of valuation one have normalized units
  congruent modulo \(p^{L_p}\); the same holds for the reciprocal
  roots of valuation minus one.

The bases can be chosen as large as desired, in particular so that
the [ordinary weak-curvature theorem](MULTIVARIATE_UNIFORM_WEAK_CURVATURE_ANGLE_SEPARATION.md)
proves coprimality at every terminal. For a fixed low-prime
precision and a fixed uniform high-prime precision, prime ranges
through a constant multiple of the logarithm of the maximal
exponent can be treated with a fixed number of rows. The proof
uses the [six-row exterior lift](MULTIVARIATE_SIX_ROW_PADIC_EXTERIOR_OBSTRUCTION.md)
with additional rows whose effects at the low primes are fixed
and whose high-prime factors are units.

## An integer monoid for low-prime valuations

Let \(p\le ND\) and \(j=n,N\). For a mirror pair of gap \(D\)
and smaller exponent \(A\), write
\[
 K_{p,j,k}(A)=
 \sum_{h=1}^{kD}v_p(jA+h)+
 \sum_{h=1}^{(j-k)D}v_p(jA+h)-
 \sum_{h=1}^{jD}v_p(jA+h).
\]
These are the differences of the factorial coefficient valuations
from the constant coefficient, with the binomial factor omitted.
It suffices to take \(1\le k\le\lfloor j/2\rfloor\), by reflection.
Use these coordinates for every low prime in an integer vector.

The following refinement of the
[fixed-prime cancellation](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md#a-simultaneous-obstruction-at-finitely-many-primes)
permits cancellation of any specified integer vector by adding a
finite multiset of actual mirror pairs. List the fractions in
increasing order and put an integer \(u\) strictly between the
adjacent fractions after multiplication by \(D\).
At a selected low prime \(p\), impose
\(A+u\equiv p^H\pmod{p^{H+1}}\). At the other low primes choose
fixed congruences which keep every sampled \(jA+h\) nonzero
modulo a sufficiently high fixed power. Such congruences exist:
each linear congruence eliminates at most \(q^{v_q(j)}\)
classes modulo \(q^b\), and the number of sampled factors is
fixed. These residues and the one at the selected prime are
compatible by CRT. Choose a positive representative.

For all sufficiently large \(H\), the resulting vector is
exactly \(H v+c\), with \(c\) fixed. The exceptional factor
\(h=ju\) has valuation \(H+v_p(j)\), and all others have
valuation independent of \(H\). In the selected prime block
the vector \(v\) is one on the fraction indices above \(u/D\)
and zero elsewhere. It is zero in all other prime blocks.
The analogous construction with \(u\) replaced by \(D-u\)
gives limiting vector \(-v\) after division by \(H\).
The positive vectors \(v\) form an integer basis: they are
the nested step vectors in each prime block.

Take the columns \(P_H\) from all positive choices at height
\(H\), \(P_{H+1}\) from the same choices one height higher, and
\(Q_H\) from the negative choices at height \(H\).
The differences of columns of \(P_{H+1}-P_H\) give the integer
basis just described. The matrices divided by \(H\) tend,
respectively, to \(V,V,-V\), with \(V\) invertible. Hence
\[
          -P_H^{-1}(P_{H+1}{\bf1}+2Q_H{\bf1})
                             \longrightarrow {\bf1}.
\]
For a sufficiently large fixed \(H\), this is a strictly
positive rational vector. Clearing denominators gives an
integer zero relation with strictly positive weights on every
column in all three lists. Every integer vector is already an
integer combination by the column differences. Adding enough
copies of the positive zero relation makes all its
multiplicities nonnegative. This proves the stated monoid
property.

## Keeping the high-prime adjustment independent

Begin with the three mirror pairs of gaps \(D,1,1\) used in the
six-row construction, all with fixed positive prototype bases.
Add finitely many gap-\(D\) mirror pairs given by the monoid
property to cancel their entire low-prime valuation vector.
After this cancellation, the normalized factorial ratios are
units at all low primes. Let \(R\) be a common multiple of
the orders of
\((\mathbb Z/p^{L_p}\mathbb Z)^\times\) over the low primes.
Take \(R\) copies of this prototype multiset.

At a low prime \(p\), fix for every physical copy a base class
modulo a power of \(p\) at least \(L_p\) beyond the largest
valuation of its sampled factors \(jA+h\) at the prototype.
Thus every physical copy keeps both every individual
valuation and every stripped factor unit modulo \(p^{L_p}\).
Its factorial ratios therefore have exactly the prototype's
valuation and stripped residues. The repeated product has
valuation zero and unit residue one. Including the binomial
coefficients gives the claimed congruence at each low prime.

Now designate exactly one physical copy of each of the original
three pairs to be active at the high primes. At a high prime
\(p>ND\), require all other pairs' bases to be zero modulo
\(p\). Every sampled factor of those pairs is a unit because
\(1\le h\le ND<p\). The active marked pair is assigned the
same class modulo \(p^2\) as in the six-row construction.
Hence its single endpoint hit is still the only hit at this
prime. The remaining pairs contribute only a nonzero fixed
factor to the ratio of the two normalized residue roots.
The two active neutral unit-gap bases can absorb this factor
by the same two Möbius images and the same invertible
two-variable lift. Their contribution is used once at a high
prime even though the different physical copies agree at
the low primes. The high-prime exterior slopes therefore
remain \(1,-1\), with the stated normalized units.

Take the low and high congruences together by CRT, fixing the
shadow bases and the active marked base before the two active
neutral adjustments. The common low-prime moduli are fixed
independently of the selected high primes. Every base has
high-prime modulus dividing \(\prod_{ND<p\le Y}p^{\max(2,L)}\)
when all high primes up to \(Y\) are selected with precision
\(L\). By taking representatives above a fixed lower bound,
the maximal exponent is at most
\(C2^{2Y\max(2,L)}\), using the same elementary least-common-
multiple estimate as in the six-row construction. This proves
the size and fixed-dimension statement.
