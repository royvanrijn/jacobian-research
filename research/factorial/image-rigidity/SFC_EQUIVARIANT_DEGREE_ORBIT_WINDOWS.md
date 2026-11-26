# Factorial windows for one coordinate orbit at each degree

Let \(H\) be a finite group of permutations of the variables.
Let a finite nonempty support of distinct nonnegative
integer exponent columns be \(H\)-stable, with each
occupied degree forming one orbit. Fix a character
\(\chi:H\to\mathbb C^\times\) with image of order \(e\).
For a sufficiently large positive integer \(D\), form
\[
                    F_D(X)=\sum_i c_i X^{Db_i}
\]
and require \(F_D(\sigma X)=\chi(\sigma)F_D(X)\) for all
\(\sigma\in H\). The coefficients can vary with \(D\).
Let \(b\) be the number of occupied degrees after dropping
zero coefficients.

**Theorem.** There is a dilation threshold independent
of the coefficients. Above it, every \(b\) consecutive
positive orders divisible by \(e\) contain an order
with nonzero factorial moment. In particular the
ordinary factorial window has length \(eb\), which is
at most the actual support size. The theorem also
holds over any characteristic-zero field, with a
character into that field's multiplicative group.

There are no face restrictions on the exponent orbits.
The surviving foreign monomials can differ between
peaks of one degree. For example, with a transposition
of \(X,Y\) and its sign character, the theorem allows
independently chosen orbit coefficients in
\[
           a(X-Y)+b_0(X^2-Y^2)+c(X^3Y-XY^3)
\]
after common dilation. At the two degree-one peaks,
different foreign degree-two terms survive. For the
trivial character the window is the number of occupied
degrees, even if their orbits contain many monomials.

## Orbit reduction and the early orders

Coordinate permutations preserve the factorial
functional, or equivalently the product-exponential
integral. For \(e\nmid t\), equivariance therefore
gives \(\Lambda(F_D^t)=0\). In a nonzero orbit the
stabilizer of any column is contained in
\(\ker\chi\), since a stabilizer fixes the coefficient
on that column. Hence its size as an orbit is a
multiple of \(e\), and the actual support size is
at least \(eb\). Dropping a zero orbit preserves
every hypothesis.

Choose peak boxes equivariantly under the
coordinate permutations. This can be done with
a common box radius small enough for separation,
using a positive interval on active coordinates
and an interval starting at zero on inactive
coordinates. At each order divisible by \(e\),
the original scaled integral over each of the
boxes in an orbit is the same complex number,
by change of variables and \(\chi^t=1\).
Thus one can use one representative per orbit
and multiply the integral by its positive orbit
size. This equality is exact and applies to
the full local phase, including foreign terms.

Use the normalized profiles \(P_i\) and raw weights
\[
   q_i(t)=c_i(Dt)^{D d_i}
                  \left(e^{-d_i}\prod_\nu b_{\nu i}^{b_{\nu i}}\right)^D
\]
in the [separated-support estimates](SFC_DILATED_SUPPORT_SADDLE_ENVELOPE.md#exact-integral-and-separation-of-the-peaks).
Here \(0^0=1\). Across one orbit the profiles are
coordinate permutations and the weights have
equal moduli. Their log moduli as functions of
\(\log t\) are one affine line per occupied degree.

Consider \(b\) consecutive sampled orders divisible
by \(e\), starting at \(m\). For small fixed \(\eta>0\)
and \(m\le e^{\eta D}\), the change-point estimate
in the [exponential-range proof](SFC_MIXED_DEGREE_TOP_BLOCK_TAIL_WINDOWS.md#all-coefficient-vectors-up-to-an-exponential-starting-order)
shows that at most \(b-1\) of these integers are
nonclean. In detail, choose \(L>r+s+2\), where \(r\)
is the number of variables and \(s\) is the size
of the given support. A nonclean integer is at
distance at most \(2L\log(tD)/(Dg)\) from an
envelope change, where \(g>0\) bounds the distinct
degree gaps below. For \(t\le e^{\eta D}+es\)
this distance is less than \(1/3\) if \(\eta\)
is small and \(D\) sufficiently large. Each change
accounts for at most one integer. When there is
just one degree all sampled orders are clean.

At a clean sampled order \(t\), only one orbit
leads. Divide by \(Q(t)^t(tD)^r\), where
\(Q(t)=\max_i|q_i(t)|\). On the real box of one
leading label the normalized integral differs
from the pure leading monomial integral by at
most \(O(t e^{-cD})\) times its magnitude. Indeed
the other terms in the box before taking the
\(t\)-th power have relative size \(O(e^{-cD})\)
by peak separation and the bound on all raw
weights. The pure integral on that box has
magnitude bounded below by a positive constant
times \((tD)^{-\rho_i}\), with
\(\rho_i=r-|S_i|/2\). This follows directly by
integrating over active widths proportional
to \((tD)^{-1/2}\) and inactive widths proportional
to \((tD)^{-1}\) about the peak. The leading
orbit sum is its orbit size times one local
integral, by the exact symmetry.

For every other box the maximum relative
raw weight is \((tD)^{-L/t}\). The same foreign
box estimate as in the exponential-range
proof bounds its normalized integral by
\(O((tD)^{-L})\). The complement of the boxes
loses \(O(e^{-ctD})\), from the common integrable
peak envelope. Taking \(\eta\) small enough
also gives \(t e^{-cD}=o(1)\). These errors
are smaller than the leading orbit sum, so
the clean sampled moment is nonzero.

## The large orders and the ordered determinant

Suppose \(m>e^{\eta D}\). Use the local expansion
for a [full, order-dependent singleton phase](SFC_MIXED_DEGREE_SINGLETON_FOREIGN_BLOCKS.md#the-shifted-expansion)
on each representative box. This local construction
requires only separated individual centers. It
includes every foreign term in the full phase
and does not use singleton status until the
linear algebra after the expansion. Apply it
at the fixed shifts \(h=e\ell\), with
\(0\le\ell<b\). Bounds on a fixed finite range
of shifts are the same as for consecutive ones.

Let \(Q=\max_i|q_i(m)|\), and admit an orbit if
\(|q_i(m)|\ge e^{-\gamma D}Q\) for a representative,
with \(\gamma\) small relative to the fixed distinct
degree gaps. Eligibility is independent of its
representative. Use the notation \(u_i,v_i,d_*\)
from the singleton phase expansion. In particular
the frequency satisfies
\[
      c e^{-\gamma D}e^{D(d_i-d_*)}\le |v_i|
              \le C e^{D(d_i-d_*)},\qquad d_*=\max_i d_i ,
\]
and \(|u_p|\ge e^{-C e^{-cD}}\) for a raw
maximizing label \(p\). On a sampled order the
mode of each admitted orbit, including its
positive orbit size, is
\[
          u_i^m(v_i^e)^\ell (mD)^{-\rho_i} B_{i,\ell},
\]
where, after including the orbit size in \(K_i\),
\[
              B_{i,\ell}=K_i+o(1),\qquad 0<c\le K_i\le C .
\]
Indeed the pure degree correction
\(\exp(Dd_i[(m+h)\log(1+h/m)-h])\) tends to one,
and the peeled full saddle amplitude at a
fixed shift tends to its nonzero pure value as
in the singleton estimate. The time correction
beyond the linear phase is already part of that
amplitude. The \(o(1)\) estimates are uniform for
the fixed support and bounded shifts. Omitted
boxes and the complement contribute
\(R_\ell=O(e^{-cmD})\) in the common normalization
of that expansion.

Let \(a\) be the number of admitted orbits. Use
the first \(a\) sampled orders. The matrix of
the admitted modes has entries
\[
                 M_{\ell i}=(v_i^e)^\ell B_{i,\ell},
                    \qquad 0\le\ell<a .
\]
Order the columns by increasing degree. Since
each degree has only one orbit, consecutive
column moduli \(|v_i^e|\) increase by a factor
at least \(\exp(cD)\). In the determinant the
term assigning exponent \(\ell=i-1\) to the
\(i\)-th column is strictly largest in its
frequency factors. Any other assignment has
an inversion; swapping inverted exponents
increases the frequency magnitude by at least
\(\exp(cD)\). Ordering all exponents therefore
shows that every other term loses this factor.
The \(B_{i,\ell}\) are bounded above and below
in modulus, so the sum of the other terms
cannot cancel the ordered term. In particular
\(|\det M|\ge e^{-CD}\), since every frequency
also has magnitude at least \(e^{-CD}\) and
the matrix size is bounded. Its entries are
uniformly bounded above, so
\(\|M^{-1}\|\le e^{CD}\).

If all sampled moments vanish, this system
gives
\[
        |u_i|^m(mD)^{-\rho_i}\le C e^{CD}e^{-cmD}
\]
for every admitted orbit. But the raw maximizing
orbit satisfies
\(|u_p|^m\ge\exp(-Cm e^{-c'D})\), contradicting
the displayed bound for large \(D\). The argument
also works for one admitted orbit. Hence a
sampled moment is nonzero in the large range.

Every \(eb\) consecutive positive integers contain
exactly \(b\) consecutive multiples of \(e\).
This proves the stated window. For a general
characteristic-zero field embed the finitely
generated field of coefficients and character
values into \(\mathbb C\). Coordinate symmetry
and the integer factorial equations are preserved.
The actual support-size convention is that of
[Edo and van den Essen, Conjecture 2.8](https://arxiv.org/abs/1304.3956).
