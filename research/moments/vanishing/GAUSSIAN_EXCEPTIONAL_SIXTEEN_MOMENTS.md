# Local sixteen-moment rigidity at the exceptional ternary cubic

Let \(V=\mathbb C[Z,W,Y]_{\leq3}\) and use the normalized Gaussian Wick
functional with \(E[ZW]=E[Y^2]=1\) and all other basic pairings zero.
Write \(\mu_j(P)=E[P^j]\). Put \(P_*=W+ZY-Z^3/2\), and let
\(\mathcal O\) be its orbit under \(\mathrm O(2ZW+Y^2)\times\mathbb G_m\),
acting by orthogonal substitutions and scalar multiplication.

**Theorem (sharp local moment prefix).**
The reduced germ of \(\{\mu_1=\cdots=\mu_{16}=0\}\) along \(\mathcal O\)
is \(\mathcal O\). The local intersection on a transverse slice at
any point of the orbit is zero-dimensional and has length \(40040\). In every Zariski neighbourhood of \(P_*\) there are
points with the first fifteen moments zero and \(\mu_{16}\ne0\).
Such witnesses can be chosen with algebraic coefficients.

In particular, sixteen is the least local set-theoretic prefix at the
exceptional orbit. The universal sixteen-moment problem also requires
control away from this orbit.

## A transverse chart and six formal parameters

The polynomial shear
\[
 Z\mapsto Z,\quad Y\mapsto Y-Z^2,\quad W\mapsto W+ZY-Z^3/2
\]
preserves the Gaussian functional and sends \(W\) to \(P_*\), so \(P_*\)
is pure. The orbit tangent vectors are
\[
 \begin{split}
 P_* &= W+ZY-Z^3/2,\\
 (Z\partial_Z-W\partial_W)P_* &= -W+ZY-3Z^3/2,\\
 (Z\partial_Y-Y\partial_W)P_* &= Z^2-Y,\\
 (W\partial_Y-Y\partial_Z)P_* &= ZW-Y^2+3Z^2Y/2.
 \end{split}
\]
Their coefficient matrix on \(W,ZY,Z^2,Z^2Y\) has determinant \(3\).
Fix these four coefficients to their values at \(P_*\). The resulting
affine slice \(S\) has dimension sixteen. The action map from the
group near the identity times \(S\) has invertible differential at
\((1,P_*)\) and gives an isomorphism on formal completions. For example,
the orthogonal coordinates can be supplied by the Cayley transform of
a skew-adjoint \(3\times3\) matrix near zero. The formal inverse follows
by successively solving each homogeneous degree using the invertible
linear part. Orthogonal substitution preserves every moment, and
scalar multiplication multiplies \(\mu_j\) by a unit on the chart.
Thus each prefix on the completion is the product of its slice scheme
with four smooth orbit directions.
[Purity of the shear and the mixed derivative cutoff](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md#3-the-exceptional-shear-exact-original-images-and-sharp-jet-order)
· [Orbit and separation calculation](GAUSSIAN_CUBIC_PREFIX_GEOMETRY.md)

Write slice increments as \(\sum_{i=1}^{10}v_i M_i + aT_a+bT_b+cT_c+dT_d+eT_e+fT_f\), with
\[
 (M_1,\ldots,M_{10})=(1,Y,Y^2,Y^3,WY,WY^2,W^2Y,W^3,Z,ZW^2)
\]
and
\[
\begin{array}{ll}
 T_a=W^2+(11Y-Y^3)/4,    & T_b=2WY-3Z+ZY^2,\\
 T_c=ZW+(Y^2-3)/2,      & T_d=ZWY+(Y^3-7Y)/2,\\
 T_e=Z^2W-WY-Z,         & T_f=Z^3+2WY^2/9+4ZW^2/9 .
\end{array}
\]
These sixteen polynomials form a basis of the allowed increments.
Indeed, each \(T\) has a distinct monomial with coefficient one
outside the pivot list, and its remaining monomials are pivots.

The ten-by-ten matrix \(D_{ji}=jE[P_*^{j-1}M_i]\), \(1\le j\le10\),
is invertible, and each \(T\) is killed by these ten differentials.
These assertions are checked below by the original Wick formula.
Hence the first ten moments solve uniquely for \(v_i=v_i(a,b,c,d,e,f)\)
as formal series with zero constant and linear terms. Equivalently,
replacing the ten \(v_i\) by the first ten moments is an invertible
formal coordinate change. After elimination let \(F_j\) denote the
remaining moment series, \(11\le j\le16\). Their linear parts before
the elimination are zero, as the original cubic marker cutoff at the
shear is nine in the mixed moment index. The finite verification
also checks this vanishing for all six required \(j\).

## Weighted initial forms and their certificate

Give \((a,b,c,d,e,f)\) weights \((3,2,1,3,2,4)\). The original Wick
calculation gives \(\operatorname{ord}_{w}(F_j)=j\) for \(11\le j\le16\).
Denote the weight-\(j\) leading parts by \(H_j\). The finite computation
checks every lower weight, including all ordinary quadratic terms,
which have weight at most eight. Here is a recurrence specifying the
six forms without expanding the full moment ideal.

For \(n\ge0\) and a spatial monomial \(Z^i W^j Y^k\), expand \(P_*^n\)
by using \(s\) copies of \(ZY\), \(r\) copies of \(-Z^3/2\) and the
remaining copies of \(W\). A term contributes only when
\(i+s+3r=j+n-s-r\) and \(k+s\) is even. In that case its contribution to
the period of the monomial is
\[
 \binom{n}{r}\binom{n-r}{s}(-1/2)^r
       (i+s+3r)!\,\frac{(k+s)!}{2^{(k+s)/2}((k+s)/2)!}.
\]
At weight \(m\), take the known lower-weight slice increment \(A_{<m}\).
The nonlinear contribution to \(\mu_n\) at that weight is
\[
 B_{n,m}=\sum_{k=2}^{\min(n,m)}\binom nk
          E\bigl[P_*^{n-k}(A_{<m}^{\,k})_{[m]}\bigr],
\]
where \([m]\) denotes the weight-\(m\) part. Multiply the first ten
entries of \(B_{\bullet,m}\) by \(-D^{-1}\) to get the ten pivot
corrections at weight \(m\). Add the prescribed \(uT_u\) terms of
weight \(m\) to obtain the full increment in that weight. This works
because nonlinear products at weight \(m\) only involve positive
weights below \(m\). For \(11\le n\le16\), the linear term vanishes, so
\(B_{n,m}=0\) for \(m<n\) and \(H_n=B_{n,n}\), as verified. The
recurrence therefore computes the required coefficients of the exact
formal elimination.

The [certificate](gaussian_exceptional_sixteen_moment_certificate.json.gz)
contains polynomial multipliers modulo `32003`. At \(c=0\) the
identities express \(a^{12},b^{18},d^{12},e^{18},f^9\) in the ideal of
the six specialized forms. The verifier checks these identities
by direct multiplication. At \(c=1\) the reduced graded reverse
lexicographic Gröbner basis of the six specialized forms is \((1)\).
The verifier recomputes that basis using msolve with exact sparse
linear algebra. A separate Singular computation with the added
generator \(c-1\) gives the same unit ideal. All input coefficients
have denominators prime to \(32003\).

If a common zero of the reduced forms over \(\overline{\mathbb F}_{32003}\)
has \(c\ne0\), rescaling every variable of weight \(w\) by \(c^{-w}\)
gives a common zero at \(c=1\). The unit identity excludes that.
At \(c=0\) the five power identities force all remaining coordinates
to be zero. Thus the modular weighted quotient is finite-dimensional.
In each sufficiently large weight its homogeneous Macaulay matrix
has full column rank modulo \(32003\). Its entries are rational with
denominators prime to \(32003\), so the corresponding matrix has full
column rank over \(\mathbb Q\). Therefore the six rational \(H_j\)
generate a homogeneous ideal of finite codimension. This is a
weight-by-weight rank lift, so it uses no inference of rational ideal
membership from membership of an individual modular polynomial.

## The local scheme and multiplicity

In the positively weighted polynomial ring the six \(H_j\) form a
regular sequence. In detail, they are a system of parameters at the
homogeneous maximal ideal, and the polynomial local ring is regular
and hence Cohen--Macaulay. Homogeneity makes regularity global,
since a nonzero homogeneous kernel element cannot be killed by a
polynomial with nonzero constant term. The commutative algebra used
here is the same system-of-parameters criterion as in the
[five-direction slice proof](GAUSSIAN_GENERIC_FIVE_DIRECTION_MOMENT_SCHEMES.md#3-complete-intersection-supplies-a-rational-lower-hilbert-bound).

The weight-graded ideal of the formal ideal \((F_{11},\ldots,F_{16})\)
is exactly \((H_{11},\ldots,H_{16})\). To check the reverse inclusion
after the evident one, take a representation of a nonzero element
as a combination of the \(F_j\). If its lowest possible weight
cancels, the lowest-weight coefficients give a homogeneous syzygy
of the \(H_j\). The Koszul relations generate the syzygies of a
regular sequence. Indeed, reducing a relation modulo all but its
last generator expresses its last coefficient in the preceding
ideal, and subtracting pair relations removes that coefficient.
Induction gives the assertion with homogeneous coefficients.
Subtracting the corresponding Koszul relations
of the \(F_j\) from the coefficients preserves the represented
element and increases the minimum possible weight. Repeat until
its actual lowest weight is reached. Each step increases an integer
bounded by that actual weight, so this argument ends.

Consequently the associated graded slice quotient has Hilbert series
\[
 \frac{\prod_{j=11}^{16}(1-z^j)}
 {(1-z)(1-z^2)^2(1-z^3)^2(1-z^4)} .
\]
The finite length follows either by division in the complete weighted
filtration or by Nakayama's lemma and separation of that filtration.
Evaluating the series at \(z=1\) gives
\[
 \frac{11\cdot12\cdot13\cdot14\cdot15\cdot16}{3\cdot2\cdot1\cdot3\cdot2\cdot4}
 =40040.
\]
Thus the slice has a single reduced point, with the stated nonreduced
length. The action chart gives the four orbit directions back.
Faithful flatness of local completion descends the radical equality
to the algebraic germ. Concretely, any element in the local orbit
ideal has a power in the completed prefix ideal, hence a power in
the prefix ideal by contraction. The other containment follows
from purity of the orbit. On the open where the orbit is closed this
gives a Zariski neighbourhood of each of its points with no additional
reduced components. Orthogonal and scalar transport gives the same
result throughout \(\mathcal O\).

On \(S\), fifteen equations in a regular local ring of dimension
sixteen leave positive local dimension. Since adding \(\mu_{16}\)
leaves dimension zero, some component through \(P_*\) is not
contained in \(\{\mu_{16}=0\}\). This gives points with
\(\mu_1=\cdots=\mu_{15}=0,\ \mu_{16}\ne0\) in every Zariski
neighbourhood of \(P_*\). The slice and equations are rational,
and the nonempty constructible set has algebraic points by the
Nullstellensatz over \(\overline{\mathbb Q}\). The same is true
inside any neighbourhood defined over \(\overline{\mathbb Q}\).

## Replay

The recorded replay used Python 3.12.3, python-flint 0.9.0, msolve 0.6.5
and Singular 4.3.2. From the repository root, run

```sh
.venv/bin/python research/moments/vanishing/verify_gaussian_exceptional_sixteen_moments.py \
  --cross-check-singular \
  --write-receipt research/moments/vanishing/gaussian_exceptional_sixteen_moment_verification.json
```

Omitting `--cross-check-singular` replays the proof using msolve alone
for the affine chart.

The [receipt](gaussian_exceptional_sixteen_moment_verification.json)
records the source and certificate hashes and the exact bounds of
the calculation. The proof establishes the local intersection and
does not classify components elsewhere in the coefficient space.
