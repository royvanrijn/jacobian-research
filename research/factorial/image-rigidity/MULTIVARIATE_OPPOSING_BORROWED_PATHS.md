# All-degree coprimality for opposing borrowed paths

Let \(J\ge4\) be an
[absolute eventual cutoff for ordinary binomials](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md).
Fix positive integers \(r,s\). There is a prime \(p\), determined
by \(J,r,s\), with the following property.

Take \(r\) paths ending at \((0,1)\) and \(s\) paths ending at
\((1,0)\). Each begins at some \((c,c)\), \(c\ge1\), and has
positive integer coordinates at all other vertices before its
endpoint. Each edge \(A\to B\) supplies the ordinary exponent row
\(p^e A-B\), with \(p^e>J\max_\nu B_\nu\). The ordinary binomial
whose exponent rows are these edges has coprime consecutive
factorial moments in every degree.

The theorem allows arbitrary internal exponents and arbitrary
lengths of the paths. The support is fully shared, and there is
no retained core. Both endpoint coefficients vanish in the
residue arrays. The interior has a strict interlacing comparison,
and direct first-edge calculations separate the two exterior
branches at adjacent orders.

At the prime chosen for \(r=s=1\), the same conclusion holds
whenever \(r,s\equiv1\bmod(p-1)\). In particular, it applies to
all ordinary unit-gap profiles with the larger exponent in each
row divisible by \(p\) and these counts in the two orientations.
Indeed, such a row is the single edge from \((c,c)\) to \((0,1)\)
or \((1,0)\), with \(e=1\). These heights may vary independently.

## Interior and endpoint coefficients

Use the [borrowed factorial identity](MULTIVARIATE_BORROW_CYCLE_ROW_TRANSFER.md#the-borrowed-factorial-identity)
at precision \(p\) with \(s=0\) in its notation. For order \(j<J\)
and coefficient index \(k\), a vertex row \(A\) gives
\(x=(j-k)A_1+kA_2\). After the rational scalar and common diagonal
rescaling, an edge with both arguments positive contributes
\((x-1)!/(y-1)!\) times a unit congruent to one modulo \(p\).
If the last argument is zero it instead contributes
\(-p^e x(x-1)!\) times such a unit, by the endpoint calculation
for [one borrowed path](MULTIVARIATE_BORROWED_JACOBI_PATHS.md#the-reduced-arrays-including-the-endpoint).
The ordinary rows have positive coordinates. Their columns are
distinct, since the signed gap of a final row is congruent to
\(1\) or \(-1\) modulo \(p\).

The products along the paths telescope. Also divide the order
\(j\) polynomial by the product of \((jc-1)!\) from all the
starting rows. Denote the resulting coefficients by \(c_{j,k}\).
At \(0<k<j\),
\[
          c_{j,k}\equiv
         \frac{\binom jk}{(k-1)!^{\,r}(j-k-1)!^{\,s}}\pmod p .
                                                               \tag{1}
\]
Take \(p>2J\), so these coefficients are all units. For each
path ending at \((0,1)\), write \(a_i\) for the first entry
of its last source vertex and \(e_i\) for its last exponent.
Put \(u=\sum_i(e_i+v_p(a_i))\) over these paths, and
\(a^\circ=\prod_i a_i p^{-v_p(a_i)}\).
At \(k=0\), these paths contribute the extra factor
\(\prod_i(-p^{e_i}j a_i)\) after telescoping.
The \(s\) other paths contribute the denominator
\((j-1)!^s\). Therefore
\[
        v_p(c_{j,0})=u,\qquad
        p^{-u}c_{j,0}\equiv
           \frac{(-1)^r j^r a^\circ}{(j-1)!^s}\pmod p.       \tag{2}
\]
The reversed formulas hold at \(k=j\), interchanging \(r,s\)
and using the second entries of the final sources of the
paths ending at \((1,0)\). This calculates exact valuations
even when some \(a_i\) or their reversed counterparts are
nonunits.

## A strict interlacing lemma for the interior

For \(m\ge0\) define
\[
 D_m(z)=\sum_{\ell=0}^m
 \frac{\binom m\ell z^\ell}
 {(\ell+1)!(\ell!)^{r-1}(m-\ell+1)!((m-\ell)!)^{s-1}} .
                                                               \tag{3}
\]
Putting \(m=j-2,\ \ell=k-1\) in (1) gives the residue
polynomial \(j(j-1)zD_{j-2}(z)\). We prove that successive
\(D_m\) have simple negative roots and strictly interlace.

For an integer \(c\ge0\), the diagonal map
\(T_c(z^k)=z^k/(k+c)!\) preserves real-rootedness and
makes every nonzero root of the result simple. For a direct
proof, the Euler map \(1+\theta/h\), with
\(\theta=z\partial_z\) and \(h>0\), preserves real-rootedness.
At a nonzero root of the input of multiplicity \(b\) it
leaves multiplicity \(b-1\). Every new nonzero root is
simple: away from the old roots and zero, its equation is
\(f'/f+h/z=0\), and this function has strictly negative
derivative. Between each two consecutive distinct poles
there is exactly one zero, since the pole residues are
positive. Together with the surviving multiple roots
and the unchanged zero-root multiplicity, these exhaust
the degree. The finite products of these Euler maps,
multiplied by the positive dilation \(z\mapsto z/N\) and
the scalar \(1/c!\), have diagonal entries tending to
\(1/(k+c)!\), since
\[
 \lim_{N\to\infty}
 \frac{\prod_{h=c+1}^{c+N}(1+k/h)}{c!\,N^k}
       =\frac1{(k+c)!}.
\]
For inputs of degree at most \(M\), set aside the first
\(M\) Euler maps. The remaining products, including
the dilation and scalar, converge coefficientwise to a
real-rootedness preserver with strictly positive diagonal
entries. This uses closure of real-rootedness at fixed
degree. The set-aside maps now remove all nonzero
multiplicities. The zero-root multiplicity is preserved.

It follows that \(T_c\) preserves strict interlacing for
adjacent-degree polynomials with simple negative roots and
positive coefficients. If the pair interlaces strictly,
every real linear combination \(f+t g\) of the monic
higher-degree and lower-degree polynomials is simple and
real-rooted, by the signs at the roots of \(g\) and at
infinity. Any zero at the origin is at most simple.
The transformed combinations are therefore all simple,
and the transformed pair cannot share a root: otherwise
one could cancel the derivative there with a choice of \(t\).
In the transformed pair the quotient of the higher-degree
polynomial by the lower-degree one has nonvanishing
derivative wherever it is finite, for the same reason.
It is thus strictly monotone in each interval between its
poles, with distinct infinite limits at the interval ends
and at the two exterior infinities. It has one zero in
each interval, which gives strict interlacing.
Reversing both polynomials preserves the strict ordering
since their roots are negative and nonzero.

The sequence \(\sum_{\ell=0}^m\binom m\ell z^\ell/(\ell+1)!\)
is \(L_m^{(1)}(-z)/(m+1)\). Its simple negative roots and
strict consecutive interlacing follow from Laguerre
orthogonality for the weight \(t e^{-t}\) on \(t>0\).
Equivalently, the coefficient expression follows by
differentiating \(e^{-t}t^{m+1}\) exactly \(m\) times
and dividing by \(m!te^{-t}\), and integration by parts
gives orthogonality. Apply \(r-1\) copies of \(T_0\).
Reverse the sequence, apply \(T_1\) and \(s-1\) copies
of \(T_0\), and reverse again. This produces exactly
(3). The strict interlacing statement follows. In
particular the rational resultants of \(D_m,D_{m+1}\)
are nonzero.

## Separation at the chosen prime

Choose \(p>2J\) avoiding the numerators of the nonzero
resultants of the relevant \(D_m,D_{m+1}\), and the
numerators of the differences
\[
       \frac{j^{r-1}}{(j-1)^s}
            -\frac{(j+1)^{r-1}}{j^s}
        \qquad (2\le j<j+1<J)
                                                               \tag{4}
\]
and their versions with \(r,s\) interchanged. The differences
are nonzero, since equality would give
\(j^{r+s-1}=(j+1)^{r-1}(j-1)^s\), incompatible with
\(j\ge2\) and coprimality of \(j\) with \(j^2-1\).
All denominators here and in (3) are units at \(p\).

Suppose an actual normalized consecutive pair at \(j,j+1\),
\(2\le j<j+1<J\), has a common root in an algebraic closure
of \(\mathbb Q_p\). At valuation zero the unit residue of
the root must be a common zero of \(D_{j-2},D_{j-1}\),
impossible by the prime choice. At positive valuation,
the only competing terms of smallest valuation are the
constant term and the term of degree one, since all
interior coefficients are units and the last coefficient
has positive valuation. Thus the root has valuation \(u\)
from (2). After scaling it by \(p^{-u}\), the necessary
residue is
\[
                    (-1)^{r+1}a^\circ
                           \frac{j^{r-1}}{(j-1)^s}.       \tag{5}
\]
This follows by using \(c_{j,1}\equiv j(j-1)^s/(j-1)!^s\)
from (1). Formula (4) prevents the residue (5) from
matching between the two orders. At negative valuation
apply the identical argument to the common reciprocal
root and the reversed arrays. The root cannot be zero
since both actual endpoints are nonzero. This proves
coprimality at all the stated pairs. At terminal two the
even moment is strictly positive on the real line,
whereas the first moment has a real root. Terminals at
least \(J\) are covered by the ordinary absolute theorem.
Rational Bézout identities transfer coprimality to any
characteristic-zero field.

Finally, if \(p\) is taken for \(r=s=1\), and the actual
path counts satisfy \(r,s\equiv1\bmod(p-1)\), the interior
coefficients in (1) agree modulo \(p\) with those for
counts \(1,1\). In (5) the profile-dependent product
\(a^\circ\) remains a unit and independent of order,
while the \(j\)-dependent factor and sign reduce to
those for \(1,1\). The same reasoning applies in the
reversed orientation. This proves the claimed unbounded
count version at one fixed prime.

The digit gap and path representation are hypotheses on
the ordinary profiles. The moment comparisons in this
proof do not supply such a representation for an
arbitrary exponent array.
