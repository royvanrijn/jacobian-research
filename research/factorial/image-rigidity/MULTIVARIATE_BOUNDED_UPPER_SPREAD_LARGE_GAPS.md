# All-degree coprimality for large gaps with bounded upper spread

In an ordinary two-monomial profile, discard equal rows and let
\(M\) be the largest remaining exponent. Fix integers \(A\ge1\)
and \(B\ge0\). Select a nonempty group of \(R\) rows with
gaps \(d_i\in[D,AD]\) and larger exponents \(M_i=M-e_i\),
where \(0\le e_i\le BD\). Let \(G\) be the sum of absolute
gaps in the remaining rows.

**Theorem.** There are computable \(\delta,D_0>0\), depending
only on \(A,B\), such that \(D\ge D_0\) and \(RG\le\delta D\)
imply ordinary SFC(2) in every degree over every
characteristic-zero field. All certificates use the original
profile. The selected factorial ramps may occupy disjoint
intervals and the selected larger exponents need not include \(M\).

The proof constructs a progression of depth marks from the
union of the selected ramps. The product estimate is the
one in the
[common-upper-exponent theorem](MULTIVARIATE_COMMON_UPPER_LARGE_GAPS_ALL_DEGREES.md#bounding-the-possible-height).
The change is that a mark may lie below the beginning of any
individual ramp.

## The progression before the first ramp ends

Take \(J\ge4\) at least the
[effective absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
The [curvature twelve theorem](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
allows us to restrict to total curvature below \(12\).
Since a selected row has gap at least \(D\) and larger
exponent at most \(M\), this gives \(M>D^2/12\).

Fix a terminal \(4\le m<J\), and put \(m'=m-1\).
For \(u\ge0\), the proposed mark is \(x_u=mM-u\). Set
\[
 U_0=m\min_i(e_i+d_i),\qquad
 L_\pm(u)=\max\left(\{-1\}\cup
    \left\{\frac{u-me_i}{d_i}:i\text{ has orientation }\pm\right\}\right)
 \quad(0\le u\le U_0).
\]
At \(u=0\) the sum \(L_+(u)+L_-(u)\) is at most zero.
At \(U_0\) a row attaining the minimum contributes \(m\) to
one of the two maxima, so the sum is at least \(m-1\).
The two maxima are nondecreasing, continuous and each is
\(1/D\)-Lipschitz. Therefore the set
\[
                         m-2\le L_+(u)+L_-(u)<m-1
\]
contains a half-open interval of length at least \(D/2\)
inside \([0,U_0)\). Take the integers \(u\equiv1\pmod m\)
in that interval. They form a progression with
\(k\ge D/(2m)-1\) terms.

To check the marked coefficient indices, a selected row
contains \(x_u\) somewhere in its full factorial ramp exactly
when \(0\le u-me_i<md_i\). The upper inequality always
holds for \(u<U_0\). A row where \(u-me_i<0\) has not
reached the mark. Thus the increasing marked indices are
\(k\ge m-\lfloor L_+(u)\rfloor\) and the decreasing
marked indices are \(k\le\lfloor L_-(u)\rfloor\).
A floor of \(-1\) gives an empty set in the relevant
orientation. The two floors sum to \(m-3\) or \(m-2\),
so the common unmarked interval has width one or zero.

## The uniform offset bounds

Define \(K=J(A+B)\) and \(H=K^2(A+B)\).
The progression lies in \(0\le u<U_0\le m(A+B)D\), so
its span \(W\) and its count are at most \(KD\). For
\(D\ge48K^2\) it has at least \(D/(4K)\) terms, all with
\(x_u>m(M-(A+B)D)\ge M\). Also \(x_u\) is coprime to \(m\).

At the marked order the other factors of selected row \(i\)
differ from \(x_u\) by the offsets \(u-me_i-h\), where
\(0\le h<md_i\). At the other order use
\[
 m(m'M_i-h)=m'x_u+(m'u-mm'e_i-mh),
                    \qquad 0\le h<m'd_i .
\]
The latter offset is nonzero since \(m\nmid u\) and
\(\gcd(m,m')=1\). The nonzero magnitudes in both offset
families are at most \(HD\), and any given magnitude
occurs at most twice per row in either family.

For each \(p\mid x_u\), let \(U,L\) bound the two cleaned
valuation arrays using these offsets and the full retained
factorial blocks, including \(m!\) and \(m'!\).
If \(v_p(x_u)>U+aL\), where \(a\) is the greatest distance
to the unmarked interval, all the offsets stabilize and the
other Newton polygon has slopes in \([-L,L]\).
The [deep-block degree bound](MULTIVARIATE_CENTRAL_DEEP_PRIME_DIVISOR_GATE.md#marks-in-arbitrary-blocks)
then excludes any common factor, since the surviving
interval has width at most one. This comparison also applies
when several copies of the mark occur at an index.
Here \(1\le a\le m\), since the unmarked interval may meet
a terminal endpoint.

Consequently a common root forces the least common
multiple bound
\[
             \operatorname{lcm}_{u}x_u
                \mid B_mB_{m'}^m((HD)!)^{2R(m+1)},
 \qquad
             B_j=j!\prod_{i\ {\rm retained}}
                       \prod_{h=1}^{j d_i}(j\ell_i+h),
\]
where \(\ell_i\) is the retained smaller exponent. If
\(p>HD\) divides a mark, it must divide \(B_mB_{m'}\),
since all selected offsets are \(p\)-units.

## Completion of the product comparison

The [height estimate](MULTIVARIATE_COMMON_UPPER_LARGE_GAPS_ALL_DEGREES.md#bounding-the-possible-height)
and [prime deletions](MULTIVARIATE_COMMON_UPPER_LARGE_GAPS_ALL_DEGREES.md#a-prime-that-excludes-contact)
now apply with that proof's \(J\) replaced by \(K\), using the
offset parameter \(H\) above. Their inputs are \(x_u>M>D^2/12\),
a progression step \(m\le K\) coprime to all terms, span and
count at most \(KD\), count at least \(D/(4K)\), and retained
factors at most \(KM\).
The larger factorial powers in the least common multiple bound fit
the same estimates: \(m+1\le K\) bounds the total order-factorial
power and \(m+m'm=m^2\le K^2\) bounds the retained height exponent
per unit of \(G\). The group factorial has power at most \(2RK\).
One may use
\[
 \begin{split}
 C&=8K(2K+4KH+2K^2),\qquad \delta=\frac{1}{128K^2(C+1)},\\
 D_0&=\max\{H,48K^2,292K,2^{1024KH}\}.
 \end{split}
\]
This excludes a common root at every terminal in the indicated
range. Terminal two follows from even positivity and terminal
three from
[quadratic nondivisibility](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients).
The absolute cutoff covers the higher terminals. A rational Bézout
identity transfers coprimality to any characteristic-zero field.
