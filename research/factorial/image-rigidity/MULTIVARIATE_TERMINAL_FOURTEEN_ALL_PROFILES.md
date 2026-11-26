# Terminal fourteen for every ordinary two-monomial profile

For distinct nonnegative integer exponent columns \(a,b\) in any
finite dimension, write
\[
            Q_j(z)=\sum_{k=0}^j\binom jk
                         \prod_i((j-k)a_i+kb_i)!\,z^k .
\]

**Theorem.** The order-thirteen and order-fourteen polynomials
\(Q_{13},Q_{14}\) have disjoint root arguments. This remains
true under independent positive rescalings of the variable.
In particular they are coprime over every characteristic-zero field.

Delete equal rows, which give only a nonzero scalar at each
order. The [retained-profile theorem](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_PROFILES.md)
covers every profile with a private row and every fully shared
profile having some reciprocal \(d/(13\ell+1)>1/8\), with
\(\ell\) the smaller row exponent and \(d\) its gap.
For the remaining fully shared profiles put
\[
       u^2=\prod_i\prod_{h=1}^{d_i}
                   \frac{13\ell_i+h}{13\ell_i+12d_i+h}.
\]
The distinct columns imply \(0<u<1\). All reciprocals in
this case are at most \(1/8\).

| Range of \(u\) | Separation proof |
| --- | --- |
| \((0,.02]\) | [Small-outer-mean rays](MULTIVARIATE_TERMINAL_FOURTEEN_SMALL_OUTER.md) |
| \([.02,.04]\) | [Paired rays and retained profiles in the low corner](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_CORNER.md) |
| \([.04,.25]\) | [Direct small-reciprocal bands](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_LOW_U.md) and their [retained complements](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_LOW.md) |
| \([.25,.50]\) | [Left flank](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID_LEFT.md), [middle](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID.md) and [right flank](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_MID_FLANKS.md) |
| \([.50,.64]\) | [Direct disks](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_LOW_MID.md) and [retained complements](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_LOW_MID.md) |
| \([.64,.76]\) | [Middle disks](MULTIVARIATE_TERMINAL_FOURTEEN_SHARED_DISKS.md) and [retained complements](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_MIDDLE.md) |
| \([.76,.855]\) | [Direct upper band](MULTIVARIATE_TERMINAL_FOURTEEN_HIGH_U.md) and [retained bridge](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_BRIDGE.md) |
| \([.855,1)\) | [Direct high-mean bands](MULTIVARIATE_TERMINAL_FOURTEEN_HIGH_U.md) and [upper retained complements](MULTIVARIATE_TERMINAL_FOURTEEN_RETAINED_UPPER.md) |

In the \([.04,.25]\) range, compare the two cited tables with
\(c_{\max}=\max_i d_i/(13\ell_i+1)\). The direct cap is \(1/8\)
except for the bands listed with smaller caps. The retained
tables begin at those same caps on \([.10,.11]\), \([.12,.16]\)
and \([.16,.21]\). On \([.21,.22]\) the intermediate strip
\(1/14<c_{\max}\le3/40\) connects the direct cap \(1/14\) to
the retained band above \(3/40\); on \([.22,.25]\) the
retained band begins at \(1/14\). Thus all reciprocal
values are covered, with equality included at the direct
caps and at the upper endpoints of the retained bands.

All bands share their closed interior endpoints. Since
both moments have nonzero endpoints, disjoint complex root
arguments imply coprimality. The integer coefficient arrays
have a rational Bézout identity, which transfers the
conclusion to every characteristic-zero field.
