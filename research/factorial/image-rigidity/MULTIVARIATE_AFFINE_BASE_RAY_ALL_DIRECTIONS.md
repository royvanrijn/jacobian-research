# All-degree coprimality along every affine ray of row bases

Prescribe finitely many unequal rows by their positive integer gaps
\(d_i\), orientations, and smaller exponents
\(\ell_i(T)=a_i+s_iT\), with \(a_i,s_i\) nonnegative integers and
at least one \(s_i>0\). Let \(T\) range over nonnegative integers.
Fixed rows can have either orientation and can include private rows.
There is no restriction on coincidences among the moving
directions or on the signed total gap.

**Theorem.** For all but finitely many \(T\), the ordinary consecutive
factorial moments on this original profile are coprime in every
degree over every characteristic-zero field. At every bounded
terminal there is an explicit nonzero integer \(C_N\) and a
primitive affine form \(B_N T+A_N\), with \(A_N,B_N>0\), such
that a common factor requires \(B_N T+A_N\mid C_N\).
The integer \(C_N\) is a resultant of two explicitly specialized
polynomials, one a monomial or binomial. The construction does
not require a prime to be primitive for the degree, or any
genericity of the affine directions.

Let \(J\ge4\) be any integer at least the
[computable absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Its numerical value is unevaluated. Write \(G=\sum_{s_i>0}d_i\).
A uniform bound for the number of exceptional parameters is
\[
              \sum_{N=3}^{J-1} (N-1)\bigl(NG-c_N\bigr),        \tag{1}
\]
where \(c_N\) is the smaller of the two orientation counts
defined below. The witness excludes physical common roots.
The [primitive affine-ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_PADIC_ANGLES.md)
gives the stronger angular separation under its stated
native primitivity condition.

## Exact specializations at the first affine factor

For \(m\) strip the factorial \((m\ell_i(T))!\) from each row.
The remaining factor at index \(k\) is
\[
 \prod_{h=1}^{d_i k}\bigl(m(a_i+s_iT)+h\bigr)
 \quad\text{or}\quad
 \prod_{h=1}^{d_i(m-k)}\bigl(m(a_i+s_iT)+h\bigr).
\]
Along with \(\binom mk\), these define a polynomial
\(F_m(T,z)\) in two variables with integer coefficients.
At every nonnegative \(T\), \(F_m\) is the ordinary moment
polynomial divided by a nonzero scalar.

Fix \(3\le N<J\), and put \(n=N-1\).
Set
\[
        \tau=\max_{s_i>0}-\frac{Na_i+1}{N s_i}
                 =-\frac{A}{B},\qquad A,B>0,\quad\gcd(A,B)=1 .
\]
This is the rightmost zero of the upper moving factors. Let
\(u,v\) count its attaining rows in the increasing and decreasing
orientations, and put \(c=\min(u,v)\). Only \(h=1\) can
vanish in those rows at \(T=\tau\). All other upper factors and
every lower factor are strictly positive at \(\tau\).
Indeed, for every moving row the lower-order first zero is
strictly smaller than the upper-order first zero, and all
later offsets have smaller zeros as well. The fixed factors
are strictly positive constants.

Make the integer parameter change \(x=BT+A\). In the two
products above replace each moving factor by the integer
linear factor
\[
           m s_i x+B(ma_i+h)-m s_i A
\]
and keep the fixed factors and binomial coefficients.
Denote the resulting polynomial by \(\mathcal F_m(x,z)\).
If \(G_+,G_-\) are the respective total moving gaps,
\[
 \mathcal F_m(x,z)
   =B^{mG_-}
      F_m\!\left((x-A)/B,\ B^{G_+-G_-}z\right).
\]
The positive argument scale is the same for the two orders.
The coefficient orders at \(x=0\) of the upper array are
\(v\) at index zero, \(u\) at index \(N\), and \(u+v\)
at every interior index, exactly. Thus \(x^{-c}\mathcal F_N\)
is an integer polynomial and specializes at \(x=0\) to
\[
 \begin{cases}
   C_0 & u>v,\\
   C_1z^N & u<v,\\
   C_0+C_1z^N & u=v,
 \end{cases}
\]
where each constant appearing is positive.
The lower specialization has degree \(n\) and positive coefficients,
including its endpoints. In the unbalanced cases the two
specializations immediately have no common finite root.

## Balanced cohorts have a positive geometric specialization

If \(u=v\), there is at least one marked row of each orientation.
At order \(n\) every marked row has first Gamma shape
\(1+n\ell_i(\tau)=1/N\). For a marked row of gap \(d\),
the Gauss grouping of its factor at distance \(k\)
contains the rising factorial \((1/(Nd))_k\):
\[
         (1/N)_{d k}
           =d^{dk}\prod_{h=0}^{d-1}
                       \bigl((h+1/N)/d\bigr)_k .
\]
The same applies to \(n-k\) in the opposite orientation.
For \(0<\gamma<1\), the ratio \((\gamma)_j/j!\) is the
\(j\)-th moment of a \(\mathrm{Beta}(\gamma,1-\gamma)\)
variable. Use one such small shape from each marked orientation.
It absorbs the two factorials in \(\binom n k\), up to
the single scalar \(n!\). Each other rising factorial,
including any from the unmarked moving rows, is the moment
of a Gamma variable with a strictly positive shape.
The powers of \(d\) and \(B\) are positive multiplicative
scales in these laws. Consequently the lower specialization
has coefficients proportional to
\(\mathbb E[U^{n-k}V^k]\) for two positive variables, or,
after constant normalization, to \(\mathbb E_n[T^k]\)
under the \(U^n\) tilt, with \(T=V/U\).
The two reserved Beta factors are independent and
nonconstant, so \(T\) is nonconstant.

The
[geometric moment angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md#the-geometric-resolvent-without-reflection-symmetry)
therefore apply. The \(j\)-th positive lower argument lies in
\[
 \left(\frac{(2j-1)\pi}{n},\frac{2j\pi}{n+1}\right).
\]
This interval lies strictly between the two consecutive
upper binomial rays \((2j-1)\pi/(n+1)\) and
\((2j+1)\pi/(n+1)\). Reflection treats the lower half-plane.
The binomial has a negative real root just when
\(N\) is odd; the lower geometric moment polynomial
has a negative real root just when \(n\) is odd.
There is no common argument and hence no common root
at this balanced specialization.
The argument depends on the small virtual shapes at the
specialization and applies even if the actual ordinary
profile has no private variables.

## Native divisor certificates and the exceptional bound

Use the fixed-degree Sylvester resultant
\[
             H_N(x)=
       \operatorname{Res}^{(n,N)}_z
            \bigl(\mathcal F_n(x,z),\,x^{-c}\mathcal F_N(x,z)\bigr).
\]
It is in \(\mathbb Z[x]\). In every case \(H_N(0)\ne0\):
the lower leading coefficient is nonzero at zero and
the upper specialized value is nonzero at each lower root.
The fixed-degree product formula remains valid when the
upper specializes to a constant, with the lower leading
coefficient still raised to the formal upper degree \(N\).
Thus the integer \(C_N=H_N(0)\) is obtained entirely from
the two specialized arrays.

For an original nonnegative integer \(T\), the affine
parameter \(x=BT+A\) is a positive integer. The stripped
array still has exactly its intended degree because the
original factorial coefficients are positive and \(x\ne0\).
A common root would make \(H_N(x)=0\), implying
\(x\mid C_N\) by reducing the integer polynomial modulo \(x\).
In particular any prime power dividing \(x\) and exceeding
the corresponding valuation of \(C_N\) is an exact
arithmetic certificate for that pair. It can be checked
using only the fixed specialized integer \(C_N\).
The necessary divisibility does not assert that its
divisors give common roots.

The \(T\)-degree of \(F_m\)'s coefficient at \(k\) is at most
\(mG_-+k(G_+-G_-)\), and the integer parameter and scalar
changes preserve this bound. The resultant has degree
\(N\) in the lower coefficients, degree \(n\) in the upper
coefficients, and total weighted index degree \(nN\).
The last equality follows by applying \(z\mapsto t z\) to
both polynomials and using its scaling on the resultant.
Thus the degree of the unstripped resultant in \(x\)
is at most \(nNG\). Removing \(x^c\) from the upper
reduces the resultant degree bound by \(nc\).
Each exceptional integer \(T\) gives a distinct integer
root of the nonzero polynomial \(H_N\).
Summing its degree bounds proves (1). Terminal two
is always separated by real positivity of the quadratic.
The absolute angular cutoff handles \(N\ge J\).
Coprimality transfers to all characteristic-zero fields
by the rational Bézout identity.

For example, take the two shared rows
\((T+1,T+2)\) and \((T+3,T+1)\), together with any fixed
rows. The first-block marker cohort has one row of each
orientation at every terminal, so the binomial
specialization applies at every relevant degree.
The moving portion has no exact mirror and its signed total
gap is \(-1\).
