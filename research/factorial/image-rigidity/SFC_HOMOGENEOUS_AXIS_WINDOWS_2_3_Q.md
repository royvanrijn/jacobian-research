# A shifted degree sieve at factorial orders two, three and \(q\)

Use a fixed homogeneous trinomial axis profile: column \(i\)
has exponents \((u_i,D-u_i)\) or \((D-u_i,u_i)\), with
\(u_i\) a fixed nonnegative integer and \(D\) the degree.
The assignment of a column to a side is fixed. Offsets on the
same side are distinct.

**Theorem.** For every integer \(q\ge4\), the ordinary degrees
at which the window of orders \(2,3,q\) can have a nonzero
common factorial zero lie in a computable finite set.
This includes the consecutive window starting at two.
The conclusion holds over every characteristic-zero field.
There is no degree bound in the input and the three fixed
axis offsets are arbitrary.

The proof uses the
[fractional degree sieve for orders \(1,2,3\)](SFC_HOMOGENEOUS_TRINOMIAL_AXIS_DEGREE_SIEVE.md)
with the pole now at \(D_*=\max_i u_i-1/q\).
When the maximal offset occurs on both sides, the pole leaves
a root-of-unity condition that the second and third moments
exclude.

## 1. The fractional specialization

Use \(M=\max_i u_i\) and the decoupled forms \(P_m(D)\)
of the initial axis sieve. A term with count vector \(\alpha\)
on one side has coefficient
\[
             \binom m\alpha\frac{U!}{(mD)_{\underline U}},
             \qquad U=\sum_i\alpha_i u_i .
 \tag{1}
\]
Mixed terms between sides are omitted. At \(D_*=M-1/q\),
there is a simple pole in \(P_q\) precisely on the pure terms
of the columns with offset \(M\): the final denominator
factor at \(U=qM\) vanishes there. Their residues are the same.
All coefficients of \(P_2\) and \(P_3\) are regular and
positive whenever present. Indeed for \(m=2,3\) the last
factor for \(U\le mM\) is at least \(1-m/q>0\).
After multiplying \(P_q\) by \(D-D_*\), its specialization
is a nonzero multiple of the sum of the maximal pure
\(q\)-th powers.

If the maximum is unique, only the other two columns can
remain. On opposite sides their second and third forms
are diagonal with positive coefficients. A ratio with
negative real square cannot have negative real cube.
On the same side the surviving columns have both entries
nonnegative at \(D_*\), since their offsets are at most
\(M-1\). By the beta integral, the two forms are positive
multiples of the ordinary product-exponential moments of
the two real-power monomials with total degree \(D_*\).

Here the pair has no common second and third moment zero.
To see this, use the
[boundary pair argument](SFC_UNIFORM_HOMOGENEOUS_TRINOMIAL_INITIAL_CORNER.md#3-boundary-clusters-and-the-final-two-orders)
with \(h_m(t)=m^{-1}\log\Gamma(mt+1)\) added over the two
exponent rows. The difference \(h_3-h_2\) is strictly
convex by the trigamma integral. The Jensen gap for a
one-third or two-thirds mixture, multiplied by three,
strictly exceeds twice the midpoint gap for \(h_2\).
Thus both normalized mixed cubic coefficients are smaller
than the normalized mixed quadratic coefficient. As in
the boundary pair proof, comparing the sum of the mixed
coefficients in a putative quadratic-times-linear
factorization gives a contradiction by the
arithmetic-geometric mean inequality. The argument uses
nonnegative real exponents and does not require them to
be integers.

## 2. Two maximal columns

If the maximum is twice attained, put the smaller offset
\(a\) and the larger \(b=M\) on one side and the equal
maximal offset on the other. Use coefficients \(X,Y,Z\).
At \(D_*\), the second and third forms have the shape
\[
 \begin{split}
 Q_2&=A X^2+2L XY+B(Y^2+Z^2),\\
 Q_3&=C X^3+3U X^2Y+3V XY^2+E(Y^3+Z^3),
 \end{split}
 \tag{2}
\]
with positive rational coefficients. The marked equation
is \(Y^q+Z^q=0\). If \(Y=0\), both \(Z\) and then \(X\)
vanish. Otherwise put \(Y=1,\ Z=\xi\), with \(\xi^q=-1\).

Write \(\delta=b-a>0\) and \(t=D_*-a=\delta-1/q\).
By (1), the ratios of monomial moments obtained by
replacing exponent \(ma\) by \(ma+h\) are
\[
                   r_m(h)=\prod_{j=1}^h
                              \frac{ma+j}{mt-j+1}.
\]
In particular \(L/A=r_2(\delta)\), \(U/C=r_3(\delta)\),
\(B/A=r_2(2\delta)\), and \(V/C=r_3(2\delta)\).
The quotients of the individual factors for \(m=3\) and
\(m=2\) are positive, less than one and strictly decreasing
with \(j\) for \(1\le j\le2\delta\). Indeed they are
\((3a+j)/(2a+j)\) times \((2t-j+1)/(3t-j+1)\).
The first factor is nonincreasing, the second is strictly
decreasing, and the difference between the product's
denominator and numerator is \(jD_*-a>0\).
We will use
\[
 \lambda=\frac{U/C}{L/A}\in(0,1),\qquad
 \mu=\frac{V/C}{B/A}\in(0,\lambda^2).                         \tag{3}
\]
Put \(x=L/\sqrt{AB}\in(0,1)\) by the strict square-integral
inequality on the two left columns, and
\(w=(E/C)(A/B)^{3/2}>0\).
In the variable \(z=X\sqrt{A/B}\), equations (2) become
\[
 z^2+2xz+S=0,\qquad z^3+3x\lambda z^2+3\mu z+wT=0,
            \qquad S=1+\xi^2,\quad T=1+\xi^3 .
 \tag{4}
\]

### Even \(q\)

If \(\xi^2=-1\), the quadratic in (4) forces \(z\) to be
real, but the cubic has nonzero imaginary part.
Otherwise \(\xi^2\) is nonreal, since the real choices
are \(1,-1\) and the first is incompatible with even \(q\).
The order of \(\xi\) is a multiple of four, so the
cyclotomic automorphism \(\xi\mapsto-\xi\) fixes
\(\xi^2\). Work in the rational normalization (2).
If \(s=X/Y\) is a common root, the cubic shows
\(\xi^3\in\mathbb Q(\xi^2,s)\), hence \(\xi\) belongs
to that field. The quadratic gives degree at most two
over \(\mathbb Q(\xi^2)\), so this field equals
\(\mathbb Q(\xi)\). The conjugate of \(s\) is the other
quadratic root.

Put \(\ell=2L/A\) and \(b_*=B/A\).
The two roots have traces of their first, second and
third powers
\(-\ell,\ \ell^2-2b_* S,\ -\ell^3+3\ell b_* S\),
respectively. The trace of the cubic equation is
therefore a rational constant plus \(3b_*(C\ell-2U)S\).
The coefficient is positive by (3), and \(S\) is
nonreal. This trace cannot vanish.

### Odd \(q\)

At \(\xi=-1\), (4) would give a shared root of the real
quadratics \(z^2+2xz+2\) and \(z^2+3x\lambda z+3\mu\).
The first has conjugate nonreal roots, so the quadratics
would coincide. This would require
\(\lambda=\mu=2/3\), contrary to (3).

For the other roots, eliminate \(z^2\) in the cubic of
(4). The remainder is \(Pz+Q\), where
\[
  k=x(2-3\lambda),\qquad p_0=2xk+3\mu,\qquad
  P=p_0-S,\qquad Q=k S+wT .
\]
Hence any common root requires
\[
                         R(\xi)=Q^2-2xP Q+S P^2=0.          \tag{5}
\]
In the original rational normalization the same
elimination gives a rational polynomial proportional
to \(R\). For example, if \(P_0 s+Q_0\) is the
unnormalized cubic remainder, its eliminant is
\(Q_0^2-\ell P_0 Q_0+b_* S P_0^2\), unchanged up to a
constant by the real scalings leading to (4).
The degree of \(R\) is six. Its coefficients at degrees
six, five and one are, respectively,
\(1+w^2,\ 2w(k+x)>0,\ 0\).

The order of \(\xi\) is \(2r\) for an odd divisor \(r\)
of \(q\). If \(r>1\), the only possibilities with
\(\varphi(r)\le6\) are \(3,5,7,9\).
For \(r=7\), the cyclotomic polynomial of degree six
has a negative degree-five coefficient relative to its
leading coefficient. For \(r=9\), that coefficient
is zero. Neither divides \(R\).

For \(r=3\), \(S=\xi\) and \(T=0\). The cubic gives
\(z^2+3x\lambda z+v=0\), with \(v=3\mu>0\), and
(4) gives \(z=(v-\xi)/k\), since \(k=0\) is immediately
impossible. The real cubic factor is quadratic after
division by \(z\), so its two roots are \(z,\bar z\).
Their product and sum give
\[
 v^2-v+1=v k^2,\qquad 1-2v=3x\lambda k .
\]
The first equation requires \(k^2\ge1\).
For \(\lambda\ge1/3\) we have \(|k|<1\).
Otherwise substituting the second equation into the
first gives a right side
\(v(1-2v)(2-3\lambda)/(3\lambda)<1/3\),
since \(0<v<3\lambda^2\). The left side is at least
\(3/4\). Thus \(r=3\) is excluded.

For \(r=5\), divisibility by
\(\xi^4-\xi^3+\xi^2-\xi+1\) and the zero coefficient
of \(\xi\) force the coefficient of \(\xi^3\) to be the
negative of the leading coefficient. In fact if the
quotient is \(u_2\xi^2+u_1\xi+u_0\), the equation at
\(\xi\) gives \(u_1=u_0\), and the coefficient of
\(\xi^3\) is then \(-u_2\). Write \(J=p_0-1\).
In (5) this necessary equality reads
\[
             2w(k+w-xJ)=-(1+w^2).
\]
It implies \(xJ-k=\tfrac32 w+1/(2w)\ge\sqrt3\).
On the other hand (3) gives
\[
 xJ-k
 \le x\bigl((2-3\lambda)(2x^2-1)+3\lambda^2-1\bigr)
 \le \max\{x(4x^2-3),\,x(3-2x^2)\}\le\sqrt2 .
\]
The second inequality uses convexity in \(\lambda\)
on \([0,1]\). This excludes \(r=5\), completing the
tied-maximum case for all \(q\ge4\).

## 3. The effective ordinary degree test

The specialization of \(P_2,P_3,(D-D_*)P_q\) has no
common projective direction by Sections 1–2.
Thus the decoupled resultant is a nonzero rational
function of \(D\). The exact ordinary moment divided by
\((mD)!\) differs from \(P_m\) only in terms with at least
one column from each side. At \(D\ge2qM\) their two
factorial arguments are at least \(D/2\). Their normalized
coefficients are bounded by \(q!2^{-\lfloor D/2\rfloor}\),
as in the [initial sieve](SFC_HOMOGENEOUS_TRINOMIAL_AXIS_DEGREE_SIEVE.md#2-returning-to-the-ordinary-mixed-axis-moments).
All coefficients stay bounded. The universal
homogeneous resultant is a fixed polynomial in them,
so the exact resultant differs from the decoupled one
by an effectively bounded multiple of this exponential.
A nonzero rational function has an effective
inverse-power lower bound past a computable threshold.
Comparison bounds all possible ordinary failures and
leaves finitely many exact moment tests. Nonvanishing
for complex projective coefficients gives nonvanishing
over every characteristic-zero field by the
rational-coefficient Nullstellensatz.

For exact verification without constructing the universal
resultant, clear each order's row denominator by
\((mD)_{\underline{mM}}\). In degree \(q+3\), form the
linear map whose columns are these forms times all
monomials of complementary degree. At \(D_*\) its
rank is full, since three homogeneous forms with
only the zero affine solution are a regular sequence
and the quotient has top degree
\((2-1)+(3-1)+(q-1)=q+2\).
Select a maximal minor nonzero there. Its determinant
is a nonzero integer polynomial in \(D\) for the
decoupled forms. The perturbation from mixed terms
is exponentially small times a power of \(D\), which
gives the same finite-degree test.

## 4. Finite height lists for fixed homogeneous supports

**Corollary.** Fix any bivariate homogeneous trinomial offset
support and any positive multiple \(q\) of four. There are only
finitely many independent nonnegative common translations
at which orders \(2,3,q\) can have a common factorial zero.
The exact list of failed integer height pairs is computable.

The collinear quadratic Gaussian forms at these orders
have no common direction. The
[even-order valuation argument](SFC_COLLINEAR_TERNARY_EVEN_START_AND_DYADIC_REDUCTION.md#even-orders-and-distinct-valuations)
uses an odd order and two even orders with distinct
two-adic valuations, as is the case here. By the
[independent-height Gaussian normalization](SFC_BIVARIATE_TRANSLATION_GAUSSIAN_JETS.md#1-exact-graded-normalization),
all possible failures are therefore on finitely many
computable bounded-height rows and columns. On each such
line one exponent row is fixed. The support there is a
one-side axis profile, and the theorem gives a computable
finite list of degrees, hence of heights on that line.
Test the original equations at those heights. In particular
there is no infinite exceptional height line for the
consecutive window starting at two.

## Exact checks

The [Macaulay checker](verify_sfc_homogeneous_axis_windows_2_3_q.py)
chooses a full-rank minor at the fractional mark. For a
one-side profile the decoupled forms are the ordinary forms,
so the checker extracts all possible ordinary integer
roots of that minor and tests the actual rank at each.
It verified all degrees for \(q=4\) on the one-side
profiles \((0,1,2),(0,1,3),(0,2,3),(0,2,4)\).
For the tied mixed profile \((0,1)\) opposite \(1\),
the fractional mark had full rank at
\(q=4,5,6,7,9,15\). In mixed-side mode this checker
tests the fractional mark only.

~~~sh
python3 research/factorial/image-rigidity/verify_sfc_homogeneous_axis_windows_2_3_q.py 0 1 2
python3 research/factorial/image-rigidity/verify_sfc_homogeneous_axis_windows_2_3_q.py 0 1 1 --opposite-third --terminal 15
~~~
