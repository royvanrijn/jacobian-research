# Effective zero-strength exclusion for full-rank affine factorial tails

Consider a fixed ordinary two-monomial core and the
[balanced affine template](MULTIVARIATE_BALANCED_AFFINE_JET_DIFFUSE_EXCLUSION.md)
with rational row slopes \(\beta_i>0\), rational offsets, and
fixed nonzero integer gaps \(d_i\). The count of a type is
\(rc_i+s_i\), where the nonnegative integers \(c_i\) are not
all zero, the \(s_i\) are fixed integers, and
\(\sum_i c_i d_i=0\). Only ordinary heights and nonnegative
integer counts are used in the conclusion. At least one moving
row must be present, and the total profile must have distinct
exponent columns. Fix adjacent orders
\(n,N=n+1\), \(n\ge2\), at which the
[affine divisor matrix](MULTIVARIATE_AFFINE_DIVISOR_RANK_DIFFUSE_EXCLUSION.md)
has full row rank.

**Theorem.** At this pair there is a computable height beyond
which all original consecutive moment pairs in this template
are coprime, uniformly over admissible integer multiplicities
\(r\). In particular the conclusion includes the
zero-strength limit even when the fixed core has a common
root. If the matrix has full row rank at every pair below
the [effective absolute cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md),
there is a single computable all-degree height threshold.
Pairs already covered by another uniform exclusion do not
require the rank test. The assertion transfers to every
characteristic-zero field.

The new input is Ax's differential theorem. It forces every
analytic branch of the common-quadratic locus at zero
strength to have constant repetition parameter. An affine
ray specialization then excludes the admissible integer
constants from these branches.

## Normalized coefficients and algebraic independence

Write \(u=1/M\), \(A=\sum_i c_i d_i^2/\beta_i>0\),
\(\tau=Aru\), and \(h=\sum_i s_i d_i\). Use any consistent
choice of the common lower endpoint scale from the affine
divisor proof. The coefficient arrays, excluding the two
constants and the lower leading coefficient from the
multiplier coordinate list, have the form
\[
             b_i\,\exp(G_i(u))\,\exp(r H_i(u)),\qquad
                      1\le i\le s=n+N-1.                   \tag{1}
\]
Here \(b_i\) are fixed nonzero rational numbers, including
the upper argument twist \((N/n)^h\). Each \(H_i,G_i\)
vanishes at \(u=0\), is analytic with rational Taylor
coefficients, and is \(1/n\) times the logarithm of a
nonzero rational function of \(u\) taking value one at zero.
For \(H_i\) these rational functions are the normalized
repeated-block coefficient products recorded by the affine
divisor matrix. For \(G_i\) they come from the fixed additive
counts. The notation \(H_i/u\) always uses the holomorphic
extension at zero. Formula (1) retains the constant and
lower leading coefficients of the normalized fixed core in
the omitted positions.

The \(H_i\) are linearly independent over \(\mathbb Q\)
modulo constants. An integral relation would make the
corresponding product of rational functions constant and
give a vector in the left kernel of the full-row-rank
divisor matrix. This argument also applies after any
nonconstant analytic substitution \(u=u(w)\) at zero,
since an identity of analytic functions after substitution
is an identity before substitution.

We use [Ax, On Schanuel's conjectures, formal
power-series statement (SP)](https://doi.org/10.2307/1970774).
For \(k\) analytic germs of one variable linearly independent
over \(\mathbb Q\) modulo constants, the field generated
by the germs and their exponentials has transcendence degree
at least \(k+1\) over \(\mathbb C\). Subtracting their
constant terms gives the formal statement with zero constant
terms. In particular, since the \(\exp H_i\) are algebraic
over \(\mathbb C(u)\), Ax gives
\[
         \operatorname{trdeg}_{\mathbb C(u)}
                  \mathbb C(u,H_1,\ldots,H_s)=s.             \tag{2}
\]
Indeed adjoining \(u\) and the algebraic exponentials
increases transcendence degree over \(\mathbb C\) by at most
one over the degree of the \(H_i\) alone. The same proof
works after substitution \(u=u(w)\). The independence in
(2) also holds over algebraic extensions of \(\mathbb C(u)\)
inside the field of meromorphic germs.

## A contact branch has constant rational repetition parameter

Regard the \(s\) normalized multipliers as coordinates in
a complex torus. The locus of coefficient arrays with a
common divisor of degree at least two has codimension at
least two there. For fixed lower coefficients there are
only finitely many monic quadratic divisors of the lower
polynomial. For each one, requiring the upper polynomial
to vanish modulo that divisor gives two independent linear
conditions in its free coefficients: the classes of \(z,z^2\)
are independent modulo a monic quadratic with nonzero
constant term. Intersecting with the multiplier torus and
taking the finite incidence over lower coefficients retains
the stated codimension. The locus is defined by the
corank-two Sylvester minors, including the resultant if
desired. A physical common root of the ordinary arrays
belongs to this locus because it is nonreal, by strict even
moment positivity, and its conjugate is also common.

Use the analytic variables \(u,\tau\) in (1), replacing
\(rH_i\) by \((\tau/A)(H_i/u)\). The reference resultant
on \(u=0\) is a nonzero analytic function of \(\tau\) by the
[Gaussian reference comparison](MULTIVARIATE_AFFINE_DIVISOR_RANK_DIFFUSE_EXCLUSION.md#algebraic-projection-of-a-common-factor).
Near a finite contact at \(u=0\), the analytic zero set of
the resultant therefore has finitely many local branches
parametrized after finite ramification by analytic germs
\(u=u(w)\ne0,\ \tau=\tau(w)\). One can obtain this description
without choosing individual multiple roots at \(u=0\):
use a small circle in \(\tau\) to construct the monic
polynomial of the enclosed roots, take its squarefree part
over the meromorphic germs in \(u\), and shrink the punctured
disk to avoid its discriminant zeros and poles. The simple
root branches on the punctured disk have finite permutation
monodromy. A finite power substitution makes them
single-valued, and boundedness of the roots gives analytic
extension to zero. Any accumulating sequence with all the
corank-two minors zero lies on such a branch of their common
zero set, by the analytic identity theorem and the finite
number of branches.

Let one of these common branches be given, and set
\(r=\tau/(Au)\), a meromorphic function of the branch
parameter. Both \(H_i\) and \(rH_i\) are analytic. The
codimension bound gives
\[
 \operatorname{trdeg}_{\mathbb C}
     \mathbb C\bigl(\exp(G_i+rH_i):1\le i\le s\bigr)\le s-2.
                                                               \tag{3}
\]
Suppose first that \(r\) is nonconstant. Any nonzero rational
linear relation modulo constants among the \(2s\) germs
\(H_i,rH_i\) has the form
\[
         a\cdot H+r\,b\cdot H=c,\qquad b\ne0,
         \qquad a,b\in\mathbb Q^s,\ c\in\mathbb C.
\]
It gives \(r=(c-a\cdot H)/(b\cdot H)\).
These two affine linear polynomials in \(H\) are
nonproportional, since \(r\) is nonconstant, and are
coprime in the polynomial ring over \(\mathbb C\).
By algebraic independence of the \(H_i\), equality with the
ratio from a second relation gives a polynomial identity.
Unique factorization forces the numerator and denominator of
the second ratio to be the same scalar multiples of the first.
The scalar is rational by comparing a nonzero coordinate of
\(b\). Hence the space of relations has dimension at most one.

If there is no relation, Ax gives transcendence degree at least
\(2s+1\) for \(H,rH,\exp H,\exp(rH)\).
On the other hand they are generated up to algebraic
extension by \(u,H,r\) and the exponentials in (3),
because \(\exp H,\exp G\) are algebraic over \(\mathbb C(u)\).
The upper bound is \(1+s+1+(s-2)=2s\), a contradiction.
If there is one relation, \(r\) already lies in \(\mathbb C(H)\).
At least \(2s-1\) of the germs are independent modulo
constants, so Ax gives a lower bound \(2s\), whereas the
upper bound is \(1+s+(s-2)=2s-1\). Thus \(r\) is constant.
If that constant is not rational, independence of the
\(H_i\) again makes the \(2s\) germs independent modulo
constants: a relation would require \(a_i+rb_i=0\) for
all \(i\). The same transcendence bounds exclude it.
Thus every common branch has constant rational \(r\).
In particular all such branches meeting \(u=0\) at finite
strength lie on lines \(\tau=Aru\) through zero strength.

If a branch has constant admissible integer \(r\), the
multiplicities \(rc_i+s_i\) are fixed. The exact arrays along
that branch are the original affine-ray factorial arrays
up to scalar and common argument normalization. Their
resultant cannot vanish identically. The
[all-directions affine-ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md)
excludes accumulation at fixed row multiplicities. For
rational row slopes and offsets, clear denominators by
restricting the integer height to any admissible residue
class and translate the integer parameter to make the
smaller-exponent intercepts nonnegative. The theorem applies
on each such class. If no admissible residue class exists,
the candidate branch cannot contain the original profiles.
A constant branch with no moving row, or with coincident
monomials, is outside the stated hypotheses. This excludes every
accumulating ordinary collision.

## Effective isolation at zero strength

Write \(R(u,\tau)\) for the normalized analytic resultant.
Expand it at \((0,0)\) and compute its first nonzero homogeneous
layer \(R_b\) of total degree \(b\). The expansion uses
rational coefficients. It terminates since \(R(0,\tau)\)
is not identically zero. Every branch with constant
\(r=c\) satisfies \(R_b(1,Ac)=0\). This is a nonzero
ordinary polynomial: a nonzero homogeneous polynomial remains
nonzero on setting its first variable to one. Its rational
roots give a computable finite candidate set.

For each such \(c\), compute the smallest vanishing order
\(e_c\) along \(\tau=Acu\) among the resultant and the
corank-two minors. To test any fixed partial derivative in
\(\tau\) on that line, substitute into (1). The derivative
becomes a polynomial in \(H\) over algebraic functions of
\(u\): the factors \(\exp(G_i+cH_i)\) are algebraic for
rational \(c\), and \(\partial_\tau\) contributes
\(H_i/(Au)\). By (2), vanishing is exactly vanishing of the
algebraic coefficient functions, tested by operations on
finite algebraic function extensions with the local branches
specified at zero. Increasing the derivative index
terminates because the original resultant is nonzero.
Thus \(e_c>0\) exactly at the genuine common lines.

Divide all the analytic equations by
\[
                           P(u,\tau)=
               \prod_{c:e_c>0}(\tau-Acu)^{e_c}.
\]
This is analytic division by common linear factors.
After division, at least one equation is not identically
zero on each removed line, by the minimum in the definition
of \(e_c\). Any other common curve of the reduced system
would be a common branch of the original system. There is
no vertical component \(u=0\) because \(R(0,\tau)\) is
not identically zero. Hence the reduced system has an
isolated common zero or an empty common zero germ at the
origin.

For clarity, an effective punctured exclusion radius can be
obtained directly from this isolation. Use a small rational
circle around zero in \(\tau\) on which the reduced
resultant at \(u=0\) is separated from zero, using its first
nonzero Taylor coefficient there. Cauchy bounds on the
normalized coefficient functions and on the divided equations
give a rational small \(u\)-disk preserving this separation.
Division by a linear factor can be bounded on smaller
polydisks by the analytic identity
\[
 \frac{f(u,\tau)}{\tau-Acu}
     =\int_0^1 f_\tau\bigl(u,Acu+t(\tau-Acu)\bigr)\,dt
     \quad\text{when }f(u,Acu)=0 .
\]
Multiply a formal weighted sum of the reduced minors over the
enclosed zeros of the reduced resultant, counting
multiplicities. The product is polynomial in the weights
and analytic in \(u\). Its coefficients are obtained by
Newton identities from contour power sums, as in the
[exact contour eliminant](MULTIVARIATE_ISOLATED_TWO_CLOCK_DIFFUSE_EXCLUSION.md#exact-elimination-of-an-accumulating-contact).
At \(u=0\) each Taylor coefficient in \(u\) is an effective
residue in \(\tau\), with algebraic coefficients. Isolation
ensures this product is not identically zero; the empty
root case gives the product one. A first nonzero coefficient
of one of its weight monomials is found by finite search.
The corresponding Cauchy remainder bound then gives a
computable punctured disk on which the reduced equations do
not vanish together. Every removed line has a repetition parameter that fails the
ordinary hypotheses, by the affine-ray argument.
This provides the desired uniform zero-strength height bound.

At strictly positive strengths the effective
[affine divisor-rank theorem](MULTIVARIATE_AFFINE_DIVISOR_RANK_DIFFUSE_EXCLUSION.md)
already gives computable neighborhoods of the finitely many
reference contacts. The absolute curvature bound bounds
strength at a collision and coefficient continuity excludes
the remaining compact set away from all these neighborhoods.
Taking maxima over the finitely many terminals below the
absolute cutoff completes the effective theorem.

In particular, the
[overlapping four-type affine family](MULTIVARIATE_OVERLAPPING_AND_SEPARATED_PAIR_DIVISOR_RANK.md)
has the required full rank at every pair. Its all-degree
threshold therefore holds for every admissible multiplicity
with an arbitrary fixed ordinary core, including cores whose
own consecutive moments are not known to be coprime.
