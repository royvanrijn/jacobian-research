# Quartic affine-source components through ordinary degree seven

**Theorem KellerQuarticAffineSourceReductionThroughDegreeSeven.**
Let \(F=(P,Q,R)\) be a noninvertible three-variable Keller map of
ordinary total degree at most seven over a field of characteristic zero.
Suppose its quartic component is \(R=a(x,y)z+b(x,y)\) with \(a\ne0\).
Over the algebraic closure, affine changes of coordinates and an
affine rescaling of \(R\) put the scalar in one of these forms:
\[
\begin{aligned}
  &g(x)+x^2(z+xy),\\
  &g(x)+x^2(z+y^2),\\
  &g(x)+x^2(xz+y^2),\\
  &x(xy-1)z+\lambda_0y+\mu_0+
          xy\big((\lambda_1-\lambda_0)y+\mu_1-\mu_0\big),
          \qquad \lambda_0\lambda_1\ne0 .
\end{aligned}
\tag{1}
\]
In the first three lines \(g'(0)\ne0\).
At degree seven there is also the possible parallel scalar
\(R=g(x)+A(x)z\) with \(A\) nonconstant; at degree at most
six that scalar is excluded. The affine degree bounds
for the other two components are unchanged.

**Theorem SexticTripleLineShearQuadraticMateExclusion.**
For the first scalar in (1), set \(v=z+xy\) and use \(y\) as the
free variable. If \(P,Q\) are of ordinary degree at most six
in the original affine coordinates, neither mate can have
\(y\)-degree at most two in the sheared coordinates.
With degree at most seven, if one of the mates is quadratic
in that variable, the free-variable degrees are \(3,2\).
After making the scalar \(r=x+x^2v\) as below, the leading
coefficients are nonzero scalar multiples of \(rx^3,rx^2\),
respectively.

At ordinary degree at most six, the
[triple-line exclusion](SEXTIC_KELLER_TRIPLE_LINE_SEXTIC_QUARTIC_MATES.md)
and the [parabolic exclusion](KELLER_PARABOLIC_SHEAR_QUADRATIC_MATES.md#the-full-sextic-exclusion)
remove the first two forms in (1). The
[tangent-conic exclusion](SEXTIC_KELLER_TANGENT_CONIC_HIGH_MATES.md#exclusion-of-the-tangent-conic-quartic-form)
removes the third form. For the fourth form the
[squarefree mate restrictions](KELLER_SQUAREFREE_QUARTIC_FIBRES_AND_MATES.md)
require both mates to depend on \(z\) and at least one to have
\(z\)-degree at least two. For sextic maps the
[high source-degree exclusion](SEXTIC_KELLER_SQUAREFREE_MATE_DEGREES.md)
also bounds both \(z\)-degrees by three and gives the
two possible leading coefficient forms at \(z\)-degree three.

**Corollary SexticKellerQuarticAffineSourceSquarefreeResidual.**
A noninvertible characteristic-zero polynomial Keller map of ordinary
degree at most six with a quartic component affine in a source
direction can be put in the fourth form in (1). Both mate degrees
in that source direction are in \([1,3]\), and at least one is at
least two. The cubic leading restriction above applies to each mate
of degree three in that direction.

## Affine normal forms and the completion bounds

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. The residue arguments and cubic coefficient
classification in the
[quartic affine-source classification](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md)
use the mates for their polynomial primitive on the fibres.
They apply to the present degree bound. The propagating
degree bounds for the coordinate-completion cases are as follows.
A reduced coefficient of degree at most two has an
[inverse completion of degree at most eight](KELLER_QUARTIC_TRIPLE_BOUNDARY_AND_LINEAR_COMPONENTS.md#invertibility-with-a-reduced-quadratic-coefficient).
The cubic graph case has inverse degree at most six.
For a parallel coefficient with no constant restriction of
\(b\) on a root line, the inverse parametrizations of the
individual fibres have degree at most thirteen, as in the
[parallel fibre construction](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md#coordinates-on-the-parallel-support-fibres).
Here it is the fibre degree that is bounded by thirteen.
Each plane Keller map on a fibre therefore has degree at most
\(7\cdot13=91\), so the
[bounded planar reduction](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
still gives inversion. For the residue method on a repeated
root line, the primary antecedent is
[Wang and Yan, Lemma 3.14](https://arxiv.org/pdf/2610.36495v1).

For a repeated parallel coefficient with a constant fibre
restriction, translate that root to zero. The
[repeated-line residue result](KELLER_LINEAR_SOURCE_COEFFICIENT_CUBIC_CLASSIFICATION.md#normal-residues-along-a-repeated-line)
gives \(b(0,y)=\mu\), \(b_x(0,y)=\lambda\ne0\).
Thus
\[
     a=x^2(a_0+a_1x),\qquad
     b=g(x)+x^2\big((u_0+u_1x)y+c y^2\big).
\]
If \(a\) has a second, distinct root, it is simple, and the
restriction of \(b\) there must be affine linear in \(y\).
This forces \(c=0\). When \(c=0\) and the two linear
polynomials \(a_0+a_1x,\ u_0+u_1x\) are dependent, a
constant transverse frame gives \(R=g(x)+A(x)v\) with
\(A\) nonconstant. This is excluded through degree six by
the [parallel scalar theorem](SEXTIC_KELLER_PARALLEL_COEFFICIENT_SCALAR.md).
At degree seven in the present quartic setting this
parallel scalar is an additional residual unless its
mates fall under a separate exclusion. When those
two linear polynomials are independent, a constant
transverse frame gives the first scalar of (1).
When \(c\ne0\), \(a\) is proportional to \(x^2\) or
\(x^3\). Complete the square in \(y\) by an affine
translation and rescale to get the second or third
scalar in (1).

For the irreducible squarefree cubic option,
an affine \(z\)-shear gives
\[
                  (x^2 y+\delta x-1)z+\lambda xy+\mu ,
                  \qquad \lambda\ne0 .
\]
Interchange the two transverse directions and rescale
the variable formerly called \(z\). The coefficient
in the complementary direction becomes proportional
to \(x(xy-1)\). On \(x=0\) and \(xy=1\) the
new \(b\) has the same nonzero linear slope, with
constants that may differ. Thus this case belongs
to the last line of (1). The other squarefree cubic
residual already has that coefficient and the two
specified affine restrictions. Interpolating between
the two components gives the displayed representative
of degree at most three. In both reductions the
difference from the original \(b\) is divisible by
the cubic coefficient and has degree at most four,
so the quotient used in the \(z\)-shear is affine.

## Support after the triple-line shear

For the first scalar of (1), first subtract \(g(0)\)
and rescale the scalar to make \(g'(0)=1\), adjusting
the original affine transverse coordinates accordingly.
Put \(v=z+xy+g_2(x)\), with \(\deg g_2\le2\), so that
\(r=x+x^2v\). This is a polynomial shear; bounds
on the mates are still measured before the shear.
By the [affine-mate exclusion](KELLER_REPEATED_PARALLEL_SCALAR_AFFINE_MATE.md),
degrees zero and one in the new free \(y\) are impossible.
For ordinary original degree at most \(D\in\{6,7\}\), the
leading coefficient of a mate of free \(y\)-degree \(m\)
has weighted degree at most \(D\) for
\({\rm wt}(x)=1,\ {\rm wt}(v)=2\). Its \(v\)-degree is
at most \(\lfloor(D-m)/2\rfloor\).

Indeed, before accounting for \(g_2(x)\), group the original
polynomial by its homogeneous degree \(d\) in \(y,z\),
and write the group as \(y^d S_d(x,z/y)\), with
\(\deg_x S_d\le D-d\) and \(\deg_s S_d(x,s)\le d\).
If \(d>m\), the absence of powers of \(y\) larger than
\(m\) in the sheared group is equivalent to divisibility
by \((s+x)^{d-m}\). This factor is monic in both \(x\)
and \(s\). The quotient has \(x\)-degree at most \(D-2d+m\)
and \(s\)-degree at most \(m\). Evaluating the quotient at
\(s=-x\) bounds the \(x\)-degree of the coefficient of
\(y^m v^{d-m}\) by \(D-2(d-m)\).
For \(d=m\) the same bound follows by evaluating \(S_d\)
directly. The extra change \(v\mapsto v-g_2(x)\)
does not affect either asserted bound.

## Quadratic mates in the sheared frame

Suppose \(Q=qy^2+B y+C\) in the sheared coordinates.
The
[degree-free quadratic-mate theorem](KELLER_REPEATED_PARALLEL_SCALAR_QUADRATIC_MATE.md)
gives
\[
  q=e(r)x^{2j},\qquad Z=x^j\left(y+\frac{B}{2q}\right),
  \quad K=C-\frac{B^2}{4q},\quad Q=e(r)Z^2+K,\qquad j\in\mathbb Z.
\]
The scalar \(e\) is polynomial in \(r\): along each fibre
there is a nonvertical component with \(x\) generically
nonzero, so a denominator factor in \(e(r)\) would be
a pole of the polynomial \(q\) on that component.
The theorem also gives
\[
  P=\Psi(Q,r)+\int_0^Z F(Q-e(r)u^2,r)\,du,\qquad
  \partial_x J(K)=\kappa x^{-j-2},\qquad J_T=F .
 \tag{2}
\]
The derivative fixes \(r\), and \(\kappa\ne0\).
In particular \(j=-1\) is impossible by the simple pole
of \(x^{-1}dx\).
The odd polynomial has degree \(2s-1\) in \(y\) where
\(s=\deg_T J=\deg_T F+1\).

In the generic scalar fibre ring, \(K\) is Laurent in \(x\).
From \(J(K)-U(r)=\lambda x^{-j-1}\) and \(j\ne-1\)
it follows that \(J(T)-U(r)\) has a unique root \(\tau\)
over \(\overline{\mathbb C(r)}\). In fact each distinct
root would give a unit \(K-\tau\) in the Laurent ring.
Two such units with nonzero constant difference would
both be constant in \(x\), contrary to \(DK\ne0\).
Consequently
\[
   K=\tau(r)+u(r)x^d,\quad J(T)-U(r)=c(r)(T-\tau(r))^s,
             \qquad sd=-j-1,\qquad c u^s=\lambda\ne0 .
 \tag{3}
\]
The coefficients \(\tau,u,c\) lie in
\(\mathbb C[r,e(r)^{-1}]\). For \(\tau,u\) this follows
by comparing the Laurent coefficients of \(K\).
For \(c\), the recursion on the odd part of (2)
uses the fact that the centered polynomial has
coefficients in \(\mathbb C[r,e(r)^{-1},x,x^{-1}]\).
Its leading odd coefficient is invariant,
and determines the leading coefficient of \(F\)
up to a nonzero rational scalar and a power of \(e\).
The same ring has invariant subring
\(\mathbb C[r,e(r)^{-1}]\), and the recursion also places
all coefficients of \(F,\Psi\) in that subring.

The [full exceptional-fibre theorem](KELLER_QUADRATIC_MATE_EXCEPTIONAL_FIBRE.md)
forces \(r\mid q\). Since \(x\) is a unit at a general point
of \(H=0\) for \(r=xH\), this implies \(e(0)=0\).

The support bound gives \(\deg e\le2\) and
\(4\deg e+2j\le D\). Polynomiality gives
\({\rm ord}_0 e+2j\ge0\).
For \(\deg e=2\), these inequalities and \(e(0)=0\)
force \(e\) proportional to \(r^2\) and \(j=-1\),
already excluded. Hence after a constant rescaling
\[
                       q=rx^{2j},\quad j\in\{0,1\}.
 \tag{4}
\]
In particular the case of two quadratic mates is
impossible by their invariant leading ratio, which
in (4) forces a constant ratio. Subtraction would
give an affine mate. Write \(m>2\) for the other mate
and put \(h=-d\); thus \(sh=j+1\).
For an odd \(m\), (2) forces \(m=2s-1\), so only
\(m=3,\ j=1,\ s=2,\ h=1\) is possible.
The leading coefficient is proportional to
\(c(r)r x^3\) by the highest odd term.
Polynomiality along the nonvertical component of
each fibre shows that \(c(r)r\) is a polynomial in \(r\).
For \(D=6\) its weighted bound forces \(c(r)r\)
constant. This makes \(c u^2=\lambda\) impossible
by parity of the \(r=0\) valuation.
For \(D=7\) it can instead be linear in \(r\).
The same parity forces its constant term to vanish,
so \(c\) is a nonzero constant and the leading
coefficient of \(P\) is proportional to \(r x^3\).

## Excluding even higher mate degrees

Suppose \(m=2a\ge4\). The odd part has degree at most
three, so the leading term comes from \(a\)-th degree
of \(\Psi\) in \(Q\). Write \(\psi_a\) for that coefficient.
Polynomiality of the leading coefficient
\(\psi_a r^a x^{jm}\) makes \(\psi_a r^a\) polynomial
in \(r\), testing on the nonvertical part of each fibre.
For \(j=1\) the weighted bound \(m+4\deg(\psi_a r^a)\le7\)
makes \(\psi_a\) have pole order \(a\) at \(r=0\).
For \(j=0\) the \(v\)-degree bound makes
\(\psi_a r^a\) at most linear in \(r\).
A constant leading coefficient is forbidden by the
[formal residue identity](SEXTIC_KELLER_CUBIC_COMPONENT.md#the-exceptional-scalar-and-three-residue-restrictions).
That identity also makes an invariant leading coefficient
vanish at \(r=0\). Thus \(\psi_a\) has pole order \(a-1\)
in this case. Let \(k_a\) denote the resulting positive
pole order.

Use the two components \(x=0\) and \(H=1+xv=0\) of \(r=0\).
Along \(H=0\), \(x\) is a unit and \(r\) has order one.
If \(K\) were regular there, (3) would make \(\tau,u\)
regular at \(r=0\). Then \(B\) vanishes along \(H=0\),
so \(Z\) is a unit-slope affine change of \(y\) over
the local ring and the even and odd parts must
separately be regular. If \(u(0)=0\), the leading
odd coefficient in \(Z\), a nonzero scalar times
\(c(r)r^{s-1}\), has negative order by \(cu^s=\lambda\).
If \(u(0)\ne0\), \(Q|_{H=0}=K|_{H=0}\) is
nonconstant in \(x\), so the coefficient poles of
\(\Psi\) cannot disappear upon substitution.
Both contradict regularity. Thus \(K\) has a simple
pole along \(H=0\), the only possible order there
from \(K=C-B^2/(4q)\). This also makes \(Q|_{H=0}\)
nonconstant as \(B|_{H=0}\ne0\).

Put \(b={\rm ord}_0 u\). The simple pole of \(K\)
gives \(b\ge-1\) and \({\rm ord}_0\tau\ge-1\).
If \(b=-1\), then \(c\) has order \(s\), \(Z\)
has at most a simple pole, and every term of the
integral in (2) is regular along \(H=0\).
For example its terms are scalar multiples of
\(c r^i(Q-\tau)^{s-1-i}Z^{2i+1}\), each of order
at least \(s+i-(s-1-i)-(2i+1)=0\).
Then the poles of \(\Psi\) cannot cancel, a
contradiction. Hence \(b\ge0\) and \(\tau\) has
a simple pole with nonzero leading coefficient.

For a formal variable \(Q\), choose one branch of
\(Z_*^2=(Q-\tau)/r\), and let \(O_*(Q,r)\) denote
the odd integral in (2) evaluated at \(Z_*\).
By substituting \(u'=Z_* t\) in the integral,
\[
                 O_*=\gamma_s\, c(r)r^{-s}
                         [r(Q-\tau)]^{\,s-1/2},
 \qquad
                 \gamma_s=s\int_0^1(1-t^2)^{s-1}dt\ne0 ,
 \tag{5}
\]
with the branch chosen accordingly.
Its pole order is \(N=s(b+1)\), and the coefficient of
\(r^{-N+l}\) has degree at most \(l\) in \(Q\), by
the binomial expansion at the nonzero leading constant
of \(r(Q-\tau)\).

Write \(U=u x^{-h}\) so \(K=\tau+U\).
With a consistent sign for \(Z=\sqrt{(Q-K)/r}\),
expansion of the odd integral at fixed \(Q\) gives
\[
 O-O_*=-\frac{cU^s}{2r Z_*}
                    \left(1+O\!\left(\frac{U}{Q-\tau}\right)\right).
 \tag{6}
\]
Indeed its derivative in \(K\) is
\(-s c(K-\tau)^{s-1}/(2r Z)\), which can be integrated
formally from \(\tau\) to \(\tau+U\) whenever
\(U/(Q-\tau)\) has positive order.
The product \(cU^s=\lambda x^{-(j+1)}\) is independent
of \(r\).

Along \(H=0\), the ratio in (6) has positive order,
and the difference is regular. Since \(Q|_{H=0}\)
is nonconstant, regularity of \(P\) there makes
all negative Laurent coefficients of \(\Psi+O_*\)
in the symbolic variable \(Q\) zero. In particular
(5) and the pole order \(k_a\) of \(\psi_a\) give
\(N\ge a+k_a\). If \(j=0\), then \(s=h=1\),
so \(N\ge2a-1\) gives \(b\ge2\).
If \(j=1\), then \(N\ge2a\) and \(s h=2\);
this gives \(b\ge h\) for both possible \(s\).
In either case \(b\ge h\) and \(N>j+1\).

Now along \(x=0\), \(r=x+x^2v\) has order one, and
\(U\) is regular because \(b\ge h\).
Again the expansion (6) is valid. If the sign of
\(Z_*\) differs from the sign selected on \(H=0\),
then \(P\) has the pole of \(2O_*\) of order \(N\);
the correction (6) has smaller pole order \(j+1\).
So the sign must agree. The expression \(\Psi+O_*\)
has no negative powers of \(r\), as checked on
\(H=0\), and \(Q\) is regular on \(x=0\).
But (6) now has a nonzero pole of order \(j+1\)
since \(cU^s=\lambda x^{-(j+1)}\) and \(rZ_*\)
is a unit. This is incompatible with polynomiality
of \(P\), completing the even case.

## The degree-seven boundary

Alpöge's example has ordinary degrees \((7,6,4)\),
constant Jacobian \(-2\), and generic sheet degree three.
Its quartic component \(2x-3x^2y-x^3z\) is the first form
in (1) after the affine choices
\(y_{\rm new}=-z/2,\ z_{\rm new}=-3y/2,\ r=R/2\).
In the sheared frame \(v=z_{\rm new}+xy_{\rm new}\),
the mates have free degrees \(3,2\) and leading
coefficients \(-16 r x^3/27,-8 r x^2/3\).
The [exact replay](verify_quartic_affine_source_sextic_reduction.py)
checks those coefficients, the Jacobian, and three distinct
rational points with the same image.
[Gao, arXiv:2608.00222v2, Theorem 3.3](https://arxiv.org/pdf/2608.00222v2).
