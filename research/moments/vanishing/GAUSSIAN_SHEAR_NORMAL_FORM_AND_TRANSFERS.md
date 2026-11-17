# Gaussian shear normal forms and exact transfer boundaries

15 November 2026. All-degree algebraic theorem in every dimension with one
circular pair and any number of real spectators. The Gaussian shear
mechanism is already retained in the [affine-W theorem](CUBIC_AFFINE_W_CLASSIFICATION.md)
and its earlier sparse-shear source. The coefficient normal form below is
an explicit deduction, not a novelty claim. Its proof uses polynomial
Gaussian integration by parts and Euclidean division; it does not depend
on the independently unreviewed logarithmic contour argument.

The [checker](verify_gaussian_shear_normal_form.py) verifies generic
coordinate, quadratic, Jacobian, integration-by-parts and operator identities
and two exact Bezout examples. The all-dimension statements have written
proofs below. No independent replay, external review or formal verification
is claimed.

## 1. Polynomial Gaussian shears in every dimension

Let E(ZW)=1, E(T_i T_j)=delta_ij, and all other pairings be zero. For an
arbitrary vector u=(u_1(Z),...,u_s(Z)) of complex polynomials define

    Phi_u:
      Z -> Z,
      T_i -> T_i+Z*u_i(Z),
      W -> W-sum_i u_i(Z)*T_i-Z*sum_i u_i(Z)^2/2.                (1)

The number s>=1 and every degree are arbitrary. These are polynomial
automorphisms with

    Phi_u^-1=Phi_(-u),       Phi_u Phi_v=Phi_(u+v),
    det Jac(Phi_u)=1,
    2ZW+sum_i T_i^2 invariant.                                (2)

More strongly, they preserve the complete polynomial Gaussian functional:

    E(F composed with Phi_u)=E(F) for every polynomial F.       (3)

To prove (3), put

    D_u=sum_i u_i(Z)*(Z*d/dT_i-T_i*d/dW).

Give T_i auxiliary degree one, W degree two and Z degree zero. D_u strictly
lowers this degree, so it is locally nilpotent and exp(D_u) is a finite
sum on every polynomial. Its coordinate values give (1). Gaussian
integration by parts gives

    E(u_i Z F_(T_i))=E(u_i F_(W T_i))=E(u_i T_i F_W),

because u_i depends only on Z. Hence E(D_u F)=0 for every F, and all
positive terms in exp(D_u) have zero expectation. This proves (3)
algebraically; no real change of integration contour is assumed. Preserving
a quadratic form and volume alone is not being used as an unproved general
criterion for Gaussian preservation.

## 2. A finite-dimensional coefficient normal form in arbitrary degree

Consider

    P=W*(a(Z)+sum_i b_i(Z)*T_i)+G(Z,T_1,...,T_s),               (4)

with arbitrary polynomial a,b_i,G and at least one b_i nonzero. Let

    g=gcd(b_1,...,b_s), chosen monic.

Under (1), expressed in the new coordinates, the W-coefficient has normal
form

    R(Z)+sum_i b_i(Z)*T_i,
    R=a mod (Z*g),             deg R<deg g+1.                 (5)

The residue class of a in C[Z]/(Z*g) is the complete orbit invariant of
the W-coefficient under this shear subgroup at fixed b. It is not a
classification of the whole polynomial P, because G changes as well.

Choose Bezout polynomials v_i with sum_i b_i v_i=g. Divide a=Z*g*q+R and
set u_i=q*v_i. In new coordinates use Phi_u^-1, so T_old=T-Zu and
W_old=W+u dot T-Z*(u dot u)/2. Therefore

    P composed with Phi_u^-1=W*(R+b dot T)+G_tilde,
    G_tilde=(u dot T-Z*(u dot u)/2)*(R+b dot T)+G(Z,T-Zu).       (6)

Everything is polynomial. Conversely a shear changes a by an element of
Z*(b_1,...,b_s)=(Z*g), proving necessity and completeness of the coefficient
invariant. A different Bezout choice can act on G through a stabilizer
b dot u=0. For example u_i=b_j*v,u_j=-b_i*v is such a stabilizer for i!=j.
No uniqueness of G_tilde or full-orbit classification is asserted.

For one spectator g is the monic scalar multiple of b, so arbitrary a
reduces modulo Z*b to degree at most deg b. In particular every nonzero
constant-b family reduces to a constant a=a(0), regardless of the original
degree of a. Several coprime b_i give the same constant normal form.

If g has order m_0 at zero, the residue class records the jet of a through
order m_0 there. At each nonzero root alpha of multiplicity m_alpha it
records the jet through order m_alpha-1. These are the Chinese-remainder
data of C[Z]/(Z*g). High-degree a is removable; these root jets and the
transformed G remain. If all b_i=0, the subgroup leaves a unchanged; do
not perform a gcd-zero division in that separate case.

## 3. GMC safety and counterexamples transfer exactly

For P_u=P composed with Phi_u, equation (3) implies

    E(P_u^m)=E(P^m),
    E(Q*P_u^m)=E((Q composed with Phi_u^-1)*P^m).                (7)

Multipliers are transported bijectively. All-order pure vanishing,
eventual pure vanishing, GMC safety and GMC failure are invariant. A
witness Q for P becomes Q composed with Phi_u for P_u and retains its
mixed sequence at every index. Dimension is preserved. Degree, homogeneity
and term counts need not be preserved.

Let delta be the maximum degree of a coordinate of Phi_u or Phi_u^-1.
For nonzero u with j=max deg u_i, delta<=max(1,2j+1). Use actual coordinate
degrees when sum_i u_i^2 has leading cancellation. If P has cutoff
m>C*deg Q for every Q, then P_u has the valid bound

    m>C*delta*deg Q.                                           (8)

Apply the old bound to Q composed with Phi_u^-1 to prove this. Optimality
of the transported bound is not automatic.

The loss can be sharp. In one spectator take u=-Z^j,j>=1. The coordinate
W, with cutoff coefficient one, becomes

    S_j=W+Z^j*T-Z^(2j+1)/2.

The retained sharp shear theorem G3SHEARCUTOFF20261114 proves

    E(W^d*S_j^((2j+1)d))=((2j+1)d)!*(-1/2)^d !=0.

Here delta=2j+1 is exactly the sharp coefficient. Nonlinear Gaussian
conjugacy preserves qualitative safety, but must not silently preserve a
sharp linear degree-Q cutoff.

## 4. The ordinary Laplacian is preserved exactly in the linear case

Write L=2*d_Z*d_W+sum_i d_(T_i)^2 and D_i=Z*d_(T_i)-T_i*d_W.
Each D_i commutes with L. Since [L,u_i(Z)]=2*u_i'(Z)*d_W,

    [L,D_u]=2*sum_i u_i'(Z)*(Z*d_(T_i W)-T_i*d_(W W)).          (9)

This vanishes as an operator iff every u_i is constant. Those are the
linear null rotations in this subgroup. Moreover the same iff holds
for the finite substitution commuting with L. Indeed, from (1),

    L((W composed with Phi_u)^2)
       =-4*sum_i u_i'(Z)*(T_i+Z*u_i(Z)).                       (10)

If the finite substitution commuted with L, this would be zero because
L(W^2)=0; its independent T_i coefficients force all u_i'=0. Constant u
is linear orthogonal and does commute with L.

For a concrete all-power obstruction take s=1,u=Z,F=W^2. Every L^m(F^m)
is zero because F depends only on W. But

    F composed with Phi_Z=(W-ZT-Z^3/2)^2,
    L(F composed with Phi_Z)=-4*(T+Z^2) !=0.                   (11)

Thus even the pure differential-vanishing hypothesis is lost at m=1.
Conjugation gives a valid order-two operator with variable coefficients;
it does not give the original constant-coefficient GVC/Laplacian problem.
This is the exact transfer boundary, consistent with Gaussian invariance.

## 5. An unbounded-degree Gaussian orbit class with a finite pure test

More generally let Psi be a finite composition of the shears above and
complex linear orthogonal transformations of Gaussian coordinates. Each
factor preserves the Wick functional, so (7) applies to Psi, with delta
its actual polynomial degree and inverse degree bound. For every base
class C, its orbit class {P composed with Psi:P in C} has exactly the
same possible pure moment sequences as C. In particular any sufficient
pure-prefix length K transfers unchanged, and the least universal prefix
is the same when the identity transformation is admitted. Before applying
a prefix proved only for a bounded-degree normal form, verify that the
underlying base polynomial belongs to that class.

Apply this to the internally closed ternary degree<=3 class, using
G3CUBICALL20261115 and the retained G3PURE9420261115. Every member of its
Gaussian orbit class is GMC-safe. Its first94 zero pure moments imply all
positive pure moments zero, although its degree is unbounded. The least
universal prefix for this orbit class equals the still-uncomputed least
cubic prefix, in[20,94]. A valid mixed cutoff is m>3*delta*deg Q. Neither
94 nor the transported coefficient is asserted optimal. These corollaries
inherit the independent analytic-review obligations of the base cubic
proof; the algebraic transfer itself needs no contour theorem.

The class has arbitrarily high degree because it contains the S_j above.
It is a proper class. The [all-degree barrier](ALL_DEGREE_STRUCTURE_AND_BARRIERS.md)
family P_d=Z+W^d has first nonzero index d+1, so for d>=94 it cannot be
Gaussian-conjugate to a cubic. Even within the pure-null locus, the
retained quartic failure

    B=W*(1+Z)-(1+3*Z/2+Z^2/2)*T^2,
    E(Z*B^m)=m! for every m>=1,

cannot be Gaussian-conjugate to a cubic, because mixed failure is invariant
and the cubic class is safe. No claim that every GMC-safe polynomial belongs
to this orbit class is made. Thus a uniform pure test on this structured
unbounded-degree family is consistent with the impossibility of a uniform
test for all degrees and all polynomials.

## 6. Consequences for the all-degree frontier

The [squarefree linear-T exclusion](SQUAREFREE_LINEAR_T_ALL_DEGREE.md)
already covers arbitrary a,G when b is nonconstant squarefree. For the
remaining constant or repeated-root b, normalize a by (5) first. Constant
b requires only h=a(0)+bT; repeated-root b requires only a mod Z*b, with
the exact G transformation (6).

This does not classify the remaining G. Even P=WT and P=WT+Z^2 have the
same W-coefficient but different pure behavior: the former is pure-null,
whereas the latter has E(P^3)=6. The
[complete h(0,T)=0 branch](NEUTRAL_FACE_ALL_DEGREE_CLASSIFICATION.md)
remains closed in three variables; this theorem does not extend its
one-spectator analytic classification to more spectators. The algebraic
shear, coefficient normal form and Gaussian transfer are uniform in dimension.

Seek a complete classification of G in this normal form or an invariant
closing another whole architecture. Do not start separate calculations for
ever larger degrees of a or transfer Laplacian claims through a nonlinear
shear. This algebraic theorem does not upgrade earlier analytic assurance
or change any least cubic prefix.
[Source audit](gaussian_shear_normal_form_source_audit.json).
