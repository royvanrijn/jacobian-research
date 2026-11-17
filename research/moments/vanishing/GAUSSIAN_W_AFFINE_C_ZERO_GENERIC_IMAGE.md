# Exact generic original images for W-affine phases with c=0

Let K be any field of characteristic zero, with one Witt pair Z,W and
an independent ordinary Gaussian Y:

    E[Z^i W^j Y^k]=delta_(i,j) i! E[Y^k],
    E[Y^(2q)]=(2q-1)!!, E[Y^(2q+1)]=0.

For arbitrary h,u,k in K[Z] and constant a in K, consider

    P=W h(Z)+aY^2+u(Z)Y+k(Z).                         (G1)

Write Phi[Q]=E[Q exp(tP)], V=Phi[K[Z,W,Y]], and M=K[t]V. This note
classifies the FULL ACTUAL GENERIC IMAGE K(t)M, including a faithful
original marked frame and an exact finite algebraic kernel test.
It does not classify the global K[t]-lattice M or the literal filtered
images of fixed spatial markers.

The proof is algebraic over K. The accompanying exact computation is
bounded internal support; it is not the proof of the unrestricted
statement. Independent replay, formal verification and external review
are not claimed.

## 1. Selected field and exact marked frame

Put F=K(t), and let zeta in tK[[t]] be the unique selected root

    zeta=t h(zeta), J=1-t h'(zeta), d=1-2at,
    R=t k(zeta)+t^2 u(zeta)^2/(2d),
    f=exp(R)/(J sqrt(d)), L=F(zeta).                  (G2)

All selected square roots have constant coefficient one. In particular
J(0)=d(0)=f(0)=1. The field L is embedded in K((t)) by the selected
formal root. Then

    F M=L f subset K((t)),
    dim_F(F M)=[L:F].                                (G3)

If h(0)=0, the selected root is identically zero and L=F: the actual
generic rank is ONE, even when the entire polynomial z-t h(z) has
higher degree or is reducible. The one-element marked frame is f.

If h(0)!=0 and D=deg h>=1, then [L:F]=D and a faithful ORIGINAL frame is

    Phi[1],Phi[Z],...,Phi[Z^(D-1)]
        =f,zeta f,...,zeta^(D-1)f.                   (G4)

For nonzero constant h, L=F and the frame is again f. Thus for h(0)!=0
the rank is max(1,deg h); h=0 belongs to the h(0)=0 case.

Here is the degree argument over the given field K. If h(0)!=0 then
gcd(z,h)=1. The polynomial z-t h(z) is primitive and linear in t over
K[z], hence irreducible in K[z,t]. It is also primitive over K[t]:
any common factor of its z-coefficients divides t because its constant
coefficient is -t h(0), and cannot divide its z-coefficient
1-t h'(0). Gauss's lemma therefore gives irreducibility in F[z].
For D>=1 its z-degree is D. In characteristic zero the resulting
field extension is separable. If h(0)=0, the factor z=0 is selected
because the complementary factor of z-t h(z) is a t-adic unit.

## 2. Full original auxiliary-source identity

Introduce marker variables A,B,C. There is a unique root

    zeta_B=t h(zeta_B)+B, zeta_B(0,B)=B,
    J_B=1-t h'(zeta_B).

The exact original auxiliary source is

    E exp(tP+AZ+BW+CY)
      = exp(A zeta_B+t k(zeta_B)
               +(t u(zeta_B)+C)^2/(2d))
          /(J_B sqrt(d)).                            (G5)

This is a joint formal identity: at every fixed t-order and fixed
marker degrees, all required Wick contractions and algebraic operations
are finite. The root has polynomial B-coefficients at each t-order.
At t=0 the expression is exp(AB+C^2/2), as required by the original
Gaussian covariance.

For completeness, with a t-adic polynomial-coefficient numerator T,
circular Wick contraction gives

    E_ZW[T(Z)exp(W(t h(Z)+B))]
      =Res_z T(z)/(z-t h(z)-B) dz
      =T(zeta_B)/J_B.                                (G6)

For the first equality expand the denominator geometrically; the
coefficient of z^j W^j is contracted by j!. For the second use the
near-identity formal change v=z-t h(z), whose inverse is zeta_v;
its Jacobian is 1/(1-t h'(zeta_v)), and extract the formal Cauchy
residue at v=B. Formal residue invariance can be checked monomial by
monomial: all powers except v^-1 dv are exact, while
dlog(v)=dlog(z)+dlog(1-t h(z)/z) has the same residue. Every t-coefficient
has a finite Laurent principal part. Apply this to the ordinary Gaussian
completion

    E_Y exp(taY^2+(tu(z)+C)Y)
       =d^(-1/2) exp((tu(z)+C)^2/(2d))

and the base source exp(Az) to obtain (G5).

This is the n=s=1, constant-A specialization of the retained
[multibase original contraction](GAUSSIAN_MULTIBASE_QUADRATIC_FIBRE_PRIMITIVITY.md#1-a-self-contained-multibase-wick-contraction),
Section 1. The scalar specialization is also retained in the
[c=0 scalar criterion](GAUSSIAN_QUADRATIC_W_FINITE_IMAGE_CLASSIFICATION.md#3-a-scalar-criterion-for-every-w-affine-phase-with-c0),
equations (C8)-(C11). The present field/frame conclusion concerns the
entire generic marked image, beyond that scalar criterion.

To prove the upper inclusion in (G3), differentiate (G5) finitely many
times in A,B,C and set the sources to zero. Implicit differentiation
gives partial_B zeta_B=1/J_B. Inductively all higher B-derivatives of
zeta_B at zero lie in F(zeta). Differentiating the remaining factors
introduces only rational functions of t,zeta,J^-1,d^-1; the common
undifferentiated factor is f. Hence EVERY original polynomial marker
has its normalized period in L.

For the reverse inclusion, base-only contraction gives
Phi[Z^j]=zeta^j f for every j>=0. The finite field basis in (G4), or
the selected rank-one case, therefore spans Lf. Since f is a nonzero
formal series, any F-linear relation among these actual marked periods
would give the same relation among the field basis after division by f.
This proves faithfulness directly in K((t)), without specialization
of an abstract tensor product.

## 3. Finite algebraic original-kernel test

There is a convenient implementation of the source differentiation.
Temporarily use z as an independent variable and put

    J_z=1-t h'(z), mu=t u(z)/d, v=1/d,
    Rhat(t,z)=t k(z)+t^2u(z)^2/(2d).

Let H_j(mu,v) be the ordinary Gaussian moment polynomials, characterized by

    H_0=1, H_1=mu,
    H_(j+1)=mu H_j+j v H_(j-1).

Define the finite differential operator on F(z)

    Bcal=J_z^-1 partial_z
             +(partial_z Rhat)/J_z+t h''(z)/J_z^2.    (G7)

Then for every i,j,k>=0,

    Phi[Z^i W^j Y^k]/f
       = (Bcal^j (z^i H_k(mu,v))) evaluated at zeta.  (G8)

Indeed differentiation in C gives H_k, and differentiation in B acts
on the normalized factor by J_z^-1 partial_z, while
partial_B log(exp(Rhat)/J_z)
=(partial_z Rhat)/J_z+t h''/J_z^2. This proves (G8) without infinite
source expansion.

For Q=sum q_ijk Z^i W^j Y^k calculate the finite rational expression
sum q_ijk Bcal^j(z^i H_k), then evaluate it in the selected field.
If h(0)=0, evaluate at z=0. If h(0)!=0, reduce numerator and
denominator modulo z-t h(z), inverting the denominator in L.
The denominators d and J_z are nonzero on the selected branch, as
their t-adic values are units. This gives the exact equivalence

    Phi[Q]=0 iff the resulting element of L is zero. (G9)

Equivalently one may apply Q(partial_A,partial_B,partial_C) to (G5).
For a monomial, differentiating is i!j!k! times its source coefficient;
the factorials must be retained. Formula (G9) is a finite algebraic
test for each fixed ORIGINAL spatial marker. It does not say that
the kernel is a spatial ideal.

## 4. Differential connection on the faithful frame

There is a unique extension delta of partial_t to the finite separable
field L. On expressions in t,zeta it is calculated by

    delta=partial_t+(h(zeta)/J) partial_z,
    tau=delta R-delta J/J-d'/(2d), f'=tau f.         (G10)

For h(0)=0 this means evaluation on the selected factor z=0, not on
the whole reducible polynomial. The full generic image is the
restriction of scalars to F of the rank-one L-connection Lf. On the
frame (G4),

    partial_t(zeta^j f)
      =[j zeta^(j-1)h(zeta)/J+tau zeta^j]f.         (G11)

For j=0 the first term is zero. Reduction in the field L gives an
explicit rational connection matrix in the original marked frame.
This does not require the unit alone to be cyclic.

## 5. Unit rank can be strictly smaller than full marked rank

Take h=1+Z^2 and a=u=k=0. Then

    Delta=1-4t^2,
    zeta=(1-sqrt(Delta))/(2t),
    J=sqrt(Delta), f=Delta^(-1/2),
    f'=4t/Delta*f.                                  (G12)

Thus the unit-cyclic differential rank is ONE. The full original
generic image has rank TWO by (G3), with frame f,Phi[Z]. More
explicitly,

    Phi[Z]=(f-1)/(2t),
    Phi[1-2tZ]=1.                                   (G13)

Since Delta is not a square in F, f is not rational, so f and 1
are F-linearly independent. The time-dependent amplitude in (G13)
is an original K[t]-linear combination of two fixed spatial markers.
It does not prove that the literal image V is K[t]-closed. A nontrivial
fixed spatial kernel element in this example is W(1-Z^2), also
detected by (G8)-(G9).

Generic rank one likewise does not mean that M is finite over K[t].
For example P=Y^2 has h=0 and generic rank one, but its original
Y^(2j) periods have unbounded powers of (1-2t)^-1 multiplying
f=(1-2t)^(-1/2). A finite set of K[t]-generators has bounded pole order
after division by f, so it cannot generate them all. Global lattices
and exact literal degree filtrations remain separate questions.

## 6. Exact bounded support

Run

    python verify_gaussian_w_affine_c_zero_generic_image.py \
      --output gaussian_w_affine_c_zero_generic_image_verification.json

from this directory. The [checker](verify_gaussian_w_affine_c_zero_generic_image.py)
compares the normalized original-marker formulas with independent
finite Wick contractions, checks selected-zero and nonzero field
branches, exact kernel/connection identities, and the rank-defect
example. Its [receipt](gaussian_w_affine_c_zero_generic_image_verification.json)
pins the proof, checker and retained source files. These finite
controls do not establish the unrestricted degree or field statements;
those follow from the proof above.
