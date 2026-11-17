# Cubic scalar operators with algebraic coefficients and every finite mate

**QSO3A.** Let F be a one-variable algebraic function field over a
characteristic-zero constant field C, with a nonzero C-derivation.
If P has differential order three and [P,Q]=kappa is a nonzero
constant, then
\[
                  Q=A(P)+Q_1,\qquad A\in C[z],\qquad
                              \operatorname{ord}Q_1=1.      \tag{1}
\]
In particular no such pair has Q order at least two coprime to
three. Arbitrarily many lower constant flows and arbitrary finite
algebraic coefficient covers are included. As usual C denotes the
actual constant field; extension of constants and descent give the
same assertion over non-algebraically closed characteristic-zero fields.

After a finite algebraic coefficient extension and a change and shift
of the differential generator normalizing P to E^3+uE+w, every
mate is
\[
 Q=A(P)+aE+b,\quad a,b\in C,\quad a\ne0,\qquad
                  u'=0,\quad w'=-\kappa/a.                  \tag{2}
\]
Thus these possibilities really occur, as witnessed by E^3+x and
-E. The general algebraic cubic obstruction does not assert that a
cubic operator has no finite mate of any order.

The proof replaces a meromorphic-reconstruction obligation by a weaker
fact sufficient for algebraic solutions: scalar Stokes data are
injective. Algebraic continuation around a branch point already exists;
injectivity forces the continued branch to return. A universal Jacobi
polynomial identity then supplies the infinity contradiction for all
partner orders simultaneously.

## A finite normal form and the time primitive

First extend constants algebraically and choose a component. Adjoin a
cube root of P's leading coefficient, rescale its differential generator,
and shift that generator to remove its quadratic term. These finite
Ore changes give
\[
                           P=E^3+uE+w.                      \tag{3}
\]
Its new derivation is denoted by a prime. Constants remain constants.
Let W=P^(1/3) be the unique formal pseudodifferential root with
W=E+O(E^-1), and set B_j=(W^j)_+ for j>=0. Its coefficients
are differential polynomials over Q in u,w: at every recursive
step the new root coefficient occurs with coefficient three.
Each specified coefficient uses only finitely many derivatives and
formal coefficients. No infinite operator is asserted to be actual.

Every [P,B_j] has order at most one, since W^j commutes with P and
its negative part has order at most minus one. If a finite operator
with commutator of order at most one has current leading term fE^j,
its commutator with P has coefficient 3f' at E^(j+2). Hence f is
constant. Subtracting fB_j lowers the order. Finite repetition gives
\[
                           Q=\sum_j c_jB_j,\qquad c_j\in C. \tag{4}
\]
Remove the terms j divisible by three, which are polynomials in P.
The remaining highest order m is not divisible by three. A nonzero
commutator requires a nonzero remainder. We will exclude m>=2.

Put R=sum_(j) c_j W^j and write its negative part as
R_-=r_1 E^-1+r_2 E^-2+... . Direct multiplication gives
\[
 [P,Q]=[R_-,P]
          =-3r_1'E-3(r_1''+r_2').                          \tag{5}
\]
The right side has no lower-order remainder because its left side is
a differential operator. Thus a scalar commutator implies
\[
                         r_1'=0,\qquad r_2'=-\kappa/3,
           x=-3r_2/\kappa\in F,\qquad x'=1.                \tag{6}
\]
This constructs the time variable before any analytic theorem is used.
It is nonconstant, so F is finite over C(x), and u,w are algebraic
functions of x. The derivation is the extension of d/dx.

## Primary-source input and its precise scope

The primary input is [Beals--Sattinger, arXiv:solv-int/9801010v1,
Section 9](https://arxiv.org/abs/solv-int/9801010). The remark on
printed page 43 includes finite combinations with arbitrary constant
coefficients. Theorem 9.1 and Lemma 9.3 supply a spectral system with
an ordinary finite spectral point and a fixed formal phase. Theorem
9.4 is scalar Stokes injectivity. We use that forward statement,
not the generic inverse-data solvability hypothesis of Theorem 9.5.
The exact inspected pages, cached PDF hash, and imported scope are in
[the source audit](all_cubic_source_audit.json).

Here are the needed admission and injectivity details, including the
ratio orientation. They also explain why an artificial resonant pole
from a singular gauge causes no gap in this application.

After a harmless nonzero scaling/sign convention for Q, the compatible
scalar equations are P f=lambda*f and f_lambda=Q f/kappa.
Writing Psi=(f,f_x,f_xx)^t gives a rank-three companion system. Reducing
each derivative in Q modulo P-lambda expresses its spectral matrix as
a polynomial in lambda. Pull back lambda=z^3. The original companion
spectral equation is ordinary at every finite z, including z=0.
The singular shearing used to diagonalize infinity need not be used
at zero.

For m coprime to three, the normalized sectorial fundamental matrices
at infinity have the form
\[
 \Psi_\nu=\Lambda_z(I+O(z^{-1}))e^{\Phi(x,z)},\qquad
 \Lambda_z=d(z)\Lambda,\quad d(z)=\operatorname{diag}(1,z,z^2),
                                                                  \tag{7}
\]
where Lambda is the fixed Vandermonde matrix of the cube roots of
unity and J=diag(1,omega,omega^2). The leading phase coefficient is
a nonzero constant multiple of J^m. Its eigenvalues are distinct
for all finite x; lower flows cannot change that leading coefficient.

The finite polynomial/logarithmic part of the phase is
\[
                         \Phi(x,z)=\Phi_0(z)+xzJ.            \tag{8}
\]
All coefficients in Phi_0 are constants on the coefficient curve.
For clarity, formal diagonalization at infinity recursively divides
only by differences of the fixed distinct leading eigenvalues. Its
finite phase coefficients therefore lie in the algebraic coefficient
field, with its lower-flow constants. The x equation, after conjugation
by Lambda_z, has the form zJ+O(z^-1). Comparing diagonal polynomial
and logarithmic parts in the simultaneous formal solution forces
Phi_x=zJ; all other phase coefficients have derivative zero. They
belong to C and are single-valued. This is the content of the
formal construction used in Section 9, and retains every lower flow.

The right Stokes jumps Psi_(nu+1)=Psi_nu*S_nu are constant in x.
The sectors may be fixed throughout finite time, since the leading
irregular eigenvalues never coalesce. No nonresonance condition at an
apparent regular singularity and no genericity of the Stokes matrices
is required here.

Suppose two local scalar solutions have the same phase and Stokes
matrices. Use the **left ratio**
\[
                        G(z)=\widetilde\Psi_\nu\Psi_\nu^{-1}.
                                                                  \tag{9}
\]
Equal right jumps cancel, so G is independent of the sector. Both
companion spectral equations are ordinary at finite z, making G
entire there. Equation (7) gives
\[
           G=d(z)\Lambda(I+O(z^{-1}))\Lambda^{-1}d(z)^{-1}.
                                                                  \tag{10}
\]
Hence G_(ij)=delta_(ij)+O(z^(i-j-1)), with indices 0,1,2. Sectorial
estimates cover infinity. Entire functions with these polynomial
bounds make G a unipotent lower triangular polynomial matrix;
its first row is exactly (1,0,0). Therefore the scalar first rows
of the two fundamental matrices agree. This holds on their common
x domain, so their x derivatives give the second and third companion
rows as well. The scalar coefficients u,w coincide.
This reconstructs the injectivity argument without inheriting the
inconsistent ratio orientation in the text extraction of Section 8.

## Algebraic continuation gives rational coefficients

Take an algebraic coefficient branch and continue it around a finite
x value, along a loop avoiding its finitely many singular values.
The formal constants in (8) and all Stokes matrices remain unchanged.
On return the continued branch has a common local domain with the
original branch. The injectivity argument forces them to agree.
Thus the functions u,w have trivial monodromy around every finite
branch value, including pole branches. No global inverse-monodromy
reconstruction, analytic Fredholm continuation, or assertion of
meromorphic scalar coordinates for arbitrary resonant systems is used.

For completeness, first work over complex constants. For arbitrary
characteristic-zero constants, choose a finitely generated constant
subfield k_0 supporting the finitely many algebraic minimal polynomials
of P,Q's coefficients over k_0(x), together with the nonzero string
constant and the finite normalization data. Adjoin those coefficients
to k_0(x), obtaining a finite extension F_0. It is stable under d/dx:
differentiating a separable minimal polynomial expresses the derivative
of its root in the same field. Embed k_0 into the complex numbers,
**keeping x transcendental**. The tensor product
F_0 tensor_(k_0(x)) C(x) is a nonzero finite separable algebra. Choose
a field component. The restriction F_0 to that component is injective,
being a unital map of fields, and the derivations agree by uniqueness
of algebraic extension. Thus the complex argument applies without
sending the constructed time x to a scalar. This uses finite algebraic
data, not an embedding of an unrestricted constant field.

The function field C(x,u,w) is unramified over every finite x.
If its degree is d, its genus is g and it has b places at infinity,
Riemann--Hurwitz gives 2g-2=-d-b. Consequently d=b=1 and g=0.
Both u and w are rational functions of x.

## One dispersionless polynomial controls every infinity balance

For m>=1 not divisible by three define
\[
 H_m(u,w)=\operatorname{res}_{\xi}
                         (\xi^3+u\xi+w)^{(m+3)/3}.         \tag{11}
\]
It is a polynomial of weight m+4 for wt(u)=2, wt(w)=3. Its two
partial derivatives are the leading, derivative-free string invariants:
\[
 \partial_w H_m=\frac{m+3}{3}[\xi^{-1}](\xi^3+u\xi+w)^{m/3},
 \partial_u H_m=\frac{m+3}{3}[\xi^{-2}](\xi^3+u\xi+w)^{m/3}.
                                                                  \tag{12}
\]
We need the universal fact that these derivatives have no common zero
away from (0,0), not a finite check of selected m.

Choose a,b in {0,1} and an integer l>=0 so that
\[
                 m+4=2a+3b+6l.
\]
These choices exist for all m under consideration. Put
nu=(2+2a)/3 and epsilon=(1-2b)/2. Direct multinomial expansion gives
\[
 H_m=C_m u^a w^{b+2l}
 {}_2F_1\left(-l,-l+\epsilon;\nu;
                           -\frac{4u^3}{27w^2}\right),      \tag{13}
\]
with C_m nonzero. This hypergeometric expression is a terminating
polynomial identity after clearing the displayed w powers.
For example the ratio of consecutive coefficients comes from the
support (a+3k,b+2l-2k); cancellation of the triplication and duplication
factors gives exactly the two numerator and one denominator rising
factorials in (13).

Let z=-4u^3/(27w^2). The finite Pfaff identity transforms its polynomial
factor into
\[
 (1-z)^l\frac{l!}{(\nu)_l}
 P_l^{\alpha,\beta}\left(1-\frac{2z}{z-1}\right),\qquad
       \alpha=\nu-1\in\{-1/3,1/3\},\quad
       \beta=-\epsilon\in\{-1/2,1/2\}.                    \tag{14}
\]
The Jacobi polynomial has l simple roots in (-1,1), with neither
endpoint a root. One elementary justification is its Rodrigues
formula: integration by parts makes it orthogonal to all lower-degree
polynomials for the positive integrable weight
(1-x)^alpha*(1+x)^beta. Fewer than l simple sign changes would give
a lower-degree test polynomial with a nonzero signed integral, a
contradiction. The endpoint values are nonzero binomial coefficients.
This proves the needed root statement for every l.

Under the change of argument in (14), these roots yield distinct,
finite, nonzero values of u^3/w^2. Thus H_m is a nonzero scalar times
u^a*w^b times distinct factors u^3-rho_i*w^2, with every rho_i nonzero.
On u*w nonzero, weighted Euler's identity says a common gradient zero
must lie on H_m=0; a simple such factor has a nonzero gradient.
On either axis the corresponding prefactor has multiplicity at most
one, and the remaining endpoint coefficient is nonzero. This proves
\[
       \partial_uH_m=\partial_wH_m=0
                              \quad\Longrightarrow\quad u=w=0.
                                                                  \tag{15}
\]

## The uniform rational infinity contradiction

Give a derivative weight one. The coefficient r_i from W^j has
weight j+i. At rational infinity, let
R=max(deg(u)/2,deg(w)/3)>0. For a differential monomial of weight k
and total derivative order d, its ordinary x degree is at most
R*k-(R+1)*d. Therefore derivative terms are strictly lower than
their derivative-free weight, and every lower flow j<m has smaller
degree than the top flow at the same invariant.

Let (u_*,w_*) be the nonzero weighted leading pair. Since r_1 is
constant, its leading degree R(m+1) must vanish, so
partial_w H_m(u_*,w_*)=0. By (15), partial_u H_m at that pair is
nonzero. Hence r_2 has degree R(m+2). But (6) makes its degree one.
If m>=2, rationality gives R>=1/3 and
\[
                          R(m+2)\ge4/3>1,
\]
a contradiction. If u,w are both bounded at infinity, all their
derivatives and every term of r_2 are bounded, also contradicting (6).
This excludes every m>=2 coprime to three, with every lower flow.

Only m=1 remains after the finite polynomial subtraction (4).
Directly comparing [E^3+uE+w,aE+b] gives a,b constant, u'=0 and
-a*w'=kappa. This proves (1)--(2), and descending successive
order-three leading-coefficient ratios proves (1) over the original
constant field without requiring descent of the normalizing cube root.

## The prescribed cubic Long repair is closed in every partner order

Suppose actual regular finite scalar S,T in the exact-R chart have
[S,T]=h, generic ord(S)=3, and classical orders two and three.
Equation (1) gives T=A(S)+T_1 of order one. Absorb A's constant term
into T_1. Every coefficient of A is h-regular: otherwise its largest
pole reduces to a nonconstant polynomial in the order-two classical
S, which cannot be cancelled by the order-one T_1. Therefore T_1
is regular too. Reduction has order at most one or a positive even
order, never the prescribed order three. Such a generic cubic
repair is impossible for every finite partner order.

This does not exclude generic scalar order at least four, arbitrary
nonpositional first images or genuinely nonscalar carriers. It supplies
neither a finite bounded Long lift nor a solution of DC1/DC2.

[The exact checker](check_all_cubic_algebraic_string.py) and
[its retained receipt](all_cubic_controls.json) verify declared finite
normal-form identities, the time coefficient formula, polynomial
spectral reduction, the left-ratio bounds, terminating Jacobi identities
and sample gradient resultants. The all-order normalization, Stokes
injectivity application, continuation, Jacobi factorization and infinity
argument use the written proof and pinned primary source, rather than
an extrapolation from those controls. No formal verification, independent
whole-proof replay, external review or novelty priority is claimed.
