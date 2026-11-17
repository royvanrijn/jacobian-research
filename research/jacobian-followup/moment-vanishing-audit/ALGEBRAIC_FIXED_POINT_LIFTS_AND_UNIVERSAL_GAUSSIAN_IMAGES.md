# Algebraic fixed-point lifts and universal regular Gaussian images

16 November2026. Internal all-degree/all-ring proofs. No independent,
formal or external assurance, and no priority claim. The previous
[origin theorem](TIME_ZERO_SATURATION_AND_PRINCIPAL_PART_RECONSTRUCTION.md)
and all sealed Oaku/formal/ordinary packets are unchanged. The
[checker](verify_algebraic_gaussian_lifts.py),
[certificate](algebraic_gaussian_lift_certificate.json) and
[source audit](algebraic_gaussian_lift_source_audit.json) retain finite
transcriptions; universal assertions use the proofs below.

## 1. A polynomial determinant compensator

Let R be ANY commutative Q-algebra and H(U,Z) any polynomial. Use seven
Gaussian variables U,A,Z,W,Y,V,T with covariance pairings

    E[U A]=E[Z W]=E[Y V]=E[T^2]=1,

all other pairs zero. This is fixed nonsingular rational split precision.
It is equivalent to the standard Gaussian after extending coefficients
by i: a circular pair is U=x+i*y,A=(x-i*y)/2. No nonzero real pure phase
is being asserted. If H is independent of U, omit U,A and work in FIVE
Gaussian variables.

Put C(Y)=4+6Y+4Y^2+Y^3 and define the POLYNOMIAL phase

    P_H=A+W H(U,Z)+V H_z(U,Z)(1+Y)+G(U,Z,Y)T^2,
    G=-(1/2) H_z(U,Z)(1+Y) C(Y).                          (1)

When U,A are omitted, omit the A term and write H=H(Z). If deg H=m>=1,
deg P_H<=m+5. No degree-minimality assertion is made.

Let g in tR[[t]] be the unique fixed point and define units D and 1+y:

    g=t H(t,g), D=1-t H_z(t,g), y=D^-1-1.                (2)

Existence/uniqueness is coefficient recursion; no field, domain,
Noetherian or convergence hypothesis enters. Circular Gaussian
contraction is formal Lagrange residue. At the chosen fixed points,

    U=t, Z=g, Y=t H_z(t,g)(1+Y)=y.

The normal-pair Jacobian determinant is D^2. Indeed the U equation has
derivative1, the Z equation derivative D, and the Y equation derivative
D; the constraint matrix is triangular in this order. The identity

    tG(t,g,y)=(1-(1+y)^4)/2

shows that the remaining T precision is (1+y)^4=D^-4. Its normalized
Gaussian expectation contributes D^2 and EXACTLY cancels the normal
Jacobian. Thus

    Phi_H[Q(U,Z,Y) T^(2k)]
      =(2k-1)!! Q(t,g,y) D^(4k), k>=0,                  (3)
    Phi_H[Q(U,Z,Y) T^(2k+1)]=0.

In particular Phi_H[1]=1 in ALL orders, over EVERY R. These are formal
identities, not contour substitutions. A proof of the circular residue
rule is to expand exp(tW h): E[Z^a W^b]=delta_ab a! gives the Lagrange
coefficient formula. Its version with a W^j marker is

    j! Res Q(z)/(z-t h(z))^(j+1) dz.

Apply it successively to the triangular normal pairs. Each coefficient
uses only finitely many polynomial Wick contractions. This is the same
classical Gaussian-Lagrange mechanism retained in
[the earlier pullback proof](RATIONAL_GAUSSIAN_PULLBACKS_AND_IRREGULAR_RANK.md),
now with an explicit polynomial inverse-Jacobian variable and compensator.

## 2. The FULL actual image, including momentum amplitudes

The full actual polynomial-Weyl image is EXACTLY the commutative ring

    M_H=R[t,g,D^-1] subset R[[t]].                        (4)

Here the right side is an actual subring of formal series, not an
assertion that the abstract implicit algebra has no other components.
The inclusion from right to left follows from(3): every polynomial in
t,g,y is the period of the ORIGINAL polynomial marker Q(U,Z,Y).
Time-polynomial coefficients are permitted in the exact source; U can
also realize t as an original spatial marker in the seven-variable lift.

For the other inclusion retain EVERY amplitude A^a W^b V^c T^j.
Introduce sources alpha,beta,gamma for A,W,V. Their fixed-point equations
and determinant are

    u=t+alpha, z=beta+t H(u,z),
    D_*=1-t H_z(u,z), y_*=(gamma+t H_z(u,z))/D_*,
    normal determinant D_*^2.

Since 1+y_*=(1+gamma)/D_*, the remaining T precision is

    B_*=(1+y_*)^4-gamma C(y_*).

For a coordinate polynomial Q, before differentiating the sources the
T^(2k) expectation is

    (2k-1)!! Q(u,z,y_*) D_*^-2 B_*^(-k-1/2).             (5)

Differentiate (5) a,b,c times and set sources zero. This gives the
ORIGINAL A^a W^b V^c marker, with no factorial ambiguity. Implicit
differentiation of z uses only D_*^-1. All derivatives of u,z,y_*,D_*
therefore lie in R[t,g,D^-1] at zero sources. There B_*=D^-4; derivatives
of its half power are rational constants times integer powers of D.
Odd T powers vanish. Consequently EVERY original polynomial amplitude
lands in the ring(4). This is a finite symbolic constructor for each
specified amplitude, not a full actual-image-kernel computation.
For example g'=(H+tH_U)(t,g)/D, so(4) is closed under differentiation.

The Jacobian-localized algebra

    R[t,z,(1-tH_z)^-1]/(z-tH(t,z))

maps onto(4) via the selected formal branch. Its Jacobian is invertible.
The relevant classical interpretation is a smooth presentation of
relative dimension zero: [Stacks00T6](https://stacks.math.columbia.edu/tag/00T6),
[00T7](https://stacks.math.columbia.edu/tag/00T7) and
[00U0](https://stacks.math.columbia.edu/tag/00U0). The equation need NOT
be monic in z; the monic definition of standard etale is not invoked
without a change of presentation. No all-ring injectivity or branch
component/base-change assertion is imported from the geometric terminology.

## 3. Two fixed linear markers and a complete all-ring safe normal form

Equation(4) gives, over EVERY R,

    P_H safe iff g in R[t] AND D^-1 in R[t].              (6)

Thus polynomiality of the TWO FIXED LINEAR markers Z and Y decides
all-original-marker safety throughout this arbitrary-degree family.
This is not a claim that two linear markers suffice for every Gaussian
phase. If either selected series is nonpolynomial it itself is a fixed
linear bad multiplier. Purity is already proved in(3).

There is an exact finite-polynomial normal form for the safe members.
They are precisely the H for which there exist p in tR[t] and B(t),
J(t,z) polynomial such that

    H(t,z)=p(t)/t+(z-p(t))B(t)+(z-p(t))^2 J(t,z),         (7)

and the ideal generated by the finitely many coefficients of B is
NILPOTENT. Indeed p is then the unique fixed point g, and
D=1-tB. A polynomial with constant coefficient1 is a unit of R[t]
exactly when all its nonconstant coefficients are nilpotent. Finitely
many nilpotent coefficients generate a nilpotent ideal, even for
non-Noetherian R. Conversely(6) gives p=g, D a polynomial unit and
the required nilpotent B=H_z(t,p); division by the monic z-p twice
gives(7), over rings with arbitrary zero divisors.

Over a characteristic-zero FIELD this reduces to B=0, hence the exact
safe normal form

    H(t,z)=p(t)/t+(z-p(t))^2 J(t,z), p in tK[t].          (8)

The constant-H boundary is safe. H=z^2 has p=0 and is safe in every
degree variant (8). H=1+z is already an unsafe rational example.
Over R=Q[e]/e^2, H=1+e*z gives p=t+e*t^2 and D^-1=1+e*t;
it is safe with nonzero nilpotent B. Thus replacing the ordinary
nilpotent condition in(7) by the field equation B=0 would be false.

## 4. EVERY algebraic formal series occurs in seven Gaussian variables

Let K be characteristic zero and f in K[[t]] algebraic over K(t),
with its actual branch specified. There is an explicit polynomial
H(U,Z) and an ORIGINAL polynomial marker Q(U,Z) in the seven-variable
pure phase(1) for which

    Phi_H[Q]=f,
    K(t) M_H=K(t,f), EXACTLY.                            (9)

Here is a direct univariate encoding proof, including a singular
plane-model branch. Choose F(t,x) in K[t,x] with F(t,f)=0 and
F_x(t,f)!=0. Write v=ord_t F_x(t,f), choose integer k>v, and retain
a(t)=f mod t^k, a polynomial of degree<k. Taylor expansion proves

    G(t,u)=t^(-k-v)F(t,a(t)+t^k*u) in K[t,u].

Its derivative G_u(0,u) is a NONZERO CONSTANT lambda: terms quadratic
in u have extra t-order because k>v. The chosen constant u0 is the
root of G(0,u). Put

    R(t,w)=w-lambda^-1 G(t,u0+w).

Then R(0,0)=0 and R_w(0,0)=0. Let c be the coefficient of t in R(t,0).
The polynomial

    H(t,z)=(R(t,t*(c+z))/t-c)/t                         (10)

is well-defined: the first division uses R(0,0)=0 and absence of its
linear w term; after substitution the constant in t is exactly c,
so the second division is also polynomial. The unique fixed point
g=tH(t,g) recovers

    f=a(t)+u0*t^k+c*t^(k+1)+t^(k+1)*g.                  (11)

Its marker is (11) with t replaced by U and g by Z. Conversely g is
rational in t,f, so the fields in(9) agree. All original amplitudes
are in that field by(4), while powers of the ORIGINAL marker Z span
it over K(t). This proves exact full actual rank [K(t,f):K(t)], not
only a marked scalar order or an annihilating-source upper bound.
Finite branch jets are essential input; no choice of algebraic conjugate
is hidden. Over Q the construction stays over Q with split precision.

Artin--Mazur implicit polynomial encoding is classical; the introduction
of [Alonso--Castro-Jimenez--Hauser1403.4104v1](https://arxiv.org/html/1403.4104v1)
states that encoding and the need to distinguish branches. We use no
division/standard-basis theorem or complexity claim from that paper.
Equations(10)-(11) supply the special one-variable encoding directly.
The Gaussian compensator is an application/construction here; historical
priority has NOT been investigated.

All generic connections in(9) have FINITE geometric monodromy and are
REGULAR SINGULAR on the entire projective time line. To see this without
testing a nonminimal equation, split the algebraic field at any geometric
place into tame Puiseux fields Kbar((s)), t-a=s^e after a unit change.
In basis 1,s,...,s^(e-1), (t-a)partial_t is diagonal with entries j/e.
There are no positive slopes; the same argument at infinity applies.
Analytic continuation permutes the finitely many field embeddings.
Punctual lattice changes preserve regular holonomicity. The field
connection, not a larger exponential source, is being classified.
This CLOSES regular NONRATIONAL ALGEBRAIC realizability in full generality
for algebraic formal field branches; regular transcendental realizations
and classification of all Gaussian connections remain separate.

## 5. Five variables already realize every rank and symmetric monodromy

Take H(z)=1+z^d, d>=2, and use the FIVE-variable version of(1). Put

    h=1+z^d, q=h-z h'=1-(d-1)z^d, t=z/h.

The exact full global image(4) becomes the explicit differential ring

    M=K[z,h^-1,q^-1] subset K[[t]],
    partial_t=(h^2/q)partial_z.                          (12)

Indeed gcd(z,h)=1 gives h^-1 in K[t,z]; D=q/h, and gcd(h,q)=1
gives q^-1 using h/q=D^-1 and Bezout. The reverse inclusion is immediate.
There is no actual-image-kernel calculation. The polynomial z-t h(z)
is irreducible over K(t): the coprime rational map z/h has degree d.
Thus the EXACT actual rank is d in a fixed Gaussian dimension FIVE.
The fixed LINEAR bad marker Z has

    [t^(d*j+1)]g=binomial(d*j+1,j)/(d*j+1), j>=0.        (13)

This coefficient is nonzero in characteristic zero.
All other positive coefficients vanish unless their index is1 modulo d.

Over C the geometric monodromy is S_d. The map z/(1+z^d) has one
unramified point z=0 and one point of ramification d-1 at z=infinity
over t=0. Its d finite critical points q=0 are simple and have DISTINCT
nonzero critical values a_i=(d-1)z_i/d. Infinity in time is unramified.
Loops around the finite critical values are transpositions. The product
relation expresses the loop at0 as their inverse product, so the full
transitive group is generated by transpositions and is S_d. This is a
direct branch-cycle argument, not an irreducibility assertion inferred
from a few computed ranks. The connection is the degree-d permutation
representation, with trivial and irreducible standard summands.

For d>=3 its RATIONAL CORE is exactly K[t], despite genuine regular
nonrational failure. A rational pole in t pulls back to a pole at EVERY
point in its fibre. At each critical value a_i only the ramified point
is deleted from(12); its d-2 unramified points are retained. They forbid
a rational pole there. Every other finite fibre retains a point too,
including z=0 at t=0. Therefore a rational element of(12) has no finite
pole. For d=2 each critical fibre is entirely ramified/deleted and the
rational core is K[t,(1-4t^2)^-1].

The global boundary is also explicit. Let J be the intermediate extension
of this finite-monodromy connection on the affine time line. Local tame
decomposition at a ramification point gives one trivial invariant direction
and nonintegral Kummer directions. The latter are already meromorphic in
their minimal extension; deleting the point adds exactly ONE delta in
the invariant direction. The deleted source point at infinity adds delta_0,
and the d deleted simple critical points add delta_(a_i). Thus

    M/J=delta_0 direct_sum (direct_sum_i delta_(a_i)).   (14)

At0 the meromorphic boundary has two invariant directions and the actual
principal-part map selects the hyperplane that is regular on z=0. This
agrees with the universal origin theorem. J is O direct_sum J_standard.
The unit is therefore NOT essential, and M is not globally semisimple,
even though its generic connection is semisimple with FINITE monodromy.
Over a geometric splitting field its composition length is d+3 and
its characteristic cycle is

    d[zero section]+(d-1)[T0*]+2 sum_i[T_(a_i)*].         (15)

These local assertions follow from the displayed tame normal forms:
at a_i rank-minus-invariant-dimension=1, at0 it is d-2;
adding the point quotients(14) gives(15). No Stokes computation is used.

## 6. A genuine Gaussian failure of the generic trace projection

For d=3, the ORIGINAL marker (1+2Y/3)^2 has period q^-2, since

    q^-1=1+(d-1)y/d.

Exact quotient-field linear algebra gives

    Tr_(K(t,z)/K(t))(q^-2)=(9t^3-4)/(27t^3-4).          (16)

Its normalized trace has finite poles, hence is NOT in M by the rational
core result. Thus generic algebraic trace fails to preserve a GENUINE
pure Gaussian actual lattice. The earlier abstract square-root gluing
control in the boundary note is no longer the only illustration.
For every d>=3 the trace of q^-2 has a nonzero simple pole at each
critical value: its two ramified local branches contribute equal leading
terms; all other branches are regular. The same trace failure is general.

## 7. Scope and credits

The family(1) has an exact all-ring image and safe normal form, not just
a collection of examples. Every algebraic formal branch is realizable
with the full corresponding generic field in seven variables; five
variables already realize all ranks>=2 and S_d monodromy, with rational
core O for d>=3. Regular transcendental images, arbitrary transcendental
generic connections, exhaustive Gaussian safe normalforms and remaining
sharp real-cutoff geometry are NOT closed by this construction.

Classical Gaussian contraction/Lagrange inversion, implicit algebraic
encoding, tame finite-cover classification and point-module calculus
receive their ordinary credit. The older ternary even-rank tower and
ternary monomial-recomposition ranks are preserved; this result adds the
REGULAR algebraic branch and its exact ring/trace/lattice structure.
Sabbah/Maisonobe retain convergent rational-slope credit; arbitrary
nonconvergent formal extension and all-ring descent priority remain
separately unassessed. All previous protected endpoints, Oaku operators,
torsion witnesses, specialization thresholds and terminal UNKNOWN jobs
are unchanged. The broad goal remains active.
