# Source-line conductors and finite primitive-order collisions

HC4MYKW0/1. Full HC4/JC2 remains OPEN. This note gives a universal
actual-source conductor constraint and an exact higher-dimensional
countercontrol to a local primitive-hypersurface argument. It constructs
no plane Keller counterexample and no four-variable unit-Hessian example.

## 1. Retained input and the universal conductor theorem: KW0

Let (P,Q):A2->A2 be an ACTUAL complex polynomial Keller map in source
coordinates (t,z), with j=Jac_(t,z)(P,Q) a nonzero constant. Choose the
linear source coordinate t so that P(t,z) is monic of positive degree M
in z, and t is primitive over C(P,Q). PF2OMIT3 ALREADY proves that such
coordinates exist outside a finite exceptional set, that the source
R=C[t,z] is finite over C[P,Q,t], and that its omitted normalization
boundary is exactly the global pole fibre of t:N->P1. These are PRIOR
results, not new consequences claimed here.

Put D=C[t,P], A=C[t,P,Q], and let

    Phi(T;p,q)=0

be the primitive irreducible inverse equation. Normalize it to be MONIC
in q, of degree M: Q is integral over the normal D, and primitivity of
t makes Q a primitive element of L/Frac(D). The same polynomial is the
primitive minimal relation of t over C(p,q); it need not be monic in T.
Set the SOURCE polynomial

    h(t,z)=Phi_T(t;P(t,z),Q(t,z)).                    (1)

Then the FULL conductor ideal, including every finite source point, is

    (A:R)={r in R:r*R subset A}=h*R.                (2)

This upgrades the divisorial inverse-conductor identity already retained
in PF2OMIT4 to an exact global ideal for the source-finite chart. The
trace-dual method and the distinction between primitive collisions and
true ramification are credited there; neither is new here.

For EVERY c in C, put p_c(z)=P(c,z), q_c(z)=Q(c,z). The polynomial
Phi(c;p,q) is reduced and irreducible, and its normalization is precisely
the finite IMMERSIVE BIRATIONAL map A1_z->H_c it defines. Its curve
conductor is exactly

    (C[p_c,q_c]:C[z])=h(c,z)*C[z].                 (3)

Thus conductor SPECIALIZATION IS EXACT in this actual-source family,
with no missing finite branch or exceptional-parameter correction. If
delta_c=dim_C(C[z]/C[p_c,q_c]) is its total AFFINE normalization defect,

    deg_z h(c,z)=2*delta_c                         (4)

for every c. The polynomial h(c,z) is never identically zero. Its zeros
are exactly the finite points whose image has ANOTHER source point on
the SAME line t=c; no cusp branch or omitted point accounts for them.

Consequently, under noninvertibility, Gwozdziewicz's ALREADY IMPORTED
one-line theorem gives

    deg_z h(c,z) is a POSITIVE EVEN integer for EVERY c. (5)

In particular h is nonconstant, has no vertical source-line factor, and
its zero curve maps surjectively and quasi-finitely to A1_t. The generic
degree is an even upper bound for every special degree; drops can only
be even and occur at finitely many c. This is an all-degree necessary
constraint, not an exclusion of every polynomial satisfying it.

## 2. Exact trace-dual proof of the conductor ideal

The elementary monogenic trace-dual formula over the polynomial base D
gives, for the trace of L/Frac(D),

    A^vee=(1/Phi_q)*A,
    R^vee=(1/P_z)*R.                               (6)

Both algebras are finite free of rank M over D: use the power bases of
Q and z respectively. The generic trace is nondegenerate. These are the
standard complementary modules already used in PF2OMIT4 and KU0.
They also follow directly from Lagrange interpolation at the split roots:
tracing a polynomial divided by the monic derivative extracts its top
remainder coefficient, and the resulting dual power basis is triangular
with unit diagonal. This identifies the entire modules over D, not merely
their height-one localizations.

Here is the conductor argument, written globally rather than only at
height-one DVRs. For any finite normal overorder R of the monogenic A,

    (A:R)=Phi_q*R^vee.                             (7)

Indeed, if r*R subset A, then (r/Phi_q)*R subset A^vee and tracing shows
r/Phi_q belongs to R^vee. Conversely R^vee is an R-module and is contained
in A^vee, so r=Phi_q*q with q in R^vee has r*R subset Phi_q*A^vee=A.
Equation(6) now gives conductor=(Phi_q/P_z)*R. The quotient is regular
by (7), including at zeros of P_z.

Differentiate the identity Phi(t;P,Q)=0 with respect to t,z. Since the
Keller Jacobian j is constant and nonzero, solving the two equations gives

    Phi_p(t;P,Q)=-h*Q_z/j,
    Phi_q(t;P,Q)= h*P_z/j.                         (8)

Hence (7) is (h/j)*R=h*R, proving (2) at EVERY source point. No assumption
that P_z itself is invertible was made. Without the whole Keller
condition the factor is h/J(t,z), and zeros of h can also reflect true
ramification rather than normalization conductor; Section5 controls this.

## 3. Every line, every parameter, and the even conductor degree

Since P is monic of z-degree M, each p_c is monic of that SAME degree.
The line parametrization is finite onto its image. Its normal image curve
is A1: its function field is rational, while a finite map from A1 makes
its coordinate units scalar and leaves only one place at infinity.
Factoring through that A1 normalization gives a polynomial s=s(z).
The invertible Keller differential makes (p_c',q_c') nowhere simultaneously
zero. Thus s'(z) cannot vanish. In characteristic zero this makes s
linear, so the parametrization is birational and immersive for EVERY c.
This is the same elementary coordinate-line argument retained in PF2RULED53.

The minimal equation of q_c over C(p_c) therefore has degree M. The
specialization Phi(c;p,q) remains monic of q-degree M and vanishes on
the image, so it IS that irreducible equation. In particular it cannot
acquire another component or a repeated factor at a special c. Applying
the global argument (6)--(8) over C[p_c] proves (3); conductor formation
has not been commuted with specialization without this proof.

For completeness, (4) also follows from exact trace lattices, without
an additional plane-curve duality import. Let E_c=C[p_c], U_c=C[p_c,q_c]
and N_c=C[z]. These are free rank-M E_c-lattices, and N_c/U_c has finite
C-length delta_c. A power-basis inclusion matrix L_c has

    deg_p det(L_c)=delta_c

by Smith form over C[p]. Their trace determinants differ by det(L_c)^2.
The monogenic determinant formulas and (8) therefore give

    Norm_(C(z)/C(p_c))(h(c,z)/j)
            = nonzero_constant * det(L_c)^2.

The degree of that norm in p is the degree of h(c,z) in z: a finite
polynomial map A1_z->A1_p pushes its finite zero divisor forward without
changing the total degree over C. Taking degrees proves (4). The standard
conductor-degree relation was already used in the specific target-curve
family PF2K1CC1; this supplies the universal source-family proof.

An immersive normalization has smooth individual branches. A singular
image point here has at least two finite preimages; conversely such
preimages make its local ring nonnormal. Thus the conductor support is
exactly the line-collision locus stated after (4). If delta_c=0, the line
map is an embedding and hence injective. Gwozdziewicz then gives the
Keller inverse. Under noninvertibility delta_c>0 for EVERY c, proving (5).
This use of the one-line theorem is credited, not a new independent proof.

## 4. Finite collisions need not be hypersurface smoothness: KW1

The following control uses the KNOWN three-dimensional whole-source
Keller map already credited and checked in KO2. It is NOT a new Jacobian
counterexample and is NOT an HC4 counterexample. Put u=1+x*y and

    A=u^3*z+y^2*u*(4+3*x*y),
    B=y+3*x*u^2*z+3*x*y^2*(4+3*x*y),
    C=2*x-3*x^2*y-x^3*z.

Its Jacobian is -2. The actual source points

    (0,0,-1/4), (1,-3/2,13/2), (-1,3/2,13/2)

all map to b=(-1/4,0,0). The source COORDINATE z is a primitive INTEGRAL
algebraic function over C(A,B,C), of degree3. An explicit compact monic
equation is as follows. With independent target variables A,B,C and w,
write

    P=B*C, Q=A*C^2,
    gamma=(P-4*w+3*w^2)/2,
    f=w^3-2*w^2+P*w-2*Q,
    g=8*(C^2*Z-5*gamma^2+3*w*gamma+gamma^3),
    H(A,B,C,Z)=Res_w(f,g)/(512*C^6).                (9)

The displayed quotient is an EXACT polynomial, monic of Z-degree3.
Its full coefficients are retained in the checker receipt. On the source,
gamma=C/x generically, w=u*gamma, and
C^2*z=5*gamma^2-3*w*gamma-gamma^3. Together with the KO2 irreducible
cubic f, these identities show H is the FIELD NORM of Z-z.

At b it specializes to

    H(b,Z)=(Z+1/4)*(Z-13/2)^2.                     (10)

If z were in the target field, its norm would be (Z-z)^3; its regular
Z^2 coefficient would make z a target polynomial, and every specialization
would be a cube. Equation(10) rules that out. Since the field degree is3,
z is primitive and H is its irreducible minimal equation. Thus (9) gives
a connected, genuine primitive order, not a split formal substitute.

At (A,B,C,Z)=(-1/4,0,0,13/2), H AND ALL FOUR PARTIAL DERIVATIVES vanish.
The primitive hypersurface is singular at the actual affine source image,
although the original Keller differential is invertible at BOTH distinct
finite preimages having z=13/2. Individual local inverse graphs are
smooth; their union is not a smooth hypersurface there. No omitted point
or infinite source coordinate is required for this collision.
In fact all three degree-three sheets are retained at b: in the completed
full finite normalization the three etale source points give three rank-one
clopen factors. They exhaust its generic rank3; a remaining factor would
have rank0, impossible for this torsion-free finite target module. Thus
the full normalization is etale over b and has no omitted point there.

This also gives a genuine whole-Keller control of KU0's order distinction.
The trace-dual algebra of the finite free monogenic ORDER C[A,B,C,z] is
its localization at H_Z. It therefore deletes the displayed double
primitive image and BOTH actual etale source points above it. The full
normalization retains them. Even an integral primitive SOURCE COORDINATE
does not make its order the right input for source-etale reconstruction.

The coordinate z also gives a concrete instance of the already known
failure of the UNRESTRICTED higher-dimensional D-algebra nonembedding
conjecture: z is not a target polynomial, yet its target derivatives
and their algebra lie in C[x,y,z], because the whole Keller inverse
Jacobian matrix has polynomial entries. No finite-generation assertion
for this three-dimensional D-algebra is needed or inferred from KV0/1,
whose source-layer hypotheses are PLANE-specific. The plane version and
HC4 remain untouched by this higher-dimensional countercontrol.

## 5. A small connected control and the exact literature boundary

For a low-cost transcription control, take the non-Keller plane map

    (t,z)->(a=t-z^2,b=-t*z).

The source coordinate t is primitive, with monic relation
t^3-a*t^2-b^2=0 and finite normalization C[t,z]. At (a,b,t)=(-1,0,0),
the two points (t,z)=(0,1),(0,-1) have nonzero source Jacobian -2,
but the irreducible primitive hypersurface is singular. For the equation
Phi=b^2-(T-a)*T^2, h=Phi_T(t;a,b)=t*J and its TRUE normalization
conductor is t*R, not h*R. The critical-Jacobian factor distinguishes
this control from the actual Keller theorem (2). Its t=0 line map is
two-to-one; the all-parameter immersive-birational proof in Section3
therefore cannot be applied to it.

Makar-Limanov's primary MPIM22-11 Section4 (manuscriptp5,PDFp7) attempts
to infer nonsingularity of the primitive hypersurface at finite source
images from the independent linear parts of the Keller outputs, then
attributes possible singularities to infinity. That LOCAL INFERENCE
does not account for several finite branches. The small connected control
refutes the local argument; (9)--(10) exhibits the issue even with a
GLOBAL constant Jacobian in dimension3. This is a correction of that
inference, not a claim that a plane Keller counterexample exists.

Section2's nonsingular-primitive-hypersurface reciprocal argument remains
valid, as does its definition of D(y). Its separate element-existence
regularity conjecture was never proved there. Smoothness of the whole
normal derivative algebra from KV0/1 still does not supply the missing
primitive hypersurface element. Section1--3 here make the finite collision
obstruction explicit in every hypothetical noninvertible plane source:
EVERY source-line member has positive even affine conductor degree.

FULL THEOREM TARGET: full HC4/JC2 proof or an ACTUAL counterexample.
CURRENT LAST UNSOLVED IMPLICATION is unchanged:

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.          OPEN

The conductor theorem is a necessary constraint, not the missing proof.
It does not make the polynomial h a unit, force one line injective, or
prove full marked canonical-energy integrality. Quartic realization
remains retired by the published Domrina theorem. External review and
outreach remain DEFERRED, with all assurance flags unchanged.
