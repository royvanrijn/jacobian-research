# Primitive-gamma source coordinates and the remaining Darboux barrier

2026-11-16. Claims HC4MYLJ0/1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
External review/outreach is DEFERRED. The preceding goal turn was PROGRESS
(LI0/1). No registered mathematical checker is replayed.

The primitive-gamma source question left open in LI1 is now CLOSED
POSITIVELY: its defining polynomial is a coordinate of A5, uniformly
at every level. The beta source remains excluded at every level and
after every affine stabilization. These are different charts.

This note gives a general arbitrary-degree coordinate construction and
explicit formulas for the known gamma family. It then identifies two
obstructions to retaining the most natural position pairs in Darboux
coordinates. No standard polynomial Darboux chart, four-variable Hessian
potential, or HC4/JC2 counterexample is constructed. The original LI
proof, context and receipt remain frozen historical evidence; the
previous absence of a gamma trivialization is superseded here.

## 1. LJ0: a general coordinate theorem for a two-generated coordinate ideal

Let k be a field, R=k[w_1,...,w_m,y,z], m>=1, and suppose

    (A,B)=(y,z)=I in R.

Suppose the residue of C in R/I is a polynomial coordinate. After a
polynomial change of the w coordinates, write

    C=a_0*w_1+d(w_2,...,w_m)+a*A+b*B,  a_0 in k*.

Then, for two independent helpers s,t,

    h=s*A+t*B-C

is a polynomial coordinate of R[s,t]. This holds in ALL DEGREES and
with arbitrarily many passive variables. No Keller hypothesis, topology,
Quillen--Suslin theorem, or cancellation conjecture is needed.

Here is a constructive proof. Choose a polynomial 2-by-2 matrix M with

    (A,B)^T=M*(y,z)^T.

On I/I^2, the two generating sets are bases. Thus M modulo I is
invertible over R/I=k[w], and det(M) modulo I is a FIELD UNIT delta.
Write det(M)-delta=r*A+s_0*B. The determinant can be corrected without
changing the two polynomials, using row syzygies as follows.
For row corrections h_0*(z,-y), k_0*(z,-y), the exact identity is

    det(M_corrected)=det(M)+h_0*B-k_0*A.

Take h_0=-s_0 and k_0=r. Then det(M_corrected)=delta.
Denote this matrix by Mhat. Its inverse is polynomial, adj(Mhat)/delta.
This rank-one determinant correction is the SAME algebraic identity
already used in retained HH0/JV2 for oriented generating pairs; its
coordinate-ideal application here is not claimed as a new orbit theorem.

Now make the polynomial helper change

    (p,q)=(s-a,t-b)*Mhat.

Its polynomial inverse is (s,t)=(p,q)*Mhat^-1+(a,b), with w,y,z fixed.
In these coordinates

    h=p*y+q*z-a_0*w_1-d.

This eliminates w_1 linearly. More strongly (w_2,...,w_m,y,z,p,q,h)
is an ambient polynomial coordinate system, with

    w_1=(p*y+q*z-d-h)/a_0.

Every level h=c is therefore A^(m+3), with an actual polynomial inverse.
The hypothesis is equality of the SCHEME IDEALS, not merely that one
selected component or the normalization of the zero set is a line.
The offset C modulo I being a coordinate is also essential.

## 2. Explicit uniform coordinates for the known primitive-gamma family

Retain the KNOWN Alpoge map (A,B,C) and its Jacobian -2, generic degree3
and three-point collision from KO2/KW. Put

    u=1+x*y,  gamma_0=2-3*x*y-x^2*z,
    A=u^3*z+y^2*u*(1+3*u),
    B=y+3*x*u^2*z+3*x*y^2*(1+3*u),
    C=x*gamma_0.

The primitive-gamma moment polynomial on A5 is

    h=2*alpha*A+beta*B-C.

The previous Hodge computation proves only that the zero set of A,B is
the coordinate line; we now supply SCHEME-IDEAL identities and a complete
ambient coordinate change. All displayed expressions are polynomials.
Define

    D=3*u*gamma_0-3,
    r=(3*x/2)*(1-D),  s=(u*D-1)/2,
    T=y*u*(1+3*u)+x*(u^2+u+1)*z.

Direct identities give

    u*B-3*x*A=u*y,
    r*A+s*B=y,
    A-T*y=z.

Thus (A,B)=(y,z) as ideals, not only as reduced sets. In particular

    z=(1-T*r)*A-T*s*B,
    C=2*x+a*A+b*B,
    a=-3*x^2*r-x^3*(1-T*r),
    b=-3*x^2*s+x^3*T*s.

For an entirely explicit corrected matrix set

    M=[[y*u*(1+3*u),             u^3],
       [1+3*x*y*(1+3*u),     3*x*u^2]],
    J=x*(u^2+u+1),  H_0=J*s,  K_0=-J*r,
    Mhat=M+[[H_0*z,-H_0*y],[K_0*z,-K_0*y]].

Here det(M)=-u^3, Mhat*(y,z)^T=(A,B)^T and

    det(Mhat)=-u^3+J*(r*A+s*B)=-1.

Define new helpers

    (p,q)=(2*alpha-a,beta-b)*Mhat.                 (1)

This is an ambient polynomial automorphism with x,y,z fixed. Its inverse
is explicitly

    (v_1,v_2)=(p,q)*Mhat^-1,
    alpha=(v_1+a)/2,  beta=v_2+b.                 (2)

The full moment polynomial becomes

    h=p*y+q*z-2*x.                               (3)

Consequently (y,z,p,q,h) are polynomial coordinates of A5. On h=c the
inverse is (2) evaluated at

    x=(p*y+q*z-c)/2.                             (4)

This proves Y_c=A4 for EVERY c, and a globally trivial polynomial
family over the level line. No localization, square root, selected open
subset, or affine stabilization is used. There is no degree-bounded
extrapolation: the maps (1)--(4) are exact identities.

## 3. LJ1: what this gives, and what it does not give, symplectically

The target gamma section has coordinates (A,B,alpha,beta), with C
eliminated by 2*alpha*A+beta*B-C=c. Its standard symplectic form is

    omega_0=d(alpha) wedge dA+d(beta) wedge dB.

Using (2)--(4), the projection

    Phi_c:A4_(y,z,p,q)->A4_(A,B,alpha,beta)

is a polynomial ETALE map, the same incidence base change of the known
map as before. It has generic degree3: the target incidence projection
to A3_(A,B,C) is dominant and meets the finite rank3 open, and its source
has now been proved irreducible. For EVERY c the retained three-point
collision over (A,B,C)=(-1/4,0,0) survives by taking alpha=-2*c,beta=0.
The known collision is reused, not recomputed.

The Jacobian in the displayed source and target coordinate orders is
1/2. This follows WITHOUT REPLAYING the old Jacobian calculation:
the ambient map (x,y,z,alpha,beta)->(A,B,alpha,beta,h) has determinant2
from the retained -2 determinant; (x,y,z,alpha,beta)->(y,z,p,q,h) has
determinant4, since the helper determinant is -2 and eliminating x
has determinant -2. Divide these constants and restrict to h=c.

The pulled-back form Omega_c=Phi_c^*omega_0 is a polynomial EXACT
symplectic form on A4, with primitive alpha*dA+beta*dB. In the displayed
coordinates Omega_c^2/2=-(1/2)*dy wedge dz wedge dp wedge dq.
This is an actual affine-space source; it is NOT a supplied polynomial
identification of Omega_c with a chosen standard source symplectic form.
A polynomial Keller map on A4 alone is not an HC4 gradient map.

### The source pair y,z cannot be retained as Darboux positions

This obstruction holds at EVERY level. On the original cotangent lift,
the target momentum gamma has torus weight+1. Invariant lifts of the
section coordinates y,z are y/gamma and z/gamma^2. If D_C is the
polynomial target-dual source derivation with D_C(C)=1,D_C(A)=D_C(B)=0,
the ordinary cotangent Poisson rule gives

    {y/gamma,z/gamma^2}
       =gamma^-4*(y*D_C(z)-2*z*D_C(y)).             (5)

These functions are moment-invariant, so their bracket descends to the
reduction. On gamma=1, this is exactly {y,z}_(Omega_c). Evaluate the
source base at (x,y,z)=(0,1,0). Here

    grad A=(7,8,1),  grad B=(12,1,0),
    D_C(y)=-6,  D_C(z)=89/2,
    {y,z}_(Omega_c)=89/2 !=0.

This base point occurs on EVERY Y_c, since in (4) take y=1,z=0,p=c
and, for example, q=0. Thus no Darboux coordinate system for Omega_c
can have y,z as its two commuting positions, regardless of how its
two momenta are changed. The computation controls a bracket of the
new reduction, not a replay of the known map determinant or collision.

### The target pair A,B cannot be retained as global Darboux positions

This is a general elementary obstruction. Suppose a polynomial
symplectic form is d(alpha) wedge dA+d(beta) wedge dB and there are
GLOBAL polynomial coordinates (A,B,U,V) with the same form equal to
dU wedge dA+dV wedge dB. Comparing coefficients of dU wedge dA,
dV wedge dA,dU wedge dB,dV wedge dB gives

    alpha=U+a_1(A,B),  beta=V+b_1(A,B).

Hence (A,B,alpha,beta) are themselves global polynomial coordinates.
The map Phi would be a polynomial automorphism. Its retained collision
rules this out for every Phi_c. Therefore standardizing Omega_c, if
possible at all, must also change its target position pair.

These two barriers do NOT exclude a Darboux transformation mixing all
four coordinates. No such transformation, no quantization, and no
four-variable Hessian potential is constructed or ruled out here.

## 4. Provenance, exploration and the full unresolved target

The determinant-correction identity is credited to the retained HH0
argument and its Guyot orbit-theorem context; the present proof is
elementary and does not import another classification theorem. The
known map, torus weights, reduced form and collision are retained from
KO2 and the archived natural-torus note. LI1's beta exclusion and gamma
axis calculation are retained byte-for-byte. No mathematical premise
from the newly browsed iterated-image literature is imported: stable
surjective images alone do not supply properness or trace regularity.

The new checker verifies only the new ideal identities, corrected
matrix, polynomial inverse schemas and moment identity, along with the
new reduced Poisson calculation and retained byte pins. The universal
coordinate theorem is a written proof, not a sampled degree claim.
No independent replay, formal verification or external review is claimed.

The positive source trivialization removes ONE obstruction to a concrete
counterconstruction route. It neither proves ordinary trace-pairing
regularity on an actual plane nor supplies the missing standard source
symplectic/gradient structure. Full HC4/JC2 remains OPEN. Published full
plane degree-four exclusion remains in force; these degree3 maps have
FOUR-dimensional sources and do not reopen any quartic plane queue.
