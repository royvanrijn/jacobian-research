# Jacobian mate rigidity and quadratic cyclic target exclusions

PF2RULED145. 16 November 2026. ALL degrees over C.
This closes an unbounded target-coordinate family on the surviving
QUADRATIC cyclic source images, without assuming target-field invariance.
The companion target function has NO degree bound. It also upgrades
PF144's field escape to target-coordinate pencil escape.
Unrestricted cyclic targets, generic first-velocity integrality and full
HC4/JC2 remain OPEN. No actual counterexample is constructed.

## 1. Rational constants of a genuine Keller coordinate

Let B=C[P,Q] subset R=C[x,y], Jac(P,Q)=1, L=C(x,y), F=C(P,Q).
For EVERY polynomial target coordinate P (after any target rechart),

    {H in L:{P,H}=0}=C(P),
    C(P) intersect R=C[P].                             (1)

Here {P,Q}=1 uses the actual constant-Jacobian volume.

We reproduce the needed primitive-coordinate argument. Let K0 be the
relative algebraic closure of C(P) in L; it is finite over C(P).
Normalize C[P] in K0. Its finite normal curve algebra A lies in R,
since its elements are integral over C[P] and R is integrally closed.
Thus A has scalar units. Restrict the resulting dominant A2->Spec A
to a source affine line on which P is nonconstant. This makes the
smooth affine curve Spec A rational. A smooth rational affine curve
with scalar units is A1 (the rational-curve argument already used in MR0).
Write A=C[h], h in R and P=f(h) with f polynomial. If deg f>1,
a root c of f' is reached by h on A2, and dP=f'(h)dh vanishes there,
contradicting the actual Keller differential. Hence K0=C(P).

The derivation D={P,-} kills K0. Conversely, a D-constant H transcendental
over C(P) would make L finite algebraic over C(P,H), and D would vanish
by separability, contradicting D(Q)=1. This proves the first equality.
The second follows from PF52's F intersect R=B, or directly from
surjectivity of the nonconstant polynomial P and denominator poles.

It follows that ALL polynomial Jacobian mates of a fixed genuine
Keller coordinate differ by a scalar multiple and a polynomial in it:

    {P,Q'}=c in C* => Q'=c*Q+h(P), h in C[P].          (2)

This is a uniqueness statement, not an existence theorem for a mate.

## 2. One target pencil already forces full field descent

Let V be one of PF144's actual cyclic source images with
C=O(V) subset R and unit volume Omega pulling back to dx wedge dy.
Let rho be a regular automorphism of V with rho^*Omega=chi*Omega,
chi in C*. This holds for its fibre scalings and all their conjugates.

Suppose an actual target pair P,Q in C exists. If rho preserves the
SINGLE coordinate pencil C(P), then rho(P)=a*P+b, a!=0:
both rho(P) and rho^-1(P) lie in R intersect C(P)=C[P], so they are
mutually inverse polynomials. The coframe equation now gives

    {rho(P),rho(Q)}=chi,
    {P,rho(Q)-(chi/a)*Q}=0.

Equation (1) and polynomial source regularity force

    rho(Q)=(chi/a)*Q+h(P), h polynomial.              (3)

Thus rho(B)=B and rho(F)=F. The constant Jacobian has supplied the
missing passage from ONE preserved target pencil to full equivariance.

More generally, for finite-order rho, if rho sends ANY polynomial
target coordinate to ANY polynomial target coordinate of B, it preserves
B. Complete rho(P) to a target pair rho(P),T in B, with bracket one.
Then rho(Q)-chi*T centralizes rho(P), so belongs to C[rho(P)] by (1).
Hence rho(B) subset B; finite order turns this inclusion into equality.

Apply PF144. No nonidentity regular fibre scaling on MX0 or on the
classical S_g can send a target coordinate to another target coordinate.
In particular NO target coordinate can be a fibre eigenfunction plus
a constant, even with a completely unrestricted companion.

The same holds for EVERY automorphism-conjugate of these fibre groups:
compose a hypothetical target with that automorphism and apply PF144
to the same whole plane source. This is stronger than requiring the
original two target generators to form an invariant span. It does NOT
say every polynomial target coordinate has a preserved pencil.

## 3. The FULL quadratic cyclic source-image algebra

Specialize MX0 to g=2. Retain every k>=2 and every positive-degree
square-free b(u), b(0)!=0. Put

    sigma(u)=-b(u)+2*(1-k)*u*b'(u),
    M(u)=u^2*b(u)^k, N(u)=u*sigma(u)^2.

The ENTIRE quotient V has the presentation

    C=O(V)=C[u,w,z]/(w^2-N(u)-M(u)*z).               (4)

This is not a partial primitive order. Indeed MX0's cover is
t^2*b(t^2)^k*z=v^2-sigma(t^2)^2 and the deck involution is
(t,v,z)->(-t,-v,z). Reduce any polynomial modulo the monic v^2 relation.
Invariance makes its part without v even in t, and its coefficient
of v odd. Thus invariants are exactly A(u,z)+w*B(u,z), u=t^2,w=t*v,
with the sole relation (4). This also proves completeness scheme-wise.

The source chart has

    u=t^2, w=t*(sigma(t^2)+t^2*b(t^2)^k*x).

MX0 supplies the whole source, surjectivity and etaleness. The reduced
multiple-fibre curve C0={u=w=0} is reached isomorphically by t=0.

Consider ANY proposed first target coordinate

    P=c*z+A(u)*w+B(u), c in C, A,B arbitrary polynomials. (5)

If c=0, P is constant on C0. In any actual Keller pair a complementary
Q then restricts to an etale polynomial on the whole source line t=0,
hence an affine-linear function. The retained one-line theorem forces
inversion, contradicting the degree-two source factor. So c=0 is excluded.

For c!=0 rescale P to c=1. Define

    E(u)=B(u)-M(u)*A(u)^2/4.

The polynomial automorphism

    tau_h(u)=u, tau_h(w)=w+M*h,
    tau_h(z)=z+2*h*w+M*h^2                         (6)

preserves (4), with inverse tau_-h. Taking h=A/2 sends z+E to P.
Thus the critical scheme of P equals that of z+E under an automorphism.

## 4. Every such moving-label coordinate has an explicit critical scheme

On D(M), the surface coordinates are u,w, and

    z+E = w^2/M + D(u), D(u)=E(u)-N(u)/M(u).

A critical point is exactly w=0,D'(u)=0. Undoing (6), its coordinates
for P are w=-M*A/2, z=M*A^2/4-N/M.

There is NO additional critical point over M=0. At a root of b,
w!=0 and u,z are local coordinates; the derivative of z+E in z is one.
At u=0, N has a simple zero and M a double zero; w,z are local coordinates
and again that derivative is one. This checks the whole surface.

Write r=deg b>=1 and put

    T(u)=u^2*b^(k+1)*E'
          +sigma^2*(b+k*u*b')-2*u*b*sigma*sigma'.   (7)

Direct differentiation gives D'=T/[u^2*b^(k+1)]. Moreover

    T(0)=b(0)^3 !=0,
    T(beta)=k*beta*b'(beta)*sigma(beta)^2 !=0
       for EVERY root beta of b.                    (8)

Thus gcd(T,M)=1 and the COMPLETE critical algebra of P is

    O(Crit P)=C[u]/(T).                              (9)

This includes multiplicities and all closed critical points; there is
no discarded bad-fibre or isolated-point case.

Its EXACT length is

    3*r,                    if E is constant;
    (k+1)*r+m+1,            if m=deg E>=1.           (10)

For the first case sigma has degree r with leading coefficient
b_lead*(-1+2*(1-k)*r), nonzero. The last two terms of T have combined
leading coefficient
b_lead*sigma_lead^2*(1+(k-2)*r), nonzero, and degree 3r.
For nonconstant E the first term of T has degree
(k+1)*r+m+1>3r, so cannot cancel.

Consequently EVERY function (5), c!=0, has a nonempty finite critical
scheme. It cannot be a coordinate of ANY etale plane target. Together
with c=0, this excludes the ENTIRE family (5), in all polynomial
A,B degrees and for every possible companion Q. Nonzero A generally
breaks the original fibre scaling; target invariance was never assumed.
For example k=2,b=u-1,A=1,E=0 gives P=z+w+M/4, with the explicit
critical point (u,w,z)=(1/3,-2/81,1/81).

The geometric mechanism is the unavoidable finite critical point of a
rational function with the fixed poles of N/M: a simple pole at zero
and poles of order k at every b-root. PF137's pole-count method is credited,
but (7)--(10) give the exact critical scheme for the MOVING labels which
passed its different exact-volume obstruction.

## 5. The canonical-trivial exceptional S2 also has no such target

For the classical surface S2={w^2=u+u^2*z}, replace M,N by u^2,u.
The same square completion reduces (5), c=1, to z+E(u), with

    E=B-u^2*A^2/4, D=E-1/u, D'=(1+u^2*E')/u^2.     (11)

If E is nonconstant of degree m, Crit(P) has exact length m+1,
away from u=0. If E is constant, P is nonsingular and conjugate to z
plus a constant. A no-critical-points test ALONE therefore does not
exclude this final case.

The actual coframe excludes it. Normalize the unit volume as
Omega=du wedge dw/u^2. All units are scalar by pullback to the whole
plane, so any etale target coframe is a nonzero constant times Omega;
rescale its companion. If a RATIONAL Q satisfied d z wedge d Q=Omega,
then on the geometric generic fibre z=c!=0 it would satisfy

    dQ=-du/(2*w), w^2=u+c*u^2.                     (12)

Adjoin a square root a of c and set
H=1+2*c*u+2*a*w. Its inverse is 1+2*c*u-2*a*w. The generic conic
has function field C(c,a)(H), and

    u=(H+H^-1-2)/(4*c), w=(H-H^-1)/(4*a),
    -du/(2*w)=-dH/(2*a*H).                         (13)

The last differential has nonzero residues at H=0 and infinity.
A rational derivative has zero residue at every point, so no rational
Q exists. This remains true after any finite field extension by trace:
an exact pullback would give an exact differential downstairs after
normalized trace. Thus even a rational mate is impossible.

The automorphisms (6) preserve Omega. Hence the nonsingular E-constant
case also has no regular Jacobian mate. The c=0 source-line argument
works as before. Therefore (5) is excluded on S2 for ALL A,B,c as well,
with no target-field symmetry hypothesis.

## 6. Scope, transport and the full normalization

The exclusion applies if ANY polynomial coordinate in ANY target
rechart has form (5). It also survives arbitrary automorphisms of V:
a hypothetical target with such a transformed coordinate can be
precomposed with the inverse surface automorphism. No claim classifies
all regular coordinates or all nonsingular functions on these surfaces.

In particular functions with nonconstant z coefficient, higher nonlinear
dependence on z, or general w*z terms have not been exhausted.
One cannot replace them by (5) without a proved reduction.

For a hypothetical original Keller map retain B,R,F,L, the FULL finite
normalization S and its full trace-dual Lambda, including EVERY ramified
boundary prime. V is the actual relative source image, not S or Spec S.
The contradictions above use critical points of a proposed actual
coordinate, or its forced coframe on a complete generic curve.
They neither delete ramified primes nor claim marked trace regularity.

PF142's mixed 2+3 minimal case stays retired. PF143's marked quadratic
integrality implication remains UNSUPPLIED. PF144's field-stable targets
are already excluded; this note removes the larger concrete quadratic
coordinate family even when the target breaks the original scaling.
The full HC4/JC2 theorem and an actual counterexample remain OPEN.

## 7. Evidence

The universal arguments are written above. New cheap exact controls
support the coframe formula, the COMPLETE quadratic invariant algebra,
square completion, critical numerator and exact length formulas,
the nonsingular S2 control, and its nonzero-residue generic conic.
Finite tests do not classify arbitrary target functions or prove HC4/JC2.

MR0's normal-curve/local-slice foundations, PF52's field intersection
and one-line theorem, PF137's pole method, MX0's actual source-image
construction, and PF142--144 are retained with their original assurance.
No predecessor checker, new large search or new external theorem import
is used. Independent whole-proof replay, formal verification,
external review and novelty priority are not claimed.
