# Every plane Keller extension admits Pascal-nonterminating HC4 normal forms

HC4MYMK0. 16 November 2026. ALL original plane Keller maps, ALL generic
sheet degrees, and two unbounded integer framing parameters.

Every actual plane Keller pair can be reframed, without changing its
embedded target ring, full finite normalization or trace-dual module, so
that its cotangent HC4 potential has a normalized gradient with exponential
iterate growth and a NONTERMINATING Pascal series for its canonical energy.
This holds for automorphisms too. It is not a counterexample to HC4 or to
the formal Legendre identity. Full HC4/JC2 remains OPEN.

## 1. Primary-source screen and the distinction being tested

Adamus--Hajto,
[*The Legendre potential and the gradient case of the differential Galois criterion for the Jacobian conjecture*, arXiv:2610.36736v1](https://arxiv.org/html/2610.36736v1),
Sections3--6, give energy membership and Nambu-tower criteria, identify
the doubled criterion with its base version, and study the Pascal series.
Proposition6.10 and Remark6.11 ALREADY give an invertible four-variable
example with nonterminating coordinate and energy Pascal differences.
That negative result is credited to them; it is not a new discovery here.
The precise pertinent statements and proofs were read. The earlier
differential-Galois foundations cited by them were NOT independently
reproduced or imported as new HC4 theorems.

The result below is a uniform framing statement for EVERY plane extension,
retaining its ACTUAL target subring. It extends the retained MG0 distinction
between target reframing and dynamical conjugacy to normalized gradient
maps and their scalar energy differences. All formulas needed for this
framing and its nontermination are proved directly.

## 2. Arrange a useful top term while keeping the first derivative identity

Let K:A2 -> A2 be any complex polynomial Keller map. Fix a source point x0
and put A0=DK(x0). Let q=deg K>=1 and let T_q be the nonzero degree-q top
vector of A0^-1*(K(x0+z)-K(x0)).

Choose v!=0 with T_q(v)!=0. There is a covector ell with

    ell(v)=1,       ell(T_q(v))!=0.

Indeed if T_q(v) is proportional to v its nonzero proportionality factor
works for every such ell; otherwise avoid one point in the affine line of
covectors with ell(v)=1. Choose S in GL2 with second column v and first
column in ker ell. Its second inverse row is ell. In the actual affine
source coordinates x=x0+S*z, write

    C(z)=(P(z),Q(z))=S^-1*A0^-1*(K(x0+S*z)-K(x0)).

Then C(0)=0, DC(0)=I, det DC=1, and

    deg Q=q,       Q_q(e2)!=0,       deg P<=q.              (1)

This is a choice of source coordinates and an invertible affine TARGET
change. The generated subring inside the intrinsic original source ring
is still B=C[K]. Source-coordinate notation does not replace the actual
source by an abstract etale surface.

For ANY integers M,N>=2 apply the two target shears

    T_(M,N)(p,q)=(p+q^M, q+(p+q^M)^N),
    G=T_(M,N) composed with C.                             (2)

Their inverse sends (u,v) to

    q=v-u^N,       p=u-(v-u^N)^M.

Thus G(0)=0, DG(0)=I, det DG=1 and C[G]=B. Put

    D=M*N*q,       H=Q_q^(M*N),       v0=e2,
    a=H(v0)!=0.

The full top term is

    G_D=v0*H.                                             (3)

This framing is NOT dynamical conjugacy of K. Its purpose is precisely
to change degree dynamics while keeping B, its fraction field, its full
normalization S_B in the source field, and its trace-dual Lambda LITERALLY
unchanged. Sheet degree d and the full original trace-pairing problem are
also unchanged. Arbitrarily large D can be chosen.

The two shears in(2) have IDENTITY derivative at zero. This matters:
postcomposing by a high-degree shear with a nonidentity linear part and
then normalizing that part can destroy its degree growth. We are not
silently applying MG0's earlier one-shear framing to a normalized gradient.

## 3. The cotangent potential and an actual gradient normal form

Let z=(z1,z2), w=(w1,w2), A(z)=DG(z), and put

    h(z,w)=w dot G(z),
    g_h=(z,w) dot grad(h)-h=w dot A(z)*z.                  (4)

Then det Hess(h)=1 and the Hessian at zero is the constant swap matrix

    J=[0 I; I 0].

As in retained KO3, grad(h) is invertible iff the original K is. Its
source change (z,w)->(z,xi=A^t*w) identifies the gradient extension with
the original base extension after adjoining two independent xi variables.
The original finite trace-pairing issue has not disappeared in the lift.

Consider the tangent-to-identity polynomial map

    F=J*grad(h)=(G(z), A(z)^t*w).                          (5)

It is LINEARLY CONJUGATE to a genuine gradient in Euclidean Hessian normal
form. Explicitly let i^2=-1 and set

    U=[I iI; I/2 -iI/2],       U^t*J*U=I,
    psi(X)=h(U*X),       Phi=grad(psi)=U^-1*F*U.            (6)

Now psi(0)=0, grad psi(0)=0, Hess psi(0)=I and det Hess psi=1. Its energy
is g_psi=g_h composed with U. This normalization IS a linear dynamical
conjugacy of F, unlike the target reframing of K in Section2. Degrees and
Pascal-difference vanishing are therefore transported correctly here.

## 4. Exact iterate and energy degrees in every order

Write z_j=G^j(z), z0=z. From(3) induction gives

    deg z_j=D^j,
    (z_j)_top=c_j*v0*H^(D^(j-1)) for j>=1,
    c_1=1,       c_(j+1)=a*c_j^D,       c_j!=0.            (7)

Also A_top=v0*(grad H)^t. Euler gives

    (grad H(v0)) dot v0=D*a!=0,
    partial_(v0) H is not identically zero.               (8)

The full iterates of(5) have

    w_n=A(z_(n-1))^t ... A(z0)^t*w,
    g_h composed with F^n
       =w^t*A(z0)*A(z1)*...*A(zn)*zn.                    (9)

The order of these matrix factors is important. Their leading rank-one
products are nonzero by(8). Thus

    deg F^n=D^n,
    deg(g_h composed with F^n)=(D+1)*D^n,       n>=0.      (10)

Here is the precise degree count, including the scalar energy rather than
only the coordinate iterates. For j>=0 the z-degree of A(z_j) is
(D-1)D^j. Summing j=0,...,n and adding deg z_n gives

    sum_(j=0)^n (D-1)D^j+D^n=(D+1)D^n-1.

The independent w adds one. For n>=1 the top in(9) is a NONZERO scalar
multiple of

    (w dot v0)*(partial_(v0) H)*H^((D+1)D^(n-1)-1).

For n=0 it is D*(w dot v0)*H. This proves the noncancellation used in
(10), not just an upper bound. Similarly w_n has z-degree D^n-1, so its
total degree is D^n, agreeing with the base z_n. All n,M,N are covered by
this written argument; no iterate campaign or numerical growth estimate
is needed.

Define Delta_F f=f composed with F-f. The binomial identity gives

    Delta_F^m(g_h)=sum_(j=0)^m (-1)^(m-j)*binom(m,j)
                                      *(g_h composed with F^j).

The term j=m has uniquely greatest degree, hence for EVERY m>=0

    deg Delta_F^m(g_h)=(D+1)*D^m,
    Delta_F^m(g_h)!=0.                                    (11)

The same holds for Phi,g_psi by(6). In particular neither the gradient
nor its scalar energy is Pascal finite. Every original plane Keller
extension admits these genuine HC4 normal forms in arbitrarily high degree.

## 5. A nonterminating formal series can still have a polynomial sum

Because Phi=X+terms of order>=2, Delta_Phi raises the order of a formal
series with zero constant term. The geometric operator series

    sum_(m>=0) (-1)^m Delta_Phi^m(g_psi)

converges X-adically and equals g_psi composed with Phi^-1, namely the
formal Legendre transform of psi. This follows directly by inverting
the operator id+Delta_Phi; it is also the Pascal identity credited in
Section1. Formula(11) says the individual terms NEVER vanish. It does
not say the formal sum cannot be a polynomial: cancellation among the
infinite terms is permitted in the adic topology.

For example start from the IDENTITY plane map. For every M,N>=2, take

    G=(z1+z2^M, z2+(z1+z2^M)^N).

It has the explicit polynomial inverse in Section2. Consequently psi has
a polynomial gradient inverse and polynomial Legendre transform, and its
energy is integral over its gradient target. Yet (11) holds in EVERY
order. This is an actual whole-plane inverse control, not a local Laurent
model or a hypothetical noninverse realization.

Pascal composition differences must also be distinguished from iterated
Nambu DERIVATIVES. The latter terminate on the energy exactly when the
gradient is invertible, by the retained criteria and the primary source.
The operator Delta_Phi does not commute with the Nambu frame in general.
Its nontermination cannot be substituted for a nonintegral energy, a
counterexample, or failure of Nambu-tower termination.

## 6. The non-Ricci-flat regime can be imposed without changing the gap

The constructed potential cannot have zero Levi-Civita Ricci tensor.
Indeed MJ0's cotangent proof would then force the plane G to be an affine
shear of the form c+L(z+v*f(ell(z))), ell(v)=0. Since G(0)=0 and DG(0)=I,
subtract the constant and linear terms of f: G itself is a shear tangent
to identity. Its iterate degrees are bounded. This contradicts(7), since
D>=4. Hence the full Hessian metric has NONZERO Ricci tensor.

This does not weaken MJ0, which is a valid sufficient subclass theorem.
It shows that an arbitrary original trace-pairing problem can also be
encoded in the complementary non-Ricci-flat regime by a legitimate target
framing. No new metric-rank classification or connection-parameter queue
is proposed.

## 7. Verification and remaining goal

The cheap checker transcribes the two-shear inverse and first jet, the
normalizing congruence/conjugacy, leading rank-one matrix multiplication,
and small symbolic energy degree controls. These are not actual noninverse
examples, and finite controls do not certify the unbounded theorem; its
proof is (1)--(11). Prior registered mathematical checkers are not replayed.

The primary paper already states that its criteria do not settle JC2/HC4.
This note does not replace those criteria by a stronger automatic claim,
adopt their differential-Galois foundations, or claim their known
nontermination example as new. B,S,Lambda and the full original goal are
retained, not replaced by a Pascal-nilpotence objective.

FULL THEOREM TARGET: prove or break full HC4/JC2.
CURRENT LAST UNSOLVED IMPLICATION remains

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.                 OPEN

Full marked-energy integrality is also OPEN. Published generic-degree-four
PLANE exclusion remains; no intermediate plane is manufactured. No
independent replay, formal verification or external review is claimed.
Review and outreach remain DEFERRED; no missing artifact is rebuilt.
