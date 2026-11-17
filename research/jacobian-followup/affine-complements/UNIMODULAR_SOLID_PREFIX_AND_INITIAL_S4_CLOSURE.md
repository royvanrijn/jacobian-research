# Unimodular solid prefixes and complete initial simple-S4 closure

15 November 2026. PF2RULED88 closes ALL trajectories starting at
PF82's initial full-S4/simple-meridian state3 on a torus knot,
in BOTH source pencils and through EVERY finite genuine cable tail.
There is no prefix eta condition. It includes the previously open
intrinsic seed ranges and the capped solid-torus exception.

The new reusable tools are a whole-chart source determinant congruence
and a solid-prefix argument that never computes a recursive cap slope.
They have no numerical bound on source orders, cable parameters or
tail length. All source smoothness/exactness/whole-chart hypotheses,
and the actual primitive peripheral maps, remain essential.

## 1. Two whole-chart complementary source obstructions

Retain PF77/79's smooth exact complementary source with TWO entire
omitted lines and their FULL charts. Its R exterior is irreducible
and boundary-incompressible, and its group is NONABELIAN: its actual
terminal disk has two genuine cones, and its group injects through
the incompressible path gluings. Even the terminal(2,2) Klein case
has nonabelian group. Therefore a target SOLID TORUS R exterior,
whose group is Z, is impossible by PF77's stable exterior-group
isomorphism. This uses no Waldhausen theorem with a compressible
target boundary. For nonsolid targets use the full PF79 peripheral
recognition, with all its stated irreducibility/asphericity/essential-
surface hypotheses and explicit Klein treatment.

A source with a WHOLE orientable negative-Euler disk/two-cone
R exterior can only be the rho=2 collar joined to a SINGLE D node.
The other whole orientable possibility is the corner fusion disk
(q,q+1,2), with THREE genuine cones. Put p=p_D, sigma>=1,
p>sigma, gcd(p,sigma)=1. The actual single-node chart has

    kappa=p+1-sigma>=2,
    S_D=(sigma-1)*p,
    n_D=(sigma+1)*p-sigma+1,
    delta=3*n_D-p^2.                                  (1)

Every nontrivial sigma conjugate is detected by the leading term,
proving the contact formula even with lower noncharacteristic terms.
Modulo p-1 the determinant satisfies the GENERAL congruence

    delta=5 mod(p-1).                                 (2)

If delta=1, then p-1 divides4, so p=2,3 or5. Substituting in(1)
gives sigma=-4/3,-1/3 or2/3, respectively, none positive integral.
Thus no such whole two-cone source has determinant1. This is a
whole SINGLE-node obstruction, not a claim about arbitrary deep
sources of determinant1 or arbitrary homology spheres.

A second useful comparison concerns an equal target pair(r,r), r>=3,
with target determinant r and R-meridian/fibre distance1. Matching
forces p=kappa=r, hence sigma=1. Then n_D=2r and
delta=r*(6-r). Equality with r forces r=5. But the actual collar
source R-meridian distance is |p-3sigma|=2, not1. Therefore this
marked target cannot be a complementary source either.

## 2. All initial full-S4/simple torus representations

Let the initial genuine torus knot have full S4 image, meridian a
single transposition and preferred longitude identity. Its central
fibre image is identity. Coprime generator orders and full S4 force
generator types4-cycle and3-cycle: transposition/3-cycle generates
full S4 only when their product is a4-cycle; a double transposition
with a3-cycle lies in A4. Thus, interchanging exponents if needed,

    e=4a, o=3t odd, gcd(4a,3t)=1.

The lifted Seifert base is planar; meridian boundary orbits are(2,1,1).
Fill the TWO fixed meridians and drill only the simple R core.
Because the fibre permutation is identity, these are primitive base
section caps. The ACTUAL remaining exterior E has base disk

    (a,t,3t), omitting orders1.                         (3)

After the R-meridian fill, the closed cover W has epsilon=-4/(eo),
Lambda=lcm(a,t,3t)=eo/4. PF72 gives regular-fibre order1 and cyclic
H1 of order t (at most two cones divisible by any prime). The
preferred R longitude becomes that fibre, since eo is even. Hence

    H1(W)=Z/t, [R]=0.                                 (4)

PF87's cyclic-pattern transfer, equivalently PF82's actual simple
tail calculation, preserves zero R through all later odd-winding
cables. PF67's nonzero common-source mark therefore excludes ALL
initial-state3 trajectories in the COMMON pencil at every length.

## 3. Nonsolid initial exteriors in the complementary pencil

In (3) both a,t>1 give THREE genuine terminal cones. With no tail
this cannot be the whole fusion(q,q+1,2): t and3t are odd>=3, so
the required2 would have to be a, leaving a nonconsecutive pair.
With a tail, it remains a terminal three-cone disk, forbidden by
the source terminal budget.

If a=1,t>=3, the exterior is disk(t,3t). With NO tail it must be
the whole collar/single-D model(1). Assigning p=t,kappa=3t makes
sigma=1-2t<0. The other assignment p=3t,kappa=t gives
sigma=2t+1 and delta=t*(9t+12), contradicting the target delta=t.

With at least one tail the initial disk stays terminal. Actual
fibre determinants are nonzero: the first is proportional to
r*eo-k, and subsequent ones to r*r_old*k_old-k, which cannot vanish
by gcd(r,k)=1. Incompressible gluing, fibre uniqueness and PF79's
recognition keep the edges canonical. Even-k tails have an annulus
with TWO equal cones(r,r), impossible at the root or internally.
Thus all tails have odd k. They add no homology, so delta=t and R=0;
the final R-fibre distance2 forces rho=2, n_R=3. PF83 then gives
3 not dividing p_D. No source D WINDING can be divisible by3.
The terminal winding must therefore be t and the free-tail order3t.
Its actual normalized target coefficient is +/-inverse4 modulo3t:
the original odd-exponent torus cone has inverse(e)=inverse4, and
its degree1 lift has fibre degree1. Source free-tail coefficient is
+/-1 (PF85), requiring3t to divide4-1 or4+1. This is impossible
for t>=3. Both orientation signs were allowed.

If t=1,a>=2, (3) is disk(a,3), with a coprime3. With NO tail the
target delta1 contradicts Section1's whole two-cone obstruction.
With any tail, the same canonical-edge/two-equal-cone argument makes
all tails odd-k. Their H1 stays trivial. Final distance2 and zero R
force the collar source congruence delta=2 mod3 (PF83), contradicting
delta1. All nonsolid initial cases are now closed.

## 4. A unimodular solid-prefix lemma avoids recursive cap slopes

The remaining initial case is a=t=1, T(4,3): (3) is an ACTUAL solid
torus and (4) is trivial. More generally, start from ANY actual
simple-state prefix whose R exterior E is solid and whose R-meridian
fill has H1=0. In the true basis m=M^2,L, its primitive solid meridian
has the form

    eta=A*m+L, A an arbitrary INTEGER.                 (5)

Indeed its coefficient of L is the intersection with m and equals
the order of the filled H1, namely1; choose the common sign so it
is +1. No value, sign, bound or recurrence for A is assumed.

Every passing next state3 cable has odd winding r>=3. Retain its
two degree1 fixed pieces and all genuine branch fills. Its UNFILLED
active double cover is an annulus with:

    k odd:  one cone r, fibre H=k*m+2r*L;
    k even: two cones(r,r), fibre H=(k/2)*m+r*L.

Capping by (5) adds the cone of order

    k odd:  J=|2Ar-k|,  an ODD positive integer;
    k even: J=|Ar-k/2|, a positive integer coprime r.   (6)

Neither is zero: the first is odd; the second would make r divide k.
Thus there is no zero-distance/reducible cap in this lemma.

If k is odd, the new R exterior is disk(r,J), with its R-meridian
distance2. J=1 makes it solid again. PF87 adds no homology when r,k
are odd, so its normal fill still has H1=0; (5) applies again with
some new integer A, whose value is irrelevant. J>=3 is the FIRST
NONSOLID outcome: a negative-Euler disk with two genuine odd cones,
closed H1=0 and R=0. If the trajectory ends there, Section1 excludes
its whole two-cone determinant1. If further tails follow, edges remain
canonical by the same fibre determinant. Even-k continuations fail
the equal-two-cone budget; all-odd continuations preserve delta1 and
end at distance2, contradicting the zero-mark collar congruence.

If k is even, the new R exterior is disk(r,r,J), orders1 omitted,
with R-meridian distance1 and filled H1=Z/r,R=0. It is NEVER solid.
If J>1 it has three genuine cones, with two equal orders r>=3.
It cannot be the whole(q,q+1,2) fusion and cannot be a terminal
piece under later tails. If J=1 it is disk(r,r). Ending there is
excluded by Section1's equal-pair determinant/distance argument.

For later tails the two order-r cones persist terminally. The original
cable cone has coefficient +/-inverse(k) mod r. The unbranched double
cover has fibre degree1/base degree2, so its TWO cone lifts keep that
coefficient. Since J=1 in(6), k=+/-2 mod r; both coefficients are
therefore +/-inverse2. A source free-tail order r has coefficient
+/-1, forcing r to divide3 or1, hence r=3. But all subsequent
passing tails must have odd k; their H1 remains Z/r, with zero R
and final distance2. The source delta=2 mod3 contradicts delta3.

Every finite trajectory either stays solid throughout (already
impossible by its group Z), or has a FIRST nonsolid outcome covered
above. This proves the all-length solid-prefix exclusion without
computing a single recursively updated cap slope.

The local cap/determinant proof is unbounded in A,r,k and tail length.
It can be reused at any global cover degree when the ACTUAL local
restriction is one cyclic2 orbit and the other attached pieces perform
the stated fixed meridian fills. It is not a theorem about arbitrary
noncyclic patterns or nonunimodular cap A*m+bL with b!=1.

## 5. Evidence and the frontier

The [checker](verify_unimodular_initial_s4.py),
[receipt](unimodular_initial_s4_receipt.json) and
[manifest](sources/unimodular-initial-s4/manifest.json) retain full
initial torus representation/residue classification, integral CW
filled/unfilled markings, actual repeated tails, exact slope and
determinant identities and solid/nonsolid controls. Finite controls
support, rather than replace, the written unrestricted proof.
Complete retained Hatcher/Waldhausen/Neumann-Raymond primary audits
are reused; no fresh visual audit or independently reproved deep
rigidity is claimed. Related primary records were consulted, with
no additional theorem imported.

With PF85-87 this removes the INITIAL-state3 frontier entirely.
Remaining quartic first arrivals are incoming2(S3) in the still-open
source sectors and complementary incoming7(A4) of initial T(2,3)
origin outside prior closed seeds. The A4 cap has coefficient2 on L,
so the unimodular lemma does NOT close it. Full longer-stage/quartic/
higher-degree/degenerate goals remain OPEN. No finite Keller map,
formal verification, external review or novelty/priority is asserted.
