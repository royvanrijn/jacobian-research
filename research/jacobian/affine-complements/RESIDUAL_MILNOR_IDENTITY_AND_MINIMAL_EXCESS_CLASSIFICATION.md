# Residual Milnor identity for every image branch

PF2RULED115, 16 November 2026. The following local theorem applies to an
ARBITRARY irreducible SINGULAR image branch, including several Puiseux pairs.
It needs no scalar-Morse, coprime-leading-form or contractibility condition.
Its global applications retain the ORIGINAL finite Keller normalization,
the WHOLE marked U=A2, and ALL omitted primes. The full goal stays ACTIVE.

Let h:(N,p)->(C2,0) be finite with just one reduced critical branch D.
Assume the source pair is a smooth-axis pair or a small cyclic-quotient-axis
pair, as supplied for actual omitted primes by the retained Borisov input.
Let n=1 in the smooth case, or let n be the cyclic uniformizing degree.
Write e>=2 for transverse index, q for LOCAL degree D->C^nu, ell for
surface rank, eta=ell-e*q, and c>=2 for the conductor of the chosen
irreducible singular image branch C. On the C2 uniformization, straighten
the invariant critical branch A={x=0} and put g=h composed with rho=(P,Q).
For a reduced local equation f of C, factor

    f(P,Q)=x^e*H,  x does not divide H.

Then H is a nonempty REDUCED residual curve and its ordinary plane Milnor
number satisfies the exact identity

    mu(H)=n*(eta-q)*(c-1).                         (1)

Consequently eta>=q. Equality eta=q holds precisely when H is smooth.
In that equality case C is the ORDINARY cusp, and the complete numerical
local possibilities are

| e | q | n | ell | eta | source |
|---|---|---|---|---|---|
| 2 | 1 | 1 | 3 | 1 | smooth |
| 3 | 2 | 1 | 8 | 2 | smooth |
| 3 | 1 | 2 | 4 | 1 | A1 quotient |

All three occur as genuine finite germs. This is a classification of the
equality-case numerical data and image type, not a claim of a unique
analytic equivalence class for every finite map germ with these data.
In particular, for EVERY e>=4 and EVERY singular image branch,

    eta>=q+1,  mu(H)>=n*(c-1)>0.                  (2)

No one-Puiseux-pair assumption is used in (1)-(2). A local unit excess
eta=1 is possible only for an ordinary cusp with either the smooth
index2/local-degree3 germ or the A1 index3/local-degree4 germ.

The PF114 linear mass inequality is not asserted for a general branch
with several Puiseux pairs; its separate torus-knot hypothesis remains.

## 1. Actual local degrees and the residual divisor

The small cyclic cover is free off its origin and on the punctured axis
has degree n. Thus g has local degree n*ell, its axis A has normalized
residue degree n*q, and

    Jac(P,Q)=epsilon*x^(e-1), epsilon(0)!=0.       (3)

This is the PF106 codifferent/holomorphic-Jacobian argument; no extra
isolated Jacobian zero is possible. Other omitted primes are disjoint,
so there is no other critical curve in this local factor. The other
components of g^-1(C) are GENERICALLY UNRAMIFIED, hence have multiplicity
one in H. They are ordinary preimage curves, not extra omitted primes.

Let H=B_1*...*B_k, ignoring a unit, and let q_j be the local normalized
degree B_j->C. Degree conservation at a general point of C gives

    sum_j q_j=n*eta.                              (4)

The retained PF106 conductor identity gives

    I(A,H)=n*q*(c-1)+1.                           (5)

Since c>=2, (5) is positive, so H is nonempty and meets A. We do not
assume H irreducible or smooth. The normalized excess eta used here is
ell-e*q, not the prime-index summand ell-e used before residue splitting.
Choose ONE local image branch C; replacing f by the entire reducible
target branch union would invalidate these formulas.

## 2. Polar orders on every residual branch

Take the first target coordinate generically, of order a on C^nu.
In characteristic zero the polar/conductor identity gives

    ord_Cnu f_Q=c+a-1.

For each B_j let c_j be its conductor and I_jl=I(B_j,B_l). For ANY
function P nonconstant on B_j, the corresponding source polar identity is

    ord_Bj Jac(P,H)=c_j+sum_(l!=j)I_jl+ord_Bj P-1. (6)

Here P need not be a generic LINEAR source coordinate. To prove this
extension, choose a generic source linear coordinate z on B_j. The
Dedekind formula gives the order of the complementary derivative of B_j
as c_j+ord_Bj z-1. Along a parametrization, the chain rule multiplies
this derivative by (dP/dt)/(dz/dt). This changes the order by
ord_Bj P-ord_Bj z. Multiplying by the other factors of H adds their
intersection numbers, proving (6).

Differentiate f(P,Q)=x^e H and wedge with dP. After dividing by
x^(e-1), (3) gives, up to an irrelevant sign,

    f_Q(P,Q)*epsilon=e*H*Jac(P,x)+x*Jac(P,H).

Restrict to B_j. Use (6), ord_Bj P=q_j*a, and the target polar order.
The result is

    q_j*(c-1)+1=I(B_j,A)+c_j+sum_(l!=j)I_jl.      (7)

Sum (7) over all k branches and use (4)-(5):

    n*eta*(c-1)+k
      =n*q*(c-1)+1+sum_j c_j+2*sum_(j<l)I_jl.

The characteristic-zero Milnor/conductor formula is

    mu(H)=sum_j c_j+2*sum_(j<l)I_jl-k+1.

This proves (1). Nonnegativity of mu, and mu=0 exactly for a smooth
reduced germ, prove eta>=q and the equality statement. Moreover the
residual Milnor number is divisible by n*(c-1). The equation measures
ACTUAL residual-curve complexity; it is not a free numerical assignment
of a different finite cover or a realization theorem.

The primary input is Garcia Barroso--Ploski, *On the Milnor Formula in
Arbitrary Characteristic* (2018), complete PDF2-4 and10, printed120-122
and128: Proposition2.1(ii),(iii),(v), its full proof, and Corollary3.4.
Their barred conductor-defined invariant equals the ordinary Milnor
number in characteristic zero; we explicitly use that specialization.
The all-branch summation (7), uniformizing factor n and residual identity
(1) are written here. No positive-characteristic version is asserted.

## 3. Smooth residual forces a one-Puiseux-pair image

Assume eta=q. On the smooth uniformization both A and H are smooth,
with intersection

    t=n*q*(c-1)+1>=2.

Their reduced union is analytically equivalent to x*(x-y^t)=0. Its link
exterior E_up is Seifert fibred (the weighted circle action has weights
(t,1)); its base is an annulus with one cone of order t. The restriction
of g makes E_up a connected UNBRANCHED finite cover of the knot exterior
of C. This uses the local-link construction, not a map after deleting
another global omitted prime while pretending the original map finite.

A finite cover of an algebraic-knot exterior with at least two essential
Puiseux pairs CANNOT be Seifert fibred. Here is the needed precise reason.
Use its standard iterated-torus decomposition, with winding p_i>=2,
coefficients alpha_i and

    alpha_i>p_i*p_(i-1)*alpha_(i-1).

The initial torus piece has disk base with coprime cone orders>=2; every
later cable piece has annulus base and cone order>=2. All these bases
have NEGATIVE orbifold Euler characteristic. Their finite covers retain
that sign, so are neither solid tori nor any product/Klein exception
to boundary-fibre uniqueness. The lifted cutting tori stay essential.

On a companion torus the two regular-fibre slopes, in its actual
meridian/Seifert-longitude basis, are

    (p_(i-1)*alpha_(i-1),1), (alpha_i,p_i).

Their determinant p_i*p_(i-1)*alpha_(i-1)-alpha_i is nonzero. Passing
to any finite torus cover applies an invertible rational linear map to
both lines, so cannot make them coincide. If E_up were Seifert fibred,
Hatcher Proposition1.11 would make the lifted essential tori vertical.
A closed horizontal torus cannot cover its base with nonempty boundary.
Restricting that putative fibration to the lifted pieces and using
boundary-fibre Lemma1.15 forces their boundary fibre lines to be the
inherited ones: all its exceptional products/Klein bundles have already
been excluded. Verticality would make those lines agree on each cutting
torus, contradicting the nonzero determinant. This proves directly that
the entire finite cover cannot be Seifert; no general descent theorem
from a Seifert finite cover is needed.

The imports here are the algebraic iterated-torus description in Lokteva
Proposition5, PDF9, the standard actual cable slopes (Ito Section2,
PDF3 uses the reversed p,q convention), and Hatcher Theorem1.9,
Propositions1.11/1.13 and Lemmas1.14-1.15, complete PDF15,17-18,20-21. These were
freshly rendered/read. The finite-cover slope argument is written above;
no blanket H1-based manifold recognition or unproved Seifert descent
is used. All exceptional zero-orbifold-Euler pieces are excluded by the
strict negative sign, not silently ignored.

It follows that C has one essential Puiseux pair. PF114 now applies to
the ORIGINAL local cyclic-axis pair. Since eta=q and q|delta|eta,
delta=q. Its inequality e*q<=5*eta-2*delta forces e<=3. This argument
only uses a torus knot AFTER deriving it, so does not assume the
conclusion in proving (1) or the equality case.

## 4. Finish the minimal-excess classification

The PF106 coefficient argument for cyclic-axis germs gives n<=e-1,
gcd(n,e)=1 and marked type1/n(1,-e), when n>1.

For e=2 it forces n=1. The Jacobian has order1, so the map has rank1.
PF105's all-index rank1 calculation gives q=1 and ell=m*(e-1)+1,
m>=2 for a singular image. Equality ell=(e+1)*q=3 forces m=2 and
image orders(2,3): the ordinary cusp.

For e=3,n=1, rank1 would give q=1 and odd ell=2m+1>=5, contradicting
ell=4. At rank0 PF105 gives ell=2m, m>=3, q<=2, and first image orders
(m,2m-2). Equality ell=4q forces q=2,m=4, leaving normalized orders
(2,3) and ell8. For n>1 the only possibility is n=2. The same retained
local theorem gives n*q<=2, so q=1,ell4, source A1 and image(2,3).
This proves exactly the displayed table. The PF105/PF106 LOCAL arguments
used here do not require a global Keller realization of these germs.

## 5. Global costs for arbitrary image arrangements

For an actual e>=4 omitted prime D, let its image normalization have k
nonimmersive finite parameters. Over them let M be the number of points
on D, counted once each. The sum of their local normalized degrees is
k*f. Applying (2) at ALL of those source points gives

    sum_(these points)eta_p >= k*f+M >= k*(f+1).   (8)

This applies even when the image has additional collisions or several
Puiseux pairs. Remaining local excesses stay nonnegative. It is a
singular-branch statement; distinct immersive branches crossing in the
target are not individually declared nonimmersive.

If an actual f>1 occurs, the retained Borisov power-map LOCATION theorem
puts its finite normalized ramification at a cyclic quotient point,
with q=f. PF106 already makes that chosen image branch singular, and
PF105/PF106 exclude f>1 at e in {1,2,3,4,6}. Thus e>=5, and (2) gives

    eta_at_quotient>=f+1,
    d>=(e+1)*f+r+1.                              (9)

For actual e=3, f=1. At any smooth source point a singular image branch
costs eta>=2; the only eta1 point could be its single A1 quotient point.
Thus k such nonimmersive parameters cost at least 2k-1, or 2k when the
prime has no quotient point. This does not exclude the allowed A1 cusp.

A useful CURRENT low-degree consequence already follows from PF114 and
is credited to it: every actual f>1 satisfies

    d>=5*f+ceil(7*f/5)+3.                         (10)

Indeed put w=e*f, e>=5. If 2w<=d, d>=10f, which is stronger than (10)
for f>=2. If 2w>d, apply PF114's same-image cusp cost and its different
RAMIFIED image cost: d>=e*f+ceil(f*(e+2)/5)+3, minimized at e=5.
In particular ALL omitted residue degrees are one through field degree15.
This is a residue statement, not inversion through degree15. No new
all-degree residue-one theorem is claimed, and the f>1 frontier stays
open with the actual singular residual H constrained by (1)-(2).

## 6. Exact controls, including a genuinely higher-characteristic branch

All equality profiles occur. The ordinary finite cubic cusp supplies
e2,q1,n1. On the smooth plane put

    P=X^2/2+X*Y^3-Y^6/6, Q=X*Y.

Then Jac=(X-Y^3)^2, ell8, q2, and

    (9*P^2-16*Q^3)=(X-Y^3)^3*(9*X-Y^3)/4.

Thus H is smooth, eta2=q. Quotient by -I gives the genuine finite A1
profile e3,q1,ell4,eta1. These are the retained PF105 finite maps; their
etale opens are C* times C, not the WHOLE Keller plane. The smooth q2
profile is allowed locally and cannot be promoted to an actual omitted
f2 prime on a smooth Keller normalization.

To guard against a hidden one-Puiseux-pair restriction, use the genuine
finite germ (the extra factor1+Y is a unit near zero)

    P=21X^2+42XY^3-7Y^6+42X^2Y+28XY^4-6Y^7, Q=XY,
    Jac=42(1+Y)(X-Y^3)^2.

On the critical axis it has parametrization
(Q,P)=(t^4,56t^6+64t^7), characteristic(4;6,7), q1,ell8,eta5.
Its chosen branch equation is

    f(U,V)=U^4-6272U^2V^3-917504UV^5
                   -16777216V^7+9834496V^6.

Exact elimination and division f(P,Q)/(X-Y^3)^3 produce the actual
reduced residual H. Local Jacobian-algebra calculations give c=mu(f)=16
and mu(H)=60, agreeing with (1): (5-1)*(16-1)=60. The higher-characteristic
germ is ALLOWED, not spuriously excluded or called an actual Keller map.

For the local algebra checks the verifier computes the lengths modulo
the gradient ideal plus consecutive powers of the maximal ideal.
Equal lengths at orders N,N+1 imply m^N M=m^(N+1)M in the localized
Jacobian algebra; Nakayama gives m^N M=0. Thus the stable length is the
EXACT local Milnor number, not an extrapolation from a jet cutoff.
The written universal proof, not these finite controls, proves (1)-(10).

The checker seals predecessors/primary pages, verifies the unbounded
summation and torus-cover determinant, checks all three genuine equality
profiles and the non-torus control, and keeps the genuine higher-index
finite quotient families. Classification/topology/polar imports are
declared. No historical priority, formal verification, external review,
automatic source-action descent, all-degree inversion or actual finite
target construction is claimed. The proportional PF114 bound remains;
the next broad gate is compatibility of the actual residual divisors and
global finite algebra. Published quartic retirement/shared priorities stay.

PF107's genuine finite A1 degree12/e9 control is also retained: on its
smooth uniformization, P^3-Q^2=-X^9*(9X^3+8Y^3)/64. The residual
ordinary triple point has mu4, exactly 2*(3-1)*(2-1). This actual finite
high-index cover remains ALLOWED locally and has the wrong Keller open.
