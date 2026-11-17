# Pole-four coordinates: odd companions excluded, even companions restricted

14 November 2026. PF2RULED9 proves two statements under the compatible
F1/contact-four hypothesis of [PF2RULED1](PROOF.md).

* No finite quartic map has coordinate poles (4,n) with n odd.
* After divisible-pole target shears, a surviving pole-four coordinate
  must have companion n=2r>=10, r odd. Its leading symbols are
  sigma^2,sigma^r, where sigma is a sextic section of L^2. It has a
  normal form P=S^2+H2 with S,H2 in V2, and sigma has no simple zero.
  Hence its root partition is (6),(4,2),(3,3) or (2,2,2).

This does **not** exclude the whole pole-four coordinate layer. The
minimum coordinate pole remains four, and the first open height-eight
pairs remain (6,8),(7,8). The arbitrary higher-pole and actual
compatible-embedding problems remain open.

Here L=O_B(B)=O_P1(3), Vm=H0(F1,O(mB)), and a function of exact pole m
has symbol in H0(L^m). We use the established full symbol restriction
and lifting sequence. Locally B=(z=0), w is a parameter on B, and

 dP wedge dQ = U*ell*z^-3 dz wedge dw,

where U is a unit and ell defines E. Away from contact ell is a unit;
at contact q it has the form unit*(z-phi(w)), ord(phi)=4.

## Weighted power criterion

For odd n, the leading differential equation gives symbols tau^4,tau^n
with a cubic tau. At a tau zero of multiplicity k<=3, give w,z positive
integer weights a,b, and let c=b-ak>0. Suppose

 P_in=t^(-4c)*H(Z),  w=t^a, z=t^b Z,
 H(Z)=(A+B Z+C Z^2+D Z^3)/Z^4, A!=0,

with nonconstant numerator. Write Q_in=t^-K*G(Z). The leading tau^n
term gives K>=cn, while the actual B-pole bound makes G a polynomial
in Z^-1 with pole at most n. If e is any lower bound for the weighted
value of ell, the initial wedge must vanish whenever

 (n+2)c+e > a(2k-1).                                      (1)

Indeed the prescribed volume has value >=a-2b+e>-(n+4)c,
whereas an initial wedge would have value -4c-K<=-(n+4)c.
The exact initial ODE is

 4c H G' - K H'G=0.

Comparison at Z=0 gives pole(G)=K/c<=n, hence K=cn and
G^4=constant*H^n. Since gcd(4,n)=1, every zero of H must have
multiplicity divisible by four. A nonconstant numerator of degree at
most three cannot satisfy this, including when roots coalesce.

The weights are local valuations, without any descended grading
assumption. Positive weights give finite initial sums in the formal
local ring. At a negative weighted value no regular z^j, j>=0 term
can occur, which justifies the denominator bound on G.

For the even case n=2r, r odd, the common leading section is instead
sigma in H0(L^2): symbols sigma^2,sigma^r. At a sigma zero of order k,
put c=2b-ak>0 and suppose P_in=t^(-2c)*H(Z), with H having pole four
at zero. The corresponding criterion is

 (r+1)c+e > a(k-1).                                      (2)

It forces K=rc and G^2=constant*H^r. Thus H must be a rational
**square**, which is a strictly weaker condition than a fourth power.

## Odd companions: complete initial atlas

Companion poles <=7 are excluded by [PF2RULED7](POLE_SEVEN_EXCLUSION.md),
which uses the published2022 degree alternative with its exception.
We therefore assume n>=9 odd. The following all-height step is intrinsic;
only this low-height base imports that degree theorem.

Lift tau to T in V1 and write P=T^4+A3, A3 in V3. Let a3 be its symbol.
If h=ord(a3)<2k at a tau zero, choose a=1,b=4k-h,c=3k-h.
The z^-3 term and tau^4/z^4 are initial, all lower terms are above them,
and H has a linear numerator. Inequality(1) holds at n=9 and increases
with n. Thus tau^2 divides a3 globally. Lift the quotient to A1 in V1:

 P=T^4+T^2*A1+H2,  H2 in V2.                             (3)

Write alpha for the cubic symbol of A1. If alpha is nonzero at a tau
zero, weights a=1,b=2k,c=k give a nonconstant numerator of degree at
most two and satisfy(1). Thus alpha vanishes on the support of tau.
If alpha is not proportional to tau, some multiple root has

 (k,h=ord(alpha))=(2,1),(3,1) or (3,2).

In the actual local Laurent expansion let j,l be the orders of the
z^-2,z^-1 coefficients of P. The z^-3 coefficient has order 2k+h;
regular terms in the lift of T cannot cancel it. Cover all possibilities:

| Condition | a,b,c | Why lower terms do not dominate |
| --- | --- | --- |
| j>=2h | 1,2k-h,k-h | l>=3h-2k is automatic on the three rows |
| j<2h | 2,4k-j,2k-j | 2l>=3j-4k is automatic on the three rows |

The first row retains a nonzero linear numerator coefficient; the
second retains a nonzero quadratic one. Ties can add a cubic term,
without making a fourth power. In each row(1) holds for n=9.
Therefore alpha is proportional to tau. Then A1=nu*T+constant.
Replacing T by T+nu/4 depresses the cubic term and absorbs the rest in
V2. We may write

 P=T^4+H2.                                               (4)

Let beta be the sextic symbol of H2. If it is not proportional to tau^2,
some tau zero has h=ord(beta)<2k. Choose a=2,b=4k-h,c=2k-h.
The numerator is quadratic, possibly with a cubic tie, unless the
z^-1 coefficient of order l satisfies 2l<3h-4k. The intermediate
z^-3 coefficient from the T lift lies strictly above the initial face.
Inequality(1) always holds for these weights at n=9. The complete
exception list where the z^-1 term dominates is

 (k,h,l)=(2,3,0),(3,5,0),(3,5,1).                         (5)

These require the actual fibre geometry, treated below. They are all
impossible. It follows that beta=mu*tau^2 globally, since both sections
have degree six. Subtracting mu*T^2 leaves a V1 function:

 P=T^4+mu*T^2+H1.                                       (6)

If its cubic symbol v is proportional to tau, P is a polynomial in T,
excluded by the [polar-centralizer lemma](LOW_COORDINATE_POLE.md).
Otherwise choose a tau root with h=ord(v)<k and weights

 a=3,b=4k-h,c=k-h.

The initial is (A+D Z^3)/Z^4 with A,D nonzero; the mu*T^2 term and
lift corrections are above it. Inequality(1) at n=9 holds except for
k=3,h=2. The last exceptional case is also treated below.

## Odd companions: exceptional fibre census

Put F=z^4(P-c0). On a normalized generic P fibre the volume identity
is, up to a nonzero unit and sign,

 dQ/dt = U*ell*z*w'(t)/F_z.                              (7)

Thus for w=t^a,z~t^b, ord(F_z)=s, the Q pole away from q is s-a-b.
A pole of degree zero would require a nonzero simple pole in dQ,
which is impossible for a derivative. A positive derivative order or
order zero makes Q regular at a missing point, contrary to finiteness.
Every generic P fibre is finite over the Q-line with total degree four;
all missing points must be poles and their pole orders sum to four,
summing over components if the fibre is reducible.

For(5) the relevant Newton models and normalized branches are:

| k,h,l | Leading terms of F | (ord(w),ord(z),ord(F_z)) | Q poles away from q |
| --- | --- | --- | --- |
| 2,3,0 | w^8+B w^3 z^2+C z^3 | (3,8,16) | 5 |
| 3,5,0 | w^12+B w^5 z^2+C z^3 | three branches (1,4,8) | 3,3,3 |
| 3,5,1 | w^12+B w^5 z^2+C w z^3+D(c0) z^4 | (3,11,25),(1,1,3) | 11,1 |

Here B,C are nonzero, and D(c0) is nonzero for generic c0. Coefficients
are units with the stated orders, not restricted numerical choices.
The omitted T-lift terms lie above all indicated Newton edges. All
edge roots used here are simple; their cusp parameters and derivative
orders persist for every permitted higher coefficient.

In the first row a single missing point already contributes five if
it is not q. If it is q, contact four has order eight along that branch
and Q becomes regular. Other roots of tau do not affect either local
contradiction. In the second row tau has a sole triple root. Away
from q the full pole total is nine. At q, ell has order at least four
on each branch, even allowing cancellation when ord(z)=4 ord(w);
Q is then regular. In the third row the whole infinity list again
lies over that triple root. Away from q its pole total is twelve.
At q, ell has orders eleven and one, making dQ have a nonzero simple
pole on each branch. This is impossible. Thus(5) is exhausted.

For the last case of(6), k=3,h=2, the Newton model is

 F=w^12+mu*w^6*z^2+B*w^2*z^3+D(c0)*z^4+(higher edge terms), B!=0.

Its two branches have (ord(w),ord(z),ord(F_z))=(3,10,26),(1,2,6).
Away from q they have Q poles thirteen and three, total sixteen.
There are no other infinity points because tau has a sole triple root.
Thus that root must be q. At q use the same local weights a=3,b=10,c=1.
The contact factor has weighted value at least min(10,12)=10, so

 (n+2)c+e >= n+12 > 15=a(2k-1),

already for n>=9. This repairs the previously failed inequality(1).
The initial numerator A+B Z^3 is not a fourth power, a contradiction.
The contact pole total alone would be 3+1=4 and would **not** exclude
this configuration; the strengthened differential inequality is needed.

All odd companion poles are therefore excluded.

## Even companions: completing a square

A companion pole divisible by four is lowered by a target shear,
repeatedly if necessary. Poles below four are already excluded.
Together with the odd-companion theorem and height-seven base, the
only remaining reduced case is n=2r>=10 with r odd>=5.
Write symbols sigma^2,sigma^r, sigma in H0(L^2), lift sigma to S in V2,
and put P=S^2+A3. At a sigma root k<=6, if a3 has order h<k, use

 a=1,b=2k-h,c=3k-2h.

The initial numerator is linear and not a square. All lower terms are
strictly above it, and inequality(2) holds at r=5 for all k<=6,h<k.
Thus sigma divides a3. Lift its quotient to A1 in V1; after
S -> S+A1/2 the exact algebraic square completion gives

 P=S^2+H2, S,H2 in V2.                                  (8)

This normal form is necessary, not sufficient for a finite map.
In particular S has pole two: its lower Laurent term must be retained
in any further weighted analysis. Treating it as a pole-one lift would
silently discard possible initial terms.

## Even companions: a simple sextic root is impossible

Suppose sigma has a simple zero. Use weights a=b=1,c=1 and write the
weighted initial of S as (A+B Z)/Z^2, A!=0. If eta is the value at
this point of the symbol of H2, then

 H(Z)=((A+B Z)^2+eta*Z^2)/Z^4.

Inequality(2) holds. Therefore H must be a rational square. The
numerator discriminant is -4*A^2*eta, so eta=0, including degenerate
quadratic and linear cases.

Now use actual regular local coordinates R=z^2*S,z. They are coordinates
because sigma has a simple zero. Let h(R,z)=z^2*H2, with h(0,0)=0.
The generic fibre is

 F=R^2+z^2*h(R,z)-c0*z^4.

Its R-critical point Rcrit(z)=O(z^2) is independent of c0. Formal
completion of the square gives, after a unit coordinate change,

 F=rho^2+f(z)-c0*z^4, ord(f)>=3.

For completeness, solve F_R=0 recursively: its R-linear coefficient
is a unit, and its value at R=0 is divisible by z^2. Taylor expansion
at Rcrit factors F-F(Rcrit,z) as (R-Rcrit)^2 times a unit. Absorbing
the formal square root of that unit gives rho and the displayed form.

For generic c0 the second term has order exactly d=3 or4. This assertion
retains all higher coefficients: if the cubic coefficient vanishes,
the free generic value c0 makes the quartic coefficient nonzero.
The volume remains a unit times ell*z^-3 dz wedge d(rho), and the fibre
formula in these coordinates is

 dQ/dt = unit*ell*z*z'(t)/F_rho.

For d=3 the single cusp has z=t^2,rho~t^3, and dQ has order at least
zero. Then Q is regular at a missing point, contradicting finiteness.
For d=4 there are two smooth branches z=t,rho~t^2. On each branch,
if ell is a unit dQ has a nonzero simple pole, impossible for a
derivative; if ell vanishes, Q is regular and finiteness fails.
This also covers q and all possible contact cancellations. A generic
P fibre cannot contain E, so the nonzero unit differential is retained.
Consequently sigma has no simple root.

The complete partitions of a degree-six divisor with all parts >=2 are

 (6), (4,2), (3,3), (2,2,2).

In all but(3,3), sigma is a nonzero constant times a cubic square.
In the remaining case sigma is a nonzero constant times L1^3*L2^3 for
distinct linear forms. Over C the constants can be absorbed. This is
a root restriction only; none of these four types is excluded here.

## Assurance and next gate

The [exact checker](verify_pole_four_companions.py) checks the universal
weighted ODE, complete finite weight/dominance lists, all exceptional
Newton edges and derivative orders, contact alternatives, square
completion and the sextic partition census. Its
[receipt](pole_four_companions_receipt.json) retains those controls.
The all-height initial coverage, formal local-coordinate argument,
generic-fibre normalization and finiteness interpretation are written
proofs, without formal verification, independent full reproof or external
review. The height<=7 base imports only the already pinned2022 theorem
through PF2RULED7; no new literature theorem or novelty priority is claimed.

The next pole-four gate is precisely the square-leading n=2 mod4 layer
with the four sextic partitions above. First height-eight pairs and
genuinely mixed-weight candidates on the marked Wright ring remain open,
as does actual ruled compatibility. No finite map has been constructed.

The subsequent [PF2RULED10 formal-jet theorem](POLE_FOUR_NONSQUARE_EXCLUSION.md)
excludes the entire(3,3) case. The subsequent
[PF2RULED11 leading-period theorem](BOUNDARY_LEADING_PERIOD.md) excludes
(2,2,2) and requires contact in the mixed(4,2) type. The four-partition
list and this note's receipt retain the necessary types established by
PF2RULED9. Current survivors are only(6),(4,2), with that contact
restriction; neither remaining type is excluded in full.

The adopted frontier(6,8),(7,8) incorporates the subsequent
[PF2RULED12 full(5,8) exclusion](POLE_FIVE_EIGHT_EXCLUSION.md); it is
not an additional exclusion proved by this earlier note.
