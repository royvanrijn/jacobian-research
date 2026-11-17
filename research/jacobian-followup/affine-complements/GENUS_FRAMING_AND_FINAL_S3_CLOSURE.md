# Genus, infinity framing and all final-stage S3 arrivals

15 November 2026. PF2RULED91 proves a degree-free canonical/genus
identity and closes EVERY incoming S3 first arrival which is the
FINAL cable stage, in BOTH source pencils and at arbitrary companion
depth. It also closes specified arbitrary later tails, stated precisely
below. A later prefix is never assumed to be its own affine packet.
The full quartic/longer-stage/higher-degree goal remains OPEN.

## 1. The all-depth canonical/genus bridge

Use ACTUAL positive iterated torus-knot data at infinity, with initial
T(p,q), p,q>=2 coprime, and genuine ordinary cables(r_j,k_j),
r_j,k_j>=2 coprime. Put B_0=pq and B_j=r_j*k_j, the actual regular
fibre coefficients in the preferred meridian/longitude bases.
Neumann--Norbury Section4 Figure1 and Section5.2 equation(7) give

    d_0=p+q,
    Delta_j=k_j-r_j*B_(j-1)<0,
    d_j=r_j*d_(j-1)+Delta_j.                           (1)

Here d=-K_X.Cbar-1 for the corresponding transverse infinity
curvette. The normal degree of the immersive rational component is
eta=d-1, exactly as in PF59/80. This uses the corrected ACTUAL keys,
not a literal binomial Wronskian or connectedness of a moduli space.
Complete primary PDF6/12/13 were re-rendered and visually read.

Let c_j=2g(J_j)-1, with g the Seifert genus of the actual knot.
Positive cable surfaces give

    c_0=pq-p-q,
    c_j=r_j*c_(j-1)+(r_j-1)*k_j.                       (2)

For completeness, the cable page consists of r_j copies of a
companion page and the connecting positive cable bands, giving
chi_new=r_j*chi_old-(r_j-1)*k_j. These pages are fibres: the pattern
map cyclically permutes the companion pages and glues along their
boundary collars. The Alexander polynomial of a fibre surface is
det(tI-monodromy) on its first homology, with nonzero extreme
coefficients and breadth twice its genus. The general Seifert-matrix
breadth bound makes these surfaces genus minimizing. This supplies
the elementary genus calculation without assuming the affine curve
itself is a generic smooth fibre or equating its geometric and knot
genera.

Subtracting(1) from B_j, and using(2), proves at EVERY depth

    B_j-d_j=c_j.

In particular, at a last(r,k) cable on companion J,

    eta=k-r*(2g(J)-1)-1.                               (3)

There is no field-degree or cable-length bound. The actual quartic
has eta<=0 at its FINAL stage by PF59; no such sign is assigned to
an arbitrary earlier prefix. The strict infinity inequality in(1)
is also independent of all source coefficients.

## 2. Solid S3 histories and the positive L-space genus threshold

PF90 already excludes every NONSOLID incoming S3 history in both
pencils at every length. Its exact necessary solid history is

    J_0=T(2,3), A_0=1,
    w_i odd>=3, k_i=2*A_(i-1)*w_i+epsilon_i,
    epsilon_i in {+1,-1}, A_i=A_(i-1)*w_i^2.            (4)

Let W=product w_i, so A=W^2 is an odd square. Define
D_i=2*A_i-c_i. Equations(2),(4) give

    D_0=1,
    D_i=w_i*D_(i-1)-(w_i-1)*epsilon_i>=1.              (5)

All D_i are odd. Equality1 holds exactly when ALL signs are plus.
A minus sign cannot be followed by another solid continuation:
after it, the ordinary fibre coefficient is B_i=2*A_i-w_i, whereas
the next solid k would be2*A_i*w_next+/-1, strictly greater than
w_next*B_i, contrary to(1). Thus a genuine infinity solid history
has only plus signs, with at most one minus sign at its END.
After such a final minus, D=2*w_last-1. The empty chain has D=1.

Every J in(4) is a POSITIVE L-space knot. The trefoil is one; and

    k_i/w_i=2*A_(i-1)+epsilon_i/w_i > c_(i-1)

by(5). Hom0912.4046v1 Theorem (the forward construction is due to
Hedden) inductively gives the assertion. Its Section2 rank formula
then gives the EXACT ordinary filling interval

    L(E(J))=[c,infinity], c=2g(J)-1=2A-D,              (6)

with c included and outside the INTERIOR. Negative finite slopes
are not L-spaces; zero filling has positive first Betti number.
These are NEW precise Floer imports, not independently reproved.
Ozsvath--Szabo math/0504404v1 Propositions9.5--9.6 were also read as
the primary rank/genus input behind this formula. Across both surgery
signs we use Hom's SIGNED-numerator formula, not the absolute-value
display in the OS v1 statement; the latter's proof itself gives the
signed kernel/cokernel computation. The manifest pins
both complete relevant HTML sources. Endpoint versus interior is
essential in the next argument.

## 3. Common pencil: all nonpositive arrival degrees

At an S3 first arrival put r=3s odd,k even. Let c=2A-D as above,
and suppose its ACTUAL arrival degree eta=k-r*c-1<=0. This is
automatic if it is the final stage. Define

    Jcap=Ar-k/2=(D*r-eta-1)/2 >= (r-1)/2.              (7)

The genuine solid cap in PF90 gives Y, disk with cones(s,Jcap),
attached to the ACTUAL ordinary E(J). Order1 cones are omitted.
The base fibre/section and the sign come from PF81's actual j0
presentation. Choose r*u-k*v=1, 1<=v<r. In Rasmussen's(h,-Y)
boundary basis the two unnormalized cone coefficients and the
ordinary meridian n have slopes

    a1=-3v/r,
    a2=(2A*v-u)/Jcap=2v/r-1/(r*Jcap),
    x(n)=v/r.                                        (8)

The actual ordinary longitude is L1=h^u*Y^(-k); thus t*n+L1 has
slope (u-t*v)/(k-t*r). The sum a1+a2 is nonzero, so the
section0 filling is a rational homology Seifert space and retained
Rasmussen Proposition3.9 applies to Y when it is nonsolid.

Assume first r>=9. If eta<0, its parity makes eta<=-2. Put
l0=(2A-1)*r-k=(D-1)*r-eta-1>0 and N=(l0*v-1)/r>=0.
In the exact floor/ceiling interval formula, denominator r gives
y_minus>=v/r. Denominator l0 gives

    ceil(l0*a1)=-3N, ceil(l0*a2)=2N+1,
    y_plus<=N/l0=x-1/(r*l0)<x<=y_minus.                (9)

Indeed 2*Jcap-l0=r, so the second fractional part is1/Jcap;
r>3 makes the first ceiling exact. The circular-interval convention
puts n OUTSIDE the L-space INTERIOR of Y. The ordinary meridian
is outside that of E(J), also by PF90's null-core lemma. The retained
published graph gluing theorem excludes their union as an L-space.

For eta=0 use a DIFFERENT matched slope: the genus endpoint
c*n+L1 of(6). Here k=r*c+1, choose v=r-1, and(8) becomes

    a1=-3+1/s, a2=2-D/Jcap,
    slope(c*n+L1)=1,
    0<D/Jcap<1/s.                                    (10)

For every l>=1, ceil(l*D/Jcap)<=floor(l/s)+1. The lower-endpoint
formula therefore gives y_minus<=1, with equality at l=1.
Thus y_minus=1, and this matched genus slope is outside the
INTERIORS of BOTH pieces. Again the actual union is non-L-space.
This avoids treating the meridian as an endpoint when it need not be.

All pieces in these r>=9 cases have negative Euler bases and
incompressible boundaries. The remaining fibre determinant is
r*B_companion-k>0 by(1). The union is an irreducible graph manifold,
with nontrivial LO group, contradicting the common source by PF74.
Its ACTUAL drilled R exterior is an annulus with cones(s,Jcap)
attached to E(J), also irreducible/incompressible graph. PF82's
non-L-space tail theorem excludes ANY later genuine simple tails
of these seeds. A prefix need not possess its own affine source.

Now r=3. The product collar BEFORE the solid cap gives the actual
closed end S3_t(J), with

    t=k/2-A,
    t<=c+(1-D)/2.                                    (11)

If eta<0 or D>1, then t<c. Formula(6) makes this non-L-space
unless t=0, which already has positive b1. In the former case it
is irreducible graph: for a nonempty solid chain, the last ordinary
fibre has B-c=w_last+1 or w_last-1>=2, so the last surgery cap
has distance |B-t|>=3; it leaves an incompressible disk two-cone
piece attached to the older companion. The initial trefoil case
has the same direct disk-cone check. Its drilled R exterior has
root annulus cone Jcap>=2 and ordinary E(J), and is incompressible.
The non-L-space or positive-b1 obstruction therefore propagates
through ANY later simple tails by PF82.

The remaining r=3,eta=0,D=1 case has t=c=2A-1 and may be an
L-space. Here the actual homology/linking calculation of PF73,
with the ACTUAL solid meridian A*m+L left symbolic, gives

    H1=Z/c, R=k*g=g, lambda(R,R)=-k^2/c=-1/c modZ.      (12)

For c=1 the mark is zero; for c>=3 its least nonnegative linking
value is1-1/c>1/2. Both contradict the common-source strict mark
gate PF67/71. This excludes the FINAL arrival, but this linking
argument is NOT asserted to persist under arbitrary later tails:
odd tail windings can change it by a square. Such tails remain open.
The retained actual cover control has solid prefix(3,7), A=9,
arrival(3,52), and tail(3,5). Its H1 is Z/17 and its R self-linking
becomes8/17. This passes that strict linking test alone; it is not
an affine packet or a full marked source match.

## 4. Complementary pencil at r>=9

Continue to assume the actual arrival eta<=0. For r>=9, both
orders(s,Jcap) in the drilled R root annulus are genuine. PF79's
whole source theorem forces a two-cone boundary annulus to have
CONSECUTIVE orders. From(7), Jcap-s>1 except

    s=3, Jcap=4, D=1, eta=0.                          (13)

At(13), the source corner merger must have rho=4,q=3,
sigma=3a,p=4a. Its remaining D path matches the WHOLE ordinary
E(J), consisting of its annulus cable pieces w_i and terminal
trefoil disk(2,3). The source terminal winding cannot be2 with
free order3: the last characteristic index sigma-1+3 is divisible
by2 and cannot break that last gcd. Equivalently the required
source coefficient inverse(3-1) modulo2 does not exist.
The terminal winding is therefore3 with free order2.
Products of the remaining characteristic gcd drops give EXACTLY

    a=3*W, A=W^2.                                    (14)

This uses the actual full canonical rooted path and cone multiplicities,
not only a terminal-order or unfilled-orbit comparison. The all-depth
contact bound of PF89 gives nD>=(sigma+1)*p-sigma+1, hence

    det Q=13*nD-9*p^2 >=12*a^2+13*a+13
                       =108*A+39*W+13.

The target has H1 order k/2-A*s=6*A-4, a contradiction. This
closes(13) at ANY number of older solid continuations.

If any later genuine simple tail follows a root with s,Jcap>=2,
that canonical annulus TWO-cone piece becomes internal. Source
internal D pieces have exactly ONE cone. The actual nonzero fibre
determinants preserve the two edges, so PF83 excludes every such
later complementary tail as well. This structural tail statement
even holds without eta<=0 whenever both orders are genuine.

## 5. Complementary pencil at r=3,Jcap>=2

Here the target root is an annulus with ONE cone Jcap, regular R
meridian distance1, attached to the whole E(J). Write Jcap=Jc
in this section. Coprimality gives3 not dividing Jc. Capping the
actual cone relation Z^Jc=h^(2A*v-u) in(8) gives

    3*(2A*v-u)=-1 mod Jc.                             (15)

Thus the target normalized cone coefficient is +/-inverse3 mod Jc.
The retained Neumann--Raymond slice convention and PF85/89 inverse
direction are used; replacing inverse3 by3 would lose information.

A simple rho=1 source has cone2 and distance3, so cannot match.
An UNMERGED higher-rho source has q=Jc and its unit-meridian cable
has normalized coefficient+/-1. Equation(15), allowing BOTH signs,
forces Jc dividing2 or4: Jc is2 or4.

Alternatively a rho=2 collar places the first D annulus at the root.
Write p=Jc*a,sigma=b*a,1<=b<Jc,gcd(b,Jc)=1. Its coefficient is
+/-inverse b and its R distance is |Jc-3b|. Matching(15) and
distance1 gives b=3 or Jc-3 and EXACTLY

    (Jc,b)=(4,1),(5,2),(8,3),(10,3).                  (16)

There is no parameter cutoff: these follow by solving two linear
absolute-value equations. Terminal Klein alternatives are excluded
by the retained full group/fibre argument PF79.

The terminal winding3 argument from Section4 holds again. For the
unmerged source, product of ALL D windings gives p=3W. For the
collar case the first D winding is Jc and the remaining product is
3W, so a=3W,p=3JcW,sigma=3bW. All fibre interfaces are canonical
by(1) and the actual source rigidity hypotheses are retained.

If the history is empty, A=1; k>=2 then gives Jc<=2, and Jc=2
has target b1=1, impossible for the source. Otherwise W>=3,A>=9,
and every candidate Jc in(15)--(16) has target determinant

    det_target=2A-Jc>0.                               (17)

For the unmerged Jc=2 source, nR=7,q=2,p=3W; determinant matching
requires7*nD=38A-2, hence A=3 mod7, impossible for a square.
For unmerged Jc=4, nR=21,q=4 gives21*nD=146A-4, forcing A=2 mod3,
again impossible. For collar Jc=4 or10, matching3*nD-p^2=2A-Jc
likewise requires A=2 mod3.
For the remaining collar Jc=5,b=2 or Jc=8,b=3, the SAME all-depth
contact lower bound gives

    det Q >= Jc*(3b-Jc)*a^2+3*(Jc-b)*a+3,

which is45A+27W+3 or72A+45W+3, respectively. Each strictly
exceeds(17). Thus EVERY candidate is excluded without finding any
new source coefficient or enumerating deeper clusters.

## 6. The last complementary product collar

For r=3,Jcap=1, equation(7) and eta<=0 force D=1,eta=0:
the solid history has ALL plus signs. The R exterior is E(J)
itself; its actual normal meridian is c*n+L1, by(11).
For the empty chain this is the already sealed PF78 trefoil
case with r=3,k=4 and full peripheral map.

Otherwise let w be its last odd winding. The boundary root of
E(J) is the ordinary annulus cable with cone w, coefficient+/-1
since its k=2*Aprev*w+1. Its fibre is B*n+L1 with B=2A+w,
so the R meridian has distance B-c=w+1.
An unmerged higher-rho source has distance1, and the rho=1
cone2 root cannot have odd order w. A rho=2 collar has coefficient
+/-inverse b, forcing b=1 or w-1, and distance |w-3b|.
Neither equals w+1 for odd w>=3: the second would force w=4.
This excludes every nonempty product-collar history.

Sections3--6 close ALL final-stage S3 arrivals in BOTH pencils.
They do not assert all later tails are closed.

## Evidence and exact remaining scope

The [checker](verify_genus_final_s3.py),
[receipt](genus_final_s3_receipt.json), and
[manifest](sources/genus-final-s3/manifest.json) retain symbolic
unbounded identities, complete divisor/residue deductions, exact
composed CW target controls, full primary HTML Floer imports and
the fresh complete NN PDF6/12/13 visual audit. Deep Floer, graph
gluing/LO, JSJ and full peripheral recognition are imported, not
independently reproved. Finite controls are not universal proofs.

Remaining S3 trajectories require AT LEAST ONE later simple tail.
COMMON: positive arrival eta, or r=3,eta=0,D=1 followed by tails,
outside prior PF81/83/90 closures. COMPLEMENTARY: r=3 at arrival,
or Jcap=1 (absolute value if arrival eta positive); otherwise a
two-cone internal annulus already excludes every later tail.
Keep all necessary infinity, sign, source and marked-peripheral
conditions. PF89's later A4 terminal pairs(2,4),(2,8),(4,10), with
their order/tail restrictions, remain separate and OPEN.
No full longer-stage/quartic/higher-degree/degenerate exclusion,
finite Keller map, formal verification, external review,
independent full reproof or novelty/priority claim is made.
