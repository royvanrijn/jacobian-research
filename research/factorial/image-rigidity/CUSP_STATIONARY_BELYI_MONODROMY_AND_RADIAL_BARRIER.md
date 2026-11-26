# The stationary Belyi cover at the ordinary Gamma cusp

This note concerns the complete algebraic stationary cover of the actual
positive-axis moment integral
\[
 I_m(w)=\int_0^\infty [S^a(w+S^d)]^m e^{-mS}\,dS.
\]
It proves its exact cusp normalization, geometric and arithmetic monodromy,
and a radial-action identity. It also gives an explicit obstruction to
treating every stationary solution as an active contribution. It does not
locate every actual moment zero or prove a new Strong Factorial case.

Throughout d>=2 is an integer,
\[
 b=(d-1)/2,\quad a=b^2,\quad D=a+d=(b+1)^2,\quad
 s_c=b(b+1)=(d^2-1)/4,\quad c=b/(b+1).
\tag{B1}
\]
The stationary geometry allows real a. Ordinary polynomial exponents on this
exact cusp require odd d, so a is an integer. The monodromy statements concern
rational maps over Q for every integer d>=2.

## 1. Exact identification with a classical rational Belyi family

**Theorem CuspStationaryBelyiCover.** On a local logarithm,
\[
 \Phi(S,w)=a\log S+\log(w+S^d)-S.
\]
Its nonzero, nonsingular stationary solutions satisfy
\[
 w=\Omega(S)=S^d\frac{D-S}{S-a}.
\]
Set
\[
 w_c=s_c^d/c,\qquad r=S/s_c,\qquad W=w/w_c.
\]
Then the stationary cover is exactly
\[
 W=B_d(r)=r^d\frac{1-cr}{r-c},\qquad
 B_d'(r)=-dc\,\frac{r^{d-1}(r-1)^2}{(r-c)^2}.
\tag{B2}
\]
It has degree d+1 and precisely three branch values 0,1,infinity. The
ramification partitions, in that order, are
\[
 (d,1),\qquad(3,1^{d-2}),\qquad(d,1).
\tag{B3}
\]

Indeed the stationary equation gives
\(w+S^d=dS^d/(S-a)\).
Substitute \(S=s_cr\), using
\(D-S=(b+1)((b+1)-br)\) and
\(S-a=b((b+1)r-b)\), to obtain B2. Differentiation gives the displayed
factorization. The zero fiber is r=0 of multiplicity d and r=1/c simple;
the pole fiber is infinity of order d and r=c simple. At r=1 the derivative
has order two and B_d(1)=1, hence the local degree is three. Every remaining
point over1 is simple, because there are no other derivative zeros.
These account for all branch points, including infinity. Numerator and
denominator are coprime since c is neither0 nor1.

The rational family itself is classical. Raimundas Vidunas,
*Belyi Maps from Zeroes of Hypergeometric Polynomials*, Mathematics13(2025),
156, Section5 Remark3, Eq98, printed15, gives
\[
 \phi(x)=
 \frac{(1-(p-q)x)^p(1-2px)^q}{(1-(p+q)x)^p}.
\]
Here the parameter q denotes the exponent called r in that source.
Taking p=1,q=d and x=(1-r)/2 gives B_d(r) exactly; the original source parameters are p=1,q_source=-1,r_source=d. No new Belyi-family
classification is claimed. Its identification with this actual stationary
map, and the direct consequences below, are proved in this note.

## 2. The cusp is the unique positive-parameter three-value collapse

Here briefly allow arbitrary real a>0, retaining integer d>=2, and put
b=(d-1)/2. The exact stationary derivative is
\[
 \Omega'(S)=-dS^{d-1}
 \frac{S^2-(2a+d-1)S+a(a+d)}{(S-a)^2}.
\tag{B13}
\]

**Theorem PositiveParameterBelyiCollapse.** This actual stationary cover
has exactly three branch values if and only if a=b^2. Otherwise it has
exactly four:0,infinity and two distinct finite nonzero values.

To prove completeness and distinctness, choose q with q^2=b^2-a.
The two nontrivial critical points and values are
\[
 S_\pm=a+b\pm q=(b\pm q)(b+1\mp q),\qquad
 w_\pm=(b\pm q)^{d-1}(b+1\mp q)^{d+1}.
\tag{B14}
\]
Neither point is0 or a: their product is a(a+d)>0, and the quadratic
in B13 equals a at S=a. Their values are finite nonzero, since a>0
also excludes S=a+d. At the cusp q0 they merge at s_c and w_c.

If0<a<b^2 choose0<q<b. Both values are positive, and, along this
q-parametrization,
\[
 \frac{d}{dq}\log(w_+/w_-)
 =\frac{4dq^2}{a(a+d)}>0.
\tag{B15}
\]
The log ratio starts at0, hence the two values are distinct.

If a>b^2 write q=i v,v>0. A continuous argument of w_+ is
\[
 \theta(v)=(d-1)\arctan(v/b)
                 -(d+1)\arctan(v/(b+1)).
\]
Its derivative is
\(-2dv^2/[a(a+d)]<0\), and its limits at0 and infinity are0 and-pi.
Thus -pi<theta<0 for every finite v>0; the two conjugate values
are nonreal and distinct. Formula B13 lists all finite critical
points, and the zero and pole fibers retain types(d,1), so there
are no further branch values. This proves both implications.
The a>0 hypothesis is essential; no boundary-a0 classification is asserted.

**Theorem OffCuspStationarySymmetricMonodromy.** For every fixed real a>0
different from b^2, the geometric monodromy is S_(d+1).

The cover is connected, and its branch cycle over0 is a d-cycle fixing
one hub. This forces primitivity directly. A block containing the hub
is invariant under that d-cycle, since its translate meets it at the
hub. If it contains any other point it contains all d other points;
otherwise it is a singleton, as are all its translates. Thus every
block system is trivial. A nontrivial critical point is simple, and
its branch value is separate by B15 and its complex counterpart, so
its branch cycle is a transposition. Form the graph of all conjugates
of that transposition. Its connected components form a block system.
A nontrivial edge rules out singleton components, hence primitivity
makes the graph connected. Transpositions along a connected graph
generate the full symmetric group, proving the assertion without a
classification theorem.

**Corollary OddGapAlternatingCuspEquivalence.** For every fixed odd integer
d>=3 and real a>0, the following are equivalent: the actual stationary
cover has exactly three branch values; a=(d-1)^2/4; and its geometric
monodromy is A_(d+1). This uses the cusp group calculation below and
the off-cusp symmetric theorem above. For even d the group is symmetric
in both regimes, so alternating monodromy is not a cusp test.

## 3. Geometric monodromy from the passport

**Theorem CuspStationaryGeometricMonodromy.** The geometric monodromy of B_d,
equivalently its generic Galois group over C(W), is
\[
 G_{\rm geom}=
 \begin{cases}
 A_{d+1},&d\ {\rm odd},\\
 S_{d+1},&d\ {\rm even}.
 \end{cases}
\tag{B4}
\]

Choose positively oriented standard branch loops, and the convention that
permutations act right to left, with
\(\sigma_0\sigma_1\sigma_\infty=1\).
A connected rational covering has a transitive monodromy action: after
deleting the finite preimages of the three branch values its source sphere
remains connected, and lifting paths connects the sheets. Equivalently,
the equation of the generic fiber is irreducible because it is linear in W
with coprime numerator and denominator in r.

Label the d-cycle over0 as \(\sigma_0=(1\,2\,\ldots\,d)\), fixing a hub h.
The generator \(\sigma_1\) is a3-cycle. Transitivity forces it to contain h;
otherwise both generators, and hence their product and inverses, fix h.
Write \(\sigma_1=(h\,i\,j)\). The product \(\sigma_0\sigma_1\) has cycle
type(d,1), the inverse type of \(\sigma_\infty\), so it has a fixed point.
Neither h nor j can be fixed, and any base label other than i or j is
moved by \(\sigma_0\). Therefore its fixed point is i and
\(i=\sigma_0(j)\).

Conjugating by powers of \(\sigma_0\) supplies hub3-cycles for every edge of
the base d-cycle. They generate every hub3-cycle: from the triples
\((h\,u\,v)\) and \((h\,v\,z)\), conjugation yields a triple on
\(\{h,u,z\}\); inverses fix its orientation if needed. Propagate along a
simple path in the connected edge graph. Products then supply
\[
 (h\,u\,v)(h\,v\,z)=(u\,v\,z).
\]
All3-cycles are obtained, and these generate the alternating group
(an even permutation is a product of pairs of transpositions, each
expressible as3-cycles). For d2 this argument already gives A3, while
\(\sigma_0\) adds a transposition. In general both generators are even
when d is odd, so the group equals A_(d+1). When d is even the d-cycle
is odd, giving the whole symmetric group. No Jordan theorem or
classification of primitive groups is used.

In particular a global fixed choice of a proper subset of saddle sheets
cannot be preserved by continuation on the full three-punctured base.
Local active subsets can still exist on cut domains; this statement does
not determine their original-cycle coefficients.

## 4. Exact discriminant, arithmetic monodromy and constants

**Theorem CuspStationaryArithmeticMonodromy.** The integer-cleared generic
fiber polynomial is
\[
 P_d(r,W)=(d-1)r^{d+1}-(d+1)r^d+(d+1)Wr-(d-1)W.
\]
Its discriminant in r is
\[
 \operatorname{Disc}_rP_d=
 (-1)^{d(d-1)/2}d^d(d-1)^{d-1}(d+1)^{d+1}
 W^{d-1}(W-1)^2.
\tag{B5}
\]
Let L be its splitting field over Q(W), and let its full constant field be
\(k=L\cap\overline{\mathbb Q}\). Then:
- if d is even, the arithmetic group is S_(d+1) and k=Q;
- if d is an odd square, both groups are A_(d+1) and k=Q;
- if d is odd and not a square, the arithmetic group is S_(d+1) and
\[
 k=\mathbb Q\!\left(\sqrt{(-1)^{(d-1)/2}d}\right).
\tag{B6}
\]

Here is a direct discriminant calculation. Put A=d-1,B=d+1,n=d+1.
At a generic root r_i,
\[
 P'_d(r_i)=dA\,r_i^{d-1}(r_i-1)^2/(r_i-c).
\]
Vieta and evaluation at1 and c give
\[
 \prod r_i=(-1)^dW,\quad
 \prod(r_i-1)=(-1)^{d+1}2(W-1)/A,\quad
 \prod(r_i-c)=(-1)^d4dc^d/(AB).
\]
Thus
\(\prod P'_d(r_i)=(-1)^d d^dB^{d+1}W^{d-1}(W-1)^2\).
The standard identity
\(\operatorname{Disc}P=(-1)^{n(n-1)/2}A^{n-2}\prod P'(r_i)\)
proves B5, with the stated sign.

For even d the geometric group is already symmetric, so the arithmetic
group is the same. Its geometric subgroup fixes every algebraic constant
in L, hence k=Q.

For odd d the arithmetic group lies between A_n and S_n. The square class
of B5 in Q(W) is \(\kappa=(-1)^{(d-1)/2}d\); every W-factor has even
exponent. The discriminant is a square exactly when the arithmetic group
acts evenly on the roots. Therefore the arithmetic group is A_n if
kappa is a rational square, and S_n otherwise. Since d>0 is odd,
kappa is a square exactly when d is an odd square (then d is1mod8).
In the nonsquare case the Vandermonde square root of B5 explicitly puts
sqrt(kappa) in L. The fixed field of the geometric A_n has degree two,
so it is precisely Q(sqrt(kappa))(W); every algebraic constant lies there.
The only elements of that rational-function field algebraic over Q are its
constant coefficients. This proves B6 and excludes extra constants.
In the square case the same index-one argument gives k=Q.
These are arithmetic generic-cover statements, not specializations at a
fixed numerical W or assertions about a specified real coefficient field.

## 5. Reciprocal geometry and a useful action ordering on the unit circle

The cover satisfies
\[
 B_d(1/r)=1/B_d(r),\qquad |B_d(r)|=1\quad(|r|=1).
\tag{B7}
\]
For r=e^(i theta), its continuous lifted phase has derivative
\[
 \frac{d}{d\theta}\arg B_d(e^{i\theta})
 =d-\frac{1-c^2}{1-2c\cos\theta+c^2}
 =\frac{2dc(1-\cos\theta)}
        {1-2c\cos\theta+c^2}.
\tag{B8}
\]
It is positive except at theta=0mod2pi, and has total winding d-1.
For W on the unit circle other than1, exactly d-1 distinct simple
preimages lie on the unit circle. The remaining two are off it and are
paired by r->1/conjugate(r). At W=1, all d+1 roots are on the unit circle,
with r1 triple and d-2 additional simple roots.

At a stationary solution the critical exponential factor is
\[
 \lambda(S)=dS^{a+d}e^{-S}/(S-a).
\tag{B9}
\]
For real a its absolute value is unambiguous even when a logarithm branch
is needed. Among unit-circle r=e^(i theta), put x=Re r. Apart from one
common constant,
\[
 \log|\lambda(s_cr)|=-s_cx-\tfrac12\log(1-2cx+c^2).
\]
Its derivative in x is
\[
 -s_c+\frac{c}{1-2cx+c^2}<0\quad(-1<=x<1),
\]
with equality only at1, since c/(1-c)^2=s_c. Thus two distinct unit-circle
saddles over the same W can have equal factor modulus only when they are
complex conjugates, and then W is real, hence W=+/-1. This is an action
ordering for that slice of the full saddle cover, not an all-plane
classification of active Stokes curves.

## 6. Exact radial criterion and a larger inactive saddle

**Lemma CuspRadialActionAndActivationBarrier.** On any cut domain with two
distinct simple nonzero stationary branches S_i(w),S_j(w), let D_ij be the
difference of their logarithmic critical actions. Then
\[
 w\,\partial_wD_{ij}=-(S_i-S_j)/d.
\tag{B10}
\]
Consequently the level curve Re(D_ij)=0 is transverse to radial dilation
at w precisely when Re(S_i-S_j) is nonzero. A nonzero imaginary action
difference is not enough for this radial test.

Indeed the envelope derivative is
\(w\Phi_w=(D-S)/d\), and root derivatives vanish by stationarity.
This proves B10. In particular conjugate saddles at a real w have
equal real action and equal real parts of S, so radial transport is
tangent to their equimodular locus. Such pairs exist on the positive
real axis away from the cusp, for example at d3,a1: there are exactly
two real roots of its quartic stationary equation and one nonreal
conjugate pair whenever w>0 and w!=16. Global radial transversality for
all formal saddle pairs is therefore false.

For odd d>=3, at W=1 there is a simple saddle S=-s_c in addition to the
physical quartic saddle S=s_c. Their moduli obey
\[
 \frac{|\lambda(-s_c)|}{\lambda(s_c)}
       =e^{2s_c}/d>1.
\tag{B11}
\]
The inequality follows already from exp(2s_c)>1+2s_c>d.
For d3,a1 the actual stationary polynomial at w16 factors as
\[
 S^4-4S^3+16S-16=(S-2)^3(S+2),
\]
while \(\lambda(2)=48e^{-2}\) and \(\lambda(-2)=-16e^2\).

Nevertheless the original positive-axis integral at positive w has a
unique positive stationary maximum. On S>0 its stationary map has
derivative
\(-dS^{d-1}(S-s_c)^2/(S-a)^2\), so at w_c its unique stationary point
in (a,D) is s_c. The logarithmic phase tends to minus infinity at0 and
infinity. Hence the elementary positive Laplace principle gives
\[
 \lim_{m\to\infty}I_m(w_c)^{1/m}=\lambda(s_c),
\]
not the larger factor modulus in B11. One can prove this without a
quartic expansion: a small interval around the unique maximum gives the
lower bound; for the upper bound, if psi(S)=S^a(w_c+S^d)e^-S and p>0
is fixed, then
\(\int\psi^m <=\|\psi\|_\infty^{m-p}\int\psi^p\).
Thus the negative saddle is not an active unique dominant packet of
the positive cycle. This is an explicit activation countercontrol,
not a numerical absence-of-zeros observation.

As an anchor on the negative real axis, put w=-t,t>0. There are exactly
two positive stationary points S_L in(0,a), S_R in(D,infinity), one
on each side of the zero t^(1/d). Their absolute critical-factor
difference H=log|lambda(S_R)|-log|lambda(S_L)| satisfies
\[
 \frac{dH}{d\log t}=-(S_R-S_L)/d<-1.
\tag{B12}
\]
As t goes from0 to infinity, H goes from+infinity to-infinity, by the
endpoint stationary asymptotics. Hence these two active positive-axis
packets have exactly one equimodular crossing and it is radially
transverse. Independently, for odd integer m the original I_m(-t)
is strictly decreasing in t, since its derivative is the negative
integral of a nonnegative even power; it has exactly one negative real
root. Even m has no real root. None of this excludes a common nonreal
root of adjacent orders.

### Distinguished physical branches have global radial transversality

**Theorem DistinguishedCuspSaddleSectorBranches.** Fix integer d>=3 at
the cusp B1. For every W in the upper half-plane there is exactly one
stationary preimage r_L in
\[
 {\cal D}_L=\{0<|r|<1,\ -\pi/d<\arg r<0\}
\]
and exactly one r_R in the same angular sector with |r|>1. They are
holomorphic simple inverse branches throughout that half-plane and
\[
 \Re r_L<1<\Re r_R.
\tag{B16}
\]
They continue the positive left and right stationary points from the
negative real W axis. Conjugation supplies the corresponding branches
in the lower half-plane.

Here are the boundary and root-count details. Write
\(K(r)=(1-cr)/(r-c)\). For Im r<0,
\[
 \Im K(r)=(c^2-1)\Im r/|r-c|^2>0.
\]
On the angular ray arg r=-pi/d, the factor r^d is negative real,
so B_d(r) lies strictly in the lower half-plane. On the unit-circle
arc with -pi/d<theta<0, B8 and its strict upper bound by d imply that
the continuous phase of B_d lies strictly between d*theta and0,
hence in(-pi,0). The positive-axis boundary maps to the real line,
with its pole mapping to infinity.

Consequently roots for W in a compact upper-half-plane set cannot
cross either domain boundary. Their counts, with multiplicity, are
constant throughout that connected half-plane. At W=iM for large M,
the simple pole gives
\(r_L=c+c^d(1-c^2)/W+O(W^{-2})\), inside the inner sector.
The d large roots satisfy \(r^d\sim-W/c\); exactly one has angle
\(-\pi/(2d)\), inside the outer sector. There are no other roots.
Thus each count is one. All roots are simple off the real branch
values, and uniqueness makes the local inverse branches global.

It remains to prove the strict right-hand inequality in B16; being
outside the unit circle alone would not suffice. On r=1+iy,y>0 put
\[
 \psi(y)=d\arctan y-\arctan(by)-\arctan((b+1)y).
\]
This is a continuous argument of B_d(1+iy), and
\[
 \psi'(y)=
 \frac{d s_c y^2[1+(s_c-1)y^2]}
 {(1+y^2)(1+b^2y^2)(1+(b+1)^2y^2)}>0.
\tag{B17}
\]
Here s_c>=2 since d>=3. If 1-iy lies in the lower angular sector,
then arctan y<pi/d, so0<psi(y)<d*arctan y<pi.
Therefore B_d(1-iy) lies in the lower half-plane, not the upper.
The right inverse branch cannot cross Re r=1. It has Re r>1 at
large iM, so that inequality holds everywhere in the upper half-plane.
The left inequality follows from |r_L|<1. This proves B16.
On negative W the only boundary preimages in these domain closures are
r_L in(0,c) and r_R in(1/c,infinity), which identifies their physical
positive-axis origin.

**Theorem DistinguishedCuspRadialEquimodularGraph.** Define the
distinguished critical-action difference on either open half-plane by
\[
 H(W)=\log|\lambda(s_c r_R(W))|
                    -\log|\lambda(s_c r_L(W))|.
\]
For every fixed0<theta<pi, H(rho*e^(i theta)) decreases strictly from
+infinity to-infinity as rho runs from0 to infinity. It has exactly
one zero on every ray, defining a smooth radial graph. Its upper and
lower graphs meet at the unique negative-axis crossing in B12, and
both end at the positive cusp W1. This statement classifies that
distinguished formal saddle pair, not the active original cycle or
the complete zero set of I_m.

The branches never hit a zero or pole of lambda in the open
half-plane. Their logarithms may therefore be chosen holomorphically.
B10 and B16 give the exact strict derivative
\[
 \frac{dH(\rho e^{i\theta})}{d\log\rho}
       =-\Re(S_R-S_L)/d<0.
\tag{B18}
\]
At rho0 the inner root tends to0 while the outer root tends to1/c,
so its two factor moduli tend respectively to0 and a positive
constant; H tends to+infinity. At rho infinity the inner root tends
to c, with |lambda_L| a positive constant times rho. The outer root
is asymptotic to(-W/c)^(1/d), with angle(theta-pi)/d; its real part
tends to+infinity. Thus its exponential e^(-s_c r_R) dominates every
power, |lambda_R| tends to0 and H tends to-infinity. Existence,
uniqueness and smoothness now follow from B18 and the real implicit
function theorem. These asymptotics are on each fixed ray; no
unjustified uniform estimate as the angle approaches an endpoint
is used.

For endpoint identification, each branch extends to every positive
real W other than1, with the same strict real-part inequalities:
a point on the potential crossing line Re r1 with Im r<0 maps strictly
below the real axis by B17, and the inner unit arc also maps strictly
below. At W1 every preimage is on the unit circle by Section5; the
only such preimage in either angular-domain closure is r1.
Thus both distinguished branches tend to1 and H tends to0.
On the positive real axis away from1 the radial derivative is
strictly negative, so H is positive on(0,1) and negative on(1,infinity).
For any rho0<1<rho1 these strict signs persist for sufficiently small
positive theta. The unique radial zero lies between them; letting
rho0,rho1 approach1 proves convergence to W1. On the negative real
axis the branches extend to the two simple positive saddles, and
B12 gives the single crossing with nonzero radial derivative.
The same argument and the implicit function theorem identify that
endpoint and join the conjugate graphs. This does not assert that
their radii exceed1, that their active coefficients are nonzero, or
that every actual polynomial zero approaches them.

## 7. Classical zero-locus theorem and the precise next gap

The actual primary Sokal2004 paper *Chromatic Roots are Dense in the Whole
Complex Plane*, Theorem1.5, printed226, assumes an EXACT finite sum
f_m=sum alpha_j(w) beta_j(w)^m with analytic fixed amplitudes and a
nondegenerate-dominance condition. It identifies limiting zero sets,
not a finite-m adjacent-order separation bound. A saddle expansion of
I_m with m-dependent prefactors and a remainder is not itself that
hypothesis. Its more general Theorem3.2, printed235, can detect limiting
zeros from a proved nonharmonic logarithmic growth limit; it still does
not supply active thimble coefficients or consecutive-order
transversality. These are classical methods, not imported conclusions
about the present moment polynomials.

The concrete remaining global problem is to decompose the positive
original cycle on cut cells of this explicit three-branch-value cover,
determine which saddle coefficients are nonzero, and classify
Re(D_ij)=0 only for co-dominant ACTIVE pairs. For the two distinguished physical continuations, B16--B18 already prove the global radial test and its single radial equimodular graph. Other formal pairs can fail transversality, as shown above. On any additional active regular cell, B10 reduces the test exactly to Re(S_i-S_j)!=0.
Any finite-order noncollision proof still needs uniform packet errors,
control near active multiple ties and branch values, and proof that
no actual root escapes those cells. The Belyi passport, full
monodromy and formal equimodular curves alone prove none of those
global SFC obligations.

## Quantitative torus-chain comparison (open)

[Long, *A Resolution-Free Proof of the Duistermaat–van der Kallen
Nonvanishing Theorem*, Theorems `thm:sublevel` and
`thm:log-integrability`, Proposition
`prop:primitive`](https://github.com/long-mathematics/duistermaat-van-der-kallen-nonvanishing/blob/20f13044f15fd445845bf33fe74e24c3a7424598/duistermaat_van_der_kallen_nonvanishing.tex)
proves an \(O(\varepsilon^{1/D})\) holomorphic-volume sublevel bound on
bounded semialgebraic families, integrability of the logarithmic form on
torus chains under Newton-polytope interiority, and an endpoint primitive
with a power-law bound. To test a transfer to the present active-cycle
and global-error problem, first exhibit one *actual* chain in the
required semialgebraic form with fixed complexity, a bounded relevant
Laurent polynomial on that chain, and stable interiority of its Newton
polytope. Then check that its form, dimension and degeneration match
the cited hypotheses. The source distinguishes smooth ODE transport of
homology classes from semialgebraic representatives; the former does not
supply the latter. It gives no numerical constants or first-nonzero-moment
bound, hence no explicit SFC cutoff. This is a targeted source reading,
not a verified transfer or a replay of its Lean development.

## Evidence and source scope

The checker reconstructs exact rational identities, literal resultants
and branch permutations; finite controls support rather than replace
the proofs for all d. The only imported canonical stationary geometry
is recorded in COALESCING_BINOMIAL_SADDLE_BRIDGE.md. The explicit
Belyi family is credited to the actual publisher Eq98, and the
monodromy and arithmetic arguments here are direct. All assurance
flags remain false; internal replay is separate from external review.

    python3 verify_cusp_stationary_belyi_monodromy.py --repo-root REPO --output RECEIPT
