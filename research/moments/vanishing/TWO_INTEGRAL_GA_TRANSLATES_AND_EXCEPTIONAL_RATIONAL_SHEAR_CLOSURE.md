# Two integral Ga translates and constant-volume birational closure

17 November 2026. Written all-degree regularity and safety theorem.
Both FULL exceptional Mathieu--Zhao statements remain **OPEN**.
External review is deferred; independent replay and formal verification
are not claimed.

A rational Ga translation between regular phases has a sharp divisorial
pole bound. In particular a reduced Ga-leading divisor and Ga degree
at least three force that translation to extend to a GLOBAL volume-
preserving automorphism. Thus some comparisons on the finite Ga chart
really do globalize while retaining the original unit and every marker.

On both exceptional quotients, every REGULAR rational Ga translate of
a single nonzero-weight phase is either an actual regular shear of that
seed, or has only one strict torus sign after subtracting a constant.
Every eventually pure polynomial in it is safe for ALL original markers.
The allowed rational translation is arbitrary, not assumed homogeneous.
The extension in Section6 covers every constant-volume birational change
preserving this Ga fibration up to a regular base automorphism.
This closes a birational construction architecture in every degree.
Its safe phases belong to the already known automorphic/nullcone or
one-sided sectors; no coverage of the full remaining pure locus is claimed.

## 1. Two integral translates over a discrete valuation ring

**Lemma.** Let R be a discrete valuation ring with fraction field K and
residue characteristic zero. Let p(V) in R[V] have degree d>=1 and
leading coefficient a_d. Suppose b in K has valuation -r<0 and

    p(V+b) belongs to R[V].

Then

    ord(a_d) >= r ceil(d/2).                            (1)

**Proof.** Put F(T)=p(bT), and let -N be the minimum valuation of its
coefficients. If N<=0, its leading coefficient already gives
ord(a_d)>=dr, which implies (1). Suppose N>0. For a uniformizer t,
the nonzero reduction q(T) of t^N F(T) has degree at most d.

Since p(V) is integral, its coefficient of T^j in F(T) has valuation
at least -jr. Thus coefficients with j<N/r disappear after reduction,
and q(T) is divisible by T^ceil(N/r). Apply exactly the same argument
to F(T+1)=p(bT+b), using integrality of p(V+b). Its reduction is q(T+1),
so q(T) is also divisible by (T-1)^ceil(N/r). Consequently

    2 ceil(N/r)<=d,       N<=r floor(d/2).

But N>=dr-ord(a_d), from the leading coefficient of F. Combining these
inequalities proves ord(a_d)>=r ceil(d/2). QED.

The bound is sharp for EVERY d and r. Put k=floor(d/2), epsilon=d-2k,
b=t^-r, and take

    p(V)=t^(r epsilon) V^(k+epsilon)(t^r V-1)^k.

Then

    p(V+b)=V^k(t^r V+1)^(k+epsilon)

is also integral, and the leading coefficient has valuation exactly
r ceil(d/2). These DVR polynomials are controls for the sharp lemma,
not asserted to be exceptional phases or Mathieu counterexamples.

## 2. The global pole bound on the ORIGINAL Ga torsor

Let X=SL2/Gamma, Gamma=2O or2I, A=C[X],

    D=z partial_x+w partial_y,       B=ker D=C[z,w]^Gamma.

Use the entire original torsor over U=Spec B minus its singular vertex,
as in [the higher-plinth model](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md).
For b in Frac(B), the finite sum

    sigma_b(f)=sum_(j>=0) b^j D^j f/j!

defines a rational Ga translation. In a rational slice it is V->V+b;
it is a function-field automorphism and fixes Frac(B). Suppose f and
g=sigma_b(f) are BOTH original regular functions. Let their Ga degree
be d>=1 and put

    ell=D^d f/d! =D^d g/d! in B.

For EVERY prime divisor E of U where b has pole order r_E, choose a
regular local Ga slice near its generic point. The torsor is locally
trivial there, so f=p(V) and g=p(V+b) have coefficients in the DVR at E,
with the SAME leading coefficient ell. Lemma1 gives

    ceil(d/2) r_E <= ord_E(ell).                        (2)

Thus poles of b can occur only on the ACTUAL Ga-leading divisor. In
particular b belongs to B[ell^-1]: normality extends a rational function
with no other codimension-one poles across the missing vertex of
codimension two. This also follows by pulling back to the punctured
plane and then applying finite-group invariance.

If every multiplicity of ell on U is smaller than ceil(d/2), there are
NO poles, so b belongs to B. This includes reduced ell with d>=3.
Then bD is a regular locally nilpotent derivation and sigma_b is an
actual global automorphism, with inverse sigma_(-b).

No component of ell=0 is dropped in this argument. It controls poles
of the proposed change of coordinates; it does not assert support
localization or vanishing of the corrected Neron component. The sharp
DVR examples explain why arbitrary multiplicities cannot be ignored.

## 3. When the canonical finite Ga model gives a genuine comparison

For f and g of the same d and ell, use the canonical slices from the
[finite Ga trace calculation](CANONICAL_FINITE_GA_TRACE_AND_LOSS_OF_THE_ORIGINAL_UNIT.md):

    v_f=D^(d-1)f/(d!ell),       v_g=D^(d-1)g/(d!ell),
    f=P_f(v_f),                g=P_g(v_g).

If P_f=P_g as polynomials over B_ell, put b=v_g-v_f. It is in Frac(B)
because D b=0, and g=sigma_b(f). Under the multiplicity condition of
Section2, this rational comparison therefore extends across ell=0 to
the GLOBAL automorphism sigma_b. This is a sufficient globalization
lemma, not a claim that all canonical finite Ga models coincide.

It is unit-preserving. Indeed D preserves the original volume Omega,
Db=0, and Lie_(bD) Omega=0. Equivalently for every original regular R,

    h(sigma_b R)=h(R),

because each nonconstant Taylor term b^j D^j R is
D(b^j D^(j-1)R), whose Haar mean is zero. Hence for ONE FIXED original
marker Q and every m,

    h(Q g^m)=h((sigma_(-b)Q) f^m).                     (3)

The inverse marker is regular and independent of m. The automorphism
is also an isomorphism of the full phase complements over the SAME
phase parameter, preserving Omega. It therefore transports the full
source, its canonical quotient and corrected obstruction before any
purity assumption. All phase values, boundary components and branches
are transported; none is removed by a scalar trace or norm.

This supplies a genuine unit-retaining bridge for these rational Ga
comparisons. The finite-chart restriction by itself still kills the
ambient unit, as previously proved. Here globalization, not contraction
of the discarded divisor, restores the full original comparison.

## 4. Positive-weight seeds force EVERY rational translation to be regular

Let f be a nonzero original phase of ONE left-torus weight w>0. Let n
be its highest Peter--Weyl degree. The lowering ladder gives

    d=(n+w)/2,       ell=D^d f/d! in B_n, ell!=0.       (4)

Only that highest isotypic degree contributes to ell. Lower isotypic
degrees have smaller Ga degree, even when there are right multiplicities.
In particular ell is a homogeneous ORIGINAL invariant binary form.

Every projective root orbit has size at least six on2O and twelve on2I.
Thus every root multiplicity M satisfies respectively

    M<=n/6 or M<=n/12,       hence M<ceil(d/2).         (5)

Here d>n/2, so ceil(d/2)>n/4. The finite quotient is etale on the
punctured plane, so these are also the multiplicities at divisors of U.

If sigma_b(f) is regular for ANY rational b in Frac(B), (2) and (5)
force b in B. No homogeneity of b was assumed. Its image is therefore
an actual regular volume-preserving shear of f and is pure SAFE against
every original marker by (3). Its safe polynomial recompositions are
exactly H(sigma_b f) with H(0)=0.

The same argument works for an invariant binary model with minimum root
orbit at least four; only the two requested quotients are needed here.
No claim about arbitrary finite groups with smaller root orbits is made.

## 5. Negative-weight seeds: regularity forces a one-sided image up to a constant

Let instead f have ONE weight w<0. If d=0, the translation fixes f,
so there is nothing to prove. Suppose d>=1 and g=sigma_b(f) is regular.
As above ell is homogeneous of binary degree n=2d-w. Section2 gives
b in B_ell, which has a finite homogeneous expansion

    b=sum_e b_e,       e is binary degree.

For b=0 the conclusion is immediate. Otherwise let e_min be its least
degree. A homogeneous b_e has LEFT weight -e, while D lowers weight by2.
If e_min+2>=0, every term in sigma_b(f) has weight at most w<0, so g is
strictly one-sided.

If e_min+2<0, the largest possible weight of g is attained uniquely at
the top Taylor term with all factors b_(e_min). Its nonzero coefficient is

    b_(e_min)^d ell,
    W_max=w-d(e_min+2)=-(n+d e_min).                  (6)

Because g is regular, its torus-weight component (6) is regular. That
component is D-invariant, hence belongs to B. Its binary degree is
n+d e_min, which cannot be negative. Thus W_max<=0. If W_max=0, the
component is a nonzero CONSTANT c; every other component is strictly
negative. Otherwise put c=0. In all cases

    g=c+g_minus,       g_minus has only strict negative weights. (7)

Regularity of the weight components follows either from the rational
torus grading or by Vandermonde projection of finitely many torus
translates of the regular g. This step does not assume that the
individual rational Taylor terms are all regular.

Consequently h(g^m)=c^m for every m>=0. For every polynomial H,

    H(g) eventually pure iff H(c)=0;
    when H(c)=0, H(g) is safe for EVERY fixed original marker.   (8)

Indeed H(c+g_minus) has only negative weights once its constant term
is zero. If H(c)!=0 its unit moments equal H(c)^m and cannot vanish
eventually. No inverse rational translation is applied to the marker
in this case; its original torus-weight support supplies the cutoff.

## 6. All constant-volume birational changes preserving the Ga fibration

The result extends from translations to the entire constant-volume class
preserving this Ga fibration up to a REGULAR automorphism of its base.
The assumptions on the base and volume are essential.

### 6.1 Base-fixed maps are affine, and small multiplicities force linear part one

Use any rational slice V with DV=1. Since Omega=dV wedge omega_B,
a base-fixed rational self-map Phi with Phi^*Omega=lambda Omega,
lambda in C^*, satisfies partial_V(Phi^*V)=lambda. Hence

    Phi^*V=lambda V+b,       b in Frac(B).               (9)

This already makes it birational. In a REGULAR local slice u_i the
formula is lambda u_i+b_i. If f and Phi^*f are regular, apply Lemma1
to p(lambda V) and its translate by b_i/lambda. The same pole bound
(2) holds for b_i, since lambda is a constant unit. If all leading
multiplicities are below ceil(d/2), every b_i is regular, giving a
GLOBAL affine torsor automorphism, with a regular inverse.

Now lambda MUST equal one. Write c_ij=u_i-u_j on overlaps. Because
the map fixes the base,

    b_i-b_j=(1-lambda)c_ij.

If lambda!=1, the local functions u_i-b_i/(1-lambda) agree. They
would give a global REGULAR slice v with Dv=1. But then

    Omega=d(v i_D Omega)

on the entire original X, contradicting its nonzero Haar cohomology
class. Therefore lambda=1 and the global map is a regular Ga shear.
This uses the ORIGINAL volume line, not a scalar rescaling of a
unit-killing comparison.

For a positive-weight seed, the multiplicity estimate (5) applies, so
EVERY regular endpoint of (9) forces this globalization and lambda=1.

### 6.2 Negative-weight seeds also allow arbitrary constant fibre scaling

Choose the homogeneous rational slice V=a_F/F, with F a minimal
homogeneous invariant and a_F=(xF_z+yF_w)/deg F. It has left weight2
and is regular over F!=0. For a negative-weight f with d>=1,
regularity of f(lambda V+b) implies b is in B[(F ell)^-1]: on F!=0
the local pole bound permits poles only on ell=0. Thus b again has
a finite homogeneous expansion.

The substitution rho_lambda:V->lambda V commutes with the left torus.
In the Taylor expansion

    f(lambda V+b)=sum_j b^j rho_lambda(D^j f)/j!,

rho_lambda(D^j f) still has weight w-2j. Therefore the entire argument
of Section5 applies unchanged. If the least binary degree of b makes
weights increase, the unique top term is again b_min^d ell. Regularity
puts it in B and forces its weight at most zero. Otherwise all weights
are already negative. The image is c plus a strict negative-weight
phase, and (8) follows for every H. The case b=0 is itself homogeneous
of weight w<0 whenever regular; d=0 is fixed by every base-fixed map.
No inverse of the rational map is applied to an original marker.

### 6.3 Every regular automorphism of the punctured Klein base lifts linearly

For completeness, every regular automorphism psi of U has a linear
lift T in GL2(C) normalizing Gamma. Here is an argument for the two
irreducible exceptional groups that does not use a classification of
acyclic curves or a plane degree-reduction theorem.

The finite etale cover C^2 minus0 -> U is universal: its domain retracts
onto the simply connected sphere S^3. Pull it back along psi and then
along that same cover. Each connected component over C^2 minus0 has
degree one, since finite etale complex maps are finite topological
coverings. It is therefore an algebraic isomorphism. A component gives
an algebraic lift of psi. Lifting the inverse and adjusting by a deck
transformation gives an inverse lift. Regular functions and the inverse
extend across the origin by Hartogs for the normal affine plane, so
the lift is a polynomial automorphism fixing0 and normalizing Gamma.

Let N be its maximum component degree and L_N its leading homogeneous
pair. Equivariance makes the linear span of the image of L_N stable
under the irreducible output Gamma representation. It must be all C^2.
If N>1, the constant nonzero Jacobian of the plane automorphism forces
Jac(L_N)=0. Two homogeneous forms of the SAME degree with zero Jacobian
are proportional: Euler's identity gives d(q/p)=0 for their quotient.
The same conclusion holds if one component has smaller degree, since
then the leading pair has only one nonzero component. Either case
makes the image span a line, a contradiction. Thus N=1 and the lift
is linear.

Compare the known normalizer description in
[Arzhantsev--Zaidenberg, MPIM11-61, Theorem5.3](https://archive.mpim-bonn.mpg.de/id/eprint/538/1/preprint_2011_61.pdf).
The argument above proves the needed punctured-base statement directly;
their acyclic-curve classification is not imported.

### 6.4 Lift the base automorphism without changing volume or seed weight

For a row-coordinate lift T, put delta=det T. The explicit map

    L_T(M)=diag(delta^-1,1) M T                         (10)

is a GLOBAL automorphism of SL2/Gamma covering psi. It descends because
T normalizes Gamma. It preserves Omega: choosing a square root of delta
writes (10) as left multiplication and right multiplication by elements
of SL2. It also commutes with the left torus, so it preserves every
seed weight. Its inverse and all transformed markers are regular.

Suppose now Phi covers ANY regular base automorphism psi and has
constant nonzero volume Jacobian. Factor Phi=L_T composed Phi0.
Then Phi0 fixes the base and has constant nonzero Jacobian, and

    Phi^*f=Phi0^*(L_T^*f).

If f is a single-weight seed, L_T^*f has the same weight. Sections6.1
and6.2 prove the safety classification for EVERY regular image Phi^*f.
For an arbitrary mixed seed satisfying the small-multiplicity condition,
Ga degree and leading multiplicities are preserved by L_T, so Section6.1
globalizes the comparison as well. In that case Phi is a global
Omega-preserving automorphism, transporting the full source and every
fixed original marker as in (3), before purity.

### 6.5 An actual Long counterexample when the exceptional root bound is absent

The distinction is realized by the retained original SL2 counterexample,
not just by abstract DVR polynomials. Use Long's matrix convention

    M=[[a,c],[b,d]],       ad-bc=1,
    D=b partial_a+d partial_c.

The seed f=cad has the single positive left weight1 and is pure SAFE.
The rational Ga translation with parameter h=1/d acts by

    a -> a+b/d,       c -> c+1,       b,d fixed.

It has volume Jacobian ONE on its regular locus, and its regular phase
image is exactly

    sigma_h(f)=(1+c)(ad+b)=g.

Here D-degree(f)=2, ell=D^2f/2=b d^2. The shift has a pole of order1
along d=0, where ell has multiplicity2. It satisfies the general pole
inequality but violates the small-multiplicity globalization hypothesis.
The trivial quotient has no exceptional root-orbit estimate to rule it
out. The resulting g is Long's actual pure unsafe phase.

There is an even stronger marker calibration. The fixed original bad
marker Q=-c has REGULAR inverse pullback sigma_(-h)Q=1-c. Nevertheless,
the [retained exact Long identities](COHOMOLOGICAL_PULLBACK_AND_A_UNIFORM_ORIGINAL_MARKER.md)
give for EVERY m>=1

    h((-c)g^m)=(-1)^(m-1)/(m+1),
    h((1-c)f^m)=0.

The second equality also follows immediately from positive weights.
Thus constant rational Jacobian, a regular phase image AND a regular
inverse marker do not establish (3). The rational map loses the Haar
line on its pole-free chart; a global Ga slice there makes Omega exact.
On the two exceptional quotients, Sections4 and6 exclude this entire
single-weight construction by forcing the map global whenever it could
produce a mixed regular phase. This is not an exceptional counterexample.

## 7. Full scope and the remaining exceptional obligation

Combining Sections4--6 proves, for EVERY regular image g of a nonzero
single-weight seed under a constant-volume birational change preserving
the Ga fibration up to a regular base automorphism, the existence of c with

    H(g) eventually pure iff H(c)=0 iff H(g) all-marker safe.

In the positive-seed case c=0 and the entire comparison is global. In
the negative-seed case g-c is already one-sided. These conclusions cover
arbitrary rational affine fibre terms and every polynomial H, in ALL degrees. They
remain valid after transport by any regular Haar-preserving automorphism.

The new input is the sharp regularity/globalization mechanism. The
resulting pure phases lie in the previously closed automorphic-nullcone
or one-sided classes. Thus this does not assert new mixed safe points
outside those classes or exhaust the remaining mixed pure locus. It
does rule out manufacturing a pure bad phase by a rational Ga shift of
an eigenphase while relying on regularity only at the final endpoint.

More generally (2)--(3) can be used for any proposed Ga normal form:
with small enough leading-divisor multiplicities, a regular rational
identification automatically transports the ORIGINAL full extension.
A birational base change that is not regular, a nonconstant volume
Jacobian, a ramified cover or general mixed seed weights remain outside
the architecture-closure theorem. Those cannot be treated as
unit-preserving comparisons without an additional argument.

The full target remains corrected coupling zero on the ENTIRE pure
locus, or an ACTUAL pure exceptional phase with ONE FIXED original bad
marker. The capacity, automorphism, primitive-fibre, relative Ga and
full family-transfer inputs are preserved. No finite-degree campaign,
extra decidability result or external review is introduced. Both FULL
exceptional Mathieu--Zhao cases remain OPEN.
