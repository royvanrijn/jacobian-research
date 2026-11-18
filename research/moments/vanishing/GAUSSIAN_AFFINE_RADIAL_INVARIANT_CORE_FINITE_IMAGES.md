# Affine radial invariant cores: scalar rigidity and finite actual images

Use one Witt Gaussian pair Z,W with E[ZW]=1 and an independent ordinary
Gaussian Y with E[Y^2]=1; all other basic pairings vanish. Over a
characteristic-zero field K put S=ZW and consider

    P=P0+R, P0=S H(Y)+G(Y), H,G in K[Y].             (A1)

Every monomial of R has strictly positive Z-minus-W angular weight, or
every monomial has strictly negative weight. The case R=0 is included.
There is no bound on H,G or on the W or Y degree of R.

For ORIGINAL fixed spatial markers define

    Phi_P[Q]=E[Q exp(tP)], V_P=Phi_P[K[Z,W,Y]],
    M_P=K[t]V_P.

The complete scalar and full-image criterion is

    Phi_P[1]=p(t)exp(psi(t)), p in K[t],p(0)=1,psi in tK[t]
      iff M_P is finite over K[t]
      iff H=0 and G=c+bY.                           (A2)

Every finite case has primitive vacuum

    p=1, psi=ct+b^2t^2/2, M_P=exp(psi)K[t].         (A3)

For nonconstant P one has literally V_P=M_P; for the constant phase c,
V_P=K exp(ct) while M_P=exp(ct)K[t]. Within the whole stated family,
purity and safety are equivalent. Both require zero core H=G=0, so
P=R. General mixed-parity invariant cores of higher S degree are not
classified here.

## 1. Exact Gamma elimination preserves the original endpoint

The core P0 has weight zero, so the actual vacuum ignores every strict
remainder factor:

    Phi_P[1]=Phi_P0[1].                              (A4)

The stronger full-image equality M_P=M_P0 is the
[angular-core theorem](GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md).
We also prove all-original-marker sufficiency directly below.

For the original radial variable E[S^n]=n!, put q=1-tH(Y). Exact formal
factorial contraction gives

    Phi_P0[S^m Y^j]=m! E_Y[Y^j exp(tG)/q^(m+1)].    (A5)

Unbalanced Z,W markers have zero core period. Formula (A5) is the full
original balanced-marker map, including m=0. In particular

    Phi_P0[1]=E_Y[exp(tG)/q].                       (A6)

The Gamma endpoint is encoded by q^(-1), not discarded. Its relative
Stein relation is (m+1)S^m-q S^(m+1), whose image in (A5) is zero.
Multiplication by q in (A6) recovers the endpoint Gaussian source:

    Phi_P0[q F(Y)]=E_Y[F(Y) exp(tG)].               (A7)

Here q is a permitted polynomial-time amplitude. Thus all powers
q^(-n), n>=1, and all polynomial Y amplitudes are realized by original
balanced markers with polynomial-time coefficients; q^0 is obtained by
(A7). After rational time extension they span the entire q-localized
Y source.

Because q has constant t coefficient one, q^(-n) has polynomial-Y
t-coefficients. The ACTUAL Gaussian functional therefore extends to
C(t)[Y,q^-1] by coefficientwise expansion, with a finite Laurent
principal part. Multiplication by exp(tG), Gaussian integration by parts
and time differentiation are locally finite in every t coefficient.
Its spatial twist and time connection are those of

    phi=-Y^2/2+tG(Y), eta=dY/q.                     (A8)

Indeed differentiating Y^j/q^(m+1) gives the term
(m+1)tH'Y^j/q^(m+2), exactly matching the original S H' term in the
Y Stein operator. Its time derivative gives
(m+1)H Y^j/q^(m+2), matching multiplication by S H in the original
connection. The remaining derivatives and G terms match directly.
Thus every spatial identity and every parameter-covariant iterate is
tracked on the FULL localized source, not only on one scalar value.
There is no leftover Gamma boundary after this exact elimination.

If Phi_P0[1]=p exp(psi), gauge phi to phi-psi and use eta/p. The
resulting genuine horizontal cyclic functional takes its seed to one.
The inverse of p is rational time and creates no spatial pole; it adds
deg p to every infinity order below.

## 2. Moving logarithmic sections and all refined component orders

Assume H nonconstant and put M=deg H>=1, r=1/t, v(r)=1, Q=H(Y)-r.
Exactly,

    eta=-r dY/Q.                                    (A9)

The polynomial Q is squarefree over C(r), since a root shared with H'
would make r equal a fixed critical value. After a finite algebraic base
cover its roots are moving sections beta(r). The volume has a simple
logarithmic pole on each generic section. At spatial infinity it is
logarithmic for M=1 and holomorphic for M>=2.

Mark ALL those sections, spatial infinity and the phase fibres
(phi-psi)=0,1,infinity. Resolve their collisions after finite base change
and retain a complete proper marked model. The following bounds apply
to EVERY generic component of every further refinement.

At an ordinary finite centre or nonconstant finite Y residue, H is a
unit and v(eta)>=1. At a root alpha of H of multiplicity m, write
N=v(Y-alpha)>0.

* If mN<1, then v(Q)=mN and
      v(eta)>=1+(1-m)N>=1/m>0.
* If mN>1, then v(Q)=1 and v(eta)>=N>0.
* If mN=1 without extra cancellation, v(eta)>=1/m.
* With extra cancellation at mN=1, v(Q)>1 and v(H')=(m-1)/m.
  The EXACT relative differential identity dQ=H' dY gives
      eta=-r dQ/(QH'),
      v(eta)>=1-v(H')=1/m>0,
  using v(dQ)>=v(Q).

At spatial infinity, v(Y)=-N<0 and H dominates r, so

    v(eta)>=1+(M-1)N>=1.                            (A10)

These cases include arbitrary multiplicity, cancellation and every
phase-dependent refinement. They do not depend on G or psi. They are
not restricted to a stationary-face list selected before the scalar
gauge. A polynomial p only improves their orders.

The moving sections themselves have the same positive normalization.
For a branch approaching an m-fold root of H, beta-alpha is a unit
times r^(1/m), and the logarithmic residue of eta is

    -r/H'(beta), of order 1/m.                       (A11)

After separating marked sections they meet smooth finite components,
away from nodes. Their local logarithmic coefficients consequently have
the displayed positive factor. Pole circles and approach segments are
controlled with these residues, including their collisions. If the
phase value at a puncture has nonzero polar action, that local endpoint
has no zero-exponential contribution; if its action is regular, (A11)
makes its seed value a positive power times at most logarithms.

## 3. Complete moving-puncture comparison and the unit contradiction

This is a moving-section argument with the new bounds in Section2.
The imported comparison is the complete proper marked-curve construction,
not a fixed-puncture statement applied outside its hypotheses.

On a finite/finite node, write the relative logarithmic volume as
f(u,v)du/u with uv=r after suitable ramification. Extract its two
positive vertical orders from Section2. The remaining coefficient is
bounded holomorphic because marked horizontal poles have been separated
from nodes. Its crossover annulus contributes a positive power times
1+|log r|. At a finite/pole node the positive finite-side factor survives
the substitution into the pole chart. Descending exponential integration
absorbs every fixed extra pole order. Rational phase and inverse branches
have only algebraic growth, dominated by the deep phase cutoff. At
phase-regular marked punctures the seed has only logarithmic poles and
the positive residue normalization (A11). Thus compact pieces, every
annulus, all tails and their gluing maps are retained.

The relevant construction and primary-source scope are recorded in the
retained [all-centres and crossover proof](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#all-centres-punctures-and-crossover-estimates),
[normalized moving-cover construction](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#normalized-cover-phase-and-horizontal-poles)
and [complete period comparison](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#infinity-and-complete-period-comparison).
The three underlying formal-nearby-cycle, finite/pole-node and ordinary
period-pairing inputs are identified and scoped in the
[primary manifest](../../../results/gaussian-boundary/sources/laurent_primary_manifest.json).
The quantitative bounded-degree path estimate is
[analytic_core.tex](../../../results/gaussian-boundary/analytic_core.tex),
Lemma `lem:curvelength`; the exact passage from a horizontal cyclic
unit quotient to an analytic period is
[seed_period_comparison.tex](../../../results/gaussian-boundary/seed_period_comparison.tex),
Lemma `lem:seedfunctional`. These are the retained complete-comparison
imports, with their existing review limitations.

All hypotheses of that seed passage are checked here. The total space
is the smooth affine open q!=0 in the t,Y plane after removing finitely
many time denominators. Both t and phi-psi are regular there and are
algebraically independent: the t^0 coefficient -Y^2/2 cannot be removed
by a scalar gauge. Section1 constructs the actual horizontal functional
on this exact de Rham source and all its covariant iterates. No inverse
of H(Y) is formally assumed at t=0, and no localization isomorphism at
the moving q-divisor is needed.

Hence every zero-exponential normalized eta/p period decays at time
infinity. The ordinary-parameter perfect period pairing would represent
the transported unit quotient by a full curve period with seed value
one, a contradiction. This excludes every nonconstant H for arbitrary G.

## 4. Constant H and the ordinary Gaussian scalar criterion

Suppose H=h!=0 constant. The volume is dY/(1-th), with scalar infinity
order one before division by p. The phase -Y^2/2+tG(Y)-psi has degree
D=max(2,deg G)>=2. Its top coefficient has order at most zero: it is
-1/2+tG_2 in degree two, and tG_D in higher degree. Therefore

    D v_infinity(1/[p(1-th)])>v_infinity(top).        (A12)

The retained absolute polynomial all-period leading theorem,
[Formal stationary balls](../../../results/gaussian-boundary/polynomial_stationary_balls.tex),
Theorem `thm:stationaryballs`, excludes its normalized unit quotient.
This theorem retains the irregular infinity pole of dY: its scalar
tail prefactor has order v(1/[p(1-th)])-v(top)/D>0. No logarithmic-only
criterion is used for this polynomial-volume case.

If H=0 and deg G>=2, the actual phase has top coefficient of order -1
and normalization 1/p of nonnegative order. The same strict criterion
excludes the scalar premise. Thus G=c+bY. Its exact original unit is
exp(ct+b^2t^2/2), so polynomial-exponential rigidity gives p=1 and the
exponent in (A3). This proves the scalar necessity in every case.

All complex-analytic necessities descend to characteristic-zero K by
embedding the finitely generated field of phase coefficients and the
hypothetical finite p,psi data into C. The sufficient identities are
algebraic over K. No arithmetic specialization or EIT theorem is used.

## 5. Original finite images, safety, recompositions and jets

If H=0,G=c+bY, let gamma be the minimum strict angular weight of R,
after reversing negative support. For R=0 retain only n=0. For every
original marker Q of degree D,

    Phi_P[Q]=exp(ct+b^2t^2/2)
       *sum_(n=0)^floor(D/gamma) t^n/n!
                  E[Q(Z,W,Y+bt)R(Z,W,Y+bt)^n].     (A13)

This proves all-original-marker finite generation and primitive vacuum.
A finite full image supplies a scalar p exp(psi) by the
[general actual-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image),
so the scalar result is equivalent to the full-image classification.
The same theorem gives literal V=M for nonconstant P, with the constant-
phase exception stated after (A3).

Purity in (A2) forces p=1,psi=0, hence c=b=0. Safety makes the actual
unit polynomial, and the scalar theorem again makes it exp(ct+b^2t^2/2),
which is polynomial only for c=b=0. Conversely zero core leaves strict
angular support, so every marker has cutoff

    E[Q P^m]=0 for m>floor(deg Q/gamma).             (A14)

For b!=0 and R!=0, normalized time degree is at most
D+(1+deg_Y R)floor(D/gamma). For b=0 it is floor(D/gamma); with R=0
use degree at most D (zero if b=0).

For a separate recomposition variable U, the
[general safe-recomposition theorem](GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md)
gives

    S_P=(U-c) if H=0 and G=c constant;
    S_P=(0) otherwise.                             (A15)

The constant phase has the same formula by direct evaluation f(c)=0.
At a safe phase, the
[mixed-cutoff finite-jet and tail theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md#1-the-finite-jet-theorem)
gives r-jet budget floor((1+d/gamma)r) for perturbations of original degree
at most d. An actual coefficient maximal-ideal certificate of index nu
leaves only intervening moments through floor((1+d/gamma)(nu-1)) to
check. For the zero phase take gamma=infinity. These are sufficient
budgets at known safe phases, not least-prefix or computed-index claims.

## 6. Evidence and remaining scope

The [focused exact checker](verify_gaussian_affine_radial_invariant_core_finite_images.py)
and [receipt](gaussian_affine_radial_invariant_core_finite_images_verification.json)
retain bounded original balanced-marker comparisons with the exact
q-localized source, unbalanced-marker zeros, the tracked Gamma endpoint
identity, localized Gaussian Stein and time identities, moving-log
residues through multiplicity six, and finite original branches with
higher W and Y degree. They check finite algebraic scope. Uniform
all-model volume orders, complete comparison and the unrestricted scalar
criterion are the written proof above.

The core is affine in S but arbitrary in Y; general higher-S mixed-parity
invariant cores remain outside this theorem. The strict remainder has no
degree restriction. Imported comparison statements retain their primary-
source audits and pending independent-review status. No independent
replay, formal verification, external mathematical review, publication
or priority is asserted, and no universal Gaussian safety-to-purity claim
is made.
