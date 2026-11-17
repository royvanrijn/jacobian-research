# Scalar time equations and actual apolar gaps

16 November 2026. Internal written proofs; no novelty, independent replay,
external review or formal verification is claimed. See the
[checker](verify_scalar_time_gaps.py), [receipt](scalar_time_verification.json)
and [source audit](scalar_time_source_audit.json).

One COMMON scalar differential equation in time, with coefficients formal
in the marker variables but polynomial in time, suffices for the conditional
linear degree bound. Joint holonomicity and a source integration b-function
are unnecessary for this implication. We apply it to a source recording
the ACTUAL apolar deficits, rather than all original marked outputs.
Existence of the required equation for arbitrary formal GVC pairs is NOT
proved. The full formal gap problem remains open.

## 1. An elementary indicial bound with formal marker coefficients

Let K be a characteristic-zero field and

    F(t,s)=sum_alpha F_alpha(t) s^alpha in K[[t,s]], F_0=1.

Suppose ONE nonzero scalar operator annihilates F:

    L=sum_(i,j) a_ij(s) t^i partial_t^j,
    a_ij in K[[s]], with FINITELY MANY pairs(i,j).           (1)

Equivalently, F satisfies an ordinary time equation over
Frac(K[[s]])(t); clear its finitely many coefficient denominators.
The distinction from K[t][[s]] is essential: time degrees and derivative
orders in(1) have common finite bounds.

**MVTIMEINDICIAL20261116.** If EVERY F_alpha is polynomial in t, then
there is one nonnegative integer C with

    deg_t F_alpha <= C|alpha| for all alpha.                (2)

No equation in a marker derivative is assumed. The equation alone does
NOT reduce the individual polynomiality premise to finitely many markers.

Here is an explicit sufficient C from a supplied equation. Introduce
independent parameters A_i and put s_i=A_i z. Divide the transformed
operator by the largest common power of z. It has an expansion

    L_0(t,partial_t)+sum_(beta>=1) z^beta L_beta(t,partial_t),
    L_0!=0,                                                (3)

over K(A). All time powers and derivative orders are still from a fixed
finite set. Let M=max(i-j) over its nonzero coefficient positions. For
L_0=sum c_ij t^i partial_t^j set

    d=max{i-j:c_ij!=0},
    p(N)=sum_(i-j=d)c_ij (N)_j,                            (4)

where(N)_j=N(N-1)...(N-j+1). The polynomial p is nonzero, since the
falling factorials indexed by the distinct j are linearly independent.
Choose B>=0 beyond every nonnegative integer root of p, allowing equality
at the largest such root. Then

    C=max(B,M-d)                                           (5)

suffices. Roots are taken as identities over K(A). This is a supplied-
equation bound, not a general algorithm to produce(1) from a formal source.

Proof. If a nonzero polynomial g has degree N>B, its leading term gives

    deg_t L_0(g)=N+d.                                      (6)

The coefficient is p(N) times the leading coefficient of g. The same
formula holds even if N is smaller than some derivative orders, because
their falling factorials are0. If N+d were negative, p(N) would be0;
therefore that case cannot occur here.

Write F(t,Az)=sum_k g_k(t)z^k, with g_0=1. For deg g_k>B, extracting z^k
in(3) and using(6) gives

    deg g_k <= max_(1<=beta<=k)(deg g_(k-beta)+M-d).         (7)

Terms that vanish are omitted. If deg g_k<=B, the desired bound Ck is
immediate for k>=1. Otherwise induction in(7), with C>=M-d, gives the
same bound. The initial coefficient g_0 has degree0. Finally g_k is a
homogeneous polynomial in the INDEPENDENT A variables, with coefficients
F_alpha(t), |alpha|=k. Hence no cancellation between different original
markers can hide a larger degree. This proves(2).

If F_0 is instead any polynomial, the same argument gives an affine bound
A+C|alpha|. For F_0=1,(5) gives the zero-offset bound without knowing the
degrees of any positive-marker coefficients in advance. All assertions
here are over fields; no inversion of a nonunit indicial coefficient in a
general nonreduced coefficient ring is assumed.

## 2. Why this is broader, and what it does not decide

The retained [MVHOLLINEAR20261115](HOLONOMIC_SOURCES_AND_CUTOFF_LOCI.md)
uses a joint holonomic source and a source b-function. Its stronger finite
marker polynomiality/decision statement is retained. Section1 needs only
the ordinary time equation, but retains ALL individual polynomiality
premises. It supplies a simpler conditional degree bound, not a replacement
for those finite decision procedures.

Indeed for every N>=1 the source

    F_N=1+s^N(exp(t)-1)

has pure seed1, satisfies the SAME equation(partial_t^2-partial_t)F_N=0,
and has polynomial coefficients at every marker degree<N. Its coefficient
at degree N is not polynomial. Thus no finite marker test depending only
on that scalar equation can imply all-marker polynomiality.

The formal coefficient field strictly enlarges the source class. Let
phi(s)=sum_(k>=2)s^(k!) and

    B(t,s)=sum_(m>=0)(t phi(s))^m/(m!)^2.

It satisfies(t partial_t^2+partial_t-phi(s))B=0. Every marker coefficient
is polynomial in time, since ord(phi)=2. The source is not jointly
D-finite over K(t,s): the Bessel eigenoperator would force phi algebraic,
contrary to the retained lacunary recurrence-gap proof. Section1 still
applies. This is a source-class boundary, not a new GVC counterexample.

It would be insufficient to replace(1) by coefficients in K[t][[s]] with
unbounded time degrees. For example

    F=1+sum_(k>=1)t^(k^2)s^k

has all marker coefficients polynomial but no linear degree bound.
It is annihilated by the tautological first-order operator
F partial_t-partial_t(F), whose coefficients lie in K[t][[s]]. This does
not satisfy the common finite time-degree requirement of(1).

## 3. The reversed source of actual apolar heights

Let lambda be an arbitrary formal constant-coefficient symbol over a
characteristic-zero field, and P a NONZERO polynomial with
T^mP^m=0 for every m>=1, T=lambda(partial). Then lambda(0)=0 and define

    h_m=max{k>=0:T^kP^m!=0}, b_m=m-h_m, m>=1;
    R_m(z)=sum_(k>=0) T^kP^m z^(m-k)/k!, R_0=1;
    Psi(t,z,x)=sum_(m>=0) R_m(z)t^m/m!.                    (8)

Each R_m is a polynomial in z with coefficients in K[x], because purity
gives h_m<=m-1. Its valuation at z=0 is EXACTLY b_m; output x remains
independent and cannot be replaced by evaluation at0. Thus Psi belongs
to K(x)[[t,z]], with Psi(t,0,x)=1. There is no undefined substitution into
an infinite exponential here: each fixed t^m coefficient is a finite sum.

The retained actual-gap reduction gives

    every fixed original marker is eventually zero
      iff b_m tends to infinity
      iff Psi belongs to K(x)[t][[z]].                     (9)

For the last equivalence, each fixed z^d receives exactly the terms with
k=m-d. If b_m>d those terms vanish. Conversely, polynomiality for d<=D
excludes h_m>=m-D at every sufficiently large m. This proves the converse
by taking the maximum of those FINITELY MANY time degrees.

**GVCREVERSEDHEIGHT20261116.** Suppose this ACTUAL reversed source has
one nonzero scalar time equation with coefficients in K(x)[[z]][t]. Under
the equivalent eventuality premises(9), section1 computes a sufficient
integer C from that equation such that

    m<=C b_m for EVERY m>=1,
    T^m(QP^m)=0 whenever m>C deg Q.                         (10)

Indeed the nonzero t^m z^(b_m) coefficient implies m<=C b_m by(2).
If delta=degQ and m>C delta, then b_m>delta, so
T^(m-delta)P^m=0. The retained formal commutator implication
T^(m-delta)f=0 => T^m(Qf)=0 proves(10). The zero marker is automatic,
and delta0 is the pure premise. No phase-degree factor is lost in(10).

A sufficient gate is any NONZERO finite recurrence for R_m with
coefficients in K(x)[[z]][m], of bounded polynomial m-degree, holding
for all sufficiently large m. For the exponential generating function,
replace m by theta_t=t partial_t and the forward shift by partial_t.
This is the injective Ore substitution S m=(m+1)S. The exceptional finite
initial indices give a polynomial-time inhomogeneous term. Applying enough
time derivatives kills it and gives an equation of the specified form.
No such recurrence is proved or produced for a general formal symbol.

This gate concerns a TWO-INDEX UNMARKED source that records actual heights.
It does not require the full original-marker source, or a recurrence for
an abstract height sequence unrelated to actual outputs. It is a sufficient
criterion, not a claimed characterization of every linearly safe pair.

## 4. Retained closures and the genuinely missing hypothesis

For algebraic formal lambda the construction in
[GVCALGFORMALSOURCE20261116](ALGEBRAIC_FORMAL_SYMBOLS_AND_GLOBAL_LINEARITY.md)
also gives a joint D-finite series

    H(u,v,x)=sum_(k,m) T^kP^m u^k v^m/(k!m!).              (11)

Use its proper factorial diagonals, retain the two exponent indices
INDEPENDENTLY, and evaluate the auxiliary Taylor variables using the bound
deg P^m<=Dm, uniform in k. Do not take the final k=m diagonal used there.
A scalar v-equation for H transforms under

    u=1/z, v=t z, partial_v=z^(-1)partial_t

to a NONZERO scalar time equation for Psi, after clearing z denominators.
The rational coordinate map is invertible(t=uv,z=1/u). The substituted
series is defined coefficientwise because every u coefficient at fixed m
is a finite polynomial, and purity removes all negative z exponents.
Thus section3 gives another proof of conditional linearity for algebraic
symbols, without a full marker b-function.

Formal units preserve the exact heights. Representatives modulo J_P,
the common annihilator of ALL phase powers, preserve every unmarked
output. Therefore a supplied representative with this reversed-source
equation transfers the resulting HEIGHT bound to the original pair;
the original formal commutator then gives the SAME sufficient C. It need
not transfer the actual marked outputs or their least cutoff. This
observation also applies to nonalgebraic representatives satisfying the
new gate, should they be supplied.

For example lambda in J_P gives Psi=exp(tzP), with the first-order
equation(partial_t-zP)Psi=0 and sufficient C1. The phase can be concise:
P=a+2wb+w^2c has R=AC-B^2/4 in J_P. Multiplying R by the formal factor
W-phi(B) produces a symbol with no nonzero algebraic multiple, by the
retained graph-divisor proof, while still lying in J_P. Its full raw source
is nonholonomic. This case was ALREADY safe by the phase-quotient theorem;
it illustrates why the reversed-source gate has no algebraicity necessity.

What remains unknown is the existence of a suitable finite time equation
or recurrence for GENERAL individually-eventually-safe formal pairs. Formal
local Noetherianity, an abstract first apolar algebra, or a joint-source
equation with arbitrarily growing time degrees does not supply that gate.
The monomial compact-face theorem is retained separately; it is not being
asserted that its entire reversed source is D-finite.

## 5. Weakening the rational raw-source hypothesis does not enlarge it

For comparison retain the FULL original-marker source

    F_raw=sum_(m,alpha) T^m(x^alpha P^m)t^m s^alpha/
                                      ((m!)^2 alpha!).     (12)

**GVCTIMERIGID20261116.** For EVERY NONZERO polynomial phase P,
F_raw has a scalar ordinary time equation over K(x,s,t) if and only if
lambda is algebraic over K(Z). This holds without purity or eventuality.
Thus replacing joint holonomicity by rational-marker time holonomicity
does NOT enlarge that raw-source method. The field of formal marker
coefficients in section1 and the DIFFERENT source(8) matter.

The necessity proof in
[GVCJOINTSOURCERIGID20261116](FORMAL_SYMBOL_SOURCE_RIGIDITY.md)
actually uses ONLY time derivatives. Gauge by exp(-s dot x), substitute
x=A/h,t->t h^D, and take the leading h coefficient. A scalar time equation
survives: clear denominators and divide its common h power before setting
h=0; no h or marker derivative occurs. The resulting Bessel source is

    H_0=sum_(m>=0)(P_D(A)t lambda(s))^m/(m!)^2.

The time operator A_t=t partial_t^2+partial_t satisfies
A_t^j H_0=(P_D(A)lambda(s))^j H_0. A scalar time equation makes all pure
time derivatives span a finite-dimensional space over K(A,s,t), so these
eigenvalue powers are linearly dependent. Cancel H_0 and descend the
independent parameters to get algebraicity of lambda over K(s).
The reverse implication is the retained algebraic full-source theorem.
Nonzero CONSTANT phases are included; P=0 must be excluded, since its raw
source is exp(s dot x) for every lambda. Unit stripping consequently has
the same algebraic-principal-divisor boundary as before.

## 6. Primary source audit and verification boundary

We read all32 physical pages of Shaoshi Chen's March28,2024 Wuhan slides,
[Degree Growth of Holonomic Polynomial Sequences and Applications](https://perso.math.u-pem.fr/liao.lingmin/SAP-WH/CHEN_Shaoshi.pdf),
joint work credited there to Bell,Chen,Wan,Wang,Yin. Logical slide6 gives
a linear lower bound on a SYNDETIC subset when degrees are unbounded;
it does NOT state a rational slope with bounded error at every index.
Logical slide9 discusses bounded-degree/separability equivalence. The
degree-growth theorem's proof is not present in the slides. No theorem
from those slides is imported here. Equations(3)-(7) are the full elementary
proof needed for our DIFFERENT conditional upper bound. The slides were
a useful research lead, not evidence for a stronger unproved statement.

Small exact controls check normal-order indicial extraction, recurrences
for polynomial time coefficients, a reversed actual-height profile and
its bound, top-output/Bessel transcription, a delayed nonpolynomial marker,
and the forbidden unbounded-time-degree completion. They do not prove
existence of a general time equation, implement a recurrence producer,
decide the full formal gap, or increase any assurance flag.

Both full exceptional Haar MZ cases remain open at actual global unit
couplings. Closure-first work and deferred external review/outreach persist.
