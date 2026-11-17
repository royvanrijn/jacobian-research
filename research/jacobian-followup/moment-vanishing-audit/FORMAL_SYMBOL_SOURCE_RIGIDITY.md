# Exact algebraicity boundary for joint formal-symbol sources

16 November 2026. Internal written proof. No independent replay, formal
verification, external review or novelty is claimed. The
[checker](verify_formal_symbol_source_rigidity.py),
[receipt](formal_symbol_source_rigidity_verification.json) and
[source audit](formal_symbol_source_rigidity_source_audit.json) contain exact
transcription controls; they do not replace the universal argument.

For EVERY nonzero polynomial phase, the full raw GVC Taylor source is
D-finite precisely when the formal symbol is algebraic. Thus the previous
algebraic-symbol source theorem is an exact boundary of this method.
Individual-marker eventuality does not force source holonomicity, even when
the principal formal symbol ideal has no algebraic generator up to a unit.
The general eventuality-to-linearity question nevertheless remains OPEN:
the nonholonomic examples below have an explicit sharp linear cutoff.

## 1. Setup and the specialization lemma

Let K be any characteristic-zero field, let lambda belong to K[[Z_1,...,Z_n]],
and let P in K[x_1,...,x_n] be NONZERO. Retain independent output x and all
marker variables s. Define the actual raw source

    F(t,s,x)=sum_(m,alpha) lambda(partial_x)^m(x^alpha P(x)^m)
                             t^m s^alpha/((m!)^2 alpha!).       (1)

Every coefficient is a finite differential action on a polynomial.
D-finiteness here means finite span of all (t,s)-derivatives over K(x,t,s).
No derivatives in the parameter x are part of that hypothesis.

We use this elementary relative specialization lemma. Suppose H belongs to
L[[h,t,s]], is D-finite in (t,s) over L(h), and h is a parameter. Then
H|_(h=0) is D-finite in (t,s) over L. For each coordinate y among (t,s),
choose a nonzero ordinary differential annihilator in partial_y, with
coefficients in L(h,t,s). Clear denominators, obtaining coefficients in
L[h,t,s], and divide their largest common power of h. The resulting
operator has a nonzero specialization at h=0. Its identity on H can be
specialized coefficientwise, because it contains no partial_h. Thus the
specialized series has a nonzero ordinary differential equation in each
coordinate. These separate equations are the standard equivalent criterion
for multivariate D-finiteness (bounded pure derivative orders reduce the
mixed derivatives to a finite rational span). This proves the lemma without
assuming any differential equation in h or any noncharacteristic hypothesis.

## 2. Gauge and top-output degeneration

Multiplying (1) by exp(-s dot x) gives the coefficientwise identity

    G(t,s,x)=sum_m t^m/(m!)^2
                  lambda(s+partial_x)^m P(x)^m.               (2)

For clarity the finite Taylor operator in each coefficient is

    lambda(s+partial_x)^m P(x)^m
      =sum_beta (partial_s^beta lambda(s)^m)/beta!
                                      *partial_x^beta P(x)^m. (3)

Only |beta|<=Dm occurs, where D=degP. Formula (2) follows either by
exponential conjugation or by summing the original markers first:
sum_alpha x^alpha s^alpha/alpha!=exp(s dot x). In particular every t^m
coefficient of G is a polynomial in x of degree<=Dm with coefficients in
K[[s]]. This finite x-degree assertion is essential; F itself has the
exponential x-dependence before the gauge is removed.

Introduce algebraically independent parameters a=(a_1,...,a_n) and h. In
(2) substitute x=a/h and replace t by t h^D. The bound just proved gives
an actual series

    H(h,t,s,a)=G(t h^D,s,a/h) in K[a][[h,t,s]].                (4)

At each fixed m its h-degree is at most Dm; no negative h-power remains.
If F is D-finite in (t,s) over K(x), so is G: multiplication by the
rank-one exponential gauge preserves the finite derivative span. The
substitution in (4) embeds K(x,t,s) into K(a,h,t,s), and rescales partial_t
by h^(-D); it therefore preserves D-finiteness in (t,s) over K(a,h).
This does not assume joint D-finiteness in the new parameters.

The relative specialization lemma now applies over L=K(a). If P_D denotes
the highest homogeneous part of P and c=P_D(a), then c is NONZERO and

    H_0(t,s,a)=H(0,t,s,a)
             =sum_m (c t lambda(s))^m/(m!)^2                 (5)

is D-finite in (t,s) over K(a). Indeed only beta=0 and the highest output
part P_D(x)^m survive in (3)-(4). Constant nonzero phases are included:
D=0, c=P, and (5) already equals G. The excluded phase P=0 would give
F=exp(s dot x) for every symbol, so the exclusion cannot be dropped.

## 3. The Bessel eigenoperator forces algebraicity

Set A=t partial_t^2+partial_t. The factorial normalization gives

    A H_0=c lambda(s) H_0,
    A^j H_0=(c lambda(s))^j H_0 for every j>=0.              (6)

The second equality uses the independence of c lambda(s) from t. If H_0
is D-finite, all the left sides lie in one finite-dimensional derivative
span over K(a,t,s). Consequently some nonzero polynomial in T over that
field vanishes at c lambda(s): cancel H_0, whose constant coefficient is1,
from a linear dependence of the right sides. Since c!=0, lambda is
algebraic over K(a,t,s).

This descends to K(s). Clear denominators in a polynomial relation for
lambda and expand its coefficients as polynomials in the independent
parameters (a,t). The resulting identity holds in K[[s]][a,t]; at least one
nonzero coefficient polynomial in (s,T) must vanish at T=lambda(s). Thus
lambda is algebraic over K(s). No algebraic closedness assumption is used.

**GVCJOINTSOURCERIGID20261116.** For every characteristic-zero field,
every formal symbol and every NONZERO polynomial phase in any dimension
and degree, the full raw source (1) is D-finite in (t,s) over K(x) if and
only if the symbol is algebraic over K(Z). The reverse implication is
exactly the retained
[GVCALGFORMALSOURCE20261116](ALGEBRAIC_FORMAL_SYMBOLS_AND_GLOBAL_LINEARITY.md),
with its proper factorial diagonals and finite rho evaluation. The new
necessity proof is (2)-(6), together with the elementary relative
specialization lemma. It applies without purity or eventuality assumptions.

In particular the constant phase P=1 alone detects the algebraicity of
every formal symbol through joint-source holonomicity. Separate
D-finiteness of finitely many marked coefficients is not substituted for
the full joint hypothesis.

## 4. Exact boundary after stripping units

**GVCUNITSOURCEMAX20261116.** For any fixed NONZERO P, there exists a formal
unit U such that the raw source for lambda/U is holonomic if and only if
the principal symbol ideal has an algebraic formal generator. This follows
directly from section3 and the retained
[GVCUNITEXACT20261116](ALGEBRAIC_FORMAL_SYMBOLS_AND_GLOBAL_LINEARITY.md).
Every zero predicate and every cutoff is invariant under these units.
Thus stripping a unit followed by a joint-holonomic-source proof covers
EXACTLY the already closed algebraic-divisor class. Choosing another
nonzero polynomial phase cannot enlarge that class by making the raw
source holonomic. This is a boundary for this specified proof method,
not a proof that no other finiteness method can close the remaining OP.

## 5. Eventuality does not force the missing finiteness

Let

    phi(v)=sum_(k>=2) v^(k!),
    lambda(X,Y,Z)=X-phi(Y),     P(x,y,z)=z^D, D>=1.          (7)

The lacunary phi is not D-finite, hence is not algebraic. Here is a direct
formal proof. A nonzero polynomial differential equation for a univariate
series induces, after coefficient extraction and a finite shift, a
nontrivial recurrence

    sum_(j=0)^r p_j(N) a_(N+j)=0, p_r!=0,                  (8)

valid for all sufficiently large N. If r=0 the coefficients eventually
vanish. If r>=1, a block of r consecutive zeros beyond every root of p_r
propagates by (8) to all subsequent coefficients. The factorial gaps give
arbitrarily late such zero blocks followed by a coefficient1. Both cases
contradict (8). Algebraic formal series are D-finite by the finite separable
field derivative-span argument already retained. The proof works over
every characteristic-zero K, with no analytic convergence theorem.

Moreover the principal ideal (lambda) has NO nonzero algebraic formal
element, even with extra unused variables. If mu!=0 is algebraic, choose
a cleared minimal polynomial B(X,Y,Z,T) with B(X,Y,Z,0)!=0. A zero constant
coefficient would make T a factor of its irreducible minimal polynomial,
contradicting mu!=0. If mu is divisible by lambda, the well-defined formal
substitution X=phi(Y) gives

    B(phi(Y),Y,Z,0)=0.                                      (9)

At least one nonzero coefficient polynomial in Z of B(X,Y,Z,0) would then
make phi algebraic over K(Y), a contradiction. In particular no formal
unit can convert this lambda into an algebraic symbol.

Nevertheless (7) is pure and has the SHARP integer cutoff1. Its operator
acts only on x,y, has no constant term, and lowers total x,y degree by at
least1. It commutes with multiplication by z. For every polynomial Q,

    lambda(partial)^m(Q z^(Dm))=0 when m>deg_(x,y)Q,

hence when m>degQ. For Q=x^q at m=q the output is q! z^(Dq)!=0, proving
sharpness. This works for arbitrary D, including the unresolved OP's
n=3,D>=3 range. It is not a counterexample to eventuality-to-linearity.
By section3 its actual full joint source is NOT D-finite.

**GVCQUALNONHOL20261116.** Even purity together with individual-marker
eventuality, and indeed an explicit sharp linear cutoff, does not imply a
joint holonomic source or an algebraic principal symbol divisor up to a
unit. Thus the proposed shortcut from eventuality to algebraic-source
finiteness is false. The example is a boundary control for a universal
method question, not a replacement research programme of easy safe phases.

## 6. Remaining work and verification boundary

OP-GVC-FORMAL-LINEAR stays OPEN for genuinely nonalgebraic principal formal
divisors in n>=3 and polynomial phase degree>=3. Section3 proves that its
general resolution needs a mechanism beyond applying the retained raw
joint holonomic source, even after stripping units. Formal-local
Noetherianity or a local Bernstein--Sato polynomial does not by itself
provide the required polynomial Weyl source or linear time/source degree
control. No general impossibility of other formal-local arguments is claimed.

The checker verifies finite coefficient instances of the exact gauge,
output-degree bound, top-output degeneration and Bessel eigenoperator;
tests both a smooth nonalgebraic-graph jet and algebraic/constant controls;
and checks sharp transverse differential actions. It cannot prove the
infinite lacunary transcendence or universal D-finite rigidity by truncation.
Those are the written arguments above. No general algebraicity recognizer,
source producer or formal-local integration algorithm is implemented.

Both complete exceptional Haar MZ cases remain open at the actual
complement zero-value logarithmic extension. Closure-first priority and
deferred external review/outreach remain in force. All protected Gaussian,
Haar, binary, cubic-prefix, ordinary-ideal and radial-transfer results remain
unchanged. No new numerical endpoint or universal positive GVC is claimed.
