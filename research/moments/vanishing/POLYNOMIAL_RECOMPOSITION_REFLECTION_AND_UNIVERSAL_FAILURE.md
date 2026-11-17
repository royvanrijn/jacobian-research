# Polynomial recomposition reflects safety and preserves every pure failure

16 November 2026. Internal written proofs; no novelty, independent replay,
external review or formal verification is claimed. See the
[checker](verify_polynomial_recomposition_reflection.py),
[receipt](polynomial_recomposition_reflection_verification.json) and
[source audit](polynomial_recomposition_reflection_source_audit.json).

This closes a GLOBAL question left beyond the
[factorial recomposition theorem](ACTUAL_JOINT_IMAGE_AND_FACTORIAL_RECOMPOSITION.md).
For EVERY pure unsafe Gaussian phase, in every dimension and degree, EVERY
nonzero zero-constant polynomial recomposition is again pure and unsafe.
The same theorem holds on both full exceptional Haar quotients, and in the
general intrinsic Haar setting below. No factorial law, support assumption,
homogeneity, coefficient arithmetic or asymptotic estimate is required.

A stronger reflection theorem explains why: H(P) is all-marker safe exactly
when P-c is all-marker safe for a unique scalar c with H(c)=0. The proof keeps
ALL original markers and the actual formal period image. A finite resolvent
identity and its rank-one quotient at time zero do the work.

Full exceptional MZ remains OPEN: this theorem propagates ANY actual failure
but neither constructs the first one nor rules every failure out.

## 1. Hypotheses and the actual formal image

Let K have characteristic zero, A be a commutative K-algebra, and
h:A->K be K-linear with h(1)=1. For a fixed phase P define

    E_Q(t)=sum_(n>=0) h(QP^n)t^n/n!,
    V_P=sum_(Q in A) K[t] E_Q(t) subset K[[t]].

It is stable under multiplication by t and differentiation, since
E_Q'=E_(PQ). Assume the following actual-image condition, also after scalar
extension to an algebraic closure:

    dim_K(V_P/tV_P)=1.                                      (1)

This is a substantive hypothesis, not a consequence of an arbitrary
functional being normalized. It holds for EVERY polynomial Gaussian phase
over EVERY characteristic-zero field, by
[GMCFORMALCOMP20261115](GAUSSIAN_FORMAL_COMPLETION_AND_POLYNOMIAL_HOM.md):
the exact Gaussian integration source has H_P/tH_P=K and maps onto V_P.
Its image is nonzero and has a series with nonzero constant term, so its
quotient dimension is exactly1. Extension fields have the same source proof.

It also holds for the full Haar setting of
[GHAARTOPSTATE20261116](INTRINSIC_HAAR_STATES_AND_NONSPLIT_PERIOD_EXTENSIONS.md):
X is smooth connected affine, its nowhere-zero algebraic top form Omega
generates H^dimX_dR(X)=C, and the normalized functional is the nonzero
compact top-cycle period. The full twisted top source modulo t is C[Omega].
Both exceptional quotients satisfy these hypotheses. Arbitrary extension
fields of C are allowed using its exact coefficient-base-change statement.
No localization, semisimplification or generic source-rank replacement occurs.

Call P **safe** if for EVERY fixed Q the sequence h(QP^n) is eventually zero.
This definition does not itself require pure vanishing. It is equivalent to
V_P subset K[t]. Call P pure when h(P^n)=0 for every n>=1.

Every nonzero K[t]-submodule L of K[[t]] has L/tL!=0: otherwise L=tL would
make every element divisible by arbitrarily high powers of t inside K[[t]],
forcing L=0. Consequently(1) makes V_P indecomposable as a Weyl module.
This is the elementary time-zero argument retained in GHAARUNITEXT20261116;
here it works equally for Gaussian formal periods, with no convergence premise.

## 2. Exact resolvent reconstruction in every degree

For any Q put C_Q^P(s)=sum_(n>=0)h(QP^n)s^(-n-1), a FORMAL Laurent series
at infinity; it need not converge. Let H(T) be nonconstant of degree d and set

    (H(s)-H(T))/(s-T)=sum_(j=0)^(d-1) r_j(s) T^j.

The polynomial identity gives

    C_Q^P(s)=sum_(j=0)^(d-1) r_j(s) C_(QP^j)^(H(P))(H(s)).  (2)

Indeed expand 1/(H(s)-H(P)) at infinity and multiply its polynomial
numerator. Each coefficient of s^-1 is a finite sum, since H(s) has degree
d>0; coefficientwise application of h is therefore legitimate. No interchange
of divergent analytic series or preservation of a weighted measure is used.

If all the d markers QP^j have eventually zero tails for H(P), their resolvents
on the right are polynomials in 1/H(s). Thus C_Q^P is rational, with finite
poles only among the roots of H. In particular safety of H(P) makes EVERY
original C_Q^P rational with this common finite pole set.

## 3. Rational ordinary marked series force a single exponential centre

**MVRECOMPREFLECT20261116.** Under(1), if EVERY C_Q^P is rational, there is a
unique c in K with

    V_P=K[t] exp(ct),                                       (3)

and hence P-c is safe. If the rational resolvents have poles in one specified
finite set S, then c belongs to S.

Proof. Over an algebraic closure, a proper rational Laurent series at infinity
is a finite sum of (s-c)^-k. Under exponential generating normalization this
term becomes t^(k-1) exp(ct)/(k-1)!. Therefore V_P is a submodule of

    direct_sum_(c in Kbar) Kbar[t] exp(ct).                  (4)

This sum is direct, and every submodule of it is a sum of some displayed
simple modules. Here is an explicit proof sufficient for this use. For a
finite exponential-polynomial sum, choose powers (partial_t-c)^N annihilating
each of its components. The distinct factors are coprime; Bezout/Chinese
remainder polynomials in partial_t isolate each component inside the same
submodule. From any nonzero p(t)exp(ct), apply (partial_t-c)^degp to obtain
a nonzero multiple of exp(ct); multiplication by t generates its full summand.
Thus V_P is a direct sum of these summands. Indecomposability from(1) forces
exactly ONE. Uniqueness and Galois invariance of this sole centre imply c in K.
For general extension fields the same argument is algebraic and formal; no
embedding in C is needed. Pole support gives the last assertion.

For an arbitrary fixed nonconstant H, (2) and(3) now prove the exact criterion

    H(P) safe
       iff there is c in K with H(c)=0 and P-c safe.         (5)

The centre is unique whenever it exists. Necessity follows from the common
root pole set in(2). Conversely, expand H(c+R) with R=P-c and zero constant
term. Its m-th power uses only R^k with k>=m, so every fixed Q tail vanishes
once m exceeds the eventual cutoff for Q,R. More sharply, if H(c+T) has
order e at T=0 and h(QR^k)=0 for k>N_Q, the sufficient cutoff is e*m>N_Q.

The rational-series premise is equivalent to existence of this unique safe
centre: if P-c is safe, every E_Q is exp(ct) times a polynomial, hence every
ordinary C_Q^P is rational. The statement concerns the ACTUAL full image,
not a larger annihilating equation admitting an irrelevant exponential solution.

## 4. Universal persistence of pure failures and finite marker blocks

**MVPUREFAILRECOMP20261116.** If P is pure and unsafe under(1), there is a
fixed Q whose ordinary C_Q^P is NONRATIONAL. For every nonconstant polynomial
H of degree d, at least one marker in the finite block

    Q, QP, ..., QP^(d-1)                                   (6)

has an infinite nonzero marked tail for H(P). In a polynomial or filtered
coordinate algebra these markers have degree at most degQ+(d-1)degP.

Otherwise all ordinary marked resolvents would be rational, so(3) would hold.
But purity puts the constant1 in V_P. A single summand in(3) can contain1
only for c=0, which would make P safe. This contradiction proves existence
of a nonrational Q. Formula(2) proves the block assertion by contraposition.
Nonrationality persists under EVERY field extension: a scalar formal series
is rational iff its Hankel matrix has finite rank; ranks of finite matrices
over K are unchanged by extending K. This also justifies the field-extension
use in section5 without choosing a new numerical or analytic witness.

If H(0)=0, H(P) is pure by expanding its powers in positive powers of P.
Hence EVERY nonzero zero-constant H(P) is again a genuine pure failure.
If H(0)!=0, its pure moments equal H(0)^n, so it is not eventually safe
already with marker1. The same nonrational-block obstruction still applies.
An eventually pure but unsafe P has the same conclusion: its unit period
is a nonzero polynomial, which again forces c=0 in(3); zero-constant
recomposition remains eventually pure.

The marker Q in(6) can be selected once for P, independent of H and d. Its
finite block depends only on d. No theorem that the SAME original Q alone
works for every recomposition, or that all sufficiently late marked indices
are nonzero, is asserted. The earlier exact factorial theorem remains
stronger in precisely those two respects, with its asymptotic constants intact.

## 5. The complete coefficient-ring safety and normalization ideal

**MVRECOMPNIL20261116.** Fix ANY pure unsafe P over K as above and d>=1.
Over ANY commutative K-algebra B let

    R=sum_(j=1)^d b_j P^j, I=(b_1,...,b_d).

The family is entirely pure. Its all-marker safe locus is exactly

    R safe over B iff I is nilpotent.                       (7)

If I^L=0 then R^L=0, proving sufficiency with a marker-independent index.
Conversely reduce a safe input modulo any prime ideal of B and pass to its
fraction field. Every original marked tail remains finite. If any coefficient
survived, section4 over that field would give a pure failure, a contradiction.
Thus every b_j lies in the nilradical. Finitely many nilpotent generators
generate a nilpotent ideal, with no Noetherian assumption on B.

The exact scalar normalization ideal is

    J_R={a in B: a*h(Q exp(tR)) is polynomial for every Q}
       =Gamma_I(B)=union_(L>=0) Ann_B(I^L).                  (8)

Indeed a belongs to J_R iff R is safe over B/Ann(a). Apply(7) there.
For Gaussian and intrinsic Haar sources their retained unique-formal-state
theorems identify(8) with the full polynomial-Hom normalization ideal.
Safety requires1 in this ideal, not merely a nonzero element of it.
For the universal domain K[b_1,...,b_d], Gamma_I=0, whereas the origin fibre
has one normalized polynomial state: failure of naive Hom base change occurs
in EVERY such pure recomposition family, not only a factorial example.

This classifies all coefficient rings and all recomposition degrees for every
pure failure; it does not compute each ordinary prescribed-cutoff ideal or its
least nilpotent thickness. No infinite coefficient-ring Noetherianity is used.

## 6. Evidence, scope and the remaining full exceptional problem

The checker verifies universal polynomial divided-difference identities,
partial-fraction to exponential conversion, Chinese remainder projectors,
and small rational reconstruction examples. The nonrational control uses
Long's original SU2 pair, F=(1+y)(1+z+yz), G=-y in determinant coordinates.
Its exact marked law is (-1)^(n-1)/(n+1), n>=1; hence

    C_G^F(s)=1/s-log(1+1/s).

Its derivative has residues1 and-1 at0 and-1. A rational function's derivative
has zero residues, proving nonrationality exactly. The all-index input is
[Long2607.19012v1, Theorems2.1 and4.2](https://arxiv.org/html/2607.19012v1);
the checker supplies only selected Haar transcriptions. This SU2 example is
NOT invariant under either exceptional group and is not transferred to them.

The universal proofs are elementary algebra plus the retained exact
time-zero source property. No new analytic comparison, asymptotic fit,
support sweep, invariant marker search or external theorem is imported.
The finite checker cannot certify the universal statements by itself.

On either FULL exceptional quotient, any actual pure failure would therefore
generate an unsafe polynomial algebra in ALL zero-constant recompositions,
with the same finite-marker-block rule and exact ring normalization(8).
This forbids a general repair strategy by taking polynomial functions of
an unsafe pure phase. It does not provide the first invariant failure or
split the remaining actual unit-period extension. The full mixed-isotypic,
image-interior, strict-leading opposite-weight cancellation problem remains
open. Closure-first work continues; external review and outreach stay deferred.
