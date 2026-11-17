# Logarithmic trace-defect spectrum and target-derivative domains

18 November 2026. Finite normal complex surface covers and actual
polynomial plane Keller sources. Ordinary codifferents are used throughout.

Let B=C[P,Q], F=C(P,Q), L/F finite and separable, and S the full
finite normal integral closure. Put

    Lambda={a in L:Tr_(L/F)(a*S) subset B},  T=Lambda/S.

Let H be a reduced equation of the ENTIRE branch divisor, an empty
product being one. Choose a constant target derivation D transverse
to every branch component. The height-one calculations in
[codifferent preimages](CODIFFERENT_PREIMAGES_UNDER_TARGET_DERIVATIONS.md)
give H*Lambda subset S, D(S) subset Lambda, and the exact derivative
pole -q-e of a pole of order q at a ramified prime of index e.

The logarithmic vector field W=H*D preserves both S and Lambda and
induces a B-LINEAR endomorphism of T. At a branch component Gamma,
after division by the nonzero function D(H)|Gamma, its characteristic
polynomial is

    product_(E over Gamma) product_(i=1,...,e_E-1)
         (Z+i/e_E)^(f_E),                                  (1)

where f_E is the residue-field degree. It is semisimple, with trace

    -rank_(C(Gamma))(T_Gamma)/2.                            (2)

Equal fractions from DIFFERENT ramified factors must be grouped:
e=2,i=1 and e=4,i=2 have the same eigenvalue. Actual prime-factor
decomposition uses the finite local algebra, not a false claim that
all displayed fractions are distinct.

For an actual Keller source, a second consequence is analytic:
a generic target derivative has a NON-DENSE natural maximal domain
on the source Bergman space of a sufficiently small target ball
meeting a generic ramified image point. The missing domain modes
are exactly the first e-1 normal jets of a canonical coefficient.
No automatic density, velocity integrality or Keller inverse is proved.

## 1. The logarithmic field acts on the actual trace defect

At a source divisor of normal index e, the ordinary codifferent test is

    v_E(lambda)>=1-e.

Multiplication by the reduced branch H adds e at every ramified
factor and a nonnegative amount elsewhere. Hence H*Lambda subset S.
For s in S, D(s) belongs to Lambda, so W(s) belongs to S.

For lambda in Lambda, either it is regular at E, when
v_E(Dlambda)>=1-e, or it has a pole q with 1<=q<=e-1,
when transversality gives v_E(Dlambda)=-q-e. Multiplication by H
then gives v_E(Wlambda)>=1-e. The unramified tests are regular.
All height-one tests together imply W(Lambda) subset Lambda.

Thus W induces an additive endomorphism Wbar of T. For b in B,

    W(b*lambda)-b*W(lambda)=H*D(b)*lambda in S.

This proves B-linearity. The endomorphism is annihilated by H,
just as T is. This is a logarithmic action on the torsion module;
D itself does NOT preserve Lambda in the noninverse case.

## 2. Compute the full generic spectrum, including residue degrees

Fix a branch component Gamma and localize B at its reduced equation.
The quotient T_Gamma is a vector space over k=C(Gamma), since H
kills it and all other branch factors become units.

After completion and a finite unramified residue extension, the
normal algebra splits into its ramified DVR factors with

    t=u^e,

where t is a target normal parameter. This is tame characteristic-zero
DVR theory, at the same exact scope as the codifferent proof. On each
factor the quotient has the k(E)-basis

    u^-1,...,u^-(e-1) modulo the regular ring.                (3)

Tangential residue coefficients contribute their usual degree f_E
over k. Normalize Wbar by the scalar

    c_Gamma=D(H) mod Gamma,  c_Gamma!=0.

Modulo regular elements this normalized operator equals
(t/D(t))*D. Indeed H=t*h with h a unit; its leading normal
coefficient is h*D(t), and any difference in the normalized
coefficient divisible by t sends all modes in (3) into the
regular ring.

In the strict coordinate t=u^e,

    (t/D(t))*D(u)=u/e.

The derivative of a tangential residue coefficient, after multiplication
by t, contributes a regular element when paired with u^-i, i<e.
Therefore the induced operator sends u^-i to -(i/e)*u^-i.
This gives (1) after splitting, and hence over k by descent.

A product of distinct linear polynomials annihilates the normalized
operator, using the UNION of all fractions occurring in (1).
Its minimal polynomial is squarefree, so it is semisimple. The
multiplicity of an eigenvalue is the SUM of f_E over all matching
fractions. Lagrange projectors recover these aggregate eigenspaces;
they do not separate primes that have the same fraction.

Summing -i/e gives -(e-1)/2 for each residue copy. Thus, with
r_Gamma=sum_E f_E*(e_E-1),

    Tr(Wbar|T_Gamma)=-(r_Gamma/2)*c_Gamma.                   (4)

In particular the induced operator is nonzero whenever the generic
defect is nonzero. Formula (4) also forbids representing this generic
operator as a commutator of two k-linear endomorphisms when r_Gamma>0,
because the trace of such a commutator is zero. NO such commutator
representation from source Ga actions is asserted.

The ordinary discriminant multiplicity along Gamma is r_Gamma:
the tame different exponent is e_E-1 and finite pushforward multiplies
it by f_E. Thus the logarithmic trace records half that multiplicity.
This is a derived description of the existing finite defect, not a
new proof that the defect vanishes.

## 3. A small target ball needs no singular-point extension theorem

Now assume the ACTUAL source R=C[x,y] and Jac(P,Q)=1, so
Omega=dx wedge dy=dP wedge dQ. Take an omitted ramified prime E.
Choose a smooth generic point p of E whose target lies on a smooth
branch point of its image Gamma. Avoid other branch curves, singular
points, ramification of residue-curve maps, and the images of other
exceptional finite sets. At every point over that target the finite
map has its tame smooth local normal form. Properness of the finite
map permits a sufficiently small target ball V whose whole finite
preimage lies in these good neighborhoods.

Write Omega_X=K^-1(V) for the original source open and

    H_V=A^2(Omega_X),  Dom(D)={f in H_V:Df in H_V}.           (5)

The domain may be disconnected; the Hilbert space is defined on all
components. No density assumption is built into the notation.

In a polydisk about p, with E={u=0}, write

    t=u^e,  Omega=g(t,z)*e*u^(e-1) du wedge dz,

where g is a holomorphic unit. The source fills this polydisk minus E.
For every f in H_V, its canonical coefficient

    F_f=e*g*u^(e-1)*f

is square integrable on the punctured polydisk and holomorphic there.
Its Laurent expansion has no negative normal powers: such a power
has divergent square integral. Thus F_f extends holomorphically over E.
Mean value and Cauchy estimates give continuous functionals

    f |-> partial_u^j F_f(p),  j>=0.                       (6)

Their local norm is bounded by the global H_V norm, because Omega
is the original unit source volume. This is an interior evaluation
at a missing divisor, not an estimate at the outer boundary of V.

No quotient-singularity L2 theorem is used here: V was deliberately
chosen at a GENERIC smooth branch point. For a larger ball seeing
all omitted primes, one must separately use the established ordinary
L2/canonical extension theorem, including its finite-point inputs.

## 4. The derivative removes exactly e-1 local pole modes

In the chosen coordinates a regular target derivative is

    D= a(t,z)/(e*u^(e-1))*partial_u + b(t,z)*partial_z,

with a a unit because D is transverse. Put F=F_f. The canonical
coefficient of Df has the singular part

    (a/e)*u^-e*(u*partial_u F-(e-1)*F),                     (7)

plus holomorphic terms from b and the regular unit g.
If F=sum_(k>=0) A_k(z)*u^k, the coefficient multiplying A_k
in (7) is k-(e-1).

For k=0,...,e-2 it is nonzero; its exponent k-e is negative.
The first such nonzero coefficient cannot cancel against any later
normal term or a holomorphic tangential/unit term. At k=e-1 the
coefficient is ZERO, and all k>=e give nonnegative exponents.
It follows that the local canonical coefficient of Df is holomorphic
if and only if

    F is divisible by u^(e-1),
    equivalently f is regular across E.                    (8)

In particular every member of the natural domain (5) satisfies

    partial_u^j F_f|E=0, j=0,...,e-2.                       (9)

This is the exact INTERIOR condition; outer-boundary behavior can
impose additional requirements, so (9) alone is not claimed to be
the global maximal-domain description.

At least one evaluation in (9), already j=0 at a suitable p, is
NONZERO on H_V. The finite canonical module has a global section
lambda*Omega not vanishing at p. Equivalently choose a global
Lambda generator attaining value 1-e at generic E. Normalization
and compactness over the small ball, together with the same tame
charts at every preimage, make lambda square integrable on Omega_X.
Thus F_lambda(p)!=0 after choosing p away from its zeros.

Hence Dom(D) lies in the kernel of a nonzero continuous functional:

    closure(Dom(D)) is a PROPER subspace of H_V.            (10)

This proves NON-DENSITY, not merely unboundedness or lack of a norm
estimate. On polynomial H_V vectors, the missing negative modes
are precisely the ordinary codifferent modes measured in Section2.

## 5. The existing all-boundary ball and the inverse criterion

For the bounded target ball seeing every omitted prime in
[ordinary Bergman trace duality](BERGMAN_TRACE_DUAL_AND_MULTIPLIER_DOMAIN_RIGIDITY.md),
the separately established identity R intersect H=Lambda supplies
the same nonzero global lambda in H. The local argument at p therefore
proves (10) on that larger Hilbert space as well whenever a ramified
prime is present. Its finite-point/quotient inputs are NOT bypassed
by the small-ball proof.

If the Keller map is invertible, pullback identifies this Hilbert
space with the ordinary Bergman space of the target ball, and D with
a constant derivative. Target polynomials lie in its domain and are
dense, as follows from the usual orthogonal monomial expansion on
a ball. Conversely a noninverse has a ramified prime in its full
finite normalization by the retained finite-normalization inverse
criterion, so the preceding non-density applies. Thus, WITH the
specified ordinary L2 input for that whole ball,

    Dom(D) dense  iff  the Keller map is invertible.         (11)

Neither density nor stability of the source polynomial space in H
is automatic. A polynomial vector field on the whole source may
have a non-dense maximal Bergman domain on this preimage domain.
The source is unbounded; no missing divisor is part of the domain.
The interior jet functionals explain exactly why a naive
differentiation/density argument would already assume inversion.

## 6. Exact support and limitations

The checker computes actual logarithmic Laurent modes, characteristic
polynomials including coincident eigenvalues from different factors,
traces and spectral projectors, the B-linearity correction, and the
normal-jet cancellation with variable units/tangential coefficients.
Exact radial exponents distinguish ordinary L2 modes from their
nonintegrable derivatives, including the logarithmic thresholds.

The finite controls support the written unrestricted proof. Tame
DVR normal forms, normal height-one tests, finite-map compactness,
Laurent removability, mean-value/Cauchy bounds and the ball Bergman
monomial expansion are the adopted elementary foundations. The
large-ball statement additionally uses the precise retained ordinary
L2 theorem; it does not identify other moment lattices with Lambda.

No automatic codifferent multiplication, dense derivative domain,
target Ga descent, trace perfection, HC4 or JC2 inversion is proved.

Replay:

    .venv/bin/python research/hessian/hc4-triangular-graphs/verify_logarithmic_trace_defect.py
