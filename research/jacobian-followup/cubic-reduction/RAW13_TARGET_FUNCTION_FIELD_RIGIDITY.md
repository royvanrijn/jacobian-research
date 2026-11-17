# Uncapped rational-target rigidity and complete first nonlinear target layer (KCUB42)

Subsequent [KCUB44](POSITIVE_SOURCE_QUADRATIC_FIELD_OBSTRUCTION.md) proves
that every nonzero positive-axis quadratic source shear has no rational
Keller repair through source degree eight. It also extends the original
target-field intersection unchanged through degree nine. The transformed
field at degree nine remains unclassified;prior positive-source-gate
proposals below are stage history.

Subsequent [KCUB43](QUADRATIC_TARGET_AFFINE_DEGREE8.md) closes the proposed
joint affine degree-eight gate:every nonzero shear has at least16 terms
in every affine source/output frame,so none improves raw13. Its nonlinear
family minimum is only bounded[16,19];the nineteen fixed-source count below
remains correctly scoped. Source-affine-open proposals below are stage history.

Fix the established three-variable degree-seven map H. For every
characteristic-zero field k,

    k(H1,H2,H3) intersect k[a,b,c]_(degree<=7)
        = span_k(1,H1,H2,H3).                         (1)

Thus any rational target function whose pullback is a polynomial of degree
at most seven is affine. There is **no target-degree or denominator bound**.
The assertion persists in every invertible affine source/target frame of H.
In particular a target-only sequence, even with temporary degree growth,
cannot end at degree<=7 with a smaller map: its total rational target
operation must be affine, and [KCUB41](RAW13_AFFINE_ORBIT_MINIMUM.md) gives
at least thirteen nonlinear occurrences.

The first nonlinear layer is also completely determined. Put

    Q=Y1^2,
    S=Y2^3+27*Y1*Y3^2-9*Y2*Y3,
    T=Y2^2-4*Y1*S.

Then

    k(H) intersect k[a,b,c]_(degree<=8)
      = span_k(1,H1,H2,H3,Q(H),S(H),T(H)).              (2)

If a rational target map R has polynomial Keller pullback R(H) of degree
at most eight, R is polynomial and, after invertible affine output
normalization, is exactly

    R(Y)=Y+(0,lambda,mu)*Y1^2.                         (3)

Its inverse subtracts the same vector times Y1^2. Every constant and tangent
coefficient is covered by the normalization. When (lambda,mu)!=0, the
pullback has degree eight and, **with the original source coordinates
fixed**, needs at least nineteen nonlinear occurrences even under arbitrary
affine output mixing. This bound is attained by (H1,H2+H1^2,H3), with row
counts (2,11,6). It is a threshold control, not a global improvement over H.

The nineteen-term bound is not asserted after a source affine change.
The remaining degree-eight target route is the explicit two-parameter
family (3) jointly with such source changes. Nonlinear source operations,
unweighted constructions and higher-degree simplicity remain open. No
headline positive bound or unrestricted 62/12 profile count changes.
The full global simplicity goal remains open.

## The full cubic target function field

The saved map, with v=1+2ab, is

    H1=a-3*a^2*b-a^3*c,
    H2=3*v^2*a*c+(9*v^2-6*v-2)*b,
    H3=v^3*c+(6*v^2+2*v)*b^2.

Write its target coordinates as (s,t,u). Its tangent map is I, so
substitution k[s,t,u]->k[a,b,c] is injective: the lowest nonzero homogeneous
part of any polynomial is unchanged under H=X+terms of order>=2. It also
embeds rational function fields. Put

    g=1-3ab-a^2c,   w=1/a+2b,   W=s*w=g*v.

Direct polynomial identities, independently checked after clearing all
source denominators, are

    W^3-W^2+s*t*W-2*s^2*u=0,
    E=3*W^2-2*W+s*t=g,
    a*E=s,
    2*s*b=W-E,
    2*s^2*c=5*E^2-3*W*E-2*E^3.                       (4)

Consequently the source field is

    k(a,b,c) = k(s,t,u)[W]/(W^3-W^2+s*t*W-2*s^2*u),

with source recovery

    a=s/E,
    b=(W-E)/(2*s),
    c=(5*E^2-3*W*E-2*E^3)/(2*s^2).                  (5)

The cubic is irreducible over every k(s,t,u). To see this, rescale W=s*w.
The equation becomes s*w^3-w^2+t*w-2u=0. In the full polynomial ring this
polynomial is irreducible, since it is linear in u with unit coefficient
and its quotient is k[s,t,w]. As a polynomial in w it is primitive
(the coefficient of w^2 is -1). Gauss's lemma therefore gives
irreducibility over k(s,t,u). Rescaling by the nonzero field element s
preserves that conclusion. The root in (4) embeds this degree-three field
in the source field, and (5) recovers every source coordinate, proving
equality. These are generic function-field assertions, not claims that
s or E are nonzero at every affine point of the original map.

The monic cubic gives a unique basis 1,W,W^2 over k(s,t,u). For a bounded
source polynomial P, substitute (5) and write its normal form as

    P=r0(s,t,u)+r1(s,t,u)*W+r2(s,t,u)*W^2.

Then P lies in the entire target field if and only if r1=r2=0. This is the
finite linear membership operator used below. It bounds the desired
source degree, never the unknown target function.

## Exact specialization minors, with no heuristic extrapolation

Source polynomials of degree<=7 and degree<=8 have dimensions 120 and165.
The corresponding membership operators have two rational-function output
coordinates. Evaluate those coordinates at finitely many integer target
points where s and multiplication by E are invertible. This gives exact
rational linear functionals on each finite source space.

The producer selects independent rows using three-coefficient polynomial
arithmetic modulo the monic cubic. The independent checker uses the
companion matrix

    C = [ 0  0  2*s^2*u ]
        [ 1  0    -s*t ]
        [ 0  1       1 ],

verifies its cubic, forms E(C)=3*C^2-2*C+s*t*I, inverts it by adjugates,
and reconstructs all three source multiplication matrices from (5).
It evaluates every source monomial as a matrix product. The coefficients
of W,W^2 in the first column reproduce the two membership functionals.
It also checks H evaluated on those matrices equals (s,t,u)*I.

The selected square minors have the following independently computed data:

| Source degree | Columns | Selected rank | Target points | Determinant modulo p |
| --- | ---: | ---: | ---: | ---: |
| <=7 | 120 | 116 | 58 | 866695135 |
| <=8 | 165 | 158 | 79 | 224657276 |

Here p=1000000007 is checked prime by trial division. Every fixed
specialization denominator is nonzero modulo p. Thus each modular matrix
is the reduction of the corresponding exact rational matrix, and a nonzero
modular determinant proves its rational determinant nonzero. The rational
minor remains nonzero over every characteristic-zero extension. Producer
row selection is not trusted: the checker reconstructs all entries using a
different representation and performs independent column elimination.

These are certificates of exact linear rank, not candidate-map sampling or
an assumption that the tested points are generic. An unknown target
function's denominator may vanish at a selected point; this is irrelevant.
Only the bounded source polynomial's normal-form coefficients are evaluated.
Their fixed denominators come from (5). If the unknown function belongs to
the base field, its two nonconstant normal-form coefficients are identically
zero, so those functionals vanish without ever evaluating the unknown
function itself.

The four independent functions 1,H1,H2,H3 lie in the degree-seven kernel.
The rank116 minor leaves dimension at most 120-116=4, proving (1).
For degree eight, Q(H),S(H),T(H) have distinct single leading monomials

    a^6*c^2,   -36*a^3*b^4*c,   36*a^4*b^2*c^2,

so they add three independent functions to the previous four. The rank158
minor leaves dimension at most 165-158=7, proving (2). As pullback embeds
rational function fields, equality of a rational function's pullback with
one of these polynomial combinations forces equality in the target field
itself. This explains why the arbitrary rational target becomes polynomial
in this low source-degree range, without a separate denominator theorem.

## Complete constant-Jacobian locus in degree eight

By (2), each component of a rational R with polynomial pullback degree<=8
is a linear combination of 1,Y1,Y2,Y3,Q,S,T. The last three have zero
constant and linear parts. If R(H) is Keller, the now-polynomial R has
nonzero constant Jacobian: injectivity of pullback forces its Jacobian to
be that constant, and its tangent matrix A is invertible. Subtract its
constant and multiply outputs by A^-1. The full normalized model is

    R_i=Y_i+u_i*Q+v_i*S+w_i*T,   i=0,1,2,

with nine free coefficients. No parameter is required nonzero. Its
Jacobian-minus-one has exactly27 source-monomial coefficient equations.
Seven explicit rational membership identities in this full ideal force

    u0=v0=v1=v2=w0=w1=w2=0.

The independent checker rebuilds the determinant from the nine-parameter
maps, reconstructs every generator, and multiplies all seven identities
exactly: **17 multiplier monomials, maximum multiplier degree one**. Thus
this is a full ideal consequence, with no saturation or omitted zero
stratum. The survivors are u1=lambda,u2=mu and give precisely (3).
The checker also verifies its generic two-sided inverse and determinant
one. The initial primary basis and subsequent lift logs are retained;
classification relies on exact replay, not a basis assertion alone.

## Sharp degree and fixed-source support consequences

Equation (1) makes eight the least source degree of a polynomial pullback
of any non-affine rational target function. The scalar Y1^2 attains it.
The triangular map (Y1,Y2+Y1^2,Y3) also shows eight is the least possible
degree of a genuine nonlinear target Keller modification of H.

For a nonzero surviving vector v=(0,lambda,mu), the polynomial H1^2 is

    a^2-6*a^3*b-2*a^4*c+9*a^4*b^2+6*a^5*b*c+a^6*c^2.

These six source monomials are disjoint from all thirteen original
nonlinear source columns of H. Under any invertible output matrix A,
each original column stays nonzero, and the six new columns are nonzero
multiples of A*v. There are therefore at least nineteen distinct nonlinear
source columns and hence at least nineteen occurrences. The displayed
single-row shear attains19 and is independently checked with determinant
one and the saved distinct collision points. There is no claim of a new
global upper bound: degree8/19 is dominated by the original degree7/13.
The nineteen lower bound uses the original source coordinates explicitly.

For degree<=7 the conclusion extends to every affine source/target frame.
If H'=B composed H composed A and q(H') is a polynomial of degree<=7,
compose that polynomial with A^-1. Ordinary degree is unchanged, and
(q composed B)(H) has degree<=7. Equation (1) makes q composed B affine,
hence q affine. For vector targets the same argument applies componentwise;
a Keller pullback then forces the affine tangent matrix invertible.
KCUB41 supplies raw13 throughout the resulting affine orbit. Arbitrarily
many target-only stages are covered, even if intermediate degrees grow,
because their total rational target functions are covered by (1).

## Reproduction, failed shortcut and open next gate

Run `python verify_raw13_target_field.py` from this directory. It checks all
polynomial identities (4), both complete field-rank certificates, the seven
known kernel functions and full determinant/membership/inverse model.
[Artifact index](raw13-target-field/README.md) and
[validation](raw13-target-field/validation.json) pin all sources and inputs.
The producers `raw13_target_field.py` and `raw13_target_field8.py` recreate
the finite specialization certificates; they use no polynomial ideal
solver. `raw13_target_locus.py` recreates the retained27-generator input
and Singular scripts without running the solver. Existing logs suffice
for independent certificate verification.

A naive leading-degree argument is invalid: the highest forms satisfy a
cubic relation, and subsequent cancellations produce S(H),T(H) of degree
eight. These discoveries and a further degree-ten expression are retained
in [discovery.json](raw13-target-field/discovery.json). A SAGBI completion
was considered and its library documentation inspected, but no completion
query was started and no finite-basis assertion is used. The full cubic
field operator provides the uncapped argument instead. The new small
nine-variable primary basis and membership jobs both completed; there are
no new solver timeouts or candidate coefficient grids.

Generic irreducibility, the field-normal-form criterion, specialization
lifting, dimension arguments and affine-frame consequences are written,
not formalized or externally reviewed. All finite polynomial and matrix
premises are independently verified. KCUB41 is imported only for the
complete affine-orbit count. The separate twelve-variable cubic retains
KCUB26's own scope; no new whole-field theorem is transferred to it.

The next concrete raw-term gate is affine source sparsification of the
explicit degree-eight family H+(0,lambda,mu)*H1^2. The fixed-source minimum
nineteen does not settle that joint family. More general nonlinear source
changes, source/target coordination and changed constructions remain
open routes to the original global objective. Another target-only search
ending at degree seven, or a fixed-source target search through degree
eight, cannot improve raw13.
