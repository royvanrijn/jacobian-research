# One zero diagonal weight: all-degree classification and sharp dimension four (KCUB39)

Subsequent [KCUB40](EQUAL_NEGATIVE_FULL_DEGREE7.md) closes full111 through
degree seven by a changed curve/Euler-Wronskian proof. The current
degree<=7 diagonal-weight raw<=12 frontier is empty:only112 can collide,
with exactly13 terms. All larger frontier and full111-unresolved descriptions
below are historical stages;original timeout records are unchanged.

Over every characteristic-zero field, every three-variable polynomial Keller
map equivariant for a nontrivial diagonal G_m action with exactly one zero
weight is a polynomial automorphism. There is no degree or support bound.
For primitive opposite-sign weights (0,p,-q), p,q positive and coprime,
output translation and diagonal normalization give exactly

    F(x,y,z)=(x+g(y^q*z^p),y,z),   g(0)=0.

The inverse subtracts the same g. Same-sign weights also give explicit
polynomial inverses below. In this additional form, allowing dimension and
degree to vary, the minimum collision dimension is **four**: dimensions at
most three are excluded here, and identity stabilization of the established
three-variable map supplies a directly checked four-variable upper control.
It is an existing construction, not a new witness.

Combining this theorem with [KCUB38](DISTINCT_NONZERO_WEIGHT_MINIMUM.md), the
only remaining weight type for a degree<=7 three-variable diagonal-weight
collision with at most twelve raw nonlinear occurrences is **(1,-1,-1)**,
up to permutation and simultaneous sign reversal. Counts refer to the
coordinates in which the action is diagonal. The full (1,-1,-1) family,
including its degree-seven 77* sector, is not closed here.

A separate degree/profile consequence is that a collision in the strict
invariant degree-flag census below 764 cannot have any nontrivial linear
G_m symmetry. That statement permits linear diagonalization and uses only
the strict part of [KCUB23](../equal-negative-completion/PROOF.md). It is
not an unrestricted exclusion of any of the twelve remaining whole flags.
Global positive bounds, normalized support interval [10,48], raw13,
cubic12/36, cubic-linear89 and whole-profile counts 62 excluded/12 remaining
are unchanged. All global simplicity optima beyond the settled degree-three
minimum across dimensions remain open.

## Complete opposite-sign modules and boundary

Work first with det JF=1 and tangent map I. This normalization is legitimate:
F(0) is fixed by the torus, and its derivative is invertible and commutes
with the action. For weights (0,p,-q) the three weights are distinct, so
subtract the constant in the first output and rescale outputs diagonally.
For a nonzero original constant determinant this reduces to the stated case
without changing equivariance or invertibility.

Set t=y^q*z^p. The integer weight equation, using gcd(p,q)=1, gives the
complete modules in every degree:

    F=(A(x,t), y*B(x,t), z*C(x,t)).

For example p*b-q*c=p forces b=1+q*j,c=p*j,j>=0; the other two outputs
follow the same way. The map K[x,t] -> K[x,y,z] is injective, so coefficient
identities in x,t are genuine polynomial identities. Chain rule gives

    J = A_x*(B*C+t*(q*C*B_t+p*B*C_t))
        - A_t*t*(q*C*B_x+p*B*C_x) = 1.                 (1)

At t=0 this is A_x(x,0)*B(x,0)*C(x,0)=1. Every factor is a unit in K[x].
Tangent and constant normalization therefore force

    A(x,0)=x,   B(x,0)=C(x,0)=1.                      (2)

This argument removes every pure-x nonlinear coefficient, without dividing
by an unknown parameter or discarding any zero branch.

## Universal positive-slope obstruction

This section is a finite-support proof in arbitrary degree, not an inference
from the degree-seven models. Expand the remaining nonzero terms in x,t.
For an A-term x^k*t^j put m=k-1; for B- or C-terms put m=k.
Every such term has j>=1 by (2), with k>=0. Let

    R=max({0} union {m/j : actual nonzero remaining terms}).

Suppose R>0, and write R=r/s with coprime positive integers r,s. Give the
source variables auxiliary integer weights

    omega=(q*s,-r,0).

The normalized source monomial relation of any remaining output term is
(m,q*j,p*j). Its omega-weight is q*(s*m-r*j)<=0. Consequently every F_i
has omega-degree at most omega_i, with equality supplied by its linear term.
Let G_i be the sum of the terms of weight exactly omega_i in F_i.

Each derivative entry dF_i/dx_j has degree at most omega_i-omega_j.
In every determinant summand these bounds sum to zero. A product containing
any strictly lower term has weight strictly negative. Thus the weight-zero
part of det JF is precisely det JG. Since det JF=1, we have det JG=1.
This argument applies to integer weights of either sign; every polynomial
has finite support, so no completion or limiting argument is needed.

For a nonlinear top term, s*m=r*j. Coprimality gives j=s*l,m=r*l with
an integer l>=1. Therefore

    G_i=x_i*f_i(Z),   Z=x^r*y^(q*s)*z^(p*s),
    f_i in K[Z],   f_i(0)=1.

At least one f_i is nonconstant, since R was an attained positive maximum.
All entries v=(r,q*s,p*s) are strictly positive. The rank-one determinant
identity, obtained by expanding the diagonal matrix plus the logarithmic
rank-one derivative term, is

    det JG = product_i f_i
             + Z*sum_i v_i*f_i'*product_(j!=i) f_j.   (3)

It is a polynomial identity; a derivation in the fraction field clears to
this identity and requires no nonvanishing evaluation of x,y,z. If d_i is
the degree of f_i, its top coefficient is product_i lc(f_i) multiplied by

    1+sum_i v_i*d_i.

This positive integer is nonzero in characteristic zero. If any d_i>0,
(3) cannot be the constant 1. This contradiction proves R=0.

It follows that A=x*f(t)+g(t), B=B(t), C=C(t), with f(0)=B(0)=C(0)=1
and g(0)=0: all A exponents have k<=1 and all B,C exponents k=0.
Equation (1) becomes

    f(t)*(B*C+t*(q*C*B'+p*B*C'))=1.                  (4)

Both factors are polynomial units, hence f=1 and the bracket equals 1.
Its leading scalar is 1+q*deg(B)+p*deg(C), which is again nonzero in
characteristic zero. Thus B=C=1, completing the exact shear classification.
All identities and inverses are over the original field; no roots or
algebraic-closure parameter choices occur in this proof.

## Same-sign weights and smaller dimensions

After weight reversal and permutation write same-sign weights as (0,p,q),
0<p<=q. If p<q the complete modules are

    F=(A(x), y*B(x), z*C(x)+E(x)*y^m),

where the E term occurs only if q=m*p for an integer m>=2. The determinant
is A'*B*C. Units and tangent normalization give A=x,B=C=1. Subtraction
of E(x)*y^m is an explicit inverse. If p=q, the moving pair is M(x)*(y,z)
for a two-by-two polynomial matrix M. Now J=A'*det(M)=1 forces A=x and
det(M)=1 after normalization. Its adjugate gives a polynomial inverse.
These are complete modules because positive weights preclude nonlinear
terms of the same weight involving a higher-weight variable.

In dimension two, a nontrivial diagonal action with exactly one zero weight
has complete form (A(x),y*B(x)); its determinant A'*B is a nonzero constant,
so the same unit argument gives an automorphism. In dimension one, exactly
one zero weight is the trivial action, and there is no eligible nontrivial
action; independently a constant nonzero derivative gives an affine map.
Together these facts exclude collision dimensions at most three in the
stated exactly-one-zero form, in all degrees.

For the matching upper control set v=1+2*a*b and use the saved raw13 map

    H1=a-3*a^2*b-a^3*c,
    H2=3*v^2*a*c+(9*v^2-6*v-2)*b,
    H3=v^3*c+(6*v^2+2*v)*b^2.

The four-variable identity stabilization (x,H1,H2,H3) has weights
(0,1,-1,-2), degree seven, identity tangent map and nonlinear row counts
(0,2,5,6). Its full four-by-four determinant is one. Distinct points
(0,0,1,-8) and (0,-1,1,4) both map to (0,0,1,0). The independent checker
expands the determinant over all 24 permutations, checks every coefficient
against the retained [raw13 input](raw13-specialization/input.json), and
checks weights, counts, degrees and the collision exactly. This establishes
the sharp additional-form dimension four, without asserting a dimension
optimum for any other form.

## Degree-seven frontier and strict linear-symmetry corollary

KCUB38's ten-type raw<=12 frontier contains full111 and nine opposite-sign
one-zero signatures. The theorem removes all nine simultaneously, including
all their parameter-zero and high-weight cases. The remaining raw<=12 type
is full111. The restriction to displayed diagonal coordinates is essential:
a general linear source change can change monomial counts.

For collisions without a count cap, KCUB36 and KCUB38 leave only111 or112
through degree seven; same-sign and two-zero signatures are covered by
their complete unit/rank-one arguments. KCUB15's full112 collision
classification gives the known representative up to diagonal transforms,
with invariant output degree flag 764: the degree-seven part occurs in
only the third output, the degree-six part of the remaining quotient in
only the second, and the surviving first output has degree four. KCUB23
excludes111 strictly below764, including all its parameter branches in
that scope. Therefore a strict lower-flag collision has no nontrivial
linear G_m symmetry. Ordinary degree and the invariant flag survive
invertible linear source and output changes, so linear diagonalization is
allowed for this consequence. It is not allowed as a count-preserving step
in the raw sparsity statement.

Neither statement treats nonlinear torus actions through degree-changing
conjugacies. Neither closes full111's 77* sector. Its old fixed-linear and
broader chart-II forty-second timeouts, and the unreplayed chart-I primary
certificate, remain untouched and unused. No old timeout is relabelled as
a successful solver query. No new ideal or coefficient-grid query is used.

## Reproducible finite controls and trust boundary

Run `python verify_zero_weight.py` from this directory. The stdlib Fraction
checker independently reconstructs all nine complete degree-seven native
models produced by [zero_weight.py](zero_weight.py): 240 atoms and 438
native determinant coefficients. The t=0 boundary removes 162 pure-x
atoms, leaving 78. Across the nine models it checks all twelve possible
positive-slope principal-face specializations, including the filtered
Jacobian identity and the exact rank-one determinant expression. It also
checks 27 rank-one leading degree cases, 16 final-tail leading cases,
nine generic two-sided shear inverses, the same-sign block identities and
the four-dimensional upper control. See [the artifact index](zero-weight/README.md)
and [validation receipt](zero-weight/validation.json).

These finite checks validate algebraic premises and implementation; the
unbounded theorem is proved by the written attained-maximum, integer
filtration and polynomial-unit arguments above. They are not formalized or
externally reviewed. The frontier and strict corollary explicitly import
KCUB36/38 and KCUB23 with their retained boundaries. The direct all-degree
one-zero theorem uses no external plane theorem or archived radical
certificate. The earlier degree-four quotient analysis was read for
provenance and to avoid duplicate bounded work; it is not an extra premise.

Next work on a degree<=7 diagonal-weight raw improvement requires new
structure in full111, preserving its old unresolved solver records. General
source-frame sparsification, unweighted maps, higher-degree support and
changed carrier constructions remain separate open routes.
