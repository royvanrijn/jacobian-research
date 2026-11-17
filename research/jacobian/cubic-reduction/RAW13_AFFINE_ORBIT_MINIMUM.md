# Every affine frame of the established map needs thirteen nonlinear terms (KCUB41)

Subsequent [KCUB43](QUADRATIC_TARGET_AFFINE_DEGREE8.md) closes the joint
affine first-nonlinear-target gate:all nonzero degree-eight shears need
at least16 terms in every affine frame. The interval[16,19] is not a sharp
minimum claim;only the zero shear can attain the complete-class minimum13.

Subsequent [KCUB42](RAW13_TARGET_FUNCTION_FIELD_RIGIDITY.md) closes all
rational target-only changes ending at degree<=7,with no target degree or
denominator cap,and classifies the first nonlinear degree-eight layer. Its
nonzero quadratic shears need19 terms with source coordinates fixed;affine
source frames of that new two-parameter family remain open.

Let H be the established degree-seven three-variable Keller map with raw
nonlinear row counts (2,5,6). Over every characteristic-zero field, every
map

    G(X)=N*H(M*X+t)+u,   M,N invertible,

has at least **thirteen nonlinear monomial occurrences**. This includes all
source translations and all independent source/output linear changes,
without requiring identity linear part. Output translation u has no effect
on this count. If determinant one is required, impose det(N)*det(M)=1;
the bound holds even without that extra condition.

Equality holds precisely when t=0 and both M and N are monomial matrices:
coordinate permutations and nonzero coordinate scalings. Thus a genuinely
mixed affine frame cannot improve or even retain thirteen. The original
map attains the bound. This is an orbit minimum and rigidity theorem, not
a new witness or a global unrestricted sparsity minimum.

Combining [KCUB40](EQUAL_NEGATIVE_FULL_DEGREE7.md), thirteen is also sharp
in **every affine coordinate frame** of every degree<=7 three-variable
collision map equivariant for a nontrivial displayed diagonal G_m action.
The earlier escape through general source affine coordinates is now closed
for that entire classified family. Unweighted maps, nonlinear coordinate
changes, other constructions and higher-degree sparsity remain open. The
global positive bounds and unrestricted 62/12 profile counts do not change.

## A summed homogeneous-rank lower bound

Write the homogeneous expansion at t as H(X+t)=sum_d H_{d,t}(X). Let r_d(t)
be the dimension of the span of its three degree-d scalar output polynomials.
Equivalently it is the rank of the three-row monomial coefficient matrix
C_d(t). Invertible source linear substitution is an invertible operator on
homogeneous polynomials, and invertible output mixing preserves their span.
Consequently r_d(t) is unchanged by M and N.

Every output frame has at least r_d(t) nonzero degree-d scalar polynomials,
and each contributes at least one monomial occurrence. Hence

    nonlinear occurrences of G >= sum_(d=2..7) r_d(t).       (1)

This is a universal linear-algebra argument, not a coefficient search. The
complete translated coefficient matrices give the exact table below, with
t=(alpha,beta,gamma):

| Centre | r2 | r3 | r4 | r5 | r6 | r7 | Sum |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| alpha=beta=0 | 2 | 3 | 3 | 2 | 2 | 1 | 13 |
| (alpha,beta)!=(0,0) | 3 | 3 | 3 | 2 | 2 | 1 | 14 |

All entries hold over every extension field, including all parameter-zero
strata. The next two sections prove the table with small explicit minors.
Thus (1) already proves the global affine-orbit bound thirteen and forces
an equality centre onto the original c-axis.

## High-degree minors and the alpha=0 branches

Source monomials below use the original variable names a,b,c; alpha,beta,
gamma are centre parameters. All determinants list columns in the displayed
order. The top output is exactly 8*a^3*b^3*c, giving r7=1. Rows one and two
are zero above their original degrees, giving the stated upper bounds.
For the lower bounds:

- Degree6, rows2,3 and columns (a^3*b^2*c,a^2*b^4): minor288.
- Degree5, rows2,3 and columns (a^2*b^3,a*b^4): minor1728*alpha.
  At alpha=0 use (a^2*b^3,a^2*b^2*c), with minor432.
- Degree4, all rows and columns (a^3*c,a*b^3,b^4): minor-1728*alpha^3.
  At alpha=0 use (a^3*c,a^2*b*c,a*b^3), with determinant-336.
- Degree3 at alpha=0, columns (a^2*b,a*b^2,a*b*c): determinant-432.
- Degree2 at alpha=0, columns (a^2,a*c,b^2): determinant-72*beta.
  If alpha=beta=0 the first row is zero and rows2,3 on (a*c,b^2)
  have minor24, so the rank is exactly two.

Every listed polynomial is reconstructed independently from the full
six-variable translation, rather than from the producer's binomial formula.
The remaining alpha!=0 degree2/3 branches are covered next.

## Two complete selected-minor units in the alpha!=0 chart

H has weights (1,-1,-2). For D_lambda=diag(lambda,lambda^-1,lambda^-2),
H(D_lambda X)=D_lambda H(X). Expanding at a centre shows

    H_{d,D_lambda*t}(D_lambda X)=D_lambda H_{d,t}(X),

so all r_d are constant on these centre orbits. If alpha!=0, choose
lambda=1/alpha and normalize the centre to (1,beta,gamma), replacing the
last two coordinates by alpha*beta,alpha^2*gamma. These new beta,gamma are
arbitrary; no zero condition on either is imposed. The checker verifies
weight covariance of every translated matrix coefficient.

For degree2 at this normalized centre put

    z=2*beta+1,
    P=6*beta^2+2*beta*gamma+beta+gamma.

The determinants on columns

    (a^2,a*c,b*c),   (a*b,a*c,b*c)

are m1=-54*z*P and m2=108*z*(2*P-1). Thus
m2+4*m1=-108*z. If the rank were at most two, z=0 and beta=-1/2.
But the determinants on (a^2,a*b,a*c) and (a^2,a*c,b^2) then equal
-27*(2*gamma-7)/2 and -36*(gamma-3). They cannot vanish simultaneously.
This proves r2=3 everywhere in this chart.

For degree3, the three determinants on columns

    (a^2*b,a*b*c,b^2*c),
    (a^2*b,b^3,b^2*c),
    (a^3,a^2*b,a^2*c)

are respectively

    -648*(2*beta+1)^2,
    288*(3*beta+gamma-1),
    144*beta^2*gamma^2.

If they all vanished, the first two force beta=-1/2,gamma=5/2, while
the third equals225. Therefore r3=3.

For direct certificate arithmetic, the producer also gives a rational
polynomial combination of the four degree2 determinants equal to one and
a combination of the three degree3 determinants equal to one. The two unit
identities use **31 multiplier monomials in total**, with maximum multiplier
degree recorded in [validation](raw13-affine-jet/validation.json). The
independent Fraction checker recomputes all seven determinants and multiplies
the units exactly. No primary ideal-membership answer is a proof premise.
Together with the alpha=0 branches, these arguments prove every rank in
the table, without a coordinate-entry bound or sampled centre grid.

## Equality forces the original centre and monomial frames

If the total is thirteen, the rank table forces alpha=beta=0, and every
degree must attain its rank bound individually. In particular there is
exactly one degree-seven occurrence. Under a source linear matrix M, the
top scalar polynomial becomes

    8*ell1(X)^3*ell2(X)^3*ell3(X),

where the three ell_i are independent linear forms. Output mixing
multiplies this polynomial by the nonzero vector N*e3. For the total number
of top-degree occurrences to be one, the scalar product itself must be
one monomial and the vector must have one nonzero coordinate. Unique
factorization then forces every ell_i to be a scalar coordinate variable;
independence makes the variables distinct. Thus M is a monomial matrix.
This conclusion uses the full top factorization, not a finite list of
source matrices.

A monomial source substitution merely permutes and rescales source columns.
At t=(0,0,gamma), gamma!=0, the translated map has **eighteen distinct
nonlinear source monomials**. The original thirteen remain nonzero; the
five extra columns are

    a^3, a^3*b^2, a*b, a^2*b^2, a^3*b^3,

with coefficients -gamma*e1,12*gamma*e2,6*gamma*e3,12*gamma*e3,
8*gamma*e3, respectively. Invertible output mixing cannot annihilate a
nonzero coefficient column. Therefore it needs at least eighteen
occurrences in this case, contradicting thirteen. Equality requires t=0.

At t=0 the original thirteen source columns are all pure nonzero multiples
of one of e1,e2,e3, and all three axes occur. After monomial source changes,
thirteen total occurrences require every output-transformed column to have
one coordinate. Thus all three columns of N have exactly one nonzero
entry; invertibility makes N monomial. Conversely any monomial M,N with
t=0 retain exactly thirteen nonlinear occurrences. Output translation u
remains arbitrary. This proves the complete equality statement.

For a genuinely mixed source frame on the c-axis the top degree already
needs at least two occurrences, so the total is at least fourteen. Off the
c-axis the rank sum itself is fourteen. These are lower bounds, not claims
that fourteen is attained. No next-value optimum is asserted.

## Scope, verification and provenance

[The retained input](raw13-specialization/input.json) supplies H. The
independent checker reconstructs its full determinant one, identity tangent,
degree seven, (2,5,6) counts and the distinct collision points (0,1,-8),
(-1,1,4) mapping to (0,1,0). It is the established witness, not a new map.

Run `python raw13_affine_jet.py`, then `python verify_raw13_affine_jet.py`
from this directory. The producer reconstructs six complete translated
jet matrices by binomial coefficients. The checker instead performs the
full substitution H(x+alpha,y+beta,z+gamma) in six polynomial variables,
collects all **348 matrix entries**, verifies covariance and every explicit
minor, multiplies both unit identities, and checks all eighteen c-axis
source columns and the thirteen pure-axis columns at the origin. See
[the artifact index](raw13-affine-jet/README.md).

Three preliminary, subsecond, two-variable SymPy Groebner calculations on
rank ideals returned unit ideals and are retained as discovery provenance.
They are not used as evidence for the theorem: the written field arguments
and independently multiplied selected-minor units replace them. Production
and verification make no solver calls. No old solver query or coefficient
grid is rerun. The universal rank invariance, weight normalization,
monomial-factor and equality proofs are written, not formalized or externally
reviewed; exact finite premises are independently checked.

KCUB40 classifies every degree<=7 diagonal-weight collision as the known112
representative up to diagonal operations and permutation. This theorem
therefore closes arbitrary affine sparsification throughout that class,
including source frames in which its action is no longer diagonal. It does
not imply a term lower bound for maps with no such affine-coordinate form,
or for nonlinear coordinate changes. The original global raw/normalized
support, cubic dimension/support and cubic-linear objectives remain open.
Future work must leave this complete affine orbit or use genuinely nonlinear
changes; another affine parameter grid cannot improve raw13.
