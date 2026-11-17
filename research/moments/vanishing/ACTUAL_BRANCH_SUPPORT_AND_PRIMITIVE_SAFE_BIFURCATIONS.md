# Actual branch support and primitive safe bifurcations

16 November 2026. Written all-degree proof with exact algebraic controls.
Both FULL exceptional MZ cases remain **OPEN**. Independent replay,
formal verification and external review are not claimed; review is
deferred. No degree sweep or new holonomic CAS construction is used.
The [checker](verify_haar_branch_support_controls.py) and
[receipt](haar_branch_support_controls.json) verify the exact algebraic
controls, not the universal analytic-continuation or primitivity proofs.

The positive mechanism concerns ACTUAL marked periods, not all geometric
bifurcations. A pure fixed-marker failure must couple the logarithm at0
to another finite branch value of that SAME marker period. Conversely
at most one finite branch value for the entire actual period packet
forces rank1 and a unique safe centre. This does not prove that every
exceptional pure packet has the required branch restriction.

## 1. A single-valued infinity germ cannot have just one branch point

For any finite Gamma, put X=SL2/Gamma and use its original compact Haar
cycle. For a nonconstant regular f and ANY fixed original regular Q,

    C_Q(s)=h(Q/(s-f))=sum_(m>=0) h(Q f^m)/s^(m+1)

is a convergent, SINGLE-VALUED germ at infinity, of order O(1/s).
Its continuations are geometric regular-singular Cauchy periods on
the full regular base; keep all finite asymptotic bifurcations as well
as ordinary critical values. These are retained source properties in
[the resolvent audit](ALGEBRAIC_RESOLVENT_SPLITTING_AND_ZERO_VALUE_LOGARITHMS.md).

A branch point here means nontrivial local analytic continuation on
SOME branch of this actual multivalued function. A meromorphic pole
alone is not a branch point. Apparent singularities of a displayed
equation are not declared branch points either.

Suppose this function has at most ONE finite branch point a, allowing
other finite meromorphic poles. Its analytic continuation orbit factors
through pi1(C minus {a}): loops about all other finite points act
trivially on every branch. The large infinity loop is the inverse of
the remaining generator. It fixes the original germ by its convergent
Laurent series, so that generator fixes it too. The entire orbit has
one element. Thus C_Q is single-valued everywhere off finitely many
points. Regular-singular moderate growth makes it meromorphic at each
such point. It is rational.

Therefore EVERY NONRATIONAL original marked Cauchy function has at
least TWO distinct FINITE branch points. No purity assumption, support
bound, degree bound or proper-volume model is used.

For an eventual-pure failure, the retained zero-log theorem supplies
ONE fixed original Q with an actual integer-exponent logarithm at0.
That same C_Q is nonrational and must have another finite branch point
a!=0. The second point need not be a second nonzero intrinsic Tate
component: it could involve noninteger or genus monodromy. We do NOT
claim that Bad_boundary(f) must have two points.

## 2. A sufficient rank-one mechanism for the full original packet

Suppose the UNION of the finite branch points of ALL original C_Q has
size at most one. Section1 makes every C_Q rational, retaining all
actual relations. The full actual scalar Cauchy connection consequently
has rank1 over C(s): it is contained in C(s), and contains its literal
unit1 via s C_1-C_f=1.

The retained [actual safe-centre theorem](CONSTANT_GENUS_TRANSLATIONS_AND_ACTUAL_RECOMPOSITION.md),
Section5, gives a UNIQUE c with

    V_f = exp(c t) C[t],
    every original marked sequence for f-c eventually vanishes.

Here V_f is the full ACTUAL exponential image, not the abstract
integration module. One can also see the image statement directly:
rational Cauchy functions have finite principal parts at finitely many
values, so their Borel inverse images lie in the sum of the corresponding
exponential polynomial modules. Distinct exponentials give distinct
simple modules. Entire-image indecomposability retains only one sector;
a nonzero polynomial in that sector and differentiation generate all
of exp(c t)C[t]. The constant-unit normalization excludes the zero image.
This uses the GLOBAL polynomial-Weyl image and indecomposability, not
a localized connection or a claim that the unit is always essential.
For the relevant indecomposability, a formal map from the original Haar
source is determined by its initial functional on A/Lie(A)=C; the
recurrence u_Q'=u_(fQ) then determines every coefficient. Two nonzero
direct-summand projections into the original formal image would give
two independent such maps, contradicting that uniqueness. This is
the retained Haar/Peter-Weyl argument, without dropping any amplitudes.

For EVERY polynomial H, including every degree and all branches,

    H(f) eventually pure iff H(c)=0,
    and in that case H(f) is safe for ALL original markers. (1)

The rank-one source gate is global. The intrinsic unit extension splits;
this also follows by applying the intrinsic criterion to the safe
eventual-pure f-c and reflecting the affine translation back to f.
No scalar unit Cauchy function alone is substituted for the packet:
the known Long failure has rational C_1 but a nonrational marked C_Q.

In particular an ACTUAL counterexample needs at least two finite
branch values in its marked period packet, as well as the retained
capacity, primitive-character and automorphism constraints.

## 3. One geometric finite bifurcation is sufficient, but not necessary

If the FULL geometric bifurcation set of f has at most one finite
value a, all original marked periods are holomorphic continuations
off a. The argument of Section1 then makes every C_Q rational with
its only possible finite POLE at a. There must be a finite bifurcation:
with none, C_1 would be entire in s and vanish at infinity, contrary
to C_1(s)=1/s+O(1/s^2).

For that one value every original C_Q has an exact finite expansion

    C_Q(s)=sum_(j=1)^N_Q b_(Q,j)/(s-a)^j.

Equivalently, for every polynomial R,

    h(Q R(f))=sum_(j=1)^N_Q b_(Q,j) R^(j-1)(a)/(j-1)!.

Thus c=a in Section2, and the entire original marked functional is a
finite jet at that centre for each fixed Q. There is no uniform bound
on N_Q asserted here, and no claim that the original unit has jet
order0. In particular eventual purity is not silently replaced by
all-positive purity. This is an all-degree sufficient geometric
criterion, independent of whether the original volume extends to a
proper phase model.

The condition requires full geometric bifurcations, including finite
values from infinity. Checking only df=0 is insufficient. Nor does
regularity of a minimal connection imply regularity of all the marked
sections at a point: apparent poles of those sections must be retained.

## 4. Actual primitive pure SAFE phases with multiple geometric bifurcations

Use the full Klein presentations and normalizations of
[the quadratic quotient proof](QUADRATIC_TRACE_MODELS_AND_HAAR_RESIDUES.md).
Take its binary forms F_n,G_m in the TOP row (x,y), and define

    f_O=F_8(x,y)+G_12(x,y),
    f_I=F_12(x,y)+G_20(x,y).                              (2)

These are original regular invariant polynomials over Q. They have
two strictly positive weights, so ALL positive pure moments vanish
and every fixed original marker eventually vanishes. A sufficient
bound is n*l>deg_matrix Q for the l-th marked power. Their actual
exponential image is C[t], and C_1(s)=1/s exactly.

The TOP-row quotient is the punctured surface U_Gamma with equation

    C^2=K(F,G),
    K_O=F^3 G/36-3G^3,
    K_I=12F^5+G^3/144.

The map X->U_Gamma is an affine-line torsor (the opposite unipotent
version of the retained construction). Fibres of f over s are torsors
over the affine curves

    C^2=p_s(F)=K(F,s-F),

with the origin removed at s=0. The generic curves are geometrically
integral, so (2) are primitive on the ORIGINAL quotients. In fact
they are primitive upstairs too: the 2I root must be invariant and
has minimum degree12, ruling out a nonlinear root for degree20.
For2O the only possible nonlinear root has degree6 and sign character.
Differentiating f=H(q) along the unipotent fixing the top row gives
q in C[x,y]. The degree6 sign space is the one-dimensional 2T-fixed
space, with square proportional to G_12. A polynomial of that root
cannot also supply the nonzero F_8 term. All other root degrees are
excluded by the minimum invariant degree8. These use the retained
[primitive-base character theorem](PRIMITIVE_BASE_CHARACTERS_AND_FIBRE_EULER_CONSTRAINTS.md).

Exact discriminants, as POLYNOMIAL IDENTITIES in the entire base s, are

    disc_F p_s = -s^9(s-16)/6912                     on2O,
    disc_F p_s = s^10(50000s^2-1)/331776            on2I. (3)

At nonzero critical values the critical quotient points are

    (F,G,C)=(48,-32,0)                              on2O,
    (F,G,C)=(u,-3u/5,0), u^2=1/8000                on2I.

Their values are16 and +/-1/(100 sqrt(5)). Each is smooth on U_Gamma:
K_G is nonzero. On local coordinates F,C the Hessian determinant
of F+G is

    -2(K_FF-2K_FG+K_GG)/K_G^2,

which is nonzero at these points. Thus these are ordinary nodal
degenerations of the curve family. The top-degree coefficient of p_s
is constant and nonzero; the smooth projective hyperelliptic models
have fixed infinity punctures. Outside the values in (3) there is a
smooth proper family with those punctures removed, hence no additional
finite asymptotic bifurcation. The affine-line torsor is smoothly
locally trivial and preserves this conclusion for X.

At s=0 the normalization is given by y=C/F:

    y^2=F(108-F)/36                                 on2O,
    y^2=F(1728F^2-1)/144                            on2I.

Removing the origin removes one further point from each normalization.
The exact genera, punctures and Euler characteristics are:

| Group | Generic genus / infinity points | Generic chi | chi at0 | chi at each nonzero critical value |
|---|---|---:|---:|---:|
| 2O | 1 / 2 | -2 | -1 | -1 |
| 2I | 2 / 1 | -3 | -2 | -2 |

At a nonzero critical value one node identifies two points of the
normalization and lowers Euler characteristic by1. The difference at0
proves that0 is also a genuine finite bifurcation, although it is not
an ordinary critical value on the punctured base. Hence the FULL finite
bifurcation sets in these controls are EXACTLY

    {0,16} on2O,
    {0,+/-1/(100 sqrt(5))} on2I.                         (4)

All these extra bifurcations coexist with pure safety and actual image
C[t]. On smooth fibres the affine-line torsor retracts to an affine
curve, so H2 of the fibre is ZERO. The unit-relevant complement H3
source is therefore only its ambient rank-one line. The nontrivial
nodal monodromy lies in fibre H1, not in that source. These controls
do not contradict the intrinsic boundary criterion.

## 5. Consequence for the uniform research obligation

Purity does NOT force a primitive phase to have one finite geometric
bifurcation, or force all its fibre cohomology monodromy to vanish.
The explicit SAFE controls (2)--(4) rule out those stronger shortcuts
on BOTH original exceptional models. They are not counterexamples.

What remains is to split the unit in the actual unit-relevant source
on EVERY pure phase, or construct one actual pure phase and fixed
original bad marker. A failure marker with an integer logarithm at0
must also have another finite branch value; geometric critical values
which are invisible to that packet do not supply such a marker.
The improved capacity and automorphism inputs are retained. No
finite-degree calculation, additional decision procedure or external
review is substituted for full pure-to-safe inclusion.
