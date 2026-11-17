# Every nonzero positive-axis quadratic source shear has no Keller repair through degree eight (KCUB44)

Let H be the saved three-variable degree-seven map. Over any
characteristic-zero field k, put

    q(b,c)=n1*b^2+n2*b*c+n3*c^2,
    E_q(a,b,c)=(a+q(b,c),b,c),   F_q=H composed E_q.

For every nonzero quadratic q,

    k(F_q) intersect k[a,b,c]_(degree<=8)
        = span_k(1,(F_q)_1).                          (1)

The same equality holds through degree seven. Every rational target
function with polynomial pullback of degree<=8 is therefore affine in
the first target coordinate alone. A vector of such functions has Jacobian
rank at most one, so **no rational target repair can produce a polynomial
Keller map of degree at most eight**. There is no target-degree or denominator
bound. Arbitrary affine source frames after E_q, affine target frames, and
lower affine terms in q are covered as explained below.

This is the full three-coefficient family, including mixed quadratics and
all coefficient-zero strata. The zero quadratic is retained separately:
it is the old map H, with the larger low-degree field intersections from
[KCUB42](RAW13_TARGET_FUNCTION_FIELD_RIGIDITY.md), and is not excluded.

For nonzero q, the least possible degree of a Keller/collision endpoint
in this repair class is bounded by **[9,10]**. Nine is not claimed attained.
The source-only map F_q is a determinant-one collision of degree ten. A
checked q=b^2 control has row degrees (7,9,10) and nonlinear counts
(8,14,15), totaling37; this is not a global improvement over degree7/raw13.

The new degree-nine base-field calculation below concerns k(H), not
k(F_q). It does **not** classify degree-nine polynomials or repairs in the
transformed field. Global positive bounds and unrestricted 62/12 profile
counts remain unchanged; the overall simplicity goal is open.

## A new gap in the original target field

Write V_d=k[a,b,c]_(degree<=d), and let K=k(H1,H2,H3). Set

    Q=Y1^2,
    S=Y2^3+27*Y1*Y3^2-9*Y2*Y3,
    T=Y2^2-4*Y1*S,
    B=span_k(1,H1,H2,H3,Q(H),S(H),T(H)).

KCUB42 proves K intersect V8=B. The new exact rank certificate proves

    K intersect V9=B subset V8.                       (2)

There are220 source monomials through degree nine. Using KCUB42's full
cubic source-field presentation, the two nonconstant field-basis
coefficients give a linear membership operator on those220 columns.
A213-square exact rational specialization minor is nonzero: its independently
computed determinant is116435059 modulo the checked prime1000000007.
It uses107 good integer target points, yielding214 available functionals;
only213 selected rows enter the minor. All fixed denominators are nonzero
modulo the prime. Thus the rational minor is nonzero and persists over
every characteristic-zero field.

The checker reconstructs every entry using companion matrices and
adjugates, independently of the producer's polynomial triples and row
selection, then uses column elimination. The seven known functions in B
are independent: four have independent constant/linear jets, and the
three new degree-eight leading terms are

    a^6*c^2,  -36*a^3*b^4*c,  36*a^4*b^2*c^2.

Hence rank213 leaves precisely seven kernel dimensions, proving (2).
The exact-specialization argument and unknown-denominator boundary are
the same as in KCUB42: only normal forms of bounded source polynomials
are evaluated, not the unknown rational target functions themselves.

In particular, every source polynomial in K of degree<=9 actually has
degree<=8. For the original map H, its entire degree-nine rational-target
class is therefore the same as the degree-eight class of KCUB42/43.
This extends that original target-only/affine-source raw13 closure through
an endpoint degree of nine. It is not the transformed-field assertion (1).

## A formal obstruction for every positive-order quadratic deformation

Consider any finite polynomial family

    q_tau(b,c)=sum_(e>=d) tau^e*q_e(b,c),

where d>0, q_d!=0 and every q_e is homogeneous quadratic. Suppose
P_tau is polynomial in tau with source degree<=8 coefficient by coefficient,
and P_tau belongs to k(tau)(H composed E_(q_tau)). Define

    Z_tau=P_tau composed E_(-q_tau)=sum_j tau^j*Z_j.

It is a polynomial in tau and the source variables, and it lies in K(tau).
Every coefficient Z_j lies in K. To justify this without specializing an
unknown target denominator, use KCUB42's source-field basis 1,W,W^2 over K.
Expand each coefficient in that basis. The two nonconstant coordinates
of Z_tau are polynomials in tau over K. They vanish because Z_tau lies in
K(tau), and basis independence persists after adjoining a transcendental
variable. Comparing tau coefficients proves the assertion. No regularity
of the unknown target function at tau=0 is assumed.

Taylor expansion in a is exact because q_tau is independent of a:

    P_tau=Z_tau composed E_(q_tau)
         =sum_(m>=0) q_tau^m/m! * partial_a^m(Z_tau).   (3)

Let the coefficients of Q(H),S(H),T(H) in a function of B be its nonlinear
triple. The only degree-nine leading terms needed in (3), with q=q_d, are

| Source term | Degree-nine leading polynomial | a exponent |
| --- | --- | ---: |
|q*partial_a(Q(H))|6*a^5*c^2*q|5|
|q*partial_a(S(H))|-108*a^2*b^4*c*q|2|
|q*partial_a(T(H))|144*a^3*b^2*c^2*q|3|
|q^2/2*partial_a^2(H3)|24*a*b^3*c*q^2|1|
|q^3/6*partial_a^3(H2)|12*b^2*c*q^3|0|

They are nonzero and linearly independent for every nonzero q(b,c), since
q has no a dependence and their a exponents differ. The coefficient
constants are nonzero in characteristic zero. The independent checker
reconstructs all five terms from the full source polynomials with all
three quadratic coefficients formal.

For orders j<d, Z_j=P_j belongs to B by (2). Induction on the next orders
in (3) gives the following facts.

1. For d<=j<2d, only first Taylor derivatives occur. Earlier coefficients
   Z_i lie in B, so the source degree of Z_j is at most nine; (2) puts it
   back in B and forbids a degree-nine part. The first three distinct
   heads force the nonlinear triple of Z_(j-d) to vanish. Higher q_e
   multiply earlier triples already killed. Thus all Z_j through order
   2d-1 lie in B, and nonlinear triples of Z_i vanish for i<d.
2. For 2d<=j<3d, second derivatives can occur, but their indices are i<d,
   where the nonlinear triples are already zero. Their source degree is
   at most nine. Again Z_j lies in B. Besides the first three heads, the
   fourth forces the H3 coefficient of Z_(j-2d) to vanish. The same
   induction kills nonlinear triples for i<2d and H3 coefficients for i<d.
3. At order3d, a third derivative can only come from Z0 with three factors
   q_d. Its nonlinear triple and H3 coefficient are already zero. The
   only new degree-nine head is the last table row, from its H2 coefficient.
   All terms still have degree at most nine, so Z_(3d) lies in B and has
   no degree-nine part. The five different a exponents force that H2
   coefficient to vanish as well.

For precision, the possible second-derivative leading H3 term at order3d
comes from Z_d; first-derivative nonlinear terms come from Z_(2d). The
remaining lower indices have been killed by the preceding induction, so
no cancellation at the same a exponent is omitted. Third derivatives of
H1 have source degree at most seven after multiplication by q_d^3.
No fourth Taylor term can occur before order4d.

Consequently

    P_0=Z_0 belongs to span_k(1,H1).                   (4)

This is a written induction valid for every P_tau and every positive-order
quadratic deformation, not an inference from a finite coefficient grid.
It uses (2) only where the induction has established source degree<=9;
no statement about all higher-degree Z_j is assumed.

## Weighted degeneration covers every nonzero coefficient triple

The saved map has weights (1,-1,-2), so H composed D_tau=D_tau composed H,
where D_tau=diag(tau,tau^-1,tau^-2). For every fixed n=(n1,n2,n3),

    E_(n(tau))=D_tau composed E_n composed D_tau^-1,
    n(tau)=(tau^3*n1,tau^4*n2,tau^5*n3),
    F_(n(tau))=D_tau composed F_n composed D_tau^-1.    (5)

Every monomial of the full three-parameter family is checked against these
augmented source/parameter weights. The first nonzero q_tau term has order
3 if n1!=0, order4 if n1=0,n2!=0, or order5 if only n3!=0. These three
cases exhaust all nonzero triples; no roots, normalization division or
nonvanishing discriminant of q is needed.

Let U_n=k(F_n) intersect V8. For any P in U_n, put

    P_tau=tau^m*P(D_tau^-1 X),

where m is its greatest source monomial weight i-j-2k. Its tau powers
m-(i-j-2k) are all nonnegative, so P_tau is polynomial in tau and every
source coefficient has degree<=8. By (5) it belongs to the deformed target
field. Its value at tau=0 is the greatest-weight initial form of P. The
formal result (4) therefore puts every such initial form in span(1,H1).

A finite-dimensional polynomial subspace and its associated graded space
for a monomial weight filtration have the same dimension: sum the
dimensions of the successive filtered quotients. Their initial forms
realize those quotients. Thus dim(U_n)<=2. Constants and H1 composed E_n
are two independent members, the latter of source degree seven for n!=0
because its top form is -q^3*c. They give equality, proving (1). They also
belong to V7, so the same two-dimensional intersection holds there.

This argument turns every fixed nonzero triple into a positive-order
formal deformation. It does not infer universal rank from the seven
fixed-parameter pilots retained in discovery. The zero triple has no
positive leading order and is handled by its different base-field result.

## Repair endpoints, affine boundaries and controls

If R is any three-component rational target map and R(F_q) is polynomial
of degree<=8, (1) makes every component a constant plus a scalar multiple
of (F_q)_1. Its Jacobian rank is at most one, so it cannot be Keller.
Pullback embeds target function fields, so the target functions themselves
are affine in Y1 alone. Neither polynomiality nor a degree/denominator
bound on R was assumed.

An invertible affine source map A after E_q does not change degree:
compose a proposed polynomial endpoint with A^-1 and apply (1). An affine
target frame can be folded into R. An arbitrary polynomial q of degree two
with nonzero homogeneous quadratic part is covered too: its lower affine
part is a source affine shear fixing b,c, which factors on the right of
E_(q2). These are affine frames **after the displayed source shear**.
We do not classify H composed A composed E_q for an arbitrary affine A
inserted between H and E_q; that can change the nonlinear axis relative
to the displayed weights and lies outside this proof.

The full source map has top forms

    (F_q)_1,7=-q^3*c,
    (F_q)_2,9=12*q^3*b^2*c,
    (F_q)_3,10=8*q^3*b^3*c.

They are nonzero for every nonzero q. The source shear has polynomial
inverse E_(-q) and determinant one, so F_q is Keller of degree ten.
The original collision transports to

    p_n=(-n1+8*n2-64*n3,1,-8),
    q_n=(-1-n1-4*n2-16*n3,1,4),

both mapping to (0,1,0). Their last coordinates differ, so they remain
distinct for every triple in characteristic zero. The checker verifies
the all-parameter inverse transport and replays the full determinant and
collision of q=b^2 with counts (8,14,15). This gives upper10 against the
proved lower9 for the least repair/collision endpoint degree. Existence
at degree nine is unproved, and no global record is improved.

## Reproduction and remaining work

Run `python raw13_target_field9.py`, `python positive_source_quadratic.py`
and `python verify_positive_source_quadratic.py` from this directory.
The first produces the220-column original-field rank certificate; the
second retains the complete nonlinear source family and Taylor heads.
The independent checker rebuilds the213-square minor via companion
matrices/adjugates and column elimination, verifies the full source
composition/inverse, all weight identities, all five heads, bookkeeping
for leading orders3,4,5 through3d, the zero case and collision controls.
See [artifact index](positive-source-quadratic/README.md) and
[validation](positive-source-quadratic/validation.json).

The seven early fixed-triple rank pilots are retained as discovery only.
They suggested the dimension drop through degree seven; they are not
universal certificates. No polynomial ideal solve or large parameter
search is used in the final proof. Coefficientwise field descent, the
formal induction and weighted-filtration dimension argument are written,
not formalized or externally reviewed. Exact finite premises are
independently checked. KCUB42 supplies the generic field model and its
specialization boundary; KCUB43 supplies the prior zero-quadratic affine
baseline, with their existing trust boundaries.

The positive-axis quadratic gate through degree eight is closed. A
concrete next low-degree route is the negative-axis source shear
(a,b,c)->(a,b+r*a^2+s*a*c+h*c^2,c) with full rational target repair.
Its parameter weights are (-3,0,3), so this positive-weight degeneration
cannot be transferred to it. In particular the resonant s=-1/3 branch
already contains a saved degree9/19-term control, which must be read and
preserved rather than rerun as a new candidate. Every mixed and zero
coefficient stratum needs a genuinely new argument. Degree-nine repairs
in the positive-axis family also remain unclassified. The global
simplicity objective, other constructions and all unrestricted minima
remain open.
