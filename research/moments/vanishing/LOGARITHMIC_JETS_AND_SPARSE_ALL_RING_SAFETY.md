# Logarithmic jets, sparse proper models and all-ring safety

16 November 2026. Internal written all-degree proof using retained Haar,
global-bundle and projective-period inputs. Independent replay, formal
verification and external review are not claimed. Both FULL exceptional
MZ cases remain **OPEN**. Review is deferred and is not the next task.

EVERY phase on SL2/2O supported on at most FOUR left-torus weights,
and EVERY phase on SL2/2I supported on at most SIX weights, satisfies
MZ in every degree and over EVERY commutative C-algebra. Zero need not
be in the support. Over C, every mixed-sign phase in these classes,
and EVERY nonconstant polynomial recomposition, is NOT eventually pure.
A pure phase has zero neutral part and only one sign of charged weights.

The proof is a general valuation inequality, not a degree/support search
or another auxiliary curve construction. Torus derivatives separate all
weight components. Their valuation loss is bounded by logarithmic
discrepancy. A transverse simple coefficient root supplies the extra
logarithmic coordinate that closes the sharp cases. A general charged-
degree/root-multiplicity criterion applies beyond the four/six thresholds.

See the [checker](verify_haar_sparse.py), [receipt](haar_sparse_verification.json)
and [source audit](haar_sparse_source_audit.json).

## 1. A logarithmic derivation loses at most log discrepancy

Work locally on a smooth threefold with coordinates R,xi,t, where R=0
is a boundary divisor and t=0 is another boundary when it meets the
chart; at the opposite fibre end replace t by 1/t. Put

    Lambda=dR/R wedge dxi wedge dt/t.

For a divisorial valuation v centred on this chart, let

    b=v(R), c=v(t) at the section, c=0 off the sections,
    A=A_Z(v), d=A-b-c.

Lambda has exact order d-1, and d>=0 because the boundary is SNC.
For either delta=t*partial_t or partial_xi, and EVERY rational function g,

    v(delta g)>=v(g)-d.                                (1)

Here is a direct proof. On a smooth model extracting v, write a local
uniformizer u and g=u^m h with h a unit. Pullbacks of dR/R and dt/t
are regular 1-forms plus a multiple of du/u; dxi is regular. In the
wedge dg wedge contraction_delta(Lambda), at most ONE du/u can occur.
Its order is therefore at least m-1, for negative m as well as positive.
Dividing by Lambda, of order d-1, proves (1). This also proves the
inequality for arbitrarily iterated applications, with loss j*d.
No chosen finite resolution or simple-critical-point assumption is used.

If A0(R,xi) restricts at R=0 to a root of multiplicity m, then
partial_xi^m A0 is a unit at the centre. Applying (1) repeatedly gives

    h=v(A0)<=m*d.                                     (2)

For m=0 the coefficient itself is a unit and h=0. Thus the base-root
estimate follows from the SAME derivation argument, without a new
inversion-of-adjunction theorem or an assumed normal crossing for A0.
It includes singular coefficient curves and all infinitely near valuations.

If m=1, A0 is smooth transverse to R=0. Coordinates R,A0,t give instead

    Lambda'=dR/R wedge dA0/A0 wedge dt/t,
    d'=A-b-h-c=d-h>=0.

The torus derivation delta=t*partial_t leaves A0 unchanged. The same
wedge argument now gives the STRONGER estimate

    v(delta^j g)>=v(g)-j*(d-h).                        (3)

This is where a simple root differs from merely substituting m=1 in (2).
The improvement is essential for the four/six thresholds below.

## 2. Vandermonde separation on the ORIGINAL global bundle

Use the retained [global-bundle model](GLOBAL_TORUS_BUNDLE_AND_UNIT_INVISIBLE_INFINITY_POLES.md):
PGL2 is the nonzero Hom-line torsor over Q2=(P1 times P1) minus diagonal B.
Its smooth P1-bundle completion Z has disjoint fibre sections D0,Dinf,
and the ORIGINAL Haar volume Theta has divisor

    div(Theta)=-2*pi^*B-D0-Dinf.

Let f=sum_(i in S) f_i have k=|S| DISTINCT effective torus charges;
the usual SL2 weights are 2i. Each f_i is nonzero. The neutral component
may contain both a constant and a nonconstant part; it is one weight.
For a charged component let 2D_i be its MINIMUM matrix-entry polynomial
representation degree, i.e. its filtration degree, and let F_i be its
nonzero leading binary form on B. In the row-line chart

    [x:y:z:w]=[t:t*xi:1:eta], R=eta-xi,

homogenizing lower terms by powers of the determinant gives EXACTLY

    f_i=t^i*A_i(R,xi)/R^D_i,
    A_i(0,xi)=F_i(1,xi) nonzero.                       (4)

Minimality of the representation degree ensures the last polynomial is
not identically zero: otherwise the homogenized numerator is divisible
by the determinant and its degree can be reduced. Its binary degree is
2D_i. For right-Gamma invariant f_i it is a right-Gamma invariant binary
form. At other base/fibre charts the same statement holds with local
units and the appropriate regular frame. Choose a chart around the
valuation centre; at Dinf replace t by 1/t and negate all charges.

The derivation delta=t*partial_t satisfies delta f_i=i f_i. Hence

    f_i=L_i(delta) f,
    L_i(T)=product_(j in S,j!=i) (T-j)/(i-j).            (5)

This has degree k-1, independently of charge magnitudes or gaps. If
v(f)>=0, (1) gives for EVERY component

    v(f_i)>=-(k-1)*d.                                 (6)

There is no genericity condition on coefficients or critical points.
In particular (5) keeps zero-free supports and arbitrary nonlinear
relations between their coefficients.

## 3. General charged-component capacity criterion

For a leading binary form F_i let M_i be its MAXIMUM projective root
multiplicity, and define its integer capacity

    cap_i=D_i                    if M_i<=1,
    cap_i=D_i-M_i                if M_i>=2.            (7)

Suppose f has BOTH signs and there is a positive component and a negative
component whose capacities are BOTH at least k. Then f has a proper
model carrying its SAME ORIGINAL holomorphic Haar volume.

To prove this, the global-bundle theorem already handles EVERY finite-
base boundary valuation over finite phase values. Only b>0 remains.
Write c for the one fibre-section order, or zero off both sections. A
negative order for the original volume would imply, by integrality,

    ord_v Theta=A-1-2b-c<0 => d=A-b-c<=b.              (8)

At D0 choose the specified NEGATIVE component, so i=-a<0. At Dinf
use the positive component and inverse frame. Off both sections choose
either specified component and put a*c=0. Formula (4) gives

    v(f_i)=h-D_i*b-a*c,
    (k-1)*d>=D_i*b+a*c-h.                             (9)

Let m be the multiplicity at the ACTUAL base point of F_i. If m=0,
h=0 and (9),(8) contradict k-1<D_i. If m>=2, (2) gives

    (k-1+m)*d>=D_i*b+a*c.                             (10)

Here k<=D_i-M_i and m<=M_i, so k-1+m<=D_i-1. Together with d<=b
this is impossible since b>0. If m=1, use (3) in (5) instead. Then

    (k-1)*(d-h)>=D_i*b+a*c-h,
    (k-1)*d>=D_i*b+a*c+(k-2)*h.                       (11)

The hypothesis cap_i>=k implies D_i>=k. Since k>=2, the last term is
nonnegative, and (8) makes (11) impossible. This also treats a simple
local root when some OTHER root of F_i is multiple.

Thus every finite-phase divisorial valuation has nonnegative original
volume order. The retained
[valuative original-model criterion](CONSTANT_VOLUME_MODELS_AND_ACTUAL_MARKER_DEFECTS.md)
constructs a smooth PROJECTIVE A1_f model containing the ORIGINAL open
PGL2 with the SAME holomorphic Theta. Finite etale quotient by the right
projective exceptional group gives the same property on SL2/Gamma by
the retained constant-volume correspondence criterion. Normalized compact
periods are the original unweighted Haar periods, with mass1.

The retained [projective original-volume theorem](SEMISTABLE_LOG_PERIODS_AND_SINGLE_ISOTYPIC_SAFETY.md)
therefore excludes eventual purity of f and of EVERY nonconstant H(f).
The same model works by composition with the finite polynomial map H.
Primitive purity is NOT reflected from a composite. This general
criterion is valid for any invariant phase on PGL2 with the stated
capacities, not just at the exceptional four/six thresholds.

## 4. Four weights on 2O and six on 2I, in EVERY degree

On 2O, D_i>=4 and every leading root has multiplicity at most floor(D_i/3).
On 2I, D_i>=6 and every leading root has multiplicity at most floor(D_i/6).
These are the retained minimum matrix degrees8/12 and right-projective
orbit sizes6/12. They apply to each original nonzero charged component,
including nonhomogeneous coefficients of arbitrary matrix degree.

On 2O: if M_i<=1, cap_i=D_i>=4. If M_i>=2, D_i>=6 and

    cap_i>=D_i-floor(D_i/3)>=ceil(2D_i/3)>=4.

On 2I: if M_i<=1, cap_i=D_i>=6. If M_i>=2, D_i>=12 and

    cap_i>=D_i-floor(D_i/6)>=ceil(5D_i/6)>=10>=6.

Consequently EVERY charged component has capacity at least4 on2O and
at least6 on2I. Section3 applies to EVERY mixed-sign phase with

    k<=4 on2O, k<=6 on2I.                              (12)

Every nonconstant H of each such mixed phase is NOT eventually pure.
This closes, in particular, ALL three-NONZERO-weight supports on both
quotients, as well as all four-weight supports on2O and all supports
through six weights on2I. The earlier zero-containing-three-weight
results remain valid; their auxiliary-model collision proofs are not
replaced or declared incorrect.

If f has only one charged sign, h(f^m)=h(p^m) for its neutral part p.
The retained [entire neutral radical theorem](VALUATIVE_BOUNDARY_VOLUMES_AND_TORUS_INVARIANT_CLOSURE.md)
then says eventual purity forces p=0. Conversely zero-neutral one-sided
phases are pure and safe, by weight selection against EVERY fixed
original marker. A purely neutral eventually-pure phase is zero.
Thus over C, within (12), the EXACT classification is

    f eventually pure iff p=0 and f has only one charged sign.       (13)

Zero is included. If a one-sided f has constant neutral part c, all
nonzero H(f) are pure precisely when H(c)=0, and those composites are
safe. If its neutral part is nonconstant, no nonzero H(f) is eventually
pure: neutral projection gives H(p), and the neutral radical is zero.
These assertions also include constant phases in the evident way.

## 5. Exact all-ring criterion and an explicit all-marker cutoff

Let R be ANY commutative C-algebra; R may be nonreduced or infinitely
generated. Work in R tensor_C A_Gamma with a fixed support satisfying
(12). In fixed finite C-bases of the phase components let

    I0=ideal_R(all neutral coefficients),
    I+=ideal_R(all positive-charge coefficients),
    I-=ideal_R(all negative-charge coefficients),
    J=I0+I+*I-.                                      (14)

These coefficient ideals are independent of the chosen bases. Then

    f eventually pure iff J is NILPOTENT.             (15)

In either case f is safe against EVERY original marker over R.
Necessity: let R0 be the finitely generated C-subalgebra of R generated
by the finitely many phase coefficients. Every closed point of SpecR0
specializes the actual eventual-purity identities to C. Classification
(13) makes I0 and I+*I- vanish at every such point. Nullstellensatz
puts J0 in the nilradical of this Noetherian algebra, so J0, and hence
its image J in R, is nilpotent. This does not require R itself to be
Noetherian, nor a coefficient-ring extension of Hodge theory.

For sufficiency and safety choose nu>=1 with J^nu=0. If charged
components occur, put

    w=min nonzero absolute phase charge,
    W=max absolute phase charge,
    B_Q=max absolute charge of the fixed marker Q.

In a term of Q*f^m let z be the number of neutral factors and n+,n-
the positive and negative counts. Its coefficient lies in
J^(z+min(n+,n-)). If it survives, z+min(n+,n-)<=nu-1.
For a weight-zero term with, say, n+>=n-, charge balance implies
w*n+<=W*n-+B_Q. Therefore

    m=z+n++n- <= (1+W/w)*(nu-1)+B_Q/w.                (16)

The opposite case is identical. Hence the explicit sufficient bound

    w*m>(W+w)*(nu-1)+B_Q => h_R(Q*f^m)=0.             (17)

This includes Q=1 and proves eventual purity as well. If there are no
charged factors, f^nu=0 already. No numerical nilpotence index is
claimed without a certificate J^nu=0.

One CANNOT replace (15) by requiring either I+ or I- nilpotent globally.
For R=C times C, the phase (U_positive,V_negative) has I+*I-=0 while
neither ideal is nilpotent. The paired-product ideal in (14) keeps the
correct componentwise orientation and all ordinary nilpotents.

## 6. Consequences for the FULL remaining target

Every genuine exceptional MZ failure must now have at least FIVE
distinct left-torus weights on2O and at least SEVEN on2I, after EVERY
algebraic automorphism and EVERY choice of maximal left torus. Haar
and all fixed markers transport exactly under those automorphisms.
This improves the former 'four weights OR three nonzero weights'
residual boundary. It is a necessary constraint, not a failure producer.

There is also an all-support degree-versus-support obstruction. Let
C+ and C- be the maximum capacities (7) on the two charged sides.
An eventually-pure mixed phase with k weights must satisfy

    min(C+,C-)<k.                                     (18)

Thus ONE entire charged side has cap_i<k for EVERY component. The orbit
bounds give cap_i>=2D_i/3 on2O and >=5D_i/6 on2I. On that side the
original charged matrix degrees consequently satisfy

    2D_i<3k on2O, 2D_i<(12/5)k on2I.                  (19)

This general obstruction applies for arbitrary support size, not only
(12). It identifies the low-degree charged side any full pure mixed
phase must retain; it does not bound the whole phase's degree.

No full unit-splitting bridge is inferred merely from finite-base volume
extension. Here ALL finite-phase original-volume poles are excluded,
and the projective-period theorem rules out the mixed phases themselves.
The full canonical boundary extension/corrected torus-component target
remains OPEN for larger supports failing the capacity test. Original
marked functions, actual global relations, all primitive branches and
all nonlinear recompositions remain in that target. Both FULL MZ cases
remain OPEN. Gaussian/formal Bernstein effectivity stays separately
owned; no bounded support/degree/parameter campaign was run.
