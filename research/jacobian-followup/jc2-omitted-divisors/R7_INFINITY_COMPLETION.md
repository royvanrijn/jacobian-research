# The complete infinity root ledger closes the single-prime r7 stratum

14 November 2026. PF2OMIT5 uses the full finite-normalization infinity
ledger, including unramified rows. Assume exactly the hypotheses of
[PF2R7INF1](../jc2-r7-universal-contact/PROOF.md): a complex degree-six
plane Keller map, whose entire affine-normalization ramification divisor
is one prime E of index two and residue degree one, with birational
polynomial image (p(t),q(t)) of degrees (6,10) and gap21. The statement
does not derive this entire-ramification hypothesis from a terminal r7
passport. Every omitted unramified prime remains allowed.

Under these hypotheses the stratum is empty, without importing source
finality or a conclusion of Borisov Theorem3.8. The earlier type1 residue,
common-boundary-point, finite coexceptional-curve and positivity inputs
remain those explicitly isolated in PF2R7INF1. This is a written geometric
proof with exact algebraic controls, not formal verification or external
review. Singleton-zero properness for unrestricted Keller maps remains open.

## 1. Retained geometry and the missing generic-root constraint

On the minimal smooth resolution Y of the finite projective normalization
W over target P2, put H=g^*O(1). The retained adjoint pencil F=K+H has
connected rational fibres and H.F=2. Its map to P1 is basepoint free.
The affine normalization N is smooth and is unchanged by Y over A2.
The two primes over the generic target infinity line have residue degree
one. After the existing exclusions, their indices are (1,5). Call them
C1,C2, with index five on the vertical C2.

The remaining cases from PF2R7V1 and PF2R7G1 are:

| Case | F.C1,F.C2 | Degree of E over the pencil base |
|---|---|---:|
| Reciprocal conics | 1,0 | 3 |
| Quartic contacts | 2,0 | 4 |

A general pencil fibre maps isomorphically onto a target conic. Its
primitive conic equation G(z;x,y) generates the degree-six function field
over C(x,y), so its degree in z is exactly six. Its coefficients in z
have target degree at most two. Homogenize it to degree six in the pencil
coordinate and degree two in the target coordinates.

**Generic infinity cluster lemma.** At the generic target infinity DVR,
this binary root order is finite flat: its quadratic target part is
nonzero, hence its coefficients are not all zero there. Its normalization
is the same function-field normalization as W over that DVR. The z-value
on C1 varies with target direction; the value z0 on C2 is constant.
They are distinct over a general target infinity point. Split the two
special root supports over a henselian extension. Each root-order factor
is free, and normalization preserves its generic rank. Its special root
multiplicity therefore equals the sum of e*f of the normalization primes
in that cluster. Consequently those two multiplicities are exactly 1,5.

Thus the homogeneous quadratic part of G has one moving root and a
fixed root of multiplicity five. A normalization-index defect cannot
change these cluster ranks. This is the same flat-rank transport used
in [PF2OMIT3](ADMISSIBLE_INVERSE.md), now at target infinity. We keep
ramification indices, root multiplicities and order indices distinct.

## 2. A regular primitive pencil has a finite projective inverse order

We will also use the pencil on the affine N. More generally let N be
normal and finite over A2, and z:N->P1 a morphism generating its function
field over C(x,y). The graph map N->A2 times P1 is finite: it is the
closed graph in N times P1 followed by a finite base change. Its image
is the primitive binary equation G=0. This image is finite over A2,
so the coefficients of G have unit ideal; otherwise it would contain
a whole P1 over a target point. N is precisely its normalization.

In particular a smooth point of this inverse hypersurface is already
a point of N with the same local morphism. This conclusion does not
require z to be a linear source coordinate or finite valued on the
whole source. It does **not** say that the omitted boundary is z=infinity.
That stronger source-chart assertion requires PF2OMIT3's source-finiteness
criterion. Here the regular pencil comes from Y and restricts to N.

## 3. The quartic case has an impossible fifth-power quadratic part

Use [PF2R7QBASE1, Section2](../jc2-r7-quartic-contact/PROOF.md).
All exceptional curves over the common infinity point are vertical in
this case. Their connected fibre has a single base value, chosen z0=0.
It contains the infinity points of E and C2 and the point of C1 over
target a=[0:1:0]. On C1 the pencil map has degree two. Accordingly

    G=h(z)*(A2(x,y)-z*B2(x,y))+l(z)*x+m(z)*y+n(z),
    deg h<=5, deg l,deg m,deg n<=6,
    A2(0,1)=0, B2(0,1)!=0,

where A2,B2 are coprime binary quadratics. PF2R7QBASE1 proves that
z(t) has its unique zero at t=infinity, of order four. No squarefree
derivative or nonimmersion assumption is needed in this step.

Section1's complete infinity cluster now gives the further equality

                         h(z)=c*z^5, c!=0.             (I1)

Indeed the moving factor A2-z*B2 is nonzero at z=0 over a generic
target direction; all five units of the other cluster are at z=0.

At t=infinity write p~t^6, q~t^10, z~lambda*t^(-4), lambda!=0.
Since A2 has no y^2 term, every term of the quadratic part in (I1)
tends to zero, with order at least four in 1/t:

| Monomial | Order |
|---|---:|
| z^5*x^2 | 8 |
| z^5*x*y | 4 |
| z^6*x^2 | 12 |
| z^6*x*y | 8 |
| z^6*y^2 | 4 |

The same lower bound holds after applying p*(partial_x+s*partial_y),
where s=q'/p'~(5/3)*t^4. Incidence and tangency both vanish on E.
If l or m has a term producing a negative order, the only possible
cancellation is between l*p and m*q. Their leading coefficient ratios
in these two identities are 1 and 5/3. Those distinct ratios forbid
such cancellation. Thus

                    ord_0 l>=2, ord_0 m>=3.             (I2)

Absent polynomials count as having infinite order. The surviving linear
terms have order at least two; the quadratic terms have order at least
four. Incidence then forces n(0)=0. Equations (I1),(I2) show that **every
coefficient of G is divisible by z**, contradicting primitivity.

This excludes the quartic case for the whole stated family. No universal
quartic-contact elimination or new bounded coefficient search is used.

## 4. The reciprocal profiles cannot supply a finite ramification point

The retained reciprocal equation and complete profiles are
[PF2R7G1](../jc2-r7-cubic-cover/PROOF.md):

    G=a(z)*(x*y-z*x^2)+b(z)*x+c(z)*y+d(z),
    deg a<=5, deg b,deg c,deg d<=6.

The map z:E^nu=P1->P1 has degree three. At t=infinity it has local
degree two over its limiting value z0, and there is one additional
finite parameter u above z0, of local degree one. Its complete possible
profiles are:

| z0 | Coefficient profile |
|---|---|
| finite, translated to0 | ord_0(a,b,c)=(5,0,2) |
| infinity | (deg a,deg c,deg b)=(0,3,5) or (1,4,6), with deg d<=deg b |

The first infinite profile gives deg_z G<=5, contradicting the actual
degree-six field. The second has quadratic part, in homogeneous base
coordinates [Z:W],

    W^4*(a1*Z+a0*W)*(W*x*y-Z*x^2), a1!=0.

Over a generic target infinity point its three distinct root supports
have multiplicities 4,1,1. That contradicts the complete 5,1 cluster
ledger in Section1. Hence only the finite profile could remain.

For that profile, a(z)=a5*z^5, b(0)=b0!=0 and c is divisible by z^2.
At z=0 the entire conic fibre is the target line

                          b0*x+d0=0.                   (I3)

In the affine inverse hypersurface,

    G_x|_(z=0)=b0,   G_y|_(z=0)=0,
    G_z|_(z=0)=b1*x+d1.

The hypersurface is smooth along that line, hence Section2 identifies
it locally with N. In coordinates (z,y) its finite projection has
Jacobian -G_z/G_x. On (I3), this equals the constant

                   -(b0*d1-b1*d0)/b0^2.                (I4)

Therefore this whole line is etale if that constant is nonzero, and
the whole line is ramified if it is zero. Yet the additional finite
point u on E lies on this fibre, so the constant must be zero. The
ramification divisor then contains the line (I3). Since its entire
support is the single prime E, E would be that line and its target
image would be a line, contradicting degree10. This excludes the last
reciprocal profile, including the special immersed multiple-root case.

The order-two vanishing of c is essential. For example
G=x+z*y+z^2 has G_x=1 and G_z|_(z=0)=y, so its ramification meets the
z=0 line at just one point. This valid finite quadratic control has
ord c=1, outside our profile. The argument does not forbid arbitrary
conic fibres from meeting a ramification curve.

## 5. Consequence for the requested omitted-divisor ledger

The balanced infinity case was excluded by PF2R7INF1, and the varying
cubic case by PF2R7V1. Sections3--4 exclude both remaining cases.
Thus **no map satisfies the entire single e2,f1 ramified-prime,
degree-six gap21 hypotheses**, without the source-finality import.

The conditional omitted-prime counts and exceptional lengths in
PF2OMIT1--4 remain correct implications and explicit controls. Under
the just excluded hypotheses they have no global Keller realization.
Their previously null target equations are not replaced by invented
equations or interpreted as empty lists for other ramification patterns.

There is a residual ledger if a degree-six application has an immersed
e2,f1 gap21 prime R but the ramification is not assumed entire. The
local smooth-branch argument of [PF2OMIT2](GLOBAL_LEDGER.md) makes every
ambient length on R equal two: other ramification primes cannot meet R.
Thus E_R=0. The all-prime identity is

    5=2+sum_(other ramified D)(e_D+E_D)+t,

where t counts omitted unramified primes. The single-prime possibility
has just been excluded. Since each other ramified e_D>=2, there is
exactly one other ramified prime, with these exhaustive possibilities:

| Other ramified index | Its excess E_D | Unramified omitted primes t | Its ambient lengths |
|---:|---:|---:|---|
| 2 | 0 | 1 | two at every finite point |
| 2 | 1 | 0 | two except at one point of length three |
| 3 | 0 | 0 | three at every finite point |

Its image differs from C. Indeed C has positive affine delta12, so an
immersive normalization has a value with at least two branches. If the
other ramified prime shared C, its finite normalization map onto C^nu
would have a point over each of those two distinct parameters. Boundary
unibranchness keeps those two points distinct. They contribute at least
2*2=4 omitted units in addition to R's four, exceeding degree six.

Thus such an application needs exactly two ramified primes with distinct
images, and at most one unramified omitted prime. The other ramified
residue degree and image equation are not determined here. Its exceptional
point and target value in the middle row remain unknown. These nulls
are recorded explicitly in the residual ledger.

If the same application has a singleton source fibre at zero outside
the finite normalization's branch locus, its five omitted points must
all be unramified. Only the first residual row has any such prime.
Consequently both ramified indices are two, both excesses vanish, and
there is exactly one unramified omitted prime D. Its image H has exactly
five normalization branches over zero, each contributing ambient length
one. Since zero lies on H but on neither ramified image, the full
nonproper locus has exactly three distinct components C,C',H. The five
immersed branches of H give affine delta at least10 at zero. This does
not produce their actual equations or prove that this global pattern exists.

The immersed birational parametrization (t^5-t,t*(t^5-t)) is a local
image control: its five distinct zero preimages have five distinct tangent
slopes, so the origin is an ordinary five-branch point of delta10. It is
not the asserted image of a Keller boundary or a global cover realization.

This completes that conditional r7 obstruction, not every normal-r7
terminal application. To exclude an intended terminal passport one
must still establish the entire single ramified-prime hypothesis or
supply its additional ramified rows. The unrestricted singleton-zero
ledger still requires the actual normal algebra, all omitted primes,
their images/e/f and individual exceptional ambient Artin lengths.

## Replay and proof boundary

[verify_infinity_completion.py](verify_infinity_completion.py) checks
the homogeneous root clusters, the six exact leading coefficient
constraints, reciprocal field-degree profiles, line-fibre Jacobians
and positive controls outside the hypotheses. Its
[receipt](infinity_receipt.json) binds the prior conic-profile data and
canonical proof inputs. [infinity_ledger.json](infinity_ledger.json)
records the excluded cases with separate indices, residue degrees,
root multiplicities, ramification and omission fields.

Finite graph normalization, henselian cluster rank transport and the
inherited compactification/pencil results are written proof steps.
The algebraic replay does not independently certify those geometric
inputs, formalize the argument or prove unrestricted properness.
