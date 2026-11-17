# Monic keys force a Jacobian factor; whole-chart necessity in every index

16 November 2026. PF2RULED99 gives two degree-free algebraic results and
completes the specified PF45/40/46 necessity audit with explicit primary
imports. A monic-key valuation can force EVERY pair of regular target
polynomials to have a common Jacobian factor. Separately, whenever a
regular source coordinate exists on a smooth affine one-line extension,
the full root planes, support range, exactness test and minimal source
additive action follow for EVERY positive generic index e.

The positive-index whole-chart theorem is conditional on coordinate
existence and smoothness. It does not assert these for arbitrary
higher-index omitted primes. The two-index coordinate-existence bridge
is audited below; higher-index all-polar packets remain open.

During this audit the shared authority imported the PUBLISHED Domrina
degree-four exclusion as HC4MYKR0. Its unrestricted main theorem and
definition of generic sheet degree were independently checked in the
[primary full text](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&option_lang=eng&paperid=273&what=fullteng),
printed1--2, together with the conclusion on printed32. Therefore
quartic plane Keller REALIZATION is already excluded by literature.
The present degree-free results remain useful; the unfinished internal
quartic chain is optional independent-proof work, not a realization
frontier. This does not certify that chain or solve the all-degree goal.

## 1. An all-degree obstruction from one monic key

Let v be a discrete valuation of C(x,y), trivial on C. Work in a DVR
presenting it, or an unramified extension with the same integer values.
Suppose, for a uniformizer u,

    v(x)=-m<0, v(y)=-n<0,
    in(u^m*x)=A in C*, in(u^n*y)=B in C*.              (1)

Put g=gcd(m,n), M=m/g. Suppose G in C[x,y] is MONIC of y-degree M,
and

    v(G)=s>0, in(u^(-s)*G)=gamma notin C.             (2)

Here in means residue of a valuation-zero element. Since C is
algebraically closed, gamma is transcendental over C. Conditions(1)--(2)
are genuine leading-residue hypotheses, not consequences merely of
the three integer values. No canonical coefficient, field degree,
ramification index, finite target map or smooth surface is needed for
the following algebraic lemma.

Monic division gives a UNIQUE finite expansion of every polynomial:

    F=sum c_(i,j,k)*x^i*y^j*G^k,
    i,k>=0, 0<=j<M, c_(i,j,k) in C.                 (3)

For example this follows by presenting C[x,y] as the free C[x,T]-module
with basis1,y,...,y^(M-1), using the monic relation T=G(x,y).
The value of a displayed monomial and its leading coefficient are

    -m*i-n*j+s*k,        A^i*B^j*gamma^k.            (4)

At a fixed value and fixed k, there is at most one pair(i,j): equality
would give (m/g)*(i-i')+(n/g)*(j-j')=0, hence M divides j-j'; the
range in(3) makes j=j' and i=i'. At a fixed value the remaining terms
have distinct k, so their nonzero leading coefficients form a nonzero
polynomial in gamma. They CANNOT cancel. Therefore v(F) is exactly
the minimum in(4) over its nonzero coefficients.

Write A_v=C[x,y] intersect O_v. Every monomial in the expansion of a
member of A_v has nonnegative value. For k=0 this forces i=j=0.
Consequently

    A_v subset C+G*C[x,y],  image(A_v mod G)=C.       (5)

If P,Q in A_v, write P=p0+G*H and Q=q0+G*K. Direct differentiation
gives the stronger certificate

    Jac(P,Q)=G*(H*Jac(G,K)+K*Jac(H,G)+G*Jac(H,K)).    (6)

Thus EVERY such Jacobian is divisible by the nonconstant polynomial G.
In particular it cannot be a nonzero constant. Equivalently, every
regular polynomial algebra at this valuation contracts the source
curve G=0 to a point. This does not depend on finite generation of A_v.

For ANY actual plane Keller normalization, both target coordinate
pullbacks P,Q are regular at EACH omitted-prime valuation. Hence no
omitted prime can satisfy(1)--(2). Other omitted primes, all their
indices, the total field degree and the target image curve are
unrestricted in this conclusion. It is a single-prime necessary
obstruction to the WHOLE finite algebra; it does not turn a deleted
component open into a finite target cover.

The requirement gamma notin C is essential. An actual smooth affine
control and its precise cancellation appear in Section5.

## 2. Whole planes for every regular-coordinate positive-index sector

Let N be a smooth affine complex surface containing marked U=A2_(x,f),
with sole complement D=A1 and

    div_N(dx wedge df)=(e-1)D, e>=1.                 (7)

Assume f is an actual source coordinate regular on N. It is constant
on D, by the retained BGV coordinate-extension lemma; subtract that
value. Put r=v_D(f)>=1. The complementary coordinate x must have a
pole. If both coordinates were regular, the normal intersection
description O(N)=C[x,f] intersect O_(v_D) would be C[x,f], and the
specified source open would already be the whole affine N.

The source coordinate line L={f=0} is closed in N and disjoint from D.
Indeed its normalized affine closure cannot add the unique point at
infinity of A1 without becoming P1, which is not affine. Thus

    div_N(f)=L+rD.

Normalize the ACTUAL surface in f=z^r. This finite cyclic cover N'
is affine. Over U it is exactly A2_(x,z), with retained L'={z=0};
near L' it is a smooth reflection, and near D it is etale because
locally f=s^r*unit. The finite etale cover over D=A1 splits into r
copies D_i=A1. The root surface is smooth everywhere.

Remove L' and all D_j except D_i. This gives an AFFINE open X_i:
it is the complement of an effective Cartier divisor in smooth affine
N'. The map z:X_i->A1 is affine, faithfully flat and of finite type.
Its generic fiber is the source affine line. Every nonzero fiber is
that full line, and its zero SCHEME fiber is precisely the reduced
line D_i, since ord_(D_i)(z)=1. Thus every fiber is geometrically
integral, including its scheme multiplicity. Flatness also follows
locally from torsion-freeness over the smooth base curve.

Miyanishi TIFR60, ChapterI Theorem4.1.1, applies with exactly these
hypotheses: locally noetherian, locally factorial integral base,
affine faithfully flat morphism, generic A1 and other geometrically
integral fibers. The complete printed24/PDF40 was freshly read, as
was the local algebra on PDF41. The theorem is imported, not
independently reproved. Its application is NOT made merely from the
reduced set underlying a multiple fiber.

The resulting A1-bundle on A1 is trivial: its line bundle has Picard
class zero and its affine torsor class has H1(A1,O)=0. Hence

    X_i=A2_(z,t_i),

and these WHOLE planes together with A2_(x,z) cover N'. Their common
overlaps are Gm_z times A1. An automorphism of the one-variable
polynomial algebra over C[z,z^-1] is affine in that variable, so

    x=a_i*z^k*t_i+h_i(z), a_i in C*, h_i Laurent.

The pulled-back volume has order e-1 along D_i, since the root cover
is etale there. Its order from this formula is r-1+k, giving

    k=e-r,  x=z^(e-r)*t_i+h_i(z).                   (8)

Rescale t_i to remove a_i. Polynomial shifts of t_i remove exactly
the powers of h_i with exponent at least e-r. Thus the unique gauge is

    h(z)=sum_(j>=r-e+1) c_j*z^(-j), finite support,
    q=max{j:c_j!=0}>=1.                            (9)

Negative or zero j ARE ALLOWED when r<e. A pole remains because x is
polar. Deck equivariance gives h_i(z)=h(zeta^(-i)*z), with an indexing
choice, and separatedness requires all h_i distinct. Otherwise two
zero lines would give different limits of the same punctured path.
Equivalently gcd(r,{j:c_j!=0})=1, with all signed support indices kept.

This proves whole-chart NECESSITY at arbitrary positive e and r.
No converse affineness theorem for arbitrary proposed data is asserted.
For e=1, (9) begins at j=r; for e=2,r>=2 it begins at j=r-1; for
e=2,r=1 it also retains the constant term. These are precisely the
necessity statements used by PF40/43/46. Their affine converses are
separate theorems and are not needed to infer these charts from N.

## 3. The residue and minimal additive direction in every index

On the root planes take polynomial primitives

    theta_src=r*z^(r-1)*x dz,
    theta_i=r*z^(e-1)*t_i dz.

Their difference on the overlap is r*z^(r-1)*h_i dz. Its only residue
term is j=r and equals r*c_r, the same on every chart. Therefore

    omega globally algebraically exact iff c_r=0.    (10)

Necessity uses the polynomial Poincare lemma on the whole planes.
For sufficiency, residue-zero Laurent differences have Laurent
antiderivatives. Choose these and their differences as a Cech cocycle
on the plane cover. Since the existing N' is affine, H1(N',O)=0;
subtracting differentials of regular chart functions glues a global
primitive. Average over the finite deck group. Etale descent applies
near each D_i; across retained L' an invariant dz coefficient is
z^(r-1) times a function in x,z^r, so it descends as a regular df
coefficient. This proves(10) on N, including its branch line.

No bound on e is used here: e>=1 is exactly what makes theta_i
polynomial. In an actual finite Keller normalization the primitive
P*dQ/c already exists, so(10) is necessary without any extra target
invariance assumption.

For a nonnegative integer a the source field f^a*partial_x has boundary
expression

    z^(r*(a+1)-e)*partial_(t_i).                     (11)

It extends regularly if and only if r*(a+1)>=e. The LEAST multiplier is

    a_min=max(0,ceil(e/r)-1).                        (12)

The field descends, is locally nilpotent on O(N) as a restriction of
the source locally nilpotent derivation, and has boundary vanishing
order r*(a+1)-e. In particular partial_x alone need not extend when
r<e. This provides an exact source additive direction in every
smooth regular-coordinate positive-index sector. It says nothing
about preservation of the actual target algebra C[P,Q].

## 4. Audit of the all-polar index-two bridge

The following source checks resolve the necessity question posed at
the PF98 checkpoint. Mondal's published sign criterion, Michigan
Math.J.66(2017), Theorem1.4(1), complete PDF3/printed229, has exactly
the sign quantifiers used by PF45. The author HTML1307.5577v4,
Definitions2.8/4.5, Theorem4.6(1), Definition7.1 and Theorem7.2 apply
to ANY divisorial semidegree, before its separate primitive-algebraic
specialization. The normal-form existence argument in AppendixD.1
and the full canonical proof in Section7.3 were also read. Negative
final key values are permitted. These are explicit primary imports,
not new proofs of the valuation classification.

Put delta=-v_D. A nonzero regular function vanishing on D has negative
delta; the sign criterion makes the last key value negative. If all
source coordinates are polar, the normal-form weighted case is
impossible, and the initial pair has

    w0=g*A, w1=g*B, A>B>=2, gcd(A,B)=1.

Let alpha_i be the successive gcd ratios. The PRIMARY general formula
has the sign

    e=1+v_D(dx wedge dy)
      =-w0-w_last+sum_(i=1)^n (alpha_i-1)*w_i.        (13)

Set C0=-w0, C_i=C_(i-1)+(alpha_i-1)w_i,
z_i=w_i-C_(i-1). The strict key inequalities give
z_(i+1)<z_i and z_last=-e. Stop at the first nonpositive key value.
The preceding positive contributions imply

    g*(A*B-A-B)=C1<=e-1.                            (14)

This reproduces the written all-index bound without assuming that
an affine partial extension is a primitive compactification. For e=1
there is no all-polar case. For e=2 the only initial pair is(3,2),
its gcd is already1, and the final key value is-1. This does not
exclude all-polar higher-index cases.

BGV extends the pole-three coordinate to P1 with infinity fiber3D.
Its normalized x-root cover is exactly the WHOLE torsor construction
of PF37: remove the other boundary copies, and apply the reduced
scheme-fiber bundle theorem over P1. With x=u^-3 the volume order
forces its infinity transition y=u^5*t+h(u). Gauging h and making
the specified affine source shear/translations gives the entire form

    x=u^-3,
    y=u^-2+a*u^2+b*u^3+c*u^4+u^5*t.                 (15)

The shear and translations change the root cover, but the whole
torsor argument applies anew to the resulting source coordinates;
a formal substitution alone is not the existence argument. The
two planes' primitive difference is

    3*(u^-6+a*u^-2+b*u^-1+c)du.

Exactness forces b=0. The MONIC polynomial

    G=y^3-x^2-3a*y-3c

then has value1 and normalized leading residue3t. It satisfies
Section1 with(m,n,M,s)=(3,2,3,1). Thus every target polynomial
regular at that prime is constant mod G, and its pairwise Jacobian
has factor G. Alternatively the normal intersection algebra contracts
G=0 and cannot contain the marked affine plane as an open. This
verifies the exact coordinate-existence bridge used by PF45 and its
componentwise application via the already verified affine opens.

Combined with Section2, this audits the source-coordinate and WHOLE
chart NECESSITY input for every smooth index1/2 component open.
No source action has been pushed to the target, and removing other
primes has not been used to preserve finiteness.

## 5. A genuine affine countercontrol when the leading residue is constant

Use the retained exact smooth affine surface

    N=Spec C[f,W,V]/(f^2*V-W*(W-1)),
    D={f=0,W=1}, U=N minus D=A2_(x,f),
    W=f^2*x, V=x*(f^2*x-1).

At D a whole chart is f=u, x=u^-2+t. Put the NEW source coordinates
X=x,Y=x+f and monic G=Y-X=f. Both X,Y have pole2, M=1, v(G)=1,
but the normalized leading residue of G is the CONSTANT1. The
regular function

    V=G^2*X^2-X

has value0 and nonconstant residue t; the two pole-order2 terms
cancel. On G=0 it equals-X, not a constant. Thus(5) fails when
gamma is constant, on an actual smooth affine plane extension.
The surface is not claimed to carry a finite Keller map. This keeps
the new obstruction from incorrectly ruling out known source models.

Similarly the local higher-index control e5,r3,h=z^-6+z has an
essential retained positive power in(9), yet c_r=0. It is not
declared a globally affine realization here. The all-index chart
and exactness laws do not justify discarding such signed support
terms or automatically applying the index1/2 contact bounds.

## 6. Published quartic closure and remaining all-degree obligations

This closes the specified PF98 key-form/whole-chart NECESSITY audit,
using the explicitly retained primary normal-form, canonical, sign
and bundle statements. It adds the all-degree monic-key Jacobian
factor and the all-index regular-coordinate chart/residue/action
theorem. It does not independently reprove the external valuation
classification or all prior source affine-converse theorems.

The shared [published-theorem adoption](../hc4-triangular-graphs/QUADRATIC_SOURCE_RECONSTRUCTION_AND_QUARTIC_RETIREMENT.md),
HC4MYKR0, excludes ALL complex plane Keller maps of geometric degree4.
Domrina's PartII theorem on printed1 has no boundary-count, monodromy,
image irreducibility or polynomial-degree restriction. Proposition1.2's
two-dicritical reduction is a proof step, not a theorem hypothesis.
This note independently checks that precise statement, not the full
33-page splice proof. It does not claim a new quartic theorem.

The application to the ORIGINAL finite-completion requirement is
immediate but must retain its marking. If pi:N->A2 is finite dominant
of degree4 and a DENSE MARKED OPEN U in N is A2, with pi|U etale,
then target coordinates restrict to polynomials P,Q on U. Their
Jacobian is nowhere zero, hence a nonzero constant over C. Density
gives C(N)=C(U), so [C(x,y):C(P,Q)]=4. Domrina excludes this actual
Keller-compatible finite completion, regardless of its other boundary
data. A completion satisfying only numerical or local necessary
conditions has not thereby been shown to give such a marked map.
Arbitrary finite quartic maps with interior ramification are not
excluded by this argument.

Accordingly the current priority is HIGHER-DEGREE boundary reduction:
use Section1 to exclude monic primitive-key prime patterns directly
from the finite target algebra, and Sections2--3 for smooth regular-
coordinate primes at any positive index. Remaining cases include
nonconstant-leading-residue hypotheses that fail, Laurent terminal
keys outside the monic-degree condition, higher-index all-polar
packets, and singular boundary charts. Apply the canonical budget
and actual finite-algebra constraints without deeper coefficient
sweeps. The local passport, normal deformation and later JSJ chain
can still be audited for all-degree extensions or an independently
verified proof; those audits are not required to settle degree4.
No full all-degree/Keller theorem or research-goal completion is
registered here. Preserve all older exact controls and failed routes.

The [checker](verify_monic_key_chart_audit.py),
[receipt](monic_key_chart_audit_receipt.json), and
[source manifest](sources/monic-key-chart-audit/manifest.json) retain
exact polynomial identities, primitive residues, monic expansion
controls and the essential constant-residue counterexample. These
controls support the written unbounded proofs; they do not infer
universal statements from a finite search or construct a Keller map.
