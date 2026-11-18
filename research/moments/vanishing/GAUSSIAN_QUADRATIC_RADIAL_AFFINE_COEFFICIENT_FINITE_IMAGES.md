# Quadratic radial Gaussian phases with affine leading coefficient

Let K be a characteristic-zero field. Use an original Witt Gaussian
pair Z,W and an independent ordinary Gaussian Y, normalized by
E[Z^iW^j]=delta_(i,j)i! and E[Y^(2j)]=(2j-1)!!. Put S=ZW and consider

    P=P0+R, P0=A(Y)S^2+H(Y)S+G0(Y), deg A<=1.       (Q1)

The polynomials H,G0 are arbitrary. Every monomial of R has strictly
positive Z-minus-W weight, or strictly negative weight after reversing
the pair. R=0 is allowed; there is no degree bound on R.

For ORIGINAL polynomial markers define

    Phi_P[Q]=E[Q exp(tP)], V_P=Phi_P[K[Z,W,Y]],
    M_P=K[t]V_P.

**Finite-image theorem.** The full actual M_P is finite over K[t] if
and only if

    A=H=0 and G0=c+bY.                              (Q2)

In this branch

    Phi_P[1]=exp(psi), psi=ct+b^2t^2/2,
    M_P=exp(psi)K[t].                                (Q3)

Thus every finite vacuum in (Q1) is primitive. If P is nonconstant,
V_P=M_P; if P=c is constant, V_P=K exp(ct). Safety, retaining every
original marker, holds precisely when P0=0. In particular safety
implies purity throughout this family.

In the ZERO-boundary subfamily G0=0 there is also an exact converse:
polynomiality of the unit, purity and safety are equivalent, and all
require A=H=0. This restricted assertion is proved in Section 9.

This is a FULL finite-image theorem. It does not assert that purity
alone forces safety for arbitrary G0, or classify scalar identities
p exp(psi) with an arbitrary exponential psi. The finite-image premise
is used to fix psi
from the entire original boundary image. The proof retains an arbitrary
polynomial factor p in the hypothetical vacuum until necessity is proved.

The main new exclusion is A!=0 and deg H>=2. It uses an exact rank-two
relative Gaussian connection, two sectorial Gaussian endpoint tails,
and a complete relative period frame. The written proof and primary
comparison imports retain their internal-review scope; independent
replay, formal verification and external review are not asserted.

## 1. Original markers, angular reduction and the boundary

For an invariant core P0(S,Y), a monomial Z^iW^jY^ell has zero period
when i!=j. When i=j it is exactly S^iY^ell. Consequently the complete
original core image is represented by K[S,Y] with E_S[S^n]=n!, not by
a smaller class of distinguished amplitudes.

The [angular-core theorem](GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md)
gives equality of actual embedded modules and preservation of vacua:

    M_P=M_P0, Phi_P[1]=Phi_P0[1].                    (Q4)

We may prove necessity on P0. The
[general finite actual-image theorem](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image)
then supplies

    M_P0=exp(psi)K[t], Phi_P0[1]=p exp(psi),
    p in K[t], p(0)=1, psi in tK[t].                 (Q5)

Neither p=1 nor universal safety implies purity is assumed.

The [boundary-image theorem](GAUSSIAN_INVARIANT_CORE_BOUNDARY_REDUCTION.md)
applies to the FULL image and gives

    G0(Y)=c+bY, psi=ct+b^2t^2/2.                    (Q6)

Its exact original endpoint amplitude is 1-tP0_S. In particular it
does not discard the Gamma endpoint to obtain (Q6).

If A=0, the [S-affine classification](GAUSSIAN_AFFINE_RADIAL_INVARIANT_CORE_FINITE_IMAGES.md)
already forces H=0. If deg H<=1, after (Q6) the phase is linear in Y
and belongs to the [quadratic-Y classification](GAUSSIAN_VARIABLE_RADIAL_QUADRATIC_FIBRE_FINITE_IMAGES.md),
which forces A=H=0. It remains to exclude

    A=alpha Y+beta !=0, deg H=k>=2.                 (Q7)

Work first over C, put r=1/t and gauge by exp(-psi). The spatial phase is

    phi=(AS^2+HS)/r-S-Y^2/2+bY/r-b^2/(2r^2).        (Q8)

## 2. The exact original relative source

Let B=C(t)[S,Y] and define

    D_S=partial_S-1+t(2AS+H),
    D_Y=partial_Y-Y+t(A'S^2+H'S+b).

The original period functional annihilates D_Y B and D_S(SB), by
ordinary Gaussian and Gamma integration by parts coefficientwise.
The second condition requires the primitive to vanish at S=0. The
full two-variable relative de Rham source is therefore

    H_rel^2=B/(D_Y B+D_S(SB)).                       (Q9)

Its complex is SB -> B dS+SB dY -> B dS wedge dY. The time connection
is partial_t+P0, or partial_t+P0-c-b^2t after the boundary gauge. All
identities are exact for the original formal Wick functional; an
analytic interpretation of a divergent original integral is not assumed.

On A!=0 set q=1-tH-2tAS. In the S-relative quotient,

    [qS^n]=n[S^(n-1)] (n>=1),                       (Q10)

because D_S(S^n)=nS^(n-1)-qS^n. The invertible leading coefficient 2tA
reduces every S monomial to 1,S. These two classes are independent:
a nonzero D_S(SF) has highest S degree at least two. Use the basis

    v=[dS], e_B=[q dS].

The class e_B is the actual endpoint class. Equations [qS]=1 and
[qS^2]=2S give

    nabla_Y e_B=(-Y+tb)e_B,
    nabla_t e_B=(bY+c)e_B.

On a finite square-root time and spatial cover put

    u=sqrt(tA)v, x=(H-r)/(2sqrt(rA)),
    G=-(Y-b/r)^2/2, F=G-x^2.                       (Q11)

Choose the cover consistently, so sqrt(r)=1/sqrt(t). After the boundary
gauge the exact connection in BOTH Y and t is

    nabla u=dF u-dx e_B, nabla e_B=dG e_B.          (Q12)

For verification, in the unscaled basis the Y coefficients on v are

    a_Y=-Y+tb-A'/(2A)
          +(tH'+A'(1-tH)/(2A))*(1-tH)/(2tA),
    k_Y=-(tH'+A'(1-tH)/(2A))/(2tA),

and the gauged time coefficients are

    a_t=bY-b^2t-H^2/(4A)+r^2/(4A)-r/2,
    k_t=-(r^2+rH)/(4A).

Scaling gives a_Y+A'/(2A)=F_Y, a_t+1/(2t)=F_t and
sqrt(tA)k_Y=-x_Y, sqrt(tA)k_t=-x_t. Thus the endpoint extension and
the original parameter connection are both present.

## 3. Spatial localization and ordinary-point comparison

For constant A the preceding reduction is global in Y. Otherwise let
D={A(Y)=0}, a smooth divisor transverse to the original boundary S=0.
The relative Gysin sequence for the polynomial exponential connection
contains

    H^0(D,D intersect {S=0};exp(phi|D))
       -> H_rel^2(A^2,{S=0};exp(phi))
       -> H_rel^2(A^2-D,{S=0}-D;exp(phi)).           (Q13)

The left term is zero. A horizontal algebraic function on the connected
S line that vanishes at S=0 is zero, by uniqueness for its first-order
linear equation. Hence the original source injects into the localized
one. The localization cokernel, which may contain a residue line, is
retained. No equality of the two source modules is asserted.

In particular A^-1 has not been treated as a t-adic unit or evaluated
by an inverse-marker version of the original Wick functional. Instead
choose an ordinary time point. A horizontal functional on the original
source extends linearly at that point to the finite-dimensional localized
source, then extends horizontally by analytic parallel transport. Its
restriction remains the original functional. This extension need not
be rational. It suffices for a period contradiction on an original seed.

The precise absolute period input is Hien--Roucairol, *Integral
representations for solutions of exponential Gauss--Manin systems*,
Bull. SMF 136(4) (2008), 505--532, Theorems 2.7 and 3.5
([published primary text](https://www.numdam.org/item/10.24033/bsmf.2564.pdf)).
Theorem 2.7 gives a perfect finite-dimensional period pairing for an
elementary exponential line on a smooth affine complex variety, in
arbitrary dimension. Theorem 3.5 represents every local holomorphic
horizontal functional away from a finite parameter set by its
rapid-decay periods. These are not quantitative decay statements and
are not directly theorems about an arbitrary rank-two curve connection.

For the pair (A^2,{S=0}), apply those two absolute pairings to the
relative de Rham cone and dualize its finite-dimensional chain cone.
The resulting relative pairing permits the finite boundary on S=0
and is horizontal. This provides ordinary-point relative periods of
the original two-variable exponential line. Reduction (Q12) is then
an explicit relative Gaussian fibre construction, not an application
of the rank-one theorem to an arbitrary rank-two system.

The Y long exact sequence gives

    0 -> H_Y^1(boundary G) -> H_Y^1(relative source)
      -> H_Y^1(absolute Gaussian F) -> 0.           (Q14)

The relevant H^0 terms vanish: neither a nonconstant quadratic-S
exponential nor exp(-G), exp(-F) can have a nonzero algebraic horizontal
coefficient. There is no H_Y^2 on the affine Y curve. For affine A the
boundary curve is punctured at A=0: its Gaussian frame has BOTH the
usual Gaussian contour and a puncture loop. The loop is never dropped.
For constant A the boundary frame has only the Gaussian contour.

## 4. All absolute critical points are bounded

Completing the S square in (Q8) gives

    phi=A/r*(S+(H-r)/(2A))^2+F,
    F=-(H-r)^2/(4Ar)-Y^2/2+bY/r-b^2/(2r^2).        (Q15)

The absolute Gaussian quotient has seed

    eta=sqrt(r)dY/sqrt(A).                          (Q16)

Its critical points are exactly the critical solutions of (Q8) after
S=(r-H)/(2A). They satisfy

    2AS+H=r, A'S^2+H'S=rY-b.                       (Q17)

A(Y)=0 is impossible on such a solution: it would fix Y in the
constant field and require the constant H(Y) to equal r. Eliminating S
gives the exact polynomial

    N=A'(H-r)^2-2AH'(H-r)-4rA^2Y+4bA^2=0.          (Q18)

If alpha!=0, its degree is 2k and its leading coefficient is
(1-2k)alpha h_k^2!=0. The competing rA^2Y term has degree at most
three, less than 2k. If alpha=0, its degree is 2k-1 with leading
coefficient -2k beta h_k^2!=0; the competing terms have degrees at
most one. Every coefficient is integral in r. Thus every Puiseux root
Y is integral: an escaping root would make the highest-degree term the
unique term of least valuation. Equivalently the coefficients after
monic normalization are uniformly bounded, giving one fixed root disk.

S is integral too. This is immediate for constant A. For affine A,
if v(S)<0, the term alpha S^2 in the second equation (Q17) has valuation
2v(S), strictly below that of H'(Y)S; its right side is integral, a
contradiction. This proof includes common roots, repeated roots and
multiple critical points; no reducedness is assumed.

If b!=0, every normalized absolute critical action has unique leading
term -b^2/(2r^2), since all other terms in (Q8) have valuation at least
-1. Hence none has zero exponential factor. If b=0, retain those
critical actions whose polar part vanishes; these actions are bounded
after finite ramification.

## 5. The full absolute seed estimate, including infinity

For constant A=beta, F is a polynomial of degree D=2k in Y with leading
coefficient -h_k^2/(4beta r). Its absolute source on A^1 has rank
D-1=2k-1, exactly the number of finite critical points counted with
multiplicity.

For affine A with alpha!=0, put V^2=A(Y). Then

    eta=(2sqrt(r)/alpha)dV,
    D=4k-2, C_D=-h_k^2/(4alpha^(2k)r).              (Q19)

F has exactly two generic phase poles: order D at V=infinity and
order two at V=0. The latter coefficient is
-(H(-beta/alpha)-r)^2/(4r), nonzero over C(r). The full exponential
source on the V punctured line has rank D+2=4k, exactly the finite
critical multiplicity. Its descended Gaussian half-density summand
has rank 2k. There is no extra regular-puncture mode. Keep the deck
summand corresponding to the original Gaussian fibre.

The finite proper phase map and the full scalar curve comparison give
a vanishing-cycle frame of these ranks. Multiple critical points are
handled by finite local morsification and their full local cycle groups,
not one selected branch. After a finite parameter cover, a sufficiently
small generic sector avoids critical/Stokes alignments; continuous
transport gives a flat frame. The retained proper-model comparison
includes finite phase poles, all finite-component necks, and monodromy.

For b=0 all critical points lie in a fixed coordinate disk. Write X=Y
or V as above; throughout that disk eta=c sqrt(r)dX, with c constant,
is regular, including at V=0. Descending inverse-phase arcs from a
zero-polar-action critical value have Re F bounded above. Their physical
length inside the disk is uniformly bounded: the inverse images of
phase rays are real algebraic curves of fixed degree, and the
curve-length bound controls their intersections with a fixed disk
independently of their coefficients. This holds for every branch in
a multiple local cycle group and for compact portions ending at the
finite pole. Thus their complete compact seed contributions are
O(sqrt(r)), up to finite logarithmic factors after ramification.

The volume dX has a double pole at infinity, so a logarithmic-volume
lemma is insufficient. Every other nonconstant coefficient of F has
valuation at least -1; for b=0 there is no scalar coefficient of
valuation -2. Increase the fixed disk so the top term dominates
uniformly on its outside. There |F|>=c/|r|. On a descending inverse-F
path from a bounded critical value, Im F stays bounded and Re F is
bounded above. Hence outside the disk Re F<=-c/|r|.

The exact inverse phase satisfies, uniformly on these tails,

    |dX/dF|<=C |r|^(1/D)|F|^(1/D-1).               (Q20)

This follows by dividing by C_D X^D and using top-term dominance for
both F and F_X. Equivalently the leading descending coordinate is
T=-C_D X^D, with seed

    (c/D)sqrt(r)(-C_D)^(-1/D)T^(1/D-1)dT.

Its prefactor has order 1/2+1/D>0. The exact tail (Q20), starting at
Re F<=-c/|r|, is exponentially small times an algebraic power. Finite
pole tails and intervening compact necks were already bounded inside
the disk; no complete Gamma integral is substituted for a truncated
finite-pole tail. This proves strictly positive order for the FULL
zero-exponential absolute seed block when b=0.

For b!=0 all absolute critical actions have nonzero polar part by
Section 4. The rank counts above account for the entire absolute
frame, including multiplicities. Thus there is no zero-exponential
absolute block. In particular no uncounted cycle at infinity can
supply one: the pole orders and leading coefficient are fixed, the
critical set is bounded, and the full finite critical multiplicity
equals the full source rank.

The precise scalar comparison inputs are the retained three facts in
[the affine-isotropic source](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md): regular
nearby cycles with proper-image compatibility; finite/pole-node and
finite/finite annulus complexes with all gluing maps; and ordinary
perfect exponential period pairing. Roucairol Section 1.2 and Lemma
1.4, and Sabbah Section 1.b, formula (1.4) and Proposition 1.5, supply
the first two comparison facts, with monodromy. They do not supply
the quantitative bounds proved above. The
[primary manifest](../../../unpublished/gaussian-boundary/sources/laurent_primary_manifest.json)
pins the exact inspected primary texts. The quantitative path input is
[analytic_core.tex](../../../unpublished/gaussian-boundary/analytic_core.tex),
Lemma lem:curvelength; ordinary seed-functional representation is
[seed_period_comparison.tex](../../../unpublished/gaussian-boundary/seed_period_comparison.tex),
Lemma lem:seedfunctional. Their existing review limits are unchanged.

## 6. Two Gaussian endpoint sectors with an exact jump

Let E(x)=integral_0^x exp(z^2)dz and define the entire functions

    J_+(x)=exp(-x^2)(i sqrt(pi)/2-E(x)),
    J_-(x)=exp(-x^2)(-i sqrt(pi)/2-E(x)).

They are the tails to +i infinity and -i infinity and satisfy

    J'=-2xJ-1,
    J_+-J_-=i sqrt(pi)exp(-x^2).                   (Q21)

On the fixed closed sectors

    Sigma_+: -pi/6<=arg x<=7pi/6,
    Sigma_-: -7pi/6<=arg x<=pi/6,

which cover the x circle with overlap, one has

    |J_+/-|<=C/(1+|x|),
    J_+/-=-1/(2x)+O(|x|^-3).                       (Q22)

Here is an elementary uniform proof. Put
I(z)=integral_0^infinity exp(-s^2-2zs)ds; then J_+(x)=iI(-ix) and
J_-(x)=-iI(ix). For |arg z|<=3pi/4-delta, rotate the s ray to
theta=-sign(arg z)(pi/4-delta/2), taking theta=0 at arg z=0. Entirety
permits the rotation, since quadratic decay dominates the linear term
on the connecting large arcs. On the rotated ray s=exp(i theta)v,

    Re(s^2)>=sin(delta)v^2,
    Re(zs)>=sin(delta/2)|z|v.

The integral is bounded by the minimum of
sqrt(pi)/(2sqrt(sin(delta))) and 1/(2sin(delta/2)|z|). One integration
by parts gives I(z)=1/(2z)-(1/z)integral s exp(-s^2-2zs)ds; the same
bound gives the O(|z|^-3) remainder. Taking delta=pi/12 proves (Q22).
These constants depend on the fixed sectors, not on r.

The exact local horizontal values of (Q12) are

    h_B=C_B exp(G), h_u=exp(F)(C_A-C_B E(x)).

On either sector this is

    h_u=C_B exp(G)J_+/-+(C_A-C_B C_+/-)exp(F),       (Q23)

where C_+/-=+/-i sqrt(pi)/2. The original seed v has endpoint density
eta J_+/- exp(G), and absolute density eta exp(F). Thus a change of
endpoint chart adds EXACTLY a constant multiple of an absolute
Gaussian period. No exponentially growing remainder is discarded.

At bounded Y, (Q22) gives the positive seed factor sqrt(r)dY/sqrt(A).
At an A-root this is regular on V^2=A. At Y=r^-N times a unit, N>0,
x is large and the endpoint amplitude is

    eta J_+/-=-r/(H-r)(1+O(x^-2))dY,               (Q24)

of volume order 1+(k-1)N>0. The sign in (Q24) and the common branch
sqrt(r)=1/sqrt(t) are retained.

## 7. A complete relative frame, including endpoint gluing

It remains to build actual lifts of the entire boundary frame in
(Q14). This keeps the extension and its Stokes jumps, instead of
semisimplifying the rank-two connection.

First let b=0. Choose a Gaussian contour between the two real decay
ends, with fixed small algebraic detours around every root of A,H.
For affine A also choose a fixed small loop around its root, avoiding
all roots on the loop itself. These paths represent the complete
boundary frame on the punctured original Y line. On them

    |H(Y)|>=c(1+|Y|)^k,
    |x|>=c|r|^-1/2,
    |eta|/(1+|x|)<=C|r||dY|/|H-r|.

Consequently their endpoint pieces have O(r) seed integrals, with
Gaussian decay on the unbounded ends. This construction works for
common and repeated roots: detour around the fixed roots, rather than
forcing a path through their moving collision regions.

Choose J_+ on upper-half-plane x pieces and J_- on lower pieces. The
switches may be placed at real x. There are only finitely many, uniformly
in small r: on each fixed straight or algebraic detour segment the
equation Im x=0, after squaring and clearing denominators, is a real
algebraic equation of fixed degree. Identically real segments need no
switch. A generic parameter ray and a small sector keep the nontrivial
switches continuous and avoid critical/Stokes alignments.

At a switch Y_j, Re F(Y_j)=Re G(Y_j)-x(Y_j)^2<=C-c/|r|.
The discontinuity (Q21) is a closed Gaussian S-thimble. Attach it over
a descending inverse-F path from Y_j to a phase pole, avoiding the
finite critical values. Small detours around those values stay in a
left half-plane Re F<=Re F(Y_j)+C. The attached surface closes the
jump and decays at its other end. Algebraic inverse-root and length
bounds contribute at worst a power of |r|^-1; the initial exp(-c/|r|)
dominates it. The double-pole infinity tail is controlled by (Q20).
Thus the corrections are exponentially small, and the resulting lift
is a genuine relative rapid-decay cycle. Its homology class varies
continuously, hence flatly, on the selected sector. No jump is set to zero.

For b!=0 choose a generic parameter ray with |Im(b/r)|>=c/|r|. The
boundary Gaussian contour Y=b/r+z, z real, stays |Y|>=c/|r| and has
G=-z^2/2. Equation (Q24) gives an O(r^(k+1)) endpoint integral. At a
switch its F action is at most -c|r|^(-(2k+1-deg A)); the same closing
surfaces are exponentially small. Here 2k+1-deg A>=4, so this starting
phase level dominates the scalar gauge term of size O(r^-2). On the
far infinity part of a closing path, exact top-term inversion therefore
again gives (Q20); its other compact parts require only polynomial
algebraic length bounds. The b=0 fixed-disk dominance argument is not
being assumed unchanged for b!=0. The extra puncture-loop class, when
A is affine, has polar boundary action

    G(Y_A)=-b^2/(2r^2)+bY_A/r-Y_A^2/2.

It has no zero-exponential boundary factor. By exactness of formal
exponential decomposition, the zero-exponential part of any lift of
this factor lies in the absolute part, already treated in Section 5.

A full absolute-F frame, together with these lifts of the FULL boundary
frame, has a block-triangular period matrix with invertible diagonal
blocks. It is therefore the complete relative frame (Q14), including
the extension, finite poles, infinity tails and all closing chains.
Every zero-exponential seed period in that frame has strictly positive
r-order, up to finite logarithmic factors.

On V^2=A, retain the appropriate deck summands throughout. The boundary
line descends invariantly and the absolute Gaussian half-density has
its matching deck sign. The identity J_-(-x)=-J_+(x), together with
the sign of 1/sqrt(A), makes the endpoint density descend. Additional
boundary cycles on the whole cover are not silently counted as
original-line cycles.

## 8. Excluding the normalized scalar unit

The exclusion just proved needs only the ACTUAL scalar identity

    Phi_P0[1]=p(t)exp(ct+b^2t^2/2), p in C[t],p(0)=1,  (Q26)

under (Q7) and G0=c+bY. Full image finiteness was needed in Section 1
to force this particular boundary normalization. It is not needed for
the cyclic seed argument that follows.

The gauged original Wick functional is defined on the full original
source (Q9), with values in C((t)); every rational-time amplitude has
a finite Laurent principal part and polynomial spatial coefficients.
Under (Q26) it takes w=v/p to one. Its restriction to the rational
cyclic differential submodule U=C(t)<partial_t>w therefore takes
values in C(t): it sends w,nabla_t w,nabla_t^2 w,... to 1,0,0,... .
Every rational relation is respected because the formal functional was
defined on the whole original source. Thus U has a horizontal rational
unit quotient, without any rationality assertion for other markers.

This is exactly the cyclic argument of Lemma lem:seedfunctional in
[seed_period_comparison.tex](../../../unpublished/gaussian-boundary/seed_period_comparison.tex).
Choose a time point where the original source, U and its quotient are
ordinary. Extend the quotient functional from that fibre of U to the
original source fibre, then through the injection (Q13) to the localized
source fibre. Parallel transport gives an analytic horizontal extension;
uniqueness preserves its restriction to U. Its value on w remains
identically one, and Section 3 supplies its complete relative period
representation. No convergence of the other original formal periods,
or rational horizontal extension to the whole source, is assumed.

Division by a nonzero polynomial p adds deg(p)>=0 to all r-orders.
By Sections 5-7 the zero-exponential part of EVERY such seed period
has strictly positive order. Distinct formal exponential factors are
independent after ramification, with monodromy and logarithms retained.
A period identically one has zero-exponential component one, giving a
contradiction. Thus (Q26) is impossible under (Q7) for every polynomial
p. Together with (Q5)-(Q6), this proves necessity in (Q2).

The argument is first over C. A hypothetical finite identity over K
uses only finitely many coefficients of A,H,G0,p. Their finitely
generated characteristic-zero field embeds in C, preserving the
nonzero leading coefficients and every formal identity. Every original
monomial marker has its normalized Taylor coefficients in that same
phase/psi field, so its polynomiality transfers coefficientwise; no
single uniform marker cutoff is needed for this transfer. Necessity
therefore holds over every characteristic-zero field. Sufficiency
below is algebraic over K.

## 9. Every original marker in the finite branch

For P0=c+bY let gamma>=1 be the least positive angular weight of R!=0,
after choosing its orientation. If an original marker Q has total
spatial degree D, angular balance and ordinary Gaussian translation give
the exact finite formula

    Phi_P[Q]=exp(ct+b^2t^2/2)
       *sum_(n=0)^floor(D/gamma) t^n/n!
          E[Q(Z,W,Y+bt)R(Z,W,Y+bt)^n].              (Q25)

For R=0 take only n=0. Every summand is polynomial; Q=1 has only n=0,
giving the primitive vacuum (Q3). This proves full finite-image
sufficiency without enlarging the original marker ring. The retained
fixed-marker theorem gives V=M for nonconstant P; the constant-phase
exception is stated at the beginning.

Safety forces full finiteness and then (Q3) to be polynomial. Over a
field this requires c=b=0. Conversely (Q25) with c=b=0 is polynomial
for every marker. Thus safety is exactly zero core and implies purity.
No converse from purity outside the finite branch is inferred.

There is a stronger conclusion when the boundary itself is zero.
For G0=0 in (Q1), preservation of the unit under R gives

    Phi_P[1] is polynomial
      iff Phi_P[1]=1
      iff P is safe
      iff A=H=0.                                    (Q27)

Indeed a polynomial unit would be (Q26) with c=b=0. Section 8 excludes
A!=0,deg H>=2. The S-affine scalar criterion excludes H!=0 when A=0;
the retained quadratic-Y scalar criterion excludes the remaining
A!=0,deg H<=1 cases. Both imports allow a polynomial factor p in the
unit, so no primitivity is being assumed. Conversely A=H=0 leaves
only the strictly one-sided R, whose unit is one and whose every
original marker is safe by (Q25). The restriction G0=0 is essential
to this argument: polynomiality of the unit alone has not been used
to force an arbitrary G0 to vanish or become affine.

There is an algebraic consequence for finite moment calculations. Fix
a support bound deg H<=k and take the independent coefficients of
A=alpha Y+beta and H=sum_(j=0)^k h_jY^j as parameters over Q. Let
m=(alpha,beta,h_0,...,h_k), let mu_n=E[P0^n], and put
I_N=(mu_1,...,mu_N), I_infinity=(mu_n:n>=1). Each mu_n is homogeneous
of coefficient degree n. Equation (Q27) says that the common field
zero set of I_infinity over an algebraic closure is exactly the origin.
The Nullstellensatz and Noetherianity therefore give

    radical(I_infinity)=m,
    m^nu subset I_N for some finite N,nu.           (Q28)

Thus this fixed zero-boundary moment ideal has finite colength. No
finite-moment extrapolation enters (Q28): finite generation chooses a
prefix generating I_infinity, and the finitely generated radical m
has a power in that ideal. No
effective N, nu or length is claimed. Once a concrete computation
certifies such a containment, homogeneity puts every moment mu_n with
n>=nu directly in I_N. To certify I_infinity=I_N only the intervening
moments N<n<=nu-1 need further checking. This is also the zero-base
C=0 case of the mixed-cutoff jet theorem; no generic cubic budget is
needed at that base point.

The same exact containment has a ring consequence. Over any commutative
Q-algebra satisfying the certified prefix, the ideal J generated by
the actual core coefficients has J^nu=0, hence P0^nu=0 as a polynomial.
With a strict R of minimum angular weight gamma and an original marker
Q of degree D, expansion of (P0+R)^m uses only j<=nu-1 surviving core
factors and at most floor(D/gamma) strict factors. Consequently

    E[Q(P0+R)^m]=0
       for m>nu-1+floor(D/gamma).                  (Q29)

For R=0 the cutoff is m>=nu. This is a coefficient-nilpotence and
marked-tail statement over rings. It does not force the coefficients
to be zero or the unit to be primitive over a nonreduced Q-algebra.

In the safe branch the valid cutoff is m>floor(D/gamma), hence C=1
is always valid. With b!=0 a valid normalized degree bound is
D+(1+deg_Y R)floor(D/gamma); for R=0 it is D. These are upper bounds,
not filtered-surjectivity assertions for arbitrary R. The
[mixed-cutoff theorem](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md) therefore
gives degree-d perturbation r-jets through moment
floor((1+d/gamma)r) in the safe branch with R!=0; (d+1)r is always a
valid looser budget. No optimality for arbitrary R is asserted. For
R=0 the safe phase is zero and the sharper mixed constant is C=0,
giving budget r. The theorem's local/global certificate requirements
remain in force.

The [safe-recomposition theorem](GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md)
also gives the precise safe-recomposition ideal: it is (U-c) when
A=H=0 and G0=c is constant, and the zero ideal otherwise. Any safe
shift would itself have finite full image, so no excluded radial core
can be rescued by nonlinear univariate recomposition.

## 10. Full invariant finite images through original degree five

Let P0(S,Y) be ANY invariant core of original spatial degree at most
five, and allow arbitrary one-sided R. Full finiteness first makes
its boundary affine. Its remaining support is exactly

    P0=S H3(Y)+S^2(alpha Y+beta)+bY+c,
    deg H3<=3.

It belongs to (Q1), so full finiteness holds exactly when P0=c+bY.
Its finite vacuum is primitive and safety holds exactly for P0=0.
This closes the full finite-image quintic frontier, including its
mixed-parity cases. Pure invariant quintics without a finite full-image
premise, and arbitrary scalar exponential normalizations, are not
classified by this corollary.

After the same boundary reduction, an invariant core of original
spatial degree at most six has support

    delta S^3+S^2 A2(Y)+S H4(Y)+bY+c,
    deg A2<=2, deg H4<=4.

The [logarithmic-cover theorem](GAUSSIAN_QUADRATIC_RADIAL_LOGARITHMIC_COVER_FINITE_IMAGES.md)
also classifies every quadratic A2 and arbitrary H. Thus outside the
retained quadratic-Y family the remaining sextic support requires
delta!=0 and deg H4>=3. The logarithmic-cover theorem additionally
excludes every scalar p exp(psi), with arbitrary polynomial psi, when
deg A>=2 and all roots of A have multiplicity at most two. It does not
settle the nonzero cubic radial term. Known parity subfamilies remain
available independently.

## Evidence and replay

Run the [focused checker](verify_gaussian_quadratic_radial_affine_coefficient_finite_images.py)
with `--output gaussian_quadratic_radial_affine_coefficient_finite_images_verification.json`
from this directory. It verifies the exact full-source connection and
flatness, bounded original Gamma/Gaussian controls, sector constants,
Gaussian jump and deck signs, leading-coefficient controls and finite
original-marker formulas. Its
[receipt](gaussian_quadratic_radial_affine_coefficient_finite_images_verification.json)
pins this proof and all indispensable retained imports, including the
primary manifest and its exact cached sources. The unrestricted
classification is the written proof, not a finite moment prefix or a
finite set of sampled multiplicities. Existing internal-review limits
remain; no independent replay, formal verification or external review
is asserted.
