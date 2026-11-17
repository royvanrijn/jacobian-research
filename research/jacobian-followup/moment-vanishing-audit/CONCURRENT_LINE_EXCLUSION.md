# Three concurrent distinct leading lines: complete affine GMC safety

The subsequent [weighted binary finite-prefix proof](SQUAREFREE_SINGULAR_FINITE_PREFIXES.md)
supplies the sufficient pure gate concurrent104. The qualitative
arguments below retain their scopes; earlier absence-of-prefix statements
are historical. No numerical arbitrary-cubic bound follows yet.

14 November 2026. Internal written formal proof with exact auxiliary
algebra. No independent replay, external review or formal verification.

## Statement and boundary

**G3CONCURRENT20261114.** Every ternary cubic whose leading form is
a product of three distinct concurrent linear factors satisfies GMC,
with arbitrary lower terms. If their common vertex is nonisotropic,
some positive pure moment is nonzero. The same pure exclusion holds
when F2(vertex)!=0, regardless of isotropy. The isotropic zero-F2 case
uses the retained W-linear GMC theorem, without a separate assertion
here that its entire pure-null locus is empty.

Independent leading triangles were closed by the
[affine triangle theorem](AFFINE_TRIANGLE_EXCLUSION.md). The remaining
possible ternary cubic GMC failures therefore have a line-smooth-conic
leading curve or a repeated leading line. Neither an explicit global
pure prefix nor the coefficient-scheme/lower-power GVC questions are
settled by this note.

The subsequent [affine line-conic theorem and squarefree-leading
corollary](LINE_CONIC_AFFINE_EXCLUSION.md) close every reduced leading
ternary cubic. Only repeated leading lines remain as possible cubic
GMC failures; the earlier arguments below retain their own scopes.

## 1. Vertex reduction and an exact binary phase

Choose coordinates in which the vertex is the T-axis and F3=F(u,v)
is a squarefree homogeneous binary cubic. Its binary gradient has
only the zero at the origin. If F2(vertex)!=0, any sequence violating
complex linear gradient coercivity must approach that vertex, so
u,v=o(T); then P_T=2*alpha*T+O(u,v,1), alpha!=0, contradicts
grad P=o(||x||). The [coercive theorem](COERCIVE_GRADIENT_MOMENTS.md)
pure-excludes this case.

For an isotropic vertex with F2(vertex)=0, the
[W-linear Gaussian frame criterion](RANK_ONE_WT2_CLASSIFICATION.md)
gives GMC safety. It includes the full affine polynomial, not only F3.

Otherwise retain the full Gaussian with m=Q_TT!=0 and det Hess Q!=0:

    P=T*H(u,v)+G(u,v), deg H<=1, G3=F,
    Q=m*T^2/2+T*J(u,v)+q2(u,v), J=Du+Ev+l_T.       (1)

Here q2 includes all transverse quadratic, linear and constant terms.
The parameter-dependent polynomial change T=T'+(tH-J)/m is exact:

    tP-Q=-m*T'^2/2+Phi_t(u,v),
    Phi_t=tG-q2+(tH-J)^2/(2m).                    (2)

It has volume Jacobian1. The Gaussian factor has rank1; the remaining
binary cubic has squarefree leading part tF and top twisted rank4.
Indeed its two quadratic leading gradients have no common projective
zero and form a regular sequence of total length4. Equivalently the
original projective gradient quadrics have a unique length4 vertex
at infinity: the T equation forces the homogenizing coordinate to0
by its nonzero -m coefficient, then the binary leading gradients give4.
Bezout leaves8-4=4. The same vertex argument gives complex gradient
coercivity for each fixed nonzero t, so the retained tame comparison
also gives full twisted rank4.

Use the [formal seed criterion](CUSPIDAL_AFFINE_EXCLUSION.md), Section2.
All pure zeros give a rational horizontal quotient of the seed submodule
with seed value1. A complete semisimple formal frame after ramification
lets that quotient extend to the full module. Nonzero exponential factors
cannot map to the trivial factor; positive valuation of every remaining
zero-exponential scalar seed then contradicts its required value1.

## 2. Constant H

For H=h, (2) becomes

    Phi_t=h^2*t^2/(2m)+t*R-q_eff,
    R=G-(h/m)*J, q_eff=q2-J^2/(2m).               (3)

The Schur Hessian is nondegenerate: det Hess Q=m*det Hess q_eff.
R retains squarefree binary cubic leading part. It has finite critical
Milnor rank4 and no escapes under R-q_eff/t. At a critical point whose
Hessian has rank at least1, the analytic Morse lemma gives the local
type alpha+x^2+y^(mu+1), mu<=4. Formal Gaussian/power recurrences give
scalar characters exp(alpha*t)*t^(-1/2-j/(mu+1)),1<=j<=mu, all with
strictly positive prefactor powers. Finite Hermite jets give a complete
frame when there are several such points.

If the Hessian has rank0, the cubic Taylor expansion is exactly alpha+F
centered there, of Milnor rank4, exhausting the full rank. Rescaling both
variables by t^(-1/3) gives a fixed homogeneous rank4 cubic phase, volume
t^(-2/3) and a Gaussian perturbation of positive parameter order. Its
leading quotient has Hilbert function1,2,1; uv(u+v), for example, has
basis1,u,v,v^2. The fixed-phase conjugacy ell(exp(F_r-F0)*a) gives four
independent formal functionals with positive scalar seed valuation.
This also proves completeness of the rank0-Hessian alternative.

Thus h=0 has only positive-valued zero-exponential seeds, and h!=0
has a common nonzero quadratic action h^2*t^2/(2m) on every factor.
The criterion proves pure exclusion in both cases. The normalized
conditional Gaussian integration in (2) changes the seed only by a
nonzero parameter-independent constant, never by a hidden t power.

## 3. Nonconstant H: two exhaustive normal forms

Take u=H as an affine coordinate. Absorb all degree<=2 terms of G
divisible by u by an affine T change. A linear v shift by a multiple
of u depresses the leading cubic, after which lower divisible terms
are again absorbed. The full Q changes and is retained, with m!=0.
If u does not divide F, the result is

    P=uT+a*u^3+b*u^2*v+d*v^3+q*v^2+fv+k,
    d!=0,27*a^2*d+4*b^3!=0.                     (4)

If u divides F, the result is

    P=uT+a*u^3+c*u*v^2+q*v^2+fv+k, a*c!=0.      (5)

In the latter case squarefreeness forces the uv^2 coefficient nonzero;
depressing its quadratic cofactor makes a nonzero as well. These cases
exhaust all squarefree F. In (2), set u=tU,v=tV. The leading binary
phase is t^4*Phi0, with Phi0=F(U,V)+U^2/(2m).

## 4. The nonfactor case (4)

Finite P-critical points have u=T=0 and v critical for
d*v^3+q*v^2+fv+k. They have total Milnor rank2, including a double
critical point. Their three-variable types are s*w+y^D,D=2 or3;
the zero-exponential scalar powers t^(-1-j/D) are positive.

For a nonzero critical point of Phi0 put V=zU. Then

    R(z)=a+bz+d*z^3, 3d*z^2+b=0,
    U=-1/(3m*R(z)), V=zU.                        (6)

The resultant is d^2*(27a^2*d+4b^3)!=0. Hence R does not vanish at
these roots, and their critical values are U^2/(6m)!=0. If b!=0,
both z roots are nonzero and simple, with Hessian determinant
-6d*z*U/m!=0. They give two escaping Morse characters of prefactor
t^-2 and nonzero quartic action. Their physical Gaussian fluctuations
are bounded; the T centers are distinct at order t^2. The
[confluent numerator frame](AFFINE_TRIANGLE_EXCLUSION.md), of degree<4
here, joins them to the finite rank2. This is a full formal frame,
even if the action values coincide. The seed criterion applies.

### The entire b=0 degenerate escape

Now a*d!=0. Phi0 has a nonzero A2 point
(U,V)=(-1/(3am),0), with value C0=1/(54a^2*m^3)!=0. Put r=1/t.
A formal parameter-dependent coordinate change near that point gives
the scaled binary phase exactly as

    C(r)+X^2+Y^3+A(r)*Y, C(0)=C0,A(0)=0.         (7)

Construct this recursively: split the nondegenerate X direction;
at each r order divide the remaining analytic Y coefficient by Y^2,
remove the divisible part by a coordinate vector field, retain constant
and linear terms, and normalize the leading cubic unit. This formal
division works to all orders. At r=0, Y is a nonzero constant multiple
of V, and the radial coordinate depends only on U.

The exponent is r^-4 times (7). If ell=ord A is1 or2, its two roots
Y~r^(ell/2) yield Morse characters. Their width is r^(2-ell/4), and
width/separation has positive power2-3ell/4. Cubic corrections shrink.
If ell>=3 or A=0, use X=r^2*x,Y=r^(4/3)*y. The linear perturbation
has positive power ell-8/3, so a fixed Gaussian/cubic quotient of
rank2 gives two formal characters. This exhausts every order, including
exact vanishing of A. Each factor retains C0*t^4 as a nonzero leading
action; a lower unfolding term cannot cancel it.

These modes and the finite characters form a complete rank4 frame.
For an explicit leading matrix, multiply the finite jet rows, of
degree<2, by [1-T/(K0*t^2)]^8,K0=-1/(3a*m^2). Their growth at escape
is at most t^2 and the cutoff is O(t^-8). For the escape rows use
[T/(K0*t^2)]^8 times1,(V-beta(r))/sigma, where beta=O(r), and sigma
is r^(ell/2) in the Morse case or r^(4/3) in the cubic case. The local
coordinate Y differs from this leading polynomial coordinate by terms
of higher valuation in these charts. At finite points its worst growth
is r^-1/3, dominated by the reverse cutoff r^16. At escape the two-row
leading matrix is respectively the two-root Gaussian matrix or the
fixed cubic dual basis. Cross blocks vanish, proving completeness.
This closes every degenerate escape in (4).

## 5. The factor case (5)

There is always one primary escaping Morse point of Phi0:

    U=-1/(3am),V=0,
    Phi0=1/(54a^2*m^3), det Hess Phi0=2c/(3a*m^2).

Its prefactor is t^-2 and its quartic action is nonzero. The other
three dimensions split exhaustively as follows.

**q!=0.** There is one finite Morse point
(0,-f/(2q),-c*f^2/(4q^2)), of value k-f^2/(4q). Two secondary escapes
have u=U+o(1),v=t^(1/2)*(V+o(1)),T=t*(W+o(1)), with

    U=-q/c,W=-q/(cm), V^2=q/(mc^2).              (8)

Their limiting critical system is W+cV^2=0,2(cU+q)V=0,U-mW=0,
with Jacobian4q!=0. Their action has leading q^2*t^2/(2mc^2)!=0.
Fluctuation scales t^-1,t^(-1/2),1 in u,v,T give Gaussian Hessian
[[0,2cV,1],[2cV,0,0],[1,0,-m]], determinant4q and prefactor t^(-3/2).
Higher cubic terms shrink. A generic physical linear form separates
the two secondary points by at least t^(1/2), and the primary by t^2;
all fluctuations remain bounded. The numerator frame gives all four
characters: finite1,secondary2,primary1. Any zero-exponential finite
seed has positive valuation.

**q=0,f!=0.** There are no finite P-critical points. Three secondary
escapes have u=t^(-1/3) U,v=t^(1/3) V,T=t^(2/3) W, where

    V^3=f/(2mc^2), U=-mcV^2,W=-cV^2.             (9)

Their critical Jacobian and limiting Gaussian Hessian determinant are
6mc^2*V^2!=0. Their leading action is (3fV/4)*t^(4/3)!=0, prefactor
t^(-4/3), and fluctuation scales t^-1,t^(-1/3),1. The three V^2 values
are distinct, giving T separation of order t^(2/3). The numerator frame
combines these with the primary into rank4. All exponential factors
are nonzero.

**q=f=0.** The critical set is the parabola u=0,T=-cv^2, with value k.
Put s=T+cv^2,u=rU,s=rS,r=t^(-1/2). Exactly,

    t(P-k)=US+a*r*U^3,
    Q(rU,v,-cv^2+rS)=Q_curve(v)+O(r).

The fixed phase US-Q_curve(v) has a Gaussian pair and a quartic of
leading coefficient -mc^2/2!=0. Its quotient has rank3 with basis1,v,v^2,
including every critical degeneracy. Fixed-phase perturbation yields
three formal characters with seed exp(k*t)*t^-1*C[[r]]. Together with
the primary this is a full rank4 frame: use the same eighth-power T
cutoffs on the three parabola rows and the primary row. At the primary
v=O(1); at the parabola T=-cv^2+O(r) is a fixed polynomial in the block
variables. Cross terms have positive parameter order, and the diagonal
blocks are the fixed quartic dual basis and a nonzero Gaussian. If k=0
all three zero-exponential scalar seeds have positive valuation;
otherwise their exponential is nonzero. This finishes (5).

## Evidence and remaining scope

The [checker](verify_concurrent_line_profiles.py) and
[receipt](concurrent_line_profiles_verification.json) verify Gaussian
completion, the Schur determinant, binary Jacobian basis, resultant,
actions/Hessians, secondary scalings, quartic restriction, parabola
normal scaling and A2 threshold arithmetic. They do not certify formal
normal forms, completeness or quotient extension. The
[source audit](concurrent_line_source_audit.json) records those boundaries.

The first checker missed one factor d in the resultant; its nonzero
locus for d!=0 was unchanged. The second compared expanded and factored
gradients structurally. Both exact inputs and terminal failures are
retained. The corrected checker passes; no failed run is promoted.

All concurrent-leading cubics are now GMC-safe. Nonhomogeneous
line-smooth-conic and repeated-line leading curves, numerical uniform
gates, coefficient schemes, lower-power GVC and external review remain open.
