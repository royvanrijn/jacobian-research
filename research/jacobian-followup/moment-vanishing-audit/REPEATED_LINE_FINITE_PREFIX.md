# An explicit pure prefix for every repeated-line leading ternary cubic

15 November 2026. Internal written quantitative corollary of the
[repeated-line exclusion](TERNARY_CUBIC_GMC_CLOSURE.md). Exact auxiliary
connection-degree checks are separate from its analytic period input.
No independent specialist replay, external review or formal verification.

**G3REPEATEDFINITE20261115.** For any complex ternary cubic whose
leading form is L^2*M or L^3, the first99 zero pure Gaussian moments
force all positive pure moments to vanish. Consequently mixed moments
vanish for m>3*deg Q, sharply by the retained cubic shear. The first99 zero pure moments
force a Gaussian W-linear frame. Every branch pure-excluded by the
repeated-line proof has a nonzero pure moment by99. The pure gate99 is
not claimed optimal or sufficient for arbitrary cubics.

## A common rational degree bound for both isotropy charts

In every pure-exclusion branch of that proof, exact quadratic elimination
leaves a degree d=3 or4 polynomial phase and determinant prefactor

    I_Gamma(t)=constant*D(t)^(-1/2)*integral exp(U(t,u)/D(t)) du,
    deg_t D<=2, deg_t U<=3, deg_u U=d.           (1)

For a nonisotropic repeated form, D=det(I-2tA). The inverse2-by2 matrix
has adjugate of degree<=1 in t, so multiplying the phase
t*g-u^2/2+(t^2/2)*B^T*(I-2tA)^(-1)*B by D has degree<=3 in t.
The prefactor is exactly the determinant in(1).

For an isotropic repeated form with nonzero W^2 coefficient lambda,
put C_T=A*t-1/2 after completing W. Completing T gives D=t*C_T,
up to a nonzero constant. Before T completion the phase has t,1,1/t
terms, its T-linear coefficient has t-degree<=1, and C_T has degree<=1.
The resulting numerator on the common denominator t*C_T has degree<=3.
The two Gaussian determinant factors combine to a constant times D^-1/2.
All coefficient cancellations and lower-degree D are allowed. D is
nonzero; neither D(0)!=0 nor convergence of a Wick series is assumed.
The eliminated phase has actual degree3 or4 in every exclusion branch;
the zero-quartic/no-cubic boundaries are handled by the independent frame
theorem rather than this exclusion calculation.

## The exact polynomial connection bound

Let a(t)=d*[u^d]U!=0. The top one-variable twisted quotient has basis
1,u,...,u^(d-2). For n>=d-1 its integration-by-parts relation is

    a*u^n = -sum_(k=0)^(d-2) ([u^k]U_u)*u^(n-d+1+k)
            -(n-d+1)*D*u^(n-d).                 (2)

The normalized connection on a basis vector has source

    C=(U_t*D-U*D')/D^2-D'/(2D).                 (3)

Its numerator has t-degree<=4 and u-degree<=d. Multiplying a basis
vector by it gives u-degree at most2d-2. Repeated use of(2) lowers
u-degree strictly, and needs at most d reductions. Each introduces
one denominator a; its other numerator coefficients have degree<=3
(the derivative term contains D of degree<=2). Therefore the full
connection matrix can be written N/H with

    H=D^2*a^d, B=4+3d,
    deg H<=B, deg N_ij<=B.                      (4)

These are upper bounds. Common factors, sparse sources, vanishing
columns, degree drops and special coefficients can only reduce them.

## A seed relation of bounded polynomial degree

At minimal seed rank s<=d-1, write the iterated connection columns

    v_j=h_j/H^j, deg h_j<=j*B,  j=0,...,s.

Indeed h0 is the constant seed and

    h_(j+1)=H*h_j'-j*H'*h_j+N*h_j.              (5)

Choose s independent rows of these columns at minimal rank. The signed
maximal minors c_j of the h-column matrix give a nonzero relation;
each has degree<=B*(s*(s+1)/2-j). Setting q_j=c_j*H^j gives

    sum_(j=0)^s q_j(t)*I_Gamma^(j)(t)=0,
    deg q_j<=C=B*s*(s+1)/2.                    (6)

The chosen independent rows determine the full column relation since
the seed rank is s. Minimal rank ensures q_s!=0. In every pure-exclusion
branch q0!=0 as well: otherwise(6) is the minimal seed equation with
the constant function as a solution. The retained seed-submodule/actual
period bridge would then give an actual nonzero constant seed period,
which the exhaustive repeated-line proof excludes. This is where that
written analytic input enters; the polynomial degree calculation alone
does not prove q0!=0. Seed cyclicity in the full ambient module is not
required.

Apply (6) to the original formal Wick series I=E(exp(tP)), which is a
horizontal functional by Gaussian integration by parts before elimination.
If M1..M_K=0, then I=1+O(t^(K+1)). For K=C+s, every term q_j*I^(j),
j>=1, has order at least C+1, while q0*I=q0+O(t^(K+1)). A nonzero
polynomial q0 of degree<=C cannot be canceled. Thus some pure moment
by C+s is nonzero. This argument is valid even when D(0)=0, because
the final relation has polynomial coefficients and the original Wick
series is formal; no eliminated t-adic contour integral is used.

The worst values are

    d=3: B=13, s<=2, C+s<=41;
    d=4: B=16, s<=3, C+s<=99.                  (7)

The nonisotropic double-line zero-quartic branch and isotropic zero-W^2
branch already give frames directly. Every other repeated-line branch
is pure-excluded by99, by the retained exhaustive case cover. Hence
99 initial zero moments force a frame. The specific frames supplied
here are more restricted than an arbitrary W-linear frame, which lets
us retain99 for the whole repeated-line class.

## The surviving frames require at most93 moments

In the nonisotropic double-line zero-quartic branch, take T along the
repeated form and W along the constructed transverse isotropic v.
The leading cubic is T^2*(alpha*T+beta*Z), independent of W. Thus P
belongs to W*(a+bZ+cT)+G(Z,T), the d=0 face of the retained
[fourteen-coefficient theorem](MIXED_PERIOD_AT_INFINITY.md), with
sufficient pure gate92.

In the isotropic zero-W^2 branch take L=Z. The leading cubic is
Z^2*M. Its W coefficient is delta*Z^2. If delta=0 the same d=0
face applies, regardless of G3. If delta!=0, a Gaussian null rotation
fixing Z removes the T coefficient from M, leaving

    P=W*(a+bZ+cT+dZ^2)+alpha*Z^3+G2(Z,T), d!=0, deg G2<=2.

When c=0 this is the retained c=0 face with sufficient pure gate35.
When c!=0, the retained full mixed-phase theorem pure-excludes it.
Its quantitative degree bound improves here because the actual
one-variable phase has degree4 rather than possibly5/6: put

    h(Z,t)=(Z-t*(a+bZ+dZ^2))/(c*t),
    Phi=t*(alpha*Z^3+G2(Z,h))-h^2/2,
    prefactor=constant/t.

The quartic coefficient is exactly d^2*(t*[T^2]G2-1/2)/c^2, never
identically zero, and t^2*Phi has t-degree<=3. The retained mixed
connection bound is B=3d_phase+3=15, seed rank<=3, giving the explicit
prefix15*3*4/2+3=93. This uses the complete retained interior exclusion
to ensure a nonzero zeroth scalar coefficient; it does not infer
exclusion from degree alone. The c=0 and d=0 bounds are35 and92.

Thus every frame branch either is pure-excluded by93 or has all-order
nullity forced by at most92 moments. Combined with the complementary99
bound, this proves the class-wide pure gate99. The retained W-linear
mixed theorem gives m>3*deg Q. Its sharp shear W+ZT-Z^3/2 has repeated
leading form Z^3, so the coefficient3 is also sharp within this class.

## Evidence and remaining boundary

The [checker](verify_repeated_line_finite_prefix.py) and
[receipt](repeated_line_finite_prefix_verification.json) verify the common
numerator degrees of both eliminations, actual connection reductions on
degree-complete rational models including D(0)=0, the common-denominator
and iterated-column degree bounds, the restricted mixed quartic identity
and the35/92/93/99 arithmetic. They do not certify analytic period
completeness or evaluate moments through99.
The [source audit](repeated_line_finite_prefix_source_audit.json) records
these distinctions.

The [least global pure prefix](CUBIC_FRAME_AND_PREFIX_COMPUTABILITY.md)
is now proved computable but numerically uncomputed. This note supplies
an explicit99 gate on the nonreduced-leading class; the arbitrary squarefree
affine class still needs a quantitative bound before any numerical global
pure-only prefix is asserted. Scheme-theoretic ideal cutoffs, term-count
optimality, independent review and lower-power GVC remain separate.
