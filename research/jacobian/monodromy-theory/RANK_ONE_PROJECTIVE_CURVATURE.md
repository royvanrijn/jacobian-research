# Nontriangular quadratic projective curvature of rank one

The triangular result extends first to every quadratic triangular
projective presentation, without a constant-diagonal assumption.
A separate negative-cofactor argument and a universal curl coefficient
then close a genuinely nontriangular rank4 class: the entire quadratic
projective nonlinearity has one-dimensional image.

These are homogeneous input-construction theorems. They do not exclude
exceptional11 absolutely or complete the faithful-transitive spectrum.

## 1. All quadratic triangular projective presentations

**KMFI130.** Let r>=2, p=(u,x1,...,x_(r-1)), and V a homogeneous
quadratic tuple with V0=u^2. Suppose v_i=V_i(1,x) depends only on
x1,...,xi. If H=g-fp=FV with any nonzero homogeneous polynomial F
is normalized and has full rankDX=r, then v is affine, its constant
matrix has distinct nonzero spectrum different from1, and V=uMp.
Thus every full-rank survivor has primitive degree1 and is classified
by KMFI121--122. No constant diagonal was assumed here.

Indeed KMFI125's determinant proof, also valid without primitivity,
gives f=-Fh, h linear and detDV=constant*h^r nonzero.
The first derivative row of V is (2u,0,...,0), so u divides detDV.
Hence h=lambda*u. On u=1 the determinant of B=Dv is a nonzero constant.
Triangularity gives detB=product partial_i(v_i), so every factor is
a nonzero constant. Thus v_i=a_i*x_i+Q_i(earlier), and KMFI129 applies.

This does not extend to every homogeneous degree by the same argument:
for d>=3 the polynomial h of degree d-1 can have other factors.
Nor does it make every quadratic v triangular.

## 2. A negative logarithmic cofactor cannot be polynomial

Let M,b,c be constant matrix/vector data and let
\[
                      W=Mx+b+cQ(x)-x.
\]

**KMFI131.** If Q is a homogeneous quadratic and its directional
derivative D_c Q is nonzero, there is no nonzero polynomial A satisfying
\[
                     W(A)=(\kappa-D_c Q)A               \tag{1}
\]
for any constant kappa. Also, if Q=R(ell*x) for a univariate polynomial
R of degree>=2 and ell*c!=0, the same nonexistence holds with
D_c Q=(ell*c)R'(ell*x), even in higher degree. No invertibility or
triangularity of M is assumed in this differential lemma.

Here are proofs including the isotropic quadratic case.
After a constant linear change put c=e1.
If Q(c)!=0, write Q=alpha*x1^2+x1*L(other)+Q0(other), alpha!=0.
Let m be the degree of A in x1, with nonzero leading coefficient
A_m(other). Comparing degree m+1 in(1) gives
\[
 D A_m=-\alpha(m+2)A_m,\qquad
                   D=\sum_{i>1}M_{i1}\partial_i .
\]
D is a constant directional derivative, hence locally nilpotent on
polynomials. It has no nonzero constant eigenvalue, a contradiction.

If Q(c)=0 but D_cQ!=0, choose x2=D_cQ and make a linear change
of x1 preserving c so that
\[
                         Q=x1*x2+Q0(x3,\ldots,xn).
\]
This is completion of the nondegenerate cross term; Q0 is independent
of x1,x2. Comparing the highest degree in x2 in(1) gives
\[
 \left[(x1+M_{12})\partial_1+
                    \sum_{i>2}M_{i2}\partial_i\right]A_m=-A_m.
\]
After shifting x1+M12, this operator is an Euler derivation plus
a commuting constant directional derivative. Its polynomial eigenvalues
are nonnegative integers: compare the highest Euler-variable coefficient
and use local nilpotence of the other summand. Eigenvalue -1 is impossible.

For the univariate case conjugate so ell*x=x1 and c is a nonzero
multiple of e1. If degR=q>2, the highest x1-degree gives m=-q,
impossible for a polynomial. If q=2, the first argument applies;
linear and constant terms are lower in x1 and do not affect it.

To apply the lemma to homogeneous normalized inputs with V0=u^d,
suppose detDv is a nonzero constant. Then the determinant gate gives
h=lambda*u^(d-1). If s=degF, the trace gate on u=1 is
\[
 (v-x)(Fbar)=
       [r\lambda-d-\operatorname{tr}M-s-D_cQ]Fbar.        \tag{2}
\]
Thus(1) is a necessary condition with a constant kappa.
Any nonlinearity covered by KMFI131 therefore admits no nonzero
homogeneous polynomial normalized multiplier at all, before full rank.

## 3. The complete rank4 quadratic rank-one image class

**KMFI132.** Let p=(u,x1,x2,x3) and let V be homogeneous quadratic,
V0=u^2, with projective dehomogenization
\[
                    v=Mx+b+cQ(x),                       \tag{3}
\]
where Q is a nonzero homogeneous quadratic form and c!=0.
Then no nonzero homogeneous polynomial F makes H=FV a normalized
input of full rank4. The constants M,b,c and the quadratic form Q
are arbitrary initially. The quadratic part has one-dimensional image;
Q is not assumed to be a square of a linear form.

The same exclusion holds for any higher displayed homogeneous degree
obtained by multiplying V by a power of u. Such a representation merely
absorbs that power into F. Affine projective inputs are the separately
classified survivors; no general quadratic nonlinearity of image
dimension>=2 is classified here.

Normalization first forces detB constant nonzero, exactly as in KMFI130.
Thus M is invertible, Delta=detM!=0, and the determinant lemma gives
\[
                  dQ(M^{-1}c)=0                         \tag{4}
\]
identically. For a quadratic form this says M^-1 c is in its radical.
In particular rankQ<=2.

If D_cQ!=0, KMFI131 already excludes every normalized polynomial F.
It remains c in the radical of Q.

### Rank-two quadratic form

If rankQ=2, its radical is one-dimensional.
Equation(4) makes M^-1 c proportional to c, so kc is M-invariant.
Choose a basis with c last and triangularize the induced two-dimensional
quotient matrix. Then M is lower triangular, and Q depends only on
the two earlier variables. Constants b are allowed in the triangular
architecture. KMFI129 excludes full rank.

### Rank-one quadratic form with an invariant subspace

If rankQ=1, absorb a constant into c and write Q=(ell*x)^2.
Here ell*c=ell*M^-1*c=0.
If also ell*M*c=0, Cayley--Hamilton gives ell*M^2*c=0.
Consequently the invariant cyclic subspace
K=span(c,Mc,M^2c) lies in kerell and has dimension1 or2.

If dimK=1, put c last and triangularize M on its two-dimensional
quotient; Q depends only on the earlier variables.
If dimK=2, K=kerell. Put this invariant plane last, refine it by
an invariant line, and triangularize M there. Q depends only on the
first variable while c lies in the last plane.
Both are exactly the constant-diagonal triangular projective class
KMFI129. The original nonzero quadratic part cannot become affine
under a constant affine conjugacy, so full rank is impossible.

### The genuinely nontriangular companion case

The remaining case is
ell*c=ell*M^-1*c=0 but a=ell*M*c!=0.
Set tau=trM, sigma=the second characteristic coefficient, Delta=detM.
The three forms
\[
            X=ell*x,\quad Y=ell*M*x,\quad
                            Z=Delta*ell*M^{-1}*x
\]
are independent. To see this, applying a relation to c kills the X,Z
terms but not Y. Proportionality of X and Z would make ell a left
eigenvector of M, contradicting ell*M*c!=0.
Cayley--Hamilton now gives the linear normal form
\[
                   v=(Y,\ Z-\sigma X+\tau Y+aX^2,\ \Delta X)+b'.
\]
A common scaling of X,Y,Z normalizes a to1. The constant vector b'
can also be removed by affine projective conjugacy.
Indeed W=v-x has a zero: its first and third equations express Y,Z
linearly in X, and the middle equation is a genuine quadratic in X.
Over algebraically closed constants it has a root. Translation by this
fixed point removes b'; it only changes sigma in the displayed formula.
Thus it suffices to treat
\[
                 v=(y,\ z-\sigma x+\tau y+x^2,\ \Delta x),
                 \qquad \Delta\ne0.                      \tag{5}
\]

## 4. One universal curl coefficient closes the companion case

For (5), the full quadratic tuple and its chart derivative are
\[
 V=(u^2,\ uy,\ uz-\sigma ux+\tau uy+x^2,\ \Delta ux),
\]
\[
 A=DV|_{u=1}=
 \begin{pmatrix}
 2&0&0&0\\
 y&0&1&0\\
 z-\sigma x+\tau y&2x-\sigma&\tau&1\\
 \Delta x&\Delta&0&0
 \end{pmatrix},\quad
 w=V|_{u=1}=(1,y,z-\sigma x+\tau y+x^2,\Delta x).
\]
Its characteristic polynomial is
\[
             \chi_A=(T-2)
                  [T^3-\tau T^2+(\sigma-2x)T-\Delta].
\]
Let C=(w,Aw,A^2w,A^3w) and D=detC/Delta.
D is polynomial, and
\[
                             [x^5z]D=-2.                 \tag{6}
\]
Hence C is invertible for every parameter specialization Delta!=0,
including every repeated or resonant spectrum.

Write mu_0,...,mu_3 for the first four coefficients at T infinity of
\[
                         1-(T-\lambda)^4/\chi_A.
\]
The unique feedback one-form is eta=(mu_0,...,mu_3)C^-1.
These are exactly KMFI125's equations, with h=lambda*u.
In a polynomial multiplier it must be dlogF; in particular its
transverse restriction must be closed. But the exact identity is
\[
 [x^2yz^4]\left\{
      D^2(\partial_x\eta_y-\partial_y\eta_x)\right\}=8.    \tag{7}
\]
The expression in braces is polynomial. The coefficient is independent
of lambda,tau,sigma,Delta and never vanishes in characteristic zero.
Thus eta cannot be closed for any parameters. No normalized polynomial
multiplier exists in (5), irrespective of rank.

Equations(6)--(7) are finite symbolic polynomial certificates, not a
sample search. The determinant and moment definitions specify every
coefficient without an omitted computation. The exact checker expands
the whole identities over Q(lambda,tau,sigma,Delta), retains their
polynomial models, and verifies both constant coefficients.
Since D has the nonzero coefficient(6), parameter specialization cannot
invalidate the denominator gate or conceal a solution. All scalar
degrees were already eliminated by the unique feedback form.

The three cases exhaust(3), proving KMFI132.

## 5. Evidence and remaining frontier

[Helper](rank_one_projective.py) ·
[checker](verify_rank_one_projective.py) ·
[model](rank_one_projective_model.json) ·
[receipt](rank_one_projective_receipt.json).

The proof is unbounded in the homogeneous multiplier degree.
The negative-cofactor lemma also handles univariate nonlinearity in
every degree and dimension. Exact controls check the universal curl,
normal-form conjugacy, triangular determinant gate, negative-cofactor
identities, and the nonlinear rank-deficient triangular survivors.
The latter show why KMFI132 is a full-rank obstruction for the entire
class, although its two nontriangular mechanisms exclude normalization
itself. No blanket absence of rank-deficient multipliers is asserted.

This closes a genuinely nontriangular primitive quadratic rank4 class,
including arbitrary linear mixing, constants, square and nonsquare
quadratic forms of one-dimensional nonlinear image.
Remaining homogeneous quadratic directions can have nonlinear image
dimension>=2, or lack a pure-square first coordinate after the allowed
matching linear changes. Higher primitive degrees and nonhomogeneous
fields remain open. None of the theorems is an absolute PSL2(11)
obstruction or a new abstract occurrence group. The full monodromy
classification remains the active unresolved objective.

Formal verification, independent whole-proof replay and external
review are not claimed.
