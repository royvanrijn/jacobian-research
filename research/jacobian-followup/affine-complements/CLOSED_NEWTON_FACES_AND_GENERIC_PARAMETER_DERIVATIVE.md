# Closed Newton faces exclude mixed quadratic cyclic target coordinates

PF2RULED148. 17 November 2026. ALL degrees over C.

On EVERY quadratic MX0 surface, there is no rational Jacobian mate for

    P=a*w+F(u,z)+g(u), a in C*,
    F=sum c_(j,i)*u^j*z^i, 1<=i<=n, j>=0, n>=2,
    c_(0,n)!=0,
    (2*n-1)*j+m*i <= n*m for EVERY term of F,       (1)

where m=deg M and g is an ARBITRARY polynomial. Coefficients on the
CLOSED Newton face are unrestricted. No simple-root, double-root,
nondegeneracy, target-symmetry or companion-degree condition is imposed.
The exclusion also holds in every finite source-field extension.

This includes genuine products in u,z, not just the displayed additive
family of [PF147](ADDITIVE_BASE_PERTURBATIONS_AND_RESONANT_CYCLIC_TARGET_EXCLUSION.md).
It does not reduce arbitrary coordinates to (1). Full first-derivative
integrality, unrestricted cyclic targets and HC4/JC2 remain OPEN.

## 1. Scope, actual coframe and retained inputs

More generally let

    V={w^2=N(u)+M(u)*z} smooth,
    m=deg M >=3, 1<=ell=deg N<m,
    Omega=du wedge dw/M.

The residue form is a regular unit on ALL V. We prove the assertion
above for V. The leading z coefficient in (1) is a nonzero CONSTANT.
Every term with i<n has u-degree at most

    jmax <= (n-1)*m/(2*n-1) < m/2.                 (2)

An absent mixed term contributes zero, not a required degree.
Pure base terms are all placed in g; thus F(u,0)=0.

On the WHOLE geometric generic curve P=p,

    G=N+M*z-(p-F-g)^2/a^2=0,
    G_u=N'+M'*z+(2*w/a)*(F_u+g'),
    G_z=M+(2*w/a)*F_z,
    Omega=dP wedge omega,
    omega=dz/(a*G_u)=-du/(a*G_z).                  (3)

Both mixed partial derivatives are actual partials of F. Omitting F_u
or using PF147's additive derivative would be incorrect.

As in PF147, this generic curve is smooth, componentwise if necessary:
on a smooth characteristic-zero surface each critical curve of P maps
to a single value, and the remaining isolated critical points are
finite. Formula (3) is regular and NONZERO on the entire affine curve,
including all M=0 patches. A mate would give dQ=omega on the fibre.

The argument below proves holomorphicness at EVERY infinity valuation.
The symbolic controls are support for this universal proof, not a
bounded census promoted to a theorem.

## 2. Finite z at infinity

A pole of z with u finite is impossible: the leading coefficient
c_(0,n) is constant, so F has a pole of order n*(-nu(z)), while
N+M*z has at most order -nu(z). Thus at every infinity point u has
a pole alpha=-nu(u)>=1.

Suppose z stays finite there. By (2), all mixed coefficients grow
strictly slower than u^(m/2). If d=deg g, with d=0 for constant g,
comparison in w^2=N+Mz gives

* 2*d<m: the residue of z is zero;
* 2*d=m: its residue is g_d^2/(a^2*M_m), nonzero;
* 2*d>m: no such branch exists.

On every finite-z branch the pole order of w is at most
max(d,jmax)*alpha. Consequently w*F_z has pole order strictly smaller
than m*alpha, since d<=m/2 and jmax<m/2. In (3), G_z is therefore led
by M, and

    ord(omega)=(m-1)*alpha-1 >=0.                  (4)

This includes finite NONZERO limiting z, zero limiting z, all leading
coefficient cancellations in N-g^2, and every ramification index.

## 3. Polar z away from the closed face

Now beta=-nu(z)>0. The surface equation gives the EXACT w pole

    W=(m*alpha+beta)/2, H=n*beta, J=d*alpha.        (5)

For every term u^j z^i of F with i<n, its pole order T=j*alpha+i*beta
satisfies

    T <= max(W,H),                                 (6)

with equality possible ONLY when W=H and its term lies on the closed
face in (1). To see this, write r=beta/alpha and t=m/(2*n-1).
For 0<=r<=t compare j+i*r to (m+r)/2: at r=0 it is strictly lower
by (2), and at r=t inequality (1) gives the desired comparison.
For r>=t compare to n*r; i<n makes it strictly lower after t.

If H>W, all other terms of F are strictly below H. The equation P=p
forces J=H>W, so the leading f_n*z^n and g cancel. In G_u, g'
strictly dominates F_u and its product with w dominates N'+M'z.
Thus the EXACT order is

    ord(omega)=W-alpha+(n-1)*beta-1 >=0.           (7)

Here W-alpha is a positive integer because
2*(W-alpha)=(m-2)*alpha+beta>0.

If W>H, equation P=p forces J=W. All mixed terms are below W.
The two possible leading pieces of G_u have combined coefficient

    (m-2*d)*w^2/u,

which is nonzero since beta=(2*d-m)*alpha>0. Equation (4) again
gives the exact differential order.

The remaining case is W=H, with J<=W. Unlike PF147, arbitrary
closed-face terms can make its leading zero TRIPLE or higher.
We keep all of them.

## 4. Generic-parameter derivative lemma

Let k be an algebraic closure of C(p). Suppose A(v,Z) is a formal
power series, with coefficients in C, regular near (0,Z0), with
Z0 in C. Let D>=1, and take ANY normalized branch of

    A(v,Z)=p*v^D, Z(0)=Z0.                        (8)

After a formal unit root and a choice of normalized parameter s,
write v=s^e, Z in k[[s]], with Z-Z0 of positive order.
The derivation partial_p of k extends coefficientwise to k((s)),
fixing s and hence v. It commutes with the substitution in A.
Since every coefficient of A is INDEPENDENT OF p,

    A_Z * partial_p Z = v^D.                     (9)

In particular partial_p Z is nonzero. Its order is at least ONE:
the constant coefficient of Z is Z0 in C, whose derivative is zero;
all other terms have positive order. Therefore

    nu_s(A_Z)=e*D-nu_s(partial_p Z) <= e*D-1.      (10)

This proof covers every branch, every multiplicity and every
ramification index. It does not assume a binomial normal form,
simple critical points, a bound on leading multiplicity, or
independence of the branch coefficients from p. Their possible
p-dependence is exactly what (9) measures.

Only generic-CURVE completed fields are involved. The parameter
derivative is not a source vector field asserted to preserve the
original finite normalization or its trace-dual.

## 5. Every closed-face branch

Put h=gcd(m,2*n-1) and

    l=(2*n-1)/h, k0=m/h, D=n*k0.
    l*m+k0=2*D.

For the valuation calculation make the finite curve-base change

    u=v^(-l), Z=v^k0*z, T=v^D*w,
    M0(v)=v^(l*m)*M(v^-l).

It is only a local cover used to examine each normalized curve
valuation; it is not an intermediate affine plane. In case W=H,
Z,T are units, with nonzero constant residues. The surface relation is

    T^2=M0(v)*Z+v^(2*D)*N(v^-l).

The last term vanishes at v=0. Choose the actual leading T on
EACH branch. Then T(v,Z) is a formal unit, with T_Z=M0/(2*T) a unit.
Every exponent in

    A(v,Z)=a*T(v,Z)
             +sum c_(j,i)*v^(D-l*j-k0*i)*Z^i
             +v^D*g(v^-l)

is nonnegative: use (1) and J<=W. Its coefficients are in C.
The generic equation is exactly (8), and Z0 is a nonzero root
of its p-independent leading function. The leading function is
a nonzero constant times sqrt(Z), plus a polynomial of degree n
and possibly a constant. It is not identically zero, so its roots
are isolated. Multiple leading roots are NOT removed.

Direct pullback of the ORIGINAL coframe (3) gives

    Omega=-l*v^(l*m-D-l-1)*(T_Z/M0)*dv wedge dZ,
    omega=l*v^(l*(m-1)-1)*(T_Z/M0)*dv/A_Z.          (11)

At v=s^e, use (10) to obtain

    ord_s(omega)
      = e*l*(m-1)-1-nu_s(A_Z)
      >= e*(l*(m-1)-D) >=0.                       (12)

The last integer is nonnegative because

    l*(m-1)-D
       = ((n-1)*m-(2*n-1))/h
       = ((n-1)*(m-2)-1)/h >=0                   (13)

for m>=3,n>=2.

This proves holomorphicness on EVERY branch of the finite local
cover. A differential pole downstairs could not disappear: its
pullback order at index E is E*ord+(E-1), still negative if
ord<=-1. Thus (12) proves the same assertion on every original
generic curve branch. No ramified branch, residue or multiplicity
is discarded.

## 6. Nonexactness and the actual cyclic application

The required omega is now a nonzero holomorphic differential on
EVERY complete projective generic component. It cannot be a
rational derivative: a rational primitive with a pole creates a
differential pole, and a regular function on a complete projective
curve is constant.

The same obstruction holds in every finite characteristic-zero
extension of C(V). Normalized FIELD trace would trace a rational
primitive to one downstairs; equivalently the full ramified
pullback of a holomorphic nonzero differential stays holomorphic
and nonzero. This is a rational exactness statement only.

For every quadratic MX0 model,

    M=u^2*b(u)^k,
    N=u*(-b+2*(1-k)*u*b')^2,
    k>=2, b square-free, b(0)!=0, r=deg b>=1,

one has m=k*r+2>=4 and ell=2*r+1<m. The full smooth quotient and
whole-plane chart, including q^*Omega=-2*dx wedge dt, are retained
from [MX0](../hc4-triangular-graphs/LOGARITHMIC_DESCENT_AND_CYCLIC_EXACT_SOURCE_IMAGES.md)
and [PF146](HOLOMORPHIC_GENERIC_FIBRES_AND_NONLINEAR_CYCLIC_TARGET_EXCLUSION.md).
Thus every P in (1), with arbitrary g, is excluded as an ACTUAL
target coordinate even if its hypothetical companion only lives
in the whole plane source field.

For example, z^2+u*z is allowed on EVERY quadratic MX0 model:
3+m<=2*m since m>=4. A cubic z^3 permits all mixed monomials
5*j+m*i<=3*m. The equality cases are included.

On the marked plane chart, every fixed-t line passes the old
necessary degree test: z is quadratic in x when a(t)!=0, and
linear with nonzero coefficient when a(t)=0; w is affine or
constant, respectively. Since c_(0,n)!=0, the degree of P there
is 2*n or n, both >=2. This asserts no test of all source curves
under all polynomial changes of coordinates.

The theorem survives surface automorphisms of these cyclic models
and polynomial target rechartings, just as PF147 does. It neither
assumes nor proves that every target coordinate has the displayed
form in some such coordinates.

## 7. A genuine triple-root MX control

Use k=4,b=u-1, so

    M=u^2*(u-1)^4, N=u*(1-7*u)^2, m=6,ell=3,
    P=-8*w+z^2+(-6*u^2+8*u-4/3)*z+g(u),
    g=-3*u^4+8*u^3-(20/3)*u^2+(32/27)*u+32/81.  (14)

Here n=2,l=1,k0=2,D=4. The term u^2*z is ON the closed
face. For the leading T=-1,Z=1,

    A0(Z)=8*sqrt(Z)+Z^2-6*Z-3.

Its value and first TWO derivatives vanish at Z=1; its third
derivative is THREE. PF147's double-root normal form cannot be
used at this point.

Let Zc=(1-v)^(4/3) and Y=Z-Zc. The exact coefficients in (14) give

    A(v,Zc)=O(v^5),
    A_Z(v,Zc)=(8/27)*v^3+(10/81)*v^4+O(v^5),
    A_ZZ(v,Zc)=O(v^5),
    A_ZZZ(v,Zc)=3+O(v).

Consequently the generic local equation has initial balance

    Y^3/2=p*v^4.

There is ONE normalized local branch with ramification e=3 and
three conjugate Puiseux expressions; gcd(3,4)=1. Its derivative
A_Z has order EIGHT in the normalized parameter and (11) gives
ord(omega)=3*(6-1)-1-8=SIX. All conjugates are accounted for,
not mistaken for three separate unramified branches.

Arbitrarily higher leading multiplicity is also permitted as n
grows: choose m=2*(2*n-1), an actual MX model with r=1,k=4*n-4,
and choose the closed-face polynomial to cancel the degree-n
Taylor polynomial of sqrt(Z) at Z=1. The leading zero then has
multiplicity n+1. The proof of (10)--(13) is unchanged.

## 8. Boundaries, original finite normalization and remaining goal

For classical S2, m=2. No term with j>=1,i>=1 satisfies (1);
this displayed closed family is already the additive family of
PF147. We do not assert (13) for m=2. Above-face S2 products,
arbitrary mixed coordinates, nonconstant w coefficients and
higher w dependence remain unclassified here.

The nonlinear hypothesis and degree gap remain material. On S2,
P=w+g(u) has the rational mate 1/u; P=w+u*z has the rational mate
1/w. These are polar mates and supply NO whole-plane counterexample.
When M=1,N=u, the plane pair (w+z^2,z) is regular with unit
Jacobian. No extension outside the hypotheses is inferred.

Retain the ORIGINAL target B, source R, full finite normalization S
and trace-dual Lambda from [PF143](MARKED_QUADRATIC_VELOCITIES_AND_CYCLIC_TRACE_POLES.md),
with EVERY ramified and unramified boundary prime. The local finite
curve cover, completed field, scaled T and parameter derivative
replace none of these objects. No Tr(R) subset B, no quadratic
marked regularity, and no normalization stability is inferred.

[PF142](PARTIAL_KUMMER_PLANE_CHARTS_AND_CYCLIC_EQUIVARIANT_FACTORS.md)'s mixed
2+3 minimal case remains retired. The generic first-derivative
integrality implication remains unsupplied. This theorem excludes
a larger surviving displayed mixed-coordinate family; full HC4/JC2
and an actual counterexample remain OPEN.

The universal proof is (2)--(13), including the coefficientwise
parameter-derivative argument. Supporting exact controls retain
the full actual coframe, closed-face weights, multiple and
ramified local branches, example (14) and boundary controls.
No new external theorem import or predecessor checker replay is
used. Independent whole-proof replay, formal verification,
external review and novelty priority are unclaimed.
