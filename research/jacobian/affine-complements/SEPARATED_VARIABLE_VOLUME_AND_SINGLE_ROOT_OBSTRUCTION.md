# Exact volume forces one base root in separated-variable plane extensions

PF2RULED137. 16 November 2026. All polynomial degrees in characteristic
zero, geometrically. This excludes a complete test family of possible
affine completions, not all affine-line fibrations or plane Keller maps.
The full Wright/HC4/JC2 and finite trace-dual gate remain OPEN.

## 1. Statement and relation to the actual problem

Let k be algebraically closed of characteristic zero, let a(x) be a
nonzero polynomial, and let p(z) be square-free of degree at least2.
Set

    S(a,p) = {a(x)*y=p(z)} subset A3,
    Omega = dx wedge dz/a(x).                                  (1)

Then S is smooth affine, contains a whole open A2, has only constant
units, and Omega is a regular nowhere-zero two-form. We prove

    Omega has a GLOBAL REGULAR primitive
      iff a is constant, or a=c*(x-b)^m with c!=0 and m>=2.       (2)

Consequently NO etale morphism S(a,p)->A2 exists if a has at least TWO
distinct roots, or if a is linear. There is no bound on either degree,
the number of roots, their multiplicities, or the proposed target degree.
An actual canonical completion U_N has an etale target and an exact unit
volume, so these surfaces cannot be such a U_N, even if its chosen plane
chart and all its other necessary conditions were available.

The one-root m>=2 family is exactly the retained higher-power Danielewski
family after an affine change and rescaling. It passes exactness; a
compatible etale target is STILL OPEN. This is not a new construction of
that family, its atlas or its higher-power primitive. We credit
UADANCHART1 and UADANVOL1 and extend their one-base-root calculation to
ARBITRARY a. HC4MYMP0/MQ0 already settle actual noncyclic deck cases and
source-image point gaps; neither result assumes all remaining surfaces
have the separated equation (1).

## 2. Whole plane charts and a rational-translation exactness test

Write the roots of p as alpha_i and p_i(z)=p(z)/(z-alpha_i). For each i
remove the curves {x=b,z=alpha_j} for every root b of a and every j!=i.
The resulting open U_i is A2_(x,w_i), by the explicit inverse charts

    z=alpha_i+a(x)*w_i,
    y=w_i*p_i(alpha_i+a(x)*w_i).                              (3)

Away from a=0 the inverse is w_i=(z-alpha_i)/a. Near the retained
curves z=alpha_i it is w_i=y/p_i(z), with p_i(alpha_i)=p'(alpha_i)!=0.
These opens cover U_i. All U_i together cover S, since at a=0 precisely
one root label occurs at each point. They establish smoothness and
regularity of (1), whose pullback is dx wedge dw_i. Equivalently the
hypersurface Jacobian is nonsingular there because p'!=0.

The surface is integral: its defining polynomial, linear in y, is
primitive over k[x,z], since a(x) and p(z) are relatively prime; it is
irreducible over k(x,z). Every U_i is dense. A unit on S restricts to
a unit of k[x,w_i], hence is constant on that dense chart and on S.

For i!=j the EXACT intersection is

    U_i intersection U_j = {a!=0} times A1_w,
    w_j=w_i+(alpha_i-alpha_j)/a(x).                            (4)

Suppose dtheta=Omega for a globally regular one-form theta. On U_i,
theta_i=-w_i*dx is a regular primitive. The polynomial Poincare lemma
on A2 in characteristic zero gives a polynomial H_i(x,w_i) with
theta-theta_i=dH_i. One may prove this lemma by integrating the first
coefficient in x and then the remaining coefficient in w; closedness
removes the other variable. Thus no topological period construction or
formal-to-algebraic extension is being assumed.

Equality of theta on (4) now gives

    d[H_j(x,w_i+(alpha_i-alpha_j)/a)-H_i(x,w_i)]
       =(alpha_i-alpha_j)*dx/a(x).                            (5)

The bracket is rational in x and polynomial in w_i. Its w_i derivative
is zero; characteristic zero makes it an element of k(x). Since
alpha_i!=alpha_j, dividing by their difference produces R in k(x) with

    R'(x)=1/a(x).                                             (6)

More generally, for any two whole plane charts with the SAME base x
and transition w_j=w_i+h_ij(x), an exact common dx wedge dw volume
requires the rational differential h_ij(x)*dx to have a rational
primitive. In particular every finite residue must vanish. This is a
NECESSARY test; we assert no general sufficiency for arbitrary atlases.
Different x-dependent root labels can pass it, so (6) is specific to
the separated-variable constant-root-label equation (1).

## 3. The all-degree rational primitive lemma

For nonconstant a, a rational R satisfying (6) exists iff

    a=c*(x-b)^m with m>=2.                                   (7)

Here is a degree proof. Factor a=c*product_(i=1)^r(x-b_i)^n_i,
with N=sum n_i. A derivative of a rational function has zero residue
at every finite pole. A simple root of a would give the nonzero
residue 1/a'(b_i), impossible. Hence all n_i>=2 and N>=2.

R has precisely the finite poles b_i of orders n_i-1: differentiation
raises each pole order by one and cannot create or cancel another
pole. There is no nonconstant polynomial part, since R'=1/a tends to
zero at infinity. Thus R has finite limit R_infinity and, as a map
P1->P1, degree

    deg R=sum_i(n_i-1)=N-r.                                 (8)

Expanding R' at infinity and integrating the leading term gives

    R-R_infinity = [1/(c*(1-N))]*x^(1-N)+lower powers of x.    (9)

The local degree at infinity is therefore N-1. It cannot exceed the
degree N-r of the rational map: equivalently, a zero of the rational
function R-R_infinity has order at most its total pole degree. Hence
N-1<=N-r, or r<=1. This contradicts r>=2 in EVERY degree, even when
all roots are repeated and some individual residues happen to vanish.

Conversely for (7), R=(x-b)^(1-m)/(c*(1-m)) satisfies (6).
Partial fractions also gives the equivalent residue formulation:
1/a has a rational primitive iff all its finite residues vanish.
The simple-pole coefficients are precisely the terms whose integrals
are not rational; every higher-order term has a rational primitive.

Together with Section2 this proves the NECESSITY in (2). It also
shows that at least one residue of dx/a is nonzero whenever a has two
distinct roots. This does not identify a particular residue without
calculation, nor claim that every residue is nonzero.

## 4. Sufficiency, credited to the retained monomial primitive

If a is constant, (x,z) are global plane coordinates and
(x/a)*dz is a regular primitive of Omega.

If a=c*(x-b)^m, m>=2, put u=x-b and v=c*y. Then u^m*v=p(z) and
Omega=c^(-1)*du wedge dz/u^m. Reuse UADANVOL1 Section3's explicit
regular primitive. Choose E(z) inverse to p'(z) modulo p(z), set

    D=(1-E*p')/p,
    H'=(D-E')/(m-1), J=H-E.

Here D and H are polynomials, using square-freeness and characteristic
zero. The GLOBAL REGULAR form

    theta=c^(-1)*(H(z)*u*dv+J(z)*v*du)                       (10)

satisfies dtheta=Omega because

    (H-J)*p'+(m*H'-J')*p=E*p'+D*p=1.

This is the prior monomial formula under an isomorphism, not a new
primitive construction or predecessor-checker replay. It proves the
remaining implication in (2), but no polynomial Keller pair follows.

## 5. Target obstruction, controls and exact scope

If e=(P,Q):S->A2 were etale, dP wedge dQ would be a nowhere-zero
regular form and hence f*Omega for a unit f of S. Section2 gives
f in k*. It is exact, with primitive P*dQ. Thus Omega would be exact,
contradicting (2) in the excluded cases. This needs neither a finite
map e, a degree bound, a deck group nor a minimal counterexample.

The hypothesis deg p>=2 is essential. For p=z, S is A2_(x,y) for
EVERY a, and Omega=dx wedge dy, even for multiple distinct roots of a.
There are no two different root charts from which to obtain (5).
Square-freeness is also retained; singular repeated-root models are
outside this statement.

The genuinely exact mixed-fibre (2,3) surface in PF2RULED130 Section5
has an x-dependent right side, and its ruling is not (1). It is not
excluded here. Neither arbitrary affine-line surfaces nor all the
actual reduced-fibre universal covers in MQ0 are classified by (1).
The surviving one-root higher-power surfaces still require compatible
plane Keller pairs on their WHOLE atlas (UADANCHART1), the original
full finite order and regular trace-dual multiplication. The trivial
deck case is also open. Original ramified primes and the residual
discriminant are never removed by this obstruction.

The new checker supplies only small exact chart, residue, degree and
scope controls; Sections2--4 establish the unbounded conclusions.
No bounded scan, new surface classification, external theorem import,
independent replay, formal verification, external review or novelty
claim is made. The literature search did not establish priority for
this elementary extension. The classical atlas and monomial primitive
are explicitly credited to the retained source notes, not claimed new.
