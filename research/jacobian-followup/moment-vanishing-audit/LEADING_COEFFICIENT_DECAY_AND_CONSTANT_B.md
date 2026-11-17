# Leading-coefficient decay and the full constant-b classification

15 November 2026. Internal written analytic/algebraic proof; independent
specialist review outstanding. The [checker](verify_constant_b_all_degree.py)
checks finite identities, reductions and boundary models. It does **not**
certify rapid-decay comparison, action separation or the contour estimate.
Older source-bound proofs remain unchanged.

## 1. A uniform inverse-ray estimate needs no saddle expansion

Let p(x)=a_d*x^d+...+a_0 be a complex polynomial, a_d!=0, d>=2.
Fix an angle theta with cos(theta)>=kappa>0 and an endpoint A in C.
Consider inverse lifts of the ray

    w=A-exp(i*theta)*u,          u>=0.

Away from its finite branch vertices the full inverse image is a real
algebraic curve. A chain using its edges with multiplicity at most B obeys

    integral exp(Re p(x))*|dx|
       <= C_(d,kappa)*B*exp(Re A)*|a_d|^(-1/d).              (1)

The constant is independent of **every lower coefficient**, of the locations
and separations of roots and critical points, and of A except for exp(Re A).
Multiple critical points at the endpoint are permitted.

Here is an elementary proof. On the slice n<=u<=n+1 one has |p-A|<=n+1.
Factor p(x)-A=a_d*product_j(x-r_j). At least one factor satisfies

    |x-r_j| <= ((n+1)/|a_d|)^(1/d).

Thus the entire slice lies in d disks of this radius. No bound on their
centers is needed. The full inverse line is contained in

    Im(exp(-i*theta)*(p(x)-A))=0,

a real algebraic curve of degree at most d. Its length inside a disk of
radius R is at most pi*d*R by Crofton/coarea: a generic real line has at
most d intersections and its offset inside the disk ranges over length2R.
Line components and isolated singularities do not affect the generic-line
estimate. Sum over the d covering disks and multiply by the chain's bounded
multiplicity. Since exp(Re p)<=exp(Re A-kappa*n), summation gives (1) with

    C_(d,kappa)=pi*d^2*sum_(n>=0) exp(-kappa*n)*(n+1)^(1/d).

The sum converges. This proof bounds the length of the complete inverse
curve, so a singular inverse derivative near a multiple critical point
does not invalidate it.

The sharper Cartan disk-cover estimate and the projection-length method
in [Eremenko--Hayman,On the length of lemniscates](https://www.math.purdue.edu/~eremenko/dvi/erdos23.pdf),
lemmas1--3 and the proof of lemma4, motivated this step. The weaker d-disk
cover above follows directly from factorization and is sufficient; no
optimal lemniscate-length or covering constant is needed.

## 2. A no-constant-period criterion for the ordinary measure dx

Let Phi(t,x) be a rational-in-t polynomial in x of fixed generic degree d>=2,
with leading coefficient a_d(t), and let rho(t) be a nonzero algebraic
normalization. Work after finite ramification in a small generic sector at
t=infinity, avoiding a finite exceptional base set. Suppose

    |rho(t)|*|a_d(t)|^(-1/d) -> 0.                            (2)

Then no flat rapid-decay period

    rho(t)*integral_(Gamma_t) exp(Phi(t,x))*dx                 (3)

is a nonzero constant.

The action filtration used here is the inverse-phase-cut construction of
[the logarithmic criterion,section3](LOGARITHMIC_CURVE_PERIOD_CRITERION.md).
For a polynomial phase on A1 there are no finite removed points or residue
circles; inverse-sheet gluing gives critical branch stars and decay ends at
infinity. Its rank is finite. Critical locations and values are algebraic
in t; label them after Puiseux ramification. Choose a fixed generic ray
angle with cos(theta)>=kappa>0 so the corresponding critical cuts do not
hit other critical values on a sufficiently small sector. Equal critical
values use one branch star with all its incident lifts. The sheet-gluing
multiplicities and flat changes of basis are bounded there.

Critical values with a nonzero positive-power Puiseux part supply the
corresponding exponential factors. The zero-exponential block is represented
by bounded-action critical cuts. Apply (1) to each such cut: all its
unnormalized periods are O(|a_d(t)|^(-1/d)). There are finitely many cuts.
Equation (2) makes their normalized periods tend to zero. Formal/sectorial
separation of the different exponential blocks prohibits a nonzero constant
from being supplied by the other blocks. This includes merging saddles:
the whole branch star is estimated before separating its inverse sheets,
not by a nonuniform Gaussian approximation.

This action-separation/period-comparison step is analytic and remains an
explicit specialist-review obligation. The finite checker does not prove
it. The key new quantitative estimate (1) is independent of lower phase
coefficients and does not require logarithmic differentials.

The normalization condition is essential. For Phi=-x^2/(2t^2), the real
Gaussian cycle has period sqrt(2*pi)*t; rho=1/t makes it constant, while
the product in (2) is constant rather than tending to zero.
The measure is also essential: for Phi=-t^2*x^2/2, multiplying dx by t
cancels the decay. Amplitudes or coordinate Jacobians cannot be discarded.

## 3. From a formal constant seed to an analytic constant period

For the affine polynomial phase use

    H=C(t)[x]/{q_x+Phi_x*q:q in C(t)[x]},
    nabla[q]=[q_t+Phi_t*q],        v=[1].

Suppose a formal horizontal functional J exists with J(v)=1.
On the cyclic submodule N generated by v its values on positive connection
iterates are zero. Thus J|N:N->C(t) is a rational horizontal quotient.
At an ordinary base point extend its fibre functional to H, solve the dual
analytic connection, and continue on a simply connected ordinary sector.
Uniqueness preserves the restriction to N, so the period of v remains1.
A rational extension of J to all H is neither claimed nor needed.

[Hien--Roucairol,0704.1739v1](https://arxiv.org/html/0704.1739v1),
theorems2.4 and3.7 with proposition3.4, represents this analytic functional
by a flat rapid-decay cycle after removing the finite exceptional set.
For a nonconstant polynomial-in-x phase, t and Phi are algebraically
independent; clear rational t-denominators on a smooth affine total space.
This supplies the required hypotheses. There is no assertion that the
original formal Wick series converges or defines a real integration contour.

Consequently (2), with rho=1, forbids a formal horizontal functional sending
v to1. This excludes **every seed rank and every possible trivial quotient
with that seed value**, rather than deducing irreducibility from a rank bound.

## 4. Polynomial Gaussian-normalized arcs are exactly translations

Let Phi in C[t,x] and assume Phi(0,x)=-x^2/2. The coefficientwise formal
Gaussian expectation is well defined because Phi+x^2/2 is divisible by t.
Then

    E_x exp(Phi(t,x)+x^2/2)=1
        iff Phi(t,x)=-x^2/2+B(t)*x-B(t)^2/2,                 (4)

where B in t*C[t]. There is no rank-one hypothesis.

Indeed its Gaussian functional satisfies both the x integration-by-parts
relations and the actual tangent connection in section3. If d=deg_x Phi>=3,
the nonzero polynomial a_d(t) vanishes at0. It is nonconstant, hence
|a_d(t)| tends to infinity on a sector at infinity. Sections2--3 give a
contradiction to the selected seed value1.

Thus d=2. Write Phi=A(t)*x^2+B(t)*x+C(t), A(0)=-1/2, B(0)=C(0)=0.
The exact formal quadratic identity gives

    I=exp(C-B^2/(4A))/sqrt(-2A).

I=1 implies exp(2C-B^2/(2A))=-2A. The derivative of a rational function
has zero residues at every finite pole, whereas A'/A has nonzero residues
at any zero or pole of a nonconstant rational A. Differentiating the
identity forces A constant. Its initial value gives A=-1/2 and then
C=-B^2/2. Conversely the Gaussian linear generating identity proves (4)
coefficientwise. The high-degree exclusion is analytic; this final
quadratic step is algebraic.

The proof uses polynomial dependence on t. For example the rational arc
Phi=-x^2/2+t*x^3/(1+t)^2 has the same Gaussian initial phase, but its
cubic leading coefficient tends to zero at infinity, so (2) does not
apply. This is a boundary of this argument, not a counterexample to (4).
No classification of all rational-time arcs is asserted.
The earlier [rank-one rigidity theorem](ACTUAL_DIRECTION_RANK_ONE_RIGIDITY.md)
retains its broader rational-time direction criterion and remains unchanged.

## 5. The complete one-spectator constant-b architecture

Let E(ZW)=1, E(T^2)=1, and consider

    P=W*(a+c*T)+G(Z,T),       a,c in C, c!=0, G polynomial.   (5)

The circular fixed point is Z=t*(a+c*T), with Jacobian factor exactly1.
Thus its formal pure series equals

    I(t)=E_T exp(t*G(t*(a+c*T),T)).

Its moving phase Phi=-T^2/2+t*G(t*(a+c*T),T) is polynomial in t,T and
has the Gaussian initial phase. Apply (4): t*G(t*(a+c*T),T) must be linear
in T, with constant term minus one half the square of its linear coefficient.

Write G=sum_i Z^i*g_i(T). Distinct i give distinct t^(i+1) coefficients.
The degree bound in T forces all g_i=0 for i>=2, g_1=k constant, and
g_0=p*T+q. Hence

    t*G(t*(a+c*T),T)=t*(p*T+q)+k*t^2*(a+c*T).

Its linear coefficient is B=p*t+k*c*t^2. Comparing its constant coefficient
with -B^2/2, first t^4 gives k=0 because c!=0, then t^2 gives p=0, then
t gives q=0. Conversely G=0 has strictly negative circular weight. Therefore

    all positive E(P^m)=0      iff G=0.                      (6)

This removes the actual seed-rank restriction from the earlier constant-b
corollary and closes the formerly open higher-rank trivial-quotient branch
for (5). For nonzero G, the formal pure series is not1; a positive pure
moment must be nonzero. No larger fixed-scale Krylov calculation is needed.

The surviving P=W*(a+c*T) are GMC-safe with sharp coefficient1:

    E(Q*P^m)=0 whenever m>deg Q.

Circular charge gives the upper bound. If a!=0, Q=Z,m=1 has value a.
If a=0, Q=Z^2,m=2 has value2*c^2. Thus every smaller real coefficient fails
for each nonzero-c normalized member, including both a branches.

This theorem has one real Gaussian spectator. It does not classify the
several-spectator constant vector-b case, the b=0 case, or repeated-root
nonconstant b. It does not assert the same classification over nonreduced
coefficient algebras: the square arguments use a field.

This scheme boundary already occurs in degree two. On the line G=p*T,
the pure series is exp(p^2*t^2/2), so its full ordinary positive-pure-moment
ideal is (p^2), not (p). Over Q[p]/(p^2), a nonzero G=p*T is pure-null.
The classification (6) describes field points and their radical locus;
it does not replace the retained ordinary pure-ideal algorithms by the
reduced equations G=0.

## 6. Explicit classification before the shear normalization

Now allow arbitrary polynomial a(Z), still with constant c!=0:

    P=W*(a(Z)+c*T)+G(Z,T).

Set

    a_0=a(0),        u(Z)=(a(Z)-a_0)/(c*Z).

The quotient is polynomial. The retained exact Gaussian shear theorem gives

    P composed with Phi_u^(-1)
      =W*(a_0+c*T)
        +(u*T-Z*u^2/2)*(a_0+c*T)+G(Z,T-Z*u).

By (6), pure nullity is equivalent to vanishing of the whole displayed
remainder. Equivalently the original polynomial has the explicit factorization

    P=(W-u*T-Z*u^2/2)*(a(Z)+c*T),                            (7)

or

    G=-(u*T+Z*u^2/2)*(a(Z)+c*T).

These are exactly the Gaussian shears Phi_u of W*(a_0+c*T).
Thus the arbitrary-degree one-spectator nonzero-constant-b pure-null locus
has a complete polynomial normal form, is safe, and has the transported
bound m>delta*deg Q, with delta the actual shear/inverse degree.
This bound is not asserted optimal before normalization. Dimension stays3;
nonlinear shears do not preserve the ordinary constant-coefficient Laplacian.

At any fixed coefficient support, the explicit equations for the remainder,
together with c!=0, describe the pure-null locus algebraically on that chart.
They add a concrete safe structural family, not a proof of the global
coefficient safe-locus stabilization conjecture. The usual finite-support
pure ideal algorithms remain available; no new least pure prefix is computed.

## 7. Updated frontier and assurance

The one-spectator constant-b higher-rank rational-trivial-quotient problem
is now closed internally by (1)--(7), with analytic comparison/action
separation still requiring independent specialist review. Its old rank-one
subcase remains a valid special result; its broader rational-time theorem
has not been superseded.

Remaining structural frontiers are nonconstant repeated-root b, the b=0
branch beyond retained cases, several-spectator normal forms, higher neutral
radial degree, genuinely mixed angular signs, and actual generic higher-rank
connections in architectures not excluded here. Fixed-connection boundary
subspaces, exceptional coefficient restrictions, global complex-field
uniformity, valid operator compression, dense SIC33, least Delta power and
prior analytic assurance remain separate. Numeric cubic refinement stays
secondary.

[Source and assurance audit](constant_b_all_degree_source_audit.json).

