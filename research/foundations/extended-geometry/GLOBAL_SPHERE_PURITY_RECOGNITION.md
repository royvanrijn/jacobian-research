# Exact recognition of global sphere purity

Let K have characteristic zero and let

    Q = Spec K[x,y,t]/(t^2+2xy-1).

For a polynomial class p define

    Lambda(p) = CT_x (1/2) integral[-1,1]
                         p(x,(1-t^2)/(2x),t) dt,
    G_p(z) = sum[m>=0] Lambda(p^m) z^m.

The integral in the definition is the rational rule on polynomials. Purity
means G_p=1. The following statements concern arbitrary polynomial classes
on Q, without any restriction on the number of height terms.

**Theorem GlobalSpherePurityRecognition.** There is a terminating exact
algorithm which, on input p over Q, decides purity. It returns either a
nonzero rational moment or a differential equation and a finite initial
vanishing certificate proving that every positive moment is zero. The
algorithm also works over an explicitly presented finitely generated field
extension of Q.

**Theorem ComputableUniformSpherePrefix.** There is a computable function
B(d) such that for every characteristic-zero field K and every class p
represented by a polynomial of total degree at most d,

    p is pure  iff  Lambda(p^m)=0 for 1<=m<=B(d).

Here computable means that a specified terminating algorithm outputs B(d).
No numerical bound on its size or practical complexity is asserted. The
argument proves equality of radicals of the infinite moment ideal and a
finite prefix ideal. It does not prove equality of the ideals themselves.

The D-finiteness and regular telescoping results used below are classical.
The contribution of this argument is their explicit connection to this
functional, an exact identity test at a possibly singular origin, and the
coefficient-space termination argument. No novelty of the general period
or creative-telescoping methods is asserted.

## 1. A closed rational period

First take K to be a finitely generated extension of Q and embed it into C.
This embedding is used only in the proof, not in the exact computation.
Write

    q=t^2+2xy-1, Omega=dx wedge dy wedge dt,
    Gamma={y=conjugate(x), t real, 2|x|^2+t^2=1}.

The affine complex quadric q=0 is smooth and Gamma is a compact oriented
real two-sphere in it. Its residue form is

    omega=Res(Omega/q).

For the radial vector field E=x d/dx+y d/dy+t d/dt, one has
omega=(i_E Omega)/2 on q=0, because dq(E)=2 there. On the chart x!=0 this
is, with the corresponding orientation, -dx wedge dt/(2x). Parametrizing
x=sqrt((1-t^2)/2) exp(i theta) proves

    integral_Gamma f omega / integral_Gamma omega = Lambda(f),
    integral_Gamma omega != 0.                                  (G1)

The poles at x=0 in the chart formula are coordinate artifacts; omega is
regular on the smooth quadric. The two exceptional poles of the angular
parametrization are measure-zero points of the real sphere.

There is an explicit closed tube, so no endpoint terms in t are discarded.
Fix 0<epsilon<1. For rho=epsilon exp(i theta) and P in Gamma put

    Tube(Gamma): (rho,P) |-> sqrt(1+rho) P,

using the square root analytic near rho=0 with value 1. This is a compact
closed real three-cycle in C^3 on which q=rho. For all sufficiently small
complex z, neither q nor 1-zp vanishes on this cycle. Locally, using rho
as the transverse coordinate, the pullback of Omega is

    d rho wedge (1+rho)^(1/2) omega.

The ordinary one-variable Cauchy formula in rho therefore gives

    I_p(z) := integral_Tube(Gamma) Omega/[q(1-zp)]
            = 2*pi*i integral_Gamma omega/(1-zp),
    I_p(z)/I_p(0)=G_p(z),  I_p(0)!=0.                           (G2)

Compactness justifies the convergent geometric expansion for small z and
its termwise integration. Equation (G2) is also the equality of the
resulting formal series. It proves that the functional is independent of
the polynomial representative modulo q.

## 2. The precise imported creative-telescoping theorem

Bostan, Lairez and Salvy, *Creative Telescoping for Rational Functions Using
the Griffiths-Dwork Method*, arXiv:1301.4313v2 (21 April 2013), Section 7,
Theorem 12 (PDF page 7), treats an arbitrary affine rational function
F in k(z)(x_1,...,x_n), over a characteristic-zero field. It gives an
algorithm computing a nonzero telescoper with a **regular certificate**.
No smoothness hypothesis on its affine denominator is required in this
statement; Section 7 supplies the deformation for singular denominators.

Regular means that, over k(z), every irreducible factor in a certificate
denominator already divides the denominator of F. Proposition 4, Section
4.1 (PDF page 4), then proves that this telescoper annihilates the integral
over a closed n-cycle avoiding the poles. The warning and example just
before Proposition 4 explain why an arbitrary certificate is insufficient.

Apply this theorem to

    F(z,x,y,t)=1/[q(x,y,t)*(1-zp(x,y,t))].                      (G3)

The cycle in Section 1 satisfies exactly the closed-cycle and pole-avoidance
hypotheses. A rational function of z occurring in certificate coefficients
may have a pole at z=0. Apply Proposition 4 first at generic sufficiently
small nonzero z and clear all coefficient denominators. Analyticity of
G_p then gives the polynomial differential identity at z=0 as well.
Consequently a computable nonzero operator

    L=sum[j=0..r] p_j(z) D_z^j,  p_j in K[z],

annihilates G_p. The identity holds over K: after the chosen injection into
C it holds coefficientwise, and injectivity returns it to K.

An alternative elementary implementation makes computability explicit
without relying on a black-box implementation of their algorithm. Put
D=q(1-zp). Search, fairly over all orders and degree bounds, for

    L(F)=d_x(U_x/D^N)+d_y(U_y/D^N)+d_t(U_t/D^N),               (G4)

where the coefficients of L and the U_i are polynomials over K, the
coefficients of L depend only on z, and N>=max(1,r). Clear D^(N+1) and
solve the resulting finite homogeneous linear system over K. Require the
projection to the coefficients of L to be nonzero; a nonzero certificate
with zero L does not qualify. Search all finite tuples of bounds, for
example in order of their sum. Such a search terminates: the theorem gives
a regular certificate, each of whose denominators divides some power of
D over K(z); clearing the remaining denominators in z puts it in (G4).
All polynomial identities in (G4) can be replayed by exact multiplication
and differentiation. This implementation is a termination proof, not an
efficiency claim.

## 3. Identity recognition at a singular origin

Here is the needed formal lemma, valid over every characteristic-zero field.
Let

    M=sum[j=0..r] a_j(z) D_z^j != 0, a_j in K[z].

For a_j!=0 write v_j=ord_0(a_j) and l_j=[z^v_j]a_j, and put

    nu=min_j(v_j-j),
    Ind_M(s)=sum[v_j-j=nu] l_j (s)_j,                         (G5)

where (s)_j=s(s-1)...(s-j+1). This is not the zero polynomial: the term
with largest participating j has strictly highest degree. If a nonzero
formal H in K[[z]] satisfies M(H)=0 and N=ord_0(H), the coefficient at
z^(N+nu) is its initial coefficient times Ind_M(N). Thus Ind_M(N)=0.
This calculation is still valid if N+nu<0, because all relevant falling
factorials are then zero. In particular a nonzero formal solution must
start at one of the finitely many nonnegative integer roots of Ind_M.

To apply this to H=G_p-1, let h=-L(1)=-p_0(z). If h=0, take M=L.
Otherwise take

    M=(h D_z-h') composed with L.                            (G6)

This operator is nonzero and M(H)=0. The product is nonzero, for example
by its nonzero leading differential coefficient. Let

    B_L=max({0} union {N>=0 integer: Ind_M(N)=0}).

Since H(0)=0, checking the moments 1,...,B_L proves H=0. Conversely a
nonzero one of these moments disproves purity. This proves a terminating
exact decision procedure once the integer roots are computed. Ordinary
point uniqueness for an ODE is not used.

For K=Q, integer roots are computable by rational polynomial factorization
or a root bound and exact tests. For a finitely generated field
K=Frac(Q[c_1,...,c_s]/P), with P prime presented by a Groebner basis, clear
all parameter denominators in Ind_M(s) and reduce its coefficients modulo
P. Expand the result in standard parameter monomials. It is a finite sum
sum_b b * I_b(s), with I_b in Q[s] and at least one I_b nonzero. Any
integer root in K must be an integer root of this one nonzero rational
polynomial. Enumerate those finitely many candidates and test them in K.
Thus the generic-point calls used below are effective too.

The exact sphere moments themselves require no period-zero oracle. For
monomials,

    Lambda(x^i y^j t^k)=0 unless i=j and k=2r,
    Lambda(x^i y^i t^(2r))
       =2^i*i!*(2r)!*(i+r)!/[r!*(2i+2r+1)!].                 (G7)

Expand a finite polynomial power and apply (G7). All initial values used
by the decision procedure are exact rational expressions in the input
coefficients.

## 4. A computable degree-only prefix

Fix a nonnegative integer d and let c be one variable for every monomial x^i y^j t^k with
i+j+k<=d. Redundancy modulo q is harmless. Let p_c be this universal
polynomial and let

    mu_m(c)=Lambda(p_c^m) in R=Q[c].

These homogeneous coefficient polynomials are computable by (G7).
Maintain an ideal I generated by finitely many mu_m, initially I=(0).
It is always proper: every positive moment vanishes at the zero coefficient
vector c=0. Thus the proper-input hypothesis of minAssGTZ is satisfied.
The following is an algorithm.

1. Compute all minimal prime ideals P of I in R.
2. For each such P form the effective field K_P=Frac(R/P) and apply the
   decision procedure of Sections 1--3 to the class of p_c over K_P.
3. If the result is nonpure, choose a returned nonzero moment mu_r mod P,
   replace I by I+(mu_r), and restart from step 1.
4. If every component is pure, stop, and output B equal to the largest
   moment index ever added, or B=0 if none was needed.

Minimal-prime decomposition and exact computations in these finitely
generated fields are standard effective Groebner and factorization
operations over Q. This is the only commutative-algebra algorithmic
primitive beyond field arithmetic required here; no complexity estimate
for it is assumed. Gianni--Trager--Zacharias (1988), *Groebner bases and
primary decomposition of polynomial ideals*, JSC 6, 149--167, is the
classical algorithmic reference. The publisher abstract explicitly states
primary decomposition, radicals and primality over algorithmic base rings;
the period argument does not depend on a numerical bound from that paper.
For the exact algorithmic primitive one can use the primary maintained
Singular 4.4.1 manual, Section D.4.27.7, `minAssGTZ`, at
https://www.singular.uni-kl.de/Manual/latest/sing_1534.htm . Its specified
return is precisely the minimal associated prime ideals of a proper input
ideal; it is designed for characteristic zero. Use a global polynomial
ordering over Q. This proof needs Q-primes, not absolute prime decomposition.
The saved primary page is [the retained official manual](sources/singular_sing_1534.html),
SHA256 48e747ffacca53148107971d39415f9217c8d9b95f690728b08faf30a35793c9.

Every unsuccessful pass strictly increases sqrt(I). Indeed mu_r is not in
one current minimal prime P, whereas sqrt(I) is contained in P. Therefore
mu_r is not in sqrt(I). Polynomial rings are Noetherian, so infinitely many
such strict increases are impossible. Each pass makes only finitely many
terminating field-level calls. The entire algorithm consequently terminates.

At termination every mu_m belongs to every minimal prime of I, since the
generic phase was pure on each component. Hence

    sqrt((mu_1,mu_2,...))=sqrt(I).                            (G8)

Because I is generated by moments of index at most B, this is also the
radical of the full prefix ideal (mu_1,...,mu_B). It follows over ANY
characteristic-zero coefficient field: a coefficient evaluation killing
that prefix has a prime kernel containing I, hence containing every mu_m
by (G8). The converse is immediate. This proves the uniform theorem.

The reverse implication also holds at every fixed B. Suppose the field-valued prefix criterion holds over every characteristic-zero field, but some mu_m does not belong to the radical of the prefix ideal. Choose a prime ideal containing that prefix but excluding mu_m. Evaluation in the fraction field of its quotient makes the first B moments zero and mu_m nonzero, a contradiction. Thus the universal field-valued criterion and reduced moment-ideal equality are equivalent. This proves no ordinary-ideal equality or nilpotence exponent.

This is an effective, potentially very expensive bound, rather than a
closed formula. A bound on telescoper order and polynomial coefficient
degree alone would not justify it: (z D_z-N) z^N=0 has fixed order and
coefficient degree and arbitrarily late first nonzero coefficient.

A small illustration is the linear family p=c+a x+b y+e t. The first two
moments give c=0 and 2ab+e^2=0. For the centered family direct expansion of
(G7) gives Lambda(p^(2r))=(2ab+e^2)^r/(2r+1) and odd moments zero. Thus
the generic component is pure after two strict radical enlargements, and
B(1)=2 is valid. For d=0, the constant-input case has the valid bound B(0)=1. The same
constant case applies when the original balanced degree d is zero. This
illustration is not a computation of B(d) for higher d.

## 5. An exact three-height fixture

Consider

    p=t*x+2*t*y*(1+t^2)=t/x*(1+x^2-t^4).                    (G9)

It has three height coefficients C=1+x^2, D=0, E=-1. Constant-term
extraction of the rational resolvent gives

    G_p(z)=J(4z^2),
    J(s)=integral[0,1] (1-s*(t^2-t^6))^(-1/2) dt.           (G10)

Indeed CT_x 1/[1-z(a x+b/x)]=(1-4z^2 ab)^(-1/2) as a formal series,
with a=t and b=t(1-t^4). Put A=t^2-t^6 and Y=1-sA. A literal fourth-order
operator annihilating J has coefficients

    c0=-12*(4s^2-84s+7),
    c1=-6*(232s^3-3612s^2-14s-147),
    c2=-12s*(232s^3-3052s^2-94s+2121),
    c3=-8s^2*(152s^3-1764s^2-378s+5103),
    c4=-16s^3*(2s-21)*(4s^2-27).

The exact sparse input global_sphere_three_height_certificates.json gives
a polynomial R with R(s,0)=R(s,1)=0 and the identity

    sum[k=0..4] ck(s)*(1/2)_k*A^k*Y^(4-k)
        = (d_t R)*Y-(7/2)*R*(d_t Y).                        (G11)

Here (1/2)_k is the rising factorial. Dividing by Y^(9/2) gives
L(Y^(-1/2))=d_t(R/Y^(7/2)). The boundary terms vanish literally, so this
fixture has a fully explicit ODE proof, independently of a generic
creative-telescoping computation.

For H=J-1 the homogenized operator in (G6) has lowest shift -2 and

    Ind_M(N)=-10584*N*(N-1)*(2N-1)*(6N-7)*(6N+1).

Its only nonnegative integer roots are 0 and 1, so a single coefficient
of J-1 decides this fixture. It is J_1=2/21, hence Lambda(p^2)=8/21!=0.
The portable checker verifies (G11), the boundary values, the homogenized
indicial polynomial, direct moments and ODE coefficients, and rejects a
corrupted certificate. The fixture is an actual computation beyond the
two-height ansatz, not an implementation of a general three-variable
creative-telescoping engine.

There is also an explicit pure control using the direct rational tube
integrand: p=x*(1+t^2+t^4). For V=x*d_x-y*d_y one has V(q)=0, V(p)=p,
and div(V)=0. Hence for F=1/[q(1-zp)],

    z*d_z F = d_x(xF)-d_y(yF).                               (G11a)

This is a regular certificate for L=zD_z. Its origin indicial is N, so
H=G_p-1, H(0)=0 must be zero without any additional moments. The checker
verifies the cleared polynomial identity behind (G11a) as well as finite
raw moments. Thus both purity outcomes have an exact certificate control.

## 6. Every bihomogeneous balanced pair number

Let n>=1, let F in K[xi_1,...,xi_n,z_1,...,z_n] be bihomogeneous of degree
(d,d), and define E_n(xi^alpha z^beta)=d_z^alpha z^beta. Every E_n(F^m) is
then a scalar. On the complex unit sphere Sigma={xi=conjugate(z),
sum |z_i|^2=1}, its normalized average satisfies

    average_Sigma F^m
      =(n-1)!*E_n(F^m)/(dm+n-1)!.                           (G12)

One direct proof uses phase invariance for alpha!=beta, and independent
complex Gaussian variables for alpha=beta: each coordinate has moment
alpha_i!, while the radial moment is (dm+n-1)!/(n-1)!. Radial and angular
variables are independent. This proves (G12) term by term.

Use q=sum xi_i z_i-1 in C^(2n). The real sphere Sigma is a compact closed
(2n-1)-cycle on this smooth quadric. Its residue form has a nonzero
constant multiple of Euclidean spherical volume; this follows by the
linear transformation z_i=X_i+iY_i, xi_i=X_i-iY_i and the radial contraction
formula used in Section 1. For completeness the scaling calculation in m real coordinates is
entirely elementary. After the indicated complex linear coordinate change,
q=sum w_i^2-1. For s in the unit real sphere set w=r(eta)*s with
r(eta)=sqrt(1+eta). The ambient volume form pulls back to
r^(m-1) dr wedge omega_S, where omega_S is the usual spherical volume.
Since dr=deta/(2r), integration around |eta|=epsilon gives

    integral_Tube Omega/[q(1-zp)]
       =pi*i integral_S omega_S/(1-zp(s)),

up to the single nonzero constant from the complex linear coordinate
change. This proves nonzero unweighted normalization directly, without
appealing to an unmentioned residue theorem. In the balanced n-pair case
m=2n and F(r*s)=(1+eta)^d F(s). The same scaling tube is a closed
2n-cycle. Apply the imported affine rational theorem in 2n integration
variables to 1/[q(1-zF)]. All arguments of Sections 2--4 now apply, and
(G12) transfers moment vanishing to E_n(F^m).

**Corollary BalancedImagePurityRecognition.** For every fixed n there is a
terminating exact purity algorithm for arbitrary rational-coefficient
bihomogeneous F of degree (d,d). There is a computable function B(n,d)
such that checking E_n(F^m)=0 for 1<=m<=B(n,d) proves vanishing for every
positive m, over every characteristic-zero field. The coefficient-space
algorithm can use exactly the (d,d) monomials of F, and its moments can be
computed directly by E_n, with the nonzero factorial normalization in
(G12).

Bihomogeneity is essential to this specific normalization: if terms of
different balanced degrees are added, their factorial factors differ.
No assertion about arbitrary nonhomogeneous Image inputs follows merely
by applying (G12).

## 7. Deciding eventual vanishing for a fixed marker

**Theorem FixedMarkerEventualVanishingRecognition.** For any rational
sphere polynomials p and g, there is a terminating exact algorithm deciding
whether Lambda(g*p^m) is zero for all sufficiently large m. For any fixed
bihomogeneous balanced F on n pairs and ANY fixed polynomial marker Q,
there is a terminating exact algorithm deciding whether E_n(Q*F^m) is
zero for all sufficiently large m. No purity assumption on p or F is
required. This theorem does not decide safety quantified over all markers.

For sphere inputs use the numerator g in (G3). Regular telescoping applies
to g/[q(1-zp)], with the same compact closed cycle. Normalize this weighted
period by I_1(0), the UNWEIGHTED volume, which is nonzero. In general the
weighted mean I_g(0) can be zero and cannot be used as a denominator. Its
series A_g(z)=sum Lambda(g*p^m) z^m is D-finite, and eventual vanishing is
exactly the assertion that A_g is a polynomial.

Polynomiality of a D-finite series with computable coefficients is also
decidable. For a verified nonzero annihilator L=sum p_j D^j put

    sigma=max_j(deg p_j-j),
    Ind_infinity(N)=sum[deg p_j-j=sigma] lc(p_j)*(N)_j.        (G13)

The polynomial Ind_infinity is nonzero, by its largest participating j.
Every NONZERO polynomial solution of L has degree at a nonnegative integer
root: compare the highest possible exponent N+sigma. This remains valid
when some derivatives kill low-degree polynomial terms, using falling
factorials exactly as in Section 3. If there are no such roots, the only
polynomial solution is zero. Otherwise let Nmax be the largest such root.
Truncate A_g through degree Nmax, obtaining P. Then A_g is polynomial if
and only if A_g=P. If the root set was empty put P=0 instead.

For an exact test of this equality set h=-L(P). Use L itself if h=0, and
otherwise (hD-h')L. This nonzero operator annihilates A_g-P. The origin
indicial test of Section 3 proves or disproves zero by finitely many exact
coefficients. Include the constant coefficient in this test; it is not
automatically zero when the infinity-root set is empty. This gives an
exact polynomiality certificate or a proof that infinitely many coefficients
are nonzero. Integer roots over generic finitely generated fields are
effective by the same rational-coordinate procedure as in Section 3.

For a bihomogeneous balanced original marker of degree (e,e), formula
(G12) already transfers this decision to E_n(G*F^m). To cover a general
original marker and all possible cancellations, decompose it into its
bihomogeneous pieces Q_(r,s). Pieces with r>s are killed by E_n identically.
For each output degree ell>=0 put

    Q_ell=sum[s-r=ell] Q_(r,s).

The polynomial E_n(Q_ell*F^m) is homogeneous of degree ell in z. For every
multi-index gamma with |gamma|=ell, the operator identity

    d_z^gamma E_n(H)=E_n(xi^gamma H)

shows that its gamma-th coefficient, times gamma!, is
E_n(xi^gamma Q_ell F^m). There are only finitely many ell and gamma to test.
For a fixed pair (ell,gamma), the markers

    G_s=xi^gamma Q_(s-ell,s)

are balanced of degree (s,s). Let e be the least occurring s, and let
A_s(z) be their normalized sphere-period series. By (G12), the scalar
sequence being tested, divided by the nonzero factorial
(dm+e+n-1)!, is the coefficient sequence of the D-finite series

    B_(ell,gamma)(z)
      =1/(n-1)! sum_s P_s(theta) A_s(z),
    theta=zD_z,
    P_s(v)=product[j=e+n .. s+n-1] (dv+j),                   (G14)

with an empty product equal to 1. Indeed P_s(m) is exactly the ratio of
the two factorials in (G12). Taking the least s makes P_s a polynomial,
so no division by a parameter-dependent sequence is needed.

For a direct algorithm, apply the finite sum of Euler operators in (G14)
to the corresponding rational tube integrands before creative telescoping.
Differentiating in z preserves rationality, and integration commutes with
these derivatives on the fixed compact tube for small z. Consequently
B_(ell,gamma) itself is one rational period to which Theorem 12 applies;
no separate closure theorem for D-finite series is needed. Its initial
coefficients are computed by E_n or by (G7). Polynomiality of every
B_(ell,gamma) is therefore decidable by (G13) and the origin test. All coefficient cancellations between different
bidegrees are preserved in (G14). Their eventual vanishing is equivalent
to eventual vanishing of E_n(QF^m), because there are only finitely many
output coefficients and every factorial normalization is nonzero.

For example p=y and g=x has A_g=z/3, a polynomial with zero constant term.
The choice p=t, g=t gives A_g=z/3+z^3/5+..., which is not polynomial. These
are literal weighted-moment fixtures in the companion formal-indicial
packet. For polynomiality, in contrast to identity with 1, one must test
both the indicial behavior at infinity and the exact residual at the origin.

**Corollary ImageFailureFiniteCertificates.** For a specified rational
bihomogeneous balanced F and a specified polynomial Q, the conjunction

    E_n(F^m)=0 for every m>=1, and
    E_n(QF^m)!=0 for infinitely many m

is equivalent to existence of finite exact certificates of the preceding
forms: a regular telescoper and origin vanishing certificate for the pure
series, together with the finite marker aggregation and an annihilator,
infinity degree bound and a nonzero origin residual certificate for at
least one aggregated output series. Soundness follows from the origin and
infinity lemmas. Existence follows from the terminating recognition
algorithms. This equivalence is about the specified pair (F,Q); it does not
decide whether some marker Q exists or whether a parameter variety has a
rational point.

## 8. Scope of the evidence and next quantitative problem

The compact-domain polynomial moment-vanishing problem and the implication
from pure vanishing to eventual mixed vanishing already appear in
Francoise--Pakovich--Yomdin--Zhao, *Moment vanishing problem and positivity:
Some examples*, Bull. Sci. math. 135 (2011), 10--32, Sections 2.1--2.2 and 3.
The actual author PDF https://www.math.bgu.ac.il/~pakovich/Publications/FPYZ.pdf
was inspected for this comparison. This is classical problem context,
not an additional premise or a claim about what that article does not prove.

The all-input and uniform termination claims are proved using the cited
regular-telescoping theorem, formal leading-term algebra, and effective
commutative algebra. The finite checker verifies the displayed three-height
certificate and monomial/normalization fixtures. It does not execute the
universal minimal-prime algorithm or compute a useful numerical B(d).

The next quantitative problem is a practical certified telescope backend
and a usable upper bound on the coefficient-space algorithm. Picard--Fuchs
D-finiteness itself is classical. Neither ordinary-point ODE uniqueness,
a rational-diagonal assertion, nor raw operator order/degree bounds
supply the missing numerical cutoff. For example p=t^2 already gives
G_p(z)=atanh(sqrt(z))/sqrt(z), with coefficients 1/(2m+1); it need not be
a diagonal of a rational function regular at the origin over Q. Indeed the coefficients of any such rational function, and hence its diagonal, belong to Z[S^-1] for a fixed finite set of primes S: clear the finitely many coefficient denominators and invert the nonzero constant denominator coefficient. The numbers 1/(2m+1) involve every odd prime in their denominators, so no such fixed S exists.

Primary period source inspected: Bostan--Lairez--Salvy arXiv1301.4313v2,
Proposition 4 and Theorem 12, including its singular-denominator deformation.
Local primary PDF: [retained Bostan--Lairez--Salvy source](sources/bostan_lairez_salvy_1301_4313v2.pdf),
SHA256 c7d666525392f58a13d1fb5b8b77805e612391c069c07f2f8c53c482d4d40886.
The closed-tube reduction and the singular leading-term lemma are proved
above rather than imported as unstated period or ODE hypotheses.

Proof support: [source ledger](global_sphere_purity_source_ledger.json), [exact certificate input](global_sphere_three_height_certificates.json), [checker](verify_global_sphere_three_height.py), [receipt](global_sphere_three_height_verification.json), and [formal indicial proofs](FORMAL_INDICIAL_PURITY_CERTIFICATES.md). All reviews and replays recorded here are internal; no formal verification, external human review, practical universal backend or historical-priority claim is asserted.
