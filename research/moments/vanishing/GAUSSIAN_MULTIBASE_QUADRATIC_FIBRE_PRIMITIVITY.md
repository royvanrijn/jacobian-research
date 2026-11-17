# Primitive vacua for multibase polynomial quadratic Gaussian fibres

For this architecture, every actual marked image finite over K[t] has
a primitive original vacuum. There is an exact all-degree branch
criterion for that finite image, and safety implies purity as a corollary.

## Finite actual images and exact safety criterion

Let K be a field of characteristic zero. Use n independent Witt pairs
(Z_i,W_i), E[Z^alpha W^beta]=delta_(alpha,beta)alpha!, and s independent
ordinary standard Gaussian variables Y. Let

    P=W^T H(Z)+Y^T A(Z)Y+B(Z)^T Y+G(Z),              (1)

where H is an arbitrary polynomial n-vector, A an arbitrary symmetric
polynomial s by s matrix, B a polynomial s-vector, and G a polynomial.
The crucial hypothesis is that H depends on Z alone, not on Y or W.
The cases n=0 or s=0 are allowed with the empty-matrix conventions.

Let r(t) in tK[[t]]^n be the unique root r=tH(r), and put

    J(t)=I-t DH(r(t)),
    D(t)=I-2t A(r(t)),
    R(t)=tG(r(t))+t^2 B(r(t))^T D(t)^(-1)B(r(t))/2.

Let V=Phi[K[Z,W,Y]] be the literal vector image of original fixed markers
and M=K[t]V its polynomial-time span. The principal finite-image theorem is

    M is finite over K[t]
      iff r in K[t]^n and detJ=detD=1.                 (F)

These conditions automatically make R polynomial: D^(-1)=adj(D),
and every remaining entry in its definition is polynomial in t.

In that case

    f=exp(R), M=exp(R)K[t].

Thus the distinguished original unit generates the whole lattice:
M=fK[t]. If P is nonconstant, V=M as literal K-linear spaces by the
retained general fixed-marker-image theorem. For a constant phase
P=lambda, the exact exception is V=K exp(lambda t), M=exp(lambda t)K[t].

Safety holds precisely when the conditions in (F) hold and R=0.
Equivalently, the following are equivalent:

1. P is safe for EVERY original polynomial marker.
2. Every base-only marked period Phi[F(Z)] is polynomial in t.
3. The scalar f=Phi[1] is polynomial, of some degree L, and, for every
   i, both Phi[Z_i] and Phi[Z_i^(L+1)] are polynomial.
4. The following exact branch conditions hold:

       r(t) in K[t]^n,
       det J(t)=1,   det D(t)=1,   R(t)=0.           (2)

Every such phase is pure: f=1. Thus safety implies purity throughout
this unrestricted multibase architecture. This is a classification by
its formal source branch, not an exhaustive normal form for H,A,B,G.
It does not prove general field safety implies purity.

When s=0, D is empty and R=tG(r). When A is constant, detD=1 says A is
nilpotent, but for variable A the condition is only along r(t). A need
not be nilpotent away from that branch.

The tests in item3 concern entire marked series, not finitely many
moment coefficients. For L=0 they reduce to the unit and n linear base
markers. The bound L+1 is sharp in the retained polynomial-seed family.

## 1. A self-contained multibase Wick contraction

Introduce sources v in K^n for W and eta in K^s for Y. Let r(t,v) be
the unique t-adic polynomial-coefficient root

    r=tH(r)+v,        r(0,v)=v.

Write J(t,v)=I-t DH(r(t,v)), D(t,v)=I-2t A(r(t,v)). Then for every
polynomial base amplitude F,

    E[F(Z)exp(tP+v^T W+eta^T Y)]
      = F(r(t,v)) exp(E(t,v,eta))
          / [det J(t,v) sqrt(det D(t,v))],           (3)

where the square root has constant coefficient one and

    E(t,v,eta)=tG(r)
                +(tB(r)+eta)^T D(t,v)^(-1)(tB(r)+eta)/2. (4)

The identity is formal in t,v,eta. Every fixed coefficient uses finitely
many Gaussian contractions. Differentiating (3) in v and eta gives ALL
original W,Y markers.

Here is a complete formal-residue proof, without invoking a multivariate
Good theorem or a theorem stated only for one variable. For a numerator
Q(t,Z) polynomial in Z at each t order, direct circular Wick contraction
in each pair gives

    E_ZW[Q(t,Z)exp(W^T(tH(Z)+v))]
      =Res_Z Q(t,Z) product_i (Z_i-tH_i(Z)-v_i)^(-1) dZ. (5)

Indeed expand each denominator geometrically: its term indexed by m_i
is (tH_i+v_i)^(m_i)/Z_i^(m_i+1), exactly the coefficient obtained from
E[Z_i^(m_i)W_i^(m_i)]=m_i!. At fixed t,v order the sums are finite.

For clarity, formal residue change of variables for U=Z-tH(Z) can be
proved directly. Its inverse Z=r(t,U) has polynomial U coefficients
at every t order. Residue is invariant under this substitution. It is
enough to prove this for a Laurent monomial top form in U, coefficientwise
in t. If some exponent k_i differs from -1, that form is exact: integrate
U_i^(k_i) as U_i^(k_i+1)/(k_i+1). Pullback preserves exactness, and the
residue of an exact Laurent top form is zero. For the exceptional form
product_i dU_i/U_i, write

    dU_i/U_i=dZ_i/Z_i+d log(1-tH_i(Z)/Z_i).

The logarithm is a legitimate t-adic Laurent series. Every difference
term in the wedge product contains an exact dlog factor and the other
factors are closed, so it too is exact and has zero residue. The sole
remaining residue is that of product_i dZ_i/Z_i, namely one. This proves
formal residue invariance for all relevant series by linearity; all
coefficientwise operations remain finite.

Apply this substitution to (5). The Jacobian of the inverse is
1/det(I-tDH(r(t,U))). The resulting residue is

    Res_U [Q(t,r(t,U))/det J(t,U)] product_i(U_i-v_i)^(-1)dU
      =Q(t,r(t,v))/det J(t,v),

by geometric expansion in v. This proves the multibase circular identity.
Ordinary finite-dimensional Gaussian square completion in Y gives
(3)-(4), including the determinant square root with value one at t=0.
That Gaussian formula also follows directly by Wick contraction and
formal differentiation in eta; it uses no convergence assumption.

In particular, at v=eta=0,

    f=exp(R)/(detJ sqrt(detD)),
    Phi[Z_i^m]=f r_i^m, EVERY m>=0.                 (6)

## 2. Base powers force a polynomial source root

Safety implies every series in (6) is polynomial. More generally,
assume the finite all-order tests in item3. Since f(0)=1,

    r_i=Phi[Z_i]/f

is rational and regular at zero. If r_i has any nonconstant denominator
q in lowest terms, polynomiality of f r_i^(L+1) forces q^(L+1) to
divide f. This is impossible because deg f=L. Hence every r_i is
polynomial. This proves the source-root conclusion from exactly the
stated base tests, without any momentum-marker argument.

## 3. Rational exponential and determinant rigidity

Once r is polynomial, J,D have polynomial entries and both determinants
have value one at t=0. Also R is rational, regular at zero, with R(0)=0.
Squaring the unit identity in (6) yields

    exp(2R)=f^2 (detJ)^2 detD in K[t].                (7)

A rational function T regular at zero cannot have a nonconstant rational
exponential. To see this algebraically, differentiate exp(T)=S in K(t):
T'=S'/S. Over an algebraic closure, a finite pole of T would give a pole
of order at least two on the left, while logarithmic derivatives have
only simple poles. Thus T is polynomial; at infinity S'/S=O(1/t), so
T'=0. Characteristic zero then makes T constant.

Apply this fact to (7). Its value at zero gives R=0 and

    f^2 (detJ)^2 detD=1.

Every factor is a polynomial. Each is therefore a polynomial unit;
normalization at zero gives f=detJ=detD=1. This proves both purity and
all the conditions in (2). No assumption of positivity, actual rank one,
formal image-kernel identification, or a bounded moment prefix is used.

## 4. The branch conditions control EVERY original marker

Assume (2). Every v derivative r_beta(t)=partial_v^beta r(t,v)|_0 is
polynomial in t. For first derivatives this is J^(-1), a polynomial
matrix by detJ=1. Differentiating r=tH(r)+v at higher v order gives
J r_beta equal to a polynomial expression in t, derivatives of H at
r(t), and lower-order r_gamma. Induction proves the assertion.

Consequently every fixed v jet of A(r),B(r),G(r), and J(t,v) is polynomial.
The inverse-matrix identity D(t,v) D(t,v)^(-1)=I recursively gives
polynomial inverse jets because D(t,0)^(-1)=adj(D(t,0)). The scalar
inverse determinant J jets are polynomial since detJ(t,0)=1. Likewise
all jets of (detD(t,v))^(-1/2) are polynomial by the finite binomial
expansion in detD(t,v)-1, which has zero v-constant coefficient.

The exponent (4) has polynomial fixed source jets and value R(t)=0
at v=eta=0. Thus every fixed derivative of its exponential is a finite
polynomial expression in those jets. The same holds after multiplying
by r(t,v)^alpha and the two determinant factors. By (3), every original
monomial marker Z^alpha W^beta Y^gamma has polynomial period. Linear
combination proves safety for every original polynomial marker.

This also supplies an exact original-marker formula and a recursive
algorithm for its period. It controls arbitrarily high marker order;
it is not an inference from finitely many tested jets.

## 5. Every finite actual image has primitive vacuum

Here is the extra argument proving (F), beyond the safe case. The
general finite-image classification is retained in the
[amplitude-ideal note, section7](GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md#7-every-finite-actual-image-and-its-original-fixed-marker-image). Its short
source/PID argument is as follows. The exact integrated Gaussian source
H specializes at t=0 to its one-dimensional unperturbed Gaussian
quotient. Its period surjection H->M gives K->M/tM. A nonzero finite
torsion-free K[t] module is free and has rank dim_K(M/tM), hence M has
rank one. Choose its generator g with g(0)=1. Derivative stability gives
g'=a(t)g for a polynomial a. Integration in characteristic zero yields

    M=exp(psi)K[t], f=p exp(psi),
    psi in tK[t], p in K[t], p(0)=1.                (9)

The t=0 source quotient used here follows directly from Gaussian
Hermite reduction, or the polynomial Stein divergence identity:
every amplitude minus its Gaussian expectation lies in the sum of the
unperturbed spatial Stein images. The expectation kills those images
and takes1 to1. No actual-image kernel identification is required.

By (6), every normalized base power satisfies p r_i^m in K[t]. The
same denominator argument as section2, with m=deg(p)+1, forces every
r_i to be polynomial. Thus J,D are polynomial matrices and R is rational.
Substitute (9) into the exact unit formula and square:

    exp(2(R-psi))=p^2(detJ)^2 detD.

The rational-exponential argument forces R=psi, and the polynomial
product is1. Therefore p=detJ=detD=1. This proves necessity of (F)
and primitive vacuum, without a prior generic-rank assumption.

Conversely assume the branch conditions in (F); R is then automatically
polynomial. The source-jet proof
of section4 applies after factoring exp(R(t)) out of (3). Every original
marker then has period exp(R) times a polynomial. Since the original
unit has period exp(R), its K[t] span is exactly M=exp(R)K[t], a finite
image. This proves (F). It also gives the safety criterion: exp(R) is
polynomial exactly when R=0.

This result rules out nonprimitive finite actual Gaussian images in
this entire architecture. It does not exclude them for phases with
nonlinear momentum dependence or nonquadratic ordinary fibres.

## 6. Optional centered-momentum identity

There is a useful independent pole diagnostic. For coordinate indices
i,j and every integer m>=0,

    Phi[(Z_i-r_i(t))^m W_j^m]
       =m! f(t) (J(t)^(-1)_(ij))^m.                 (10)

The amplitude has polynomial t coefficients once r is polynomial, so
safety applies to it by finite K[t]-linear combination of original
markers. To derive (10), differentiate the source expression for
(Z_i-r_i(t))^m exactly m times in v_j. Every factor vanishes at v=0,
so the only surviving term differentiates each factor once. All
remaining scalar/source derivatives are multiplied by a zero factor.
This identity forces inverse-Jacobian pole ladders under safety. The
shorter necessity proof above already obtains detJ=1 from (7).

## 7. Exact examples and rejection boundaries

A nonlinear triangular example uses n=2,

    H_1=1+Z_2^2, H_2=1,
    r=(t+t^3,t),
    F(Z)=Z_1-Z_2-Z_2^3.

Here detJ=1 and F(r)=0. With s=0 and G=F, (2) proves pure safety.
The source root is explicit:

    r_2=t+v_2, r_1=v_1+t(1+(t+v_2)^2),
    F(r)=v_1-v_2-v_2(t+v_2)^2.

Every marked period follows from polynomial source differentiation.
The exact cutoff deg_t Phi[Q]<=3 deg Q holds, and Phi[Z_1^d]=(t+t^3)^d
attains it. This is an all-order sharp control, not a finite experiment.

Variable nonnilpotent fibre matrices are allowed. For the same H choose
s=1, A(Z)=F(Z), B(Z)=1, G(Z)=F(Z)-Z_2/2. Along the source branch,
D=1 and R=-t^2/2+t^2/2=0. Thus the phase is pure and safe even though
A(Z) is nonzero and not globally nilpotent.

A second exact matrix example has one pair Z,W and two ordinary
variables, with H=1, B=G=0 and

    A(Z)=[[0,-1/2],[-1/2,-Z/2]].

Then r=t,J=1 and

    D(t)=[[1,t],[t,1+t^2]],
    detD=1, D^(-1)=[[1+t^2,-t],[-t,1]].

Thus P=W-Y_1Y_2-ZY_2^2/2 is pure and safe. The matrix A(Z) is not
pointwise nilpotent. This checks why the criterion concerns the branch
determinant rather than a global nilpotence condition on A.

Three rejection controls distinguish the conditions:

* H=1+Z,G=A=B=0 has r=t/(1-t), hence a nonpolynomial root.
* H=Z,G=A=B=0 has r=0 but detJ=1-t and inverse-Jacobian poles.
* H=1,G=Z^2,A=B=0 has r=t,J=1 but R=t^3 and seed exp(t^3).

The retained Long phase is an additional important boundary:

    H=1+Z, A(Z)=-(1+Z)(2+Z)/2, B=G=0.

It has f=1, but r=t/(1-t); its linear base marker already proves it
unsafe. Thus the theorem does not infer safety from purity alone.

For the retained polynomial-seed compensator, f=(1-t)^M and its base
root contains r_Y=t/(1-t). Then f r_Y^j=t^j(1-t)^(M-j) is polynomial
for j<=M, but the base marker of degree M+1 is nonpolynomial. This
proves sharpness of L+1 in item3 for every L=M>=1. That compensator's
phase is itself of the architecture (1): the three isotropic momenta
are A,W,V, and its remaining ordinary variable occurs quadratically.

## 8. Spectators, classical mechanisms, and scope

Adding unused Gaussian variables, including an unused Witt pair, preserves
safety and purity by tensor-factorization of every polynomial marker.
Adding a used pair with its new H component equal to1 is already within
(1), and can encode polynomial dependence on the corresponding root t.
The theorem allows an arbitrary polynomial base potential G and arbitrary
polynomial A,B. It does not permit H to depend on the ordinary Y fibres,
or higher Y degree, or nonlinear dependence on the momenta W.

Over C, every original marked period in this architecture converges
near t=0. The root r(t,v) is analytic by the implicit-function theorem
at the identity Jacobian, both determinant factors are analytic units,
and any fixed source derivative in (3) is analytic. This gives an exact
obstruction to a universal finite auxiliary-Gaussian encoding preserving
the original unit: for P=(ZW)^2 the unit is
sum_m (2m)!t^m/m!, whose successive coefficient ratio is 2(2m+1), so
its radius is zero. Such a phase cannot have its original unit realized
by any phase of (1), even after adjoining finitely many variables.
This does not obstruct an encoding restricted to finite actual images.

Formal residues, Gaussian square completion and the PID argument are
classical mechanisms. The proof above derives the multibase identity
directly; a one-variable Lagrange theorem is not used as a multibase
result. No external novelty or priority claim is made. The finite-image
argument uses the retained general source/PID classification, and the
explicit original-marker formula retains the distinguished vacuum.

The [checker](verify_gaussian_multibase_quadratic_fibre_primitivity.py) and
[receipt](gaussian_multibase_quadratic_fibre_primitivity_verification.json)
retain the finite exact controls. Replay from the repository root with
`.venv/bin/python research/moments/vanishing/verify_gaussian_multibase_quadratic_fibre_primitivity.py --output research/moments/vanishing/gaussian_multibase_quadratic_fibre_primitivity_verification.json`.
They independently compare the source formula to original
Wick contractions, check the nonlinear triangular and variable-fibre
examples, and retain the three unsafe branches plus Long's pure unsafe
boundary. Those checks support the identities; sections1-5 are the
unrestricted proof. Independent replay, external review and formal
verification are not claimed.

## 9. A separate mixed-source extension with zero root

The following family is outside (1), because its isotropic-momentum
coefficient depends on the ordinary fibre. Let K have characteristic zero,
use one Witt pair Z,W and ONE ordinary standard Gaussian Y, and put

    P=W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),
    h(0)=c(0)=0.

Here a is constant and h,c,u,k are arbitrary polynomials. There is no
degree or root-multiplicity bound. For the FULL actual marked image,

    M is finite over K[t]
      iff a=0, h'(0)=c'(0)=0.                       (Z1)

Thus h,c must lie in z^2 K[z], while u,k are unrestricted. In that case,
writing u0=u(0), k0=k(0),

    psi=k0 t+u0^2 t^2/2,
    f=exp(psi), M=fK[t].                            (Z2)

The original vacuum is primitive. For nonconstant P, also V=M by the
general fixed-marker theorem; for constant P=k0, V=K exp(k0t).
Safety, purity and polynomiality of f are equivalent here. They hold
exactly when (Z1) holds and u0=k0=0.

Here is a self-contained scalar derivation. Set eta=h'(0), gamma=c'(0).
The angular-weight-zero part of P is

    Pzero=ZW(eta+gamma Y)+aY^2+u0Y+k0,               (Z3)

and its remainder has strictly positive Z-minus-W weight. Every unit
power containing a positive-weight factor has zero circular contraction.
Since E[(ZW)^n]=n!, summing the zero-weight part gives the actual identity

    f=E_Y[exp(t(aY^2+u0Y+k0))/(1-t(eta+gamma Y))].

Put

    delta=1-2at, mu=t u0/delta,
    Dlin=1-t eta-t gamma mu,
    q=t gamma/[sqrt(delta)Dlin],
    R0=t k0+t^2u0^2/(2delta).

Ordinary Gaussian completion, with square-root constant term one, yields

    f=exp(R0)/[sqrt(delta)Dlin] * M0(q^2),
    M0(s)=sum_(n>=0)(2n-1)!!s^n, (-1)!!=1.          (Z4)

Indeed the shifted denominator is Dlin(1-qS); its geometric expansion
has only finitely many terms at each t order. This is a formal Wick
identity, not a convergence assertion about the original integral.

Suppose M is finite. The source/PID argument of Section 5 gives

    f=p exp(psi), psi in tK[t], p in K[t], p(0)=1,   (Z5)

without assuming vacuum primitivity. Over C, the scalar prefactor in
(Z4) is analytic and nonzero near zero. If gamma!=0, q is analytic with
q(0)=0 and q'(0)=gamma, so it has an analytic inverse t=t(q). Equation
(Z5) would make M0(q^2) analytic. This is impossible: its q^(2n)
coefficient (2n-1)!! has unbounded 2n-th root, since it is at least n!.
Thus gamma=0.

For an arbitrary characteristic-zero K, adjoin to Q the finitely many
coefficients of P,p,psi. The resulting finitely generated field embeds
in C. Every coefficient of (Z4)-(Z5) lies in that field, so the entire
formal identities specialize under this embedding, with gamma remaining
nonzero if it was nonzero. The embedded p exp(psi) is an analytic germ,
contradicting the preceding argument. This proves gamma=0 over K without
assuming that K itself has a complex embedding.

With gamma=0, square (Z4)-(Z5):

    exp(2(R0-psi))=p^2(1-t eta)^2(1-2at).           (Z6)

The rational-exponential rigidity of Section 3 forces its exponent to
vanish and its polynomial right side to be one. Normalization at zero
therefore gives p=1, eta=a=0 and psi=R0=k0t+u0^2t^2/2. This proves
necessity of (Z1) and primitive vacuum. If only f is polynomial, the
same proof takes psi=0,p=f and additionally forces u0=k0=0; thus a
polynomial unit period is already pure in this family.

For sufficiency and EVERY original marker, assume (Z1) and split

    P=P0+Pplus,
    P0=u0Y+k0,
    Pplus=W(h+cY)+(u-u0)Y+(k-k0).

Every monomial of Pplus has angular weight at least one. For an original
marker Q of degree D, terms Q Pplus^n with n>D cannot contract. Gaussian
translation in Y consequently gives the exact finite formula

    Phi[Q]=exp(psi) sum_(n=0)^D t^n/n!
       E[Q(Z,W,Y+t u0)Pplus(Z,W,Y+t u0)^n].          (Z7)

Thus all original marked periods are exp(psi) times polynomials. The
unit is exp(psi), proving (Z2) and sufficiency of (Z1).

This retains a uniform original-marker degree bound:

    deg_t(exp(-psi)Phi[Q])<=2 deg Q,
    deg_t(exp(-psi)Phi[Q])<=deg Q if u0=0.

For a monomial Q=Z^i W^j Y^l, a nonzero contraction requires n<=j-i.
Its shifted Y degree is at most l+n, so the resulting time degree is
at most l+2n<=l+2j<=2 deg Q. Without the mean shift it is at most n.
The first coefficient two is sharp: for u0!=0 and P=u0Y+Z^2WY,

    exp(-u0^2t^2/2)Phi[W^D]
       =(2D)!/D! * t^D E[(Y+t u0)^D],               (Z8)

of exact degree 2D with leading coefficient (2D)!u0^D/D!. In the safe
case P=Z^2WY, even D gives exact degree D, so the safe bound is sharp
as well. These are all-D identities, not finite-prefix controls.

Finally exp(psi) is polynomial exactly when psi=0, which over a field
means u0=k0=0. Equations (Z2),(Z7) prove the asserted safety/purity/unit-
polynomial equivalence. The single ordinary fibre is an explicit scope
condition: with several ordinary fibres a nonzero source vector can be
isotropic, so the implication gamma=0 would require a different argument.
This extension does not allow arbitrary Y-dependent H in (1), nonlinear
W dependence, or remove the hypotheses h(0)=c(0)=0.

## 10. A quadratic zero-weight part with arbitrary positive remainder

Let K have characteristic zero. Use one Witt pair Z,W and s ordinary
standard Gaussian coordinates Y, and put S=ZW. Consider

    P=Pzero+Pplus,
    Pzero=lambda S^2+eta S+S gamma^T Y+Y^T A Y+u^T Y+k, (Q1)

where A is a constant symmetric matrix and all other displayed
coefficients are constant. EVERY monomial of Pplus has strictly positive
Z-minus-W weight; apart from this, Pplus is an arbitrary polynomial.
Its W and Y degrees are unrestricted. Exchanging Z,W covers negative
angular support.

The full actual M is finite over K[t] exactly when

    lambda=eta=0, A is nilpotent,
    gamma^T A^j gamma=gamma^T A^j u=0,
                      j=0,...,s-1.                 (Q2)

Under these conditions let nu>=1 satisfy A^nu=0 and set

    C(t)=(I-2tA)^(-1)=sum_(j=0)^(nu-1)(2tA)^j,
    psi=kt+t^2 u^T C(t)u/2.                        (Q3)

Then the original vacuum is primitive:

    f=exp(psi), M=fK[t].                            (Q4)

For nonconstant P, also V=M literally; for constant P=k, V=K exp(kt).
Safety holds exactly when (Q2) holds and

    k=0, u^T A^j u=0, j=0,...,s-1.                 (Q5)

Every safe member is pure. Purity alone need not imply safety in this
vector extension, as the exact boundary below shows. Empty-matrix
conventions cover s=0. This theorem concerns a quadratic zero-weight
part and arbitrary strictly positive remainder, not arbitrary mixed
angular support.

### The actual unit and the divergent radial series

Unit and Y-only markers ignore Pplus: every term containing a positive-
weight factor has zero circular contraction. Put

    D=I-2tA, Delta=det D, C=D^(-1),
    sigma2=gamma^T C gamma,
    R=kt+t^2 u^T C u/2,
    L=1-t eta-t^2 gamma^T C u,
    g=(2t lambda+t^2 sigma2)/L^2.

Ordinary Y completion leaves the radial functional E_S[S^m]=m!. Its
exact identity

    E_S[exp(aS^2+bS)]=(1-b)^(-1)M0(2a/(1-b)^2),
    M0(x)=sum_(n>=0)(2n-1)!!x^n

follows by expanding in a, summing
sum_m(2n+m)!b^m/m!=(2n)!(1-b)^(-2n-1), and using
(2n)!/n!=2^n(2n-1)!!. Consequently the ACTUAL unit is

    f=exp(R)/[sqrt(Delta)L] * M0(g).                (Q6)

Every operation is coefficientwise finite; square roots have value one
at zero. If M is finite, the general source/PID result gives

    f=p exp(psi), p in K[t], p(0)=1, psi in tK[t],   (Q7)

without assuming primitive vacuum.

Over C, all scalar factors outside M0 in (Q6) are analytic units at
zero. If g is nonzero, write g=b t^r w(t), b!=0, r>=1, w(0)=1. The
analytic change tau=t w(t)^(1/r) has derivative one and g=b tau^r.
Equation (Q7) would make M0(b tau^r) analytic. Its coefficient at degree
rn is b^n(2n-1)!!, of unbounded rn-th root. Thus g=0. Its first t
coefficient is 2lambda, so lambda=0 and sigma2=0 identically.

For arbitrary K, adjoin to Q the finitely many coefficients of the phase,
p and psi, and embed that finitely generated field in C. All exact
formal identities specialize, while a nonzero rational g stays nonzero.
The embedded p exp(psi) is analytic, giving the same contradiction.
This proves the conclusion over every characteristic-zero field.

### A bounded-pole lemma using every Y marker

Let p in K[t] be nonzero, m,v rational vectors, and Sigma a symmetric
rational matrix. If EVERY xi coefficient of

    p exp(m^T xi+xi^T Sigma xi/2)/(1-v^T xi)          (Q8)

is polynomial in t, then all entries of m,v,Sigma are polynomial.

Work at a finite point over an algebraic closure, with valuation val.
All nonzero constants and factorials are units. If N=val(p), the fraction
in (Q8) has coefficients uniformly bounded below by -N. In one variable,
write A(x)=exp(mx+sigma x^2/2)/(1-vx). Multiplication by 1-vx preserves
some uniform lower coefficient bound, so E(x)=exp(mx+sigma x^2/2) and
E(-x) are uniformly bounded. Their product exp(sigma x^2) is uniformly
bounded. Its x^(2n) coefficient has valuation n val(sigma), forcing
val(sigma)>=0. Now exp(-sigma x^2/2) is integral, so multiplying E by
it gives bounded coefficients of exp(mx), forcing val(m)>=0. The inverse
E is then integral; multiplying A by it leaves the bounded geometric
series (1-vx)^(-1), forcing val(v)>=0. Every coefficient convolution is
finite; no convergence or finite truncation is used.

For several variables restrict xi to x e_i and x(e_i+e_j). The first
restrictions make m_i,v_i,Sigma_ii integral; the second make
Sigma_ii+2Sigma_ij+Sigma_jj integral and hence Sigma_ij integral. Thus
there are no finite poles, proving the lemma.

With lambda=sigma2=0, adding the ordinary source xi^T Y gives exactly

    Phi[exp(xi^T Y)]
      =f exp(m^T xi+xi^T C xi/2)/(1-v^T xi),
    m=t C u, v=t C gamma/L.                        (Q9)

The conditional linear source has zero variance gamma^T C gamma, so its
geometric Gaussian series evaluates at its mean. Equivalently, complete
the square and replace L by L-t gamma^T C xi in (Q6).

Every xi coefficient of exp(-psi) times (Q9) is polynomial, by finite M.
The bounded-pole lemma with prefactor p forces C,m,v polynomial. Since
C=D^(-1) is polynomial, Delta is a polynomial unit and equals one; hence
A is nilpotent. Now L and R are polynomials. Squaring the unit identity
against (Q7) gives exp(2(R-psi))=p^2 Delta L^2. The field rigidity in
Section 3 forces R=psi and p=L=1. From L=1-t eta-t^2 gamma^T C u one
gets eta=0 and gamma^T C u=0. Expanding the finite Neumann series for C
turns this and sigma2=0 into the finite chains (Q2). This proves necessity
and primitive vacuum from ALL original Y markers, rather than the scalar
unit alone.

### Sufficiency for all original markers

Assume (Q2). At fixed S, the ordinary Y mean is mu_S=t C(u+S gamma)
and covariance is C, both polynomial in t,S. Its scalar exponent is

    R+S(t eta+t^2 gamma^T C u)
       +S^2(t lambda+t^2 gamma^T C gamma/2)=R,

and the determinant factor is one. Denote the resulting finite Wick
polynomial of an amplitude F by

    W_(C,mu_S)(F)
      =[exp((1/2)sum C_ij partial_Yi partial_Yj)F]_(Y=mu_S).

It preserves angular weight because S has weight zero. For an original
marker Q of degree d, at most d positive-weight factors can survive.
Thus the exact all-marker formula is

    Phi[Q]=exp(R) sum_(n=0)^d t^n/n!
             E_ZW[W_(C,mu_S)(Q Pplus^n)]_(S=ZW).    (Q10)

Every term is polynomial in t after removing exp(R); all S powers are
finite and circular contraction uses only finitely many factorials.
The unit itself is exp(R), proving (Q4) and sufficiency.

A sufficient original-marker bound is also explicit. Let b be the maximum
Y degree of Pplus, taking b=0 if Pplus=0. In Z^i W^j Y^alpha one has
n<=j-i, and the Wick time degree is at most nu(|alpha|+bn). Hence

    deg_t(exp(-R)Phi[Q])
       <=max(nu,1+nu b) deg Q.                     (Q11)

No optimality is asserted for every phase. Finally exp(R) is polynomial
exactly when R=0; its coefficients in (Q3) give precisely (Q5). This
proves the safety criterion and safety implies purity.

### A pure unsafe boundary with its full original image

Over C take s=2, A=I_2, gamma=(1,i)^T, eta=-2, u=(-2,2i)^T and
lambda=k=Pplus=0. With d=1-2t one gets Delta=d^2, sigma2=R=0, L=1/d,
so f=1. But Phi[Y1]=-t(1+2t)/d=1+t-d^-1 is not polynomial.
Its entire original image is

    V=M=K[t,d^-1].                                  (Q12)

Here is the all-marker proof. For every a>=0,

    Phi[S^a exp(xi^T Y)]
      =a!d^a exp(m^T xi+(xi_1^2+xi_2^2)/(2d))
                        /(1-t(xi_1+i xi_2))^(a+1),
    m=t(-2,2i)^T/d.

Each fixed xi coefficient has poles only at d=0. Unbalanced Z,W markers
have zero period, so V is contained in K[t,d^-1]. Conversely
Phi[S^a]=a!d^a spans every polynomial in t. The fixed marker
B=(3-S)/2-Y1 has Phi[B]=d^-1, and actual phase multiplication gives

    Phi[P^n B]=partial_t^n(d^-1)=2^n n!d^(-n-1).

These fixed original markers realize every pole order, proving (Q12).
This example is not a counterexample to safety implies purity. It shows
why the all-Y-marker pole lemma is essential in the vector theorem; the
unit formula alone does not force A nilpotent.

## 11. Integer-Gamma radial shapes and an exact cubic delay family

Section 10 extends to any POSITIVE INTEGER r of independent Witt pairs.
Set

    S=sum_(i=1)^r Z_i W_i,

and use the same quadratic Pzero in (Q1). Every monomial of Pplus must
have strictly positive TOTAL weight sum_i(deg_Zi-deg_Wi), with no further
restriction. The finite-image criterion is still exactly (Q2), its
primitive image is (Q3)-(Q4), and safety additionally requires (Q5).
This is an original Gaussian realization in 2r+s variables; no assertion
about arbitrary Gamma shapes without such a realization is made.

Here are the changes in the proof. Independent circular contraction gives
E[exp(xS)]=(1-x)^(-r) and E[S^m]=(r)_m. With the notation preceding (Q6),
the actual unit is

    f=exp(R)/[sqrt(Delta)L^r] * G_r(g),
    G_r(x)=sum_(n>=0)(r)_(2n)x^n/(2^n n!).          (G1)

Indeed E[S^(2n)exp(bS)]=(r)_(2n)(1-b)^(-r-2n). The ratio of successive
G_r coefficients is (r+2n)(r+2n+1)/(2(n+1)), tending to infinity.
The same order-straightening and characteristic-zero field descent as
Section 10 therefore force g=0 under finite M, hence lambda=sigma2=0.
The ordinary-Y source then becomes

    f exp(m^T xi+xi^T C xi/2)/(1-v^T xi)^r,
    m=tCu, v=tCgamma/L.                            (G2)

The bounded-pole proof multiplies by the fixed polynomial (1-vx)^r to
obtain the uniformly bounded exponential. The product with its reflection
still forces the variance, then mean, integral. The remaining geometric
series has coefficients binom(r+n-1,n)v^n, with nonzero characteristic-
zero scalar factors, so v is integral too. Thus C is polynomial and
Delta=1. Finally exp(2(R-psi))=p^2L^(2r) forces p=L=1 and gives (Q2).

For sufficiency, the polynomial ordinary-Y covariance C and mean
 tC(u+Sgamma) again have scalar exponent independent of S. Only finitely
many positive-TOTAL-weight factors can contract with any fixed original
marker. The remaining operation is finite circular contraction in all r
pairs. This proves all-original-marker sufficiency and the same bound
(Q11), without requiring positivity in each individual pair separately.

### Polynomial units of arbitrary degree from cubic phases

Over C use two ordinary fibres and, for every r>=1, take

    A=I_2, gamma=(1,i)^T, eta=-2, u=(-2,2i)^T,
    lambda=k=Pplus=0,
    P_r=-2S+S(Y1+iY2)+Y1^2+Y2^2-2Y1+2iY2.

Every P_r has total degree THREE in 2r+2 Gaussian variables. Put d=1-2t.
Then Delta=d^2, sigma2=R=0 and L=1/d. Thus

    Phi_(P_r)[1]=d^(r-1).                           (G3)

It is pure only for r=1. None of these phases is safe. The complete
original vector and polynomial-time images are

    V=M=K[t,d^-1].                                  (G4)

There is also an exact original O-module marker filtration. For
F_D=K[t]-span{Phi[Q]:deg Q<=D}, one has for EVERY D>=0

    F_D=d^(r-1-D)K[t].                              (G5)

In particular every original marker of degree at most r-1 has polynomial
period, while the least bad original marker degree is exactly r.

To prove all statements, unbalanced markers in any Witt pair have zero
period. For a multi-index alpha let a=|alpha|. The full balanced source is

    Phi[(prod_i(Z_iW_i)^alpha_i)exp(xi^T Y)]
      =(prod_i alpha_i!) d^(r+a-1)
          exp(m^T xi+(xi_1^2+xi_2^2)/(2d))
                    /(1-t(xi_1+i xi_2))^(r+a),
    m=t(-2,2i)^T/d.                                (G6)

This follows by differentiating the product radial geometric series
before Y completion. Every fixed Y jet is rational with poles only at
d=0. For Y degree q, its period lies in d^(r+a-1-q)K[t]. The total
original degree is 2a+q, giving the inclusion in (G5).

For the reverse inclusion write

    Phi[Y1^D]=d^(r-1-D)q_D(t).

The polynomial q_D satisfies q_D(1/2)=(-1)^D!=0. Its highest pole comes
uniquely from taking all D derivatives of the linear mean -2t/d;
covariance factors use two derivatives for only one d denominator, and
the geometric factor has no d denominator. Now F_D contains both this
period and Phi[1]=d^(r-1). After factoring d^(r-1-D), their remaining
factors are q_D and d^D, which are coprime. Bezout over K[t] gives the
reverse inclusion, for negative and nonnegative r-1-D alike.

For the literal image (G4), radial markers give
Phi[S^a]=(r)_a d^(r+a-1), polynomials of unbounded degree. The polynomial
part of V is differentiation-stable by actual phase multiplication;
the triangular derivative argument of the original-image theorem gives
K[t] contained in V. The marker Y1^r has nonzero simple principal part
(-1)^r/d. Subtracting its polynomial part, already realized by fixed
markers, gives d^-1 in V. Derivatives then realize every d^-n by fixed
markers obtained through multiplication by powers of P_r. Formula (G6)
gives the opposite containment, proving (G4).

This is an all-order, fixed-cubic-degree and growing-dimension polynomial-
seed family with a complete original image and exact marker delay. It
is not a counterexample to safety implies purity, and no novelty or
external-review claim is made.
