# Polynomial unit periods, exact marker filtrations and pure Airy delays

These are internal all-degree proofs with finite exact controls.
Independent replay, formal verification, external review and priority are
not claimed. The [constructor](gaussian_seed_filtration.py),
[checker](verify_gaussian_seed_filtration.py),
[certificate](gaussian_seed_filtration_certificate.json) and
[receipt](gaussian_seed_filtration_verification.json) transcribe the formulas;
the universal assertions use the proofs below.

## 1. Arbitrary polynomial unit periods in seven variables

Let R be ANY commutative Q-algebra, f(t) in R[t] with f(0)=1, and M>=0
an integer. Put

    B(U)=(1-f(U))/U, N=2M+4,
    C_N(Y)=((1+Y)^N-1)/Y
          =sum_(j=0)^(N-1) binomial(N,j+1)Y^j.

Use the fixed nonsingular split Gaussian covariance UA=ZW=YV=T^2=1,
with all other pairings zero. Define the ORIGINAL polynomial phase

    P_(f,M)=A+B(U)WZ+B(U)V(1+Y)
              -(1/2)B(U)(1+Y)C_N(Y)T^2.                 (1)

Write Phi[Q]=E[Q exp(tP)]. All statements are coefficientwise Wick
identities in R[[t]], with no convergence, domain or reducedness premise.
The normal variables satisfy U=t,Z=0,Y=f^-1-1. Their determinant is f^2.
The remaining T precision is f^-N. Thus, with r=N/2=M+2,

    Phi[1]=f^(r-2)=f^M,
    Phi[(1+Y)^d]=f^(M-d), EVERY d>=0.                  (2)

The half determinant contributes f^r with constant term one; it is not
an unspecified square-root branch. If f=1, phase(1) is simply A.
If deg f=L>=1, deg P<=L+2M+5. No degree-optimality is asserted.

This is the usual formal Gaussian/Lagrange elimination. For published
provenance see [Gessel, Lagrange Inversion, Theorem2.1.1,
equation(2.1.4), and Section4.1](https://arxiv.org/html/1609.05988v1).
Gessel's stated coefficient ring is a Q-DOMAIN. Our all-ring assertion
follows by coefficientwise polynomial identities in universal rational
coefficient variables, then arbitrary base change. Every time coefficient
uses finitely many Wick contractions. No domain theorem is silently
applied to a ring with zero divisors.

## 2. Every original marker and the EXACT degree filtration

Let F_d be the R[t]-span of periods of ORIGINAL spatial polynomial
markers of total degree at most d. Time-polynomial coefficients are
external and unrestricted in this definition. Then

    F_d=f^(M-d)R[t], EVERY d>=0,
    full actual image M_P=R[t,f^-1] inside R[[t]].      (3)

This is a statement about the original-marker O-module filtration.
It is NOT a bound on the time degree of each period, an ordinary last
nonzero-moment formula, or an assertion about a good D-order filtration.
Those distinctions matter for nonsquarefree f as well as for general f.

Here is the full momentum proof. For the monomial

    U^p A^q Z^a W^b Y^c V^d T^(2j),

introduce sources alpha,beta,gamma for A,W,V. Put

    u=t+alpha, D_*=1-tB(u),
    z=beta/D_*, y=(1-D_*+gamma)/D_*,
    R_N(D,gamma)=D^N+(1-D)(1+gamma)
                   *sum_(i=0)^(N-1)(1+gamma)^(N-1-i)D^i.

The exact residual T precision is D_*^-N R_N(D_*,gamma), and
R_N(D,0)=1 as a POLYNOMIAL identity. Source differentiation therefore gives

    Phi[U^p A^q Z^a W^b Y^c V^d T^(2j)]
      =delta_ab a! (2j-1)!!
         * [partial_alpha^q {
              u^p D_*^(M-a-c+Nj)
                * [partial_gamma^d ((1-D_*+gamma)^c
                     R_N(D_*,gamma)^(-j-1/2))]_(gamma=0)
           }]_(alpha=0).                              (4)

Every odd T period is zero. The inner gamma derivative in(4) is a
POLYNOMIAL in D_* because R_N(D,0)=1. The alpha derivatives lower the
displayed D_* exponent by at most q; derivatives of the inner polynomial
and of B(u) introduce polynomial factors only. Hence the result belongs to

    f^(M-a-c-q)R[t] subset f^(M-delta)R[t],
    delta=p+q+2a+c+d+2j,

where a=b; the a!=b case is zero. This proves the inclusion in(3) for
every marker, including all momenta. The reverse inclusion follows from
the ORIGINAL degree-d polynomial (1+Y)^d in(2). Taking the union over d
gives the full image. Differentiating a period multiplies its amplitude
by P, so this image is already a D_R[t]-module. No image kernel is computed.

## 3. Exact all-ring safety, normalization and purity boundaries

Let I be the ideal of the finitely many nonconstant coefficients of f.
From(3), over EVERY R,

    P_(f,M) safe iff f^-1 in R[t]
                 iff f is a unit in R[t]
                 iff I is nilpotent.                   (5)

Indeed a polynomial with constant term one is a polynomial unit precisely
when all its nonconstant coefficients are nilpotent; finitely many such
coefficients generate a nilpotent ideal. This argument works without
Noetherianity. If f is NOT a polynomial unit, the LEAST degree of a bad
original polynomial marker is EXACTLY M+1. All lower periods are
polynomials by(3), and (1+Y)^(M+1) has period f^-1.

The scalar normalization ideal is also exact:

    {c in R : c Phi[Q] polynomial for EVERY original Q}
       =Gamma_I(R)={c : I^L c=0 for some L}.             (6)

For the converse, c/f polynomial implies f^-1 polynomial in
(R/Ann(c))[[t]], so the polynomial-unit criterion there gives I^L c=0.
The forward implication follows from the finite geometric expansion
modulo Ann(c). No radical replacement loses ordinary torsion here.

Purity is a different condition. At M=0 every phase(1) is pure. At M>0
it is pure iff f=1: if f^M=1, the first nonzero positive coefficient
would survive multiplication by the invertible rational integer M.
For example over Q[e]/(e^2), f=1-et,M=1 is safe but nonpure. Over a
field this particular family is safe only when f=1. This does NOT
decide whether every safe phase over a characteristic-zero field is pure.

Any finite prescribed pure-moment packet mu_0=1,mu_1,...,mu_L is realized
exactly with zero tail in seven variables: take M=1 and
f=1+sum_(n=1)^L mu_n t^n/n!. Thus eventual vanishing of pure moments
does not force all-positive purity, even when the entire tail is zero.
These phases generally fail safety at a higher marker, as(5) specifies.
Polynomiality of the unit period alone is not an all-marker safety test.

## 4. Five variables, all-positive purity and arbitrary marker delay

Take f=1-s, so B=1 is constant. Omit U,A entirely. In FIVE variables set

    P_M=WZ+V(1+Y)-(1/2)(1+Y)C_(2M+4)(Y)T^2.           (7)

The source calculation above, with p=q=0, gives

    E[exp(s P_M)]=(1-s)^M,
    F_d(P_M)=(1-s)^(M-d)K[s].                          (8)

Let K be a characteristic-zero field. For any INTEGER k>M, form the
ORIGINAL polynomial phase Q_(M,k)=P_M^k. Then

    E[Q_(M,k)^n]=0, EVERY n>=1,                        (9)
    least bad original marker degree = M+1,            (10)
    E[(1+Y)^(M+1) Q_(M,k)^n]=(kn)!, EVERY n>=0.        (11)

Proof: equation(9) extracts moments of P_M at indices kn>M from the
degree-M polynomial in(8). Every marker of degree at most M has a
polynomial old period by(8), hence its subsequence of moments at indices
kn is eventually zero too. The old period of (1+Y)^(M+1) is(1-s)^-1;
its moment at EVERY old index m is m!, giving(11). Factorials do not
vanish over a characteristic-zero field. This proves minimality in(10)
for the complete degree packet, not just selected low markers.

Taking M=k-1 proves: NO bound depending only on Gaussian dimension
can bound the least bad marker degree, EVEN UNDER ALL-POSITIVE PURITY.
The dimension stays five; deg Q=k(2M+6), so this does not contradict
the retained fixed-DIMENSION-AND-DEGREE uniformity theorems. It also
does not contradict two-linear-marker sufficiency in the older
specific compensator architecture; taking powers changes that architecture.

## 5. IDENTICAL full actual images and seeds have different marker filtrations

For every k>=1 define the actual formal series

    h_k(t)=sum_(n>=0) (kn)!/n! *t^n,
    theta=t partial_t,
    A_k=1-t*k*prod_(j=1)^(k-1)(k theta+j).

Its exact recurrence gives A_k h_k=1. For EVERY 0<=M<k,

    M_(Q_(M,k))=D_K[t] h_k inside K[[t]],
    Phi_(Q_(M,k))[1]=1.                               (12)

To prove BOTH inclusions, retain all old markers using(4). Their old
periods are rational functions with their only possible finite pole at
s=1. Polynomial parts transfer to time polynomials by the linear map

    sum c_m s^m -> sum c_(kn) (kn)!/n! t^n.

Each pole transfers by the EXACT Euler formula

    (1-s)^-p -> [prod_(j=1)^(p-1)(k theta+j)/(p-1)!] h_k.
                                                               (13)

Since 1=A_k h_k, every transferred period belongs to D h_k. Conversely
h_k itself is the period of the original marker in(11), and the actual
image is D-stable. This proves global equality as EMBEDDED modules,
not just equality of generic rank, scalar annihilator or formal type.

Consequently, for fixed k the very SAME full actual module and SAME
unit element have least bad marker degrees 1,2,...,k as M varies.
Even the full actual image with its distinguished unit does not recover
the ORIGINAL marker filtration. Any classification intended to predict
least bad degree must retain the source-to-image marker filtration.
Neither an ambient kernel nor its omission explains this loss of data.

## 6. Exact Airy descent and actual rank

For k>=2, put ell=k-1 and make the finite cover and rational gauge

    t=(-1)^ell/(k x^k), theta_t=-(1/k)theta_x, y=x u.

Then the identity on arbitrary functions u is EXACT:

    A_k[x u]=-(partial_x^ell-x)u.                       (14)

Indeed prod_(j=1)^ell(k theta_t+j)=(-1)^ell prod(theta_x-j),
and prod_(j=1)^ell(theta_x-j)(x u)=x^(ell+1)partial_x^ell u.
This proves(14) with its sign and without any asymptotic substitution.

For ell>=2 the right side is the generalized Airy operator. The
published primary source [Sabbah--Yu, Hodge properties of Airy moments,
Section6.a, Lemma6.2 and the following paragraphs](https://arxiv.org/html/2112.13405v2#S6.SS1)
identifies this operator as a Fourier transform of a rank-one exponential,
proves irreducibility and restates Katz's Galois theorem: SL_ell for odd
ell and Sp_ell for even ell. The source also records infinity slope
(ell+1)/ell. These are classical imports, not computed Galois certificates.
Only this differential-equation material is used; no Hodge or arithmetic
moment theorem is imported. The paper is published in Tunisian Journal
of Mathematics5(2)(2023),215-271, DOI10.2140/tunis.2023.5.215.

It follows that A_k is irreducible over C(t), hence over every
characteristic-zero field: a factor over any such field uses finitely
many coefficients, which embed into C. For ell=1 this is immediate.
Finite algebraic base extension preserves the connected differential
Galois group, so the homogeneous quotient has the displayed SL/Sp
connected group; at k=2 the first-order exponential has G_m instead.
We do NOT specify the entire finite descent group or the full actual
unit-extension Galois group in this packet.

The actual rank of(12) is EXACTLY k. The coefficients (kn)!/n! have
zero complex convergence radius when k>=2, so h_k is not rational.
This nonrationality over arbitrary characteristic-zero fields follows
by the same finite-coefficient embedding argument for any putative
rational expression. In the generic actual module, the class of h_k
modulo K(t)*1 is therefore nonzero and is killed by irreducible A_k.
It generates a rank-ell quotient, giving rank ell+1=k, rather than
merely an order-k upper bound. At k=1, h_1=(1-t)^-1 and the rank is1.

The quotient has positive origin slope 1/ell, obtained by dividing
the Airy infinity slope (ell+1)/ell by the cover degree k=ell+1.
This is NOT a regular-transcendental Gaussian realization. No conclusion
about computed Stokes matrices is drawn.

## 7. Consequences for the still-open general programme

Closed here: all-ring polynomial-seed compensators with exact original
marker O-filtration, full image, safety and scalar normalization; arbitrary
delay of the first bad marker under all-positive purity in dimension five;
and equality of full actual images/unit periods across different marker
filtrations, with exact Airy identification and actual rank.

Still OPEN: exhaustive safe phase normal forms; whether all-marker safety
over a characteristic-zero FIELD forces all-positive purity; realization
of regular transcendental periods in an ORIGINAL Gaussian full image or
a valid general obstruction; and a classification retaining the original
source filtration. Our unsafe families do not settle the safety-to-purity
question. Existing finite-primary cutoff algorithms do not by themselves
provide those structural theorems.

The focused checker compares30 complete original marker formulas against
Wick contractions,12 filtration generators,all36 degree-at-most-two
monomials in seven variables,three pure power cases,three coefficient
transfers/recurrences/Airy gauges,one nilpotent and one arbitrary seed
control,and six rejection cases. These finite calculations transcribe
the proofs; they do not prove their universal quantifiers or compute
Galois groups.
