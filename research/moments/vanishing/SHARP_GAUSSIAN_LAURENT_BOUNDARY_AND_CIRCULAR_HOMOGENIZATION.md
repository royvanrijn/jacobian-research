# The sharp Gaussian--Laurent boundary and circular homogenization

Let K be any characteristic-zero field. All moment functionals below are defined algebraically by their rational monomial values. Thus a real probability realization is useful for terminology, but is not required for the field statements.

Write gamma(X^{2j})=(2j-1)!!, gamma(X^{2j+1})=0, gamma(1)=1. Write C(Z^aW^b)=a! delta_ab for a circular pair, and use tensor products for independent coordinates. For a functional lambda, a phase f is pure if lambda(f^m)=0 for every m>=1, and safe if for every fixed algebra marker h, lambda(hf^m)=0 for all sufficiently large m.

The positive classification and structured Gaussian assertions below are COROLLARIES of Christopher D. Long's arbitrary-torus factorial theorem. The negative direction has a direct rational central-binomial circuit proof. These facts are recorded with their exact structural and field hypotheses; no priority or general Jacobian-conjecture claim is made.

## 1. The sharp boundary when at least one Laurent variable is present

For n>=0 and r>=1 define
A_{n,r}=K[x_1,...,x_n,t_1^{+-1},...,t_r^{+-1}]
and Lambda_{n,r}=gamma_n tensor CT_t, where gamma_n is the standard n-dimensional Gaussian moment functional.

**Theorem 1 (sharp Gaussian--Laurent boundary).** Over EVERY characteristic-zero field K, ker Lambda_{n,r} is Mathieu--Zhao if and only if n<=2.

For n>=3 there is already a rational-coefficient phase of Gaussian degree TWO, with T=1+t_1 and t=t_1,

\[
 P=\frac{T}{t}\,[x_1^2+x_2^2-t(t+2)x_3^2],\qquad H=t,
 \tag{A}
\]

such that for every m>=1,

\[
 \Lambda_{n,r}(P^m)=0,\qquad
 \Lambda_{n,r}(HP^m)=2^m m!\ne0.
 \tag{B}
\]

The marker is a UNIT. The phase is absolutely nonpower even up to multiplication by a unit. Its t-support contains adjacent exponents, so it is not a Laurent polynomial in t^b after a monomial shift for any integer b>=2.

For n=0, the positive direction is the Duistermaat--van der Kallen Laurent theorem, transferred to arbitrary characteristic zero using the finite coefficient field. For n=1 and n=2 it follows from the exact Long embeddings in Sections 3 and 4. The negative phase is also the exact algebra-intertwining image of Long's known ternary quartic, as explained below. A complete local all-order proof of (B) is included, so no assertion rests on initial moment checks.

### 1.1 Exact central-binomial calculation at every order

Use only x=x_1,y=x_2,z=x_3,t=t_1, and set U=x^2+y^2. The formal Gaussian generating function gives
gamma(U^a)=2^a a!, while gamma(z^{2j})=2^j(1/2)_j. Indeed gamma_1(exp(sx^2))=(1-2s)^{-1/2}; multiplying the independent x,y series gives (1-2s)^{-1}. These are identities of formal series, not analytic integrability claims. Here (b)_j=b(b+1)...(b+j-1) and (b)_0=1.

Put a(t)=t(t+2). Expanding the bracket to power m gives

\[
 \gamma_3(P^m)
  =2^m m!\,T^m t^{-m}
       \sum_{j=0}^{m}\frac{(1/2)_j}{j!}[-a(t)]^j .
 \tag{C}
\]

The term j has Gaussian coefficient
binom(m,j) 2^{m-j}(m-j)! 2^j(1/2)_j
=2^m m!(1/2)_j/j!.

In K[[t]], the tail beyond j=m has order at least m+1, since a(t) has order one. The full binomial series is

\[
 \sum_{j>=0}\frac{(1/2)_j}{j!}[-a(t)]^j
  =(1+a(t))^{-1/2}
  =((1+t)^2)^{-1/2}=T^{-1},
 \tag{D}
\]

using the formal root with constant term one. Thus (C) differs from
2^m m! T^{m-1}t^{-m} by a series of order at least one. For every integer ell>=0,

\[
 \boxed{\quad
 \Lambda_{3,1}(t^\ell P^m)
  =2^m m!\binom{m-1}{m-\ell}
  =2^m m!\binom{m-1}{\ell-1}.
 \quad}
 \tag{E}
\]

The binomial coefficient is zero outside 0<=m-ell<=m-1. This proves (B), and positivity for 1<=ell<=m in the real realization. No extension of (E) to negative ell is asserted: for example Lambda(t^{-1}P)=-3.

Extra Gaussian and Laurent variables are ignored; their functional on 1 is one. Therefore this SAME rational phase and unit marker work for every n>=3 and r>=1 over the native field K.

### 1.2 Exact relation to Long's known counterexample

Long's paired ternary quartic is

\[
 P_3=(1+Z)\bigl(W-(2+Z)z^2/2\bigr),\qquad Q_3=Z,
\]

with circular covariance C(Z^aW^b)=a! delta_ab and independent standard Gaussian z. Define

\[
 Z\longmapsto t,\qquad W\longmapsto U/(2t),\qquad z\longmapsto z.
 \tag{F}
\]

This is an injective algebra map: the independent variables U,z,t recover the monomial exponents in Z,W,z, and U=x^2+y^2 is transcendental over K(z,t). It intertwines EVERY polynomial moment, because

\[
 (\gamma_2\otimes\operatorname{CT}_t)
      \bigl(t^a[U/(2t)]^b\bigr)=b!\,{\bf1}_{a=b}.
\]

Its images are P_3->P/2 and Q_3->t. Hence (B) is also the precise pullback of Long's known marked identity, with the scaling factor 2^m. This is a direct classical-corollary relation, stronger than just a similar cancellation method. The independently displayed proof (C)--(E) establishes the same facts locally over Q.

### 1.3 Actual divisor and adjacent Laurent support

The ideal (t+1) is prime in
Kbar[x_1,...,x_n,t_1^{+-1},...,t_r^{+-1}].
The bracket in (A), restricted to t=-1, is x^2+y^2+z^2, a nonzero polynomial in the residue domain. Its (t+1)-valuation is therefore zero, while T/t has valuation one. Thus ord_{t+1} P=1. Units have valuation zero, so P cannot be a unit times an a-th power for any a>=2, even after algebraic closure.

The coefficients of t-exponents -1 and 0 are U and U-2z^2, both nonzero. A monomial shift preserves their difference one, proving the stated Laurent-support primitivity.

### 1.4 Genuine product instability over formally real fields

If K is formally real, the Gaussian factor alone has only the zero pure phase in every finite dimension. The monic rational Hermite polynomials have generating series exp(sX-s^2/2). The Gaussian exponential identity gives
gamma_n(H_alpha H_beta)=alpha! delta_{alpha,beta}. Since these monic polynomials form a basis, for nonzero f=sum c_alpha H_alpha,

\[
 \gamma_n(f^2)=\sum_\alpha \alpha!\,c_\alpha^2\ne0.
\]

Each positive integer weight is a finite sum of identical squares, and formal reality forbids a nontrivial sum of squares zero. Thus the Gaussian kernel is Mathieu--Zhao. The Laurent constant-term kernel is Mathieu--Zhao over every characteristic-zero field by Duistermaat--van der Kallen. Formula (B) gives genuine failure after taking their product when n>=3.

The condition r>=1 is part of Theorem 1. No arbitrary-field boundary for r=0 is inferred; real-coefficient Gaussian positivity persists in every dimension.

## 2. The exact classical input and characteristic-zero transfer

For s>=1 let
B_s=K[u,z_1^{+-1},...,z_s^{+-1}] and define

\[
 \Gamma_s(u^a z^\nu)=a!\,{\bf1}_{\nu=0}.
 \tag{1}
\]

For F in B_s let Omega(F) be its angular exponent support, collecting all nonzero radial coefficients. Use Omega(0)=empty.

**Classical input (Long).** The following are equivalent: all positive pure Gamma_s moments vanish; all sufficiently large pure moments vanish; zero is outside conv Omega(F); and every fixed marker has eventually vanishing mixed moments. If zero belongs to conv Omega(F), there is b>=1 for which Gamma_s(F^{bp}) is nonzero for EVERY sufficiently large rational prime p.

The pinned primary states this over C. For each phase, embed only the field generated over Q by its finitely many coefficients into C. The embedding preserves support and all exact moments. It therefore transfers the angular criterion and the prime-dilation statement to K. Strict separation then proves the marker statement for every K-marker directly. No embedding of all K is used.

For completeness, if ell strictly separates Omega(F) from zero, set delta=min ell(Omega(F))>0. For a nonzero fixed marker G put B=-min ell(Omega(G)). Every exponent of GF^m has ell-value at least m delta-B, so its angular constant term vanishes when m delta>B. For zero F or zero G the conclusions are immediate.

The exact existing repository input is the one-radial factorial theorem with UNRESTRICTED angular rank, together with its prime-index nonvanishing. Increasing angular rank below does not increase the number of radial factorial variables.

## 3. One real Gaussian coordinate and an arbitrary Laurent torus

On A_r=K[X,t_1^{+-1},...,t_r^{+-1}], r>=0, put
Lambda_r=gamma_X tensor CT_t. If F=sum_mu F_mu(X)t^mu is nonzero, its original angular support is the set of mu with F_mu!=0. When r=0 this is the singleton zero vector.

Define an algebra map into B_{r+1} by

\[
 \iota(X)=z_0+\frac{u}{2}z_0^{-1},\qquad
 \iota(t_i)=z_i .
 \tag{2}
\]

It is injective: the highest X-degree a_mu of each nonzero F_mu produces a uniquely largest z_0-exponent a_mu, with its nonzero leading coefficient. The binomial theorem gives

\[
 \Gamma_{r+1}(\iota(X^{2j}t^\mu))
  ={\bf1}_{\mu=0}\binom{2j}{j}2^{-j}j!
  ={\bf1}_{\mu=0}(2j-1)!!,\qquad
 \Gamma_{r+1}(\iota(X^{2j+1}t^\mu))=0.
 \tag{3}
\]

Thus Gamma composed with iota is exactly Lambda_r, simultaneously on every power and every marker.

For each mu, every new angular exponent has first coordinate between -a_mu and a_mu. Both endpoints (a_mu,mu) and (-a_mu,mu) occur without cancellation. At the positive endpoint the unique contribution has X-degree a_mu and u-degree zero; at the negative endpoint the unique contribution has that X-degree and u-degree a_mu. Smaller X-degrees cannot reach either endpoint. Therefore

\[
 \operatorname{Newt}_{ang}(\iota F)
 =\operatorname{conv}\{(a_\mu,\mu),(-a_\mu,\mu):F_\mu\ne0\}.
 \tag{4}
\]

Projection and the symmetric vertical midpoints (0,mu) prove

\[
 0\in\operatorname{Newt}_{ang}(\iota F)
 \ \Longleftrightarrow\
 0\in\operatorname{conv}\{\mu:F_\mu\ne0\}.
 \tag{5}
\]

**Theorem 2 (one Gaussian and Laurent tori).** For F in A_r, all positive pure vanishing, eventual pure vanishing, safety against every fixed A_r-marker, and omission of zero from the ORIGINAL angular convex hull are equivalent. If the angular hull contains zero, one b>=1 gives Lambda_r(F^{bp})!=0 for every sufficiently large prime p. In particular ker Lambda_r is Mathieu--Zhao. For r=1 the pure phases have angular support of one strict sign. For r=0 no nonzero phase is pure.

**Proof.** Apply the classical input to iota F and use (3)--(5). The algebra-homomorphism identity iota(GF^m)=iota(G)iota(F)^m handles every fixed marker. Conversely safety implies eventual pure vanishing by taking G=1. The angular separation bound of Section 2 gives the other easy directions directly. The zero phase is included by empty support.

## 4. Two Gaussian coordinates and Laurent tori

On C_r=K[Z_0,W_0,t_1^{+-1},...,t_r^{+-1}], r>=0, put
Psi_r=C_{Z_0,W_0} tensor CT_t. Define

\[
 \jmath(Z_0)=z_0,\qquad
 \jmath(W_0)=u z_0^{-1},\qquad
 \jmath(t_i)=z_i .
 \tag{6}
\]

The monomial Z_0^a W_0^b t^mu maps to
u^b z_0^{a-b} z^mu. The data (b,a-b,mu) recover the original monomial, proving injectivity. Also

\[
 \Gamma_{r+1}(\jmath(Z_0^aW_0^bt^\mu))
   =b!\,{\bf1}_{a=b,\mu=0}
   =\Psi_r(Z_0^aW_0^bt^\mu).
 \tag{7}
\]

No angular-support cancellation occurs: terms with the same (a-b,mu) have distinct powers of u unless they were the same monomial. Hence the expanded angular support is exactly

\[
 \Omega_{\rm rot,tor}(F)
  =\{(a-b,\mu):[Z_0^aW_0^bt^\mu]F\ne0\}.
 \tag{8}
\]

**Theorem 3 (a circular pair and Laurent tori).** In C_r, positive pure vanishing, eventual pure vanishing, safety against every fixed polynomial/Laurent marker, and zero outside conv Omega_rot,tor(F) are equivalent. If zero is inside, the same fixed-power prime-dilation nonvanishing conclusion holds.

**Proof.** Apply Long through (6)--(8), exactly as for Theorem 2.

The angular condition here is the EXPANDED rotational-plus-Laurent condition (8). It cannot be replaced by the original t-exponents alone: F=Z_0 has only t-exponent zero but is pure and safe.

To realize an actual standard Gaussian pair Y_1,Y_2, one may work over K(i) and set
Z_0=Y_1+iY_2, W_0=(Y_1-iY_2)/2. These are invertible linear coordinates with covariance C. The usual symmetric normalization (Y_1+-iY_2)/sqrt(2) is equivalent after adjoining sqrt(2). The paired functional itself and formulas (6)--(8) are defined over K without either extension. A statement about standard diagonal Gaussian coordinates descends from K(i), as detailed in Section 8.

## 5. Exact circular homogenization with arbitrary coefficient markers

Let A be any commutative K-algebra carrying a K-linear functional lambda. If P in A[t,t^{-1}] has support in [-L,L], L>=0, define

\[
 {\cal H}_L(P)=(ZW)^L P(Z/W)
  =\sum_{j=-L}^{L}P_j Z^{L+j}W^{L-j}.
 \tag{9}
\]

This is a genuine polynomial homogeneous of total degree 2L in the PAIR (Z,W). It is not asserted to be separately homogeneous in Z and W; the coefficients may have arbitrary degrees in A.

For h=sum_{|j|<=e}h_j t^j put
Q_e(h)=(ZW)^e h(Z/W). For EVERY m>=0,

\[
 (\lambda\otimes C)\bigl(Q_e(h){\cal H}_L(P)^m\bigr)
  =(Lm+e)!\,(\lambda\otimes\operatorname{CT}_t)(hP^m).
 \tag{10}
\]

Indeed a product with angular exponents j from h and k from P^m has circular exponents Lm+e+j+k and Lm+e-j-k. Balance is exactly j+k=0, and the factorial is always (Lm+e)!. In particular the marker is h(t), not h(t^{-1}).

Conversely every polynomial of circular total degree 2L is uniquely (9). Every arbitrary polynomial marker splits into finitely many circular total degrees. Odd-degree components have zero expectation against every power of H_L(P); an even component of degree 2e is uniquely Q_e(h_e). Since factorials are nonzero units in K, (10) proves:

* given-phase purity is an iff, with identical zero sets of pure orders;
* safety against ALL source markers is an iff with safety against ALL target markers;
* unsafety is an iff: infinitely many detections by a finite sum force infinitely many detections by at least one fixed even component;
* every-order nonzero marked values transfer for each matching pair h,Q_e(h).

For the last assertion, no reverse existence claim is made for an arbitrary inhomogeneous target marker detecting every order. A finite collection can cover all orders without an individually chosen component doing so. The safety equivalence needs only a common maximum of finitely many eventual thresholds.

## 6. Several fixed circular degrees

Apply (9)--(10) in each constrained pair (Z_i,W_i). If the phase H has circular degree 2L_i in each such pair, write

\[
 H=\left(\prod_i(Z_iW_i)^{L_i}\right)
       P(\text{base coordinates},Z_1/W_1,\ldots,Z_s/W_s).
 \tag{11}
\]

For a source marker with angular support coordinatewise bounded by e_i, its lift gives exact marked multiplier

\[
 \prod_i(L_i m+e_i)! .
 \tag{12}
\]

Every target marker splits into finitely many multidegrees in the constrained pairs. A component odd in at least one pair contributes zero; the all-even components recover precisely the source markers. This proves all-marker equivalence in both directions.

There are two base choices:

(a) ONE unpaired standard Gaussian X, using Theorem 2. Then H lies in K[X,Z_1,W_1,...,Z_s,W_s] and has fixed degree in EVERY circular pair.

(b) ONE unrestricted circular pair (Z_0,W_0), using Theorem 3. Then H lies in K[Z_0,W_0,...,Z_s,W_s] and has fixed degree only in pairs 1,...,s; it may have arbitrary polynomial degrees in the remaining pair 0.

For even constrained degrees, the original source angular weights become the target rotational weights by multiplying the constrained coordinates by 2; the unrestricted pair's rotational coordinate is unchanged. This invertible linear rescaling preserves whether the origin lies in the convex hull.

If one or more constrained degrees D_i are odd, odd powers of H have zero pure expectation. Thus H is pure iff H squared is pure. Safety of ANY phase H is equivalent to safety of H squared: even orders use H squared; odd orders use the fixed marker QH. The squared phase has even degree in every constrained pair and fits (11). Its angular Newton polytope is twice that of H: an exposed vertex's coefficient survives squaring in the relevant polynomial coefficient domain. Eventual purity also passes to the square and then implies the full criterion.

**Theorem 4 (structured Gaussian classification).** In either base case (a) or (b), assume the phase has a fixed total degree D_i>=0 in EACH constrained circular pair. Positive purity, eventual purity, safety, and omission of zero from the convex hull of its FULL rotational-weight vectors are equivalent. In case (a) X has weight zero. In case (b) the unrestricted pair's rotational weight MUST be included. If zero belongs to the hull, one fixed power has nonzero pure moments at every sufficiently large prime dilation.

When an odd constrained degree requires the square, the fixed power in the last statement may be chosen even. When all constrained degrees are even, (12) transfers prime nonvanishing directly. Degree-zero constrained pairs cause no difficulty.

This covers a structured class in every odd dimension 2s+1, and a class in every even dimension 2s+2 with all but one circular pair constrained. It is not a theorem for arbitrary Gaussian phases in those dimensions.

## 7. Explicit marker cutoff and the orthogonal nullcone

If zero is outside the rotational Newton polytope, choose a rational strict separating form and clear denominators to an integer vector c. Let

\[
 \delta=\min_{\omega\in\operatorname{supp}_{rot}(H)}
                  \langle c,\omega\rangle>0,\qquad
 B=-\min_{\nu\in\operatorname{supp}_{rot}(Q)}
                  \langle c,\nu\rangle .
\]

Then every fixed marker Q has zero marked moment whenever m delta>B. This argument holds directly in the polynomial ring and does not rely on a large-order bound derived from factorial sizes.

The cocharacter Z_i->a^{c_i}Z_i, W_i->a^{-c_i}W_i preserves the paired covariance form (and fixes the unpaired X, when present). Every term of H acquires a positive power of a, so H tends algebraically to zero. Thus these pure phases lie in the corresponding split orthogonal nullcone after passing to circular coordinates. Over a nonsplit original field this cocharacter statement is made over a field where those circular coordinates exist.

The converse geometric implication also holds. Fix an ordinary degree bound D containing H, and let V_D be the finite-dimensional vector space of polynomial phases of degree at most D. For every m>=1, I_m(F)=E(F^m) is a polynomial in the coefficients of F and is invariant under the orthogonal group of the Gaussian covariance form. This follows locally from the formal identity E exp(v dot X)=exp(q_cov(v)/2): a covariance-preserving linear substitution leaves every coefficient, hence every polynomial moment, unchanged. The same argument covers the split paired coordinates. Moreover I_m(0)=0.

If zero is in the geometric orthogonal orbit closure of H in V_D over an algebraic closure, invariance makes I_m constant on that orbit, and polynomial continuity gives I_m(H)=I_m(0)=0 for every m. Combined with the explicit cocharacter above, this proves WITHIN the structured classes of Theorem 4 the exact equivalence: purity iff eventual purity iff safety iff rotational Newton separation iff membership in the geometric orthogonal nullcone. The cocharacter even lies in the special orthogonal subgroup. Only the nullcone-implies-purity direction was used without structural hypotheses; no converse for arbitrary Gaussian phases is asserted. This geometric statement does not assert a cocharacter over a nonsplit original field.

For a single pair one can choose c=1 or -1. The negative rotational weight of Q is at most its circular degree, so m delta>deg_circ Q is a sufficient cutoff. In higher angular rank the displayed chosen-c bound remains exact; no uniform c-independent degree estimate is asserted.

## 8. Scalar extension and the limits of the construction

For a GIVEN phase, purity and safety are preserved and reflected under every field extension L/K for these rational moment functionals. Reflection uses the original markers. For preservation, every L-marker is a finite L-linear combination of original coordinate monomials, so finitely many original eventual thresholds have a maximum. Pure moments are exact rational-coefficient linear expressions in the phase coefficients.

Thus working over K(i) to form standard Gaussian circular coordinates does not require a global embedding K->C and does not silently descend an arbitrary new upstairs phase. The field statement in Theorems 2--4 concerns each specified phase and the native scalar-compatible functional.

The structural hypotheses distinguish this result from known Gaussian failures. Long's ternary quartic
(1+Z)(W-(2+Z)X^2/2) has mixed circular degrees and is outside case (a). Long's four-variable cubic with two circular pairs also has mixed degree in BOTH pairs, so is outside case (b). No general ternary or higher-dimensional Gaussian conjecture, no general Mathieu classification and no Jacobian-conjecture conclusion follows.

There is also a useful complex-coefficient guard. For nonzero f in K[X], the radial ansatz
p=(1+t)t^{-1}f((1+t)X)^2 has angular exponents on both sides of or at zero, hence is not pure by Theorem 2. For real nonzero f, the separate Gaussian obstruction proves eventual positive moments, a stronger conclusion. Its modified-family claim cannot be extended indiscriminately to complex constants: if f is constant with f squared=-1, then p=-(1+t^{-1}), and p^a(1+p) is a polynomial of STRICTLY NEGATIVE angular support, hence pure and safe. This agrees with, rather than contradicts, the exact angular criterion.

## 9. Classical sources and proof scope

Christopher D. Long, A Factorially Weighted Constant-Term Theorem on Algebraic Tori and the Gaussian Moments Conjecture in Dimension Two, pinned source commit
2f1ba8ae5d3d282141aca2285326a2d31ce20cd2:
https://raw.githubusercontent.com/octonion/mathematics/2f1ba8ae5d3d282141aca2285326a2d31ce20cd2/gmc/gmc2_stronger_arbitrary_torus.tex .
The imported statements are thm:main and thm:prime-nonvanishing, with the established binary embedding in prop:intertwining. The classical proof uses balanced-face extraction, cofinite good reduction and normalized factorial isolation.

The repository already states these results and the binary Gaussian embedding in
[Laurent and binary Gaussian vanishing equivalences](LAURENT_AND_BINARY_GAUSSIAN_VANISHING_EQUIVALENCES.md).
The linked exposition retains the original Duistermaat--van der Kallen and cofinite-good-reduction inputs; the present proof uses Long's theorem without re-proving those foundational inputs.

Formula (6) is that known binary embedding tensored with more angular variables. Formula (2) is its restriction to one standard Gaussian linear coordinate, and the symmetric vertical-segment argument supplies the stronger original-angular criterion. Formula (10) and the full marker decomposition explicitly record the fixed-circular-degree transfer.

The negative antecedent is Christopher D. Long, Small Counterexamples to the Gaussian Moments Conjecture, arXiv:2607.18186v1, Section 5, Theorem 5.1 and its complete proof:
https://arxiv.org/html/2607.18186v1#S5 .
Its exact five-term paired quartic and the identity E[A(Z)P_3^m]=m![z^m]A(z)(1+z)^{m-1} are the only source scopes used in the comparison (F). Section 1 supplies its own direct normalized proof.

The conclusions are the stated algebraic embeddings, native rational counterpair and structured Gaussian corollaries of the cited results. The fixed-degree hypotheses and the condition of positive Laurent rank remain part of their respective statements.
