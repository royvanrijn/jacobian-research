# Stepanov auxiliary determinants for two-support factorial moments

Research investigation, 21 November 2026. Inspected base revision:
`19eb9912fb60142f8b620dac24896931d4ff3e46`.

## Outcome and assurance

**Continuation update.** The subsequent investigation has closed all four
formerly remaining lower-three gaps, so **every ordinary phase
`X^3(A+BX^d)`, d>=1, now satisfies all-order SFC**. It also proves the infinite
support family `{3^r,3^(r+1)}`, closes lower-four gaps eight and eleven, and
establishes new multivariate transfer families. The complete current synthesis,
proof map, and remaining scope are in
[STEPANOV_SFC_PROGRESS_AND_LOWER_THREE_CLOSURE.md](STEPANOV_SFC_PROGRESS_AND_LOWER_THREE_CLOSURE.md).
The earlier reductions below remain valid; their formerly surviving classes
are explicitly identified as intermediate stages.

This investigation supplies a complete finite-certificate extension of the
ordinary univariate two-support theorem to **every terminal order through 12**,
an algebraic marked-moment recovery theorem valid even over commutative
Q-algebras, and precise obstructions to several proposed determinant
amplifications. A new valuation step in the retained native-prime machinery
improves its ordinary residual congruences from modulus p^2 to p^3.
An additional harmonic-number extraction gives explicit all-order families;
that extraction, unlike the extra valuation step, is already a consequence
of the retained pointwise theorem.

No unrestricted all-order SFC(2) theorem, explicit global eventual cutoff, or
general multivariate unmarked theorem is proved here. The algebraic determinant
arguments assume no coefficient positivity, ordered coefficient field, or
Hermitian inner product. No Alpöge or Riemann-hypothesis claim is a premise.

The new finite certificate contains every modular Bézout witness, not just a
digest or a sampled prefix. Two separately written implementations checked all
35,061 new cases. The written universal algebraic arguments have internal
cross-review and exact controls. These are not formal verification or external
peer review. The repository-wide navigation check fails on an existing stale
supporting-checker pin. The global registry and generated status pages are
unchanged; the new results are recorded in these research notes only.
The precise replay boundary is recorded in Section 9.

## 1. The actual problem and current proof boundary

For a characteristic-zero field K, extend the factorial functional K-linearly:

\[
\mathcal L(X^\alpha)=\alpha!:=\prod_i\alpha_i!.
\]

Here SFC(2) means **two monomials**, not two variables. Its ordinary univariate
case is

\[
f=X^a(A+BX^d),\qquad a\ge0,\quad d\ge1,
\]

with integer exponents. For terminal N >= 2, it asks whether the two moments
of orders N-1,N can both be zero for nonzero f. Define

\[
H_{m;a,d}(x)=\sum_{j=0}^{m}\binom mj(am+1)_{dj}x^{m-j}.
\tag{1}
\]

For B nonzero, the exact normalization is

\[
\mathcal L(f^m)=B^m(am)!H_{m;a,d}(A/B).
\tag{2}
\]

These are monic integer polynomials with nonzero constant terms. The lower
and upper polynomials use the **different** shapes a(N-1)+1 and aN+1. A
fixed-shape Appell discriminant cannot replace this actual pair when a > 0.

For fixed a,d,N, coprimality over Q is equivalent to no common complex root.
A rational Bézout identity then persists over every characteristic-zero
field. An embedding of that entire field into C is neither needed nor
assumed. If A or B is zero, the phase is a monomial and its positive moments
are nonzero. Conversely a nonconstant rational gcd supplies a complex
counterexample via (2).

### 1.1 Latest recorded theorem, with its quantifiers

The [support-uniform synthesis](SUPPORT_UNIFORM_EVENTUAL_ORDINARY_BINOMIAL_SYNTHESIS.md)
records

\[
\exists N_*\ \forall a\ge0\ \forall d\ge1\ \forall N\ge N_*:
\gcd(H_{N;a,d},H_{N-1;a,d})=1.
\tag{3}
\]

N_* is unevaluated. Together with fixed-terminal support finiteness this
implies a finite set of bad triples (a,d,N), **not that the set is empty**.
It is stronger than the earlier fixed-support quantifier order
`for every (a,d), there exists N_(a,d)`.

The elementary synthesis was reconstructed. Put c=a/(d-1)^2, s=(d-1)/2,
and r=(a-s^2)/(a+s), for a >= 1,d >= 2. The identity

\[
r=\frac{c-1/4}{c+1/[2(d-1)]}
\]

allows a finite covering by a small-c sector, a compact subcritical interval,
a fixed critical strip, a compact supercritical interval, and a large-c
sector. The large-thermal sequential exclusion gives fixed delta,U,N_l by
a diagonal contradiction; the bounded-thermal theorem is then applied at
that **one fixed U**. With h=min(1/8,delta/16), |c-1/4| <= h implies
|r| < delta. At c >= 4*10^24, the existing a >= 10^24 d^2 theorem applies
because (d-1)^2 >= d^2/4. Taking a maximum of the finitely many regional
cutoffs, with the a=0 and d=1 boundaries, gives (3).

No quantifier error was found in that covering. The complete contour estimates
in its regional inputs remain essential mathematical premises; this pass
does not supply a second full derivation of that transitive analytic proof.
Moreover the current global and two regional source guards do not replay:
see Section 9. Thus (3) is not used to prove any new result in Sections 2–7.

### 1.2 Reproduced prefix and precise remaining scope

The [terminal-through-ten theorem](ORDINARY_BINOMIAL_TERMINALS_THROUGH_TEN.md)
was reproduced in full: 24,159 supports, 24,160 full modular Bézout identities,
and a byte-identical receipt. Its sole first-prime exception
(N,a,d)=(9,54,8) succeeds at 1000003 after failing at 65537.
The small-ratio regional checker also reproduced its recorded receipt.

The recorded all-order ordinary families include a=0,1,2; d=1; degree at
most seven; d >= 3a for a >= 2; and a >= 10^24 d^2 for d >= 2.
Their large finite computations were not all repeated in this investigation.
They are explicitly retained proof inputs when used to narrow future targets.

Let

\[
e_N=\min\{e:2^e\ge4N(N-1)\},\qquad D_N=2(N-1)e_N.
\]

After Section 2, a convenient **containing** set of unresolved ordinary
triples, using those recorded families, is

\[
\begin{gathered}
N\ge13,\quad a\ge3,\quad d\ge2,\quad a+d\ge8,\\
d<\min(D_N,3a),\qquad
a<\min\{d(N-2),10^{24}d^2\}.
\end{gathered}
\tag{4}
\]

Special arithmetic results remove further points. Equation (4) is not a
minimal list of open triples. Accepting all inputs of (3) also gives N < N_*,
but without an evaluated N_* this is not an executable finite global box.

For multivariate two-support phases the moment is instead

\[
\mathcal L((AX^\alpha+BX^\beta)^m)
=\sum_{j=0}^{m}\binom mj A^{m-j}B^j
  \prod_i((m-j)\alpha_i+j\beta_i)!.
\tag{5}
\]

The ordinary one-variable theorem does not cover (5) in general.

## 2. A proved extension through terminal order twelve

**Theorem T12.** For every characteristic-zero field K, every nonzero ordinary
univariate binomial f=X^a(A+BX^d), and every integer 2 <= N <= 12,

\[
\bigl(\mathcal L(f^{N-1}),\mathcal L(f^N)\bigr)\ne(0,0).
\tag{6}
\]

This is a bounded-terminal theorem with **unbounded exponents and arbitrary
field coefficients**. It is not pairwise coprimality of every two moments
among the first twelve.

### 2.1 Universal complements

Use exactly these existing proof inputs:

1. `thm:rational-annulus-cutoff` and `lem:monotone-centers` in
   [lower_exponent_initial_and_boundary.tex](lower_exponent_initial_and_boundary.tex).
   For C_(N-1),C_N > 0 satisfying (1+1/C_j)^j < 2, and
   B_j=((d_0+1)j)!/j!, the inequality

   \[
   B_N^{N-1}>(C_{N-1}C_N)^{N(N-1)}B_{N-1}^{N}
   \tag{7}
   \]

   excludes a common root for every a >= 1,d >= d_0. Endpoint-normalized
   factorial log-convexity puts the roots in rational annuli. The ratio of
   their centers increases with both a and d, so the verified endpoint
   a=1,d=d_0 suffices for the unbounded region.

2. `ImprovedLinearLowerExponentCriterion` in
   [GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md](GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md):
   a >= d(N-2) excludes the actual consecutive pair, **including equality**.
   The upper-order inner annulus is strict. This input uses the written
   Gamma-factor upper disk proof and closed-disc Grace composition with
   reversal; both dependencies are pinned by the new replay.

3. `thm:large` in
   [ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex](ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex):
   when a=0 and 2^d >= 4N^2, the degree-N constant-family polynomial is
   squarefree. Only at a=0 one has H_N'=N H_(N-1), giving coprimality.

All numerical hypotheses of these complements are checked with integers:

| N | C_(N-1) | C_N | d_0 | First automatic a=0 gap |
|---|---:|---:|---:|---:|
| 11 | 14 | 16 | 57 | 9 |
| 12 | 16 | 17 | 65 | 10 |

For example 2^9 >= 4*11^2, but 2^9 < 4*12^2, so the constant support
a=0,d=9 must be checked at N=12. It is included.
The inequalities (7) have positive integer differences of bit lengths
50,047 and 71,073, respectively; both implementations reproduce their exact
digests and the strict radius inequalities.

These retained complex root bounds concern positive numerical Gamma shape
parameters, which here are positive integers. They do not impose positivity
on A,B or an ordering on K. They are separate from the algebraic marked-rank
arguments below.

### 2.2 Complete finite certificates

After applying only these three complements, the remaining domains are:

| N | Positive-a domain | Positive cases | Constant cases | Total |
|---|---|---:|---:|---:|
| 11 | 1 <= d <= 56; 1 <= a < 9d | 14,308 | a=0, d=1,...,8: 8 | 14,316 |
| 12 | 1 <= d <= 64; 1 <= a < 10d | 20,736 | a=0, d=1,...,9: 9 | 20,745 |

For each of these **35,061** triples the retained certificate records p and
all ascending coefficients of U,V satisfying

\[
U(x)H_{N;a,d}(x)+V(x)H_{N-1;a,d}(x)=1\pmod p,
\quad \deg U<N-1,\quad\deg V<N.
\tag{8}
\]

Every case succeeds at p=65537. The largest factorial argument is
8436 < p, so all normalizing factorials are units; both actual degrees and
their monic leading coefficients are preserved. Thus (8) proves that the
integer Sylvester resultant is nonzero modulo p and hence nonzero over Q.

The compressed complete witness file is
[stepanov_terminal_11_12_bezout.jsonl.gz](stepanov_terminal_11_12_bezout.jsonl.gz).
Its SHA256 is
`401fe2b233dd45dff604aa955b0b6134662ccfd894675b2155375f082d699ae1`.
It contains actual identities rather than just successful-gcd assertions.

Generation used Python-FLINT 0.9.0. The first checker reconstructs coefficients
from factorial tables and verifies products with standard-library arithmetic.
The separate verifier imports neither that generator nor FLINT nor SymPy:
it enumerates the domain in a different order, builds all 70,122 coefficient
vectors from binomial recurrences and rising-factorial blocks at each own-order
shape, and checks all 777,771 coefficients of (8) by dot products. It also
checks 17 Appell boundary identities and rejects 15 malformed, altered, or
incorrect-domain controls, including the wrong use of the upper Gamma shape
for the lower polynomial. All five written complement inputs are strictly
source-pinned.

Together the complements and complete finite certificate prove the N=11,12
cases. The reproduced retained theorem supplies N <= 10. Finally (2) and
rational Bézout transfer prove (6) over every characteristic-zero field.
This proof does not use the support-uniform eventual synthesis. QED.

## 3. A universally invertible marked-moment matrix

The following argument works in any number of variables, even when the two
exponent vectors are incomparable.

**Theorem M.** Let R be a commutative Q-algebra, let
f=AX^alpha+BX^beta with distinct nonnegative integer exponent vectors, fix
n >= 1, and choose i with delta=beta_i-alpha_i nonzero. Set

\[
\gamma_k=(n-k)\alpha+k\beta,\quad e_k=\gamma_{k,i},\qquad
V_{rk}=\mathcal L(X_i^rX^{\gamma_k}),\quad0\le r,k\le n.
\]

Then

\[
\boxed{\det V=
 \left(\prod_{k=0}^n\gamma_k!\right)
 \delta^{n(n+1)/2}\prod_{k=1}^n k!.}
\tag{9}
\]

**Proof.** The entry is gamma_k! (e_k+1)_r. Factor the factorial from each
column. Changing the monic rising-factorial basis (z+1)_r to ordinary powers
has determinant one. The remaining determinant is the Vandermonde
product over j<k of e_k-e_j. Since e_k-e_j=(k-j)delta, (9) follows,
including its sign for negative delta. Every factor is a nonzero rational
number and hence a unit in R. No positivity is used. QED.

### 3.1 Explicit inverse and coefficient-ring consequence

Expand the Lagrange interpolation polynomial in the rising basis:

\[
p_k(z)=\prod_{\ell\ne k}\frac{z-e_\ell}{e_k-e_\ell}
      =\sum_{r=0}^{n}c_{kr}(z+1)_r,
\qquad Q_k(X_i)=\sum_{r=0}^{n}c_{kr}X_i^r.
\]

All c_kr are rational. Then

\[
\mathcal L(Q_kX^{\gamma_j})=\gamma_j!\,\delta_{kj},
\qquad
\mathcal L(Q_kf^n)=\binom nk A^{n-k}B^k\gamma_k!.
\tag{10}
\]

Consequently

\[
\boxed{\mathcal L(X_i^rf^n)=0\ (0\le r\le n)
\iff f^n=0\iff(A,B)^n=0.}
\tag{11}
\]

The exponents gamma_k are distinct, so (10) recovers every coefficient of
f^n; the mixed monomials generate (A,B)^n. The reverse implications are
immediate. Over a field, (11) forces A=B=0. Over a nonreduced Q-algebra it
allows nonzero nilpotent coefficients, as it should.

For f=AX+BX^2,n=2, explicit duals are

\[
Q_0=10-5X+X^2/2,\quad Q_1=-15+9X-X^2,\quad
Q_2=6-4X+X^2/2,
\]

with contractions 2A^2,12AB,24B^2. The marker degree n is independent of
the support exponents and variable count. The n+1 measurements are minimal
for linear reconstruction on the full (n+1)-dimensional support span;
no minimality for merely detecting the smaller family of binomial powers
is asserted.

### 3.2 Exact operator bridge and the missing implication

Write M_n(A,B)=L(f^n), theta=B partial_B and delta_j=beta_j-alpha_j.
Coefficientwise,

\[
\mathcal L(X^\eta f^n)
=D_{n,\eta}M_n,
\qquad
D_{n,\eta}=\prod_j(n\alpha_j+1+\delta_j\theta)_{\eta_j},
\tag{12}
\]

and

\[
M_{n+1}=A D_{n,\alpha}M_n+B D_{n,\beta}M_n.
\tag{13}
\]

On A^(n-k)B^k the operator multiplies by the exact marked-to-unmarked
factorial ratio, proving (12); multiplication by f proves (13). Coefficients
A,B stand **outside** the operators.

The two scalar zeros M_n=M_(n+1)=0 give one zero marked form and one
specific combination of other marked forms. They do not supply all n+1
zero forms required by (11). Differentiating a polynomial identity that
holds globally is legitimate; differentiating a single zero value as if it
were an identity is not.

In fact direct same-degree ideal membership cannot provide the bridge.
In Q[A,B], put I=(M_n,M_(n+1)). Its degree-n part is Q M_n, whereas (9)
shows that the n+1 marked forms span all binary forms of degree n. They
cannot all belong to I. Any ideal-based amplification must increase
coefficient degree, use a justified localization, or establish an extra
actual-family relation. This is a degree obstruction even when the
particular pair is already coprime; it is not a counterexample to SFC.

The raising relation alone is insufficient too. For alpha=1,beta=2,n=2,
dehomogenize A=1,B=t and set theta=t partial_t. The operator
T=(3+theta)+t(3+theta)(4+theta) sends
p=11-20t+9t^2 to 33+52t-355t^2+270t^3. Both vanish at t=1, while
p'(1)=-2. This is an abstract seed, **not** an actual factorial polynomial.
An argument using only the operator and degree would therefore be invalid.

### 3.3 A general factorial Gram factorization

For a finite downward-closed S in N^d, order monomials by total degree and
put G_(alpha,beta)=(alpha+beta)!. The binomial identity

\[
\frac{(\alpha+\beta)!}{\alpha!\beta!}
=\sum_{\gamma\le\alpha,\beta}
  \binom\alpha\gamma\binom\beta\gamma
\]

gives G=DTT^tD, with D_(alpha,alpha)=alpha! and unit triangular
T_(alpha,gamma)=binom(alpha,gamma). Thus

\[
\det G=\prod_{\alpha\in S}(\alpha!)^2.
\tag{14}
\]

This is an algebraic invertibility statement over any Q-algebra. It does
not turn the factorial functional on complex-coefficient polynomials into
a positive inner product. Theorem M uses many fewer markers when the support
lies on the binomial progression.

## 4. Hankel determinants: exact failures and the scope of nonvanishing

Let mu_j=L(f^j) for f=X+tX^2. Exact computation gives

\[
\det\begin{pmatrix}\mu_0&\mu_1\\\mu_1&\mu_2\end{pmatrix}
=1+8t+20t^2,
\tag{15}
\]

and

\[
\mu_1\mu_3-\mu_2^2
=2(6t+1)(72t^3+60t^2+12t+1).
\tag{16}
\]

Equation (15) vanishes at t=(-2+i)/10 and its conjugate. More strongly,
at the rational binomial f=X-X^2/6 all of mu_1,mu_2,mu_3 equal 2/3,
so (16) is zero although none of its entries is zero.

For the spatial marked Hankel family

\[
\mathsf H_{r,n}(A,B)=\det[
\mathcal L(X^{u+v}(AX^a+BX^b)^n)]_{0\le u,v<r},\qquad a<b,
\]

each entry is homogeneous of coefficient degree n. Its endpoint values are

\[
\mathsf H_{r,n}(A,0)=A^{nr}\prod_{j=0}^{r-1}j!(an+j)!,\quad
\mathsf H_{r,n}(0,B)=B^{nr}\prod_{j=0}^{r-1}j!(bn+j)!.
\tag{17}
\]

To verify the scalar determinant, factor (s+u)! from row u of
[(s+u+v)!]; the remaining columns are monic rising polynomials at nodes
0,...,r-1, giving the product of j!. Thus (17) is algebraic.

It follows that H_(r,n)(1,t) has positive degree nr and nonzero constant
term. It has a nonzero root over an algebraic closure. **No single such
determinant is nonzero at every genuine binomial over every characteristic-zero
field.** This does not rule out several minors together or a nonvanishing
theorem specifically on the hypothetical common-moment locus.

A separate control shows why the support assumption matters. Over Q(i),

\[
g=X^2+(-4+2i)X+2-2i
\]

has L(g)=L(g^2)=0 and L(g^3)=32+80i, but

\[
\begin{pmatrix}L(g)&L(Xg)\\L(Xg)&L(X^2g)\end{pmatrix}
=\begin{pmatrix}0&2i\\2i&4+8i\end{pmatrix},\qquad\det=4.
\tag{18}
\]

This phase has **three** terms. It refutes support-free propagation from
two zero scalar moments to spatial marked rank loss. It is not an SFC(2)
counterexample, and it does not violate the length-three SFC window.

## 5. What Stepanov amplification would have to add

The inspected primary sources are Stepanov's 1969 paper, especially printed
pages 1110–1113, and Bombieri's Bourbaki exposition, printed pages 236–238:

- S. A. Stepanov, *On the number of points of a hyperelliptic curve over a
  finite prime field*, Math. USSR-Izv. 3(5) (1969), 1103–1114,
  [primary record and DOI](https://www.mathnet.ru/php/archive.phtml?jrnid=im&option_lang=eng&paperid=2197&wshow=paper).
- E. Bombieri, *Counting points on curves over finite fields*, Bourbaki
  exposé 430 (1973; published 1974), 234–241,
  [primary PDF](https://www.numdam.org/item/SB_1972-1973__15__234_0.pdf).

The classical construction requires a nonzero auxiliary function, forced
high-order zeros, and a smaller pole/degree budget. Its linear dimension count
is accompanied by a separate injectivity proof. In Bombieri's formulation,
Frobenius supplies a linear space of p-power functions; separated pole orders
prove the required multiplication map injective. Those p-powers provide
multiplicity without independent ordinary-derivative conditions at every
point. Characteristic zero has no corresponding Frobenius-fixed evaluation
identity here. A sum of high powers need not itself have high-order roots:
x^2+1 already has simple roots. An SFC analogue must prove its own relation
among vanishing conditions and its own nonzero determinant at every intended
specialization.

### 5.1 Raw interpolation gives no gain

For distinct nodes alpha_i and multiplicities e_i, let M=sum e_i. The jet map
on polynomials of degree at most D,

\[
H\longmapsto(H^{(j)}(\alpha_i)/j!)_{i,\,j<e_i},
\]

has kernel the multiples of product_i(x-alpha_i)^(e_i). Its rank is
min(D+1,M). A nonzero kernel requires D >= M; a root-count contradiction
would require D < M. A special auxiliary space must supply a real dependence
between conditions without losing the corresponding dimension in its
polynomial image. Counting formal coefficients alone does not do that.

### 5.2 Full shifted jets are exactly the original resultant

For nonzero P,Q of degrees m,n, use columns

\[
P,xP,\ldots,x^{n-1}P,\quad Q,xQ,\ldots,x^{m-1}Q
\]

and rows D_x^r/r!, r=0,...,m+n-1. The determinant is

\[
(-1)^{mn}\operatorname{Res}(P,Q).
\tag{19}
\]

Indeed the normalized Taylor-translation matrix from ascending monomial
coefficients has determinant one. The coefficient matrix is the Sylvester
multiplication matrix with the displayed sign convention. If e=deg gcd(P,Q),
its kernel consists of (U,V)=(TQ/G,-TP/G), deg T<e, and its rank is m+n-e.

For the actual pair P=H_N,Q=H_(N-1), mn=N(N-1) is even: (19) is precisely
the original resultant. Unnormalized derivatives multiply it by the known
nonzero product of r!. This is an excellent finite certificate, as Section 2
illustrates. It supplies no independent universal rank theorem.

### 5.3 Product Wronskians amplify degree and multiplicity equally

Let W denote the ordinary, unnormalized Wronskian. For s >= 1,

\[
\boxed{W(P^s,P^{s-1}Q,\ldots,Q^s)
=\left(\prod_{j=0}^{s}j!\right)
 (PQ'-P'Q)^{s(s+1)/2}.}
\tag{20}
\]

Proof: W(Gh_0,...,Gh_s)=G^(s+1)W(h_0,...,h_s) by the triangular Leibniz
matrix. Set z=Q/P in K(x), factor G=P^s, and change derivatives from x
to z. The diagonal factors are 1,z',...,(z')^s; z'=(PQ'-P'Q)/P^2.
All powers of P cancel. This rational calculation is a polynomial identity,
including dependent cases.

For the actual monic degrees N,N-1, PQ'-P'Q has degree 2N-2 and leading
coefficient -1. Thus (20) is already a nonzero polynomial even if P,Q have
a common root. If G=gcd(P,Q), then G^2 divides PQ'-P'Q, and the amplified
degree comparison gives only deg G <= N-1. The amplification parameter s
cancels completely. This rules out an improved multiplicity/degree bound
from this whole product-Wronskian family.

For a simple abstract countercontrol, P=(x+1)(x+2),Q=x+1 have the nonzero
Wronskian -(x+1)^2 and still share x=-1. Nonzero as a polynomial is different
from nonzero at the alleged common root.

## 6. A precise derivative-closure certificate, with its circularity boundary

For a characteristic-zero field, let I=(P,Q)=(G) and let sigma be a
polynomial. Then

\[
\sigma P',\ \sigma Q'\in I
\iff G\mid\sigma G'
\iff\operatorname{rad}(G)\mid\sigma.
\tag{21}
\]

For the first implication apply sigma D_x to G=UP+VQ. For the reverse,
differentiate P=GP_0 and Q=GQ_0. Characteristic zero makes each irreducible
factor's multiplicity drop by exactly one on differentiation, proving the
second equivalence.

For the actual pair, whose constant terms are nonzero, explicit identities

\[
x^kP'=UP+VQ,\qquad x^kQ'=RP+TQ
\tag{22}
\]

would prove coprimality. Conversely coprimality manufactures (22) by Bézout.
Thus (22) is **a certificate form equivalent to the target**, not an
independent weaker sufficient condition. Solving it by dividing by the
unknown resultant would be circular. It becomes useful only if its
coefficients and harmless denominators are constructed independently.

The actual shape identities explain the obstacle. Put n=N-1,
gamma=an+1, E=x D_x, and

\[
T_\ell(E)=(\gamma+dn-dE)_\ell.
\]

Then

\[
(\gamma)_aP=[xT_a(E)+T_{a+d}(E)]Q,
\qquad
(\gamma)_aP'=(n+1)T_a(E)Q.
\tag{23}
\]

These follow coefficientwise from factorial shifting and Pascal's identity.
In xT_a(E), apply T_a first, then multiply by x. At a common root they are
coupled jet equations, not the simultaneous derivative closure (22). For
a=1, the second equation reads

\[
P'+dxQ'=[N+d(N-1)]Q,
\]

which permits both derivatives to be nonzero at a common root. At a=0 it
becomes P'=NQ, but one still needs an independent squarefreeness argument.

## 7. Native-prime determinants: use the existing mechanism first

The repository already contains a genuine algebraic auxiliary determinant
method in [PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex](PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex),
`thm:criterion`, built on the full actual concomitant of
[LOWER_TWO_PRIME_GAP_SFC.tex](LOWER_TWO_PRIME_GAP_SFC.tex).
Its independent reduced cyclic rank and affine nonzero-slope determinant
are extra structure beyond the generic Wronskians of Section 5.

The pointwise criterion, its harmonic evaluation, a further ordinary
valuation refinement, and their exact certificates are written in
[STEPANOV_NATIVE_PRIME_OFFSET_COROLLARIES.md](STEPANOV_NATIVE_PRIME_OFFSET_COROLLARIES.md).
The following identifies exactly which additional step is obtained here.

### 7.1 An extra ordinary p-adic digit

Let p be prime, 1 <= a < p, n=N-1 >= 1, gamma=an+1, and

\[
C=(a+p)n+p+1,\qquad \kappa_a=\frac{p(n+1)}{(C)_a},
\qquad \sigma_a(x)=[Y^2](x+Y)_a.
\]

The retained theorem excludes a common root if sigma_a(gamma mod p) is
nonzero or v_p(kappa_a) >= 0. The new ordinary boundary argument also
excludes the case sigma_a(gamma mod p)=0 and v_p(kappa_a)=-1.

**Proof of the new step.** A common lower root t gives the reciprocal
u=epsilon/t, epsilon=(-1)^(p+1). The coefficient of t^k in
T_(n,gamma,p) has valuation at least k: each of the k complete blocks of p
consecutive integer factors contains a multiple of p. In the monic
reciprocal, v_p(u)<1 would make its leading term uniquely smallest in
valuation. Hence

\[
v_p(u)\ge1.
\tag{24}
\]

The retained reduced-kernel lemma, applied to the actual primitive rising
jet w, gives a unit last coordinate w_p and

\[
\overline{v\cdot w}=c\,\sigma_a(\bar\gamma),\qquad c\ne0,
\qquad w_p=\epsilon\kappa_a u(v\cdot w).
\tag{25}
\]

Here v is the integral row of rising-basis coefficients of
V_a(E)=((E)_a-(C)_a)/(E-C), with its constant term omitted because the
lower moment is zero. Both equations in (25) use the same primitive jet
normalization; the last equality is the actual upper-moment equation.

If sigma_a(bar gamma)=0, then v_p(v dot w)>0. This valuation may be
fractional, which causes no difficulty. If v_p(kappa_a) >= -1, (24)
makes epsilon kappa_a u integral, so the right side of (25) has strictly
positive valuation, contradicting the unit w_p. This proves the new
exclusion. It uses no positivity principle. QED.

Combining this with the retained nonzero-sigma branch, any common root
requires v_p(kappa_a) <= -2. The valuation is integral because kappa_a is
rational. Exact cancellation gives

\[
\kappa_a=\frac{p}{(a+p)\prod_{j=0}^{a-2}(C+j)}.
\tag{26}
\]

Since a<p, at most one of the remaining denominator factors is divisible
by p. Thus survival requires p^3 | C+j for some j, not just p^2 | C+j.
At that residue bar gamma=-j. Writing H_m for the harmonic number only
in the following formula, the exact integer Hasse value is

\[
s_{a,j}=\sigma_a(-j)=(-1)^j j!(a-1-j)!
 (H_{a-1-j}-H_j).
\tag{27}
\]

Consequently the new necessary condition for a bad terminal is

\[
\boxed{(a+p)N\equiv a-1-j\pmod{p^3},\quad
0\le j\le a-2,\quad p\mid s_{a,j}.}
\tag{28}
\]

The companion note proves (27) by differentiating the remaining factors
of (x)_a at x=-j, and gives exact residue certificates. The same new
valuation step works for positive rational p-integral Gamma shape, since
each p-block still contains a factor in p Z_p. It must not be inferred
from integrality in an arbitrary ramified extension. For example, at p=3,
gamma=1+sqrt(3), n=1, the reciprocal root has valuation 1/2 because
(gamma)_3=24+14sqrt(3). This is a scope countercontrol, not an SFC failure.

The retained nonzero-sigma branch uses a proved real-Gamma even-moment
property after showing that a hypothetical common factor is linear.
That dependency is declared in the companion note; it is not an assumption
that a complex Hankel determinant is positive. Every ordinary conclusion
descends to rational coprimality and then extends to arbitrary
characteristic-zero coefficient fields.

### 7.2 Explicit all-order corollaries and remaining targets

For even a, all integers in (27) are nonzero. The retained pointwise
criterion therefore excludes every terminal at every prime p>a not
dividing their product. For a=6 the exact list is
274,-26,4,4,-26,274, with exceptional primes above six only 13,137.
The relative-gap theorem already excludes 137. Thus all prime gaps p>6
except possibly 13 work at every terminal, with gaps 7,11,17 the useful
additional explicit consequences inside p<3a. Lower exponent eight
likewise works at all primes p>8 except possibly 11. These harmonic
corollaries package consequences of the retained theorem; no new mechanism
or literature-priority claim is attached to them. The bare claim of
all sufficiently large prime gaps was already implied by relative-gap.

At the end of the extra-digit stage, the useful first target was lower
exponent a=3. The recorded degree-seven and relative-gap results reduced its
completion to gaps 5,6,7,8. The following table retains those intermediate
domains and records their subsequent closure:

| Support | Gap | Domain left by the extra-digit stage | Current result |
|---|---:|---|---|
| {3,8} | 5 | N = 47 + 125k, k >= 0 | [Closed at every order](STEPANOV_LOWER_THREE_PRIME_CLOSURE.md) |
| {3,9} | 6 | N >= 13 | [Closed, including all adjacent powers of three](STEPANOV_ADJACENT_POWERS_OF_THREE_SFC.md) |
| {3,10} | 7 | N = 103 + 343k, k >= 0 | [Closed at every order](STEPANOV_LOWER_THREE_PRIME_CLOSURE.md) |
| {3,11} | 8 | N >= 13 | [Closed at every order](STEPANOV_LOWER_THREE_EIGHTH_GAP_CLOSURE.md) |

The two prime-gap rows arose from the central Hasse-derivative failure:
sigma_3(-1)=0.
The affine determinant's slope is no longer a unit. The valuation argument
above gains one digit, but cannot be iterated merely from the qualitative
fact v_p(v dot w)>0. The continuation supplies precisely the missing control:
a fixed unit minor produces an integral analytic normalized jet, whose first
two scalar coefficients are proved uniformly. Injectivity excludes deep
valuation classes; a complete Hensel certificate handles the last gap-seven
quadratic branch.

## 8. Failed routes and concrete next work

| Attempt | Exact obstruction or outcome |
|---|---|
| Infer nonzero pure Hankel determinants from a nonzero binomial | Rational counterexample (16); complex counterexample (15) |
| Use one spatial Hankel determinant uniformly | Positive homogeneous degree with nonzero endpoint terms forces a nonzero algebraic zero, (17) |
| Infer all marked zeros from two scalar zeros without using support | Three-term countercontrol (18) |
| Express same-power coordinate marks directly in the two-moment ideal | Its degree-n part has dimension one, whereas the marked span has dimension n+1 |
| Obtain multiplicity solely from the raising identity and degree | Explicit abstract seed p(1)=Tp(1)=0 but p'(1) nonzero |
| Count auxiliary coefficients and use ordinary interpolation | Jet rank min(D+1,M) supplies no degree contradiction |
| Enlarge to the full shifted-jet determinant | It is exactly the Sylvester resultant (19) |
| Take all products P^(s-j)Q^j | Identity (20) gives only deg gcd <= N-1 at every s |
| Solve derivative closure by generic inversion | The missing denominator is the same unknown resultant |
| Derive the kappa-integral prime-gap criterion again | Valid, but already contained in the stronger retained native-prime theorem |
| Turn one bad univariate pair into infinitely many by new-variable monomials | Preserves zeros but leaves the univariate theorem's domain |

The last distinction matters even if (3) is accepted. Multiplying f(X) by
Y^r in a new variable multiplies its m-th moment by (rm)!, but creates a
multivariate phase. Multiplying by X^r or replacing X by X^q changes the
factorial weights and does not preserve scalar moment zeros in general.
For example L(X-X^2/2)=0, while multiplying by X gives moment -1 and
substituting X^2 gives moment -10. Taking f^q for q>1 generally produces
q+1 supports. Thus finite univariate exceptions cannot be declared empty
using any of these operations.

The four targets in Section 7 are now closed. Their proofs show why a larger
generic Hankel determinant was not needed: the effective invariants were the
actual normalized jet, simple residual roots with an unramified lift, and
scalar equations retaining every descended seed. The current next targets
are in the [continuation synthesis](STEPANOV_SFC_PROGRESS_AND_LOWER_THREE_CLOSURE.md),
including the remaining lower-four composite gaps and a precise polynomial-
observable contact bound toward broader support. A finite prefix alone does
not settle those unbounded domains.

## 9. Reproducibility, source failures, and commands

### 9.1 Existing current-source failures

At the inspected base revision:

| Existing check | Observed outcome |
|---|---|
| `verify_ordinary_binomial_terminal_ten.py` | PASS; complete 24,159 cases, byte-identical receipt |
| `verify_unrestricted_small_ratio.py` | PASS; byte-identical receipt |
| `verify_support_uniform_ordinary_completion.py` | FAIL: `Missing/stale source:large_thermal_checker` |
| `verify_large_thermal_quartic_fold_matching.py` | FAIL: `Missing/stale source:exact_checker` |
| `verify_compact_supercritical_ratio_uniform.py` | FAIL: `Missing/stale source:main_proof` |
| `make check-navigation` | FAIL: `changed supporting checker: FACTORIAL-UNIFORM-BOUNDED-QUARTIC-CUSP-EXCLUSION` |

The support-uniform source ledger has these two direct mismatches:

| File | Expected SHA256 | Current SHA256 |
|---|---|---|
| `verify_large_thermal_quartic_fold_matching.py` | `bbe51ac5537be9e6c98b6a983b71f5a16acc3bf9e5cacf2344042ae738d42549` | `ac81edd0e28421b4579e3822e422989b4e69bc0d76955a7b9912ca345ae1edfb` |
| `COMPACT_SUPERCRITICAL_RATIO_UNIFORM_EVENTUAL_SFC.md` | `0a4945c8a491c20a628db6e1fd37af555b70fa5e4b850d4202087c6bc5c4219f` | `f50a193ff8d583d817cfeb5218933a5ddcddbe6cd6762049d8485db9ce33d417` |

The relevant records are in
[support_uniform_ordinary_completion_source_ledger.json](support_uniform_ordinary_completion_source_ledger.json),
[large_thermal_quartic_fold_source_ledger.json](large_thermal_quartic_fold_source_ledger.json),
and [compact_supercritical_ratio_source_manifest.json](compact_supercritical_ratio_source_manifest.json).
These are evidence-identity failures, not mathematical counterexamples.
Blindly replacing the hashes would not reproduce the cited proof versions.
The original ledgers and receipts were preserved.

### 9.2 New exact reproduction

Run from the repository root. The complete independent terminal replay uses
only the Python standard library:

```sh
python -B research/factorial/image-rigidity/verify_stepanov_terminal_extension_independent.py \
  --repo-root . \
  --certificate research/factorial/image-rigidity/stepanov_terminal_11_12_bezout.jsonl.gz \
  --output /tmp/stepanov-terminal-independent.json
```

The alternate standard-library replay is
`verify_stepanov_terminal_extension.py`, with the same arguments. Adding
`--generate` regenerates all witnesses and requires Python-FLINT; use a
new output certificate path when reproducing generation. Its coefficient
construction differs from the independent replay.

The marked-rank controls use only the standard library:

```sh
python -B research/factorial/image-rigidity/verify_stepanov_marked_determinants.py \
  --output /tmp/stepanov-marked.json
```

They check 50 Vandermonde determinants, 100 inverse markers, 200 independently
expanded original contractions, 200 marker-operator identities, 50 raising
recurrences, 11 factorial Gram determinants, and the explicit negative
controls. These finite checks support the written universal proofs; their
bounded ranges do not establish universal quantifiers by extrapolation.

The auxiliary identities use exact SymPy arithmetic (tested with 1.14.0):

```sh
python -B research/factorial/image-rigidity/verify_stepanov_auxiliary_identities.py \
  --output /tmp/stepanov-auxiliary.json
```

They check the jet ranks, symbolic shift-resultant signs, product-Wronskian
identities, 28 actual shape-identity cases, four exact actual resultants,
and positive and negative derivative-closure examples. The native-prime
corollary checker is documented in its companion proof.

Retained new receipts are
[stepanov_terminal_extension_verification.json](stepanov_terminal_extension_verification.json),
[stepanov_terminal_extension_independent_verification.json](stepanov_terminal_extension_independent_verification.json),
[stepanov_marked_determinants_verification.json](stepanov_marked_determinants_verification.json),
and [stepanov_auxiliary_identities_verification.json](stepanov_auxiliary_identities_verification.json).
The [native-prime receipt](stepanov_native_prime_corollaries_verification.json)
includes 418 primitive-kernel/Hasse evaluations, 110 excluded valuation-minus-one
lifts, 90 rational-shape coefficient guards, and the ramified countercontrol,
in addition to the harmonic, congruence, and cyclic-rank controls.
The [package manifest](stepanov_investigation_manifest.json)
records the inspected revision, exact new file hashes,
and the separation between successful mathematical replays and the unrelated
baseline navigation/source failures.
