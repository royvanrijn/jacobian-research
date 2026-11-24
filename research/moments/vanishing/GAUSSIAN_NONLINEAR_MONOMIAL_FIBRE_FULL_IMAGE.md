# Actual Gaussian lattices for nonlinear monomial fibres

Let \(K\) be a characteristic-zero field. Use a normalized Wick pair
\(Z,W\) with \(\mathbb E[ZW]=1\), and \(s\ge1\) independent ordinary
Gaussian coordinates \(Y_1,\ldots,Y_s\), independent of that pair.
Fix \(m\ge2\) and integers \(k_i\ge1\), and put
\[
P=W(1+Z^m)\prod_{i=1}^sY_i^{k_i},\qquad
\Phi[A]=\mathbb E[Ae^{tP}],\qquad
M=K[t]\Phi[K[Z,W,Y_1,\ldots,Y_s]],\qquad V=K(t)M.
\]
All the series in this definition are formal original Wick periods.
This extends the [explicit linear-fibre result](GAUSSIAN_W_AFFINE_ROTATIONAL_DIVERGENCE.md)
when \(s=1,k_1=1\). In particular \(s=1,k_1\ge2\) is a nonlinear
\(Y\)-coefficient of \(W\), outside the earlier class
\(W(h(Z)+c(Z)Y)+G(Z,Y)\).

Let \(\ell=m\) if all \(mk_i\) are even and \(\ell=2m\) otherwise.
Set
\[
d=\ell/m,\quad e=\ell-d,\quad r_i=k_i\ell/2,\quad
\rho=\sum_i r_i,\quad x=t^\ell,\quad \vartheta=(t/\ell)\partial_t.
\]
For a positive integer \(v\) write \(G_v=\{1/v,2/v,\ldots,1\}\).
Define
\[
\mathcal N=G_e\ \cap\ \bigcup_i(G_{2r_i}\setminus G_{r_i}),\qquad
p=\ell+\rho-d-|\mathcal N|.
\]
Then the original unit is \(f(t)=g(x)\), where
\[
g(x)=\sum_{q\ge0} a_q x^q,\qquad
a_q=\binom{\ell q}{dq}\prod_i\frac{(2r_iq)!}{2^{r_iq}(r_iq)!}.
\]
Its cyclic differential subspace over \(K(t)\) has rank \(p\).
For each \(\beta\in\mathcal N\), define
\(\kappa_\beta(x)=\sum_{q\ge0}a_q x^q/(q+\beta)\).
A basis of the full actual space is
\[
\mathcal B=\bigl(1,g,\vartheta g,\ldots,\vartheta^{p-1}g,
                  (\kappa_\beta)_{\beta\in\mathcal N}\bigr).
\]
In particular \(\dim_{K(t)}V=\ell+\rho-d+1\). With
\(R=K[t,t^{-1}]\), the exact polynomial-time lattice is
\[
M[1/t]=R\mathcal B,\qquad M=(R\mathcal B)\cap K[[t]].
\]
Over \(\mathbb C\), the convergent Laurent germs in \(V\) form
exactly \(\mathbb C(t)\), the convergent elements of \(M\) are
exactly \(\mathbb C[t]\), and the maximal globally regular actual
differential submodule is the constant line.

For a single \(Y^k\), the extra-parameter set is
\(G_{\gcd(e,2r_1)}\setminus G_{\gcd(e,r_1)}\).
For example \(W(1+Z^5)Y^4\) has \(p=12\) and full rank \(15\).
Its two extra periods have denominators \(q+1/4\) and \(q+3/4\).
The obstruction applies to this whole monomial-fibre family
within \(\texttt{GAP-GAUSSIAN-REGULAR-TRANSCENDENTAL}\).
The unrestricted higher-rank realization question remains open.

## The unit differential order

Original Wick evaluation gives \(a_q\): the only nonzero unit
coefficients have \(m\mid n\) and all \(k_i n\) even, so \(n=\ell q\).
The coefficient ratio is \(a_{q+1}/a_q=C P_0(q)/Q_0(q)\), with \(C\ne0\),
and with monic coprime polynomials obtained by cancelling the
factor grids \(G_\ell+\sum_iG_{2r_i}\) and \(G_d+G_e+\sum_iG_{r_i}\).
A parameter \(\alpha\) denotes the factor \(q+\alpha\).
First cancel each \(G_{r_i}\) against its own \(G_{2r_i}\), and
\(G_d\) against \(G_\ell\). Since \(G_e\cap G_\ell=G_d\), the only
remaining cancellations of \(G_e\) are precisely the single
occurrences of the parameters in \(\mathcal N\).
Thus \(Q_0\) has simple parameters \(G_e\setminus\mathcal N\),
\(\deg P_0=p\), and \(\deg P_0-\deg Q_0=\rho>0\).
Also \(1\notin\mathcal N\), so \(Q_0(-1)=0\). In the multiple
Gaussian case a parameter in \(\mathcal N\) can still occur in
\(P_0\), since only one occurrence was cancelled. For example
\(m=3,s=2,k_1=k_2=1\) retains \(1/2\) in \(P_0\).
No reduced numerator parameter and reduced denominator parameter
differ by an integer, since they lie in \((0,1]\) and have
disjoint supports.

The recurrence gives the formal differential equation
\[
[Q_0(\vartheta-1)-CxP_0(\vartheta)]g=0.
\]
Its coefficient at \(x^0\) vanishes by \(Q_0(-1)=0\).
This is the Euler-operator form of
[DLMF equation 16.8.3](https://dlmf.nist.gov/16.8.E3), with
the parameter \(1\) accounting for the first Euler factor.
The displayed recurrence proves it formally in the divergent
case. The order assertion needed here is a minimality statement
for the actual unit, so it requires a lower bound beyond the form
of the hypergeometric equation.

Indeed, clear denominators in an equation with rational right side
and left side \(\sum_{j=0}^J x^jL_j(\vartheta)g\), \(L_J\ne0\).
At large coefficient index \(n\), dividing by \(a_{n-J}\) gives a
rational identity in \(n\). At \(n=J-\alpha\) for any numerator
parameter \(\alpha\), all terms with \(j<J\) vanish with at least
its numerator multiplicity, by the first factor
\(P_0(n-J)/Q_0(n-J)\) in the forward ratio. No denominator factor
can vanish there because the reduced parameter grids are disjoint
modulo integers. Thus \(P_0\) divides \(L_J\), which proves the order
lower bound \(p\), including for an equation with rational
inhomogeneity. For equations over \(K(t)\), clear denominators
and split the \(t\) exponents modulo \(\ell\); all operators on the
unit split into equations over \(K(x)\) of no larger order.
The unit derivatives are actual marked periods because
\(\partial_t^j\Phi[1]=\Phi[P^j]\). This proves their exact cyclic
rank and the independence of the constant line.

## Original markers and the additional periods

Write \(A=Z^aW^b\prod_iY_i^{c_i}\). At coefficient index \(n\), a
term survives only if
\[
j=(n+b-a)/m\in\mathbb Z,\quad 0\le j\le n,\quad
k_i n+c_i\text{ is even for every }i.
\]
Its coefficient, when it survives, is
\[
\frac{(n+b)!}{j!(n-j)!}\prod_i (k_i n+c_i-1)!!.
\]
The congruences either conflict or select one residue \(u\) modulo
\(\ell\), \(0\le u<\ell\). For \(n=\ell q+u\) and
\(v=(u+b-a)/m\), the coefficient for large \(q\) is \(a_q H(q)\)
for a rational \(H\). The ordinary Gaussian factorial shifts only
give the polynomials
\(\prod_{j'=1}^{(k_i u+c_i)/2}(2r_iq+2j'-1)\).
The first numerator factorial also has a nonnegative shift \(u+b\).
If \(v>0\), the factors \(dq+j'\) for \(1\le j'\le v\) from the
first denominator factorial are cancelled by the numerator factors
\(\ell q+mj'\), because \(mv=u+b-a\le u+b\). If \(v\le0\) this
denominator shift contributes a polynomial. Therefore all
remaining poles of \(H\) are simple and on the grid
\(-\beta-j'\), \(\beta\in G_e\), \(j'\in\mathbb Z_{\ge0}\), from
the second denominator shift.

The constant is an actual time-polynomial period:
\[
\Phi\!\left[1-tmZ^{m-1}\prod_iY_i^{k_i}\right]=1.
\]
Conditional on the ordinary coordinates, the amplitude is the
Jacobian of the selected circular inverse
\(z=t(1+z^m)\prod Y_i^{k_i}\), so the identity follows from the
[original circular Wick residue](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#exact-circular-elimination).
Equivalently, at order \(n>0\) the coefficient identity is
\([z^n](1+z^m)^n=[z^{n-1}]m z^{m-1}(1+z^m)^{n-1}\),
obtained by differentiating the polynomial.

All simple poles on the \(Q_0\)-grid are generated in the unit
cyclic space modulo finite Laurent corrections. Indeed \(x^{-J}g\)
has tail multiplier
\(C^J\prod_{h=0}^{J-1}P_0(q+h)/Q_0(q+h)\), with coprime numerator
\(N_J\) and denominator \(D_J\). For any rational \(H\) supported
only on the \(Q_0\)-grid, choose \(J\) large enough that \(D_JH\)
is polynomial. Find polynomials \(U,V_0\) with
\(H=U N_J/D_J+V_0\) by inverting \(N_J\) modulo \(D_J\).
The actual Laurent expressions
\(U(\vartheta)x^{-J}g+V_0(\vartheta)g\) then have the desired tail.

For \(\beta\in\mathcal N\), put \(b_\beta=e\beta\in\mathbb Z_{>0}\).
The spatial marker \(Z^{m b_\beta}\) has residue \(u=0\) and tail
multiplier
\[
\frac{\prod_{h=0}^{b_\beta-1}(dq-h)}
     {\prod_{h=1}^{b_\beta}(e q+h)} .
\]
Its residue at \(q=-\beta\) is nonzero, and all other poles are
at parameters smaller than \(\beta\). Starting at the least
element of \(\mathcal N\), subtract the terms on the \(Q_0\)-grid
and the already obtained \(\kappa\) terms, along with any
polynomial part of the rational multiplier. The available
constant line corrects the finite Laurent difference. Dividing
by the nonzero residue then realizes the exact \(\kappa_\beta\)
in \(M[1/t]\). For shifted poles use \(x^{-J}\kappa_\beta\).
Its only pole outside the \(Q_0\)-grid is at \(q=-J-\beta\),
with nonzero residue. The argument for noncancellation even
if \(\beta\) occurs in \(P_0\) is that all sampled arguments
\(q+h\), \(0\le h<J\), are strictly less than \(-1\) at this
point, while every root of \(P_0,Q_0\) is in \([-1,0)\).
Removing all such poles and using the cyclic construction
proves generation of every marked tail, hence of \(V\).

For independence, after clearing denominators and taking a
single residue of the time exponents, suppose that a relation
has a nonzero term \(A_{\beta,k}x^k\kappa_\beta\). For large
indices \(n\), divide by \(a_{n-J}\), where \(J\) is the maximal
shift in that relation. The term with coefficient \(A_{\beta,k}\)
then gives a nonzero simple pole at \(n=k-\beta\); its prior
forward-ratio factors, if any, are evaluated strictly to the
left of \(-1\). The other \(\kappa\) terms cannot have this pole,
as their \(\beta\)'s lie in \((0,1)\) and the shifts are integers.
The factorial ratios in other terms have no pole in that residue
modulo one because \(Q_0\) has no parameter \(\beta\). This
contradicts the rational identity. Together with the unit lower
bound this proves the asserted basis. Each additional period
satisfies \((\vartheta+\beta)\kappa_\beta=g\), so its quotient
line is regular of local exponent \(-\ell\beta\) in the \(t\)
coordinate. This exponent can be nonintegral, as in the example
\(m=5,k_1=4,\beta=1/4\).

## Polynomial time and the regular core

The entire generation proof uses Laurent time shifts, polynomial
Euler operators and finite Laurent corrections. The unit operator's
leading coefficient \(-Cx\) is a Laurent unit, so high Euler
derivatives reduce into the displayed basis over \(R\). This proves
\(M[1/t]=R\mathcal B\). The [original Wick Stein backward shift](GAUSSIAN_W_AFFINE_ROTATIONAL_DIVERGENCE.md#the-polynomial-time-lattice)
preserves \(M\), since it applies to every fixed polynomial phase
and uses only triangular centered-polynomial reduction. Thus
\(h\mapsto(h-h(0)f)/t\) proves the origin saturation
\(M=(R\mathcal B)\cap K[[t]]\). In particular, given Laurent
coordinates \(h=\sum r_b(t)b(t)\) with \(t^Jr_b(t)\in K[t]\),
the \(J\) coefficients of \(t^Jh\) of indices \(0,\ldots,J-1\)
give the complete finite membership test.

For any polynomial-time combination the tails of each residue
class modulo \(\ell\) are \(a_q\) times a rational function of \(q\).
This holds also after finite shifts of the coefficients, by the
unit ratio. A nonzero rational multiplier has polynomial-size
asymptotics on the positive integers, whereas the unit ratio is
asymptotic to \(Cq^\rho\) with \(\rho>0\). That residue forces
divergence. Convergence in \(M\) therefore requires vanishing of
every tail and yields a polynomial. Clearing rational
denominators gives the statement for \(V\). The regular-singular
Frobenius implication used in the
[Gaussian regular-core audit](GAUSSIAN_G_FUNCTIONS_AND_GLOBAL_REGULAR_CORE.md)
then confines every globally regular actual submodule to the
constant line. No \(G\)-operator or hypergeometric irreducibility
theorem is invoked for the rank or the convergence argument.

## Unit divergence for arbitrary nonconstant spatial coefficient

Let \(c\in\mathbb C[Z]\) be nonconstant, with \(c(0)\ne0\).
For the same independent ordinary Gaussian variables and positive
integer powers \(k_i\), consider
\[
P=Wc(Z)\prod_iY_i^{k_i}+a_0+\sum_i a_iY_i^2 ,
\qquad a_i\in\mathbb C .
\]
Then the actual unit diverges and the original marked periods of
\(1\) and \(P\) are linearly independent over \(\mathbb C(t)\).
The proof allows coincident dominant moduli and requires neither
coefficient positivity nor rotational symmetry of \(c\). This is
the slice with no added spatial perturbations of the
[tied-singularity problem](GAUSSIAN_W_AFFINE_SELECTED_DOMINANT_DIVERGENCE.md#4-the-concrete-bottleneck-beyond-these-criteria).
The full marked rank and lattice above concern only the
binomial spatial coefficient.

First set \(a_0=a_1=\cdots=a_s=0\). For \(X\) near zero, the selected
inverse \(z=Xc(z)\), \(z(0)=0\), gives the algebraic kernel
\[
K(X)=\frac{1}{1-Xc'(z(X))}.
\]
Exact circular Wick contraction gives
\([X^n]K(X)=[z^n]c(z)^n\), and hence
\[
[t^n]\Phi[1]=[X^n]K(X)\prod_i\mathbb E[Y_i^{k_i n}] .
\]
The unit only needs the even part
\(K_{\mathrm{ev}}(X)=(K(X)+K(-X))/2\) for this argument.
It is a nonpolynomial algebraic germ. Indeed, on every branch
of \(X=z/c(z)\) above sufficiently large \(|X|\), \(z\) lies
near some finite root \(z_*\ne0\) of \(c\). To see this uniformly
over the branches, the rational function \(z/c(z)\) has a
finite limit at spatial infinity because \(\deg c\ge1\).
The preimage of \(|X|\) large therefore lies near the finitely
many poles, exactly the roots since \(c(0)\ne0\). At a root of
multiplicity \(\mu\),
\[
1-\frac{z c'(z)}{c(z)}
   =-\frac{\mu z_*}{z-z_*}+O(1).
\]
Thus \(K(X)\to0\) at infinity on every branch, and the same
holds for \(K(-X)\). Both selected germs can be continued along
a common path to infinity avoiding a finite set of branch
values. A polynomial \(K_{\mathrm{ev}}\) would continue there
as the same polynomial and tend to zero, contrary to
\(K_{\mathrm{ev}}(0)=1\). Equivalently, an entire continuation
of the germ would be bounded at infinity by the same
branch-uniform estimate. This proves that its power series
has finite positive convergence radius.

By the coefficient formula for the radius, there is an infinite
subsequence of even \(n\) with
\(\left|[X^n]K(X)\right|\ge c_0^n\) for some \(c_0>0\).
On this subsequence every \(k_i n\) is even, so each ordinary
Gaussian moment is the positive double factorial
\((k_i n-1)!!\). These factors grow faster than every fixed
exponential along the subsequence. There is no sign cancellation
because the angular coefficient is already the exact full
coefficient \([X^n]K(X)\). This proves actual-unit divergence.

For the diagonal quadratic additions, put \(D_i(t)=1-2a_i t\).
Exact formal Gaussian completion gives for every integer \(j\ge0\)
\[
\mathbb E[Y_i^j e^{t a_iY_i^2}]
   = D_i(t)^{-(j+1)/2}\,\mathbb E[Y_i^j],
\]
with the unit branch of the power at zero. Both sides vanish
for odd \(j\). For even \(j=2h\), the coefficientwise identity
follows by factoring \((2h-1)!!\) out of
\((2h+2n-1)!!\) in each coefficient and summing the binomial
series. If \(f_0\) denotes the unperturbed actual unit, the full
unit is therefore
\[
e^{a_0t}\prod_i D_i(t)^{-1/2}\,
 f_0\!\left(t\prod_iD_i(t)^{-k_i/2}\right).
\]
The prefactor is an analytic unit and the time substitution
has derivative one at zero, hence analytic inverse. Divergence
persists. Finally \(\Phi[P]=\partial_t\Phi[1]\). A rational
logarithmic derivative would be regular at zero since
\(\Phi[1](0)=1\), and its first-order equation would make the
unit analytic. This proves the claimed marked independence.

## Positive moment law and its determinacy

The raw coefficients \(a_q\) on the stride \(t^\ell\) form a
Stieltjes moment sequence, despite the divergence just proved.
Their representing law on \([0,\infty)\) is determinate exactly
when \(\rho\le2\).

First, \(\binom{\ell q}{d q}=\mathbb E[B^q]\) for a compact
positive random variable \(B\). Order the numerator grid \(G_\ell\)
as \((\alpha_j)\) and the denominator multiset \(G_d+G_e\) as
\((\beta_j)\). The inequality
\(\lfloor\ell u\rfloor\ge\lfloor du\rfloor+\lfloor eu\rfloor\)
for \(0\le u\le1\) implies \(\alpha_j\le\beta_j\) term by term
by comparing the numbers of grid elements below each level.
Take \(B=(\ell^\ell/(d^d e^e))\prod_j B_j\) with independent
\(B_j\sim\operatorname{Beta}(\alpha_j,\beta_j-\alpha_j)\),
taking \(B_j=1\) at equality. The ordinary beta integral and
factorial grouping give the claimed moments. This is the
zero-shift integer-grid instance of
[Młotkowski–Penson, Lemmas 3.1–3.2 and Theorem 3.3](https://arxiv.org/pdf/1309.0595#page=8).
Their proof uses Gauss multiplication, matches the ordered
upper and lower gamma parameters, and normalizes at moment
zero before taking a Mellin product. The floor argument above
checks the matching independently for these grids.

Take \(B\) independent of all \(Y_i\). Then
\(T=B\prod_i |Y_i|^{2r_i}\) has moments \(a_q\).
If \(\rho\le2\), the ratio or Stirling's formula gives
\(a_q\le D^q(2q)!\) for a fixed \(D\).
For any representing law, the symmetrized square root \(S\)
has even moments \(a_q\) and zero odd moments.
By the power series for \(\cosh\) and monotone convergence,
\(\mathbb E\cosh(\delta S)<\infty\) for sufficiently small
\(\delta>0\). Its characteristic function is analytic in a
strip and determined by those moments, so \(S^2\) has only
one possible law.

If \(\rho>2\), the gamma-power variables \(X_i=|Y_i|^{2r_i}\)
have positive densities of the form
\(f_i(u)=c_i u^{\gamma_i-1}\exp(-u^{1/r_i}/2)\),
\(\gamma_i=1/(2r_i)\). Let \(f_T\) denote the density obtained by
conditioning on \(B\prod_{i<s} X_i\) and scaling \(f_s\).
For \(x\ge1\) restrict \(B\) to any compact interval in
\((0,\infty)\) with positive mass and \(X_i\) for \(i<s\) to
\([x^{r_i/\rho},2x^{r_i/\rho}]\). The last density is then
evaluated at a point comparable to \(x^{r_s/\rho}\).
The interval measures and the final density divided by the
scaling factor give
\[
f_T(x)\ \ge\ c x^{(\sum_i r_i\gamma_i)/\rho-1}e^{-C x^{1/\rho}} .
\]
For \(0<x\le1\) condition instead on the scale in one fixed compact
positive interval of positive probability. This gives
\(f_T(x)\ge c' x^{\gamma_s-1}\). The case \(s=1\) uses just the
compact interval for \(B\) in both arguments.

Choose \(\lambda\ge\gamma_s\) and \(D_0>C\), and a sufficiently
small positive \(c_0\). Then
\(h(x)=c_0x^{\lambda-1}e^{-D_0x^{1/\rho}}\le f_T(x)\) for \(x>0\).
Set \(\theta_0=\pi/\rho<\pi/2\). For all integers \(q\ge0\),
\[
\int_0^\infty x^q h(x)
 \sin\bigl(D_0\tan\theta_0\,x^{1/\rho}-\pi\lambda\bigr)\,dx=0.
\]
Indeed, use \(v=x^{1/\rho}\) and the gamma integral with complex
rate \(D_0(1-i\tan\theta_0)\). After the fixed phase
\(-\pi\lambda\), its imaginary part has phase
\(\rho(q+\lambda)\theta_0-\pi\lambda=q\pi\).
Adding or subtracting one half of the resulting nonzero signed
density gives distinct positive probability measures with the
same moments.

The large-\(x\) bound is the explicit gamma-power instance of the
density-window step in
[Lin–Stoyanov, proof of Theorem 7 and Corollary 1](https://arxiv.org/pdf/1406.1654#page=8).
Their proof obtains a tail bound for the product density and
applies Krein's criterion to the gamma product. Here the bounded
additional scale \(B\) is retained, and the explicit zero-moment
perturbation proves indeterminacy of this particular mixed
product directly. Neither the positive moment law nor its
determinacy supplies convergence of the original Wick series.

### All-order Hankel minors and the order-moment source

Both \(b_q=\binom{\ell q}{d q}\) and \(a_q\) have strictly positive
Hankel minors on increasing row and column indices, including the
indeterminate cases \(\rho>2\) for \(a_q\). More generally their
interpolations \(z\mapsto\mathbb E[B^z]\) and
\(z\mapsto\mathbb E[T^z]\) have strictly totally positive additive
kernels for nonnegative real orders. Both measures have infinite
positive support. For \(B\), the least numerator parameter
\(1/\ell\) is smaller than the least denominator parameter
\(1/\max(d,e)\), so at least one beta factor is nonconstant and
absolutely continuous. For \(T\), at least one additional factor
has a positive density on \((0,\infty)\).

This applies the exact Mellin transport in
[Salazar, arXiv:2606.31647v2, Theorems 2.1–2.3](https://arxiv.org/pdf/2606.31647v2#page=5).
For either positive measure \(\mu\), the determinant at orders
\(u_i+v_j\) is, by expansion and symmetrization,
\[
\frac1{r!}\int
 \det[t_k^{u_i}]_{i,k=1}^r
 \det[t_k^{v_j}]_{j,k=1}^r\,d\mu(t_1)\cdots d\mu(t_r).
\]
Fubini follows by expanding the finite products and using the
moments at the specified sums. On an ordered chamber the two
generalized Vandermonde determinants are positive: pass to
\(\log t_k\) and use the exponential Wronskians, which are a
positive exponential times an ordinary Vandermonde in the order
parameters. The Chebyshev induction on initial Wronskians, or
successive Rolle's theorem after division by the first
exponential, gives nonvanishing. Coalescence fixes the sign.
All chambers give nonnegative contributions, since simultaneous
permutation of the nodes multiplies both determinants by the
same sign. Strict positivity follows by taking a product of
disjoint ordered neighborhoods of \(r\) support points with
positive mass. Thus the determinant property is compatible
with a nonunique representing law and with a divergent actual
unit. It gives no analytic Gaussian realization of a Markov
generating function for the unnormalized coefficients.

The same proof in
[Salazar, Theorem 4.1, equations (4.5)–(4.9)](https://arxiv.org/pdf/2606.31647v2#page=10)
uses a discrete Mellin measure for positive Jackson
\(\Gamma_q\), \(0<q<1\). The reciprocal product generating
function \(F(\xi)=1/(\xi;q)_\infty\) satisfies
\(F(\xi/q)(1-\xi/q)=F(\xi)\) for \(|\xi|<q\), equivalently
\(F(q\xi)=(1-\xi)F(\xi)\) for \(|\xi|<1\).
Its coefficients therefore satisfy
\((1-q^n)c_n=c_{n-1}\), \(c_0=1\). Substituting \(\xi=q^z\)
gives the distinct atoms \(q^n/(1-q)\) with weights
\((1-q)(q;q)_\infty/(q;q)_n\). The weights are bounded, so
the geometric tail has finite moments for \(z>0\), even
though the unweighted measure has infinite mass. This verifies
the actual domain required in the Mellin determinant argument.
The reciprocal \(1/\Gamma_q\) and the ordinary gamma quotient
with numerator shift larger than the denominator shift are
separate problems. The Gaussian application here uses its
own positive measure rather than a Jackson specialization.
