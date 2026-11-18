# Pair-primitive binary sums: a quaternary local-global deduction

Fix two primitive integral positive-definite binary quadratic forms
\(Q_i(z)=\tfrac12 z^{\mathsf T}A_i z\), with symmetric integral Gram
matrices \(A_i\) of even diagonal. The forms are fixed throughout; no
class-number-one or maximal-order hypothesis is imposed.

This note proves eventual existence and a power-saving asymptotic for
**coordinate representations**, with each binary pair primitive. It does
not prove an asymptotic for the number of distinct value pairs. The argument
is a deduction from classical quaternary theta-series estimates, using
two elementary Möbius inversions. Priority for this particular deduction
has not been established. This is a research proof, without external review.

## Statement and local criterion

Let
\[
 R^*(N)=\#\{(z,w)\in\mathbb Z^2\times\mathbb Z^2:
 Q_1(z)+Q_2(w)=N,\ \gcd(z_1,z_2)=\gcd(w_1,w_2)=1\}.
\]
For each prime let
\[
 \beta_p^*(N)=\lim_{h\to\infty}p^{-3h}
 \#\{(z,w)\bmod p^h:F(z,w)=N\bmod p^h,\
              z\not\equiv0\bmod p,\ w\not\equiv0\bmod p\},
 \qquad F=Q_1\oplus Q_2 .
 \tag{1}
\]
Write
\[
 J=\frac{4\pi^2}{\sqrt{\det A_1\det A_2}}.
\]

**Theorem.** There is an absolute positive \(\eta\) (for example the
deliberately nonoptimized value \(1/1000\)) such that
\[
 R^*(N)=JN\prod_p\beta_p^*(N)+O_{Q_1,Q_2}(N^{1-\eta}).
 \tag{2}
\]
For every fixed \(\epsilon>0\), on the locally admissible targets,
\[
 N^{-\epsilon}\ll_{Q_1,Q_2,\epsilon}
       \prod_p\beta_p^*(N)
       \ll_{Q_1,Q_2,\epsilon}N^\epsilon .
 \tag{3}
\]
Consequently every sufficiently large positive integer admitting a
pair-primitive solution over every \(\mathbb Z_p\) is globally represented
with both pairs primitive. This holds for forms in any ideal class and
for norm forms of nonmaximal quadratic orders.

Local admissibility is a finite, explicitly enumerable congruence test.
Let
\[
 S=\{p:p\mid 2\det A_1\det A_2\},\qquad
 t_p=\max\{v_p(\det A_1),v_p(\det A_2)\}.
\]
It is equivalent to soluble pair-primitive congruences modulo
\(p^{2t_p+1}\) for all \(p\in S\). At primes outside \(S\), every target
has a pair-primitive solution. Thus an admissibility modulus is
\[
 M=\prod_{p\in S}p^{2t_p+1}.
 \tag{4}
\]
This modulus need not be minimal.

A local obstruction gives exact nonrepresentation, including below the
eventual cutoff. Formula (2) does not determine the finite list of
admissible exceptions or a numerical cutoff. The separate specialized
\(-8\) sieve below supplies an explicit cutoff for that form.

## Exact good-prime factors

For \(p\notin S\), set
\[
 \varepsilon_i=\left(\frac{-\det A_i}{p}\right)\in\{-1,1\}.
\]
The local factors are
\[
 \beta_p^*(N)=
 \begin{cases}
 1-\dfrac{2+\varepsilon_1\varepsilon_2}{p^2}
     +\dfrac{\varepsilon_1+\varepsilon_2}{p^3},&p\nmid N,\\[6pt]
 1+\dfrac{\varepsilon_1\varepsilon_2}{p}
  -\dfrac{2+\varepsilon_1+\varepsilon_2+
                  \varepsilon_1\varepsilon_2}{p^2}
  +\dfrac{1+\varepsilon_1+\varepsilon_2}{p^3},&p\mid N .
 \end{cases}
 \tag{5}
\]
They depend only on whether \(p\mid N\), even when \(v_p(N)\) is large.
In particular they are positive for every \(p\ge3\).

To verify (5), the number of vectors representing \(a\) under a
nondegenerate binary form over \(\mathbb F_p\) is \(p-\varepsilon_i\)
when \(a\ne0\), and \(p+(p-1)\varepsilon_i\) when \(a=0\).
Subtract the zero vector at \(a=0\), convolve the two resulting functions,
and divide by \(p^3\). Every surviving root is nonsingular, so every
root modulo \(p\) has exactly \(p^{3(h-1)}\) lifts modulo \(p^h\).
This also establishes the existence assertion outside \(S\).

On admissible targets the resulting singular series also satisfies
\[
 \prod_p\beta_p^*(N)\asymp_{Q_1,Q_2}
 \prod_{\substack{p\mid N\\p\notin S}}
       \left(1+\frac{\left(\frac{\det A_1\det A_2}{p}\right)}p\right).
 \tag{5a}
\]
The ratios of the good-prime factors to the displayed target factors
are \(1+O(p^{-2})\), and the bad-prime factors have uniform positive
bounds. In the self-sum case the target symbol is always \(+1\).

## Bad primes: finite testing and uniform density

If \(v\in\mathbb Z_p^2\) is primitive, the identity
\(\operatorname{adj}(A_i)A_i v=(\det A_i)v\) shows that some component of
\(A_i v\) has valuation at most \(v_p(\det A_i)\). At a pair-primitive
solution the gradient of \(F\) therefore has a component of valuation
\(s\le t_p\).

A solution modulo \(p^{2t_p+1}\) consequently lifts by univariate Hensel:
fix three variables and solve for a coordinate whose derivative has
valuation \(s\). This proves both directions of the criterion (4).
Conversely, given an exact solution, fix those three coordinates modulo
\(p^{2s+1}\). Hensel supplies \(p^s\) choices of the remaining coordinate
modulo \(p^h\), for every \(h\ge2s+1\). All these lifts keep the prescribed
nonzero blocks. Hence
\[
 \beta_p^*(N)\ge p^{-5s-3}\ge p^{-5t_p-3}.
 \tag{6}
\]
The limit exists: partition the primitive vectors into the finitely many
residue neighborhoods on which the valuation of a selected derivative
is fixed and apply the same Hensel parameterization.

The same argument gives a uniform upper bound for the density with only
the *whole* four-vector required primitive. Its derivative has valuation
bounded in terms of \(\det A_1\det A_2\). Splitting off vectors divisible
by \(p\) then gives the recursion
\[
 \beta_p(F;N)=\beta_p^{\rm vec}(F;N)
       +\mathbf1_{p^2\mid N}p^{-2}\beta_p(F;N/p^2).
 \tag{7}
\]
It follows that \(\sup_{N\ge1}\beta_p(F;N)<\infty\), for each fixed \(p\).
Equations (5)--(7) imply (3): the nontarget factors converge absolutely,
the finitely many bad factors are bounded above and below on admissible
targets, and the remaining target factors are \(1+O(1/p)\).
The elementary divisor estimate bounds their product and its reciprocal
by \(N^\epsilon\). No bound on target valuations is being assumed.

## The imported theta-series estimate

The analytic input is Jeremy Rouse,
*Integers represented by positive-definite quadratic forms and Petersson
inner products*, [arXiv:1802.03437v1](https://arxiv.org/abs/1802.03437),
9 February 2018. The precise parts needed are equation (1), PDF page 5
(the genus density formula), Theorem 3, PDF page 10, and the **general
all-coefficient** bound at the start of Section 5, PDF page 15, as
combined on PDF page 16:
\[
 r_T(n)=E_T(n)+c_T(n),\qquad
 |c_T(n)|\ll_\epsilon
 \max\{L(T)^{5/4+\epsilon}D(T)^{1/8+\epsilon},
                         L(T)^{3/2+\epsilon}\}
       \tau(n)\sqrt n .
 \tag{8}
\]
Here \(T\) is a primitive positive-definite integral quaternary form,
\(D(T)\) its Gram determinant, and \(L(T)\) its level.
The genus formula, with the real factor computed independently in the
Gram convention of this note, gives
\[
 E_T(n)=\frac{4\pi^2 n}{\sqrt{D(T)}}\prod_p\beta_p(T;n).
 \tag{9}
\]
The real factor follows by differentiating the ellipsoid volume
\(2\pi^2n^2/\sqrt{D(T)}\). It is not transcribed from the
prefactor displayed in Rouse's Section 3: that line displays
\(\pi^2n/\sqrt{D(T)}\), despite the stated Gram convention.
The independent real-volume calculation resolves this apparent
factor-of-four discrepancy for the convention used here. For four squares, \(A=2I_4\) and
the real density is \(\pi^2n\), fixing the factor 4 above.

These statements apply to all positive \(n\); the sharper bound that
requires \((n,L(T))=1\) is not used. The general estimate is effective.
A copy of the exact v1 PDF and its SHA-256 are recorded in
[the source ledger](primitive_norm_quaternary_source_ledger.json).
Rouse derives the general coefficient estimate from
Schulze-Pillot--Yenirce, Theorem 10; this deduction uses (8) as a single
imported bound, rather than the existential conclusion of Rouse's
Theorem 1.

## Möbius truncation

For \(d,e\ge1\), let
\[
 T_{d,e}=d^2Q_1\oplus e^2Q_2,\qquad
 r'_{d,e}(N)=\#\{z,w\ne0:T_{d,e}(z,w)=N\}.
\]
Positive definiteness bounds all divisors relevant below by
\(O_{Q_1,Q_2}(\sqrt N)\). Exact Möbius inversion gives
\[
 R^*(N)=\sum_{d,e\ge1}\mu(d)\mu(e)r'_{d,e}(N).
 \tag{10}
\]
For a fixed positive-definite binary form \(Q=ax^2+bxy+cy^2\),
\[
 r_Q(m)\ll_{Q,\epsilon}m^\epsilon\qquad(m\ge1).
 \tag{11}
\]
For completeness, associate a solution to the algebraic integer
\(2ax+by+y\sqrt{b^2-4ac}\), of norm \(4am\), in the maximal order of
its imaginary quadratic field. There are at most \(\tau(4am)\)
integral ideals of that norm, and finitely many generators per principal
ideal. This proves (11), without a class-number-one assumption.

If \(d>D\), choose \(m=Q_1(z)\ge1\). There are \(O(N/d^2)\) possible
\(m\), and (11) bounds the choices in both blocks by \(N^\epsilon\).
The unrestricted sum over the divisors \(e\) of the second coordinate
gcd costs another \(N^\epsilon\). The same argument with the blocks
reversed therefore proves
\[
 R^*(N)=\sum_{d,e\le D}\mu(d)\mu(e)r'_{d,e}(N)
                 +O_{Q_1,Q_2,\epsilon}(N^{1+\epsilon}/D).
 \tag{12}
\]
Neither block is zero here; this exclusion is essential to the estimate.

Put \(g=(d,e)\). If \(g^2\nmid N\), the corresponding representation
and genus coefficients vanish. Otherwise divide \(T_{d,e}\) and \(N\)
by \(g^2\). The resulting form is primitive because each \(Q_i\) is
primitive and \((d/g,e/g)=1\). Its determinant is
\(D(F)((d/g)(e/g))^4\), and its level divides
\(L_0\operatorname{lcm}((d/g)^2,(e/g)^2)\), where \(L_0\) is any
common multiple of the two fixed binary levels. In particular (8)
gives, uniformly,
\[
 r_{d,e}(N)=E_{d,e}(N)
             +O_{Q_1,Q_2,\epsilon}((de)^{3+\epsilon}N^{1/2+\epsilon}).
 \tag{13}
\]
Allowing the whole form to have content \(g^2\) thus causes no gap.

The theta coefficient includes zero blocks. Their removal in the
truncated sum costs \(O(D^2N^\epsilon)\), by (11). Combining (12),(13)
and enlarging the harmless power \(3+\epsilon\) to 4 gives
\[
 R^*(N)=\sum_{d,e\le D}\mu(d)\mu(e)E_{d,e}(N)
 +O\left(N^{1+\epsilon}/D+D^{10}N^{1/2+\epsilon}
                                      +D^2N^\epsilon\right).
 \tag{14}
\]

## Absolute completion of the genus main term

For \(a,b\in\{0,1\}\), define \(\lambda_p(a,b;N)\) to be the local
surface density of \(F=N\) with the first block divisible by \(p\)
if \(a=1\), and the second block divisible by \(p\) if \(b=1\).
Thus
\[
 \lambda_p(a,b;N)=p^{-2a-2b}
       \beta_p(p^{2a}Q_1\oplus p^{2b}Q_2;N).
 \tag{15}
\]
The coefficient-change Jacobian in (15) is essential. Formula (9)
becomes, for squarefree \(d,e\),
\[
 E_{d,e}(N)=JN\prod_p
       \lambda_p(\mathbf1_{p\mid d},\mathbf1_{p\mid e};N).
 \tag{16}
\]
All these quantities are nonnegative.

At a good prime,
\[
 \lambda_p(1,1;N)
 =\mathbf1_{p^2\mid N}p^{-2}\beta_p(F;N/p^2)\ll p^{-2}.
 \tag{17}
\]
For \(\lambda_p(1,0)\), separate solutions whose second block is
primitive from those whose second block is divisible by \(p\).
The latter contribute (17). In the former, fixing the first block,
which lies in \(p\mathbb Z_p^2\), leaves a binary equation for a
primitive second block. At a good prime its binary primitive density
is at most 2: a split binary form becomes a unit times \(uv\), with
one coordinate a unit, while an inert one represents only units
primitively. Integrating over the first block gives at most \(2p^{-2}\).
Consequently, uniformly in \(N\),
\[
 \lambda_p(1,0),\lambda_p(0,1),\lambda_p(1,1)\ll p^{-2}.
 \tag{18}
\]
At bad primes all four \(\lambda_p\)'s are bounded by the uniformly
bounded unrestricted density in (7).

For \(p\nmid N\) good, \(\lambda_p(0,0)=1+O(p^{-2})\);
for \(p\mid N\) good, it is \(1+O(p^{-1})\), uniformly in the valuation.
Taking \(0<\alpha<1/2\), (16)--(18) and multiplicativity imply
\[
 \sum_{d,e\ge1}\mu(d)^2\mu(e)^2(de)^\alpha E_{d,e}(N)
              \ll_{Q_1,Q_2,\alpha,\epsilon}N^{1+\epsilon}.
 \tag{19}
\]
Indeed the good-prime Euler factor is
\[
 \lambda_p(0,0)+p^\alpha
       (\lambda_p(1,0)+\lambda_p(0,1))
                 +p^{2\alpha}\lambda_p(1,1).
\]
Away from target primes this is \(1+O(p^{-2+2\alpha})\).
At target primes the extra \(O(1/p)\) product is \(O_\epsilon(N^\epsilon)\).
Bad primes are a fixed finite set. The argument still applies if some
unrestricted bad density vanishes, since then every term vanishes;
it does not divide by a zero density.

It follows from (19) that the main term in (14) can be completed with
error \(O(N^{1+\epsilon}D^{-\alpha})\). Local inclusion-exclusion gives
\[
 \sum_{d,e\ge1}\mu(d)\mu(e)E_{d,e}(N)
 =JN\prod_p(\lambda_p(0,0)-\lambda_p(1,0)
                  -\lambda_p(0,1)+\lambda_p(1,1))
 =JN\prod_p\beta_p^*(N).
 \tag{20}
\]
Choosing \(\alpha=1/4\), \(D=N^{1/100}\), and then sufficiently small
\(\epsilon\), proves (2), with \(\eta=1/1000\).
Equations (3) and (2) prove eventual existence.

## Quadratic-kernel corollary of the existing splitting theorem

There is also a direct application of
[SubgroupSemigroupScale](ASYMMETRIC_SUBGROUP_SPLITTING.md), without the
quaternary argument. For any two fixed nonprincipal quadratic literal
characters modulo \(L\), take the first subgroup to be the kernel of the
first character and the second quadratic kernel to be that of the second.
The source then gives (22) for their prime-factor semigroups.

For the norm interpretation, let \(K_i\) have class number one and fundamental
discriminant \(\Delta_i<0\). Choose a fixed even \(L\) divisible by
\(|\Delta_1|\) and \(|\Delta_2|\). Take \(Q_i\) to be an integral-basis norm
form of the maximal order of \(K_i\), inflate both quadratic characters
literally modulo \(L\), and put
\[
 H=\ker\chi_{\Delta_1},\qquad K=\ker\chi_{\Delta_2}
       \quad\text{in }(\mathbb Z/L)^* .
\]
Define the unramified, \(L\)-coprime primitive norm cores
\[
 P_i^0=\{n\ge1:(n,L)=1,\ n=Q_i(x,y),\ (x,y)=1\}.
\]
Then
\[
 P_i^0=\{n\ge1:p\mid n\Longrightarrow\chi_{\Delta_i}(p)=1\}.
 \tag{21}
\]
Indeed an inert prime dividing a norm divides both coordinates; the
converse follows by choosing one oriented prime ideal above each split
prime, then using class number one to make their product principal.
The resulting generator cannot have both coordinates divisible by a
rational prime, since its ideal contains at most one of the two
conjugate prime factors.

Set
\[
 \Gamma(N)=\frac{L}{|H||K|}
                 \#\{(a,b)\in H\times K:a+b=N\bmod L\},\qquad
 J_{12}(N)=\prod_{\substack{p\mid N,\ p\nmid L\\
                     \chi_{\Delta_1}(p)=\chi_{\Delta_2}(p)=-1}}
                            \frac{p-1}{p-2}.
\]
Specializing the subgroup theorem to index \(\kappa=2\) gives
\[
 \#\{a+b=N:a\in P_1^0,\ b\in P_2^0\}
       \asymp_{K_1,K_2,L}\frac{N J_{12}(N)}{\log N}
 \tag{22}
\]
for every sufficiently large target with \(\Gamma(N)>0\), and exact
zero when \(\Gamma(N)=0\). The kernels need not be nested and the fields
may coincide. The weighted asymptotic in that source specializes to
the two weights \(\tau(n)\mathbf1_{P_i^0}(n)\).

Equation (22) is an existing-theorem corollary, not a new analytic
theorem. It gives the sharper logarithmic order for distinct values
on these cores, whereas (2) gives an explicit leading main term for
coordinate counts with full local ramification. Multipliers supported
on \(p\mid L\) are outside the allowed coefficient hypotheses of the
subgroup theorem. Passing from (22) to all ramified norm strata is not
justified by calling such multipliers fixed; the quaternary deduction
above supplies eventual existence for those strata independently.

## Angular equidistribution and fixed macroscopic balance

The full primitive count also supplies the large-prime tail needed to
pass fixed-lattice harmonic equidistribution through the sieve. This
extends the fixed-lattice argument already used for
[Eisenstein shells](EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md).

Let \(T\) be the real block diagonal linear map with
\(F(v)=|Tv|^2\), and give \(S^3\) its probability surface measure
\(\sigma\). Write \(S^*(N)=\prod_p\beta_p^*(N)\).

**Corollary PairPrimitiveAngle.** For every continuous function
\(f:S^3\to\mathbb C\), uniformly as locally admissible \(N\to\infty\),
\[
 \frac1{R^*(N)}
 \sum_{\substack{F(v)=N\\\text{both blocks primitive}}}
             f(Tv/\sqrt N)\longrightarrow\int_{S^3}f\,d\sigma.
 \tag{23}
\]
In particular, for every fixed interval \(I\subset[0,1]\) of positive
length,
\[
 \#\{(z,w)\text{ counted by }R^*(N):Q_1(z)/N\in I\}
       =|I|JN S^*(N)+o_{Q_1,Q_2,I}(N S^*(N)).
 \tag{24}
\]
Every sufficiently large locally admissible target therefore has
representations with both summand values in \([N/3,2N/3]\).
Each fixed positive-length balance interval has
\(\gg_{Q_1,Q_2,I,\epsilon}N^{1-\epsilon}\) distinct ordered value pairs.

Here uniformity means that the displayed limit holds through **all**
locally admissible targets, without an exceptional target subsequence.
This is a qualitative angular assertion: no power-saving error for
(23)--(24), moving intervals or shrinking windows is asserted.

**Proof.** Fix \(Y\ge\max S\), with \(Y\ge2\). Let \(R_Y(N)\) count
shell points for which neither block is zero modulo any \(p\le Y\).
Put
\[
 S_Y(N)=\prod_{p\le Y}\beta_p^*(N)
                          \prod_{p>Y}\beta_p(F;N).
\]
Finite inclusion-exclusion over the full rank sublattices
\(d\mathbb Z^2\oplus e\mathbb Z^2\), where
\(d,e\mid\prod_{p\le Y}p\), and (13) give
\[
 R_Y(N)=JN S_Y(N)+O_{Q_1,Q_2,Y,\epsilon}(N^{1/2+\epsilon}).
 \tag{25}
\]
The fixed prime set is chosen before \(N\) grows. The finite sieve
excludes zero blocks automatically.

For every good prime, (17)--(18) give
\[
 0\le\beta_p(F;N)-\beta_p^*(N)\ll p^{-2}.
\]
The factors (5) are bounded below by an absolute positive constant
for \(p\ge3\). Thus, uniformly in \(N\) and its valuations,
\[
 1\le\frac{S_Y(N)}{S^*(N)}
      =\prod_{p>Y}\frac{\beta_p(F;N)}{\beta_p^*(N)}
      =1+O_{Q_1,Q_2}(Y^{-1}).
 \tag{26}
\]
All bad primes were included, so division by a zero bad density
does not occur on the admissible targets. Combining (25),(2),(3)
gives the positive tail estimate
\[
 0\le\frac{R_Y(N)-R^*(N)}{JN S^*(N)}
          \ll_{Q_1,Q_2}Y^{-1}+o_{Q_1,Q_2,Y}(1)
 \tag{27}
\]
as admissible \(N\to\infty\). This derives the needed tail from the
full count, rather than assuming a separate uniform Hecke estimate.

Now let \(H\) be the restriction to \(S^3\) of a homogeneous Euclidean
harmonic polynomial of degree \(d>0\). On a fixed full rank lattice
\(L\subset\mathbb Z^4\), the theta series
\[
 \sum_{v\in L}H(Tv)q^{F(v)}
\]
is a cusp form of weight \(2+d\). Odd degrees vanish by sign symmetry.
For the even degrees, Walling's Theorems 1 and 2 establish modularity
and cuspidality. Use her bilinear form \(B(v,v)=2F(v)\): her
\(\exp(\pi i B(v,v)\tau)\) then equals \(q^{F(v)}\).
The quadratic nebentypus need not be trivial when
\(\det A_1\det A_2\) is nonsquare.

Theorem 12 of Schulze-Pillot--Yenirce,
[arXiv:1602.01803v6](https://arxiv.org/abs/1602.01803), supplies the
all-\(n\) bound for arbitrary fixed integral weight, level and character.
After multiplying its coefficient estimate by \(N^{-d/2}\), one gets
\[
 \sum_{\substack{v\in L\\F(v)=N}}H(Tv/\sqrt N)
             \ll_{L,H,\epsilon}N^{1/2+\epsilon}.
 \tag{28}
\]
In particular no condition that \(N\) be coprime to the lattice level
is required. The finite inclusion-exclusion giving \(R_Y\) gives
(28) for its finitely sifted points as well.

Since the fully primitive points form a subset of these points,
their harmonic sum differs by at most
\(\|H\|_\infty(R_Y-R^*)\). Equations (3),(27),(28) imply
\[
 \limsup_{\substack{N\to\infty\\N\ {\rm admissible}}}
 \frac{\left|\sum_{\text{both primitive},\,F(v)=N}
                          H(Tv/\sqrt N)\right|}
                         {JN S^*(N)}
             \ll_{Q_1,Q_2,H}Y^{-1}.
\]
Let \(Y\to\infty\). All nonconstant harmonic moments vanish. Polynomial
approximation on the sphere, with the decomposition of polynomial
restrictions into spherical harmonics, proves (23).
Under uniform surface measure on \(S^3\subset\mathbb R^2\oplus\mathbb R^2\),
the squared length of the first block is uniform on \([0,1]\).
Indeed writing the two radii as \(\sqrt t,\sqrt{1-t}\) gives constant
surface density in \(t\). Interval boundaries have measure zero, so
(23) implies (24). Finally (11) removes coordinate multiplicities.
\(\square\)

The exact additional primary sources are Walling,
*The Eichler Commutation Relation for theta series with spherical
harmonics*, [Theorems 1 and 2](https://matwbn.icm.edu.pl/ksiazki/aa/aa63/aa6335.pdf),
printed pages 238 and 241, and Schulze-Pillot--Yenirce, v6, Theorem 12,
PDF page 12. Walling's
[1994 correction](https://matwbn.icm.edu.pl/ksiazki/aa/aa66/aa6647.pdf)
withdraws the corollary to Theorem 3 and Section 4; neither is used.
Exact copies and hashes are in the source ledger. This angular
deduction uses classical fixed-level inputs and makes no priority claim.

## An exact dyadic shell formula for discriminant \(-8\)

Here \(Q(x,y)=x^2+2y^2\). Let \(r(N)\) count all signed ordered
coordinate representations by \(Q\oplus Q\), and let \(G_2(N)\) count
those for which neither binary pair has both coordinates even.
**Only primitivity at 2 is imposed in \(G_2\)**; odd inert prime
obstructions remain to be removed.

Put \(\sigma(x)=0\) when \(x\) is not a positive integer, and define
\[
 f_{32}(\tau)=\eta(4\tau)^2\eta(8\tau)^2
       =\sum_{N\ge1}a_Nq^N
       =q\prod_{j\ge1}(1-q^{4j})^2(1-q^{8j})^2 .
\]
Then, for every \(N\ge1\),
\[
 \begin{split}
 r(N)&=4\sigma(N)-4\sigma(N/2)+8\sigma(N/4)-32\sigma(N/8),\\
 s(N)&:=r_{Q\oplus4Q}(N)
   =\sigma(N)-\sigma(N/2)+8\sigma(N/16)-32\sigma(N/32)+a_N,\\
 G_2(N)&=r(N)-2s(N)+r(N/4).
 \end{split}
 \tag{29}
\]
In the last term \(r(N/4)\) is zero if \(4\nmid N\).
The last equality is ordinary inclusion-exclusion of the two
coordinate-divisibility conditions.

Consequently, if \(N=2^e m\) with \(m\) odd,
\[
 G_2(N)=
 \begin{cases}
  2\sigma(m)-2a_m,&e=0,\\
  4\sigma(m),&e=1,\\
  20\sigma(m),&e=2,\\
  16\sigma(m),&e=3,\\
  0,&e\ge4 .
 \end{cases}
 \tag{30}
\]
Moreover \(a_N=0\) for even \(N\), and
\[
 |a_N|\le\tau(N)\sqrt N.
 \tag{31}
\]
Thus the dyadic part of the \(-8\) sieve has an exact coefficient
formula and an explicit square-root error on odd targets. Formula
(30) is not the full primitive norm sumset classification.

**Modular proof and finite certificate.** The two theta products in
(29) lie respectively in \(M_2(\Gamma_0(8))\) and
\(M_2(\Gamma_0(32))\), with trivial character: their integral even
Gram matrices have levels 8 and 32 and square determinants 64 and
1024. This is the degree-zero theta transformation law, with the same
Gram convention as above.

Write \(E_2(\tau)=1-24\sum_{n\ge1}\sigma(n)q^n\).
The proposed Eisenstein series are
\[
 E_r=-\frac1{24}\{4E_2(\tau)-4E_2(2\tau)
                         +8E_2(4\tau)-32E_2(8\tau)\},
\]
\[
 E_s=-\frac1{24}\{E_2(\tau)-E_2(2\tau)
                         +8E_2(16\tau)-32E_2(32\tau)\}.
\]
They are holomorphic modular forms at the stated levels. Indeed
\(dE_2(d\tau)-E_2(\tau)\) is holomorphic of weight two on
\(\Gamma_0(d)\), as follows directly from the transformation law of
\(E_2\); each displayed combination is a linear combination of these
forms, since its coefficients \(c_d\) satisfy \(\sum c_d/d=0\).
Both displayed constant terms equal 1.

Martin--Ono, *Eta-quotients and elliptic curves*, Theorem 1,
level-32 row, identifies \(f_{32}\) as a normalized weight-two newform
of level 32; Theorems 2 and 3 identify its Gaussian CM curve and
prime coefficients. The finite Sturm bounds are
\[
 \frac2{12}[\mathrm{SL}_2(\mathbb Z):\Gamma_0(8)]=2,\qquad
 \frac2{12}[\mathrm{SL}_2(\mathbb Z):\Gamma_0(32)]=8 .
\]
The coefficient checks, including the constant terms, are exactly
\[
 \begin{array}{c|rrrrrrrrr}
 N&0&1&2&3&4&5&6&7&8\\ \hline
 [q^N]\Theta_{Q\oplus Q}&1&4&8&&&&&&\\
 [q^N]E_r&1&4&8&&&&&&\\
 [q^N]\Theta_{Q\oplus4Q}&1&2&2&4&4&4&8&8&8\\
 [q^N](E_s+f_{32})&1&2&2&4&4&4&8&8&8 .
 \end{array}
\]
These integer equalities prove the two identities by Sturm's
Theorem 1, applied modulo every prime: the difference has all
coefficients divisible by every rational prime, hence is zero.
This is a finite proof certificate with a stated completeness bound,
not an extrapolation of a target search.

For (31), Martin--Ono's explicit CM formula gives
\(|a_p|\le2\sqrt p\) at every odd prime, including \(a_p=0\) for
\(p\equiv3\pmod4\). The roots of
\(1-a_pX+pX^2\) therefore have absolute value \(\sqrt p\).
The newform Euler recurrence and multiplicativity give (31) at every
odd \(N\). The eta product is supported on \(N\equiv1\pmod4\), so
the even coefficients vanish.

The checker verifies the full finite Sturm coefficient lists by
independent lattice enumeration, divisor sums and multiplication of
truncated eta factors, and also checks (29)--(30) against direct
dyadic coordinate enumeration through its recorded bound.
The source ledger retains Martin--Ono's author PDF and the two-page
primary Sturm preview containing Theorem 1 (printed page 276);
only that theorem statement, not the unavailable remainder of its
proof, is imported here. The ordinary \(r(N)\) formula is classical;
no priority claim is made for this explicit dyadic calculation.

## Explicit \(-8\) sieve cutoffs and a finite exception reduction

The dyadic formula can be combined with the finite Hecke-quotient
argument from the Eisenstein proof. The order and the normalization
must be changed explicitly. The result here is an effective cutoff
and a finite reduction; it does not assert that the remaining finite
exception set is empty.

### The quaternion order and the unramified operator

Let \(H\) be the Hurwitz order in Hamilton's quaternion algebra, with
standard generators \(I,J,K\). Put
\[
 \iota=I+J,\quad j=K,\qquad
 O=\mathbb Z[\iota]\oplus\mathbb Z[\iota]j
   =\{a+bI+cJ+dK:a,b,c,d\in\mathbb Z,\ b\equiv c\bmod2\}.
\]
Thus \(\iota^2=-2\), \(j^2=-1\), \(j\iota=-\iota j\), and the norm on
\(O\) is \(Q\oplus Q\). The order has index 4 in \(H\).

Here is a direct verification that \(O\) has one locally principal
right-ideal class. The Hurwitz order is norm Euclidean. After rounding
a real quaternion to the nearest integer vector, write the four
coordinate errors as \(r_i\), \(|r_i|\le1/2\), and put \(s=\sum r_i^2\).
The squared distance to a suitable half-integer vector is
\(s-\sum|r_i|+1\). Since \(s\le\frac12\sum|r_i|\), one of these two
distances is at most \(1/2\). Euclidean division makes every Hurwitz
right ideal principal.

Away from 2, \(O_\ell=H_\ell\). Normalize the additive volume of \(H_2\)
to one. Its integer and half-integer cosets each have volume \(1/2\).
Every element of the half-integer coset has odd norm, as do half the
integer elements. Hence \(\operatorname{vol}(H_2^\times)=3/4\).
The additive volume of \(O_2\) is \(1/4\), and half its elements have
odd norm, so \(\operatorname{vol}(O_2^\times)=1/8\). The unit index is
therefore 6. The 24 Hurwitz units have stabilizer
\(H^\times\cap O_2^\times=O^\times=\{\pm1,\pm K\}\), of size 4.
Their orbit already has all six cosets of
\(H_2^\times/O_2^\times\). Consequently
\[
 H^\times\backslash H_2^\times/O_2^\times
\]
has one element. Principalizing in \(H\), then using this unit
transitivity at 2, principalizes every locally principal \(O\)-ideal.
This is a double-coset argument, not an inference from a mass alone.

For an odd inert prime \(p\), \(({-8}/p)=-1\), reduction identifies
\(O/pO\) with \(M_2(\mathbb F_p)\). Set
\[
 G_p=\mathrm{PGL}_2(\mathbb F_p),\quad
 \Gamma_p=\rho_p(O^\times),\quad
 Y_p=\Gamma_p\backslash G_p .
\]
Then \(|\Gamma_p|=2\). For \(p\nmid n\), use the concrete operator
\[
 (T_nf)(\Gamma_p g)=\frac14
       \sum_{\substack{\alpha\in O\\\operatorname{nrd}\alpha=n}}
                           f(\Gamma_p\rho_p(\alpha)g).
 \tag{32}
\]
For odd \(m\) with \(p\nmid m\), class number one identifies this with
the all-integral-right-ideal correspondence of norm \(m\).
It has constant eigenvalue \(\sigma(m)\); its determinant-character
eigenvalue is \((m/p)\sigma(m)\). On the orthogonal complement of those
two characters,
\[
 \|T_m\|\le\tau(m)\sqrt m .
 \tag{33}
\]
The proof of the
[Eisenstein spectral lemma](../../unpublished/eisenstein-restricted-prime-sums/sections/prime-bounds.tex)
applies with the following checked replacements: the definite algebra
is ramified at \(2,\infty\), the fixed local level is \(O_2^\times\),
the unit count is 4, and the prime being reduced is any inert \(p\ge5\).
Class number one gives the same finite adelic quotient \(Y_p\).
The only character constituents on \(\mathrm{PGL}_2(\mathbb F_p)\)
are the trivial and quadratic determinant characters. All other
constituents transfer by the same Jacquet--Langlands theorem to
weight-two cusp forms. Every prime dividing \(m\) is a split,
maximal-level place, so its two Hecke roots have size \(\sqrt\ell\).
The all-ideal prime-power recurrence then gives (33), uniformly in
\(p,m\). No unramified estimate at 2 is being assumed.

For \(N=2^e m\), \(m\) odd, \(0\le e\le3\), let
\(A_e=T_{2^e}\) as defined by (32). Its row sum is
\[
 c_e=r(2^e)/4\in\{1,2,6,6\}.
\]
It has nonnegative entries, is self-adjoint by quaternion conjugation,
and has norm at most \(c_e\) by the row-sum bound. Coprime ideal
factorization gives
\[
 T_N=A_eT_m=T_mA_e .
 \tag{34}
\]
For clarity, if \(q\in O\) has norm \(N\), the unique norm-\(m\)
superideal of \(qO\) has localization \(O_2\) at 2 and \(qO_\ell\)
at every odd \(\ell\). It is integral and locally principal.
Write it as \(\beta O\), with \(\operatorname{nrd}\beta=m\).
Then \(q=\beta\alpha\), \(\operatorname{nrd}\alpha=2^e\).
Exactly four choices differ by an \(O\)-unit, giving the factor \(1/4\)
in (34). Reversing the two coprime local factors gives commutation.
Thus the complementary norm of \(T_N\) is at most
\(c_e\tau(m)\sqrt m\).

### The actual sieve inequality

Put
\[
 \mathcal I=\{p\text{ prime}:p\equiv5,7\bmod8\},\qquad
 E(m)=\frac{\tau(m)\sqrt m}{\sigma(m)},\qquad
 M(X)=\max_{1\le t\le X}
       \prod_{\substack{q\mid t\\q\equiv1,3(8)}}(v_q(t)+1).
\]
Only primes in \(\mathcal I\) and the dyadic condition need to be
removed to obtain **primitive norm values**. Split prime factors may
remain in a coordinate gcd: their norm values admit separate
primitive binary representations.

The image of \((\mathbb Z[\iota]/p)^\times/\mathbb F_p^\times\) is a
nonsplit torus, with normalizer of order \(2(p+1)\). It contains
\(\Gamma_p\), and an invertible quaternion has a zero binary component
precisely when its projective image lies in this normalizer.
Its orbit fraction is
\(\mu_p=2/(p(p-1))\). Averaging over this orbit exactly as in the cited
spectral proof, and using (33)--(34), shows that, for \(p\nmid N\),
\[
 \frac{b_p(N)}{r(N)}
       \le\mu_p+E(m)(1-\mu_p),
 \tag{35}
\]
where \(b_p\) counts all shell points with a binary component divisible
by \(p\). If \(p\mid N\), anisotropy of \(Q\bmod p\) forces both
components to be divisible by \(p\); the spoiled count is zero unless
\(p^2\mid N\), and otherwise its ratio is
\(\sigma(m/p^2)/\sigma(m)<p^{-2}\le\mu_p\).

For a large prime retain only the \(G_2\) shell. Both its components
are nonzero. Unique ideal factorization in \(\mathbb Z[\sqrt{-2}]\)
gives \(r_Q(t)\le2M(X)\) for \(1\le t\le X\). Also
\[
 \#\{(x,y)\ne(0,0):x^2+2y^2\le X\}\le4X .
\]
To see this, let \(g(y)=\sqrt{\max(X-2y^2,0)}\). It is even and
decreasing on \(y\ge0\), so
\(\sum_{y\in\mathbb Z}g(y)\le g(0)+\int_{\mathbb R}g(y)\,dy\).
Counting rows and deleting the origin therefore gives
\[
 S(X)\le2\sum_y g(y)+2\lfloor\sqrt{X/2}\rfloor
      \le\frac{\pi X}{\sqrt2}+(2+\sqrt2)\sqrt X\le4X\quad(X\ge4).
\]
For the last inequality use
\(\pi/\sqrt2+(2+\sqrt2)/2<110/49+12/7=194/49<4\).
The three smaller integer ranges have respectively 2, 4 and 8 points,
so the bound holds for every real \(X\ge1\).
Choosing the divisible component first, and counting both
orientations, gives
\[
 \frac{b_p^{(2)}(N)}{r(N)}
       \le L_e\,\frac{M(N)}{p^2},\qquad
       L_e=\frac{4\,2^e}{c_e}\le\frac{16}{3}.
 \tag{36}
\]
The exclusion of zero components is essential here.

For every real \(P>0\), write
\[
 B=\sum_{p\in\mathcal I}\mu_p,\qquad
 T(P)=\sum_{\substack{k>P\\k\equiv5,7(8)}}k^{-2}.
\]
Elementary progression bounds give
\[
 \#\{p\in\mathcal I:p\le P\}\le P/4+1,\qquad
 T(P)\le\frac1{4P}+\frac2{P^2},\qquad B<\frac3{16}.
 \tag{37}
\]
For the last bound sum the inert primes through 64 exactly, then
majorize the rest by
\[
 \frac52\left(\frac1{69^2}+\frac1{8\cdot69}
                       +\frac1{71^2}+\frac1{8\cdot71}\right).
\]
The resulting rational upper bound is
\(49287299505759574293/269836951018231737140<3/16\).

By (30)--(31), \(G_2(N)/r(N)\ge1/2-E(m)/2\) in all four strata.
A sufficient, fully explicit criterion for a primitive **value**
decomposition is therefore
\[
 B+E(m)(P/4+3/2)+L_eM(N)
                           \left(\frac1{4P}+\frac2{P^2}\right)<\frac12.
 \tag{38}
\]

### Certified constants, global onset and the smaller residual set

The finite Euler-factor certificate in the checker proves
\[
 \begin{split}
 M(X)&\le6X^{1/6},\\
 \tau(m)&\le7m^{1/6}
       &&\text{if every prime of odd }m\text{ is in }\mathcal I,\\
 \tau(m)&\le42m^{1/6}
       &&\text{for every odd }m,\\
 M(X)&\le2X^{1/4}.
 \end{split}
 \tag{39}
\]
Here is the completeness argument for these finite checks.
For \(k=6\), maximize \((a+1)^6/p^a\) for each odd prime \(p<64\)
and \(0\le a\le6\). For \(a\ge6\) this expression decreases; for
\(p\ge64\) it is at most one because \(a+1\le2^a\).
The split-prime product is below \(6^6\), and the inert-prime product
is below \(7^6\). Their product gives \(42^6\) for all odd primes.
The \(k=4\) split calculation uses only primes 3 and 11:
the product is \(4096/297<2^4\). These are finite exact rational
certificates for (39), with no bound on \(X\) or \(m\).

Since \(N\le8m\), (39) gives \(M(N)\le9m^{1/6}\).
For arbitrary odd \(m\), take \(P=m^{1/4}\) in (38).
The terms after \(B\) are at most
\[
 \frac{45}{2}m^{-1/12}+159m^{-1/3}.
\]
If \(m\ge2^{75}\), then \(m^{1/12}>224/3\) and \(m^{1/3}\ge2^{25}\);
\[
 \frac3{16}+\frac{135}{448}+\frac{159}{2^{25}}<\frac12 .
\]
Thus every \(N\ge2^{78}\) with \(16\nmid N\) is a sum of two primitive
\(-8\) norm values. This is an unconditional explicit global onset,
not a claim about the exact smaller exceptions.

The same inequalities give a fully numerical lower bound for distinct
ordered value pairs. The surviving shell proportion exceeds `1/100`,
since

```
1/2-3/16-135/448-159/2^25 = 2620327/234881024 > 1/100.
```

For every admissible dyadic stratum, `r(N)>=3N`. Each value pair has
at most `(2M(N))^2<=144N^(1/3)` coordinate witnesses. Hence for every
`N>=2^78` with `16` not dividing `N`,

```
#{(a,b) in P_-8^2 : a+b=N} >= N^(2/3)/4800.
```

This explicit lower bound has a numerical onset. It does not identify
the sharp logarithmic scale or an unweighted leading constant.

There is a sharper reduction for the exact classification. If \(m\)
has only inert prime factors, use \(P=2m^{1/4}\) and the second bound
in (39). The terms after \(B\) are at most
\[
 \frac{19}{2}m^{-1/12}+\frac{69}{2}m^{-1/3}.
\]
They prove existence for all \(0\le e\le3\) once \(m\ge2^{60}\), since
\[
 \frac3{16}+\frac{19}{64}+\frac{69}{2^{21}}<\frac12 .
\]
For an odd prime \(N\), use \(M(N)\le2N^{1/4}\),
\(E(N)\le2N^{-1/2}\), \(e=0\), \(L_0=4\), and \(P=2N^{3/8}\).
The corresponding terms are at most
\(2N^{-1/8}+7N^{-1/2}\).
All odd primes \(N\ge2^{22}\) are therefore represented, because
\[
 \frac3{16}+\frac3{10}+\frac7{2048}<\frac12 .
\]

The companion checker exhausts all 295947 primes \(N\le2^{22}\),
with no exception. It constructs the unramified core by deleting
multiples of 2 and every prime \(5,7\bmod8\), then adds twice that
core. The maximal-order primitive criterion proves that these are
exactly the primitive norm values; direct coordinate enumeration
independently matches the flags through 10000. A packed bitset then
forms the union \(a+P_{-8}\) over the 1064 primitive values
\(a\le5000\). Every prime through the cutoff belongs to this union;
no fallback search was needed. The retained receipt contains the
exact range, target count and flag/coverage hashes.

Together with the analytic prime threshold, this proves the
unrestricted statement
\[
 \boxed{\text{Every prime }N\ge2\text{ belongs to }P_{-8}+P_{-8}.}
 \tag{40a}
\]
The cutoff is what makes the finite certificate a complete proof of
this prime-family statement.

Finally multiplication of a primitive \(-8\) norm value by a split
prime preserves membership in that norm set. A least counterexample
to the proposed exact identity, if divisible by a split prime \(q\)
with \(N/q\ge2\), would descend to a smaller counterexample.
The boundary \(N=q\) is ruled out by (40a). Consequently the exact
identity is equivalent to coverage of the **finite** residual family
\[
 N=2^e m\ge2,\quad 0\le e\le3,\quad m<2^{60},\qquad
       q\mid m\Longrightarrow q\equiv5,7\bmod8 .
 \tag{40}
\]
This family has not been exhausted here. A least counterexample, if
one exists, is therefore smaller than \(2^{63}\) and has no odd
split prime factor. Equation (40) is a least-counterexample reduction;
the reduced core cutoff \(2^{60}\) alone is **not** a global onset
for arbitrary targets with split prime factors. The unconditional
all-target onset proved above is \(2^{78}\).

## Consequences and limits

* The local-global existence conclusion applies to mixed quadratic
  fields, nonprincipal classes and nonmaximal quadratic orders.
  Ramified and conductor primes enter the finite local test (4).
  They do not need to be fitted into the coefficient hypotheses of
  the splitting-semigroup theorem.
* The theorem permits every target valuation compatible with the
  primitive local conditions. It is stronger than applying a theorem
  about one primitive four-vector: both blocks are sieved separately.
* Distinct value pairs have at least \(N^{1-\epsilon}\) representations
  on every sufficiently large admissible target, using (11) to remove
  coordinate multiplicities from (2). This does not identify their
  logarithmic order or leading constant.
* The exact \(-8\) identity is reduced to the explicit finite family
  (40), with an unconditional global onset \(2^{78}\). The residual
  inert-core family has not been exhausted. No analogous numerical
  cutoff for \(-11\) is computed here.
* Simultaneous Gaussian--Eisenstein primitive norms in *each* summand
  impose two independent norm conditions on each value. Their
  \(C_{12}+C_{12}\) problem is not a single binary-form sum and does
  not follow from this theorem.
* Fixed macroscopic balance and qualitative angular equidistribution
  follow from (23)--(28). No power-saving angular discrepancy, moving
  intervals or shrinking-window result is proved here. Unweighted
  value-pair leading constants remain a separate obligation.
* The imported general coefficient bounds are effective, and the
  other estimates above admit effective constants for fixed forms.
  For general forms this gives an effective route in principle.
  The separate \(-8\) sieve supplies the numerical onset \(2^{78}\)
  and the smaller finite reduction (40). No novelty claim is made
  before a priority comparison with work on coprimality restrictions
  in quaternary representations.

For comparison, Jörg Brüdern's
[*Sums of Squares Coprime in Pairs*](https://doi.org/10.4064/ba62-3-3),
2014, Theorem 2 and Sections 2--4, already uses a truncated coprimality
Möbius expansion, a large-common-factor tail, and a completed singular
series for a stronger collection of pairwise conditions on four squares.
The present proof uses the easier two-block pattern and general
quaternary modular bounds. This is evidence that the method is classical,
not a priority claim for the exact statement here.

The companion checker
[verify_primitive_norm_quaternary.py](verify_primitive_norm_quaternary.py)
verifies the good-prime counts, bad-prime finite-lifting criterion on
sample forms, and the exact Möbius identity on finite ranges. These
checks challenge normalizations and examples; they do not prove the
uniform analytic estimates or the unrestricted theorem.
