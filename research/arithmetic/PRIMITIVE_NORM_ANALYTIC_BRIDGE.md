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

**Theorem.** For every fixed \(0<\eta<1/8\), in particular
\(\eta=1/9\),
\[
 R^*(N)=JN\prod_p\beta_p^*(N)+O_{Q_1,Q_2,\eta}(N^{1-\eta}).
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

Dependence of the onset on the forms is necessary. For every integer
\(r\ge1\), the primitive principal norm form of the order
\(\mathbb Z[17^r i]\) and the target
\[
 Q_r(x,y)=x^2+17^{2r}y^2,\qquad N_r=17^{2r}-2
\]
give a locally pair-primitive exception. Indeed a positive primitive
value below \(17^{2r}\) must have \(y=0\), hence \(x=\pm1\), so
the only such value is 1, and \(N_r>2\) is not a sum of two positive
primitive values.

For completeness, put \(c=17^r\). At 17 take the two vectors
\((0,1)\) and \((t,0)\), where \(t^2=-2\) in \(\mathbb Z_{17}\);
the root \(t=7\bmod17\) lifts by Hensel. At two take
\((1,c^{-1})\) and \((s,2c^{-1})\), where
\(s^2=c^2-8\) in \(\mathbb Z_2\). This odd square root exists because
\(c^2-8\equiv1\bmod8\). These pairs are primitive and their norms
sum to \(N_r\). At each odd \(q\ne17\), every unit is represented
primitively by \(Q_r\) over \(\mathbb Z_q\): after scaling \(y\) by
the unit \(c\), the equation \(x^2+y^2=a\ne0\bmod q\) has
\(q-\left(\frac{-1}{q}\right)>0\) nonsingular solutions and lifts.
Choose a residue \(a\) different from both 0 and \(N_r\bmod q\);
then both \(a\) and \(N_r-a\) have primitive local representations.
Thus \(N_r\to\infty\) gives an infinite obstruction to a cutoff
independent of the order conductor. The fixed-form quantifiers in
(2)--(4) retain precisely this necessary dependence.

## Additive disappearance of class obstructions

Here is the class-group consequence with no hypothesis on the class
number. It concerns primitive elements in the rational-coordinate sense:
an element of an order is primitive if it is outside
\(m\mathcal O\) for every integer \(m>1\).

**Corollary AdditiveClassDisappearance.** Fix any imaginary quadratic
order \(\mathcal O_D\) of discriminant \(D<0\). In an integral basis
\(1,\omega_D\), put
\[
 Q_D(x,y)=\operatorname{Norm}(x+y\omega_D),\qquad
 P_D=\{Q_D(x,y)>0:\gcd(x,y)=1\}.
\]
Let \(\mathcal A_D\) be the set of integers \(N\ge2\) for which
\(Q_D(z)+Q_D(w)=N\) has a solution over every \(\mathbb Z_p\), with
both blocks outside \(p\mathbb Z_p^2\). Then there is a finite set
\(E_D\subset\mathcal A_D\) such that
\[
                   P_D+P_D=\mathcal A_D\setminus E_D.       \tag{CG.1}
\]
The set \(\mathcal A_D\) is a union of residue classes, restricted to
\(N\ge2\), modulo the explicit sufficient modulus
\[
             M_D=\prod_{p\mid2D}p^{\,2v_p(|D|)+1}.          \tag{CG.2}
\]
Each class is decided by finite enumeration of the pair-primitive
congruences at the indicated prime powers. Thus ramification and the
conductor are included in the criterion.

For maximal orders, [PL3All](PRIMITIVE_NORM_LOCAL_STRUCTURE.md#the-complete-maximal-order-self-sum-criterion)
gives the complete local set by a smaller closed congruence table,
including the special residue condition at a ramified three.

More precisely, the number \(R_D^*(N)\) of ordered coordinate
representations with both pairs primitive satisfies, for every fixed
\(\epsilon>0\),
\[
 R_D^*(N)=\frac{4\pi^2}{|D|}N\prod_p\beta_{D,p}^*(N)
                       +O_{D,\epsilon}(N^{7/8+\epsilon}),  \tag{CG.3}
\]
and \(R_D^*(N)\gg_D N\) for every sufficiently large
\(N\in\mathcal A_D\). In particular, for every fixed \(\epsilon>0\)
there are \(\gg_{D,\epsilon}N^{1-\epsilon}\) distinct ordered value
pairs on these targets. These are lower bounds; the coordinate
asymptotic is not an unweighted value-pair asymptotic. No numerical
cutoff or list of \(E_D\) is asserted by this general corollary.

**Corollary GenusInvariantConvolution.** Fix primitive positive-definite
integral binary forms \(Q_1,Q_2,Q'_1,Q'_2\), with \(Q_i\) and \(Q'_i\)
in the same genus for \(i=1,2\). Then
\[
 R_{Q_1,Q_2}^*(N)-R_{Q'_1,Q'_2}^*(N)
                    =O_{\mathbf Q,\epsilon}(N^{7/8+\epsilon}). \tag{CG.4}
\]
Their local admissible target sets are identical, and their primitive
value sumsets differ by at most finitely many integers. For any fixed
two genera, a single eventual threshold can be chosen for all pairs
of form classes in those genera.

**Proof.** A quadratic order has a basis giving
\[
 Q_D(x,y)=
 \begin{cases}
 x^2+xy+(1-D)y^2/4,&D\equiv1\pmod4,\\
 x^2-(D/4)y^2,&D\equiv0\pmod4.
 \end{cases}
\]
Both are primitive positive integral forms, and their Gram determinant
is \(|D|\). Primitivity of the element is exactly primitivity of
its two coordinates; a change of integral basis preserves it.
The theorem (2), with its Waibel error proved below, gives (CG.3).
The local criterion (4) gives (CG.2) and exact nonrepresentation
outside \(\mathcal A_D\). For the self-sum all good-prime target
ratios in (20b) exceed one, while admissible bad-prime densities have
a uniform positive lower bound. Hence the singular series is bounded
below by a positive constant depending only on \(D\). Taking, for
example, \(\epsilon<1/8\) in (CG.3) proves the linear lower bound
and the finiteness in (CG.1). The binary multiplicity bound (11),
applied with exponent \(\epsilon/2\) in each block, converts this to
the stated distinct-value lower bound.

For (CG.4), local integral equivalence gives matrices
\(U_{i,p}\in\mathrm{GL}_2(\mathbb Z_p)\) carrying \(Q_i\) to \(Q'_i\).
Their reductions are bijections at every prime power and preserve
nonzero blocks modulo \(p\). Thus every local density, including the
bad-prime densities, agrees exactly. The Gram determinants agree as
well. Subtract the two instances of (2); their entire main terms
cancel. The common local criterion and eventual existence give the
finite symmetric difference of sumsets. There are finitely many
classes in each fixed positive-definite genus, so taking the maximum
of their finitely many thresholds proves the last assertion.
\(\square\)

There is also a direct proof of (CG.4) before singular-series
completion. For every Möbius divisor pair \(d,e\), the scaled
quaternary forms \(d^2Q_1\oplus e^2Q_2\) and
\(d^2Q'_1\oplus e^2Q'_2\) remain locally isometric, so their genus
coefficients cancel exactly. Truncate the two Möbius sums at
\(d,e\le H\). The large-gcd tails cost \(O(N^{1+\epsilon}/H)\),
the removed zero blocks \(O(H^2N^\epsilon)\), and the total cusp
difference costs
\(O(N^{1/2+\epsilon}H^{3+\epsilon}
      +N^{3/4+\epsilon}H^{1+\epsilon})\).
Taking \(H=N^{1/8}\) proves (CG.4) directly. Thus genus invariance
does not require a separate class-group equidistribution estimate.

This proves that finer ideal-class obstructions within a genus do
not persist in these additive existence problems, and do not change
the leading coordinate count. Genus restrictions remain in the
local factors. The conclusion allows arbitrary class groups, mixed
form classes and nonmaximal orders; it does not require converting
principal norm values into a prime-factor semigroup. It is a
corollary of the established quaternary deduction, rather than a new
distribution theorem for class-group characters. The precise
classical attribution and limits of the literature search are
recorded in [the novelty comparison](PRIMITIVE_NORM_NOVELTY.md).

For a concrete class obstruction, the
[discriminant -23 criterion](PRIMITIVE_NORM_CLASS_NUMBER_THREE.md)
shows that \(2^{3a+1}\), \(a\ge0\), is locally a primitive norm
everywhere but is never a global primitive principal norm: the
prime above \(2\) has ideal class of order three, and a primitive
power must use a single orientation. Nevertheless every
sufficiently large member of this sequence is a sum of two
primitive principal norms. Its only additive local condition is
\(23^2\nmid N\), by PL3All, and powers of two always satisfy it.

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
Waibel's genus formula, equations (1)--(2) of
[arXiv:2012.04311v2](https://arxiv.org/pdf/2012.04311v2), independently
displays this same factor in this Gram convention.
The real factor also follows by differentiating the ellipsoid volume
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
Schulze-Pillot--Yenirce, Theorem 10 in the edition he cites. In the retained
[arXiv:1602.01803v6](https://arxiv.org/pdf/1602.01803v6), it is
Theorem 12, PDF page 12. It covers every positive coefficient index,
arbitrary level and arbitrary character, including oldforms. Its bound is
\(\ll_\epsilon\langle C,C\rangle^{1/2}L^{1+\epsilon}
\tau(n)\sqrt n\), using \(\dim S_2(\Gamma_0(L),\chi)\ll_\epsilon
L^{1+\epsilon}\). Rouse's inner product is \(3/\pi\) times that
source's normalized inner product, a fixed factor. Combining this bound
with Rouse's Theorem 3 proves (8), with constants depending only on
\(\epsilon\). Local solubility and coprimality of \(n\) with the level
are not hypotheses of this coefficient estimate. This deduction uses
the coefficient estimate, rather than the existential conclusion of
Rouse's Theorem 1.

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
and retaining the power \(3+\epsilon\) gives
\[
 R^*(N)=\sum_{d,e\le D}\mu(d)\mu(e)E_{d,e}(N)
 +O_{Q_1,Q_2,\epsilon}\left(
 N^{1+\epsilon}/D+D^{8+2\epsilon}N^{1/2+\epsilon}
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
Taking any fixed \(0<\alpha<1\), (16)--(18) and multiplicativity imply
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
If \(p\nmid N\), then \(\lambda_p(1,1;N)=0\) exactly, so this is
\(1+O(p^{-2+\alpha})\). The product over these primes converges absolutely
for every \(\alpha<1\). At target primes it is at most
\[
 1+C_\alpha(p^{-1}+p^{-2+\alpha}+p^{-2+2\alpha})
 \le 1+C'_\alpha p^{-\delta},
 \qquad \delta=\min(1,2-2\alpha)>0.
\]
For any \(\epsilon>0\), choose a fixed prime threshold after which
\(1+C'_\alpha p^{-\delta}\le p^\epsilon\). Their product over \(p\mid N\)
is at most a fixed constant times \(\operatorname{rad}(N)^\epsilon\),
hence \(O_{\alpha,\epsilon}(N^\epsilon)\). This proves uniformity even
when \(v_p(N)\) is arbitrarily large. For \(\alpha>1/2\) the extra
target-prime term need not be \(O(1/p)\); only its product needs to be
\(N^{o(1)}\). Bad primes are a fixed finite set. The argument still
applies if some unrestricted bad density vanishes, since then every
term vanishes; it does not divide by a zero density. The endpoint
\(\alpha=1\) is not claimed: the positive nontarget Euler factors
can then have nonsummable \(1/p\) terms.

It follows from (19) that the main term in (14) can be completed with
error \(O(N^{1+\epsilon}D^{-\alpha})\). Local inclusion-exclusion gives
\[
 \sum_{d,e\ge1}\mu(d)\mu(e)E_{d,e}(N)
 =JN\prod_p(\lambda_p(0,0)-\lambda_p(1,0)
                  -\lambda_p(0,1)+\lambda_p(1,1))
 =JN\prod_p\beta_p^*(N).
 \tag{20}
\]
The resulting error budget is
\[
 O_{Q_1,Q_2,\alpha,\epsilon}\left(
 N^{1+\epsilon}D^{-1}
 +N^{1/2+\epsilon}D^{8+2\epsilon}
 +D^2N^\epsilon+N^{1+\epsilon}D^{-\alpha}\right).
 \tag{20a}
\]
For fixed \(\alpha<1\), take \(D=N^{1/(16+2\alpha)}\).
Before the arbitrarily small \(\epsilon\) losses, the limiting saving
is \(\alpha/(16+2\alpha)\). This Rouse-based route proves (2)
for every fixed \(\eta<1/18\), including \(\eta=1/20\) by
\(\alpha=15/16\), \(D=N^{1/18}\), \(\epsilon=1/1000\).
A stronger existing theta estimate gives the theorem's stated range.

## Sharper power saving from Waibel

Use Waibel, *Uniform bounds for norms of theta series and arithmetic
applications*, [arXiv:2012.04311v2](https://arxiv.org/pdf/2012.04311v2),
Theorem 1 and definition (6), PDF p.3, and Lemma 15, PDF p.16.
The paper was subsequently published in 2022. The imported edition
here is the inspected v2. For a primitive positive quaternary form of
Gram matrix \(A\), level \(L\), and determinant \(\Delta\), its
all-index estimate is
\[
 |c_A(n)|\ll_\epsilon n^{1/2}
 \left(\frac{L}{\sqrt{\mathcal F(A,2)}}+
                \frac{\sqrt L}{\Delta^{1/8}}\right)
 \min\left(\sqrt L,\
       1+\frac{n^{1/4}(n,L)^{1/4}}{\sqrt L}\right)(nL)^\epsilon.
\]
Here \(\mathcal F\) is Waibel's local Jordan invariant, not the
quadratic form \(F\) of this note. This estimate requires neither
coprimality nor local solubility of \(n\).
Waibel uses the unnormalized Petersson norm; the displayed coefficient
estimate already incorporates that convention.

For the dyadic polynomial normal form in the cited v2 equation (5),
allow a scalar coefficient exponent \(\nu_i(2)=0\). The printed
condition \(\nu_i(2)\ge1\) cannot hold for four squares: its polar form
is zero modulo two, forcing all binary cross-term blocks to have even
scale, whereas the polynomial represents 1. With scalar exponents
allowed to be zero, definition (6) gives
\(\mathcal F(2I_4,2)=4\) and \(\mathcal F(2I_4,4)=16=\det(2I_4)\).
Only comparison up to fixed-form constants is used below; this
dyadic convention does not alter the claimed power saving.

Only squarefree \(d,e\) occur. Put \(g=(d,e)\),
\(u=d/g\), \(v=e/g\), \(t=uv\), and \(n=N/g^2\), discarding
\(g^2\nmid N\). For \(A=u^2A_1\oplus v^2A_2\),
\[
 \Delta=\Delta_0t^4,\qquad L\asymp_{Q_1,Q_2}t^2,\qquad
 \mathcal F(A,2)\asymp_{Q_1,Q_2}t^4.
\]
Indeed outside the fixed bad primes a prime of \(t\) has Jordan
exponents \(0,0,2,2\), giving level \(p^2\) and
\(\mathcal F\)-factor \(p^4\). At bad primes the squarefree coprime
pair \((u,v)\) has just three scaling patterns; their factors differ
from these expressions by fixed constants. This proves both the
upper and lower level comparisons needed in the minimum above.

Waibel's estimate, using the second entry of that minimum, therefore
gives
\[
 |c_{d,e}(N)|\ll_{Q_1,Q_2,\epsilon}N^\epsilon(de)^\epsilon
 \left\{N^{1/2}(de)^{1/2}
   +N^{3/4}(de)^{-1/2}(N,L_0(de)^2)^{1/4}\right\}.
\]
The discarded powers of \(g\) are \(g^{-2}\) and \(g^{-1/2}\), so both are favorable.
The target gcd can be averaged, rather than maximized. For squarefree
\(d\), set \(G_N(d)=(N,d^2)^{1/4}\). Expanding its prime factors gives
\[
 \sum_{\substack{d\le D\\d\ {\rm squarefree}}}
       d^{-1/2}G_N(d)
 \le 2D^{1/2}
       \prod_{p\mid N}\left(1+
           \frac{p^{\min(v_p(N),2)/4}-1}{p}\right)
 \ll_\epsilon D^{1/2}N^\epsilon.
\]
The expansion uses
\(\sum_{d\le D,\ k\mid d}d^{-1/2}\le2D^{1/2}/k\);
the last target-prime product is \(N^{o(1)}\).
Also
\[
 (N,L_0(de)^2)^{1/4}
 \le L_0^{1/4}(N,(de)^2)^{1/4}
 \le L_0^{1/4}G_N(d)G_N(e),
\]
including common primes.
Consequently the total truncated cusp error is
\[
 O_{Q_1,Q_2,\epsilon}
       (N^{1/2+\epsilon}D^{3+\epsilon}
          +N^{3/4+\epsilon}D^{1+\epsilon}).
\]
With the unchanged Möbius tail, zero blocks, and genus completion,
the full budget is
\[
 O_{Q_1,Q_2,\alpha,\epsilon}\left(
 N^{1+\epsilon}D^{-1}
 +N^{1/2+\epsilon}D^{3+\epsilon}
 +N^{3/4+\epsilon}D^{1+\epsilon}
 +D^2N^\epsilon+N^{1+\epsilon}D^{-\alpha}\right).
\]
Take \(D=N^{1/8}\) and \(\alpha<1\) arbitrarily close to 1.
After allocating arbitrarily small epsilon losses, this proves
\[
 R^*(N)=JN\prod_p\beta_p^*(N)
       +O_{Q_1,Q_2,\epsilon}(N^{7/8+\epsilon}),
\]
equivalently every fixed \(\eta<1/8\). For an explicit
\(\eta=1/9\), take \(\alpha=15/16\) and choose all combined epsilon
losses below \(1/1000\). The largest pre-loss exponent is
\(1-15/128<8/9-1/1000\). The literal endpoint
\(O(N^{7/8})\) is not asserted.
Equations (3) and (2) prove eventual existence.

The same proof extends to any two fixed positive-definite integral
binary forms with fixed coefficient contents \(c_1,c_2\).
After dividing by \(g^2\), the additional content is
\(h=(c_1u^2,c_2v^2)\), which divides
\(\operatorname{lcm}(c_1,c_2)\). Divide the form and target by \(h\)
when possible. These finitely many extra content factors change only
fixed-form constants in the determinant, level and local-invariant
comparisons. The local derivative and Möbius-tail arguments use the
original fixed Gram matrices unchanged. In particular fixed positive
coefficients on either binary form are permitted by this extension.

## Stronger lower bounds and distinct value pairs

Let \(b_p\) be the good-prime factor in (5) when \(p\nmid N\).
Factoring the two lines of (5) gives the exact target ratio
\[
 T_p=\frac{\text{target-prime factor in (5)}}{b_p}
 =\begin{cases}
 (p+3)/(p+2),&\varepsilon_1=\varepsilon_2=1,\\
 (p-1)/(p-2),&\varepsilon_1=\varepsilon_2=-1,\\
 1-1/p,&\varepsilon_1=-\varepsilon_2.
 \end{cases}
 \tag{20b}
\]
Thus, writing \(S^*(N)=\prod_p\beta_p^*(N)\),
\[
 S^*(N)=\left(\prod_{p\notin S}b_p\right)
        \left(\prod_{p\in S}\beta_p^*(N)\right)
        \prod_{\substack{p\mid N\\p\notin S}}T_p .
\]
The first product converges to a positive constant. If
\(\det A_1\det A_2\) is a square, all the \(T_p\) exceed one, and
\[
 c_{\rm loc}:=\left(\prod_{p\notin S}b_p\right)
                 \prod_{p\in S}p^{-5t_p-3}>0,\qquad
 R^*(N)\ge \tfrac12Jc_{\rm loc}N
 \tag{20c}
\]
for every sufficiently large admissible \(N\). In particular
\(R_Q^*(N)\gg_Q N\) for self-sums. More generally this applies to
two forms from the same imaginary quadratic field, even with different
ideal classes or conductors. The displayed lower coefficient is
explicit; a numerical onset still requires an explicit analytic error
constant.

If the determinant product is nonsquare, let \(\psi\) be its associated
nonprincipal quadratic character. The exact ratios imply
\[
 S^*(N)\gg_{Q_1,Q_2}
 \prod_{\substack{p\mid N\\p\notin S,\ \psi(p)=-1}}(1-1/p)
 \gg_{Q_1,Q_2}(\log\log(3N))^{-1/2},
 \qquad
 R^*(N)\gg_{Q_1,Q_2}\frac{N}{\sqrt{\log\log(3N)}}.
 \tag{20d}
\]
For the second inequality use the fixed-modulus Mertens theorem,
[Williams, Theorem 1, printed p.356](https://people.math.carleton.ca/~williams/papers/pdf/057.pdf),
on the \(\varphi(q)/2\) reduced classes where \(\psi=-1\). Its full
product through \(x\) is comparable to \((\log x)^{-1/2}\).
Split the target primes at \(x=\log N\); above \(x\),
\(\sum_{p\mid N,p>x}1/p\le\log N/(x\log x)\), so their product is
bounded below. All constants depend on the fixed forms, not on target
valuations. A linear bound is in general false in this nonsquare case:
multiply the product of the good negative-character primes through
\(x\) by a bounded factor chosen to put it in any fixed admissible
residue class modulo (4). The singular series tends to zero along
these admissible targets. The error in (2) is negligible.

There is also a logarithmic improvement on the distinct-value lower
bound obtained from the pointwise estimate (11). Let \(f_K(n)\) count
integral ideals of norm \(n\) in an imaginary quadratic field \(K\).
The norm embedding used for (11) gives
\(r_Q(n)\le w_Kf_K(4an)\) for \(Q=ax^2+bxy+cy^2\), where \(w_K\)
is the number of units. At split primes
\((j+1)^3\le\binom{j+7}{7}\); at inert and ramified primes the
corresponding coefficient domination is immediate. Hence
\(f_K(n)^3\le f_K^{*4}(n)\), where the right side is the fourfold
Dirichlet convolution. The elementary ideal-count estimate
\(\sum_{n\le X}f_K(n)\ll_K X\) follows from \(f_K=1*\chi_K\)
by Dirichlet's hyperbola method. Partial summation and convolution give
\[
 \sum_{n\le X}f_K^{*4}(n)\ll_K X(1+\log X)^3,\qquad
 \sum_{n\le X}r_Q(n)^3\ll_Q X(1+\log X)^3.
 \tag{20e}
\]
The second inequality bounds the subsequence \(4an\) by the full
sum through \(4aX\).

Let \(r_Q^*(n)\) count primitive binary coordinate representations of
\(n\), and let \(V(N)\) count distinct ordered pairs \((a,N-a)\)
with both values primitively represented by their respective forms.
Put \(w(a)=r_{Q_1}^*(a)r_{Q_2}^*(N-a)\). Hölder and Cauchy give
\[
 R^*(N)\le V(N)^{1/3}\left(\sum_a w(a)^{3/2}\right)^{2/3},
 \qquad \sum_a w(a)^{3/2}\ll_{Q_1,Q_2}N(1+\log N)^3.
\]
Consequently, on sufficiently large admissible targets,
\[
 V(N)\gg_{Q_1,Q_2}
 \frac{N S^*(N)^3}{(1+\log N)^6}.
 \tag{20f}
\]
This yields \(V(N)\gg N/(\log N)^6\) for same-field pairs, and
\(V(N)\gg N/((\log N)^6(\log\log(3N))^{3/2})\) for distinct-field
pairs. The same bounds, with interval-dependent constants, hold in
each fixed positive-length balance interval by (24).
They are lower bounds, not unweighted value-pair asymptotics.

### Fractional moments and the sharp upper-sieve scale

Classical discriminant-uniform multiplicative-function bounds improve
(20f). For every fixed \(\epsilon>0\), all sufficiently large
admissible targets satisfy
\[
 V(N)\gg_{Q_1,Q_2,\epsilon}
       \frac{N}{(\log N)^{\log 4+\epsilon}}.
 \tag{20g}
\]
Here and below \(\log\) is natural. For all sufficiently large targets,
\[
 V(N)\ll_{Q_1,Q_2}\frac{N}{\log N}\,
 \prod_{\substack{p\mid N,\ p\notin S\\
                  \varepsilon_1(p)=\varepsilon_2(p)=-1}}
                        \frac{p-1}{p-2}.
 \tag{20h}
\]
No matching logarithmic-scale lower bound for arbitrary form classes
is claimed.

The exact input is Henriot's **New Theorem 5**, printed p.377 of the
[2014 erratum](https://doi.org/10.1017/S0305004114000280), together with
the growth-class definitions in
[arXiv:1102.1643v1](https://arxiv.org/pdf/1102.1643v1).
The discriminant-uniform constant is essential for the factors
\(X,N-X\). Only this upper-bound theorem is used; its lower-bound
counterpart requires positive weights and does not apply to functions
vanishing at inert primes.

Here are the majorant and local calculations. For
\(Q=ax^2+bxy+cy^2\), let \(K=\mathbb Q(\sqrt{b^2-4ac})\).
At a fixed finite set \(B_Q\) containing the primes of \(4aD_K\), and
at split primes, set \(h_Q(p^j)=j+1\); at all remaining inert primes
set \(h_Q(p^j)=\mathbf1_{2\mid j}\).
The norm embedding used in (11) proves
\[
 r_Q^*(n)\le r_Q(n)\le w_K\tau(4a)h_Q(n).
\]
This accommodates every ideal class and conductor. It avoids the
false inequality \(f_K(4an)\le f_K(4a)f_K(n)\), which can fail at
inert primes.

Fix \(t>1\), put \(c=2^t\), and set
\(c_i(p)=c\,\mathbf1_{\varepsilon_i(p)=1}\) outside a fixed finite
set \(B\supset S\cup B_{Q_1}\cup B_{Q_2}\). Enlarge \(B\) so that the factors below are
positive. The corrected Henriot theorem gives
\[
 \sum_{1\le a<N}
   \bigl(r_{Q_1}^*(a)r_{Q_2}^*(N-a)\bigr)^t
 \ll_{Q_1,Q_2,t}
 N(\log N)^{2^t-2} W_t(N),
\]
\[
 W_t(N)=\prod_{\substack{p\mid N\\p\notin B}}
       \left(1+\frac{(c_1(p)-1)(c_2(p)-1)}p\right).
 \tag{20i}
\]
To verify the uniformity hypotheses, remove the endpoints
\(\min(a,N-a)\le N^{1/2}\), whose divisor bound costs
\(O_\epsilon(N^{1/2+\epsilon})\). Cover the remainder using \(a\)
on the left and \(N-a\) on the right,
by dyadic intervals with \(x\gg N^{1/2}\), \(y=x\), and \(2x<N\),
of total length \(O(N)\).
Use Henriot's \(\alpha=\delta=1/4\), total degree 2, and growth
parameter below \(1/2400\). Since
\(\|X(N-X)\|\asymp N\), all height and interval conditions hold.
The functions \(h_1^t h_2^t\) belong to the prescribed growth class,
bounded by fixed divisor powers.

In the corrected Euler sum, a nontarget prime has first coefficient
\(c_1(p)+c_2(p)\), with root factor \(1-2/p\).
A target prime has first coefficient \(c_1(p)c_2(p)\), with root
factor \(1-1/p\); the corrected cross-coprimality condition excludes
a positive valuation in only one factor. Higher valuation terms are
\(O_t(p^{-2})\), uniformly in \(v_p(N)\), since their tails are
bounded by
\(\sum_{r\ge2}(2r+1)(r+1)^{2t}p^{-r}\).
The target ratio has first coefficient
\((c_1(p)-1)(c_2(p)-1)\), proving (20i).
Fixed-character Mertens products give the baseline exponent \(2^t-2\).
The product \(W_t(N)\) and its reciprocal are bounded by fixed powers
of \(\log\log(3N)\).

Hölder and (20i) give
\[
 V(N)\gg_{Q_1,Q_2,t}
 \frac{N}{(\log N)^{(2^t-2)/(t-1)}}
       \left(\frac{S^*(N)^t}{W_t(N)}\right)^{1/(t-1)}.
\]
The parenthesized expression is bounded below by an inverse fixed
power of \(\log\log(3N)\). Let fixed \(t>1\) approach 1:
\((2^t-2)/(t-1)\to\log4\). Absorbing the remaining log-log factor
proves (20g). The constants may depend strongly on \(\epsilon\).

For (20h), majorize primitive value membership by the indicator
excluding all good inert primes. Apply the same upper theorem,
now with first coefficients \(1\) at split primes and \(0\) at
inert primes. The baseline exponent is \(-1\). The target correction
has a first-order contribution only at primes inert in both fields;
there its exact factor is \((p-1)/(p-2)\). All other ratios are
\(1+O(p^{-2})\). This proves (20h), uniformly over target valuations.
Both bounds also apply on each fixed positive-length balance interval,
using (24) for the lower bound.

The [novelty comparison](PRIMITIVE_NORM_NOVELTY.md) gives the exact
classical attribution and the hypotheses preventing stronger
unweighted correlation theorems from being specialized to one fixed
target.

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

The unrestricted \(Q\oplus Q\) shell formula is classical:
Liouville, *Sur la forme \(x^2+y^2+2(z^2+t^2)\)* (1860),
[pp.269--272](https://www.numdam.org/item/JMPA_1860_2_5__269_0.pdf),
and Deutsch, [arXiv:math/0406434v1](https://arxiv.org/pdf/math/0406434v1),
give direct precedents, including parity and whole-vector primitive
counts. Those counts do not impose the two designated binary gcds.
The calculation below records the additional dyadic filter in the
normalization used by the sieve.

The eventual existence statement in this specialization also has a
classical prime-variable proof. Plaksin's theorem, used with
\(x^2+2y^2\), gives the odd strong core after regrouping and deleting
prime-coordinate gcd failures. Its Gaussian specialization, followed by
\(u=(x+y)/2,\ v=(x-y)/2\) for even targets, gives
\(N=Q(p,u)+Q(q,v)\) throughout \(N\equiv2,4,6\pmod8\).
The deletion costs \(O_\epsilon(N^{1/2+\epsilon})\) against a uniformly
positive main term of order \(N/\log^2N\).
Doubling covers \(v_2(N)=3\). Thus all eventual admissible -8 existence
is a classical deduction; the explicit spectral bound and complete
finite bridge supply additional quantitative and exact information.
See the [novelty comparison](PRIMITIVE_NORM_NOVELTY.md) for the
source hypotheses, local factors, parity bijection, and prime-coordinate
count. This alternative does not supply a numerical onset here.

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
They are holomorphic modular forms at the stated levels. To check
every cusp explicitly, put
\[
 D_d(\tau)=dE_2(d\tau)-E_2(\tau)
   =\frac1{2\pi i}\frac{d}{d\tau}
                   \log\frac{\Delta(d\tau)}{\Delta(\tau)},
 \qquad \Delta(\tau)=\eta(\tau)^{24}.
\]
The quotient is a modular function on \(\Gamma_0(d)\), with no
zero or pole in the upper half-plane and meromorphic at each cusp.
At a cusp of width \(h\), write its expansion as
\(q_c^\nu U(q_c)\), where \(q_c=e^{2\pi i\tau/h}\) and
\(U(0)\ne0\). The weight-two transform of its logarithmic derivative
is
\[
                 \frac1h\left(\nu+q_c\frac{U'(q_c)}{U(q_c)}\right),
\]
which has no pole at that cusp. Thus \(D_d\) is holomorphic of
weight two, with trivial character. Explicitly,
\[
 E_r=\frac{D_2}{12}-\frac{D_4}{12}+\frac{D_8}{6},
 \qquad
 E_s=\frac{D_2}{48}-\frac{D_{16}}{48}+\frac{D_{32}}{24}.
\]
These identities prove membership at levels 8 and 32, including all
cusps. Both constant terms are 1, and their positive coefficients
are integers by the displayed divisor-sum formulas.

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
and an intermediate finite reduction. The
[exact completion](PRIMITIVE_NORM_MINUS8_EXACT.md) proves that its remaining
finite exception set is empty using sharper estimates and complete certificates.

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
Here is the transfer argument directly for this order. Let \(V_\ell\)
be the image of \(O_\ell^\times\) in \(PB^\times(\mathbb Q_\ell)\),
and replace \(V_p\) by its reduction kernel to obtain \(U(p)\).
Class number one identifies
\(PB^\times(\mathbb Q)\backslash PB^\times(\mathbb A_f)/U(p)\)
with \(Y_p\). Its stabilizer is exactly \(\Gamma_p\): local scalar
corrections have even norm valuations and combine into one rational
scalar, leaving an element of \(O^\times\). The transitive right
\(G_p\)-action makes the adelic measure proportional to counting
measure. Principalizing odd-norm ideals acts by \(\alpha^{-1}\)
at \(p\); quaternion conjugation changes this to the shell sum in
(32), including its factor \(1/4\).

To apply Jacquet--Langlands one must exclude constituents with a
finite-dimensional local factor at a split place. Such a smooth
irreducible local factor is a determinant character: normality of its
open kernel and conjugation by diagonal matrices kill both unipotent
subgroups. Strong approximation for \(B^1\), away from infinity and
that split place, then makes the automorphic constituent invariant
under \(B^1(\mathbb A)\), so it is a global norm character.
The hypothesis holds because \(B^1\) is noncompact at the chosen
split place; see Voight, *Quaternion Algebras*, Main Theorem 28.5.3.
On the transitive quotient \(Y_p\), norm characters restrict to
characters of \(G_p\). These are exactly \(1\) and the quadratic
determinant character, since the elementary unipotents generate
\(\mathrm{SL}_2(\mathbb F_p)\). Merely counting finite-group characters
without this strong-approximation step would not justify transfer.

Thus the complement of those two characters satisfies the
infinite-dimensional-at-split-places hypothesis of
[Jacquet--Langlands, Theorem 14.4](https://publications.ias.edu/sites/default/files/automorphic-forms-on-gl2_rpl_9.pdf).
It transfers to cuspidal \(\mathrm{GL}_2\) representations, and the
trivial definite infinite type corresponds to weight two (Section 5,
the \(n=0\) case on printed pp.89--90 and Proposition 5.20, pp.105--106). Every \(\ell\mid m\) is split and has maximal local
level, so the transferred factor is unramified there. The weight-two
Ramanujan theorem gives two Hecke roots of magnitude \(\sqrt\ell\).
The all-integral-ideal recurrence
\(T_{\ell^{a+2}}=T_\ell T_{\ell^{a+1}}-\ell T_{\ell^a}\)
therefore bounds the \(\ell^a\) eigenvalue by
\((a+1)\ell^{a/2}\). Multiplicativity proves (33) with constant 1,
uniformly in \(p,m\). The level at 2 affects the conductor of the
transfer; no unramified estimate at 2 is used. The exact primary source
editions and theorem locators are pinned in the source ledger.
The [earlier Eisenstein lemma](../../unpublished/eisenstein-restricted-prime-sums/sections/prime-bounds.tex)
is a comparison, rather than an additional unexpanded spectral input.

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
identity is equivalent to coverage of a **finite** residual family.
One dyadic stratum can also be removed by leastness. Write
\(C=\{n\text{ odd}:q\mid n\Rightarrow q\equiv1,3\bmod8\}\), so
\(P_{-8}=C\cup2C\). If a least counterexample had \(N=4m\), \(m\)
odd, its smaller admissible target \(2m\) would be represented.
Both summands of \(2m\equiv2\bmod4\) must be odd, since every even
primitive norm is \(2\bmod4\). Doubling those two odd values preserves
primitive norm membership and represents \(N\), a contradiction.
Thus the smaller residual family is
\[
 N=2^e m\ge2,\quad e\in\{0,1,3\},\quad m<2^{60},\qquad
       q\mid m\Longrightarrow q\equiv5,7\bmod8 .
 \tag{40}
\]
The argument through (40) is a finite reduction. For the stratum
\(e=3\), both summands must be even: elements of \(C\) are
\(1,3\bmod8\), so two cannot sum to zero modulo eight. Its existence
problem is exactly \(4m\in C+C\). The separate
[exact -8 proof](PRIMITIVE_NORM_MINUS8_EXACT.md) now covers this family:
107 rational interval certificates, complete strong-core coverage
through \(2^{27}\), and split-prime descent prove
\(P_{-8}+P_{-8}=\{N\ge2:16\nmid N\}\). This does not rely on
exhausting all inert cores below \(2^{60}\) directly.
A least counterexample, if one exists, is therefore smaller than \(2^{63}\) and has no odd
split prime factor. Equation (40) is a least-counterexample reduction;
the reduced core cutoff \(2^{60}\) alone is **not** a global onset
for arbitrary targets with split prime factors. The unconditional
all-target onset proved above is \(2^{78}\).

The [discriminant -11 specialization](PRIMITIVE_NORM_MINUS11_EFFECTIVE.md)
now supplies the numerical global onset `2^71` for every even target
not divisible by121, and at least `N^(2/3)/100000` distinct ordered value
pairs. Its maximal order requires retaining all ideal classes: the
principal mass is3/5, and the ramified reduction has ten equally weighted
points. Sharper lattice/divisor estimates, 35 certified intervals and a
complete finite core computation now prove the exact -11 identity:
`P_-11+P_-11={N>=2:N even,121 not dividingN}`. These are additional
explicit deductions, separate from the general nonnumerical cutoff above.

## Consequences and limits

* The local-global existence conclusion applies to mixed quadratic
  fields, nonprincipal classes and nonmaximal quadratic orders.
  Ramified and conductor primes enter the finite local test (4).
  They do not need to be fitted into the coefficient hypotheses of
  the splitting-semigroup theorem.
* The theorem permits every target valuation compatible with the
  primitive local conditions. It is stronger than applying a theorem
  about one primitive four-vector: both blocks are sieved separately.
* Equations (20c)--(20i) give linear coordinate counts for same-field
  pairs and \(N/\sqrt{\log\log(3N)}\) for distinct-field pairs.
  Fractional moments improve the general distinct-value lower bound
  to \(N/(\log N)^{\log4+\epsilon}\); the upper sieve gives
  \(N J_{12}(N)/\log N\). No value-pair leading asymptotic is claimed.
* The exact \(-8\) identity is reduced to the explicit finite family
  (40), with an unconditional global onset \(2^{78}\). The residual
  inert-core family is now covered by the separate exact -8 proof.
  Both -8 and -11 binary classifications are proved, including ramification.
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
Certificate verification requires active Python assertions. The checker rejects optimized execution (-O, -OO, or optimization enabled through PYTHONOPTIMIZE) before generating a receipt.


The [proof audit](PRIMITIVE_NORM_ANALYTIC_AUDIT.md) records the reviewed
revision, the exact analytic hypotheses, source normalization issues,
the independent spectral specialization, the exponent improvement and
the exact-classification completion and remaining priority obligations. It is an
internal mathematical audit, not external human review.

## Explicit first-power adaptation of Baker

This gives a second proof of the full coordinate asymptotic with a
small positive power saving. It does not recover the stronger
\(7/8+\epsilon\) error above. The input is the analytic estimates
and displayed proofs in R.C.Baker, The values of a quadratic form
at square-free points, Acta Arith.124(2006),101--137:
[retained primary edition](sources/baker-squarefree-points-2006.pdf).
The relevant locations are p.102 (weight scope), (3.1)--(3.4),
Lemmas 2--10, (4.7)--(4.12), and Section 6.
His Theorem 1 imposes individual coordinate squarefreeness.
The first-power lattice changes and local completion below are
an adaptation, not an invocation of that theorem.

### Fixed shell weight and scaled lattice

Put \(F=Q_1\oplus Q_2\). To avoid a half-integral matrix take
\(f=2F=\mathbf x^{\mathsf T}M\mathbf x\), where
\(M=\operatorname{diag}(A_1,A_2)\) is fixed integral nonsingular.
Let \(m=2N>0\) and \(P=\sqrt m\).
Choose a fixed smooth \(\phi\), supported in \((1/2,3/2)\)
and equal to one near 1, and put \(w(\mathbf z)=\phi(f(\mathbf z))\).
Its support is compact and avoids zero; on \(f(\mathbf x)=m\)
we have \(w(\mathbf x/P)=1\) exactly.

This weight includes coordinate axes. Baker's p.102 first states
the underlying Heath-Brown analytic weight class, smooth compact
support avoiding zero, and then adds coordinate-plane vanishing
for his scalar-squarefree inversion. Only that inversion and its
tail use the extra restriction. The analytic integral and
exponential-sum estimates used below remain in the original
weight class, and we replace the inversion.

For positive integers \(\mathbf s=(s_1,\ldots,s_4)\), let
\(D_{\mathbf s}=\operatorname{diag}(s_i)\), \(B=\prod_i s_i\), and
\[
 f_{\mathbf s}(\mathbf x)=f(D_{\mathbf s}\mathbf x),\qquad
 w_{\mathbf s}(\mathbf z)=w(D_{\mathbf s}\mathbf z).
\]
Define the complete sums, with \(e_q(t)=\exp(2\pi it/q)\), by
\[
 S_q(\mathbf s,\mathbf c)=
 \sum_{\substack{a\bmod q\\(a,q)=1}}\ \sum_{\mathbf b\bmod q}
 e_q\!\left(a(f_{\mathbf s}(\mathbf b)-m)+\mathbf c\cdot\mathbf b\right).
\]
The delta-symbol and Poisson identity (Baker (3.4)) is
\[
 r_{f_{\mathbf s}}(m)
 =c_P P^{-2}\sum_{\mathbf c\in\mathbb Z^4}\sum_{q\ge1}
 q^{-4}S_q(\mathbf s,\mathbf c)I_q(\mathbf s,\mathbf c),
 \qquad c_P=1+O_L(P^{-L}).
 \tag{B.1}
\]
The support of the kernel in (3.1) restricts \(q\) to \(q\le P\)
for this radial weight. All counts on the left are raw coordinate
counts, because \(w_{\mathbf s}(\mathbf x/P)=1\) on the shell.

### First-power gradient and complete-sum estimates

The source uses \(s_i=d_i^2\). For first powers, write
\(g_i=(s_i,q)\). The gradient congruence is
\[
                      2D_{\mathbf s}MD_{\mathbf s}\mathbf u=0\bmod q.
\]
For \(\mathbf v=MD_{\mathbf s}\mathbf u\), each condition
\(2s_iv_i=0\bmod q\) permits at most \(2g_i\) values.
The fixed matrix \(M\) has at most \(|\det M|\) preimages
of any vector modulo \(q\), by its Smith form.
Each possible \(D_{\mathbf s}\mathbf u\) then has at most
\(\prod_i g_i\) preimages \(\mathbf u\).
Thus the gradient kernel contains at most
\(C_M\prod_i g_i^2\) vectors. Applying Baker's Lemma 7 proves
the first-power version of Lemma 9:
\[
                  |S_q(\mathbf s,\mathbf c)|
                       \ll_M q^3\prod_i(s_i,q).
 \tag{B.2}
\]
This includes common primes among the \(s_i\); neither their
squarefreeness nor pairwise coprimality was used.

Here is the corresponding target-uniform version of Lemma 10:
\[
 \sum_{q\le X}|S_q(\mathbf s,\mathbf c)|
 \ll_{f,\epsilon}
 B^{2+\epsilon}X^{7/2+\epsilon}(m+|\mathbf c|+1)^\epsilon.
 \tag{B.3}
\]
To check it, split \(q=uv\), where \(u\) is the product of the
primes occurring exactly once in \(q\) and outside
\(2\det(M)B\); \(v\) is squarefull away from that prime set.
Chinese remaindering is legitimate even when the entries of
\(\mathbf s\) overlap, since \((u,v)=1\).
At primes dividing \(u\), the entries \(s_i\) are units.
Changing variables there and applying Baker's Lemma 8 bounds
the good part by
\[
                  C^{\omega(u)}u^{5/2}(u,m)^{1/2}.
\]
The transformed dual value need not equal the unscaled dual
value; bounding the exceptional gcd by \((u,m)\) removes that
issue. Since \(m>0\),
\[
 \sum_{u\le Y}(u,m)^{1/2}
 \le\sum_{u\le Y}(u,m)\le Y\tau(m).
\]
For the bad part, (B.2) gives at most \(Bv^3\).
The sum of \(v^{-1/2}\) is
\(O_{f,\epsilon}(B^\epsilon X^\epsilon)\):
primes dividing \(2\det(M)B\) contribute their elementary
Euler factors, and a squarefull number is \(r^2t^3\) with
\(t\) squarefree, giving a logarithmic bound for the remainder.
These estimates prove (B.3), even with a smaller power of \(B\).
No condition on \((m,B)\) or on any \(v_p(m)\) was imposed.
The zero-target/zero-dual exception in Baker's Lemma 10 cannot
occur when \(m>0\).

### First-power oscillatory integrals

The exact change of variables behind (4.7) now gives
\[
 I_q(\mathbf s,\mathbf c)
       =B^{-1}I_q(c_1/s_1,\ldots,c_4/s_4).
 \tag{B.4}
\]
In particular the Jacobian is \(B^{-1}\). In Baker\'s
parameterization \(s_i=d_i^2\), this same Jacobian is
\((\prod_i d_i)^{-2}\).
Let \(J_f=\pi^2/\sqrt{\det M}\).
Baker's Lemmas 2--5 imply the following deliberately loose bounds:
\[
 \begin{aligned}
 I_q(\mathbf s,0)
   &=\frac{P^4}{B}\left(J_f+O_L((q/P)^L)\right),\\
 |I_q(\mathbf s,\mathbf c)|
   &\ll_L B^L P^5q^{-1}|\mathbf c|^{-L}
                  &&(\mathbf c\ne0),\\
 |I_q(\mathbf s,\mathbf c)|
   &\ll_\epsilon B^2P^{3+\epsilon}q
                         (1+|\mathbf c|)^\epsilon
                  &&(\mathbf c\ne0,\ q\le P).
 \end{aligned}                                                   \tag{B.5}
\]
The fixed-form integral estimates apply to real frequencies;
integrality is not needed after (B.4), just as the source itself
uses the rational frequencies \(c_i/d_i^2\).
Since \(|(c_i/s_i)_i|\ge|\mathbf c|/\max s_i\) and
\(\max s_i\le B\), enlarging the powers of \(B\) gives (B.5).
All implied constants depend only on the fixed forms, weight,
chosen decay order and epsilon.

### A complete growing-divisor error budget

Take
\[
 \mathbf s=(d,d,e,e),\quad d,e\le H,\quad
 H=P^{1/1024},\quad C=P^{1/64},\quad Q_0=P^{3/4}.
\]
Use weights \(\mu(d)\mu(e)\); they are bounded by one in all
absolute analytic estimates. Now \(B=d^2e^2\), and
\[
 B\le H^4,\qquad
 \sum_{d,e\le H}B^2\ll H^{10},\qquad
 \sum_{d,e\le H}B^4\ll H^{18}.
\]
From (B.3), a dyadic interval \(R<q\le2R\) has
\(q^{-4}|S_q|\) sum
\(\ll B^{2+\epsilon}R^{-1/2+\epsilon}
(m+|\mathbf c|+1)^\epsilon\).
Combining this with (B.5) gives the following bounds,
all including the \(P^{-2}\) in (B.1):

| Terms | Bound before arbitrarily small epsilon losses | Exponent of \(P\) |
|---|---|---|
| \(|\mathbf c|>C\), decay order \(L=512\) | \(P^3H^{4L+10}C^{4-L}\) | \(-2.927734375\) |
| \(0<|\mathbf c|\le C\) | \(P^{3/2}H^{18}C^4\) | \(1.580078125\) |
| \(\mathbf c=0,\ q>Q_0\), including series extension | \(P^2H^{10}Q_0^{-1/2}\) | \(1.634765625\) |
| Integral-expansion error at \(\mathbf c=0,\ q\le Q_0\), order \(16\) | \(P^2H^{10}(Q_0/P)^{16}\) | less than \(-1.99\) |

For the first row, sum
\(\sum_{|\mathbf c|>C}|\mathbf c|^{-L+\epsilon}
\ll C^{4-L+\epsilon}\); the additional \(q^{-1}\) makes
the dyadic modulus sum convergent. For the second row,
the factor \(q\) in (B.5) changes the dyadic sum to
\(O(P^{1/2+\epsilon})\), and there are \(O(C^4)\) frequencies.
For the third row, sum the dyadic \(R^{-1/2+\epsilon}\) tail.
For the last row, bounding \(q/P\) by \(Q_0/P\) already suffices.
Thus all analytic errors are \(O(P^{1.7})\) after choosing
epsilon losses sufficiently small. The factor \(c_P-1\)
is negligible by taking its fixed decay order sufficiently large.

The retained main term is
\[
 J_fm\sum_{d,e\le H}
  \frac{\mu(d)\mu(e)}{d^2e^2}
          \sum_{q\ge1}q^{-4}S_q((d,d,e,e),0).
 \tag{B.6}
\]
This is an absolute remainder budget for the new lattice.
It does not use Baker's scalar-squarefree inclusion-exclusion
or his completed squarefree-point singular series.

### Local normalization and the new Euler factors

Let \(\beta_p(f_{\mathbf s};m)\) be the usual unrestricted
density of the scaled equation. By (B.3), its additive-character
series is absolutely convergent and
\[
              \sum_{q\ge1}q^{-4}S_q(\mathbf s,0)
                    =\prod_p\beta_p(f_{\mathbf s};m).
\]
This identity is the unrestricted local-count identity in
Baker's Section 5; the scalar-squarefree Euler product in
that section is not used.

The lattice Jacobian must be tracked locally as well.
For \(h>\max_i v_p(s_i)\), the map
\(\mathbf b\mapsto D_{\mathbf s}\mathbf b\bmod p^h\)
has kernel \(p^{\sum_i v_p(s_i)}\). Hence the density in
the original shell with the prescribed coordinate divisibilities is
\[
        \gamma_{p,\mathbf s}
          =p^{-\sum_i v_p(s_i)}\beta_p(f_{\mathbf s};m).
 \tag{B.7}
\]
Multiplying these identities over the primes recovers exactly
the factor \(B^{-1}\) in (B.6).

For (B.8)--(B.10) use the original \(F,N\) convention for both
\(\gamma\) and \(\beta\); the conversion from \(f=2F,m=2N\)
is accounted for explicitly below. For squarefree \(d,e\), a prime
occurs with exponents \(a,b\in\{0,1\}\) in the two blocks.
Write \(\gamma_p^{ab}\) for the ambient surface mass with the first block divisible
by \(p^a\) and the second by \(p^b\). Signed local inversion gives
\[
             \beta_p^*(N)=
          \gamma_p^{00}-\gamma_p^{10}
                         -\gamma_p^{01}+\gamma_p^{11}.
 \tag{B.8}
\]
This includes the overlap \(p\mid(d,e)\).
In the original \(F,N\) convention that overlap is
\[
 \gamma_p^{11}=
 \begin{cases}
 p^{-2}\beta_p(F;N/p^2),&p^2\mid N,\\
 0,&p^2\nmid N.
 \end{cases}                                                     \tag{B.9}
\]
It cannot be discarded by assuming \((d,e)=1\).
Equations (B.7)--(B.9) hold at every prime, including two,
ramified primes and conductor primes, with unrestricted target
valuations.

For completeness, the good-prime factors are explicitly
recomputed. Put \(\varepsilon_i=(-\det A_i/p)\) for
\(p\nmid2\det A_1\det A_2\). A nondegenerate binary form has
\(p-\varepsilon_i\) vectors of norm \(a\ne0\), and
\(p+(p-1)\varepsilon_i-1\) nonzero vectors of norm zero.
Convolving these two functions over \(\mathbb F_p\), and
dividing by \(p^3\), gives
\[
 \beta_p^*(N)=
 \begin{cases}
 1-\dfrac{2+\varepsilon_1\varepsilon_2}{p^2}
        +\dfrac{\varepsilon_1+\varepsilon_2}{p^3},&p\nmid N,\\[5pt]
 1+\dfrac{\varepsilon_1\varepsilon_2}{p}
  -\dfrac{2+\varepsilon_1+\varepsilon_2+\varepsilon_1\varepsilon_2}{p^2}
  +\dfrac{1+\varepsilon_1+\varepsilon_2}{p^3},&p\mid N.
 \end{cases}                                                     \tag{B.10}
\]
Every retained root is nonsingular, so each has \(p^3\) lifts
per further power. These factors therefore depend on divisibility
by \(p\), not the size of \(v_p(N)\). At bad primes the exact
factor is the primitive-block limit (1), equivalently (B.8).
The finite lifting criterion and uniform density bounds (6)--(7)
apply without alteration.

### Möbius tails, zero blocks and absolute completion

Perform exact double Möbius inversion only on nonzero binary
blocks, as in (10). The fixed-binary bound (11) gives
\[
                       O_\epsilon(N^{1+\epsilon}/H)
 \tag{B.11}
\]
for the large-block-gcd tails, by summing
\(O(N/g^2)\) vectors in a block with gcd multiple of \(g>H\)
and applying (11) to the other positive norm.

The unrestricted circle counts include a zero block.
For first block zero the \(d\)-sum costs at most \(H\), and
\(\sum_{e^2\mid N}r_{Q_2}(N/e^2)\ll_\epsilon N^\epsilon\).
Consequently removing both possible zero blocks costs
\(O_\epsilon(HN^\epsilon)\); both cannot vanish when \(N>0\).
Individual-coordinate axes remain included.

The absolute genus/local moment argument in (15)--(19), with
moment \(1/2\), depends on local lattice masses, not on the
chosen theta estimate. It therefore applies to (B.6):
\[
 \sum_{d,e\ge1}\mu(d)^2\mu(e)^2(de)^{1/2}E_{d,e}(N)
                    \ll_\epsilon N^{1+\epsilon}.
\]
Here \(E_{d,e}(N)\) is the nonnegative real-times-local main
coefficient for the scaled lattice in the original convention.
Thus completing either divisor sum costs
\(O_\epsilon(N^{1+\epsilon}H^{-1/2})\).
This absolute moment includes the overlap (B.9);
no division by a possibly zero local density is used.
Signed completion gives (B.8), rather than Baker's
coordinate-squarefree factors.

Finally clear the initial factor two. At odd primes the
equations \(2F=2N\) and \(F=N\) have the same density.
At two the former has twice the density: reduction from
modulo \(2^h\) to modulo \(2^{h-1}\) gives 16 lifts,
against the normalization factor 8.
This applies with every divisibility/primitive-block condition.
Together with
\(J_f=\pi^2/\sqrt{\det A_1\det A_2}\) and \(m=2N\),
the total main term becomes
\[
 \frac{4\pi^2N}{\sqrt{\det A_1\det A_2}}\prod_p\beta_p^*(N).
\]
The largest pre-epsilon error is the completion
\(N^{1-1/4096+\epsilon}\), since \(H\asymp N^{1/2048}\).
Choosing all epsilon losses sufficiently small proves the
conservative full-shell bound
\[
 R^*(N)=JN\prod_p\beta_p^*(N)+O_{Q_1,Q_2}(N^{1-1/16384}).
 \tag{B.12}
\]
All target valuations and all overlaps between the two divisors
are included.

### What this adaptation establishes historically

The new lattice is not a literal instance of Baker's squared
coordinate lattice. Equations (B.2)--(B.5) give the needed
first-power reconstruction, and (B.7)--(B.10) replace his
squarefree-point Euler factors. The analytic input and the
proof strategy are classical; the modification is short and
uses no new coefficient estimate. It is reasonable to describe
this as a short adaptation of a published convergent-sieve proof,
rather than an immediate application of his stated Theorem 1.
Whether someone previously printed this exact adaptation remains
a separate bibliographic question.

In particular, an unspecified-power-saving general two-block
theorem, qualitative class-obstruction disappearance and genus
invariance alone do not support a strong claim of a new analytic
breakthrough. The sharper \(7/8+\epsilon\) error, explicit local
classification and exact exception theorems must be evaluated
on their own statements and arithmetic content.
