# A pointwise asymmetric Gaussian–Eisenstein convolution

One summand can be restricted to primes congruent to 1 modulo 12 while
all prime factors of the other are congruent to 1 modulo 3. A positive
weighted convolution has an asymptotic for every sufficiently large
congruence-compatible target. The proof combines a four-character
expansion on the first side with a two-character divisor expansion on
the second. It does not prove the corresponding two-sided modulus-12
statement.

This is a research proof with internal independent arithmetic and analytic
reviews. External review and a priority comparison have not been completed.
One-sided existence may already follow from Linnik-type prime plus
primitive binary quadratic norm results; that specialization has not
been audited here. The concrete result proved below is the weighted
pointwise asymptotic.

## Theorem AsymPointwise

Include 1 in the multiplicative sets
\[
 C=\{n\ge1:p\mid n\Longrightarrow p\equiv1\pmod{12}\},\qquad
 A=\{n\ge1:p\mid n\Longrightarrow p\equiv1\pmod3\}.
\]
Let \(\tau\) be the divisor function, and put
\[
 U(n)=\tau(n)^2\mathbf1_C(n),\qquad
 V(n)=\tau(n)\mathbf1_A(n),\qquad
 \mathcal C(N)=\sum_{\substack{a,b\ge1\\a+b=N}}U(a)V(b).
 \tag{1}
\]
There is an absolute constant \(\eta>0\) such that, for every sufficiently
large integer \(N\equiv2\pmod6\),
\[
 \mathcal C(N)=N\mathfrak S(N)+O(N^{1-\eta}).
 \tag{2}
\]
There is no exceptional set of targets. All valuations of every prime
dividing \(N\) are allowed.

Here is an explicit definition of the main term. Let
\[
 \rho_U=\operatorname*{Res}_{s=1}\sum_{n\ge1}\frac{U(n)}{n^s}>0,
 \qquad
 \rho_V=\operatorname*{Res}_{s=1}\sum_{n\ge1}\frac{V(n)}{n^s}>0.
 \tag{3}
\]
For each prime \(p\ge5\), set \(x=p^{-1}\), \(q_p=1-x\), and
\[
 A_e=U(p^e),\quad B_e=V(p^e),\quad
 M_{A,p}=q_p\sum_{e\ge0}A_ex^e,\quad
 M_{B,p}=q_p\sum_{e\ge0}B_ex^e.
 \tag{4}
\]
Thus \(A_0=B_0=1\). For \(e\ge1\),
\[
 A_e=\begin{cases}(e+1)^2,&p\equiv1\ (12),\\0,&\text{otherwise},\end{cases}
 \qquad
 B_e=\begin{cases}e+1,&p\equiv1\ (3),\\0,&\text{otherwise}.\end{cases}
\]
For \(r\ge0\), define
\[
 \begin{split}
 I_{p,r}={}&q_p\sum_{i<r}A_iB_ix^i+x^r(1-2x)A_rB_r\\
 &+q_p\left(B_r\sum_{i>r}A_ix^i+
                   A_r\sum_{j>r}B_jx^j\right).
 \end{split}
 \tag{5}
\]
An empty sum is zero. Then
\[
 \boxed{\displaystyle
 \mathfrak S(N)=6\rho_U\rho_V
    \prod_{p\ge5}\frac{I_{p,v_p(N)}}{M_{A,p}M_{B,p}}.}
 \tag{6}
\]
The factors at primes not dividing \(N\) are \(1+O(p^{-2})\), so (6)
is an absolutely convergent product with finitely many changed factors.
All factors are positive. With absolute positive comparison constants,
\[
 \mathfrak S(N)\asymp T_{\rm asym}(N),
 \tag{7}
\]
where
\[
 T_{\rm asym}(N)=
 \prod_{\substack{p\mid N\\p\equiv1\ (12)}}(1+3/p)
 \prod_{\substack{p\mid N\\p\equiv5,11\ (12)}}(1+1/p)
 \prod_{\substack{p\mid N\\p\equiv7\ (12)}}(1-1/p).
 \tag{8}
\]
In particular, for every \(\varepsilon>0\),
\[
 N^{-\varepsilon}\ll_\varepsilon\mathfrak S(N)
                 \ll_\varepsilon N^\varepsilon.
 \tag{9}
\]
The decreasing factor at primes congruent to 7 modulo 12 belongs to the
natural target normalization. A fixed positive lower bound is not asserted.

## Corollary AsymEveryTarget

Every sufficiently large \(N\equiv2\pmod6\) is a sum
\[
 N=a+b,\qquad a\in C,\quad b\in A.
 \tag{10}
\]
For every fixed \(\varepsilon>0\), the number of these ordered value pairs
is \(\gg_\varepsilon N^{1-\varepsilon}\).

Indeed, choose the exponent in the lower bound (9) smaller than both
\(\eta/2\) and \(\varepsilon/4\). The main term in (2) dominates
its error. The elementary bound
\(U(a)V(b)\ll_\varepsilon N^{\varepsilon/2}\), after reducing the
divisor-bound exponent if necessary, gives the asserted number of pairs.
The order is fixed: the first summand satisfies the modulus-12 restriction.
No balance, explicit cutoff, or unweighted \(N/(\log N)^{5/4}\) estimate
is proved here. The second summand is not required to belong to \(C\).

## 1. Primary analytic inputs and dependencies

The imported dispersion estimate is Proposition 5.1 of
V. Blomer, L. Grimmelt, J. Li and S. Rydin Myerson,
*Additive problems with almost prime squares*, Geometric and Functional
Analysis **33** (2023), 1173–1242, printed pp.1192–1193:
[DOI](https://doi.org/10.1007/s00039-023-00635-w),
[primary published PDF](https://pure.mpg.de/pubman/item/item_3521894_4/component/file_3548627/Blomer-Grimmelt-Li-Myerson_Additive%2Bproblems%2Bwith%2Balmost%2Bprime%2Bsquares_2023.pdf),
[cached arXiv v5](sources/blomer-grimmelt-li-rydin-myerson-2111.01601v5.pdf).
Proposition 5.1 is PDF p.15 in that version. Its printed conditions,
including the location of its absolute value and the coprimality
conditions, are used here. We abbreviate this input as BGLM.

The following specialization records the aspects needed below. In a
bilinear box of product scale \(X\), there is a sufficiently small
\(\eta_0>0\), and an associated \(\delta_0>0\), such that the short
factor may have length between \(X^{\eta_0}\) and
\(X^{1/4+\eta_0}\). Both coefficients may be divisor-bounded.
The moduli go up to \(X^{1/2+\delta_0}\); auxiliary moduli and the
truncated-character conductor \(R\) may be at most \(X^{\delta_0}\).
The main residue may have size at most \(X^{1+\delta_0}\), and the
inverse multiplier size at most \(X^{\delta_0}\). The signed averaged
error from the model with primitive conductors at most \(R\) is
\(\ll X(\log X)^{O(1)}R^{-1/3}\). The source permits its stated small
auxiliary-modulus average. A fixed auxiliary modulus may be selected by
an indicator weight, retaining the same full error bound.

The source requires that the short variable be coprime to the inverse
multiplier and that the varying modulus be coprime to the residue and
specified auxiliary factors. These conditions are checked below.
One may select a reduced class of the varying modulus modulo 12.
Crucially, the absolute value is **outside its modulus sum**. We use
neither an \(L^1\) sum of individual progression errors nor arbitrary
modulus-dependent phases.

The second input is
[Theorem FixedPeriodicBoxes](FIXED_PERIODIC_DIVISOR_BOXES.md). For one,
two or three interval factors multiplied by arbitrary functions of
period 12 and modulus at most 1, their Dirichlet convolution at product
scale \(Y\) has reduced-progression error
\[
 \ll Y^{1-\sigma_0}/Q,\qquad
 \sigma_0=1/10000,\quad (Q,12)=(a,Q)=1,
 \quad Q\le Y^{1/2+\sigma_0}.
 \tag{11}
\]
Pad with singleton factors at 1 when there are fewer than three factors;
bounded-variation cutoffs follow by partial summation. Its full proof
transfers the estimates in Sections 2–6 of D. R. Heath-Brown,
*The divisor function \(d_3(n)\) in arithmetic progressions*, Acta
Arithmetica **47** (1986), 29–56:
[DOI](https://doi.org/10.4064/aa-47-1-29-56),
[primary scan](https://matwbn.icm.edu.pl/ksiazki/aa/aa47/aa4713.pdf).
Completion modulo \(12Q\) introduces separate bounded local Gauss
factors, removed before the deep incomplete correlation is estimated.
Thus (11) is the proved periodic extension, not an application of an
indicator-only theorem to character-weighted intervals. No estimate for
moduli with unbounded 2- or 3-parts is imported.

The dependency graph is
\[
 \begin{gathered}
 \text{BGLM Proposition 5.1}+\text{FixedPeriodicBoxes}
       +\text{elementary character means and kernel bounds}\\
 \Downarrow\quad\text{signed progression estimate of Section 4}\\
 \Downarrow\quad\text{gcd strata and two-character hyperbola}\\
 \Downarrow\quad\text{AsymPointwise}\quad\Longrightarrow\quad
 \text{AsymEveryTarget}.
 \end{gathered}
 \tag{12}
\]

## 2. Dirichlet series and literal character means

All four characters \(\psi\bmod12\) have their literal values, including
zero on nonunits. Euler factors give
\[
 D_U(s):=\sum_{n\ge1}\frac{U(n)}{n^s}
       =E_U(s)\prod_{\psi\bmod12}L(s,\psi),
 \tag{13}
\]
where, for \(z=p^{-s}\) and \(p\ge5\),
\[
 E_{U,p}(s)=
 \begin{cases}1-z^2,&p\equiv1\ (12),\\
               (1-z^2)^2,&p\equiv5,7,11\ (12).
 \end{cases}
 \tag{14}
\]
There are no local factors at 2 or 3. For the second weight, use the
principal and nonprincipal characters literally modulo 6:
\[
 D_V(s)=E_V(s)L(s,\psi_0^{(6)})L(s,\chi_3^{(6)}),
 \qquad
 E_V(s)=\prod_{\substack{p\ge5\\p\equiv2\ (3)}}(1-p^{-2s}).
 \tag{15}
\]
Here \(\chi_3\) is the quadratic character of discriminant \(-3\),
with its value at even integers removed when viewed modulo 6.
All correction products converge absolutely for \(\Re s>1/2\).
Their coefficients have a fixed divisor-function bound and finite
absolute norm at \(\Re s=3/5\).

The literal principal character modulo 6 or 12 has residue \(1/3\)
at 1. Each nonprincipal factor is a fixed real Dirichlet \(L\)-function
with positive value at 1, and every correction factor in (14)–(15) is
positive there. This proves the positivity and existence of (3).

For \(B\in C\), normalize the shifted weights by
\[
 F_B(n)=\frac{U(Bn)}{U(B)},\qquad
 G_B(n)=\frac{V(Bn)}{V(B)}.
 \tag{16}
\]
They retain the four- and two-character factorizations. At \(p^v\Vert B\),
the replacement correction for the first series is
\[
 (1-z)^4\sum_{i\ge0}
       \left(\frac{v+i+1}{v+1}\right)^2z^i,
 \tag{17}
\]
a polynomial of degree at most 3 with uniformly bounded coefficients.
For the second it is
\[
 (1-z)^2\sum_{i\ge0}\frac{v+i+1}{v+1}z^i,
 \tag{18}
\]
a polynomial of degree at most 1, again uniformly bounded in \(v\).
If the first sequence is additionally restricted by \((n,T)=1\),
replace its whole local correction at each \(p\mid T\) by the polynomial
cancelling the four corresponding \(L\)-factors. This replacement has
uniform degree and coefficient bounds. It does not require a bound for
the reciprocal of an arbitrary shifted local series.

Consequently every kernel \(e\) used here satisfies
\[
 \sum_{j\ge1}|e(j)|j^{-3/5}
       \ll_\varepsilon (BT)^\varepsilon,
 \qquad |e(j)|\le \tau(j)^{O(1)}.
 \tag{19}
\]
The coefficients' fixed exponent is independent of \(B,T\). The product
of bounded local norms over \(p\mid BT\) is
\(\ll_\varepsilon(BT)^\varepsilon\), by separating finitely many small
primes. All kernel indices in the actual families are 1 modulo 12:
these families are the target-deleted \(F_B\) and the **undeleted**
\(G_B\), with \(B\in C\). Base corrections are supported on squares,
and exceptional nonsquare indices occur only at primes of \(B\), or
at allowed primes deleted from \(F_B\), all congruent to 1 modulo 12.
No such support assertion is made for arbitrary Euler deletions of
\(G_B\).

Here are the elementary character estimates used in the model. If
\(\chi\) has primitive conductor \(d\nmid12\), every one of the four
literal mod-12 characters twisted by \(\chi\) is nonprincipal, with
period dividing \(12d\). In a dyadic four-variable box, sum last over
the largest variable; the product of the other three scales is at most
\(Y^{3/4}\), and a nonprincipal periodic sum is \(O(d)\). Summing the
boxes and convolving the absolutely convergent kernel gives
\[
 \sum_{n\le Y}F_B(n)\mathbf1_{(n,T)=1}\chi(n)
    \ll_\varepsilon d(BTY)^\varepsilon Y^{3/4}.
 \tag{20}
\]
The factor \(j^{-3/4}\) in the kernel convolution is retained, so (19)
suffices. Imprimitive twists add only unit-removal divisor losses.

The untwisted mean has the same \(Y^{3/4+o(1)}\) error. After extracting
the principal factor, the convolution of the three remaining fixed
nonprincipal characters has partial sums
\(O(Y^{2/3}(\log Y)^{O(1)})\), by the largest-variable argument.
Hyperbola splitting its convolution with 1 at \(Y^{3/4}\) gives
\(O(Y^{3/4}(\log Y)^{O(1)})\): the long side uses those partial sums,
and the short side has at most \(Y^{3/4}\) terms. Its kernel preserves
the error. Thus, with
\(H_{B,p}(z)=\sum_{i\ge0}F_B(p^i)z^i\),
\[
 \begin{split}
 \sum_{n\le Y}F_B(n)\mathbf1_{(n,T)=1}
   &=\rho_{B,T}Y+O_\varepsilon((BTY)^\varepsilon Y^{3/4}),\\
 \rho_{B,T}&=\rho_{F_B}\prod_{p\mid T}H_{B,p}(1/p)^{-1}.
 \end{split}
 \tag{21}
\]
For an additional modulus \(L\) prime to \(6T\), let
\[
 E_B(L)=\prod_{p\mid L}H_{B,p}(1/p)^{-1}.
\]
The principal-character density in a reduced class modulo \(L\) is
\[
 \frac{\rho_{B,T}E_B(L)}{\varphi(L)}.
 \tag{22}
\]
This specifies the principal model, not an individual high-modulus
progression estimate. Extra unit removal at \(L\) has the same
polynomial-kernel proof and contributes only a divisor loss.

Orthogonality, (20), and the reciprocal \(\varphi\)-weights show that
the nonprincipal conductors at most \(R\) in the BGLM model contribute
at most
\[
 O_\varepsilon(R^6N^{3/4+\varepsilon})
 \tag{23}
\]
after the relevant modulus sums, with an additional fine-box count if
boxes are used. This generous bound absorbs character counts and small
auxiliary divisor factors. Because the moduli are prime to 6, a
nonprincipal conductor in that model cannot divide 12. Thus (20)
applies to every nonprincipal model character.

## 3. Absolute truncations and a parameter hierarchy

Against a divisor-bounded companion, a kernel index \(j>K\) contributes
in absolute value at most
\[
 N^{1+\varepsilon}\sum_{j>K}\frac{|e(j)|}{j}
       \ll_\varepsilon N^{1+2\varepsilon}K^{-2/5}.
 \tag{24}
\]
Count multiples of \(j\) and use (19); signed kernels cause no difficulty.
The bound also applies to the second kernel. An original summand at most
\(N^{1-\lambda}\) contributes \(O(N^{1-\lambda+\varepsilon})\).

Choose \(\eta_0\) permitted by BGLM so that
\(\eta_0<\sigma_0/1000\), and let \(\delta_0\) be its BGLM output.
Choose fixed positive parameters in the order
\[
 \begin{gathered}
 0<r<\frac1{10000}\min(\eta_0,\delta_0,\sigma_0,1/100),\\
 0<\kappa<r/10000,\qquad 0<\varepsilon<\kappa/10^6,\\
 R=N^r,\quad K=N^\kappa,\quad
 \lambda_0=10\kappa,\quad\zeta=\kappa/100,
 \quad\Delta=N^{-\zeta}.
 \end{gathered}
 \tag{25}
\]
Decrease \(\kappa\), chosen last, to meet the finite list of strict
inequalities below. These choices give no numerical cutoff. A sufficiently
small final saving \(\eta<\kappa/10000\) will suffice.

Discard the original summand \(a\le N^{1-\kappa}\). In the application
\(a=gm\), with \(g\le K\), so \(m\ge N^{1-2\kappa}\). After a
first kernel index \(j\le K\), the four-character product scale \(X\)
is at least \(N^{1-3\kappa}\). After extracting a scalar
\(c_0\le N^{10\kappa}\), the BGLM scale \(x\) is at least
\(N^{1-13\kappa}\).

An explicit sufficient list of source inequalities is
\[
 \begin{gathered}
 11\kappa<(1-13\kappa)\delta_0,\qquad
 1<(1-13\kappa)(1+\delta_0),\\
 \tfrac12<(1-13\kappa)(\tfrac12+\delta_0),\qquad
 r<(1-13\kappa)\delta_0,\\
 \eta_0(1-3\kappa)\ge10\kappa(1-\eta_0),\\
 \tfrac12+\kappa<
 (1-3\kappa)(1-8\eta_0)(\tfrac12+\sigma_0).
 \end{gathered}
 \tag{26}
\]
All hold with sufficiently small \(\kappa\). The last inequality places
the **full** modulus, including the second kernel, in the periodic-box
range.

## 4. The signed progression estimate on the hyperbola

The estimate needed below is restricted to the actual application.
Take \(B\in C\), \(B\le K\), and target residue \(A\ll N\).
The first sequence is \(F_B(n)\mathbf1_{(n,A)=1}\). Thus the deletion
modulus in (21) is **\(T=A\)**. No independent arbitrary \(T\) and
extra varying-modulus filter are introduced. The inverse multipliers
and extra unit congruences are polynomially small, and the varying
moduli are prime to \(6A\). In each reduced class modulo 12, the signed
sum of progression errors over the hyperbola region is
\(O(N^{1-\eta_1})\), for some fixed \(\eta_1>0\), with principal
density (22). The following proves precisely this assertion.

### Medium variables and short-factor coprimality

Truncate the first kernel at \(K\). In a four-character box choose a
variable whose scale belongs to
\[
 [X^{2\eta_0},\ X^{1/4+\eta_0}]
 \tag{27}
\]
if one exists. Combine the other three factors into the arbitrary longer
coefficient, and apply BGLM to this bilinear convolution.

The chosen variable need not initially be coprime to the inverse
multiplier. In this application that multiplier is the first kernel
index \(j\). Write the variable uniquely as \(c_0n'\), where
\(c_0\mid j^\infty\) is its entire part supported on primes of \(j\),
and \((n',j)=1\). The portion \(c_0>N^{10\kappa}\) satisfies
\[
 \sum_{\substack{c_0>N^{10\kappa}\\c_0\mid j^\infty}}\frac1{c_0}
 \le N^{-5\kappa}
       \prod_{p\mid j}(1-p^{-1/2})^{-1}
 \ll_\varepsilon N^{-5\kappa}j^\varepsilon.
 \tag{28}
\]
For any fixed \(t>0\), the analogous product over prime divisors of an
integer is \(\ll_{t,\varepsilon}\) that integer to \(\varepsilon\).
Counting multiples with \(\lfloor Y/c_0\rfloor\le Y/c_0\), and
reordering the actual divisor choices, turns (28) into an actual-count
error \(O(N^{1-5\kappa+\varepsilon})\). There is no number-of-scalars
endpoint loss. The corresponding nonprincipal model uses (20) and
the smooth Euler product at \(3/4\), not an \(R^2\) multiple of this
tail.

For the remaining \(c_0\), use inverse multiplier \(jc_0\). If the
original short scale is \(S\), its lower bound
\(S/c_0\ge(X/c_0)^{\eta_0}\) follows from (26). Division by \(c_0\)
only decreases the upper relative exponent, preserving the upper
BGLM short bound. The new short variable is coprime to the full inverse
multiplier, and its coefficient remains divisor-bounded.

### Boxes with no medium variable

If (27) contains no variable scale, every scale is either below
\(X^{2\eta_0}\) or above \(X^{1/4+\eta_0}\). There are at most three
long variables. Fix the short variables, whose total scale is at most
\(X^{8\eta_0}\). The product \(Y\) of the long scales is at least
\(X^{1-8\eta_0}\). Apply (11) to the one, two or three remaining
character factors, padding with singletons if necessary.

The full modulus is at most \(kq\le N^{1/2+\kappa}\); (26) puts it
below \(Y^{1/2+\sigma_0}\). It is prime to 12. A fixed scalar sharing
this modulus gives no actual solution for a unit residue, and gives
zero in the principal unit model. Otherwise invert that scalar in the
residue; (11) is uniform there. Summing the short values gives, crudely,
\[
 O(N^{1-\sigma_0+8\eta_0+O(\kappa)+\varepsilon}),
 \tag{29}
\]
which saves a power. This three-factor step is not supplied by the
bilinear dispersion proposition alone.

### Shared auxiliary factors and exact BGLM parameters

Write \(q=q_kq'\), where \(q_k\mid k^\infty\) is the entire part
supported on the second kernel index \(k\), and \((q',k)=1\).
Discard \(q_k>N^{10\kappa}\) by the smooth Rankin estimate. In an
actual progression, sum first
\(N^\varepsilon(N/(kq)+1)\): the reciprocal part has the Rankin
saving, and the endpoint is \(O(N^{1/2+\varepsilon})\) times a smooth
summable factor. The second-kernel square-root endpoint costs at most
\(K^{1/10}\), by (19). Model tails use (20) and smooth summation at
\(3/4\).

In the application below, the BGLM parameters are
\[
 a_1=M,\quad a_2=jc_0,\quad a_3=1,\quad
 d=kq_k,\quad c=12.
 \tag{30}
\]
Terms with \((j,kq)>1\) have no actual solution because
\((M,kq)=1\); their principal unit model is also zero. In every
remaining term, \(q'\) is coprime to \(a_1a_2d\), and the auxiliary
residue modulo \(d\) is a unit. Also \((n',a_2)=1\), as required.
The sizes are
\[
 |a_2|,d\le N^{11\kappa},\quad |a_1|\le N,
 \quad Q'\le N^{1/2},\quad x\ge N^{1-13\kappa}.
\]
Thus (26) verifies the exact source conditions. A fixed \(d\) is
selected by an indicator in the auxiliary average. All moduli and
auxiliary factors are prime to 6. No unbounded 2- or 3-parts are
introduced.

The two kernels have total absolute mass at most
\(N^{6\kappa/5+\varepsilon}\), by (19). Counting the \(c_0,q_k\)
choices wastefully by \(N^{20\kappa}\), and using fine boxes in five
coordinates, is safely covered by the aggregate BGLM error
\[
 O(N^{1-r/3+100\kappa+\varepsilon}).
 \tag{31}
\]
The hierarchy (25) leaves a strict saving. Within each reduced class
of \(q'\bmod12\), \(\chi_3(q')\) is constant. The other character
factors are fixed in that application. Hence the signed source estimate
is exactly the kind required: no arbitrary \(q'\)-dependent coefficient
is inserted. A fixed subinterval of a dyadic modulus range is obtained
by differences of cumulative sums and dyadic decomposition, costing
only logarithms.

### Moving boundaries and the truncated-character model

Use relative-width \(\Delta=N^{-\zeta}\) boxes in the four character
variables and in \(q\), at most \(N^{5\zeta+\varepsilon}\) boxes.
For boxes crossing \(m+kq^2=M\), restore the actual divisibility
conditions and write \(M-m=kqt\). A crossing box forces
\[
 |kq(t-q)|\ll\Delta N.
\]
The number of such \((q,t)\) is
\[
 O\left(\frac{\Delta N\log N}{k}+\sqrt{\frac Nk}\right).
 \tag{32}
\]
At any original \(m\), recombine all first-kernel indices:
\[
 \sum_{j\mid m}|e(j)|\tau_4(m/j)\ll_\varepsilon N^\varepsilon.
\]
There is no \(K^{3/5}\) loss from the first kernel. The second kernel
uses \(\sum|e(k)|/k<\infty\), and its endpoint uses
\(\sum_{k\le K}|e(k)|k^{-1/2}\ll N^\varepsilon K^{1/10}\).
The actual boundary error is
\[
 O_\varepsilon\left(N^\varepsilon
        (\Delta N\log N+\sqrt N K^{1/10})\right).
 \tag{33}
\]
The same reasoning handles cumulative cuts in the original summand.

The principal model's boundary volume at kernel index \(j\) retains
its factor \(1/j\), giving an absolutely convergent sum. Bound every
nonprincipal model box with (20), retaining \(j^{-3/4}\). Formula
(23), the fine-box count, and even an additional \(N^{100\kappa}\)
loss save a power under (25). Unit removal at the full modulus has
only a divisor loss. One must not instead bound a model boundary or
scalar tail by \(R^2\) times its actual size: the smooth scalar sum at
\(3/4\) is finite and is used separately.

Only the principal model remains, with exactly (22). Equations
(24), (28), (29), (31), and (33), together with (23), prove the signed
progression assertion. Boundaries save \(\zeta\), kernels save
\(2\kappa/5\), and the other terms save larger fixed powers under
the hierarchy. A final saving below \(\kappa/10000\) is permissible.

## 5. Gcd strata and the strict two-character hyperbola

In (1), let \(g=(a,N)\). Since \(U(a)>0\) implies \(a\in C\),
every nonzero stratum has \(g\in C\) and \(g\mid N\). The strata
\(g>K\) contribute \(O(N^{1-\kappa+\varepsilon})\), by divisor bounds
and \(\sum_{g\mid N,g>K}N/g\). For \(g\le K\), write
\[
 M=N/g,\qquad a=gm,\qquad b=g(M-m),\qquad(m,M)=1.
\]
The weight is exactly \(U(g)V(g)F_g(m)G_g(M-m)\). Since
\(g\equiv1\pmod6\), we still have \(M\equiv2\pmod6\).

Put
\[
 g_{E,6}(t)=(1*\chi_{-3})(t)\mathbf1_{(t,6)=1}.
\]
The second kernel gives
\[
 G_g(t)=\sum_{k\mid t}e_g(k)g_{E,6}(t/k).
 \tag{34}
\]
Truncate \(k\le K\) by (24). If \((k,M)>1\), no term survives,
since \(m=M-kt\) would share that prime with \(M\). Otherwise every
nonzero residual \(t\) is 1 modulo 6, because \(M\equiv2\pmod6\),
\(m\equiv1\pmod{12}\), and the actual kernel support has
\(k\equiv1\pmod{12}\).

For such \(t\), divisor pairing is exact:
\[
 g_{E,6}(t)=
 2\sum_{\substack{q\mid t\\q<\sqrt t}}\chi_3(q)
 +\mathbf1_{t=\square}\chi_3(\sqrt t).
 \tag{35}
\]
Indeed \(\chi_3(t)=1\), so paired divisors have the same character.
The strict inequality matters. Replacing it by a weak inequality in
the subsequent cumulative expression changes only a square term.
That square contribution is
\[
 \ll_\varepsilon N^\varepsilon\sqrt M
        \sum_{k\le K}|e_g(k)|k^{-1/2}
 \ll_\varepsilon N^{1/2+\kappa/10+\varepsilon}
 \ll N^{3/5+\varepsilon}.
 \tag{36}
\]
The outer \(g\)-strata add only a divisor loss \(N^\varepsilon\).

Consequently the normalized fixed-\(g\) count, apart from the reserved
power errors, is
\[
 2\sum_{\substack{k\le K\\(k,M)=1}}e_g(k)
 \sum_{\substack{q\le\sqrt{M/k}\\(q,6M)=1}}\chi_3(q)
 \sum_{\substack{m\le M-kq^2\\m\equiv M\ (kq)}}
       F_g(m)\mathbf1_{(m,M)=1}.
 \tag{37}
\]
The restriction \((q,M)=1\) follows from \((m,M)=(k,M)=1\).
The first kernel in Section 4 is that of
\(F_g\mathbf1_{(\cdot,M)=1}\), so its target residue and deletion
modulus are both \(M\), exactly the case \(T=A\) proved there.

Apply Section 4, separating the reduced classes of \(q\bmod12\).
The main term is
\[
 2\rho_{g,M}
 \sum_{\substack{k\le K\\(k,M)=1}}e_g(k)
 \sum_{\substack{q\le\sqrt{M/k}\\(q,6M)=1}}
  \chi_3(q)\frac{E_g(kq)(M-kq^2)}{\varphi(kq)}.
 \tag{38}
\]
We now evaluate this expression rather than invoke an unspecified
singular series.

## 6. Completion of the main term and the factor 6

For fixed \(g\), write \(H_A(z)=\sum_{i\ge0}F_g(p^i)z^i\), and set
\[
 W_p=[(1-1/p)H_A(1/p)]^{-1}=1+O(1/p),
\]
uniformly in \(g\). For fixed \(k\),
\[
 \frac{E_g(kq)}{\varphi(kq)}=\frac{w_k(q)}{kq},\qquad
 w_k(q)=\prod_{p\mid k}W_p
          \prod_{\substack{p\mid q\\p\nmid k}}W_p.
 \tag{39}
\]
The series of \(\chi_3(q)w_k(q)\), with \((q,6M)=1\), is
\(L(s,\chi_3^{(6)})\) times a correction absolutely convergent for
\(\Re s>0\), with finite Euler modifications at \(kM\). The product
of their norms is \(\ll_\varepsilon(gMk)^\varepsilon\). Convolution
with the bounded partial sums of that fixed nonprincipal character gives
\[
 \sum_{q\le Q,(q,6M)=1}\chi_3(q)w_k(q)
       \ll_\varepsilon(gMkQ)^\varepsilon.
 \tag{40}
\]
Partial summation therefore completes the reciprocal \(q\)-sum in
(38) to infinity and removes its \(q^2\) term with error
\(O_\varepsilon(N^\varepsilon\sqrt{M/k})\). After summing \(k\),
(19) bounds this by \(N^{3/5+\varepsilon}\). Complete the \(k\)-sum
with a power error by (24). The outer divisor strata retain these
power savings.

For \(p\nmid6M\), put \(x=1/p\),
\(H_B(z)=\sum_{i\ge0}G_g(p^i)z^i\), and let \(e_i\) be its
local kernel coefficients relative to \(g_{E,6}\). Then
\[
 \sum_{i\ge0}e_ix^i=(1-x)(1-\chi_3(p)x)H_B(x).
\]
The completed local double sum from (38) is exactly
\[
 \begin{split}
 1+\frac1{(1-x)H_A(x)}
   \left(\frac{\sum_{i\ge0}e_ix^i}{1-\chi_3(p)x}-1\right)
 &=\frac{M_A+M_B-1}{M_A},\\
 M_A=(1-x)H_A(x),\qquad M_B&=(1-x)H_B(x).
 \end{split}
 \tag{41}
\]
These \(M_A,M_B\) refer to normalized shifted weights in this fixed
stratum. At \(p\mid M\), both \(k\) and \(q\) are excluded, so the
double sum has local factor 1; the first mean contributes
\(H_A(x)^{-1}\).

The literal principal character residue gives
\[
 \rho_{F_g}=\frac13\prod_{p\ge5}M_A.
\]
Here and in the next display, products not explicitly absolutely
convergent are in increasing prime order. The fixed nonprincipal
\(L\)-values and absolute corrections in (13)–(18) justify this
interpretation. The fixed-\(g\) main coefficient, including the factor
2 in (38), is therefore
\[
 \frac23\prod_{\substack{p\mid M\\p\ge5}}(1-1/p)
          \prod_{p\nmid6M}(M_A+M_B-1),
 \tag{42}
\]
multiplied by \((N/g)U(g)V(g)\). In completing the \(q\)-sum the only
*new* conditional factor is \(L(1,\chi_3^{(6)})\); the first mean
\(\rho_{F_g}\) already contains the other fixed \(L\)-values.

Every local factor in (42) is positive. Complete the positive divisor
sum over \(g\in C\), \(g\mid N\); its omitted main tail is
\(O(N^{1-\kappa+\varepsilon})\). All shifted local means and their
finitely modified products have at most \(N^\varepsilon\) growth,
uniformly in \(g\), and
\(\sum_{g\mid N,g>K}1/g\ll N^{-\kappa+\varepsilon}\).

Recover (5) without a restriction on \(r=v_p(N)\). Choices
\(v_p(g)=i<r\) contribute
\[
 q_p\sum_{i<r}A_iB_ix^i.
\]
The choice \(i=r\) contributes
\[
 x^r\left((1-2x)A_rB_r+
 q_p\left[B_r\sum_{j>r}A_jx^{j-r}
              +A_r\sum_{j>r}B_jx^{j-r}\right]\right).
\]
This is exactly the rest of (5). At a prime not congruent to 1 modulo
12, only \(i=0\) is possible, since \(A_i=0\) for \(i\ge1\);
the formula still applies. Hence (42) and the divisor sum give
\((2/3)\prod_{p\ge5}I_{p,v_p(N)}\) in the stated prime order.
Finally,
\[
 \rho_U=\frac13\prod_{p\ge5}M_{A,p},\qquad
 \rho_V=\frac13\prod_{p\ge5}M_{B,p}.
\]
The coefficient is exactly (6). The factor 6 records the restrictions
at 2 and 3 with literal-character normalization. Combining the power
errors proves (2).

## 7. Local positivity and uniform target comparison

One may also derive (5) by integrating the two local valuation weights
on \(a+b=N\). Below \(r=v_p(N)\), both valuations equal \(i\),
of measure \(q_px^i\). When both equal \(r\), the measure is
\((1-2x)x^r\). If one exceeds \(r\), the other equals \(r\), giving
the two remaining tails. This accounts for every singular target class.
At \(r=0\),
\[
 I_{p,0}=M_{A,p}+M_{B,p}-1.
 \tag{43}
\]
Every local factor is positive: the unit–unit stratum alone has positive
measure for \(p\ge5\).

For \(p\equiv5,7,11\pmod{12}\), \(A_e=0\) for \(e\ge1\),
so \(I_{p,r}=q_p\) for all \(r\ge1\). In classes 5 and 11,
\[
 M_{A,p}=M_{B,p}=q_p,\qquad I_{p,0}=1-2x,
\]
giving \(I_{p,r}/I_{p,0}=1+x+O(x^2)\). In class 7,
\[
 M_{A,p}=q_p,\quad M_{B,p}=q_p^{-1},\quad
 I_{p,0}=q_p+q_p^{-1}-1=\frac{1-x+x^2}{1-x},
\]
so \(I_{p,r}/I_{p,0}=1-x+O(x^2)\).

For \(p\equiv1\pmod{12}\),
\[
 M_{A,p}=\frac{1+x}{(1-x)^2},\quad
 M_{B,p}=\frac1{1-x},\quad I_{p,0}=1+4x+O(x^2).
\]
Substituting \(r=1\) in (5) gives \(I_{p,1}=1+7x+O(x^2)\).
For \(r\ge2\), the first two terms of
\(q_p\sum_{i<r}(i+1)^3x^i\) give \(1+7x\), and all remaining
terms of (5) are \(O(x^2)\), uniformly in \(r\). Sum the
polynomially weighted geometric tails; here \(p\ge13\).
Consequently, uniformly for every \(r\ge1\),
\[
 \frac{I_{p,r}}{I_{p,0}}=1+\frac{c_p}p+O(p^{-2}),\qquad
 c_p=\begin{cases}3,&p\equiv1\ (12),\\
                  1,&p\equiv5,11\ (12),\\
                 -1,&p\equiv7\ (12).
      \end{cases}
 \tag{44}
\]
The unit factor in (6) is \(1+O(p^{-2})\), since
\(M_{A,p}=1+(A_1-1)x+O(x^2)\) and likewise for \(M_{B,p}\).
Every implied constant here is absolute and uniform in the target
valuation. Thus the ratio of (6) to (8) is bounded above and below by
positive convergent products. Finitely many small primes also have
uniform positive bounds, from (5) and its unit–unit stratum. This proves
(7).

For (9), absorb finitely many small primes into the constant. At every
other prime \(1-1/p\ge p^{-\varepsilon}\), while all positive target
factors are at least 1. The negative product in (8) is therefore at
least a constant times \(N^{-\varepsilon}\). The upper bound follows
from \(1+3/p\le p^\varepsilon\) for sufficiently large \(p\).
This completes the proof.

## 8. Scope and review record

The argument works because the second weight has a two-character
expansion. Its divisor pairing produces one family of square-root
moduli to which BGLM and FixedPeriodicBoxes apply. Replacing \(V\)
by a weight supported on \(C\) introduces four character factors on
both sides. No such pointwise estimate is proved here. The every-target
\(C+C\) theorem, its raw simultaneous-norm lower bound, and its
pointwise bad-prime marginal estimates remain open.

There is no balance restriction in (1). A uniform shifted-coefficient
interface for a later roughness/cofactor argument is not asserted.
Removing the divisor weights at the expected \((\log N)^{5/4}\) scale
would need additional work. The unweighted conclusion proved here is
the elementary \(N^{1-\varepsilon}\) corollary. No novelty, effective
cutoff, or finite computational completion is claimed.

Independent internal arithmetic and analytic campaigns reviewed this
proof on 18 November 2026. The arithmetic review reconstructed the
all-valuation local formula, gcd/kernel/hyperbola identities, and factor
6. Small exact rational diagnostics supported those checks but did not
replace the proof. The analytic review inspected the printed BGLM
Proposition 5.1 and reconstructed (25)–(26), its signed-modulus
restriction, scalar tails, and boundary estimates. A separate independent
review checked the primary-source transfer proving FixedPeriodicBoxes.
These are scoped internal mathematical reviews, not formal verification
or external peer review. Publication priority, including comparison with
prime plus primitive quadratic norm methods for one-sided existence,
remains unresolved.
