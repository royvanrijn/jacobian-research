# Pair correlations for quadratic character indicators

**Scope.** This note proves a fixed-shift pair asymptotic for a
multiplicative indicator with even valuations at inert primes and
arbitrary fixed valuation masks at finitely many primes. It uses the
same explicit ZF and weighted-dilation inputs as the
[Connors--Keating argument](CONNORS_KEATING_INDICATOR_CORRELATIONS.md),
through the [shared arithmetic-model comparison theorem](SPARSE_ARITHMETIC_MODEL_COMPARISON.md).
The theorem is about the stated Euler indicator. Identifying it with
values of a particular quadratic form is a separate arithmetic step;
no arbitrary principal-form or class-number assertion is implicit.

The sums-of-two-squares indicator direction was explicitly proposed
by Walsh, Section 1.2.4; see the exact
[source comparison](CONNORS_KEATING_INDICATOR_CORRELATIONS.md#1-statement-and-imported-inputs).
This note completes the arithmetic model checks and local constants
for the stated quadratic-indicator family within that framework.
The existing ZF/G inputs remain part of its proof. Section 9 gives
the mixed fixed-shift correlation between the square and hexagonal forms.

## 1. The indicator and its coupled bad-prime factor

Fix a nonprincipal primitive real quadratic Dirichlet character
\(\chi\) of conductor \(q\). Fix a finite set \(B\) of primes
containing all primes dividing \(q\). For every \(p\in B\),
fix a mask \(a_p:\mathbb Z_{\ge0}\to\{0,1\}\), with
\(a_p(0)=1\). Define a multiplicative indicator \(f\) by
\[
 f(p^e)=
 \begin{cases}
 a_p(e),&p\in B,\\
 1,&p\notin B,\ \chi(p)=1,\\
 \mathbf1_{2\mid e},&p\notin B,\ \chi(p)=-1.
 \end{cases}                                             \tag{QP.1}
\]
Set \(f(n)=0\) for nonpositive integers. Write
\[
 F_p(s)=\sum_{e\ge0}f(p^e)p^{-es},\quad
 F(s)=\sum_{n\ge1}f(n)n^{-s},\quad
 t_p=(1-p^{-1})F_p(1),\quad A_B=\prod_{p\in B}t_p.
                                                               \tag{QP.2}
\]
In particular \(A_B>0\). The masks may have infinite exponent
support; bounded coefficients make every \(F_p\) holomorphic on
\(\Re s>0\).

Write \(q=\prod_{p\mid q}p^{k_p}\), and decompose \(\chi\)
by CRT as the product of its local primitive unit characters
\(\chi_p\) modulo \(p^{k_p}\). On
\(Z_B=\prod_{p\in B}\mathbb Z_p\), outside the measure-zero
set with a zero coordinate, put
\[
 e_p=v_p(z_p),\qquad r(\boldsymbol e)=\prod_{p\in B}p^{e_p},
\]
and define the unit character
\[
 \chi_B(\boldsymbol z)=
 \prod_{p\mid q}\chi_p\!\left(
       \frac{z_p}{p^{e_p}}
       \prod_{\substack{\ell\in B\\\ell\ne p}}\ell^{-e_\ell}
                              \right).
                                                               \tag{QP.3}
\]
Every factor is evaluated on a \(p\)-adic unit. Define
\[
 P_B(\boldsymbol z)=
       \left(\prod_{p\in B}a_p(e_p)\right)
                    \frac{1+\chi_B(\boldsymbol z)}2,
 \qquad
 W_B(h)=\int_{Z_B}P_B(\boldsymbol z)P_B(\boldsymbol z+h)\,d\boldsymbol z.
                                                               \tag{QP.4}
\]
Here \(h\) is added in every coordinate and the Haar measure is
normalized. Values on the zero-coordinate set do not affect the
integral. For a positive integer \(n\), with
\(n_B=\prod_{p\in B}p^{v_p(n)}\), the character in (QP.3)
is exactly \(\chi(n/n_B)\).

This bad-prime profile is generally coupled. Dividing out a power
of one bad prime changes the unit at the other bad primes. Replacing
(QP.4) by the product of separate ramified-prime norm densities is
not valid for an arbitrary character with several ramified primes.
Conditional on all valuations, at least one nontrivial local
\(\chi_p\) averages to zero on its unit group. Since the masks
depend only on valuations, the marginal identity is
\[
              \int_{Z_B}P_B\,d\boldsymbol z=A_B/2.        \tag{QP.5}
\]

## 2. Statement and normalization

Let
\[
 K_f=\frac1{\Gamma(1/2)}\lim_{s\downarrow1}(s-1)^{1/2}F(s).
                                                               \tag{QP.6}
\]
Its positive real determination is fixed by \(F(s)>0\) for
real \(s>1\). The exact Euler identity is
\[
 F(s)^2=\zeta(s)L(s,\chi)
 \prod_{\substack{p\notin B\\\chi(p)=-1}}(1-p^{-2s})^{-1}
 \prod_{p\in B}F_p(s)^2(1-p^{-s})(1-\chi(p)p^{-s}).       \tag{QP.7}
\]
In particular
\[
 K_f^2=\frac{L(1,\chi)}\pi
 \prod_{\substack{p\notin B\\\chi(p)=-1}}(1-p^{-2})^{-1}
 \prod_{p\in B}F_p(1)^2(1-p^{-1})(1-\chi(p)/p)>0.
                                                               \tag{QP.8}
\]
For a fixed nonzero integer \(h\), define
\[
 C_f(h)=\frac{4L(1,\chi)}\pi\,W_B(h)
       \prod_{p\in B}\frac{1-\chi(p)/p}{1-p^{-1}}
       \prod_{\substack{p\notin B,\ \chi(p)=-1\\p\mid h}}
             \frac{1-p^{-v_p(h)-1}}{1-p^{-1}}.            \tag{QP.9}
\]
Valuations of a negative shift mean those of its absolute value.

**Theorem QuadraticCharacterPairs.** Assume the scoped ZF and G
inputs stated in [the analytic input record](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md).
For every fixed data set (QP.1), every fixed nonzero \(h\), and
every \(\Phi\in C_c^\infty((0,\infty))\),
\[
 \sum_n\Phi(n/X)f(n)f(n+h)
   =C_f(h)\frac X{\log X}\int_0^\infty\Phi(u)\,du
                       +o_{f,h,\Phi}(X/\log X).         \tag{QP.10}
\]
Consequently
\[
 \sum_{1\le n\le X}f(n)f(n+h)
   =C_f(h)\frac X{\log X}+o_{f,h}(X/\log X).             \tag{QP.11}
\]
If \(W_B(h)>0\), (QP.11) is a relative asymptotic. If
\(W_B(h)=0\), the pair count is identically zero.
There is no uniformity in varying conductors, masks or shifts.

The hypotheses retain one common zero-free half-plane for all
Dirichlet conductors and heights and the specified weighted-dilation
root moment with its original ideal-kernel conditions. The theorem
does not independently revalidate either imported preprint. Its
proof below supplies the indicator-specific conditions and the
separate arithmetic model evaluation.

## 3. Rough density with exactly two principal branches

Put \(L=\log X\), \(y=\exp(L^{2/5})\),
\(S=\{p:p\le y\}\), \(V=\prod_{p\le y}(1-p^{-1})\),
and \(E=\exp(\sqrt L)\). Eventually \(B\subset S\).
For \(n=a v\), \(a=n_S\), write the rough series as
\[
 F_R(s)=F(s)\prod_{p\le y}F_p(s)^{-1}
                =(s-1)^{-1/2}H_R(s).                    \tag{QP.12}
\]
The left side is defined by its rough Euler factors; any zeros
of a removed finite factor are removable in this expression.
Once \(B\subset S\), these factors are independent of the
bad-prime masks: at \(p>y\) they are
\((1-p^{-s})^{-1}\) at split primes and
\((1-p^{-2s})^{-1}\) at inert primes.
One has
\[
 \frac{H_R(1)}{\Gamma(1/2)}
     =K_f\prod_{p\le y}F_p(1)^{-1},\qquad
 \frac{|H_R^{(j)}(1)|}{j!}
                 \ll_j H_R(1)(C\log y)^j.              \tag{QP.13}
\]
The latter bound follows by differentiating the removed Euler
product in a disk of radius \(c/\log y\); the inert-square
correction in (QP.7) is absolutely convergent to the right of
one half.

For a sufficiently large fixed depth \(J\), take
\[
 M_J(t)=\sum_{j=0}^{J}\frac{H_R^{(j)}(1)}{j!}
                    \frac{(\log t)^{-1/2-j}}{\Gamma(1/2-j)},
 \quad
 D(n)=f(a)(1+\chi(v))\frac{M_J(v)}V,\quad \xi=f-D.       \tag{QP.14}
\]
Use the actual Hankel formula throughout
\(X^{1/4}\le t\le X^3\). On this fixed power range,
the leading term of \(M_J\) is positive and its following
terms have relative size \(O_j(L^{-3j/5})\). Extend it smoothly
and nonnegatively outside that bulk so that
\(D\ll_J L^{-3/10}\) on \(1\le n\le X^3\).
This makes \(0\le D\le1\) eventually for each fixed depth.

For any rough integer with \(f(v)=1\), one has \(\chi(v)=1\).
Thus the principal twist and the \(\chi\) twist of the rough
indicator coincide coefficientwise. These are exactly the two
principal branches, and both must be centered.

Here are the uniform one-variable estimates required for that
assertion. For a test of period at most \(L^{B_0}\), refine its
modulus to a multiple of \(q\). Twisting the square of (QP.12)
by a character \(\eta\) gives
\(L(s,\eta)L(s,\eta\chi)\), analytic inert-square
corrections, and explicit removed-prime factors. Only an inflated
principal character or an inflated \(\chi\) gives a pole at one;
all primes in its modulus are at most \(y\).

Use smooth Mellin inversion at height \(\exp(L^{1/10})\).
The classical zero-free region, together with Siegel's theorem
for a possible exceptional real character of polylogarithmic
conductor, permits the line \(\Re s=1-c/\log y\).
Its \(X\)-saving is \(\exp(-cL^{3/5})\), whereas the
Euler products and character sums cost only the allowed fixed
logarithmic factors and smaller exponential factors.
At either principal branch, Hankel expansion to depth \(J\)
has relative remainder \(O_J(L^{-3(J+1)/5})\).
The unweighted rough population has principal density \(V\);
its projection by \(1+\chi(v)\) therefore has the same two
principal means as the actual rough indicator.

It follows, for every requested fixed accuracy \(D_0\) and
period/smoothness bound \(B_0\), that a sufficiently large fixed
depth gives
\[
 \left|\sum_{v\ S\text{-rough}}\xi(rv)a(v)W(rv/Y)\right|
                 \ll_{B_0,D_0,J}(Y/r)L^{-D_0},           \tag{QP.15}
\]
uniformly for \(S\)-smooth \(r\le E\),
\(X^{1/2}\le Y\le X^2\), \(|a|\le1\) periodic of
period at most \(L^{B_0}\), and the smooth weights of
condition (C3) in the shared theorem. The rough variable is
still in a fixed positive power range because \(r\le E=X^{o(1)}\).
The factor \(f(r)\) is outside the rough sum. All structured
bank weights likewise depend on that fixed record, not arbitrarily
on \(v\). Polynomial-logarithmic Mellin phases are included in
the smooth derivative budget before choosing \(J\).
No binary correlation estimate has entered (QP.15).

## 4. Verifying the shared comparison hypotheses

Choose a fixed even comparison modulus \(m\) divisible by
\(q\) and every prime in \(B\). Use its identity prime
class and the scale \(\mathcal S_X=X/L\) on the fixed-shift
line. The required checks are as follows.

**(C1): boundedness.** Section 3 gives a nonnegative bounded
model, eventually at most one for every fixed depth. Its endpoint
moment exponents and its envelope exponents are independent of
the final depth; only constants and thresholds change.

**(C2): common exact invariance.** A prime \(p\le y\) with
\(p\equiv1\pmod m\) is split and outside \(B\). Therefore
\(f(pn)=f(n)\), including repeated factors. It changes only
the smooth record in (QP.14), leaving \(v\) unchanged and
\(f(a)\) unchanged. Hence \(D(pn)=D(n)\) and
\(\xi(pn)=\xi(n)\) exactly. The chosen class has fixed
positive prime density. Deleting primes dividing the fixed
shift preserves the required bank masses and bank counts.

**(C3): every recordwise principal mean.** This is (QP.15).
The fixed comparison modulus can have more character classes
than \(\{1,\chi\}\), because its identity prime banks form
a smaller class than all split primes. Those extra characters
have no additional rough pole; their means are among the
nonprincipal estimates already covered by arbitrary periodic
tests. On gcd profiles use primitive characters whose conductor
divides the actual residual modulus, as in the shared theorem.
No inflated character is evaluated on a nonunit bad-prime record.

**(C4): the positive sparse envelope.** Let \(\mathcal P\)
be a permitted union of the small and big identity-class banks.
Each bank prime avoids \(h\), so its joint valuation law is
\[
 \mu_p(0,0)=1-2/p,\qquad
 \mu_p(e,0)=\mu_p(0,e)=(1-1/p)p^{-e}\quad(e\ge1),        \tag{QP.16}
\]
with no mass on two positive valuations.
The tilted smooth-part moment removes records exceeding \(E\),
with error \(XL^{O(1)}\exp(-cL^{1/10})\). There are at most
\(E^2\) retained pairs \(r,s\). They are coprime, and their
divisibilities specify one progression in the interval
\(I\subset[X,2X]\), with \(|I|/(rs)+O(1)\) terms.
Exactness at the remaining bank primes yields the factors in
(QP.16), including the single extra-root exclusion at a prime
dividing a record.

At an inert prime outside \(B\), a necessary condition at an
actual endpoint is that its valuation is not exactly one.
This excludes \(p-1\) classes modulo \(p^2\). For
\(p\nmid h\), the two forms exclude \(2(p-1)\) classes.
Their slopes stay units there because record primes are split.
Drop the finite primes of \(B\), inert primes dividing \(h\),
and two if necessary, for this upper bound.
The prime-square upper sieve thus has the same dimension as in
the ordinary Gaussian-indicator argument, with only a convergent
\(1+O(p^{-2})\) correction to its first-order factors.

Take a Selberg level \(z=X^{1/100}\), with nonbank cutoff
\(w=z^{1/16}\). Its denominator can be bounded directly for
prime-square atoms: if the forbidden probability is \(\rho_p\),
put \(g(p)=\rho_p/(1-\rho_p)\). Under the Euler-product
subset law, \(\mathbb E\log d=\sum\rho_p\log p\ll\log w\).
Markov therefore leaves a fixed positive fraction of its mass in
\(d\le z\). This gives the usual lower bound for the Selberg
denominator. The same argument includes bank atoms and the
different cutoffs on actual and model endpoints.

For \(ff\) the nonbank denominator has order \(\log w\).
For \(Df\) it has order \((\log y\log w)^{1/2}\), and
for \(DD\) it has order \(\log y\). Multiplying by the
corresponding model bounds \(L^{-3/10}\) and \(L^{-3/5}\)
gives \(L^{-1}\) in every case. The CRT endpoint budget is
\(z^{4+o(1)}\); summing over \(E^2=X^{o(1)}\) records
is power smaller than \(X^{1/2}\). Consequently, for each
nonnegative record test \(F_0\le L^{B_0}\),
\[
 \sum_{n\in I}F_0(n_{\mathcal P},(n+h)_{\mathcal P})
                         w_1(n)w_2(n+h)
 \ll_{f,h,B_0,J}\frac{|I|}L\mathbb E_\mu F_0
       +X^{1/2}L^{C_{B_0}}
       +XL^{C_{B_0}}e^{-cL^{1/10}},                     \tag{QP.17}
\]
for all \(w_i\in\{f,D\}\). Its constants may depend on
the masks but its endpoint and logarithmic exponents do not
increase with the final fixed depth. This is (C4), uniformly
in both interval endpoints, and uses only an upper sieve.

The fixed-shift clause of the shared theorem now gives, for
some fixed \(\delta_m>0\),
\[
 \sum_{n\in I}\{f(n)f(n+h)-D(n)D(n+h)\}
                        \ll_{f,h,J}XL^{-1-\delta_m}.    \tag{QP.18}
\]
It uses the fixed-shift low/high Mellin argument, rather than
applying reflected target-grid energy outside its frequency
range. In that argument the high prime sum expands the fixed
identity class in characters modulo \(m\); its finitely many
principal integrals decay above a polylogarithmic frequency,
and ZF supplies the smooth prime error at polynomial heights.
The two free big products give the required cofactor mean-square
length. Those analytic steps are the shared theorem's proof,
not new assumptions of this application.
Partial summation and a finite dyadic partition insert a fixed
smooth positive-bulk weight into (QP.18).

## 5. The full model product

Let \(b_{I,y}(n)\) impose even valuations at all inert primes
outside \(B\) up to \(y\). Since an allowed good-prime smooth
record has character one, (QP.13)--(QP.14), outside the negligible
large-record set, give
\[
 D(n)=\lambda_y P_B(n)b_{I,y}(n)(1+O_J(L^{-1/2})),\qquad
 \lambda_y=\frac{2K_f}{A_B\sqrt L}
       \prod_{\substack{p\le y,\ p\notin B\\\chi(p)=-1}}(1+1/p).
                                                               \tag{QP.19}
\]
Here \(P_B(n)\) means the diagonal embedding of \(n\) in
\(Z_B\). To verify the normalization, each prime in \(B\)
contributes \(F_p(1)^{-1}/(1-p^{-1})=t_p^{-1}\) to
the ratio in (QP.13). A good split prime contributes one and a
good inert prime contributes \(1+1/p\).
The retained record has \(\log(n/n_S)=\log X+O(\sqrt L)\)
on a fixed positive bulk, producing the displayed relative error.

We evaluate the support of the leading model without assuming
binary cancellation. Truncate every bad-prime valuation so that
\(p^{v_p(n)},p^{v_p(n+h)}\le E\). The discarded count is
\(O_B(X/E+1)\), and its Haar probability is \(O_B(E^{-1})\).
Fixing all these valuations fixes their cross-prime character
phases in (QP.3). Its remaining unit conditions are a finite
union of residue cells modulo the fixed prime powers
\(p^{k_p^*}\), where \(k_p^*=\max\{1,v_p(q)\}\)
for every \(p\in B\). The extra digit at a prime outside
the conductor fixes its valuation exactly. Thus each bad-prime stratum is a fixed finite
union of ordinary congruence classes with modulus
\(E^{O_B(1)}\).

Write the permitted good inert parts of the two endpoints as
\(d^2,e^2\), with \(d,e\) supported on inert primes outside
\(B\) at most \(y\). Discard \(d>E\) or \(e>E\), at
cost \(O_h(X/E+\sqrt X)\). On the retained strata the joint
divisibilities and bad unit conditions give either empty cells
or a fixed finite union of progressions with moduli
\(E^{O_B(1)}=X^{o(1)}\). The residual quotient forms must
avoid each inert prime up to \(y\). They exclude zero, one,
or two residue roots; a cell excluding the full residue field
is empty and is removed first.

There are at most \(O_{B,q}(E^2(\log E)^{2|B|})\)
retained cells, including the fixed number of unit choices.
Their moduli are bounded by \(O_{B,q}(E^{4+2|B|})\).
Use upper and lower fundamental-lemma weights of level
\(D_0=X^{1/100}\).
The sieve dimension is at most one, uniformly in the prefixes,
and its parameter is
\(\log D_0/\log y\asymp L^{3/5}\to\infty\).
Its relative error is superlogarithmic. The CRT remainder summed
over the cells is at most
\(E^{O_{B,q}(1)}D_0^2(\log D_0)^{O(1)}=X^{1/50+o(1)}\),
using the squarefree sieve root bound \(2^{\omega(d)}\).
It is negligible compared with
\(X/L^A\) for each fixed \(A\).

The total main mass of the prefix sum is bounded. Compatible
inert square prefixes satisfy \((d,e)^2\mid h\); their
masses are bounded by constant multiples of
\(1/\operatorname{lcm}(d^2,e^2)\), with a summable total.
At a bad prime the two valuations cannot both exceed
\(v_p(h)\), so its compatible joint masses also have a
summable geometric tail. Completing these disjoint strata
therefore gives
\[
 \begin{split}
 &\sum_n\Phi(n/X)P_B(n)P_B(n+h)b_{I,y}(n)b_{I,y}(n+h)\\
 &\quad=XW_B(h)
    \prod_{\substack{p\le y,\ p\notin B\\\chi(p)=-1}}
                  \delta_p(h)\int\Phi
                              +o_{f,h,\Phi}(X/\log y),
 \end{split}                                             \tag{QP.20}
\]
where the even-valuation local density is
\[
            \delta_p(h)=\frac p{p+1}(1-p^{-v_p(h)-1}).    \tag{QP.21}
\]
This formula holds also at two when it is an inert prime
outside \(B\); empty quotient cells are removed as above.
It follows by partitioning the valuations below, equal to,
and above \(v_p(h)\), and summing the even geometric series.

Because \(\lambda_y^2\asymp_f(\log y)/L\), the absolute
error in (QP.20) becomes \(o(X/L)\) after multiplication.
For \(p\nmid h\), the factor
\((1+1/p)^2\delta_p(h)\) is \(1-p^{-2}\). For
\(p\mid h\), its ratio to that factor is
\((1-p^{-v_p(h)-1})/(1-p^{-1})\). Finally (QP.8) gives
\[
 \frac{K_f^2}{A_B^2}
 \prod_{\substack{p\notin B\\\chi(p)=-1}}(1-p^{-2})
       =\frac{L(1,\chi)}\pi
                 \prod_{p\in B}\frac{1-\chi(p)/p}{1-p^{-1}}.
                                                               \tag{QP.22}
\]
The remaining Euler tail is absolutely convergent. Equations
(QP.18)--(QP.22) prove (QP.10), with precisely the constant
(QP.9). In particular the individual finite-prime mask masses
cancel from the normalization; their pair dependence remains
in the coupled integral \(W_B(h)\).

The coefficient is unchanged if a good prime is moved into
\(B\) with its old valuation rule. At a split prime its new
bad factor and its prefactor are both one. At an inert prime
\(\ell\), the even mask makes its cross-prime phase
\(\chi(\ell)^{v_\ell(n)}\) equal to one. Thus the new
bad density is \(W_B(h)\delta_\ell(h)\), and the new
prefactor gains \((\ell+1)/(\ell-1)\). Their product is
\[
 \frac{\ell+1}{\ell-1}\delta_\ell(h)
      =\frac{1-\ell^{-v_\ell(h)-1}}{1-\ell^{-1}},
\]
the previous good-prime correction, equal to one when
\(\ell\nmid h\). This also checks the normalization when
using a common enlarged bad set for several applications.

## 6. Zero constants, sharp cutoffs and assurance scope

If an integer pair has both indicators equal to one, its good
prime parts have character one. Thus its diagonal embeddings
in \(Z_B\) have both profiles in (QP.4) equal to one.
All valuations are finite, and each unit character condition is
constant on an open finite-level cylinder after those valuations
are fixed. A positive-integer solution would consequently give
\(W_B(h)>0\). This proves the identically-zero assertion when
\(W_B(h)=0\).

For a sharp cutoff, use nonnegative smooth upper and lower
approximations on fixed dyadic intervals. The interval envelope
(QP.17), with the unweighted record test, bounds the part below
\(\epsilon X\) by \(O_{f,h}(\epsilon X/\log X)+o(X/\log X)\).
First let \(X\) grow, then shrink the smooth transition widths
and \(\epsilon\). The finitely many excluded nonpositive
endpoints for negative \(h\) do not affect the limit. This
proves (QP.11).

The analytic argument is an application proof under the explicitly
recorded imports. Its arithmetic steps are the two-branch rough
model, the prime-square upper envelope, and the coupled finite-place
model evaluation. A finite local checker can validate factors and
form classifications; it cannot establish the imported infinite
zero-free or graph estimates. External human review, formal
verification, varying-parameter uniformity and an effective
numerical onset remain separate.

## 7. Smilansky's five ordinary quadratic-form indicators

Let \(b_d(n)\) indicate the positive integer values of
\(x^2+dy^2\), with \(x,y\in\mathbb Z\), and set it to
zero at nonpositive integers. Coordinates may be zero, and no
primitivity condition is imposed. Smilansky's
[Theorem 10 and Section 6.3](https://arxiv.org/html/1208.2419v2)
give the exact value criteria and the conjectured constants.
With primitive fundamental characters, their criteria are
\[
 \begin{array}{c|c|c|c}
 d&\chi&B&\text{valuation masks in }B\\\hline
 1&\chi_{-4}&\{2\}&\text{all exponents}\\
 2&\chi_{-8}&\{2\}&\text{all exponents}\\
 3&\chi_{-3}&\{3\}&\text{all exponents}\\
 4&\chi_{-4}&\{2\}&a_2(e)=\mathbf1_{e\ne1}\\
 7&\chi_{-7}&\{2,7\}&a_2(e)=\mathbf1_{e\ne1},\quad a_7(e)=1.
 \end{array}                                             \tag{QP.23}
\]
Every prime outside \(B\) follows (QP.1). In particular two
is inert for \(\chi_{-3}\), so its valuation must be even
in the \(d=3\) case. The diagonal forms for \(d=3,4,7\)
belong to conductor-two orders, with discriminants \(-12,-16,-28\).
Their ordinary value sets are exactly (QP.23); the conductor
does not introduce a further restriction for \(d=3\).
For example the Eisenstein units move an odd norm into
\(\mathbb Z[\sqrt{-3}]\), and even norms have even
two-valuation. The finite restriction at two is essential for
\(d=4,7\).

Put \(Q_d=\{p:\chi_D(p)=-1\}\). These sets and the
exceptional-place sets \(R_d=B\) are

| \(d\) | \(Q_d\) | \(R_d\) | \(c_d\) | \(C_d(1)\) |
|---|---|---|---|---|
| 1 | \(p\equiv3\pmod4\) | \(\{2\}\) | \(2\) | \(1/2\) |
| 2 | \(p\equiv5,7\pmod8\) | \(\{2\}\) | \(2\sqrt2\) | \(1/\sqrt2\) |
| 3 | \(p\equiv2\pmod3\), including 2 | \(\{3\}\) | \(2/\sqrt3\) | \(1/(3\sqrt3)\) |
| 4 | \(p\equiv3\pmod4\) | \(\{2\}\) | \(2\) | \(1/4\) |
| 7 | \(p\equiv3,5,6\pmod7\) | \(\{2,7\}\) | \(2\sqrt7/3\) | \(1/(2\sqrt7)\) |

The bad profile (QP.4) is a product of the following local
predicates in these five cases. Write \(u_p(n)=n/p^{v_p(n)}\).

| Place | Local predicate |
|---|---|
| \(d=1,p=2\) | \(u_2(n)\equiv1\pmod4\) |
| \(d=2,p=2\) | \(u_2(n)\equiv1,3\pmod8\) |
| \(d=3,p=3\) | \(u_3(n)\equiv1\pmod3\) |
| \(d=4,p=2\) | \(u_2(n)\equiv1\pmod4\) and \(v_2(n)\ne1\) |
| \(d=7,p=2\) | \(v_2(n)\ne1\) |
| \(d=7,p=7\) | \(u_7(n)\) is a nonzero square modulo 7 |

For \(d=7\), the removal of a two-power in (QP.3) does not
couple the two places, because \(\chi_{-7}(2)=1\).
If \(a=v_p(h)\), the local pair probabilities are
\[
 \begin{aligned}
 W_{1,2}(h)&=\begin{cases}1/4,&a=0,\\1/2-3/2^{a+2},&a\ge1,\end{cases}\\
 W_{2,2}(h)&=\begin{cases}1/4,&a=0,1,\\1/2-3/2^{a+1},&a\ge2,\end{cases}\\
 W_{3,3}(h)&=1/2-1/3^{a+1},\\
 W_{4,2}(h)&=\begin{cases}
       1/8,&a=0,\\0,&a=1,\\5/16,&a=2,\\
       3/8-3/2^{a+2},&a\ge3,
       \end{cases}\\
 W_{7,2}(h)&=\begin{cases}1/2,&a=0,1,\\3/4,&a\ge2,\end{cases}\\
 W_{7,7}(h)&=1/2-2/7^{a+1}.
 \end{aligned}                                           \tag{QP.24}
\]

**Corollary SmilanskyFiveForms.** Under the same ZF and G inputs,
for each \(d\in\{1,2,3,4,7\}\) and every fixed nonzero
integer \(h\),
\[
 \sum_{1\le n\le X}b_d(n)b_d(n+h)
       =C_d(h)\frac X{\log X}+o_{d,h}(X/\log X),
 \quad
 C_d(h)=c_d\prod_{p\in R_d}W_{d,p}(h)
       \prod_{\substack{p\in Q_d\\p\mid h}}
              \frac{1-p^{-v_p(h)-1}}{1-p^{-1}}.           \tag{QP.25}
\]
The constant is positive except for \(d=4\) and
\(v_2(h)=1\), where the count is identically zero. In
particular
\[
 \boxed{\displaystyle
 \sum_{n\le X}b_2(n)b_2(n+1)
               \sim\frac1{\sqrt2}\frac X{\log X}.}       \tag{QP.26}
\]

**Proof.** Apply Theorem QuadraticCharacterPairs using (QP.23).
Its bad density is \(\prod_{p\in R_d}W_{d,p}(h)\).
The primitive-character values are
\[
 L(1,\chi_{-4})=\frac\pi4,\quad
 L(1,\chi_{-8})=\frac\pi{2\sqrt2},\quad
 L(1,\chi_{-3})=\frac\pi{3\sqrt3},\quad
 L(1,\chi_{-7})=\frac\pi{\sqrt7}.
\]
Substitution into (QP.9) gives exactly the table's \(c_d\).
The factor at the split bad prime two for \(d=7\) is one.
The masks at two in \(d=4,7\) affect the pair probability,
while their marginal masses have already canceled in (QP.22).
This is Smilansky's constant in equations (6.2)--(6.4).
At shift one the inert divisor product is empty, and (QP.24)
gives the last column, including (QP.26). The zero case also
follows directly from \(x^2+4y^2\equiv0,1\pmod4\). \(\square\)

For completeness, the factors in (QP.24) have short independent
derivations. At a ramified odd prime \(p=3,7\), let the
unit part be a quadratic residue. At \(v_p(h)=0\), the
nonzero residue pairs number \((p-3)/4\), since
\(\chi_p(-1)=-1\) and
\(\sum_u\chi_p(u)\chi_p(u+h)=-1\). Exactly one of the
two zero-residue neighborhoods has the other endpoint a square;
its conditional norm mass is \(1/(2p)\). The initial pair
probability is therefore \((p-1)/(4p)\). For higher
valuations the recurrence is
\[
 W_p(a)=\frac{p-1}{2p}+\frac1pW_p(a-1),
 \qquad W_p(0)=\frac{p-1}{4p},
\]
giving \(1/2-(p+1)/(4p^{a+1})\).

At two the first two unit predicates are invariant under
multiplication by two. Their odd-shift base probability is
\(1/4\). For the \(d=1\) predicate the odd-endpoint
contribution is zero at \(a=1\), and \(1/4\) at
\(a\ge2\); add half the probability at \(a-1\).
For \(d=2\), counting the intersection of \(\{1,3\}\)
with its translate modulo eight gives odd contributions
\(1/8,0,1/4\) at \(a=1,2,\ge3\), respectively;
again add half the preceding probability. This proves the
first two dyadic formulas. For \(d=4\), odd shifts give
\(1/8\), shifts with \(a=1\) are impossible, and for
\(a\ge2\) the odd contribution is \(1/4\) while the
even endpoints are both divisible by four. Thus
\(W_{4,2}(h)=1/4+W_{1,2}(h/4)/4\). Finally the
\(d=7\) two-adic predicate excludes just one residue,
\(2\pmod4\); its translate excludes a distinct residue
unless \(4\mid h\). This gives its two displayed values.

Two scope distinctions are essential. First, Smilansky's
\(L_3,L_7\) in equation (6.3) are series over odd integers;
they omit the two-Euler factor. Their values
\(\pi/(2\sqrt3)\) and \(\pi/(2\sqrt7)\) must not be
substituted for the primitive values in (QP.8)--(QP.9).
Second, a general quadratic Euler indicator need not be an
ordinary principal-form indicator. For example \(\chi_{-20}\)
allows three, but three is not \(x^2+5y^2\). The general
theorem remains valid for its explicitly defined Euler set;
the five-form corollary uses the exact separate classifications.

## 8. Exact verification and source scope

The [finite checker](verify_quadratic_indicator_pairs.py) passed
327,680 direct representation/classification comparisons through
65,536, 121,656 certified prime-power residue memberships,
3,584 rational pair-density enclosures for shifts 1 through 128,
28 singleton enclosures and five exact shift-one normalization
checks. It also checks the character/order distinctions above,
the \(\chi_{-20}\) scope example, and 128 exact identities
for enlarging the bad set by inert two in the \(d=3\) case.
Unknown high-valuation residues are treated by lower and upper
rational bounds. Run
`python3 research/arithmetic/verify_quadratic_indicator_pairs.py`.

These are finite arithmetic checks. The analytic implication
uses the shared comparison theorem and the same scoped ZF/G
imports as before. The source comparisons, one-variable
centering, positive envelope and local model evaluation have
written internal application review. No independent validation
of the imported preprints, external human review, formal
verification, growing-shift uniformity or numerical onset is
inferred from these checks.

## 9. Mixed square and hexagonal indicator correlations

Let \(u_p(t)=t/p^{v_p(t)}\) on \(\mathbb Z_p\setminus\{0\}\)
and let \(\nu_2(t)\), \(\nu_3(t)\) be the indicators of
\(u_2(t)\equiv1\pmod4\), \(u_3(t)\equiv1\pmod3\), respectively.
The two mixed local masses, with normalized Haar measure, are
\[
 \begin{split}
 W_2(h)&=\int_{\mathbb Z_2}\nu_2(t)
                    \mathbf1_{2\mid v_2(t+h)}\,dt,\\
 W_3(h)&=\int_{\mathbb Z_3}\mathbf1_{2\mid v_3(t)}
                    \nu_3(t+h)\,dt .
 \end{split}                                                   \tag{QP.25}
\]
Null sets of infinite valuations can be discarded. These definitions
keep the square form at the unshifted endpoint.

**Theorem MixedSquareHexagonalPairs.** Assume the same scoped ZF and G
inputs. For each fixed nonzero integer \(h\), with \(\delta=1/32768\),
\[
 \sum_{1\le n\le X} b_1(n)b_3(n+h)
      = C_{1,3}(h)\frac X{\log X}
                  +O_h\!\left(\frac{X}{(\log X)^{1+\delta}}\right),
                                                               \tag{QP.26}
\]
where
\[
 C_{1,3}(h)=
 \frac{4\,3^{1/4}}{\sqrt2}\,W_2(h)W_3(h)
 \prod_{p\equiv5,7\ (12)}(1-p^{-2})^{-1/2}
 \prod_{\substack{p\mid h\\p\equiv11\ (12)}}
            \frac{1-p^{-v_p(h)-1}}{1-p^{-1}}.
                                                               \tag{QP.27}
\]
The finite-place masses are positive and satisfy
\[
 \begin{aligned}
 W_2(h)&=
 \begin{cases}
  1/4,& h\equiv1\pmod4,\\
  5/12,& h\equiv3\pmod4,\\
  1/2-W_2(h/2)/2,&2\mid h,
 \end{cases}
 &
 W_3(h)&=
 \begin{cases}
  1/4,& h\equiv1\pmod3,\\
  1/2,& h\equiv2\pmod3,\\
  1/2-W_3(h/3)/3,&3\mid h.
 \end{cases}
 \end{aligned}                                                 \tag{QP.28}
\]
In particular
\[
 C_{1,3}(1)
     =\frac{3^{1/4}}{4\sqrt2}
                    \prod_{p\equiv5,7\ (12)}(1-p^{-2})^{-1/2}.
                                                               \tag{QP.29}
\]

### Joint sieve envelope

Use (QP.23) for each indicator, enlarging both finite sets to
\(B=\{2,3\}\) with the original valuation rules. On the square
side, two has all exponents and three has even exponents.
On the hexagonal side, two has even exponents and three has
all exponents. Put \(\chi_1=\chi_{-4}\), \(\chi_3=\chi_{-3}\)
and \(I_d=\{p\notin B:\chi_d(p)=-1\}\) for \(d=1,3\).
Both sets of inert primes have density \(1/2\).
The two models \(D_d\) from (QP.14) satisfy (C1) and (C3)
separately with the same outer parameters. Choose a fixed depth
sufficient for both. For (C2), the identity class modulo 12
consists of split primes for both characters. A bank prime
changes only the small-prime factor and preserves both the
indicator and its model exactly.

For (C4), use the bank law (QP.16) after deleting primes
dividing \(h\), and the same exact-record progressions and
large-record bound as in Section 4. For an endpoint \(d\), use
a cutoff \(w\) when its weight is \(b_d\) and \(y\) when its
weight is \(D_d\). At a prime \(p\in I_d\) below that cutoff,
the endpoint's valuation cannot be one. Bank-record slopes
are units at these inert primes since the bank class is
\(1\bmod12\). Omit the finitely many good primes dividing \(h\)
for the upper bound. If \(k_p\in\{0,1,2\}\) endpoint tests
are active at \(p\), their forbidden classes are disjoint and
occupy a fraction \(\rho_p=k_p(p-1)/p^2\) modulo \(p^2\).
Since \(1-\rho_p=(1-1/p)^{k_p}(1+O(p^{-2}))\),
Mertens' theorem applied to the two inert sets gives the
nonbank sieve product
\[
 \prod_p(1-\rho_p)\asymp_h
                 \frac{1}{\sqrt{\log Y_1\log Y_3}},
 \qquad Y_d=\begin{cases}w,& w_d=b_d,\\ y,&w_d=D_d.\end{cases}
\]
Here \(w_d\) denotes the chosen endpoint weight. A model
weight supplies a factor \(O(L^{-3/10})\). With
\(z=X^{1/100}\) and \(w=z^{1/16}\), every choice of endpoint
weights therefore has sieve scale \(X/L\).
The Selberg denominator argument of Section 4 applies since the
forbidden fractions, including the bank atoms, are at most
\(2/p\). At a squarefree least common multiple \(m\) of the
sieve labels the CRT root count is at most \(2^{\omega(m)}m\).
Thus the endpoint error is still \(O(z^{4+o(1)})\) per retained
record, and summing over the \(E^2\) records costs less than
\(X^{1/2}\). The retained bank factors give exactly
\(\mu_{\mathcal P}(r,s)\), including the reciprocal record modulus.
Together with the large-record bound this yields (C4) for each
\(w_d\in\{b_d,D_d\}\), with absolute interval dependence through
the length. The shared comparison theorem gives a mixed
indicator-to-model replacement error \(O_h(XL^{-1-\delta})\)
on the dyadic interval, using the fixed-shift branch and modulus 12.

### Mixed model product and normalization

The cross-prime unit phases in (QP.3) equal one on the allowed
bad-prime masks. In fact \(\chi_1(3)=-1\) and \(\chi_3(2)=-1\)
are raised to even valuations. Thus the product of the two
bad profiles has Haar mass \(W_2(h)W_3(h)\).
Write \(K_d=K_{b_d}\) from (QP.6) and let \(A_d\) be the
product \(A_B\) from (QP.2) for \(b_d\) with enlarged \(B\).
By (QP.19), outside large smooth records the prefactor
on side \(d\) is
\[
 \lambda_{d,y}=\frac{2K_d}{A_d\sqrt L}
                           \prod_{\substack{p\le y\\p\in I_d}}(1+1/p)
\]
with relative error \(O(L^{-1/2})\).

To count the joint model support, use the square-prefix argument
of Section 5 with the prefix on side \(d\) supported on
\(I_d\cap[1,y]\). Retain both prefixes up to \(E\) and the
four bad-prime powers up to \(E\). The omitted integer count
is \(O_h(X/E+\sqrt X)\) by summing divisibility by squares and
bad-prime powers, and the omitted Haar mass is \(O(E^{-1})\).
In each retained progression the bad unit conditions are fixed.
The two quotient forms at good primes exclude at most one
root per inert endpoint, with an empty cell discarded if a
required quotient is identically zero there. The upper
dimension bound is one since the root fraction is at most
\((\mathbf1_{p\in I_1}+\mathbf1_{p\in I_3})/p\).
The fundamental lemma with level \(X^{1/100}\) has
\(\log(X^{1/100})/\log y\asymp L^{3/5}\).
The relative sieve error is smaller than every fixed log power.
The \(O(E^2(1+\log E)^4)\) retained cells, with at most
\(2^{\omega(m)}\) CRT roots at a squarefree sieve divisor \(m\),
have total endpoint error \(X^{1/50+o(1)}\).
Their main coefficients are the exact Haar masses of disjoint
cells and have sum at most one. Completing the omitted Haar
tails gives the joint mass
\[
 W_2(h)W_3(h)\prod_{5\le p\le y}\delta_{1,3;p}(h),
\]
where a prime inert on neither side has density one, a prime
inert on exactly one side has density \(p/(p+1)\), and a
prime inert on both sides has the pair density (QP.21).
The full smooth-record tails removed in using (QP.19) cost
\(XL^{O(1)}\exp(-cL^{1/10})\) on the two endpoints.

Multiply the densities by the two prefactors. At a singly
inert prime, \((1+1/p)\delta_{1,3;p}=1\). At a doubly inert
prime the two factors give \((1-p^{-2})\) times the
shift factor in (QP.27), which equals one for \(p\nmid h\).
Equation (QP.22) supplies the remaining constants
\[
 \begin{split}
 \frac{K_1^2}{A_1^2}\prod_{p\in I_1}(1-p^{-2})
 &=\frac{L(1,\chi_1)}{\pi}\,
   \frac1{1-1/2}\frac{1+1/3}{1-1/3}=1,\\
 \frac{K_3^2}{A_3^2}\prod_{p\in I_3}(1-p^{-2})
 &=\frac{L(1,\chi_3)}{\pi}\,
   \frac{1+1/2}{1-1/2}\frac1{1-1/3}=\frac{\sqrt3}{2}.
 \end{split}
\]
The Euler normalizations cancel the doubly inert factors
\(1-p^{-2}\). The primes inert at just one endpoint are
\(5,7\bmod12\), and those inert at both are \(11\bmod12\).
Thus the remaining absolutely convergent product and the
bad-place factors give (QP.27). Their omitted good-prime tail
as \(y\) grows has relative size \(O(y^{-1})\) once every prime
dividing \(h\) is included. The model sum on \(X<n\le2X\)
therefore has main term \(C_{1,3}(h)X/L\) and error
\(O_h(XL^{-3/2})\). Combining with the replacement and summing
over dyads down to \(\sqrt X\) proves (QP.26). The omitted
\(O(\sqrt X)\) terms and the \(O(X/L^2)\) error from summing
the dyadic main terms fit the asserted error.

### The two finite-place factors

For odd \(h\), the even \(t\) contribution at two has mass
\(1/4\), since \(t+h\) is odd. An odd \(t\) must be \(1\bmod4\);
this contributes nothing when \(h\equiv1\bmod4\), and contributes
the mass \(1/6\) of positive even valuations of \(t+h\) when
\(h\equiv3\bmod4\). For even \(h\), the odd \(t\) contribution
is \(1/4\). On the even part, division by two preserves
\(\nu_2\) and reverses the valuation parity test at \(t+h\).
This gives \(W_2(h)=1/4+(1/2)(1/2-W_2(h/2))\).

At three, for \(h\equiv1\bmod3\), the classes \(t\equiv0,2\)
contribute \(1/12\) and \(1/6\) respectively. The first
is the mass of positive even valuations of \(t\). The second
is the mass of positive valuations of \(t+h\) with unit
part \(1\bmod3\). For \(h\equiv2\bmod3\), the classes
\(t\equiv1,2\) contribute \(1/6\) and \(1/3\) respectively.
If \(3\mid h\), the \(t\equiv1\) class contributes \(1/3\).
Division by three on \(t\equiv0\) reverses the valuation parity
test and preserves \(\nu_3\) at the translated endpoint, giving
\(W_3(h)=1/3+(1/3)(1/2-W_3(h/3))\).
These computations prove (QP.28). Both integrals are at most
\(1/2\) and both recurrences and their initial values give a
lower bound \(1/4\). At \(h=1\) both masses equal \(1/4\).

In Iyer's notation the mixed count corresponds to
\(S(\triangle,\square_2,-h)\), by translating the starting
variable to \(n+h\). At \(h=1\), the hypothesis
\(\chi_{-3}(-1)=-1\) in Lemma 1.3 gives an unconditional
global lower bound of order \(X/(\log X)^{3/2}\).
Equation (QP.26) gives the global leading term under ZF and G.
Iyer's Corollary 1.4 gives a short-interval lower bound on
intervals of length \(H_a x^{5/6}(\log x)^{19}\).
[Iyer, Gaps between hexagonal and square lattices](https://doi.org/10.1017/S0004972726101191).
