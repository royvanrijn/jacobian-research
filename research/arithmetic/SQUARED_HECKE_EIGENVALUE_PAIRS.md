# Squared Hecke-eigenvalue pairs and the prescribed shift one

**Mathematical status.** The asymptotic for the full squared-eigenvalue
correlation at the prescribed shift one remains **open in this work**.
This record gives an unconditional truncated-pair asymptotic, its explicit
positive constant, an exact signed remainder criterion, and the scope of a
weighted marked-graph extension. The truncated theorem follows from
classical Hecke identities, Deligne's bound and Miller's cancellation
estimate; no historical novelty is claimed for that consequence.
The full target is not a consequence of the existing Dirichlet-only input.

The first target is the normalized Ramanujan form
\[
 \Delta(z)=\sum_{n\ge1}\tau(n)e(nz),\qquad
 \lambda_\Delta(n)=\tau(n)n^{-11/2},\qquad
 a(n)=|\lambda_\Delta(n)|^2.
\]
The arithmetic statements below also apply to two fixed normalized
holomorphic Hecke cusp forms of full level, with real eigenvalues and
trivial central character. Write \(a_f(n)=\lambda_f(n)^2\), and set
every sequence to zero on nonpositive integers. No ramified-level,
Maass-form or general varying-shift extension is asserted here.

## 1. Exact relation to the inspected literature

Walsh's [12 November revision, Section 1.2.4](https://arxiv.org/html/2610.40439v2#S1.SS2.SSS4),
printed pp.5–6, identifies squared holomorphic Hecke-eigenvalue
correlations as an extension to pursue. That passage does not prove a
weighted extension of its graph argument.

Kim's [arXiv:2011.06142v5](https://arxiv.org/abs/2011.06142v5),
Theorem 1.1, printed p.3, treats a fixed full-level holomorphic form.
For \(X^{2/3+\varepsilon}\le H\le X^{1-\varepsilon}\), its
dyadic correlation has error \(O(X(\log X)^{-A})\) outside
\(O(H(\log X)^{-3A})\) shifts in \([1,H]\). The exceptional
set may contain one at every scale. Equation (1.9) gives its predicted
coefficient. We do not claim its range is the best currently available:
a later Huang–Lü paper, [JNT 278 (2026), 128–177](https://doi.org/10.1016/j.jnt.2025.04.017),
also concerns almost-all shifts, but its full theorem ranges were not
verified in this audit.

Mangerel's [arXiv:2108.11401v2](https://arxiv.org/abs/2108.11401v2),
Section 1.3.1 and Corollary 1.1, gives one-variable short-interval
variance estimates. Lemma 5.6 supplies the harmonic prime mass for
\(|\lambda_f(p)|^2\), with the holomorphic Ramanujan hypotheses
satisfied. These are useful inputs, with a different conclusion from a
prescribed-shift binary asymptotic.

The exact cached primary editions and scoped imported statements are
registered with this note. Claims of exhaustive priority or an exhaustive
classification of the literature are separate from this source comparison.

## 2. The Hecke convolution and its one-variable cancellation

The Hecke relation gives the exact identities
\[
 c_f(n):=\lambda_f(n^2),\qquad
 a_f(n)=\sum_{d\mid n}c_f(d),\qquad c_f=\mu*a_f.       \tag{H.1}
\]
Deligne's bound implies
\[
 |\lambda_f(p)|\le2,\quad |c_f(n)|\le d(n^2)\le d_3(n),
 \quad \sum_{n\le Z}|c_f(n)|\ll Z\log^2(2Z).           \tag{H.2}
\]
Put
\[
 R_f=\frac{L(1,\operatorname{sym}^2f)}{\zeta(2)}>0.
\]
For \(\Re s>1\),
\[
 \sum_n\frac{a_f(n)}{n^s}
   =\zeta(s)\frac{L(s,\operatorname{sym}^2f)}{\zeta(2s)},
 \qquad
 C_f(s):=\sum_n\frac{c_f(n)}{n^s}
   =\frac{L(s,\operatorname{sym}^2f)}{\zeta(2s)}.        \tag{H.3}
\]
Here \(R_f\) is the usual Rankin–Selberg mean of \(a_f\).
Kim's actual Lemma 1.3 records the consequence of
[Miller, Theorem 1.1](https://arxiv.org/abs/math/0404521v3)
\[
 \sum_{n\le T}c_f(n)e(\alpha n)
                 \ll_{f,\varepsilon}T^{3/4+\varepsilon}
 \quad(\alpha\in\mathbb R).                            \tag{H.4}
\]
The conversion from symmetric-square coefficients uses the square-Möbius
factor in (H.3); its absolute budget is
\(\sum_{d\ge1}d^{-3/2-2\varepsilon}<\infty\).
At \(\alpha=0\), partial summation and Abelian passage from
\(\Re s>1\) therefore give
\[
 B_f(T):=\sum_{n\le T}\frac{c_f(n)}n
       =R_f+O_{f,\varepsilon}(T^{-1/4+\varepsilon}).    \tag{H.5}
\]
In particular these partial sums are uniformly bounded. None of
(H.1)–(H.5) uses the recent ZF or G preprints.

## 3. The positive normalized constant

Write \(A_{f,p}=\lambda_f(p)^2\). The actual local Euler factor is
\[
 F_{a_f,p}(s):=\sum_{k\ge0}a_f(p^k)p^{-ks}
   =\frac{1+p^{-s}}
    {(1-p^{-s})(1-(A_{f,p}-2)p^{-s}+p^{-2s})}.
\]
For a twisted local factor, replace \(p^{-s}\) throughout by
\(\eta(p)p^{-s}\). Define
\[
 m_{f,p}=(1-p^{-1})\sum_{k\ge0}\frac{a_f(p^k)}{p^k}
       =\frac{1+p^{-1}}
              {1-(A_{f,p}-2)p^{-1}+p^{-2}}.            \tag{H.6}
\]
The second formula follows from the prime-power Hecke recurrence. For
the shift one, at most one endpoint is divisible by \(p\). Hence
the raw \(p\)-adic pair mean is
\[
 H_{f,g,p}=1-2/p+(1-1/p)
       \sum_{k\ge1}\frac{a_f(p^k)+a_g(p^k)}{p^k}
             =m_{f,p}+m_{g,p}-1.                       \tag{H.7}
\]

**Proposition HeckePairConstant.** For fixed full-level holomorphic
forms \(f,g\), the normalized expression
\[
 C_{f,g}(1)=R_fR_g\prod_p J_{f,g,p},\qquad
 J_{f,g,p}=\frac{m_{f,p}+m_{g,p}-1}{m_{f,p}m_{g,p}}     \tag{H.8}
\]
defines a finite strictly positive constant; the displayed product is
absolutely convergent. This is an arithmetic normalization and is not
an assertion of the full shifted-pair asymptotic.

Indeed \(m_{f,p}=1+(A_{f,p}-1)/p+O(p^{-2})\), uniformly
in the Deligne range. Thus \(J_{f,g,p}=1+O(p^{-2})\).
Each raw mean in (H.7) is positive. For \(p>2\) its unit-unit
term is already positive. At two, the terms \(a_f(2^k)\)
cannot all vanish, by the Hecke recurrence. Absolute convergence and
positivity of the individual factors prove the assertion.

For \(f=g=\Delta\), (H.8) becomes the particularly simple formula
\[
 C_{\Delta,\Delta}(1)=R_\Delta^2
 \prod_p\left[1-
  \left(\frac{\lambda_\Delta(p)^2-1-p^{-1}}{p+1}\right)^2
       \right].                                       \tag{H.9}
\]
For example, \(\tau(2)=-24\) gives
\(A_{\Delta,2}=9/32\), \(m_{\Delta,2}=32/45\), and the
normalized factor at two is exactly \(855/1024\).

The local interpretation can also be made without a growing-prime
limit inside an integer sum. For fixed \(y\), put
\(D_{f,y}(n)=R_f\prod_{p\le y}a_f(p^{v_p(n)})/m_{f,p}\).
CRT, truncation of the finitely many valuations and
\(a_f(p^k)\le(k+1)^2\) give
\[
 \lim_{X\to\infty}\frac1X\sum_{n\le X}
       D_{f,y}(n)D_{g,y}(n+1)
       =R_fR_g\prod_{p\le y}J_{f,g,p}.                 \tag{H.10}
\]
The subsequent \(y\to\infty\) limit is (H.8). Transferring
this local model to the full coefficients is a separate obligation.

A graph application needs a model preserving every active bank record.
For a bank prime greater than \(y\), the multiplier of the model
in (H.10) is one, whereas the actual coprime multiplier is \(a_f(p)\).
Thus its centered remainder cannot be substituted directly into the
weighted marked comparison. The finite-prime model here identifies the
arithmetic coefficient; the graph model must retain the bank factors.

## 4. An unconditional truncated-pair asymptotic

For \(1\le Z\le X\), define the generally signed truncation
\[
 a_{f,Z}(n)=\sum_{\substack{d\mid n\\d\le Z}}c_f(d).
\]

**Theorem HeckeTruncatedPairs.** For fixed forms \(f,g\) as above,
\[
 \sum_{n\le X}a_{f,Z}(n)a_{g,Z}(n+1)
   =C_{f,g}(1)X+
       O_{f,g,\varepsilon}\bigl(XZ^{-1/4+\varepsilon}
                                  +Z^2\log^4(2Z)\bigr).       \tag{H.11}
\]
In particular \(Z=X^{4/9}\) gives error
\(O_{f,g,\varepsilon}(X^{8/9+\varepsilon})\).

To prove this, expand the two finite divisor sums. CRT gives
\[
 \sum_{n\le X}a_{f,Z}(n)a_{g,Z}(n+1)
       =XM_Z+O(Z^2\log^4(2Z)),\quad
 M_Z=\sum_{\substack{d,e\le Z\\(d,e)=1}}
                         \frac{c_f(d)c_g(e)}{de},       \tag{H.12}
\]
using (H.2) for the sum of the \(O(1)\) progression errors.
The convergence of this signed double sum requires an argument.

Let
\[
 C_{f,p}(z)=\sum_{k\ge0}c_f(p^k)z^k
       =\frac{1+z}{1-(A_{f,p}-2)z+z^2},\qquad
 U_{f,p}(z)=1-C_{f,p}(z)^{-1}
       =\frac{z(A_{f,p}-1-z)}{1+z}.
\]
All positive-degree coefficients of \(U_{f,p}\) have absolute
value at most four. The coprime double Dirichlet series factors as
\[
 \sum_{(d,e)=1}\frac{c_f(d)c_g(e)}{d^se^t}
       =C_f(s)C_g(t)\mathcal H_{f,g}(s,t),\quad
 \mathcal H_{f,g}(s,t)=\prod_p
            \{1-U_{f,p}(p^{-s})U_{g,p}(p^{-t})\}.       \tag{H.13}
\]
Its double Dirichlet coefficients \(h(k,l)\) satisfy
\[
 \sum_{k,l}\frac{|h(k,l)|}{k^\sigma l^\tau}<\infty
 \quad(\sigma,\tau>0,\ \sigma+\tau>1).                 \tag{H.14}
\]
This follows directly by bounding the local positive-degree double
sum by \(16p^{-\sigma-\tau}/
((1-p^{-\sigma})(1-p^{-\tau}))\).

Set \(B_f(T)=0\) for \(T<1\). Convolution in (H.13) gives
\[
 M_Z=\sum_{k,l}\frac{h(k,l)}{kl}
                    B_f(Z/k)B_g(Z/l).
\]
For \(0<\delta<1/4\), (H.5), uniform boundedness of the
partial sums, and (H.14) imply
\[
 M_Z=R_fR_g\mathcal H_{f,g}(1,1)+O(Z^{-\delta}).
\]
Explicitly, the error is bounded by a constant times
\(Z^{-\delta}\sum|h(k,l)|(k^\delta+l^\delta)/(kl)\),
which converges. At \((1,1)\), the factor in (H.13) is exactly
\(J_{f,g,p}\). Taking \(\delta=1/4-\varepsilon\) proves
(H.11). The argument uses one-variable cancellation and elementary
progression counting; it does not estimate the full binary remainder.

## 5. The exact remaining signed estimate

Define
\[
 r_{f,Z}(n)=\sum_{\substack{d\mid n\\d>Z}}c_f(d),
 \qquad a_f=a_{f,Z}+r_{f,Z},
\]
and, with \(Z=X^{4/9}\),
\[
 \begin{split}
 \mathcal R_{f,g}(X)=\sum_{n\le X}\{&
       r_{f,Z}(n)a_{g,Z}(n+1)
       +a_{f,Z}(n)r_{g,Z}(n+1)\\
       &+r_{f,Z}(n)r_{g,Z}(n+1)\}.
 \end{split}                                           \tag{H.15}
\]

**Proposition HeckeLargeDivisorReduction.** For these fixed forms,
the desired shift-one asymptotic with coefficient (H.8) holds if and
only if \(\mathcal R_{f,g}(X)=o(X)\).

This follows from the exact expansion and (H.11). The case
\(f=g=\Delta\) is the current open estimate. Its summands are
signed, even though the full \(a_f\) is nonnegative. Replacing
(H.15) by termwise absolute values is an additional, stronger request.
The triangle bound for the full divisor expansion is on the scale
\((\sum_{d\le X}|c_f(d)|)(\sum_{e\le X+1}|c_g(e)|)\),
which is much too large; (H.4) alone does not control the coupled
divisibility constraints in (H.15). This reduction makes the unresolved
estimate explicit rather than including it in a model hypothesis.

## 6. What changes in the graph argument

**Proposition HeckePrimeInvarianceObstruction.** The whole sequence
\(a(n)=\tau(n)^2/n^{11}\) fails the current common-prime
invariance condition at every prime.

Indeed \(a(p)=1\) would imply \(\tau(p)^2=p^{11}\),
an equality between an integer square and a nonsquare. Thus
\(a(p\cdot1)\ne a(1)\) for every prime. Also
\(a(p^2)=(a(p)-1)^2\), so repeated powers cannot be replaced by
complete multiplicativity. No assumption that \(\tau(p)\ne0\)
is made; division by individual eigenvalues is unnecessary.

For a finite bank union \(\mathcal P\), however, the exact
factorization \(a(n)=a(n_\mathcal P)a(n/n_\mathcal P)\)
leaves an invariant rough factor. The varying record coefficient can
therefore be treated separately. The
[weighted marked-comparison calculation](HECKE_WEIGHTED_BANK_COMPARISON.md)
uses the original graph label law and normalizes the physical vertex
weights. For a self-correlation, set
\[
 b(p)=a(p),\quad V_g=\sum_{p\in\mathcal P_g}1/p,\quad
 W_g=\sum_{p\in\mathcal P_g}b(p)^2/p,\quad R_g=W_g/V_g.
\]
With bounded invariant factors, normalization by \(R_g^{-M/2}\)
per endpoint keeps the physical squared norm bounded independently of
the eventual fixed mark length \(M\). Shared-label compression
leaves the explicit comparison factors \(\prod_gR_g^{1-t_g}\).
They cost only a fixed logarithmic power when the number of banks is
\(O(\log\log X)\). This result retains the stated G and ideal
comparison imports; it is not an application of the whole-weight
invariance theorem to \(a\).

Positive bank mass already follows from known input: if
\(\sum_{p\in\mathcal P_g}b(p)/p\asymp V_g\), then
Cauchy–Schwarz and \(0\le b(p)\le4\) imply
\[
 \frac{(\sum b(p)/p)^2}{V_g}\le W_g\le4\sum b(p)/p.
\]
Individual small or zero prime coefficients cause no singularity.

The remaining analytic prime sum is
\[
 \sum_{p\in I}b(p)\eta(p)p^{it}
     =\sum_{p\in I}\eta(p)p^{it}
       +\sum_{p\in I}\lambda_\Delta(p^2)\eta(p)p^{it}. \tag{H.16}
\]
The second term belongs to twisted symmetric-square automorphic
theory. Harcos–Thorner's
[arXiv:2303.16889v6, Theorem 2.1](https://arxiv.org/abs/2303.16889v6)
supplies a genuine fixed-representation Siegel–Walfisz input for
Rankin–Selberg prime coefficients in polylogarithmic progressions.
Its precise hypotheses and pole terms matter. A low-frequency prime
progression theorem does not settle the full height range in (H.16).
A standard region of width \(c/\log(q(|t|+3))\) alone produces
no fixed logarithmic saving for subpower bank primes at polynomial
\(X\)-height. Such a bound, or a sufficient integrated substitute,
must be established at the exact graph comparison scales.

There is a useful distinction between rough means and prime sums.
For a primitive Dirichlet character \(\eta\), the twisted
coefficient series is
\[
 \sum_n a_\Delta(n)\eta(n)n^{-s}
    =\frac{L(s,\eta)L(s,\operatorname{sym}^2\Delta\otimes\eta)}
            {L(2s,\eta^2)}\qquad(\Re s>1).
\]
For a character induced to a larger modulus, delete the Euler factors
at every additional prime of that modulus on both sides. This convention
is necessary when expanding periodic tests into imprimitive characters.
The principal pole has order one. Higher Taylor terms of its Hankel
density carry \(1/\Gamma(1-j)=0\) for \(j\ge1\); its
principal rough density is therefore constant. Zeros of the entire
symmetric-square numerator do not obstruct a coefficient-mean contour
shift. The prime estimate instead concerns a logarithmic derivative,
where those zeros matter. Bounds uniform in the required twists and
heights still have to be proved in either application.

The current obligations are weighted tails for the unbounded invariant
factors, a positive bank-record envelope, the required weighted frequency
estimate, and attachment of the rough-mean calculation below to the full
record-preserving graph model. The marked
normalization and the arithmetic constant do not imply those estimates.
The general \(f,g\) case also needs a separate audit of shared
cross-masses; the self-correlation normalization alone supplies no such
extension.

## 7. A rough-mean contour calculation with an explicit analytic input

This paragraph isolates the automorphic assertion used in the coefficient
mean; it is separate from the prime estimate (H.16). Let **(AUT)** mean:
for every primitive character \(\eta\) of conductor \(q\), the
coefficient product
\[
 Q_\eta(s)=\frac{L(s,\eta)L(s,\operatorname{sym}^2 f\otimes\eta)}
                         {L(2s,\eta^2)}
\]
is meromorphic on \(\Re s\ge3/4\), has only the simple pole
at one for the principal primitive character, of residue \(R_f\),
and satisfies, for some fixed \(C_f\),
\[
 |(s-1)^{\mathbf1_{\eta=1}}Q_\eta(s)|
             \ll_f \{q(3+|\Im s|)\}^{C_f}
 \quad(3/4\le\Re s\le2).                              \tag{H.17}
\]
Use the corresponding polynomial bound for \(L(s,\eta)\).
These are the standard automorphic continuation and growth inputs;
Kim's proof of Lemma 3.7 uses the associated twisted convexity bound.
The implication proved here retains (AUT) explicitly, rather than treating
it as part of the imported Dirichlet-only ZF assertion.

Put \(L=\log X\), \(y=e^{L^{2/5}}\),
\(V_y=\prod_{p\le y}(1-1/p)\),
\(\rho_y=R_f\prod_{p\le y}F_{a_f,p}(1)^{-1}\), and
\(\kappa_y=\rho_y/V_y\). For a bounded periodic test \(b\)
of period at most \(L^B\), and a smooth weight \(W\) supported
on \([T,2T]\), with fixed derivative budgets bounded by powers of
\(L\), (AUT) implies, uniformly for \(X^{2/5}\le T\le X^2\),
\[
 \sum_{\substack{n\ge1\\(n,\prod_{p\le y}p)=1}}
       (a_f(n)-\kappa_y)b(n)W(n)
       \ll T\exp\{-c\log T/\log y\}L^{C'}.            \tag{H.18}
\]
The constants may depend on the fixed form, period exponent and test
budgets. In particular the bound is smaller than every fixed negative
power of \(L\) times \(T\).

To check this implication, expand the periodic test into characters on
the units modulo its period; all nonunits are absent from the rough
sum for sufficiently large \(X\). Inflate primitive characters with
the Euler-factor convention above. On
\(\sigma=1-1/\log y\), the inverse small-prime factors obey
\[
 \prod_{p\le y}|F_{a_f,\eta,p}(\sigma+it)^{-1}|
 \le\exp\{O(\sum_{p\le y}p^{-\sigma})\}
                                      \ll(\log y)^{O(1)}.
\]
This uses the exact rational local factor and \(A_{f,p}\in[0,4]\).
The same bound holds for deleting the small Euler factors of ordinary
\(L(s,\eta)\). Mellin inversion and (H.17), with sufficiently
many integrations by parts in \(\widehat W\), permit moving the
contour to \(\sigma\). The new integral is bounded by
\(T^\sigma L^{C'}\). The only residue is principal: \(\rho_y\)
for the actual coefficients and \(V_y\) for the rough constant
sequence. Their residues agree after multiplication by \(\kappa_y\).
The principal coefficient of the periodic test is its average on unit
residue classes, \(\varphi(q)^{-1}\sum_{(a,q)=1}b(a)\), in
both sums. This proves (H.18).

For exact small records \(r\le e^{\sqrt L}\), their multiplier
has the crude bound \(a_f(r)\le d(r)^2\le r^2\). This costs
at most \(e^{2\sqrt L}\), still dominated by the exponential
saving in (H.18). Thus the one-variable estimate has the required
recordwise margin. The remaining truncation and bank-envelope arguments,
and the full weighted frequency estimate, are separate from this lemma.

## 8. Verification scope

The companion [exact checker](verify_squared_hecke_pairs.py) checks finite Hecke,
convolution, local-density and CRT identities using independently generated
integer coefficients of \(\Delta\). Its finite checks do not test
or prove (H.15), a weighted frequency estimate, or a full pair asymptotic.
The written arguments and source comparisons have internal AI review.
No external human review, independent whole-proof replay, formal
verification, prescribed numerical onset or historical priority is claimed.
