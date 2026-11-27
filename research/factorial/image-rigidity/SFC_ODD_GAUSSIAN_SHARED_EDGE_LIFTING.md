# Higher compatibility at coincident odd Gaussian Hermite edges

The quadratic Gaussian window is zero-free at **every odd prime-power
start**. More generally, the dominant initial prime-power theorem in
[the odd-start note, Section 7](SFC_ODD_GAUSSIAN_LOCAL_WINDOW_TESTS.md#7-large-prime-powers-dividing-the-initial-order)
does not need a quadratic-character restriction on the prime.

Retain the same prospective root in the two unit-order equations.
Their leading binomials can coincide. A short polynomial combination
cancels that layer and the next one, leaving a single nonzero monomial.
The quadratic coefficient forced by the first order has too large a
valuation to cancel it.

## 1. Scope and existing reduction

Let \(G_m(A,B,C)=\mathbb E(A+BX_m+CX_m^2)^m\), with
\(X_m\sim\mathcal N(0,1/m)\). Fix \(q=(2k+1)p^e\), where
\(k\ge0,e\ge1\), and \(p\) is an odd prime with \(p>4k+2\).
Define
\[
 S_k(Z)=\sum_{j=0}^k\frac{k_{\underline j}}{(2j+1)!}Z^j,\qquad
 T_k=\operatorname{Res}_Z((2k+1)!S_k(Z),(2k+1)!S_k(2Z)),\quad T_0=1.
\]

**Prime-power window theorem.** For every \(p>4k+2\), the forms
\(G_q,G_{q+1},G_{q+2}\) have no common projective zero over any
characteristic-zero field. In particular this holds at every odd
prime-power start \(q=p^e\).

The [final-residue refinement](SFC_ODD_GAUSSIAN_FINAL_RESIDUE_LIFTING.md)
proves the final-slope conclusion for every exponent and both
quadratic characters, without the resultant condition. The argument below proves the
internal-edge exclusion independently of that refinement.

Use the initial-order reduction proved in Section 7 of the linked
note. A hypothetical algebraic zero, embedded in
\(\overline{\mathbb Q}_p\), has the chart and valuations
\[
 B=1,\quad x=A/B,\quad h=C/B,\quad v(p)=1,\quad v(x)>0,
 \qquad v(h)\ge b=e/2-L/2,\quad L=\sum_{j=1}^e p^{-j}.       \tag{1}
\]
The even-order polygon forces \(v(x)<v(h)\). The two unit-order
pure Hermite equations have distinct first slopes. Their final common
slope, when \(k>0\), is excluded by \(p\nmid T_k\). The remaining
possible common valuations are
\[
 \alpha=\frac{J+1}{2J(p-1)},\qquad J=p^a,\quad1\le a<e.       \tag{2}
\]
A quadratic nonresidue \(2\) already separates their leading
binomials. The case requiring a new argument therefore has
\(e\ge2,(2/p)=1\), and \(p\ge7\). The preceding reduction and
the optional final-slope resultant test are unchanged inputs to
the earlier shortcut; the final-residue refinement now closes
that slope even when the resultant vanishes modulo \(p\).

Put \(N=(q+1)/2\). The normalized pure even and odd polynomials are
\[
 E(Y)=\sum_{i=0}^NE_iY^i,\quad E_i=\frac{(2(q+1))^iN_{\underline i}}{(2i)!},
 \qquad O(Y)=\sum_{i=0}^NO_iY^i,\quad O_i=\frac{(2(q+2))^iN_{\underline i}}{(2i+1)!}. \tag{3}
\]
The pure moments are nonzero row constants times \(E(x^2)\) and
\(xO(x^2)\), respectively.

## 2. The complete next-layer inventory

Fix (2) and set
\[
 \lambda=2\alpha,\quad r=(J+1)/2,\quad r'=(pJ+1)/2,\quad s=r+J,
 \quad\delta=2/(p-1),
\]
\[
 \gamma=v(E_r)+r\lambda=-a/2+\frac{3J+1}{2J(p-1)},\qquad
 \omega=\gamma+\delta.                                      \tag{4}
\]
Give a coefficient \(c_iY^i\) the weight \(v(c_i)+i\lambda\).

**Near-face lemma.** The following are the complete lists of terms
of weight at most \(\omega\). Include the constant in the first row
only when \(0\le\omega\).

| Polynomial | Indices and weights |
|---|---|
| \(E\) | \(r,r':\gamma;\ s:\omega;\ 0:0\) |
| \(YO\) | \(r,r':\gamma;\ r+1,r'+1:\gamma+\lambda;\ s:\omega\) |
| \(YE\) | \(r+1,r'+1:\gamma+\lambda\) |

Legendre's formula and the doubling carry in the factorial quotient give
\[
 v(E_i)=-v(i!)-\sum_{j=1}^e\mathbf1_{i\bmod p^j=(p^j+1)/2},
 \qquad v(O_i)=v(E_i)-v(2i+1).                              \tag{5}
\]
Indeed the carry from \(\binom Ni\) at level \(p^j\) is
\(\mathbf1_{i\bmod p^j>(p^j+1)/2}\), and the carry from
\((2i)!/(i!)^2\) is the same indicator with \(\ge\).
No level above \(e\) occurs because \(2N=q+1<p^{e+1}\).

Write \(c=(p-1)/2\). The sum of indicators is the length \(t\),
truncated at \(e\), of the initial base-\(p\) digit string
\((c+1,c,\ldots,c)\), read from the least significant digit.
With \(s_p(i)\) the digit sum, the even weight is
\[
                       (s_p(i)+i/J)/(p-1)-t.                \tag{6}
\]
For a fixed positive prefix length, its minimum occurs at
\(N_t=(p^t+1)/2\); adding higher digits increases both terms in
the numerator. The exact gap there is
\[
 W(N_t)-\gamma=(a-t)/2+\frac{p^t-J}{2J(p-1)}.
\]
For \(1\le t<a\), its minimum is at \((p^t+1)/2\), and its
gap above \(\gamma\) is at least \((p-1)/(2p)>\delta\).
For \(t=a\), write \(i=r+Ju\). The gap is
\((s_p(u)+u)/(p-1)\), retaining only \(u=0,1\).
For \(t=a+1\), write \(i=r'+pJu\); the gap is
\((s_p(u)+pu)/(p-1)\), retaining only \(u=0\).
For \(t\ge a+2\), the minimum gap is at least
\((p-1)/2>\delta\). For \(t=0\), every positive index has weight
at least \(\lambda>\omega\). Only the constant may remain.
Here \(\lambda-\omega=a/2-(5J-1)/(2J(p-1))>0\) for \(p\ge7\).
This proves the first row, including coefficients away from the hull.

For \(YO\), put \(j=i-1\). If \(v(2j+1)=0\), its weight is
the even weight at \(j\) plus \(\lambda\), retaining just \(j=r,r'\).
Otherwise write \(j=(p^t-1)/2+p^tu\), with \(t=v(2j+1)\).
The weight becomes
\[
 -t/2+\frac{s_p(u)+(p^t-1)/(2J)+p^tu/J}{p-1}+\lambda.
\]
At \(u=0\) this is the even minimum for prefix length \(t\).
The same four cases retain \(i=r,s,r'\). This proves the second
row. Shifting the first row by one proves the third.

## 3. Two cancellations leave a monomial

Define rational units
\[
 \kappa=E_r/O_{r-1},\qquad c_*=E_{r+1}/E_r-O_r/O_{r-1},
 \qquad R(Y)=(1-c_*Y)E(Y)-\kappa YO(Y).                     \tag{7}
\]
The coefficient at \(r+1\) is exactly zero. The coefficient ratios
\[
 \frac{E_i}{O_{i-1}}=\frac{(q+1)(N-i+1)}i
       \left(\frac{q+1}{q+2}\right)^{i-1},\qquad
 \frac{O_i}{E_i}=\frac1{2i+1}\left(\frac{q+2}{q+1}\right)^i \tag{8}
\]
give \(\bar\kappa=2^{2-r}\) and \(\bar c_*=-2\). Here
\(E_{r+1}/E_r\) has positive valuation and \(O_r/O_{r-1}\equiv2\).

Since \((2/p)=1\), \(2^{r'-r}=2^{(p-1)J/2}\equiv1\pmod p\).
Equations (8) cancel the weight-\(\gamma\) coefficients at \(r,r'\)
modulo one further power of \(p\). They also cancel the layer at
\(r+1,r'+1\), the left cancellation being exact. One valuation
step suffices because \(\delta<1\).

At \(s=r+J\), instead
\(\kappa O_{s-1}/E_s\equiv2^J\equiv2\pmod p\). The \(YE\)
term is above the retained layer there. Hence \(R_s/E_s\equiv-1\),
with weight exactly \(\omega\). These are all possible low terms
by the near-face lemma. Thus the least weight of \(R\) is
\(\min(0,\omega)\), achieved by a single monomial: the constant
\(1\) or the term of index \(s\). A tie is impossible since
\[
 2J(p-1)\omega=J(7-a(p-1))+1\ne0                            \tag{9}
\]
by reduction modulo \(p\).

## 4. The quadratic perturbation cannot cancel it

At a unit order \(m\), let \(a_l\) be the coefficient of \(x^l\)
in \(\mathbb E(x+X_m)^m\). A nonzero coefficient of \(x^ih^j\)
in \(\mathbb E(x+X_m+hX_m^2)^m\) is
\[
 a_{i+j}\binom{i+j}{j}\prod_{u=1}^j\frac{m-i-j+2u-1}{m}.    \tag{10}
\]
The ratio is \(p\)-integral. For \(v(x)=\alpha<v(h)=\beta\),
every term with \(j\ge1\) lies at least \(\beta-\alpha\) above
the minimum pure weight. This remains true after row normalization.

At a common quadratic zero, the two unit-order equations imply
\[
 v(E(x^2))\ge\gamma+\beta-\alpha,\qquad
 v(x^2O(x^2))\ge\gamma+\beta-\alpha.
\]
Since \(1-c_*x^2\) is a unit, (7) gives
\(v(R(x^2))\ge\gamma+\beta-\alpha\). Its unique least monomial
instead gives \(v(R(x^2))=\min(0,\omega)\). The bounds conflict:
if the minimum is zero, (1) and (4) give
\[
 \gamma+b-\alpha=(e-a)/2+1/(p-1)-L/2>0.
\]
If the minimum is \(\omega\), the gap is at least
\[
 b-\alpha-\delta=e/2-L/2-\frac{5+1/J}{2(p-1)}>0,
\]
using \(p\ge7,e\ge2,J\ge p,L<1/(p-1)\). Both margins exceed
\(5/12\). A zero quadratic coefficient is excluded by the same
argument with infinite valuation. All shared edges (2) are now
excluded, proving the prime-power window theorem.

## 5. Example, evidence and remaining scope

**Factorial consequence.** For every fixed homogeneous bivariate
trinomial support and every start covered by the prime-power window
theorem, its factorial moments at \(q,q+1,q+2\) cannot all vanish
once both common heights exceed a computable support- and
start-dependent threshold. This follows from the
[exact Gaussian transfer](SFC_BIVARIATE_TRANSLATION_GAUSSIAN_JETS.md#1-exact-graded-normalization).
It includes all odd prime-power starts and makes no assertion of a
cutoff uniform over starts or supports.

At start \(49\), prime \(7\), the shared valuation is
\(\alpha=2/21\), with \(b=45/49,\gamma=-5/21\). Here
\(\kappa=34375000/132651\) and \(c_*=-47/6\). The polynomial
(7) has constant \(1\) as its unique least term. A quadratic
correction has weight at least \(86/147>0\), so it cannot vanish.
The old leading binomials coincide at this same prime.

The [checker](verify_sfc_odd_gaussian_shared_edges.py) checks the
complete layer inventories and cancellations in twelve selected
cases, including prime-power depth five and nontrivial dominant
cofactors. It checks 20 shared edges and 1498 coefficients in (10),
including both unit orders at start 49, and verifies that example's
full syzygy separately over the rationals. It rejects a corrupted \(\kappa\),
and retains the old leading collision as a control. The
[receipt](sfc_odd_gaussian_shared_edge_verification.json) records
the inputs and outcomes. Run from the repository root:

```sh
python3 research/factorial/image-rigidity/verify_sfc_odd_gaussian_shared_edges.py
```

These are internal exact controls. Unrestricted primes and exponents
are covered by the digit proof, not by extrapolating the finite cases.
There is no independent review, external replay or formal verification.

The pure forms agree with [DLMF 18.12.16](https://dlmf.nist.gov/18.12.E16)
after imaginary-argument and variance scaling. That classical
generating function does not exclude coincident scaled initial forms.
The valuation method is classical; no historical priority claim is
made. The final-slope resultant argument retains its separate Schur
input as documented in the original note.

The unrestricted odd-start problem remains open outside the stated
dominant initial and middle prime-power families and verified finite
domain. The character restriction on the initial shared edges is
removed and is no longer a remaining obstruction.

The final-slope family \(q=5\cdot13^e\) is excluded for every
\(e\ge1\) by the [same-residue compatibility proof](SFC_ODD_GAUSSIAN_FINAL_RESIDUE_LIFTING.md#5-the-requested-family-and-the-next-boundary).
At the common residue \(Z=2\), the next corrections disagree.
A profile requiring the second secondary layer is
\(q=9\cdot89^e\), \(e\ge2\):
the reductions of \(S_4(Z),S_4(2Z)\) have gcd \(Z+25\).
Here \((2/89)=1\), and the first secondary compatibility vanishes
at \(Z=64\). The [next layer](SFC_ODD_GAUSSIAN_FINAL_RESIDUE_LIFTING.md#6-the-second-secondary-layer-at-positive-character-primes)
has normalized odd/even ratio two, forcing Wronskian valuation
\(46/89<1\) against a lower bound of one. This excludes every
exponent in that family and removes the remaining dominant
initial-prime final-residue exceptions. Unrestricted odd starts
remain separate.
