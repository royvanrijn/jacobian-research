# Well-factorable weights and the critical modulus-12 remainder

This note proves a support obstruction for pure well-factorable lower
weights, gives the exact large-prime branch of the half-dimensional
Rosser remainder, and records the applicable range of the primary
Fouvry theorem. It does not prove the critical signed remainder or the
every-target \(C_{12}+C_{12}\) assertion.

## 1. Sequence and the unresolved estimate

Include 1 in
\[
 C=\{n:p\mid n\Rightarrow p\equiv1\pmod{12}\},\qquad
 A=\{n:p\mid n\Rightarrow p\equiv1\pmod3\}.
\]
For \(N\equiv2\pmod{12}\), put
\[
 U=\tau^2\mathbf1_C,\quad V=\tau\mathbf1_A,\quad
 w_N(b)=U(N-b)V(b),\quad X_N=N\mathfrak S_{\rm asym}(N).
\]
Use the exact local density
\[
 g_N(p)=
 \begin{cases}
 (2p-1)/(p^2-p+1),&p\equiv7\ (12),\ p\nmid N,\\
 0,&\text{otherwise},
 \end{cases}
\]
extended multiplicatively to squarefree integers, and write
\[
 A_d(N)=\sum_{\substack{1\le b<N\\d\mid b}}w_N(b),\qquad
 r_d(N)=A_d(N)-X_Ng_N(d).
 \tag{1}
\]
The canonical definitions and all-valuation local calculation are in
[AsymSemilinear](ASYMMETRIC_SEMILINEAR_SIEVE.md).
In particular \(r_d=0\) when a prime \(7\bmod12\) divides \((d,N)\).

Fix \(0<\epsilon<1/100\), \(B\ge0\), and set
\[
 z=N^{1/2-\epsilon},\qquad D=N^{1/2}(\log N)^{-B}.
 \tag{2}
\]
All statements at these critical scales are for \(N\) sufficiently large
after \(\epsilon,B\) are fixed; in particular \(1<\sqrt D<z<D\).
Let \(\mathcal P_7\) be the primes \(7\bmod12\), and
\(P_7(z)=\prod_{p<z,\ p\in\mathcal P_7}p\).
The lower half-dimensional weight is
\[
 \lambda_d^-(D;z)=
 \mu(d)\mathbf1_{d\mid P_7(z)}
 \mathbf1_{d\in\mathcal D_D^-},
 \tag{3}
\]
where, for \(d=p_1\cdots p_r\) with \(p_1>\cdots>p_r\),
\[
 \mathcal D_D^-=
 \{d:p_1\cdots p_{2j-1}p_{2j}^2<D
                  \text{ for every }2j\le r\}.
 \tag{4}
\]
At the present scales \(z<D\), this weight is supported on \(d<D\).
For \(r=1\), the condition in (4) is empty; for odd \(r\ge3\),
the last even-index condition implies \(d<D\) since \(p_r<p_{r-1}\).
For even \(r\) the same conclusion is immediate.

A signed estimate
\[
 \sum_d\lambda_d^-(D;z)r_d(N)
       \ge-o\!\left(X_N\mathcal V_N(z)\right),
 \qquad
 \mathcal V_N(z)=\prod_{p<z}(1-g_N(p)),
 \tag{5}
\]
with the appropriate uniformity for one sufficiently small fixed
\(\epsilon\), would complete the critical switching criterion.
The lower main is of size
\(\sqrt\epsilon\,X_N\mathcal V_N(z)\); the already proved bad-survivor
upper bound is \(O(\epsilon^{3/2}X_N\mathcal V_N(z))\).
See [AsymCriticalCriterion](ASYMMETRIC_CRITICAL_SWITCHING.md).
Equation (5) remains open.

## 2. Balanced factorization excludes every large prime factor

A sequence \(\lambda\) is well factorable of level \(Q\) if, for
every real factorization \(Q=Q_1Q_2\), \(Q_i\ge1\), one can write
\(\lambda=\lambda_1*\lambda_2\) with
\(\operatorname{supp}\lambda_i\subset[1,Q_i]\), and with the stipulated
fixed divisor bounds on the factors. The following obstruction uses
only the support requirement.

**Proposition WFSupport.** If \(\lambda\) is well factorable of
level \(Q\), then
\[
 \lambda(n)=0\quad\text{whenever }P^+(n)>\sqrt Q.
 \tag{6}
\]
In particular \(\lambda(p)=0\) for every prime \(p>\sqrt Q\).
The same conclusion holds for a finite linear combination of
well-factorable sequences whose levels are at most \(Q\), regardless
of their orders or scalar coefficients.

**Proof.**
Choose \(Q_1=Q_2=\sqrt Q\). If \(p\mid n\) and \(p>\sqrt Q\), every
factorization \(n=ab\) has \(p\mid a\) or \(p\mid b\). The corresponding
factor exceeds its support bound. Thus every term in
\((\lambda_1*\lambda_2)(n)\) is zero. Apply this separately to each
member of a finite linear combination. \(\square\)

In contrast, (3) gives
\[
 \lambda_p^-(D;z)=-1
       \qquad(\sqrt D<p<z,\ p\in\mathcal P_7).
 \tag{7}
\]
Thus the actual half-dimensional lower weight is not well factorable
of level \(D\), nor a finite linear combination of such weights. This
is a support contradiction, independent of estimates for exponential
sums. Padding the well-factorable level to \(Q\) can remove this
particular obstruction only when \(Q\) reaches at least \(z^2\),
up to the endpoints of the available primes.

## 3. A pure well-factorable lower polynomial has no positive main

For a finite prime set \(\mathcal P\), let \(\lambda\) be real and supported on
squarefree divisors of \(\prod_{p\in\mathcal P}p\). A pointwise divisor
lower weight satisfies
\[
 L_\lambda(n):=\sum_{d\mid n}\lambda_d
        \le \mathbf1_{(n,\prod_{p\in\mathcal P}p)=1}
                  \qquad(n\ge1).
 \tag{8}
\]

**Proposition WFLowerObstruction.**
Suppose (8) holds, \(\mathcal P\) contains a prime \(p>\sqrt Q\), and
\(\lambda\) is a finite linear combination of well-factorable
sequences of levels at most \(Q\). Then
\[
 \lambda_1\le0,\qquad L_\lambda(n)\le0\quad(n\ge1).
 \tag{9}
\]
For every collection \(0\le g(p)\le1\), with multiplicative extension
to squarefree integers,
\[
 \sum_d\lambda_d g(d)\le0.
 \tag{10}
\]

**Proof.**
By (6), \(\lambda_p=0\). Evaluating (8) at \(n=p\) gives
\(\lambda_1\le0\). On a sifted integer \(L_\lambda(n)=\lambda_1\);
on every other integer (8) gives \(L_\lambda(n)\le0\). This proves
(9). For (10), average the same inequality over the finite probability
space in which each prime divides independently with probability
\(g(p)\). The expectation of \(\mathbf1_{d\mid n}\) is \(g(d)\).
\(\square\)

Consequently such a weight cannot produce a positive sieve lower
bound when \(z>\sqrt Q\). The argument is independent of the sieve
dimension. It applies to pure divisor lower polynomials; it does not
exclude a switched inequality containing other sequence-dependent
terms.

### 3.1 The obstruction persists on actual even-parity pairs

The abstract test \(n=p\) is not an actual second summand for this
problem, because every \(b\) with \(w_N(b)>0\) has even
\(\Omega_7(b)\). There is nevertheless an actual-pair obstruction.

**Proposition WFActualPairObstruction.**
Fix \(0<\epsilon<1/100\), and \(0<\theta<1-2\epsilon\). There are
infinitely many \(N\equiv2\pmod{12}\) with a positive-weight bad pair
\(1+b=N\) such that
\[
 b=cp^2,\quad c\in C,\quad p\in\mathcal P_7,\qquad
 N^{\theta/2}<p<N^{1/2-\epsilon}.
 \tag{11}
\]
For each such \(N\), any pure divisor lower weight on the
\(\mathcal P_7\)-divisors, which is a finite linear combination of
well-factorable weights of levels at most \(Q=N^\theta\), and which
is a valid lower bound on the actual positive-weight support, has
\(\lambda_1\le0\). Its lower polynomial is then nonpositive on
every actual positive-weight pair.

**Proof.**
Choose
\[
 \frac{4\epsilon}{1-2\epsilon}<\alpha<\frac2\theta-2,
 \qquad
 c=13^{\lceil\alpha\log p/\log13\rceil},\qquad N=cp^2+1,
 \tag{12}
\]
along an unbounded sequence of primes \(p\equiv7\pmod{12}\).
The interval for \(\alpha\) is nonempty exactly because
\(\theta<1-2\epsilon\). Since \(p^\alpha\le c<13p^\alpha\),
\[
 \frac{\log p}{\log N}\longrightarrow\frac1{2+\alpha}
       \in(\theta/2,\,1/2-\epsilon).
\]
This proves (11). Also \(c\equiv p^2\equiv1\pmod{12}\), so
\(N\equiv2\pmod{12}\), \(1\in C\), \(b\in A\), and
\(w_N(b)=V(cp^2)>0\). Here \(\Omega_7(b)=2\), so the exact parity
constraint is respected.

The only squarefree \(\mathcal P_7\)-divisors of \(b\) are \(1,p\).
By (6), its lower polynomial is therefore \(\lambda_1\).
Since \(p<z\), this is a bad, unsifted pair, forcing
\(\lambda_1\le0\). Every sifted actual pair also has polynomial
\(\lambda_1\), and the defining lower inequality bounds the
polynomial by zero on the other actual pairs. \(\square\)

For the particular critical level \(Q=D\) in (2), one can use the
fixed choice \(\alpha=1/10\) in (12):
\[
 c=13^{\lceil\log p/(10\log13)\rceil},\qquad N=cp^2+1.
 \tag{13}
\]
Then \(\sqrt D<p<z\) eventually, because
\(1/(2+1/10)=10/21<1/2-\epsilon\) for \(\epsilon<1/42\), while
\(\sqrt D\le N^{1/4}\). The permitted \(\epsilon<1/100\) suffices.
The same pair is a valid \(C_4+C_3\) pair, \(1+cp^2\), so the
support obstruction is not removed by starting from that larger
unweighted sequence.

### 3.2 Failure of a negligible local absolute approximation

Let \(Q=N^\theta\), where \(0<\theta<1-2\epsilon\), and let \(h_d\)
be any finite linear combination of well-factorable sequences of
levels at most \(Q\). Then, at (2),
\[
 \begin{split}
 \sum_{d<D}|\lambda_d^-(D;z)-h_d|g_N(d)
 &\ge \sum_{\substack{\sqrt Q<p<z\\p\equiv7\ (12)\\p\nmid N}}g_N(p)\\
 &=\frac12\log\frac{1-2\epsilon}{\theta}+o(1)>0.
 \end{split}
 \tag{14}
\]
Indeed \(p<z<D\), (6)--(7) apply, and
\(g_N(p)=2/p+O(p^{-2})\) off the target. The prime number theorem
in the fixed class \(7\bmod12\), followed by partial summation,
gives the displayed limit before deleting target primes.
There are at most \(2/\theta+O(1)\) prime factors of \(N\) above
\(\sqrt Q\), so deleting them changes the sum by \(O_\theta(N^{-\theta/2})\).
For \(Q=D\), the same argument gives
\[
 \sum_{d<D}|\lambda_d^--h_d|g_N(d)
       \ge\frac12\log(2-4\epsilon)+o(1).
 \tag{15}
\]
This rules out a negligible local absolute replacement by such
weights. It does not rule out cancellation between omitted prime
and composite terms, or an argument using the actual sequence.

## 4. Exact switching of the large-prime lower branch

For squarefree \(s=q_1\cdots q_m\), \(q_1>\cdots>q_m\), supported on
\(\mathcal P_7\), define the upper half-dimensional support and weight by
\[
 \begin{split}
 \mathcal D_y^+
 &=\{s:q_1\cdots q_{2j-2}q_{2j-1}^2<y
                           \text{ for every }2j-1\le m\},\\
 \lambda_s^+(y)&=\mu(s)\mathbf1_{s\in\mathcal D_y^+},
 \qquad \lambda_1^+(y)=1.
 \end{split}
 \tag{16}
\]
The primary support definitions are Iwaniec, *The half dimensional
sieve*, Acta Arith. 29 (1976), printed p.75, (3.5)--(3.6),
[cached source](sources/iwaniec-half-dimensional-sieve-1976.pdf).
The finite pairing cutoff can be chosen above every possible
factor count for the fixed support bound.

**Proposition RosserLargePrimeSwitch.**
For \(p\in\mathcal P_7\), \(\sqrt D<p<z\), one has exactly
\[
 \lambda_{ps}^-(D;z)=-\lambda_s^+(D/p)
 \tag{17}
\]
on the possible squarefree cofactors \(s\); the nonzero upper side
automatically has
\[
 P^+(s)<\sqrt{D/p}<D^{1/4}<p,\qquad s<D/p.
 \tag{18}
\]
Consequently the entire large-prime branch of the signed remainder is
\[
 \boxed{\displaystyle
 \sum_{\substack{d<D\\P^+(d)>\sqrt D}}
       \lambda_d^-(D;z)r_d(N)
   =-\sum_{\substack{\sqrt D<p<z\\p\equiv7\ (12)}}
        \sum_s\lambda_s^+(D/p)\,r_{ps}(N).}
 \tag{19}
\]

**Proof.**
The first even-index lower condition for \(ps\) is \(pq_1^2<D\).
The subsequent conditions are
\[
 p q_1\cdots q_{2j-2}q_{2j-1}^2<D,
\]
exactly the upper conditions in (16) with \(y=D/p\).
The first one gives (18); the product bound for \(s\) follows from
the last odd-index condition and the decreasing prime ordering.
Thus all cofactor primes are below \(z\), \(p\nmid s\), and
\(\mu(ps)=-\mu(s)\), proving (17). An integer \(d<D\) has at most
one prime factor above \(\sqrt D\), so this parametrization is
unique. This proves (19). Terms with \(p\mid N\) vanish by (1).
\(\square\)

The exact switched quantity on the right of (19) is a genuine
remaining obligation. No claim that \(\lambda^+(D/p)\) is well
factorable follows from its definition or from (18). Any use of a
well-factorable theorem for these coefficients needs its own
factorization statement or a different switching argument.

## 5. Exact scope of Fouvry's primary theorem

The relevant primary source is E. Fouvry, *Autour du théorème de
Bombieri--Vinogradov. II*, Ann. Sci. École Norm. Sup. (4) 20 (1987),
617--640:
[Numdam primary PDF](https://www.numdam.org/item/ASENS_1987_4_20_4_617_0.pdf).

Let \(\alpha_m,\beta_n\) have fixed divisor bounds and supports
\([M,2M]\), \([Y,2Y]\), and set \(x=4MY\). The source uses \(N\)
for the present \(Y\). Its discrepancy is
\[
 \Delta(f;x;q,a)
 =\sum_{\substack{n\le x\\n\equiv a\ (q)}}f(n)
  -\frac1{\varphi(q)}
       \sum_{\substack{n\le x\\(n,q)=1}}f(n),
 \quad
 E(c,f;x,Q;a)=\sum_{\substack{q\le Q\\(q,a)=1}}c_q\Delta(f;x;q,a).
 \tag{20}
\]
The absolute value is outside the signed weighted modulus sum.

The essential hypothesis on \(\beta\) is the source's (1.3):
for every \(B_1>0\), uniformly in \(d,k\ge1\) and \((a_0,k)=1\),
\[
 \sum_{\substack{n\equiv a_0\ (k)\\(n,d)=1}}\beta_n
 =\frac1{\varphi(k)}\sum_{(n,dk)=1}\beta_n
   +O_{B_1,K}\!\left(
       \frac{Y\tau_K(d)}{(\log(2Y))^{B_1}}\right).
 \tag{21}
\]
For \(Y\ge x^\delta\), the printed p.619 main theorem gives
\(E(\gamma*\xi,\alpha*\beta;x,4RS;a)\ll x(\log x)^{-A_0}\),
with \(\gamma,\xi\) supported at scales \(R,S\), under one of five
sets of inequalities. Only its alternatives (C.1) and (C.2) are
stated uniformly for the large residue \(|a|\le x\):
\[
 \begin{array}{ll}
 {\rm(C.1)}&
 \max(Y^{-1}R^2S,\ Y^{3/2}R^{7/4}S^2)\le x^{1-\delta},\\[2mm]
 {\rm(C.2)}&
 \max(R^2S,\ Y^{5/4}R^{7/4}S^2)\le x^{1-\delta}.
 \end{array}
 \tag{22}
\]
Alternatives (C.3)--(C.5) are stated only for
\(|a|\le(\log x)^{A_0}\), so they cannot be applied to the target
residue of size comparable to the convolution scale.

Printed p.620 Corollary 2(i) gives a well-factorable weight of
order 1 (the convention defined on p.618), with level
\[
 Q=x^{5(1-\nu)/9-\delta},\qquad
 \delta\le\nu=\frac{\log Y}{\log x}\le\frac1{10},
 \qquad |a|\le x.
 \tag{23}
\]
For a literal choice with both factors at least 1, take
\[
 s=\max\!\left(0,\frac{1-10\nu}{9}-\frac\delta2\right),\qquad
 r=\frac{5(1-\nu)}9-\delta-s,\qquad R=x^r,\ S=x^s.
 \tag{24}
\]
When the maximum is attained by its second entry, this is the
printed choice \(r=(4+5\nu)/9-\delta/2\),
\(s=(1-10\nu)/9-\delta/2\), which satisfies (C.2).
In the clamped range \(s=0\), one has
\(\nu\ge1/10-9\delta/20\), and directly
\[
 2r\le1-\tfrac32\delta,\qquad
 \tfrac74r+\tfrac54\nu\le1-\tfrac74\delta.
\]
Thus (C.2) still holds, including \(\nu=1/10\), and
\(r+s=5(1-\nu)/9-\delta\) is unchanged.

### 5.1 A range barrier in the full large-residue theorem

Write \(R=x^r\), \(S=x^s\), \(Y=x^\nu\), and discard the favorable
positive losses in (22). In case (C.1),
\[
 2r+s-\nu\le1,\qquad
 \tfrac74r+2s+\tfrac32\nu\le1.
\]
Multiply the first inequality by \(1/9\) and the second by \(4/9\):
\[
 r+s\le\frac{5(1-\nu)}9.
 \tag{25}
\]
In case (C.2), the same combination of
\[
 2r+s\le1,\qquad
 \tfrac74r+2s+\tfrac54\nu\le1
\]
gives exactly (25). Hence even the full two-factor modulus form of
the primary theorem cannot exceed the \(5/9\) combined-modulus
exponent in its large-residue range. This is a limitation of these
stated hypotheses, not a general impossibility theorem.

### 5.2 A real beyond-half range for arbitrary modulus weights

The same printed p.619 Corollary 1, taking \(S=1\) in (C.1),
requires no factorability of the modulus weight. For any
divisor-bounded \(c_q\), it gives
\[
 E(c,\alpha*\beta;x,Q;a)\ll x(\log x)^{-A_0},
 \quad |a|\le x,\qquad
 Qx^\delta\le
 \min(x^{1/2}Y^{1/2},\ x^{4/7}Y^{-6/7}).
 \tag{26}
\]
The exponent range is strictly beyond one half for
\(0<\nu<1/12\), away from endpoints by a fixed amount.
Its maximum is \(10/19\), at \(\nu=1/19\).
This gives an actual usable bilinear-box estimate whenever (21)
and the exact residue/scale hypotheses hold. It does not cover all
boxes of the asymmetric critical remainder.

A literal principal factor
\(\beta_n=\mathbf1_{[Y,2Y]}(n)\chi_0^{(12)}(n)\)
satisfies (21), after enlarging the fixed divisor order to at least 2.
Indeed inclusion-exclusion for \((n,12d)=1\) gives identical
progression and averaged main terms, with error \(O(\tau(12d))\);
in the averaged error one uses
\(\tau(k)/\varphi(k)\ll1\). This is absorbed by (21) since
\(Y(\log(2Y))^{-B_1}\) tends to infinity.
By contrast, a nonprincipal literal character
\(\beta_n=\psi(n)\mathbf1_{[Y,2Y]}(n)\) modulo 12 does not satisfy
(21): take \(d=1,k=12\) and a unit \(a_0\). Its progression sum is
\(\psi(a_0)Y/12+O(1)\), whereas the averaged sum is \(O(1)\).
A fixed-period extension subtracting all resonant fixed characters
would have to be proved; it is not part of (21).

## 6. Relation to Wu's formulation

J. Wu, *Chen's double sieve, Goldbach's conjecture and the twin prime
problem*, Acta Arith. 114 (2004), 215--273,
[official primary PDF](https://www.impan.pl/shop/publication/transaction/download/product/82961),
printed p.222 Lemma 2.6, explicitly cites Fouvry's Corollary 2(i).
It gives (23) for \(|a|\le X\).
Wu's common hypotheses on p.221 include both a condition of the
form (21) and a roughness condition excluding prime factors below
\(\exp((\log\log n)^2)\).
Fouvry's original theorem assumes (21) and the divisor bounds,
without that additional roughness hypothesis. Thus one should use
the original statement for a literal principal interval factor.

Wu's p.222 Lemma 2.5 has more favorable levels in several factor
ranges, but is stated only for \(|a|\le(\log X)^A\).
His p.223 Lemma 2.8 similarly has constants depending on the fixed
residue \(a\). Neither statement supplies the required uniformity
when the target residue grows like \(X\).
Wu's well-factorable linear-sieve formula, p.220 Lemma 2.2, assumes
\(z\le Q^{1/2}\), consistent with Proposition WFLowerObstruction.
It is not a factorization theorem for the present half-dimensional
weights with \(z\) almost \(D\).

The full 1993 Iwaniec--Pomykała paper, *Sums and differences of
quartic norms*, Mathematika 40 (1993), 233--245, has not been audited
here. Wu's reproduction of its auxiliary Proposition in Lemma 2.5
does not constitute an audit of that paper's main theorem. The
scope and priority comparison with the full 1993 main theorem
remain unresolved. No novelty or priority conclusion is drawn
from the auxiliary estimate alone.

## 7. Why the critical application still lies beyond these inputs

The exact compression already proved in
[AsymAuxiliaryAverage](ASYMMETRIC_AUXILIARY_AVERAGE.md)
expresses a lower-weight remainder, after a common-square cutoff
\(k\le K\), as
\[
 \sum_{\substack{k\le K\\k\ {\rm squarefree},\ p\mid k\Rightarrow p\in\mathcal P_7\\(k,N)=1}}
 \mu(k)
 \sum_{\substack{\ell<D/k\\\ell\ {\rm squarefree},\ p\mid\ell\Rightarrow p\in\mathcal P_7\\(\ell,kN)=1}}
 2^{\omega(\ell)}\lambda^-_{k\ell}(D;z)
 \left[C_{k^2\ell}(N)
          -\frac{X_NH(k\ell)}{k^2\ell}\right]
   +O_\delta(N^{1+\delta}/K).
 \tag{27}
\]
Here \(C_c(N)=\sum_{a+ct=N}U(a)V(t)\) and
\(H(c)=\prod_{p\mid c}(1-p^{-1}+p^{-2})^{-1}\).
Already the \(k=1\) slice allows \(c=\ell\) almost \(D\).
Expanding the second divisor convolution gives a divisor variable
\(q\) as large as \((N/c)^{1/2}\), and combined progression modulus
\[
 cq\asymp (Nc)^{1/2},
 \tag{28}
\]
which can be almost \(N^{3/4}\). The first product scale can be
comparable to \(N\); hence even the \(5/9\) ceiling (25) is
insufficient. A very small principal first factor can place some
boxes in (26), but no such factor of the required length is
guaranteed across the whole decomposition.

The \(q\)-sum phases allowed by BGLM Proposition 5.1 are not
arbitrary modulus weights, and its outer auxiliary modulus has a
small-power range. Therefore neither an unproved factorability
assertion for (3), nor a relabeling of \(cq\) as the old varying
modulus, repairs the missing critical estimate.
The exact additive resonance described in
[AsymSignTransfer](ASYMMETRIC_SIGN_TRANSFER.md)
also prevents replacing this remainder by a zero-main periodic
twist.

One concrete remaining route is to estimate the switched
prime-cofactor expression (19), retaining its nonzero local main,
together with the complementary branch having
\(P^+(d)\le\sqrt D\). Such an estimate would use information beyond
pure well-factorable divisor lower weights. No bound of size (5)
for those two pieces is established here.

## 8. Primary-source pins and verification scope

The primary statements above were read in the actual PDFs;
printed pp.619--620 of Fouvry were also visually inspected.
The retained primary sources are:

| File | SHA-256 | Exact locators |
|---|---|---|
| [Fouvry 1987](sources/fouvry-bombieri-vinogradov-II-1987.pdf) | `13dec04a3f215c0b5809dcefe77cd11fcad8a8982ddb9c471f1b6b08fad61c7f` | printed pp.618--620; PDF pp.3--5 including the Numdam cover |
| [Wu 2004](sources/wu-chens-double-sieve-2004.pdf) | `ebe75ef36478213bde937dce7d163bdd14b5b0352be54e085c9f696afc6e8a30` | printed pp.219--223; PDF pp.5--9 |

The support propositions, the actual-pair family, the local
absolute-approximation lower bound, the large-prime switch, and the
linear combination proving (25) are elementary arguments supplied
here. Their scope is a rigorous obstruction and an exact
reformulation. They do not establish (5), a full critical
distribution level, or an every-target two-sided representation.

## Dependencies and finite controls

Propositions WFSupport, WFLowerObstruction and
WFActualPairObstruction use only the stated factorization support,
the defining lower inequality, and the elementary target family
(12). The independent-model conclusion is finite probability.
The local absolute-approximation asymptotic uses the prime number
theorem in the fixed residue class \(7\bmod12\).
Proposition RosserLargePrimeSwitch uses precisely Iwaniec's support
definitions (4) and (16). The source-range ceiling (25) is the
displayed linear combination of Fouvry's two large-residue
hypotheses. No dispersion theorem is needed for these elementary
barriers or the exact switch.

The critical significance of (5) uses
[AsymCriticalCriterion](ASYMMETRIC_CRITICAL_SWITCHING.md);
the exact coefficient compression (27) and its currently available
small range use
[AsymAuxiliaryAverage](ASYMMETRIC_AUXILIARY_AVERAGE.md).
The resonance limitation uses
[AsymSignTransfer](ASYMMETRIC_SIGN_TRANSFER.md).
The exact density and original sieve setting are retained in
[AsymSemilinear](ASYMMETRIC_SEMILINEAR_SIEVE.md).

The [standard-library checker](eisenstein-restricted-prime-sums/verify_asymmetric_well_factorable_switch.py)
and its [compact receipt](eisenstein-restricted-prime-sums/asymmetric_well_factorable_switch_checks.json)
retain 1,426 finite large-prime switch cases, balanced-support and
actual-pair controls, and exact rational parameter checks with
deliberately invalid alternatives. The receipt pins this proof and
its canonical and primary dependencies. Run from any working directory:

```sh
python3 research/arithmetic/eisenstein-restricted-prime-sums/verify_asymmetric_well_factorable_switch.py
```

from the repository root, or use the script's absolute path elsewhere.
An explicit \(\texttt{--output PATH}\) writes a new receipt.

These finite controls do not prove the infinite support
propositions, the prime-distribution asymptotics, the imported
dispersion theorem, or the critical remainder. The proofs supplied
above establish the first of those claims; the precise primary
statements supply the cited inputs. No finite binary coverage scan,
external specialist review, formal verification, or priority
claim is used.
