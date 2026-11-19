# A weighted marked comparison for Hecke-square self-correlations

The bounded invariant-factor marked comparison below is proved under the
existing scoped G and ideal-comparison inputs. Its prime weights are bounded
and nonnegative, may vanish, and have fixed positive quadratic bank masses.
The fixed-shift Hecke-square asymptotic remains **OPEN**. In particular this
comparison does not prove the required automorphic prime estimate or the
unbounded-weight arithmetic reduction.

## 1. Source statements actually inspected

The distinction from [the existing bounded-indicator comparison](SPARSE_ARITHMETIC_MODEL_COMPARISON.md)
is visible in its (C1)--(C4) assumptions and the marked
comparison calculation (5.15)--(5.19) in
[the graph proof](../../results/pointwise-binary-splitting-semigroups/graph.tex).
The imported [weighted-dilation preprint](sources/weighted-dilation-graphs-sep24-2026.pdf), Definitions3.2--3.3, Lemma3.4,
Theorem3.5 and Corollary3.11 (printed pp.11--12 and26).

G uses label laws `mu_g(p)=1/(V_g*p)`, physical state mass `V_g^(-M)`
per ordered distinct list, and coefficients K_nu that ignore small-bank
labels. Corollary3.11 permits a mark-independent first endpoint with
supremum `exp(O(sqrt(log X)))` and Hilbert-space norms
`X^(1/2)*(log X)^C3`. Only the exponent C3, not the supremum threshold,
enters the required fixed saving. The ideal comparison coefficients and
their total logarithmic cost must be fixed independently of mark length M.

Walsh v2 §1.2.4, printed pp.5--6, identifies the normalized holomorphic
Hecke-square correlation and sums-of-two-squares indicator extension.
His earlier §1.2.1 explains the arithmetic difference: on one-prime
bank strata the divisor factor is the *constant* `k*l`. For Hecke
weights it is `a(p)*a(q)` and varies with the prime labels. These are
future extensions in that source, not results claimed by its Theorem1.1.
Primary source: https://arxiv.org/abs/2610.40439v2.

## 2. Exact whole-weight obstruction, not an impossibility of factoring

For Delta put `a(n)=tau(n)^2/n^11`. Exact invariance of the whole
weight fails for every prime: `a(p)=1` would make the integer square
`tau(p)^2` equal the nonsquare `p^11`. Thus (C2) already fails at n=1.
Also `a(11)=285809990544/285311670611>1`, violating the theorem's
literal `[0,1]` bound. The coefficient tau(11)=534612 can be verified
directly from `q*prod_(m>=1)(1-q^m)^24`, so no numerical approximation
enters these observations.

The Hecke recurrence is `lambda(p^2)=lambda(p)^2-1`, hence
`a(p^2)=(a(p)-1)^2`, rather than `a(p)^2`. Repeated prime factors
therefore require their actual recurrence or a proved negligible tail.

This does not stop Walsh's factorization. For any active prime union P,
`a(n)=a(n_P)*a(n_(P^c))`, and the second factor is invariant under every
prime in P, including repeated multiplication. A record-preserving model
gives the same factorization for its centered remainder. The issue is
the varying record factor, and not the lack of a formal invariant rough
factor.

## 3. Weighted marked-comparison lemma for a self-correlation

Fix the original disjoint banks, with at least `c_0*log log X` of
each small and big type, and total `s<=C_0*log log X`. The small-bank
primes lie between `exp((log X)^(1/8))` and `exp((log X)^(1/4))`,
and the big-bank primes between `exp((log X)^(1/3))` and
`exp((log X)^(2/5))`. Use the same fixed reduced prime class as in
the scoped ideal comparison; for Delta one may take all odd primes.
They have their original
reciprocal masses `v_-<=V_g<=v_+` and minimum prime
`P_*=exp((log X)^c)`, c>0. Let b be any fixed nonnegative prime weight
with `0<=b(p)<=A`, for fixed A. Assume

    w_- <= W_g := sum_(p in P_g) b(p)^2/p <= w_+,
    R_g := W_g/V_g,

with fixed positive bounds independent of X, g and mark length. Let
alpha,beta be bounded by one and invariant under multiplication by all
bank primes. The usual fixed-shift sign/cutoff conditions are retained.

Delete all bank primes dividing a fixed nonzero shift h. Put
`L=log X`, let `omega_g(n)` count distinct bank-prime divisors of
`|n|`, and take a fixed smooth function Psi compactly supported in
`(0,infinity)`. Extend alpha and beta by zero on nonpositive integers
if only positive endpoints are wanted; this preserves their prime
invariance. They may depend on X, but their entire definition stays
tied to that outer scale under every dilation.

For `1<=t_g<=m+1` define the weighted list weight

    W_t^(b,u)(n)
      = prod_g [1_(omega_g(n)=t_g) V_g^(-t_g)]
        * sum_(ordered distinct t-lists dividing n)
                  u(list) prod_(all endpoint labels p) b(p).

Here u is a usual bounded list coefficient; it ignores small-bank
labels, as in the original comparison. The *additional displayed b*
factors may depend on those small labels. This distinction is important.

Define the harmonic raw correlation by

\[
 C_0^{b}=\sum_{n>0}\frac{\Psi(n/X)}n
          \alpha(n)\beta(n+h)
          \mathcal W_{\mathbf1}^{b,1}(n)
          \mathcal W_{\mathbf1}^{b,1}(n+h).
\]

**Proposition HeckeWeightedMarkedComparison.** Under the stated bank, weight,
invariant-factor and cutoff hypotheses, and the scoped G and ideal-comparison
inputs, for every fixed accuracy A_0 there is a comparison family such that

\[
 C_0^b+\sum_\nu c_\nu
       \left(\prod_gR_g^{1-t_{\nu,g}}\right)C_\nu^b
                                             =O_{A_0}(L^{-A_0}),
\]

where

\[
\begin{split}
 C_\nu^b={}&\mathbb E_{D_1,D_2} A_\nu(D_1)B_\nu(D_2)
  \sum_{n>0}\frac{\Psi(n/(XD_1D_2))}{n}
       \alpha(n)\beta(n+hD_1D_2)\\
 &\hspace{12mm}\cdot
       \mathcal W_{\mathbf t_\nu}^{b,u_\nu}(n)
       \mathcal W_{\mathbf t_\nu}^{b,v_\nu}(n+hD_1D_2).
\end{split}
\]

The products D_1,D_2 are independent products of labels drawn with the
original laws `1/(V_g*p)`, use at most m slots per big bank and at
least one big-bank prime each. The functions A_nu, B_nu, u_nu, v_nu
have modulus at most one, and u_nu,v_nu ignore small-bank labels.
In each small bank `t_(nu,g)=1`. There are fixed m and C such that
`#{nu}+sum|c_nu|<=L^C`; both are fixed before choosing the eventual
mark length M. The additional displayed scalar is

                    prod_g R_g^(1-t_g).                    (H.1)

The number of comparisons and their total logarithmic coefficient cost
are independent of the eventual fixed mark length. As before, for every
requested fixed accuracy, the error is smaller than that inverse power
of log X. No inverse b(p) is used, and b(p)=0 is permitted.

### Proof: normalized vertices

Let `M=J+1`, the original physical number of marks per bank. Write
`B(n)=prod_(p in P,p|n)b(p)` for the radical record factor. Compress
the physical graph to exact-count vertices `omega_g(n)=M`, as in the
existing proof, and put

    F(n)=alpha(n) B(n) prod_g R_g^(-M/2),
    G(n)=beta(n) B(n) prod_g R_g^(-M/2).

These are functions of the integer only. They are independent of the
ordering or choice of its marks. On exact-count vertices, every bank's
ordered M-list contains all its M prime divisors. Their supremum is at
most `exp(O_M(s))=(log X)^(O_M(1))`, and therefore satisfies the source
`exp(O(sqrt(log X)))` condition for each fixed M.

Their physical L2 norms have no M-dependent logarithmic exponent.
For each full ordered list with product Q, count its multiples in a
containing interval of length O(X). Because M is fixed and the banks
are subpower, `Q=X^o(1)<X^(1/2)` eventually. The count is at most
`C X/Q`, with C independent of the list and of M (the threshold may
depend on M). Drop the exact-count restrictions for this upper bound.
For one bank the normalized list sum is at most

    V_g^(-M) R_g^(-M)
      * sum_(distinct ordered p_1,...,p_M) prod_j b(p_j)^2/p_j
      <= (V_g R_g)^(-M) W_g^M = 1.               (H.2)

Thus `||F||^2,||G||^2<=C X`. With bounded invariant factors at most
`(log X)^B`, the corresponding bound is `C X(log X)^(2B)`;
the exponent remains independent of M. Corollary3.11 now applies
without increasing ideal accuracy after choosing M. The original
kernel and physical operators are unchanged. In particular this proof
does not assert G for a new label law.

### Proof: shared-label normalization

For a comparison with reduced counts t_g, each bank has M-t_g shared
labels. Away from the same collision sets treated in the original
compression, writing `n=D_R*w`, `n'=D_R*w'` gives

    B(n)B(n') = B(D_R)^2 B(w)B(w').

The shared normalization formerly gave independent law
`mu_g(p)=1/(V_g*p)`. Inserting the shared factors b(p)^2 and averaging
therefore gives

    R_g^(M-t_g)
      * [independent draws with law b(p)^2/(W_g*p)].

Combined with the two vertex normalizers `R_g^(-M)`, this leaves
`R_g^(-t_g)`. The raw pattern has t_g=1 in every bank. Multiply the
entire relation by `prod_g R_g` to obtain (H.1).

Both R_g and R_g^(-1) have fixed bounds. Since `0<=t_g-1<=m` and
`s=O(log log X)`, the multiplier (H.1) costs only
`exp(O_m(s))=(log X)^(O_m(1))`, independent of M. Weighted endpoint
list factors also have size at most `exp(O_m(s))`. The chosen graph
accuracy can absorb these costs before M is finally selected.

### Proof: square/collision errors

The argument uses no lower bound on any individual b(p). The tilted
shared draw has maximal atom at most `A^2/(w_-*P_*)`. There are
`O_M(s)` shared draws. For fixed reduced endpoints in a polynomial
range and a fixed auxiliary product, the number of excluded primes
dividing `w*w'*D_C` is O(log X). Repetitions and intersections with
those divisors therefore have probability

    O_(M,A,w_-)(s*log X/P_*).

The unweighted auxiliary draws have the same exponential saving.
Products of permitted comparison costs and remaining list weights
contribute only fixed logarithmic powers, independently of M; even a
fixed M-dependent logarithmic cost is beaten by `P_*`. Hence the total
harmonic-weight error is `O_M((log X)^C/P_*)`, smaller than every fixed
inverse logarithmic power. The cases in which a shared prime also
divides a reduced endpoint are exactly among these excluded events.
Thus the radical factorization used above is legitimate at the stated
accuracy. This does not yet prove that higher bank powers have
negligible weight in the original *unbounded a(n)* correlation; that
belongs to the separate arithmetic tail/envelope reduction.

These calculations reproduce the original proof's shared-label,
auxiliary-label and forbidden-list cleanup with changed maximal atoms
and explicit scalar masses. The ideal coefficients still ignore small
labels, so all their original structural hypotheses remain true.

## 4. Specialization to Delta and what it does not supply

Use `b(p)=a(p)=|lambda_Delta(p)|^2`. The holomorphic Hecke bound gives
`0<=b(p)<=4`. Mangerel's Lemma5.6, printedpp34--35, gives
`sum_(p in P_g)b(p)/p=V_g+O(1/log P_*)` for the bank intervals.
No fourth-moment theorem is needed to obtain the quadratic mass bounds.
Indeed, Cauchy--Schwarz and `b^2<=4b` give

    (sum_(p in P_g)b(p)/p)^2 / V_g
        <= W_g <= 4 sum_(p in P_g)b(p)/p.

Consequently `R_g` lies between `1+o(1)` and `4+o(1)`, uniformly in
the actual banks. Fixed positive bounds, for example `1/2<=R_g<=5`
for sufficiently large X, suffice throughout Section3. The marked
comparison is therefore compatible with the Delta self-case without
using Sato--Tate or symmetric-power functoriality beyond this fixed
untwisted prime first moment.

There is an important two-form limitation. With different weights b_f
and b_g, shared labels contribute `sum b_f(p)b_g(p)/p`. Normalizing the
two physical L2 norms instead uses `sum b_f^2/p` and `sum b_g^2/p`.
Their ratio can be strictly below one. In this normalization route a
resulting M-dependent deficit cannot be dismissed as a fixed coefficient
cost. This is a limitation of the present proof, not an impossibility
theorem for other graphs or normalizations. The self-case avoids the
mismatch; the proposition makes no general f,g claim.

The later analytic factorization now extracts the explicit small-slot
weight b(p). In place of the existing unweighted prime sum it needs

    sum_(p in I) b(p) eta(p) p^(it)
      = sum_(p in I) eta(p)p^(it)
        + sum_(p in I) lambda_Delta(p^2) eta(p)p^(it).       (H.3)

The first term is the current Dirichlet problem. The second is the
prime trace for the twisted symmetric-square GL3 L-function. In
particular, ordinary fixed-modulus Sato--Tate or a prime number theorem
at bounded height does not imply the uniform conductor/height/cell
estimates required by the comparison-energy argument. Standard
zero-free regions narrowing like `1-c/log(q(|t|+3))` do not yield a
fixed-power saving on `P=exp((log X)^a)` at polynomial X-height.

There remain separate obligations: invariant-factor truncation with
negligible weighted tails; a positive bank-record envelope for all
actual/model choices; recordwise automorphic rough centering; the
weighted version of the major-arc/high-frequency estimate represented
by (H.3); and a rigorous model-product evaluation. The marked graph
calculation supplies none of those estimates by itself.

No new fixed-shift asymptotic has been proved here.

The inspected basic automorphic source is Mangerel,
[*Divisor-bounded multiplicative functions in short intervals*, v2](https://arxiv.org/abs/2108.11401v2),
Section1.3.1 (printedp4) for the holomorphic Hecke bound and Lemma5.6
(printedpp34--35) for prime first moments and harmonic interval masses.
The [cached primary edition](sources/mangerel-divisor-bounded-2108.11401v2.pdf)
states a fixed-representation untwisted result. It does not supply (H.3)
uniformly over the required Dirichlet twists and heights.
