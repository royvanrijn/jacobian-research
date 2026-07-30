# Keller/Tschirnhaus bridge card in ranks five through seven

This is an exact regression and comparison card, not a second canonical
proof.  The all-rank statement is the
[generic Tschirnhaus non-descent theorem](GENERIC_TSCHIRNHAUS_NON_DESCENT.md).
Its inputs live in:

- the
  [quadratic-gauge stable-moduli theorem](QUADRATIC_GAUGE_STABLE_MODULI.md#41-explicit-quotient-coordinates-on-the-compiler-slice),
  which proves that the compiler seed quotient has the saturated coordinates
  \(I_5,J_6,\ldots,J_N\);
- the
  [all-rank collision-projective theorem](ALL_RANK_COLLISION_PROJECTIVE_DESCENT.md),
  which proves that framed projective transport has \(N-3\) independent
  residuals; and
- the
  [universal relative Keller map](UNIVERSAL_RELATIVE_KELLER_MAP.md#41-quantitative-seed-descent-defect),
  which combines them into the quantitative Tschirnhaus-descent obstruction.

The purpose here is only to record the exact \(N=5,6,7\) specialization and
to compare four arithmetic/geometric complexity measurements without
conflating their status.

## 1. Exact low-rank comparison

For the \(S_N\)-extension problem over a characteristic-zero field \(k\),
let:

1. \(\operatorname{ed}_k(S_N)\) be essential dimension;
2. \(\operatorname{RD}_k(S_N)\) be resolvent degree;
3. \(\operatorname{ktdim}_k(S_N)\) be the proposed least target dimension
   of one fixed Keller self-map carrying a field-versal family of rank-\(N\)
   finite-etale fibres; and
4. \(\operatorname{kdeg}_{k,N}(m)\) be the proposed least maximum coordinate
   degree at fixed Keller target dimension \(m\).

Only (1)--(2) are established external invariants.  Items (3)--(4) are
organizing definitions for this research programme.  They satisfy the
formal bounds

\[
 \operatorname{RD}_k(S_N)
 \le \operatorname{ed}_k(S_N)
 \le \operatorname{ktdim}_k(S_N)
 \le N,                                                \tag{1.1}
\]

where the last inequality is supplied by the promoted universal map
\({\cal U}_N:\mathbb A^N\to\mathbb A^N\).

\[
\begin{array}{c|c|c|c|c|c}
N&\operatorname{ed}_k(S_N)&\operatorname{RD}_{\mathbb C}(S_N)
&\operatorname{ktdim}_k(S_N)&
\deg_{\rm vertical}&\deg {\cal U}_N\\ \hline
5&2&1&2\le\operatorname{ktdim}\le5&32&33\\
6&3&1\le\operatorname{RD}\le2&3\le\operatorname{ktdim}\le6&38&39\\
7&4&1\le\operatorname{RD}\le3&4\le\operatorname{ktdim}\le7&44&45
\end{array}                                             \tag{1.2}
\]

The vertical coordinate degree is \(6N+2\).  Promoting the top seed
coefficient \(u_N\) adds one to total degree, giving \(6N+3\).  These are
upper bounds from the present construction, not minimality results.

The essential dimensions in (1.2) are exact.  The equalities
\(\operatorname{RD}_{\mathbb C}(S_6)=2\) and
\(\operatorname{RD}_{\mathbb C}(S_7)=3\) remain respectively Hilbert's
sextic conjecture and the algebraic form of Hilbert's thirteenth problem.

References:

- [Edens--Reichstein, essential dimension of symmetric groups](https://arxiv.org/abs/2308.10096);
- [Duncan, essential dimensions of \(A_7\) and \(S_7\)](https://arxiv.org/abs/0908.3220);
- [Farb--Wolfson, resolvent degree and Hilbert's thirteenth problem](https://arxiv.org/abs/1803.04063);
- [Edens--Reichstein, current resolvent-degree status](https://arxiv.org/abs/2406.15954).

## 2. The exact common-algebra witness

Take the split algebra \(A_N=\mathbb Q^N\) with primitive coordinate

\[
 r=(1,2,\ldots,N)
\]

and the quadratic Tschirnhaus coordinate

\[
 u=r+r^2=(2,6,12,20,\ldots,N+N^2).                   \tag{2.1}
\]

Put

\[
 P_r(T)=\prod_{i=1}^N(T-i),\qquad
 P_u(U)=\prod_{i=1}^N(U-i-i^2).
\]

Exact interpolation produces a polynomial \(h_N(U)\), of degree less than
\(N\), such that

\[
 h_N(T+T^2)\equiv T\pmod{P_r(T)},\qquad
 h_N(U)+h_N(U)^2\equiv U\pmod{P_u(U)}.                \tag{2.2}
\]

Hence

\[
 \mathbb Q[U]/(P_u)\simeq\mathbb Q[T]/(P_r),           \tag{2.3}
\]

and the universal compiler realizes each presentation as a complete Keller
fibre with the same abstract finite-etale algebra.  This does not assert a
source--target equivalence carrying one marked fibre to the other.

The coordinate change is nevertheless nonprojective.  After the first
three root pairs determine the unique Möbius transformation, the residual
at \(i=4,\ldots,N\) is

\[
 R_i=-2(i-1)(i-2)(i-3).                               \tag{2.4}
\]

Thus all \(N-3\) framed projective residuals are nonzero.

## 3. Exact stable-boundary separation

On the compiler slice, one residual scaling removes one of the \(N-3\)
seed coordinates.  The canonical stable-boundary coordinates specialize to

\[
 I_5=\frac{u_5^5}{u_4^6},\qquad
 J_6=\frac{u_4u_6}{u_5^2},\qquad
 J_7=\frac{u_5u_7}{u_6^2}.                            \tag{3.1}
\]

The exact fingerprints of the two presentations are:

\[
\begin{array}{c|c|c}
N&\Phi(P_r)&\Phi(P_u)\\ \hline
5&
\displaystyle\frac{75076}{968203125}&
\displaystyle\frac{1296}{50236123}\\[4pt]
6&
\displaystyle\left(\frac{734832}{1220703125},\frac{25}{63}\right)&
\displaystyle\left(\frac{13436928000}{70946061021073},\frac{83}{224}\right)\\[6pt]
7&
\displaystyle\left(
\frac{68696948942127}{44596500536000000},\frac{280}{529},\frac{23}{56}
\right)&
\displaystyle\left(
\frac{254014812517232151861328125}
     {548382681573171461681967824896},
\frac{43696}{88725},\frac{65}{168}
\right).
\end{array}                                             \tag{3.2}
\]

Every applicable coordinate changes.  Equations (2.3) and (3.2) isolate
the obstruction:

\[
\boxed{
\text{the finite-etale fibre descends, but the intrinsic
reconstruction-boundary Fitting divisor does not.}}    \tag{3.3}
\]

For \(N=5,6,7\), the respective counts are

\[
\begin{array}{c|ccc}
N&5&6&7\\ \hline
\text{compiler seed parameters}&2&3&4\\
\text{projective residuals}&2&3&4\\
\text{residual scaling redundancies}&1&1&1\\
\text{stable boundary coordinates}&1&2&3.
\end{array}                                             \tag{3.4}
\]

This is an obstruction to descent of the quadratic-gauge atlas.  It proves
no lower bound for a different Keller construction and no new
resolvent-degree equality.

For every `N>=6`, the same split witness is separated uniformly by the top
stable coordinate:

\[
 J_N(P_u)-J_N(P_r)
 =
 -\frac{(N-1)(7N+11)}
        {30N(N+1)(N+2)}\ne0.
\]

That symbolic all-rank calculation belongs to the canonical generic
non-descent theorem.  The table above remains a low-rank regression because
it checks every applicable fingerprint coordinate and the promoted
coordinate degrees.

## 4. Exact regression

Run

```bash
.venv/bin/python scripts/verify_keller_tschirnhaus_descent_567.py
```

Refresh the pinned certificate only intentionally:

```bash
.venv/bin/python scripts/verify_keller_tschirnhaus_descent_567.py --write
```

The generated witness is
[`artifacts/generated-results/keller_tschirnhaus_descent_567.json`](../artifacts/generated-results/keller_tschirnhaus_descent_567.json).
The checker verifies the quotient-algebra isomorphisms, compiler targets,
projective residuals, low-rank stable fingerprints, and coordinate degrees.
Stable invariance itself is imported from the canonical stable-moduli
theorem linked above.
