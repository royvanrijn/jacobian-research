# Integer Gamma weights and an ordinary radial Strong Factorial family

Internal research consequence, 21 November 2026. Retained base inspected:
`19eb9912fb60142f8b620dac24896931d4ff3e46`.

## Result and scope

Let `q,p_1,...,p_r` be arbitrary **positive integers**; no primality assumption
is made. Define the normalized product-Gamma functional by

\[
 \mathcal L_{\mathbf p}(Y^\alpha)=\prod_i(p_i)_{\alpha_i},
 \qquad (p)_k=p(p+1)\cdots(p+k-1).
\]

For every characteristic-zero field and every nonzero phase

\[
 F(Y)=\left(\prod_{i=1}^{r}Y_i^{q p_i}\right)\sum_{i=1}^{r}c_iY_i,
\tag{1}
\]

every block of `h` consecutive positive moments under `L_p` contains a nonzero
moment, where `h` is the number of **distinct nonzero values** among the `c_i`.
In particular the window of length `s=#{i:c_i!=0}` works.

There is an exact ordinary-factorial consequence. For pairwise disjoint blocks
of variables put `S_i=X_(i,1)+...+X_(i,p_i)`. Then the ordinary polynomial

\[
 G(X)=\left(\prod_{i=1}^{r}S_i^{q p_i}\right)\sum_{i=1}^{r}c_iS_i
\tag{2}
\]

satisfies the same `h`-moment window, and hence its actual support-size SFC
window. This is a family with unbounded dimension, degree, and monomial count.

**The ordinary unmarked phase (1) with unequal exponents remains uncovered.**
The weights in `L_p` cannot be omitted. The construction (2) uses an exact
block-sum moment identity, not identification of variables or a general
moment-preserving change of coordinates.

## 1. Excluding the zero critical value in the contact bound

For a nonconstant polynomial `P`, put

\[
 \operatorname{CV}^{\times}(P)
   =\{P(z):P'(z)=0\}\setminus\{0\}.
\]

**Lemma 1.** If `P(0)=0`, `P'(0)=1`, and `P!=z`, then

\[
 \operatorname{ord}_0(P-z)-1\leq\#\operatorname{CV}^{\times}(P).
\tag{3}
\]

**Proof.** Use the retained perturbation argument in
[RIGIDITY_CRITICAL_VALUES.tex](RIGIDITY_CRITICAL_VALUES.tex), `lem:basins` and
`lem:contact`, lines 29–91. If `P=z+c z^(k+1)+O(z^(k+2))`, the polynomial
`P_delta=(1+delta^k)P` has `k` distinct nonzero attracting fixed points for
sufficiently small positive `delta`. The origin is a fixed point with multiplier
`1+delta^k>1`. Each attracting basin contains a critical value, and distinct
basins give distinct values. None of those values can be zero: the forward
orbit of zero is the fixed point zero, so it cannot converge to any of these
nonzero attracting points. Scalar multiplication preserves the number of
nonzero critical values. This proves (3). QED.

The underlying argument is credited to liqsweep, *Polynomial composition
rigidity via critical values*, September 2026,
[pinned source, Sections 2–5](https://github.com/blueberryvertigo/polynomial-composition-rigidity/blob/051c74aad9ef9bf8804c371964803aee099c5b63/rigidity.typ).
The deletion of zero is the additional observation used here. The retained
basin and perturbation arguments are complete unconditional written proofs;
no rigidity conjecture or SFC statement is assumed.

If `a(0)=b(0)=0`, the chain rule implies

\[
 \operatorname{CV}^{\times}(a\circ b)
 \subseteq\operatorname{CV}^{\times}(a)
       \cup\bigl(a(\operatorname{CV}^{\times}(b))\setminus\{0\}\bigr).
\tag{4}
\]

Indeed a critical value zero of `b` maps to `a(0)=0` and is excluded.
Consequently, for normalized `a,b` and nonidentity `a composed with b`,

\[
 \operatorname{ord}_0(a\circ b-z)
 \leq1+\#\operatorname{CV}^{\times}(a)
         +\#\operatorname{CV}^{\times}(b).
\tag{5}
\]

**Corollary 2.** Let `a=z+O(z^2)` be a nonidentity polynomial and put
`c=#CV^times(a)`. Its inverse has no `c` consecutive zero coefficients
starting at degree at least two.

**Proof.** First `c>=1` by Lemma 1 applied to `a`. If a block starts at degree
`n+1`, truncate `a^{-1}` at degree `n`. The resulting normalized polynomial
`b` gives contact at least `n+c+1`. But (5), together with
`#CV^times(b)<=deg b-1<=n-1`, bounds contact by `n+c`. If the composition is
the identity, multiplicativity of degree forces `a=z`, again a contradiction.
QED.

As in the retained proof, every characteristic-zero consequence descends using
an embedding into `C` of the finitely generated coefficient field. An embedding
of an arbitrary entire coefficient field is not required.

## 2. The weighted coefficient identity

For (1), set

\[
 a(z)=z\prod_{i=1}^{r}(1-c_i z^q)^{p_i},\qquad b(z)=a^{-1}(z).
\tag{6}
\]

For every integer `m>=0`, direct expansion gives

\[
 \mathcal L_{\mathbf p}(F^m)
 =m!\sum_{k_1+\cdots+k_r=m}
      \prod_i\frac{(p_i)_{q p_i m+k_i}}{k_i!}\,c_i^{k_i}.
\tag{7}
\]

Lagrange inversion gives

\[
\begin{aligned}
 [z^{qm+1}]b
 &=\frac1{qm+1}[u^{qm}]
       \prod_i(1-c_i u^q)^{-p_i(qm+1)}\\
 &=\frac1{qm+1}\sum_{|k|=m}
       \prod_i\frac{(p_i(qm+1))_{k_i}}{k_i!}\,c_i^{k_i}.
\end{aligned}
\tag{8}
\]

The rising factorials factor as

\[
 (p_i)_{q p_i m+k_i}
 =(p_i)_{q p_i m}\,(p_i(qm+1))_{k_i}.
\]

Consequently equations (7)–(8) give

\[
 \boxed{\mathcal L_{\mathbf p}(F^m)
 =m!(qm+1)\prod_i(p_i)_{q p_i m}\,[z^{qm+1}]b.}
\tag{9}
\]

Every scalar factor in (9) is a nonzero integer. At
`m=0`, both sides equal one. Every coefficient of `b` outside degrees
congruent to one modulo `q` is zero, by the same Lagrange formula.

The corresponding ordinary **marked** identity is

\[
 \mathcal L\!\left(\left(\prod_iY_i^{p_i-1}\right)F^m\right)
 =\left(\prod_i(p_i-1)!\right)\mathcal L_{\mathbf p}(F^m).
\tag{10}
\]

This marker is fixed as `m` varies. It is not an independent monomial factor
of `F^m` and cannot be deleted from the statement.

## 3. Nonzero critical values and the support-size window

Group the active indices by their distinct nonzero coefficient values
`d_1,...,d_h`, and put `P_j=sum_(c_i=d_j) p_i`. Thus

\[
 a(z)=z\prod_{j=1}^{h}(1-d_jz^q)^{P_j}.
\]

Its derivative factors exactly as

\[
\begin{aligned}
 a'(z)&=\prod_{j=1}^{h}(1-d_jz^q)^{P_j-1}\,D(z),\\
 D(z)&=\prod_{j=1}^{h}(1-d_jz^q)
       -qz^q\sum_{j=1}^{h}P_jd_j
                     \prod_{\ell\ne j}(1-d_\ell z^q).
\end{aligned}
\tag{11}
\]

Every critical point supplied by the first product maps to zero under `a`.
The remaining polynomial `D` has degree `qh`; in particular it has at most
`qh` distinct roots. Therefore

\[
 \#\operatorname{CV}^{\times}(a)\leq qh.
\tag{12}
\]

Suppose the `h` moments at `n,...,n+h-1`, with `n>=1`, all vanish. By (9)
and the structural zero classes, all inverse coefficients in

\[
 q(n-1)+2\leq\ell\leq q(n+h)
\]

vanish. This interval starts at degree at least two and has
`q(h+1)-1>=qh` terms. Equations (12) and Corollary 2 force `a=z`.
But a nonzero `F` has at least one active factor and `deg a>=q+1`, a
contradiction. This proves the asserted all-order `h`-window under `L_p`.
Repeated coefficient values and inactive indices are handled by the grouping
above, not by assuming generic coefficients. Inactive indices merely contribute
the nonzero scalar factors `(p_i)_(q p_i m)` in (9).

## 4. Exact ordinary block-sum transfer

For `S=X_1+...+X_p` in `p` independent ordinary variables,

\[
 \mathcal L(S^k)
 =\sum_{k_1+\cdots+k_p=k}\frac{k!}{k_1!\cdots k_p!}
       k_1!\cdots k_p!
 =k!\binom{k+p-1}{p-1}=(p)_k.
\tag{13}
\]

Disjoint variable blocks make their factorial moments multiply. Applying (13)
coefficientwise proves the exact identity

\[
 \mathcal L(G^m)=\mathcal L_{\mathbf p}(F^m)
 \qquad(m\geq0),
\tag{14}
\]

where `G` is (2). Hence `G` has the proved `h`-window.

The summand indexed by an active `i` in (2) has block-total-degree vector
`(q p_1,...,q p_r)+e_i`. These vectors are distinct for different `i`.
Within each summand, expansion coefficients are positive integers times the
single nonzero scalar `c_i`; thus no cancellation occurs within a summand
or between different active summands. In particular `N(G)>=s>=h`.
Every actual support-size window contains an initial `h`-window, so `G`
satisfies ordinary SFC. This conclusion permits arbitrary complex coefficients
and arbitrary characteristic-zero coefficient fields; positivity of those
coefficients is not assumed.

For an explicit check, the actual number of monomials is

\[
 \mathcal N(G)=\sum_{i:c_i\ne0}
    \binom{(q+1)p_i}{p_i-1}
       \prod_{j\ne i}\binom{(q+1)p_j-1}{p_j-1}.
\tag{15}
\]

When all `p_i=1`, this recovers
[POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md](POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md).
Independent linear additions and independent monomial factors can then be
adjoined by the already established transfer identities. No general nonlinear
substitution or identification of two existing variables is used.

## 5. The unmarked unequal-exponent gap is real

Take `q=2`, `(p_1,p_2)=(1,2)`, and

\[
 F=X^2Y^4(5X-3Y).
\]

Directly,

\[
 \mathcal L(F)=5\cdot3!4!-3\cdot2!5!=0,
 \qquad
 \mathcal L(YF)=5\cdot3!5!-3\cdot2!6!=-720.
\tag{16}
\]

Here `L_(1,2)(F)=L(YF)`. Thus even a single scalar zero for the ordinary
functional need not be a zero for the weighted functional. This is a scope
countercontrol, not an SFC counterexample. The proof above supplies no
propagation from ordinary zero windows to the marked windows in (10).

Unrestricted SFC still requires projective regularity of every support-size
moment system for every distinct integer support. The present exponent patterns
and their block-sum lifts form a special class. General unit-row inflation,
collision-free row identification, and the scalar-to-marked zero-window bridge
remain unproved here.

## 6. Verification and source boundary

The retained critical-value argument, the new exclusion of zero, the coefficient
identity, and the block-sum transfer have written proofs above and internal
cross-review. No literature-priority claim, external peer review, or formal
verification is asserted. The result does not rely on the support-uniform
eventual ordinary-binomial theorem or its stale source receipts. No existing
registry or status file is changed by this note.

Run the standard-library checker from the repository root:

```sh
python -B research/factorial/image-rigidity/verify_weighted_radial_linear_sfc_transfer.py \
  --repo-root . --output /tmp/weighted-radial-transfer.json
```

The checker independently expands the original weighted phases and ordinary
block-sum phases, inverts (6) by formal composition, verifies the complete
derivative factorization (11), and checks active/duplicate coefficients, the
window arithmetic, `m=0`, and the explicit wrong-functional control (16).
The [receipt](weighted_radial_linear_sfc_transfer_verification.json) records
the exact finite domains and source hashes. These checks support the written
proof; they do not establish its unrestricted quantifiers by extrapolation.
