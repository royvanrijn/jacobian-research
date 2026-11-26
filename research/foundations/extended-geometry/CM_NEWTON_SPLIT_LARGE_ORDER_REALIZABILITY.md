# Split traces with large residue orders

Finite-field examples with one nonordinary factor persist when
both ordinary residue eigenvalues, their ratio, and the integral
trace cofactor are required to have large multiplicative order.
The count is uniform in the failed coefficient prime and in a
growing auxiliary congruence. The examples at different residue
characteristics are finite-field isogeny classes.

## Simultaneous large orders in each trace box

Let \(F\) be a fixed totally real cubic field. At an odd rational
prime \(p\) split completely in \(F\), prescribe a prime
\(\mathfrak q\mid p\). Define traces by
\[
  0\ne a\in\mathfrak q,\qquad a\equiv a_0\pmod{M O_F},\qquad
  (\sigma_i(a)/\sqrt p)_i\in B\subset(-2,2)^3,
\]
where \(M\) is prime to \(p\) and \(B\) is an axis-parallel cube
of side \(w\). Write \(u_1,u_2\in\mathbb F_p^\times\) for the
residues of \(a\) at the other two coefficient primes, and let
\(b=N_{F/\mathbb Q}(a)/p\). For sufficiently large \(p\), one has
\(v_{\mathfrak q}(a)=1\), \(v_{\mathfrak q_i}(a)=0\) at the other two
primes, and \(0<|b|\le8\sqrt p<p\). These follow from
\(0<|N(a)|\le8p^{3/2}<p^2\).

**Theorem.** Fix \(\eta,\kappa\ge0\) with \(\eta+\kappa<1/6\),
and fix \(0<\epsilon<1/2-3(\eta+\kappa)\).
Uniformly for \(M\le p^\eta\), \(w\ge p^{-\kappa}\), and every
choice of \(\mathfrak q,a_0,B\), the number of these traces is
\[
  (1+o(1))\frac{w^3\sqrt p}{M^3\sqrt{|\operatorname{disc}(F)|}}.
\]
A proportion \(1-o(1)\) of them have all four orders
\[
 \operatorname{ord}(u_1),\quad \operatorname{ord}(u_2),\quad
 \operatorname{ord}(u_2/u_1),\quad
 \operatorname{ord}(b\bmod p)
\]
greater than \(p^{1/2-3(\eta+\kappa)-\epsilon}\).
Each trace gives a simple abelian threefold over \(\mathbb F_p\)
with real multiplication by \(F\), relative Frobenius polynomial
\(X^2-aX+p\), and exactly one nonordinary coefficient factor.

### The trace lattice and the individual residues

By the
[balanced ideal basis](CM_NEWTON_SPLIT_COFACTOR_NORM_AND_TRACE_COUNT.md#counting-candidate-traces-with-a-moving-auxiliary-level),
the embedding matrix \(C\) of a basis of \(\mathfrak q\) and its inverse
satisfy \(\|C\|\ll_F p^{1/3}\) and \(\|C^{-1}\|\ll_F p^{-1/3}\).
Put \(N=p^{1/6}\). In coordinates \(a\leftrightarrow x\in\mathbb Z^3\)
for this basis, \(\|x\|\ll_F N\). The congruence is one residue class
of \(x\) mod \(M\), since \(M\) is prime to \(p\). The real trace
region in those coordinates has diameter \(O_F(wN)\).
Counting fundamental parallelepipeds gives the stated main term
with error \(O_F((1+wN/M)^2)\), including removal of zero.
The error is uniform in the center of \(B\), since the boundary
layer estimate uses only its side length and the diameter of a
fundamental cell. Write \(L=wN/M\), which tends uniformly to infinity.

For either ordinary prime \(\mathfrak q_i\), the map
\(a\mapsto u_i\) is injective on the trace box even before imposing
the auxiliary congruence. Indeed two traces with the same residue
have difference in \(\mathfrak q\mathfrak q_i\), while
\(|N(a-a')|\le64p^{3/2}<p^2\) for sufficiently large \(p\).

For a fixed nonzero integer \(b\), the principal ideal \((a)\)
has norm \(p|b|\le8p^{3/2}\). The number of ideals of that norm is
at most the ternary divisor function at \(p|b|\). Each such ideal
has \(O_F((1+\log p)^2)\) generators in the trace box. To see this,
each absolute embedding of a qualifying generator is between
\(|b|/4\) and \(2\sqrt p\), by bounding the other two embeddings.
Relative to one such generator, the possible units have logarithms
in a box of side \(O(1+\log p)\) in the rank-two unit lattice.
Thus the multiplicity of each \(b\) is \(p^{o(1)}\), uniformly in
\(b\). Also, distinct integers \(|b|\le8\sqrt p\) reduce to distinct
residues for sufficiently large \(p\).

### The ratio and planes in the trace lattice

For a given ratio \(r\in\mathbb F_p^\times\), the condition
\(u_2=ru_1\) defines a two-dimensional plane \(W_r\) in
\(\mathfrak q/p\mathfrak q\). It contains the line
\(\mathcal K=pO_F/p\mathfrak q\), the kernel of the map to the two
ordinary residues. Suppose two linearly independent coordinate
vectors \(x,y\) for this ratio occur in the trace box. Their
cross product is nonzero modulo \(p\), since its entries are
\(O_F(N^2)\). Its primitive integer normal \(n\) also annihilates
\(W_r\) modulo \(p\). Every other qualifying coordinate vector
lies in their rational plane, since its determinant with \(x,y\)
is divisible by \(p\) and has size \(O_F(N^3)<p\).

The functional defined by \(n\) corresponds under the trace
pairing to some \(\lambda\in\mathfrak q^\vee\). Since \(n\) kills
\(\mathcal K\) modulo \(p\), this functional is integral on \(O_F\):
for \(u\in O_F\), its value at \(pu\in\mathfrak q\) is divisible
by \(p\). Therefore \(\lambda\) belongs to the fixed codifferent of
\(O_F\). The embedding vector of a nonzero element of this lattice
has norm bounded below in terms of \(F\). The relation
\(\sigma(\lambda)=C^{-t}n\) now gives \(\|n\|\gg_F N^2\).
For any nonzero integer vector in the rational plane, its cross
product with at least one of \(x,y\) is a nonzero integer multiple
of \(n\). Since \(\|x\|,\|y\|\ll_F N\), every nonzero integer
vector in the plane has length \(\gg_F N\). A packing bound gives
\(O_F(1)\) qualifying traces for such a ratio.

Otherwise all vectors for that ratio lie on an integer line
with primitive direction \(z\). Different ratios with nonempty
line fibers have different lines, since nonzero multipliers in
the trace box have absolute value \(O_F(N)<p\). On such a line
there are \(O_F(1+L/\|z\|)\) qualifying points. The multipliers
are in one class mod \(M\) because \(z\) is primitive, and the
intersection with the trace region has length \(O_F(wN)\).
Among \(Z\) distinct integer lines the sum of the inverse norms
of the primitive directions is \(O(Z^{2/3})\), by the bound
\(O(T^3)\) for the number of directions of length at most \(T\).
Hence any set of \(Z\) ratios accounts for
\[
                         O_F(Z+LZ^{2/3})
\]
qualifying traces.

### Applying the fiber bounds

The number of units in \(\mathbb F_p\) with multiplicative order
at most \(H\) is at most \(H\tau(p-1)\), by summing \(\varphi(d)\)
over divisors \(d\) of \(p-1\) with \(d\le H\). Take
\(H=p^{1/2-3(\eta+\kappa)-\epsilon}\). The individual residues and
the cofactor account for at most \(Hp^{o(1)}=o(L^3)\) exceptions.
For the ratio, \(Z\le H\tau(p-1)=o(L^3)\), so the last fiber bound
is also \(o(L^3)\). This proves the simultaneous assertion.

Finally \(a\) generates \(F\) because a rational integer in
\(\mathfrak q\) is divisible by \(p\) and any such nonzero integer
is outside the trace box. A root \(\pi\) of \(X^2-aX+p\) is a Weil
\(p\)-integer and generates the CM quadratic extension
\(F(\sqrt{a^2-4p})\). At \(\mathfrak q\) the discriminant has
valuation one, yielding one degree-two local place of slope
\(1/2\). At each ordinary coefficient prime the local quadratic
splits with slopes zero and one by the Newton polygon. The
Honda–Tate local invariants are all zero. The division algebra
is therefore the CM field itself and the simple isogeny class
has dimension three. [Waterhouse, Chapter 2, pp. 527–528](https://www.numdam.org/article/ASENS_1969_4_2_4_521_0.pdf).

For a fixed characteristic-zero threefold the Frobenius trace
is a single choice at each prime. The theorem counts the
available finite-field traces. Together with the
[density-zero exclusions for small residue orders](CM_NEWTON_GROWING_RESIDUE_ORDER_EXCLUSIONS.md),
this shows that simultaneous large orders leave a substantial class
of finite-field candidates for the global problem.
