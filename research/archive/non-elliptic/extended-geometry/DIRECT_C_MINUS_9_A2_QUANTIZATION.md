# Direct quantization of the `c=-9` rank-two completion

## Result

Let `(R,T,D,S)` be the exact polynomial symplectic map in
[Quadratic ladder and rank-two Poisson audit](QUADRATIC_LADDER_AND_POISSON_AUDIT.md),
with the pole-cancelling shear `W=Z-9Q^2`.  The smallest direct
parity-preserving filtered quantization test fails by an explicit rational
obstruction.

This is a result for the actual `c=-9` completion.  It is distinct from the
degree-five sample studied in
[Filtered `A_2` quantization obstruction](RANK_TWO_FILTERED_QUANTIZATION_OBSTRUCTION.md).

## Ore/Weyl reduction

On `B=Q[X,Q,Z]`, put

\[
 \delta=3X^2\partial_X+(2-6XQ)\partial_Q,\qquad
 [Z,a]=\hbar\delta(a).
\]

Then \(R=2X-3X^2Q\) is central.  Since `delta` and `partial_Z` commute, the
Weyl-symbol commutator is governed by

\[
 \Pi=\partial_Z\otimes\delta-\delta\otimes\partial_Z
\]

and

\[
 \frac{[\operatorname{Weyl}(f),\operatorname{Weyl}(g)]}{\hbar}
 =
 \Pi(f,g)+\frac{\hbar^2}{24}\Pi^3(f,g)
 +\frac{\hbar^4}{1920}\Pi^5(f,g)+\cdots .             \tag{1}
\]

For the `c=-9` pair, set

\[
 W=Z-9Q^2,\qquad Y=Q-\frac{XW}{3},\qquad U=1+XY
\]

and

\[
\begin{aligned}
S&=\frac12\left(U^3W+Y^2U(4+3XY)\right),\\
T&=Y+3XU^2W+3XY^2(4+3XY).
\end{aligned}                                         \tag{2}
\]

Exact expansion gives

\[
\deg_Z(S,T)=(3,2),\qquad
\deg_B(S,T)=(15,11),\qquad \Pi(S,T)=1,                \tag{3}
\]

where `deg_B(X)=deg_B(Q)=1` and `deg_B(Z)=3`.

## The complete bounded correction space

Take the standard parity-preserving, order-lowering ansatz

\[
 S_\hbar=S+\hbar^2S_2,\qquad
 T_\hbar=T+\hbar^2T_2                              \tag{4}
\]

in Weyl symbols, with the complete natural spaces

\[
\begin{array}{c|cc}
 &\max\deg_Z&\max\deg_B\\ \hline
S_2&1&11\\
T_2&0&7.
\end{array}                                           \tag{5}
\]

There is no `hbar^4` correction in this order-lowering class: its permitted
`Z`-order would be negative.

The `hbar^3` commutator equation is

\[
 \{S_2,T\}+\{S,T_2\}
 =-\frac1{24}\Pi^3(S,T).                              \tag{6}
\]

The domain in (5) has 159 monomials.  Exact row reduction over `Q` gives

\[
 \operatorname{rank}(d_2)=149,\qquad
 \dim\ker(d_2)=10,                                    \tag{7}
\]

and (6) is soluble.  Retain all ten affine parameters.

At the next order the obstruction is

\[
\mathcal O_5=
\{S_2,T_2\}
+\frac1{24}\Pi^3(S_2,T)
+\frac1{24}\Pi^3(S,T_2)
+\frac1{1920}\Pi^5(S,T).                              \tag{8}
\]

Substitution of the complete ten-parameter solution of (6) gives the exact
identity

\[
 \boxed{[X^{12}]\,\mathcal O_5=-49.}                  \tag{9}
\]

The left side is independent of all ten parameters.  Thus coefficient
extraction

\[
 \Lambda(F)=[X^{12}]F
\]

is an explicit obstruction functional on the `hbar^3` solution torsor, and
\(\Lambda(\mathcal O_5)=-49\ne0\).  For the free-zero representative, the
entire surviving residue is

\[
 -49X^{12}+147X^{13}Q-\frac{441}{4}X^{14}Q^2
 =
 -49X^{12}\left(1-\frac32XQ\right)^2.                 \tag{10}
\]

This is a cokernel certificate, not a Gröbner-exhaustion report.

## Consequence and exact scope

Equation (9) proves:

\[
\boxed{\text{The natural parity-preserving filtered Weyl-symbol lift of
the exact `c=-9` pair is obstructed at }\hbar^5.}
\]

Because `[S_hbar,T_hbar]=hbar` already fails, no four-operator lift of
`(R,T,D,S)` exists in this bounded class; the connection equations involving
`D` need not be reached.

This does **not** prove that the polynomial symplectic map is unquantizable in
`A_2`.  In particular it does not exclude:

- odd powers of `hbar`;
- corrections that enlarge the Bernstein-degree or differential-order
  bounds;
- a different filtered presentation not preserving the Ore polarization.

The exact certificate is
[`verify_c_minus_9_parity_quantization_obstruction.py`](../scripts/verify_c_minus_9_parity_quantization_obstruction.py).
