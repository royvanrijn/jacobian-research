# Filtered `A_2` quantization obstruction for rank-two descent

This note isolates the direct quantization test formerly appended to the all-degree classical descent theorem.  It concerns one filtered Weyl-symbol ansatz at an exact degree-five sample; it neither weakens nor extends the classical rank-two completion in [Rank-two symplectic descent](RANK_TWO_SYMPLECTIC_DESCENT.md).

Once a polynomial exact symplectic map

\[
 G=(R,T,D,S):\mathbb A^4\to\mathbb A^4
\]

is available, a direct lift to the second Weyl algebra must solve more than
the associated-graded problem.  Choose Weyl generators with

\[
 [D,R]=[S,T]=\hbar
\]

and all mixed commutators zero.  A filtered quantization seeks ordered lifts
`R_hbar,T_hbar,D_hbar,S_hbar` with the exact same relations.  The principal
symbols give `G`, but at filtration order `m-2` the commutator differs from
the Poisson bracket by an ordering cocycle.  Corrections must therefore be
found recursively from linear equations of the form

\[
 d_G(\Delta_m)=-\mathcal O_m,
\]

where `d_G` is the linearized commutator map and `O_m` is determined by lower
orders.  Nonvanishing classes in the corresponding cokernel are genuine
ordering obstructions.

The practical first test is the completed degree-five family: use symmetric
(Weyl) ordering, compute the six commutator defects through the first two
lower filtration levels, and solve simultaneously for corrections.  Success
would give a filtered endomorphism of `A_2`; injectivity follows from the
injective associated-graded map, while non-surjectivity would still require a
separate filtered argument.

## Ore reduction and exact localized Darboux coordinates

The adapted coordinates make three of the six quantum relations automatic.
Put

\[
 \delta=3X^2\partial_X+(2-6XQ)\partial_Q
\]

and quantize `B=K[X,Q,Z]` as the Ore algebra

\[
 [Z,a]=\hbar\delta(a),\qquad a\in K[X,Q].
\]

Then `R=2X-3X^2Q` is central in this Ore algebra.  Consequently PBW-ordered
quantizations of `S,T,f` satisfy

\[
 [R,S]=[R,T]=0,\qquad [E+f,R]=\hbar
\]

exactly.  Symmetric ordering in the original four Darboux variables hides
this simplification and creates avoidable mixed commutator defects.

After localizing at `X`, set

\[
 v=X^{-1},\qquad P=-Z/3,
 \qquad U=E+\frac{(3-Rv)v^2}{2}P,
\]

with the coefficient placed to the left of `P`.  Direct first-order
differential-operator calculation gives

\[
 [P,v]=[U,R]=\hbar,
 \qquad [P,R]=[U,v]=[U,P]=0.                         \tag{Q1}
\]

Thus the localized problem is a rank-one fiber Weyl problem in `(v,P)` over
the parameter `R`, followed by a quantum Hamiltonian connection in the
`R`-direction.  If corrected fiber operators satisfy `[S_h,T_h]=hbar` and
`D_h=U+A_h`, the remaining equations are

\[
 [A_h,S_h]=-\hbar\partial_RS_h,
 \qquad [A_h,T_h]=-\hbar\partial_RT_h,               \tag{Q2}
\]

whose compatibility follows by differentiating the canonical commutator.
The exact splitting certificate is
[`verify_rank_two_quantum_darboux.py`](../scripts/verify_rank_two_quantum_darboux.py).

## A degree-five parity obstruction

The Ore reduction makes the first direct test finite.  At the admissible
rational point

\[
 (\kappa,\tau)=(0,1),\qquad s_2=69/7,
\]

the classical fiber symbols have differential orders `(5,4)` and leading
Bernstein degrees `(29,25)`.  Fiber Weyl ordering is the finite PBW operator

\[
 \operatorname{Weyl}_{v,P}(F)
 =\exp\!\left(\frac{\hbar}{2}\delta\partial_Z\right)F.
\]

For parity-preserving symbols

\[
 S_h=S+\hbar^2S_2+\hbar^4S_4,
 \qquad T_h=T+\hbar^2T_2+\hbar^4T_4,                \tag{Q3}
\]

the `hbar^3` equation is

\[
 \{S_2,T\}+\{S,T_2\}=-\frac1{24}\Pi^3(S,T),        \tag{Q4}
\]

where `Pi=partial_Z tensor delta-delta tensor partial_Z`.  In the complete
filtered spaces

\[
 \begin{array}{c|cc}
       &\text{maximum }Z\text{-order}&\text{maximum Bernstein degree}\\ \hline
 S_2&3&25\\
 T_2&2&21\\
 S_4&1&21\\
 T_4&0&17
 \end{array}
\]

equation (Q4) has rank `1527`, nullity `42`, and exact rational solutions.
Keeping the full 42-parameter solution space, the `hbar^5` equation projects
to a 53-dimensional cokernel.  Exact row reduction over `Q` contains the
constant equation `1=0`, independently of all 42 parameters.  Hence:

\[
 \boxed{\text{The standard parity-preserving filtered Weyl-symbol lift is
 obstructed at }\hbar^5\text{ at }(\kappa,\tau)=(0,1).}
\]

The exact calculation is
[`explore_degree_five_a2_subprincipal.py`](../scripts/explore_degree_five_a2_subprincipal.py).
It also repeats the ranks modulo `32003` and obtains the same constant
obstruction.  This closes the originally proposed symmetric-ordering test;
it does **not** rule out lifts with odd powers of `hbar`, a nonstandard wider
filtration, or quantization as a Weyl bimodule rather than an algebra
endomorphism.  Because the fiber relation already fails in the declared
ansatz, the connection equations (Q2) and quantum pole-descent test are not
reached there.  The same exact audit finds that allowing odd powers begins
with a rank-`2075` linear equation on `2132` coefficients, hence a
57-dimensional first-order kernel.  Determining its nonlinear higher-order
branches is the remaining direct-endomorphism problem; it is not hidden in
the now-closed symmetric-ordering calculation.

