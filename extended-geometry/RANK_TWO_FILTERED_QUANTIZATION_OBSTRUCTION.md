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
obstruction.  This closes the originally proposed symmetric-ordering test.
Because the fiber relation already fails in the declared ansatz, the
connection equations (Q2) and quantum pole-descent test are not reached
there.

## Odd first corrections and the gauge quotient

The full standard filtered symbol ansatz at this sample is

\[
 S_\hbar=\sum_{n=0}^{5}\hbar^nS_n,\qquad
 T_\hbar=\sum_{n=0}^{4}\hbar^nT_n,                 \tag{Q5}
\]

with

\[
 \begin{aligned}
 \operatorname{ord}_Z(S_n)&\leq5-n,&
 \deg_B(S_n)&\leq29-2n,\\
 \operatorname{ord}_Z(T_n)&\leq4-n,&
 \deg_B(T_n)&\leq25-2n .                           \tag{Q6}
 \end{aligned}
\]

Here `deg_B(X)=deg_B(Q)=1` and `deg_B(Z)=3`.  No parity condition is imposed.
The coefficient of `hbar` in the normalized commutator is

\[
 d(S_1,T_1)=\{S_1,T\}+\{S,T_1\}=0.                \tag{Q7}
\]

In the complete spaces prescribed by (Q6), `d` has `2132` columns, exact
rank `2075`, and kernel dimension `57` over `Q`.

Nineteen of these directions are target-Hamiltonian gauge.  Since `R` is
fiber-central, the admissible Hamiltonians are

\[
 H=\sum_{a=0}^{9}\alpha_aR^aT+\frac{\gamma}{2}T^2
   +\sum_{a=0}^{7}\beta_aR^aS,                     \tag{Q8}
\]

and their first variations are

\[
 (\{H,S\},\{H,T\})=
 \left(-\sum_{a=0}^{9}\alpha_aR^a-\gamma T,
             \sum_{a=0}^{7}\beta_aR^a\right).     \tag{Q9}
\]

These 19 vectors are independent over `Q`.  Thus the essential first-order
space has dimension `38`.  For a reproducible coordinate slice, order the
filtered monomials by increasing `Z`-order, then `X`- and `Q`-exponent, put
the `S_1` columns before the `T_1` columns, and normalize each null vector at
its free column.  Row reduction of (Q8) has pivot indices

```text
0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,18,20,22,24
```

inside the 57-vector null basis.  Setting those coordinates to zero gives
the quotient basis with indices

```text
13,15,17,19,21,23,25,26,27,28,29,30,31,32,33,34,35,36,37,
38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56.
```

This is an exact parameterization over `Q`, not just a dimension count.

## Parked wider-ansatz calculation (`OP-QO`)

The continuation beyond the proved bounded ansatz is tracked only as `OP-QO`
in [STATUS.md](../STATUS.md).  At the next order one must solve

\[
 d(S_2,T_2)=
 -\{S_1,T_1\}-\frac1{24}\Pi^3(S,T).                \tag{Q10}
\]

The last term already lies in the image of `d` by (Q4).  Consequently the
well-defined quadratic obstruction on the gauge quotient is

\[
 \kappa_2([S_1,T_1])=
 \bigl[\{S_1,T_1\}\bigr]\in\operatorname{coker}(d_2). \tag{Q11}
\]

Exact rational projection of its `741` quadratic coefficient vectors gives
`934` nonzero cokernel coordinates spanning a 41-dimensional space.  The
same rank is obtained modulo `32003`.  The zero locus is not the origin:
in the quotient coordinates above, the ten axes

```text
0,1,2,3,8,9,10,23,24,35
```

survive (Q10), and there are mixed-support solutions as well.  Hence the odd
problem cannot be dismissed by linear rank or by the parity obstruction.

For each of the ten rational axis branches, the next coefficient can still
be tested without choosing a second-order gauge.  If
`(S_1,T_1)=lambda(s_1,t_1)`, put `u=lambda^2`, divide the `hbar^3` equation by
`lambda`, and solve it simultaneously with (Q10).  The resulting system is
affine linear in `(S_2,T_2)`, `(S_3/lambda,T_3/lambda)`, and `u`.  Exact row
reduction over `Q` is inconsistent for all ten axes.  Thus none of these
simple branches reaches the third correction.  This does not yet classify
all mixed components of the 41-quadric locus, so it is not a global
nonexistence theorem for a polynomial filtered Weyl endomorphism.

The exact kernel, gauge slice, quadratic map, and axis continuation are
computed by
[`explore_rank_two_odd_quantization.py`](../scripts/explore_rank_two_odd_quantization.py).

## A surviving formal-local quantization

The localized Darboux reduction does, however, construct an all-orders lift
once the coefficient ring is changed from filtered polynomials to the
formal etale completion.  Fix `R` and complete at any point of a fiber.  As
`{S,T}=1`, the formal inverse-function theorem makes `(S,T)` formal
coordinates.  In these coordinates the linearized normalized commutator is

\[
 d(f,g)=\partial_Sf+\partial_Tg.                    \tag{Q12}
\]

It has the explicit contracting homotopy

\[
 \mathsf h(F)=\left(-\int F\,dS,0\right),
 \qquad d\mathsf h(F)=-F.                           \tag{Q13}
\]

Suppose `S_hbar,T_hbar` have been constructed through order `hbar^(n-1)`
and let `F_n` be the coefficient of the remaining normalized commutator
defect.  Adding `hbar^n h(F_n)` kills that coefficient.  Formal integration
is coefficientwise and is always defined in characteristic zero.  Induction
therefore gives

\[
 [S_\hbar,T_\hbar]=\hbar
\]

to every order in the completed local Ore/Weyl algebra.

The construction need not start on the parity branch.  Any of the 57
solutions of (Q7), including any of the 38 quotient representatives, is a
divergence-free formal vector field and hence is formally Hamiltonian after
completion.  Prescribing it as `(S_1,T_1)` and then applying (Q13) at every
higher order produces an odd-in-`hbar` formal lift.  Thus every first-order
branch survives in the formal-local category, even though the ten simplest
polynomial filtered branches above do not survive their next bounded
equation.

The `R`-connection is unobstructed there as well.  Differentiating the exact
canonical commutator says that the `R`-variation is a divergence-free formal
vector field in `(S_hbar,T_hbar)`.  The formal Poincare lemma supplies a
Hamiltonian `A_hbar`; then `D_hbar=U+A_hbar` solves (Q2).  This yields a
homomorphism from the two-pair formal Weyl algebra to the completed localized
source algebra.  Equivalently, its graph gives a rank-one formal Weyl
bimodule.

The contracting-homotopy and connection primitives are checked in
[`verify_rank_two_formal_local_quantization.py`](../scripts/verify_rank_two_formal_local_quantization.py).
This construction is the first surviving formal quantization beyond the
parity ansatz.  It deliberately does **not** claim pole descent, bounded
Bernstein degree, or a global filtered endomorphism of `A_2`; those are the
global conditions not controlled by the formal contracting homotopy.
