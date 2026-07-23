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

## Low-support mixed odd branches

The coordinate axes do not describe the quadratic Maurer--Cartan locus.
Write \(x_0,\ldots,x_{37}\) for the ordered quotient coordinates displayed
above.  Exact restriction of the 41 quadratic obstruction equations to
coordinate subspaces gives the following sharper picture.

First, monomial equations force

\[
 x_6=x_7=x_{19}=x_{34}=0
\]

on the reduced locus.  Exactly ten coordinate axes survive, as already
recorded.  Among those ten axes, the coordinate-line graph has two maximal
linear subspaces:

\[
\begin{aligned}
 L_1&=\langle x_0,x_1,x_2,x_3,x_8,x_9,x_{10}\rangle,\\
 L_2&=\langle x_0,x_1,x_9,x_{10},x_{23},x_{24},x_{35}\rangle .
\end{aligned}
\tag{Q11a}
\]

Both are seven-dimensional vector spaces contained in the quadratic locus.
Outside their 36 coordinate lines, the complete exact-support-two audit
finds seven isolated rational directions and nine quadratic closed points
(nine conjugate pairs geometrically).  Exact row reduction over each
quadratic residue field shows that all nine quadratic classes, as well as
all seven rational directions, are already inconsistent at \(u=0\) in the
next simultaneous correction equation.  Every case has nullity zero; the
quadratic classes all have augmented rank \(2574\), while the rational ranks
range from \(2571\) to \(2576\).

There is a useful uniform relaxation on \(L_1\) and \(L_2\).  Allow the
second-order kernel parameters and every quadratic particular solution to
vary independently before imposing the factorization of their coefficients.
Exact row reduction over \(\mathbb Q\) gives

\[
\begin{array}{c|c|c}
&\text{relaxed rank jump}&\text{necessary first-order subspace}\\ \hline
L_1&1071\longrightarrow1076&
\left\langle x_9+\frac95x_0,\,
x_{10}+\frac{26}{15}x_0\right\rangle,\\[1mm]
L_2&1056\longrightarrow1058&
\left\langle
x_9+\frac95x_0,\,
x_{10}+\frac{26}{15}x_0,\,
x_{23}-\frac{37696}{1575}x_0,\,
x_{24}-\frac{142216}{4725}x_0,\,
x_{35}+\frac{205276828}{165375}x_0
\right\rangle .
\end{array}
\tag{Q11b}
\]

The first two mixed directions in (Q11b) genuinely survive the equation
which kills every coordinate axis.  For either direction, the simultaneous
second/third correction matrix has rank \(2572\), nullity \(63\), and leaves
\(u=\lambda^2\) free.  Thus these are not artifacts of the linear
relaxation.

They nevertheless fail at the following bounded level.  For either fixed
mixed direction, retain \(u\), all 63 compatible lower-lift parameters, all
2079 linear and quadratic coefficients of the next obstruction, and every
bounded fourth correction.  In the normalized commutator this is order four,
equivalently the \(\hbar^5\) coefficient of the unnormalized commutator.  The
span rank is

\[
 626\longrightarrow627
\tag{Q11c}
\]

when the constant obstruction is appended, over each of
\(\mathbf F_{31991}\), \(\mathbf F_{32003}\), and
\(\mathbf F_{65521}\).  All rational denominators are units at these primes.
The rank jump after good reduction proves that the characteristic-zero
constant is outside the span as well.  Since the actual lower-lift image is a
subset of this deliberately enlarged coefficient span, no amplitude and no
choice of compatible lower corrections continues either of these two mixed
directions.

The same calculation closes the rest of exact support two.  Intersecting the
five-dimensional residual space for \(L_2\) with every coordinate line
leaves only

\[
 x_9-\frac{27}{26}x_{10}.
\]

It survives the simultaneous equation with the same nullity \(63\), but its
all-scale next-obstruction rank is again \(626\to627\) at all three primes.
Every other residual coordinate-line intersection forces \(u=0\).
Together with the axes and isolated points above, this proves that **every
exact-support-two branch of this symbol fails by normalized order four**.
The calculation is
[`explore_rank_two_odd_mixed_quantization.py`](../scripts/explore_rank_two_odd_mixed_quantization.py).

The whole projective line spanned by

\[
 d_1=x_9+\frac95x_0,\qquad
 d_2=x_{10}+\frac{26}{15}x_0.
\]

is eliminated.  Over \(\mathbb Q(r)\), a fixed three-term functional

\[
 \mathcal R(F)=
 \frac{[X^{18}Z^3]F}{211680}
 +\frac{19[X^{14}Z]F}{1270080}
 +\frac{[X^{16}Z^2]F}{105840}
\tag{Q11d}
\]

annihilates all 614 fourth-correction columns and all 2079 enlarged
lower-lift coefficients, while \(\mathcal R\) of the constant obstruction
is \(1\).  The chosen lower-lift basis has only one parameter pole,
\(4r+3=0\), with multiplicity two.  At \(r=-3/4\) the lower-lift nullity
jumps from \(63\) to \(74\); a direct exact calculation there gives
\(646\to647\).  The two endpoints are already among the fixed directions
above.  Hence there are no exceptional projective parameters.  The exact
checker is
[`verify_rank_two_odd_mixed_function_field.py`](../scripts/verify_rank_two_odd_mixed_function_field.py).

## Exact support-three atlas and complete continuation

There are \(\binom{38}{3}=8436\) coordinate projective planes.  Saturating
the exact restricted quadratic ideals by all three coordinates gives:

\[
\begin{array}{c|r}
\text{chart type}&\text{number}\\ \hline
\text{plane contained in }L_1\cup L_2&66\\
\text{empty by an exact monomial equation}&7924\\
\text{empty after exact saturation}&268\\
\text{positive-dimensional}&149\\
\text{zero-dimensional}&29.
\end{array}
\tag{Q11e}
\]

The positive-dimensional charts comprise 52 line and 97 conic charts.  The
zero-dimensional schemes have length profile
\(1\cdot2+19\cdot3+9\cdot4=95\); after factoring their residue polynomials,
they give 38 rational and 28 quadratic closed-point classes.  This count is
characteristic-zero complete.  A finite-field pass at \(31\) is used only
as an independent profile: it misses 28 of the conic charts and therefore
is not used as a completeness certificate.

Of the 66 isolated closed-point classes, 65 are inconsistent already at
\(u=0\).  The sole survivor is the rational double point

\[
 x_4+\frac{126}{1963}x_{17}-\frac{63}{3926}x_{18}.
\tag{Q11f}
\]

Its simultaneous system leaves \(u\) free with rank \(2561\) and nullity
\(74\).  At the following bounded equation, retaining all 74 lower-lift
parameters, all 2849 enlarged obstruction coefficients, and every fourth
correction gives

\[
 646\longrightarrow647
\]

over each of the same three good primes.  Thus every isolated
exact-support-three branch is eliminated.  The atlas and continuation are
checked independently by
[`explore_rank_two_odd_support_three.py`](../scripts/explore_rank_two_odd_support_three.py)
and
[`verify_rank_two_odd_support_three_points.py`](../scripts/verify_rank_two_odd_support_three_points.py).

The 149 positive-dimensional charts can also be closed.  Quotient the
third-order equation by its fixed third differential once.  Of the resulting
2295 cokernel coordinates, only 920 occur in the 42 kernel couplings and
affine right sides.  Exact row reduction on each coordinate plane gives the
following dimensions for the necessary first-correction subspace:

\[
\begin{array}{c|ccc}
\text{necessary vector dimension}&0&1&2\\ \hline
\text{number of line/conic charts}&18&110&21.
\end{array}
\tag{Q11g}
\]

No necessary projective line is contained in its quadratic curve.  After
exact torus saturation, the 149 curves therefore leave only 23 closed-point
classes: four rational and nineteen quadratic.  All nineteen quadratic
classes force \(u=0\), with rank \(2571\) and nullity \(64\).  The four
rational directions are

\[
\begin{aligned}
&\frac{331}{126}x_1+x_{11}-\frac12x_{12},\\
&-\frac{2948}{189}x_2-\frac83x_{13}+x_{14},\\
&-\frac{305}{9}x_3-\frac{10}{3}x_{15}+x_{16},\\
&\frac23x_8+x_{21}-\frac34x_{22}.
\end{aligned}
\tag{Q11h}
\]

Each leaves \(u\) free with nullity \(74\), and each then has exact enlarged
rank \(646\to647\).  Consequently every exact-support-three branch outside
the 66 coordinate planes contained in \(L_1\cup L_2\) is eliminated.  The
\(L_1\) residual line was closed above, and the \(L_2\) planes are closed by
the residual audit below.  The compressed atlas and continuation are
[`explore_rank_two_odd_support_three_curves.py`](../scripts/explore_rank_two_odd_support_three_curves.py)
and
[`verify_rank_two_odd_support_three_curves.py`](../scripts/verify_rank_two_odd_support_three_curves.py).

## The residual five-space

The remaining necessary space inside \(L_2\) is the five-space in (Q11b).
Write its ordered basis as \(e_0,\ldots,e_4\), in the order displayed there,
and write a direction as \(z_0e_0+\cdots+z_4e_4\).  After quotienting by the
third differential, the 42 kernel-coupling columns depend on only the two
linear forms

\[
\begin{aligned}
 a={}&z_0+\frac43z_1+\frac{2648}{189}z_3
       -\frac{580504}{735}z_4,\\
 b={}&z_2+2z_3-\frac{9838}{105}z_4.
\end{aligned}
\tag{Q11i}
\]

Their union has rank \(14\), while a nonzero pencil member has rank \(11\).
The cubic column carrying the nonzero scale has 49 nonzero coordinates
modulo that union, and all 49 are constant multiples of one polynomial

\[
 b\,Q(z_0,\ldots,z_4),
\tag{Q11j}
\]

where \(Q\) is the explicit primitive quadric printed by
[`explore_rank_two_odd_residual_five_space.py`](../scripts/explore_rank_two_odd_residual_five_space.py).
The three remaining pencil-cokernel functionals vanish identically after
putting \(a/b\) into the generic pencil.  On \(b=0\) the corresponding fixed
pencil functionals also vanish, and on \(a=b=0\) the cubic column itself is
zero.  Thus this is an equality, not merely a necessary factor:

\[
 \{\text{nonzero-scale directions in the residual }\mathbf P^4\}
 =V(b)\cup V(Q).
\tag{Q11k}
\]

This factorization also finishes the contained support-three charts.  An
exact-support-three direction in this residual space either has two nonzero
\(z_i\) and nonzero \(x_0\), or three nonzero \(z_i\) with their \(x_0\)
coefficients cancelling.  Intersecting all twenty such projective lines with
\(V(bQ)\), then saturating by the three support coordinates, gives:

- the line \(\langle e_0,e_1\rangle\), already eliminated by the complete
  residual-line certificate;
- three rational and nine quadratic closed-point classes, all on residual
  coordinate pairs;
- no class on an \(x_0\)-cancelled residual triple.

All twelve closed points leave \(u\) free with nullity \(63\), and exact
number-field row reduction gives \(626\to627\) for each.  This is checked by
[`verify_rank_two_odd_residual_support_three.py`](../scripts/verify_rank_two_odd_residual_support_three.py).
Together with the atlas above, this proves that **every exact-support-three
branch of the displayed symbol is eliminated**.

The pointwise fourth-order gap on the residual five-space can in fact be
closed uniformly.  First, the projected affine third-order equation has one
second correction on the whole residual \(\mathbf P^4\): in the ordered
42-dimensional second-kernel basis it has only

\[
 k_{18}=-\frac{710122059}{392},
 \qquad
 k_{33}=-\frac{922185}{56}.
\tag{Q11l}
\]

The resulting zero-scale fourth defect is independent of \(z\).  The fixed
three-monomial residue from the complete residual line evaluates to \(1\)
on it.  Over the generic coupling pencil, and also on the \(b=0\) chart,
this residue kills the linear and quadratic fourth-order variation on all
31 compatible second-kernel directions.

There is a further structural collapse:

\[
 Q=
 1093955625a^2-42674688000ab+416287859400b^2,
\tag{Q11m}
\]

whose discriminant is
\(-2(486202500)^2\).  Thus over \(\mathbf Q(\sqrt{-2})\), \(V(Q)\) is the
union of the two conjugate hyperplanes

\[
 \frac ab=\frac{2048}{105}\mathbin{\pm}\frac29\sqrt{-2}.
\tag{Q11n}
\]

On \(b=0,a\ne0\), and on either hyperplane in (Q11n), normalize the nonzero
linear form.  The coupling becomes constant and the remaining three
coordinates parametrize an affine three-space.  Exact polynomial
right-inverse calculation shows that the same three-term residue kills the
unique nonzero-scale lower direction, its diagonal term, all 31 cross terms
with compatible second-kernel directions, and all 31 third-gauge terms,
identically in those three parameters.

The intersection \(a=b=0\) is a separate rank-zero projective plane: all 42
second-kernel directions become compatible and the lower nullity rises from
63 to 74.  A second, explicitly recorded 16-monomial residue has value \(1\)
on the same constant defect and annihilates the 614 fourth-correction
columns, all 42 second-kernel freedoms, all 31 third gauges, the scale
direction, and every cross term, identically in three homogeneous plane
coordinates.  These identities are checked by
[`verify_rank_two_odd_residual_fourth_identity.py`](../scripts/verify_rank_two_odd_residual_fourth_identity.py).
Consequently **every nonzero-scale branch in the residual
\(\mathbf P^4\) is obstructed at fourth order**, including all higher-support
points of that residual family.

The remaining unrestricted-odd problem for this symbol is now solely to
classify components of the original 41-quadric locus whose general point
has essential support at least four outside \(L_1\cup L_2\), then continue
any such components through the filtered equations.  Even eliminating all
of them would prove only that this classical symbol does not furnish a
counterexample; it would not prove \((DC_2)\).

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

The low-support calculation above has exactly the same logical scope.  It
eliminates additional branches of this one classical symbol in the declared
global filtered spaces.  Even a complete nonexistence proof for this symbol
would not prove `(DC_2)`.  Conversely, one genuine global filtered lift,
together with the required nonsurjectivity argument, would disprove
`(DC_2)`.

The follow-up [quantum-residue comparison audit](QUANTUM_RESIDUE_OBSTRUCTION.md)
defines the relevant connecting torsor and tests the full degree-five seed
constructor.  It also records an essential correction to the naive residue
proposal: formal-étale solvability does not imply solvability after
localizing only at `X`.  An explicit 16-term exact rational functional now
proves that the known `hbar^5` obstruction survives the full localization
`Q[X,X^-1,Q,Z]`, not merely finitely many pole bands.
