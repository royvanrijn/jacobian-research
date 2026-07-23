# Quantum residue obstruction: definition, comparison maps, and first audit

## 1. Outcome

The residue idea has a valid obstruction-theoretic core, but the coefficient
category in the first formulation is too small.

Let

\[
 A=K[X,Q,Z],\qquad A_X=K[X,X^{-1},Q,Z].
\]

Formal solvability in an étale completion of a branch does **not** imply
solvability over \(A_X\).  The completion permits infinite series in inverse
branch coordinates; \(A_X\) permits only finite Laurent polynomials.  Thus
the sentence

\[
 \text{formal-local solvability}\Longrightarrow
 \text{vanishing in }C^\bullet_{\mathrm{loc}}
\]

is false when \(C^\bullet_{\mathrm{loc}}\) means localization only at \(X\).
At the known degree-five sample, an exact 16-term functional proves that the
full parity obstruction at \(\hbar^5\) remains nonzero over
\(K[X,X^{-1},Q,Z]\), including every Laurent \(\hbar^2\) gauge direction and
every Laurent \(\hbar^4\) correction.  This is the opposite of what would be
required to define its residue by an unavoidable finite negative-\(X\)
principal part.

The corrected invariant is a connecting class for a *specified comparison
algebra*.  There are at least two different candidates:

1. the finite-Laurent comparison \(A\subset A_X\), which detects genuine
   pole descent but need not kill the obstruction; and
2. the formal-étale comparison \(A\subset\widehat A_\xi^{\,\mathrm{et}}\),
   which kills every local defect but records an infinite formal tail rather
   than a Laurent principal part.

The exact experiment is
[`explore_degree_five_quantum_residue.py`](../scripts/explore_degree_five_quantum_residue.py).
The characteristic-zero function-field certificate at \(\hbar^3\) is
[`verify_degree_five_third_order_function_field.py`](../scripts/verify_degree_five_third_order_function_field.py).
The all-pole exact-\(\mathbb Q\) certificate is
[`verify_degree_five_laurent_quantum_obstruction.py`](../scripts/verify_degree_five_laurent_quantum_obstruction.py).

## 2. The correction complex

Fix a classical canonical pair \(S,T\in A\), with

\[
 \{S,T\}=1.
\]

At one filtered order, suppressing the order-dependent degree bounds, use

\[
 C^0(A)=A,\qquad C^1(A)=A\oplus A,\qquad C^2(A)=A
\]

and

\[
\begin{aligned}
 d_0H&=(\{H,S\},\{H,T\}),\\
 d_1(s,t)&=\{s,T\}+\{S,t\}.
\end{aligned}
\tag{QR1}
\]

Jacobi and \(\{S,T\}=1\) give \(d_1d_0=0\).  The actual filtered complex at
order \(n\) is the finite subquotient obtained from the declared
\(Z\)-orders and Bernstein bounds.  For the parity branch in degree five:

\[
\begin{array}{c|cc}
 &\operatorname{ord}_Z&\deg_B\\ \hline
 S_2&3&25\\
 T_2&2&21\\
 S_4&1&21\\
 T_4&0&17 .
\end{array}
\tag{QR2}
\]

Here \(\deg_BX=\deg_BQ=1\) and \(\deg_BZ=3\).

For Weyl symbols, write

\[
 \Pi=\partial_Z\otimes\delta-\delta\otimes\partial_Z,\qquad
 \delta=3X^2\partial_X+(2-6XQ)\partial_Q .
\]

The first two defects are

\[
\begin{aligned}
 O_3&=\frac1{24}\Pi^3(S,T),\\
 O_5(S_2,T_2)
 &=\{S_2,T_2\}
 +\frac1{24}\Pi^3(S_2,T)
 +\frac1{24}\Pi^3(S,T_2)
 +\frac1{1920}\Pi^5(S,T).
\end{aligned}
\tag{QR3}
\]

Thus \(d_1(S_2,T_2)=-O_3\), and then
\(d_1(S_4,T_4)=-O_5(S_2,T_2)\).

## 3. What the connecting class really is

For a chosen comparison algebra \(B\supset A\), suppose the filtered
complexes form a short exact sequence

\[
 0\longrightarrow C^\bullet(A)
 \longrightarrow C^\bullet(B)
 \longrightarrow C^\bullet(B)/C^\bullet(A)
 \longrightarrow0.
\tag{QR4}
\]

If a polynomial defect \(O\) is killed by a correction
\(\Delta_B\in C^1(B)\), then

\[
 [\overline{\Delta_B}]
 \in H^1\!\left(C^\bullet(B)/C^\bullet(A)\right)
\tag{QR5}
\]

maps under the connecting homomorphism to \([O]\in H^2(C^\bullet(A))\).
This is the precise residue construction.

There are two qualifications.

First, (QR5) is not generally one canonical element.  Changing the local
primitive changes it by the image of \(H^1(C^\bullet(B))\).  The intrinsic
object is the torsor

\[
 \operatorname{QRes}_B(O)=
 \delta^{-1}([O])
 \subset H^1\!\left(C^\bullet(B)/C^\bullet(A)\right),
\tag{QR6}
\]

or its image after an explicit gauge normalization.  At higher order it is
a Massey-type obstruction set because \(O_5\) depends quadratically on the
chosen solution of the \(O_3\) equation.

Second, the set (QR6) exists only if \([O]\) dies in \(H^2(C^\bullet(B))\).
Formal-étale exactness proves this for \(B=\widehat A_\xi^{\,\mathrm{et}}\).
It says nothing by itself for \(B=A_X\).

Consequently the safe conjecture is:

> **Comparison-quantization criterion.** A global filtered quantization
> exists exactly when there is a compatible tower of lower corrections for
> which every polynomial obstruction class vanishes.  After choosing a
> comparison algebra \(B\) in which those classes die, this is equivalent to
> the compatible connecting torsors (QR6) containing zero.

This is tautologically correct as obstruction theory.  Its value depends on
finding a comparison algebra and a finite decomposition that make the
classes computable.

## 4. Why the formal primitive is not a Laurent principal part

At a completed étale branch, \(S,T\) are formal coordinates and

\[
 d_1(f,g)=\partial_Sf+\partial_Tg.
\]

The homotopy

\[
 \mathsf h(F)=\left(-\int F\,dS,0\right)
\tag{QR7}
\]

kills every defect.  Pulling (QR7) back requires the formal inverse branch of
\((S,T)\).  There is no reason for that inverse, or the resulting primitive,
to be a finite Laurent polynomial in \(X\).

In fact the degree-five family is especially transparent at \(X=0\):

\[
 T|_{X=0}=-\frac32Q,\qquad
 S|_{X=0}=-\frac13Z+c(a,\tau)Q^2.
\tag{QR8}
\]

This is already a polynomial coordinate pair in \(Q,Z\).  The inverse
function theorem therefore produces an \(X\)-adic **regular power series**
near the divisor \(X=0\), not a negative-\(X\) principal part.  The missing
global property can be an infinite positive-\(X\) tail.  Localizing at \(X\)
does not capture that failure.

A comparison quotient such as

\[
 K[Q,Z][[X]]/K[X,Q,Z]
\tag{QR9}
\]

is therefore at least as natural as \(A_X/A\) for this chart.  It records
non-polynomial tails, but it is infinite-dimensional until an additional
weight/jet argument is supplied.

## 5. Full degree-five seed constructor

On the generic chart put

\[
 a=-\frac{1+\kappa}{2+\kappa},\qquad
 \kappa=-\frac{1+2a}{1+a},
\]

and

\[
 s_2=
 \frac{
 -12a(a+1)^2\tau^2
 +18a(a+1)(4a+5)\tau
 +216a^3+648a^2+738a+315
 }{28(a+1)^2}.
\tag{QR10}
\]

The new script constructs \(S,T\) over an arbitrary exact coefficient field
without first expanding a large rational SymPy expression.  On
\(\kappa\ne-1,-2\) it uses

\[
\begin{aligned}
 H(w)&=w^2(w-1)(\tau w^2+Aw+B),\\
 H'(w)&=w\,p_1(w),\\
 wH'(w)-H(w)&=w^2q_2(w)
\end{aligned}
\]

to cancel the apparent powers of \(X\) before sparse substitution.  This
gives the complete two-parameter family, verifies \(\{S,T\}=1\), and agrees
coefficient-for-coefficient over \(\mathbb Q\) with the old
\((\kappa,\tau)=(0,1)\) constructor.

The divisor \(\kappa=-1\) is not discarded.  The script separately implements
the replacement chart

\[
 R=2X+2X^3Q,\qquad
 \delta_{-1}=-2X^3\partial_X+(2+6X^2Q)\partial_Q,
\]

with the completing shear

\[
 W=Z+\frac{2(2\tau^2-15\tau-18)}{105}XQ.
\]

In this chart the adapted operator \(Z\) has Bernstein weight four.  The
parity correction bounds are therefore

\[
 (S_2,T_2,S_4,T_4):
 \quad(28,3),(24,2),(24,1),(20,0),
\tag{QR10a}
\]

where each pair is `(maximum Bernstein degree, maximum Z-order)`.  Direct
exact computation verifies the exceptional pair is canonical for
\(\delta_{-1}\).  Thus the constructor covers every normalized admissible
degree-five parameter \(\kappa\ne-2\), using the correct chart on
\(\kappa=-1\).

Every rank below is an exact rank in the displayed field, not a floating
point computation.

## 6. The full first-order quotient across the charts

The unrestricted first correction uses

\[
 (\deg_BS_1,\operatorname{ord}_ZS_1)=(27,4),\qquad
 (\deg_BT_1,\operatorname{ord}_ZT_1)=(23,3)
\]

on the generic chart.  At every tested \(\tau\ne0\) generic-chart seed,
including both \(\kappa=0\) and \(\kappa\ne0\), the exact calculation is

\[
 2132\text{ columns},\quad
 \operatorname{rank}d_1=2075,\quad
 \dim\ker d_1=57.
\]

The 19 target-Hamiltonian directions have the previously displayed basis,
so the quotient remains

\[
 \boxed{57-19=38}
\tag{QR10b}
\]

across the sampled generic chart, not only at \((0,1)\).

On the \(\kappa=-1\) replacement chart, the correct bounds are `(30,4)` and
`(26,3)`.  There are 2424 columns, rank 2374, and kernel dimension 50.  Since
\(\deg_BR=4\), the bounded target-Hamiltonian gauge now has

\[
 8+1+7=16
\]

directions, giving quotient dimension

\[
 \boxed{50-16=34.}
\tag{QR10c}
\]

This number is chart/filtration sensitive.  It is another warning that a
global quantum-residue invariant must specify how the filtered complexes on
the two classical descent charts are compared.

## 7. Polynomial \(QRes_3\) and \(QRes_5\)

The \(\hbar^3\) assertion is now uniform, rather than sampled.  Over
\(\mathbb Q(a,\tau)\), the 1569-variable generic-chart equation has a
solution on a fixed support of only 42 monomials: 27 in \(S_2\) and 15 in
\(T_2\).  The corresponding 42 columns have rank 42, the right-hand side
has 89 nonzero output equations, and exact substitution leaves zero
residual in all 89 equations.  Every coefficient denominator divides

\[
 a^4(a+1)^{10}.
\tag{QR10d}
\]

There are consequently no extra parameter poles: the solution specializes
at every point of the generic chart \(a\ne0,-1\), for every \(\tau\).
On the \(\kappa=-1\) replacement chart, a parallel computation over
\(\mathbb Q(\tau)\) uses 34 monomials, gives rank 34 against 57 output
equations, and has polynomial coefficients in \(\tau\).  It also has zero
residual.  Hence

\[
 \boxed{QRes_3=0}
\tag{QR10e}
\]

on the full normalized admissible degree-five parameter surface
\(\kappa\ne-2\), including the degree-drop divisor \(\tau=0\).

At \(\hbar^5\), the basic exact fiber strata for \(\tau\ne0\) are:

\[
\begin{array}{c|ccc|ccc}
 &\#C^1_3&\operatorname{rank}d_3&\dim\ker d_3
 &\operatorname{rank}d_5&
 \operatorname{rank}\langle d_5,O_5^{>0}\rangle&
 \operatorname{rank}\langle d_5,O_5^{>0},O_5^0\rangle\\ \hline
 \text{generic open}&1569&1528&41&594&642&643\\
 \kappa=0&1569&1527&42&594&646&647\\
 \kappa=-1\text{ replacement}&1853&1817&36&769&809&810
\end{array}
\tag{QR11}
\]

Here \(O_5^0\) is the constant term obtained from one affine solution of the
\(\hbar^3\) equation, and \(O_5^{>0}\) denotes all linear and quadratic
coefficient vectors in the complete \(\ker d_3\)-parameter family.  The
rank jump in the last column proves that the constant term is outside the
span of the correction image and every nonconstant gauge coefficient.
Therefore no value of the \(\hbar^3\) parameters can kill the
\(\hbar^5\) defect.

The \(\kappa=0\) row is the previously studied sample stratum; its nullity
42 is not the generic-chart nullity.  The replacement chart has a third
nullity, 36, and again has the final rank jump.  This is why both the full
parameter family and the missing chart had to be implemented.

The label “generic open” is essential.  Four parameterized periods found
below leave a cubic-field locus on which the last two generic ranks become
\(636\to636\).  There the span test is inconclusive, and the genuine
quadratic equations must be solved.  They are solved explicitly in
(QR19l)--(QR19m).

At \(\tau=0\), where the seed drops below degree five, the observed row is

\[
 \operatorname{rank}d_3=1500,\quad
 \dim\ker d_3=69,\quad
 594,\ 915,\ 915.
\tag{QR12}
\]

The linear-span certificate is inconclusive there, as it should be: this is
the degree-drop divisor, not the degree-five problem.

Interpretation:

* \(QRes_3=0\) in the polynomial comparison at every admissible seed, by
  the two characteristic-zero function-field certificates above.
* The complete polynomial \(\hbar^5\) obstruction is nonzero on both the
  generic sampled stratum, the special \(\kappa=0\) sampled stratum, and the
  exact rational sample on the \(\kappa=-1\) replacement chart.
* (QR11) is a fiberwise exact computation and a generic-open certificate,
  not yet a symbolic determination of every exceptional divisor.  A
  function-field row reduction or determinantal interpolation is still
  required before claiming the result for every \((\kappa,\tau)\).

## 8. Finite Laurent bands

For pole order \(N\), define

\[
 V(D,r;N)=
 \left\langle X^iQ^jZ^k:
 i\ge-N,\ 0\le k\le r,\ i+j+3k\le D\right\rangle .
\tag{QR13}
\]

This is exactly \(X^{-N}\) times a bounded polynomial numerator space.  The
script solves the \(\hbar^3\) equation in
\(V(25,3;N)\oplus V(21,2;N)\).  Since the quadratic term
\(\{S_2,T_2\}\) can have pole \(2N\), it tests the \(\hbar^5\) correction in
\(V(21,1;2N)\oplus V(17,0;2N)\).

At \((\kappa,\tau)=(0,1)\) over \(\mathbf F_{32003}\):

\[
\begin{array}{c|ccc|ccc}
N&\#C^1_3&\operatorname{rank}d_3&\dim\ker d_3&
\#C^1_5&\operatorname{rank}d_5&
\operatorname{rank}(B_N)\to\operatorname{rank}(B_N|O_5^0)\\ \hline
0&1569&1527&42&614&594&646\to647\\
1&1719&1675&44&741&721&773\to774\\
2&1876&1831&45&880&860&912\to913\\
3&2040&1994&46&1031&1011&1063\to1064
\end{array}
\tag{QR14}
\]

Here \(B_N\) contains the Laurent \(d_5\)-columns and every nonconstant
linear/quadratic coefficient vector from the full Laurent
\(\hbar^3\)-solution space.  Each final rank jump is a rigorous
nonexistence certificate in that pole band.  The rows \(N\le2\) repeat
unchanged over the good primes \(31991,32003,65521\).

On the replacement chart at \((\kappa,\tau)=(-1,1)\), the parallel rows for
\(N=0,1,2\) are

\[
\begin{array}{c|ccc|ccc}
N&\#C^1_3&\operatorname{rank}d_3&\dim\ker d_3&
\#C^1_5&\operatorname{rank}d_5&
\operatorname{rank}(B_N)\to\operatorname{rank}(B_N|O_5^0)\\ \hline
0&1853&1817&36&787&769&809\to810\\
1&2015&1976&39&930&912&956\to957\\
2&2184&2144&40&1085&1067&1111\to1112 .
\end{array}
\tag{QR15}
\]

The pole-band pattern at the known \((0,1)\) sample admits a finite all-pole
certificate.  Define a 16-term functional \(\Lambda\), supported on

\[
\begin{gathered}
(19,7,0),(18,2,2),(14,0,1),(18,4,1),\\
(17,1,2),(18,6,0),(16,0,2),(17,3,1),\\
(19,1,3),(17,5,0),(16,2,1),(19,3,2),\\
(15,1,1),(18,0,3),(19,5,1),(16,4,0),
\end{gathered}
\tag{QR16}
\]

with the exact rational coefficients recorded in the checker.  Its
\(X\)-support is the interval \(14\le i\le19\).  The two columns of the
\(\hbar^5\) linearized operator shift input \(X\)-degree by at most 11 and
12.  Hence every negative-\(X\) correction monomial is automatically
invisible to \(\Lambda\); checking the polynomial columns proves

\[
 \Lambda\circ d_5=0
 \quad\hbox{on the full finite Laurent union.}
\tag{QR17}
\]

For the lower lift, only \(S_2\)-terms of \(X\)-degree at least \(-8\) and
\(T_2\)-terms of degree at least \(-12\) can contribute to
\(\Lambda(O_5)\).  Inputs below degree \(-13\) cannot affect the projection
of the \(\hbar^3\) equation to output degree at least \(-1\).  The resulting
finite projected affine system has

\[
 4065\text{ variables},\qquad
 \operatorname{rank}=2990,\qquad
 \dim=1075.
\tag{QR18}
\]

On this deliberately larger affine space, exact rational calculation shows
that the quadratic matrix of \(\Lambda(O_5)\) is skew, all 1075 linear
derivatives vanish, and

\[
 \boxed{\Lambda(O_5)=-\frac{47547660815739}{190658}\ne0.}
\tag{QR19}
\]

Every finite Laurent \(\hbar^3\) solution restricts to a point of (QR18), so
(QR19) proves:

> **All-pole obstruction at the known seed.** At
> \((\kappa,\tau)=(0,1)\), the parity-preserving \(\hbar^5\) obstruction does
> not vanish after localization at \(X\).

Therefore the originally proposed \(QRes_5\in
H^1(C^\bullet_X/C^\bullet_{\rm poly})\) is not merely nonzero at this seed:
its connecting preimage is empty because the obstruction class remains
nonzero in \(H^2(C^\bullet_X)\).

The same finite-cutoff argument gives exact all-pole certificates on the
other two parameter strata:

\[
\begin{array}{c|c|c}
(\kappa,\tau)&\#\operatorname{supp}\Lambda&
\Lambda(O_5)\\ \hline
(1,1)&16&
\displaystyle\frac{18408906052532094299482611}
{3338052203018715136}\\[2mm]
 (0,-3)&15&1\\
(-1,1)&7&\displaystyle\frac{40998496}{305}.
\end{array}
\tag{QR19a}
\]

For the generic point \((1,1)\), the projected Laurent affine space again
has dimension 1075.  On the replacement chart at \((-1,1)\), the exact
finite superset has 5559 variables, rank 3627, and dimension 1932; the
quadratic functional matrix has eight entries and is skew, and every linear
variation vanishes.  These certificates prove that the all-pole obstruction
holds on a nonempty Zariski-open subset of the generic chart, on a nonempty
open of the special \(\kappa=0\) divisor, and on a nonempty open of the
\(\kappa=-1\) replacement divisor.  The vertical 15-term row additionally
proves the isolated four-period common point \((0,-3)\).  They do not yet
determine every complement of those opens.

The replacement-chart certificate is
[`verify_degree_five_exceptional_laurent_quantum_obstruction.py`](../scripts/verify_degree_five_exceptional_laurent_quantum_obstruction.py).

The invariant which actually survives these calculations is dual rather than
a connecting residue.  If

\[
 \Lambda\in H^2(C^\bullet_X)^\vee
\]

annihilates the linearized correction image and the variation of the
higher obstruction under every lower lift, define the **quantum Laurent
period**

\[
 \operatorname{QPer}_{n,\Lambda}
 :=\Lambda(O_n).
\tag{QR19b}
\]

It is gauge-independent by construction, and nonvanishing obstructs both
polynomial and \(X\)-Laurent quantization.  Equations (QR19) and (QR19a) are
explicit quantum Laurent periods.  Unlike the proposed `QRes`, they live
where the calculation says the obstruction actually remains: in the dual
of localized \(H^2\), not in a connecting \(H^1\)-torsor.

There is a parameter-uniform formulation which includes the vertical
phenomenon.  Let \(B\) be a bounded seed chart, let \(E_n\) be the finite
free module of order-\(n\) defect symbols modulo the order-\(n\) correction
image, and let

\[
 p:\mathcal L_{<n}\longrightarrow B
\]

be the affine scheme of lower-order lifts.  The order-\(n\) defect is a
section

\[
 o_n\in\Gamma(\mathcal L_{<n},p^*E_n).
\]

Define the **quantum-period sheaf**

\[
 \mathcal P_n=
 \ker\!\left(
 E_n^\vee\longrightarrow
 p_*\mathcal O_{\mathcal L_{<n}}/\mathcal O_B
 \right),
\qquad
 \Lambda\longmapsto \Lambda(o_n)\bmod\mathcal O_B.
\tag{QR19b'}
\]

Its sections are precisely the dual functionals whose value is independent
of the lower lift.  Evaluation defines the **linear quantum residue
section**

\[
 \operatorname{QRes}^{\mathrm{lin}}_n
 \in\Gamma(B,\mathcal P_n^\vee),
\qquad
 \operatorname{QRes}^{\mathrm{lin}}_n(\Lambda)=\Lambda(o_n).
\tag{QR19b''}
\]

This definition is finite and coherent in the bounded model, is
gauge-independent, and specializes to (QR19b) on a Laurent chart after the
all-pole cutoff has been proved.  Nonvanishing is a rigorous obstruction.
Vanishing is only necessary in general: \(o_n\) is quadratic in the lower
lift at order five, so linear functionals need not separate its image from
zero.  This is why the cubic component requires the actual quadratic solve
(QR19l)--(QR19m).

Formation of \(\mathcal P_n\) need not commute with base change.  The
natural map

\[
 \mathcal P_n\otimes k(b)\longrightarrow \mathcal P_{n,b}
\]

can fail to be surjective on a Fitting stratum.  The 15-term functional at
\((a,\tau)=(-1/2,-3)\) is exactly such a vertical class.  Thus a complete
computational invariant consists of

1. the coherent sheaf and evaluation section
   \((\mathcal P_n,\operatorname{QRes}^{\mathrm{lin}}_n)\);
2. its flattening/Fitting stratification and the extra fiberwise residue
   sections there; and
3. on the common zero scheme of all linear residues, the nonlinear
   zero-locus of \(o_n\).

This replaces the original single class in
\(H^1(C^\bullet_X/C^\bullet_{\rm poly})\).  It also refines the conjecture:
all linear residues vanish at every genuinely quantizable point, but the
converse requires the nonlinear zero-locus calculation and cannot be
asserted from dual periods alone.

There is also a useful exact reduction of the pending parameter-uniform
calculation.  Restrict a candidate functional to the fixed 16 monomials in
(QR16).  At each of the exact generic-chart samples

\[
(\kappa,\tau)=(1,1),\quad(2,3/2),\quad(-1/2,2),
\]

the complete bounded calculation produces only 41 nonzero conditions:
14 from the \(\hbar^5\) correction image, three first variations, and 24
quadratic lower-lift variations.  Their exact rank is 15.  Thus the bounded
period is the one-dimensional kernel of a \(15\times16\) subsystem.  Ten
image conditions and five lower-lift conditions suffice; the same pivot
pattern occurs at all three samples.  At \(\kappa=0\) the rank remains 15
and only the label of the last kernel direction changes because
\(\dim\ker d_3\) jumps from 41 to 42.

This reduction is reproduced by
[`explore_degree_five_fifth_order_period_constraints.py`](../scripts/explore_degree_five_fifth_order_period_constraints.py).
It replaces the bounded part of the parameter calculation by a small,
explicit symbolic system.

All 16 bounded lower-kernel directions which can meet (QR16) lift over
\(\mathbb Q(a,\tau)\), in nine sparse blocks.  Evaluating every variation
recovers exactly the 41 nonzero conditions seen at the samples.  Their
function-field rank is 15, so the five-condition pivot reduction did not
discard a parameter-dependent condition.  Rational row reduction suffers
severe coefficient swell, but a fraction-free Singular syzygy computation
gives exactly one kernel generator.  Contracting it with \(O_5\) and
factoring the numerator gives, up to a nonzero rational normalization,

\[
 \operatorname{QPer}^{\mathrm{bd}}_5(a,\tau)
 =
 \frac{\tau(a+1)^5}{a^5}\,P_{41}(a,\tau),
\tag{QR19c}
\]

where \(P_{41}\in\mathbb Q[a,\tau]\) is irreducible, has total degree 41,
and has 381 terms.  Before factorization the numerator has 451 terms and
total degree 47; the chosen syzygy normalization has denominator \(4a^5\).
The computation is implemented by
[`compute_degree_five_fifth_order_function_field.py`](../scripts/compute_degree_five_fifth_order_function_field.py).
Its `--full-factors` option prints \(P_{41}\).

After saturating by the generic-chart and degree-drop factors
\(a(a+1)\tau\), the bounded period therefore leaves the genuine candidate
exceptional curve

\[
 P_{41}(a,\tau)=0.
\tag{QR19d}
\]

This does not assert that quantization exists on that curve.  It says that
this one 16-support period ceases to certify nonexistence there.  A second
period was extracted at a finite-field point of \(P_{41}=0\), by replacing
the support monomial \(X^{19}Q^7\) with \(X^{20}Z^4\).  Its complete
41-condition function-field calculation gives

\[
 \operatorname{QPer}^{\mathrm{bd},2}_5(a,\tau)
 =
 \frac{\tau^2(a+1)^4}{a^5}\,Q_{32}(a,\tau)
\tag{QR19e}
\]

up to a nonzero rational normalization, where \(Q_{32}\) is irreducible,
has total degree 32, and has 258 terms.  In particular
\(P_{41}\mathrel{\nmid}Q_{32}\); the two candidate exceptional curves have
only a finite intersection.

Two further supports resolve almost all of that finite intersection.  The
third support contains both \(X^{19}Q^7\) and \(X^{20}Z^4\), dropping
\(X^{19}Q^5Z\).  The fourth replaces \(X^{19}Q^7\) by
\(X^{20}Q^2Z^3\).  Their complete 41-condition calculations give

\[
\begin{aligned}
 \operatorname{QPer}^{\mathrm{bd},3}_5(a,\tau)
 &\doteq \frac{\tau^2(a+1)^4}{8a^5}\,R_{38}(a,\tau),\\
 \operatorname{QPer}^{\mathrm{bd},4}_5(a,\tau)
 &\doteq \frac{\tau(a+1)^4}{8a^5}\,S_{39}(a,\tau),
\end{aligned}
\tag{QR19h}
\]

where \(\doteq\) suppresses a nonzero rational scalar.  The irreducible
polynomials \(R_{38}\) and \(S_{39}\) have respectively 356 and 371 terms
and the indicated total degrees.

The ideal generated by \(P_{41},Q_{32},R_{38},S_{39}\) has dimension zero
over each of
\(\mathbf F_{31991},\mathbf F_{32003},\mathbf F_{65521}\), vector-space
dimension 44, and elimination degree 24.  Its \(a\)-eliminant has the same
factor pattern at all three primes:

\[
 a^{10}(2a+3)^3(2a+1)^3
 (11a^2+28a+18)
 (94a^3+335a^2+400a+160)^2.
\tag{QR19i}
\]

The exact characteristic-zero fiber gcds are

\[
\begin{aligned}
 \gcd_i F_i(a,0)
   &=a^8(2a+3)^3(11a^2+28a+18),\\
 \gcd_i F_i(0,\tau)
   &=\tau^3(6\tau^2+4\tau+9),\\
 \gcd_i F_i(-3/2,\tau)&=\tau^2,\\
 \gcd_i F_i(-1/2,\tau)&=(\tau+3)^2,
\end{aligned}
\tag{QR19j}
\]

where \(F_i=P_{41},Q_{32},R_{38},S_{39}\).  Thus \(a=0\) is the replacement
boundary, and \(\tau=0\), \(a=-3/2\), and
\(11a^2+28a+18=0\) belong to the degree-drop boundary calculation.  At the
remaining rational point

\[
 (a,\tau)=(-1/2,-3)
\]

the exact bounded ranks are \(646\to647\), so the fifth-order obstruction
is nonzero even though all four displayed periods vanish.

Replacing \(X^{19}Q^7\) by \(X^{20}Q^4Z^2\) gives a fifth complete
41-condition function-field chart,

\[
 \operatorname{QPer}^{\mathrm{bd},5}_5(a,\tau)
 \doteq \frac{\tau(a+1)^4}{8a^5}\,T_{45}(a,\tau),
\]

where \(T_{45}\) is irreducible with 485 terms.  It still vanishes at
\((-1/2,-3)\), but lowers the local common-scheme length there from five to
four.  This is not a contradiction with the rank obstruction.  At the
rational point the dual kernel acquires a vertical generator: an exact
15-term functional, supported on

\[
\begin{gathered}
(14,0,1),(15,1,1),(16,0,2),(16,2,1),(16,4,0),\\
(17,1,2),(17,3,1),(17,5,0),(18,2,2),(18,4,1),\\
(18,6,0),(19,1,3),(19,3,2),(19,5,1),(20,4,2),
\end{gathered}
\]

annihilates the complete bounded span and has normalized value
\(\Lambda(O_5)=1\).  The same finite-cutoff calculation used in (QR18)
shows that it annihilates every Laurent correction and is independent of
all 1075 Laurent lower-lift directions.  Hence the rational point is
obstructed even after localization at \(X\).  The vertical generator is not
the specialization of the generic one-dimensional syzygy; a scalar period
chart alone therefore misses base-change failure in the dual obstruction
sheaf.

There are no hidden characteristic-zero closed points.  Exact local standard
bases give the following residue-field degrees and local lengths:

\[
\begin{array}{c|c|c|c}
\text{component}&[k(p):\mathbb Q]&
\operatorname{length}\mathcal O_{V(F_1,\ldots,F_4),p}&
\text{contribution}\\ \hline
(0,0)&1&18&18\\
a=0,\ 6\tau^2+4\tau+9=0&2&3&6\\
(-3/2,0)&1&4&4\\
\tau=0,\ 11a^2+28a+18=0&2&1&2\\
(-1/2,-3)&1&5&5\\
\text{cubic component (QR19k)}&3&3&9
\end{array}
\]

The contributions sum to \(44\).  Reduction modulo 32003 gives an upper
bound of 44 for the characteristic-zero quotient length, so equality proves
that this list exhausts the complete zero-dimensional support, including
scheme length.  Consequently, after removing the replacement and
degree-drop boundaries, the only four-period common zeros are
\((-1/2,-3)\) and the cubic component (QR19k).

For the five-period ideal the rational local length is four and every other
entry is unchanged, so the exact and modular lengths both become 43.  The
closed support is unchanged, exactly as predicted by the vertical
functional.

The other surviving reduced component is

\[
\begin{aligned}
 94a^3+335a^2+400a+160&=0,\\
 8\tau+658a^2+1593a+976&=0.
\end{aligned}
\tag{QR19k}
\]

Every \(F_i\) reduces exactly to zero modulo this ideal over \(\mathbb Q\).
Put

\[
 K=\mathbb Q[a]/(94a^3+335a^2+400a+160)
\]

and use the second equation of (QR19k) for \(\tau\).  Projection of the
complete bounded fifth-order defect modulo \(\operatorname{im}d_5\) gives
680 nonzero coordinate equations in the 41 lower-lift parameters.  Their
quadratic coefficient rank is 34; cancellation of the quadratic terms
gives eight consistent independent linear consequences.  Singular computes
a proper 19-generator standard basis of dimension 27.

More strongly, in the exact row-reduced kernel basis \(u_1,\ldots,u_{41}\),
set every parameter to zero except

\[
\begin{aligned}
u_{34}={}&-\frac{101668771215}{2097152}a^2
          -\frac{487549466415}{4194304}a
          -\frac{18474132105}{262144},\\
u_{19}={}&\frac{243}{10485760}
 \left(122116896574a^2+292049112895a+176583624080\right).
\end{aligned}
\tag{QR19l}
\]

All 680 projected equations then vanish exactly.  Solving the unreduced
correction equation gives

\[
 \operatorname{rank}d_5=594,\qquad
 \dim\ker d_5=20,
\]

with a particular correction having 14 nonzero coefficients, and direct
re-expansion verifies

\[
 d_5(S_4,T_4)+O_5=0.
\tag{QR19m}
\]

Thus the bounded \(\hbar^5\) obstruction genuinely vanishes on the cubic
component; the failure of four periods there is not a defect of their
supports.

The explicit branch does stop at the next order.  The filtration allows no
\(S_6,T_6\) symbols: their permitted \(Z\)-orders would be negative.  The
order-seven defect can therefore vary only through the 20-dimensional
kernel of the fifth-order correction.  For (QR19l), those variations have
rank six, while adjoining the 27-term constant defect raises the rank to
seven.  In fact coefficient extraction at \(X^{18}\) annihilates every
variation and gives

\[
 [X^{18}]O_7=
 \frac{2189187}{83886080}
 \left(587583566a^2+1388701707a+831388850\right)\ne0
 \quad\text{in }K.
\tag{QR19n}
\]

Hence no choice among the 20 fifth-order corrections extends this branch
through \(\hbar^7\).  The fifth-order solution scheme is much larger:
its standard basis has dimension 27 in the 41 lower-lift coordinates.  On
the explicit line

\[
 u_{36}=r,\qquad u_{35}=-\frac43r,\qquad
 u_{25}=-\frac{32a+24}{9}r
\]

with (QR19l) unchanged, all 680 fifth-order equations still vanish.  The
point \(r=1\) is again obstructed at order seven with ranks \(6\to7\), now
by a two-term period.  This is evidence for an order-seven obstruction on
the component, but it is not yet a proof over the entire
27-dimensional fifth-order solution scheme.

The extra Laurent lower-lift directions remain to be controlled before any
of (QR19c), (QR19e), or (QR19h) is called a parameter-uniform Laurent
period.

The replacement chart is cleaner.  The seven-term support in the exact
\((\kappa,\tau)=(-1,1)\) certificate has a complete
\(\mathbb Q(\tau)\) bounded period with primitive numerator

\[
 p_5(\tau)=
 744\tau^5+13020\tau^4+17000\tau^3+
 3240\tau^2-5670\tau-2187.
\tag{QR19f}
\]

It is irreducible.  At roots of \(p_5\) modulo 32003 the full bounded rank
still jumps, and replacing \(X^{17}QZ^3\) by \(X^{18}QZ^4\) produces a
second parameterized period with primitive numerator

\[
 q_4(\tau)=
 34\tau^4+472\tau^3+341\tau^2-162\tau-81.
\tag{QR19g}
\]

After primitive polynomial normalization the two actual period values have
greatest common divisor \(\tau^2\).  Therefore, on the degree-five locus
\(\tau\ne0\), at least one bounded period is nonzero at every point of the
entire \(\kappa=-1\) replacement chart.  This is computed by
[`compute_degree_five_exceptional_fifth_order_function_field.py`](../scripts/compute_degree_five_exceptional_fifth_order_function_field.py).

## 9. A canonical formal-tail representative

The formal comparison does kill the defect, but its first failure is not a
negative-\(X\) term.  Fix the free-zero polynomial solution
\((S_2,T_2)\), put \(T_4=0\), and solve

\[
 \{S_4,T\}=-O_5
\tag{QR20}
\]

in \(K[Q,Z][[X]]\), using zero \(Z\)-integration constants.  Since

\[
 T|_{X=0}=-\frac32Q,\qquad
 \{-,T|_{X=0}\}=-3\partial_Z,
\]

(QR20) has a unique recursive primitive in this gauge.  It is regular and
begins

\[
\begin{aligned}
S_4={}&-\frac{17302696248868}{7203}X^{12}Z
-\frac{137721444722906}{2401}X^{13}QZ\\
&+\frac{69014023267862}{2401}X^{14}Z^2
+\frac{33989945008104}{16807}X^{14}Q^2Z+\cdots .
\end{aligned}
\tag{QR21}
\]

The allowed filtered space has \(\operatorname{ord}_Z S_4\le1\) and
\(\deg_BS_4\le21\).  Hence the first failure in this canonical formal gauge
is already

\[
 \boxed{[X^{14}Z^2]S_4=\frac{69014023267862}{2401}\ne0,}
\tag{QR22}
\]

a differential-order failure.  Continuing the recursion through \(X^{30}\)
also exhibits the growing positive-\(X\) formal tail.  There is no
negative-\(X\) principal part.

This representative depends on the free-zero lower lift, \(T_4=0\), and the
zero integration constants; it is not yet a gauge-independent invariant.
It is computed by
[`explore_degree_five_quantum_formal_tail.py`](../scripts/explore_degree_five_quantum_formal_tail.py).

## 10. Reproduction

Exact rational polynomial calculation at a generic seed:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_degree_five_quantum_residue.py \
  --kappa 1 --tau 1 --max-pole 0 --exact
```

Known sample through pole order three:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_degree_five_quantum_residue.py \
  --kappa 0 --tau 1 --max-pole 3
```

All-pole exact rational certificate:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_five_laurent_quantum_obstruction.py
```

Generic-chart all-pole point and replacement-chart all-pole point:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_five_laurent_quantum_obstruction.py --generic-point
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_five_exceptional_laurent_quantum_obstruction.py
```

Canonical formal tail through \(X^{30}\):

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_degree_five_quantum_formal_tail.py
```

Uniform characteristic-zero \(\hbar^3\) certificate on both charts:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_five_third_order_function_field.py
```

Exact reduction of the bounded \(\hbar^5\) period conditions:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/explore_degree_five_fifth_order_period_constraints.py
```

Generic function-field bounded period and factorization:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/compute_degree_five_fifth_order_function_field.py
```

Remaining generic periods and the complete replacement-chart pair:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/compute_degree_five_fifth_order_function_field.py --second-support
PYTHONPATH=scripts .venv/bin/python \
  scripts/compute_degree_five_fifth_order_function_field.py --third-support
PYTHONPATH=scripts .venv/bin/python \
  scripts/compute_degree_five_fifth_order_function_field.py --fourth-support
PYTHONPATH=scripts .venv/bin/python \
  scripts/compute_degree_five_fifth_order_function_field.py --fifth-support
PYTHONPATH=scripts .venv/bin/python \
  scripts/compute_degree_five_exceptional_fifth_order_function_field.py
```

Write the four residual factors, verify their exact cubic component, and
analyze their modular intersection:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/compute_degree_five_fifth_order_function_field.py \
  --factor-output /tmp/degree_five_P41.singpoly
PYTHONPATH=scripts .venv/bin/python \
  scripts/compute_degree_five_fifth_order_function_field.py \
  --second-support --factor-output /tmp/degree_five_Q32.singpoly
PYTHONPATH=scripts .venv/bin/python \
  scripts/compute_degree_five_fifth_order_function_field.py \
  --third-support --factor-output /tmp/degree_five_R38.singpoly
PYTHONPATH=scripts .venv/bin/python \
  scripts/compute_degree_five_fifth_order_function_field.py \
  --fourth-support --factor-output /tmp/degree_five_S39.singpoly
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_period_intersection.py \
  --verify-cubic-component /tmp/degree_five_{P41,Q32,R38,S39}.singpoly
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_period_intersection.py --exact-decomposition \
  /tmp/degree_five_{P41,Q32,R38,S39}.singpoly
```

Exact solution of the genuine quadratic \(\hbar^5\) system on the cubic
component:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --exact-cubic --groebner --seventh-order
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --exact-cubic --seventh-order --u36 1
```

Extract and verify the vertical all-pole period at the rational common
point:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/extract_degree_five_fifth_order_functional.py \
  --kappa 0 --tau=-3
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_five_laurent_quantum_obstruction.py \
  --rational-common-point
```

Regression test:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/test_degree_five_quantum_residue.py
```

## 11. Next finite calculation

The all-pole certificate rules out the proposed Laurent connecting class at
the known seed, (QR22) gives a canonical-gauge formal representative, and
(QR19l)--(QR19m) close the genuine bounded fifth-order equations on the
four-period cubic locus.  The remaining useful directions are:

1. **Formal-tail invariantization.** Determine the quotient of (QR22) by
   changes of the lower lift, formal Hamiltonian gauge, and homotopy.  The
   displayed coefficient is a representative, not yet a cohomology class.
2. **Dual-sheaf formulation.** The exact length comparison exhausts the
   period support, the rational point is obstructed, and the cubic component
   is solvable at order five.  The rational vertical functional proves that
   the dual kernel does not commute with specialization.  Formulate the
   period invariant as a coherent cokernel/dual-sheaf section with its
   Fitting strata rather than as one rational scalar.
3. **Parameter-uniform Laurent periods.** Enlarge the 41 bounded conditions
   to the finite all-pole lower-lift superspace used in (QR18).  This is the
   missing descent step between the bounded period stratification and a
   full statement in localized \(H^2\).
4. **Uniform seventh order on the cubic component.** The branch (QR19l) and
   a second point on its \(u_{36}\)-line are obstructed by \(6\to7\) rank
   jumps, with the one-term period (QR19n) at the first point.  Compute the
   order-seven residue on the complete 27-dimensional fifth-order solution
   scheme, including its other components and vertical strata.

Only after one of these steps should a scalar or finite vector be called
the quantum residue of the degree-five family.
