# A count saving for small square kernels on all three letters

For a positive integer \(T\), let
\(\operatorname{sf}_2(T)=\prod_{q:\ v_q(T)\ {\rm odd}}q\)
be its squarefree kernel. There are absolute
\(f_*,\delta,\lambda_*>0\) such that for every fixed
\(0<\lambda<\lambda_*\) the number of primitive \(a+b=c\le X\)
satisfying
\[
 \operatorname{rad}(abc)<c^\lambda,\qquad
 \prod_{T\in\{a,b,c\}}\operatorname{sf}_2(T)\le c^{f_*\lambda}
 \tag{K1}
\]
is \(O_{\lambda,\varepsilon}(X^{(1/3-\delta)\lambda+\varepsilon})\).
There is no short-letter condition and no exact valuation divisibility
requirement beyond the prescribed parity defect. In particular even
actual valuation masses may be distributed over many different
degrees. The argument uses the actual three-letter height and radical
constraints. Larger parity defects require additional weight, as in the selector below.

## A conic of controlled coefficient height

Work on a dyadic shell \(c\asymp X\). Write
\[
            a=A x^2,\qquad b=B y^2,\qquad c=C z^2
\]
with \(A,B,C\) squarefree. There are
\(X^{f_*\lambda+o_X(1)}\) triples of kernels in (K1), since
their product and its ordered prime allocation determine them.
Fix one such triple. The roots have no common divisor, and
their rational projective point lies on the nonsingular conic
\(Q(x,y,z)=A x^2+B y^2-C z^2=0\).

Set \(M=\max(A,B,C)\). Cassels' small-zero theorem gives
a nonzero integral zero \(\mathbf u\) of \(Q\) with height
\(O(M)\). Complete \(\mathbf u\) to a rational basis by two
standard vectors \(\mathbf v,\mathbf w\), and let \(B_Q\) be
the symmetric bilinear form for which \(Q(\mathbf v)=B_Q(\mathbf v,\mathbf v)\).
For \(\mathbf r=s\mathbf v+t\mathbf w\), define
\[
        (F_1,F_2,F_3)(s,t)
                = Q(\mathbf r)\mathbf u-2B_Q(\mathbf u,\mathbf r)\mathbf r .
\]
These are integral quadratic forms of coefficient height \(O(M^2)\).
They parametrize the conic by secants through the point of
\(\mathbf u\). They have no common geometric zero:
such a zero would give an isotropic two-dimensional space
spanned by \(\mathbf u,\mathbf r\) in a nondegenerate
three-dimensional quadratic space over the algebraic closure,
which is impossible. The parametrization is birational, with
inverse given off the base point by the line through \(\mathbf u\).
It is therefore an isomorphism of the projective line onto the
conic. The three coordinate lines meet the conic transversely
at distinct pairs of points, and no two such lines meet on
the conic. Hence each \(F_i\) has two simple geometric roots
and the pairs are disjoint.
[Cassels, Theorem on small zeros of an isotropic integral quadratic form](https://www.cambridge.org/core/journals/mathematical-proceedings-of-the-cambridge-philosophical-society/article/abs/bounds-for-the-least-solutions-of-homogeneous-quadratic-equations/9C3C57BC5549819DF3DBF5F3BC42DCE8)

## Numerical and finite-place bounds on the conic parameter

For primitive integer parameters \((s,t)\), put
\(H=\max(|s|,|t|)\) and \(g=\gcd(F_1(s,t),F_2(s,t),F_3(s,t))\).
The resultant identities in both affine charts show
\[
 g\le M^{O(1)},\qquad
 \max_i|F_i(s,t)|=H^2 M^{O(1)}
\]
where the second equality means upper and lower bounds with
fixed polynomial factors in \(M\). The formulas for the same
root choice are \(|F_1|=g x,|F_2|=g y,|F_3|=g z\), since the
original roots are primitive. There are two complete values
of height comparable to \(X\), so \(H=X^{1/4+O(f_*\lambda)+o_X(1)}\)
and all three coordinate logarithmic heights are
\[
       \log_X |F_i(s,t)|
         =\tfrac12\log_X T_i+O(f_*\lambda)+o_X(1).
 \tag{K2}
\]
Here and below the suppressed constants are absolute.

For an outer sequence of small-budget boxes, let \(h\)
be its limiting smallest complete height exponent.
One may partition the real parameter into finitely many
charts in which every cross determinant between two
primitive parameters in a chart satisfies
\[
        |s t'-s't|\le X^{h/2+O(f_*\lambda)+o(1)}.
 \tag{K3}
\]
The \(o(1)\) includes convergence of the box height exponents.
If \(h=1\), use the usual \(2H^2\). If \(h<1\), choose one
of two real affine charts in which the denominator has
absolute value \(H\). The short coordinate polynomial is
quadratic, or linear in an affine chart with a root at
infinity. Its distinct finite roots are separated by
\(M^{-O(1)}\), from the coefficient and discriminant bounds.
Equation (K2), divided by \(H^2\), forces the chart coordinate
to lie within \(X^{-(1-h)/2+O(f_*\lambda)+o(1)}\) of one
real root. Indeed a first root proximity follows by taking
the square root of the bound on the polynomial value;
after this is much smaller than the separation, only that
root is close, and the full bound applies after dividing by
the other factor. A nonreal root and its conjugate would
contradict the separation. Record that root in the chart,
and multiply the interval length by the two denominators
to prove (K3).
[Height and real correction for a projective carrier](ABC_RATIONAL_PROJECTIVE_CARRIERS.md#height-local-branch-weights-and-the-real-correction)

Fix a bad-prime product \(J\ne0\) containing \(2ABC\), the
leading coefficients and each nonzero homogeneous pairwise
resultant of the \(F_i\), and the derivative resultants of
their nonconstant dehomogenizations in both affine charts.
It is \(M^{O(1)}\) up to an absolute factor. If a prescribed
whole prime power \(q^e\Vert T_i\) is in a lock, work in a
\(q\)-adic projective chart whose denominator is a unit.
For \(q\nmid J\), the coordinate polynomial has a simple
root and its value has valuation \(e/2\). Recording the chart
and the branch therefore fixes the projective congruence
to depth \(e/2\).
For \(q\mid J\) the value valuation is at least
\((e-1)/2\) even accounting for the evaluated gcd. If its
valuation exceeds \(2v_q(J)\), the derivative resultant
bounds the derivative valuation by \(v_q(J)\) and Hensel
gives a root class with depth at least
\((e-1)/2-v_q(J)\). At lesser valuations discard the
precision at cost at most \(O(v_q(J)+1)\). The total lost
log weight across the selected bad support is \(O(\log J)\).
A chart and one of at most two roots per prime cost
\(X^{o_X(1)}\) combined records, by the divisor estimate.
Thus locks on any collection of letters of combined
numerical weight \(w\) give a divisibility of the cross
determinant of any two parameters in the same record by a
modulus of size
\[
                         X^{w/2-O(f_*\lambda)-o_X(1)} .
 \tag{K4}
\]
The selected prime sets on the three letters are disjoint
by primitivity. The weights are the valuations of the primitive quotient throughout.
[Exceptional support in projective residue records](ABC_RATIONAL_PROJECTIVE_CARRIERS.md#projective-residue-records-including-exceptional-primes)

## Using all actual endpoint locks

The affordable whole-coordinate pair-lock result gives an
absolute \(\eta_0>0\) usable at all balanced endpoint
sequences. It supplies full divisors of some two letters
of combined normalized radical cost at most \(1/3-\eta_0\)
and numerical weight at least the shorter height of that
pair plus \(\eta_0\). Its height threshold is at least the
global shortest height \(h\), and it permits macroscopic
individual primes. Use these locks on a hypothetical
sequence of boxes near saturation before conditioning on
the square kernels. After conditioning, the old bound on
the possible lock values still holds, and the kernel
record adds at most \(f_*\lambda\) to the exponent.
[Affordable whole-coordinate pair locks at every balanced endpoint](ABC_WIDE_COFACTOR_CONGRUENCES.md#affordable-whole-coordinate-pair-locks-at-every-balanced-endpoint)

Choose \(f_*<\eta_0/4\). If the claimed bound failed for
arbitrarily small fixed budgets and arbitrarily small absolute
gaps, form a sequence tending to the balanced \(1/3\)
endpoint as in the existing small-budget box reduction.
The total radical exponent \(L\) satisfies \(L/\lambda\to1\).
The available divisor-record count, including all kernel
choices, has exponent at most
\((1/3-\eta_0+f_*+o(1))\lambda\).
For such a record (K4) beats (K3) by an absolute numerical
margin (multiplied by \(1/2\)), since \(f_*\lambda\to0\).
It fixes the primitive parameter, and hence the actual
triple. This contradicts saturation and yields the absolute
\(\delta,\lambda_*\). Every kernel, resultant, branch chart,
and finite truncation has been paid at its fixed budget
before the outer limit.

In the balanced minimally repaired sector with
\(\operatorname{rad}(c)\ge c^{(1/3-\xi)\lambda}\) and repair
height \(H_N(c)\le Y\) for a fixed order \(N\), the height
cutoff implies \(c\le Y^{1/(N(1/3-\xi)\lambda)}\).
Applying the original count gives repaired exponent at most
\((1/3-\delta)/(N(1/3-\xi))+\varepsilon\), which is below
\(1/N\) for \(\xi<\delta\). This formulation takes the
original total radical budget to be \(\lambda\), with fixed
small enlargements absorbed when using the alternative
central-band conventions.
[One-target minimal repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## A quantitative charged-kernel selector

The fixed-record proof also yields a direct sector at
non-negligible kernel mass. Let \(L=\log_X\operatorname{rad}(abc)\).
Define the normalized kernel mass \(f=\log_X(ABC)/L\) and the
numeric weight
\[
          W=\sum_{q\mid ABC} v_q(abc)\log_X q .
\]
The sum uses entire odd-valuation prime powers, not just their
squarefree parts. Once \(A,B,C\) are paid, all exact prime
powers on their support have only subpower many possibilities,
by the fixed-support divisor estimate. Include all of them
in (K4), gaining weight \(W/2\) with no additional radical cost.
If \(W>h\) by a fixed margin, this already proves a count
exponent \(fL+o(L)\) whenever \(f<1/3\) by a fixed margin.

More generally, for the remaining even-valuation support,
a whole-prime selector of normalized mass \(s\) and numerical
weight \(w\) yields exponent \((f+s)L+o(L)\) whenever
\(W+w>h\) by a fixed margin. If its largest individual
normalized prime mass is at most \(\tau\), one can use the
explicit bound
\[
       f+\frac{(h-W)_+(1-f)}{2+h-W}+\tau
 \tag{K5}
\]
for the coefficient, up to arbitrary fixed allowances and
vanishing small-budget losses, in the range \(W\le h\).
Indeed the remaining support has mass \(1-f\) and numerical
weight \(2+h-W+o(1)\), since two complete letters are tall.
Sort its primes by their scaled valuation \(L v_q(abc)\).
To achieve a desired numerical weight \(h-W+\gamma\)
with an arbitrarily small fixed margin \(\gamma>0\), the
fractional top tail costs at most that target times
\((1-f)/(2+h-W+o(1))\). Round up the last prime to obtain
an actual full divisor, increasing the cost by at most
\(\tau\). Enumerating the selected product by its radical
cap costs exactly this normalized amount up to subpower
losses. Neither the order nor the subset need be common
across the triples. The proof applies in boxes with
uniform values of these parameters up to arbitrarily
small slack; dyadic shells and fixed parameter bins
handle their variation uniformly.

Thus with \(\tau=o(1)\) and \(W\le h\), the explicit
sufficient condition is \(W>h-1+3f\) by a fixed margin,
since this places (K5) strictly below \(1/3\).
In particular a positive ordinary short-height gap
permits a fixed positive kernel cap without any lower
bound on \(W\). With no height gap, a positive proportion
of odd radical at sufficiently high degrees can supply
the needed weight even when \(f\) is a non-negligible
fraction of the total budget. This calculation does
not guarantee a saving on every balanced parity profile.

## Two-square-kernel version and its rate obstruction

If only two letters \(U,V\) have fixed squarefree kernels,
write them as \(A x^2,B y^2\). The third letter \(T\) is
\(A x^2\pm B y^2\) up to an overall sign, so the primitive
parameter is simply \([x:y]\). It has height
\(X^{1/2+O((\log_X A)+(\log_X B))+o(1)}\).
There are now one contact of multiplicity two on
each kernel letter and two simple contacts on \(T\).
Consequently a full power divisor contributes half its
numerical weight on \(U,V\), and all its weight on \(T\).
The same finite-place precision proof applies with
coefficient error \(O(\log_X(AB))\). The cross determinant
threshold in ordinary numerical units is
\[
 \begin{cases}
   1, & \text{all three complete heights tend to one},\\
   (1+h)/2, & \text{one of \(U,V\) is short of height \(h<1\)},\\
   h, & \text{\(T\) is short of height \(h<1\)} .
 \end{cases}
\]
Indeed the real chart widths have respective exponents
zero, \(-(1-h)/2\) near the double zero, and \(-(1-h)\)
near a simple zero. Thus if a kernel letter has a
fixed positive short-height gap and the tall remaining
letter \(T\) has a prescribed positive-mass bounded
scaled-degree diffuse block, a strict saving follows
for sufficiently small joint kernel cost on **only**
\(U,V\). Start with the full lock on \(T\) of mass
\(1/3+o(1)\) and numerical weight \(1+o(1)\), then omit
a small positive block from that diffuse band: its
numerical loss can be made smaller than the fixed
determinant margin \((1-h)/2\). Pay the two kernel
radicals from the omitted radical margin. Equivalently
use actual bounded-degree-band primes small enough to
round the omitted block. If \(T\) is short by a fixed
margin below height \(1/2\), the analogous count starts
with one tall kernel letter of numerical branch
weight \(1/2\) and omits a diffuse block there if
available.

This specific two-kernel argument has a height-one
residual obstruction at the level of profiles. In the
balanced case with all heights tending to one and
each scaled degree supported on a narrow diffuse
interval about \(3\), all \(U,V\) half-rates are less
than every \(T\) full rate. The entire radical budget
\(1/3\) on \(T\) gives exactly the threshold weight one,
and replacing any of it by equally cheap locks from
\(U,V\) decreases weight. Thus no such threshold lock
gives a saving on this profile, even with trivial
kernels on \(U,V\). This does not assert that the
profile is realized by primitive triples, or rule out
a different use of the conic and the third radical cap.
