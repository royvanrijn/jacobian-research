# Two high power orders with an arbitrary third letter

Two exact power orders are enough for a small-budget saving
when the coefficient radical on those letters is small. The orders
need not have a common divisor. The third letter obeys a radical cap
and a lower numerical height cap, with arbitrary valuations.

## Count and range

Fix \(h\in(0,1]\), \(r>0\), \(0<t<r\) and \(\eta>0\).
There are constants \(A\) and \(\lambda_0>0\), depending only on
these parameters, with the following property. For
\(0<\lambda<\lambda_0\) fixed, count primitive positive triples
with \(c\le X\) which, after permuting the letters and fixing the
signs, have the form
\[
 s_1 x^P \ \pm\ s_2 y^Q=\pm T,\qquad
 \operatorname{rad}(s_1s_2)\le X^{f\lambda},\quad
 X^h\le T\le X,\quad
 \operatorname{rad}(T)\le X^{r\lambda},
\]
where \(P,Q\) are positive integers with
\[
                    P,Q\ge A,\qquad \lambda P Q\ge A .
\]
Then
\[
 \#\{(a,b,c)\}\ll_{h,r,t,\eta,\lambda,\varepsilon}
        X^{(f+\max(t,r-t)+\eta)\lambda+\varepsilon}.
 \tag{TG1}
\]
This is uniform in \(f\ge0\). The orders can vary over the triples
and with \(X\). Taking \(t=r/2\) makes the coefficient
\(f+r/2+\eta\). If every sufficiently high-order prime on \(T\) has
root mass at most \(t\lambda\), only the first alternative is used,
and \(f+t+\eta\) suffices. The exact condition on such primes is
specified in the proof.

In particular, the order hypotheses hold for all sufficiently small
\(\lambda\) if \(P\ge\lambda^{-\alpha}\) and
\(Q\ge\lambda^{-\beta}\) with \(\alpha,\beta>0\) and
\(\alpha+\beta>1\). They also admit, for example, a large bounded
order on one letter and an order of size \(1/\lambda\) on the other
with a suitable constant. The count always takes the height limit
at fixed \(\lambda\) first.

## A rectangular determinant over a full modulus

Let \(C_x,C_y\ge1\) be real and \(k\) an integer with \(k\ge B C_x C_y\), where \(B\) is
a sufficiently large absolute constant. Let \(m\ge2\) be \(k\)-full,
so \(p^e\Vert m\) implies \(e\ge k\). For fixed positive integers
\(P,Q\), count coprime positive pairs \((x,y)\) with
\[
 x\le m^{C_x/k},\qquad y\le m^{C_y/k},\qquad
 x^P y^{-Q}\equiv K\pmod m ,
\]
both unit modulo \(m\). The count is
\[
       O_{P,Q,k,C_x,C_y,\varepsilon}(m^{1/k+\varepsilon}),
 \tag{TG2}
\]
uniformly in the unit \(K\). The lower bound on \(k\) is independent
of \(P,Q\). The rectangle and its coefficient-height estimate are
important when the root height constants are allowed to diverge as
the budget decreases.

At each \(p^e\Vert m\) record both coordinates modulo
\(p^{b_p}\), where \(b_p=\lfloor e/k\rfloor\). At 2 use depth
\(\max(2,b_2)\). The condition modulo this shallow power gives at most
\(O_{P,Q}(p^{b_p})\) records. Fixing one coordinate leaves a power
congruence on the other, and its unit kernel has size at most
\(2\max(P,Q)\) by the prime-power unit groups. The exceptional factor
at 2 costs an absolute constant. Multiplying gives
\(O_{P,Q,k,\varepsilon}(m^{1/k+\varepsilon})\) charts, using
\(C^{\omega(m)}\ll_{C,\varepsilon}m^\varepsilon\).

On one chart choose a pair \((x_0,y_0)\) as reference at each prime.
Put \(v=\min(v_p(P),v_p(Q))\). For any other pair in the chart,
\((x/x_0)^{P/p^v}(y/y_0)^{-Q/p^v}\) is a principal unit,
in \(1+4\mathbb Z_2\) at 2, and its \(p^v\)-th power is
\(1\bmod p^e\). The principal-unit valuation formula gives the
relation modulo \(p^{e-v}\) before raising to \(p^v\).
One of \(P/p^v,Q/p^v\) is a \(p\)-adic unit, so one coordinate
is a binomial series in the relative increment of the other, modulo
\(p^{e-v}\). Its coefficients are \(p\)-adically integral, as
binomial coefficients with exponent in \(\mathbb Z_p\) are integral.
The depth of the free increment is at least \(b_p\).
Discard \(p^e\) with \(e\le 2v+4\) here. For \(k>4\) their product
is bounded in terms of \(P,Q\), allowing the other primes to retain
\(1-o(1)\) of the log weight of \(m\). The parameter choice and
orientation are allowed to vary across primes.

Take monomials \(x^i y^j\) with \(0\le i\le D_x\),
\(0\le j\le D_y\), where \(D_x=\lceil B_0 C_y\rceil\),
\(D_y=\lceil B_0 C_x\rceil\) for a sufficiently large absolute
\(B_0\). Write \(N=(D_x+1)(D_y+1)\), \(R=\lfloor N/2\rfloor\)
and \(H_0=D_x C_x+D_y C_y\). A minor of \(R\) columns and rows
has Archimedean size \(O(m^{R H_0/k})\), with constants fixed
before the height limit. Taylor expansion at a retained prime gives
divisibility by \(p^{b_p R(R-1)/2}\). To see this, truncate every
series at order \(R-1\). Distinct Taylor columns have distinct
orders, and an error column has order at least \(R\).
The same bound holds for errors in the congruence because
\(e-v\ge b_p R\), once \(B\) is large compared with \(B_0^2\).
Summing over primes and using \(b_p\ge e/(2k)\) gives
divisibility of log weight at least
\((1-o(1))R(R-1)\log m/(4k)\).
Choose \(B_0\) so that \((R-1)/4>2H_0\), possible since
\(R\gg B_0^2 C_x C_y\) and \(H_0\ll B_0 C_x C_y\).
All those minors vanish for sufficiently large \(m\).

The evaluation rank on a chart is thus less than \(N/2\).
Use at most \(N/2\) independent rows and apply the pigeonhole
principle to coefficient vectors with coordinates in
\([0,\lceil 4N m^{H_0/k}\rceil]\cap\mathbb Z\).
There are more vectors than their possible images under the selected
rows. A nonzero kernel polynomial follows with coefficient height
\(O_N(m^{H_0/k})\). This avoids the factor of the rank in the
height obtained from a single cofactor minor. The polynomial has
bounded degree at fixed \(C_x,C_y\), and each of its primitive
irreducible factors has coefficient height \(O_{D_x,D_y}(m^{H_0/k})\)
and bidegree at most \((D_x,D_y)\). For the coefficient estimate,
normalize the two factors to coefficient sup norm one in the fixed
finite-dimensional real coefficient spaces. The product has positive
minimum norm on the product of these compact spheres, since the
polynomial ring is a domain. Undoing the normalization gives a bound
on the product of the two norms by a fixed multiple of the original
polynomial's norm. A nonzero integral cofactor has norm at least one.

For an irreducible factor \(G\), set
\(J=QxG_x+PyG_y\). If \(G\mid J\), then \(x^P/y^Q\) is constant
on each geometric component of \(G=0\) in the torus. Indeed the
logarithmic differential of this quotient is dependent with \(dG\)
at smooth points, and characteristic zero identifies the elements
with zero differential as the constants on a geometric curve.
An exact quotient determines at most one coprime positive pair
because its reduced numerator and denominator are \(x^P,y^Q\).
If \(G\nmid J\), the points \(G=J=0\) are bounded in number by
Bézout. Monomial factors such as \(x\) or \(y\) give no unit pairs.

On the remaining points \(0<|J|\ll m^{2H_0/k}\) with constants
depending on the fixed parameters. Discard primes where
\(v_p(J)>e/5\). The lost log-weight is at most \(5\log |J|\),
which is less than \((1/10)\log m\) for large \(m\) if \(B\)
is sufficiently large. At other primes choose an integral
representative of \(K\) and apply the multivariate Hensel lemma
to \(G=0\) and \(F=x^P-Ky^Q=0\). On the original unit point
the Jacobian determinant is a unit multiple of \(J\)
modulo \(p^e\), so its valuation is at most \(e/5\).
A lifted solution is isolated and congruent to the point modulo
\(p^{\lfloor 4e/5\rfloor}\). There are a bounded number of such
isolated solutions per prime by Bézout, uniformly in coefficients.
Record the primes discarded and the lifted point at every other prime.
The number of choices is \(m^{o(1)}\) at fixed \(P,Q,k,C_x,C_y\).
Two original integer pairs in the same record agree modulo a product
of powers greater than \(m^{1/2}\), which is larger than either
coordinate bound for \(B\) large. Thus each record gives at most
one original pair. Summing the charts proves (TG2).

The fixed-order conditions here are local to the lemma. A
large-order separation will keep them legitimate in (TG1).

## How the third letter supplies the modulus or a high power

First choose \(s_1,s_2\), at a cost \(X^{f\lambda+\varepsilon}\).
Choose their joint radical and partition its prime support between
the two coprime coefficient letters. The number of partitions is
\(X^{o(1)}\) uniformly at polynomial height. The fixed-support
Rankin estimate counts the exact coefficients at \(X^{o(1)}\)
cost in the exponent sense. More precisely, for any fixed
\(\varepsilon>0\), there are at most
\(X^\varepsilon\prod_{p\mid R}(1-p^{-\varepsilon})^{-1}
=X^{\varepsilon+o(1)}\) integers up to \(X\) supported on \(R\),
uniformly for \(R\le X^2\). Small exponent losses are chosen
before the inner height limit. This also gives the claimed
uniformity in \(f\). Primitive solutions make \(x,y,T\)
pairwise coprime.

Choose a fixed \(\kappa\in(0,h/(2r))\) and put
\(k=\lfloor\kappa/\lambda\rfloor\). The primes of valuation at
least \(k\) on \(T\) carry numerical weight at least \(h/2\):
all other primes contribute at most \(kr\lambda<h/2\).
There are two cases.

If none of these high-order primes has \(\log_X p>t\lambda\),
partition the high-order primes into at most \(\lceil 2r/t\rceil+1\)
blocks, each of total prime log mass at most \(t\lambda\).
Take a maximal consecutive block in any order and continue.
Every pair of consecutive blocks has combined mass greater than the
block cap, which proves the number bound. Some block has numerical
weight at least a fixed \(w_0>0\), depending only on \(h,r,t\).
Use its full prime powers as a modulus \(m\mid T\).
There are at most \(X^{t\lambda+\varepsilon}\) choices of \(m\)
by the same fixed-support argument. The equation gives
\(x^P/y^Q\equiv K\bmod m\) with \(K\) a fixed unit depending on
the coefficients and signs.

First suppose \(P,Q\le D/\lambda\), with \(D\) a fixed parameter
specified below. At fixed \(\lambda\) these are finitely many
orders. Use a fixed positive integer amplification order \(\ell\).
Color primes in two classes and keep the pairs with \(x\) supported
on the first class and \(y\) on the second. Some coloring keeps
\(X^{-o(1)}\) of all pairs. This follows by averaging the uniform
coloring probability \(2^{-\omega(xy)}\), since all pairs are
coprime. On taking \(\ell\)-fold coordinatewise products, the
congruence has right side \(K^\ell\), and the two coordinates are
still coprime. The number of presentations of any product
is \(X^{o(1)}\) at fixed \(\ell\), by the divisor bound.
Apply (TG2) to the products with
\[
      C_x=\max(1,\ell k/(P w_0)),\qquad
      C_y=\max(1,\ell k/(Q w_0)).
\]
The size condition on \(k\) follows for large \(A\) and small
\(\lambda\), since
\[
 \frac{C_x C_y}{k}
 \le \frac1k+\frac{\ell}{P w_0}+\frac{\ell}{Q w_0}
            +\frac{\ell^2 k}{P Q w_0^2}.
\]
After taking \(\ell\)-th roots, the analytic cost is
\(m^{1/(\ell k)+\varepsilon}\le
X^{\lambda(2/(\ell\kappa))+\varepsilon}\).
Take \(\ell\) large in terms of \(\eta,\kappa\), and then choose \(A\)
and \(\lambda_0\). For fixed coefficients and moduli the pairs
determine the entire triple.

If \(P>D/\lambda\), enumerate \(x^P\) at cost at most
\(X^{\lambda/D+o(1)}\). For \(x=1\) there is just one value,
and otherwise \(P\le\log_2 X\), so all order choices here cost
a logarithm. If \(Q\le D/\lambda\), the root \(y\) is fixed in
divisor-many unit classes modulo \(m\) by its power congruence.
Each class contains at most one admissible \(y\) since
\(y\le X^{1/Q}\le X^{1/A}<m\). If both orders are above
\(D/\lambda\), enumerate both powered terms at cost
\(X^{2\lambda/D+o(1)}\). The same argument applies with the
orders interchanged. Take \(D\) large in terms of \(\eta\).
This completes the first case at exponent \((f+t+\eta)\lambda\).

In the other case some \(p\) has \(\log_X p>t\lambda\) and
valuation \(e\ge k\) on \(T\). Write \(T=s_3 p^e\) and
enumerate \(s_3\) at cost \(X^{(r-t)\lambda+\varepsilon}\)
by its radical. The order satisfies
\(\kappa/(2\lambda)\le e\le 1/(t\lambda)\), so there are
finitely many orders at fixed budget. When both \(P,Q\le D/\lambda\),
use Browning–Verzobio's coefficient-uniform Theorem 1.2 and
distinguish the smaller of \(P,Q\).
If \(P_+=\max(P,Q)\), \(P_-=\min(P,Q)\), the two ordinary
root boxes have exponents at most \(1/P_+\) and \(1/e\).
The theorem's \(\log_X W\) is at most \(1/\sqrt{P Q e}\).
Using \(P Q\ge A/\lambda\), \(P_+\ge\sqrt{A/\lambda}\) and
\(P_-\ge A\), its three exponents are bounded, up to their
arbitrarily small \(\varepsilon\)-losses, by terms of size
\[
 O_{\kappa}\!\left(\frac{\lambda}{\sqrt A}
            +\lambda^{3/2}\right).
\]
Indeed \(W\) has exponent \(O_\kappa(\lambda/\sqrt A)\);
the second increment is at most
\(2\max(1/P_+,1/e)/\sqrt e
=O_\kappa(\lambda/\sqrt A+\lambda^{3/2})\);
the third increment is at most
\(1/(P Q)+1/(e P_-)=O_\kappa(\lambda/A)\).
The coefficients \(s_1,s_2,s_3\) are nonzero, with signs
absorbed in the equation, and the root triple is primitive.
[Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

If one order exceeds \(D/\lambda\), enumerate its powered
term as in the first case. If the other order does not exceed
\(D/\lambda\), the remaining binomial plane curve with a nonzero
constant has degree at least \(e\) and is absolutely irreducible.
The proof uses a simple zero of the nonconstant side over
\(\mathbb C\), then Eisenstein over its rational function field.
The coefficient-uniform plane bound in the larger root box gives
exponent at most \(\max(1/Q,1/e)/e\) if the large order was
\(P\). This is \(O_\kappa(\lambda/A+\lambda^2)\).
If both orders exceed \(D/\lambda\), their powered terms determine
the third letter. This completes the second case at exponent
\((f+r-t+\eta)\lambda\), and hence proves (TG1).
[Browning–Verzobio, Lemma 2.2](https://arxiv.org/html/2608.24512v1#S2.SS1)

## Implication for the balanced endpoint

For \(r\) near \(1/3\), the saving is uniform whenever
\(f+r/2\) is bounded below \(1/3\). For example, on exceptional
triples with the third radical cap \(r=3/8\) and joint coefficient
radical cap \(f=1/8\) on the two power letters, take
\(t=3/16\) and \(\eta<1/48\). At heights comparable to \(X\),
the existing short-height saving handles the case where one
letter fails the lower cap \(X^h\), for an absolute \(h>0\).
Apply (TG1) in the other case and sum dyadic height boxes.
The resulting exponent is strictly below \(\lambda/3\) on
this two-order sector. The two divisibility assumptions are
properties of the actual integer valuations.
[Small numerical short height](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small)

For a wholly diffuse limiting \(T\) profile, every individual
high-order prime has vanishing normalized mass. Indeed a
nonvanishing prime mass forces its scaled degree to remain
bounded by the numerical height, and then supplies an atom.
On such sequences fix any small \(t>0\) in the first case.
The cost approaches \(f\lambda\), with all margins fixed before
the height limit. This allows substantially more coefficient
radical on the two power letters if the third is diffuse.

The result extends the fine-grid sector to independent orders on
only two letters. It also includes instances in which the third
letter prevents three-letter use of a small residual coefficient
radical. A broad weak degree profile can coexist with either the
two divisibility conditions or their complete failure, because
a grid of large integer orders can have scaled mesh tending to
zero. Thus an endpoint family with no such pair remains outside
this count.
[Three-letter fine grids](ABC_FINE_VALUATION_DIVISIBILITY.md)
[Fixed-height omitted-root congruences](ABC_OMITTED_COORDINATE_GRAPH.md#signed-monomial-congruences-over-full-moduli)


## A restriction on actual divisibility masses at saturation

In a balanced limiting sequence of boxes let \(D_A,D_B\) be
the limiting normalized root masses on two letters of prime
valuations divisible by their respective orders \(P,Q\).
Here the orders may depend on the box. Assume
\(\min(P,Q)\to\infty\) and \(L P Q\to\infty\), taking the usual
fixed-\(L\) height limit first. Entire prime powers at the other
primes can be assigned to the two coefficients. Their joint
normalized radical mass is \(2/3-D_A-D_B+o(1)\).

On the third letter \(T\) let \(\gamma\in[0,1/3]\) be the
limiting greatest normalized prime log mass among primes whose
valuations are at least \(\lfloor\kappa/L\rfloor\). Here
\(\kappa\) is any fixed positive constant small enough for the
lower numerical height condition. First refine the bins for the
greatest prime mass to fixed small widths, then pass to a
subsequence and let the widths tend to zero. This preserves the
limiting exponent of the selected boxes. Use (TG1) with its budget
scale \(L\) and third cap arbitrarily close to \(1/3\).
The first case with \(t>\gamma\) costs \(t+o(1)\) in that
scale, and the second case with \(t<\gamma\) costs
\(r-t+o(1)\). Let \(t\) approach \(\gamma\) from the appropriate
side after fixing each proof margin. One obtains
\[
               D_A+D_B
       \ \le\ \frac13+\min\!\left(\gamma,\frac13-\gamma\right).
 \tag{TG3}
\]
When \(\gamma=0\) or \(1/3\), the one applicable inequality
gives the result by the same limiting argument. In particular
\(D_A+D_B\le1/2\) always.

If \(T\) has any positive nonatomic scaled-degree mass, then
\[
                              D_A+D_B\le1/3.                     \tag{TG4}
\]
Indeed choose a compact positive band carrying positive nonatomic
mass. It supplies whole root positions of arbitrarily small fixed
positive total mass \(t\), all of scaled degree at least a fixed
\(\kappa'>0\). Fix this block, of numerical weight at least
\(\kappa't+o(1)\), before using the first case's determinant and
amplification. Other primes on \(T\) need satisfy no size condition.
The proof of (TG1) conditional on a cheap full modulus applies
verbatim. Choose the desired count margin and \(t\) before the
outer small-budget limit. The supported-coordinate selection is the
same as in the affordable pair lock argument. This constraint uses
actual prime valuations at fixed budgets. Replacing them by the
weak limits of the scaled exponents loses the divisibility information.
[Selecting a small full block in a diffuse band](ABC_WIDE_COFACTOR_CONGRUENCES.md#small-full-blocks-in-diffuse-bands)

For the general \(abc\) problem, first use the one-letter
quarter inequality to confine an unsaved sequence to balancing
radical caps, and use the short-term bound to keep all numerical
heights positive. Then (TG3) gives a necessary condition for every
pair of orders as above and every choice of its complementary
letter. It leaves the case of spectra that have small exact
divisibility masses on every such pair.
[One-letter quarter bound](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-quarter-inequality-for-a-letter-of-small-radical-mass)


## A saving with a joint coefficient radical below five twelfths

Fix \(\zeta>0\). There exist \(A,\delta,\lambda_0>0\), depending
only on \(\zeta\), such that for \(0<\lambda<\lambda_0\) the
primitive exceptional triples having exact power decompositions on
all three letters with
\[
 \min_i P_i\ge A,\qquad
 \lambda \min_{i\ne j} P_iP_j\ge A,\qquad
 \operatorname{rad}(s_1s_2s_3)\le c^{(5/12-\zeta)\lambda}
\]
satisfy the bound \(O_{\lambda,\varepsilon,\zeta}
(X^{(1/3-\delta)\lambda+\varepsilon})\). Here the allowable
orders and factorizations may vary with the triple.
This includes symmetric orders just above a large fixed
multiple of \(\lambda^{-1/2}\).

To prove this, suppose there are actual boxes tending to
the one-third endpoint despite the displayed coefficient radical
cap, taking \(A\to\infty\). Such a diagonal sequence would follow
from a failure of the uniform assertion by choosing \(A=j\),
\(\lambda_0=1/j\) and a prospective saving tending to zero.
Use the anatomic reduction and dyadic radical intervals for the
three coefficients, whose count at each fixed budget is subpower.
The limiting normalized masses are measured in the actual box
budget \(L\). The repository's unrestricted one-third limsup
bound forces \(L/\lambda\to1\) on such a sequence, and the
one-letter quarter bound gives letter masses \(1/3+o(1)\).
The short-term bound keeps their numerical heights bounded away
from zero. Thus the hypotheses of (TG3)–(TG4) hold for each
pair. Write \(f_i\) for the limiting normalized coefficient
radical masses and \(f=\sum_i f_i\le5/12-\zeta\). One may
use the canonical coefficients supported exactly on primes
whose valuation is not divisible by their letter's order,
because this can only decrease their radicals. Refining the
coefficient radical ranges gives a limiting subsequence even
when the admissible orders vary with the height.

If all three letters have positive nonatomic mass, (TG4) gives
\(f_i+f_j\ge1/3\) for all pairs, so \(f\ge1/2\).
Otherwise at least two letters have positive nonatomic mass
by the pairwise endpoint diffusion bound. Let \(z\) be the
coefficient radical mass of the remaining letter. The two
inequalities (TG4) give \(f+z\ge2/3\). Inequality (TG3) with the
remaining letter as \(T\) gives \(f-z\ge1/6\), since
\(D_A+D_B=2/3-(f-z)\le1/2\).
Adding yields \(f\ge5/12\), a contradiction.
The fixed absolute savings outside the balanced positive-height
sector and this diagonal contradiction prove the assertion.
[The one-third limiting count and endpoint conditions](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
· [Pairwise endpoint diffusion](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)


## A thin-band modulus raises the joint coefficient cap to four ninths

Fix \(\alpha>2/3\). If the three orders satisfy
\(P_i\ge\lambda^{-\alpha}\), the preceding uniform conclusion
holds with joint coefficient radical cap
\(c^{(4/9-\zeta)\lambda}\) in place of
\(c^{(5/12-\zeta)\lambda}\), for every fixed \(\zeta>0\).
The constants may now also depend on \(\alpha\). The new
input is a modulus of possibly vanishing numerical weight on
the nondiffuse letter.

Work in the endpoint boxes at budget \(L\). By
\(L/\lambda\to1\), after reducing \(\alpha\) slightly we can use
order lower bounds \(L^{-\alpha}\), with \(2/3<\alpha<1\).
Suppose \(T\) has a normalized root mass \(g_T>0\) on valuations
divisible by \(P_T\). All such valuations are at least
\(L^{-\alpha}\). At a fixed large \(B_1\), those above \(B_1/L\)
have total normalized mass at most \(1/B_1+o(1)\) by the height
cap. Choose \(B_1\) large in terms of \(g_T\) and partition
the remaining indices into dyadic bands \([E,2E)\), with integer
lower bounds comparable to at least \(L^{-\alpha}\).
Some band has normalized mass \(\mu\gg_{g_T}1/\log(1/L)\).
Refine by the selected band and a dyadic interval for its root
mass, so that all remaining triples have mass between
\(\mu\) and \(2\mu\). At each fixed \(L\) there are finitely
many such refinements, and the constants in the mass lower
bound can be taken uniform for \(g_T\) in a fixed positive range.

Freeze the product \(m\) of all selected full powers in that band.
Its root cost is at most \(g_T+o(1)\) as \(g_T\) is the
total normalized divisible mass. For its actual numerical weight
use bounds \(w_-\ge E L\mu\) and \(w_+\le4E L\mu\).
It is \(k\)-full with \(k=E\). For fixed coefficients on the other
two letters and fixed orders there, their power-root bounds
are \(X^{\xi_i}\) with \(\xi_i\le L^\alpha\).
After fixed amplification \(\ell\), put
\(C_i=\max(1,\ell k\xi_i/w_-)\). The condition of (TG2) holds
since
\[
 \frac{C_1 C_2}{k}
  \ \ll_{\ell,g_T}\
 L^\alpha+L^{2\alpha-1}\log(1/L)
              +L^{3\alpha-2}\log^2(1/L)=o(1).
\]
The chart exponent after amplification is at most
\(w_+/(\ell k)\le 4\mu L/\ell\). Taking \(\ell\) large
reduces it to an arbitrarily small fixed normalized allowance.
The same finite-order cutoff at \(D/L\) used in (TG1) is valid.
When one very high order is enumerated, the other root lies
in subpower many unit classes modulo \(m\), with one root per
class since \(w_-/L^\alpha\to\infty\). Thus this method gives
the limiting inequality
\[
                            d\le f_A+f_B+g_T.                   \tag{TG5}
\]
The bound pays the actual normalized radical mass of the band.
It does not use a uniform positive lower bound for its numerical
height across budgets.

In a hypothetical saturated sequence with \(f<4/9\), all
three letters cannot have positive nonatomic mass, since that
would give \(f\ge1/2\) by (TG4). Let \(z=f_T\) on the remaining
letter. The other two letters have positive nonatomic mass,
so \(f+z\ge2/3\) and \(f-z\ge1/6\) as before. In particular
\(z\le f-1/6<1/3\), and \(g_T=1/3-z>0\).
Now (TG5) gives \(1/3\le 1/3+f-2z\), so \(z\le f/2\).
Together with \(f+z\ge2/3\), this forces \(f\ge4/9\).
For a fixed margin \(\zeta\) below \(4/9\) the same
diagonal contradiction excludes every actual endpoint,
giving a uniform saving on the stated sector.

An exact single valuation grid is unnecessary if the
positions on the two free letters cluster in a bounded or
sufficiently slowly growing list of fine phases. A signed
two-basis representation controls the product of the two
rational root heights, and amplification counts distinct
exact outputs even when their two input bases have common
prime support.
[Signed valuation phases](ABC_SIGNED_VALUATION_PHASES.md)
