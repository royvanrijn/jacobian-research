# Compressing sparse signed valuation grids to two bases

A fixed number of unrelated valuation grids can be combined if
the product of their coefficient widths is smaller than the
range of actual integer valuations. The compression uses a
short kernel relation chosen in the intersection of the coefficient
box with the actual valuation strip. It can therefore handle
arbitrarily large signed grid steps whose contributions cancel.

## Rank reduction inside a bounded valuation strip

Fix \(d\ge2\). Let \(J\ge2\), \(R_i\ge1\), \(P_i\in\mathbb Z\),
and \(\mathcal V=\prod_{i=1}^d R_i\).
Given any integer set
\[
 E\subseteq[-J,J]\cap
 \left\{\sum_{i=1}^d P_i u_i:
          u_i\in\mathbb Z,\ |u_i|\le R_i\right\},
\]
there are two integer steps \(P,Q\) and two widths \(U,V\ge1\)
such that \(E\subseteq\{P u+Qv:\ |u|\le U,\ |v|\le V\}\) and,
provided \(\mathcal V/J\) is sufficiently small in terms of \(d\),
\[
              UV\ll_d J(\mathcal V/J)^{1/(d-1)}.             \tag{VG1}
\]
The implicit constant is independent of the steps \(P_i\).
An unused step may be taken nonzero with its digit identically zero.

To see this, set \(\phi(\mathbf u)=\sum_i P_i u_i\) and
intersect the coefficient box with \(|\phi|\le J\).
Call the resulting symmetric convex body \(K\).
Restrict to the real span of \(K\cap\mathbb Z^d\),
with its full induced lattice \(\Gamma\).
The case of zero span or \(\phi=0\) on this span is trivial.
In dimension \(r>0\) the relative volume
\[
                        B=\operatorname{vol}_r(K)/\det\Gamma
\]
is \(O_d(\mathcal V)\). Indeed an independent set
of \(r\) original coordinate projections identifies this span
with \(\mathbb R^r\), sends \(\Gamma\) into \(\mathbb Z^r\), and
bounds the projected volume by \(\prod (2R_i)\) over those
coordinates. Relative volume is invariant under the linear
identification. The lattice points of \(K\) span the space by
construction.

If \(r\ge3\), let \(\Gamma_0=\Gamma\cap\ker\phi\), and
let \(g\ge1\) generate the integer group \(\phi(\Gamma)\).
In the ambient Euclidean metric,
\(\det\Gamma_0=\det\Gamma\,\|\phi\|/g\), since the spacing
of the levels of \(\Gamma\) normal to \(\ker\phi\) is
\(g/\|\phi\|\). By the central-section inequality for
symmetric convex bodies (Brunn–Minkowski), the relative
volume of the middle section in its kernel lattice is at least
\[
  \frac{\operatorname{vol}_{r-1}(K\cap\ker\phi)}{\det\Gamma_0}
                 \ \ge\ \frac{gB}{2J}.
\]
Minkowski's convex body theorem therefore gives a nonzero
kernel vector in \(sK\cap\Gamma\) with
\(s\ll_r(J/B)^{1/(r-1)}\). Divide by its lattice content to
make it primitive, and denote the resulting vector by \(z\).

Project orthogonally along \(z\), obtaining \(\overline K\) and
the quotient lattice \(\overline\Gamma\). This is a lattice of rank
\(r-1\) with \(\det\overline\Gamma=\det\Gamma/|z|\);
for these facts complete \(z\) to a basis of \(\Gamma\).
The functional \(\phi\) descends to an integer-valued
functional on \(\overline\Gamma\) and still has absolute value
at most \(J\) on \(\overline K\). Every point in the original
coefficient sample projects to a lattice point of
\(\overline K\) with the same output value. The projected
lattice points still span. For each point \(y\in
\frac12\overline K\), choose a lift \(x\in\frac12 K\).
The segment \([x-z/(2s),x+z/(2s)]\) is contained in \(K\) by
convexity and has length \(|z|/s\). Integrating the chord bound
gives
\[
  \frac{\operatorname{vol}_{r-1}(\overline K)}
            {\det\overline\Gamma}
     \ll_r sB
     \ll_r J(B/J)^{(r-2)/(r-1)}.
\]
These inequalities apply throughout the iteration when
\(\mathcal V/J\) is sufficiently small. The exponent products
telescope to \(1/(r-1)\) from rank \(r\) through rank two.

At the final rank at most two, the number of lattice points of
the current convex body is \(O_d\) times its relative volume.
Indeed choose independent vectors in the body spanning its
lattice over \(\mathbb R\). Their generated sublattice has a
fundamental parallelepiped contained in a fixed dilation of
the body. A measurable fundamental domain for the full
lattice can be chosen as a subset of that parallelepiped by
selecting one representative of each full-lattice coset.
The translates of this domain by the body's lattice points
are disjoint almost everywhere and lie in a further fixed
dilation. Compare volumes. Tao–Vu's discrete John's theorem
then covers the lattice sample by a symmetric progression
of rank at most two, whose coefficient-width product is at
most a constant times the number of lattice points in the
body. Project its generators by the descended functional
\(\phi\). This proves (VG1), padding by an unused direction
if the rank is below two and using \(r\le d\).
[Tao–Vu, Theorem 1.6](https://arxiv.org/pdf/math/0701005#page=3)

## A two-letter release inequality with several unrelated grids

Use the \(abc\) box scale \(n=X^L\). Suppose two letters \(U_0,V_0\)
have retained prime-resolved factors of total normalized radical
mass \(t\), and freeze the rest on those two letters at total
normalized radical cost \(f\). Assume the third letter supplies
\(n^{\mu+o_X(1)}\) records of actual divisors \(m\), with
\(m\ge X^{w_0}\) and \(m\) at least \(\lfloor\kappa/L\rfloor\)-full,
for fixed \(w_0,\kappa>0\). Suppose also that the retained integer
valuations are at most \(J\le B_0/L\) and lie in a finite union of
translated signed progressions of total rank bounded independently
of \(L\). For these progressions let \(\mathcal V(L)\) be the product
of all their coefficient widths (omit any width below one).
Then if \(L\mathcal V(L)=o(1)\), for any fixed \(\eta>0\) the count is
\[
           O_{L,\eta,\varepsilon}
                (X^{(f+\mu+\eta)L+\varepsilon})              \tag{VG2}
\]
for all sufficiently small fixed budgets in this family.

Indeed concatenate the bases of the progressions, adding the
translate of each as a basis of width one, and use only the block
that contains the valuation concerned. For applying (VG1) one can
take \(J=B_0/L\) regardless of the lower actual cap. It gives a
common two-basis signed representation, with digit widths
\(U,V\ge1\) and \(UV=o(1/L)\). Write the ratio of the two
retained factors in these bases as \(z^P w^Q\). The two rational
coordinate height exponents satisfy
\(\xi_z\le tLU,\ \xi_w\le tLV\), so each is \(o(1)\) and their
product is \(o(L)\). Signs and overlaps of the basis coordinates
are allowed. The frozen factors and the actual full modulus
give the congruence on the exact signed monomial. Apply the
two-basis amplified rational bound (SP3), whose hypotheses are
just these height products and the fixed record parameters.
The exact ratio recovers the two complete coprime labels after
multiplication by their frozen ratio, and then recovers the
third label. All inner-height exponents are taken at fixed
\(L\).
[Amplified two-basis counting bound](ABC_SIGNED_VALUATION_PHASES.md#from-shifted-phase-concentration-to-small-height-product)

The same conclusion holds if the progressions describe only a
weighted signed digit assignment \(e=\sum_i P_i I_i\) and
the coordinate absolute-log budgets \(\xi_i\) satisfy
\(\prod_i(\xi_i+L)=o(L^{d-1})\). To use a fixed retained-mass
margin \(\zeta>0\), remove at most \(\zeta\) normalized mass
by Markov limiting each \(|I_i|\) to \(O_{d,\zeta}(1+\xi_i/L)\),
and remove valuations larger than \(B_0/L\) at further
cost at most \(2/B_0\) by the two numerical height moments.
Run (VG2) with the actual post-removal freeze cost and choose
\(B_0\) large in terms of the margin. This variant assumes
the given \(\xi_i\) bound the absolute-prime sums, or the
reduced numerator and denominator heights of the associated
basis coordinates: for primitive ratios those latter heights
bound half the absolute-prime sums. Any chosen clipped
integer coefficients may be replaced, if desired, by
bounded-in-\(L\) steps representing the same clipped
degree set; the step system is a finite list of integer
linear equations in a fixed number of variables with
bounded integer coefficients and right sides at each
fixed \(L\). Each solvable system has some fixed
integer solution that can be chosen before the height limit.

In a balanced endpoint sequence with positive nonatomic mass
on the complementary letter, a compact positive band on that
letter supplies a full modulus of any sufficiently small fixed
positive radical cost with a fixed positive lower numerical
weight. If the two free letters retain \(t>1/3\) by a fixed
margin, their freeze cost is \(f=2/3-t+o(1)\).
Choose \(\mu,\eta\) and, in the weighted form, the removals
small enough to leave a fixed gap below \(1/3\) in (VG2).
This proves a saving on this sector of actual valuation
profiles. On the one-target minimal repair formulation,
the target may be the same complementary letter and the
same full divisor is present before repair.
[Selection of an affordable full modulus](ABC_WIDE_COFACTOR_CONGRUENCES.md#small-full-blocks-in-diffuse-bands)

The same release also holds at a sufficiently small fixed
\(L\mathcal V(L)\le\tau\) rather than requiring it to tend to
zero. By (VG1),
\(\xi_z,\xi_w\ll \tau^{1/(d-1)}\) and
\(\xi_z\xi_w\ll L\tau^{1/(d-1)}\), with constants allowed
to depend on the fixed retained-mass cap and on \(B_0,d\).
In (SP3) first fix its amplification order \(\ell\) for
the desired chart allowance \(\eta\). Then choose \(\tau\)
so small that both the scaled individual-height contribution
\(\ell(\xi_z+\xi_w)/w_0\) and the product contribution
\(\ell^2 k\xi_z\xi_w/w_0^2\) meet the rational determinant
threshold. This is a fixed saving with all margins chosen
before the inner height limit.

## Mismatched shifted grids and small additive doubling

For instance, each free letter may retain a positive-mass
range \(a_i+h_i[1,M_i]\) of original valuation positions,
with arbitrary and distinct offsets \(a_i\) and steps \(h_i\),
provided \(LM_1M_2=o(1)\) and the valuation positions are
\(O(1/L)\). No common grid or common phase alignment is
assumed. To give a diffuse schematic example, let
\(J=1/L\), choose \(M_i\asymp J^{0.48}\), and take distinct
\(h_i\asymp J^{0.52}\) with \(\gcd(a_i,h_i)=1\) and
\(a_i\asymp J\). Put equal radical weights on the retained
positions, for each free letter separately. The scaled
profiles can converge to diffuse interval laws. For
every candidate growing exact power order \(D\), the
retained fraction of positions divisible by \(D\) is at most
\(1/D+O(1/M_i)\): if \(\gcd(D,h_i)>1\) the congruence has
no solutions, and otherwise it selects one residue class
of the index. Moreover no family of \(O(\log J)\) degree
windows of radius \(J^\beta\), \(\beta<1/2\), can retain a
positive fraction because the steps exceed those widths
by a power. Nevertheless \(LM_1M_2\to0\), so the two
shifted grids satisfy the new release criterion. This is
a profile illustration rather than an existence claim.

Another consequence uses degree sets \(E_A,E_B\) with bounded
additive doubling at the actual valuation scale. If their
retained root masses add to more than \(1/3\) at the
endpoint, their degrees lie in \(O(1/L)\), and
\(L|E_A||E_B|\to0\), then the sector is saved whenever
the complementary letter admits the cheap full modulus.
Indeed apply the integer Freiman theorem separately to the
degree sets, obtaining translated progressions of bounded
rank and sizes \(O_K(|E_i|)\) when the doubling constants are
at most fixed \(K\). This gives (VG2) by concatenation.
If instead the union degree set has bounded doubling, apply
the theorem once and require only \(L|E_A\cup E_B|\to0\).
[Tao–Vu, statement of the Freiman theorem](https://arxiv.org/pdf/math/0701005#page=7)

The compression hypothesis concerns actual discrete valuation
positions. A generic dense diffuse profile spread across
\(\asymp1/L\) degree positions cannot fit these width bounds
at positive retained mass. In particular this proof does
not remove the full balanced endpoint.


## Sparse diffuse spectra outside all small-width grids

For any fixed \(0<\sigma<1\), there are symmetric sets of degrees
\(E_J\subseteq[2J,4J]\cap\mathbb Z\) of size
\(J^{\sigma+o(1)}\), with uniform weak empirical distribution
after dividing degrees by \(J\), such that for every fixed
rank \(d\), any signed progression \(P\) of total coefficient
width \(\mathcal V\le J\) satisfies
\[
             |E_J\cap P| \ll_d J^\sigma \mathcal V/J+\log J .
 \tag{VG3}
\]
Equivalently one may require this for each prescribed finite
set of ranks and select a single sequence by a diagonal
choice. In particular every rank fixed in the small-budget
limit and every width product \(o(J)\) captures vanishing
normalized mass.

Here is the uniform family count needed for the probabilistic
construction. Replace each coefficient width by its integer
part, and pad by widths at least one, so the widths are at most
\(J\). For a progression intersecting \([2J,4J]\), choose
one digit vector in its box for each covered output degree.
Choose a maximal linearly independent set of the corresponding
equations for the steps. Put \(A\) for a nonsingular minor in
the constrained columns, with
\(1\le|\det A|\le d!J^d\). Reduce the free integer steps
modulo \(|\det A|\), and recover the other steps by Cramer's
rule. The change of free steps by multiples of \(\det A\)
changes the constrained steps by integers. The new steps are
all \(J^{O_d(1)}\) and satisfy the selected equations, hence
all the original equations by consistency and maximal rank.
Thus the signed progressions with bounded steps,
of which there are at most \(J^{O_d(1)}\), suffice for upper
bounds on intersections; replacing steps can add covered
outputs but may not enlarge the coefficient widths.

Choose each pair \(\{3J-i,3J+i\}\) for \(1\le i\le J\)
independently with probability \(p=J^{\sigma-1}\), and put
both degrees of every chosen pair into \(E_J\).
For a fixed progression with at most \(3^d\mathcal V\)
outputs, the intersection count is a sum of independent
Bernoulli variables of parameter \(p\) with weights \(0,1,2\).
Its exponential moment at \(1\) is at most
\(\exp((e^2-1)p\,3^d\mathcal V)\).
Markov at a sufficiently large multiple of
\(p\mathcal V+\log J\) and a union bound over the bounded-step
family prove (VG3) with probability tending to one.
Size concentration and weak empirical uniformity follow
by the same independent bounded summands, first on a
finite partition of the fixed scaled interval and then on
arbitrarily fine fixed partitions. These events are
compatible with (VG3). Symmetry ensures that the
empirical mean is precisely \(3J\). With \(L=1/J\) and
normalized letter mass \(1/3\), this gives a schematic
diffuse profile of numerical height precisely one.
It describes possible boxed mass data, and asserts no
actual solutions of the additive equation.

A translate can be included by adding its offset as a single
width-one basis. A growing exact divisor of the original
valuations gives a
width-\(O(J/k)\) progression of rank one, so these models
also have vanishing mass on every growing-order divisibility
grid. The same example rules out the smaller digit-height
product condition in a fixed number of bases on a retained
mass bounded below: applying Markov to the bounded total
digit heights would produce a width-\(o(J)\) shifted grid
containing a fixed positive fraction of the sample, contrary
to (VG3).
