# Exclusive-boundary release for nondegenerate power collisions

Consider a fixed number \(s\ge 3\) of nonzero weighted
monomials in \(m\) positive integer coordinates. Write
\[
 {\cal F}(\mathbf x)=\sum_{\nu=1}^s
       A_\nu\prod_{j=1}^m x_j^{e_{\nu j}},
       \qquad A_\nu\in\mathbb Z\setminus\{0\},
       \qquad e_{\nu j}\ge0 .
\]
Let \(E\ge1\). Assume for every \(j\) that
\(\max_\nu e_{\nu j}-\min_\nu e_{\nu j}\ge E\).
There is a bound for its nondegenerate, pairwise coprime
positive solutions \(x_j\le H\):
\[
 \#\{\mathbf x:{\cal F}(\mathbf x)=0,\ 
                   \text{no proper vanishing subsum}\}
 \ll_{m,s,\mathbf e,\varepsilon}
 H^{\,C_{m,s} E^{-1/(2(m-1))}+\varepsilon}
                        \quad (m\ge 2).                    \tag{EB1}
\]
The implied constant is uniform in all coefficients
\(A_\nu\). For \(m=1\) the count is bounded in terms of
the exponents. In particular a variable can appear
to the same high order in several of the monomials,
provided it is absent from at least one.
The geometric criterion improves the exponent in the
original bounded-root release (30) and lets it act
on shared high factors in four-term collisions.
[Bounded-root release (30)](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#releasing-a-bounded-number-of-high-degree-coordinates)

## Small-degree components and exclusive affine boundaries

If the entire polynomial is identically zero, group
monomials of equal exponent vector. The spread condition
guarantees at least two groups, and each group is then a
vanishing proper subsum. Thus the nondegenerate sample
is empty. Otherwise its ambient algebraic closure has
dimension at most \(m-1\), with degree bounded in
terms of the fixed exponents. Cover it by geometrically
integral components, and use dimension descent.

At any branch of bounded geometric degree \(D\), we may
replace the component by a locus defined by integral
polynomials of bounded degree and height \(H^{O_{m,D}(1)}\).
Take the rational nullspace of evaluation on the
actual integer sample of all affine polynomials of
degree at most \(D\). Maximal minors give a basis of
that height. A geometrically irreducible variety of degree \(D\)
is cut out set-theoretically by forms of degree at most \(D\).
For each exterior point choose a general linear space through
it disjoint from the variety, of projective dimension
\(m-\dim Z-1\). Project to \(\mathbb P^{\dim Z+1}\) from a
hyperplane of that space omitting the point. The projected
variety is a hypersurface of degree at most \(D\) and misses
the projected exterior point. Pulling back its equation
separates that point. Such forms, even when defined over an
extension, belong to the scalar extension of this
nullspace. Hence the new locus is contained in the
old component and still contains the sample; if its
dimension has not dropped then it is the very same
component. In particular a component not defined over
\(\mathbb Q\) can only give a dimension drop in that step.
All degree bounds remain uniform in the \(A_\nu\).
[Narrow-locus preparation](ABC_WIDE_COFACTOR_CONGRUENCES.md#preparing-the-short-coordinate-locus)

Let \(Z\) be a positive-dimensional geometric integral
component containing nondegenerate points, and
suppose that its affine boundary \(x_j=0\) meets an
open locus where all the other coordinates are
nonzero. Since the sample has \(x_j>0\), this boundary
has a divisor component on \(Z\) meeting that open
set. Take a general linear curve section of the
projective closure of \(Z\). It is integral, of
degree \(D\), and meets this boundary divisor in
the indicated affine open set. It is not contained in
any proper subsum locus, since none contains \(Z\).
On its normalization, at a point over that boundary
only \(x_j\) has a zero, while all other coordinates
are regular units. Choose two terms whose exponents
in \(x_j\) differ by at least \(E\); their ratio has a
valuation of absolute value at least \(E\) there.
On the other hand, divide the full relation by one
of these terms and apply Brownawell–Masser to the
remaining unit ratios. The genus is at most
\((D-1)(D-2)/2\), and the set of coordinate zeros
and poles has size at most \(2mD\). Thus the degree
of any such ratio as a rational function is
\(O_{m,s}(D^2)\). Consequently \(D\gg_{m,s}\sqrt E\).
[Brownawell–Masser, Vanishing sums in function fields](https://doi.org/10.1017/S0305004100066184)
· [Function-degree statement, Fuchs–Zannier Lemma 2.2](https://ems.press/content/serial-article-files/31810)

If \(D<c_{m,s}\sqrt E\) for a sufficiently small
constant, the exclusive affine boundary is absent
for *every* variable. Once the defining equations
\(G_r\) on the sample have the height cap above,
the polynomials
\[
          (G_r)_r,\quad x_j,\quad y\prod_{i\ne j}x_i-1
\]
have no common geometric zero, with \(y\) a new
variable. The arithmetic Nullstellensatz supplies
a nonzero integer \(B_j\) of height \(H^{O_{m,D}(1)}\)
in their integer polynomial ideal. For any pairwise
coprime integral sample point, every prime dividing
\(x_j\) leaves \(\prod_{i\ne j}x_i\) invertible in
the residue field. Evaluating the ideal relation
at the inverse choice of \(y\) shows that this prime
divides \(B_j\). Therefore \(x_j\le H\) is supported
on \(\operatorname{rad}(B_j)\), with
\(H^{o_H(1)}\) possibilities in the usual exponent
sense, uniformly in \(B_j\) under its fixed polynomial
height cap. One may see this via the Rankin estimate
at a fixed exponent and the elementary
\(\omega(B_j)=o(\log H)\). Applying this on all
coordinates counts the low-degree component at
cost \(H^\varepsilon\).
[Krick–Pardo–Sombra, Theorem 1](https://arxiv.org/abs/math/9911094)

## High-degree covers

For a geometrically integral rational component of
positive dimension \(j\le m-1\) and degree
\(D\ge c_{m,s}\sqrt E\), project finitely and
birationally to the hypersurface of degree \(D\)
in \(\mathbb P^{j+1}\) from the standard
bounded-coefficient projection. The projective
height of an integral affine point of height
\(H\) is at most \(H\). Heath-Brown's covering theorem, as recorded by Salberger, gives proper child loci at cost
\(H^{(j+1)D^{-1/j}+\varepsilon}\), with
bounded degrees in the fixed parameters.
The pullback remains proper. Dimension drops at most
\(m-1\) times, with the \(\varepsilon\)-allowance divided among the steps.
Zero-dimensional samples have bounded
cardinality. The worst exponent is at most
\(C_{m,s}E^{-1/(2(m-1))}\) after bounding the sum
of the cover costs. This proves (EB1), and the \(m=1\)
case is the ordinary univariate root bound.
[Castryck–Cluckers–Dittmann–Nguyen, Proposition 4.3.2](https://msp.org/ant/2020/14-8/ant-v14-n8-p09-s.pdf)
· [Salberger, Theorem 0.10](https://doi.org/10.1112/plms.12508)

## Primitive mixed collision estimate

Let the two complete-value dictionaries \(A,B\) have root
caps \(l_A,l_B\) in the scale \(n=X^L\), with fixed
anatomic depth at each small budget and fixed numerical
heights at most \(X^{O(1)}\).
In each choose a bounded number of original
squarefree-root slots (the number is bounded before
the outer budget limit) whose power orders are
at least \(E(L)\to\infty\). Denote their total root
caps by \(\alpha_A,\alpha_B\). All remaining root
tuples, including coefficient records if needed,
must be enumerable with cost
\(n^{l_A-\alpha_A+o_X(1)}\) and
\(n^{l_B-\alpha_B+o_X(1)}\), respectively.

Write \(E_{\rm off}^{\rm prim}(A,B)\) for the count
of primitive ordered pairs \((a_1,b_1),(a_2,b_2)\)
with both diagonal permutations removed and with
\[
                    a_1+b_1=a_2+b_2.
\]
Then in the fixed-budget-then-small-budget limit
\[
 E_{\rm off}^{\rm prim}(A,B)
     \ \ll_{L,\varepsilon}\ 
 X^\varepsilon
 n^{\,2(l_A+l_B-\alpha_A-\alpha_B)+o_L(1)}.             \tag{EB2}
\]
There is no common power-order condition and no common
prime coloring is necessary. As usual, the small-budget
error incorporates \(o(1)\) errors in the prescribed
root caps.

Freeze the unselected root tuples of the four
occurrences at the exponent in (EB2), and split the
selected squarefree roots by their exact prime
incidence patterns in these occurrences and in their
original slots. Within any one occurrence, its
selected roots have disjoint supports. Since
\((a_i,b_i)\) is primitive for each \(i\), a prime can
occur in at most two of the four occurrences.
For every pattern the corresponding product of
primes is a positive integer variable, pairwise
coprime with every other pattern variable. Each
variable is absent from at least two terms and
has an exponent at least \(E(L)\) on each
term where it occurs, even when those exponents
are different. A prime of a selected variable
can also divide a frozen coefficient; this
does not affect pairwise coprimality among the
variables and (EB1) has no condition on the
coefficients' supports. The relation after removal of the two diagonal
permutations has no proper vanishing subsum:
a two-term cancellation between opposite signs
would force one of exactly those permutations. If the dictionaries
use disjoint prime pools, the permutation is
absent apart from the common value \(1\), which
does not give distinct pairs on its own.

There are only boundedly many pattern variables
in terms of the chosen-slot bound; each is
at most \(n^{O(1)}\) since it divides one chosen
root. Apply (EB1) to the nondegenerate collision.
Its normalized exponent tends to zero since
\(E(L)\to\infty\). Reconstruction from the
patterns and frozen data is exact, which
proves (EB2) with the permutation diagonal
excluded. The same argument for fixed signs
counts the corresponding nondegenerate
equal differences. In a globally
letter-colored primitive box this is exactly
the off-diagonal mixed additive energy.

In particular, if
\(\alpha_A+\alpha_B>(l_A+l_B)/2\) by a fixed
margin, the primitive mixed energy has a
power-saving off the exponent \(l_A+l_B\) of the full dictionary-product cap, even
with different large valuation indices and
shared selected prime factors. If the selected
bounded lists account for all root mass in the
limit, the off-diagonal exponent tends to zero
in the \(n\)-scale. This requires selected
positions of substantial total normalized mass;
it does not cover arbitrary wholly diffuse
profiles where every bounded list has
vanishing mass.


## Moving numbers of exposed slots

The same argument allows the number of selected high
slots to grow in the outer limit. More precisely the
Brownawell–Masser step gives
\[
                    E\ll_s D^2+mD .
\]
When \(m\ll\sqrt E\), the threshold \(D\gg_s\sqrt E\)
can therefore be chosen independent of \(m\).
The sum of the high-degree covering exponents is then
at most
\[
             O_s\!\left(m^2 E^{-1/(2(m-1))}\right).
\]
All other factors occur in the implied constants at a
fixed budget. This tends to zero if
\(m\le c\log E/\log\log E\) for a sufficiently small
absolute \(c>0\). Thus (EB2) permits this many pattern
variables as \(L\to0\), with the height limit taken
first. In globally letter-and-slot-colored
dictionaries there are at most three pattern
variables per selected slot (the two exclusive
factors and the shared gcd), so the same order of
growth is allowed directly for the selected-slot
count. Without slot coloring a pair of primitive
pairs has at most two occurrences per prime, giving
a quadratic bound for the number of patterns in the
number of original selected slots.

The analogous release on the original three-term
primitive equation permits this many selected high
slots taken across all three letters, at once: freeze
the other root factors and use (EB1) directly.
Weak nonatomic degree profiles are compatible with
a growing, slowly sparse collection of occupied
original valuation positions. Keeping the count of
positions fixed throughout the outer limit would
lose that sector.

## A uniform small-budget count for sparse valuation positions

Fix \(\kappa>0\) and \(0<\eta<1/3\), and let
\(E=\lceil\kappa/\lambda\rceil\). There is an absolute \(c_*>0\)
such that the following holds for sufficiently small fixed
\(\lambda>0\). Suppose that, in a primitive \(a+b=c\le X\),
all but a prime support of radical at most
\(X^{(1/3-\eta)\lambda}\) is carried by at most
\[
                           c_*\frac{\log E}{\log\log E}
\]
letter and valuation pairs, with each chosen valuation at least
\(E\). The chosen indices and the defect support may depend on the
triple. Then the count is
\[
                    O_{\kappa,\eta,\lambda,\varepsilon}
                          (X^{(1/3-\eta/2)\lambda+\varepsilon}).
 \tag{EB3}
\]
This uses all three primitive letters. With the original radical
budget, the hypothesis includes retaining more than two thirds of
the budget by a fixed margin on the chosen positions. No spacing
or integer divisibility relation between the positions is imposed.

To prove it, move all primes of original valuation greater than
\(K=\lceil 24/(\eta\lambda)\rceil\) into the defect support.
Their total radical has height at most \(X^{3/K}\), from the
product of the three complete letters. Thus this enlarges the
normalized defect cost by at most \(\eta/8\). There are finitely many
choices of surviving positions at each fixed \(\lambda\). Enumerate
the joint defect radical, its letter assignments, and the three
integer values supported on it at cost
\(X^{(1/3-7\eta/8)\lambda+\varepsilon}\), using the fixed-support
count at height \(X\) and resizing \(\varepsilon\).
Fixing them gives a nondegenerate three-term relation with
pairwise coprime high roots, by primitivity and their disjoint
valuation layers. Each root has height at most \(X^{1/E}\).
Apply (EB1) with the moving-width estimate and choose \(c_*\)
small enough that its exponent in this root height tends to
zero. Dividing that exponent by \(E\lambda\ge\kappa\), it is less
than \(3\eta/8\) for sufficiently small \(\lambda\). Cases with
fewer than two roots have bounded counts at each fixed choice of
orders and coefficients. Summing proves (EB3), with the
height limit taken at fixed \(\lambda,K\).

An increasing list this short can still have an arbitrary diffuse
weak limit on a bounded interval of allowed scaled valuation
degrees, by approximating the measure by increasingly many
point masses and rounding each scaled degree to the nearest
available integer position. Thus the hypothesis controls
microscopic support size and is compatible with all three
limiting measures being nonatomic. The [earlier congruence release](ABC_ANISOTROPIC_RATIONAL_CONGRUENCES.md#an-exponential-dimension-constant)
uses an actual full modulus and two free letters. The direct count
above requires the specified radical mass on the original three
letters and also works when no such modulus has been selected.
Uniform consecutive or prime laws with order \(\lambda^{-1+o(1)}\) positions can leave vanishing
mass on each such list. Those profiles need a different count.

## The diagonal and the remaining third-letter obstruction

Let \(C\) be the actual third-letter dictionary,
and let \(r_c\) denote the number of primitive
pairs \((a,b)\) in the two free dictionaries
with \(a+b=c\). The count \(T=\sum_{c\in C}r_c\)
has exact diagonal contribution \(T\) to
\(\sum_{c\in C}r_c^2\). Off-diagonal control yields
\[
 T\le |C|+\sqrt{|C|\sum_{c\in C}r_c(r_c-1)} .
\]
If the two free dictionaries overlap, a permutation
diagonal can be treated separately; it is absent
for globally distinct letter colors.
Even a zero off-diagonal term allows one pair
per target. Consequently an off-diagonal bound
alone leaves the distinct-target occupancy at the
same exponent as the cap on the third dictionary.
More precisely, if \(O_C\) is the number of
targets with at least one primitive representation,
then \(0\le T-O_C\le\frac12\sum_{c\in C}r_c(r_c-1)\).
Small off-diagonal energy makes the occupancy
accurate; it does not bound it.

For bounded lists at the balanced endpoint,
each free letter has radical mass \(1/3\), and
the two free nonatomic masses sum to at least
\(1/3\). Bounded selected lists recover only
atomic mass in that limit, so
\(\alpha_A+\alpha_B\le1/3+o(1)\).
The exponent in (EB2) is therefore at least
\(2/3+o(1)\) on these endpoint hypotheses.
If a whole-coordinate lock of cost \(t\) is
already charged in this pair, the free residual
count sought is \(1/3-t\), and the same
degree-mass calculation gives an off-diagonal
exponent at least \(2(1/3-t)+o(1)\).
The bound has a substantive arithmetic
gain in supported packets, but does not
remove unrestricted balanced saturation.


For a fixed-target powerful-number application, let \(u+v=N\)
with \((u,v)=1\), and suppose the prime powers on \(u,v\)
belong to at most \(k\) prescribed valuation layers on each
letter, all of order at least \(E\). For fixed \(k\ge1\) and
fixed choices of the layers, the ordered count is
\[
       O_{k,\mathbf e,\varepsilon}
             (N^{\,C_k E^{-1-1/(2(2k-1))}+\varepsilon}).
\]
The disjoint squarefree layer roots are pairwise coprime, and
each is at most \(N^{1/E}\). Use (EB1) with the constant term
\(-N\) and at most \(2k\) free variables. All counted positive
pairs give nondegenerate relations. The bound concerns the
prescribed valuation layers of the primitive summands. It
does not count arbitrary \(E\)-full summands using finitely
many layers independently of \(E\). The different coefficients,
orders and overlaps in actual paths of Erdős #943 need their
own constraints. The mixed collision statement extends the
prime-incidence separation behind those paths to high
valuation layers.
[Fixed-target kernel-gcd recovery in Erdős #943](ERDOS_943_PELL_RECOVERED_KERNEL_GCD_WEIGHT.md)
