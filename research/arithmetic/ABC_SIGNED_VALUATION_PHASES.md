# Releasing shifted fine valuation phases with two signed bases

The exact power requirement can be replaced by concentration in
shifted fine valuation phases. The argument uses an actual full
divisor of a primitive letter, with no precision from repair
multipliers.

## A rational congruence bound for distinct exact outputs

Let \(m\) be \(k\)-full and \(x,y\) positive rational units
at the primes of \(m\). Assume \(H(x)\le m^{C_1/k}\),
\(H(y)\le m^{C_2/k}\), \(C_i\ge1\), and
\(k\ge B C_1 C_2\) for a sufficiently large absolute \(B\).
Let \(P,Q\) be fixed nonzero integers. Count distinct exact
rational values \(W=x^P y^Q\) satisfying
\(W\equiv K\pmod m\), for a fixed unit \(K\).
Then the count is \(O(m^{1/k+\varepsilon})\), with constant
allowed to depend on the fixed heights, orders, and \(\varepsilon\),
uniformly in \(K\).

Choose one input pair for each exact output. The rectangular
two-root determinant extends directly to these rational rows.
For local charts use unit residues of the reduced fractions.
Clear each evaluation row by the product of the respective
denominators to degrees \(D_x,D_y\). This is a local unit and
the row height is at most \(m^{(D_x C_1+D_y C_2)/k}\).
The rank-half construction gives a polynomial of this same
coefficient-height exponent. For the logarithmic Jacobian
use the numerator after row clearing. Its valuation controls
the local loss just as for the integer points. On a component
where the logarithmic Jacobian vanishes identically, the exact
output is constant, which costs one chosen pair per component.
On the other components apply Hensel to the cleared monomial
equation and the auxiliary polynomial. The good modulus in
terminal reconstruction has weight exceeding \(m^{1/2}\) if
\(B\) is sufficiently large. Equality of rational coordinates
then follows on cross-multiplying, since each height is
less than \(m^{1/4}\) with a fixed margin. These adaptations
leave the same bound and do not require disjoint supports in
the two input coordinates.
[The integer rectangular congruence proof](ABC_TWO_VALUATION_GRIDS.md#a-rectangular-determinant-over-a-full-modulus)

For a fixed amplification order \(\ell\) and full exact outputs
of polynomial height in \(X\), color primes by the sign of
valuation in the exact rational outputs. A common coloring
keeps \(X^{-o_X(1)}\) of the outputs by averaging, since the
number of prime divisors of any numerator and denominator
is \(o(\log X)\). With this prime coloring, in an \(\ell\)-fold
product there is no cancellation between numerators and
denominators. The product has \(X^{o_X(1)}\) presentations
by the divisor bound. Choose one actual coordinatewise tuple
product for each distinct output product. Its height caps
are multiplied by at most \(\ell\) and its congruence target is
\(K^\ell\). Taking \(\ell\)-th roots gives the amplified
exponent \(1/(\ell k)\) in base \(m\), subject to the same
determinant condition at the amplified parameters. Coloring
the signed exact outputs avoids any assumption on overlap of
the two basis variables.

## From shifted phase concentration to small height product

Use the anatomic box scale \(n=X^L\) and any primitive signed
additive equation among \(U,V,T\). Retain a prescribed collection
of prime-resolved root positions on \(U,V\); freeze the rest on
those two letters at total normalized root cost at most \(f\).
Take a positive integer \(P=P(L)\le \Lambda L^{-1/2}\), where
\(\Lambda>0\) is fixed. Suppose a positive integer \(k_1\) makes
every retained degree \(e\) satisfy \(\|k_1 e/P\|\le\theta\).
If the retained positions are specified by box parameters,
the fixed number of choices at each \(L\) can be summed before
the height limit. Define
\[
 d=\gcd(k_1,P),\qquad P'=P/d,\qquad k_1'=k_1/d .
\]
Require \(P'\to\infty\) as \(L\to0\). Take \(Q\in[P',2P')\)
with \(Q k_1'\equiv1\bmod P'\). Take \(r_e\) to be the
centered representative of \(k_1'e\bmod P'\) and put
\(i_e=(e-Qr_e)/P'\). Then
\[
                   e=P'i_e+Qr_e,\qquad |r_e|\le\theta P'.
 \tag{SP1}
\]
Apply the exponents \(i_e,r_e\) to the retained root factors
with opposite signs for opposite letters. This gives rational
coordinates \(z,w\) whose signed monomial \(z^{P'}w^Q\) is
exactly the ratio of the retained factors. Given the frozen
factors, this ratio recovers the complete primitive pair
and therefore the triple.

Write \(\xi_z,\xi_w\) for the numerical logarithmic height
caps of the rational coordinates in base \(X\).
For the displayed estimates, the total retained numerical
weight is at most \(2+o(1)\) and their root mass is at most
\(1+o(1)\) in normalized units. Taking absolute contributions
on the two letters gives
\[
 \xi_z\le(2+o(1))/P'+(2+o(1))L\theta P',\qquad
 \xi_w\le(1+o(1))L\theta P'.
\]
Both tend to zero and
\[
                 \xi_z\xi_w\le
                 (2+2\Lambda^2\theta+o(1))L\theta .
 \tag{SP2}
\]
No retained prime needs to have valuation divisible by \(P\).

Suppose the third letter supplies at most
\(n^{\mu+o_X(1)}\) records of full divisors \(m\) with
\(\log_X m\ge w_0>0\) and fullness threshold at least
\(k=\lfloor \kappa/L\rfloor\), for fixed positive \(w_0,\kappa\).
The record and the frozen factors give a unit congruence on
the exact signed monomial. Fix a desired normalized chart
allowance \(\eta>0\) and choose \(\ell\) large enough
that \(1/(\ell k L)<\eta\) for small \(L\), using \(m\le X\).
Take amplified height constants
\(C_i=\max(1,\ell k\xi_i/w_0)\). They satisfy
\[
 \frac{C_1C_2}{k}\le
   \frac1k+\frac{\ell(\xi_z+\xi_w)}{w_0}
      +\frac{\ell^2 k\xi_z\xi_w}{w_0^2}.
\]
The right side is below \(1/B\) if a sufficiently small
fixed \(\theta>0\) is chosen after \(\ell\) and before \(L\),
by (SP2). The rational result therefore gives the count
\[
              O_{L,\eta,\varepsilon}
                    (X^{(f+\mu+\eta)L+\varepsilon})
 \tag{SP3}
\]
for all sufficiently small fixed budgets in this range.
All choices of exponent bases and depths are fixed in the
height limit. The same proof works whenever there are two
signed bases of numerical heights tending to zero with product
\(o(L)\), by taking \(\ell\) fixed before decreasing \(L\).
In particular on a common growing grid with deviations
\(|r_e|\le D=o(P)\) one may use \(Q=1\) directly and bound
the product by \(O(LD(1+LD)/P)=o(L)\) if \(LD=o(1)\).
The \(P\le\Lambda L^{-1/2}\) assumption is not needed for
this latter common-grid conclusion.

## Arbitrary shifted clusters and logarithmically many windows

Suppose the phases \(e/P\bmod1\) of the retained positions
on the two letters lie in intervals of radius \(\rho\) about
\(s\) arbitrary centers. For an integer \(H>2/\theta\),
simultaneous Dirichlet approximation gives one integer
\(1\le k_1\le H^s\) taking every center within \(1/H\) of
an integer. Indeed partition \([0,1)^s\) into \(H^s\)
cubes and compare the \(H^s+1\) integer multiples of the
vector of centers modulo one. It suffices that
\[
                  H^s\le P^{1/2},
              \qquad H^s\rho\le\theta/2 .
 \tag{SP4}
\]
The first inequality ensures \(P'\ge P/k_1\ge P^{1/2}\),
even for composite \(P\). The second controls the changes
within the intervals. The conditions \(D H^s\le P^{1/2}\) and
\(D H^s
ho\le	heta/2\) work instead if the centers are bounded
integer combinations of \(s\) generating phases plus
rational offsets of denominator dividing a fixed \(D\).
First multiply by \(D\), and approximate the generating
phases to accuracy small enough for the stated coefficient
bound.

For a concrete application choose \(P=\lfloor L^{-1/2}\rfloor\).
Fix any \(0\le\beta<1/2\). Then for some \(c_0>0\), depending
on the prescribed counting margins and \(\beta\), it is enough
that the retained degree positions occupy at most
\(c_0\log(1/L)\) real intervals of radius \(L^{-\beta}\)
in the original integer degree coordinate. Indeed their
phase radii are at most \(O(L^{1/2-\beta})\); take
\(c_0\log H<\min(1/4,\,1/2-\beta)\) in (SP4).
The intervals can be spread across the entire numerical
degree band. Thus the conclusion allows a diverging number
of occupied windows and a nonatomic weak scaled-degree limit.

## Consequences for balanced boxes and target repair

At a balanced endpoint, suppose the third letter has positive
nonatomic mass in its scaled degrees. Choose a small fixed
root mass \(\mu>0\) from a compact positive band of that
nonatomic part. Its whole root positions give the required
positive \(w_0,\kappa\) at cost \(\mu\), even if other primes on
the letter are macroscopic. If the two free letters can retain
root mass \(t>1/3\) in phase clusters as above, the frozen cost
is \(f=2/3-t+o(1)\). Choose \(\mu,\eta\) before the height
limit to leave a fixed saving from (SP3). This also gives
a positive-width box exclusion whenever the stated margins
are uniform.
[Selecting whole root positions on a diffuse band](ABC_WIDE_COFACTOR_CONGRUENCES.md#small-full-blocks-in-diffuse-bands)

More precisely, take any sequence \(P(L)\to\infty\) with
\(P(L)\ll L^{-1/2}\), and take simultaneous weak limits of the
root measures pushed forward by \(e\mapsto e/P\bmod1\) on the
circle. In the two free phase limits the sum of the atomic
masses is at most \(1/3\) if the third letter has positive
nonatomic scaled-degree mass. Otherwise a fixed finite list
of phase atoms has combined mass greater than \(1/3\), so
choose continuity intervals of sufficiently small fixed
width around it and apply (SP4). This inequality applies to
the pair opposite any such third letter. If all three
scaled-degree measures have positive nonatomic mass it
implies total phase-atomic mass at most \(1/2\). In general
at least two letters have positive nonatomic mass, giving
total phase-atomic mass at most \(2/3\).
[Pairwise endpoint diffusion](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

For the target-repair sector impose \(T=c\) and
\(\operatorname{rad}(c)\ge c^{(1/3-\xi)\lambda}\).
In an original budget \(L\), (SP3) and
\(c\le Y^{1/(M(1/3-\xi)\lambda)}\) give repaired count exponent
\((f+\mu+\eta)L/[M(1/3-\xi)\lambda]+\varepsilon\).
For \(L/\lambda\) close to one and \(f+\mu+\eta\) below
\(1/3-\xi\) by a fixed margin this is strictly below \(1/M\).
All supplied moduli were already present on the primitive target.
[Minimal target-repair equivalence](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## A microscopic obstruction to the small-product criterion

A formal diffuse profile with consecutive equally weighted
integer degree positions on \(eL\in[2,4]\) satisfies
an actual discrete cap \(C L\) on each degree mass.
Whenever a retained mass \(s>0\) satisfies such an atom cap,
every signed two-basis assignment \(e=P i+Q r\) obeys
\[
       (\mathcal I+Ls/8)(\mathcal R+Ls/8)
                  \ \ge\ \frac{s^3}{128C}\,L,\qquad
       \mathcal I=L\int |i|,\quad \mathcal R=L\int |r|.
 \tag{SP5}
\]
Here both integrals count the retained root masses on the
two free letters. In fact Markov's inequality leaves at
least \(s/2\) of the mass with
\(|i|\le4\mathcal I/(Ls)\) and \(|r|\le4\mathcal R/(Ls)\).
That mass requires at least \(s/(2CL)\) distinct degrees.
Each integer pair determines just one degree, so the product
of the two interval lattice counts is at least this large.
Rearranging gives (SP5). By prime-resolved primitivity the
two rational height exponents are at least half the
corresponding absolute sums. Thus when both exponents tend
to zero their product cannot be \(o(L)\) on this formal
profile. The result only concerns the small-product
sufficient condition. It leaves open a saving at a fixed
product of order \(L\) through additional arithmetic.


### An entropy form for any number of signed bases

Choose a valuation randomly from the retained normalized root
mass \(s>0\), divided by \(s\), and write \(E\) for the
integer valuation drawn. Let \(E=\sum_{j=1}^d P_j I_j\)
be any signed integer basis representation, including
different assignments at different primes. If the resulting
rational height exponents are \(\xi_j\), then
\[
                \prod_{j=1}^d(\xi_j+Ls)
                  \gg_d (Ls)^d \exp H(E).                      \tag{SP6}
\]
Indeed the rational height sees at least half the absolute
primewise sum, so
\(\mathbb E|I_j|\le2\xi_j/(Ls)\).
An integer variable with mean absolute value at most \(v\)
has entropy at most \(\log(C(1+v))\). Compare its law with
the discrete exponential law proportional to
\(\exp(-|i|/(1+v))\) and use nonnegativity of relative entropy.
Since \(E\) is determined by the integer vector,
\(H(E)\le\sum_jH(I_j)\), which proves (SP6).
A mass cap \(CL\) at each actual degree gives
\(\exp H(E)\ge s/(CL)\), recovering the product-height
lower bound of order \(L^{d-1}\) at fixed retained mass.

For example, if the two letters instead have separate single
bases of height exponents \(\xi_1,\xi_2\) and a shared
correction of height exponent \(E_0\), apply (SP6) on each
letter separately. Under positive retained masses and the
same microscopic atom caps it gives
\((\xi_i+L)(E_0+L)\gg L\) for \(i=1,2\).
When all three heights tend to zero, each is large relative
to \(L\), and \(\xi_1\xi_2 E_0\gg L^2/E_0\).
Thus independent single-base rounding has a larger
combined height product than the shared two-basis
representation can attain in the clustered sector.


The small-product criterion also covers several unrelated
signed valuation grids of bounded total rank when the
product of their coefficient widths is \(o(1/L)\).
A short kernel relation in the actual valuation strip
reduces the rank while preserving every degree, allowing
the two-basis estimate above. A sparse diffuse profile
alone does not imply the width condition.
[Rank compression and sparse-profile obstruction](ABC_SPARSE_VALUATION_GRID_COMPRESSION.md)
