# Consecutive valuation profiles with simultaneous local power tests

Consecutive valuation degrees can pass every test in a fixed finite
collection of local generalized Fermat equations. For the quadratic
signature, they also occur on actual primitive sums at polynomial
height. The additional primes in the quadratic realizations all have
even valuations. The limiting small-budget count gives a lower bound
on the radical of those additional primes for almost every prescribed
support.

## Prescribed profiles and simultaneous local points

Fix \(J\ge2\) and a finite set of signatures \((q_1,q_2,q_3)\) with
\(q_i\ge2\), including \((2,2,2)\). Let \(W\) be the least common
multiple of their orders. Set
\[
 I_J=\{J,\ldots,2J-1\},\qquad S=\sum_{e\in I_J}e,\qquad
 X_0=(2Y)^{WS},\qquad L=3J/S.
\]
The parameters \(J,W\) are fixed in the height limit \(Y\to\infty\).
The list may depend on \(J\), for example by using all signatures
whose orders are at most \(2J\).

For all sufficiently large \(Y\) there are \(X_0^{L+o_Y(1)}\)
prescribed sets of \(3WJ\) primes in \([Y,2Y]\). For each set
assign \(W\) primes to each pair \((i,e)\) of a letter and a degree.
Write \(P_i=\prod_{p\text{ assigned to }(i,e)}p^e\), with the
product over all assigned degrees. The three \(P_i\) are pairwise
coprime and
\[
 Y^{WS}\le P_i\le X_0,\qquad
 Y^{WJ}\le\operatorname{rad}(P_i)\le X_0^{L/3}.                 \tag{CV1}
\]
For each tested signature put
\(D_i=\prod_{p^e\parallel P_i}p^{e\bmod q_i}\). The equation
\[
 D_1u_1^{q_1}+D_2u_2^{q_2}=D_3u_3^{q_3}                       \tag{CV2}
\]
has a nonzero point over every completion of \(\mathbb Q\).
At each prescribed \(p^e\parallel P_i\), its local point can
have \(v_p(u_i)=\lfloor e/q_i\rfloor\) and both other coordinates
units, realizing the prescribed complete valuations.

Choose \(C_W\ge2W\) so large that
\(p+1-(W-1)(W-2)\sqrt p-3W>0\) for \(p>C_W\).
Restrict the available primes to \(1\bmod A_0\), where
\[
 A_0=2W\prod_{p\le C_W}p^{2v_p(W)+1}.
\]
For each such prime \(p\) fix a character \(\chi_p\) of order
\(W\) with values in a common group of \(W\)-th roots of unity.
Colour a pair \(p<p'\) by \((\chi_p(p'),\chi_{p'}(p))\).
Finite Ramsey gives a monochromatic set of \(m=3WJ\) vertices
among every sufficiently large fixed number \(R\) of primes.
For \(n\) available primes, counting the \(R\)-sets containing
each monochromatic \(m\)-set gives at least
\(\binom nR/\binom{n-m}{R-m}\gg_{J,W} n^m\) such \(m\)-sets.
The prime number theorem in fixed progressions gives
\(n\asymp_{A_0}Y/\log Y\), and the trivial upper count \(n^m\)
gives the same height exponent.
[Ramsey, On a Problem of Formal Logic](https://doi.org/10.1112/plms/s2-30.1.264) · [de la Vallée Poussin, Recherches analytiques sur la théorie des nombres premiers](https://circolomatematico.unipa.it/download/Charles-Jean_de_La_Valle_e_Poussin-Collected_Works-Oeuvres_Scientifiques_Vol.1.pdf)

Split the increasing primes in each monochromatic set into
consecutive blocks of length \(W\), one for each pair \((i,e)\)
in a fixed order. Relative to any prescribed prime outside a block,
all its members lie on the same side in this ordering. The character
at that prime takes a constant value on the block, so its product
is a \(W\)-th power residue. At \(p^e\parallel P_i\), both other
\(D_j\), and also \(-1\), are \(W\)-th power residues modulo \(p\).
Put \(u_i=p^{\lfloor e/q_i\rfloor}\), put one other coordinate
equal to \(1\), and solve for the last by Hensel. The needed root
modulo \(p\) is a unit and \(p\nmid W\).

At \(p\le C_W\), all \(D_i\equiv1\bmod p^{2v_p(W)+1}\).
Set \(u_2=0,u_1=1\) and lift a \(q_3\)-th root of \(D_1/D_3\)
near \(1\). At any other unprescribed prime, substitute
\(u_i=t_i^{W/q_i}\). The smooth degree-\(W\) diagonal plane curve
over \(\mathbb F_p\) has a point off the coordinate axes by the
Weil bound, and the point lifts by nonsingular Hensel. A positive
real point is immediate.
[Weil, Numbers of Solutions of Equations in Finite Fields](https://www.ams.org/bull/1949-55-05/S0002-9904-1949-09219-4/S0002-9904-1949-09219-4.pdf)

With normalization \(L\log X_0\), each prescribed radical measure
has mass \(1/3+o_Y(1)\) and numerical first moment \(1+o_Y(1)\).
Each exact degree in \(I_J\) has mass \((1+o_Y(1))/(3J)\).
After the height limit, \(J\to\infty\) gives
\((1/3)\operatorname{Unif}[2,4]\) in scaled degree \(Le\) and
odd-valuation mass \(1/6\) on each letter. A divisibility grid
\(k\mathbb N\) carries only \(O(1/k+1/J)\) of the prescribed mass.
A signed grid of fixed rank and digit-width product \(o(J)\)
hits only \(o(J)\) degrees by counting its coefficient vectors,
uniformly in its integer steps. Thus growing divisibility orders
and fixed-rank signed grids of sublinear width product both
capture a vanishing fraction.

## Quadratic realization and height

For every prescribed set a positive primitive triple \(a+b=c\)
has the specified valuations on the chosen primes, with
\[
 X_0^{1-o_Y(1)}\le c\ll X_0^{6+17L/3}.                        \tag{CV3}
\]
At every other prime all valuations are even, so the odd-valuation
kernels of the complete terms equal the \(D_i\) from the quadratic
signature. These fresh even-valuation primes still contribute to the
complete radical.

Write \(Q(\mathbf z)=D_1z_1^2+D_2z_2^2-D_3z_3^2\) and let \(B\)
be its associated bilinear form. By the local points and
Hasse–Minkowski the form is isotropic over \(\mathbb Q\).
Cassels's small-zero theorem gives an integral primitive zero
\(\mathbf z_0\) of size \(\ll\max_iD_i\le X_0^{L/3}\).
[Cassels, the local-global principle for conics](https://people.math.ethz.ch/~halorenz/4students/ECC/Cassels3.pdf) · [Cassels, Bounds for the least solutions of homogeneous quadratic equations](https://www.cambridge.org/core/journals/mathematical-proceedings-of-the-cambridge-philosophical-society/article/abs/bounds-for-the-least-solutions-of-homogeneous-quadratic-equations/9C3C57BC5549819DF3DBF5F3BC42DCE8)

Fix \(p^e\parallel P_i\), put \(r_p=\lfloor e/2\rfloor\), and
let \(\mathbf s_p\) be the corresponding local zero. At least
one coordinate \(z_{0j}\) with \(j\ne i\) is a unit at \(p\).
If both others were divisible by \(p\), primitivity would make
\(z_{0i}\) a unit and give the three terms of \(Q(\mathbf z_0)\)
valuations \(0\) or \(1\), at least \(2\), and at least \(2\),
respectively. Flip the sign of \(s_{pj}\) if needed to make
\(B(\mathbf z_0,\mathbf s_p)\) a unit, since the two choices of
that bilinear value differ by a unit. The integral quadratic map
\[
 F(\mathbf z)=Q(\mathbf z)\mathbf z_0-2B(\mathbf z_0,\mathbf z)\mathbf z
\]
satisfies \(Q(F(\mathbf z))=0\) and
\(F(\mathbf z_0+\mathbf s_p)=-2B(\mathbf z_0,\mathbf s_p)\mathbf s_p\).
Impose \(\mathbf z\equiv\mathbf z_0+\mathbf s_p\bmod p^{r_p+1}\)
for every prescribed prime. The reflected vector then has exactly
the required root valuations. The total CRT modulus is
\[
 \mathcal M=\prod_{p\mid P_1P_2P_3}p^{r_p+1}
   \le (P_1P_2P_3)^{1/2}\operatorname{rad}(P_1P_2P_3)
   \le X_0^{3/2+L}.
\]
Choose \(\mathbf z=O(\mathcal M)\) in this CRT class with all
three \(F_i(\mathbf z)\ne0\). Each coordinate polynomial is
nonzero. If \(z_{0i}\ne0\), its coefficient of \(z_j^2\) is
nonzero for \(j\ne i\). If \(z_{0i}=0\), it is the nonzero product
\(-2B(\mathbf z_0,\mathbf z)z_i\). Their product has degree at most
six, so it is nonzero at some point on the seven-point grid in
each coordinate obtained by adding \(0,\ldots,6\) times
\(\mathcal M\) to a CRT representative.

Divide \(F(\mathbf z)\) by its coordinate gcd, a unit on all
prescribed primes, to obtain \(\mathbf u\). The values \(D_i u_i^2\)
are pairwise coprime by the prescribed valuations and because an
outside prime dividing any two forces all three coordinates of
\(\mathbf u\) divisible by it. They satisfy
\(D_i u_i^2\ll(\max D_i)^5\mathcal M^4\), proving the upper
height bound. The lower bound follows from \(P_3\mid c\).
The exact quadratic norm and ray-class identities hold in each
orientation for these actual triples, including ray conditions
from root divisors disjoint from the parity kernels.
[Quadratic norm and ray-class identities](ABC_WEIGHTED_SQUARE_RAY_CLASS_MOMENTS.md)

## Fresh radical at bounded height

Fix \(B_0>6\) and \(0\le b<2\). For sufficiently large fixed
\(J\), at most
\[
 O_{J,W,B_0,b}\bigl(X_0^{(1-\delta(b))L}\bigr)                 \tag{CV4}
\]
of the prescribed sets have a primitive realization \(c\le
X_0^{B_0}\) with fresh radical
\(\prod_{p\mid abc,\ p\nmid P_1P_2P_3}p\le X_0^{bL}\), for
some \(\delta(b)>0\). Thus almost all of the \(X_0^{L+o_Y(1)}\)
sets have realizations within the height range by (CV3) and
require a fresh radical above \(X_0^{bL}\) in every such
realization.

Indeed any triple counted here has \(\operatorname{rad}(abc)\le
X_0^{(1+b)L}\). Its height is at least \(X_0^{1/2}\) for large
\(Y\), and it contains only \(O_{J,W,B_0}(1)\) primes of size at
least \(Y\), so it can realize only \(O_{J,W,B_0}(1)\) prescribed
sets. Cover \(\log_{X_0}c\in[1/2,B_0]\) by a finite geometric mesh
of intervals \([h,(1+\tau)h]\). In each apply the limiting
small-budget theorem with fixed
\(\lambda=(1+\tau)(1+b)L/h\), which strictly exceeds the needed
relative cap. Choose its allowance \(\eta>0\) and then
\(\tau>0\) so that
\((1/3+\eta)(1+\tau)^2(1+b)<1\). Take \(J\) sufficiently large
to meet the resulting budget threshold in every interval, and
then choose the height exponent allowance sufficiently small
at that fixed \(J\). Summing over the intervals proves (CV4).
[Limiting small-budget coefficient](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

The construction rules out deriving a valuation grid from the
prescribed profile and any finite list of these local power tests,
even with global quadratic solvability. An endpoint proof still
needs the caps on the complete radicals or additional global
control of realized heights and fibres.
