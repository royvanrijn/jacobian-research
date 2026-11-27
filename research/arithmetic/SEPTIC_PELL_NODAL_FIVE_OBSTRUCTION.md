# Ordinary nodes at five do not support five septic Pell directions

**Theorem SN1.** Let \(P\in\mathbb Z_5[t]\) have degree seven
and unit leading coefficient, and suppose \(P^2+1\) is squarefree
over \(\mathbb Q_5\). If every geometric root of the reduction
\(\bar P^2+1\) has multiplicity at most two, then
\[
 \dim_{\mathbb F_3}
 \operatorname{Jac}(R^2=P^2+1)(\mathbb Q_5)[3]\le4 .
\]
Thus a degree-seven rational square-plus-one model with five
rational \(3\)-torsion directions cannot have such a reduction.
For an integral degree-preserving model at five, at least one
root of \(\bar P^2+1\) must have multiplicity at least three.
The smooth case is the obstruction of
[SP2](SEPTIC_PELL_FIVE_DIRECTION_CONSTRUCTION.md).

The last possible nodal residue class fails because one
three-torsion class on its normalization does not lift to the
generalized Jacobian. Counting the normalization, torus and
component groups separately would leave this class open.

## The local rank budget

Write \(\bar F=\bar P^2+1=g^2 F_0\) with \(g\) monic and
\(g,F_0\) squarefree and coprime. Since the derivative of
\(\bar P\) has degree six, there are at most six double roots
over the algebraic closure. In particular \(F_0\) has positive
degree. The integral double cover of the \(t\)-line, including
its infinity chart, has geometrically irreducible nodal special
fiber with normalization \(C_0:y^2=F_0\). There are two distinct
smooth points at infinity, because the leading coefficient
of \(\bar F\) is a nonzero square. The special fiber is reduced.
Each singularity is an ordinary node, and its local smoothing
has thickness \(n\ge1\) since the generic fiber is smooth.

For an irreducible factor \(h\) of \(g\) with degree \(d\), put
\(e_h=1\) if \(F_0\bmod h\) is a square in \(\mathbb F_{5^d}\)
and \(e_h=-1\) otherwise. This is the sign on the two branches
after \(d\) Frobenius steps. Let \(T\) be the torus in the
identity component of the special Néron fiber and let \(\Phi\)
be its component group. The contribution of this node orbit to
\(\#T(\mathbb F_5)\) is \(5^d-e_h\), and to
\(\dim T(\mathbb F_5)[3]\) it is one exactly when
\(5^d\equiv e_h\bmod3\).
For \(\Phi\), the geometric graph has one loop of length \(n\)
at each geometric node. The pairing on the loops is diagonal
with these lengths. Therefore the component group is a direct
sum of the \(\mathbb Z/n\) of the geometric nodes, with the
same signed permutation. An orbit contributes at most one to
\(\dim\Phi(\mathbb F_5)[3]\), and only if \(e_h=1\) and
its thickness is divisible by three. The identity component
sequence and the Frobenius-equivariant component group formula are
[DDMM, equations (2.3)–(2.7), reduction (2.18), and Lemmas 2.21–2.22](sources/dokchitser-dokchitser-maistret-morgan-local-hyperelliptic.pdf).
In the present model the graph pairing follows by resolving each
ordinary double point, which subdivides its loop into \(n\) edges.

The formal reduction kernel has no \(3\)-torsion because
the residue characteristic is five. Also the connected
group is an extension of \(\operatorname{Jac}(C_0)\) by
\(T\). We therefore have the initial upper bound
\[
 \dim J(\mathbb Q_5)[3]\le
 b_0 + b_T + b_\Phi,\qquad
 \begin{cases}
 b_T=\#\{h:5^{\deg h}\equiv e_h\bmod3\},\\
 b_\Phi=\#\{h:e_h=1\},
 \end{cases}
 \tag{1}
\]
where \(b_0\) is the minimum of the \(3\)-valuation of
\(\#\operatorname{Jac}(C_0)(\mathbb F_5)\) and the
multiplicity of \(X-1\) in its Frobenius polynomial modulo
three. This bound permits every split thickness to be divisible
by three and permits every rational normalization class of
order three to lift to the torus extension.

## Finite elimination and the remaining gluing condition

An affine substitution makes \(\bar P\) monic and kills the
degree-six coefficient over \(\mathbb F_5\), exactly as in SP2.
Enumerate \(t^7+\sum_{j=0}^5 c_jt^j\) up to
\(\bar P(t)\mapsto-\bar P(-t)\). The [checker](verify_septic_pell_nodal_five.py)
covers all \(7875\) representatives. There are \(616\)
representatives with a root of multiplicity greater than two.
Of the remaining \(7259\), exactly six have the bound (1)
at least five, and in all six the bound is five.
As a control, the computed order of the connected reduction is
divisible by seven throughout the nodal and smooth census.
Indeed the Pell function on the integral special fiber has
divisor seven times the infinity difference. This difference
cannot be principal: a function with one simple pole at a
smooth infinity would give a degree-one finite morphism from
that integral curve to \(\mathbb P^1\), forcing an isomorphism.

The six cases are in the table. Split orbits are indicated by
their degrees. In the last column for five cases, the point
\(Q\) gives a \(3\)-primary class of order at least nine in the
Jacobian of \(C_0\): take \([Q-\infty_+]\) and multiply by
the prime-to-three part of the indicated Jacobian order.
It follows that the actual dimension of its \(3\)-torsion is at
most \(v_3(\#\operatorname{Jac}(C_0)(\mathbb F_5))-1\).
In each of these five cases this improves (1) to at most four.

| Ascending coefficients \((c_0,\ldots,c_5)\) | Split orbit degrees | \(\#\operatorname{Jac}(C_0)(\mathbb F_5)\) | \(b_0\) | \(Q\) and order after projection |
| --- | --- | ---: | ---: | --- |
| \((0,0,1,4,3,2)\) | \(1\) | 4536 | 4 | \((0,1)\), 81 |
| \((0,2,0,0,0,1)\) | \(2,2\) | 21 | 1 | Gluing calculation below |
| \((0,4,2,2,0,2)\) | \(2\) | 756 | 3 | \((0,1)\), 9 |
| \((1,0,0,4,2,0)\) | \(1\) | 4536 | 4 | \((2,2)\), 27 |
| \((2,2,1,1,4,1)\) | \(1,1\) | 1512 | 3 | \((1,2)\), 27 |
| \((2,2,1,2,1,4)\) | \(1,2\) | 252 | 2 | \((2,1)\), 9 |

In the remaining row,
\[
 \begin{split}
 \bar P&=t^7+t^5+2t,\\
 g&=(t^2+2t+4)(t^2+3t+4),\\
 F_0&=(t^3+t^2+2)(t^3+4t^2+3).
 \end{split}
\]
The two irreducible quadratic factors give split nodes over
\(\mathbb F_{25}\). Put \(D_0=(0,1)-\infty_+\) on \(C_0\).
The class of \(7D_0\) generates its rational \(3\)-torsion.
There is a function \(f=A+B y\) of divisor \(21D_0\).
Its certificate is
\[
 \begin{split}
 A&=4t^{21}+3t^{20}+4t^{19}+3t^{18}+t^{16}+4t^{15}
     +t^{14}+4t^{13}+2t^{12}+t^{11}+3t^{10}+t^8\\
   &\qquad+2t^7+2t^6+4t^4+t^3+2t+1,\\
 B&=4t^{18}+3t^{17}+t^{16}+2t^{15}+4t^{14}+4t^{13}
     +t^{12}+3t^{11}+2t^{10}+t^8+4t^7\\
   &\qquad+2t^4+3t^2+3t+4,\\
 A^2-F_0 B^2&=2t^{21}.
 \end{split}
\]
The nonzero branch value at \((0,-1)\) and the pole of
order 21 at \(\infty_+\) specify the divisor from the norm
identity.

The divisor \(7D_0\) is disjoint from the nodes, so it lifts
to a divisor class \(E\) on the nodal curve. In the torus
coordinates of \(3E\), the first coordinate is the ratio
of the values of \(f\) at the two branches over
\(h=t^2+2t+4\), up to inversion. Indeed the torus in the
normalization sequence records precisely the gluing scalars
between the two fibers of a line bundle. Reducing \(f\) on
these branches yields the ratio \(a+2\) in
\(\mathbb F_{25}=\mathbb F_5[a]/(a^2+2a+4)\).
The exact identity \((a+2)^8=3a\ne1\) shows this ratio is
not a cube. Changing the lift \(E\) by a torus class changes
\(3E\) by a cube there. Thus the normalization's order-three
line has no order-three lift in the connected reduction.
The connected \(3\)-torsion has dimension two from the two
split quadratic torus orbits. The component group adds at
most two. This proves the bound in the last row and finishes
the census.

## Scope and replay

Run the [checker](verify_septic_pell_nodal_five.py) with SageMath's
Python. It checks the complete residue class census using exact
finite-field polynomial and Frobenius arithmetic, the five
explicit Jacobian order witnesses, and the last node's gluing
certificate. The passage from the residue group to the local
Jacobian uses the written semistable reduction argument.
SN1 concerns characteristic-zero models with a unit leading
coefficient at five and only smooth or ordinary nodal reduction
in the given equation. SN2 below treats the deeper collisions.

## All multiplicities at five

**Theorem SN2.** Let \(P\in\mathbb Z_5[t]\) have degree seven
with unit leading coefficient and \(P^2+1\) squarefree over
\(\mathbb Q_5\). Then
\[
 \dim_{\mathbb F_3}\operatorname{Jac}(y^2=P^2+1)(\mathbb Q_5)[3]\le4
\]
with no restriction on the multiplicities of its reduction.
In particular, having a geometric multiple root modulo five
cannot repair the septic rank-five obstruction as long as
the model stays integral and degree-preserving at five.

### A rank budget over a field of semistable reduction

For the reduction of the integral equation write
\(\bar P^2+1=G^2F_0\), with \(F_0\) squarefree and \(G\) monic.
Choose a finite Galois extension over which all branch roots
split and the curve is semistable. Such an extension exists
by the semistability criterion in
[DDMM, Theorem 1.8](sources/dokchitser-dokchitser-maistret-morgan-local-hyperelliptic.pdf):
after splitting the roots a further even ramification factor
makes all the depth and parity conditions hold. Work on the
geometric special fiber of its minimal regular model and retain
the action of an arithmetic Frobenius \(\sigma\) for \(\mathbb Q_5\).
There is no assumption that the semistable extension has
residue field \(\mathbb F_5\).

The outer normalization component is
\(C_0:y^2=F_0\), with its usual \(\mathbb F_5\)-Frobenius.
Indeed the full branch cluster has depth zero, its children
are the distinct root residues, and \(\bar P-2,\bar P+2\) are
coprime of odd degree seven. Each contributes an odd child and
there are at least three children altogether. If there were just
two, both would have multiplicity seven, contradicting the degree
six derivative. Thus the outer cluster is principal and has
one normalization component. Removing the pairs in the child
clusters is exactly the equation \(F_0\).
A further even ramification factor, if needed, subdivides all self-loops\nso that the regular model has smooth components.\nA root of multiplicity two gives a separate loop based on that
component. Its sign under the Frobenius of its residue field
is the branch sign used in (1). These assertions follow from
[DDMM, Theorems 8.5 and 8.7](sources/dokchitser-dokchitser-maistret-morgan-local-hyperelliptic.pdf),
which describe the normalization, chains and action of the
original Galois group on the semistable fiber.

Let \(r_2\) be the number of geometric roots of multiplicity
exactly two, and \(g_0\) the genus of \(C_0\). Use the quantities
\(b_0,b_T,b_\Phi\) from (1) on \(C_0\) and the double-root loops.
The remaining components and graph contribute at most
\(6-g_0-r_2\) to the rational torsion dimension. To check this,
the genus formula for the semistable fiber expresses
\(6-g_0-r_2\) as the sum of the genera of its remaining
normalizations and the rank of its remaining graph homology.
For an orbit of length \(d\) of normalization components of
genus \(g'\), the invariant torsion dimension is bounded by
\(dg'\). When \(d\) is odd, the Weil pairing shows that the
\(\sigma^d\)-fixed subspace in one component is isotropic,
because \(\sigma^d\) acts on \(\mu_3\) by \(5^d=-1\bmod3\).
It has dimension at most \(g'\). When \(d\) is even the bound
\(2g'\le dg'\) suffices.

For the remaining graph let \(H\) be its integral homology
lattice. The double-root loops are edge-disjoint from this
part, so the length pairing and component-group calculation
split off their terms. The torus \(3\)-torsion is
\(\operatorname{Hom}(H,\mu_3)\); its fixed dimension is the
multiplicity of the eigenspace of eigenvalue \(-1\) on \(H/3H\).
The monodromy pairing \(H\hookrightarrow H^\vee\) gives a
Frobenius-equivariant injection \(\Phi[3]\hookrightarrow H/3H\)
by tensoring its cokernel sequence. Its fixed dimension is
therefore at most the dimension of the \(+1\) eigenspace.
The sum of those dimensions is at most \(\operatorname{rank} H\),
even when the action is not semisimple. Hence the total budget is
\[
 b_0+b_T+b_\Phi+(6-g_0-r_2).
 \tag{2}
\]
The specialization of a rational \(3\)-torsion point is
\(\sigma\)-fixed, and the formal kernel still has no such torsion.

### Divisibility forced by a split component class

In a regular semistable model, suppose a rational split
double-root loop of length \(n\) is based at two points \(P_0,Q_0\)
of the outer normalization and \(3\mid n\).
Retain any subset of the other double-root node orbits on that
component, gluing those pairs of branches to form a partial nodal
curve \(Y\). If a rational \(3\)-torsion point maps to an order-three
generator on the chosen loop and to zero in all other graph
components, then
\[
 [P_0-Q_0]\in 3\operatorname{Pic}^0(Y)(\mathbb F_5)
 \tag{3}
\]
up to replacing the difference by its negative. Keeping no
other nodes gives the same necessary condition on \(C_0\).

For the proof, a rational point at infinity makes a rational
Jacobian class representable by a line bundle. Extend that
bundle to the regular model. Its component-group image is
represented by the multidegree vector modulo the graph
Laplacian by [Raynaud’s description in Baker, Appendix A, (A.3)–(A.5)](sources/baker-raynaud-specialization-pp641-644.pdf).
For a pure generator on the chosen loop, adjust
by a vertical divisor so the multidegree is \(+1\) at the
vertex \(n/3\) along the loop from \(P_0\), \(-1\) on the
outer component, and zero elsewhere. This is the divisor
representative of \(n/3\) in the length-\(n\) cycle, using
the identification of the cycle component group with the
sum of the vertex indices times the multidegrees modulo \(n\).
Such an adjustment is \(\sigma\)-invariant. The degree
correction is invariant, and its solution in the lattice of
vertical divisors is unique up to an integer multiple of
the whole fiber. Setting the outer coefficient to zero
makes the solution invariant.

Cubing the adjusted bundle gives a vertical divisor, which
can also be normalized to have outer coefficient zero.
Its coefficients have Laplacian equal to three times the
multidegree, with the intersection convention
\(a_{j-1}-2a_j+a_{j+1}\) on the loop.
On the chosen interval the slopes are \(-2\) and \(1\);
the two neighboring coefficients at the outer vertex are
\(-2\) and \(-1\). On the rest of the graph the solution
is harmonic and zero, by connectedness and the outer value
zero. Restriction to the partial curve along the retained
loops thus gives the line bundle of \(-2P_0-Q_0\).
For this restriction one can first keep the rational chains
of the retained loops and then contract them. The restricted
bundle has degree zero on each chain and descends, with the
gluing scalar given by the identifications along the chain.
After tensoring the original restriction by
\(\mathcal O(P_0)\), its cube is \(\mathcal O(P_0-Q_0)\).
The entire restriction is \(\sigma\)-invariant and gives a
rational point of the partial Picard group. This proves (3).

### Seven remaining charts

Among the \(616\) normalized representatives left out of
SN1, (2) is at most four except for the following seven.
The arithmetic here is over \(\mathbb F_5\). The Jacobian
orders in this section refer to the outer normalization.

| Coefficients \((c_0,\ldots,c_5)\) | Largest multiplicity | Jacobian order | Budget (2) | Final test |
| --- | ---: | ---: | ---: | --- |
| \((0,1,0,1,1,4)\) | 4 | 1134 | 6 | Point of \(3\)-primary order 81 |
| \((1,2,4,2,3,0)\) | 4 | 126 | 5 | Split component fails (3) on \(C_0\) |
| \((2,0,0,0,0,0)\) | 7 | 126 | 5 | Point of \(3\)-primary order 9 |
| \((2,0,0,1,0,2)\) | 3 | 21 | 5 | Connected normalization class fails gluing |
| \((2,0,1,0,0,4)\) | 6 | 42 | 5 | Split component fails (3) on \(C_0\) |
| \((2,0,3,3,2,2)\) | 3 | 7 | 5 | Split component fails (3) with a quadratic orbit retained |
| \((2,4,3,4,4,4)\) | 4 | 126 | 5 | Point of \(3\)-primary order 9 |

For the first, third and last rows, the points are respectively
\((0,1),(2,2),(1,2)\), each taken relative to \(\infty_+\)
and projected to the \(3\)-primary group as in SN1.
These exact orders lower the corresponding rank budgets to
at most four.

For the second row, the outer Jacobian has order \(126\)
and \(3\)-primary group \(C_3^2\). Two independent projected
points verify this. The rational split double root is at \(t=1\)
and the difference of its branches with \(y=\pm2\) has
nonzero \(3\)-primary projection. Thus it is not three
times a rational point. The same argument works in the fifth
row, with outer order \(42\), branches at \(t=0\), \(y=\pm2\),
and nonzero projection to its group of order three.
In each case reaching the entire budget would require
surjectivity from rational \(3\)-torsion onto the fixed
component \(3\)-torsion, including the pure generator on
this split loop. This contradicts (3). The implication to
surjectivity follows because (2) bounds the kernel and the
image separately. Equality requires both bounds and the
total fixed component dimension to be attained.

For the fourth row, the triple cluster at \(t=0\) attaches
to the outer component at a separating node. Thus its
connected Picard group contribution and that of the outer
partial nodal curve are separate factors.
Here
\[
 F_0=t(t^5+3t^4+2t^2+2t+4),\qquad h=t^2+t+2 .
\]
The degree-two point with \(h(t)=0,\ y=t+4\), minus
\(2\infty_+\), generates the rational \(3\)-torsion on
the outer normalization. A function of three times that
divisor is \(A+B y\) with
\[
 A=t^6+4t^4+3t^2+1,\qquad B=t^3+t^2+3t+3,\qquad
 A^2-F_0B^2=2h^3.
\]
Only the indicated branch vanishes over \(h\), and the
leading terms give the prescribed infinity. At the nonsplit
double root \(t=2\) the two values of the function are
\(1\pm\sqrt{2}\). Their ratio has square \(2+2\sqrt{2}\ne1\)
in the norm-one group of order six, hence is not a cube.
This excludes lifting the normalization's order-three line
to the connected reduction of the outer nodal piece.
Its bound drops by one, as in SN1.

For the sixth row the outer normalization is
\[
 y^2=t^4+3t^3+3t+3 .
\]
Let \(P_0,Q_0\) be the two branches at \(t=2\), \(y=\pm2\).
The split orbit over \(t^2+4t+1=0\) has degree two.
Retain that orbit in (3). The following function has
divisor \(7(Q_0-P_0)\):
\[
 \frac{A+B y}{(t-2)^7},\qquad
 A=3t^7+3t^6+3t^5+4t^4+4t^2+4t+1,\qquad
 B=t^4+t^3+t^2.
\]
Its numerator norm is \(4(t-2)^{14}\), and its
numerator at \(P_0\) is nonzero. Over the retained
quadratic node, with \(a^2+4a+1=0\) and branch
\(y=2a+2\), the branch ratio is \(2a\).
Since \((2a)^8=a+4\ne1\), it is noncubic.
The outer Jacobian group has order seven. If
\([P_0-Q_0]\) were three times a rational class on the
partial curve, multiplying that class by seven would give
a torus class, and the calculated ratio for \(7(P_0-Q_0)\)
would be a cube. This is impossible. The necessary lift
of the pure split component generator therefore fails also
in this row. Together with SN1 this covers every reduction.

The checker includes the full multiplicity census, the
additional point orders and split differences, and both
further gluing identities. The finite check only uses the
outer reductions. The rank budget and necessary divisibility
condition in the written argument handle all possible inner
depths and Frobenius actions, including the wild multiplicity
six and seven cases. This still leaves rational septic
models with no integral degree-preserving chart at five,
and leaves higher degrees outside the statement.

## Rational models outside the integral chart

**Proposition SN3.** For every \(a\ge1\) there are infinitely
many rational septic polynomials \(P\) with \(P^2+1\)
squarefree such that
\[
 \operatorname{Jac}(y^2=P^2+1)(\mathbb Q_5)[3^a]
             \simeq(\mathbb Z/3^a)^6,\qquad
 \operatorname{Jac}(y^2=P^2+1)(\mathbb Q)[3]=0 .
\]
None of these polynomials can be put in an integral
degree-preserving chart by a rational affine change
of \(t\) over \(\mathbb Q_5\).

Choose a positive integer \(k\equiv3\bmod5\) divisible by
\(3^a\). Put
\[
 Q(t)=t(t-1)(t-2)(t-3)(t-4)(t-5^k)(t-1-5^k),
 \qquad P(t)=5^{-kN}Q(t)
\]
with \(N\ge3,\ N\equiv4\bmod5\). These are distinct
polynomials as \(N\) varies. Over \(\mathbb Q_5\) the branch
equations are \(Q(t)=\pm i5^{kN}\), where \(i^2=-1\).
At each center \(c\) in the list of seven roots of \(Q\),
the derivative valuation \(s_c\) is \(k\) at the four
centers congruent to \(0\) or \(1\bmod5\), and zero at the
three others. Since \(kN>2s_c\), Hensel lifting puts one
root of each sign at distance of valuation \(kN-s_c\)
from \(c\). The leading Hensel corrections have opposite
signs, so the two roots have mutual distance of that same
valuation. Between different center discs the branch
distances have valuations \(k\) or zero according to the
center differences. This gives fourteen distinct rational
branch roots over \(\mathbb Q_5\).

Use the equation \(Y^2=Q(t)^2+5^{2kN}\).
Its outer branch cluster has depth zero. There are two
clusters of size four and depth \(k\), each consisting
of two twins, together with the three outer twins.
The principal clusters are all even with only even
children. The parity condition for semistability holds
since the leading coefficient is one and the sum of
depths at a size-four cluster is \(4k\).
The two normalizations of each principal cluster are
rational and individually defined over \(\mathbb F_5\).
At the outer cluster the reduction is the square
of \(Q\). At each depth-\(k\) cluster use \(t=c+5^k x\)
and divide \(Y\) by \(5^{2k}\); again the reduction is a
square, since \(2kN-4k>0\).
The ends of each twin chain lie on the two distinguished
sheets of its parent. Frobenius fixes these sheets and
every center, so the graph and all its edges are split.

The graph has two outer vertices, two vertices for each
size-four cluster, four chains of length \(k\) from the
outer vertices to those four vertices, and seven twin
chains. The three outer twin lengths are \(2kN\) and
the four inner twin lengths are \(2k(N-2)\).
Thus its homology has rank \(11-6+1=6\) and every
entry of the length pairing is divisible by \(3^a\).
The component group has trivial Frobenius action and
its \(3^a\)-torsion is therefore \((\mathbb Z/3^a)^6\).
The identity component contributes a split torus with
\(\mathbb F_5\)-points of order \(4^6\). Its preimage under
reduction is uniquely \(3^a\)-divisible, since the
formal kernel is pro-five. Smooth lifting and the
finite-field torus torsor property give surjectivity
of \(J(\mathbb Q_5)\) onto the rational component group.
Every \(3^a\)-torsion component class then has a unique
\(3^a\)-torsion lift. This proves the asserted local group.

At eleven the residue data are identical for all the
allowed \(k,N\), since five has order five in
\(\mathbb F_{11}^\times\). The equation
\(Y^2=Q^2+5^{2kN}\) has smooth degree-preserving reduction
there, and the reduced Jacobian has order \(2643158\),
prime to three. Thus the rational group has no
three-torsion by specialization at eleven.
Finally an integral degree-preserving chart over
\(\mathbb Q_5\) would contradict SN2. This construction
shows that a local obstruction at five alone cannot
exclude the remaining rational-polynomial route.

## Polynomial norm sections when the denominator is large

**Proposition SN4.** Let \(K\) be a discretely valued local
field of residue characteristic prime to six whose residue
field \(k\) has no nontrivial cube roots of unity. Let
\(P\in K[t]\) have degree \(m\ge2\), with \(3\nmid m\)
and \(P^2+1\) squarefree. Suppose its Gauss valuation is
\(-u<0\), that \(Q=\pi^uP\) has unit leading coefficient,
and that \(\bar Q\) splits over \(k\). Write \(e_r\) for
the multiplicity of each of its roots. If
\[
 \sum_{r}(-e_r\bmod3)>m-1
 \tag{4}
\]
where each summand is in \(\{0,1,2\}\), then every solution
\[
 2PZ+Z^2-1=H^3,\qquad Z\ne0,\qquad\deg Z<m
\]
has \(v(Z)=-u,\ 3\mid u\), and its primitive reductions obey
\[
 \bar H_0^3=Z_0(2\bar Q+Z_0),\qquad
 Z_0=\overline{\pi^u Z},\qquad
 \bar H_0=\overline{\pi^{2u/3}H}.
 \tag{5}
\]
Here \(v\) on each polynomial is its Gauss valuation.
In particular, (5) is necessary in this nonintegral chart,
even if the local Jacobian has maximal rational \(3\)-torsion.

Put \(z=v(Z)\). The degree identity
\(3\deg H=m+\deg Z\) gives \(\deg H<m\).
For \(-u<z<u\) the product term \(2PZ\) dominates.
The reduced cube is a nonzero scalar times \(\bar Q\)
times the primitive reduction of \(Z\). Each root
therefore forces its copy of the deficit in (4).
This contradicts the degree bound.
If \(z=u\), reducing gives \(\bar H^3+1=2\bar Q Z_0'\)
where \(Z_0'=\overline{\pi^{-u}Z}\ne0\).
Every root of \(\bar Q\) forces \(\bar H=-1\).
The other factor \(\bar H^2-\bar H+1\) takes the value
three there, so \(\bar Q\) divides \(\bar H+1\)
with full multiplicities. This again contradicts the degree
bound. If \(z>u\), then \(\bar H=-1\) identically.
Writing \(H^3+1=(H+1)(H^2-H+1)\), the latter factor
reduces to three, and the primitive reduction of \(H+1\)
must be a nonzero scalar times \(\bar Q\) times the primitive
reduction of \(Z\). Their degrees are incompatible.

If \(z<-u\), the leading relation is a cube equal to
the square of the primitive reduction \(Z_*\) of \(Z\).
Thus \(\deg Z_*\equiv0\bmod3\). On the curve,
\[
 \alpha=1-Z(P+y),\qquad y^2=P^2+1
\]
has divisor divisible by three by the argument of
[SP1](SEPTIC_PELL_FIVE_DIRECTION_CONSTRUCTION.md).
On the integral model \(Y^2=Q^2+\pi^{2u}\), the
component \(Y=-\bar Q\) has the primitive reduction
\(\alpha_*=Z_* /(2\bar Q)\) for this function, since
\(P+y=-\pi^u/(Q-Y)\). At the smooth infinity on that
component, \(\alpha_*\) has order \(m-\deg Z_*\),
prime to three. But the order of the primitive reduction
at any smooth point of the special fiber is divisible
by three when the horizontal divisor is. Indeed the
local ring is regular there, and removing the order
on the special fiber leaves a horizontal divisor
all of whose multiplicities are divisible by three.
Its intersection with the special fiber has the
same property. This excludes \(z<-u\).
For \(z=-u\), the leading polynomial \(\bar Z(2\bar Q+\bar Z)\)
is nonzero because \(\deg \bar Q>\deg\bar Z\).
Its valuation before normalization is \(-2u\).
The equality to \(H^3\) gives (5).

**Theorem SN5.** For \(K=\mathbb Q_5\) and \(m=7\),
suppose \(\bar Q=\overline{5^uP}\) as in SN4 has
all five rational residues as roots. If two residues
are double, no nontrivial polynomial norm section
of degree \(\deg Z<7\) exists. If one residue is triple,
the span in \(J(\mathbb Q_5)[3]\) of the classes from
all these sections has dimension at most three.
In particular the SN3 family has no such sections
already over \(\mathbb Q_5\), regardless of the separate
global test at eleven.

In either case the deficit in (4) is eight. The affine group
over \(\mathbb F_5\) is transitive on pairs of distinct roots
and on single roots. After translation and scaling of the
coordinates, write \(W=t^5-t\). The respective cases of
\(\bar Q\) are \(t(t-1)W\) and \(t^2W\), up to a
nonzero scalar. A scalar can be divided out of both
\(\bar Q\) and \(Z_0\), since every nonzero scalar of
\(\mathbb F_5\) is a cube. In (5) the product has degree
\(7+\deg Z_0\), so \(\deg Z_0=2\) or \(5\).
The checker enumerates these two degrees with every nonzero
leading coefficient. For \(t(t-1)W\) no product is a cube.
For \(t^2W\) there are exactly two possibilities:
\[
 Z_0=t^5+4t^3,\qquad Z_0=4t^5+4t^3.
\]
At the residues \(0,1,2,3,4\), the orders of
\((2\bar Q+Z_0)/(2\bar Q)\) modulo three are respectively
\((0,1,2,2,1)\) and its negative.

To translate these residues into a rank bound, work over
a field of semistable reduction retaining the original
arithmetic Frobenius, as in SN2. Since the initial reduction
is \(Y^2=\bar Q^2\), there are two rational outer
normalization components. Four twin chains join them at
the simple root residues. The six-root cluster over the
triple root supplies one more connecting piece. It has
a link from each of the two outer components. Its internal
graph is connected between those links: it either has a
single outer normalization component of its own, or its
two sheets are joined through descendants, ending at a
twin. This is the cluster description of
[DDMM, Theorems 8.5 and 8.7](sources/dokchitser-dokchitser-maistret-morgan-local-hyperelliptic.pdf).
Let \(H_{\rm in}\) be the graph homology inside that
connecting piece. Replacing it by one path gives a graph
of two vertices with five paths joining them, whose homology
\(U\) consists of five integer flows of sum zero.
There is a Frobenius-equivariant sequence
\(0\to H_{\rm in}\to H_1\to U\to0\), and Frobenius acts
trivially on \(U\). The sum of the rank of \(H_{\rm in}\)
and the genera of the nonouter normalizations is \(6-4=2\).

For a three-torsion divisor \(D\) with principal function
\(\alpha\), the image of its component class under
\(\Phi[3]\hookrightarrow H_1/3H_1\) can be read from the
orders of the primitive reductions of \(\alpha\) at the
attaching points. Here is the integral graph calculation.
Let \(a_v\) be the coefficients of the vertical divisor
\(3\bar D-\operatorname{div}(\alpha)\). Orient the graph
edges, take differences \(da\) of these coefficients
along them, and let \(c\) be an integral chain with
boundary the multidegree of \(\bar D\). Then
\(\partial(da)=-3\partial c\) and \(da+3c\) is a
cycle. Its length pairing against integral cycles is
three times the pairing of \(c\). This is exactly the
injection of the three-torsion of the component group
into \(H_1/3H_1\), up to the common sign convention.
Locally at a node, the order of the primitive reduction
of \(\alpha\) on one component is the adjacent vertical
valuation difference plus the horizontal contribution.
The latter is a multiple of three. Thus the differences
along the five paths at their outer ends are precisely
the five residue orders modulo three, with orientation
signs uniformly chosen on one outer sheet. On the sheet
\(Y=-\bar Q\) that reduction is
\((2\bar Q+Z_0)/(2\bar Q)\). Hence the image of all the
marked classes in \(U/3U\) has dimension at most one.

The fixed part of the torus \(3\)-torsion, from the dual
lattice twisted by \(\mu_3\), comes from the \(-1\)
eigenspace of \(H_{\rm in}/3H_{\rm in}\) at most.
Indeed \(U\) has trivial action, so its twisted dual
has no fixed vectors. The component \(3\)-torsion in the
kernel of projection to \(U/3U\) comes from the \(+1\)
eigenspace at most. The sum, with the contributions of
the normalization groups bounded as in SN2, gives
kernel dimension at most two. Thus the total span of the
marked section classes has dimension at most \(1+2=3\).
The proof leaves other residue-root profiles, including
reductions of the primitive \(Q\) missing one or more
\(\mathbb F_5\) residues, outside this bound.

## The fully split residue profiles with fewer distinct roots

**Corollary SN6.** Suppose \(P\in\mathbb Q_5[t]\) is a
septic with \(P^2+1\) squarefree, \(v(P)=-u<0\) and
\(Q=5^uP\) has unit leading coefficient. Suppose \(\bar Q\)
splits over \(\mathbb F_5\) into at least three distinct
linear factors. If polynomial norm sections as in SN5
span a five-dimensional space of local Jacobian
three-torsion, then, up to affine change over
\(\mathbb F_5\), the multiplicities of the roots of
\(\bar Q\) at residues \(0,1,2,3,4\) are one of
\[
 (0,0,1,2,4),\quad (0,0,1,5,1),\quad (0,1,1,3,2).
 \tag{6}
\]
Every nontrivial section in those three profiles is
quintic in \(Z\). The three profiles are necessary
possibilities from this test.

Let \(s\) be the number of distinct roots. On the outer
minus sheet, any function \(\alpha\) from a row with
\(v(Z)>-u\) reduces to \(1\). It therefore has zero
outer flow regardless of the deficit (4).
The smooth infinity argument in SN4 excludes
\(v(Z)<-u\) independently of that deficit. The only
remaining sources of outer flow obey (5).
There are two outer components, and each root residue
gives a connected piece with one link to each outer
component. The justification is the one in SN5 for
an even-sized cluster. Since \(s\ge3\), the largest such
cluster has size at most ten and the outer cluster
is principal. Thus the outer flow quotient has rank
\(s-1\) and trivial Frobenius action. Repeating the
kernel bound of SN5 gives
\[
 \dim\langle\text{marked section classes}\rangle
       \le b+7-s,
\]
where \(b\) is the dimension of the span of the
residue-order vectors from (5).

There are sixteen multiplicity tuples up to the affine
action, taking the lexicographically least tuple in each
orbit. Only the following five yield any critical rows.
The checker covers both allowed residue degrees of \(Z_0\).

| Multiplicity tuple | Number of \(Z_0\) | \(b\) | Bound \(b+7-s\) |
| --- | ---: | ---: | ---: |
| \((0,0,1,2,4)\) | 4 | 2 | 6 |
| \((0,0,1,5,1)\) | 4 | 2 | 6 |
| \((0,1,1,1,4)\) | 2 | 1 | 4 |
| \((0,1,1,3,2)\) | 5 | 2 | 5 |
| \((1,1,1,1,3)\) | 2 | 1 | 3 |

For the other eleven profiles \(b=0\) and \(7-s\le4\).
In each of (6) all critical rows have residue degree five.
If \(-u<v(Z)<u\), the cube in the dominant product
forces total root deficit five in each of these three
profiles, so \(\deg Z\ge5\) again. The arguments in SN4
for \(v(Z)\ge u\) use only that \(ar Q\) splits
and still exclude those cases. This gives the degree
assertion for every row in (6).

The remaining cases for this approach are (6), profiles
with fewer than three distinct roots, profiles with
irreducible residue factors of higher degree, and reductions
in which the primitive polynomial loses degree.
