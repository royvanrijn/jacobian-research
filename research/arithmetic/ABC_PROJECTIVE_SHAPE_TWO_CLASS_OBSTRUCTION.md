# Projective-shape completion with two anchors and its parity obstruction

The two split covers for an affine shape use the marked point
at infinity of the shape. A fixed projective shape of four
actual primitive rows admits a power-sized family on which
all within-role ratios of normalized summands are squares.
In particular, the two parity classes of the fourfold
products are separately trivial. The obstruction does not
require cancellation between these two classes.
[Affine-shape completion (AS1–AS4)](ABC_AFFINE_SHAPE_JOINT_POWER_CLASSES.md#two-power-class-records-determine-an-affine-shape)

## Positive primitive quadruples from odd multiples

**PSO1.** There are constants \(\kappa,c>0\) such that for
all sufficiently large \(X\) there are at least \(cX^\kappa\)
ordered quadruples of distinct positive primitive rows
\(A_i+B_i=C_i\le X\) with the following properties. The
four ratios \(b_i=B_i/C_i\), in their given order, have
one fixed projective shape. For \(a_i=A_i/C_i=1-b_i\),
every \(a_i/a_j\) and every \(b_i/b_j\) is a rational
square. Thus for any four integers \(e_i\) of even sum,
both \(\prod_i a_i^{e_i}\) and \(\prod_i b_i^{e_i}\)
are rational squares.

Put \(f(t)=t(t^2-1)\) and work first over
\(\mathbb Q(t)(\sqrt{f(t)})\). On \(E:y^2=f(x)\) let
\(P=(t,\sqrt{f(t)})\) and put \(z_j=x(jP)\). The quadratic
involution sends \(P\) to \(-P\), so \(z_j\in\mathbb Q(t)\).
For fixed \(j\ge1\) and \(t\to+\infty\), with the positive
square root chosen at \(P\),
\[
                 x(jP)\sim t/j^2,\qquad y(jP)\sim t^{3/2}/j^3 .
                                                               \tag{PSO1a}
\]
The tangent slope at \(P\) is asymptotic to
\(3\sqrt{t}/2\). For addition of \(P\) to \(jP\) with \(j>1\),
the induction hypotheses give the slope divided by \(\sqrt t\)
tending to
\[
       \frac{j^{-3}-1}{j^{-2}-1}=\frac{j^2+j+1}{j(j+1)} .
\]
Treating the tangent constant \(3/2\) as the same formula at
\(j=1\), the chord rule gives the next leading \(x\)-coefficient as
\(\bigl((j^2+j+1)/(j(j+1))\bigr)^2-1-j^{-2}=(j+1)^{-2}\).
For the next \(y\)-coefficient use the line through \(P\) and
negate its value at \(x((j+1)P)\). The result is
\((j^2+j+1)(1-(j+1)^{-2})/(j(j+1))-1=(j+1)^{-3}\).
In particular all multiples used below and the chord or tangent
formulas are defined for sufficiently large \(t>0\).

Write \(D=f(t)\). After the change
\((x,y)\mapsto(Dx,D\sqrt D\,y)\) to
\[
                         Y^2=X(X-D)(X+D),
\]
the point \(P\) becomes the rational point \((Dt,D^2)\)
over \(\mathbb Q(t)\). Fix \(r\in\{-1,0,1\}\).
If a nonvertical chord or tangent \(Y=\ell(X)\) meets
the twist at three abscissas \(X_h\), with multiplicities,
evaluating its monic intersection cubic at \(Dr\) gives
\[
                         \prod_{h=1}^3 (X_h-Dr)=\ell(Dr)^2 .
\]
The three factors are nonzero for the addition steps
used to construct the first five multiples, by (PSO1a).
The reflected third point has the same abscissa. Therefore
on those multiples the square class of \(X-Dr\) at the
\(j\)-th multiple is the \(j\)-th power of the class of
\(D(t-r)\). For odd \(j\) this gives
\[
                   (z_j-r)/(t-r)\in\mathbb Q(t)^{\times2}.
                                                               \tag{PSO1b}
\]
This is the coordinate calculation behind the usual full
two-torsion descent.
[Silverman, The Arithmetic of Elliptic Curves, Proposition X.1.4](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf#page=328)

Take the four ordered marks \(r_i=-1,0,1,\infty\), and
write
\[
 k=\frac{z_5-z_3}{z_5-z_1},\qquad
 b(r)=k\frac{r-z_1}{r-z_3},\qquad
 a(r)=(1-k)\frac{r-z_5}{r-z_3}=1-b(r).
\]
At \(r=\infty\) use \(b(\infty)=k\), \(a(\infty)=1-k\).
The map \(b\) is a nonsingular projective transformation,
since \(z_1,z_3,z_5\) are distinct. Its four values therefore
have the projective shape of the fixed marks. By (PSO1a),
\(k\to2/27\) and \(b(r)\to2/3\) at each finite mark.
Consequently all four values are in \((0,1)\) for large
positive \(t\). By (PSO1b) at \(j=1,3,5\) the values
of \(b(r)/k\) and \(a(r)/(1-k)\) are squares for each
finite mark. They equal one at infinity. Thus all the
asserted square-class identities hold and specialize
away from a finite set of poles.

For a sufficiently large integer parameter \(t\), write
each \(b(r_i)\) in lowest terms as \(B_i/C_i\) and set
\(A_i=C_i-B_i\). These are positive primitive rows and
they are pairwise distinct because \(b\) is one-to-one
on the marks. Every ratio is a fixed rational function
of \(t\) of finite degree. Hence each \(C_i\) is
\(O(t^J)\) for one fixed \(J\ge1\). At least one of the
finite-mark functions is nonconstant, since the three
finite values are distinct at every sufficiently large
\(t\) and have the same limit. Its fibres over a given
value have bounded cardinality by the degree bound.
Restricting \(t\) to large integers up to a sufficiently
small constant times \(X^{1/J}\) proves the lower bound
in (PSO1).

## Scope of the obstruction

The conclusion applies even to all parity class tests
on products of normalized summands whose weight sum
vanishes modulo two. Both complementary fourfold
products are exact squares and the joint defect
support is empty at every prime. Thus the two-class
subpower completion on an affine shape does not extend
to a projective shape alone. A pole class is an
additional arithmetic record. In this construction
the square class of \(f(z_3)\) is the class of \(f(t)\),
by multiplication of (PSO1b) over the three finite
marks. It is invisible to all the displayed relative
tests. In particular, at \(t=p>2\) prime it has odd
valuation at \(p\), so its support is not confined to
a fixed finite set.
[Projective shapes with a secant-corrected pole class](ABC_PROJECTIVE_SHAPE_SECANT_CLASS_COUNTS.md)

No upper bound of small-budget size is asserted for the original
three radical products of these rows. The obstruction concerns
completion from the shape and the indicated class data.

## An uncorrected class after fixing two original rows

**PSO2.** Fix a prime \(q\), \(s\ge4\), and weights
\(e_i\in\{1,\ldots,q-1\}\). Consider distinct positive
primitive rows \(A_i+B_i=C_i\le X^{B_0}\), with
\(a_i=A_i/C_i\) and \(b_i=B_i/C_i\). Prescribe a projective
shape given by distinct rational marks of height at most
\(X^{B_1}\), and fix two original rows in distinct slots.
There is a choice, determined by these fixed data, of one of
\[
                \prod_i b_i^{e_i},\qquad
                \prod_i a_i^{e_i},\qquad
                \prod_i(a_i b_i)^{e_i}.                       \tag{PSO2}
\]
If the outside defect radical of this chosen product at
order \(q\) is at most \(K\), then there are
\(O_\varepsilon(KX^\varepsilon)\) completions. The prescribed
prime pool may depend on the shape and the two fixed rows,
and must have prime product at most \(X^{B_2}\).
The constants are uniform at fixed \(q,s,B_0,B_1,B_2\).
The defect radical here uses primes outside the pool for
which the rational valuation of the product is nonzero
modulo \(q\). No secant class or row content is included.

Label the two fixed slots \(1,2\). Change the shape
coordinate to put their marks at \(0,\infty\), respectively.
The other marks \(r_i\) are finite, nonzero and distinct,
and their heights remain polynomial in \(X\).
For an actual tuple the map from the shape to \(b_i\),
followed by \((b-b_1)/(b-b_2)\), fixes \(0\) and \(\infty\).
It is multiplication by \(t\ne0\), so
\[
        b_i=\frac{b_1-b_2t r_i}{1-t r_i},
        \qquad a_i=\frac{a_1-a_2t r_i}{1-t r_i}
                        \qquad(i\ge3).
\]
The parameter \(t\) has polynomial rational height by
any of these nonanchor rows. Write \(S=\{1/r_i:i\ge3\}\),
with weight \(e_i\) at \(1/r_i\), and put \(u=b_1/b_2\),
\(v=a_1/a_2\). The first product in (PSO2), as a rational
function of \(t\), has zeros at \(uS\) and poles at \(S\)
before cancellations. The second instead has zeros at
\(vS\) and the same poles. The weight at a zero or pole
is the weight of its original element of \(S\), with
negative sign for poles. The anchors change only the
constant. The products have no net order at infinity.

For a positive dilation \(w\ne1\), join \(z\) to \(wz\)
whenever both belong to \(S\). This is a finite collection
of chains because every element of \(S\) is nonzero.
At the first element of each chain there is an uncancelled
pole, and at \(w\) times the last element there is an
uncancelled zero. These endpoints remain distinct across
chains, and both are branches of the \(q\)-cover because
their weights are nonzero modulo \(q\). At an interior
point the residue modulo \(q\) is the difference of two
successive weights. Thus there are at least two branches
for each product. There are exactly two branches only
when \(S\) is one consecutive chain with a constant
weight along it.

If the first or second product has at least three branches
for odd \(q\), or at least four for \(q=2\), choose that
product. The corresponding split cover has positive
genus by Riemann–Hurwitz. A parity cover has an even number
of branches, as all nonzero residues are one and their
sum is even, so the only remaining case for every \(q\)
has exactly two branches in each product. After swapping
the anchors if needed, \(u<1<v\). In the remaining
case \(S\) is a single chain for both \(u\) and \(v\), so
the ratio of the largest and smallest absolute values in
\(S\) is both \(u^{-(|S|-1)}\) and \(v^{|S|-1}\).
Since \(|S|\ge2\), this forces \(v=u^{-1}\). Write
\(S=\{z,zu,\ldots,zu^m\}\) with \(m=|S|-1\ge1\).
The two branch sets are then
\[
                         \{z,zu^{m+1}\},\qquad
                         \{z/u,zu^m\}.
\]
They are disjoint. Hence the third product in (PSO2)
has four branches and gives a positive genus cover.
This proves that a suitable choice always exists and
is determined before the remaining rows are counted.

Enumerate the outside support of the chosen defect and
its valuation residues modulo \(q\), together with the
residues on the fixed pool. This gives \(O(KX^\varepsilon)\)
exact classes, with an arbitrarily small inner allowance.
Indeed, all occurring support primes divide a
polynomial-height product of original row letters, and
the finite exponent choices on the support and pool
cost a subpower factor at fixed \(q\). The products
are positive on the original rows, so no sign choice
is needed. At each prescribed class, the chosen
rational function of \(t\) lies in that class exactly
when \(t\) lifts to a rational point on the corresponding
split \(q\)-cover, with a prescribed nonzero coefficient.
Remove \(q\)-th powers of linear factors to get exponents
between one and \(q-1\) at the actual branches. Their
sum is divisible by \(q\) since infinity is unramified.
The coefficients and point heights in this model
are polynomial in \(X\).

For genus at least two the split cyclic descent and
uniform Mordell–Lang argument of (IC2) give
\(O(X^\varepsilon)\) choices of \(t\). The genus-one
case \(q=2\) has four rational branches and a monic
split cubic model after sending a branch to infinity.
Its rational \(2\)-torsion gives the same bound.
The other genus-one case has \(q=3\) and three
rational branches. The three nonzero residues agree
modulo \(3\), since their sum is zero. Invert the
cubic class if necessary, move the branches to
\(0,1,\infty\), and use \(y^3=d z(z-1)\).
With \(h=4/d\), the point \((hy,h(2z-1))\) lies on
\(Y^2=X^3+h^2\), and \((0,h)\) is a rational
point of order three. Its bounded-height count is
again subpower. These changes of coordinate preserve
polynomial height, and can omit only a bounded number
of parameters. Finally each \(t\) recovers all the
reduced fractions \(b_i\) and the original primitive
rows. This proves (PSO2).
[Split covers and their genus-one models (IC2)](ABC_SPLIT_INDEX_POWER_CORRESPONDENCES.md#a-joint-count-from-the-branch-set)
· [Dujella, Theorem 1.1](https://arxiv.org/html/2312.03655v1#S1)

At fixed \(L>0\), put \(n=X^L\). A family of \(T\ge n^d\)
original rows has at most
\(O_\varepsilon(T^2 n^{\alpha+\beta}X^\varepsilon)\)
distinct ordered \(s\)-tuples with a projective-shape
dictionary of size \(n^\alpha\) and chosen defect
radical at most \(n^\beta\). If at least a fraction
\(n^{-\tau}\) of its tuples qualifies, then
\[
                       (s-2)d-\tau\le\alpha+\beta
\]
after taking the inner allowance arbitrarily small.
In particular, two actual anchors avoid the corrected
pole-class price when the uncorrected shape and chosen
class are jointly affordable. The shape and class tests
are invariant under separate common scaling of each
row and apply to its least-repaired primitive quotient.
[Secant-corrected single-class comparison (PM8)](ABC_PROJECTIVE_SHAPE_SECANT_CLASS_COUNTS.md#two-classes-with-small-original-overlap)
