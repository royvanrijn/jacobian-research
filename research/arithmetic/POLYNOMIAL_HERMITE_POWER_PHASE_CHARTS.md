# Rational charts for a positive-dimensional power-phase family

Let
\[
 a\ge2,\qquad n>2a,\qquad s=n-a-1,\qquad N=2n.
                                                               \tag{1}
\]
For full braids and fixed product \(P_{2n}\), consider the
polynomial passport
\[
                  (2^n),\quad(2^a1^{2n-2a}),\quad T^s.
                                                               \tag{2}
\]
There are \(\lfloor a/2\rfloor+1\) components, indexed by
the two phase-bin counts \(\{j,a-j\}\).
Every component is a rational variety over \(\mathbb Q\),
of complex dimension \(s\), after source and target affine
quotients. Explicit nonempty coefficient charts and exact
exclusion conditions are given below. All have the same
embedded monodromy \(S_n\wr C_2\), block-count divisors
\(\{1,2,2n\}\), and trivial deck group.

The component count and permutation data already specialize
the [involution-family classification](POLYNOMIAL_INVOLUTION_COMPONENT_CLASSIFICATION.md).
The contribution here is a direct coefficient-space proof,
rational charts, and exact algebraic realizations of the
positive-dimensional diagonals in
[G1--G3](POLYNOMIAL_PERIOD_DIAGONALS_AND_POWER_PHASES.md).
It does not implement the general critical-graph construction.

## Existing primary result and the new proof route

Bogatyrëv--Gendron,
[The space of solvable Pell--Abel equations](sources/bogatyrev-gendron-pell-abel-2025.pdf),
published Theorem 1.2 and Corollary 2, pp1484--1486, already
contain the component classification of (2). Their parameters
are \(A=n,B=a,g=s\), original degree \(2n\). Their degree
partitions here are
\[
                       \{n-2j,n-2(a-j)\}.                   \tag{3}
\]
The actual primary statement was checked. We do not assert
that the count or its prior classification is new.

The arguments below prove connectedness independently of
that component theorem. They use polynomial Hermite
interpolation, a nonempty algebraic open locus, and the
already proved transposition-normal-closure lemma L6.
Exact sizes use the elementary class-count formula B4
and the two-fold power interpretation, not a new count
derived from finite experiments.

## Theorem J1. Exact Hermite coefficient charts

Choose \(0\le j\le\lfloor a/2\rfloor\).
Let \(P_+,P_-\) be monic, squarefree, coprime polynomials
of degrees \(j,a-j\); put \(P=P_+P_-\). Let \(v\ne0\).
There is a unique polynomial \(Q\) of degree less than
\(2a\) such that
\[
                Q\equiv1\pmod{P_+^2},\qquad
                Q\equiv-1\pmod{P_-^2}.                      \tag{4}
\]
An empty bin uses \(P_+=1\), so \(Q=-1\).
Choose a monic polynomial \(R\) of degree \(n-2a\), and set
\[
                       H=vQ+P^2R,\qquad f=H^2.              \tag{5}
\]
All coefficients are rational functions of these parameters.
The open conditions below are necessary and sufficient for
(5) to have exactly passport (2).

Put
\[
          A(z)=\frac{H'(z)}{nP(z)},\qquad
          C(Y)=\operatorname{Norm}_{\mathbb C[z]/(A)}
                                      \bigl(Y-H(z)^2\bigr).
                                                               \tag{6}
\]
The polynomial \(A\) is monic of degree \(s\).
The exact open conditions are
\[
 \begin{gathered}
 \gcd(P,P')=\gcd(A,A')=\gcd(P,A)=1,\qquad v\ne0,\\
             C(0)\ne0,\qquad C(v^2)\ne0,\qquad
                              \operatorname{disc}(C)\ne0.
 \end{gathered}                                                \tag{7}
\]
Every polynomial in (2), up to source affine change and
target affine change, arises in this way. For
\(j<a/2\) the bin with fewer critical points fixes the
choice of sign for \(v\); for \(2j=a\), also quotient by
\[
                    (P_+,P_-,v,R)\mapsto(P_-,P_+,-v,R).
                                                               \tag{8}
\]

**Proof.** The Chinese remainder theorem gives (4),
with denominators powers of the resultant of \(P_+,P_-\).
The derivative of \(Q\) vanishes at every root of \(P\),
so \(P\mid H'\); monicity gives (6).
At roots of \(P_+\) the critical value of \(H\) is \(v\);
at roots of \(P_-\) it is \(-v\).

The first line of (7) says that these are simple critical
points and that the other \(s\) critical points are
simple and distinct from them. Their \(f\)-critical
values are precisely the roots of \(C\).
The second line says that none is zero, none equals
the bundled value \(v^2\), and they are mutually distinct.
Since \(v\ne0\), no critical value of \(H\) is zero.
Thus \(H\) has \(n\) simple roots, supplying exactly the
fiber \(2^n\) of \(f\) over zero. The identity
\(f'=2HH'\) accounts for all other ramification.

Conversely the distinguished fiber forces \(f=H^2\),
with \(H\) monic after scaling. All critical points
of \(H\) are simple, and the \(a\) points over the
bundled value of \(f\) have \(H\)-values \(v\) or \(-v\).
Their root polynomials are \(P_\pm\). Subtracting
\(vQ\) from \(H\) gives a polynomial divisible by
both squares, hence by \(P^2\). Its quotient is
the monic \(R\) in (5). This proves exhaustiveness
and the stated sign quotient. \(\square\)

## Theorem J2. Nonemptiness, components, and rationality

For every (1) and every indicated \(j\), the locus (7)
is nonempty and irreducible. Its image is one full
Hurwitz component. Different \(j\) give different
components. Each affine-quotient component has
dimension \(s\) and is rational over \(\mathbb Q\).

**Proof of nonemptiness.** Begin with
\(H_0=P^2R\), choosing all roots of \(P,R\) real and
pairwise distinct. The derivative is
\[
                  H_0'=P(2P'R+PR').
 \]
Rolle's theorem gives one further critical point in
each gap between consecutive distinct roots of \(H_0\).
There are \(s\) such gaps; counting the degree of the
derivative shows these critical points are simple,
and they have nonzero values.

Their values can be made pairwise distinct up to sign
by a small deformation preserving the double roots.
Here is a direct independence argument. If \(t_1,\ldots,t_s\)
are these critical points, tangent variations of the
monic polynomial preserving its \(a\) double zeros
are exactly polynomials \(P L\) of degree at most
\(n-1\). Indeed solve
\(L=2R\,\delta P+P\,\delta R\) by reducing modulo
\(P\), where \(R\) is invertible; the resulting
\(\delta P,\delta R\) have degrees below those of
\(P,R\), respectively. Evaluating at the \(t_i\) is surjective,
because \(P(t_i)\ne0\) and polynomial interpolation
assigns arbitrary values at \(s\) distinct points.
At a critical point its value varies by that same
evaluation. Therefore the finitely many equalities
of values or their negatives can all be avoided.

Now assign any desired signs to the roots of \(P\)
and use their Hermite interpolant \(Q\). For sufficiently
small nonzero \(v\), \(H_v=H_0+vQ\) has the prescribed
critical points with values \(\pm v\). The other
critical points and their nonzero, distinct-up-to-sign
values persist. All conditions (7) follow.
This is a uniform existence proof; finite witnesses
in the checker do not replace it.

**Proof of the component statement.** The coefficients
of \(P_\pm\), \(v\) and \(R\) form an affine rational
parameter space of dimension \(n-a+1\).
Conditions (7), squarefreeness and coprimality remove
proper algebraic subsets. The resulting nonempty open
set is irreducible and connected. It covers all
polynomials with that phase count by J1.
The unordered bin count is locally constant, by G3.
It follows that its images are exactly the components.
Source and target affine changes are connected
and do not create additional components.

**Proof of rationality.** Quotient translations by
centering the roots of \(P\); its coefficient of
\(z^{a-1}\) becomes zero. This uses a rational
translation by the root barycenter. The remaining
source scaling acts diagonally: the coefficient
of order \(i\) in each monic polynomial has weight
\(i\), and \(v\) has weight \(n\), up to reversing
all weights.

When \(2j=a\), the extra involution (8) commutes
with this action. Replace coefficients of \(P_\pm\)
by their sums and differences. It is then diagonal
with signs \(+1,-1\); \(v\) has sign \(-1\) and
the coefficients of \(R\) have sign \(+1\).
The centering constraint is invariant and linear.
Thus the function field is the invariant field
of a diagonal multiplicative group, and possibly
a diagonal sign involution, on independent rational
coordinates. On the dense coordinate torus the
invariant monomials form a free exponent lattice.
A basis gives algebraically independent Laurent
monomial coordinates for the invariant field.
It is purely transcendental over \(\mathbb Q\).
Subtracting the two source dimensions gives
\(n-a-1=s\), proving the claim. \(\square\)

## Theorem J3. Actual groups, every block, and exact sizes

Every component above has actual embedded group
\[
                            S_n\wr C_2
 \]
on the standard two residue blocks, group order
\(2(n!)^2\), block-count divisors \(\{1,2,2n\}\),
and deck order one.

Let
\[
 \begin{split}
 L&=(s+2)(s+1),\\
 C_j&=
 n^{s+1}
 \frac{(n-j-1)!(n-a+j-1)!}
 {j!(a-j)!(n-2j)!(n-2a+2j)!}.
 \end{split}
\]
The literal full orbit size and absolute full orbit
size of the component \(j\) are
\[
       M_j=L\,2^s(2-\mathbf1_{2j=a})C_j,\qquad
                         \frac{M_j}{2n}.                    \tag{9}
\]
Their total is the independent full class-character
count. For \(a=2\), the two sizes simplify to
\[
 \begin{split}
 M_0&=(n-1)(n-2)2^{n-3}(n-3)n^{n-3},\\
 M_1&=(n-1)(n-2)2^{n-3}n^{n-2}.
 \end{split}                                                  \tag{10}
\]

**Proof of the group statement.** L6's transposition
normal-closure blocks have some size \(b\).
They refine every proper block system. G3 supplies
the squaring-map blocks, so \(b\le n\).
The \(s\) simple factors each consume one return
index, giving \(b-1\ge s\), or \(b\ge n-a\).

Let \(m=2n/b\) be the number of canonical blocks.
The quotient is generated by the two bundled
involutions with full-cycle product. Hence it is
the standard dihedral group on its \(m\) sheets.
If \(m\ge4\) is even, \(b\le n/2<n-a\), impossible.
If \(m\ge3\) is odd, the fixed-point-free original
factor has a fixed quotient block, on which its
return is a fixed-point-free involution. That
contributes \(b/2\) return index in addition to
the \(s\) simple returns. Thus
\(b-1\ge s+b/2\), or \(b\ge2(n-a)>n\), again
impossible. Consequently \(m=2,b=n\).
L6 gives the full symmetric base and the stated
wreath group. Its block and centralizer assertions
give the lattice and deck count.

**Proof of the sizes.** The lower polynomial \(H\)
has the two prescribed involution values, with
\(j,a-j\) transpositions respectively, and \(s\)
other simple values. B4's elementary ordered
class formula gives exactly \(C_j\) literal
tuples, including an identity position when \(j=0\).
Every such lower tuple has group \(S_n\).
For a proper block size \(b\le n/2\), its simple
factors would require \(b-1\ge s\), contradicting
\(s+1=n-a>n/2\). Thus it is primitive and its
transposition normal closure is \(S_n\).
Dividing by \(n\) gives the absolute
degree for these fixed value roles.

For each of the \(s\) remaining critical values
of \(f\), choose either square root as the value
of \(H\), giving \(2^s\) choices. If the bin sizes
differ, the smaller bin fixes the overall sign.
If they agree, \(H\) and \(-H\) count the same
square and give a factor of two. Thus the absolute
degree for a fixed assignment of branch types is
\[
                     \frac{2^sC_j}{n(1+\mathbf1_{2j=a})}.
 \]
There are \(L\) type placements. Each absolute
class lifts to \(2n\) literal tuples because its
deck is trivial, and the full twist puts those
tuples in one orbit. This gives (9).
Substitution \(a=2\) gives (10). \(\square\)

## Proposition J4. Two explicit rational surfaces in degree ten

For \(n=5,a=2\), normalize the two bundled critical
points to \(0,1\). With \(P=z(z-1)\), the two
charts are
\[
 \begin{split}
 H_0&=-v+P^2(z+t),\\
 H_1&=v(1-6z^2+4z^3)+P^2(z+t),\qquad f_i=H_i^2.
 \end{split}                                                  \tag{11}
\]
Remove exactly the algebraic exclusions (7).
The two charts have parameter involutions
\[
       (t,v)\mapsto(-1-t,-v),\qquad
       (t,v)\mapsto(-1-t,v),                                 \tag{12}
\]
respectively. Rational quotient coordinates are
\[
       \left(\frac{2t+1}{v},v^2\right),\qquad
                         ((2t+1)^2,v).                       \tag{13}
\]
They give the two rational components of the
passport \((2^5),(2^2 1^6),T,T\), with literal
sizes \(2400,6000\), and absolute sizes \(240,600\).
For a fixed assignment of the distinguished
and bundled values and unordered remaining
two values, the degrees are \(20,50\).

**Proof.** The two-point Hermite interpolants
are \(-1\) and \(1-6z^2+4z^3\).
The monic residual factor is \(z+t\), proving
(11). Reversing \(0,1\) acts by
\(H(z)\mapsto-H(1-z)\), giving (12);
the invariant fields are exactly (13).
The critical quotient numerators in (6) are
\[
  5z^2+(4t-3)z-2t,\qquad
             5z^2+(4t-3)z+12v-2t.
 \]
All remaining open conditions are therefore
computed by a two-dimensional quotient algebra.
The sizes and degrees follow from J3.
\(\square\)

These components differ already in their intrinsic
\(S_5\wr C_2\) branch classes: the smaller bundled
factor has return types \((\mathrm{DT},1)\) or
\((T,T)\). This is a positive-dimensional phase
split, not a new genuine lifting split.
The equal-phase cubic pair in G4 is a separate
obstruction.

## Evidence and remaining scope

The [checker](verify_polynomial_phase_interpolation.py),
[receipt](polynomial_phase_interpolation_verification.json)
and [source ledger](polynomial_phase_interpolation_source_ledger.json)
retain an exact rational Hermite chart generator,
quotient-algebra critical-value polynomials and
all open-condition tests. Finite controls include
several values of \(a,n,j\), both degree-ten
surfaces and their involutions, complete \(8400\)-tuple
degree-ten coverage, exact size/character comparisons
and deliberately degenerate inputs.
The all-parameter nonemptiness and rationality claims
use J1--J3's written arguments.

The strict range \(n>2a\) is necessary for the
stated uniform group and open-chart conclusion.
The boundary range can have empty bins or
additional composition. General colored period
graphs and arbitrary disconnected passports
remain unclassified. Current audit and replay
are internal; independent replay, external review
and formal verification are not asserted.
