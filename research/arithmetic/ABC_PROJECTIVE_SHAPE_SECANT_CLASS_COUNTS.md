# Projective shapes and secant-corrected power classes

A projective shape fixes an ordered tuple of fractions up to a
fractional linear transformation. Joint power classes on the two
complementary fractions determine ratios between three binary-form
values. A weighted product of the actual secants determines the
remaining class. After this correction, three split covers recover
the three linear forms of the transformation separately. For three
rows the projective shape requires no record.

## The secant correction

Fix \(s\ge3\), a prime \(q\), and integers \(1\le e_i<q\)
with \(\sum_{i=1}^s e_i\equiv0\bmod q\). Take an ordered
tuple of distinct primitive positive triples
\(A_i+B_i=C_i\le X^{B_0}\) and write \(b_i=B_i/C_i\).
The \(b_i\) are distinct because every reduced fraction \(b_i\)
recovers its original row. Prescribe the projective shape by
distinct rational marks \(r_1,\dots,r_s\) of height at most
\(X^{B_1}\), with \(b_i=m(r_i)\) for a common fractional linear
map \(m\) over \(\mathbb Q\). We can use a polynomial-height
coordinate in which all \(r_i\) are finite. Indeed, among
\(s+1\) bounded rational pole candidates choose one away from
the marks and send it to infinity. Coordinate height bounds
below may be enlarged by fixed powers.

Choose integer edge weights \(f_{ij}\) for \(1\le i<j\le s\)
subject to
\[
        \sum_{j<i}f_{ji}+\sum_{j>i}f_{ij}\equiv e_i\bmod q,
        \qquad
        \sum_{i<j} f_{ij}\equiv0\bmod q .
                                                                \tag{PM1}
\]
For odd \(q\) take \(f_{1j}=e_j\) for \(j>3\), put all
other edges meeting those vertices to zero, and solve the
remaining vertex sums on the triangle \(1,2,3\). Its unsigned
vertex-edge incidence matrix has determinant of absolute value
two. Summing the vertex congruences gives the total edge
congruence. For \(q=2\) all \(e_i=1\) and \(s\) is even.
Use a perfect matching, then add the weights of any triangle
modulo two if the total edge weight is odd. This leaves
the vertex sums unchanged. Fix integer representatives of the
weights, of absolute value at most \(q\).

Define the signed secants and weighted products by
\[
  S_{ij}=B_j C_i-B_i C_j,\quad
  S_{\boldsymbol f}=\prod_{i<j} S_{ij}^{f_{ij}},\qquad
  T_U=\prod_{i=1}^s U_i^{e_i}\quad(U=A,B,C).
\]
The \(S_{ij}\) are nonzero. The three records in
\(\mathbb Q^\times/(\mathbb Q^\times)^q\) are
\[
                    [T_A/S_{\boldsymbol f}],\qquad
                    [T_B/S_{\boldsymbol f}],\qquad
                    [T_C/S_{\boldsymbol f}].                 \tag{PM2}
\]
They are invariant under any rational positive scaling of
each row separately. Scaling row \(i\) changes the valuation
of each representative by a multiple of \(q\), because the
exponent of its scale is the difference in the \(i\)-th vertex
congruence of (PM1). In particular these are records of the
primitive quotients of minimally repaired triples.

To see the correction directly, write \(m=L_B/L_C\) in
linear polynomials and put \(L_A=L_C-L_B\).
For \(U=A,B,C\) one has \(L_U(r_i)=\rho_i U_i\) with
a common nonzero rational \(\rho_i=L_C(r_i)/C_i\).
If \(L_B(z)=u z+v\) and \(L_C(z)=w z+x\), then
\[
 \rho_i\rho_j S_{ij}
    =L_B(r_j)L_C(r_i)-L_B(r_i)L_C(r_j)
    =(ux-vw)(r_j-r_i).
\]
The determinant \(ux-vw\) is nonzero. Set
\(\Theta=\prod_{i<j}(r_j-r_i)^{f_{ij}}\). Taking classes and
applying both congruences in (PM1) gives
\[
          \left[\prod_i L_U(r_i)^{e_i}\right]
                         =[\Theta\, T_U/S_{\boldsymbol f}]
                    \qquad(U=A,B,C).                         \tag{PM3}
\]
There is no extra charge for the determinant of the moving map.
Equivalently, two records are the joint quotient classes of
\(\prod_i(A_i/C_i)^{e_i}\), \(\prod_i(B_i/C_i)^{e_i}\),
and the third is the weighted class of
\(\prod_{i<j}(b_j-b_i)^{f_{ij}}\). These two descriptions
are related by an invertible change of basis modulo \(q\),
so they have the same union of defect supports.

## Counting from the corrected classes

**PM4.** Fix \(q,s,B_0,B_1,B_2\). For any fixed distinct
marks of the prescribed heights, fix a prime pool \(P\) of
product at most \(X^{B_2}\), possibly depending on the shape.
The defect support of a class consists of primes outside
\(P\) where its valuation is nonzero modulo \(q\).
If the product of all distinct primes in the union of the
three defect supports in (PM2) is at most \(K\), the number of
ordered tuples of that shape is
\[
                        O_{\varepsilon}(K X^\varepsilon).
                                                                \tag{PM4}
\]
The constant depends only on the fixed parameters and
\(\varepsilon>0\). For \(q=2\) the weights force \(s\ge4\) even.
For odd \(q\) the assertion includes \(s=3\).
The exact classes version gives \(O(X^\varepsilon)\) tuples.

First enumerate the union support and fix the three classes.
All outside support primes divide an actual polynomial-height
product of row and secant magnitudes. Their valuation residues,
and those on the prescribed pool, cost subpower choices at
fixed \(q\). In the quadratic case record the three signs as
well. Representatives of the three classes have polynomial
height, by choosing reduced exponents modulo \(q\).
The coefficients of \(m\), up to a common scalar, have
polynomial rational height. The maps sending respectively
the first three \(r_i\) and the first three \(b_i\) to
\(0,1,\infty\) have polynomial-height coefficient ratios,
by their formulas in three point differences, and their
composition determines \(m\).

For one of the linear forms \(L=L_U\), the class of
\(\prod_i L(r_i)^{e_i}\) is well-defined on its projective
coefficient vector because the total degree is divisible
by \(q\). A nonzero constant form gives one such vector.
For a nonconstant form write \(L(z)=a(z-t)\). The specified
class means \(t\) lifts to
\[
                           y^q=d\prod_{i=1}^s(t-r_i)^{e_i}
\]
with fixed nonzero rational \(d\) of polynomial height.
The distinct \(r_i\) are all fully ramified, and infinity
is unramified. The genus is \((q-1)(s-2)/2\).
In genus at least two apply the split cyclic-cover argument
of (IC2). Its cyclotomic isogeny descent uses the rational
branch divisors over \(\mathbb Q(\zeta_q)\) and confines Kummer
classes to the primes over \(q\), the coefficient numerator
and denominator, and the branch denominators and differences.
These form a polynomial-height product, so the Jacobian rank
is \(O(1+\omega(S))\) for such a product \(S\), at fixed
\(s,q\). Uniform Mordell–Lang gives a subpower count.
In genus one, use a branch change over \(\mathbb Q\).
For \(q=2,s=4\) this gives a split cubic with rational
\(2\)-torsion. For \(q=3,s=3\) all three weights agree, and
inverting the class if necessary gives \(y^3=d' z(z-1)\).
With \(a'=4/d'\), the point \((a'y,a'(2z-1))\) lies on
\(Y^2=X^3+(a')^2\) with the rational flex \((0,a')\) of order
three. Both models and the transformed point abscissas have
polynomial height, and Dujella's theorem gives a subpower
count in each. Thus each of the three projective forms has
\(O(X^\varepsilon)\) possibilities with a suitably smaller
allowance at each step.
[Split covers and the descent argument (IC2)](ABC_SPLIT_INDEX_POWER_CORRESPONDENCES.md#a-joint-count-from-the-branch-set)
· [Gao–Ge–Kühne, Theorem 1.1′](https://arxiv.org/html/2105.15085v4)
· [Dujella, Theorem 1.1](https://arxiv.org/html/2312.03655v1#S1)

Finally \(L_A,L_B,L_C\) are pairwise nonproportional because
\(m\) is nonsingular. After choosing a representative for
each of the three projective forms, the relation
\(L_A+L_B=L_C\) fixes the relative scales. For example expand
the chosen \(C\)-form in the basis of the chosen \(A,B\)-forms.
Both expansion coefficients are nonzero. They are the
required two scales relative to the scale of the \(C\)-form.
Thus the three forms determine one \(m\), and each of its
values \(b_i\) determines the original primitive row. This
proves (PM4).

## The cubic case with no shape record

**PM5.** For \(s=q=3\) use \(r=(0,1,2)\) and \(e_i=1\).
Every ordered triple of distinct \(b_i\) has this
projective shape. Taking \(f_{ij}=-1\) on the triangle,
the three records are simply the cube classes of
\[
       \left(\prod_i A_i\right)  S_{12}S_{13}S_{23},\qquad
       \left(\prod_i B_i\right)  S_{12}S_{13}S_{23},\qquad
       \left(\prod_i C_i\right)  S_{12}S_{13}S_{23}.
                                                                \tag{PM5}
\]
Their union defect cap counts tuples by (PM4) without any
prescribed shape or shift. The estimate has a second proof
using the uniform primitive diagonal cubic count. For
each projective form use primitive integral coefficients,
so \(L(0),L(1),L(2)\) are integers. At a prime outside
the shape determinant pool, here just \(\{2\}\), at most
one of these three values is divisible by the prime. The
specified class of their product therefore contains every
prime needed in their three cube-free parts outside the pool.
Fix these parts at subpower cost on the specified support
and the pool, and write the values as three signed cubes
times those fixed coefficients. Their relation is
\(L(0)+L(2)=2L(1)\), and the three roots are primitive since
\(\gcd(L(0),L(1))=1\). Heath-Brown's theorem counts the
roots at subpower cost, uniformly in the fixed coefficients.
They recover each projective form, and then the tuple.
[Heath-Brown, Theorem 1.3](https://arxiv.org/html/2601.07817v1#S1)

For parity quadruples take \(e_i=1\) and any of the three
perfect matchings for the edge weights. Formula (PM2) then
divides each of the original fourfold letter products by the
product of the two chosen secants. Only one projective shape
parameter is needed, for example the ordered cross ratio
of the four fractions. A moving map with just the two
uncorrected quotient classes can have a power-sized family
at one fixed projective shape.
[Explicit two-class obstruction](ABC_PROJECTIVE_SHAPE_TWO_CLASS_OBSTRUCTION.md)
· [Two-class count on affine shapes (AS1)](ABC_AFFINE_SHAPE_JOINT_POWER_CLASSES.md#two-power-class-records-determine-an-affine-shape)

## Scope at the small-budget endpoint

At fixed \(L>0\), put \(n=X^L\). If a family has
\(T\ge n^d\) primitive rows and a fraction at least \(n^{-\tau}\)
of their distinct ordered \(s\)-tuples has projective shape
dictionary of size at most \(n^\alpha\) and corrected union
defect cap \(n^\beta\), then necessarily
\[
                           sd-\tau\le\alpha+\beta.            \tag{PM6}
\]
Sum (PM4) with its inner allowance smaller than any fixed
gap times \(L\). In the cubic three-row case \(\alpha=0\).
In particular any fixed saving below the full tuple cost in
these corrected prime records would rule out the corresponding
full-exponent cohort, before taking the small-budget limit.
This test can count distinct least-repaired targets by their
primitive quotients.

The correction still prices fresh primes. At a prime dividing
one active secant and no other secant or original letter,
all three records have valuation residue \(-f_{ij}v_p(S_{ij})\).
If this is nonzero, the prime contributes once to the union.
Thus private secant mass remains an obstacle even when the
original quotient classes are cheap. At fixed shape with
\(s\ge4\), an exclusive prime must divide one of the fixed
shape determinants. Putting them in the pool removes that
charge, while the shape dictionary itself has a power cost
on typical affordable-fibre tuples. Outside the shape pool,
the remaining secant correction is the weighted class of
the primitive row contents. A saving using (PM6) needs
control of both the shape and the resulting corrected support.
[Private primes and projective shape growth (PG1–PG3)](ABC_AFFORDABLE_FIBRE_PROJECTIVE_SHAPE_GROWTH.md)
[Private prime contributions of actual secants (CR1)](ABC_CROSS_RATIO_FRESH_PRIME_OBSTRUCTION.md#exclusive-prime-packets-on-bounded-tuples)

## Two classes with small original overlap

**PM7.** In (PM4), fix two distinct roles \(U,V\), and
impose only the union defect cap \(K\) on their two corrected
classes. Impose also
\[
 \operatorname{rad}_{\overline P}\!\left(
      \gcd_{1\le i\le s} U_i\ \gcd_{1\le i\le s} V_i\right)
                                     \le G .
\]
Then the count is \(O_\varepsilon(KG X^\varepsilon)\).
In a full-exponent balanced endpoint cohort, the common-gcd
charge can be made arbitrarily small in the normalized
exponent.

After the two corrected classes are fixed, (PM3) and the
split-cover argument give subpower choices for the two
projective forms. Choose a polynomial-height representative
for each and denote the evaluated ratio of these
representatives at \(r_i\) by \(v_i\ne0\). There is a
common scalar \(k\in\mathbb Q^\times\) with
\[
                              U_i/V_i=k v_i \qquad(1\le i\le s).
\]
It has polynomial height by any actual row. Outside the
prime pool of all \(v_i\), every prime in the numerator of
\(k\) occurs in every \(U_i\), since the ratio is reduced,
and every prime of its denominator occurs in every \(V_i\).
Those two support sets are disjoint. The outside support
cap and supported-power enumeration therefore fix \(k\)
at cost \(O(G X^\varepsilon)\). The resulting ratios of
two distinct original letters determine their primitive
rows by the additive relation. In particular no third
corrected class is priced.

For the endpoint application, fix one of the rows as a
reference and use its original letter supports in the
prescribed-pool overlap estimate (DRR4). A second row with
a fixed positive normalized mass on that pool loses a fixed
power, at each sufficiently small fixed budget. On a
full-exponent cohort of \(s\)-tuples discard those pairs
and freely count the other \(s-2\) rows. This lets \(G\)
have arbitrarily small positive normalized exponent up to
an exceptional set with a fixed count loss. Thus the
comparison (PM6) at full exponent can use just the two
corrected classes and the shape. In the cubic case there
is no shape cost. This application still requires a bound
on the two corrected defects.
[Prescribed-pool overlap (DRR4)](ABC_DIVISOR_RECORD_RECIPROCITY.md#anisotropic-intersections-for-an-arbitrary-record)
· [Uniform supported-power enumeration (EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

**PM8.** With any two original rows in distinct slots fixed,
one corrected class with defect cap \(K\) suffices for
\(O_\varepsilon(K X^\varepsilon)\) completions at fixed
projective shape. Indeed, a fixed corrected class recovers
one projective linear form at subpower cost. If this is
\(L_C\), the two fixed ratios \(b_i\) specify \(L_B(r_i)\)
at the two anchor marks and so determine \(L_B\). If
the recovered form is \(L_A\) or \(L_B\), first fix its
scale and recover \(L_C\) by prescribing its values at the
two anchors using \(1-b_i\) or \(b_i\). Interpolation at the
distinct marks determines the other form. The same bound
then follows with the union record replaced by the single
defect record.

For example, with the shape-free cubic three-row
construction, a family of \(T\) rows has at most
\(O_\varepsilon(T^2 K X^\varepsilon)\) triples with that
single corrected defect at most \(K\). This imposes an
individual class cost of at least the family exponent for
any full-exponent cohort, whereas (PM7) imposes the higher
two-class joint cost. A private secant prime is still
detected by either corrected class under the condition
following (PM6).
