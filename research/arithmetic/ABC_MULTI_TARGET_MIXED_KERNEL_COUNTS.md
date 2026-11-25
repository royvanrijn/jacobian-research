# Mixed squarefree kernels of multiple actual targets

A fixed vector of complete cross ratios makes the targets in a
common-neighbor star into linear forms in two primitive integers.
A low squarefree kernel of the product of four distinct forms
then puts the star on an elliptic curve with rational
\(2\)-torsion. Six or more forms of even total degree also
give a higher-genus count. With four targets, the count
can use their kernels alone.

## Fixed cross ratios and the common factors

Fix an integer \(m\ge2\), signs \(\sigma,\tau\in\{-1,1\}\)
and an ordered family of positive primitive triples
\[
                 T_i=\sigma U_i+\tau V_i,\qquad
                 (U_i,V_i)=1,\qquad U_i,V_i,T_i\le X
                   \quad(1\le i\le m).
\]
The signs are common to the family. Prescribe the positive
reduced rational numbers
\[
 R_i=\frac{U_iV_1}{U_1V_i}=\frac{r_i}{s_i},\qquad R_1=1 .
\]
Write \(x=\gcd_i U_i\), \(y=\gcd_i V_i\) and
\(U_i=d_i x,\ V_i=e_i y\). Define
\[
 P_{\boldsymbol R}
  =\left(\prod_i r_i s_i\right)
     \prod_{\substack{i<j\\R_i\ne R_j}} |r_i s_j-r_j s_i|.
\]
Only its prime support will be used. If there is a candidate,
\(r_i,s_i\le X^2\), so \(P_{\boldsymbol R}\le X^{O_m(1)}\)
uniformly. We have \((x,y)=1\).

For a prime \(p\nmid P_{\boldsymbol R}\),
\(v_p(U_i)-v_p(V_i)=v_p(U_1)-v_p(V_1)\) for every \(i\).
Each pair of nonnegative valuations has minimum zero by
primitivity, so the difference determines that pair.
Consequently \(p\) divides none of the \(d_i,e_i\).
Fixing all of these coefficients costs \(X^{o_X(1)}\),
uniformly for the given ratios, by the fixed-support
count at polynomial height.
[Fixed-support count](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound) After they are fixed the actual
targets are
\[
                     L_i(x,y)=\sigma d_i x+\tau e_i y .
\]
The forms \(L_i,L_j\) are proportional exactly when \(R_i=R_j\).
Each is distinct from both coordinate forms \(x,y\).
Furthermore the determinant of the coefficient vectors of
\(L_i,L_j\), when nonzero, has prime support in
\(P_{\boldsymbol R}\). Indeed \(d_i e_j/(d_j e_i)=R_i/R_j\),
and outside the first product in \(P_{\boldsymbol R}\) all
\(d_i,e_i\) are units. The same support statement holds for the
determinants involving the coordinate forms. In particular the
values of any two nonproportional forms in this collection
are coprime outside \(P_{\boldsymbol R}\), since a
prime with nonzero coefficient determinant cannot divide both
values at a primitive \((x,y)\).

## The mixed-kernel count

Select \(2g+2\) pairwise nonproportional forms from
\(x,y,L_1,\ldots,L_m\), for some \(g\ge1\), and put
\(F(x,y)\) equal to their product. For four forms one may use
two distinct target directions together with \(x,y\), or
four distinct target directions alone. Let \(D=\operatorname{sf}_2(F(x,y))\)
be the positive squarefree kernel of the product of the actual
selected values. Let \(S_0\) be any fixed prime set with
\(\prod_{p\in S_0}p\le X^C\), for fixed \(C\).
With all the preceding data fixed except for the
coefficients \(d_i,e_i\) and the actual triples, impose the cap
\[
              \prod_{\substack{p\mid D\\
                       p\notin S_0,\ p\nmid P_{\boldsymbol R}}}p
                        \le H_0 .
\]
For \(H_0\ge1\), the ordered tuple count is
\[
                  O_{m,C,\varepsilon}(H_0 X^\varepsilon).
 \tag{MTK1}
\]
The selection of the forms can depend on the tuple, at
a loss bounded in terms of \(m\). The original radical caps on
all actual letters can be kept. Outside \(S_0\) and the ratio
support \(P_{\boldsymbol R}\), the cap is exactly the product
of the odd-valuation prime contributions of the selected
forms, without overlaps. Here the coordinate forms, if used,
are the actual common factors \(x,y\), with radicals bounded
by the respective original letter caps. The target forms
are the complete original targets.

To prove (MTK1), fix the coefficients as above. The portion
of \(D\) outside the specified support has at most \(H_0\)
integer values. Its part on the specified support has at most
\(2^{\omega((\prod_{p\in S_0}p)P_{\boldsymbol R})}=X^{o_X(1)}\)
values. Fix \(D\). Since \(F(x,y)=D w^2\), the tuple gives
a rational point on the double cover
\[
          D z^2=F(t,1),\qquad t=x/y,\qquad z=w/y^{g+1}.
\]
The branch points are the \(2g+2\) distinct projective zeroes of
\(F\), all rational. This includes infinity exactly when the
dehomogenized degree is \(2g+1\). The exponent \(g+1\)
in the denominator of \(z\) uses the even homogeneous degree.

### Four forms and elliptic points of bounded height

Write \(F=L_0L_1L_2L_3\), using temporary indices for the
selected forms. Choose a coordinate form \(M\) independent
of \(L_0\). Write \(L_j=\alpha_j M+\beta_j L_0\) for
\(1\le j\le3\). All \(\alpha_j\) are nonzero because the
forms are pairwise nonproportional. The rational coefficients
have polynomial height by Cramer's rule. Since every
selected form has a positive value at the actual pair,
\(L_0(x,y)\ne0\). Moving its zero to infinity, put
\[
             u=\frac{M(x,y)}{L_0(x,y)},\qquad
             v=\frac{w}{L_0(x,y)^2}.
\]
Then \(D v^2=\prod_{j=1}^3(\alpha_j u+\beta_j)\).
After division by \(D\), write the right side as
\(a u^3+b u^2+c u+d_0\), with \(a\ne0\). The coordinates
\(\xi=a u,\ \eta=a v\) satisfy
\[
                 E:\quad \eta^2=\xi^3+b\xi^2+ac\,\xi+a^2d_0.
\]
The distinct rational roots of the cubic give rational
\(2\)-torsion on \(E\). The defining coefficients of \(E\)
and the actual coordinates \(\xi,\eta\) all have polynomial
height in \(X\), uniformly. Indeed the selected values
are at most \(X\), so \(D\le X^4\), \(w\le X^2\), and
the transformations use only a bounded number of products,
sums and quotients of polynomial-height numbers.

Take \(B_1=X^{O(1)}\) sufficiently large to bound the
exponential height of \(\xi\) and exceed
\(e^{\max(h(E),e)}\), with \(h(E)\) the logarithmic
height of the displayed coefficients. Over the fixed
field \(\mathbb Q\), Dujella's bounded-height estimate
for elliptic curves with rational \(2\)-torsion bounds
these points by \(B_1^{O(1)/\log\log B_1}=X^{o_X(1)}\)
uniformly in the fixed records. The coordinate \(u\)
determines the primitive pair \(x,y\).
[Dujella, Theorem 1.1 and Corollary 1.3](https://arxiv.org/html/2312.03655v1)

### Six or more forms and a height-independent point count

Suppose \(g\ge2\). Let \(J\) be the Jacobian of the smooth
proper curve.
The branch points give full rational \(J[2]\).
For example move one branch point to infinity and form the
\(2g+1\) differences from the finite branches to that point.
Each is killed by two. In a relation among those differences
the defining function of the principal divisor is a deck
eigenfunction, since its divisor is invariant. A rational
function of the line parameter has even orders at the
branches, whereas the double-cover coordinate times such a
function has odd orders at all finite branches. The only
nonzero relation among subsets is consequently the sum of
all \(2g+1\) differences, which is the divisor of the
double-cover coordinate. This gives \(2g\) independent
classes and hence all of \(J[2]\).
[Poonen–Schaefer, Proposition 6.2](https://math.mit.edu/~poonen/papers/descent.pdf)

Outside the prime set \(S=\{p:p\mid 2D P_{\boldsymbol R}\}\),
the projective zeroes remain distinct and both \(2\) and
\(D\) are units. The degree-\(2g+2\) double cover glues its
two affine line charts by dividing the double-cover
coordinate by the \((g+1)\)-st power of the line parameter.
The resulting finite cover of the projective line is smooth
over \(\mathbb Z_p\): at a branch the derivative in the line
parameter is a unit, and away from branches the derivative
in the double-cover coordinate is a unit. Thus \(J\) has
good reduction at \(p\notin S\). Its \(2\)-Kummer
classes are unramified there, because a rational point
extends to the abelian scheme over \(\mathbb Z_p\) and its
inverse image under multiplication by \(2\) is finite étale.

The global Kummer map embeds \(J(\mathbb Q)/2J(\mathbb Q)\)
in \(H^1(\mathbb Q,J[2])\). With a rational basis of \(J[2]\)
the \(2g\) components are rational squareclasses. In each
component, unramifiedness outside \(S\) permits only
\(\pm\prod_{p\in S}p^{b_p}\), \(b_p\in\{0,1\}\).
Consequently
\[
                  \operatorname{rank}J(\mathbb Q)
                            \le 2g(1+|S|).
\]
Both \(D\le X^{2g+2}\) and \(P_{\boldsymbol R}\le X^{O_m(1)}\)
have polynomial height, uniformly. Therefore
\(|S|=O_m(\log X/\log\log X)\) for large \(X\).
Embed the curve in its Jacobian using a rational
branch point. Its degree with respect to the canonical
polarization is bounded in terms of \(g\), and its genus
is at least two, so it contains no positive-dimensional
abelian translate. Uniform Mordell–Lang bounds the number
of rational points by \(O_g(1)^{1+\operatorname{rank}J(\mathbb Q)}
=X^{o_X(1)}\), uniformly in the coefficients.
Since \(t\) recovers the coprime positive \(x,y\), each point
determines the tuple in the fixed coefficient record.
Summing the subpower completion bound over the fixed records
proves (MTK1) in both cases.
[Gao–Ge–Kühne, Theorem 1.1′](https://arxiv.org/html/2105.15085v4)

## A common-neighbor bound

Consider a simple outside-projection graph with one actual
primitive triple per edge. Assume that the complete triple
determines its outside labels. For a fixed ordered tuple of
\(m\) distinct rows, suppose the common columns have at most
\(K\) different complete cross-ratio vectors. Assume that
each vector has the required number of distinct target
directions. In any common column for which the selected
forms satisfy the mixed-kernel cap \(H_0\), the associated
tuple recovers the column. Thus at most
\(K H_0 X^{\varepsilon}\) common columns satisfy the cap.

If these properties hold for every ordered \(m\)-tuple
of rows, with all common columns satisfying the cap,
write \(V_U,V_V\) for the two part sizes and \(E\)
for the number of edges. Counting ordered \(m\)-stars
over the columns gives
\[
          E\ll_m V_V+
                   V_U V_V^{1-1/m}(K H_0 X^\varepsilon)^{1/m}.
 \tag{MTK2}
\]
For \(E>4m V_V\), columns of degree at least \(2m\)
carry at least \(E/2\) edges. Their ordered star count
is \(\gg_m E^m/V_V^{m-1}\), by the power mean and the
falling-factorial bound at that degree. Counting by
ordered row tuples instead bounds the same number by
\(V_U^m K H_0 X^\varepsilon\).

In the outside-projection setup of (TR2) and (OA1), a
fixed cheap full target modulus gives at most \(n^{o(1)}\)
cross-ratio vectors for each fixed row tuple, at
\(n=X^L\), in the stated iterated endpoint limits.
The single-anchor thinning of (OA1) can also fix those
ratios outright for the retained rows.
Distinct rows give distinct ratios within the column,
since equal complete odds of primitive triples imply equal
complete letters. In either setup, (MTK1) allows
\(n^{o(1)}\) columns for each individual mixed-kernel
record. It makes no bound on the size of the record
dictionary without the actual kernel cap.
[Complete cross ratios and one-anchor separation](ABC_DIFFUSE_FACTOR_GRAPH.md#a-single-anchor-in-one-projection)

## A saving near the minimum nonatomic mass with only a target parity cap

Let a putative balanced endpoint box have normalized count
exponent \(d=1/3\). Write \(c_U,c_V,c_T\) for its nonatomic
masses and \(\Gamma=c_U+c_V+c_T\). Suppose
\(\Gamma\le1/2+\gamma\), where \(0\le\gamma<1/6\).
A zero nonatomic mass at the balanced endpoint would force
both other masses to be \(1/3\). Thus all three are positive,
and the pairwise diffusion inequalities give
\[
          1/6-\gamma\le c_i\le1/6+\gamma\qquad(i=U,V,T).
\]
In particular \(c_T>0\). Assume one actual target letter
in this box satisfies \(\operatorname{sf}_2(T)\le n^f\)
up to a vanishing normalized tolerance. The argument
also allows the same cap only outside a fixed
polynomial-height prime pool. There is no parity-kernel
condition on either summand. The target designation may
be any of the three letters after fixing the signs.

Use the \(U,V\) projection and freeze a positive-weight
full target block at arbitrarily small normalized root
cost. The pair-uniqueness thinning and cross-ratio bounds
from (TR2) and (OA1) apply with arbitrarily small
normalized losses. Use four rows and choose their four
actual targets as the four forms, so their product's
outside squarefree kernel costs at most \(n^{4f+o(1)}\).
Then (MTK2), with the errors and freezing cost tending to
zero in the endpoint construction, implies
\[
           d\le\max\!\left(c_V,\ c_U+\frac34c_V+f\right).
 \tag{MTK3}
\]
Thus no such balanced endpoint can satisfy
\[
                 \frac74\gamma+f<\frac1{24}.
\]
Indeed \(c_U+\tfrac34c_V+f\le 7/24+\tfrac74\gamma+f\).
For example \(\gamma=f=1/100\) gives an upper bound
\(383/1200<1/3\) for the right side of (MTK3). This excludes
the corresponding near-minimum nonatomic sector with a parity cap on only the complete
target. It applies to original primitive valuations
before any common multiplier. Without the low mixed
defect the proof does not control general outside-graph
occupancy.
[Pairwise diffusion and nonatomic masses](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
