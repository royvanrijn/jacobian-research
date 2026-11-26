# Projective feature counts using the three original branch points

The primitive ratio \(a/c\) has three distinguished projective
values modulo a full divisor on one of the original letters.
A bounded-degree map of that ratio preserves these congruences,
including at poles. This gives a strict small-budget count when
the mapped ratio has a projective description in a bounded number
of integer features with small minors. The count keeps the original
three radical caps and allows arbitrary complete dictionary sizes.

## A projective feature count at all balanced profiles

**PF1.** For each fixed positive integer \(R\) there are
\(\theta_R,\delta_R,\kappa_R,L_R>0\) such that the following holds.
Fix \(0<L<L_R\), put \(n=X^L\), and count positive primitive
\(a+b=c\) in \(X/2<c\le X\) with
\[
                         \operatorname{rad}(W)\le n^{1/3+\kappa_R}
                    \quad(W=a,b,c).
\]
Take a prescribed nonconstant map
\(\Phi:\mathbb P^1\longrightarrow\mathbb P^1\) over \(\mathbb Q\),
represented by coprime integer binary forms of equal degree at
most \(D\) and coefficient height at most \(X^B\). The map may
vary with \(X\). The bounds \(D,B\) are fixed before the height
limit at each fixed \(L\) and may depend on \(L\).

Assign each row a vector \(z\in\mathbb Z^R\), and let \(r\) be
the affine dimension of all assigned vectors. Suppose the
nonzero anchored \(r\)-minors of the feature differences have
absolute value at most \(X^{\theta_R}\). For \(r=0\), take their
bound to be one. Suppose there are common real affine functions
\(U,V\) such that
\[
                    \Phi([a:c])=[U(z):V(z)]
 \tag{PF1}
\]
and the pair on the right is nonzero on every assigned vector.
Then the number of rows is
\[
                 O_{R,L,D,B,\varepsilon}
                 \left(X^{(1/3-\delta_R)L+\varepsilon}\right).
 \tag{PF2}
\]
The exponent and small-budget threshold are independent of \(D,B\).
The height threshold and implied constant are uniform under their
bounds at fixed \(L\). Thus the projective map may include a power
whose order grows as the budget tends to zero, with that order
fixed for each inner height limit.

### Rational coefficients of polynomial height

If \(r=0\), the observable is constant and yields at most \(D\)
rows. Indeed a primitive positive row is determined by its ratio,
and every fibre of \(\Phi\) has degree at most \(D\).
Assume \(r>0\). Normalize the features by choosing a witness
anchor and basis of \(r\) differences and a nonzero minor
\(\Delta_0\) on \(r\) coordinates. Give each anchored difference
its coordinates in that basis, multiplied by \(\Delta_0\).
The resulting features \(Z\in\mathbb Z^r\) have each coordinate
of absolute value at most \(X^{\theta_R}\), by Cramer's rule.
This is an affine isomorphism on the old feature span.
In every subfamily, the \(s\)-minors of anchored differences in
these coordinates have absolute value at most
\[
                            s!(2X^{\theta_R})^s .
 \tag{PF3}
\]
[Integer normalization of small feature minors (CF2)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#normalization-and-a-rank-drop-criterion)

We can choose \(U,V\) on the new coordinates with integer
coefficients of polynomial height in \(X\). To see this, write
each actual \(\Phi([a:c])=[p_i:q_i]\) primitively. Its height is at
most \(O_D(X^{B+D})\). Impose the homogeneous rational equations
\[
                         q_i U(Z_i)-p_i V(Z_i)=0
\]
on the \(2(r+1)\) affine coefficients. A basis of this kernel
consisting of integral vectors of polynomial height follows by
solving an independent subsystem by minors. Each condition
\(U(Z_i)=V(Z_i)=0\) excludes a proper subspace of the kernel,
because the original real pair is nonzero at every row.
There are at most \(X^2\) rows. Among combinations of the
integral basis with coefficients from
\(\{0,\ldots,\lfloor X^2\rfloor+1\}\), each excluded subspace
removes at most one coefficient choice after the other
coefficients have been specified along a transverse direction.
The grid therefore contains a combination avoiding all exclusions.
The polynomial height bound for the new forms, and for their
evaluations at all the features, is uniform in the assignments
at fixed \(R,D,B\).

### Affine rank loss from an original prime packet

Let \(P,Q\in\mathbb Z[s,t]\) denote the two binary forms of
\(\Phi=[P:Q]\). The branches for the original roles \(a,b,c\)
are \(b_a=[0:1]\), \(b_b=[1:1]\), and \(b_c=[1:0]\), respectively.
For each role \(W\) put
\[
             L_W(Z)=Q(b_W)U(Z)-P(b_W)V(Z),
\]
where the displayed branch pairs are used to evaluate the
binary forms. This is an integer affine function.
If \(p^j\mid W\), the primitive pair \([a:c]\) reduces to
\(b_W\) modulo \(p^j\). In each case one original coordinate
is a unit, so \((a,c)\) is a unit scalar times the branch
pair modulo \(p^j\). Homogeneity gives the same assertion
for the image pair. Outside the prime support of
\(\operatorname{Res}(P,Q)\), this image pair has a unit
coordinate and reducing it to lowest terms does not change
the projective congruence. Thus
\[
                            p^j\mid L_W(Z).
 \tag{PF4}
\]
This implication holds at a pole of \(\Phi\), by using
the homogeneous branch pair. It continues to hold when the
presenting pair \(U(Z),V(Z)\) has nontrivial common content.

On a current family, choose one actual feature \(Z_0\) for
which \(L_W(Z_0)\ne0\). If there is no such feature, the
entire family lies above \(\Phi(b_W)\) and has at most \(D\)
rows. Let \(A_0\) be the nonzero integer
\(\operatorname{Res}(P,Q)L_W(Z_0)\). Suppose \(M\) is a fixed
full divisor of \(W\) on a subfamily and \((M,A_0)=1\).
Then \(M\) divides \(L_W\) at every point of that subfamily.
If its affine dimension equals that of the current family,
it spans the same affine space. Take an anchored basis from
the subfamily and any nonzero \(s\)-minor \(\Delta\) on its
difference vectors. Cramer's rule expresses \(\Delta\)
times the affine evaluation at \(Z_0\) as an integer linear
combination of evaluations at the basis points, by applying
the rule to the homogeneous coordinates \((1,Z)\).
In particular \(M\mid\Delta L_W(Z_0)\), hence \(M\mid\Delta\).
If \(M\) exceeds the bound in (PF3), the dimension must drop.
The witness is selected before the divisor record and the
same argument applies after arbitrary thinning to that record.

### Uniform saving and its order of limits

Suppose no uniform constants in (PF2) exist for some \(R\).
As in (CF1), select budgets \(L\to0\) with cap tolerances,
feature-minor exponents and normalized count deficits tending
to zero, taking height limits at each fixed budget. The
upper bound \(1/3+o(1)\) follows from the limiting one-third
count at the total cap \((1+o(1))L\). Take the balanced
anatomic and valuation-profile refinements of such an endpoint
sequence. The total nonatomic radical mass in the actual
scaled valuation degrees is at least \(1/2\). In particular
some role \(W\) has positive nonatomic mass in a fixed
compact continuity band \([a_0,b_0]\subset(0,\infty)\).
After removing neighborhoods of the finitely many large atoms,
a fixed positive mass remains in the band and the individual
normalized prime costs in the retained part can be made
arbitrarily small by taking late profiles. The refining
cohorts retain the endpoint exponent, and the retained
band lower bound holds on each row. The same argument and
height conventions are used to select full divisor packets
in (CF1).
[Limiting one-third count](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)
· [Mass diffusion at the balanced endpoint](ABC_WIDE_COFACTOR_CONGRUENCES.md#necessary-conditions-at-normalized-exponent-13)
· [Cheap full packets in positive nonatomic bands](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#limiting-setup-and-cost-of-fixed-prime-supports)

Choose a fixed \(\xi>0\) smaller than one quarter of the
remaining band mass. For any witness \(A_0\) selected above,
delete rows whose \(W\) has normalized radical mass at least
\(\xi\) on the prime support of \(A_0\). The uniform prescribed
pool count (EB5) bounds their normalized exponent by
\(1/3-\xi/4\) at sufficiently late budgets. The polynomial
height exponent of \(A_0\) is bounded at each fixed budget,
even when \(D,B\) grow in the outer limit. More explicitly,
at fixed \(L\) its prime support has
\(O_{R,D,B}(\log X/\log\log X)\) elements. The partition sum
for the zero-price pool in (EB3) is therefore bounded by
\(\log X/\sqrt{\log\log X}\) for sufficiently large \(X\)
at that \(L\), uniformly in these witnesses. Consequently
the small-budget threshold used in the prescribed-pool
count does not require a uniform bound on \(D,B\) in the
outer limit.
[Effective count for prescribed prime pools (EB3–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

Choose a fixed \(p_0>0\) so small that
\(Rp_0<\min(\xi/8,1/12)\). Refine the retained automatic
band selections so that every individual prime cost is
less than \(p_0/2\). On each remaining row a product of
whole prime powers coprime to \(A_0\) can then be collected
in the band with normalized radical cost between
\(p_0/2\) and \(p_0\). Its numerical weight is at least
\(a_0p_0/2\). Fixing the product costs at most
\(n^{p_0+o_X(1)}\), by enumerating its squarefree support
and then its exponents at polynomial height. Retain a
largest record. Since the feature exponent tends to zero,
its fixed divisor exceeds every bound in (PF3) for large
\(X\) in all sufficiently late budgets. The affine
dimension drops by (PF4).

Repeat with the same branch and a new witness on each
stage. The support deletions cost a strict exponent margin
relative to the retained count, since there are at most
\(R\) prior locks. A family of constant image has at most
\(D\) elements, so it cannot have the retained exponent
\(1/3-Rp_0-o(1)>0\). After at most \(R\) drops the
features are constant. This contradiction gives (PF2).
The argument uses full original prime powers from a
nonatomic band and does not require a prime of positive
normalized individual cost, complete dictionary bounds,
or a prescribed minimum nonatomic mass profile.

## Rational descriptions and minimally repaired quotients

**PF5.** In (PF1) the right side can be any prescribed
pair of real polynomials of bounded degree \(d\) in
\(k\) integer parameters of absolute value at most \(X^u\),
provided \(u>0\) is sufficiently small in terms of \(k,d\).
Indeed use all nonconstant monomials of degree at most \(d\)
as the integer features. With \(R=\binom{k+d}{d}-1\)
each affine difference minor is at most
\(R!2^R X^{dRu}\). The coefficients on the feature side
may be arbitrary real numbers. At fixed \(R\), a
dictionary of at most \(n^{\delta_R/2}\) allowable feature
descriptions can be summed with a strict saving.
The projective feature condition is imposed on the actual
primitive ratio or its image. Merely prescribing a short
feature model on the repaired points does not provide this
condition.

Let \(M=\lceil K/L\rceil\) for a fixed \(K>0\). Count in
the minimally repaired primitive quotient model, with
least repaired target \(H_M(c)=\prod_{p\mid c}
p^{\max\{v_p(c),M\}}\le Y\). Assume the three separate
upper radical caps \(c^{(1/3+\xi)L}\), the target lower
radical cap \(c^{(1/3-\xi)L}\), and a feature description
of the primitive ratio as above, now with feature minors
at most \(Y^{\theta'}\) and map coefficient height at
most \(Y^{B'}\). For sufficiently small \(\xi,\theta'>0\)
in terms of \(R,K\), the number of these primitive quotients is
\[
                      O_{K,R,L,D,B',\varepsilon}
                      \left(Y^{(1-\sigma)/M+\varepsilon}\right)
 \tag{PF6}
\]
for some \(\sigma>0\) independent of \(L\).
There is no multiplicity from arbitrary common multipliers
because the count is in primitive originals and the least
repair is unique.

In fact \(H_M(c)\ge\operatorname{rad}(c)^M\) forces
\(c\le Y^{1/[M(1/3-\xi)L]}\).
For \(c\le Y^{1/[4M(1+3\xi)L]}\), the elementary radical
enumeration bound \(O(Z^{(1+3\xi)L+\varepsilon})\) up to
height \(Z\) gives at most \(Y^{1/(4M)+\varepsilon}\)
rows. Normalize the features of the entire \(Y\)-family before
passing to dyadic original height shells. By (PF3), every
\(j\)-minor in any shell is at most \(j!(2Y^{\theta'})^j\), even
if restricting to the shell lowers the affine rank. For
\(X/2<c\le X\) above the cutoff, these bounds and the map
coefficient powers of \(Y\) fit (PF1), since
\(X\ge Y^{1/[4M(1+3\xi)L]}\) and \(ML\le K+L\).
Taking \(\xi\) below the cap tolerance and the exponent
saving gives (PF6) by summing shells from (PF2), with
\((1/3-\delta_R)/(1/3-\xi)<1\). The shell height allowances
may be taken as small fixed multiples of the requested
\(Y\)-allowance, since \(ML\ge K\).
[Least target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

The branch calculation uses the same three distinguished points
as the Campana formulation of powerful triples. The counting input
is the effective radical budget on the original rows and the
determinants arise from exact linear dependencies in the lifted
features. No fixed generalized Fermat orders are required.
[Browning–Verzobio, Campana formulation](https://arxiv.org/html/2608.24512v1#S1)
· [Rational features inside a fixed divisor record](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md#lifting-short-rational-features-from-two-congruences)
