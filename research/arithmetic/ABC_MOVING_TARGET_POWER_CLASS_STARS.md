# Moving-target stars from a joint power class

Several primitive rows with an omitted-power representation give
a split cyclic cover in one rational parameter after their reduced
root ratios are fixed. Only one weighted product of actual target
values is charged. The completion works at arbitrarily large
omitted orders divisible by the class order, and keeps all original
representations even when left labels repeat.

## Actual rows and relative right kernels

Fix \(m\ge2\), an integer \(e\ge2\), a prime \(q\mid e\), and
a sign \(\sigma\in\{-1,1\}\). Use distinct actual positive
primitive rows
\[
        U_i=P_i y_i^e,\qquad U_i+\sigma B_i=T_i>0,\qquad
                        \gcd(U_i,B_i)=1 .
\]
The \(P_i\) are prescribed positive integers and may repeat.
All complete values are at most \(X^{B_0}\) for fixed \(B_0>0\).
Define
\[
       g=\gcd_i y_i,\qquad y_i=g u_i,\qquad
       G=\gcd_i B_i,\qquad B_i=G R_i x_i^e
\]
with \(R_i\) canonically \(e\)-free. Then \(\gcd_i u_i=1\)
and \(\gcd_i x_i=1\). For each \(i\) the same-row roots
\(u_i,x_i\) are coprime. Write \(y_i\le Y_i,\ x_i\le Z_i\);
the latter bounds are on the normalized right roots.
All bounds can be taken at most \(X^{B_0/e}\).

Prescribe \(D\ge1\) of polynomial height such that all
\(U_i,B_i\) are units modulo \(D\) and their ratios \(U_i/B_i\)
are congruent to one another. In particular \(D\mid T_i\)
for all \(i\) supplies this situation. The pool \(P\) is
prescribed of polynomial prime product and may depend on
the left labels and the modulus. For a nonzero rational \(t\)
write \(d_q(t;P)\) for its outside defect radical at order
\(q\), the product of \(p\notin P\) with \(v_p(t)\not\equiv0\bmod q\).

Let \(1\le w_i<q\) for \(1\le i\le m\), and \(0\le w_0<q\).
Put
\[
              \mathcal H=\operatorname{lcm}\left(
          \bigl\{d_e(B_i/B_1;P):2\le i\le m\bigr\},
                 d_q\left(G^{w_0}\prod_{i=1}^m T_i^{w_i};P\right)
                                               \right).
\]
The root ratio records use order \(e\), while the target-class
record can use the smaller prime order \(q\).

## A uniform star completion

**MS1.** Set
\[
 b=m+\mathbf1_{w_0\ne0}
              +\mathbf1_{w_0+\sum_i w_i\not\equiv0\bmod q}.
\]
If \((q-1)(b-2)/2\ge1\), the number of original ordered
\(m\)-tuples with \(\mathcal H\le H\) is
\[
 O_\varepsilon\!\left(
 H X^\varepsilon
       \prod_{i=2}^m \left(1+\frac{Y_iY_1 Z_iZ_1}{D}\right)
                                                          \right).
                                                               \tag{MS1}
\]
The constant is uniform in the numerical data at fixed
\(B_0,m,e,q\) and fixed polynomial-height bounds for the pool
and modulus. Further original radical restrictions are allowed.

Fix the union radical \(\mathcal H\), then each \(R_i\)
and the \(q\)-free representative of the weighted class at
subpower cost on prescribed prime supports. For the first
assertion, outside the pool and the relative right defect
every valuation \(v_p(B_i)\) is congruent modulo \(e\).
Subtracting their minimum shows that \(p\) does not occur in
any \(R_i\). The class representative is supported in the
specified union and pool by definition. All powers have
polynomial height at fixed \(m,e,q\).

The original projective congruence gives for \(i\ge2\)
\[
      \left(\frac{u_i x_1}{u_1 x_i}\right)^e
                   \equiv \frac{P_1 R_i}{P_i R_1}\pmod D .
\]
The unit root classes cost \(X^\varepsilon\) at fixed \(e\)
by the structure of units at prime powers. The reduced
positive fractions have numerator at most \(Y_i Z_1\)
and denominator at most \(Y_1 Z_i\). In one residue
their cross determinants are nonzero multiples of \(D\).
Rescale by the two bounds and use angular spacing in the
first quadrant. This gives
\(O(1+Y_iY_1 Z_iZ_1/D)\) fractions in that residue.
No assumption of pairwise coprimality between different
rows' reduced roots is needed.
[Angular residue count in MT1](ABC_MOVING_RIGHT_AND_TARGET_POWER_CLASS_CODEGREES.md#recovery-and-the-role-of-shared-target-primes)

Write the recovered fractions as \(r_i\), and set \(r_1=1\).
The root ratios \(u_i/x_i=\alpha r_i\) have one common
rational multiplier \(\alpha=u_1/x_1\). At a prime outside
the numerator and denominator supports of the \(r_i\),
every ratio has the same valuation \(v_p(\alpha)\).
If that valuation is positive it contradicts \(\gcd_i u_i=1\);
if negative it contradicts \(\gcd_i x_i=1\), since the
individual ratios are reduced. Thus \(\alpha\) has prescribed
support and polynomial height, and has only subpower many
possibilities. Fixing it recovers all \(u_i,x_i\) in lowest
terms. This handles any unknown gcds in the cross products.

Put \(a_i=P_i u_i^e,\ b_i=R_i x_i^e\) and
\(t=G/g^e\). The ratios \(a_i/b_i\) are pairwise distinct:
they differ from \(U_i/B_i\) by a common factor, and equality
of two actual primitive summand ratios recovers the same
original row. They are also nonzero. The recovered target
class is the class of
\[
                     t^{w_0}\prod_{i=1}^m(a_i+\sigma b_i t)^{w_i}.
                                                               \tag{MS2}
\]
Indeed \(T_i=g^e(a_i+\sigma b_i t)\), \(G=g^e t\), and
\(q\mid e\). At the fixed class coefficient \(c\) the
parameter \(t\) lies under a rational point of
\(z^q=c\,t^{w_0}\prod_i(a_i+\sigma b_i t)^{w_i}\).
This has \(b\) distinct rational branches, all completely
ramified. Its genus is \((q-1)(b-2)/2\). The usual split
cyclic descent and uniform Mordell–Lang argument applies
in genus at least two. If infinity is a branch, move it
to a finite point before applying the split model. If
necessary a single omitted parameter can be counted
directly. In genus one, the four-branch parity cover has
rational \(2\)-torsion and the three-branch cubic cover
has rational \(3\)-torsion in its split model. Their
coefficient and rational point heights are polynomial in
\(X\). Dujella's bound gives subpower completion of both.
These are exactly the bounded-degree split-cover
completions used in (IC1).
[Uniform split-cover count and elliptic models](ABC_SPLIT_INDEX_POWER_CORRESPONDENCES.md#a-joint-count-from-the-branch-set)
· [Dujella, Theorem 1.1](https://arxiv.org/html/2312.03655v1#S1)

Finally \(\gcd(G,g)=1\) by original primitivity. The reduced
rational value \(t=G/g^e\) fixes \(G,g\), and all exact
original rows follow. This proves (MS1). The proof has no
loss proportional to \(X^{O(e^{-3/2})}\). It trades the
relative target class of the pair argument for a weighted
product class with sufficiently many rational branches.

For any even \(e\), two rows and parity can use
\(w_0=w_1=w_2=1\), taking the class of \(G T_1T_2\).
The three finite roots \(0,-a_1/(\sigma b_1),-a_2/(\sigma b_2)\)
are distinct. For \(m=3\) and \(q=2\), the simple target
product already suffices with \(w_0=0\). For any odd
prime \(q\mid e\), two targets with \(w_1=w_2=1\) give
two finite branches plus infinity, with no right-gcd
class in the joint target record.

## An actual all-row graph estimate

Suppose there are \(A\) left labels \(P_i\) and \(B_1\) right
kernel labels \(Q\), with the right labels distinct modulo
rational \(e\)-th powers. The original rows have
\(B_i=Qz_i^e\) within each right column, and may have arbitrary
multiplicity at a kernel pair. In every such column tuple
the \(R_i\) above are one and \(G=Q(\gcd_i z_i)^e\). The
weighted class in (MS2) can consequently be charged by
\(d_q(Q^{w_0}\prod_i T_i^{w_i};P)\).

**MS3.** Suppose at least a fraction \(\rho>0\) of the
ordered distinct \(m\)-tuples of actual rows in common
right columns satisfy the hypotheses of (MS1) at uniform
cap \(H\) and with a prescribed common modulus \(D\).
Take \(Y_i=Y,Z_i=Z\) for uniform bounds and put
\[
                 \Delta=H X^\varepsilon
                          \left(1+\frac{(YZ)^2}{D}\right)^{m-1}.
\]
Then the total actual row count \(W\) satisfies
\[
                    W\ll_m B_1+
                        A B_1^{1-1/m}(\Delta/\rho)^{1/m}.   \tag{MS3}
\]
The constants in (MS1) are included. Every full ordered
left tuple, including repeated left labels, has at most
\(O(\Delta)\) completions over all columns. Indeed the
actual right values determine their unique kernel class
in the dictionary. There are at most \(A^m\) choices of
left labels. If \(W\ge4mB_1\), the columns with at least
\(2m\) original rows contain at least half the rows. For
their sizes \(d_j\), the sum of \((d_j)_m\) is
\(\gg_m W^m/B_1^{m-1}\) by the power mean. Apply the stated
eligible fraction. When \(W<4mB_1\) the first term suffices.

For \(n=X^L\), \(A,B_1\le n^{1/6}\), \(D\ge(YZ)^2\),
\(H\le n^\beta\) and \(\rho\ge n^{-\tau}\), the exponent
in base \(n\) is at most
\[
          \max\left\{\frac16,\,
             \frac13-\frac{1/6-\beta-\tau}{m}\right\}
\]
up to an arbitrarily small inner allowance at each fixed
budget. This saves a fixed power when \(\beta+\tau<1/6\).
The many-row theorem applies when a joint class is cheap
even if no comparable bound is supplied on the individual
target defects. It requires actual kernel dictionaries and
eligibility on the original stars. Extracting those conditions
from broad valuation profiles is a separate arithmetic step.
[Moving-target relative defects and all-pair bound](ABC_MOVING_RIGHT_AND_TARGET_POWER_CLASS_CODEGREES.md)
