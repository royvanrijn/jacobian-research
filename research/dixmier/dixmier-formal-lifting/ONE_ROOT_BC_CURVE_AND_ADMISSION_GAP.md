# A one-place BC-type curve retaining the homogeneous parameter

The homogeneous parameter in the
[one-root ODE control](ONE_ROOT_CHOPPING_ODE_HOMOGENEOUS_BRANCH.md)
also survives an explicit monic algebraic-curve and Puiseux-trace test.
This supplies curve data, not a reconstruction of commuting ordinary
differential operators or an admitted polynomial Weyl string pair.
It therefore narrows the remaining constraint without resolving the
source's main theorem.

The primary comparison is Zheglov
[arXiv:2410.06959v5](https://arxiv.org/html/2410.06959v5), especially
Definition 5.2, Remark 5.3, and Section 6.2's scalar curve/trace step.
The exact [cached PDF](sources/zheglov-2410.06959v5.pdf), dated
19 January 2026, has 78 pages and SHA256
`897913ef95e011332e86e348803c6c9da5dbd196db17a04f8f951110e9b095e3`.
No earlier source statement is adopted to produce actual operators.

## The explicit curve and simple formal branch

Let \(k\) have characteristic zero and let \(c\in k^\times\).
In commuting variables, put
\[
 Z=Y^2-X^3,\qquad
       F_c(X,Y)=(Z^2-c^2Y)^2-cX^2Z.                       \tag{1}
\]
The polynomial \(F_c\) is monic of degree eight in \(Y\), has zero
\(Y^7\)-coefficient, and has leading weighted form
\((Y^2-X^3)^4\) for weights \(\operatorname{wt}(X,Y)=(8,12)\).
Thus its displayed leading shape is of BC type with orders \((8,12)\).

There is an exact formal branch with \(X=t^8\) and
\[
 Y=t^{12}+\frac c2t^{-6}+\frac14t^{-7}
                   +0t^{-8}-\frac1{128c^2}t^{-9}+O(t^{-10}).
                                                               \tag{2}
\]
This is an all-orders existence assertion for a formal algebraic
branch, not a claim that the finite truncation in (2) solves (1).

To prove it, set \(z=t^{-1}\) and substitute
\[
             X=z^{-8},\qquad
             Y=z^{-12}+\frac c2z^6+Rz^7.
\]
Then \(G(z,R)=z^{22}F_c(X,Y)\) is a polynomial. Its first
coefficients are exactly
\[
 G(z,R)=c^2(16R^2-1)
       +z(32cR^3-2cR)+16R^4z^2+O(z^3).                 \tag{3}
\]
At \((z,R)=(0,1/4)\),
\(\partial_RG=8c^2\ne0\). Solving successive coefficients gives
a unique \(R(z)\in k[[z]]\) with \(R(0)=1/4\). This elementary
formal implicit-function argument works by dividing at each stage
by the fixed nonzero derivative. The first two new coefficients are
\(R_1=0\), \(R_2=-1/(128c^2)\), which proves (2).

## One place at infinity and the trace coefficient

After adjoining eighth roots of unity to the constants, the extension
\(k((t^{-1}))/k((t^{-8}))\) is cyclic of degree eight. The coefficient
of \(t^{-7}\) in (2) is nonzero. Since \(\gcd(7,8)=1\), the
eight series \(Y(\zeta t)\) are distinct. Each solves the same
polynomial (1). Therefore they are all of its roots, and (1) is
irreducible over the completed base field as well as over \(k(X)\).

Equivalently the degree-eight finite projection to the \(X\)-line
has one place above infinity on the normalized curve; \(t^{-1}\)
is a local uniformizer after completion. The constant extension
causes no obstruction to these conclusions over the original field.

Vieta's formula and the zero \(Y^7\)-coefficient give
\[
                         \sum_{\zeta^8=1}Y(\zeta t)=0.
\]
The coefficient of \(t^{-8}\) in this sum is eight times that
coefficient of \(Y\). Hence
\[
                            [t^{-8}]Y=0.                 \tag{4}
\]
This verifies the trace constraint directly for the entire branch,
in addition to its explicit coefficient check in (2).

## Exact chopping heads

On this branch, define ordinary commuting Laurent series
\[
               Q_1=Y^2-X^3,\qquad Q_2=Q_1^4-c^4X^3.
\]
Their heads and the next relevant coefficients are
\[
 \begin{split}
 Q_1&=ct^6+\tfrac12t^5+0t^4+\text{lower},\\
 Q_2&=2c^3t^{23}+\tfrac32c^2t^{22}+\text{lower}.
 \end{split}                                                   \tag{5}
\]
Thus the same nonzero homogeneous parameter occurs with orders
\((p,q,q_1,q_2)=(8,12,6,23)\), sequential gcds \(4,2,1\), and
chopping distances \(18,1\), whose sum is \(p+q-1=19\).
These are actual algebraic Laurent series on (1), rather than
independently prescribed coefficient tuples.

## The remaining admission boundary

Equations (1)--(5) show that the monic BC-type shape, one place at
infinity, the displayed head sequence, and the zero trace coefficient
do not by themselves force \(c=0\). The curve is a valid concrete
test for any proposed argument using only those data.

No commuting ordinary differential operators with this spectral data
are reconstructed here. A chosen admissible sheaf/Sato datum, its
ordinary-operator reconstruction, compatibility with the required
string normal form, and simultaneous global polynomial Weyl admission
remain separate obligations. In particular (1) is not a DC-pair or
a counterexample to a source theorem under all of its operator
hypotheses. The main DC1 claim remains unverified.

[The independent exact checker](check_one_root_bc_curve.py) verifies
the polynomial normalization, simple Hensel root and first coefficients,
monicity/weighted top, and the heads (5), with \(c\) indeterminate.
The all-orders branch and trace statements follow from the written
proof. [The retained receipt](one_root_bc_curve_controls.json) pins
this note, checker and exact primary PDF. No formal verification or external review is claimed.
