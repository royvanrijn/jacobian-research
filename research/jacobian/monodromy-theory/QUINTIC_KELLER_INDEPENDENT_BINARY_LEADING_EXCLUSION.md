# Independent binary leading parts of quintic Keller maps

**Theorem QuinticKellerIndependentBinaryLeadingExclusion.** A
three-variable polynomial Keller map of degree five in characteristic
zero cannot have three linearly independent degree-five leading forms
in two linear source coordinates.

Over \(\mathbb C\), put the leading part \(P=P(x,y)\) in binary
coordinates and remove its common factor \(g\), writing \(P=g\Phi\)
with \(\Phi\) a coprime homogeneous tuple of degree \(e\).
Independence implies \(e\ge2\). For \(e=2\), a constant target
change makes \(\Phi=(x^2,xy,y^2)\), which is excluded by the
[complete Veronese theorem](QUINTIC_VERONESE_BINARY_RANK_EXCLUSIONS.md).
Below \(e\ge3\). The characteristic-zero assertion follows by
embedding the finitely generated coefficient field in \(\mathbb C\).

## Mutual tangent lines

Set \(n=P_x\times P_y\). There is no nonzero constant vector \(v\)
with \(n\cdot v=0\). Indeed a constant target change making \(v\)
the third target direction would make the planar Jacobian of two
independent homogeneous binary forms of degree five vanish. By
Euler's identity a vanishing Jacobian would make their ratio
constant. This is Wang and Yan's
[Lemma 3.4](https://arxiv.org/pdf/2610.36495v1).

Let \(W\ne0\) be a homogeneous binary tuple of positive degree.
Assume both
\[
                  n\cdot W=0,\qquad P\cdot(W_x\times W_y)=0.
                  \tag{1}
\]
Then \(W\) is a rational multiple of \(P\). For a proof, work on
\(x\ne0\) using \(t=y/x\) and dehomogenize by \(x=1\). The first
equation says that \(W\) belongs to \(\operatorname{span}(P,P')\)
over \(\mathbb C(t)\), and the second that \(P\) belongs to
\(\operatorname{span}(W,W')\). Here primes are \(t\) derivatives.
Neither pair of a column and its derivative can be dependent.
For \(P\) this follows from independence of its polynomial
components. For \(W\), dependence would make its projective
direction a constant vector, contrary to \(n\cdot W=0\).
If \(P,W\) themselves were independent over \(\mathbb C(t)\), the
two planes would coincide. Then \(P'\) and \(W'\) both belong to
their span, so \((P\times W)'\) is proportional to \(P\times W\).
The ratios of the entries of that cross product are constant.
Thus \(P\) lies in a constant target plane, contrary to the
hypothesis. This proves the assertion. Since \(\Phi\) is
coprime, any polynomial tuple proportional to \(P\) is a
polynomial multiple of \(\Phi\).

## The high determinant coefficients

Write \(H_j\) for the ordinary homogeneous pieces of the map
and \(\Delta_d\) for the degree-\(d\) piece of its Jacobian
determinant. At \(\Delta_{11}\) the mixed determinant equals
\(n\cdot\partial_z H_4\), so each positive-\(z\) coefficient
of \(H_4\) is tangent to \(P\). In particular \([z^4]H_4=0\).
For \(W\) homogeneous binary of degree \(l>0\), Euler's
identity gives
\[
 \det(P_x,W_y,W)+\det(W_x,P_y,W)
                       =-\frac{5}{l} P\cdot(W_x\times W_y).
                       \tag{2}
\]
Indeed substitute \(W=(xW_x+yW_y)/l\) and then
\(xP_x+yP_y=5P\).

The \(z^5\) coefficient of \(\Delta_{10}\) is three times (2)
for \(W=[z^3]H_4\), of degree one. A rank-one \(W\) contradicts
the no-constant-annihilator condition, and a rank-two \(W\)
has nonzero constant \(W_x\times W_y\). The vanishing of (2)
would make the components of \(P\) dependent. Thus \(W=0\).
Now the \(z^3\) coefficient of \(\Delta_{10}\) is twice (2)
for \(W=[z^2]H_4\), of degree two. In both comparisons the
\((5,5,3)\) term has \(z\)-degree at most two. By tangency and
(1), the degree-two \(W\) is a polynomial multiple of
\(\Phi\), which is impossible unless it is zero. With
\(H_4\) affine in \(z\), the \(z^2\) coefficient of
\(\Delta_{10}\) is \(3n\cdot[z^3]H_3\), so the latter constant
column is zero. Consequently the full map has the form
\[
                         F=C+zB+z^2A,\qquad
                         \deg B\le3,\quad \deg A\le1,
                         \quad C_5=P.
\]

## Completing the affine-coefficient analysis

The
[binary Veronese reduction](QUINTIC_VERONESE_BINARY_Z_REDUCTION.md#the-affine-in-z-alternatives)
excludes \(A=0\) under precisely the properties used here:
independence of \(C_5\), degree at most three for \(B\), and
tangency of \(B_3\) to a normal \(n\) without a constant
annihilator. For rank three of \(B\), that proof uses the
degree-three homogeneous factorization and the absence of a
common affine zero to make \(B\) univariate; tangency then drops
its degree to at most two. The arguments for ranks one and two
do not use tangency.

For nonzero \(A\), the same reduction uses only independence
of \(C_5\) and the affine degree of \(A\) to leave \(A=(x,1,0)\)
or \(A=(0,0,h)\) up to the indicated coordinate changes.
The
[rank-two exclusion](QUINTIC_VERONESE_BINARY_RANK_EXCLUSIONS.md#the-rank-two-coefficient)
uses only these degree constraints and the resulting
\(B_3\) coordinate's independence of \(y\) (here the subscript
3 indexes the third component, not degree). For nonconstant
\(h\), the
[rank-one fibre argument](QUINTIC_VERONESE_BINARY_RANK_EXCLUSIONS.md#leading-linear-coefficients-in-the-remaining-rank-one-case)
and the binary exclusion that follows it use only the degree
bounds and independence of \(C_5\) once \(A=(0,0,x)\).
Their determinant and rational-pencil steps apply unchanged.

Finally if \(A\) is constant, let \(W=B_3\) denote the
homogeneous cubic part. The \(z\) coefficient of
\(\Delta_{10}\) gives (2) equal to zero, because
\([z^2]H_3=0\). Thus (1) makes \(B_3\) a polynomial
multiple of \(\Phi\). If \(B_3\ne0\), its degree is three,
\(e=3\) and its entries are independent. After making
\(A=(0,0,h)\) with nonzero constant \(h\), the \(z^3\)
determinant coefficient \(2h\{B_1,B_2\}\), with subscripts
indexing entries, excludes independent leading cubics of
\(B_1,B_2\). Hence \(\deg B\le2\). The final
[constant rank-one argument](QUINTIC_VERONESE_BINARY_RANK_EXCLUSIONS.md#constant-rank-one-coefficient)
then applies: depressing the third component preserves the
independent quintic leading forms while the \(z\)-linear
determinant coefficient has nonzero degree-eight top.
This finishes all ranks.
