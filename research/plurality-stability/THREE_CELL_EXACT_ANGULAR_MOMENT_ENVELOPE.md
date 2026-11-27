# Exact first angular moment comparison for arbitrary circular masses

For two ordinary circular partitions of prescribed mass
vectors $p,q\in\Delta_3$, the maximum first Fourier pairing
is exactly the maximum under the three moment norm caps
and the closure constraints. It can be computed from two
positive $2\times2$ matrices. An optimal pair can be
realized by two single arcs and their complement in each
circle. The complement can be disconnected. This strengthens
the [first-mode angular correction](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md#an-explicit-correction-for-angular-moment-closure)
at finite amplitudes. It concerns the first angular mode.

## Exact support functions for closing triples

Use circular probability measure and put
$a_i=\sin(\pi p_i)/\pi$, $b_i=\sin(\pi q_i)/\pi$.
For each moment vector one has $|u_i|\le a_i$: project
in its own direction and select the largest cosine values.
Three moment vectors from a partition sum to zero.
For putative lengths $\ell$ of a closing triple let
$G(\ell)$ be its Gram matrix, as in (DG26), restricted
to $H=\mathbf1^\perp$. Write
$\mathcal G_a=\{A\succeq0:\widetilde A_{ii}\le a_i^2\}$,
where $\widetilde A$ is the extension by zero on
$\mathbf1$. This represents exactly the closing triples
with the norm caps in a plane, allowing singular Gram
matrices. The cap lengths themselves form a triangle,
as seen by taking the moments of consecutive arcs.

For a positive definite $X$ on $H$ let
$s_a(X)=\max_{A\in\mathcal G_a}\operatorname{tr}(XA)$.
For a closing triple write
$\operatorname{tr}(XG(\ell))=\sum_i\lambda_i(X)\ell_i^2$,
using the coefficients of (DG32). If the matrix of $X$ in
the basis $(1,-1,0)/\sqrt2,(1,1,-2)/\sqrt6$ is
$\left(\begin{smallmatrix}x&y\\y&z\end{smallmatrix}\right)$,
the coefficients are $x+\sqrt3y$, $x-\sqrt3y$ and
$(3z-x)/2$. At most one is nonpositive: the first
two sum to $2x>0$, and if $x±\sqrt3y\le0$, then
$z>y^2/x\ge x/3$, so the last is positive.
The exact support is
\[
 s_a(X)=
 \begin{cases}
   \sum_i\lambda_i(X)a_i^2,
                      &\lambda_i(X)\ge0\text{ for all }i,\\
   \lambda_j(X)a_j^2+\lambda_k(X)a_k^2+
                 \lambda_i(X)(a_j-a_k)^2,
                      &\lambda_i(X)<0 .
 \end{cases}                                                  \tag{AM1}
\]
For the second case put $\lambda_i=-\gamma$. At fixed
other lengths the best choice makes the third length
$|\ell_j-\ell_k|$. The resulting expression
$(\lambda_j-\gamma)\ell_j^2+(\lambda_k-\gamma)\ell_k^2+
2\gamma\ell_j\ell_k$ is increasing in both nonnegative
lengths, since $\lambda_j-\gamma$ and $\lambda_k-\gamma$
are positive by testing $X$ on $e_i-e_j$ and $e_i-e_k$.
Saturating the two caps is feasible because
$|a_j-a_k|\le a_i$. This proves (AM1).

In particular, for interior masses every SPD support
face is physically realized by circular partitions.
With strictly positive coefficients all caps are
saturated and consecutive arcs realize the Gram matrix.
With exactly one zero coefficient, the other two caps
are saturated and the remaining length may range from
$|a_j-a_k|$ to $a_i$. Take the two saturated arcs of
the prescribed masses and vary the distance of their
centers from $\pi(p_j+p_k)$ to $\pi$. They remain
disjoint. The first moment length of their complement
ranges continuously from $a_i$ to $|a_j-a_k|$.
The negative coefficient case uses the antipodal
endpoint of the same construction.

## The exact angular maximum and its two-dimensional dual

For Gram matrices $A,B$ on $H$ set
$F(A,B)=\operatorname{tr}[(A^{1/2}BA^{1/2})^{1/2}]$.
At fixed Gram matrices, maximizing the cross moment
over relative orthogonal orientations gives $F(A,B)$:
factor the two moment matrices and use the singular value
decomposition of their cross product. Thus the relaxed
closing-triple maximum is $\max_{\mathcal G_a\times\mathcal G_b} F$.
For interior $p,q$ this relaxed maximum is attained
by actual circular partitions, as follows.

Add $\epsilon(\log\det A+\log\det B)$ to $F$ and
maximize on the two compact cap regions, giving positive
definite maximizers $A_\epsilon,B_\epsilon$ since
the cap Gram matrices are positive definite. The elementary
matrix variational identity
\[
 F(A,B)=\tfrac12\inf_{X\succ0}
                  [\operatorname{tr}(XA)+\operatorname{tr}(X^{-1}B)]
\]
follows on positive definite matrices by completing the
square, with unique optimizer $XAX=B$. It extends
to singular matrices by approximation. At the regularized
pair the partial derivatives of $F$ are $X/2,X^{-1}/2$,
as follows by differentiation of the variational value.
The first-order support inequalities thus put the two
regularized Gram matrices in the support faces with
positive definite normals $X+2\epsilon A_\epsilon^{-1}$
and $X^{-1}+2\epsilon B_\epsilon^{-1}$. Both are
realizable by the arc constructions of (AM1).
As $\epsilon\downarrow0$, their unregularized values
approach the maximum of $F$: any feasible maximizer
can first be mixed slightly with the full-rank cap
Grams before comparing the regularized values.
The arc constructions and relative orientations are
compact, so a limiting pair attains the relaxed maximum.

Let $\mathfrak F(p,q)$ denote the actual maximum
pairing. The same argument shows
\[
 \mathfrak F(p,q)=\max_{\mathcal G_a\times\mathcal G_b}F
                  =\tfrac12\inf_{X\succ0}[s_a(X)+s_b(X^{-1})].
                                                               \tag{AM2}
\]
For the last equality, the upper-bound direction is
immediate from the variational identity. At $X$ for the
regularized pair, its support optimality gives
$s_a(X)-\operatorname{tr}(XA_\epsilon)\le4\epsilon$,
by moving the regularizing trace term to the other side
and using $\operatorname{tr}(A_\epsilon^{-1}A)\ge0$
for all feasible $A$. The inverse term has the same
bound. Hence the infimum differs from the relaxed value
by at most $4\epsilon$, proving equality.

These formulas also hold at boundary masses. For example
replace $p,q$ by their $(1-\epsilon)$ mixtures with the
uniform masses. Concavity of sine shows the new caps
contain $(1-\epsilon)a,(1-\epsilon)b$. This gives
approximation of every feasible boundary Gram and
upper semicontinuity follows by compactness. The arc
realizers converge after a subsequence, even when some
arc lengths shrink to zero. For the dual, use the
same cap inclusion to bound the boundary infimum by
$(1-\epsilon)^{-2}$ times the interior infimum
and take the limit. The variational bound in the
other direction holds directly.

For interior masses, scaling $X$ first reduces (AM2)
to the following minimization on a disk:
\[
 \mathfrak F(p,q)=\min_{x^2+y^2<1}
 \sqrt{\frac{s_a(X)s_b(\operatorname{adj}X)}{1-x^2-y^2}},
 \qquad X=\begin{pmatrix}1+x&y\\y&1-x\end{pmatrix}.           \tag{AM3}
\]
Indeed $X^{-1}=\operatorname{adj}X/\det X$ and the
scalar multiple of $X$ is optimized explicitly.
The objective diverges at the circle boundary since
the two cap Grams are positive definite and both
$X$ and $\operatorname{adj}X$ have trace two.
The six affine sign lines for
$\lambda_i(X),\lambda_i(\operatorname{adj}X)$ give
piecewise products of affine numerators. There is
therefore no residual optimization over measurable
circular sets in the exact first-mode penalty.

The penalty $\sum_i a_i b_i-\mathfrak F(p,q)$ is
globally valid in place of $\Gamma$ or $\Gamma_\sharp$
in (DG29) and (DG39). Those proofs separate a cosine
term with a nonnegative coefficient and use aligned
arcs for the remainder. The exact angular penalty
does not change the quadratic coefficient near
equal masses, where all norm coefficients of the
full-cap dual are positive. It strengthens the
finite-amplitude comparison.

### An interior range with a closed cap formula

Put $p_{\min}=\arccos(4/5)/\pi$. If all
coordinates of $p,q$ belong to $[p_{\min},1/2]$,
then $\mathfrak F(p,q)=F(G(a),G(b))$. Indeed for
the cap Gram $A$,
$\lambda_i(A)=3a_i^2-\tfrac12\sum_j a_j^2$ and
$\lambda_i(\operatorname{adj}A)=\tfrac32(\sum_j a_j^2-2a_i^2)$
by the $2\times2$ coefficient formulas. For the first,
the two other $\sin^2$ terms sum to at most
$1+\cos(\pi p_i)$ by the cosine addition identity.
The inequality $5\sin^2(\pi p_i)\ge1+\cos(\pi p_i)$
follows from $\cos(\pi p_i)\le4/5$. The second is
nonnegative by the cosine law
in the sine-length triangle, whose angles are at
most $\pi/2$. The same holds for $B$. The full-cap
dual solution $X=(B+hA^{-1})/F(A,B)$, with
$h=\sqrt{\det A\det B}$, and its inverse have
nonnegative coefficients, so (AM1) proves the cap
formula throughout this range.

In particular if $|p_i-1/3|,|q_i-1/3|\le1/20$,
one has
\[
 \sum_i a_i b_i-\mathfrak F(p,q)\ge
                      \frac{3}{40}|p-q|^2.                 \tag{AM5}
\]
To prove this, use consecutive arcs of the same
cyclic order in both circles. The best relative
orthogonal alignment is then a rotation, since
the closing moment triples have the same orientation.
If $w_i=a_i b_i$, the optimized first-mode pairing
is the magnitude of the sum of complex phases
weighted by $w_i$. The differences of the matched
center phases for labels $i,j$ are
$\pm\pi(p_k-q_k)$ modulo $2\pi$ with $k$ the third
label, because the two arc centers are separated
by $\pi(p_i+p_j)$ along their common endpoint.
Thus the deficit is at least
\[
 \frac{\sum_{i<j}w_iw_j(1-\cos(\pi(p_k-q_k)))}{\sum_i w_i}.
\]
Here $w_i\ge\sin^2(17\pi/60)/\pi^2>.061$ and
$\sum_i w_i\le9/(4\pi^2)<.229$ by Cauchy–Schwarz
and (JA7). For $|d|\le1/10$,
$1-\cos(\pi d)\ge4.88\,d^2$ by the quartic Taylor
bound. These elementary inequalities give (AM5).

### The optimum need not saturate all the caps

Take $p=(5/12,5/12,1/6)$ and $q=(1/6,1/3,1/2)$.
In the basis used in (AM1), and multiplying the
Gram matrices by $\pi^2$, the full-cap matrices are
\[
 A=\begin{pmatrix}(7+4\sqrt3)/8&0\\0&3/8\end{pmatrix},
 \qquad
 B=\begin{pmatrix}1/2&-\sqrt3/4\\-\sqrt3/4&3/2\end{pmatrix}.
\]
The inverse Procrustes matrix is
$X^{-1}=(A+\sqrt{\det A\det B}\,B^{-1})/F(A,B)$.
Here $\lambda_3(B^{-1})=0$ and
$\lambda_3(A)=(1-2\sqrt3)/8<0$. Thus reducing the
third squared length of $B$ slightly, with the
other two at their caps, increases $F$ by the
derivative formula above. This is feasible since $B$
is positive definite. It is realized by two
disjoint saturated arcs and a split third cell.
Any pair in which each circle is divided into
three consecutive arcs has the original full-cap
Gram matrices and strictly loses in this example.

## Compatibility among three radii

The exact first-mode envelope still discards a constraint:
the same slice participates in every pairing. Let the
mass vectors of a collection of slices lie in a sufficiently
small neighborhood of a fixed interior profile. Then
the pair optimum is the full-cap value
$F_{rs}=F(A_r,A_s)$, where $A_r=G(a(p_r))$,
by positivity of the six dual coefficients around the
diagonal. Let $X_{rs}A_rX_{rs}=A_s$ be the positive
definite transport. For an actual slice let $U_r$ be
its first moment matrix on $H$, and put
$E_{rs}=F_{rs}-\operatorname{tr}(U_r^{\mathsf T}U_s)$.
There is a uniform bound in that neighborhood
\[
 E_{rs}+E_{st}+E_{rt}\ge
             c\|X_{rt}-X_{rs}X_{st}\|_F^2.                 \tag{AM4}
\]
It follows from the exact pair remainder,
\[
\begin{split}
 2E_{rs}={}&
 \|U_r X_{rs}^{1/2}-U_s X_{rs}^{-1/2}\|_F^2\\
 &+\operatorname{tr}(X_{rs}(A_r-U_r^{\mathsf T}U_r))
   +\operatorname{tr}(X_{rs}^{-1}(A_s-U_s^{\mathsf T}U_s)).
\end{split}
\]
The trace terms bound the sum of individual cap
deficits below by a positive constant multiple.
The squared remainder controls
$\|U_s-U_rX_{rs}\|_F^2$. Compose two such relations
and subtract the third. If $U_r$ is well conditioned
this yields (AM4). Otherwise its Gram matrix stays
a fixed distance from the uniformly positive $A_r$;
the individual cap deficits then give a fixed
positive lower bound on the left. The transports
are uniformly bounded, so (AM4) also holds in
this case.

Since $X_{rt}$ is symmetric, the right side also
controls the square of the commutator of $X_{rs}$
and $X_{st}$. For a compact radial annulus on which
the profiles stay in the neighborhood, integrating
the three-point bound controls a further nonnegative
correction to the polar comparison. In fact the
weight $V_\rho$ of the separated first mode is
bounded below on that annulus. Integration in a
third radius converts the left side of (AM4) into
three times the annular measure times the pair loss.
For $p_r=\mathbf1/3+\epsilon b_r$ with bounded
vectors $b_r\in H$, the commutator has leading
order $\epsilon^2$ proportional to the oriented
area of the triangle $(b_r,b_s,b_t)$. To check
this, differentiate $X_{rs}A_rX_{rs}=A_s$ at the
scalar equal-cap Gram $A_*$: the first derivative
is $(A'_s-A'_r)/(2A_*)$ in scalar normalization.
The derivative of the cap Gram identifies the
two-dimensional mass plane with the traceless
symmetric matrices on $H$, since its trace
derivative is zero and the diagonal cap derivatives
determine the Gram. The commutator of two such
matrices is a nonzero alternating bilinear form.
Thus profiles with a noncollinear triple of
directions pay a quartic compatibility cost beyond
the independently optimized cosine pairs.

The simultaneous equality condition in this
neighborhood can also be stated exactly. If all
pair losses vanish almost everywhere, fix a reference
radius by Fubini. Its transport maps to the other
radii commute, because their ratios are the symmetric
pair transports. They can be diagonalized in
one orthogonal basis of $H$. In this basis the cap
Grams are diagonal congruences of one fixed Gram
matrix. Conversely, a family of this form can
attain all pair cosine optima: use one full-cap
moment matrix, multiply by the corresponding positive
diagonal maps in that basis, and choose actual
consecutive arcs with those moments. This restricts
simultaneously optimal mass profiles to a
one-dimensional family up to the overall Gram scale.
