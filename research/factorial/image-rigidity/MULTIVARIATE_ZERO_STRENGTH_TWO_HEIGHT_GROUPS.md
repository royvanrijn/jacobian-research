# Zero-strength exclusion for two affine factorial height groups

Fix an ordinary two-monomial core and adjacent orders \(n,N=n+1\)
with \(n\ge2\). Append two groups of ordinary exponent rows.
In group \(b=1,2\), row type \(i\) has smaller exponent
\(a_i+\alpha_i L_b\), positive integer gap \(d_i\),
orientation \(\epsilon_i\in\{1,-1\}\), and multiplicity
\(c_i r+s_i\). Here \(a_i,\alpha_i\) are rational with
\(\alpha_i>0\), \(c_i\ge0\) and \(s_i\) are fixed integers,
and \(L_b,r\) are integers on the ordinary family. All row
exponents and counts must be admissible. Assume
\[
 \sum_i\epsilon_i c_i d_i=0,\qquad
 h=\sum_i\epsilon_i s_i d_i,\qquad
 A_b=\sum_{i\in b}c_i d_i^2/\alpha_i>0,\qquad
 B_b=\sum_{i\in b}\epsilon_i c_i d_i/\alpha_i ,
\]
and \(A_1B_2-A_2B_1\ne0\).

**Theorem.** There is no sequence of exact ordinary collisions
with \(r\to\infty\), \(L_1,L_2\to\infty\), and the added
curvature tending to zero. In particular this covers two
independent opposed heights with any fixed gaps and fixed
signed gap. The core may have a consecutive collision at
the zero-strength physical limit. There is a uniform
exclusion neighborhood of zero strength for sufficiently
large heights and counts, with the template fixed. No
numerical cutoff is claimed.

The proof combines a number-field auxiliary polynomial at
zero strength with the
[Ax rank inequality](MULTIVARIATE_OPPOSED_AX_BRANCH_RANK.md#a-bound-from-differential-transcendence).
The first forces both height functions to be algebraic over
the count function on any accumulating branch. The second
then forces the common-factor equations to hold on the
entire torus specified by the branch's factorial characters.
Their leading logarithmic direction would make the core
collide along a nonconstant Gaussian clock, which is impossible.

## Finite branches over a number field

Put \(u=1/r\) and \(p_b=r/L_b\). For \(j=n,N\) and
\(0\le k\le j\), set
\[
 R_{i,j,k}(x)=
 \begin{cases}
   \displaystyle\sum_{\ell=1}^{kd_i}
       \log\bigl(1+(a_i+\ell/j)x/\alpha_i\bigr),
                   &\epsilon_i=1,\\
   \displaystyle-\sum_{\ell=(j-k)d_i+1}^{jd_i}
       \log\bigl(1+(a_i+\ell/j)x/\alpha_i\bigr),
                   &\epsilon_i=-1 .
 \end{cases}
\]
Use the logarithm at one. After stripping the index-zero
factorials and changing the argument in both moments by
\(\prod_i(\alpha_i L_{b(i)})^{-\epsilon_i d_i(c_i r+s_i)}\),
the resulting coefficient at \((j,k)\) is
\[
 a_{j,k} j^{hk}
    \exp\!\left[\sum_i(c_i/u+s_i)R_{i,j,k}(u p_{b(i)})\right],
                                                               \tag{1}
\]
where \(a_{j,k}\) are the rational constant-normalized core
coefficients. The signed balance of the \(c_i\) gives the
fixed factor \(j^{hk}\). The apparent \(1/u\) is removable.
Calculation of the linear terms of the logarithms gives the
Gaussian central pair at
\[
                T=A_1p_1+A_2p_2,\qquad \lambda=B_1p_1+B_2p_2.
\]
The term contributed by each offset \(a_i\) is linear in
\(k\) and independent of the order. It cancels from the
adjacent endpoint scale comparison. This is the grouped
case in the
[local height-branch theorem](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md#two-groups-of-affine-row-types).

At a hypothesized zero-strength sequence, \(p_b\to0\) for
both groups, since each \(A_b>0\) and the curvature contribution
from group \(b\) is asymptotic to \(A_b r/L_b\).
The clock variables \(w=\exp(-T/(2nN))\) and
\(v=\exp(-\lambda/(2nN))\) thus approach \((1,1)\).
The determinant assumption makes them valid holomorphic
coordinates in \(p_1,p_2\) near zero. The simultaneous
angular and physical resultant equations at \(u=0\) have
an isolated zero there if there is a contact, by the same
central quotient-resultant argument as at positive strength:
the angular eliminant is a nonzero polynomial in \(w\),
and at a fixed \(w\) the physical resultant has finitely
many possible scale roots \(v\). In the absence of a central
common zero the original continuous family already excludes
the proposed sequence.

Consequently the two Weierstrass preparations and the norm
from the [Puiseux branch proof](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md#finitely-many-puiseux-branches)
put an accumulating sequence on one of finitely many
analytic branches in \(s=u^{1/M}\), with some positive
integer \(M\). In this instance the Taylor coefficients in
(1) are rational at the center \((u,p_1,p_2)=(0,0,0)\).
The clock coordinates and their inverses have rational
Taylor coefficients at this center. The two formal
preparations and the norm therefore have rational coefficients. Successive formal Puiseux expansions
of the two prepared polynomials give branch coefficients
in a fixed number field. To justify finiteness of the field,
take the squarefree part and truncate each Puiseux root
past the order of its derivative. The remaining equation
is a simple-root equation after division by a power of
\(s\). The finitely many truncation coefficients are
algebraic, and the recursion for the simple root uses
only their field. Apply this first to the outer prepared
polynomial and then to the inner one after substitution.
Equivalently this is the field assertion of the
[formal root lemma](MULTIVARIATE_OPPOSED_AUXILIARY_POLYNOMIAL_RATIONAL_BIAS.md#a-pole-bound-for-formal-root-branches),
with a number field as its base coefficient field.
Convergent expansions for \(p_1,p_2\) follow by locally
inverting the clocks. Thus \(L_b=(u p_b)^{-1}\) and
\(r=u^{-1}\) are meromorphic in \(s\), with coefficients in
a fixed number field \(K_0\). Branches with \(p_b\) identically
zero cannot contain points with both heights finite.

## Both heights are algebraic over the count

For one height \(L_b\), choose an integer \(C\) at least
the largest pole order of \(L_b,r\) in \(s\). Evaluate an
unknown rational polynomial \(P(L_b,r)\) of degree at
most \(D\), multiply by \(s^{CD}\), and require its first
\((C+1)D\) Taylor coefficients to vanish. Each coefficient
is in \(K_0\), so these are at most
\([K_0:\mathbb Q](C+1)D\) rational linear conditions.
The number \(\binom{D+2}{2}\) of unknowns is larger for
sufficiently large \(D\). Thus a nonzero \(P\) exists, and
its value on the branch is \(O(s^D)\). At ordinary points
\(L_b,r\) are integers. The values of \(P\) therefore lie in
a fixed rational lattice and are zero for all sufficiently
small \(s\) on the ordinary sequence. Analyticity after
multiplication by \(s^{CD}\) forces \(P(L_b,r)=0\) on
the branch. Since \(r=s^{-M}\) is nonconstant, this makes
\(L_b\) algebraic over \(\mathbb C(r)\). Doing this
separately for both heights proves
\[
             \operatorname{trdeg}_{\mathbb C}
                         \mathbb C(L_1,L_2,r)=1.              \tag{2}
\]

## The leading Gaussian direction in the factorial torus

On the branch write the coefficient multipliers in (1) in
the form \(G_{j,k}\exp(r H_{j,k})\), where
\[
 H_{j,k}=\sum_i c_i R_{i,j,k}(1/L_{b(i)}),\qquad
 G_{j,k}=a_{j,k}j^{hk}
       \exp\!\left(\sum_i s_i R_{i,j,k}(1/L_{b(i)})\right).
\]
The functions \(F_{j,k}=\exp H_{j,k}\) and \(G_{j,k}\)
are rational in the heights, since the row counts \(c_i,s_i\)
are integers. The constant coefficients are fixed at one.
Each \(F_{j,k}\) tends to one as \(s\to0\).

Use the fixed lattice of all integral multiplicative
relations among the \(F_{j,k}\) as functions on the entire
branch. The constants in these relations are one, by the
limit at zero. This defines a subtorus \(\mathcal T\) in
the coefficient multiplier torus: the generated character
group is torsion-free modulo constants and a multiplicatively
independent basis gives its monomial parametrization.
The logarithms are chosen near one, so a multiplicative
relation of the \(F_{j,k}\) is an additive relation of the
\(H_{j,k}\) with constant zero. If the rank is zero, then
all the \(H_{j,k}\) vanish identically.

Along the proposed branch, the common-factor minors of degree
at least two vanish by the identity theorem. They are
holomorphic after normalization, and the integral samples
have a nonreal conjugate pair by strict even-order positivity.
Regard these minors as Laurent polynomials on \(\mathcal T\),
with coefficients in \(\mathbb C(L_1,L_2,r)\) after the
\(G_{j,k}\) scales. By (2), the codimension form of the
[Ax rank inequality](MULTIVARIATE_OPPOSED_AX_BRANCH_RANK.md#a-bound-from-differential-transcendence)
gives codimension at most zero for the algebraic target of
the basis exponentials. If the rank is positive, every such
minor must therefore be the zero Laurent polynomial on
\(\mathcal T\). In rank zero the same conclusion holds because
the ordinary branch is already the only multiplier point
of \(\mathcal T\). Thus, for each small nonzero \(s\), the
entire translate by \(G(s)\) has common factor degree at
least two. This assertion uses the one fixed lattice on the
branch.

Let \(a\) be the smaller of the positive orders of
\(1/L_1,1/L_2\) in \(s\). Then
\(\omega_b=\lim_{s\to0}s^{-a}/L_b\) exists, is nonnegative,
and is positive for at least one group. On the real branch
this follows from positivity of the heights at ordinary
points. The leading logarithmic vector is
\[
                 \mathcal V_{j,k}=\lim_{s\to0}s^{-a}H_{j,k}
                      =\sum_i c_i\omega_{b(i)}R'_{i,j,k}(0).
\]
Every fixed lattice relation on \(H\) also annihilates
\(\mathcal V\), by taking this limit. Hence
\((\exp(x\mathcal V_{j,k}))_{j,k}\in\mathcal T\) for every
complex \(x\). In particular no new lattice is chosen when
\(s\) varies. Its Gaussian curvature is nonzero:
\(\sum_b A_b\omega_b>0\). As \(s\to0\),
\(G_{j,k}(s)\to a_{j,k}j^{hk}\). Evaluating along this
direction and letting \(s\to0\) now forces every Gaussian
pair with
\[
                  T=x\sum_b A_b\omega_b,\qquad
                  \lambda=x\sum_b B_b\omega_b
\]
to have a common factor of degree at least two. Its angular
quotient resultant, a nonzero polynomial in the nonconstant
clock \(w=\exp[-T/(2nN)]\), would then vanish identically.
This contradiction proves the theorem.

The argument uses unbounded counts. Bounded counts give fixed-count height families outside this
zero-strength exclusion. For fixed core and template,
the absolute terminal cutoff restricts the argument to
finitely many adjacent pairs.
