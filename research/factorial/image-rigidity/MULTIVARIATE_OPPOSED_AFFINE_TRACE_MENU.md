# A uniform polynomial bound for opposed affine angular exceptions

Fix any ordinary two-monomial core, with arbitrary private and
shared rows, and append two unequal shared rows with opposite
orientations and a common smaller exponent \(A\):
\[
                  (A,A+d_+),\qquad (A+d_-,A),\qquad d_\pm\ge1.
\]
The following bound is independent of the retained core. It also
applies to the more general affine moves below.

**Theorem.** Let \(J\ge4\) be a computable
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Set \(d_0=\min(d_+,d_-)\), \(d_1=\max(d_+,d_-)\), and
\[
 q_j=d_1j-d_0\lceil j/2\rceil,\qquad
 B_J=\sum_{n=2}^{J-2}(n-1)n\bigl((n+1)q_n+nq_{n+1}\bigr).
 \tag{1}
\]
At most \(B_J\) positive integer bases \(A\) fail to give disjoint
consecutive root arguments in every degree. Thus any \(B_J+1\)
distinct positive integer bases give a list of two-row extensions
containing an all-degree angular separation witness for every
ordinary core. The choice of which member works may depend
on the core. Both gaps can be one. With equal gaps the extensions
preserve the original signed degree difference and every exact
endpoint ratio. They can be taken at arbitrarily small added
curvature by choosing a large starting base.

More generally, take a fixed ordinary core with moving affine
rows satisfying the
[opposed-slot condition](MULTIVARIATE_MATCHED_TANH_OPPOSED_AFFINE_ANGLES.md)
at every \(3\le N<J\). If the stripped degree-\(j\) moment
array, after any common factor is removed, has parameter
degree at most \(q_j\) in each coefficient, with nonzero
endpoint coefficients on the parameter ray and throughout
the real interval used for that terminal's endpoint proof,
then the bound (1) applies on the ordinary part of the ray.
One can always take \(q_j=jG\) for \(G\) the sum of all
moving gaps. This gives a bound
\(2G\sum_{n=2}^{J-2}n^2(n^2-1)\). The cutoff \(J\) is
unevaluated.

## The trace polynomial and its degree

For one adjacent pair write \(N=n+1\) and take polynomials
\[
                 F_j(A,z)=\sum_{k=0}^j c_{j,k}(A)z^k,
                        \qquad j=n,N,
\]
with all coefficients of parameter degree at most \(q_j\).
The scalar normalizations can use rational coefficients.
On the real interval in question the two endpoints remain
nonzero. Write \(z_1,\ldots,z_j\) for the roots with
multiplicity, and \(m_j=j(j-1)/2\). Define
\[
 B_j(A,u)=c_{j,j}^{\,2(j-1)}
             \prod_{r<s}\left(u z_rz_s-z_r^2-z_s^2\right).
 \tag{2}
\]
Its ordered input roots can be chosen arbitrarily.
The output roots in \(u\) are the unordered quotient traces
\(z_r/z_s+z_s/z_r\). Its leading coefficient is
\((c_{j,0}c_{j,j})^{j-1}\), by the root product in \(z\).

Each coefficient of the product in (2) is symmetric in the
\(z_r\) and has degree at most \(2(j-1)\) in any one \(z_r\).
A symmetric polynomial with degree at most \(D\) in each
variable is a polynomial in the elementary symmetric functions
of ordinary total degree at most \(D\). To see this, use the
lexicographically leading sorted exponent tuple
\(\alpha_1\ge\cdots\ge\alpha_j\ge0\). The product
\(e_1^{\alpha_1-\alpha_2}\cdots e_j^{\alpha_j}\)
has that leading monomial, and its total degree in the \(e_i\)
is \(\alpha_1\). Its other terms have per-variable degree at
most \(\alpha_1\le D\), so subtracting the leading term
allows finite lexicographic elimination under the same cap.
Substituting the signed coefficient ratios for the \(e_i\)
and multiplying by \(c_{j,j}^{2(j-1)}\) therefore makes
\(B_j\) polynomial in \(A,u\), with
\(\deg_A B_j\le2(j-1)q_j\).

## Counting a new coincidence in the shifted resultant

Make a shifted resultant with a new variable \(t\):
\[
                    R(A,t)=\operatorname{Res}_u
                       (B_n(A,u),B_N(A,u+t)).
\]
Its leading coefficient in \(t\) is nonzero on the real
interval, and its degree in \(A\) is at most
\[
           2(n-1)q_n m_N+2(N-1)q_N m_n
            =(n-1)(N-1)(Nq_n+nq_N).                           \tag{3}
\]
This is the Sylvester multidegree bound. Let \(C(A)\) be the
coefficient of the lowest \(t\)-power that is a nonzero
polynomial, and denote that power by \(\nu\).
The product formula for \(R(A,t)\), including nonzero endpoint
factors, is the product of \(t+\xi_i-\zeta_l\) over every pair of
trace roots \(\xi_i\) of \(B_n\) and \(\zeta_l\) of \(B_N\).
Thus at a real \(A\) with \(C(A)\ne0\), exactly \(\nu\) trace
pairs match, counted with multiplicity. The use of this first
nonzero coefficient allows persistent trace matches anywhere
in the full unordered arrays.

Call a trace in one degree opposite if it comes from a root in
the open upper half-plane and one in the open lower half-plane.
The endpoint argument in the
[opposed affine proof](MULTIVARIATE_MATCHED_TANH_OPPOSED_AFFINE_ANGLES.md#the-exact-adjacent-endpoint)
gives disjoint limiting opposite trace sets. In fact after
positive scaling the upper endpoint is a binomial and its
quotients are \(N\)-th roots of unity. The angle boxes (1) in
that proof put the argument of each oriented opposite lower
quotient strictly between successive multiples of \(2\pi/N\).
Equality of traces identifies their underlying quotients or
their reciprocals, and the \(N\)-th roots of unity are closed
under inversion. On the real interval to the right of that
endpoint the product of pairwise differences of opposite
traces is real-analytic. One obtains this by taking analytic
symmetric functions of the roots inside each open half-plane,
as in the
[half-plane quotient argument](MULTIVARIATE_MATCHED_TANH_OPPOSED_AFFINE_ANGLES.md#carrying-the-half-plane-restriction-through-the-ray).
Even moments have no real roots and odd moments have just one
simple real root there, so the half-plane multisets retain
their sizes. The analytic product has nonzero limit at the
endpoint. Its zeros have empty interior. This reasoning concerns
only opposite traces, allowing identities between other traces.

Now suppose the two moments have a common nonreal argument
at \(A_0\) on the ordinary ray. The two conjugate pairs give an
opposite trace match at \(A_0\). Take real parameters approaching
\(A_0\) with \(C(A)\ne0\) and with no opposite trace matches.
Such parameters exist by the two nonvanishing statements.
Pass to a subsequence and label the roots in each order so
they converge to the roots at \(A_0\) with multiplicity.
This follows from coefficient continuity, nonzero endpoints,
and local disks separating the distinct limiting roots.
The indices of the chosen conjugate pairs are opposite in
each order for all sufficiently close parameters. Their trace
match is therefore absent along the sequence. After another
subsequence, the label pairs constituting the \(\nu\) matches
along the sequence are fixed. All those matches survive at
\(A_0\), and the chosen opposite trace pair supplies an
additional match. The product formula gives an order of
vanishing in \(t\) of at least \(\nu+1\) there, so \(C(A_0)=0\).
This proves that the number of angular exceptions at the
terminal is bounded by (3). This argument keeps all multiplicities
and needs no simple-root hypothesis at \(A_0\).

For ordinary parameters the even moment is strictly positive
on the real axis, excluding common real arguments. This also
settles the first pair of orders directly. The absolute
cutoff handles \(N\ge J\). Summing (3) proves (1).

## Two ordinary opposed rows

Strip \((jA)!\) from each of the two appended rows. At index
\(k\) their coefficient multiplier is
\[
                  (jA+1)_{d_+k}(jA+1)_{d_-(j-k)} .
\]
Here \((x)_\ell\) is the rising factorial. The lengths
\(d_+k,d_-(j-k)\) have maximum at least
\(d_0\max(k,j-k)\ge d_0\lceil j/2\rceil\).
Consequently one common factor
\((jA+1)_{d_0\lceil j/2\rceil}\) can be removed from every
coefficient. It is nonzero throughout \(A>-1/N\) at both
orders \(j=n,N\). The original degree in \(A\) at that index
was \(d_+k+d_-(j-k)\le d_1 j\). This proves the sharper
degree used in (1). The two marked first slots are at the
same affine endpoint and have opposite orientations, so
the required endpoint condition holds at every terminal.

For \(d_+=d_-=d\) the degree is \(q_j=d\lfloor j/2\rfloor\).
The left and right endpoint factors in the moment agree
before normalization, so appending the rows preserves the
endpoint coefficient ratio in every order. The added
curvature is \(2d^2/(A+d)\). Thus the exact finite menu gives
arbitrarily weak ordinary augmentations with these invariants.
Disjoint arguments imply physical coprimality over every
characteristic-zero field by a rational Bézout identity.
The result applies to the augmented supports.
