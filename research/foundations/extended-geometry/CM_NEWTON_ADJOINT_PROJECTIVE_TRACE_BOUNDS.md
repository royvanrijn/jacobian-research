# Projective trace bounds on rational adjoint factors

For a polarized abelian variety \(A/K\), write \(G\) for its
Betti Mumford–Tate group. Assume all its geometric simple
adjoint factors have type \(A_1\). Write its rational adjoint
decomposition as
\[
 G^{\mathrm{ad}}=\prod_i \operatorname{Res}_{R_i/\mathbb Q} H_i,
\]
with each \(H_i\) an absolutely simple adjoint \(R_i\)-form of
\(\mathrm{PGL}_2\). Call a label \(\tau:R_i\hookrightarrow\overline{\mathbb Q}\)
active if the projection of the Hodge cocharacter on that
simple factor is nontrivial, and let \(b_i\) count its active labels.
The cocharacter is minuscule on the adjoint groups. Thus its
nontrivial projections have root difference one, with
adjoint Hodge weights \((-1,0,1)\). Each \(b_i\) is positive:
otherwise the Hodge homomorphism lands in the rational kernel
of that adjoint quotient, contrary to Mumford–Tate minimality.

**Theorem.** After a finite extension of \(K\) defining the
tensors and a fixed level, the following holds at a set of
primes of natural density one. For each \(i\), write
\(\mathcal F_i(v)\) for the primes \(\mathfrak p\) of \(R_i\) over
the residue characteristic at which the adjoint Hodge labels
are active and the adjoint Newton slopes are flat. Then
\[
 \sum_{\mathfrak p\in\mathcal F_i(v)}
       f_{\mathfrak p}\bigl(1+\mathbf 1_{f_{\mathfrak p}\ \mathrm{odd}}\bigr)
       \le b_i-1,
 \qquad f_{\mathfrak p}=[(R_i)_{\mathfrak p}:\mathbb Q_p].
\]
At a prime split completely in \(R_i\), the number \(r_i(v)\)
of flat active labels therefore satisfies
\[
                     r_i(v)\le \left\lfloor\frac{b_i-1}{2}\right\rfloor
                     \qquad\text{for all }i.
\]
In particular \(A\) has density-one ordinary reduction over
this extension if every \(b_i\le2\). The extension need not split
the \(R_i\). There is no restriction on the
endomorphism algebra or on the dimension of the representation
of \(G\) on \(H^1(A)\).

The invariant is the relative adjoint trace plus one. Its
pole at an ordinary active label is simple. At a flat active
label of odd residue degree it has a zero from determinant
parity of the crystalline projective matrix. The resulting
norm argument works on each rational adjoint quotient, even
when \(H^1(A)\) only sees that quotient in a higher tensor.

## Simultaneous local and global Frobenius classes

Enlarge \(K\) so the absolute Hodge tensors defining \(G\) are
defined and the abelian variety with
its tensors and a fixed level away from finitely many primes
gives a point on its Hodge-type Shimura variety over \(K\).
The degree-one primes of this field have
density one, with \(K_v=\mathbb Q_p\) outside ramification.
We use good primes of this kind, with \(p>2\) and outside
the finite exceptional set for the chosen model and unramified
in the \(R_i\). At such a prime the crystalline adjoint Lie
algebra and its filtered de Rham realization decompose over
\(R_i\otimes\mathbb Q_p\). Frobenius commutes with this action
on the linear degree-one crystalline fiber. For each local
field \(L_{\mathfrak p}=(R_i)_{\mathfrak p}\), the filtration is
\(L_{\mathfrak p}\)-linear. Its positive part is either zero or
an \(L_{\mathfrak p}\)-line in the rank-three factor. Thus all
geometric labels over \(\mathfrak p\) are active together or
inactive together.

At these primes the hyperspecial model applies. Indeed one
can choose the level to be hyperspecial at all but finitely
many primes. At a good prime the point extends to the Siegel
model using the abelian scheme and its prime-to-\(p\) level.
The normalization of its closure that defines the Hodge-type
model is finite there, so the point extends to that model
by the valuative criterion. The tensors on its crystalline
cohomology identify it with a tensor fiber for \(G\) over the
completed maximal unramified extension of \(\mathbb Q_p\).
Write its Frobenius as \(b\sigma\) in such a trivialization.
The class \([b]\) has Newton cocharacter bounded by the
dominant Hodge cocharacter and Kottwitz invariant the image
of the Hodge cocharacter. We use the Hodge sign giving Newton
slopes in \([0,1]\) on \(H^1\). Inverting both conventions
has no effect on the adjoint root difference or determinant parity.
[Kisin–Madapusi Pera–Shin, (1.3.4.2), §1.3.5, Proposition 1.3.6, §1.3.7 and Lemma 1.3.8, with (1.1.5.1)–(1.1.5.2) in the linked version](https://people.math.harvard.edu/~kisin/dvifiles/newton2.pdf).

For the same good prime, use the rational semisimple
\(G\)-conjugacy class of Frobenius supplied by Kisin–Zhou,
including the crystalline realization. At degree one this
compares the linear Frobenius coming from the actual
\(\mathbb Q_p\)-crystalline tensor fiber with the auxiliary
étale Frobenius, up to tensor-preserving trivialization over
algebraic closures. For each quotient, the invariant
\[
 T_i=1+\operatorname{tr}(\operatorname{Ad}_{H_i})
\]
is a regular function over \(R_i\). Evaluated at the relevant
projection of the class it gives \(T_i(v)\in R_i\). To see the
descent, the conjugacy class is Galois-stable over \(\mathbb Q\),
and its evaluations at the labels of the rational restriction
of scalars form a Galois-equivariant tuple. This uses the
compatibility of the \(G\)-class and does not require a lift
of that class to a rational point of \(G\).
[Kisin–Zhou, Theorem 1.2 and Theorem 5.3.10 in the preprint](https://arxiv.org/abs/2505.02165).

At every complex embedding, if \(u\) and \(u^{-1}\) are the
nontrivial adjoint eigenvalues, then \(T_i=2+u+u^{-1}\).
Since the derived Lie algebra is an invariant subspace of
\(\operatorname{End}(H^1)\), these eigenvalues are ratios of
Weil numbers for the reduction. The conjugates have complex
absolute value one and are units at all finite places away
from \(p\). Consequently \(0\le\tau(T_i(v))\le4\) at every
complex embedding \(\tau\), and \(T_i(v)\) is integral away
from \(p\). The invariant is unchanged on inverting Frobenius.

## The local trace at a flat active factor

At an inactive \(L_{\mathfrak p}\)-factor the adjoint Newton
slopes are all zero by the group Newton inequality: the
averaged Hodge cocharacter on this orbit is trivial. Hence
the local valuation of \(T_i\) there is nonnegative.

At an active \(L_{\mathfrak p}\)-factor, the corresponding rank-three Lie algebra
on crystalline cohomology is a split \(L_{\mathfrak p}\)-form of
\(\mathfrak{sl}_2\). In fact the one-dimensional positive
Hodge filtration is defined over \(L_{\mathfrak p}\) and becomes
a nilpotent root line over the algebraic closure. The
adjoint Hodge filtration has this description under the
filtered comparison for the tensors. A quaternionic inner
form of \(\mathfrak{sl}_2\) has a nonzero rational nilpotent
only when the quaternion algebra is split. Thus we may
represent linear adjoint Frobenius here by a projective
matrix with a lift \(B\in\mathrm{GL}_2(L_{\mathfrak p})\).
This local description does not require a rational
trivialization of the *whole* crystalline tensor fiber.

Over the completed maximal unramified field this choice of
isomorphism of the adjoint factor changes the group Frobenius
by \(\sigma\)-conjugacy. At a split \(\mathrm{PGL}_2\) factor
the Kottwitz map is the valuation of the determinant of a
lift modulo two, independent of the scalar in the lift. By
functoriality of the Kottwitz invariant, for the unramified
restriction of scalars factor the coinvariant map sums that
determinant parity over the \(f_{\mathfrak p}\) labels. The
matrix \(B\) has equal determinant valuation at all embeddings
of \(L_{\mathfrak p}\), and the Hodge projection is minuscule
at every label of this factor. Consequently
\[
             f_{\mathfrak p}v_p(\det B)\equiv f_{\mathfrak p}\pmod2.
\]
The adjoint Newton gap is the difference of the two
eigenvalue valuations of \(B\), taken nonnegative. It is
at most one by the group Newton inequality. If it is
positive, the two valuations are integers by the Newton
polygon of a quadratic over the unramified
\(L_{\mathfrak p}\), so the gap is one. In that case
\(v_p((\operatorname{tr}B)^2/\det B)=-1\). If the gap is zero,
the two valuations equal \(v_p(\det B)/2\). The trace lies
in \(L_{\mathfrak p}\), so its valuation is at least the
ceiling of this value. When \(f_{\mathfrak p}\) is odd the
determinant valuation is odd and this bound is strict. Thus
\[
                   v_p(T_i(v))=
                   v_p\bigl((\operatorname{tr}B)^2/\det B\bigr)
                   \ge \mathbf 1_{f_{\mathfrak p}\ \mathrm{odd}} .
\]
This includes the possibility of trace zero. Using the
linear Frobenius over \(\mathbb Q_p\) on the
\(L_{\mathfrak p}\)-module gives the gap before passing to
any power, even if \(\mathfrak p\) has even residue degree.

## Norm and reduction density

For nonzero \(T_i(v)\), the local estimates give
\[
 v_p\bigl(N_{R_i/\mathbb Q}T_i(v)\bigr)\ge -b_i+
   \sum_{\mathfrak p\in\mathcal F_i(v)}
     f_{\mathfrak p}(1+\mathbf 1_{f_{\mathfrak p}\ \mathrm{odd}}).
\]
If the sum reaches \(b_i\), the norm is an integer of absolute
value at most \(4^{[R_i:\mathbb Q]}\). The same finite set
contains a zero norm. This norm is the Frobenius value of
the regular class function \(\prod_\tau T_{i,\tau}\) on \(G\).
It is nonconstant: on a torus of one geometric adjoint
factor, holding the other factors at the identity gives
\(4^{[R_i:\mathbb Q]-1}(2+z+z^{-1})\). Vasiu's Mumford–Tate
theorem in type \(A_1\) identifies the connected auxiliary
monodromy with \(G\), and the base extension above makes the
whole auxiliary image lie in \(G\). Algebraic Chebotarev for
the finite set of norm levels gives the asserted exception
of density zero, simultaneously for the finitely many \(i\).
[Vasiu, Theorem 1.3.4(a) and Example 1.3.5(iv)](https://arxiv.org/pdf/math/0212066v7),
[Rajan, Theorem 3](https://repository.ias.ac.in/38317/1/2-PUB.pdf).

Every failed active prime would contribute at least two to
the sum, so \(b_i\le2\) forces all active gaps to attain one.
In that case the dominant Newton
cocharacter has the same projection to every adjoint
factor as the dominant Hodge cocharacter. Their images in
the rational cocharacter space of the abelianization agree
by the Kottwitz equality and Galois averaging. The Hodge
class is rational over \(\mathbb Q_p\) because the tensors
and its filtration descend there. The adjoint and
abelianization maps together give an injection on rational
cocharacters of a geometric maximal torus, so the dominant
cocharacters coincide.
The Newton slopes on \(H^1\) are then all zero or one.

When a rank-two coefficient block already supplies a
relative Frobenius polynomial, the local trace-square
argument is the familiar one for that block. The result
here applies it to the adjoint quotients separately, including
tensor cases with no rank-two block in \(H^1\).
[Rank-two projective trace comparison](CM_NEWTON_PROJECTIVE_TRACE_FIELD.md).
Pink's density-one theorem over the connected field treats
geometric endomorphism ring \(\mathbb Z\), and its proof
already uses finite ranges of intrinsic Frobenius traces.
[Pink, Theorem 7.1, Corollary 7.2 and §7](https://people.math.ethz.ch/~pink/ftp/AMT-v3.pdf).

For exact real multiplication by a totally real cubic field,
the three independent active labels belong to one rational
simple adjoint quotient. At completely split primes the
bound is one flat label.
Eliminating that single failure still requires a further
argument for the split real-cubic case.

## Recovering the slopes on first cohomology

**Proposition.** At the good absolute degree-one primes of the
tensor-defining extension used for the local calculation, all Newton
slopes on \(H^1\) lie in \(\{0,1/2,1\}\). Let \(2c_i\) be the sum
of the dimensions, with multiplicity, of the irreducible geometric
summands on which any given label of the rational adjoint factor \(i\)
acts nontrivially. The sum is independent of the label. If \(r_i(v)\)
counts flat active labels, the half-slope multiplicity is
\(2\sum_i c_i r_i(v)\) and the \(p\)-rank is \(g-\sum_i c_i r_i(v)\).
Over an extension that splits the label fields, the weighted bound
makes the \(p\)-rank strictly greater than \(g/2\) on a set of natural
density one for positive-dimensional \(A\).

The proof is in [the manuscript](../../../unpublished/cm-newton-polygon-constraints/adjoint.tex).
The two Hodge weights allow at most one active standard factor in an
irreducible summand. The Newton and Hodge central weights agree at
these degree-one places. Equal adjoint slopes at that active factor
therefore replace the pair of Hodge weights by their midpoint, while
the other summands retain their weights. After the label fields split,
each flat active label costs two in the norm bound and the total
dimension of the active summands is at most \(2g\).

Fit\'e's degree-\(g\) embedded Galois endomorphism field criterion gives
positive-density \(p\)-rank at least \(\lceil g/2\rceil\). The adjoint
criterion uses a group-type hypothesis and tracks multiplicities in
tensor representations. Sheng and Zuo obtain the binomial and
supersingular shapes for Mumford's type with one active embedding;
their nonsplit case uses an orbit average that does not reduce to the
degree-one tensor-defining situation above.
[Fit\'e, Theorem 2](https://arxiv.org/html/2005.10185v3#S1)
· [Sheng–Zuo, Theorem 1.2](https://arxiv.org/pdf/1106.3505).
