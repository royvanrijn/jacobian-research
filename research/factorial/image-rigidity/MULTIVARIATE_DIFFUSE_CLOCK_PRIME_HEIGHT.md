# Prime-height restrictions at physical diffuse contacts

Fix any ordinary two-monomial core \(C\) and adjacent orders
\(n,N=n+1\), \(n\ge2\). Append arbitrary ordinary shared rows, with
signed gap \(h\), added curvature \(\tau\), signed bias \(\lambda\)
and largest relative gap \(\eta\) as in the
[physical diffuse comparison](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md).
Let
\[
 q=\exp[-\tau/(2nN)],\qquad v=\exp[-\lambda/(2nN)].
\]
In a sequence of tails take \(\eta\to0\), \(\tau\to\tau_*>0\),
\(\lambda\to\lambda_*\), and fix the integer \(h\). Write \(q_*,v_*\)
for the limiting clocks. At exact collisions the bound on the total
signed gap ensures a fixed \(h\) on a subsequence.

**Theorem.** There is a finite computable list of positive integers
\(\mathscr L_{n,N}\) depending only on the orders. Suppose the
ordinary consecutive moments have a physical common root along
the sequence. If the rational coefficient ratios defined below
lie in the units of one fixed finite prime set \(\mathcal P\), then
some \(d\in\mathscr L_{n,N}\) satisfies
\[
        q_*^d,\ v_*^d\ \in \mathbb Z[\mathcal P^{-1}]^\times .
                                                               \tag{1}
\]
Thus a positive physical reference contact that fails this
condition has an ordinary exclusion neighborhood within that
prime class. The strict prime-support assumption admits a
height relaxation. In that relaxation the conclusion is
\(q_*^d,v_*^d\in\mathbb Q_{>0}\). The subspace theorem gives
existence of the neighborhoods.

The coefficient ratios use exact rational numbers and a single
common argument normalization. For a tail \(E\) put
\[
 R_{j,k}=\prod_{(a,b)\in E}
    \frac{((j-k)a+kb)!}{(ja)!},\qquad \Lambda=N/n,
 \qquad
 \begin{cases}
 X_{n,k}=R_{n,k}/R_{n,1}^{\,k},\\
 X_{N,k}=R_{N,k}/(\Lambda^h R_{n,1})^{\,k}.
 \end{cases}                                                   \tag{2}
\]
The free coordinates of \(X\) are \(X_{n,k}\) for \(2\le k\le n\)
and \(X_{N,k}\) for \(1\le k\le N\). The others equal one.
This allows cancellations of prime factors across rows and
degrees. The arithmetic hypothesis on \(X\) permits arbitrary row counts,
gaps and heights.

## The two eliminants in rational coordinates

Write \(b_{j,k}=[z^k]Q_j^C\) for the ordinary core. Let
\[
 G_j(X,z)=\sum_{k=0}^j b_{j,k}X_{j,k}z^k,\qquad
 \mathcal D(X)=\operatorname{Res}_z(G_n(X,z),G_N(X,\Lambda^h z)).
\]
The original ordinary pair differs from this displayed pair
only by nonzero scalars and the common argument change in (2).
For a degree \(j\) polynomial \(F\), put
\[
 U_j(F;u)=\operatorname{Res}_z(F(z),F(uz))/(u-1)^j,\qquad
 \mathcal A(X)=\operatorname{Res}_u(U_n(G_n;u),U_N(G_N;u)).
\]
Division here is exact in the polynomial coefficient ring. The
off-diagonal root quotients, with multiplicity, are exactly the
roots of \(U_j\) when endpoints are nonzero. Every physical
ordinary common root is nonreal because the even moment is
strictly positive on the real axis. The common root and its
conjugate therefore force \(\mathcal A(X)=0\), as well as
\(\mathcal D(X)=0\). All coordinates of \(X\) are positive.

The simultaneous coefficient and center comparison for diffuse
tails gives \(X\to\Psi(q_*,v_*)\), where the Laurent monomial map is
\[
 \Psi(q,v)_{n,k}=q^{Nk(1-k)},\qquad
 \Psi(q,v)_{N,k}=q^{k(N-nk)}v^k .                              \tag{3}
\]
Indeed the separate tail geometric endpoint scales \(D_j\)
satisfy \(R_{j,k}/D_j^k=q^{(nN/j)k(j-k)}(1+o(1))\), and
\(D_N/D_n=\Lambda^h v(1+o(1))\). Dividing by the actual
first coefficient factor \(R_{n,1}/D_n=q^{N(n-1)}(1+o(1))\)
gives (3).

The simultaneous zero set of
\(\mathcal A(\Psi(q,v)),\mathcal D(\Psi(q,v))\) in
\((\mathbb C^\times)^2\) is finite. To see this, use the raw
Gaussian core polynomials
\[
 H_j(q,z)=\sum_k b_{j,k}q^{(nN/j)k(j-k)}z^k .
\]
Equations (3) give \(G_n(\Psi,z)=H_n(q,q^{-N(n-1)}z)\) and
\(G_N(\Psi,z)=H_N(q,q^{-N(n-1)}vz)\). Separate nonzero
argument scales leave root quotients unchanged. Thus the
angular equation implies the raw one-variable quotient
resultant in \(q\) vanishes. This resultant is a nonzero
polynomial since at \(q=0\) the endpoint binomials have
disjoint nontrivial root-of-unity quotient sets. At each
nonzero \(q\), the raw physical resultant as a function of
\(v\) is nonzero at \(v=0\). In the root product its value
there is \(b_{N,0}^{\,n}b_{n,n}^{\,N}\), by continuation from
\(v\ne0\). This proves finiteness even at complex clocks and
with multiple reference roots.

## Fixed monomial ratios at ordinary collisions

Use a finite-support expansion of each polynomial \(P=\mathcal A\)
or \(\mathcal D\), with like monomials combined and zero
coefficients removed. At any zero \(X\) partition its terms
\(c_aX^a\) into minimal nonempty blocks of sum zero. Every block
contains at least two terms. Along a subsequence fix both
partitions, and choose one anchor \(a_0\) in each block.
The ratios \(X^{a-a_0}\) of the other terms satisfy a fixed
nondegenerate equation with right side one. If the coordinates
of \(X\) are units of \(\mathbb Z[\mathcal P^{-1}]\) the ratio
tuples belong to a group of finite rank. Finiteness follows
from the [unit equation theorem of Evertse, Schlickewei and Schmidt, Theorem 1.1](https://arxiv.org/pdf/math/0409604).
For a two-term block the ratio is already determined by its
two coefficients. There are finitely many blocks. On a further
subsequence all these ratios are fixed positive rational
numbers, and remain units of \(\mathbb Z[\mathcal P^{-1}]\).

For each vector \(u\) of monomial exponents define its
two-dimensional weight \(\omega(u)\) by
\(\Psi(q,v)^u=q^{\omega(u)_1}v^{\omega(u)_2}\).
At the limit, every ratio from the preceding paragraph gives
\[
      q_*^{\omega(a-a_0)_1}v_*^{\omega(a-a_0)_2}=r_{a,a_0},
                                                                  \tag{4}
\]
where \(r_{a,a_0}\) is its fixed rational value. These weights
must span rank two over \(\mathbb Q\). Otherwise a nonzero
integer vector \((e,f)\) is orthogonal to them all, so
\((q_*s^e,v_*s^f)\) has the same ratios in every block for
all \(s\in\mathbb C^\times\). Dividing each block by its
anchor monomial then gives the same zero sum at every point
of this infinite curve, contradicting the simultaneous
finiteness proved above.

Two independent rows of (4), with determinant of absolute
value \(d\), give (1) by integer elimination. To compute a
universal finite list \(\mathscr L_{n,N}\), use the supports
of the two formal resultants in indeterminate coefficients
before specializing to core coefficients. Drop the exponents
of the three coordinates fixed to one in (2), and retain
all nonzero absolute determinants of pairs of weights of
within-support monomial differences. Allow either resultant
for each difference. Rational coefficient specialization and
the integer twist can only remove monomials from these
supports. The list therefore works for every core and twist.

## A height version allowing arbitrarily many outside primes

For a nonzero rational number \(x\), let \(\mathrm h(x)\) be the
logarithm of the maximum of the absolute reduced numerator and
the positive reduced denominator. Write \(\mathrm h_{\mathcal P^c}(x)\)
for the same height of the factor of \(x\) supported on primes
outside \(\mathcal P\). The signs of rationals have zero height.

Let \(T\ge3\) bound the number of monomials in each eliminant,
and define a positive constant
\[
 a_T=\min_{2\le k\le T} k^{-1}\exp[-(4k)^{3k}],\qquad
                     \epsilon_T=a_T/(4T).
\]
It suffices to assume, for some fixed finite \(\mathcal P\) and
fixed real \(B\), that every relevant ratio along the collision
sequence satisfies
\[
 \mathrm h_{\mathcal P^c}(X^{a-b})
              \le\epsilon_T\,\mathrm h(X^{a-b})+B.               \tag{5}
\]
The ratios needed are only differences of two monomials within
the same support. Under (5) the minimal-block ratios can still
be made constant on a subsequence, though the constants may
use primes outside \(\mathcal P\).

Here is the deduction of that constancy. First pass to a
subsequence where each of the finitely many ratios is constant
or has height tending to infinity. Bounded rational height
gives finiteness by counting possible numerators and denominators.
In a minimal zero block group terms whose ratios are constant.
If more than one class remains, collapse each class to its
representative monomial times a rational coefficient. The
coefficient is nonzero by minimality since the class is a
proper subset of the block. Thus one obtains an exact linear
relation among \(t\ge2\) monomials, all of whose pairwise ratios
have height tending to infinity and satisfy (5). We show
inductively that such a relation is impossible for \(t\le T\).

For \(t=2\) the ratio must be constant. For \(t\ge3\), choose
one monomial as anchor and divide by its term, obtaining
\(\sum_{i=1}^{t-1}y_i=1\). Enlarge \(\mathcal P\) by the
coefficient primes, and split the vector coordinatewise as
\(y=xz\), with \(x\) supported on \(\mathcal P\) and \(z\)
on its complement. Let \(H\) be the largest height of
an anchor ratio before multiplying by the fixed coefficient
ratios. Then \(H\to\infty\), and for the logarithmic vector
height of \((1:y_1:\cdots:y_{t-1})\) one has
\[
 \mathrm h(z)\le (t-1)(\epsilon_T H+B),\qquad
 \mathrm h(y)\ge H-O(1),\qquad
 \mathrm h(x)\ge\mathrm h(y)-\mathrm h(z).
\]
The vector-height inequalities follow by taking local
maxima in the definition of height. In particular
\(\mathrm h(x)\ge H/2\) eventually, and the bound on
\(\epsilon_T\) gives \(\mathrm h(z)\le a_T(1+\mathrm h(x))\)
eventually. Apply
[Evertse, Schlickewei and Schmidt, Theorem 2.1](https://arxiv.org/pdf/math/0409604)
in dimension \(t-1\). The \(y\)'s lie in finitely many proper
linear subspaces, since the \(\mathcal P\)-supported vectors
form a finite-rank group. On an infinite subsequence there
is a nonzero linear relation among the \(y_i\) alone with
fixed rational coefficients. For example, take a linear
form vanishing on that subspace and take a nonzero rational
coordinate of its coefficients in a number field basis.
This gives a relation among at most \(t-1\) of the original
monomials, without the anchor. A one-term relation is
impossible since the monomials are nonzero. Otherwise the
pairwise ratios still satisfy the hypotheses, contradicting
induction.

Therefore every original minimal block consists of just one
constant-ratio class. Equations (4) and the rank-two argument
apply and give the claimed rational-power conclusion. In
particular, at a contact failing every rational-power test
from \(\mathscr L_{n,N}\), any ordinary collision sequence
approaching the contact must violate (5) for each finite
\(\mathcal P\), for every fixed \(B\). A monomial ratio from
one of the two actual resultants must carry a positive
proportion of its height outside \(\mathcal P\) on a subsequence.
The theorem leaves open sequences with such escaping height.

## A finite radical-clock filter for each core

The candidate rational powers can also be restricted by the
prime factors of the core. For an algebraic physical reference
contact use the [prime-adic clock bounds](MULTIVARIATE_GAUSSIAN_CLOCK_PADIC_EXCEPTIONS.md)
for \(q_*\), writing \(U_p=\max(P_p,V_p)\) in the notation of
that result. Then \(|v_p(q_*)|\le U_p\) in every embedding.
Put \(a_{j,k}=v_p(b_{j,k})\) and
\(B_j=\max_{0\le k\le j}(nN/j)k(j-k)\). The diameter of the
coefficient valuation ordinates of \(H_j(q_*,z)\) is at most
\[
 T_{j,p}=\max_k a_{j,k}-\min_k a_{j,k}+2B_j U_p.
\]
Each root has valuation of magnitude at most \(T_{j,p}\)
by the Newton polygon rule, since each edge has integer width
at least one. At a physical reference contact the upper root
of the raw \(H_N\) is \(\Lambda^h v_*\) times a lower root.
Consequently
\[
                    |v_p(v_*)|\le T_{n,p}+T_{N,p}
                                             +|h v_p(\Lambda)|.
\]
All of these bounds vanish at primes \(p>\max(N,N\max_i(a_i,b_i))\),
where the maximum is over the core rows and can be taken to
be zero for the empty core. Thus for each proposed
\(d\in\mathscr L_{n,N}\), the positive rationals \(q_*^d,v_*^d\)
have an explicit finite list of possible prime factorizations.
Strict prime support further restricts the list to primes in
\(\mathcal P\).

The determinant screen can retain the actual two weight rows
in (4) as well. In particular the joint field degree
\([\mathbb Q(q_*,v_*):\mathbb Q]\) is at most their absolute
determinant \(d\). The two rational monomial equations in
(4) have \(d\) complex solutions on the torus, as the
integer Smith normal form reduces them to two independent
power equations whose exponents have product \(d\).
Every embedding of the joint field gives one of these
solutions.

There is a direct exclusion when every coordinate \(X_{j,k}\)
is a unit at all non-good primes of the core, including the
primes of \(nN\). This holds for either the exact finite prime
support condition or the height version (5) with the same
unit requirement. The constant ratios in (4) are then units
at these primes, so both rational powers in (1) are units
there. The core prime-adic bounds already make them units
at every other prime. Positivity would give \(q_*^d=1\),
contrary to \(\tau_*>0\). Thus every compact positive
strength interval admits a uniform exclusion in this class.

For a fixed core and \(\tau\) in a prescribed compact positive
interval, this gives a finite algebraic screen for all
candidate limiting points. Form the unique positive \(d\)-th
roots of the candidate tuples and test the simultaneous
resultants, including the physical cone and signed gap bound.
The [absolute cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
leaves finitely many order pairs and the absolute signed gap
bound leaves finitely many \(h\). When the screen is empty,
compactness gives a single positive diffuseness tolerance for
the whole interval and all tail row counts obeying the stated
prime or height hypothesis. The height bounds and the clock
screen are explicit, while the tolerance derived from the
subspace theorem has no computed value here.
