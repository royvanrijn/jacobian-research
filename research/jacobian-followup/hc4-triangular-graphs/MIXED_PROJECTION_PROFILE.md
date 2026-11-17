# Exact generator counts for every mixed gradient projection

14 November 2026. The stable hypersurface obstruction has a sharp
boundary under mixing of gradient outputs. The boundary can be
classified in every degree, rather than tested one projection at a time.

For a hypothetical nonautomorphic HC4 gradient, let U be its intrinsic
two-dimensional cotangent output space. If V is any k-dimensional
space of linear combinations of the four gradient components, put

\[
                            \ell=\dim(V\cap U).
\]

The minimum number of additional polynomial algebra generators over
C[V] is exactly

\[
                              4-k+\ell.                      \tag{1}
\]

Equivalently the smallest ambient affine space admitting a closed
embedding that keeps these k outputs as coordinate restrictions has
dimension **4+ell**. The same formula persists after adjoining q free
source variables, with 4 replaced by 4+q. This is a relative statement
over the chosen output algebra, not an absolute embedding dimension.

In particular, for two independent gradient combinations there are
three exhaustive cases:

| Dimension of their intersection with U | Smallest ambient dimension | An A5 completion |
|---|---:|---|
| 0 | 4 | Exists by an explicit source coordinate change |
| 1 | 5 | Exists by an explicit graph embedding |
| 2 | 6 | Impossible by HC4MYHJ1 |

Thus mixed projections outside the last row cannot be excluded by
coordinate-algebra generation alone. Their actual compatibility with
the fixed Meng--Yang potential remains an additional obligation.

## 1. A Keller component cannot become a coordinate after stabilization (HC4MYHK0)

Let K=(K1,K2) be a complex plane Keller map. For a nonzero constant
row a, put P=a1*K1+a2*K2. If P becomes a coordinate after adjoining
any finite number r of independent variables, then K is an
automorphism. Conversely, an automorphic K makes every such P a
plane coordinate already.

Here is the needed special case without importing a general stable
coordinate classification. If P is a coordinate of C[z1,z2,eta1,...,eta_r],
its zero fibre has an isomorphism

\[
                    C_P\times A^r\simeq A^{r+1},\qquad
                    C_P=\{P=0\}\subset A^2.                  \tag{2}
\]

Every fibre of P is smooth because a nonzero row of the invertible
Keller Jacobian never vanishes. Write c for the number of connected
components of C_P, G for the sum of their smooth completion genera,
and b for the total number of omitted completion points. Each affine
component omits at least one point. The retained compactly supported
Hodge--Deligne formulas give, with t=uv,

\[
       E(C_P)=ct-G(u+v)+c-b,\qquad E(C_P)t^r=t^{r+1}.         \tag{3}
\]

Equality of these polynomials forces c=1, G=0 and b=1. Hence C_P
is an actual affine line: a smooth genus-zero complex completion is
P1, with exactly one point removed. Complete the constant row a to
an invertible linear target change of K. Then C_P is an affine-line
component of a plane Keller fibre. Retained
[HC4MYGR0](CRITICAL_COMPONENTS.md#1-a-plane-keller-fibre-with-an-affine-line-component)
makes K invertible. This proves the implication for every r.

The underlying plane stable-coordinate fact is classical, and is not
asserted new. Related primary context is
[Makar-Limanov, van Rossum, Shpilrain and Yu, Theorem 1.1](https://ems.press/content/serial-article-files/43009),
which proves plane stable equivalence implies plane equivalence.
The [context record](sources/projection-profile-context.json) identifies
that theorem as context, rather than a new unproved classification
input to (3). No higher-dimensional stable-coordinate conjecture or
general affine-space cancellation theorem is assumed.

## 2. The exact stable projection profile (HC4MYHK1)

We prove a slightly more general algebra statement. Take a
nonautomorphic plane Keller map K, r>=0 independent source variables
eta, N=2+r, and

\[
 B=C[z_1,z_2,\eta_1,...,\eta_r],\qquad
 G=\operatorname{span}_C(K_1,K_2,\eta_1,...,\eta_r),\qquad
 U=\operatorname{span}_C(K_1,K_2).
\]

These N polynomials are algebraically independent, so G is an
N-dimensional vector space of actual polynomials. Let V be any
k-dimensional linear subspace of G, let S be a basis of V, and set

\[
 \ell=\dim(V\cap U),\qquad d=k-\ell.
\]

Then the minimum number of additional generators of B as a C[S]-algebra
is exactly

\[
                       \mu(B/C[S])=N-k+\ell=N-d.             \tag{4}
\]

Constants added to the S_i do not change their algebra or this count.
The assertion is conditional on K being nonautomorphic; for an
automorphic K every V is a coordinate subspace and the answer is
N-k instead.

### Upper bound with polynomial source recovery

Write S=A*K+B_eta*eta as a full-row-rank constant matrix. The rank
of B_eta is d. By invertible constant changes of the S rows and eta
coordinates, put S into the form

\[
 S_i=\eta_i+a_iK\ (1\le i\le d),\qquad
 S_{d+\alpha}=p_\alpha K\ (1\le\alpha\le\ell),               \tag{5}
\]

where the p_alpha are independent constant rows on the two K outputs.
The following N-d functions generate all of B over C[S]:

\[
                         z_1,z_2,\eta_{d+1},...,\eta_r.      \tag{6}
\]

Knowing z gives K(z), after which eta_i=S_i-a_iK for i<=d.
This is an actual polynomial recovery identity, with no inverse to
K assumed. It gives a closed embedding into A^(N+ell) whose first
k coordinate functions are S.

### Lower bound when ell=0

Algebraic independence makes the relative transcendence degree N-k.
At least N-k additional generators are necessary. Here d=k, so
(6) attains that number and supplies a genuine polynomial coordinate
system on A^N.

### Lower bound when ell=1

Suppose only N-k additional generators suffice. Together with S
they would be N generators of a polynomial algebra of transcendence
degree N, hence algebraically independent and an actual coordinate
system. The nonzero pure Keller combination p_1*K in (5) would then
be a coordinate of B. This contradicts HC4MYHK0. Thus at least
N-k+1 generators are necessary, exactly the number in (6).

### Lower bound when ell=2

Both K outputs can be recovered by an invertible constant change
among the S_i. If N-k+1 or fewer further generators sufficed, padding
by zeros would give a closed embedding

\[
                         A^N\longrightarrow A^{N+1}
\]

over the first two outputs K(z). The arbitrary-r stable hypersurface
obstruction [HC4MYHJ1](STABLE_HYPERSURFACE_OBSTRUCTION.md)
forbids this for nonautomorphic K. Thus at least N-k+2 helpers are
necessary. Again (6) attains that number. This completes the proof
of (4), for every r,k,V and every polynomial degree.

Adjoining q further independent source variables leaves V and ell
unchanged and replaces r by r+q. Hence the helper minimum increases
by exactly q, while the excess ell over the relative transcendence
degree remains unchanged. No finite rank/degree sample proves this
stable quantifier; it is part of the written argument.

## 3. Intrinsic HC4 consequences and the limit of the embedding route

Under nonautomorphy, Tang gives f=D(z)+w^T K(z). The actual source
coordinate change

\[
                     \eta=\nabla D+DK^T w
\]

has polynomial inverse w=DK^(-T)(eta-grad D), and the gradient
output space is span(eta1,eta2,K1,K2). Thus (4) applies with r=2,
independently of the degrees of D,K or a proposed graph embedding.
We use this coordinate change for the source algebra only, and do
not assert Hessian invariance under it.

An ambient coordinate completion in A5 of a prescribed set V of
gradient combinations exists **if and only if V does not contain U**.
Necessity follows from ell=2. For ell<=1, (6) constructs it, padding
with zero when fewer than five coordinates are required. These
positive constructions give even a coordinate hypersurface: when
ell=1 the single relation is S_(d+1)-p_1*K(z)=0, linear in its own
ambient S coordinate. For ell=0 one simply pads an actual source
coordinate system. They are not constructions with the fixed
Meng--Yang potential, and do not transfer a marked gradient collision.

For two outputs, write S=A*K+B_eta*eta. Then

\[
                 \ell=2-\operatorname{rank}(B_\eta).
\]

The three ranks 2,1,0 give the three rows of the opening table.
A generic mixed two-output projection has rank two and is already
a source coordinate pair. Thus a universal prohibition of mixed
physical bases cannot follow from closed-embedding geometry alone.
The potential identity f=Psi composed with i, actual Hessian
constancy in the chosen source coordinates and any marked covector
transfer must be imposed separately.

There is also an algebraic recognition of the unique cotangent plane.
For a nonzero constant source direction v,

\[
 D_v f\text{ is a polynomial coordinate}
       \quad\Longleftrightarrow\quad v\notin W,              \tag{7}
\]

where W is the intrinsic jointly affine two-plane of HC4MYHF0/1.
If v is outside W, the derivative has a nonzero eta coefficient,
and a shear as in (5) completes it to coordinates. Inside W it is
a nonzero linear combination of K, and HC4MYHK0 says it cannot
become a coordinate even after any polynomial stabilization.
For an automorphic gradient every nonzero directional derivative
is a coordinate. Thus under unit Hessian the nonzero noncoordinate
directional derivatives either form exactly W minus zero, or there
are none. The noncoordinate directions cannot span three dimensions.

The [exact checker](verify_projection_profile.py) and
[receipt](receipt-projection-profile.json) verify rank-stratum controls,
source-recovery identities, the table and curve Hodge coefficient
separation, while byte-pinning the retained theorem inputs. They do
not construct a nonautomorphic Keller map or certify the lower bounds
by finite examples.

The previous turn's fixed-output-base obstruction is therefore sharp
within all linear gradient projections. A remaining general graph
argument must use the fixed ambient potential or another restriction
beyond the existence of an embedding for an allowed projection.
The universal SL2 row equation, unit-ideal plane surjectivity, distinct
singleton-zero properness and unrestricted HC4/JC2 remain open.
No formal verification, independent replay, external review or
novelty priority is claimed.
