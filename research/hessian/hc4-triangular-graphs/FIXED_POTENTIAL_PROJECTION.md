# The fixed potential closes the rank-one mixed projection

14 November 2026. The exact mixed-projection generator profile leaves
two kinds of A5 completions: eta rank two and eta rank one. The fixed
Meng--Yang potential rules out the latter for a nonautomorphic HC4
gradient. This uses only its common-zero-base boundary, rather than
an increase in polynomial degree or a classification of embeddings.

Write the scaled potential as

\[
 \Psi=L A^2+M A+N B,\qquad LN\ne0,
\]

with the retained native linear forms A,B. On physical x=y=0,

\[
                              \Psi=Lp^2+Mp.                  \tag{1}
\]

All new statements below are over C; the elementary cotangent
calculation itself works over characteristic-zero fields.

## 1. A composite boundary closes the mixed rank-one case (HC4MYHL0)

Let

\[
 f=D(s,U)+w_1P(s,U)+w_2Q(s,U),\qquad
 J(P,Q)=c\in C^\times,
\]

and set eta=grad D+DK^T w, K=(P,Q). Suppose a physical base pair,
after invertible constant row changes preserving its common zero
scheme, has pullbacks

\[
                b_1=P(s,U),\qquad b_2=\eta_1+v(s,U),         \tag{2}
\]

for any polynomial v. If the restriction of f to b1=b2=0 is a
univariate polynomial H(h) of degree at least two in one polynomial
h on that scheme, then **K is an automorphism**.

More generally the eta coefficient in b2 may be any nonzero constant
row: choose constant quotient coordinates (s,U) making it eta1.
The pure Keller component b1 must have zero offset in the chosen
Keller pair; no arbitrary translation P-alpha is included.

Here is the all-degree argument. The polynomial vector

\[
 b=(DK)^{-1}K,\qquad
 b_2^{\rm vec}=\frac{-Q_sP+P_sQ}{c},                         \tag{3}
\]

gives the exact source-algebra identity

\[
                  f=D-b\cdot\nabla D+b\cdot\eta.             \tag{4}
\]

The notation b2_vec in (3) refers to the second vector component,
not the physical base function b2 in (2). On the scheme (2),
eta1=-v and eta2 is a free polynomial variable. Thus f is affine
in eta2, with coefficient b2_vec modulo P.

The polynomial P is squarefree. A repeated irreducible factor would
divide both P derivatives and hence the unit Jacobian c. For the
same reason gcd(P,Q)=1. Work on each irreducible component of P=0.
Its coordinate ring is a domain. If h has positive eta2 degree m,
then H(h) has degree (deg H)*m>=2 in eta2, contradicting affineness.
Therefore h is independent of eta2 on every component and

\[
             P\mid b_2^{\rm vec},\qquad P\mid P_s Q,
             \qquad P\mid P_s.                              \tag{5}
\]

The last implication uses gcd(P,Q)=1. A nonzero derivative has
smaller total degree than P, so (5) forces P_s=0. The unit Jacobian
then gives -P_U Q_s=c; both factors must be nonzero constants.
Hence

\[
 P=\alpha U+\beta,\qquad Q=\gamma s+q_0(U),
                         \qquad \alpha\gamma\ne0,            \tag{6}
\]

with an explicit triangular polynomial inverse. This proves the
claim. The source coordinate eta is used for ring identities only;
no Hessian invariance under that nonlinear coordinate change is used.

The argument needs neither a closed embedding nor Tang, once the
cotangent form and physical base data are given. It allows all degrees
of D,v,h,H and all mixed terms away from the common-zero-base scheme.
It does not apply to a variable eta coefficient row or to a shifted
pure fibre P=alpha without an additional proof.

## 2. The allowed gradient projection must now be transverse (HC4MYHL1)

Let i:A4->A5 be an actual polynomial closed embedding, f=Psi composed
with i, and assume the actual four-variable Hessian of f is a unit.
Suppose the two physical base pullbacks x composed with i and
y composed with i are independent **homogeneous constant linear
combinations of grad f**. No additive offsets are included.

Under hypothetical nonautomorphy, Tang supplies a constant cotangent
form as in Section 1. Write its gradient outputs as (eta,K), and let
U_out=span(K1,K2) be the intrinsic cotangent output plane. If V is
the span of the two physical base pullbacks, the necessary condition is

\[
                           V\cap U_{\rm out}=0.              \tag{7}
\]

Indeed, if this intersection has dimension two, the physical bases
are an invertible linear target change of K. The stable hypersurface
obstruction HC4MYHJ1 makes K invertible. If it has dimension one,
constant row operations put the common zero scheme into (2), with
v a constant linear combination of K. Equation (1) supplies the
forbidden nonlinear composite boundary, and HC4MYHL0 again makes
K invertible. The latter implication does not require the embedding
hypothesis at all; an arbitrary substitution with those identities
already has an invertible gradient under unit Hessian.

Thus the full pair must have eta rank two. It is then an actual
source coordinate pair by HC4MYHK1. The previous generator profile
remains correct: rank-one projections do admit abstract A5 embeddings,
but none with this fixed potential identity and a nonautomorphic
unit-Hessian gradient. This closes a genuinely additional stratum.

The more general base data (2), with arbitrary v(s,U), are also
closed whenever a fixed-potential realization has that common zero
scheme. General physical bases need not be gradient projections,
and rank-two eta projections are not claimed solved by (7).

## 3. The ambient critical value transfers to at most one point (HC4MYHL2)

Put

\[
                 p_*=-M/(2L),\qquad \lambda_*=-M^2/(4L),
                 \qquad o=(0,0,p_*,0,0).                    \tag{8}
\]

For any closed embedding i and unit-Hessian f=Psi composed with i,
suppose that at every full critical point of f the physical base
pullbacks vanish, and their differentials have rank two. Then

\[
 \#\{\nabla f=0,\ f=\lambda_*\}
          =\begin{cases}1&o\text{ belongs to the image of }i,\\
                         0&o\text{ does not.}\end{cases}     \tag{9}
\]

In particular the count is at most one. A sufficient condition on
the physical bases is b=phi(grad f), where phi:C4->C2 is polynomial,
phi(0)=0 and Dphi(0) has rank two. This allows arbitrary nonlinear
phi, as well as the linear projections of Section 2. The claim is
not an unrestricted singleton-zero properness statement for K.

At a source critical point with value lambda_star, (1) forces p=p*.
The native first derivatives on x=y=0 are

\[
 \Psi_p=2Lp+M,\quad \Psi_q=\Psi_r=0,\quad
 \Psi_x=3(2Lp+M)q+2Nr,\quad \Psi_y=Nq.                       \tag{10}
\]

Thus at p=p*, df=N(2r db1+q db2). Rank two of the base
differentials forces q=r=0. The image point is o, and injectivity
of the closed embedding gives at most one source point. Conversely,
all ambient first derivatives vanish at o, so a preimage of o is a
critical point of f with value lambda_star. This proves (9).
All source critical points are reduced under unit Hessian, so this
also identifies their reduced schemes at that value.

## 4. A global necessary square equation for the remaining transverse case

For eta rank two write the actual physical base pair as

\[
                           (x,y)^T=A_0K+C_0\eta,
                           \qquad C_0\in GL_2(C).
\]

Let T=C0^(-1)A0 and retain b=(DK)^(-1)K from (3). The full
zero-base fibre is eta=-TK and has coordinates z=(s,U). Its
restriction of (4) must satisfy

\[
       (1-b\cdot\nabla)D-b^T T K-\lambda_*
                              =L h(z)^2                     \tag{11}
\]

for an actual polynomial h, namely the restricted physical p minus
p*. This is a global polynomial identity, not just equality modulo
the zero-fibre ideal of K or a condition on its finitely many points.
No assertion solves (11) for arbitrary Keller data.

For each z_i in K^(-1)(0), b(z_i)=0 and Db(z_i)=I. Differentiating
the left side before subtracting lambda_star therefore gives zero
at z_i. Hence

\[
 D(z_i)-\lambda_*=L h(z_i)^2,\qquad
                           h(z_i)\nabla h(z_i)=0.            \tag{12}
\]

When D(z_i)!=lambda_star, z_i must be a critical point of h.
When D(z_i)=lambda_star, (9) says at most one such z_i can occur.
These necessary conditions do not establish a global embedding,
terminate a formal series or prove properness at another value.

The nonlinear eta coefficient case outside Section 1, translated
physical bases not satisfying the centered hypotheses, general bases
not given by gradient projections, and the remaining transverse
fixed-potential problem stay distinct. The universal SL2 row equation,
unit-ideal plane surjectivity, singleton-zero properness and unrestricted
HC4/JC2 are not resolved.

The [exact checker](verify_fixed_projection.py) and
[receipt](receipt-fixed-projection.json) check the cotangent identities,
boundary degree mechanism, native critical jets, explicit triangular
inverses and the remaining square equation on finite controls. The
arbitrary-degree conclusions are the written proofs above. The only
retained structural input beyond elementary algebra is Tang for
passing from an arbitrary HC4 potential to its cotangent alternative,
and the already proved stable hypersurface theorem for the pure pair.
No new external classification, Hessian invariance under nonlinear
source changes, formal verification, independent replay, external
review or novelty priority is claimed.
