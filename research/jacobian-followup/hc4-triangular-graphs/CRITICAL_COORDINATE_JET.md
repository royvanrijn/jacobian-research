# A critical-coordinate first jet closes all q reparameterizations

14 November 2026. **HC4MYGL2** closes the entire q-graph chart after
every polynomial base automorphism and every relative polynomial dual
automorphism, in all degrees and word lengths. It also gives another
proof of that conclusion for the r chart. Only a coordinate-pair
specialization at one base point is needed; global dual invertibility
is stronger than the hypothesis used here.

The new input **HC4MYGL0** is Gwozdziewicz's
[Injectivity on one line](https://arxiv.org/abs/alg-geom/9305008),
Theorem 1.1: a plane Keller map injective on an affine line is an
automorphism over an algebraically closed characteristic-zero field.
Its Lemma 2.3 includes the case of a map fixing a line pointwise.
The [import record](sources/critical-jet-import.json) pins the actual
[primary manuscript](sources/gwozdziewicz-injectivity-on-one-line-v1.pdf).
This theorem is imported, not independently proved by our replay.

## 1. A first-jet criterion

Write the source as (s,t,P,Q), with base origin s=t=0. Suppose

\[
 f(0,0,P,Q)=H(a(P,Q)),\qquad \deg H\geq2,          \tag{1}
\]

where (a,b) is a polynomial coordinate pair of the dual plane. Suppose
some nonzero constant base direction d satisfies

\[
 (\partial_d f)(0,0,P,Q)
   =H'(a)c(P,Q)+\lambda b(P,Q)+\mu(a),\qquad
 \lambda\ne0,                                   \tag{2}
\]

with c and mu polynomial. The hypotheses concern only the full dual
trace and first base jet at this one base point; no assertion is made
about a whole divisor or formal completion.

**HC4MYGL1.** If f has nonzero constant Hessian determinant, its
gradient is a polynomial automorphism.

Work over C, then descend the polynomial inverse to any characteristic-zero
coefficient field. The reproduced [Tang alternative HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md)
supplies a jointly affine constant two-plane under hypothetical
nonautomorphy. We exhaust its projection rank to the base.

### Rank zero

The plane is the pure-dual plane, so f is globally affine in P,Q.
This contradicts (1), because a is nonconstant and deg H is at least two.

### Rank one

Constant linear base and dual changes, with a constant base-linear
translation of the duals, put the plane in the span of partial_S and
partial_V. The base origin remains S=T=0. Joint affineness gives

\[
 f=B(T,U)+S A(T,U)+V C(T,U),\qquad
 \det\operatorname{Hess}f=\operatorname{Jac}_{T,U}(A,C)^2. \tag{3}
\]

At the base origin this is affine in V and also equals H(a0(U,V)).
Polynomial V degree in (1) forces a0 independent of V and C(0,U)=0.
Thus C=T*k(T,U). On T=0 the unit Jacobian reads

\[
 -A_U(0,U)k(0,U)\in C^\times.                    \tag{4}
\]

Both factors are nonzero constants. Hence (A,C) restricts on T=0
to (alpha U+beta,0), alpha nonzero, and is injective on that line.
HC4MYGL0 makes (A,C) an automorphism. The full gradient then has a
polynomial inverse: its S,V outputs recover T,U through (A,C), and
its T,U outputs recover S,V by the inverse transpose of D(A,C),
after subtracting grad B. This rank does not need (2).

### Rank two

A constant base-linear translation of the duals makes the plane the
base plane. In the adapted dual variables U,V,

\[
 f=H(a(U,V))+s F(U,V)+t G(U,V),\qquad
 \operatorname{Jac}(F,G)\in C^\times.             \tag{5}
\]

The translation agrees with the original dual variables at base zero.
Its correction to a first base derivative of f is a constant dual
directional derivative of H(a), and therefore is divisible by H'(a).
Consequently (2) is preserved, with a modified c, and says

\[
 d_1 F+d_2 G=H'(a)\widetilde c+\lambda b+\mu(a).  \tag{6}
\]

Choose a root xi of H'; it exists over C because deg H>=2. Compose
the **plane pair** (F,G) with the inverse of its source coordinate
pair (a,b). Its Jacobian remains a nonzero constant. On the straight
line a=xi, the linear combination in (6) equals

\[
 \lambda b+\mu(\xi),                            \tag{7}
\]

so the pair is injective on that entire line. HC4MYGL0 makes the
pair, and hence the original (F,G), an automorphism. Equation (5)
then gives the polynomial gradient inverse by the cotangent argument
used in rank one. This proves HC4MYGL1.

The nonlinear dual normalization above is used only to test the
two-variable Keller map. The four-variable Hessian argument uses
constant affine source changes throughout. We do not assume that
a nonlinear source change preserves a Hessian determinant.

## 2. The actual Meng--Yang first base jet

Use the original ambient coordinates (x,y,p,q,r), u=1+xy and

\[
 A=u^3p+3xu^2q-x^3r,
\]
\[
 B=y^2u(4+3xy)p+(y+3xy^2(4+3xy))q+(2x-3x^2y)r,
 \qquad \Psi=LA^2+MA+NB,\quad LN\ne0.
\]

For arbitrary polynomial substitutions, take a parameter base point
mapping to old (x,y)=(0,0). Write p0,q0,r0 for the specialized dual
polynomials there. The full chain-rule identity is

\[
 f_0=Lp_0^2+Mp_0,\qquad
 df|_0=(2Lp_0+M)(dp|_0+3q_0\,dx|_0)
              +Nq_0\,dy|_0+2Nr_0\,dx|_0.       \tag{8}
\]

All first jets of the substitutions are retained in dp. The terms
involving dq and dr vanish at this base point; they are not discarded
by an approximation.

### Graph eliminating q

Let (x,y)=beta(s,t) be any polynomial base automorphism, let (p,r)
be any relative polynomial dual automorphism, and let q be an arbitrary
polynomial graph function in the parameters. Translate beta^(-1)(0,0)
to the parameter origin by a constant affine base change.

The specialized pair (p0,r0) is an automorphism of the dual plane.
Choose the constant parameter base direction d with dx(d)=1,dy(d)=0
at the origin. This exists because D beta there is invertible. Then

\[
 (\partial_d f)_0=(2Lp_0+M)((\partial_d p)_0+3q_0)
                       +2Nr_0.                 \tag{9}
\]

This is exactly (1)-(2) with H(z)=Lz^2+Mz, a=p0, b=r0, lambda=2N,
mu=0. HC4MYGL1 closes the entire q chart. No straight inverse-base
image of y=0 is needed, and the retained r need not be affine on it.

### Graph eliminating r

With arbitrary base and retained (p,q) dual automorphisms, specialize
the pair (p0,q0) and choose d with dx(d)=0,dy(d)=1. Equation (8) gives

\[
 (\partial_d f)_0=(2Lp_0+M)(\partial_d p)_0+Nq_0. \tag{10}
\]

Apply HC4MYGL1 with b=q0 and lambda=N. This gives a different proof
of the all-base/all-dual r obstruction in HC4MYGC2. It uses the line
theorem instead of a rational coordinate trace and Rentschler.

**HC4MYGL2** comprises both conclusions. In fact the proof only needs
beta to have a point over old (0,0) with invertible first base Jacobian,
and the corresponding specialized retained dual pair to be an actual
polynomial automorphism. The full base and retained-dual maps may be
arbitrary polynomial substitutions away from that slice. A merely
Keller specialized pair is not substituted for a coordinate pair.

## 3. What this changes and what remains

The q and r graph charts are now closed under all polynomial base and
relative polynomial dual automorphisms. Longer fibre-preserving words
and higher degrees do not reopen either chart. The elementary
[square-axis proof](Q_GRAPH_SQUARE_AXIS.md) remains useful for its
stronger restriction on external inputs: on its smaller class it
needs neither Rentschler nor the line-injectivity theorem.

The p chart differs at the very first slice: its p0 is a graph
function, not a retained coordinate. The
[fixed-graph theorem HC4MYGP1](P_GRAPH_SINGLETON_ENDPOINT.md) still
identifies its full nonlinear-dual branch exactly with the unresolved
singleton-zero plane Keller properness problem. No theorem here makes
an arbitrary singleton-zero map injective on a line. The missing
finite-normalization fibre length d-1 remains a global obligation.

The later [HC4MYGB2](BASE_ELIMINATED_GRAPHS.md) closes both native
base-eliminated charts. Nonlinear reparameterizations of those charts,
genuinely nonlinear base/dual mixing and uncovered embeddings remain
separate. The determinant hypothesis
always concerns the actual pulled-back potential, and a marked
collision still needs the actual pulled-back covectors to match.

Condition (2) is sufficient, not automatic from (1) and unit Hessian.
For example f=P^2+s*(Q+(P+Q^2)^2)+t*(P+Q^2) has unit Hessian and
an invertible gradient. At the critical line P=0 its two first base
coefficients are Q+Q^4 and Q^2; no nonzero constant linear combination
is affine nonconstant in Q. A coordinate mate of P is cQ+k(P), so
changing that mate does not restore (2) for this example. Its
invertibility instead follows from the explicitly automorphic
coefficient pair. This control prevents treating the new jet
condition as a proof of unrestricted JC2.

## Reproduction and assurance

The [checker](verify_critical_jet.py) verifies (8)-(10), both normal
form determinants, the first-jet correction under a constant dual
translation, critical-line injectivity identities, nonlinear base
and dual specialization controls, and explicit gradient inverses.
The all-rank argument is a written proof; the exact controls do not
formally certify its generality or prove the imported line theorem.
[Actual receipt](receipt-critical-jet.json).

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_critical_jet.py
```

Tang and Gwozdziewicz are explicit imports. No new properness theorem,
unrestricted HC4 or JC2 solution, formal verification, independent
external review or novelty priority is claimed.
