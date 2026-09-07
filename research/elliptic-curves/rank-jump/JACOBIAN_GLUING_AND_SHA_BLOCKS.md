# One Jacobian carries the whole inherited Sha obstruction

For each of the six fixed-cubic controls, the inherited locally soluble
classes are killed in **one explicit genus-two Jacobian**. The existing CT
matrices certify independent Sha subspaces of dimensions
\(16,12,16,12,12,14\), respectively, inside this common kernel.
This supplies a shared low-degree arithmetic object for a large block of
obstructions; it does not supply rational points on the obstructed covers.

The same construction identifies a positive specialization mechanism:
at the anchor \(u=0\), the associated abelian-variety extension splits and
its connecting obstruction vanishes on the entire anchor group at once.
All six nonzero controls have nonsplit extensions, certified by small-prime
nonisogeny witnesses. Nonsplitting alone does not rule out isolated soluble
classes or high rank.

These are **solubility** and **solubility-obstruction** statements. In the
literature, killing Sha in an abelian variety is called “visibility.”
That term here has no relation to the project's **point-search visibility**
layer. The computation uses no rational point search.

## An explicit genus-two model with the correct labels

Keep \(f(a)=a^3+Aa+B\), \(D=1+Au^2+Bu^3\), and the
[fixed-cubic pencil](LINEAR_TWIST_SOLUBLE_BLOCKS.md). For \(uD\ne0\),
write its auxiliary curve as
\[
C_u:\quad b^2=f(a),\qquad z^2=1-ua.
\]
With \(Y=u^2b\), an even sextic model is
\[
Y^2=g_u(z)=uD-(3u+Au^3)z^2+3uz^4-uz^6.
\]
Its polynomial discriminant is
\[
\operatorname{disc}(g_u)=64u^{22}D\,\operatorname{disc}(f)^2\ne0.
\]
The two degree-two maps are
\[
\pi_0(z,Y)=\left(a,\frac{Y}{u^2}\right),\quad
\pi_u(z,Y)=
\left(\frac{a+(Aa+B)u^2}{z^2},\frac{YD}{u^2z^3}\right),
\qquad a=\frac{1-z^2}{u}.
\]
Their involutions are \((z,Y)\mapsto(-z,Y)\) and
\((z,Y)\mapsto(-z,-Y)\). They extend over the displayed poles on the
smooth projective curves.

The certificate checks both quotient equations, the sextic discriminant,
and the following labelled 2-torsion identity in \(K=\mathbb Q(\theta)\):
\[
\frac{\theta+(A\theta+B)u^2}{1-u\theta}
 =\theta+u\theta^2=\alpha_u.
\]
Thus the gluing below uses exactly the cohomology identification of the
retained local and CT calculations, not an unverified permutation.

Let \(J_u=\operatorname{Jac}(C_u)\). The complementary double covers give
dual degree-four isogenies
\[
\Phi:E_0\times E_u\longrightarrow J_u,\quad
(P,Q)\longmapsto\pi_0^*P+\pi_u^*Q,\qquad
\Psi=(\pi_{0*},\pi_{u*}):J_u\longrightarrow E_0\times E_u,
\]
with both composites equal to multiplication by two.
This is the classical bielliptic construction; see
[Wetherell, Chapter 2, Theorem 4.2](https://swc-math.github.io/aws/1999/99WetherellThesis.pdf).

The kernel of \(\Phi\) is the graph of
\(\varphi:E_0[2]\to E_u[2]\) specified by \(\theta\mapsto\alpha_u\).
To see the labels directly, take the pair of Weierstrass points with
\(z^2=1-u\theta\). Their divisor represents the pullback of the
corresponding 2-torsion point from either elliptic quotient. The two
origin fibres are the divisors above \(z=\infty\) and \(z=0\);
their difference is the principal divisor of \(z\).
The two pullback classes therefore agree. This gives four points in the
kernel, its full order. Consequently
\[
J_u\simeq(E_0\times E_u)/\operatorname{graph}(\varphi).
\]

## The connecting map is exactly the inherited torsor

There is an exact sequence of abelian varieties
\[
0\longrightarrow E_u\xrightarrow{i_u}J_u
\xrightarrow{q_u}E_0\longrightarrow0,\qquad
q_u([P,Q])=2P.
\]
Its rational-point connecting homomorphism is
\[
\partial_u:E_0(\mathbb Q)/2E_0(\mathbb Q)\longrightarrow H^1(\mathbb Q,E_u).
\]
For \(P\in E_0(\mathbb Q)\), choose \(R\) with \(2R=P\).
The cocycle of the lift \([R,0]\) is
\(\varphi(\sigma R-R)\) in the embedded \(E_u\). Hence
\(\partial_u(P)\) is exactly the image of the anchor Kummer class under
\[
H^1(\mathbb Q,E_0[2])
\xrightarrow{\varphi_*}H^1(\mathbb Q,E_u[2])
\longrightarrow H^1(\mathbb Q,E_u).
\]
In particular, for \(\beta\in W_u=W\cap\operatorname{Sel}_2(E_u)\),
its image lies in \(\Sha(E_u)[2]\) and is killed by \(i_{u*}\).
This is also the elementary gluing argument in
[Fisher, *Visible 2-torsion in the Tate–Shafarevich group*, Lemma 2.1](https://www.dpmms.cam.ac.uk/~taf1000/papers/visible2torsion.pdf).
The application here is simultaneous to the whole inherited subspace.

Exactness gives a useful equivalence for every such anchor point:
\[
\boxed{\text{its inherited cover has a rational point}
\ \Longleftrightarrow\
P\in q_u(J_u(\mathbb Q)).}
\]
This is a rational divisor-class lifting condition on one Jacobian.
It is broader than finding individual points on \(C_u\) in the earlier
square-condition test. Merely knowing that \(J_u\) is isogenous to
\(E_0\times E_u\) does not make \(q_u\) surjective on rational points.
Indeed \(\pi_0^*(P)\) always projects to \(2P\), supplying only the
trivial class in \(E_0(\mathbb Q)/2E_0(\mathbb Q)\).

## Certified block sizes and projection obstructions

For each retained CT matrix on \(W_u\), exact symplectic elimination
selects a nondegenerate subspace \(N_u\) of dimension equal to its rank.
The map \(N_u\to\Sha(E_u)[2]\) is injective: if a class mapped to zero,
it would be rational and pair trivially with the whole Selmer group,
contradicting nondegeneracy on \(N_u\).
All its images are killed in the same \(J_u\).

Let \(K_u=\ker(\partial_u|_{W_u})\). It is the space of inherited
classes with rational projection lifts. Then \(K_u\subseteq\operatorname{rad}B_u\).

| \(u\) | \(\dim W_u\) | Independent Sha block killed in \(J_u\) | \(\dim K_u\) at most | \([W_u:K_u]\) at least |
|---:|---:|---:|---:|---:|
| -3 | 17 | 16 | 1 | 65536 |
| -2 | 13 | 12 | 1 | 4096 |
| -1 | 18 | 16 | 2 | 65536 |
| 1 | 13 | 12 | 1 | 4096 |
| 2 | 13 | 12 | 1 | 4096 |
| 3 | 15 | 14 | 1 | 16384 |

These bounds are unconditional and confined to the inherited subspaces;
they do not need a complete Selmer group or finiteness of Sha.
All selected masks and symplectic pairs are retained. For example, at
\(u=-1\) the eight paired anchor masks are
\[
(1,10),(6,24),(88,149),(413,1150),
(850,17709),(2322,82425),(5253,196762),(12191,39308).
\]
Their CT matrix is eight standard alternating \(2\times2\) blocks,
and their sixteen images are independent in \(\Sha(E_{-1})[2]\).
This chosen symplectic basis is not a canonical decomposition of the
arithmetic object. The common Jacobian is the structural statement.

## A splitting event that removes the entire obstruction

More generally, the gluing extension splits over \(\mathbb Q\) if and
only if there is a homomorphism \(h:E_0\to E_u\) whose restriction to
2-torsion is \(\varphi\).
For the forward direction, a section \(s:E_0\to J_u\) gives
\(\Psi s=(1,h)\); applying \(\Phi\Psi=[2]\) to \(E_0[2]\) forces
\((P,hP)\in\ker\Phi\). Conversely, if such \(h\) exists, then
\[
s(P)=[R,hR],\qquad 2R=P,
\]
is independent of the choice of \(R\), descends over \(\mathbb Q\),
and satisfies \(q_us=1\). Such \(h\) is necessarily an odd-degree
isogeny. Thus nonisogeny rules out extension splitting.

The bounded test compared Frobenius counts at the frozen good primes.
The first distinguishing counts are:

| \(u\) | Good prime | \(\#E_0(\mathbb F_p)\) | \(\#E_u(\mathbb F_p)\) |
|---:|---:|---:|---:|
| -3 | 11 | 18 | 14 |
| -2 | 11 | 18 | 12 |
| -1 | 11 | 18 | 8 |
| 1 | 11 | 18 | 10 |
| 2 | 23 | 32 | 24 |
| 3 | 11 | 18 | 14 |

The curves therefore are not \(\mathbb Q\)-isogenous, and none of these
six extensions splits. At \(u=2\), the preceding tested primes 11 and 19
had equal counts; these equalities alone would prove nothing.
The entire retained count record has an independent Sage replay.

At \(u=0\), use the compactified double cover rather than the singular
sextic coordinates. Its two components \(z=1\) and \(z=-1\), both \(E_0\),
meet at their origins. This is a compact-type genus-two curve with
Jacobian \(E_0\times E_0\). The maps become
\[
q_0(P,Q)=P+Q,\qquad i_0(R)=(R,-R).
\]
The section \(P\mapsto(P,0)\) splits the extension. Equivalently,
\(h=1\) lifts the identity on \(E_0[2]\). The connecting map vanishes
on all of \(E_0(\mathbb Q)\), so the whole known 20-dimensional
anchor Kummer space is rationally soluble at once.

Together with the already proved generic arithmetic rank zero, this
gives an exact controlled chain:
\[
u=0\ \Rightarrow\
\text{split gluing extension}\ \Rightarrow\
\partial_0=0\ \Rightarrow\
\text{twenty independent rational directions}.
\]
It explains this engineered high-jump control before supplying its
individual exceptional points, although their known independence is
needed to count the released directions. It does not show that
unrelated MW17/MW16 high fibres satisfy the same condition.

## Mechanisms, exclusions and the next gap

1. **Solubility: extension splitting is a sufficient block mechanism.**
   It makes a whole connecting homomorphism vanish, rather than making
   isolated covers accidentally soluble. The anchor degeneration realizes
   it exactly.
2. **Solubility obstruction: a common Jacobian can carry a large Sha block.**
   All six controls demonstrate this. A shared auxiliary curve or a
   product isogeny by itself is therefore weak evidence for rational
   Mordell–Weil gain.
3. **Incidence: the relative full-Selmer theorem remains complementary.**
   At \(u=-1\), full Selmer dimension drops by just one from the anchor,
   while the common Jacobian carries the certified 16-dimensional Sha
   block. The large difference in rational solubility is not a point-search
   effect.
4. **Missing implication:** a nonsplit extension may still have many
   rational projection lifts. Explaining a large kernel of \(\partial_u\)
   without already knowing the exceptional points is the remaining
   arithmetic problem. The present computation does not derive the CT
   matrix entries from the gluing geometry.

The next useful experiment is to test whether the known CT variation
can be expressed through a concrete obstruction to rational divisor-class
lifts in this extension. A genus-two point-search expansion would test
a narrower visibility endpoint. For Agent 1, a certified lifting map or
equation-level splitting criterion would be a **solubility** feature;
the mere existence of the Jacobian is not a rank selector.

## Certificate and replay

The [frozen protocol](JACOBIAN_SHA_BLOCK_PROTOCOL.json) and
[certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_jacobian_sha_blocks_v1.json)
pin the six models, exact maps, branch labels, symplectic masks and finite
counts. From the repository root:

```sh
python3 elliptic-curves/rank-jump/jacobian_sha_blocks.py check
sage -python elliptic-curves/rank-jump/jacobian_sha_blocks.py verify
```

No new CT pairing, parameter, rational point, class group or live-search
output was computed or changed. The initial certificate before adding the
bounded nonisogeny test is retained under
`artifacts/local/rank-jump-jacobian-sha-block-v1/before-isogeny-extension/`.
