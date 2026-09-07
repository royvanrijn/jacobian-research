# The transporter is part of the CT difference extension

Follow-up: [the compatible Kummer diagram](JACOBIAN_LOCAL_CONDITIONS_AND_CT.md)
now identifies the local Baer-sum condition below with the ordinary Jacobian
condition and certifies the local connecting ranks at \(u=2\). The gap
recorded here is retained as the previous stage of the argument.

The genus-two Jacobian's 2-torsion gives a concrete bridge between three
previously separate observations: the changing 4-division fields, the
affine transporter class, and the changing CT pairings.

For every retained nonzero fixed-cubic control,
\[
0\longrightarrow V\longrightarrow J_u[2]\longrightarrow V\longrightarrow0,
\qquad V=E_0[2]\simeq E_u[2],
\]
is the difference of the two elliptic 4-torsion extensions. The
six-branch-point model determines its Galois action explicitly. It is an
indecomposable four-dimensional module, with a degree-48 splitting field.
The same abstract module type occurs in all six controls.

The associated CT comparison is obtained by retaining the **local
conditions** in that difference construction. This identifies which
higher-descent object produces the pairing difference; it does not yet
compute the CT entries directly from the new representation.
The result belongs to **solubility-obstruction structure**, not
point-search visibility or a numerical rank predictor.

## The difference of the two four-torsion extensions

Use the [labelled gluing](JACOBIAN_GLUING_AND_SHA_BLOCKS.md)
\[
J_u=(E_0\times E_u)/\operatorname{graph}(\varphi).
\]
A class \([P,Q]\) is in \(J_u[2]\) exactly when
\(2P=\varphi^{-1}(2Q)\in V\). Therefore
\[
J_u[2]\simeq
\{(P,Q)\in E_0[4]\times E_u[4]:2P=\varphi^{-1}(2Q)\}
/\operatorname{graph}(\varphi|_V).
\]
The map to \(V\) sends \([P,Q]\) to \(2P\); its kernel is
\((V\times V)/\operatorname{diag}V\simeq V\).
This is the usual pullback/pushout description of the Baer sum.
Sum and difference agree here because the end modules are killed by two.
The middle group is also killed by two, although each input middle group
is \((\mathbb Z/4)^2\).

For a coordinate check, let \(A\) be a lift of the common mod-two action
\(g\), and write the other action as \(B=(1+2c)A\).
Represent a pair of equal-parity vectors by
\[
v=x\bmod2,\qquad k=(y-x)/2\bmod2.
\]
Then the induced action on the quotient is
\[
(k,v)\longmapsto(gk+cgv,\ gv).
\]
The certificate checks this formula on all 64 pair representatives for
each of the 48 signed root permutations. It also checks the same action
independently in the six-Weierstrass-point model below.

## Three sign characters and the common shift

The branch points of \(C_u\) are
\[
z=\pm\sqrt{\gamma_i},\qquad
\gamma_i=1-u\theta_i,\qquad i=1,2,3.
\]
Recall \(\prod_i\gamma_i=D(u)\).
Model \(J_u[2]\) as even subsets of these six points modulo complements.
The three full pairs give the submodule \(V\). Parity in each pair gives
the quotient \(V\), with labels \(T_1,T_2,T_3\) satisfying
\(T_1+T_2+T_3=0\).

Choose the positive root in each pair to define a linear section of this
quotient. If a Galois element permutes the pairs by \(g\) and flips the
target signs by \(\epsilon=(\epsilon_1,\epsilon_2,\epsilon_3)\), then
\[
c(\epsilon)=\sum_i\epsilon_iN_i,\qquad
N_i(v)=e_2(T_i,v)\,T_i.
\]
Here the Weil pairing is encoded additively in \(\mathbb F_2\).
For \(T_1=(1,0),T_2=(0,1),T_3=(1,1)\), the matrices are
\[
N_1=\begin{pmatrix}0&1\\0&0\end{pmatrix},\quad
N_2=\begin{pmatrix}0&0\\1&0\end{pmatrix},\quad
N_3=\begin{pmatrix}1&1\\1&1\end{pmatrix}.
\]
They are linearly independent, have trace zero, and sum to \(I\).
Conjugation by \(S_3\) permutes them as it permutes the roots.
For example, flipping just the first branch pair acts by
\((k,v)\mapsto(k+N_1v,v)\), fixing the submodule and quotient
individually but not the extension.

Put \(\chi_D=\epsilon_1+\epsilon_2+\epsilon_3\) and
\(\epsilon_{\eta,i}=\epsilon_i+\chi_D\). Then
\[
\boxed{c(\epsilon)=\chi_DI+\sum_i\epsilon_{\eta,i}N_i,\qquad
\sum_i\epsilon_{\eta,i}=0.}
\]
In cubic squareclass coordinates, the even component is precisely
\[
\eta_u=D(u)(1-u\theta),\qquad N(\eta_u)=D(u)^4.
\]
Thus the affine class discovered by the rational-point transporter is
also the non-scalar component of this explicit torsion-difference
cocycle. The identity is verified for all eight sign vectors; the norm
identities are replayed for all six controls.

This is a decomposition in the stated root-section normalization. It
must not be promoted to two independent CT scores. Changing the
4-torsion bases can change cocycle representatives. There is even a
useful warning in dimension two: for an order-three matrix \(H\),
\(gHg^{-1}-H\) is zero on even root permutations and \(I\) on odd
permutations. Hence the scalar cocycle of the cubic discriminant is a
coboundary in \(\operatorname{End}(V)\). The full extension, rather than
an unqualified scalar-component label, is the invariant object.

## The common module is genuinely indecomposable

The earlier simple-collision witnesses prove the three \(\gamma_i\)
independent modulo squares even over \(\mathbb Q(E_0[4])\), hence over
the smaller splitting field \(L=\mathbb Q(E_0[2])\). Thus
\[
\mathbb Q(J_u[2])=L(\sqrt{\gamma_1},\sqrt{\gamma_2},\sqrt{\gamma_3}),
\qquad[\mathbb Q(J_u[2]):\mathbb Q]=8\cdot6=48.
\]
The six-point model is faithful: distinct signed permutations give
distinct actions on \(J_u[2]\), as the certificate verifies. Its image
is the full group \((C_2)^3\rtimes S_3\).

Exact linear equations for the matrices commuting with this group give
\[
\operatorname{End}_{G_{\mathbb Q}}(J_u[2])
=\{aI+bN:a,b\in\mathbb F_2\},\qquad
N(k,v)=(v,0),\quad N^2=0.
\]
This is the dual-number algebra \(\mathbb F_2[N]/(N^2)\). Its only
idempotents are zero and identity, proving that \(J_u[2]\) has no
nontrivial direct-sum decomposition as a Galois module.
This is an arithmetic block assertion, unlike separately putting a CT
matrix into symplectic normal form.

The bound does not vary across the six controls. In particular neither
the degree 48 nor this module type distinguishes their different CT
matrices or detects rational solubility. A four-dimensional coefficient
module can govern a pairing on many global cohomology classes; its
dimension does not bound CT rank by four.

## The precise CT comparison, with local data retained

Write \(L_{0,v},L_{u,v}\) for the local 2-Kummer images, and set
\[
C_v=L_{0,v}\cap L_{u,v},\qquad D_v=L_{0,v}+L_{u,v}.
\]
Start with the two exact 4-torsion sequences decorated with their local
4-Kummer conditions. Pull each right endpoint back to \(C_v\) and push
each left endpoint out to \(D_v\). Their Baer sum is an exact sequence
of Selmer objects
\[
0\longrightarrow(V,D)\longrightarrow(J_u[2],\mathcal W_{\rm sum})
\longrightarrow(V,C)\longrightarrow0.
\]
The underlying module is the one proved above. Local self-duality gives
\((V,D)^\vee=(V,C)\). Naturality and Baer-sum additivity give
\[
\boxed{\operatorname{CT}_{\rm sum}(x,y)
=\operatorname{CT}_{E_0}(x,y)+\operatorname{CT}_{E_u}(x,y)}
\quad(x,y\in S_0\cap S_u).
\]
The relevant general results are
[Morgan–Smith, *The Cassels–Tate pairing for finite Galois modules*, Theorem 1.3, Definition 4.3 and Proposition 4.4](https://arxiv.org/pdf/2103.08530).
On \(W_u\), the anchor pairing vanishes because its classes are rational.
The new decorated extension therefore produces exactly the retained
restricted \(E_u\) pairing.

Crucially, \(\mathcal W_{\rm sum}\) includes local lifting choices from
the 4-Kummer conditions. The underlying module and the two endpoint
spaces alone have not been proved sufficient to reconstruct it.
This note does not identify \(\mathcal W_{\rm sum}\) with the ordinary
local Kummer condition on \(J_u[2]\); such an identification needs its
own compatible Kummer diagram. Consequently no ordinary genus-two
Selmer dimension or new CT entry is asserted by this computation.

This also resolves an apparent compression paradox. The previous
10+2 simultaneous CT decomposition cannot factor through three linear
class features. Here the three sign characters define an extension
coefficient module, followed by global cohomology and local lifting.
They are not three linear functions on the common Selmer space.
No contradiction with the earlier rank obstruction arises.

## What changed, and the next falsifiable experiment

The chain now has an explicit intermediate object:
\[
\text{branch signs of }1-u\theta_i
\Rightarrow\text{torsion-difference extension containing }\eta_u
\Rightarrow\text{CT difference after local decoration}.
\]
This strengthens the link between geometry and simultaneous solubility
obstructions. It still stops before rational solubility: surviving CT
classes need not represent points.

The next experiment should reconstruct the local decoration at one fixed
control, preferably \(u=2\), and compare a small preselected CT submatrix
with the retained Fisher values. A successful independent reconstruction
would test the proposed explanation of CT variation. If it fails, the
failure should distinguish an incorrect local-condition identification
from insufficient global cochain data. No new parameter or point search
is needed for that test.

For Agent 1, an abstract division-field degree or common module type is
still weak evidence. A fully specified extension with computable local
lifts could provide a **solubility** feature measuring obstructions, provided its
candidate class space is constructed before exceptional points are
supplied. Neither feature is presently a prospective rank selector.

## Replay

The [protocol](TORSION_DIFFERENCE_PROTOCOL.json) and
[certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_torsion_difference_v1.json)
retain all 48 matrices, the scalar/even decomposition, commuting
algebra, six field witnesses and source hashes.

```sh
python3 elliptic-curves/rank-jump/torsion_difference.py check
```

The replay checks 768 root actions, 3072 mod-four pair actions, closure
of the 48-element image, and all commuting-algebra equations. The
simple-collision field hypotheses are the previously certified inputs;
no degree-48 number field is constructed. No active search file or
status entry was changed.
