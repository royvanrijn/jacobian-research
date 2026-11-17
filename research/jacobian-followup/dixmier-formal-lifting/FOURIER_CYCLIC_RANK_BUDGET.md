# Fourier rank constraints for cyclic branch extraction

14 November 2026. **QKR1.** Higher generic inverse-Fourier rank can
escape the rational-line obstruction in
[QFT1](RATIONAL_LINE_FOURIER_OBSTRUCTION.md), but branch extraction
has a quantitative rank cost. The cost is independent of coefficient
degrees, finite pole numbers and orders, and the chosen cyclic vector.

For a rank-N irreducible connection that becomes reducible on a cyclic
power cover, its inverse-Fourier rank r satisfies N<=2r. With at most
one finite singularity, N>2r implies irreducibility over every finite cover. More precisely,
a Kummer self-twist of order e costs at least (1-1/e)N plus the
irregularity at the branch point. A finite pushforward of a rank-m
connection along q=v^d/d costs at least (d-1)m.

An explicit all-degree family attains algebraic splitting while
retaining the classical phase equation of QPH1 on q!=0. For a polynomial
g of degree m its Fourier connection has rank2m and inverse-Fourier
rank m+1. For Long, m=2 and the inverse rank is three. This is a
positive branch-splitting control, not a full quantum carrier:
the extracted quadratic operator fails the finite-mate requirement.
In particular the rank-two proposal is insufficient in the bounded
exact-R scalar Long category already constrained by QSO2.

## Define the rank before choosing an operator presentation

Work over an algebraically closed characteristic-zero constant field k.
The central h is a nonzero constant in this generic-fiber discussion.
Let E be an algebraic connection of rank N on an open subset of A1_q.
Let j!*E be its minimal holonomic extension, and define
\[
 r(E)=\operatorname{rank}_{k(a)}
                   \mathfrak F^{-1}(j_{!*}E),\qquad
             a\mapsto h\partial_q,\quad h\partial_a\mapsto-q.
                                                               \tag{1}
\]
Only the generic connection rank of the inverse Fourier transform is
counted. It is distinct from left Weyl-module rank, completed rank,
and the order of a scalar presentation before removal of boundary
constituents.

For the formal connection V_x at a finite point x write
\[
 h_x=\dim_k H^0_{\rm dR}(k((q-x)),V_x),\quad
 \delta_x=\operatorname{irr}(V_x)+N-h_x .
\]
The Fourier rank formula is
\[
 r(E)=\sum_{x\in\mathbb A^1}\delta_x+\rho_\infty,\qquad
 \rho_\infty=
 \operatorname{irr}(V_\infty^{>1})
                         -\operatorname{rank}(V_\infty^{>1})\ge0 .
                                                               \tag{2}
\]
This is the standard Malgrange formula as stated in
[Arinkin, Proposition3.5 and definition(3.1)](https://arxiv.org/html/0808.0742).
It applies to minimal extensions with arbitrary irregular singularities.
The inverse-transform sign changes no rank or irregularity.
The formula and formal slope theory are external mathematical inputs,
not independently replayed theorems in the finite checker.

## A self twist uses up inverse Fourier rank

Let K be the Kummer line with residue 1/e at q=0, with e>=2, and
suppose
\[
                         E\simeq E\otimes K .                 \tag{3}
\]
No irreducibility assumption is needed for this part.
Put i_0=irr(V_0), and let N_reg be the rank of its slope-zero part.
Then
\[
           r(E)\ \ge\ N+i_0-\frac{N_{\rm reg}}e
                    \ \ge\ (1-1/e)N+i_0 .                   \tag{4}
\]
To prove this, the e pairwise nonisomorphic formal Kummer lines
1,K,...,K^(e-1) have the same Hom dimension into V_0, by (3).
That dimension is h_0. Their images give a direct sum of e*h_0
regular rank-one submodules. Distinct simple lines have no common
simple submodule; the corresponding evaluation maps are injective.
Thus e*h_0<=N_reg. Substitute in delta_0=N+i_0-h_0 and use (2).

Also e divides N, by taking determinants in (3). Every possible
self-twist order is covered by the same inequality. Neither local
regularity nor a bound on Stokes matrices is assumed.

For the quadratic twist, N=4 and r=2 is a rigid equality case.
It forces i_0=0, h_0=2, and the complete local module
\[
                       V_0\simeq 1^{\oplus2}\oplus
                                           K_{1/2}^{\oplus2}.
                                                               \tag{5}
\]
There are no nontrivial Jordan extensions: the four independent
line submodules already exhaust its rank. Every other finite
delta_x and rho_infinity must vanish. Thus there are no other
genuine finite singularities, and every slope at infinity is at
most one. This necessary condition does not construct such a
Long-compatible family.

## Reducibility on a cyclic cover forces a self twist

Assume now that E is irreducible and let pi(v)=v^d/d.
If pi^*E is reducible, then (3) holds for some nontrivial Kummer
character of order e dividing d. Here is the descent argument.

Restriction of an irreducible differential Galois representation
to a finite-index normal subgroup is semisimple: an irreducible
submodule exists, its finitely many conjugates span the original
representation, and that span is a sum of simple submodules.
Consequently a reducible restriction has more than scalar
horizontal endomorphisms.

The projection formula for the cyclic cover gives
\[
 \operatorname{End}_{\nabla}(\pi^*E)=
       \bigoplus_{j=0}^{d-1}
          \operatorname{Hom}_{\nabla}(E,E\otimes K_{j/d}).
                                                               \tag{6}
\]
The trivial character contributes only scalars. An additional
endomorphism supplies a nonzero Hom for a nontrivial character;
irreducibility makes that Hom an isomorphism. Equation (4) follows.
In particular
\[
 N>2r(E)\quad\Longrightarrow\quad
 \pi^*E\text{ is irreducible for every cyclic power cover }\pi .
                                                               \tag{7}
\]

If E has at most one finite singularity, N>2r(E) also excludes
reducibility over every finite algebraic function-field extension.
Its Picard--Vessiot component cover is unramified away from that
singularity and infinity; Riemann--Hurwitz makes it cyclic, as in
[QFT1](RATIONAL_LINE_FOURIER_OBSTRUCTION.md). Translate the puncture
and apply (4). Translation does not change inverse-Fourier rank,
since it corresponds to a rank-one exponential twist on the other side.
This extension permits an irregular finite puncture.

With several finite singularities, arbitrary covers with additional
branch points need not give a Kummer self-twist. In that generality
(7) concerns cyclic power covers, not all finite extensions. This
distinction must be retained when testing higher-rank models.

## Boundary extensions cannot lower this cost

Suppose H is any holonomic extension of E, allowing additional
point-supported subquotients on the finite q-line. Its composition
factors comprise those of j!*E and finitely many delta modules.
Fourier transformation is exact and a delta module has a rank-one
exponential Fourier transform. Additivity of generic rank gives
\[
 \operatorname{rank}\mathfrak F^{-1}(H)
                         =r(E)+b,\qquad b\ge0,                \tag{8}
\]
where b is the total point-supported composition length.

Thus adding nonminimal boundary constituents cannot make (4)
cheaper. This does not exclude boundary constructions of a different
classical graph or declare every proposed carrier a minimal extension.
It supplies a necessary numerical test for the connection actually
being extended.

The standard formal vanishing-cycle interpretation of this accounting
is [Arinkin, Lemma3.5 and Proposition5.1](https://arxiv.org/html/0808.0699):
minimal vanishing cycles remove the maximal horizontal submodule,
and boundary delta factors add to Fourier rank.

## Finite pushforward gives a direct test for an extracted branch

Let L be any rank-m connection over k(v), and put E=pi_*L on the
unramified locus of q=v^d/d. Restriction of scalars gives rank N=d*m.
The projection formula gives E tensor K_(1/d) isomorphic to E,
since pi^*K_(1/d) is trivial. Therefore
\[
                    r(\pi_*L)\ge(d-1)m .                    \tag{9}
\]
If L is nonsingular at v=0, E has the d Kummer characters each
with multiplicity m there, h_0=m and irr_0=0. Its contribution
at q=0 is exactly (d-1)m. All other finite singularities and
slopes greater than one at infinity add nonnegative cost.

This result holds even when E is reducible. It does not infer
generic rank from a formal classical rank. The numerical rank
must be computed for the actual generic algebraic connection.

## An explicit family really does split in every polynomial degree

Take g(a)=sum_(j=0)^m g_j*a^j with g_m nonzero and m>=1.
Put q=v^2/2, a=(h/v)*partial_v and
\[
 S_{\rm naive}=g(a)+2v,\qquad
 L_g=k(v)\langle a\rangle/
                  k(v)\langle a\rangle(S_{\rm naive}-s).
                                                               \tag{10}
\]
The source eigenvalue s is included in the constant field.
This is a rank-m connection in v. In the basis
1,a,...,a^(m-1), the ordinary v derivative has matrix
(v/h)*C(v), where C(v)=C0+v*C1 is the companion matrix of
g(a)+2v-s. This matrix is polynomial in v.
Consequently L_g extends nonsingularly across every finite v,
including v=0.

Its pushforward E_g has rank2m, has no finite singularity other
than q=0, and has local module
1^m direct_sum K_(1/2)^m there. At infinity its single slope is
\[
                      1+\frac1{2m}.
\]
Indeed g_m*a^m balances -2v, and a=(h/v)*partial_v.
The exponential phase in v has leading power 2+1/m;
on the q-line this becomes 1+1/(2m).
Equivalently the scalar operator below has a single Newton edge
balancing a^(2m) against q. Its 2m distinct leading exponential
factors form one ramification orbit. Standard formal Newton-polygon
theory therefore also makes E_g formally irreducible at infinity,
and hence generically irreducible. No bounded coefficient search
is used here.

Its infinity contribution in (2) is one. Thus the exact rank is
\[
                   r(E_g)=m+1,\qquad \operatorname{rank}E_g=2m .
                                                               \tag{11}
\]
After the quadratic cover, E_g splits into the two connections
g(a)+2v-s and g(a)-2v-s. This gives genuine algebraic extraction,
rather than merely a formal branch expansion.

For an explicit horizontal projector, use the restriction-of-scalars
matrix
\[
 H=\begin{pmatrix}
       C_0&2qC_1\\ C_1&C_0+hI/(2q)
    \end{pmatrix},\qquad
 V=\begin{pmatrix}0&2qI\\I&0\end{pmatrix}.
\]
One has V^2=2qI and
h*partial_q V+[H,V]=hV/(2q). Since partial_q v=1/v,
\[
                e_\pm=\tfrac12(I\pm V/v)                    \tag{12}
\]
are idempotent horizontal projectors after q=v^2/2.
The image frames (I,+v^-1 I)^t and (I,-v^-1 I)^t
carry the matrices C0+v*C1 and C0-v*C1, respectively.

There is also an explicit cyclic operator. Put b=s-g(a).
Conjugation gives v*b(a)*v^-1=b(a-h/(2q)), so
\[
 P_m=q^m\bigl(b(a-h/(2q))b(a)-8q\bigr).                       \tag{13}
\]
This is polynomial in q,a. It annihilates both pulled-back
branches, has a order2m and q degree m+1. Its inverse Fourier
operator has order m+1. The two branches have disjoint formal
exponential factors at infinity, making its cyclic vector span
their entire induced connection. The rank calculation (11) and
the additive boundary formula (8) show that this cyclic presentation
has no extra finite boundary composition factors.

For Long b=s-4a+3ra^2, the formula becomes
\[
 P_2=q^2 b^2-\frac h2q\,b'b+\frac{9rh^2}{4}b-8q^3 .
                                                               \tag{14}
\]
It is a finite higher-rank inverse-Fourier model with an algebraic
branch projector. Modulo h its equation is q^2*(b^2-8q)=0.
On q!=0 it recovers exactly the classical phase relation
q=(s-4a+3ra^2)^2/8.

The extra q=0 factor is essential provenance: this model has not
been supplied with an h-flat global lattice whose entire special
fiber is only the required Long graph. Its inverse-Fourier generic
rank three must not be called a rational-line repair of QFT1.

## The full right Weyl action remains a separate obstruction

For m=2 the extracted positive operator is precisely
S_naive=4a-3ra^2+2v. It has differential order two, unit leading
coefficient over the generic r field, and Long's classical S.
A full Long realization would also need a regular finite operator
T with [S,T]=h and classical T=(ra^3-a^2-va)/2, of order three.
[QSO2](QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md) excludes this
for every finite partner order, on every finite algebraic coefficient
curve. In particular no extra partner terms, finite pole enlargement
or algebraic cover supply that missing action while S_naive is fixed.

More generally, QSO2 forces every bounded full scalar Long lift with
exact undeformed R into generic fiber order m>=3 for S.
For its spectral connection L=(S-s), the actual generic v rank is m.
Equation (9) therefore forces inverse-Fourier rank at least three.
A proposed rank-two model in this category fails before further
quantum corrections are attempted.

This last application assumes the bounded exact-R scalar reduction
used in QSO2. It does not normalize every deformed first image, supply
a scalar frame for every projective carrier, or exclude a nonpositional
full carrier.

The linear case is a positive full-pair control. For g=lambda*a+mu,
\[
 S=\lambda a+\mu+2v,\qquad
 T=-\lambda a^2/8-va/2-h/(4v)
\]
satisfy [S,T]=h exactly. This control and the QFT1 rank-one/rank-two
split keep the numerical inequality separate from the additional
Long mate obstruction.

The next useful carrier construction must satisfy the rank budget,
the full right Weyl action, global admission and a finite-left lattice
together. Branch splitting alone is now explicitly achievable and
explicitly insufficient. The finite-left Long carrier in
[QMC1](FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md) and the bounded
scalar target remain open.

[Exact controls](check_fourier_rank_budget.py) verify the projector
and its branch frames, polynomial Fourier elimination in declared
degrees, the quadratic indicial roots, the classical boundary factor,
the naive mate defect and the linear exact pair. The universal rank
inequality, cyclic descent, all-degree slope and rank calculations
use the written proofs and the cited standard Fourier-rank and
formal-connection theorems. No independent theorem replay, formal
verification, external review, novelty priority, actual finite-left
Long carrier, bounded Long lift, scalar or matrix counterexample,
all-representative nonlifting or DC2 solution is claimed.
