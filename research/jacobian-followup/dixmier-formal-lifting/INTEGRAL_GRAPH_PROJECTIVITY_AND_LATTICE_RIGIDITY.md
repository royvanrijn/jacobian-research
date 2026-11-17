# Integral projectivity and rigidity of finite quantum graph carriers

15 November 2026. **QIR1.** A finite, parameter-flat full quantum graph
carrier is projective over the parameter Weyl algebra itself. It therefore
has a finite idempotent certificate with all right Weyl relations and
special idempotent e11. Generic projectivity need not be added as an
integral hypothesis.

For a fixed rational full right action, there is at most one full finite
left lattice over the generic Weyl algebra. Two parameter lattices with
the same full symplectic classical graph differ by a central parameter
power. This holds for every ramification order: rank-one special graph
modules automatically admit finite scalar corrections up to quantum
order. The statements hold in every Weyl rank and impose no spatial
degree bound.

These results control lattices for a fixed rational action. They neither
construct that action for Long nor compare arbitrary formal gauges.

## Parameter algebra and full graph condition

Let K have characteristic zero, R=K[[t]], L=K((t)), and
a(t)=t^e u(t), where e>=1 and u is a unit. Put
\[
 A_R=R\langle x_1,\ldots,x_n,p_1,\ldots,p_n\rangle/
 ([p_i,x_j]-a(t)\delta_{ij},[x_i,x_j],[p_i,p_j]),\qquad
 A=A_R[t^{-1}],\quad D=\operatorname{Frac}A,\quad
 S=A_R/tA_R=K[x,p].                                      \tag{1}
\]
Only finite spatial PBW sums belong to A_R. The scalar coefficients can
be arbitrary elements of R. After rescaling momenta, A is a Weyl algebra
over L.

A full graph carrier is an R-central A_R-bimodule M, finite on the left,
R-flat, with an identified special bimodule M/tM=S_psi. Here the left
module is S and each right generator acts by multiplication by the
corresponding coordinate of a full polynomial symplectic map psi.
All right Weyl relations must hold on M. These are the
[QMC1 hypotheses](FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md).

## Projectivity across the parameter fiber

First use the following general homological lemma. Suppose B is a
left/right Noetherian R-flat algebra of finite global dimension, and N
is a finite left B-module. If
\[
 N\text{ is R-flat},\qquad N/tN\text{ is projective over }B/tB,\qquad
 N[t^{-1}]\text{ is projective over }B[t^{-1}],
 \quad\text{then }N\text{ is projective over }B.             \tag{2}
\]

Take a finite resolution P_bullet of N by finite projectives and set
C_bullet=Hom_B(P_bullet,B). This is a bounded complex of finite right
projectives, and its cohomology E^j=Ext_B^j(N,B) is finite on the right.
Reduction of P_bullet remains a resolution: projectives are R-flat and
N is R-flat. The exact sequence of complexes given by multiplication by
t gives, for j>0,
\[
 E^j/tE^j\ \lhook\joinrel\longrightarrow\
 \operatorname{Ext}_{B/tB}^j(N/tN,B/tB)=0,
 \qquad E^j[t^{-1}]=0.                                   \tag{3}
\]
The second equality follows by localization of the finite projective
resolution. Each E^j is finite, so its t-torsion is killed by one power
t^b. The first equality says E^j=tE^j. Iterating b times gives E^j=0.
This uses bounded central torsion, rather than Nakayama with t in a
Jacobson radical.

A finite projective resolution with all positive Ext into B zero
shortens to length zero. Indeed the last dual differential is
surjective, splits because its target is projective, and dualizing the
splitting removes the last projective term. Induct on the length. This
proves (2).

The algebra A_R meets the ring hypotheses. It is an iterated polynomial
and differential Ore extension of the regular DVR R. It is Noetherian,
R-flat and has global dimension at most 2n+1. The global-dimension bound
follows successively from the induced two-term resolution for an Ore
extension, which raises a finite global-dimension bound by at most one.
A finite module consequently has a finite resolution by finite
projectives.

For a full graph carrier, M/tM=S is free. The generic M_L is projective
by [Kanel-Belov--Fedorov, Lemma3.1](https://arxiv.org/html/2610.33919v1):
a bimodule over a characteristic-zero Weyl algebra which is finite on
the left is projective there. As in QMC1, its rank is one. Thus (2)
applies. The generic projectivity theorem is an external input; the
central-fiber gluing argument above is self-contained.

## A necessary and sufficient finite corner certificate

Integral projectivity supplies a finite s and an idempotent E with
M=A_R^s E in row convention. The full right action is given by matrices
r_i=E r_i E. Conversely such a corner is finite and R-flat on the left.

The special module is identified with S, so there are a row v and a
column w over S with vw=1 and wv=E_0. This is part of the given special
frame; no new projective-freeness theorem is required. Add one zero
coordinate to E and lift v,w arbitrarily to A_R. Define
\[
 U(\widetilde w)=
 \begin{pmatrix}I_s&\widetilde w\\0&1\end{pmatrix},\qquad
 V(-\widetilde v)=
 \begin{pmatrix}I_s&0\\-\widetilde v&1\end{pmatrix},\qquad
 G=U(\widetilde w)V(-\widetilde v)U(\widetilde w).           \tag{4}
\]
All three factors and their inverses are polynomial matrices, even if
the arbitrary lifts do not satisfy tilde(v)tilde(w)=1. At t=0,
\[
 G_0=\begin{pmatrix}I_s-E_0&w\\-v&0\end{pmatrix},\qquad
 G_0\operatorname{diag}(E_0,0)G_0^{-1}
       =\operatorname{diag}(0_s,1).                        \tag{5}
\]
A final constant permutation moves the last coordinate to the first.
Conjugate E and every r_i by the same actual invertible matrix. This
produces, in some finite matrix size,
\[
 E^2=E,\quad E_0=e_{11},\quad r_i=E r_i E,\quad
 [r_i,r_j]=a(t)J_{ij}E,\quad (r_i)_0=\psi_i e_{11}.         \tag{6}
\]
Thus existence of a full finite graph carrier is equivalent to existence
of a finite certificate (6). Its size and degrees are not bounded here.
The stabilization adds one coordinate to an already obtained projective
presentation. It does not assert that the carrier is a free scalar
left module.

## Generic full lattices are unique for a fixed right action

Let V be a finite-dimensional left D-space with a commuting,
L-central full right A-action. A full finite lattice means an
A-subbimodule P in V, finite on the left, with D P=V.

A finite left A-bimodule Q which is torsion on the left must be zero.
If Q were nonzero, a good left filtration would give
\[
 \operatorname{GKdim}_A Q\le 2n-1,\qquad
 \operatorname{GKdim}_{A\otimes_L A^{op}}Q
      =\operatorname{GKdim}_A Q\ge 2n,                     \tag{7}
\]
a contradiction. The equality follows because finitely many right
generator actions raise a good left filtration by bounded amounts.
The final inequality is Bernstein's inequality for A_(2n).
For the first inequality, each of finitely many left generators of a
torsion module is a quotient of A by a nonzero left ideal, whose symbol
support has dimension at most 2n-1.

If P,Q are two full finite lattices in V, then P intersect Q is finite
on the left, right-stable and full. Fullness follows by expressing a
D-basis chosen in Q in terms of generators of P and clearing finitely
many left Ore denominators. The quotients by the intersection are
finite left torsion bimodules and hence zero. Therefore
\[
                         P=Q.                            \tag{8}
\]
This argument works for every finite D-rank. It requires the full right
Weyl algebra, not one preserving operator or only a commuting tuple.

## Rank-one graph congruence upgrades to quantum order

Let M be a full graph carrier and let P_i be any polynomial left lifts
of psi_i. Write R_i for the corresponding right generator action.
The identified graph gives (L_(P_i)-R_i)M subset tM.

Suppose, for some 1<=j<e, that this difference lies in t^j M.
The operator
\[
       T_{i,j}=(L_{P_i}-R_i)/t^j:M\longrightarrow M        \tag{9}
\]
reduces to an S-linear endomorphism of M/tM=S. Indeed its commutator
with any left coefficient b is left multiplication by [P_i,b]/t^j,
which vanishes modulo t since every algebra commutator is divisible
by a(t). Thus its reduction is multiplication by c_(i,j) in S.
Choose a finite polynomial left lift C_(i,j) and replace
P_i by P_i-t^j C_(i,j). The new difference lies in t^(j+1)M.
After at most e-1 such steps,
\[
             (L_{P_i}-R_i)M\subset a(t)M.                 \tag{10}
\]
This automatic upgrade uses the free rank-one special left module.
For a higher-rank special module the sub-quantum corrections can be
matrices, and the scalar conclusion is not automatic.

## Integral graph lattices differ only by a central power

Now fix the rational right action on V, and suppose M,N are full
finite A_R-lattices with identified special graph S_psi. By (8),
M_L=N_L. Central localization and finite generation make them
commensurable: t^b M subset N for some b in Z.

Choose the least b with this property. Such a least integer exists
because N is projective over A_R, hence t-separated as a direct summand
of a finite free module. If arbitrarily negative b worked, each element
of M would lie in every t^j N and would be zero. The induced map
t^b M/t^(b+1)M -> N/tN is nonzero, since otherwise b-1 would work.
Under the special graph frames it is multiplication by a nonzero
polynomial f in S.

Perform the preceding sub-quantum corrections simultaneously on
t^b M and N. At each j<e, their reduced operators are multiplication
by c_M and c_N. The inclusion intertwines them, giving
f c_M=c_N f. Since S is a domain and f is nonzero, c_M=c_N.
The same polynomial correction therefore works for both lattices.

The divided operators now preserve both modules:
\[
 \nabla_i=(L_{P_i}-R_i)/a(t),\qquad
 \overline{\nabla_i^M}=X_i+\alpha_i,\quad
 \overline{\nabla_i^N}=X_i+\beta_i,\quad
 X_i=\{\psi_i,-\},\quad \alpha_i,\beta_i\in S.              \tag{11}
\]
They need not be flat. Intertwining the inclusion gives
\[
            X_i(f)=(\alpha_i-\beta_i)f.                   \tag{12}
\]
A full polynomial symplectic frame spans all coordinate derivatives
with polynomial coefficients. Thus f divides every partial derivative
of f. In characteristic zero this forces f to be a nonzero constant:
a nonzero derivative of a nonconstant polynomial has smaller degree.

The inclusion is consequently an isomorphism modulo t. Its cokernel C
satisfies C=tC. Its generic fiber is zero by (8); finite generation and
centrality of t give t^k C=0 for some k. Hence C=0, and
\[
                       N=t^b M.                          \tag{13}
\]
All ramification orders e are covered by the finite simultaneous
correction argument. Mere agreement of support, special torsion, or a
higher-rank special module would not justify this argument.

## Polynomial right actions and the remaining scalar frame

Over the generic field, suppose a full finite lattice P exists for
a rational right action on D^r. If some D-linear change of frame puts
all right generator matrices in Mat_r(A), the standard lattice A^r
is stable in that frame. Uniqueness (8) then forces it to be P.
Conversely a free left frame of P makes every right action polynomial.
Consequently
\[
 P\text{ is free of rank }r
 \quad\Longleftrightarrow\quad
 \text{the full right action is polynomial in some D-frame}. \tag{14}
\]
In particular no proper nonzero left ideal of A is stable under a full
polynomial Weyl tuple. The quotient would be the forbidden torsion
bimodule of (7). This excludes all polynomial tuples in such an idealizer,
in every degree; arbitrary rational idealizer tuples remain separate.

The [QPF1 family](PROJECTIVE_FRAME_BOUNDARY.md) illustrates why formal
freeness does not settle existence. Its inclusion into A_R is an
isomorphism in every parameter jet and its generic left module is
nonfree. It cannot carry the regular full right action on that inclusion.
If it carries some other full rational action, changing only its finite
lattice cannot free it, by (8) and (14).

For Long, (6) is now a necessary and sufficient finite certificate,
rather than only a sufficient construction. The direct carrier remains
unconstructed. A genuinely new rational action or a proof that every
such certificate is impossible is still required. A formal gauge need
not be rational: the [QNI1 controls](IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md)
already prevent a universal promotion from formal equivalence to a
bounded algebraic comparison in higher rank.

The [QHI1 symbol criterion](HAMILTONIAN_IDEALIZER_CRITERION.md) remains a
possible scalar-frame gate after a carrier is found. No full generic
symbol frame is inferred from its special classical graph. The
[QGC1 closure obstruction](RAMIFICATION_AND_FINITE_CLOSURE_CARRIERS.md)
concerns extra residue covers; the rank-one scalar correction here does
not silently apply to its degree-two or higher special modules.

## Exact controls and research boundary

Run check_integral_graph_rigidity.py. It checks an explicit elementary
stabilization whose arbitrary lifts have product 1+t*p, all six full
corner Weyl relations with an independent quantum parameter, and
ramified scalar-correction controls. It also checks the Hamiltonian
frame inversion and a proper lattice preserved by a single right
coordinate which fails the full Weyl action.

The projectivity, uniqueness, scalar correction and homothety theorems
are written proofs for all ranks, degrees and parameter orders. The
listed controls do not constitute independent theorem replay, formal
verification or external review. No actual Long carrier, bounded scalar
lift, all-representative nonlifting theorem, DC2/JC2 solution or novelty
priority is claimed.
