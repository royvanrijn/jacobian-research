# The clean quadratic-gauge decorated receiver

This note constructs the enhanced receiver suggested by generic
Tschirnhaus non-descent.  Its main purpose is to keep three different
objects separate:

1. an unmarked stable ambient Keller map;
2. a marked complete fibre of such a map; and
3. the abstract finite-etale algebra carried by that fibre.

Only the second object has a natural morphism to `BS_N`.  The resulting
diagram explains both the `N-4` stable-map count and the larger marked
incidence space without assigning a nonexistent map to the unmarked
quotient.

Work over a characteristic-zero field and fix `N>=5`.  Restrict throughout
to the clean quadratic-gauge chart: all compiler seeds are nonzero, the
first target coordinate is nonzero, and the inverse polynomial is
squarefree.

## 1. The marked incidence chart

Let

\[
 K_N^\times=\mathbb G_m^{N-3}
\]

have seed coordinates `(u_4,\ldots,u_N)`.  For target coordinates
`(\pi,b,c)` put

\[
 E_{u,\pi,b,c}(S)
 =
 S+bS^2+\pi S^3+
 \sum_{j=4}^N u_j\pi^jS^j-\frac c2.                  \tag{1.1}
\]

Let `I_N^\circ` be the discriminant complement in

\[
 K_N^\times\times\mathbb G_{m,\pi}\times\mathbb A^2_{b,c}.
                                                               \tag{1.2}
\]

It has dimension `N`.  The universal relative Keller map has over this
chart the marked complete fibre

\[
 X_{u,\pi,b,c}
 =
 \operatorname{Spec}
 \frac{{\cal O}_{I_N^\circ}[S]}{(E_{u,\pi,b,c})},      \tag{1.3}
\]

which is finite etale of rank `N`.

## 2. Residual equivalence and descent of the fibre

The residual coefficient-torus equivalence is the `G_m` action

\[
\begin{aligned}
 u_j&\longmapsto\alpha^{j+1}u_j,\\
 \pi&\longmapsto\alpha^{-2}\pi,\\
 b&\longmapsto\alpha^{-1}b,\\
 c&\longmapsto\alpha c.
                                                               \tag{2.1}
\end{aligned}
\]

The inverse polynomial transforms by the exact identity

\[
 E_{\alpha\cdot(u,\pi,b,c)}(S)
 =
 \alpha E_{u,\pi,b,c}(S/\alpha).                     \tag{2.2}
\]

Thus `S -> alpha*S` identifies the two finite-etale fibres.  Consequently
the quotient stack

\[
\boxed{
 {\mathscr D}_N^{\rm quad,\circ}
 =
 [I_N^\circ/\mathbb G_m]
}                                                       \tag{2.3}
\]

has a canonical fibre-forgetting morphism

\[
\boxed{
 \rho_N:
 {\mathscr D}_N^{\rm quad,\circ}\longrightarrow BS_N.
}                                                       \tag{2.4}
\]

This is the clean quadratic-gauge decorated receiver.  An object remembers
a marked Keller fibre modulo the explicit ambient equivalence, so both the
finite-etale algebra and its ambient quadratic-gauge decoration remain
available.

The action is actually globally split and free, not merely faithful.  Put

\[
 \lambda=\frac{u_5}{u_4}.                              \tag{2.5}
\]

Then `lambda -> alpha*lambda`.  Let

\[
 Q_N^\circ=I_N^\circ\cap\{\lambda=1\}
           =I_N^\circ\cap\{u_5=u_4\}.                  \tag{2.6}
\]

The action map gives an equivariant isomorphism

\[
 \mathbb G_m\times Q_N^\circ\longrightarrow I_N^\circ,
 \qquad(\alpha,q)\longmapsto\alpha\cdot q,             \tag{2.7}
\]

whose inverse sends `x` to
`(lambda(x),lambda(x)^(-1)*x)`.  In coordinates, the invariant target
variables on the slice are

\[
 \widehat\pi=\lambda^2\pi,\qquad
 \widehat b=\lambda b,\qquad
 \widehat c=c/\lambda,                                 \tag{2.8}
\]

while invariant seed characters can be taken to be

\[
 q_4=u_4/\lambda^5,\qquad
 q_j=u_j/\lambda^{j+1}\quad(6\le j\le N).
\]

Thus the quotient stack is represented by the explicit
discriminant-open scheme

\[
 \boxed{{\mathscr D}_N^{\rm quad,\circ}\simeq Q_N^\circ.} \tag{2.9}
\]

In particular, it has no residual stack inertia from coefficient scaling,
the quotient map is a trivial `G_m`-torsor, and `rho_N` is a representable
morphism to `BS_N`.  Also

\[
 \dim{\mathscr D}_N^{\rm quad,\circ}=N-1.              \tag{2.10}
\]

This slice is defined over the ground field; unlike the normalization
`u_4=1`, it requires no root extraction.

## 3. The unmarked quotient is a different object

Forgetting the target gives the equivariant projection

\[
 I_N^\circ\longrightarrow K_N^\times.
\]

Its quotient is the stable ambient-map morphism

\[
 p_N:
 {\mathscr D}_N^{\rm quad,\circ}
 \longrightarrow
 {\mathcal M}_N^{\rm quad}
 =
 K_N^\times/\mathbb G_m
 \simeq\mathbb G_m^{N-4}.                              \tag{3.1}
\]

The geometric fibres of `p_N` are nonempty open subsets of a
three-dimensional target space.  Thus

\[
 \dim{\mathcal M}_N^{\rm quad}=N-4,\qquad
 \dim{\mathscr D}_N^{\rm quad,\circ}
 =(N-4)+3=N-1.                                        \tag{3.2}
\]

There is no corresponding morphism

\[
 {\mathcal M}_N^{\rm quad}\longrightarrow BS_N.        \tag{3.3}
\]

An unmarked polynomial map has many fibres, and no target has been selected
on the left side of (3.3).  The correct structure is the span

\[
\boxed{
\begin{array}{ccc}
 &{\mathscr D}_N^{\rm quad,\circ}&\\
 {}^{\rho_N}\swarrow&&\searrow^{p_N}\\
 BS_N&&{\mathcal M}_N^{\rm quad}.
\end{array}
}                                                       \tag{3.4}
\]

The combined map

\[
 (\rho_N,p_N):
 {\mathscr D}_N^{\rm quad,\circ}
 \longrightarrow
 BS_N\times{\mathcal M}_N^{\rm quad}                   \tag{3.5}
\]

is the precise “finite-etale algebra plus intrinsic quadratic-gauge
decoration” receiver on this chart.

## 4. What the dimension counts mean

After base change to an algebraic closure, every rank-`N` finite-etale
algebra is split.  Therefore the geometric fibre of `rho_N` has dimension
`N-1`.  Its image under `p_N` is the full `N-4` dimensional stable-map
quotient: every clean seed admits a nonempty open set of squarefree
targets.

These are compatible statements:

\[
\begin{array}{c|c}
\text{quantity}&\text{dimension}\\ \hline
\text{marked decorated receiver over a geometric algebra}&N-1\\
\text{target choices over a fixed ambient class}&3\\
\text{ambient classes possessing such a fibre}&N-4.
\end{array}                                             \tag{4.1}
\]

Thus `N-4` counts ambient stable-map realizations, not objects in the
literal fibre of (2.4).  It is also not a global lower bound for `ktdim` or
`kdeg`.

## 5. Tschirnhaus non-descent in the receiver

A primitive presentation determines a point of `I_N^\circ`, hence a point
of the decorated receiver.  Two presentations of the same abstract
finite-etale algebra have the same image under `rho_N`, but generically
have different images under `p_N`, detected by

\[
 \Phi_N=(I_5,J_6,\ldots,J_N).                          \tag{5.1}
\]

Therefore the compiler atlas does not factor through `BS_N`:

\[
\boxed{
 \text{the fibre algebra descends through }\rho_N,
 \quad
 \text{the ambient decoration }p_N\text{ does not.}
}                                                       \tag{5.2}
\]

The rank-five transition-locus calculation sharpens this further.  Equal
ambient decoration is a four-dimensional hypersurface in the ordered-root
chart, whereas the canonical equivalence carries the distinguished marked
fibre only on the one-dimensional root-scaling locus.  Any larger marked
identification must come from the target orbit of the stable
self-equivalence group of one fixed Keller map.

## 6. Physical inertia modulo vertical gauge

The
[stable intruder descent criterion](STABLE_INTRUDER_DESCENT_CRITERION.md)
and its all-rank quadratic-gauge application show more than freeness of
the explicit coefficient action.  For every geometric point of the clean
receiver, every stable marked self-equivalence is the identity on the
physical source and target coordinates.  If `Aut_vert` is the kernel of
restriction to those coordinates, then

\[
 \operatorname{Aut}_{\rm st}^{\rm marked}/
 \operatorname{Aut}_{\rm vert}=1.                     \tag{6.1}
\]

Thus the clean receiver is an explicit scheme with trivial physical
inertia.  In any future stable-map groupoid containing this chart, the only
pointwise inertia left on the chart is vertical stabilization gauge.

This does not classify that vertical group, nor does it prove that a
global quotient by vertical automorphisms is algebraic, separated, or of
finite type.

## 7. Consequences for the attack chain

The receiver changes the next attacks as follows.

1. Do not seek a section from `BS_N` directly to the unmarked quotient
   `M_N^{quad}`; that target has forgotten the datum needed to define the
   fibre.
2. Study sections, correspondences, or torsors for `rho_N`, while measuring
   their failure to have constant `p_N`-coordinate.
3. Do not spend the next attack on higher-degree fixed-map stabilizer
   calculations: the all-rank intruder theorem already makes the physical
   orbit a point in arbitrary degree and stabilization dimension.
4. For an alternative compiler, test whether its analogue of `p_N` is
   constant on Tschirnhaus arrows.  That is the actual descent criterion.
5. Use the represented receiver as a test locus for the
   [degreewise stable Keller-moduli problem](../extended-geometry/STABLE_KELLER_MODULI_PROBLEM.md):
   construct its morphism to the stable groupoid, determine whether its
   image is locally closed, and compute its normal deformation theory.
   This still requires defining and quotienting the vertical stabilization
   gauge; the pointwise physical-inertia result does not construct the
   global groupoid.

This constructs and represents the clean quadratic-gauge receiver.  It
does not construct a finite-type global moduli stack of all stable
polynomial maps.

## Exact identities

The action (2.1), identity (2.2), global slice (2.7), and stable quotient
coordinates are checked by

```bash
.venv/bin/python scripts/verify_quadratic_gauge_stable_moduli.py
.venv/bin/python scripts/verify_universal_relative_keller_map.py
```

The rank-five marked specialization is checked by

```bash
.venv/bin/python scripts/verify_rank_five_tschirnhaus_transition_locus.py
```
