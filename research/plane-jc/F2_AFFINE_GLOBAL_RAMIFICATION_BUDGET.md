# Global affine-ramification budget and the cubic E8 saturation row

> **Status.**  Exact conditional global identity and exact degree-six
> permutation/local-model theorem.  For affine nonproperness components
> whose normalizations have genera `g_j`, puncture counts `s_j`, and total
> moved-sheet degrees `A_j`, subtracting the complete logarithmic determinant
> square and all immersed conormal kernels leaves
> `d-1+sum_j A_j(2g_j-2+s_j)`.  For the rational one-puncture curves forced
> in F2 this is `d-1-sum_j A_j`.  Equivalently, with `r` components and
> meridian fixed-sheet counts `u_j`, it is
> `sum_j u_j-(r-1)d-1`.  This is a global equality, not an automatic
> module-by-module positivity statement: contracted determinant packets
> contribute their kernel degrees, and generically noncyclic packets require
> a separate filtration.  Orevkov's Euler identity plus Riemann--Hurwitz now
> proves unconditionally that the **total numerical remainder** is
> nonnegative and equals the sum of residue-normalized local multiplicity
> excesses; see
> [`OREVKOV_RESIDUE_DEGREE_BUDGET.md`](OREVKOV_RESIDUE_DEGREE_BUDGET.md).
>
> At geometric degree six, the complete `(3,5)` cusp-group enumeration with
> meridian type `3+1+1+1` consists of one conjugacy class.  Its image is
> `A_6`, its peripheral row is `(e,f)=(3,1)`, and its point budget `u-1=2`
> exactly equals the E8 cusp lower bound `2f=2`.  Its passport is precisely
> the certified F2 terminal passport
> `(5,1)|(3,3)|(3,1,1,1)`.  Thus the global `ch_2` route does not by itself
> exclude all higher inertia: it isolates one degree-six equality row.
> An exact SNC local model realizes the cubic transverse contact but creates
> a contracted component with generic Smith form `diag(t,t)`.  This locates
> the next gap at a positive-dimensional `Fitt_1` packet, outside the
> isolated-defect positivity theorem.

The algebra, complete `S_6` enumeration, and local matrix are replayed by
[`verify_f2_affine_global_ramification_budget.py`](../scripts/verify_f2_affine_global_ramification_budget.py).

## 1. The multi-component cancellation

Let

\[
 f:(X,D_X)\longrightarrow(Y,D_Y),\qquad
 L_X=K_X+D_X,\quad L_Y=K_Y+D_Y
\]

be a completed logarithmic model of geometric degree `d`.  Put

\[
 D_{\log}=L_X-f^*L_Y.                              \tag{1.1}
\]

Let `C_1,...,C_r` be the affine nonproperness components.  Write

\[
 \beta_j=2g_j-2+s_j                              \tag{1.2}
\]

for the logarithmic canonical degree of the smooth normalization of `C_j`.
If the rows over `C_j` have transverse and residue degrees `(e_{ji},f_{ji})`,
put

\[
 A_j=\sum_i e_{ji}f_{ji}.                       \tag{1.3}
\]

Finally set `c_j=L_Y.C_j`.  Projection gives

\[
 f^*L_Y.D_{\log}=\sum_j A_jc_j,                 \tag{1.4}
\]

because determinant components contracted to target points have zero
pairing with `f^*L_Y`.  Hence

\[
 D_{\log}^2=L_X^2-dL_Y^2-2\sum_jA_jc_j.         \tag{1.5}
\]

The immersed logarithmic conormal formula gives the affine kernel degree

\[
 \deg K_{\rm aff}=\sum_j A_j(c_j-\beta_j).       \tag{1.6}
\]

Subtract (1.5)/2 and (1.6) from the global budget

\[
 B_f=\frac12\bigl(L_X^2-dL_Y^2+2(d-1)\bigr).
\]

Every log square and contact cancels:

\[
\boxed{
 B_f-\left(\deg K_{\rm aff}+\frac12D_{\log}^2\right)
 =d-1+\sum_jA_j\beta_j.}                        \tag{1.7}
\]

This is the general complete-chain identity.  If a cyclic contracted packet
has total kernel degree `kappa_0` and the remaining quotient has finite
length `P`, then the exact filtration reads

\[
 \boxed{P+\kappa_0=d-1+\sum_jA_j\beta_j.}        \tag{1.8}
\]

Equation (1.8) is the safe form: `kappa_0` need not vanish on a new
contracted packet, and a generically noncyclic component cannot be inserted
as a finite point correction.

## 2. F2 rational one-puncture specialization

Every F2 affine target component has normalization `A^1`, hence
`(g_j,s_j)=(0,1)` and `beta_j=-1`.  Therefore

\[
 \boxed{P+\kappa_0=d-1-\sum_jA_j.}               \tag{2.1}
\]

If `u_j` is the number of sheets fixed by a geometric meridian of `C_j`,
the generic fiber identity is

\[
 A_j=d-u_j.                                     \tag{2.2}
\]

For `r` affine components, (2.1) becomes

\[
 \boxed{
 P+\kappa_0=\sum_{j=1}^r u_j-(r-1)d-1.}         \tag{2.3}
\]

Independently of any logarithmic cyclicity hypothesis, Orevkov's global
Euler identity and Riemann--Hurwitz on the one-puncture residue covers give
the useful necessary inequality

\[
 \sum_jA_j\le d-1.                              \tag{2.4}
\]

More precisely, the difference in (2.4) is the sum of the nonnegative local
excesses `mu_x-e*q_x`.  What remains conditional is the refinement of that
number into individual logarithmic kernel and quotient modules.  Higher cusp
inertia can create a contracted divisor on which the cokernel is already
noncyclic at the generic point; Section 4 gives an exact example.

## 3. The unique degree-six cubic-inertia action

For the E8 cusp group

\[
 G=\langle a,b\mid a^3=b^5\rangle,
 \qquad m=a^{-1}b^2,                            \tag{3.1}
\]

enumerate all pairs `(A,B)` in `S_6^2`, impose (3.1), transitivity, and

\[
 \operatorname{type}(M)=3+1+1+1.               \tag{3.2}
\]

There are exactly `720` labeled solutions.  They form one simultaneous
conjugacy orbit.  In every solution

\[
 A^3=B^5=1,\quad
 \operatorname{type}(A)=3+3,\quad
 \operatorname{type}(B)=5+1,                   \tag{3.3}
\]

and the generated group has order `360`.  Thus its image is `A_6` in the
natural degree-six action.  One representative is

\[
 \boxed{
 A=(1\ 2\ 3)(4\ 5\ 6),\quad
 B=(2\ 4\ 3\ 5\ 6),\quad
 M=(1\ 3\ 5).}                                 \tag{3.4}
\]

The preferred longitude is `ell=z m^(-15)`, with `z=a^3=b^5`.  Here
`z=1` and `M^3=1`, so `ell=1`.  The ramified peripheral orbit is therefore
one three-cycle:

\[
 \boxed{(e,f)=(3,1),\qquad A=3,qquad u=3.}      \tag{3.5}
\]

For one rational one-puncture component, (2.1) gives point budget

\[
 u-1=2.                                         \tag{3.6}
\]

The multiplicity-three E8 cusp lower is `2f=2`, so this row is exact
equality, not a deficit.  Moreover (3.3)--(3.4) give the passport

\[
 \boxed{(5,1)\mid(3,3)\mid(3,1,1,1),}           \tag{3.7}
\]

which is exactly the passport of the already certified degree-six terminal
residue cover.  The terminal and affine-cusp monodromies are therefore
compatible at the level of group, action, and passport.  This coincidence
does not construct their global gluing, but it rules out a contradiction
based only on those three invariants.

## 4. Why isolated positivity stops at the equality row

The cubic transverse contact has the following exact completed local model:

\[
 x=t^3+rt,
 \qquad
 y=t^5+\frac53rt^3+\frac59r^2t.                 \tag{4.1}
\]

Direct calculation gives

\[
 J_{(r,t)}(x,y)=-\frac59r^2t                   \tag{4.2}
\]

and

\[
 y^3-x^5=\frac{r^3t^3}{729}
 \left(125r^3+396r^2t^2+405rt^4+135t^6\right). \tag{4.3}
\]

Thus `E=(r=0)` has transverse cusp multiplicity three.  Keller boundary
purity also forces `T=(t=0)` into the source boundary.  It is contracted to
the affine cusp point.

At the SNC node use source logarithmic derivations `(r d/dr,t d/dt)`.  The
logarithmic matrix is

\[
 \Theta=
 t\begin{pmatrix}
 r&\frac59r(2r+3t^2)\\
 r+3t^2&\frac59(r^2+9rt^2+9t^4)
 \end{pmatrix},                                \tag{4.4}
\]

with

\[
 \det\Theta=-\frac59r^3t^2.                    \tag{4.5}
\]

After dividing by `t`, the determinant is `-5r^3/9`, a unit at the generic
point of `T`.  Hence

\[
 \boxed{
 \operatorname{Smith}_{\eta_T}(\Theta)=\operatorname{diag}(t,t).} \tag{4.6}
\]

The cokernel is generically split, not cyclic, along the whole contracted
component.  Its `Fitt_1` has a height-one factor `t`; it is not an isolated
point defect.  The cyclic-submodule positivity theorem therefore does not
apply to this packet, exactly as its claim boundary predicts.

The example does not prove that every cubic-inertia completion has local
form (4.1), nor that the split packet can be glued to the terminal cover.
It proves something more limited but decisive for the present route: neither
the global scalar determinant budget, isolated `Fitt_1` positivity, nor the
`A_6` passport alone excludes the first higher-inertia row.

## 5. Reduced completion target

For the degree-six floor, the former four-way frontier reduces to one exact
packet:

1. one E8 affine component;
2. cubic meridian inertia `(3,1,1,1)`;
3. natural `A_6` monodromy, matching the terminal cover;
4. equality in the global point budget; and
5. a possible generically split contracted component over the cusp point.

An exclusion now needs one genuinely new statement: classify the
positive-dimensional logarithmic cokernel on every contracted cusp packet
and prove that none can glue to the certified terminal `A_6` cover on the
compiled F2 boundary graph.  The calculation is no longer a 927-variable
Laurent descent and no longer an unspecified sign problem.  It is a
rank-two gluing problem with fixed degree, group, passport, transverse
index, and generic Smith type.

The first half of that statement is now proved by the
[`contracted-divisor Smith classification`](LOG_CONTRACTED_DIVISOR_SMITH_CLASSIFICATION.md).
The cubic packet factors as a universal rank-two quotient on the contracted
component plus the saturated normal form `[[r,0],[t^2,r^2]]`.  Its isolated
cyclic-submodule quotient has length four, and the remaining global gate is
the single incidence equation `P_other=2I-v-3n>=0`, with
`I=(D_log-2T).T`, `v=valency(T)`, and `T^2=-n`.

<!-- status-consumer: LCDSC1 07dcd994b4faf092 -->

Higher geometric degrees, additional affine components, non-E8 `k=1`
strata, and `k=2,...,24` remain separate claim boundaries.  Consequently
this note does not exclude all `(75,125)` or prove `JC(2)`.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_global_ramification_budget.py
```

The command checks the abstract genus/puncture cancellation over a grid of
exact intersection inputs, all `S_6` cusp-group pairs satisfying (3.1)--(3.2),
their single conjugacy orbit and order-`360` image, the peripheral longitude,
the saturated budget, the terminal-passport match, and every polynomial and
matrix identity in Section 4.
