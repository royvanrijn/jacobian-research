# Kummer-orbit transfer for the F2 `(75,125)` common edge

## Result and claim boundary

The earlier F2 boundary handoff treated the roots of

\[
C_0(u)=(u-1)^5(u^5-1)^2R(u^5)
\]

as unrelated centers on the selected `X^5=x` cover. That loses a rigid piece
of source structure. Every Laurent band comes from `k[X^5,y]`, so its
coefficient has one Kummer character and its exact vanishing order transfers
unchanged around every nonzero `mu_5` orbit.

For the F2 row this gives three exact improvements.

1. The five roots over the selected factor `u^5=1` have identical complete
   Newton data in their natural Puiseux charts. In particular, the final
   principal block and its bracket `X^4` transfer exactly; the four
   multiplicity-two contacts in the selected chart are not four independent
   unknown scales.
2. The printed order vertex `A'_0=(1,0)`, equivalently `v'_1=(3,0)` for
   `P=phi_0^3+...`, forces `R(0) != 0`. The two zero-root strata retained by
   the previous contact census are not compatible with the F2 normal-form row.
3. The Newton-step inequality forces the selected cofactor-root multiplicity
   to be at least two. Hence the simple roots of `R` do not supply additional
   above-bisectrix F2 continuations. If `R` has a nonzero double root,
   relabeling that squared factor gives a second copy of the same principal
   chain.

Thus the four previous algebraic contact strata reduce to two genuine
normal-form rows:

- `R(0) != 0`, `disc(R) != 0`: one admissible principal F2 chain;
- `R(0) != 0`, `disc(R) = 0`: two conjugacy-independent squared factors,
  hence two copies of the same principal F2 chain.

This is real target-input progress, but it does **not** exclude `(75,125)`.
The subsequent terminal-residue theorem identifies the selected target row
as `(e,f)=(1,6)` with geometric monodromy `A_6`; the remaining problem is to
glue one or two copies into the global log boundary. Simple cofactor roots can
still occur as spectator branches; the result only proves that they are not
additional F2 continuations of the published above-bisectrix chain.

The exact checker is
[`cas/verify_f2_kummer_orbit_transfer.py`](cas/verify_f2_kummer_orbit_transfer.py).
<!-- status-consumer: PF2GC1 33dbc5ff48b5d064 -->

## 1. One Kummer character on every Laurent band

After `x=X^5`, choose the nonzero root `1` and put

\[
u=Xy_{\mathrm{old}}=1+t,\qquad
z=(y_{\mathrm{old}}-X^{-1})^{-1}=\frac Xt.
\]

An original monomial `x^i y_old^j` on Laurent band

\[
\ell=5i-j
\]

becomes

\[
t^\ell u^j z^\ell.
\]

Since `j = -ell (mod 5)`, every coefficient on band `ell` has the unique form

\[
f_\ell(t)=t^\ell u^{k_\ell}A_\ell(u^5),
\qquad
0\le k_\ell<5,
\qquad
k_\ell\equiv-\ell\pmod5. \tag{1}
\]

Let `h_ell` be its exact tangential order at `u=1`. Then

\[
r_\ell=h_\ell-\ell
       =\operatorname{ord}_{v=1}A_\ell(v). \tag{2}
\]

The terminal halfspace gives the lower bounds used by the B0 checker,

\[
h_\ell\ge
\max\!\left(\ell,\left\lceil\frac{17\ell-3}{12}\right\rceil\right)
\quad(P),
\]

\[
h_\ell\ge
\max\!\left(\ell,\left\lceil\frac{17\ell-5}{12}\right\rceil\right)
\quad(Q),
\]

but the transfer below applies to the **actual** order, whether or not it is
strictly larger than this lower bound.

## 2. Exact transfer to every conjugate center

Let `mu^5=1`, put `s=u-mu`, and use the natural translated normal coordinate

\[
z_\mu=(y_{\mathrm{old}}-\mu X^{-1})^{-1}=\frac Xs.
\]

Because `z=(s/t)z_mu`, equation (1) gives the exact identity

\[
t^\ell u^{k_\ell}A_\ell(u^5)z^\ell
=
s^\ell u^{k_\ell}A_\ell(u^5)z_\mu^\ell. \tag{3}
\]

The map `u -> u^5` is etale at every nonzero `mu`, so
`ord_(s=0) A_ell(u^5)=r_ell`. Hence the natural coefficient order at `mu` is
again

\[
\ell+r_\ell=h_\ell. \tag{4}
\]

Therefore every actual Newton point, every edge, and every vertex
nonvanishing condition at the selected root transfers to all five conjugate
roots. The roots are five charts of one Kummer orbit, not five independent
lower-band problems.

More generally, the same statement holds over every nonzero fiber
`u^5=rho`: the five centers have identical natural Newton polygons and descend
to one orbit representative.

## 3. Transfer of the complete terminal block

The selected terminal vertices in `(t,z)` coordinates are

\[
P:(4,3),(21,15),\qquad
Q:(1,1),(18,13),(35,25).
\]

Their band characters and vanishing orders are

| term | `ell` | `k_ell` | `r_ell` | selected coefficient |
|---|---:|---:|---:|---:|
| `P_3` | 3 | 2 | 1 | `1` |
| `P_15` | 15 | 0 | 6 | `1` |
| `Q_1` | 1 | 4 | 0 | `-1` |
| `Q_13` | 13 | 2 | 5 | `-3` |
| `Q_25` | 25 | 0 | 10 | `-9/5` |

At a conjugate root `mu`, the leading coefficient is multiplied by

\[
\mu^{k_\ell+4r_\ell}.
\]

Modulo `mu^5=1`, the transferred block is therefore

\[
P_\mu
=\mu s^4z_\mu^3+\mu^4s^{21}z_\mu^{15},
\]

\[
Q_\mu
=-\mu^4s z_\mu
-3\mu^2s^{18}z_\mu^{13}
-\frac95s^{35}z_\mu^{25}. \tag{5}
\]

Using

\[
[p_i(s)z^i,q_j(s)z^j]
=\bigl(i p_iq_j'-j p_i'q_j\bigr)z^{i+j},
\]

all non-target terms cancel modulo `mu^5-1`, and

\[
[P_\mu,Q_\mu]=s^4z_\mu^4=X^4. \tag{6}
\]

Thus the complete terminal principal edge, not just its top contact, is
forced on the full selected Kummer orbit.

## 4. The zero-root strata are incompatible with `A'_0=(1,0)`

The Newton-resolution description of the `D=75`, `lambda_0=5/3` row writes

\[
\varphi_0=x(xy^5-r_1)^2q_0(xy^5),
\qquad \deg q_0=2,
\]

and has order vertex `(1,0)` for `phi_0`, equivalently `(3,0)` for
`P(e'_0)=phi_0^3`.

In the repository normalization `r_1=1` and `q_0=R`. If `R(0)=0`, then
`R(v)=v R_tilde(v)`, so

\[
\varphi_0=x^2y^5(xy^5-1)^2\widetilde R(xy^5).
\]

Its lowest `y` order is five, and the lowest `y` order of `phi_0^3` is
fifteen. The required vertices `(1,0)` and `(3,0)` disappear. Consequently

\[
\boxed{R(0)\ne0.} \tag{7}
\]

The contact rows previously labelled `R(0)=0`, both simple and double, must
be deleted before any boundary transfer.

## 5. Simple cofactor roots are not additional F2 continuations

For this row the Newton-resolution parameters are

\[
d_0=3,\quad v_0=(15,60),\quad v'_1=(3,0),
\quad(\beta_0,\gamma_0)=(1,5).
\]

The exact bound for the chosen root multiplicity `t_2` is

\[
\frac{\mu'_1\gamma_0-\nu'_1\beta_0}
     {(\gamma_0-\beta_0)d_0}
<t_2
\le
\frac{\nu_0-\nu'_1}{d_0\gamma_0},
\]

hence

\[
\frac54<t_2\le4. \tag{8}
\]

Thus `t_2=1`, corresponding to a simple cofactor root, is impossible. The
desired row has `t_2=2` and

\[
v_1=\left(\frac{21}{5},6\right)
=3\left(\frac75,2\right).
\]

The other allowed multiplicities `3,4` belong to different leading-factor
patterns/rows; the `D=75`, `lambda_0=5/3` table has a unique `t_2=2`
principal row, ending at

\[
v'_2=\left(\frac45,1\right).
\]

Therefore the simple roots of `R` are not additional choices for the
published F2 continuation. If `R` has a nonzero double root `rho`, then

\[
\varphi_0=x(v-1)^2\cdot a(v-\rho)^2
\]

can be relabelled with either squared factor selected. Both selections have
the same unique `t_2=2` chain and terminal block.

## 6. Independent principal-endpoint arithmetic

There is also a direct local check. Let the natural common root at a nonzero
orbit have leading monomial `s^e z^5`. Then the top points of `P,Q` are

\[
(3e,15),\qquad(5e,25).
\]

Suppose a principal pair ends at `(a,b)` and `(c,d)`, and its low-low bracket
is the target `s^4z^4`. Since every natural band coefficient has order at
least its band,

\[
r_P=a-b\ge0,\qquad r_Q=c-d\ge0.
\]

The target exponents imply

\[
a+c=5,\qquad b+d=4,\qquad r_P+r_Q=1.
\]

Writing `r_P` in `{0,1}`, parallelism of the two principal edges reduces to

\[
b(8e-41)=12e+36r_P-75. \tag{9}
\]

For the F2 orbit orders

\[
e=6\quad\text{(simple nonzero cofactor root)},\qquad
e=7\quad\text{(squared nonzero factor)},
\]

equation (9) has no integral solution for `e=6`, while for `e=7` it has the
unique solution

\[
r_P=1,\quad
(a,b)=(4,3),\quad
(c,d)=(1,1). \tag{10}
\]

This independently recovers the terminal endpoints and shows why only a
squared nonzero factor can carry the direct principal target block.

## 7. Revised boundary handoff

The cover-level contact census remains a correct factorization identity, but
it is no longer the correct list of independent branch-scale unknowns.

| `R` stratum | cover centers | nonzero Kummer orbits | admissible principal F2 chains |
|---|---:|---:|---:|
| `R(0) != 0`, `Delta != 0` | 15 | 3 | 1 |
| `R(0) != 0`, `Delta = 0` | 10 | 2 | 2 |
| `R(0) = 0` | — | — | excluded by (7) |

In the first row, the two simple `R` orbits are spectator orbits rather than
additional above-bisectrix F2 continuations. In the second row, the two
squared factors give two copies of the same known principal chain.

The next honest calculation is therefore no longer “recover the first normal
order at fifteen centers.”  The terminal calculation has also completed the
local target row, so the live tasks are:

1. attach the certified source ray `(12,-17)` and target ray `(5,2)` to the
   global completions for one principal chain;
2. decide whether the two packets for `Delta=0` lie over the same target
   divisor, which would force geometric degree at least twelve, or over
   distinct target divisors, which only forces degree at least six;
3. prove how the simple spectator orbits meet—or avoid—the distinguished
   finite-normalization boundary;
4. attach source-boundary branches at the three interior preimages of the
   target toric nodes, carrying different contributions `(4,2,2)`, and place
   the last contribution `2` at the source endpoint over the smooth branch
   value;
5. place the separate purity-forced affine ramification row;
6. run the global dicritical, class-group, unit, and meridian filters.

The current terminal row is centered at target infinity, so the affine-sheet
increment over a nonproperness curve cannot be added to its degree-six
contribution.  This is a substantial reduction of the boundary input, but the
global steps are not completed here.

## Reproduction

```bash
.venv/bin/python plane-jc/cas/verify_f2_kummer_orbit_transfer.py
```

Expected markers:

```text
F2_KUMMER_BAND_TRANSFER_PASS
F2_FIXED_ORBIT_TERMINAL_BLOCK_PASS
F2_ORDER_VERTEX_ZERO_ROOT_EXCLUSION_PASS
F2_NEWTON_MULTIPLICITY_GATE_PASS
F2_PRINCIPAL_ENDPOINT_FILTER_PASS
F2_KUMMER_ORBIT_REDUCTION_PASS
```
