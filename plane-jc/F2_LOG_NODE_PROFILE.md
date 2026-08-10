# F2 terminal logarithmic node profile

> **Status.**  This note proves an exact local theorem for the terminal F2
> `(75,125)` row.  Each of its three interior residue points over a target
> node has transverse coefficient order `-2`, not `0`.  Consequently the
> fixed target chart requires two successive source blowups at each point.
> After those blowups, the full logarithmic differential is invertible at
> every source node lying over either target node.  This corrects the earlier
> `16/25`-component one-blowup attachment skeleton to `19/31` components and
> shows that these terminal node slots cannot carry the sought normalization
> defect.  A boundary-support lemma also closes the remaining smooth target
> endpoint: its cokernel is `R/(w^3)` on a smooth reduced support, so its
> normalization defect is zero as well.  The result does not compute the full
> global log matrix or exclude `(75,125)`.

The rational orders and local matrices are independently checked by
[`verify_f2_log_node_profiles.py`](../scripts/verify_f2_log_node_profiles.py).
The corrected global trees are reproduced by
[`cas/verify_f2_75_125_global_attachment.py`](cas/verify_f2_75_125_global_attachment.py).

## 1. A local node-profile theorem

Let `(tau,w)` be regular parameters on a smooth source surface, with terminal
boundary `tau=0`, and let `(pi,xi)` be toroidal parameters at a target SNC
node.  Suppose that, up to units and terms of larger `(tau,w)`-order,

\[
 \pi=\tau w^p\alpha,\qquad
 \xi=w^e\beta,\qquad \alpha(0,0)\beta(0,0)\ne0.       \tag{1.1}
\]

Here `e>0`, and the characteristic is zero.

### Theorem 1.1 -- transverse-pole regularization

If `p>=0`, the logarithmic differential at the node has reduction

\[
 \overline\Theta=
 \begin{pmatrix}1&p\\0&e\end{pmatrix}.             \tag{1.2}
\]

If `p=-m<0`, fewer than `m` successive blowups in the terminal direction do
not make `pi` regular in this target chart.  After `m` blowups, the chart

\[
 \tau=u w^m                                             \tag{1.3}
\]

gives

\[
 \pi=u\alpha',\qquad \xi=w^e\beta',\qquad
 \overline\Theta=
 \begin{pmatrix}1&0\\0&e\end{pmatrix}.             \tag{1.4}
\]

In either case the determinant is the nonzero scalar `e`.  Hence the
logarithmic differential is an isomorphism in a neighbourhood of the
resolved node,

\[
 \operatorname{Fitt}_0(\mathcal T_f^{\log})
 =\operatorname{Fitt}_1(\mathcal T_f^{\log})=\mathcal O,
 \qquad \mathcal T_f^{\log}=0.                       \tag{1.5}
\]

#### Proof

In logarithmic bases, write

\[
 d\log\pi=d\log\tau+p\,d\log w+d\log\alpha,
 \qquad
 d\log\xi=e\,d\log w+d\log\beta.                 \tag{1.6}
\]

The logarithmic derivatives of a unit vanish modulo the maximal ideal, which
gives (1.2).  If `p=-m`, after `j` terminal blowups the order of `pi` along
the new exceptional coordinate is `-m+j`; it remains negative for `j<m`.
Substitution (1.3) gives (1.4).  Its determinant is a unit in characteristic
zero, proving (1.5).  Higher-order terms cannot change an invertible residue
matrix.  \(\square\)

There is a useful asymmetric companion.  If the target has only the smooth
boundary `pi=0`, with ordinary transverse coordinate `z`, and

\[
 \pi=\tau w^p\alpha,\qquad z|_{\tau=0}=w^e\beta,
                                                               \tag{1.7}
\]

then the coefficient of `dlog(tau)` in `dlog(pi)` is a unit.  Thus
`Fitt_1(T_f^log)=O` and the cokernel is cyclic.  The restriction of its
determinant to `tau=0` has leading term `e*w^e`.  Transverse terms in `z`,
however, can change the reduced determinant curve, so (1.7) alone does not
determine its normalization defect.

For a Keller compactification there is an additional closure.  The
logarithmic determinant is nonzero away from the source boundary because the
affine differential is invertible.  At a source node with boundary
`tau*w=0`, unique factorization therefore gives

\[
 \det\Theta=\tau^a w^b\cdot\text{unit}.             \tag{1.8}
\]

If the log differential is generically invertible along `tau=0`, then
`a=0`.  If (1.7) restricts the determinant to `w^e*unit` on `tau=0`, then
`b=e`.  The unit `Fitt_1` entry reduces the matrix by invertible row and
column operations to

\[
 \operatorname{diag}(1,w^e),\qquad
 \mathcal T_f^{\log}\simeq R/(w^e).                \tag{1.9}
\]

Its reduced support `w=0` is smooth, so its normalization defect is zero.

## 2. Exact target-node coordinates

For the terminal F2 block put

\[
 s=X^{17}y^5,\quad
 P=X^4y(1+s),\quad
 R=-Q=X A(s),\quad
 A(s)=1+3s+\frac95s^2.                              \tag{2.1}
\]

Near the `Q`-dominant target point set

\[
 a=R^{-1},\qquad b=P/R,qquad
 h=\frac{b^5}{a^2}
   =\frac{s(1+s)^5}{A(s)^3}
   =\frac{125s(1+s)^5}{(9s^2+15s+5)^3}.             \tag{2.2}
\]

The terminal source valuation has `nu(X)=-5`, `nu(y)=17`.  On the target
node `h=0`, regular toric coordinates are

\[
 \pi_0=\frac a{b^2}=\frac R{P^2}
      =\underbrace{X^{-7}y^{-2}}_{\tau_0}
       \frac{A(s)}{(1+s)^2},qquad \xi_0=h.          \tag{2.3}
\]

On the node `h=infinity`, use

\[
 \pi_\infty=\frac{b^3}{a}=\frac{P^3}{R^2}
      =\underbrace{X^{10}y^3}_{\tau_\infty}
       \frac{(1+s)^3}{A(s)^2},qquad
 \xi_\infty=\eta=h^{-1}
      =\frac{A(s)^3}{s(1+s)^5}.                    \tag{2.4}
\]

Both `tau_0` and `tau_infinity` have terminal order one.  The exact local
orders are therefore

| source point | target | order of transverse coefficient | residue index |
| --- | --- | ---: | ---: |
| `s=0` | node `h=0` | `0` | `1` |
| `s=-1` | node `h=0` | `-2` | `5` |
| either root of `A` | node `h=infinity` | `-2` | `3` |
| `s=infinity` | smooth point `h=125/729` | `+1` | `3` |

The pole order `-2` is the datum missed by the earlier incidence-only
attachment compiler.  Residue ramification alone records the second column
of the local exponent matrix but not the regularity of its first column.

This conclusion is stable under all unexposed lower Laurent bands.  At
`s=-1`, the pole-normalized source section defining `P` has a simple zero on
the terminal divisor; at a root of `A`, the corresponding section defining
`-Q` has a simple zero.  In either case the regular local implicit-function
theorem makes the full section itself a transverse parameter
`W=w+tau*gamma`.  Replacing `w` by `W` absorbs every higher normal term and
puts the **full** target monomials in the form

\[
 (\pi,\xi)=(\tau W^{-2}\cdot\text{unit},
             W^e\cdot\text{unit}).                \tag{2.5}
\]

Thus the two-blowup count and the invertible residue matrices do not assume
that the displayed terminal block is the whole polynomial map.  This
stability uses the simple zero; it is not available automatically at the
smooth endpoint `s=infinity`.

## 3. Corrected source attachment chains

At each of the three interior points the first blowup gives
`tau=w*u_1`, so `pi=u_1*w^-1*unit` is still not regular in the fixed target
node chart.  Blowing up the intersection of the strict terminal component
with that exceptional curve gives

\[
 \tau=w^2u,qquad
 (\pi,\xi)=(u\cdot\text{unit},w^e\cdot\text{unit}). \tag{3.1}
\]

Thus every interior slot contributes a two-component chain

```text
terminal(-7 after all three slots) -- inner(-1) -- outer(-2).
```

The labels show the final weights after the second blowup at that slot.  The
second blowup is again centered on the strict terminal component.  Hence the
terminal self-intersection falls from `-1` to `-7`, not `-4`.  The terminal
still has valency five.

The corrected principal source lower bounds are:

| F2 case | components | leaves | terminal weights | determinant | inertia |
| --- | ---: | ---: | --- | ---: | --- |
| squarefree, one packet | 19 | 6 | `(-7)` | 1 | `(1,18,0)` |
| double root, two packets | 31 | 10 | `(-7,-7)` | 1 | `(1,30,0)` |

The canonical coefficients remain integral.  On each new chain, from the
terminal outwards, they are `(-16,-8)`.  These are still lower bounds: simple
spectators, the purity row, and other global resolution centers may add
components.

## 4. The two nodes in each interior chain

The terminal--inner node maps to the original target node.  Its exponent
matrix is `diag(1,e)` and its log cokernel vanishes by Theorem 1.1.

The inner--outer node maps to the next target node in the already extracted
minimal fan.  No extra target blowup is required.  In local target fan bases,
the exponent matrices are

\[
 \begin{array}{c|c|c}
 \text{residue point}&\text{matrix}&\det\\ \hline
 s=-1&\begin{pmatrix}5&2\\0&1\end{pmatrix}&5\\[3pt]
 A(s)=0&\begin{pmatrix}3&1\\0&1\end{pmatrix}&3.
 \end{array}                                           \tag{4.1}
\]

They are invertible in characteristic zero as well.  Consequently all six
new source nodes per terminal packet are logarithmically etale.  Together
with the unramified endpoint `s=0`, all four terminal preimages of target
nodes have zero logarithmic cotangent cokernel and zero normalization defect.

At `s=infinity`, one instead has a source node over a smooth target-boundary
point.  Equations (2.4) give

\[
 \operatorname{ord}_{w=1/s}
  \frac{(1+s)^3}{A(s)^2}=1,qquad
 \operatorname{ord}_{w=1/s}
  \left(\eta-\frac{729}{125}\right)=3.             \tag{4.2}
\]

The full matrix is cyclic there because `Fitt_1=O`.  The Keller determinant
is supported on the boundary, and its generic order on the terminal
component is zero because the residue map is generically etale.  Applying
(1.8)--(1.9) to (4.2) gives the exact local presentation

\[
 \boxed{\mathcal T_f^{\log}\simeq R/(w^3).}         \tag{4.3}
\]

The reduced support `w=0` is smooth.  Hence the normalization defect and
nodal excess are zero at this endpoint too.  The multiplicity-three
divisorial contribution remains part of the global Chern ledger, but it is
not a conductor mismatch.

## 5. Consequence for the universal programme

The nodal `Fitt_1`/local-cohomology programme remains the correct universal
replacement for scalar determinant matching, but the flagship terminal
`A_6` packet does not itself supply the desired nonzero node class.  Its
target-node ramification is tame Kummer geometry; once the transverse base
points are actually resolved, logarithmic differentials absorb the residue
ramification and become invertible.

The next exact target is therefore narrower:

1. compile the same two-column node profile at the carrier/arm and spectator
   nodes, where the exponent determinant may vanish;
2. compute a nonzero normalization-defect class only at the surviving
   non-log-etale nodes; and
3. separately prove the Keller descent identity that would force that class
   to vanish.

This is real pruning: a universal theorem cannot close `JC(2)` by assigning a
positive nodal defect to any of the five marked slots of the terminal Belyi
cover: four are log-etale over target nodes, and the fifth is a cyclic
thickening of a smooth boundary branch with zero normalization defect.
