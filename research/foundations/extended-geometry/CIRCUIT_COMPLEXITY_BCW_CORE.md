# Circuit complexity, BCW output rank, and the feedback core

## Result

The proposed circuit-complexity principle has a rigorous downstream half.
Once a determinant-preserving degree-lowering circuit has produced

\[
K=X+Q+C\colon \mathbb A^N\longrightarrow\mathbb A^N,
\]

the cubic dimension, Jacobian rank, and nilpotency index are controlled by
constant output spaces.  The nilpotency invariant is not the raw
multiplication DAG.  It is a smaller **output-feedback core** obtained after
the circuit outputs have been placed back into map coordinates.

This gives a precise theorem, an exact foundational-map test, and a corrected
form of the general conjecture.

The classical determinant-preserving degree-lowering step is the BCW
suspension

\[
(f_1,\ldots,f_n,y,z)\sim
\bigl(f_1-(y+a)(z+b),f_2,\ldots,f_n,y+a,z+b\bigr).
\]

It removes a product \(ab\) while preserving stable left--right equivalence.
The formula and the subsequent Segre/cubic-homogeneous steps are recorded in
[Bass--Connell--Wright](https://doi.org/10.1090/S0273-0979-1982-15032-7)
and explicitly summarized in
[Campbell, Theorem 5](https://doi.org/10.4064/ap110-1-1).

## 1. Output-feedback theorem

Let \(k\) have characteristic zero and let

\[
\mathcal F=I+\mathcal H\colon\mathbb A^d\longrightarrow\mathbb A^d
\]

be cubic homogeneous and Keller.  Suppose the coefficient vectors of
\(\mathcal H\) span a constant \(s\)-dimensional output space.  Equivalently,
choose a full-column-rank constant matrix \(B\in k^{d\times s}\) and a
polynomial vector \(p\in k[X]^s\) such that

\[
\mathcal H=Bp.
\]

Define the \(s\times s\) feedback core

\[
R=Jp\,B.
\]

Then

\[
J\mathcal H=BJp,\qquad
(J\mathcal H)^k=B\,R^{k-1}Jp\quad(k\ge1).             \tag{1.1}
\]

The proof is one-line induction.  It yields

\[
\boxed{\operatorname{rank}J\mathcal H\le s},\qquad
\boxed{\nu(J\mathcal H)\le \nu(R)+1\le s+1}.          \tag{1.2}
\]

The final inequality also follows from the usual
\(\nu\le\operatorname{rank}+1\) bound, but (1.1) is stronger: it moves every
power calculation from the ambient Jacobian to the output-feedback core.

There is a graph refinement.  Form the directed block graph of \(R\), with an
edge \(j\to i\) when the block \(R_{ij}\) is nonzero.  Collapse strongly
connected components.  If the diagonal block on a component \(\alpha\) has
nilpotency index \(\mu_\alpha\), block multiplication gives

\[
\nu(R)\le
\max_{\text{directed paths }P}
\sum_{\alpha\in P}\mu_\alpha.                        \tag{1.3}
\]

For an acyclic scalar feedback graph, this is the number of vertices on a
longest directed path.  Thus a longest-path theorem is valid for the
**feedback graph**.  It is not valid merely from the syntax DAG of the
arithmetic circuit.

## 2. Rank-compressed BCW corollary

Write \(C=B_Cc\), where the components of \(c\) are independent and
\(t=\dim\operatorname{span}\{C_i\}\).  The repository's rank-compressed
homogenization is

\[
V=(X,Y,T)+\bigl(TQ(X)+T^2B_CY,\,-c(X),\,0\bigr)       \tag{2.1}
\]

in dimension

\[
\boxed{N+t+1}.                                        \tag{2.2}
\]

Let

\[
q=\dim\operatorname{span}\{Q_i\},\qquad
s=\dim\bigl(\operatorname{out}(Q)+\operatorname{out}(C)\bigr).
\]

The values of the homogeneous correction in (2.1) lie in a constant output
space of dimension \(s+t\).  Hence

\[
\boxed{\operatorname{rank}J\mathcal H\le s+t\le q+2t},\qquad
\boxed{\nu(J\mathcal H)\le \nu(R)+1\le s+t+1}.         \tag{2.3}
\]

This proves the rank/index part of the proposed theorem after replacing a
vague “output tensor rank” by the two exact widths \(q,t\), or by their joint
width \(s\).

Suppose a circuit compiler introduces \(a(\mathcal C)\) variables before
homogenization and produces widths \(q(\mathcal C),t(\mathcal C)\).  Then the
fully proved bounds are

\[
\begin{aligned}
n_{\rm cubic}
 &\le n+a(\mathcal C)+t(\mathcal C)+1,\\
\operatorname{rank}J\mathcal H
 &\le q(\mathcal C)+2t(\mathcal C),\\
\nu(J\mathcal H)
 &\le q(\mathcal C)+2t(\mathcal C)+1,
\end{aligned}                                         \tag{2.4}
\]

with the sharper graph bound (1.3) whenever the feedback core is known.
Consequently the original \(O(m+t)\) conjecture is reduced to one concrete
compiler theorem:

\[
a(\mathcal C),q(\mathcal C),t(\mathcal C)=O(m+t_{\rm out}). \tag{2.5}
\]

No separate raw-degree argument for nilpotency is needed after (2.5).

## 3. Exact foundational endpoint test

The verifier
[`verify_circuit_complexity_bcw_core.py`](../scripts/verify_circuit_complexity_bcw_core.py)
factors the final corrections of the two current endpoint witnesses and
computes their cores over \(\mathbb Q[X]\).

| endpoint | ambient dimension | generic rank \(JH\) | output-core dimension | exact \(\nu(R)\) | exact \(\nu(JH)\) |
|---|---:|---:|---:|---:|---:|
| rank winner | 24 | 17 | 22 | 17 | 18 |
| index winner | 22 | 18 | 21 | 17 | 18 |

In both cases \(R^{16}\ne0\) and \(R^{17}=0\) as polynomial matrices.  Formula
(1.1) therefore gives \((JH)^{18}=0\).  The stored exact nonzero seventeenth
powers give the reverse inequality.  This recovers the endpoint
nilpotency certificates through matrices smaller than the ambient
Jacobians.

Each feedback graph is one strongly connected component.  This is decisive:
the small multiplication depth of the original map does not survive target
placement as an acyclic dependency graph.  Gate sharing and grouped target
shears change this strongly connected feedback block, explaining why they
change the final rank and index profiles.

## 4. A ten-gate circuit and the determinant obstruction

The same verifier records the factorization

\[
\begin{aligned}
a&=xy,& r&=xz,& p_0&=x^2(3y+r),\\
v&=1+2a,& b&=2+3a,&
w&=v^2z+4y^2b,
\end{aligned}
\]

followed by

\[
F=(x-p_0,\ y+3xw,\ vw).                              \tag{4.1}
\]

Counting products with optimal sharing gives ten multiplication gates and
multiplication depth four.  The three stored rational points still have one
common image, and \(\det DF=1\).

There is a tempting but false shortcut: introduce one variable for each gate,
output the gate residuals, and replace the terminal expressions by the gate
variables.  That graph map agrees with \(F\) on the circuit graph, but for
(4.1) its \(13\times13\) affine Jacobian has rank \(12\) and determinant
zero.  Equality of determinants holds only after restricting to the circuit
graph.  This exact failure explains why determinant preservation must be
proved gate-by-gate using BCW source/target suspensions.

Thus the ten-gate circuit does **not** by itself prove a 14-variable cubic
endpoint.  It is instead the smallest useful compiler test.

## 5. Corrected conjecture

An ordinary arithmetic circuit is missing target ownership.  The compiler
state must be a **placed circuit** containing:

1. the hash-consed multiplication DAG;
2. the output covector or target block in which each live gate occurs;
3. exposed-coordinate ownership and live intervals;
4. the quadratic and cubic output spaces;
5. the feedback core \(R=JpB\) after every placement.

Define \(a_{\rm BCW}(\mathcal C)\) as the minimum number of auxiliary
coordinates in a determinant-preserving suspension schedule of a placed
circuit.  The sharpened conjecture is:

> There are absolute constants \(A_i\) such that every normalized
> division-free placed circuit with \(m\) multiplication gates and output
> width \(t_{\rm out}\) admits a cubic terminal schedule with
> \[
> a_{\rm BCW}\le A_1m+A_2t_{\rm out},\quad
> q\le A_3m+A_4t_{\rm out},\quad
> t\le A_5m+A_6t_{\rm out}.
> \]
> Its nilpotency index is bounded by the weighted longest path in the
> resulting feedback-core condensation graph, plus one.

The output-core theorem proves every claimed geometric consequence of this
compiler statement.

The suggested counterexample strategy also needs one correction.  Two
circuits computing the same map cannot have different “unavoidable” cubic
rank if unavoidable means the minimum over all reductions of the map.
What can differ is the rank forced by a specified compiler or schedule.
The correct falsification target is therefore two equal-size, equal-depth
placed circuits for which every ownership-respecting schedule has different
cost or feedback-core width.  That would prove that \(m,h,t_{\rm out}\) must
be supplemented by a placement or sharing-width invariant.

## 6. SAT/SMT formulation

The theorem changes the optimization target.  A bounded exact search should
not minimize only terminal dimension.  At schedule step \(j\), introduce
Boolean or finite-domain variables for:

- whether gate \(g\) is exposed;
- the coordinate owning its exposed output;
- its live interval and reuse count;
- the target block of each grouped cancellation;
- whether a cancellation is admissible (its factor outputs do not own the
  modified target block).

Polynomial identities for grouped shears are exact linear constraints on
their rational coefficients once the combinatorial schedule is fixed.
The terminal objective is lexicographic or Pareto:

\[
\bigl(a_{\rm BCW}+t,\ q+2t,\ 
 \text{largest feedback SCC},\
 \text{weighted feedback-path bound}\bigr).           \tag{6.1}
\]

Rank conditions can be handled by an outer matroid loop: guess bases for the
quadratic and cubic coefficient-row spaces, impose exact reconstruction, and
reject the schedule when the corresponding minors vanish.  This separates
the discrete ownership problem (SAT) from rational coefficient identities
(SMT over \(\mathbb Q\)) and from exact rank certification.

## Reproduction

```bash
.venv/bin/python scripts/verify_circuit_complexity_bcw_core.py
```

The generated exact record is
[`circuit_complexity_bcw_core.json`](../artifacts/generated-results/circuit_complexity_bcw_core.json).
