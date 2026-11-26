# Full power-class rank on connected abc boxes and their least repairs

An actual shared-value box has an additive presentation on a
tripartite graph. Its multiplicative power classes can still be
independent on all edges. A relation with two active edges at one
vertex is tested by recovering that vertex. If the active edges form
a matching, recover two vertices using one grid letter in the actual
additive triangle. This also tests the individual least-repair
multipliers, with a loss independent of the number of vertices.

## Counting all relations, including matchings

Use the [actual shared-value box setup (LB1)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#a-finite-difference-array-theorem).
There are \(m\ge2\) indexed potentials on each of three axes. Write
\(L_e\le X\) for the positive label of each edge \(e\) of the complete
tripartite graph, \(E_0=3m^2\) for the number of edges,
\(R_1,R_2,R_3\) for the residue alphabet bounds, and \(R=R_1R_2R_3\).
For every pair of prescribed opposite residue indices each role has
at most \(K\) possible complete labels. Fix a prime pool \(P\) of prime
product at most \(X^B\), with \(B\) fixed. For \(q\ge2\) and a positive
rational \(z\) put
\[
 {\cal D}_{q,P}(z)=
       \prod_{\substack{\ell\notin P\\q\nmid v_\ell(z)}}\ell .
\]

There are two choices of observables \(z_e\). For the original
labels take \(z_e=L_e\), all \(q\ge2\) are allowed, and put \(\Delta=1\).
For least repairs fix \(M\ge2\), allow \(q\mid M\) with \(q\ge2\), and put
\[
       z_e=Z_M(L_e)=\prod_{\ell\mid L_e}
                       \ell^{(M-v_\ell(L_e))_+},
       \qquad \Delta=X^{1/M}.
\]
Here \(Z_M(L_e)L_e\) is the least repaired value of the complete
letter. The \(z_e\) are the individual multipliers; a repair of any
one original row uses only the multiplier belonging to its
designated target. [Least repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

For \(Z\ge1\), the number of boxes for which some allowed \(q\) and
some integer weight vector \(\mathbf b\) nonzero modulo \(q\) have
\({\cal D}_{q,P}(\prod_e z_e^{b_e})\le Z\) is at most
\[
 \ll_{m,B,\varepsilon} X^\varepsilon R^m K^{3m-1}
 \left[
  \frac{Z\Delta^2}{K}\sum_{i=1}^3\frac1{R_i}
  +\sum_{1\le i<j\le3}\min\left\{
       \frac1{KR_j}+\frac{Z\Delta}{\sqrt{KR_i}},
       \frac1{KR_i}+\frac{Z\Delta}{\sqrt{KR_j}}
                        \right\}
 \right].
 \tag{BX1}
\]
The constants are uniform in \(M\) and \(Z\). There is no circulation
condition on \(\mathbf b\).

Choose a prime \(p\mid q\) for which some \(b_e\) is nonzero
modulo \(p^{v_p(q)}\). Set \(h=\min_e v_p(b_e)<v_p(q)\), treating
zero entries as having infinite valuation. Divide the weights by
\(p^h\). Every prime outside the original defect support still
gives a vanishing valuation sum modulo \(p\). Call an edge active
when its reduced weight is a unit modulo \(p\).

Suppose first that two active edges meet at one vertex \(x\) of axis
\(i\). Fix the other vertices up to translation. Their residue
choices cost at most \(R^m/R_i\), and a tree among them has
\(3m-2\) edges with at most \(K\) choices each. The known graph is
connected because it includes two full other axes. Include in the
pool the supports of every fixed label and every nonzero difference
of neighbors of \(x\). These differences are of polynomial height.
Outside that pool every prime divides at most one form incident to
\(x\), since the forms have slope \(1\) or \(-1\) and distinct
known roots. Also a prime of \(z_e\) is a prime of \(L_e\).
For original labels, any prime on a selected active form off the
pool and the recorded defect support has valuation divisible by
\(p\). For multipliers the same holds whenever its valuation on
\(L_e\) is below \(M\), because \(p\mid M\) and that prime's valuation
in \(z_e\) is \(M-v_\ell(L_e)\). The high valuations contribute to
a squarefree product at most \(X^{1/M}=\Delta\) per selected form.
Thus a single additional support of product at most \(Z\Delta^2\)
contains the remaining canonical \(p\)-free coefficient primes of
both forms.

Enumerate this support and the two coefficients on the resulting
fixed pool. After the support is chosen, the coefficients cost
\(X^{o_X(1)}\) uniformly since they are at most \(X\). At fixed
coefficients the two forms are \(C_1 y_1^p,C_2 y_2^p\). Eliminating
\(x\) gives a fixed nonzero signed sum of powers, because the
two affine roots are distinct. For bounded \(p\) the uniform
binary norm count bounds its solutions by \(X^\varepsilon\);
for larger \(p\le\log_2 X\) choose one positive root in
\([1,X^{1/p}]\). The latter cost is \(X^\varepsilon\) after
fixing the bounded-order cutoff in terms of the allowance.
For \(p>\log_2 X\) both forms are supported on the recorded pool,
so enumerate their bounded values directly. The orders below
this last threshold contribute a logarithmic factor. This yields
the first term of (BX1) by summing over the possible omitted axis
and the bounded choices of the two active forms.
[One-vertex recovery (LB3)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#recovering-one-vertex-in-a-weighted-power-relation)
· [Fixed-value binary norm counts](ABC_RECTANGLE_POWER_CLASS_DEFECT.md#the-fixed-value-bound-with-moving-coefficients)

Otherwise the active edges form a matching. Choose one active
edge \(xy\) between axes \(i,j\). Fix all vertices other than
\(x,y\) up to translation, at cost at most
\(R^m K^{3m-3}/(R_iR_j)\). The known graph is connected since
\(m\ge2\) and both touched axes retain a vertex. Every other active
edge is already fixed. Include all known labels in the pool.
This leaves \({\cal D}_{p,P'}(L_{xy})\le Z\Delta\), with \(P'\)
the enlarged fixed pool. Indeed, outside the pool and the compressed
defect support only the weight at \(xy\) can contribute; for
multipliers pay also the high-valuation tail of \(L_{xy}\).

Fix any known vertex of the third axis. The three labels on this
actual triangle determine \(x,y\) relative to the known potentials,
and are related by the original signed additive equation. The
edge from \(x\) to the known vertex takes at most \(KR_i\)
values, and the edge from \(y\) at most \(KR_j\). This uses one
known residue index and the variable residue of \(x\) or \(y\),
respectively. Apply the arbitrary-order one-grid count (OG6)
to the edge \(xy\), with defect cap \(Z\Delta\) and with those two
dictionaries. The bound includes all feasible \(p\) simultaneously.
Multiplying by the known-label count and summing over \(i,j\)
gives the second term. The choices of indexed vertices are bounded
in terms of \(m\).
[Arbitrary-order additive graph count (OG6)](ABC_ONE_POWER_GRID_ADDITIVE_GRAPHS.md#an-arbitrary-order-defect-on-a-single-letter)

## Consequences at the packet scale

In a surviving simultaneous packet tuple put
\(n=X^L\), \(d=s_1+s_2+s_3\), and \(t=1/6-d>0\).
Write \(s_{\min}=\min_i s_i\) and \(s_{\rm med}\) for the
middle of the three prices. Use the height-first convention
and the fixed positive margins of (LB4). Thus
\(R_i=n^{s_i+O(\kappa+\eta)}\), \(K=n^{t+O(\eta)}\), and the
forced box count has exponent \(md+(3m-1)t\) up to the fixed
tolerance losses. Take \(a=0\) for labels. For multipliers take
\(M\ge K_0/L\) with \(K_0>0\) fixed and put \(a=1/K_0\).
For each fixed \(\rho>0\) satisfying
\[
       \rho<\min\left\{
           t+s_{\min}-2a,\ \frac{t+s_{\rm med}}2-a
                    \right\}
 \tag{BX2}
\]
with a fixed positive margin, all but a power-small fraction of
the forced boxes satisfy
\[
  {\cal D}_{q,P}\!\left(\prod_e z_e^{b_e}\right)>n^\rho
     \quad\text{for every allowed \(q\) and \(\mathbf b\not\equiv0\bmod q\)}.
 \tag{BX3}
\]
Indeed the first term in (BX1) saves by the first margin. For
a matched pair use the denominator \(\sqrt{K\max(R_i,R_j)}\),
and \(\max(s_i,s_j)\ge s_{\rm med}\). Its diagonal term has
denominator at least \(K\), so also saves. The packet tolerances
and \(\varepsilon/L\) are taken sufficiently small in terms of
the margins and \(m\). Each bound is applied as \(X\to\infty\)
at fixed \(L\). No defect allowance grows with the number of
edges. In particular one sufficiently large \(K_0\) works for
any fixed \(m\) with the same positive price margins, allowing
the packet tolerances to depend on \(m\).
[Forced shared-value box count (LB4)](ABC_SHARED_VALUE_RESIDUE_BOXES.md#consequences-for-simultaneous-packet-strata)

Thus at prime orders the valuation map outside the prescribed
pool has full rank \(3m^2\) on the entire edge space, even for
the least multipliers at prime orders dividing \(M\). More
generally the subgroup generated by the \(z_e\) in
\(\mathbb Q_{>0}^{\times}/(\mathbb Q_{>0}^{\times})^q\) has order
\(q^{3m^2}\) for all allowed \(q\). Consequently exact power classes
on a spanning tree cannot determine any chord by a multiplicative
relation. If all valuations of a label are at most \(M\), its
repaired value is \(\operatorname{rad}(L_e)^M\) and can be a
\(q\)-th power even though its least multiplier contributes an
independent power class.

## Stability under base change

Fix a box satisfying (BX3), an allowed \(q\), and \(0\le v<\rho\)
such that the prime product of \(q\) outside \(P\) is at most \(n^v\).
Let \(F\) be any field, allowed to depend on the box. If
\[
       \operatorname{rd}(F)\le n^{(\rho-v)/2},
\]
then
\[
    [F(\mu_q, z_e^{1/q}: e):F(\mu_q)]=q^{3m^2}.
 \tag{BX4}
\]
For least multipliers with fixed \(M\) at a budget, take \(P\)
to include the primes of \(M\). This gives \(v=0\) simultaneously
for all \(q\mid M\).

If an integer vector \(\mathbf b\not\equiv0\bmod q\) gives a
\(q\)-th-power relation in \(F(\mu_q)\), choose \(p,h\) as in
the proof of (BX1). By (BX3) at order \(p\), the reduced
weights \(\mathbf b/p^h\) have defect prime product greater than
\(n^{\rho-v}\) outside \(P\) and the primes dividing \(q\).
At each such rational prime the valuation of the rational
product with the original weights has exact \(p\)-adic order
\(h\). Hence at every place of \(F(\mu_q)\) above it the
absolute ramification index is divisible by \(p\). The
cyclotomic base change is unramified away from rational primes
dividing \(q\), so \(p\) divides every ramification index in
\(F\) at each defect prime \(\ell\). The exponent of \(\ell\)
in the root discriminant is then at least \(1-1/p\ge1/2\).
The product of these contributions contradicts the stated bound.
Thus the \(z_e\) generate a subgroup of order \(q^{3m^2}\) in
\(F(\mu_q)^\times/F(\mu_q)^{\times q}\); Kummer theory gives (BX4).
[Local discriminant inequality (RL1)](ABC_RAMIFIED_IDEAL_LIFT_COST.md#the-local-cost-and-a-field-attaining-it)
· [Kummer theory, Milne Theorem 5.30](https://www.jmilne.org/math/CourseNotes/FT.pdf)

For prime \(q\), no restriction on its prime product is needed:
\(F(\mu_q)/F\) is Galois of degree dividing \(q-1\), so cannot
make a ramification index divisible by \(q\) at any prime. This
gives (BX4) with the larger root-discriminant bound
\(n^{(1-1/q)\rho}\), using all defect primes outside \(P\).
In particular, adjoining radicals from any proper subset of edges
over the specified base field leaves a further extension of degree
\(q^{3m^2-r}\), where \(r\) is the number already adjoined.
These estimates rule out constructing all the root data from
fewer radical generators at small root discriminant. They
do not force this construction or a cheap relation on endpoint boxes.
