# Trace discriminants at primes with unequal root valuations

Fix adjacent ordinary moments \(P_n,P_N\), \(N=n+1\), \(n\ge2\).
Suppose they have a rational common factor with nonzero roots
\(\alpha_i\). A reciprocal trace is
\(t_{ij}=\alpha_i/\alpha_j+\alpha_j/\alpha_i\).
Choose a Galois-invariant set of distinct such traces. Use the
monic rational polynomial \(T\) with these roots
\(t_1,\ldots,t_m\), and let \(D\) be its least positive
denominator. The selection can also be the squarefree part
of an invariant trace multiset. The argument below uses all
primes, including those dividing \(D\).

For a prime \(p\), fix embeddings in an algebraic closure of
\(\mathbb Q_p\) and put \(u_i=(-v_p(t_i))_+\), with \(u_i=0\)
if \(t_i=0\). When \(t_i\) comes from two roots of valuations
\(s,s'\), one has \(u_i=|s-s'|\). Indeed unequal root
valuations give one strictly dominant summand in the trace,
and equal valuations give two integral summands. In particular
\[
                         v_p(D)=\sum_i u_i.
\]
This is the Gauss norm identity from the
[trace height bound](MULTIVARIATE_TRACE_HERMITE_JET_HEIGHT.md).
For two distinct traces define
\[
             e_{p,ij}=v_p(t_i-t_j)+\max(u_i,u_j)\ge0.
\]
Different depths \(u_i\ne u_j\) imply \(e_{p,ij}=0\).
Equal negative trace valuations permit a positive excess
\(e_{p,ij}\), as do two integral traces.

## A discriminant height bound with a denominator-prime debit

**Theorem.** Assume \(m\ge2\). At each prime sort
\(u_1\ge\cdots\ge u_m\). Let \(E_p=\sum_{i<j}e_{p,ij}\),
which is unaffected by this ordering. The primitive integer
polynomial \(DT\) has discriminant valuation
\[
 \frac12v_p\bigl(\operatorname{disc}(DT)\bigr)
              =\sum_{i=1}^m(i-1)u_i+E_p.                     \tag{1}
\]
Suppose the real or complex traces are bounded at the
Archimedean place and \(|t_i-t_j|=O(\eta^{\gamma_{ij}})\)
for some nonnegative exponents, as \(\eta\downarrow0\).
Then
\[
 \sum_p\left[\sum_{i=1}^m(m-i)u_i-E_p\right]\log p
             \ \ge\ \sum_{i<j}\gamma_{ij}\log(1/\eta)-O(1).
                                                                  \tag{2}
\]
The implied constant depends on the finitely many
Archimedean comparison constants. No algebraicity of the
limiting traces is required. At a prime, the left side can
use an upper bound for \(\sum_i(m-i)u_i\) and a lower bound
for \(E_p\), valid jointly over every allowed common-factor
allocation. The prime sets can vary with the profile.

To prove (1), use the discriminant formula
\(\operatorname{disc}(DT)=D^{2m-2}\prod_{i<j}(t_i-t_j)^2\).
In its valuation, \(\sum_{i<j}\max(u_i,u_j)\) equals
\(\sum_{i=1}^m(m-i)u_i\). Thus the exact debit against
\((m-1)v_p(D)\) is
\(\sum_{i<j}\min(u_i,u_j)+E_p\), which is the right side
of (1). Both sides are half a nonnegative integer.
The same formula over \(\mathbb C\), along with the
factorization of the nonzero integral discriminant into
primes, gives (2). In particular, for \(a\) traces at the
same positive depth \(u\), the term
\(\binom a2 u\) is available at that denominator prime even
before testing any residual congruences.

For a cluster at one limiting trace with individual rates
\(\gamma_i\), one can take
\(\gamma_{ij}=\min(\gamma_i,\gamma_j)\) within the cluster.
Between different clusters use rate zero. The resulting bound
refines the [ordinary discriminant comparison](MULTIVARIATE_TRACE_HERMITE_JET_HEIGHT.md#comparison-with-the-elementary-discriminant-bound)
by retaining the exact non-Archimedean contribution. The
refinement also covers a transcendental reference shape, as
can occur in a compactified infinite core.

## The local excess from two normalized quotients

Choose a uniformizer scale \(p^s\) coherently for the finitely
many rational valuations of the roots at a given prime. When
two traces have the same positive depth \(u\), orient both
root quotients \(q,q'\) to have valuation \(-u\) and put
\(U=p^u q\), \(V=p^u q'\). These are units. Directly,
\[
   p^u\bigl((q+q^{-1})-(q'+{q'}^{-1})\bigr)
                    =(U-V)\bigl(1-p^{2u}/(UV)\bigr).
\]
Since the second factor is a unit,
\[
                             e_{p,ij}=v_p(U-V).                \tag{3}
\]
For two depth-zero traces, both root quotients \(U,V\) are
units and the same factorization gives
\[
                       e_{p,ij}=v_p(U-V)+v_p(U-V^{-1}).        \tag{4}
\]
These formulas hold in every residue characteristic. Equation
(3) measures a coincidence at a prime of the trace denominator
without losing precision at depth \(2u\). A sufficient ordinary
certificate for a positive bound in (3) or (4) is any bound
on the indicated difference of normalized units obtained from
the two [physical truncated Newton faces](MULTIVARIATE_PHYSICAL_TRUNCATED_NEWTON_FACES.md).
For example, if both normalized quotients agree to depth \(L\)
with the same algebraic unit, then (3) is at least \(L\).
The estimate has no requirement that all roots lie on one
slope, or that all reductions of the traces are integral.

## A selective product for Galois-invariant collections of pairs

One can replace the complete graph on the distinct traces by any
Galois-invariant graph of pairs of traces, with nonnegative
integer edge multiplicities. Here the graph lies on the actual
values and loops are excluded. The proof works even with additional
fixed algebraic vertices, provided their full Galois orbits are
included and edges join unequal values. For an edge \(ij\),
use \(\gamma_{ij}\) and \(e_{p,ij}\) as above. Then
\[
 \sum_p\sum_{ij\ {\rm edges}}
       \bigl(\max(u_i,u_j)-e_{p,ij}\bigr)\log p
       \ \ge\ \sum_{ij\ {\rm edges}}\gamma_{ij}\log(1/\eta)-O(1).
                                                                  \tag{5}
\]
Indeed \(\prod_{ij\ {\rm edges}}(t_i-t_j)^2\) is nonzero and
rational since squaring removes the sign of every edge reversal
under Galois. Its exact prime factorization gives (5). The bounded
error in (5) uses fixed multiplicities, or uniformly bounded
weighted sums of the logarithms of the comparison constants.
Equivalently, clear denominators with the square of the
least denominator of the monic polynomial
\(\prod_i(Y-t_i)^{\deg(i)}\). At each prime the resulting
integer has valuation
\(2\sum_{ij\ {\rm edges}}(\min(u_i,u_j)+e_{p,ij})\).
This follows from the Gauss norm identity for that polynomial
and the valuation of each difference. The graph keeps distinct
Galois orbits of pairs separate. In particular a graph between
two disjoint invariant trace selections uses their resultant,
while an orbital contained in one selection can use fewer
edges than a discriminant.

## A strict example from traces of an irreducible quartic

Let \(p\) be an odd prime, \(h\ge1\), and
\(t_\pm=1\pm\sqrt2/p^h\). The two traces occur from the
conjugate pairs of the irreducible rational quartic
\[
                  H(z)=\prod_{\pm}
                         \bigl(z^2+(t_\pm+2)z+(t_\pm+2)\bigr).
\]
Indeed the quadratic discriminants
\((t_\pm+2)(t_\pm-2)\) are negative at the respective real
embeddings. A root \(z\) of either quadratic generates its
\(t_\pm\) via \(t_\pm+2=-z^2/(z+1)\), so its degree over
\(\mathbb Q\) is four. The two roots of that quadratic have
trace \(t_\pm\). The conjugate-pair selection is invariant
because it is exactly the pair of roots over each value of
\(t_\pm\).

Here \(T(Y)=(Y-1)^2-2/p^{2h}\), \(D=p^{2h}\), and
\(\operatorname{disc}(DT)=8p^{2h}\). Both traces have
\(u_i=h\) and \(e_{p,12}=0\). With \(\eta=p^{-h}\), their
real difference has rate one. The bound using only
\(\log D\) gives an available budget \(2h\log p\).
The denominator-prime debit removes \(h\log p\), leaving
\(h\log p+O(1)\), sharp to a constant. This example concerns
rational factor geometry and requires no existence claim for
a common factorial factor.

In applications to factorial moments, the two parent Newton
polygons and local common-factor tests supply allowable
root valuations and bounds on normalized unit differences.
Inequalities (2) and (5) are exact necessary bounds even
with moving denominator primes and arbitrarily many diffuse
rows. A uniform bound on their left sides at all remaining
contacts is a separate arithmetic problem. A single trace
without a second invariant value gives no discriminant debit.

