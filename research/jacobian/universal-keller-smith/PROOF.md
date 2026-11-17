# Tame fixed loci and nonlinear canonical Keller costs

14 November 2026. This joins the retained [general source-deck rigidity](../monodromy-theory/GENERAL_DECK_RIGIDITY.md), [fixed-locus gate](../universal-keller-symmetry/PROOF.md) and [linear canonical costs](../universal-keller-costs/PROOF.md). Smith theory makes the low-dimensional fixed geometry automatic for tame prime groups, so the new dimension bounds require **arbitrary polynomial actions**, with no linearizability hypothesis. The characteristic-zero construction problem remains open.

## 1. Setup

Let F:X=A^D -> Y=A^D be a polynomial Keller map over k. Let a finite constant group G act on both spaces by polynomial automorphisms and make F equivariant. At a G-fixed geometric target y suppose the **whole** fibre is a finite etale faithful G set Omega. No equality between its rank N and the geometric generic degree M is assumed. Work over an algebraic closure when taking tangent spaces and fixed schemes.

A canonical marked affine-space realization over BG has this split-object form. The costs below concern that marked category. An unmarked fixed universal map need not carry these G actions or a fixed target, so its target-dimension invariant is a different question.

Write O_p(G) for the largest normal p-subgroup, and use O_0(G)=1. Smith theory is invoked only for an ell-group with ell different from the field characteristic.

## 2. UADECKWILD1: only wild symmetry can disappear

The monoid of polynomial source endomorphisms sigma satisfying F sigma=F is a finite group D(F). In characteristic zero it is trivial (KMFI47). In characteristic p it is a p-group, and
\[
             |D(F)|\mid M.                         \tag{1}
\]
The assertions hold after every constant extension and identity stabilization. In positive characteristic this does not exclude nontrivial or nonabelian p-groups, and it does not identify regular deck maps with all rational generic deck maps.

If G acts faithfully on X, then
\[
 \ker(G\to\operatorname{Aut}(Y))\subseteq O_p(G),
 \qquad
 \ker(G\to\operatorname{GL}(T_yY))\subseteq O_p(G). \tag{2}
\]
The first kernel has order dividing M. The degree divisibility is **not** asserted for the second kernel. In characteristic zero both actions on the target and on its tangent space are faithful. For the natural faithful S_N fibre, N>=5, O_p(S_N)=1 for every p, so target tangent faithfulness follows in every characteristic, without linearizing either global action.

### Finite group, free action, tame fixed points

The function-field argument in KMFI47 is characteristic-free. Dominance and F sigma=F give an injective K-linear endomorphism of the finite separable extension L/K, K=k(F). It is an automorphism of finite order; hence sigma has a polynomial inverse. All such sigma inject into Aut_K(L). The fixed-field degree formula gives |D(F)| dividing [L:K]=M.

Any nonidentity regular deck map has no fixed geometric point. Indeed F is etale and its diagonal is open and closed; the equalizer of sigma and identity over F is consequently open and closed in the connected source. A fixed point would make sigma the identity. Equivalently one may use the etale local inverse.

If D(F) had an element of prime order ell different from the characteristic, its cyclic ell-group would have a fixed point on affine space by tame Smith theory. This contradicts freeness. Cauchy's theorem therefore makes D(F) a p-group in characteristic p and trivial in characteristic zero. This extends the conclusion of KMFI47 to the exact positive-characteristic torsion restriction; it does not claim a new Smith fixed-point theorem.

For (2), the target kernel injects into D(F). At a fixed point the kernel of a faithful finite action's tangent representation is also a p-group. To see this, an element of prime order ell invertible in k cannot have identity linear part unless it is the identity: its first nonzero homogeneous term of degree j>=2 would be multiplied by ell in its ell-th iterate. If there is no first term, the action is formally the identity and then globally the identity. Cauchy's theorem finishes the assertion about that kernel. The inverse image in G of the target tangent kernel is an extension of two p-groups, hence a normal p-group. This proves (2).

The small symmetric exceptions match the modular constructions: O_2(S2)=S2, O_3(S3)=A3 and O_2(S4)=V4. Other characteristic p cores in these ranks are trivial. The canonical Frobenius line and plane models can kill exactly these translation kernels. High-rank S_N has no such kernel.

## 3. UASMITHFIX1: the whole low-dimensional fixed scheme is automatic

For every prime ell different from char k and every ell-subgroup H of G, the whole schemes X^H and Y^H are smooth, geometrically connected and mod-ell etale-cohomologically acyclic. In particular they are nonempty. Their common dimension is
\[
             t_H=\dim(T_yY)^H.                    \tag{3}
\]
The restriction F^H:X^H -> Y^H is etale, and its complete marked fibre is Omega^H. Put r_H=|Omega^H|. Then:
\[
\begin{array}{c|c|c}
 &\operatorname{char}k=0&\operatorname{char}k=p>0\\ \hline
 t_H=0&r_H=1&r_H=1\\
 t_H=1&r_H=1&r_H=1\text{ or }r_H\in p\mathbb Z_{>0}.
\end{array}                                                   \tag{4}
\]
Thus in characteristic zero every H with r_H!=1 requires t_H>=2. In characteristic p the same holds when r_H=0 or r_H>1 is not divisible by p. These conclusions do not assume that higher-dimensional fixed loci are affine spaces.

### Published input and the curve argument

The input is [Haution, Proposition3.2.1](https://arxiv.org/html/1507.04582): a tame ell-group acting on affine space has smooth geometrically connected fixed scheme, with trivial positive-degree etale cohomology. Its proof invokes classical Smith theory with Z/ell coefficients and then passes to ell-adic coefficients. This applies in every base characteristic different from ell; no topological embedding of the base into C is needed.

At any fixed source point DF and its inverse identify invariant tangent spaces. Since the whole fixed schemes are smooth and connected, they have equal constant dimension and the restricted differential is invertible. The target contains y, giving (3). The complete fibre is the fixed subscheme of the marked finite-etale fibre, hence exactly Omega^H.

A zero-dimensional smooth connected fixed scheme over the algebraic closure is a point. A one-dimensional one is A1: its smooth compactification has genus g and r>=1 points at infinity, and vanishing H1 gives 2g+r-1=0. Thus g=0,r=1. Consequently for t_H=1 the restriction is a polynomial A1->A1 with everywhere nonzero derivative, hence constant nonzero derivative. Its degree is one in characteristic zero; in characteristic p it is one or a positive multiple of p. Every whole geometric fibre has that positive cardinality. This proves (4), including the empty-fibre obstruction.

In characteristic zero, (2)--(4) give a finite representation obstruction for arbitrary polynomial actions: the target tangent representation W is faithful and
\[
      \dim W^H\ge2\quad\text{whenever }|\Omega^H|\ne1
      \quad(H\text{ any prime-power subgroup}).    \tag{5}
\]
It is stronger than merely assuming a nonlinear action might avoid the linear fixed-space test. In positive characteristic only tame subgroups and the exact alternatives in (4) may be used.

## 4. A classical low-degree representation lemma

Let V_n be the standard (sum-zero) representation of S_n in characteristic zero. For n>=7, the only irreducible representations of dimension at most n+1 are
\[
       1,\quad\mathrm{sgn},\quad V_n,\quad V_n\otimes\mathrm{sgn}.
                                                               \tag{6}
\]
We give a short all-n argument from classical Young branching and hook lengths, rather than rely on a bounded character-table search or import the full minimal-degree classification.

Call the four corresponding partitions exceptional. For n=7 and n=8, the hook formula gives minimum nonexceptional degree14, greater than n+1. Induct on n>=9. If a nonexceptional Young diagram has at least two removable corners, one child is nonexceptional: otherwise a row or column child makes the parent exceptional; a standard child leaves only parents (n-2,2), (n-2,1,1), or their transposes, each having a nonexceptional child. Induction gives this child degree greater than n, and another positive child degree makes the parent degree greater than n+1.

A rectangle with at least two rows and columns has one child and exactly two grandchildren. For n>=9 both grandchildren are nonexceptional partitions of n-2: shorten the last two rows, or remove the last two boxes of the last row. Induction gives each degree greater than n-1. Branching twice gives parent degree at least 2n, again greater than n+1. This proves (6). The irreducible dimension formula and branching rule are classical; their [primary source](https://www.math.ucla.edu/~pak/papers/whl.pdf) is pinned.

For n=2,...,6 the needed complete lists follow directly from the hook formula. In particular:

- S2 has the two characters. S3 has the standard2 plus the characters.
- S4 has characters, the2-dimensional representation through S3, and the two faithful standard3 representations. The2 kills V4.
- S5 has, besides characters, dimensions4,4,5,5,6. On C5 their fixed dimensions are0,0,1,1,2. For a5-dimensional representation use the permutation character on2-subsets minus the point permutation character; on C5 this gives2-1=1. The6-dimensional one is wedge^2 V5, whose fixed dimension is2 by pairing reciprocal nontrivial C5 characters. Sign twists agree on C5.
- S6 has four5-dimensional irreducibles and no other non-character irreducible of dimension at most7. On H=C3 times C3, acting on two3-blocks, all four have fixed dimension1. For the standard5 it is2-1; for the other5 it is the difference between the numbers of H orbits on3-subsets and2-subsets, namely4-3. Sign twists agree on H.

The subset-character decompositions are the usual multiplicity-free decompositions of the symmetric group's subset permutation representations. They can also be checked directly from the small character tables. No claim is made that arbitrary modular representations satisfy (6).

## 5. UANONLINEARCOST1: dimension costs without linearizability

For the natural marked S_N fibre of rank N>=2 in characteristic zero, with arbitrary polynomial actions, define
\[
 c(N)=\begin{cases}
 N+1,&N=\ell^a\ (\ell\text{ prime},a\ge1),\\
 N,&N=\ell^a+\ell^b\ (a,b\ge1),\text{ and the first case fails},\\
 N-1,&\text{otherwise}.
 \end{cases}                                                   \tag{7}
\]
Then
\[
                         \boxed{D\ge c(N).}         \tag{8}
\]
In particular **no characteristic-zero canonical marked S_N construction exists on A3 for any N>=3**, even with genuinely nonlinearizable actions. A generic degree larger than N does not evade the bound.

### Faithful tangent representation

The label action is faithful on X. UADECKWILD1 makes W=T_yY faithful. The smallest faithful characteristic-zero S_N representation has dimension N-1: (6) proves this for N>=7, and the displayed small lists prove it for N<=6, including S4's nonfaithful2. This proves the baseline N-1.

### Prime-power rank

For N=ell^a choose the transitive regular elementary abelian ell-group H on the labels. Apart from N=2 it lies in A_N: odd-order permutations have positive sign, and for N=2^a>=4 a nonzero translation has N/2 transpositions, an even number. It has no fixed label, so (5) needs dim W^H>=2.

For N>=7, if D<=N, (6) forces one standard or sign-standard summand, and at most one character. The standard H-fixed dimension is zero because H is transitive; sign is trivial on H. Hence dim W^H<=1, a contradiction. For N=3,4,5 the small lists give the same conclusion; the5-dimensional S5 exceptions have fixed dimension one and are thereby excluded. For N=2, faithfulness requires a sign summand, and D<=2 leaves at most one trivial direction. This again contradicts (5). Thus D>=N+1 in every prime-power rank.

### Two positive prime-power orbits

Suppose N=ell^a+ell^b with a,b>=1. Choose an ell-subgroup H contained in A_N, transitive on each block. For odd ell, the product of the two cyclic block groups works. For ell=2, take the subgroup of the two cyclic block groups with equal parity; both projections are surjective, so both orbits stay transitive. This explicitly handles blocks of size two as well.

There are exactly two orbits, no fixed labels, and dim V_N^H=1; sign twists agree. For N>=7, D<=N-1 would force W to be just one standard or sign-standard summand, contradicting (5). N=4 is already prime-power. For N=6 use H=C3 times C3 and the four5-dimensional exceptions from Section4; each again has fixed dimension one. Hence D>=N. The argument does not require a conjecture about representing every even integer by prime powers.

### What equality would require: actual acyclic surfaces

In the prime-power case D=N+1, the same H has
\[
                         \dim W^H=2.               \tag{9}
\]
For N>=7 use (6), now up to N+1. The target has one standard summand and two characters. For N=2,3,4,5 the small lists and (5) give the same conclusion: the possible faithful degree-N+1 representations surviving (5) all have H-fixed dimension exactly two. For instance two standard S3 summands would have fixed dimension zero and are excluded; S4's extra2 contributes two fixed directions; S5's wedge-square6 has two.

If the second case of (7) applies and D=N, the constructed two-orbit H also has fixed dimension two. For N>=7 the target is standard or sign-standard plus one character; N=6 uses5 plus one character.

Consequently every candidate attaining either strengthened lower bound must induce
\[
       F^H:X^H\longrightarrow Y^H                 \tag{10}
\]
as an etale morphism of **smooth connected mod-ell acyclic affine surfaces**, missing the marked target. The whole surfaces are not automatically A2. If both are A2, (10) is a plane Keller counterexample; otherwise their actual non-plane geometry is a new obligation. Mod-ell acyclicity is not contractibility or an affine-plane recognition theorem.

This relocates the equality obstruction to precise surface geometry and does not settle it. The published nonlinearizable S3 action with full fixed locus A1 is still excluded by UAFIXGEOM1; passing its C3 surface test alone would not repair that separate full-group fixed-line failure.

### A strict fixed-versus-canonical gap in every rank

The retained explicit fixed universal **absolutely stably atomic** maps have dimensions
\[
 u_N=\begin{cases}
 N,&N=3,4\quad(\mathrm{UAM1}),\\
 4,&N=5\quad(\mathrm{UACOMP1}),\\
 4,&N=6\quad(\mathrm{UASEXT1}),\\
 N-2,&N\ge7\quad(\mathrm{UATRACE1}).
 \end{cases}                                                   \tag{11}
\]
These are available over every characteristic-zero base and realize every rank-N etale algebra over every extension. They are upper bounds, not asserted minima. The authorities are the [original absolute theorem, Section3.1](../../foundations/verified/UNIVERSAL_RELATIVE_KELLER_MAP.md), [compression](../universal-keller-compression/PROOF.md), [sextic](../universal-keller-sextic/PROOF.md), and [trace](../universal-keller-trace/PROOF.md).

For every N>=3, (7) gives c(N)>u_N. Thus **every actual canonical natural marked affine-space realization has strictly larger dimension than the minimum of fixed universal maps**, even allowing arbitrary nonlinear polynomial actions and without demanding canonical atomicity. Quantitatively, D-u_N is at least one in ranks3,4; two in ranks5,6; and, for N>=7, three in prime-power rank, two in the second arithmetic case, and one otherwise.

In particular no dimension-minimizing fixed universal map can acquire the required canonical marked symmetry in its original dimension. The explicit map of dimension u_N would need at least c(N)-u_N identity coordinates before such a marked action could exist; this is a necessary cost and not a proof that stabilization ever suffices. The fixed maps already being absolutely stably atomic makes this separation independent of decomposability.

[UAK7MIN1](../universal-keller-quadrics/PROOF.md) gives ktdim_k(S7)=4 in characteristic zero when -14 is a square in k. In exactly the same bases, (8) says any canonical natural marked affine-space construction, with arbitrary polynomial actions, needs dimension at least8. Thus the marked canonical cost exceeds the unmarked fixed-map minimum by **at least four coordinates**, with no existence assertion in dimension8. The bound survives identity stabilization; m original coordinates need at least c(N)-m added coordinates whenever the resulting marked actions exist.

## 6. Exact limitations and hypothesis controls

For N>=7 outside both arithmetic cases in (7), the ordinary standard representation V_N of dimension N-1 passes **every** prime-subgroup test (5). To verify this without enumerating subgroups, any ell-group orbit has prime-power size. A failure of (5) would require at most two orbits and a fixed-label count different from one. One orbit gives N=ell^a; two nontrivial orbits give N=ell^a+ell^b. An orbit of size one plus a nontrivial orbit has exactly one fixed label, so it does not fail. These exhaust the possibilities.

Thus the fixed-tangent criterion itself cannot raise the N-1 bound in the remaining ranks, e.g. N=15 or28. This is a representation-level admissibility statement, **not** an action or Keller-map construction. At D=N-1 a sign-standard target would in fact fail the transposition test: its fixed dimension is one while the marked fixed fibre has N-2>=2 points. Therefore any surviving target tangent module at that dimension must be the ordinary standard one.

In characteristic p, x->x^p-x has regular translation deck group C_p; the group and all nontrivial deck elements are genuinely wild. For F=id on A2, the order-p action (x,y)->(x+y^p,y) is faithful but has identity differential at the origin. This shows why tangent kernels need not satisfy the deck-degree divisibility. The same action has nonreduced fixed scheme y^p=0, so treating a wild fixed locus as its smooth reduced line would invalidate the proof. Finally x->x+x^p has degree p and constant derivative one; its whole line fibres have p points, the permitted characteristic-p exception in (4).

These results neither construct the remaining characteristic-zero canonical receiver nor improve any unmarked ktdim minimum. They do close the general nonlinear three-space route and extend the N+1 cost to every prime-power rank without assuming linearizability, while identifying the exact surface and representation limitations of the method.
