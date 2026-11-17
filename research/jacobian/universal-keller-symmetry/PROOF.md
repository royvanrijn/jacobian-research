# Canonical line classification, affine cyclic covers and nonlinear fixed-locus gates

14 November 2026. This extends the [cost theorem](../universal-keller-costs/PROOF.md). The results classify an entire canonical category, construct fixed framed universal maps for every cyclic-complement affine group, decide atomicity throughout that family, and transfer the obstruction to nonlinear actions using actual fixed loci and published linearization results.

## 1. Every canonical affine-line realization

A marked realization for a finite group G with a faithful transitive N-point set Omega means a G-equivariant polynomial Keller map between affine lines, a G-fixed target, and an identification of its complete fibre with Omega. This is the split-object description of an affine-line-bundle construction over BG. Coordinates and normal forms below are geometric; no canonical coordinate or origin is being selected over BG.

**UA1CANCLASS1.** Let N>=2. Over an algebraic closure, such a realization exists exactly for the permutation groups
\[
 G=({\mathbb F}_q,+)\rtimes C_m\ \subseteq\operatorname{AGL}_1({\mathbb F}_q),
 \qquad q=p^e=N,\quad m\mid q-1,                    \tag{1}
\]
with the natural affine action, in characteristic p. Its polynomial, after affine coordinate changes, is a nonzero scalar multiple of
\[
 P_W(X)=\prod_{w\in W}(X-w),                        \tag{2}
\]
where W is a finite additive subgroup of size q and C_m acts by scalar multiplication preserving W. Geometric monodromy is regular W. Absolute stable atomicity holds exactly when q is prime.

For the natural symmetric group S_N, a canonical marked affine-line construction exists exactly in the two cases
\[
 (N,\operatorname{char} k)=(2,2),\ (3,3).           \tag{3}
\]
These exceptions work over every field of the stated characteristic, using X-X^N. The classification does not say that a single fixed univariate map is field-versal for S3: its coefficients after trivializing twisted bundles depend on the input algebra.

### Finite affine groups and the derivative character

Every polynomial automorphism of a line is affine. The source multiplier character has finite cyclic image C_m, of order prime to p, and translation kernel W. In positive characteristic W is an F_p vector space of size q=p^e. In characteristic zero it is trivial.

If m>1, an affine element with primitive multiplier a has a fixed point; center there to make it X->aX. It has order m and splits the multiplier map. Its action on W minus zero is free, so m divides q-1 when W is nontrivial. If m=1 the group is W itself. These are elementary facts about the affine group, not a new group-classification claim.

Translate the marked target to zero. Its action is linear. Differentiating f(gX)=g f(X), using f'=c nonzero constant, identifies the target character with the source multiplier. Consequently W acts trivially on the target and f is W-invariant.

The polynomial P_W is additive with nonzero constant derivative. Induct on dim_Fp W: if W=W0 plus F_p*a, then
\[
 P_W=P_{W0}^p-P_{W0}(a)^{p-1}P_{W0}.               \tag{4}
\]
The coefficient P_{W0}(a) is nonzero. Also k[X]^W=k[P_W]. For example, the field extension k(X)/k(P_W) is separable of degree q with the q translations as its full group; its invariant field is k(P_W). An invariant polynomial is a polynomial in P_W: a rational function with a finite pole would still have a pole after pulling back along the surjective polynomial P_W. Thus
\[
 f=R(P_W),\qquad R'\ \text{a nonzero constant}.      \tag{5}
\]
For W=0 the same formula uses P_W=X and q=1.

### Transitivity rules out every larger orbit

One translation orbit is a fibre of P_W, of size q. The multiplier a acts on its quotient coordinate by P_W(aX)=a^q P_W(X)=a P_W(X); the equality a^q=a follows from m dividing q-1, and is also true for q=1.

The roots of R form one C_m orbit, because the complete marked fibre is transitive. That orbit has size one if it is zero, or size m otherwise. Hence deg R is one or m. But a polynomial with constant nonzero derivative in characteristic p has degree one or degree divisible by p; m>1 is prime to p. In characteristic zero such a polynomial is already linear. Therefore R is linear in every possible case, N=q, and W is nontrivial. This proves the prime-power rank statement and (2).

To identify the permutation group abstractly, put d=[F_p(a):F_p], with d=1 for m=1. Since W is stable under a, it is a vector space over F_p(a)=F_(p^d). Therefore d divides e; an F_(p^d)-linear identification W with F_q identifies the action with (1). This relabeling need not be an affine coordinate change of the original W inside k. Different subspace polynomials can remain.

Over a specified field k the field-of-definition criterion is also exact. Write d=ord_m(p), with d=1 for m=1. The group in (1) admits a canonical marked affine-line realization over k if and only if k contains a primitive m-th root of unity and dim_(F_(p^d)) k>=e/d. Necessity follows from its actual scalar multiplier and translation subgroup W inside k. For sufficiency choose an e/d-dimensional F_(p^d) subspace W of k and use P_W. Thus over an infinite characteristic-p field the root-of-unity condition suffices; over F_(p^f) the criterion is m dividing p^f-1 and e<=f. The hypothesis k containing F_q is a convenient sufficient condition, not a necessary one.

Conversely X-X^q is equivariant for every group in (1), from X->aX+b on the source to t->at on the target, and its zero fibre is the natural F_q set. It descends on affine line bundles by the same construction as UAFROBCOST1.

For S_N the point stabilizer would be cyclic of order (N-1)!, bounded by N-1. Thus N<=3; the rank must equal a power of the characteristic, giving exactly (3).

For composite q, factoring P_W through any subgroup F_p*a gives two nontrivial polynomial Keller factors, as in UAFROBCOST1. For prime q the regular degree-p geometric action is primitive; the retained primitive-monodromy argument gives absolute stable atomicity, including identity stabilizations. This proves the atomicity assertion without assuming that every imprimitive cover admits an affine-space factorization.

For every group in (1), over a base containing F_q, the canonical affine-bundle minimum is therefore one. With absolute stable atomicity it is one for prime q and two for composite q: restrict the full H_q-equivariant plane construction UAFROB1 to G. Its marked zero fibre has the required G action, while its generic geometric group remains H_q. A canonical marked family does not require its generic frame group to equal the group acting on the one specified fibre. The fixed G-framed requirement below does impose that extra condition.

## 2. One fixed universal map for every cyclic affine complement

Assume k contains F_q, q=p^e, and m divides q-1. Put
\[
 G_{q,m}=F_q\rtimes\mu_m,\qquad r=(q-1)/m.
\]
Use the G-framed, open finite-chart definition from UAFINITE1: a specified G torsor on a nonempty target open must frame the entire rank-q cover and realize every G torsor over every extension field, including disconnected torsors and finite fields.

**UAAFFSUB1.** The single plane map
\[
 \boxed{K_{q,m}(u,x)=(u,x-u^r x^q)}                 \tag{6}
\]
has Jacobian one, arithmetic and geometric group G_(q,m), and is field-versal in this framed sense. Its minimum possible fixed framed Keller dimension is one when m=1 and exactly two when m>1.

Let d be the multiplicative order of p modulo m, with d=1 when m=1. Then
\[
 \boxed{K_{q,m}\text{ is absolutely stably atomic}
             \quad\Longleftrightarrow\quad d=e.}  \tag{7}
\]
If d<e the decomposition is global, polynomial, Keller and already defined over k; it persists after identity stabilization. This is an exact classification of this family, not a claim that generic monodromy classifies atomicity for all Keller maps.

### The specified torsor, not just the splitting field

On u!=0 define the finite etale rank-mq scheme of frames by
\[
 \lambda^m=1/u,\qquad \mu-u^r\mu^q=v.              \tag{8}
\]
Its labels are the complete q roots mu+lambda*s for s in F_q. Changing a frame acts by
\[
 (\lambda,\mu)\longmapsto(a\lambda,\mu+b\lambda),
 \qquad a\in\mu_m,\quad b\in F_q.                 \tag{9}
\]
The relation u^r*lambda^(q-1)=1 proves the second equation is preserved. This action is simply transitive on geometric frames, so (8) is an actual G_(q,m) torsor, including its chosen subgroup reduction.

For arbitrary L/k and a continuous G cocycle (a_gamma,b_gamma), multiplicative Hilbert90 supplies nonzero lambda with gamma(lambda)=a_gamma*lambda. The quantities b_gamma*lambda form an additive cocycle, so additive Hilbert90 supplies mu with gamma(mu)-mu=b_gamma*lambda. Then
\[
 u=\lambda^{-m}\in L^*,\qquad
 v=\mu-u^r\mu^q\in L.                              \tag{10}
\]
The frame at this target has the prescribed cocycle. This proves versality for every torsor, not only connected extensions or a dense set over one infinite field.

Generically lambda^m=1/u has degree m by the u-adic valuation. Over k(lambda), putting mu=lambda*z turns the remaining equation into z-z^q=v/lambda. This is a separable degree-q rational function extension with translation group F_q. The splitting field contains lambda as a difference of roots, and its full group is exactly G_(q,m). The same argument over an algebraic closure of k proves the geometric assertion.

Every one-variable Keller map is globally finite. For m>1, the cyclic quotient C_m is excluded by UAFINITE1, so dimension at least two is necessary. For m=1 the single map X-X^q realizes every translation torsor by additive Hilbert90, and dimension one is attained.

### Primitive scalar action or an explicit polynomial decomposition

Blocks for a permutation group containing the regular translation group F_q correspond to additive subgroups stable under its point stabilizer. Here they are precisely F_(p^d)-linear subspaces of F_q. Thus the action is primitive exactly when d=e. In that case UAALLATOM1's primitive-monodromy argument proves absolute stable atomicity.

If d<e, put B=p^d, s=e/d>1 and a=(B-1)/m. Define
\[
 h(u,x)=x-u^a x^B,\qquad
 Q(u,z)=\sum_{i=0}^{s-1}u^{a(1+B+\cdots+B^{i-1})}z^{B^i},       \tag{11}
\]
where the empty exponent sum is zero. Frobenius additivity telescopes:
\[
 Q(u,h(u,x))=x-u^{a(1+B+\cdots+B^{s-1})}x^{B^s}
            =x-u^r x^q.                            \tag{12}
\]
Both (u,h) and (u,Q) have Jacobian one; their generic degrees are B and B^(s-1), both greater than one. This proves the converse of (7) by an actual decomposition.

For example q=16,m=5 is already absolutely stably atomic although its complement is smaller than F16*, whereas q=16,m=3 factors through x-u*x^4. Full affine-group monodromy is sufficient but is not necessary within this family.

### Essential dimension does not give the Keller minimum

The finite affine group acts generically freely on the line. Twisting that line by any G torsor gives an affine line again, by the two Hilbert90 arguments. Over an infinite extension its free open has a rational point. Thus this one-dimensional generically free G variety is versal; the classical essential dimension is one. This agrees with Kemper--Mattig's published one-parameter generic polynomials for these groups.

Consequently, over every infinite base containing F_q and m>1,
\[
 \operatorname{ed}_k(G_{q,m})=1
       <\text{fixed framed Keller dimension}=2.    \tag{13}
\]
The extra cost comes from requiring a whole affine-space source and a constant Jacobian. A generic polynomial can have parameter-dependent derivative or omit source points. Its one parameter is not a one-dimensional Keller self-map.

For (6), the nonproper locus is exactly u=0: rank q off it and rank one on it. When m>1 the varying quotient lambda^m=1/u cannot extend across the divisor. When m=1 this particular plane still has a boundary, but the one-variable translation map is finite everywhere. The geometric source of the necessary boundary is therefore distinguished from an artifact of a chosen plane formula.

## 3. Nonlinear actions must pass their actual fixed-locus geometry

**UAFIXGEOM1.** Let a finite group act by polynomial automorphisms on two affine spaces and let F be an equivariant Keller map. If the complete fixed subschemes are each isomorphic to A^t, the restriction is a Keller map A^t->A^t, in every characteristic. A marked transitive fibre of rank greater than one has empty fixed fibre, so t>=2; in characteristic zero t=2 would contradict JC(2).

Polynomially conjugating the source and target actions preserves Keller maps and their complete marked fibres. Thus UACANCOST1 applies whenever both S_N actions are polynomially linearizable, not merely presented linearly. More generally, if adding r trivial identity coordinates makes both actions linearizable, then
\[
 d+r\ge N+1.                                       \tag{14}
\]
In characteristic zero the conditional JC(2) bound is N+2; the unconditional torus bound N+2 applies if the conjugacies simultaneously linearize a retained nontrivial commuting torus. Separate linearizations of the finite group and the torus do not suffice for that stronger assertion.

### Proof of the geometric fixed-locus restriction

At a fixed source point, DF and its inverse commute with the group, so they induce inverse maps on invariant tangent vectors. The tangent space to the fixed subscheme is exactly the invariant tangent space. The stipulated fixed schemes are smooth affine spaces, so the restricted map has invertible differential everywhere and hence constant nonzero Jacobian. This uses no averaging. The fixed fibre is precisely the fixed part of the complete marked fibre, and a nontrivial transitive set has no full-group-fixed label. Surjectivity of a univariate Keller map over an algebraic closure gives t>=2. These hypotheses concern the whole fixed schemes; a disconnected or nonaffine fixed locus is not silently replaced by one tangent space.

For the conjugacy assertion, any polynomial automorphism has constant nonzero determinant; source and target coordinate changes preserve the Keller condition. Identity stabilization preserves the fibre at added zero coordinates. Apply UACANCOST1 to the resulting linear actions. This proves (14).

### Published linearization theorems close two characteristic-zero routes

Finite polynomial actions on A2 are linearizable. Finite subgroups of the three-variable tame automorphism group are also linearizable (Lamy--Przytycki,2021; Lamy's2026 survey, Propositions2.5 and4.11).

It follows that:

- No characteristic-zero canonical marked plane construction exists for any nontrivial faithful transitive finite G set. After linearization, the source fixed space needs dimension at least two, which makes the whole two-dimensional action trivial, a contradiction.
- For S_N, N>=3, no characteristic-zero three-dimensional construction can have both actions contained in the tame automorphism groups, even after separate polynomial coordinate changes. Such a construction would require 3>=N+1. At least one action must be genuinely nonlinearizable and cannot be conjugated into the tame subgroup.

The imported linearization theorems are established literature. These consequences do not prove that all finite actions in dimension three are linearizable, nor that tame finite actions in every dimension are. The 2026 survey explicitly records nonlinearizable examples in higher dimension.

### A real nonlinearizable candidate still fails the fixed-line gate

Over characteristic zero with a primitive cube root omega, consider the published Freudenburg--Moser-Jauslin S3 action E=A4:
\[
 \sigma(a,b,x,y)=(\omega a,\omega^{-1}b,x,y),
\]
\[
 \tau(a,b,x,y)=
 (b,a,-b^3x+(1+ab+a^2b^2)y,(1-ab)x+a^3y).           \tag{15}
\]
Its algebraic nonlinearizability is published, not proved by our replay. Its actual fixed scheme is elementary:
\[
 E^{S3}=\{a=b=0,\ x=y\}\simeq A1.                 \tag{16}
\]
For source and target carrying this action, no equivariant Keller map can have a natural transitive cubic fibre at a fixed target. The fixed-line restriction is surjective. This excludes the tempting published nonlinearizable candidate without assuming JC(2). The same exclusion holds if one side carries this action and the other side is any linear representation: at a fixed point the differential forces the linear fixed space to have dimension one, and the whole fixed schemes on both sides are affine lines.

More generally, stabilizing both sides by the same representation W adds only dim(W^S3) fixed directions. If that number is zero, the exclusion persists, irrespective of how many moving directions were added. If it is one, the required map would yield a JC(2) counterexample. The original paper proves E becomes linearizable after adding the standard two-dimensional representation, but that representation adds no fixed directions and therefore does not repair this obstruction. Trivial stabilization and representation stabilization must be distinguished.

Other nonlinearizable actions with suitable higher-dimensional or nonaffine fixed loci are not excluded. No characteristic-zero canonical construction is produced here.

## 4. A finite-field bound applying to every universal map

**UAFINCOUNT1.** Let k=F_Q. Any one fixed map A^d->A^d whose complete finite-etale fibres realize every rank-N etale k-algebra must satisfy
\[
 \boxed{Q^d\ge B_N:=\sum_{j=0}^{N-1}p(j),\qquad
 d\ge\lceil\log_Q B_N\rceil,}                      \tag{17}
\]
where p(j) is the partition number, p(0)=1. This includes every field-versal Keller map, in every characteristic and every compiler architecture. It does not require generic S_N monodromy or a globally finite map.

More generally, for a specified field-versal G torsor frame with finite label set Omega,
\[
 Q^d\ge\max\left(c(G),\sum_{[g]\subset G}|\Omega^g|\right),       \tag{18}
\]
where c(G) counts conjugacy classes and the sum has one representative per class, without class-size weights.

### Disjoint whole fibres consume source points

The rank-N etale algebras over F_Q are exactly the products of finite fields with degree multiset a partition lambda of N. Such an algebra has exactly m1(lambda) rational points, the number of parts equal to one. Distinct algebras require distinct target points. Their whole fibres are disjoint subsets of A^d(F_Q), so
\[
 Q^d\ge\sum_{\lambda\vdash N}m_1(\lambda).
\]
Counting a marked occurrence of the part one gives
\[
 \sum_{\lambda\vdash N}m_1(\lambda)
   =\sum_{s=1}^N\#\{\lambda:m_1(\lambda)\ge s\}
   =\sum_{s=1}^Np(N-s)=B_N.
\]
This proves (17). Counting only target types gives Q^d>=p(N), a generally weaker bound. It is indeed implied: remove the smallest part of a partition of N; its remaining partition has size less than N and determines the removed part, giving an injection into the partitions counted by B_N.

For a constant finite G, torsors over F_Q are conjugacy classes of possible Frobenius elements g. A specified frame at one target has one such class, so realizing every class uses at least c(G) target points. Its associated label cover has |Omega^g| rational points. Summing over distinct targets proves (18). The specified frame is essential when different G classes induce isomorphic unstructured label algebras.

For example B_3=4, B_4=7 and B_7=30. Thus every fixed universal S4 map over F2 needs at least three coordinates; every fixed universal S7 map over F2 needs at least five. For rank seven there are only p(7)=15 target types, so target counting alone would have missed the source obstruction in dimension four. No matching upper bound is asserted for these small-field symmetric cases.

There is an elementary all-rank growth rule without an asymptotic partition theorem. If h(h+1)/2<=N-1, every subset of {1,...,h} is a different distinct-parts partition counted by B_N. Hence B_N>=2^h and d>=ceil(h*log_Q(2)). A fixed finite base and ambient dimension cannot support all ranks.

### Exact translation minima over every finite base

For the constant group E=(C_p)^e with its regular p^e-point action over F_(p^f), both its fixed framed universal Keller minimum and its canonical marked affine-bundle minimum are exactly
\[
 \boxed{\lceil e/f\rceil.}                         \tag{19}
\]
The split fibre alone gives the lower bound p^(fd)>=p^e. For the upper bound, decompose the abstract F_p space E into d=ceil(e/f) nonzero subspaces W_i, each of dimension at most f, and embed each W_i additively in F_(p^f). Form their subspace polynomials P_(W_i). The coordinate product
\[
 (x_1,\ldots,x_d)\longmapsto
       (P_{W_1}(x_1),\ldots,P_{W_d}(x_d))           \tag{20}
\]
has constant nonzero Jacobian (rescale outputs for determinant one), is globally finite etale of degree p^e, and has its natural translation E torsor frame everywhere.

For any extension L and any E cocycle, apply additive Hilbert90 to each W_i-valued coordinate. It is gamma(mu_i)-mu_i for some mu_i, and P_(W_i)(mu_i) belongs to L. Thus (20) realizes every framed torsor, over every extension including finite fields. Translation equivariance, with trivial target action, also twists to the canonical marked affine-bundle realization. This proves both minima, without choosing one primitive generator for a split algebra with too many rational factors.

Over every infinite characteristic-p base, both translation minima are instead exactly one: choose an e-dimensional F_p subspace W of k, use P_W in one variable and apply the same Hilbert90 argument. Infinite k has enough F_p-linearly independent constants even when it does not contain F_(p^e). Thus the translation minima are determined over every characteristic-p field, not only bases containing F_q.

These are minima without an atomicity requirement. Nontrivial products and composite-degree subspace polynomials in (20) decompose; no absolute-atomicity minimum is asserted here.

For G_(q,m) over its minimal constant field F_q, the sum in (18) is q+m-1: the identity fixes q labels, nonzero translations fix none, and each of the m-1 nonidentity multipliers supplies one class fixing one label. It gives the same lower bound two for m>1 as UAFINITE1, by an independent finite-field argument.

Finite-base field-versality is stronger than the usual generic-polynomial condition over infinite extensions. These point counts do not transfer to an infinite base merely because it contains a finite subfield, and they give no characteristic-zero dimension bound.

## 5. Remaining gates

The one-variable canonical category is now classified for every faithful transitive finite group; the old Frobenius examples are its full scalar-complement case. Fixed framed affine universality has an explicit two-dimensional map for every complement and an exact all-degree atomicity criterion, with a real polynomial decomposition in every nonprimitive case.

For characteristic-zero canonical construction, a nonlinear formula for an action is insufficient: polynomially linearizable actions inherit the ambient bound, tame actions in three variables are excluded, and the standard nonlinearizable S3 example fails its actual fixed-line test. A viable new receiver must supply both enough fibre symmetry and fixed-locus geometry supporting the required omission.

For d<e, another fixed G_(q,m)-framed universal map with absolute stable atomicity is not ruled out; our explicit family factors. Its existence is a separate gap, requiring an obstruction to polynomial affine-space factors rather than primitive monodromy. The retained [KMFI30 examples](../monodromy-theory/SIGNED_ALTERNATING.md) show in characteristic zero that imprimitive monodromy can coexist with absolute stable atomicity.

Global unmarked characteristic-zero target-dimension minima, general canonical receivers and new affine completions remain open. See the [literature and exact scope](LITERATURE_AND_GAPS.md). The [checker](verify.py) gives bounded exact controls for these written universal proofs, not formalization or independent review.
