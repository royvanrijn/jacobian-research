# Split target primes in the discriminant -23 shell

This note proves a uniform local-to-global coordinate-gcd estimate
at split primes dividing the target. It does not claim the complete
primitive -23 sumset or an effective all-prime sieve. Counts below
are signed ordered quaternion-coordinate counts, not numerical-value
pair counts; zero binary blocks are included in the bad locus.

Let \(R=\mathbb Z[(1+\sqrt{-23})/2]\),
\(Q=x^2+xy+6y^2\), and
\(O=R\oplus Rj\), \(j^2=-1\), \(ja=\bar a j\).
Its norm is \(Q\oplus Q\). It is a maximal order in the definite
quaternion algebra of discriminant 23: its trace discriminant is
\(23^2\), while the algebra ramifies at 23 and infinity. Its units
are \(\{\pm1,\pm j\}\). The classical mass formula gives total
projective mass \(11/6\), so the principal class has probability
\(w=(1/2)/(11/6)=3/11\). All ideal classes are retained.

Let \(p\ne23\) be a split rational prime in \(R\), including
\(p=2\). At \(p\), choose the integral diagonalization
\(R_p=\mathbb Z_p\times\mathbb Z_p\). Then a quaternion
\(\alpha=a+bj\) has matrix
\[
 \begin{pmatrix}a_+&b_+\\-b_-&a_-\end{pmatrix},\qquad
 \det\alpha=a_+a_-+b_+b_-.
\]
Put \(N=p^tm\), \(t\ge1\), \((m,23p)=1\).
Let \(b_p^{\rm rank1}(N)\) count the signed ordered norm-N
quaternions which are not in \(pO\), but whose first or second
binary coordinate pair is divisible by \(p\).

**Theorem CN3SplitTargetGcd.** Under the definitions and hypotheses
above, for every such split prime, every \(t\ge1\) and every
positive \(m\) prime to \(23p\),
\[
 \left|b_p^{\rm rank1}(N)
      -\frac{16w}{p+1}p^{t-1}\sigma(m)\right|
      \le16t p^{t/2}\tau(m)\sqrt m.                       \tag{S1}
\]
For \(t=1\), the right side can be improved to
\(8\tau(m)\sqrt m\). The main term divided by the primitive
local norm-shell baseline
\(4w(p+1)p^{t-1}\sigma(m)\) is exactly
\(4/(p+1)^2\). These are absolute shell counts, not a conditional
asymptotic asserted for small indices.

## 1. Exact local bad locus and its Haar mass

Modulo \(p\), a nonzero determinant-zero matrix has rank one.
If its two diagonal entries vanish or its two off-diagonal entries
vanish, rank one forces exactly one remaining nonzero entry.
Thus the bad locus is the union of the four nonzero coordinate
axes in \(M_2(\mathbb F_p)\), each of size \(p-1\).
The complete rank-one locus has size
\((p-1)(p+1)^2\). Consequently its bad fraction is
\[
                          4/(p+1)^2.                       \tag{S2}
\]
Every rank-one matrix is nonsingular for the determinant equation
(some cofactor is a unit), so every fixed target divisible by p has
the same lifting proportion among p-primitive matrices. The statement
therefore holds at every target valuation \(t\ge1\).

For comparison, at a unit target the invertible bad locus is the
split-torus normalizer and its fraction is \(2/[p(p+1)]\).
At an inert target-dividing prime, a bad binary block forces both
blocks divisible by p, so the p-primitive bad locus is empty.
This is the exact distinction missing in an inert-only argument.

## 2. Endpoint lattice chains and the positive operator

Let \(\ell_0,\ell_1\) be the two coordinate lines of
\(\mathbb F_p^2\), and fix the edge
\(L_0\supset L_1\supset pL_0\), with
\(L_0=\mathbb Z_p^2\),
\(L_1=\mathbb Z_pe_1+\mathbb Z_ppe_2\).
Its oriented stabilizer is the Iwahori compact \(I_p\).

A p-primitive matrix of determinant valuation t maps \(L_0\)
to an index-\(p^t\) lattice \(M\) with cyclic quotient. Its
rank-one image modulo p is the first edge of the length-t
geodesic from \([L_0]\) to \([M]\). Requiring that image to
be \(\ell_0\) fixes that first edge.

The condition \(\ker(\alpha\bmod p)=\ell_0\) says that its
first column belongs to \(pL_0\). Then
\(\alpha L_1=M\cap pL_0\), the unique index-p sublattice
of \(M\) lying in \(pL_0\). Projectively this is the predecessor
of \([M]\) along the same geodesic. Thus the image oriented edge
is the reversal of the last edge of that geodesic.

Let \(\mathcal B\) be forward nonbacktracking adjacency on
oriented edges of the \((p+1)\)-regular tree, of degree p, and
let \(F\) reverse an edge. The preceding two endpoint conditions
are therefore represented exactly by
\[
              V_t=\mathcal B^{t-1}F,\qquad V_t1=p^{t-1}.
                                                               \tag{S3}
\]
This is a nonnegative-kernel correspondence. Since
\(\mathcal B^*=F\mathcal BF\), \(V_t\) is self-adjoint.
It uses genuine nonbacktracking paths, not an edge-ball degree.

## 3. Full complementary operator bound

Write \(H_0\) for the orthogonal complement of all global norm
characters (here the constants). All operator norms in this section
are restrictions to \(H_0\); the constant eigenvalue
\(p^{t-1}\) of \(V_t\) is retained separately in the main term.
Let \(P_0,P_1\) be the endpoint conditional projections, with
\(P_1=FP_0F\). The forward operator can be written
\[
                 \mathcal B=((p+1)P_1-1)F.                 \tag{S4}
\]
On \(\ker P_0\cap\ker P_1\), it equals \(-F\), so its
every power has norm one. On
\(\operatorname{im}P_0+\operatorname{im}P_1\), diagonalize
\(P_0FP_0=T_p/(p+1)\) at the coarser maximal p-level.
After removing global norm characters, the shared classical Hecke
bound gives \(|\lambda|\le2\sqrt p\).

For an eigenvector \(u\) of this compression, put \(v=Fu\).
The two-frame overlap is \(\lambda/(p+1)\), of modulus less
than one. From (S4),
\(\mathcal Bu=pv\),
\(\mathcal Bv=-u+\lambda v\), so
\[
                \mathcal B^2-\lambda\mathcal B+p=0.
\]
Also \(\|\mathcal B\|\le p\), since (S4) is a projection
combination followed by a unitary. If
\(U_0=1,U_1=\lambda,U_{k+1}=\lambda U_k-pU_{k-1}\), then
\(|U_k|\le(k+1)p^{k/2}\), and for \(k\ge1\)
\[
 \|\mathcal B^k\|
 \le k p^{(k+1)/2}+(k-1)p^{k/2}.
\]
Together with the joint-projection kernel, this proves
\[
 \|V_t|_{H_0}\|\le2t p^{t/2}\quad(t\ge1),\qquad
 \|V_1|_{H_0}\|\le1.
                                                               \tag{S5}
\]
Thus the bound covers the nonspherical sector as well. The compact
unit norms at Iwahori level and at every other place are all local
units; over Q the only global norm character is the constant.
Passing to maximal p-level removes no additional character here.

## 4. The global return factor

Let \(x_A\) be the principal-class point at Iwahori level obtained
by fixing the axis flag \(\ell_0\). The global projective unit
\(j\) exchanges the two axes, so this single point represents
their two-element orbit. Its mass is
\[
                       \nu=2w/(p+1).
\]
Its projective stabilizer is trivial; its actual unit stabilizer
is \(\{\pm1\}\). Let \(f=1_{\{x_A\}}\).

There are four equally large matrix-axis cases in (S2): left
multiplication by j exchanges the image axes and right multiplication
by j exchanges the kernel axes, independently. Fix the case
image=kernel=\(\ell_0\). The lattice correspondence (S3), with
the coprime all-integral-ideal factor \(T_m\), counts exactly
this case divided by the two generators \(\pm1\) fixing the
chosen edge. Hence
\[
           b_p^{\rm rank1}(N)=8(V_tT_m f)(x_A).            \tag{S6}
\]
The factor eight is four axis cases times the two stabilizer units;
it is not the full maximal-order unit count applied twice.
The independent finite Brandt construction below checks this return
normalization against actual signed coordinate counts.

As f is a point indicator of mass nu, self-adjointness gives
\[
 (V_tT_m f)(x_A)=\frac{\langle f,V_tT_m f\rangle}{\nu}.
\]
Its constant term is \(p^{t-1}\sigma(m)\nu\), and the
complementary deviation is at most
\(2t p^{t/2}\tau(m)\sqrt m(1-\nu)\), by (S5) and the
shared unramified bound for \(T_m\). Equation (S6) yields
(S1); using \(V_1=F\) gives the sharper t=1 bound.

This local result would not alone complete the global all-prime sieve.
At split primes dividing the target, the imprimitive-quadruple branch
\(r_O(N/p^2)\) must be handled separately. At small split primes
two and three a joint primitive selection is needed to preserve a
positive total sieve margin. None of those global completion steps
is claimed by (S1).

## 5. Independent maximal-order and Iwahori return check at two

The following finite ideal calculation provides an independent check
of the return normalization, beyond the four-axis symmetry test.
At \(\ell=2,3\), identify \(\omega\bmod\ell\) with the diagonal
roots zero and one. For \(s\in\mathbb F_\ell\), the right ideal
whose reduction has image in the line \((1,s)\) has the basis
\[
 I_{\ell,s}=\langle(1,-1,0,-s),(0,s,1,-1),
                         (0,\ell,0,0),(0,0,0,\ell)\rangle
\]
in \((1,\omega,j,\omega j)\) coordinates. It contains
\(\ell O\), has norm \(\ell\) and is locally principal.
Take \(I_2=I_{2,0}\), \(I_3=I_{3,1}\).

If \(\alpha I_\ell\subset I_\ell\), then
\(\alpha\in\ell^{-1}O\), since \(\ell O\subset I_\ell\subset O\).
Thus all left-order norm-k elements are found by completely enumerating
\(\beta\in O\) of norm \(k\ell^2\), and checking
\((\beta/\ell)I_\ell\subset I_\ell\) on the displayed basis.
Norm one inclusion is equality by covolume, so it counts the full unit
group. The exact enumeration gives:

| Class representative | Left-order units | Left-order norm-two elements |
|---|---:|---:|
| \(O\) | 4 | 4 |
| \(I_3\) | 6 | 0 |
| \(I_2\) | 2 | 2 |

For example the six units of the second row include
\((-1-\omega+j)/3\) and its signed powers. The third row has
only \(\pm1\) as units and its norm-two elements are
\(\pm(1+\omega)j/2\). The finite checker lists every numerator.
Different unit-group orders prove that these ideal classes are distinct;
their projective masses sum to
\(1/2+1/3+1=11/6\), the total maximal-order mass. Hence the
three-class list is complete, with weights \((3,2,6)/11\).

The all-ideal Brandt operator at two has row sum three. Its diagonal
entries are the norm-two counts divided by the unit counts, namely
\((1,0,1)\). Weighted self-adjointness and those row sums determine
the unique nonnegative integer matrix
\[
                       T_2=\begin{pmatrix}1&0&2\\0&0&3\\1&1&1\end{pmatrix}.
\]
The Iwahori oriented-edge classes can then be labelled
\[
 (0,0),(0,2),(1,2),(2,0),(2,1),(2,2),
\]
with weights \((1,2,2,2,2,2)/11\). The principal axis point is
\((0,2)\), of weight \(2/11=2w/3\). The principal projective
order-two unit exchanges the two axis flags and fixes the third;
the projective order-three unit group in class \(I_3\) acts
transitively on the three flags (its reduction is nontrivial because
its trace is odd). The last class has trivial projective unit group.
These orbit counts and the coarse matrix give the full edge list.

Edge reversal fixes \((0,0),(2,2)\) and swaps
\((0,2)\leftrightarrow(2,0)\),
\((1,2)\leftrightarrow(2,1)\). Forward nonbacktracking adjacency
uses the coarse outgoing multiplicities at the terminal vertex,
minus one occurrence of the reverse edge. This uniquely determines
the six-by-six matrices \(F,\mathcal B\) listed in the receipt.

Their returns \(8(\mathcal B^{t-1}F)_{(0,2),(0,2)}\), for
\(1\le t\le12\), are
\[
             0,0,8,16,32,24,80,224,392,720,1440,2904.
\]
Complete signed coordinate enumeration in the independent Q-norm
shells \(N=2^t\) gives exactly the same rank-one bad counts.
This validates the actual Brandt return factor (S6) through 4096
using a complete small ideal-class computation, rather than only
checking divisibility of coordinate counts by eight. It remains a
finite check supporting the unrestricted written correspondence.

## Dependencies and verification scope

The infinite spectral input is the same
[classical all-ideal Hecke transfer](../../results/eisenstein-restricted-prime-sums/sections/hecke-transfer.tex)
used by the exact companion. It is applied to \(T_m\) away from
\(23p\), and to \(T_p\) only on the coarser maximal p-level.
The nonmaximal Iwahori operator is bounded separately by (S3)--(S5),
including the joint-projection kernel. There is no unramified
assertion made at a prime where the level is Iwahori.

Primary structural references are Voight, *Quaternion Algebras*,
author edition 1.0.7u, Theorem 25.1.1 for maximal projective mass,
Section 23.5 for lattice chains and the regular tree, and
Propositions 41.3.1/41.3.6 for the coprime all-integral-ideal
correspondences. Their statements were inspected in the retained
primary edition; the deep Jacquet--Langlands and Deligne inputs
remain imported through the pinned shared proof.

The [exact checker](verify_primitive_norm_minus23_split_targets.py)
and [receipt](primitive_norm_minus23_split_target_checks.json)
separate local matrix masses, higher-digit lifts, primitive HNF
degrees, axis symmetry, complete small left-order calculations,
the six-state p=2 Iwahori operator, and its actual return checks.
The written infinite correspondence and spectral estimate are not
proved by running those finite tests. Internal mathematical review
and the exact computations do not establish external human review,
formal verification or historical priority.
