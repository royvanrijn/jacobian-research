# Curve302: the strict ten-dimensional block has no unit kernel

The displayed strict block has **exactly ten independent ordinary half-ideal
classes and zero unit-squareclass kernel**. The same ten ideals give an
elementary direct factor in the narrow class group. This closes the
unit-versus-ideal question for this specified block, without computing
fundamental units or a complete class group.

The full Selmer boundary also closes:

\[
\boxed{\dim\operatorname{Sel}_2(E_{302})=21+c_S
       =31+(c_S-10),\qquad c_S\ge10,}
\]

where \(c_S=\dim\operatorname{Cl}(\mathcal O_{K,S_K})/2\) remains **UNKNOWN**.
This is an exact identity with an unknown class-group term, not a numerical
Selmer upper bound or an exact-rank result.

The complete four-experiment objective remains open. These are retrospective
theorems on the certified rank31 subgroup, not evidence that class-group
features predict a new large jump.

## The subgroup and cubic marking

Use the literal public rank31 group \(D\), the specialized recovered MW17
group \(M\), and their cubic Kummer images \(W,G\). The
[local-filtration certificate](CURVE302_RECOVERED_QUOTIENT_LOCAL_FILTRATION_2026-09-07.md)
proves \(\dim W=31\), \(\dim G=17\), and

\[
\dim\operatorname{loc}_S W=21,\quad
\dim\operatorname{loc}_S G=17,\quad
V=\ker(\operatorname{loc}_S|W)\cong\mathbf F_2^{10}.
\]

Here \(S\) contains infinity, 2 and all twenty retained finite bad/model
primes. In particular \(G\cap V=0\), and the existing quotient filtration is

\[
0\longrightarrow V\longrightarrow W/G
 \longrightarrow\mathbf F_2^4\longrightarrow0.
\]

The arithmetic marking is \(X=4x\), \(Y=8y+4x+4\), with
\(Y^2=f(X)\), \(K=\mathbf Q(\theta)\), and

```text
f(z) = z^3 + 5*z^2
 - 20555644225817083652508762176252556301867322902653588194515587713000*z
 + 35863572573072725712438997689722579989297270955836506038914807628765229483647299767251594888128395600.
```

The totally real cubic uses the certified maximal order, not its power
order of index `9923258571576383661788160`. All ten words are the exact
`strict_kernel_public_words` of the earlier filtration, in their original
order. No new representative selection uses Artin outcomes.

## Ten exact half ideals and a rank-nine first pairing

For each public point write \(X_i=a_i/d_i^2\), \(Y_i=b_i/d_i^3\), and put
\(\gamma_i=a_i-d_i^2\theta\). For the frozen binary word \(w_j\), set

\[
\beta_j=\prod_{i:w_{ji}=1}\gamma_i,
\qquad (\beta_j)=J_j^2,\qquad j=0,\ldots,9.
\]

The \(J_j\) are explicit integral ideals in HNF, with positive generators
\(\beta_j\). Construction uses gcd ideals \((b_i,\gamma_i)\), corrections
at the already known bad prime ideals, and exact final square identities.
It does not factor the point-coordinate numerators.

Each strict class also gives a quadratic character
\(\chi_j\) of the ordinary \(S\)-class group. These are different objects:
\(\chi_j\) is an unramified extension character; \([J_j]\) is ideal-class
2-torsion. Form

\[
A=(\chi_i(J_j))_{0\le i,j<10}.
\]

All 100 entries are exact and **rank \(A=9\)**. Its right kernel is

\[
\ker A=\langle(1,0,1,1,0,1,0,0,0,0)\rangle.
\]

Reduction of the ten ideals gives cyclic good residue rings
\(\mathcal O_K/I\simeq\mathbf Z/N\), so integer Jacobi symbols evaluate
95 entries without norm factorization. The five nonunit exceptions are
resolved only at 47, 67 and 89 by local normalization. The independent
checker uses exact Hensel lifts through valuation two, rather than the
producer's local-power routine, to check those bits.

This already supplies a split \((\mathbf Z/2)^9\) factor of the ordinary
\(S\)-class group. It does not show that its full 2-rank is nine:
the ten independent characters themselves prove \(c_S\ge10\).

## The final ideal is detected by inherited characters

Set \(J_*=J_0J_2J_3J_5\). It is invisible to the ten strict characters.
Test it against the **six previously certified ordinary unramified
characters from the generic MW17 subgroup**, with generic masks

```text
46473, 81930, 62788, 59728, 85344, 117248.
```

The result is exactly

\[
(\eta_1(J_*),\ldots,\eta_6(J_*))=(1,1,0,0,0,0).
\]

These characters need not be split over \(S\). Accordingly the calculation
restores the Artin contribution of every removed bad prime ideal; for
\(J_*\) the only one is the retained degree-one ideal above 7, with exponent
one. All six of its contributions are zero, checked separately by a
simple-root Hensel calculation. Omitting such contributions without this
check would be invalid.

Consider the map from the ten-dimensional formal ideal space to
\(\operatorname{Cl}(K)/2\). Any kernel vector is killed by the ten
\(\chi_i\), hence is zero or the displayed word for \(J_*\). The latter
is detected by \(\eta_1\). Thus all ten images in \(\operatorname{Cl}(K)/2\)
are independent. Since \(J_j^2=(\beta_j)\), they have order two and generate
a **split elementary direct factor \((\mathbf Z/2)^{10}\)** of
\(\operatorname{Cl}(K)\).

The half-ideal homomorphism on \(V\) therefore has rank ten and kernel zero.
Its kernel is exactly the intersection of \(V\) with ordinary unit
squareclasses. Consequently

\[
\boxed{V\cap\bigl(\mathcal O_K^\times/\mathcal O_K^{\times2}\bigr)=0.}
\]

Because every \(\beta_j\) is totally positive, the same ideals have order
two in \(\operatorname{Cl}^+(K)\). Projection to the ordinary group
preserves their independence; the above characters give a retraction there
as well. This proves the same split ten-dimensional factor in the narrow
group. It computes neither the full narrow group nor its signature kernel.

The ten strict characters and the six inherited ordinary unramified
characters are independent: their intersection is contained in
\(G\cap V=0\). In particular
\(\dim\operatorname{Cl}(K)/2\ge16\). This is a lower bound from known
characters, not a completed class-group measurement.

## No remaining Selmer condition inside the displayed block

All ten classes already lie in the rational Kummer image, so their covers
have rational points and zero Sha image. Unit, ideal and local/Selmer
pieces are not three unrelated summands to add together: the local
conditions filter the squareclass space, and the half-ideal map then has a
unit kernel. On this \(V\), that kernel is zero and the ideal image has
dimension ten. This says nothing about unseen Selmer classes.

To determine the remaining boundary, apply the
[strict class-field identification](../rank-jump/STRICT_SELMER_AND_ARTIN_BLOCKS.md)
and [derivative reciprocity argument](../rank-jump/DERIVATIVE_RECIPROCITY_AND_COMPLETE_BOUNDARY.md).
They are applications of classical cubic descent and global class field
theory; see [Barrera Salazar–Pacetti–Tornaría, §2](https://arxiv.org/html/2001.02263v3)
and [Milne, Class Field Theory, Chapter V](https://www.jmilne.org/math/CourseNotes/CFT.pdf).

The strict Selmer group is the full ordinary \(S\)-class character group:
\(\operatorname{Sel}_2^S(E)\simeq\operatorname{Hom}(\operatorname{Cl}_S(K),\mathbf F_2)\).
Let \(\delta=\operatorname{disc}(f)>0\) and
\(\alpha=-\delta f'(\theta)\). Its norm is \(\delta^4\); it is unramified
outside \(S\). At the ordered real roots its sign bits are \((1,0,1)\),
outside the local point images \((0,0,0),(0,1,1)\).
Local self-duality and global reciprocity therefore supply a nonzero
linear constraint on the 22-dimensional product of local point images.
The known 21-dimensional joint image already fills that hyperplane.
Hence the full Selmer localization image has dimension exactly 21.

Taking the kernel and image proves the opening formula. Equivalently,

\[
\dim(\operatorname{Sel}_2(E)/G)=14+(c_S-10),\qquad
\operatorname{rank}E+\dim\Sha(E)[2]=31+(c_S-10).
\]

Thus all unaccounted Selmer dimensions are additional strict characters.
The statement \(c_S=10\), if independently proved, would give exact rank31
and \(\Sha(E)[2]=0\). It is **not** proved here.

## Replay, scope and next decision

The [complete matrix packet](../../artifacts/generated-results/elliptic-curves/curve302_descent_anatomy_v1.json)
embeds all half ideals and the initial incomplete Artin calculation. The
[remaining-class packet](../../artifacts/generated-results/elliptic-curves/curve302_descent_remaining_class_v1.json)
retains the second character test. The
[manifest](../../artifacts/generated-results/elliptic-curves/curve302_descent_anatomy_manifest_v1.json)
pins producers, checker, helper implementations and inputs.

From the repository root:

```sh
timeout 120 sage -python research/elliptic-curves/cas/verify_curve302_descent_anatomy.sage --negative-controls
```

Replay checks all ten ideal squares by integral-lattice containment and
norm equality, all reduction transports by unimodular lattice changes,
all cyclic ring multiplication identities, all Artin entries, the final
kernel witness, positivity, complete bad support, and the real derivative
obstruction. It also reruns the existing local-filtration and ordinary
unramified-character verifiers. Two altered half ideals, including a
determinant-preserving alteration, are rejected. Sage10.9/PARI2.17.3 was
used. Each producer and the independent replay finished below one second
on this host; this is not a runtime prediction for other fibres.

There is no BNF/unit-group call, parameter search, point search, or
production-search change. Initial incomplete entries remain preserved.
The checker needed a coefficient-list coercion correction during development;
the mathematical inputs and producer outputs did not change.

The four requested experiments now stand as follows:

| Experiment | Established here | Remaining gate |
|---|---|---|
| Descent anatomy and matched controls | Exact strict unit kernel0, ideal image10; full Selmer boundary21 | Full inherited everywhere-even ideal image and its intersection with this block; matched rank-lower-bound17/21/26/27 fibres; independent class/Selmer upper bounds |
| Class1↔class6 transport | Existing maps and the common lattice are identified as distinct inputs | Exact transport of the fourteen points, preserving each new base parameter; an orbit/section comparison |
| Shared carriers | No new carrier assertion | Construct actual generic carriers through directions, then compare squareclasses and joint normalization genera |
| Self-propagation | The 13/14 finite recovery result remains intact | A point-producing implication with certified independence, beyond recovery on one fibre |

**Next priority:** finish the inherited ideal-image comparison before
replicating the anatomy on matched fibres. The new ten-dimensional image
is independent modulo the generic **strict** image (which is zero), but
may intersect the larger image of generic classes that have even finite
valuations without being strict. This distinction is needed before calling
ten dimensions new class-group contributions or a predictor.

For the geometric experiments, a point on the original302 fibre is a point
of the surface, not a generic A-section. Its B-coordinate can vary among
the fourteen points; the birational map is not a group homomorphism between
one pair of fixed elliptic curves. Likewise a specialized half-lattice
chart is not automatically a quadratic cover of the parent parameter.
These type checks are prerequisites to a valid orbit or shared-carrier test.
