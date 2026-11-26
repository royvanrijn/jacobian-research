# Residue-prime Artin conductors and spin classes

Let \(A/K\) have exact real-cubic multiplication as in the
[residue-prime ramification calculation](CM_NEWTON_RESIDUE_PRIME_TORSION_RAMIFICATION.md).
At a good split absolute degree-one place \(v\) of odd characteristic
\(p\), assume the full mod-\(p\) image is
\[
 G=\{(g_i)\in\operatorname{GL}_2(\mathbb F_p)^3:
                        \det g_1=\det g_2=\det g_3\}.
\]
This is the image for all sufficiently large split \(p\) by
[Lombardo, Theorem 1.4](https://arxiv.org/pdf/1411.1703).
Write \(D:G\to\mathbb F_p^\times\) for the common determinant and
\(N=\ker D=\operatorname{SL}_2(\mathbb F_p)^3\).

## Normalized Artin conductors

**Theorem (Uniform normalized residual conductors).** Let \(V\) be a
nonzero complex representation of \(G\) and use the actual action on
\(A[p]\) to define its local Artin conductor exponent \(a_v(V)\).
Then
\[
 \left|\frac{a_v(V)}{\dim V}
       -\left(1-\frac{\dim V^G}{\dim V}\right)\right|
              \leq\frac{6}{p-1}.
\]
In particular, for every intermediate field \(E\) of \(K(A[p])/K\),
\[
 \frac{v(\mathfrak d_{E/K})}{[E:K]}
             =1-\frac1{[E:K]}+O(p^{-1})
\]
with an absolute implied constant uniform in the subfield.

It suffices to bound the normalized conductor on each irreducible
summand. The characters trivial on \(N\) factor through \(D\);
the trivial character has conductor zero, and the other determinant
characters have conductor one. Indeed the residual cyclotomic
determinant surjects on \(\mathbb F_p^\times\) on tame inertia at
\(K_v=\mathbb Q_p\).

For an irreducible \(V\) nontrivial on \(N\), decompose its restriction
to \(N\) into tensor products of irreducibles of its three factors.
Each tensor product has at least one nontrivial factor since
\(V^N\) is a \(G\)-invariant subspace. Let \(T\) be a cyclic complement
to the wild inertia image in the inertia image and \(T^+\) its
index-two subgroup on which the common determinant is square. Such
a complement exists by Schur–Zassenhaus. The projective image of
\(T\) in each factor is a torus of order \(p-1\) or \(p+1\),
depending on the Newton type of that factor. The determinant
squareclass is surjective on either torus. For the nonsplit torus
this follows from the norm \(\mathbb F_{p^2}^\times\to\mathbb F_p^\times\),
which sends a scalar to a square. Thus the image of \(T^+\) in
each factor has order \((p-1)/2\) or \((p+1)/2\).

For \(g\in T^+\), choose \(s\in\mathbb F_p^\times\) with
\(s^2=D(g)\). Then \(g\) is the product of the central common
scalar \(sI\) and an element of \(N\). The scalar acts on \(V\) by
a complex scalar of absolute value one. On a tensor constituent
of \(V|_N\), choose a nontrivial factor. The fraction of \(g\)
that are projectively scalar in that factor is at most \(2/(p-1)\).
For the other \(g\), the component in \(\operatorname{SL}_2(\mathbb F_p)\)
is regular semisimple. Every nontrivial irreducible character there
has degree at least \((p-1)/2\) and has absolute value at most two
on a regular semisimple element. These follow from the two torus
columns of the character table, including the two pairs of
exceptional characters. The normalized character of the whole
tensor then has absolute value at most \(4/(p-1)\).
[Bonnafé, Table 5.4](https://link.springer.com/book/10.1007/978-0-85729-157-8)
· [Braun–Nebe, Theorem 4.1](https://www.math.rwth-aachen.de/~Gabriele.Nebe/papers/sl2.pdf).

Averaging these bounds over \(T^+\) and over the tensor constituents
gives \(\dim V^{T^+}/\dim V\le 6/(p-1)\). The inertia invariants
are contained in these invariants. The Swan term is at most
\(\dim V/(p-1)\) by the sole possible positive break \(1/(p-1)\)
in the [full torsion field](CM_NEWTON_RESIDUE_PRIME_TORSION_RAMIFICATION.md#the-full-torsion-field).
Consequently \(a_v(V)/\dim V\) lies between \(1-6/(p-1)\) and
\(1+1/(p-1)\) for the nontrivial-on-\(N\) irreducibles. The
conductor-discriminant formula and the single trivial summand
of a transitive permutation representation give the field
statement. The order \(p^{-1}\) is necessary for this uniform
estimate. At a mixed place the normalized discriminant exponents
of the [degree-\(2(p+1)\) fields](CM_NEWTON_RESIDUE_PRIME_TORSION_RAMIFICATION.md)
for an ordinary and a nonordinary factor differ by \(1/(p+1)\).
[Serre, Local Fields, Chapter VI §§2–3](https://link.springer.com/book/10.1007/978-1-4757-5673-9).

## Spin classes on the actual residual group

**Proposition (Determinant inflation of spin obstructions).**
Inflation along \(D\) identifies
\[
 H^2(\mathbb F_p^\times,\mathbb F_2)
                \;\xrightarrow{\ \sim\ }\;H^2(G,\mathbb F_2).
\]
Thus the first and second topological Stiefel–Whitney classes of
any real orthogonal representation of \(G\) depend on the residual
cyclotomic determinant after pullback to the local or global
Galois group. For the permutation representation of the
degree-\(2(p+1)\) field for one factor both classes vanish.

For prime \(p\ge5\), \(\operatorname{SL}_2(\mathbb F_p)\) is perfect
and centrally closed by Steinberg, so has trivial integral first
and second homology. The product kernel \(N\) has the same
vanishing by Künneth. Universal coefficients and the
Hochschild–Serre spectral sequence for \(N\to G\to\mathbb F_p^\times\)
give the isomorphism, and abelianization gives the degree-one
claim. At any higher \(p\)-power level with full image, reduction
has \(p\)-group kernel, which has no higher cohomology with
\(\mathbb F_2\) coefficients. Therefore the same mod-two obstruction
class calculation applies there. This inflation statement concerns
the cohomology of the full image group.
[Steinberg, Lectures on Chevalley Groups, Chapter 7, Theorem 10 and Corollary 2(b)](https://doi.org/10.1090/ULECT/066).

Let \(z\) be the quadratic determinant character and \(\Pi\)
the real line-permutation representation. The determinant of
\(\Pi\) is \(z\), as \(\operatorname{diag}(t,1)\) for a generator
\(t\) of \(\mathbb F_p^\times\) acts by a single \((p-1)\)-cycle
on the nonzero affine points. The representation for the larger
field is \(R=\Pi\oplus(\Pi\otimes z)\), where \(z\) in the tensor
denotes the sign line. With \(n=\dim\Pi=p+1\), the splitting
principle gives
\[
 w_2(\Pi\otimes z)=w_2(\Pi)+(n-1)w_1(\Pi)z+\binom n2 z^2.
\]
Since \(n\) is even, the Whitney sum formula gives
\(w_1(R)=0\) and \(w_2(R)=\binom{p+1}{2}z^2\).
For \(p\equiv3\bmod4\) the binomial coefficient is even. For
\(p\equiv1\bmod4\) the character \(z\) lifts to a character with
values in \(\mathbb Z/4\) through \(D\), so its Bockstein \(z^2\)
vanishes. These classes leave the integer conductor exponents of
the field free to distinguish the three local types.
[Ganguly–Joshi, §7.1 on the cyclic Bockstein](https://arxiv.org/pdf/2011.13245).

By Deligne's formula the normalized orthogonal root number, divided
by that of the determinant line, is the local invariant of
\(w_2\). Hence global Brauer reciprocity applied to these classes
only recovers cyclotomic information. For the line-and-determinant
field \(R\), this root number is one. Its trace form also has
trivial Hasse–Witt invariant by Serre's trace-form formula since
both \(w_2(R)\) and the discriminant squareclass vanish. These
are statements about signs and Brauer classes. The Artin
conductor enters the magnitude of an unnormalized epsilon factor
and retains the integer local distinction.
[Deligne, Theorem (1.5)](https://publications.ias.edu/sites/default/files/Number26.pdf)
· [Serre, Theorem 1](https://www.college-de-france.fr/media/jean-pierre-serre/UPL4088637057849025347_Serre_CMH_Trx2.pdf).

There are classes in \(H^3(G,\mathbb F_2)\) that restrict
nontrivially to \(N\). The first nonzero positive cohomology
of \(\operatorname{SL}_2(\mathbb F_p)\) has degree three and
dimension one, so \(H^3(N,\mathbb F_2)\) is the sum of the three
labeled degree-three classes. Determinant conjugation preserves
each one. In the split Hochschild–Serre sequence their only
potential outgoing differential targets bottom-row degree four,
where inflation is injective by the section. They therefore lift.
Pullback to the absolute Galois group of a finite completion gives
zero in degrees at least three, since that group has cohomological
dimension two. This leaves secondary invariants with local
cochain trivializations as a different possibility.
[Malik–Spallone, Proposition 4.1](https://arxiv.org/pdf/2208.05813).

In arithmetic Chern–Simons theory for
ramified representations, the set of local trivializations of a
pulled-back three-cocycle, modulo coboundaries, is a torsor under
local \(H^2\). The decomposition formula for everywhere unramified
representations selects unramified local trivializations. It does
not supply those for the ramified residue-prime torsion
representations considered here, or a density estimate for their
Newton types.
[Serre, Galois Cohomology, Chapter II §5.3, Proposition 15](https://link.springer.com/book/10.1007/978-3-642-59141-9)
· [Hirano–Kim–Morishita, §3.2](https://arxiv.org/pdf/2106.02308)
· [Chung–Kim–Kim–Park–Yoo, Theorem 1.1](https://arxiv.org/pdf/1609.03012).
