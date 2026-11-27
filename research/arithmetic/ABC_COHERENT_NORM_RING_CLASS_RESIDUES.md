# Ring class groups and coherent norm residues

A coherent prime path has more local structure than the nested determinant
partitions alone require. At each modulus the vertices carry classes in a
quadratic order. Divisibility of an edge coefficient is equality of the
powers of these classes. The possible powers over a prescribed maximal-order
class form one coset of a subgroup. This gives both a capacity for each block
of a partition and a test for which vertices may be in that block.

## Ideal class lifts and exact residue patterns

Let \(K=\mathbb Q(\sqrt{-d})\) for positive squarefree \(d\), with field
discriminant \(D<-4\). Then its units are \(\{\pm1\}\).
Let \(k\ge1\), let \(f\ge1\) be a conductor, and write
\[
 C_f=\operatorname{Pic}(\mathbb Z+f\mathcal O_K),\qquad
 A_f=(\mathcal O_K/f\mathcal O_K)^\times/(\mathbb Z/f\mathbb Z)^\times,
 \qquad C=C_1.
\]
Take \(A_1=0\) and use additive notation for these groups. The exact sequence
\(0\to A_f\to C_f\to C\to0\) sends a residue represented by a local unit
\(\alpha\) to the class of \((\alpha)\) in the prime-to-\(f\) ideal model of
\(C_f\). In this model the trivial principal ideals are exactly those with
a generator congruent to a rational unit at the modulus. Every class of \(C\)
has an ideal representative prime to \(f\), by approximation. A generator
of a prime-to-\(f\) principal ideal has a residue after clearing denominators
prime to \(f\). The residue is well defined in \(A_f\), since the global
units are rational and clearing denominators multiplies by a rational unit.
This gives the kernel identification and its compatibility with reduction to
smaller conductors. See [Stevenhagen, §§2 and 9, equations (2.2), (2.3), and
(9.1)–(9.2)](https://pmb.centre-mersenne.org/item/10.5802/pmb.a-56.pdf).

Define the homomorphism
\[
  \theta_{k,f}:C[k]\longrightarrow A_f/kA_f,\qquad
  t\longmapsto k\tilde t\bmod kA_f ,
\]
where \(\tilde t\) is any lift to \(C_f\). Its kernel is the image of
\(C_f[k]\to C[k]\). In terms of an ideal \(I\) prime to \(f\) representing
\(t\), compute the map by reducing a generator of \(I^k\) in \(A_f/kA_f\).

Fix \(s\ge2\) prescribed classes \(c_i\) in one coset of \(C[k]\), with repetitions
allowed. Take distinct odd unramified split primes \(p_i\nmid f\), with oriented
ideals \(\mathfrak p_i\) in the prescribed classes. As in the
[coherent prime paths](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#witness-selection),
choose generators \(\beta_i\) of
\((\overline{\mathfrak p}_1\mathfrak p_i)^k\), including \(\beta_1=p_1^k\),
and take \(\alpha_{ij}=\bar\beta_i\beta_j/p_1^k\).
Then \((\alpha_{ij})=(\overline{\mathfrak p}_i\mathfrak p_j)^k\) and
\(\alpha_{ij}\alpha_{jh}=p_j^k\alpha_{ih}\). Let \(b_{ij}\) denote the
coefficient of the second element in the integral basis \(1,\omega\), where
\(\omega=(1+\sqrt{-d})/2\) for \(d\equiv3\bmod4\) and
\(\omega=\sqrt{-d}\) otherwise. In the earlier notation \(b_{ij}=v_{ij}/\epsilon_d\).

**Proposition (exact residue criterion).** An array \(g_{ij}\in A_f\) is realizable as the
projective residues \(g_{ij}=[\alpha_{ij}\bmod f]\) if and only if
\[
 g_{ii}=0,\qquad g_{ij}+g_{jh}=g_{ih},\qquad
 g_{ij}\bmod kA_f=\theta_{k,f}(c_j-c_i).                         \tag{1}
\]
For fixed \(K,k,f,s,(c_i),(g_{ij})\) satisfying (1), the array occurs on
distinct primes in every sufficiently large dyadic interval. Among ordered
oriented tuples in \([Z,2Z)\) with the prescribed classes, its limiting
proportion as \(Z\to\infty\) is \(|kA_f|^{1-s}\).

Indeed, let \(x_i\in C_f\) be the lift given by the oriented prime \(\mathfrak p_i\).
Then \(g_{ij}=k(x_j-x_i)\) in the kernel \(A_f\). Fix a reference lift \(x_*\)
of \(c_1\) and put \(y_i=k(x_i-x_*)\). As \(x_i\) runs over the lifts of \(c_i\),
\(y_i\) runs over the entire coset of \(kA_f\) of color
\(\theta_{k,f}(c_i-c_1)\), with every value attained \(|A_f[k]|\) times.
For any array satisfying (1), choose \(y_1\in kA_f\) and take
\(y_i=y_1+g_{1i}\). This gives precisely \(|kA_f|\) compatible assignments.
All lifts \(x_i\) can then be chosen independently. The ring class field
has Galois group \(C_f\) over \(K\), with the Artin map sending an unramified
prime ideal to its lift. Chebotarev over \(K\) gives equal asymptotic counts
for its classes in \([Z,2Z)\). Prime ideals of absolute residue degree two
in this interval have norms \(p^2\) and contribute \(O(\sqrt Z)\) in total.
Consequently the oriented split primes in each class of \(C_f\) have the
same leading count. Tuples with a repeated rational prime contribute
\(O_s((Z/\log Z)^{s-1})\), since a rational prime fixes at most two choices of
orientation at each index. The limit follows by dividing the \(|kA_f|\) compatible
assignments by the \(|kA_f|^s\) equally likely ones.
[Lenstra, Theorem 3.1](sources/lenstra-chebotarev-density-notes.pdf#page=3)

## Exact capacities of nested divisibility partitions

For each \(\ell^E\Vert f\) prescribe partitions \(\mathcal P_{\ell,e}\) of the
vertices for \(1\le e\le E\), each refining the preceding one. Start each
sequence with the one-block partition \(\mathcal P_{\ell,0}\). Put
\[
  B_{\ell,e}=kA_{\ell^e},\qquad
  b_{\ell,e}=\big|\ker(B_{\ell,e}\longrightarrow B_{\ell,e-1})\big|,
  \qquad
  \lambda_{\ell,e}(i)=\theta_{k,\ell^e}(c_i-c_1).
\]
Here \(b_{\ell,e}\) is a capacity, with two subscripts for a prime and an
exponent. The coefficient \(b_{ij}\) has vertex subscripts.

**Proposition (partition criterion).** These partitions occur simultaneously with the
prescribed maximal-order classes, as the equivalences
\(i=j\) or \(\ell^e\mid b_{ij}\), exactly under the following conditions.

1. Each block of \(\mathcal P_{\ell,e}\) has one color \(\lambda_{\ell,e}\).
2. For each parent block of \(\mathcal P_{\ell,e-1}\) and each color at level
   \(e\), at most \(b_{\ell,e}\) of its child blocks have that color.

The equality \(\ell^e\mid b_{ij}\) means the residue of \(\alpha_{ij}\),
a unit at \(\ell\), is rational modulo \(\ell^e\). Thus it is equivalent
to \(y_i=y_j\) after reduction to \(A_{\ell^e}\). Reduction of the local
unit groups is surjective. It follows that \(B_{\ell,e}\to B_{\ell,e-1}\)
is surjective, and a coset of \(B_{\ell,e}\) maps onto its compatible coset
at level \(e-1\). Over a fixed parent label there are exactly \(b_{\ell,e}\)
available labels of each compatible child color. Distinct children of one
parent must use different labels. Children with different colors already
have different labels. Starting at level zero and assigning labels
recursively proves necessity and sufficiency at each prime. The Chinese
remainder theorem identifies \(A_f\) with the product of \(A_{\ell^E}\).
Combining the labels, taking their pairwise differences, and applying (1)
realizes all the partitions on the same tuple.

In fact the fixed-class limiting proportion of a partition pattern can be
read off from the recursion. Multiply \((b_{\ell,e})_u\) for every parent
and every child color, with \(u\) the number of its children of that color
and \((b)_u=b(b-1)\cdots(b-u+1)\), including \((b)_0=1\). Divide the product
by \(|kA_f|^s\). The numerator counts all compatible label assignments.
Within one coset of \(C[k]\), if the base classes are also left free,
the labels \(y_i\) are instead uniformly distributed on the subgroup
\[
 \{a\in A_f:a\bmod kA_f\in\operatorname{im}\theta_{k,f}\}.
\]
Indeed the eligible lifts form a coset of the inverse image of \(C[k]\)
and \(k\) times that subgroup has this image in \(A_f\). For example the
collision probability for all \(s\) labels modulo \(f\) in that model is
\((|kA_f|\,|\operatorname{im}\theta_{k,f}|)^{1-s}\). When base classes
vary, their colors at different primes are restricted to the image of the
*single* homomorphism \(\theta_{k,f}\).

For \(k=3\) and \(\ell\ne3\), let \(\chi_D(\ell)\) be the Kronecker
symbol, including \(0\) at ramified primes. The first group is cyclic
of order \(\ell-\chi_D(\ell)\): it is the quotient of the two split unit
groups by the diagonal scalars, the inert finite-field unit group modulo
scalars, or the additive radical in the ramified dual-number ring.
At later levels the quotient \(|A_{\ell^e}|/|A_{\ell^{e-1}}|\) is \(\ell\),
from the two-dimensional unit lift divided by the scalar unit lift.
Since cubing is invertible on the kernel of reduction for \(\ell\ne3\),
the first level has capacity
\((\ell-\chi_D(\ell))/\gcd(3,\ell-\chi_D(\ell))\) per color and every later
level has capacity \(\ell\). The first color determines each later one.
For unramified \(3\), \(|A_3|=3-\chi_D(3)\) is prime to \(3\) and
\(\ker(A_9\to A_3)\) has order \(3\). The order-three part of \(A_9\) is
this kernel, so \(3A_9\to3A_3\) is an isomorphism. Level two then has
capacity one per color. For example two vertices with equal color at
modulus \(9\) and \(3\mid b_{ij}\) must have \(9\mid b_{ij}\). All conductors,
including higher powers of \(3\) and ramified \(3\), use the group formula
in the partition criterion.

## Cubic fields detect the first color

For \(k=3\), each nonzero \(t\in C[3]\) gives a totally real \(S_3\) splitting
field \(S_t\). For an integral representative \(I\) with \(I^3=(\alpha)\)
and \(n=NI\), take the splitting field of
\[
                     F_t(T)=T^3-3nT-\operatorname{Tr}(\alpha).
\]
It depends only on the line \(\langle t\rangle\). Its quadratic subfield is
\(K^\vee=\mathbb Q(\sqrt{3d})\). The discriminant of each of its cubic
subfields is \(\operatorname{disc}(K^\vee)3^{2a}\) with \(0\le a\le2\).
The trinomial construction and the unramified-outside-\(3\) property are
classical. [Kishi, Proposition 2.1 and Remark 2.2](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/460B48BA93CA629154706CF7BFE7430D/S0017089510000431a.pdf/on_the_ideal_class_group_of_certain_quadratic_fields.pdf)

Here is the relation with the colors. The map \(t\mapsto[\alpha]\) is an
injection into \(K^\times/K^{\times3}\). Changing representatives multiplies
\(\alpha\) by a cube, and if \(\alpha\) is a cube, unique ideal factorization
makes \(I\) principal. The map remains injective over \(M=K(\zeta_3)\), by
taking norms in the degree-two extension. For \(\beta^3=\alpha\), the three
roots are \(\zeta_3^j\beta+\zeta_3^{-j} n/\beta\), \(0\le j<3\). They are
real because \(|\beta|^2=n\). The discriminant is
\(27(4n^3-\operatorname{Tr}(\alpha)^2)\), a nonzero nonsquare of squareclass
\(3d\). The three roots recover \(\beta\) over \(M\) by taking their
linear combination with coefficients \(\zeta_3^{-j}/3\). Thus
\(S_tM=M(\beta)\), a nontrivial cyclic cubic extension of \(M\) permuting
the three roots cyclically. Its real subfield is \(S_t\), since
\([S_tM:S_t]=2\). These identities give the claimed irreducibility and
independence of the representative. Distinct lines give distinct fields:
in a cyclic cubic Kummer extension over \(M\), an element with cube in
\(M\) is a Galois eigenvector, hence is in one of the three lines over \(M\)
spanned by \(1,\beta,\beta^2\).

At every prime of \(M\) away from \(3\), the valuation of \(\alpha\) is a
multiple of three. After removing a local cube, the polynomial for
\(\beta\) is a unit cube-root equation with separable reduction, so the
extension there is unramified. Ramification of degree three in
\(S_t/K^\vee\) would persist after the base change of degree two to \(M\).
Thus outside \(3\), the only possible nontrivial inertia in \(S_t/\mathbb Q\)
has order two and occurs at primes ramified in \(K^\vee\). Its decomposition
group normalizes this transposition, so is itself of order two. The cubic
algebra there is \(\mathbb Q_\ell\) times the local quadratic subfield and
has the same discriminant exponent as \(K^\vee\), including at \(2\). At
\(3\), Kishi's Remark 2.2 bounds the cubic discriminant exponent by five.
The field discriminant has squareclass \(\operatorname{disc}(K^\vee)\),
whose exponent at \(3\) is zero or one. This gives the displayed exponent.

Let \(\ell\nmid 6D\) split in \(K^\vee\), equivalently
\(3\mid\ell-\chi_D(\ell)\). Choose \(I\) prime to \(\ell\).
Then \(\theta_{3,\ell}(t)=0\) exactly when \(\ell\) splits completely
in \(S_t\). To see this, if \(\ell\) splits in \(K\), it is \(1\bmod3\).
Write the two residues of \(\alpha\) as \(a,b\in\mathbb F_\ell^\times\).
Since \(ab\) is a cube, \(a/b\) is a cube exactly when both residues are cubes.
If \(\ell\) is inert in \(K\), it is \(2\bmod3\). Rational units are then
cubes, so the image in \(A_\ell\) is a cube exactly when \(\alpha\) is a cube
in \(\mathbb F_{\ell^2}^\times\). The residue field at a prime of \(M\) is
respectively \(\mathbb F_\ell\) or \(\mathbb F_{\ell^2}\). In either case the
test is the splitting test of the unit cube-root equation. The residue degree
of \(M/K^\vee\) is at most two, so base change detects the identity in the
cyclic group of order three. This argument works when \(\ell\) divides the
trinomial discriminant because of its index. In particular an inert prime
in the associated cubic field cannot divide \(b_{ij}\) for any edge with
difference \(t\), provided \(\ell\nmid 6Dp_ip_j\). Primes dividing the incident
norms are already excluded by primitivity.

For a two-dimensional span \(U\subset C[3]\), the four lines give four
distinct fields. Two fields for a basis have intersection \(K^\vee\), and
their compositum has group \(C_3^2\) over \(K^\vee\). It contains all four fields
by the Kummer description after adjoining \(\zeta_3\) and taking the real
subfield. At any good \(\ell\) split in \(K^\vee\), either all four lines
split completely or exactly one does, according as
\(\theta_{3,\ell}|_U\) is zero or has a one-dimensional kernel.
Among all rational primes the first case has density \(1/18\); the case
with one specified splitting line has density \(1/9\).
Chebotarev counts the identity and the two nonzero elements of each
Frobenius line of \(C_3^2\) in the compositum of degree \(18\).
The two elements form one conjugacy class, since the outer involution
acts by inversion. In particular the three edge fields of an independent
triangle can have zero, two or three inert edge fields at such an \(\ell\).
A good \(\ell\) splitting in \(K^\vee\) and dividing all three coefficients
must split in all four fields.

## Relation to the classical ring spaces

For an odd prime exponent \(k\), Mayer's \(k\)-ring space consists of classes
of generators of \(k\)th powers of ideals that become rational modulo the
conductor after multiplication by a local \(k\)th power of a global element.
In the present unit range these generator classes identify with \(C[k]\).
Being rational after the multiplication is exactly
\(\theta_{k,f}(t)=0\). Mayer uses the spaces to count dihedral extensions
by conductor. The result here uses the individual cosets of \(kA_f\) and
their fibers under reduction to characterize prime-tuple residue and
partition patterns. Mayer also gives ring-space data for \(d=4027\).
[Mayer, Definition 2.2, §2.5 and Table 10](https://arxiv.org/pdf/1403.3906v1)

## Exact fixtures

In \(K=\mathbb Q(\sqrt{-4027})\), take primes \(13,17,19\) and signed
\((u,v)\) pairs \((6326,28),(3265,111),(7549,139)\) in order \(12,13,23\).
The three norm equations \(u^2+4027v^2=4(p_i p_j)^3\), incident units,
and \(\bar\alpha_{12}\alpha_{13}=13^3\alpha_{23}\) certify a coherent
triangle. Write \(t_1,t_2\) for the classes of edges \(12,13\). The
trinomial root counts in order \(12,13,23\) are \((0,3,0)\) modulo \(5\)
and \((3,0,0)\) modulo \(59\). Both primes are away from all three
polynomial discriminants. The first comparison shows \(\theta_{3,5}\)
vanishes on \(t_2\) and is nonzero on \(t_1\); the second shows \(t_2\ne0\).
Thus the triangle is independent. Edges \(12\) and \(23\) cannot have
coefficient divisible by \(5\) for any primes in these same vertex
classes away from \(5\). Edge \(13\) meets the color test, but \(5\) does
not divide its actual coefficient \(111\). Modulo \(5\) the group \(A_5\)
has order six, so there are two labels per color. Three singleton blocks
with all three vertices prescribed the same base class are impossible
at this modulus even though they form a partition with constant color
on each block.

For the modulus \(9\), use \(d=23\), primes \(13,31\), and \((u,v)=(13556,1842)\).
The norm and incident unit conditions again give an oriented edge. The
coefficient has valuation one at \(3\). The two roots of \(-23\) modulo
\(9\) are \(\pm2\), and the residues of \(\alpha\) have ratio \(4\bmod9\)
in that order. The cubes of units modulo \(9\) are \(\pm1\); this edge has
nonzero level-two color. So \(3\mid b\) and \(9\nmid b\) is compatible
here, while it is impossible for a principal edge for this \(K\).
These fixtures use only the displayed integer norm, product and residue tests.

## Use in the prime-path count

The path selection can retain independent triangles without changing its
power of \(X,Z,H\). For \(K=\mathbb Q(\sqrt{-d})\) and primes in \([Z,2Z)\)
with \(16Z^2<d\), the oriented maximal-order classes of distinct primes
are distinct. A class equality would give a generator of
\(\overline{\mathfrak p}_i\mathfrak p_j\) of norm \(p_i p_j<d/4\).
Every nonrational integral element of \(K\) has norm at least \(d/4\),
and a rational generator has a conjugation-invariant ideal.
In one coset of \(C[3]\) with \(s\) distinct classes, each pair determines
at most one third vertex on its affine line over \(\mathbb F_3\), so
only \(O(s^2)\) ordered paths have dependent differences. The
[Fermat selection](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#witness-selection)
gives \(\gg s^3\) differently colored paths, with an absolute constant.
Raise the fixed threshold for a large coset to absorb the \(O(s^2)\)
loss. In the tail range, \(16Z^2<d\) for sufficiently small auxiliary
\(\eta\) and large \(X\). The same Hölder sum then proves the existing
per-field lower bound (2) for paths satisfying both restrictions.

This adds a low-discriminant Frobenius test to the
[rank-one occupancy problem](ABC_COHERENT_NORM_PENCILS.md#exact-examples-and-the-unresolved-count).
For any reconstructed independent path in the field window, the three
trace trinomials are attached to one plane of \(3\)-classes and give
three distinct fields of discriminant \(O(X)\) with the same resolvent.
The same color exclusions hold for divisors of the normalized
coefficients, since these divide the original coefficients. The common
scale can only use primes passing all the corresponding color tests.
An aggregate upper bound for the occupied rank-one records of
\(O_\varepsilon(X^\varepsilon V^{3-\sigma})\), for \(0<\sigma\le1/9\)
in the range of that problem, would feed into the existing tail and its
[elliptic-conductor transfer](ELLIPTIC_CONDUCTOR_MOMENTS_FROM_QUADRATIC_TORSION.md).
At \(\sigma=1/9\) the exponent would be \(29/24\) for both the quadratic
second moment and the elliptic-conductor \(2\gamma\)-moment. The new test supplies a
restriction on the records, with the required aggregate bound still open.

For fixed \(K\) and plane, a count using only allowed prime supports of the
three coefficients has no fixed power saving. Choose all three coefficients
to be distinct primes in \([V,2V)\) splitting completely in the fixed common
splitting field. Chebotarev gives \(\gg_{K,U} V^3/(\log V)^3\) such coefficient
triples. They satisfy the support tests and pairwise coprimality. The
full-modulus congruences and the field window are additional restrictions.
The fixed-field pattern densities above also have errors depending on the
ring class field and do not give a summation over moving discriminants or
growing conductors. A remaining question is whether the full-modulus lattice
test correlates with the Frobenius restrictions strongly enough to save a
power across the field window.

## Cost of the Frobenius conditions on class labels

Let \(f\) be prime to \(6D\), and let \(a\) be the number of its prime
divisors split in \(K^\vee\). For an imaginary quadratic field as above,
\(\theta_{3,f}\) has image of size at most \(3^a\), even if \(f\) has
prime-power factors. In particular the \(\mathrm{Cl}(K)[3]\)-space of
classes passing *all* the first-color tests has size at least
\(|\mathrm{Cl}(K)[3]|/3^a\). If its size is \(h_0\ge9\), there are
\((h_0-1)(h_0-3)\) ordered independent pairs in it. This stronger
simultaneous test ensures that each of the four lines of every pair
passes all the prescribed prime tests. It therefore bounds from below
the number of independent class labels surviving any assignment of the
prime divisors of \(f\) to particular edges. For \(f\le X^B\) and \(|\mathrm{Cl}(K)[3]|\ge X^\tau\) with fixed
\(B,\tau>0\), the loss \(3^{2a}\) is subpower in \(X\). The first-color
tests on the class labels alone therefore leave
\(\gg_{B,\varepsilon}|\mathrm{Cl}(K)[3]|^2 X^{-\varepsilon}\) ordered pairs
for sufficiently large \(X\) depending on \(B,\varepsilon,\tau\). This is a statement about class
labels. Matching those labels with norm generators of specified size
and residue can impose additional restrictions.

## Real quadratic fields and the common unit image

Let \(K\) now be real quadratic. Use totally positive generators and write
\(C^+\) and \(C_f^+\) for the narrow class groups of the maximal order and
the order of conductor \(f\). The group \(A_f\) is the same projective unit
residue group as above. Write \(U_f\) for the image in \(A_f\) of the totally
positive units of \(\mathcal O_K\). The narrow ring class sequence is
\[
             0\longrightarrow A_f/U_f\longrightarrow C_f^+
                  \longrightarrow C^+\longrightarrow0 .
\]
Indeed a local unit residue at \(f\) can be represented by a totally
positive element by adding a sufficiently large positive rational
multiple of \(f\). Two totally positive generators of the same ideal
differ by a totally positive unit. The narrow principal criterion uses
a positive norm, equivalent to total positivity after multiplying by
\(-1\) when both embeddings are negative.
[Stevenhagen, equations (9.1)–(9.2)](https://pmb.centre-mersenne.org/item/10.5802/pmb.a-56.pdf)

Define \(\theta^+_{k,f}\) on \(C^+[k]\) by taking \(k\) times a lift to
\(C_f^+\) and reducing in \(A_f/(kA_f+U_f)\).
Fix classes \(c_i\) in one coset of \(C^+[k]\).

**Proposition (real quadratic residue criterion).** Prescribed raw projective
residues \(g_{ij}\in A_f\) occur for coherent totally positive \(k\)-power
generators if and only if they satisfy the cocycle in (1) and
\[
       g_{ij}\bmod(kA_f+U_f)=\theta^+_{k,f}(c_j-c_i).              \tag{6}
\]
For fixed admissible data the primes can be chosen distinct and in
every sufficiently large dyadic interval. Nested coefficient-divisibility
partitions occur exactly when there are vertex labels \(y_i\in A_f\) with
\(y_i\bmod(kA_f+U_f)=\theta^+_{k,f}(c_i-c_1)\), constant on each prescribed
block and distinct between blocks after the corresponding local reduction.

Reduce first to \(A_f/U_f\), where the lifting and Chebotarev argument of
(1) applies word for word to the narrow Artin map. Choose any coherent
totally positive generators at the primes obtained. Their projective
residues differ from the prescribed array by a cocycle in \(U_f\).
Every such cocycle is \(e_j-e_i\) for vertex labels \(e_i\in U_f\).
Multiply the anchor generators by totally positive units representing
these labels. Conjugation inverts a totally positive unit because its norm
is \(1\), so this changes the edge residues by exactly \(e_j-e_i\) and
preserves coherence. For the partition claim take \(g_{ij}=y_j-y_i\),
or recover \(y_i=g_{1i}\) from a realized tuple. This allows arbitrary unit
powers. The existence claim has no restriction on archimedean generator size.

The Chinese remainder theorem splits \(A_f\). It need not split the unit
image \(U_f\): one unit power changes residues at all the primes of \(f\).
The following obstruction applies to partitions even when only divisibility
is specified and no nonzero residue is prescribed.

**Proposition (a simultaneous partition obstruction).** For all large \(X\)
there are \(\gg\sqrt X\) real quadratic fields of discriminant at most \(X\)
for which the following partitions of three
totally positive coherent cube generators in the principal narrow coset
are individually realizable on the same split prime tuple but cannot
coexist on any such tuple of primes away from \(35\).

- At \(5\) the three vertices are singletons, so no edge coefficient is
  divisible by \(5\).
- At \(7\) all three vertices are in one block, so all edge coefficients
  are divisible by \(7\).

To construct the fields, take \(m\equiv3\bmod35\), \(m\ge3\), and write
\(m^2-1=d b^2\) with \(d>1\) squarefree and \(b>0\).
The totally positive unit \(\eta=m+b\sqrt d\) has norm \(1\). Both \(d\)
and \(b\) are units at \(5,7\). The prime \(5\) is inert in
\(\mathbb Q(\sqrt d)\), and \(7\) splits there.
Indeed the squareclass of \(d\) at each is that of \(8\).
Thus \(A_5,A_7\) each have order six. The nonrational coefficient of \(\eta\)
is nonzero in both residue algebras, but
\(\eta^3=m(4m^2-3)+(4m^2-1)b\sqrt d\) is rational at both, since
\(4m^2-1\equiv0\bmod35\). The rational residue is a unit. So the unit has
order three in both projective groups. The image of the totally positive
units in \(A_{35}/3A_{35}\) is cyclic, hence a line with nonzero projection
in each of its two order-three factors.

For vertices in the principal narrow class, their raw labels lie in
\(3A_{35}+U_{35}\) by (6). Equality of all their labels at \(7\) makes
their images equal in \(A_7/3A_7\), and the unit line then makes their
images equal in \(A_5/3A_5\). They have at most two distinct labels at \(5\),
the size of one coset of \(3A_5\), so cannot form three singleton blocks.
To realize the two local requirements separately, choose three primes
with trivial narrow ring class at \(35\). Their totally positive
generators can be taken rational modulo \(35\). Cubing gives one block
at \(7\). For the other pattern, multiply the three vertex cube generators
by \(1,\eta,\eta^2\), giving three different projective residues at \(5\).
For \(d=2\), the narrow prime generators \(11+7\sqrt2,13+7\sqrt2,15+7\sqrt2\)
of norms \(23,71,127\) also witness the two separate requirements:
they are rational modulo \(7\), and their cube residues at \(5\) have
color zero before multiplication by the three different unit powers.
The resulting three colors at \(5\) are distinct. In edge order \(12,13,23\),
the coefficients after these unit changes are \((1,2,2)\bmod5\).

Restrict \(m\) to be even and \(m\le\sqrt X/2\), and require \(m^2-1\)
squarefree. The factors \(m-1,m+1\) are odd and coprime. For each odd
prime \(p\nmid35\) the squarefree condition excludes exactly two residues
modulo \(p^2\), and \(5,7\) already divide neither factor. The Chinese remainder
theorem for primes \(p\le L\) gives a positive main density in the fixed
progression, bounded below independently of \(L\) by the convergent
product of \(1-2/p^2\). Among \(m\le M\), the number lost from
\(p>L\) is \(O(M/L+\sqrt M)\), summing \(O(M/p^2+1)\) over
\(L<p\le\sqrt{M+1}\). Choose \(L\) fixed large enough, then let \(M\) grow.
This gives \(\gg M\) squarefree choices and thus distinct \(d=m^2-1\)
with \(D=4d\le X\) at the stated cutoff.


The real quadratic unit correction concerns direct residue prescriptions.
The signed \(3\)-torsion tail uses Scholz reflection to transfer the
imaginary-field estimate, where the global units of the relevant
fields are rational.
