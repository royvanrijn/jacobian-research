# Two involution factors and simple branches: a complete component classification

Let \(N,A,B,s\) satisfy
\[
 s\ge1,\qquad 2\le B\le A\le\lfloor N/2\rfloor,\qquad
                         A+B+s=N-1.                         \tag{1}
\]
Consider polynomial tuples with fixed product \(P_N\) and passport
\[
       (2^A1^{N-2A}),\quad(2^B1^{N-2B}),\quad
                         s\text{ transpositions}.           \tag{2}
\]
The covering source has genus zero. The parameter \(s\) will also
be the genus of an auxiliary hyperelliptic curve; these two genera
must not be confused.

All components of (2) are classified below, including their literal
monodromy groups, full block lattices, deck groups and orbit sizes.
Component completeness uses a published theorem of Bogatyrëv and
Gendron. The group identification, permutation-level label extraction
and counting formulas are derived here. This is one family inside
the [arbitrary-passport problem](POLYNOMIAL_PASSPORT_COMPONENTS.md);
three or more bundled involutions and other arbitrary disconnected
cycle types remain outside its scope.

The braid group is the full \(B_{s+2}\), permitting all placements
of the prescribed types, with the convention
\(Q_i(X,Y)=(Y,Y^{-1}XY)\). A fixed ordered pattern is used for
counting only; no assertion about pure or individually marked
braid components is made.

## Primary theorem and conventions

Bogatyrëv--Gendron,
[*The space of solvable Pell--Abel equations*](sources/bogatyrev-gendron-pell-abel-2025.pdf),
Compositio Mathematica **161** (2025), 1483--1511, classify precisely
the passports (2) in published Corollary 2, printed pp1485--1486.
Their Theorem 1.2, printed p1484, supplies the finer degree-partition
labels. Their \(g\) equals \(s\) here. The proof uses graph normal
forms in Theorem 5.4 and their separation by degree partitions in
Section 6.1. These exact statements and the indicated proof route
were inspected in the retained published version.

After sending the two bundled values to \(+1,-1\), extract the
square-free part \(D\) of \(P^2-1\). Then \(\deg D=2s+2\), and
the Pell--Abel equation has a least-degree nonconstant solution
\(U\), of degree \(n\). The original polynomial has the form
\[
                         P=\pm T_m\circ U,\qquad N=mn,        \tag{3}
\]
where \(T_m\) is the classical Chebyshev polynomial. This fact is
part of their Section 2, Theorem 2.1 and Remark 1(2).
The primitive degree partition is the unordered pair
\(\{e,f\}\) counting roots of \(D\) where \(U=+1,-1\).

Pell primitivity in (3) means least degree for the fixed equation.
Its equivalence to permutation primitivity in this family is proved
in B2; it is not an assumption about all polynomial compositions.

Their Hurwitz components are absolute full-braid components in our
conventions. The [full-twist calculation L3](POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md)
puts all simultaneous \(P_N\)-conjugates of a tuple in its literal
braid orbit. Thus each absolute component corresponds to exactly
one literal full component. The quotient by connected affine
coordinate changes does not change these component counts.

## Theorem B1. Complete imported component labels

Write \(E=N-2A\) and \(F=N-2B\), so \(0\le E\le F\) and
\(E+F=2s+2\).

If \(E>0\), the components are indexed by
\[
       n\mid N,\qquad m=N/n\text{ odd},\qquad n\ge F,         \tag{4}
\]
with one component for each such \(n\). Its primitive degree
partition is \(\{E,F\}\).

If \(E=0\), the components are indexed by
\[
 \begin{gathered}
 n\mid N,\qquad m=N/n\text{ even},\\
 1\le e\le f\le n,\qquad e+f=2s+2,\qquad e\equiv f\equiv n
                                                    \pmod2. \tag{5}
 \end{gathered}
\]
There is one component for each pair \((n,\{e,f\})\).
In particular, for each \(n\) in the second case their number is
\[
 a(s,n)=
 \begin{cases}
 0,&n<s+1,\\
 \lfloor t/2\rfloor+1,&n\ge s+1,\ n+s\text{ odd},\\
 \lfloor(t+1)/2\rfloor,&n\ge s+1,\ n+s\text{ even},
 \end{cases}
 \qquad t=\min(s,n-s-1).                                     \tag{6}
\]

**Justification.** This is the component statement of the published
Corollary 2, using Theorem 1.2 to retain the individual partition
labels in its second case. Equations (4)--(6) are the same
hypotheses and counts, with \(g=s\) and \(m=N/n\).
For odd \(m\), the primitive partition is \(\{N-2A,N-2B\}\);
for even \(m\), one original bundled factor is fixed-point-free,
and all the admissible primitive partitions occur. This distinction
is explicit immediately before that corollary. No claim of a new
component-connectedness proof is made. \(\square\)

## Theorem B2. Monodromy, every block system and deck symmetry

Let \(\Delta_m\) be the natural degree-\(m\) Chebyshev monodromy:
the trivial group for \(m=1\), \(C_2\) for \(m=2\), and the standard
dihedral group of order \(2m\) for \(m\ge3\).
The component with least Pell degree \(n\) and multiplier \(m=N/n\)
has the same standard embedded group
\[
                  G=S_n\wr\Delta_m                          \tag{7}
\]
for every one of its partition labels.
Its block-count divisors are exactly
\[
                       \{d:d\mid m\}\ \cup\ \{N\}.           \tag{8}
\]
Its deck group is trivial for \(n\ge3\) and has order two for
\(n=2\). In particular, the original polynomial is permutation
primitive exactly when \(n=N\), and then its group is \(S_N\).

The canonical degree-\(n\) blocks are the connected components of
the support graph of conjugates of any original transposition.
They refine every proper block system, and are independent of
which transposition is chosen.

**Proof.** The transposition-normal-closure lemma L6 gives
\(G=S_a\wr K\) on its canonical blocks, for some \(a\ge2\).
Every transposition belongs to this full symmetric base. The
quotient is therefore generated by the two involutions coming
from the bundled factors, with product \(P_d\), where \(d=N/a\).
If they are \(R,S\), then \(R P_d R=P_d^{-1}\).
Consequently their group is the standard \(\Delta_d\): for
\(d\ge3\), \(R\) has form \(x\mapsto c-x\) in the cyclic labels;
degrees one and two give the stated small cases.

These canonical blocks refine every proper block system. Indeed,
if an edge transposition \((u\,v)\) crosses the boundary of a
block containing some \(w\ne u\), its image intersects that block
at \(w\), and hence equals it, a contradiction. Every edge is
therefore inside one block. Its connected components refine the
partition. As \(\Delta_d\) preserves every congruence partition
of its \(d\)-cycle, all block-count divisors of \(d\) occur and
there are no other proper block systems.

It remains to identify \(a\) with the least Pell degree \(n\).
A polynomial quotient with at most two finite branch values and
involutory branch permutations has an alternating-path plane
tree. It is the Chebyshev quotient \(T_d\), up to the coordinate
choices fixing the two values. This includes \(d=2\), where one
of the two quotient entries may be the identity. The usual
block/decomposition correspondence (Müller's
[Appendix, retained p12](sources/mueller-primitive-polynomial-monodromy-author.pdf),
via Lüroth's theorem), together with H4's plane-tree interpretation,
thus writes \(P=\pm T_d\circ V\) for the canonical quotient.
The identity
\[
          T_d(z)^2-1=(z^2-1)\,\mathcal U_{d-1}(z)^2           \tag{9}
\]
shows that \(V^2-1\) has the same square-free part \(D\).
Here \(\mathcal U_{d-1}\) is the Chebyshev polynomial of the
second kind. The Pell solution theorem therefore gives \(n\mid a\).

Conversely (3) supplies a block system of size \(n\) whenever
\(m>1\). The canonical transposition blocks refine it, so \(a\le n\).
Thus \(a=n\). If \(m=1\), a proper canonical quotient would give
a smaller Pell solution by (9), contradicting minimality; again
\(a=n=N\).
Equation (7) and the full lattice (8) follow. Full symmetric base
groups permit every fibre relabeling, so (7) is equality in the
fixed congruence labels, not just abstract group isomorphism.
Finally L6 identifies the natural centralizer: trivial for \(n\ge3\)
and the diagonal \(C_2\) for \(n=2\). Since it commutes with \(P_N\),
this is exactly the deck group used in the absolute quotient.
\(\square\)

## Theorem B3. Reading the complete component label from permutations

The labels in B1 can be recovered without a polynomial equation.
Compute the transposition-normal-closure blocks, giving \(n,m\)
as in B2. If \(m\) is odd, the label is \(n\).
If \(m\) is even, project to the \(m\) blocks and identify the
bundled factor with two fixed quotient blocks. Let its two return
involutions have \(a,b\) transpositions on their \(n\) sheets.
Then its extra label is
\[
                         \{e,f\}=\{n-2a,n-2b\}.              \tag{10}
\]
Together \(n\) and (10) are a complete component invariant for (2).

**Proof.** Two involutions of full-cycle product and total quotient
index \(m-1\) have one fixed point each when \(m\) is odd. For
even \(m\), one is fixed-point-free and the other has two fixed
points. A nontrivial return on a quotient 2-cycle would give an
original cycle of length at least four, which is excluded in (2).
All remaining heavy returns are involutions at the indicated
fixed quotient points. Each simple factor has one transposition
return in one block.

In odd degree the original fixed-point counts give \(E,F\)
directly. In even degree the two fixed points are the two
unramified endpoints over the bundled Chebyshev value. Their
returns are the primitive inner monodromies over \(+1,-1\);
their fixed-point counts are (10). Their sum is \(2s+2\) by
the return budget. Normalized tree charts change neither count.
B1 now proves sufficiency as well as necessity. \(\square\)

## Lemma B4. An exact ordered class count

Allow \(a\) or \(b\) to be zero in this lemma, retaining the
corresponding identity as a tuple position. Let \(s\ge1\),
\(0\le a,b\le\lfloor n/2\rfloor\), and \(a+b+s=n-1\).
The number of ordered tuples of types
\((2^a1^{n-2a}),(2^b1^{n-2b}),T,\ldots,T\), with product \(P_n\),
is
\[
 C_n(a,b;s)=
 n^{s+1}\frac{(n-a-1)!(n-b-1)!}
                  {a!\,b!\,(n-2a)!\,(n-2b)!}.               \tag{11}
\]
Set \(C_n=0\) for invalid indices. This is a literal count;
no division by \(n\) or rounding of a rational class weight occurs.

**Proof.** Put \(h=n-1\), and define
\[
 K_j(a)=[x^j](1+x)^{h-a}(1-x)^a.
\]
These are the hook characters of the involution class. For a
transposition,
\(\chi_j(T)=\binom hj(h-2j)/h\).
Use the exact long-cycle Frobenius formula in the
[general hook-count derivation](POLYNOMIAL_PASSPORT_COMPONENTS.md#exact-counting-and-retained-evidence).

For completeness the needed finite sum can be evaluated directly.
Binomial expansion gives
\[
 \frac{K_j(a)}{\binom hj}=\frac{K_a(j)}{\binom ha},\qquad
 (-1)^jK_j(b)=K_j(h-b).
\]
Moreover
\[
 \sum_{j=0}^h\binom hj K_u(j)K_v(j)
                         =2^h\binom hu\,\delta_{uv},         \tag{12}
\]
by extracting coefficients from
\(\sum_j\binom hj(1+x)^{h-j}(1-x)^j
 (1+y)^{h-j}(1-y)^j=2^h(1+xy)^h\).
The coefficient recurrence is
\[
 (h-2j)K_k(j)=(k+1)K_{k+1}(j)+(h-k+1)K_{k-1}(j).
\]
Since \(a+s=h-b\), the coefficient of \(K_{h-b}\) in
\((h-2j)^sK_a(j)\) is \((h-b)!/a!\); every step must raise
the index. Orthogonality eliminates every other term.
Thus the hook sum is
\[
 \sum_{j=0}^h
 \frac{(-1)^jK_j(a)K_j(b)}{\binom hj}
                   \left(\frac{h-2j}{h}\right)^s
              =\frac{2^h(h-a)!(h-b)!}{h!\,h^s}.
\]
Multiplying by the class sizes
\(n!/(2^a a!(n-2a)!)\),
\(n!/(2^b b!(n-2b)!)\), and \((nh/2)^s\), and dividing by
\(n!\), proves (11). Identity positions satisfy the same calculation.
\(\square\)

## Theorem B5. Exact literal and absolute component sizes

For an admissible positive partition \(e+f=2s+2\) on degree \(n\),
define the primitive ordered count
\[
 P_n(e,f;s)=
 \sum_{\substack{k\mid n,\ k\ {\rm odd}\\n/k\ge\max(e,f)}}
 \mu(k)\,k^{s+1}
 C_{n/k}\left(\frac{n/k-e}{2},\frac{n/k-f}{2};s\right).        \tag{13}
\]
All indices here are integers because \(e\equiv f\equiv n\pmod2\).
For a component in B1 put
\[
                  L=\frac{(s+2)(s+1)}{1+\mathbf1_{A=B}}.
\]
Its literal full orbit size is
\[
 \begin{cases}
 L\,m^{s+1}P_n(E,F;s),&m\text{ odd},\\
 L\,\dfrac m2\,m^s(2-\mathbf1_{e=f})P_n(e,f;s),
                                                &m\text{ even}.
 \end{cases}                                                 \tag{14}
\]
The size in one fixed ordered original pattern is (14) divided
by \(L\). Its absolute full orbit size is (14) divided by \(N\)
when \(n\ge3\), and twice (14) divided by \(N\) when \(n=2\).

**Proof.** First count tuples in a fixed original ordered pattern
that preserve a specified proper quotient with \(k\) blocks,
when both original fixed-point counts \(e,f\) are positive.
The quotient degree \(k\) is odd. Its two involutions can be
chosen in exactly \(k\) ways with product \(P_k\).
Each of the \(s\) simple factors can attach its unique return
to any of the \(k\) blocks. The two heavy factors have one
fixed quotient block each, with return types forced by \(e,f\).
The literal bijection L2 therefore gives
\[
                k^{s+1}
 C_{n/k}\left(\frac{n/k-e}{2},\frac{n/k-f}{2};s\right).
\]
The order of the active slots may change, but Hurwitz moves give
bijections between their type orders and preserve primitivity.

By B2, a tuple with maximal Chebyshev quotient degree \(M\)
has exactly the proper quotient degrees dividing \(M\).
Here \(M\) is odd. Summing \(\mu(k)\) over these divisors selects
exactly \(M=1\), proving (13). The singleton-block partition is
not included: every term has inner degree at least \(s+1\ge2\).
This argument also applies when one retained inner heavy position
is the identity or a transposition; it uses its specified position,
not an assumption that its type differs from all simple factors.

Now take the canonical quotient degree \(m=N/n\) of the original
component. For odd \(m\), there are \(m\) quotient pairs,
\(m^s\) choices of simple-factor blocks, and \(P_n(E,F;s)\)
primitive return tuples in each case.

For even \(m\), the prescribed fixed-point-free quotient reflection
gives \(m/2\) quotient pairs. The other reflection has two distinct
fixed blocks. The unordered primitive partition can be assigned
to them in two ways when \(e\ne f\) and one way when \(e=f\).
There are again \(m^s\) simple-factor attachments, and the
primitive return count in each assignment is \(P_n(e,f;s)\).
L2 is a literal bijection, so none of these choices is divided
by an unstated relabeling symmetry.

B1 puts each resulting canonical label into one component.
There are exactly \(L\) original type placements, with the same
count at every placement by braid bijections. This proves (14).
Finally B2 gives the exact cyclic stabilizer orders used in
the absolute count. \(\square\)

## Corollary B6. Primitive connectivity and a same-group infinite family

For every passport (2) with both \(N-2A,N-2B>0\), there is exactly
one primitive component, and its monodromy is \(S_N\).
If one bundled factor is fixed-point-free, there is no primitive
component. This settles the primitive \(S_N\) part of this entire
family, with any positive number of simple factors.

For every odd prime \(p\ge5\), the degree-\(2p\) passport
\[
           (2^p),\quad(2^{p-3}1^6),\quad T,\quad T            \tag{15}
\]
has exactly two full components. Both generate the identical
standard embedded \(S_p\wr C_2\), have block-count divisors
\(\{1,2,2p\}\), and have trivial deck group.
Their primitive degree partitions are \(\{1,5\}\) and
\(\{3,3\}\), and their literal sizes are respectively
\[
       \frac{p^3(p^2-9)(p^2-1)}{20},
               \qquad \frac{p^3(p^2-1)^2}{12}.               \tag{16}
\]
For \(p=5\), these are \(2400\) and \(6000\), with group order
\(28800\) and absolute sizes \(240\) and \(600\).

**Proof.** In the positive-fixed-point case, \(n=N\) is always
an admissible and unique primitive label in (4). The zero case
permits only even \(m\), hence no primitive label.

For (15), an even multiplier requires the inner degree to divide
the odd prime \(p\). The admissible inner degree is therefore
only \(n=p\). The positive odd partitions of \(6\) are exactly
\(\{1,5\},\{3,3\}\). B1--B2 prove the component and group
statements. In (13), every nontrivial prime divisor would leave
degree one and contribute zero. Thus \(P_p=C_p\).
Substituting \(m=2,s=2,L=12\) in (11)--(14) and simplifying
gives (16). \(\square\)

The partition in (10) is an explicit extra component invariant
within the same literal monodromy group. This result does not
identify it with a reduced Schur-multiplier or spin invariant.

## Exact controls and remaining scope

The [checker](verify_polynomial_involution_components.py),
[receipt](polynomial_involution_components_verification.json), and
[source ledger](polynomial_involution_components_source_ledger.json)
separate the imported all-degree completeness theorem, the written
counting and group arguments, and their exact bounded controls.
Ten passports are exhausted by disjoint ordinary Hurwitz closures
whose union equals an independent hook-character total. They
include both parities of \(m\), a non-square-free divisor case,
both components of (15) for \(p=5,7\), and identity inner returns.

The all-degree statement is not inferred from these finite runs.
Published Corollary 2 supplies component completeness in exactly
the family (1)--(2). General primitive \(S_N/A_N\) passports,
remaining larger-fibre return actions and reduced lifting
classification remain open. Independent review, independent
replay and formal verification of the present derivations are
not asserted.
