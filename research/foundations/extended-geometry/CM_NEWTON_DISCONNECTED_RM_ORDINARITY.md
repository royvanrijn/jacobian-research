# Ordinarity on cyclic endomorphism Frobenius classes

The coefficient field of real multiplication carries two different
Galois actions. One is the arithmetic action on an embedded number
field. The other is conjugation on endomorphisms of a fixed abelian
variety. A rational prime can split in the field for the first action
while Frobenius cyclically permutes the relative Tate-module factors
for the second action. On these Frobenius classes the permutation
forces the same Newton polygon in each factor. This works even when
the connected projective representations are independent.

## A split-prime density theorem

**Theorem.** Let \(B/k\) be a polarized abelian variety of dimension
\(d\), where \(k\) is a number field. Suppose that
\(\operatorname{End}^0(B_{\bar k})=F\) is a cyclic totally real field
of degree \(d\), that \(k\) contains a copy of \(F\) as a number
field, and that the conjugation action of \(G_k\) on
\(\operatorname{End}^0(B_{\bar k})\) surjects onto
\(\operatorname{Gal}(F/\mathbb Q)\). Let \(L/k\) be the field of
definition of the endomorphisms.

The set of good primes \(v\) of \(k\) of absolute degree one
such that \(\operatorname{Frob}_v\) generates
\(\operatorname{Gal}(L/k)\) and \(B_v\) is ordinary has natural
density \(\varphi(d)/d\). In particular the lower natural density
of ordinary primes is at least \(2/3\) for a cyclic cubic field.
At every good unramified prime \(v\) in these classes whose
rational residue characteristic splits completely in \(F\),
the reduction is either ordinary or supersingular.

**Connected monodromy.** Fix an auxiliary prime \(\ell\). Over
\(\overline{\mathbb Q}_\ell\) the Tate representation has \(d\)
two-dimensional factors \(V_\sigma\), indexed by embeddings of \(F\).
The identity component \(G^\circ\) of its algebraic monodromy
centralizes \(F\). Each factor has determinant the common cyclotomic
similitude \(\chi\). Indeed Rosati fixes the totally real field and
the alternating form is nondegenerate on each two-dimensional
factor. After passing to the connectedness field, Faltings's theorem
identifies the commutant with \(F\otimes\overline{\mathbb Q}_\ell\).
Thus the factors are irreducible and pairwise inequivalent by
semisimplicity. Their derived projections are \(\mathrm{SL}_2\).
Indeed a connected irreducible reductive subgroup of
\(\operatorname{GL}_2\) has nontrivial semisimple derived group,
which must be \(\mathrm{SL}_2\). A link between projective derived
projections would identify their projective \(G^\circ\)-actions,
since the connected center acts by scalars. The \(V_\sigma\)
would then differ by a scalar character, whose square is one
by equality of determinants. Connectedness would make the
character trivial, contrary to Faltings. By the subdirect
product lemma for simple Lie algebras, the derived group is
the product of \(d\) independent copies of \(\mathrm{SL}_2\).
[Faltings, Finiteness theorems for abelian varieties over number
fields, Theorems 3–4](https://doi.org/10.1007/978-1-4613-8655-1_2)
· [Projective trace field criterion](CM_NEWTON_PROJECTIVE_TRACE_FIELD.md).

**Cyclic Frobenius.** The action of \(G\) on \(F\) is constant
on each connected component and coincides on the image of
\(G_k\) with its conjugation action. Write \(q\) for the residue
cardinality of \(v\), and \(g=\rho_\ell(\operatorname{Frob}_v)\)
with the convention \(\chi(g)=q\). On a generator class it sends
the \(V_\sigma\) around one cycle of length \(d\). The restrictions
of \(g^d\) to these factors are mutually conjugate by powers
of \(g\). Write \(a_v=\operatorname{Tr}(g^d)/d\). Multiplying
the \(d\) off-diagonal block isomorphisms around the cycle
gives the characteristic polynomial
\[
             \det(X-g)=X^{2d}-a_vX^d+q^d.
\]
In particular \(a_v\) is a rational integer, and the Weil bound
on either \(g^d\)-block gives \(\lvert a_v\rvert\le 2q^{d/2}\).

Over the residue extension of degree \(d\) the endomorphisms act
on the reduction and commute with Frobenius. If \(p\) splits
completely in \(F\), the rational \(p\)-divisible group decomposes
by the \(d\) projectors in \(F\otimes\mathbb Q_p\). Each summand
has height two. In fact, the characteristic polynomial of a
generator of \(F\) agrees on étale and crystalline \(H^1\).
For an endomorphism \(u\) both characteristic polynomials
evaluate at an integer \(n\) to the degree of \(n-u\), by the
exterior algebra structure of abelian-variety cohomology and
the degree formula on top cohomology. This applies to a rational
endomorphism after clearing denominators. Thus each eigenvalue
has the same multiplicity two as on the Tate module. The Rosati
self-duality of each summand then forces dimension one. Hence
each slope pair is either \((0,1)\) or \((1/2,1/2)\).
The total characteristic polynomial of the residue Frobenius
after the degree-\(d\) extension is
\((X^2-a_v X+q^d)^d\). Consequently there are at most two
distinct slopes. A mixture of ordinary and supersingular
summands would have three, so all summands have the same
Newton polygon.
[Katz–Messing, Theorem 1](https://www.math.mcgill.ca/goren/SeminarOnCohomology/Katz-Messing.pdf).

**Density.** A degree-one prime of \(k\) has \(q=p\), and,
away from ramification, \(p\) splits completely in the Galois
field \(F\subset k\). In the supersingular case its common trace
is divisible by \(p^{\lceil d/2\rceil}\). If \(d\) is odd this
forces \(a_v=0\) for \(p>4\). If \(d\) is even it forces
\(a_v^2/p^d\in\{0,1,4\}\). In both cases these Frobenius
elements lie in the closed locus
\[
  \left(\frac{\operatorname{Tr}(g^d)}{d}\right)^2
                   \ \in\ \{0,\chi(g)^d,4\chi(g)^d\}.
\]
The locus is proper in each connected component mapping to a
generator. With \(g\) fixed in that component, multiply it by
an element \(h\) of one of the independent \(\mathrm{SL}_2\)
factors in \(G^\circ\). The similitude stays fixed and the trace
of the cyclic block product becomes \(\operatorname{Tr}(M h)\)
for an invertible two-by-two matrix \(M\). This trace is
nonconstant on \(\mathrm{SL}_2\), as is seen by multiplying
\(M^{-1}\) by a scalar and then by any diagonal determinant-one
matrix. Algebraic Chebotarev gives density zero for the closed
locus on the generator components. The generator classes have
density \(\varphi(d)/d\) by finite Chebotarev, and the omission
of higher absolute degree primes preserves it.
\(\square\)
[Serre, Quelques applications du théorème de densité de
Chebotarev, §5.2](https://www.numdam.org/item/PMIHES_1981__54__323_0/)
· [Sawin, Theorem 2.1](https://www.numdam.org/item/10.1016/j.crma.2016.01.025.pdf)
· [Cantoral Farfán et al., Proposition 4.8](https://arxiv.org/html/2508.11174v1#S4.SS2).

## Intersecting coefficient actions on a fixed variety

Suppose \(A\) is not potentially CM and admits a
Rosati-compatible rank-two CM field action. Write
\(D=\operatorname{End}^0(A_{\bar k})\) and \(Z\) for its center.
The rank-two isotypic argument makes \(D\) simple, and each
split two-dimensional coefficient summand of connected monodromy
is irreducible. For any rank-two CM field \(E\) acting on \(A\),
the centralizer of \(E\) in \(D\) is therefore \(E\). Indeed,
after extending coefficients it is the sum of the scalar
centralizers in these summands, by Faltings. Thus \(Z\subset E\)
and \(E\) is a maximal subfield of the central simple algebra
\(D/Z\). Put \(n=\sqrt{[D:Z]}=[E:Z]\).
[Rank-two isotypic argument](CM_NEWTON_THREE_BALANCED_GLOBAL.md).

At a good prime of a field defining all of \(D\), Frobenius
acts on each \(Z\)-component of the geometric-coefficient Tate
module as \(1\otimes M_\zeta\) on \(W_\zeta\otimes U_\zeta\),
where \(\zeta\) runs through the embeddings of \(Z\),
\(D\) acts on \(W_\zeta\) through a full \(n\)-by-\(n\)
matrix algebra, and \(U_\zeta\) has dimension two.
The dimension follows from the rank-two action of any maximal
subfield \(E\), which selects one copy of \(U_\zeta\) at
each of its \(n\) extending embeddings. In particular all
rank-two relative Frobenius polynomials over \(E\) have the
form
\[
       X^2-a_ZX+\delta_Z,\qquad a_Z,\delta_Z\in Z,
\]
where the coefficients are the same for every choice of \(E\).
Indeed, the relative characteristic polynomial over \(E\)
is defined over \(E\), and its image under any embedding of
\(E\) extending \(\zeta\) is the characteristic polynomial of
\(M_\zeta\). Equal conjugates over \(Z\) place both coefficients
in \(Z\). In particular all the relative trace norms are the
same integer \(N_{Z/\mathbb Q}(a_Z)^n\). Their analogous
signature and slope blocks are also the same, with repetition
indexed by embeddings of \(E\) over \(Z\). For the signature
this follows from the same matrix decomposition of rational
Hodge structures. For the slopes it follows from the compatible
relative Frobenius polynomial at every split place.

For a further coefficient action \(K\subset D\) of different
relative rank, let \(m_{\tau,\zeta}\) denote the multiplicity
of the embedding \(\tau:K\to\overline{\mathbb Q}_\ell\) on
\(W_\zeta\). The relative polynomial at \(\tau\) is exactly
\[
        \prod_{\zeta}\det(X-M_\zeta)^{m_{\tau,\zeta}}.
\]
Each such field therefore gives products of the same
underlying quadratic Frobenius polynomials. Changing the
rank-two CM field leaves the norms and split Newton data
unchanged. Using a higher-rank presentation can select
different polynomial combinations, and any such argument
still has to obtain a restriction on the sequence of the
underlying center blocks. This description concerns Frobenius
over the field where the full endomorphism algebra is defined.
On the disconnected classes of the preceding theorem the
additional operation is the cyclic permutation of the center
blocks themselves.

## Weight-one tensor constructions from the independent system

**Proposition.** Let \(B\) be a number-field abelian variety with
full geometric endomorphism field \(F\) totally real of degree
\(\dim B\). In rational polarizable Hodge structures, consider
the tensor category generated by \(H^1(B)\), the Tate structures,
and CM Hodge structures. Every abelian variety whose \(H^1\)
lies in this category is isogenous over \(\mathbb C\) to
a product of a power of \(B\) and a CM abelian variety.

**Proof.** Any object in question uses only a finite sum \(C\)
of CM Hodge structures. Write \(G\) for the Mumford–Tate
group of \(H^1(B)\oplus C\). Its derived group maps
isomorphically to the derived group for \(B\): projection onto
that group is surjective and the other projection is into a
torus. The derived group for \(B\) is
\(\operatorname{Res}_{F/\mathbb Q}\mathrm{SL}_2\).
Indeed it is contained there by the endomorphisms and
polarization. It contains that derived group after a coefficient
extension by Lombardo's large image theorem and Deligne's
inclusion of connected monodromy in the Mumford–Tate group.
The standard two-dimensional representation of each factor
in \(H^1(B)\) has Hodge cocharacter weights \(0,1\).
[Deligne, §3 and Main Theorem 2.11](https://jmilne.org/math/Documents/Deligne82.pdf)
· [Lombardo, Theorem 1.4](https://arxiv.org/pdf/1411.1703).

An irreducible constituent \(W\) over \(\overline{\mathbb Q}\)
of a weight-one object has a highest weight \(n_\sigma\ge0\)
on each standard \(\mathrm{SL}_2\) factor. On this constituent
the difference of extreme weights of the square of the Hodge
cocharacter is \(2\sum_\sigma n_\sigma\): modulo the center
the squared cocharacter comes from \(\operatorname{diag}(t,t^{-1})\)
on every factor. The difference is at most two for an
abelian \(H^1\). Hence \(W\) is trivial on the derived group
or is the standard representation of exactly one factor.

In the latter case, \(W\simeq V_\sigma\otimes\xi\), where
\(V_\sigma\) is the corresponding rank-two constituent of
\(H^1(B)\), and \(\xi\) is the character on the one-dimensional
derived intertwiner space. Both \(W\) and \(V_\sigma\) have
Hodge weights \(0,1\), so \(\xi\) has weight zero under the
Hodge cocharacter. The same is true for every conjugate
of \(\xi\), since both weight-one objects are rational.
The intersection of the kernels of these conjugate characters
is a rational subgroup containing both complex factors of
the Hodge morphism. By the defining minimality of \(G\),
this subgroup is \(G\). Thus \(W\simeq V_\sigma\) as
\(G\)-representations.

An occurrence of a \(V_\sigma\) gives a nonzero rational
Hodge morphism from \(H^1(B)\), by base change of the
intertwiner space. Polarizable Hodge structures are semisimple,
and \(H^1(B)\) is simple since its endomorphism algebra is
a field. Splitting off such copies leaves a structure whose
Mumford–Tate group is a torus. The correspondence between
polarizable weight-one structures and complex abelian varieties
up to isogeny, together with the torus criterion for CM,
gives the asserted isogeny. \(\square\)

In particular a new auxiliary abelian variety obtained this
way has no new non-CM slope constraints after finite base
extension. Among independent real-cubic systems the Hodge
cocharacter is noncentral on all three factors. This is the
weight obstruction to bringing the odd tensor representation
in Noot's argument into weight one here.
[Noot, Theorem 2.2](https://www.numdam.org/item/CM_1995__97_1-2_161_0.pdf).

## A fixed cubic example

Cunningham and Dembélé construct an abelian threefold \(B_1/\mathbb Q\)
with \(\operatorname{End}^0_{\mathbb Q}(B_1)=\mathbb Q\). Over the
cyclic cubic subfield \(L_{43}\subset\mathbb Q(\zeta_{43})\) it has
real multiplication by the cyclic cubic subfield \(F_{19}\)
of \(\mathbb Q(\zeta_{19})\). The latter field is generated by a root
\(e\) of \(x^3-x^2-6x+7\). These two cyclic fields are linearly
disjoint because their conductors differ. The proposition gives the
statement over \(L_{43}\), and the following reduction check rules
out extra geometric endomorphisms.
[Cunningham–Dembélé, Proposition 6.3 and Table 2](https://arxiv.org/pdf/1705.03054#page=21).

At primes of \(L_{43}\) of norm \(2\) and \(11\), Table 2 gives
relative traces \(e\) and \(-e^2\), respectively. Both reductions
are ordinary. Their rational characteristic polynomials have
middle coefficients \(3\) and \(335\), which are units at the
respective primes. The two relative discriminants are
\[
        D_2=e^2-8,\qquad D_{11}=e^4-44.
\]
Their norms to \(\mathbb Q\) are \(-31\) and \(-7\cdot449\).
Modulo \(31\), the three roots of the field polynomial are
\(4,12,16\), and \(D_2\) vanishes only at \(16\). Modulo \(7\)
they are \(0,3,5\), and \(D_{11}\) vanishes only at \(5\).
The valuations at the selected places are one by the norm
computations. In either case the three conjugate discriminants
therefore have independent square classes, witnessed by the
three valuations at the conjugate places. The splitting field
over \(F_{19}\) contains independent swaps of either root of
each relative quadratic. The ratio of the two roots in a pair
is nontorsion since the slopes are ordinary and their product
is \(p\), so one is a unit and the other has valuation one
at each place inducing such a pair. The ratio of roots from
different pairs cannot be torsion either, since swapping only
one of these pairs would then make its within-pair ratio
torsion. Thus all six Frobenius roots have distinct \(m\)th
powers for every \(m\geq1\).

For \(p=2,11\), the trace \(e\) or \(-e^2\) generates \(F_{19}\),
and for a relative root \(\pi\) the equality
\(\pi+p/\pi=a_v\) puts \(F_{19}\) in \(\mathbb Q(\pi)\).
The nontrivial discriminant extension gives degree six, and
the distinct-powers check proves degree six also for
\(\mathbb Q(\pi^m)\) for every \(m\). The finite-field theorem
of Tate gives a commutative geometric endomorphism field
\(F_{19}(\sqrt{D_p})\) for each reduction. The field from \(p=2\)
ramifies at \(31\). The field from \(p=11\) is unramified above
\(31\) since \(31\) is unramified in \(F_{19}\) and \(D_{11}\)
is a \(31\)-adic unit. Specialization of endomorphisms embeds
the geometric endomorphism algebra of \(B_1\) into each of
these degree-six fields. Since it contains \(F_{19}\) and the
two fields are nonisomorphic, this algebra must be \(F_{19}\).
[Tate, Théorème 1](https://www.numdam.org/item/SB_1968-1969__11__95_0/).

Now base change \(B_1\) to \(k=F_{19}\) as a number field.
Its endomorphism field over \(k\) is \(F_{19}L_{43}\), since
\(L_{43}/\mathbb Q\) is cubic and the endomorphisms over
\(\mathbb Q\) were only rational scalars. The theorem gives an
ordinary subset of natural density \(2/3\), all of absolute
degree one, over \(F_{19}\). This statement identifies a set
of split ordinary primes even though the usual inert coefficient
prime argument is unavailable over \(F_{19}\).

## The identity class and comparison with ordinary-reduction results

For comparison, write \(a_w\in F\) for the \(F\)-relative trace
at a good prime \(w\) above a prime of \(k\) in the endomorphism
field \(L\). Its semilinear symmetry is
\[
             a_{\sigma w}=\sigma(a_w),\qquad
                    \sigma\in\operatorname{Gal}(L/k),
\]
using the identified conjugation action on \(F\). It follows
by conjugating Frobenius by a lift of \(\sigma\), and is the
compatibility in Cunningham–Dembélé's Theorem 4.1. If
\(\operatorname{Frob}_v\) generates \(\operatorname{Gal}(L/k)\),
all \(\sigma\) fix \(w\), which is exactly the rational trace
calculation in the theorem. If \(v\) splits completely in
\(L\), the same formula relates the traces at *distinct*
primes above \(v\). For a trace divisible by one place of
\(F\) over \(p\), the relation transfers the divisibility to
the conjugate place of \(F\) at the conjugate prime of \(L\).
It supplies no additional divisibility at the original prime
of \(L\). The connected projective functions remain
independent throughout.
[Cunningham–Dembélé, Theorem 4.1](https://arxiv.org/pdf/1705.03054#page=13).

The fixed threefold \(B_1\) realizes this distinction at \(p=563\).
The prime splits in \(F_{19}\) and \(L_{43}\). In the former field
the residues of \(e\) are \(72,113,379\). In the latter, take \(w\)
with \(w^3-w^2-14w-8=0\). Its residues are \(145,157,262\).
The level-one Hilbert form in Cunningham–Dembélé's Proposition 6.3
is identified by the unique dimension-three Hecke constituent
over \(L_{43}\), with coefficient polynomial
\(e^3-e^2-6e+7\). Its Hecke data give the following three relative
Frobenius traces. The last three columns evaluate the trace at
the three coefficient places over \(563\).

| \(w\bmod 563\) | Relative trace | \(e=72\) | \(e=113\) | \(e=379\) |
| ---: | --- | ---: | ---: | ---: |
| \(262\) | \(13e^2+2e-57\) | \(482\) | \(81\) | \(0\) |
| \(157\) | \(-2e^2+11e+5\) | \(0\) | \(482\) | \(81\) |
| \(145\) | \(-11e^2-13e+52\) | \(81\) | \(0\) | \(482\) |

The three prime generators in the data download are, in table
order, \(-9w^2+35w+23\), \(13w^2-21w-171\), and
\(-4w^2-14w-9\), each vanishing at exactly the specified root
of the base field polynomial modulo \(563\). The eigenvalue
\(a=13e^2+2e-57\) has first, second and third elementary symmetric
sums \(0,-931,-10697\), by multiplication in
\(\mathbb Q[e]/(e^3-e^2-6e+7)\). The last sum is
\(-19\cdot563\), so the unique bad coefficient place has
valuation exactly one. Expanding the norm gives
\[
\begin{aligned}
 N_{F_{19}/\mathbb Q}(X^2-aX+563)
 &=X^6+758X^4+10697X^3+426754X^2+178453547\\
 &\equiv X^4(X^2+195)\pmod {563}.
\end{aligned}
\]
This is a split degree-one good reduction of \(B_1\) with
\(p\)-rank two. The relative trace table shows the unique nonordinary
factor cycling through the three coefficient places as the base
prime cycles. The displayed calculation can be replayed by forming
the characteristic polynomial of multiplication by \(a\) in
the basis \(1,e,e^2\), and substituting \(X^2+563\) for its indeterminate
after homogenizing with \(X\). The [exact replay](verify_cm_newton_split_rm_example.py)
also checks the two ordinary reductions used to identify the geometric
endomorphism field. The checked LMFDB download has SHA-256
7b106138d75b819e24accd920bbc7417885bd0653274f5cb3769a180658f4e18.
[Cunningham–Dembélé, Proposition 6.3](https://arxiv.org/pdf/1705.03054#page=21)
· [LMFDB eigenvalue data for 3.3.1849.1-1.1-c, Sage download](https://www.lmfdb.org/ModularForm/GL2/TotallyReal/3.3.1849.1/holomorphic/3.3.1849.1-1.1-c/download/sage).

For a rank-two CM example on the same variety up to taking a
square, let \(A=B_1^2\) over \(L_{43}F_{19}(\sqrt{-2})\) and let
\(E=F_{19}(\sqrt{-2})\) act by diagonal real multiplication and
by the matrix \(J=\left(\begin{smallmatrix}0&-2\\1&0\end{smallmatrix}\right)\)
for \(\sqrt{-2}\). The matrix is Rosati-skew for the product
polarization with factors of weights one and two. The rank-two
signature has three balanced pairs, one for each real embedding
of \(F_{19}\). The prime \(563\equiv3\pmod8\) splits in
\(\mathbb Q(\sqrt{-2})\), hence splits in the full base and in
\(E\). On every such reduction \(A\) has \(p\)-rank four,
and the one nonordinary balanced pair is the one selected by
the corresponding trace of \(B_1\). Its three connected projective
factors are independent by the threefold calculation. This is a
fixed characteristic-zero variety exhibiting the one-failure
Newton polygon allowed by the [three-pair bound](CM_NEWTON_THREE_BALANCED_GLOBAL.md).

Fité's Theorem 3 and Wang–Zhang's Theorem 1.3 use an inert
coefficient prime or, in Wang–Zhang, an alternative involving
two distinct equal-degree coefficient places. The cubic
calculation here uses a split prime instead. In the
Cunningham–Dembélé example the theorem of Fité already gives
positive density over \(\mathbb Q\), since the two cubic fields
are disjoint. Its noncontainment hypothesis cannot hold over
\(F_{19}\) or an extension of \(F_{19}\). For Wang–Zhang's
degree-one inert or equal-degree mechanisms in this example
over the field of definition of the action, the only
degree-one coefficient splitting is complete splitting into
three places. Suh's Theorem 4.1.1 concerns the arithmetic
action on embeddings of a Hecke field for representations over
the field of definition. It does not identify the cyclic
endomorphism-action classes above. Cantoral Farfán et al.
provide the componentwise Chebotarev method under their
CM simple-signature hypotheses. The argument here uses full
geometric endomorphism algebra a totally real field and the
cyclic block determinant. Priority for this specific
application is not established by the sources checked.
[Fité, Theorem 3 and §6](https://arxiv.org/html/2005.10185v3)
· [Wang–Zhang, Theorem 1.3](https://doi.org/10.1112/mtk.70081)
· [Suh, Theorem 4.1.1](https://arxiv.org/html/2410.01182v2)
· [Cantoral Farfán et al., Theorem 4.18](https://arxiv.org/html/2508.11174v1#S4.SS5).

The positive-density theorem is over the field where the
endomorphism action is nontrivial. It leaves the identity
Frobenius class open. On passing to the endomorphism field,
the primes detected by generator Frobenius classes have
absolute degree divisible by \(d\) and have density zero if
\(d>1\). Thus the theorem yields no potential density-one
statement for the independent split cubic case.
