# Orbit-weighted Newton inequalities in arbitrary relative rank

Let \(A/K\) be a polarized abelian variety of dimension \(eh\) over a
number field. Suppose the geometric endomorphism algebra contains a CM
field \(E\) of degree \(2e\), preserved as a subalgebra by \(G_K\), with
Rosati acting by complex conjugation. Write \(F=E^+\), and suppose
\(h\ge2\) and \(A\) is not potentially CM. The action of \(E\) need not
be defined over \(K\).

At an absolute degree-one good place \(v\) of residue characteristic
\(p\), exclude ramification in \(E\), primes dividing the order index
and the finite exceptional places of the action. A coefficient place
\(\lambda\mid p\) of \(E\) has residue degree \(f_\lambda\). Write
\(s_{\lambda,1}\le\cdots\le s_{\lambda,h}\) for the relative slopes of the
geometric \(E_\lambda\)-summand, normalized by \(v_p(p)=1\). Each slope
is repeated \(f_\lambda\) times in the underlying isocrystal. Set
\(S_{\lambda,j}=\sum_{k=1}^j s_{\lambda,k}\).

**Theorem.** There is an integer \(M\ge1\), independent of \(v\), for
which, outside a set of natural density zero over the original field
\(K\),
\[
 \sum_{\lambda\mid p} f_\lambda
             \left\lceil M S_{\lambda,j}\right\rceil
 \ \le\ eMj-1 \qquad(1\le j<h).
\]
Any positive multiple of the \(M\) constructed below also works, with
its own density-zero exception. One may take \(M=1\) on the primes
whose auxiliary algebraic Frobenius lies in the identity component.

## A Frobenius power on each monodromy component

Fix an auxiliary prime and write \(G\) for the algebraic monodromy group
on \(V_\ell(A)\), working over \(\overline{\mathbb Q}_\ell\). Faltings's
semisimplicity gives connected reductive \(G^\circ\). Its derived group
\(D\) is nontrivial. Indeed, if \(G^\circ\) were a torus, then after
passing to the connectedness field the Tate centralizer would contain
an element with distinct eigenvalues. Over the algebraic closure it
contains the diagonal matrices in a weight basis, so this is a
nonempty open condition in the centralizer. Faltings's endomorphism
theorem identifies that centralizer with geometric rational
endomorphisms after extension of scalars. A rational endomorphism in
this open condition would generate a commutative semisimple algebra
of degree \(2\dim A\), giving potential CM.

There is a finite-order representative \(g_0\) of each component of
\(G\) whose conjugation on \(D\) preserves a pinning. Start with any
representative \(g\) and adjust it by an element of \(D\) to make
\(\operatorname{Int}(g)|_D\) preserve the chosen pinning. This is the
decomposition of semisimple-group automorphisms into inner and pinned
automorphisms. The pinned action permutes the simple coroots. In
particular it fixes a nonzero cocharacter \(\nu\) obtained by summing
the coroots in an orbit. [SGA 3, Exposé XXIV, §1](https://jcreinhold.github.io/sga/iii/24-automorphisms-of-reductive-groups.html).

To arrange finite order, let \(n\) be the order of the component of the
adjusted \(g\). Its \(n\)-th power belongs to \(G^\circ\). The induced
pinned automorphism of \(D\) is then inner, hence trivial, because an
inner automorphism preserving a pinning is trivial. Thus \(g^n\)
centralizes \(D\), and in \(G^\circ\) it belongs to the center. It is
also fixed by conjugation by \(g\). The fixed subgroup
\(Z(G^\circ)^g\) has finitely many components and torus identity
component. Choose a multiple \(m\) of \(n\) with \(g^m\) in this torus.
An \(m\)-th root \(z\) of \(g^{-m}\) in the same torus exists over the
algebraic closure. Then \(g_0=gz\) has \(g_0^m=1\), induces the same
action on \(D\), and centralizes \(\nu(\mathbb G_m)\).

Take \(M\) a common multiple of the orders of the representatives.
In particular \(g^M\in G^\circ\) for every \(g\in G\), and it
centralizes \(E\). Let \(V_\sigma\) be the relative block for an embedding
of \(E\), and \(a_{j,\sigma},\delta_\sigma\) its \(j\)-th elementary
coefficient and determinant. The function
\[
 H_j(g)=\prod_{\sigma:E\hookrightarrow\overline{\mathbb Q}_\ell}
          \frac{a_{j,\sigma}(g^M)^h}{\delta_\sigma(g^M)^j}
\]
is regular and invariant under conjugation. The denominators are units
and conjugation just permutes the relative blocks.

On the component of \(g_0\), vary \(g=g_0\nu(z)\). Then
\(g^M=\nu(z^M)\). The determinants in the relative blocks are one on
\(D\), and their eigenvalues along this cocharacter are integral
powers of \(z\). In at least one block the weights are not all equal.
Otherwise the determinant-one condition makes every weight zero,
contrary to faithfulness of \(V_\ell(A)\). For \(1\le j<h\) the
elementary \(j\)-th coefficient in that block has at least two distinct
powers of \(z\), since one can exchange two unequal weights in a subset
of size \(j\). Evaluate these Laurent polynomials at a generic complex
number \(z\) of modulus one. Every coefficient has absolute value at
most \(\binom hj\), and at least one has strictly smaller value by the
triangle inequality. At \(z=1\) every coefficient equals \(\binom hj\).
This proves nonconstancy on every component.

## Valuations and coefficient places

At \(v\), the power \(\pi_v^M\) of Frobenius centralizes the specialized
action of \(E\). Let \(a_j(v),\delta(v)\in E\) be its relative elementary
coefficient and determinant. The degree polynomial of
\(a\pi_v^{Mk}\) for \(a\in E\) recovers the relative power traces,
and Newton identities give the same relative characteristic
polynomial in the Tate and crystalline realizations. As in the
[simultaneous Newton inequality](../../../unpublished/cm-newton-polygon-constraints/main.tex),
the coefficients are integral, since the roots are algebraic integers.
Every complex conjugate root has absolute value \(p^{M/2}\), and the
determinant norm is \(p^{ehM}\). Consequently
\[
 H_j(\operatorname{Frob}_v)
             =\frac{N_{E/\mathbb Q}(a_j(v))^h}{p^{ehMj}},
 \qquad
 \left|H_j(\operatorname{Frob}_v)\right|\le \binom hj^{2eh}.
\]
If \(v_p(N_E(a_j(v)))\ge eMj\), this bounded rational number is an
integer, since its denominator is a power of \(p\). Its finitely many
integer levels, including zero, contain no component. Algebraic
Chebotarev excludes them on a density-zero set. There are only \(h-1\)
indices, so the resulting strict inequality holds simultaneously.
[Rajan, Theorem 3](https://repository.ias.ac.in/38317/1/2-PUB.pdf).

The \(M\)-fold crystalline Frobenius on the residue extension where the
action is defined is the linear action of \(\pi_v^M\). If necessary
extend the unramified coefficients to split \(E_\lambda\). The
semilinear Frobenius cycles the \(f_\lambda\) idempotents, which makes
their slope lists identical after normalization. The root valuations in
each relative block of \(\pi_v^M\) are \(M s_{\lambda,k}\). Its
elementary coefficient thus has valuation at least \(M S_{\lambda,j}\)
and has integral valuation over unramified \(E_\lambda\). Finally,
\(v_p(N_E(a_j))=\sum_\lambda f_\lambda v_\lambda(a_j)\). Combining
this with the strict norm inequality proves the theorem.

On \(G^\circ\) one can use \(g_0=1\) directly and take \(M=1\).
A second setting for \(M=1\) is an \(E\)-action defined over \(K\)
and a projection of \(G^\circ\) onto the full product of
\(\mathrm{PGL}_h\) over one oriented embedding in each conjugate pair.
Then on each component multiplication by \(G^\circ\) can make all
these projective blocks scalar simultaneously. Polarization gives the
same conclusion on the conjugate blocks. Varying the projective
blocks from this scalar point proves nonconstancy without raising
Frobenius.

**Proposition.** Suppose \(h=3\), the action of \(E\) is defined over
\(K\), and \(E=\operatorname{End}^0(A_{\overline K})\). If some
conjugate pair is mixed, one can take \(M=1\) on every monodromy
component, including when different pairs are projectively linked.

Take an auxiliary prime split in a normal closure of \(E\). The
geometric endomorphism hypothesis and Faltings's theorem make every
three-dimensional \(V_\sigma\) absolutely irreducible under
\(G^\circ\). The connected center acts by scalars, so the image of
\(D\) in this block is an irreducible connected semisimple subgroup
of \(\mathrm{SL}_3\). Its possible Lie types are \(A_1\) and \(A_2\):
rank is at most two and dimension at most eight, and a nontrivial
irreducible representation of two \(A_1\) factors needs dimension
at least four. In type \(A_1\) the representation is the symmetric
square of the two-dimensional representation, and is self-dual
on \(D\). In that case
\(\operatorname{Sym}^3 V_\sigma\otimes\det(V_\sigma)^{-1}\)
is self-dual under \(G^\circ\), since its connected central character
is trivial.

For each fixed \(\sigma\), extend the base to contain the normal
closure \(N\) of \(E\). Choose a mixed embedding \(\tau:E\hookrightarrow
N\), then a place given by an embedding \(\iota:K'\hookrightarrow
\overline{\mathbb Q}_\ell\) of the extended base with
\(\iota\circ\tau=\sigma\). One obtains \(\iota\) by extending the
prescribed map on \(\tau(E)\) to \(K'\). The comparison at this place
gives the chosen mixed Hodge type on the same block \(V_\sigma\).
Connected monodromy is unchanged by the finite extension.
At a mixed type the extreme weights of the displayed
representation are \(-m_\sigma\) and \(3-m_\sigma\), for
\(m_\sigma\in\{1,2\}\). They are asymmetric about zero, excluding
self-duality. Thus every relative projective block has full type
\(A_2\) projection.

Each block sees precisely one adjoint simple factor of \(D\), and
every factor is seen since the representation is faithful. A
component representative preserves every block because the action is
defined over \(K\). It therefore preserves the factor seen by each
block. On such a factor a nontrivial diagram automorphism exchanges
standard and dual, which are distinct in rank three. Conjugation by
the representative acts within each block, so it induces an inner
automorphism on every factor. The adjoint conjugators lift together
through \(D\to D_{\rm ad}\). An automorphism trivial on the adjoint
quotient is trivial on connected \(D\), since its difference from the
identity is valued in the finite center. The adjusted representative
therefore centralizes \(D\) and acts by scalars on every block. Varying a
derived cocharacter through this scalar representative proves
nonconstancy of the normalized coefficient products with \(M=1\),
as in the first part of the proof. The norm argument then applies on
all components.

## Averaged Hodge polygons and compression

Let \(m_\sigma\) count the Hodge weights equal to one, and let
\(\theta\) be a Frobenius permutation on the \(E\)-embeddings in a
normal closure. A coefficient place corresponds to a \(\theta\)-orbit.
The prefix of its averaged Hodge polygon is
\[
 B_{\lambda,j}=\frac{1}{f_\lambda}
       \sum_{\sigma\ {\rm in\ the\ orbit}}\max(0,j-h+m_\sigma).
\]
In one step of the corresponding Dieudonné cycle, \(\bigwedge^j\)
of Frobenius on the lattice is divisible by
\(p^{\max(0,j-h+m_\sigma)}\). Compose around the orbit. The lowest
Frobenius slope on \(\bigwedge^j\) is the sum of the \(j\) lowest
slopes, so this divisibility gives \(S_{\lambda,j}\ge B_{\lambda,j}\).
The same computation on the determinant gives equality of the
endpoints. Write \(\Delta_{\lambda,j}=S_{\lambda,j}-B_{\lambda,j}\).
This is the comparison with the usual \(\mu\)-ordinary reference
polygon. [Cantoral Farfán–Li–Mantovan–Pries–Tang, Proposition 3.6](https://arxiv.org/html/2508.11174v1).

Let \(E_0\) be the connected projective spectrum field,
\(F_0=E_0\cap F\), and \(d=[F:F_0]\), as in the
[projective-spectrum compression](CM_NEWTON_PROJECTIVE_SPECTRUM_FIELDS.md).
Oriented embeddings with the same restriction to \(E_0\) have the
same connected projective spectra. For any two such embeddings, the
Hodge–Tate argument with \(\operatorname{Sym}^h V\otimes\det^{-1}V\)
shows that the two Hodge multiplicities agree if mixed, or that both
are definite. Thus their sorted Hodge vectors differ by a constant.
Their simultaneous \(\theta\)-translates also restrict identically
to \(E_0\). Averaging along the two orbits gives Hodge vectors
still differing by a constant, even when their orbit lengths differ.
As \(\operatorname{Frob}_v^M\in G^\circ\), their actual Newton slope
vectors differ by a constant by the projective spectrum comparison.
Equality of Hodge and Newton endpoints identifies these two constants.
It follows that the excess \(\Delta_{\lambda,j}\) is the same on each
oriented projective class. The conjugate excess at \(j\) is the
original excess at \(h-j\).

For any conjugate projective class \(C\) choose an orientation and
write its excess as \(\Delta_{C,j}\). Write \(m_C\) for the smaller
original Hodge multiplicity, constant on the class. Assume \(M\) is
chosen so that every \(MB_{\lambda,j}\) is integral. This can always
be arranged by making \(M\) a multiple of the possible residue
degrees. It holds with \(M=1\) when the oriented Hodge multiplicity
is constant on each coefficient place. The theorem then gives
\[
 \sum_C
     \left(\left\lceil M\Delta_{C,j}\right\rceil
              +\left\lceil M\Delta_{C,h-j}\right\rceil\right)
 \ \le\ M D_j-1,\qquad
 D_j=\sum_C\min(j,h-j,m_C).
\]
In fact the sum of the averaged Hodge prefixes over all embeddings
equals the sum before averaging. The contribution of a conjugate pair
to that sum is \(j-\min(j,h-j,m_C)\). Subtract this baseline in the
norm inequality and divide by \(d\). The \(-1\) after division follows
from integrality.

Let \(\gamma\) be the endomorphism permutation of the original
Frobenius, acting by precomposition on the labels. The automorphism
\(\gamma\) preserves connected projective equivalence by conjugation
on \(G^\circ\). It commutes with \(\theta\), which acts by
postcomposition. Frobenius on the Hodge space of the geometric special
fiber preserves multiplicities under the composite
\(\theta\gamma^{-1}\), by applying it to the coefficient eigenvector
equation. Hence \(\gamma\) preserves the averaged Hodge vectors along
coefficient orbits. The quasi-isogeny \(\pi_v\) commutes with
\(\pi_v^M\) and permutes relative blocks, so \(\gamma\) preserves the
actual Newton data. The pair summands in the compressed inequality
are therefore constant on joint \(\theta,\gamma\)-orbits. One can sum
over their representatives with coefficients equal to the orbit
lengths.

## A numerical ordinarity criterion at nonsplit places

Consider \(M=1\) at an \(E\)-action-defined degree-one place. The
oriented multiplicities are constant along \(\theta\)-orbits because
the action is already defined over the prime field. The relative
Frobenius is linear over \(E_\lambda\), so its polygons have integral
slope-change vertices. At each class the reference polygon is
ordinary. Let \(a\) be the length of a \(\theta\)-orbit of mixed
conjugate classes and \(m>0\) its common smaller multiplicity.

The smallest numerical penalty from failing one such orbit is
obtained by replacing one zero and one unit slope by two halves in
each component and its conjugate. It is \(a\) at both \(j=m\) and
\(j=h-m\) if \(m<h/2\), and \(2a\) at \(j=h/2\) in the balanced
case. This is the pointwise minimum from the
[extremal polygon calculation](CM_NEWTON_EXTREMAL_POLYGONS.md).
It is compatible with reversal when the orbit is self-conjugate,
since homogeneity then forces the balanced signature. In particular,
if every mixed orbit satisfies
\[
 \begin{cases}
 a\ge D_m,&m<h/2,\\
 2a\ge D_{h/2},&m=h/2,
 \end{cases}
\]
the compressed inequalities force ordinary reduction. Conversely,
an orbit failing the displayed test admits the stated numerical
perturbation by itself, so this is the exact criterion for these
inequalities at the given homogeneous orbit type.

Equivalently the forcing orbit types are the following.

- There are no mixed orbits.
- There is one mixed orbit of multiplicity one, in any relative rank.
- The relative rank is two and there are two equal-length mixed orbits.
- The relative rank is four and there is one balanced orbit.

An unbalanced orbit with \(m\ge2\) fails the test because \(D_m\ge am>a\).
An unbalanced orbit with \(m=1\) needs \(a\ge D_1\), the number of
all mixed classes, so it must be the only mixed orbit. If every mixed
orbit is balanced, write \(h=2H\) and let \(B\) be the number of mixed
classes. Each mixed orbit must have \(2a\ge BH\), and their lengths
sum to \(B\). For \(H=1\) this permits one orbit or two equal-length
orbits. For \(H=2\) it permits only one orbit. For \(H>2\) no
nonempty mixed signature satisfies the test.

For example, take relative rank three and two reduced mixed classes
of multiplicity one in a single nonsplit \(\theta\)-orbit of length
two. Here \(D_1=D_2=2\). A failure on the orbit costs two at both
indices, and the bound permits at most one. At the split identity
type, the corresponding single-class failure costs one and survives
both inequalities. The following construction realizes this type
and its longer-orbit counterparts on fixed abelian varieties.

## Absolutely simple examples with a nonsplit ordinary class

**Corollary.** For each \(e\ge2\) and \(h\ge3\) there is an absolutely
simple abelian variety \(A/K\) of dimension \(eh\) whose full geometric
rational endomorphism algebra is \(E=F\mathbb Q(i)\), with \(F\)
totally real of degree \(e\). The action is defined over \(K\).
Every relative rank-\(h\) conjugate pair has signature
\((1,h-1)\), and its projective class is distinct from the other
\(e-1\) classes. The set of absolute degree-one primes of \(K\)
whose rational characteristic is inert in \(F\) has positive
natural density. At a relative density-one set of these primes,
\(A\) has ordinary reduction although the characteristic is nonsplit
in \(E\). In rank four there are also examples with signature
\((2,2)\) at every pair, with the same conclusion. If \(e\) is prime,
one can choose a cyclic \(F\) for which the inert class has density
\(1-1/e\). In that case the lower natural density of ordinary reduction
is at least \(1-1/e\).

**Proof.** Write \(m=1\) for the first signatures and \(m=2\) for the
second. Take \(I=\mathbb Q(i)\) and \(B:y^2=x^3-x\) over \(I\). Fix
full level \(N\ge3\), put \(k=I(B[N])\), and choose a rational prime
\(p\nmid 2N\) split completely in \(k\).

A totally real \(F\) of degree \(e\) in which \(p\) is inert and
the primes dividing \(2N\) are unramified can be chosen as follows.
Prescribe a monic irreducible degree-\(e\) reduction at \(p\), and
a monic squarefree degree-\(e\) reduction at every prime dividing
\(2N\). For large positive integer \(T\), perturb the lower
coefficients of \(\prod_{j=1}^e(X-jT)\) by bounded integers to hit
the prescribed residues. This is possible by the Chinese remainder
theorem. After replacing \(X\) by \(TX\) and dividing by \(T^e\),
the polynomial approaches \(\prod_{j=1}^e(X-j)\), which has simple
real roots. For large \(T\) the perturbed polynomial has \(e\)
distinct real roots. The reductions give the required inertness
and unramified conditions. In particular \(F\) has odd discriminant,
so \(\mathcal O_E=\mathcal O_F\mathcal O_I\).

In each of \(h\) blocks \(B^e\) let \(\mathcal O_F\) act by its
regular action in an integral basis. Its positive definite trace
form makes this action self-adjoint, so the trace Gram matrix
defines a polarization using the product polarization. Let
\(\mathcal O_I\) act via its given action on \(m\) blocks and
via conjugation on the other \(h-m\). This gives \(A_0=B^{eh}\)
with the required polarized \(E\)-action and signature, defined
over \(k\). The degree of the polarization is prime to \(N\)
because the trace form has determinant prime to \(N\).

Take the characteristic-zero PEL moduli component \(S/k\) through
this point with full level \(N\). The polarization type and
component may depend on \(F\). The level at \(A_0\) is defined
over \(k\), since \(A_0[N]\) is rational there. The signature
reflex field is contained in \(I\), because the multiplicities are
constant over each embedding of \(I\). It is \(I\) if \(m\ne h-m\)
and \(\mathbb Q\) in the balanced case. The geometrically connected
component through the \(k\)-point descends to \(k\), is smooth,
and carries the universal abelian scheme at fine level. Its
complex uniformization is an arithmetic quotient of a product of
\(e\) unitary symmetric domains of signature \((m,h-m)\).
The arithmetic level subgroup is Zariski dense in the derived
group by Borel density, since its simple real factors are
noncompact. Comparison with étale monodromy gives the \(e\)
independent \(\mathrm{SL}_h\) factors in geometric algebraic
monodromy after extension of scalars.
[Kottwitz, §5](https://doi.org/10.1090/S0894-0347-1992-1124982-1).
[Noot, §1.8](https://www.numdam.org/item/CM_1995__97_1-2_161_0/).

To specialize while preserving both this monodromy and a degree-one
place over \(p\), let \(\Gamma\) be the compact arithmetic
\(\ell\)-adic image of the family for an auxiliary \(\ell\ne p\).
Take an open normal pro-\(\ell\) subgroup \(U\) of \(\Gamma\).
Its Frattini subgroup \(\Phi(U)\) is open and characteristic, since
\(U\) is a finitely generated pro-\(\ell\) group. A closed subgroup
surjecting onto \(\Gamma/\Phi(U)\) intersects \(U\) in a subgroup
surjecting onto \(U/\Phi(U)\). Its intersection is therefore \(U\)
by the pro-\(\ell\) Burnside basis theorem, and the closed subgroup
is all of \(\Gamma\). Use the connected finite étale cover
corresponding to \(\Gamma/\Phi(U)\).

Choose an étale chart in \(S\) over affine space and shrink the
target so that both the chart and the composed cover are finite
étale over it. A nonempty Zariski open of the chart meets a small
\(p\)-adic neighborhood of the original smooth \(k_p=\mathbb Q_p\)
point. At such a point the chart has a local section over a
\(p\)-adic open ball in its target. Hilbert irreducibility with
weak approximation gives a \(k\)-rational parameter in that ball
whose fiber in the composed cover is connected. The fiber in the
chart is then a single closed point \(\operatorname{Spec} K\).
Its lift from the local section gives \(K\hookrightarrow\mathbb Q_p\).
The connected cover fiber gives the full image modulo
\(\Phi(U)\), hence the full \(\ell\)-adic image.
In detail, the disconnected fibers of the finite separable composed
cover are contained in a thin set, by considering the proper
subgroups of the permutation group of its Galois closure that act
intransitively on the generic fiber. Theorem 3.5.3 excludes the closure
of this thin set on a product of places away from the chosen \(p\)-adic
constraint. Weak approximation in affine space then supplies the
required parameter in the original ball.
[Noot, Proposition 1.3](https://www.numdam.org/item/CM_1995__97_1-2_161_0/).
[Serre, *Topics in Galois Theory*, Theorem 3.5.3](https://indico.math.cnrs.fr/event/11410/attachments/4760/7351/SerreTopicsGaloisTheory.pdf).

The resulting \(A/K\) has one standard \(\mathrm{SL}_h\) block
and one dual block for each derived factor. Standard and dual
are distinct for \(h\ge3\), and different conjugate pairs
belong to different factors. The centralizer of the derived
product therefore consists of scalars on the \(2e\) oriented
blocks. By Faltings, the geometric endomorphism algebra is
\(E\), since it contains \(E\) and has dimension at most \(2e\).
The field algebra gives absolute simplicity. The independent
projective factors distinguish the \(e\) projective classes.

Let \(L\) be a normal closure of \(F\). The chosen degree-one
place of \(K\) over \(p\) has Frobenius in
\(\operatorname{Gal}(LK/K)\) acting as an \(e\)-cycle on the
embeddings of \(F\). The set of elements with this cycle type
is thus nonempty and conjugation invariant. Chebotarev gives
positive density of this type, and the degree-one condition
discards only a density-zero set. Its density is \(1/2\) when
\(e=2\). At these degree-one places the rational prime splits
in \(I\), which is contained in \(K\), and is inert in \(F\).
Thus \(\theta\) is transitive on the projective classes.
The action is defined over \(K\) and the projective factors
are independent, so the \(M=1\) criterion applies. For \(m=1\)
the orbit costs \(e=D_1\). For the rank-four balanced case
it costs \(2e=D_2\). In both cases the strict inequality
forces ordinary reduction outside a density-zero set.

For prime \(e\), one may instead choose \(F\) before \(p\): take the
degree-\(e\) subfield of \(\mathbb Q(\zeta_r)\) for a rational prime
\(r\equiv1\bmod 2e\) unramified in \(k\) and prime to \(2N\). This
subfield is totally real and unramified at \(2N\); its ramification
at \(r\) prevents it from being contained in \(k\). The Galois
extension \(F/\mathbb Q\) has prime degree, so \(F\) and \(k\) are
linearly disjoint. Chebotarev supplies \(p\nmid 2Nr\) completely
split in \(k\) and inert in \(F\). Proceed with the same moduli
construction. The embedding \(K\hookrightarrow\mathbb Q_p\) again
excludes \(F\subseteq K\), so \(FK/K\) is cyclic of degree \(e\).
Each nonidentity Frobenius acts as an \(e\)-cycle, giving the stated
density.

For a concrete initial PEL point in the degree-two construction, take
\(N=3\), \(F=\mathbb Q(\sqrt{17})\), and \(p=61\). The rational
prime splits in \(I\), and the sum of quadratic characters of
\(x^3-x\) over \(\mathbb F_{61}\) is 10. Hence the elliptic curve has
trace \(-10\) and determinant 61. On its 3-torsion the Frobenius at
each prime of \(I\) over 61 lies in the multiplicative group of
\(\mathcal O_I/3\mathcal O_I\simeq\mathbb F_9\), since it commutes
with the CM action. Its polynomial modulo 3 is \((X-1)^2\), so this
field element is 1. Thus 61 is completely split in \(k\). It is inert
in \(F\) since \(17^{30}\equiv-1\bmod61\).
\(\square\)

Wedhorn's theorem gives an open dense ordinary locus in a good-reduction
PEL fiber when the local reflex field is \(\mathbb Q_p\). Here the reflex
field is contained in \(I\), which splits at the indicated rational
primes. The corollary controls the reductions of the same fixed
characteristic-zero fiber as those primes vary. Noot's Proposition 1.3
is the precedent for preserving monodromy under specialization. His
ordinary-reduction Theorem 2.2 assumes the odd tensor product of
rank-two standard factors. Here the semisimple factors are of type
\(A_{h-1}\), with standard and dual summands.
[Wedhorn, Theorem 1.6.3](https://www.numdam.org/item/ASENS_1999_4_32_5_575_0.pdf).
[Noot, Proposition 1.3 and Theorem 2.2](https://www.numdam.org/item/CM_1995__97_1-2_161_0/).

The use of a sufficiently divisible \(M\) in the general theorem
can weaken the rounded penalties. The
[rank-two semilinear bound](CM_NEWTON_SEMILINEAR_PROJECTIVE_BOUNDS.md)
uses the order of the action on the coefficient field and local
gcd constraints, which can give sharper costs. The general proof
gives all exterior-coefficient constraints without a full
endomorphism-algebra hypothesis. Cantoral Farfán, Li, Mantovan,
Pries and Tang obtain ordinary and \(\mu\)-ordinary density theorems
for simple signatures using a more specialized normalized trace
and their endomorphism-algebra assumptions. The criteria here use
the actual orbit type and do not supply a fixed variety realizing
every numerical type. [Cantoral Farfán–Li–Mantovan–Pries–Tang,
Proposition 4.8 and Theorems 4.18–4.19](https://arxiv.org/html/2508.11174v1).
