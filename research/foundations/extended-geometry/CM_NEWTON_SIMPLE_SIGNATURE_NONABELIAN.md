# Density-one \(\mu\)-ordinary reduction without an abelian coefficient field

**Theorem.** Let \(A/K\) be a polarized abelian variety over a number
field and suppose its full geometric rational endomorphism algebra is
a CM field \(E\). Put \([E:\mathbb Q]=2e\) and
\(\dim A=eh\), with \(h\ge3\). Assume exactly one pair of conjugate
embeddings has Hodge multiplicities \(1,h-1\), and all the other
pairs have multiplicities \(0,h\). Then \(A\) has \(\mu\)-ordinary
reduction at a set of places of natural density one over \(K\).
If \(L/K\) is the minimal field defining all the geometric
endomorphisms, the natural density of ordinary reduction over
\(K\) is \(1/[L:K]\). Writing \(\tau:E\hookrightarrow\overline K\)
for the unique oriented embedding with multiplicity one, one has
\(L=K\tau(E)\).

The reference polygon at a coefficient place is the average of the
sorted Hodge vectors on its Frobenius orbit. The theorem permits
non-Galois \(E/\mathbb Q\), or Galois \(E/\mathbb Q\) with nonabelian
Galois group, and permits the endomorphisms to be defined only over
an extension of \(K\). Theorem 4.19 of Cantoral
Farfán–Li–Mantovan–Pries–Tang gives the same density-one conclusion
under Assumptions (F), (S) and (A), where (A) requires \(E/\mathbb Q\)
abelian. Their Theorem 4.18 already gives potentially density-one
ordinary reduction under (F) and (S), and determines the ordinary
density over \(K\) when \(E/\mathbb Q\) is Galois. At the original
primes here the reference polygon can have fractional slopes.
[Cantoral Farfán–Li–Mantovan–Pries–Tang, Assumptions 2.1, 2.5 and 2.6 and Theorems 4.18–4.19, arXiv v1](https://arxiv.org/html/2508.11174v1).

The key is to bound every relative exterior coefficient using the
order of the action of Frobenius on \(E\). This power is small enough
that the unique nondefinite orbit exhausts the available rounding
allowance if its reduction fails to attain its reference polygon.
The argument for the Frobenius power only needs one Hodge
multiplicity one, so it also yields a criterion with several mixed pairs.

## Full relative projective monodromy on each block

Suppose \(E\) is the full geometric rational endomorphism algebra
and some pair has smaller multiplicity one. Fix an auxiliary
prime \(\ell\), work over an algebraic closure of \(\mathbb Q_\ell\),
and let \(G\) be the algebraic monodromy group on \(V_\ell(A)\).
Set \(D=(G^\circ)^{\mathrm{der}}\). After a finite extension defining
the endomorphisms and killing the component group, Faltings's
semisimplicity and endomorphism theorems imply that \(G^\circ\) is
reductive with centralizer \(E\otimes_{\mathbb Q}\overline{\mathbb Q}_\ell\).
Thus each \(E\)-embedding block \(V_\sigma\) is irreducible under
\(G^\circ\) and under \(D\). For the latter assertion, the connected
center acts by scalars, and together with \(D\) generates \(G^\circ\).
[Faltings, Theorems 3 and 4](https://math.uchicago.edu/~drinfeld/Deligne%27s_conjecture_Manin_conf/Faltings_argument/Faltings.pdf).

For any chosen \(\sigma\), choose a local embedding of a finite
extension of \(K\) containing a normal closure of \(E\), to send the
Hodge label of multiplicity one to \(\sigma\). This is possible by
extending an embedding of that normal closure. The auxiliary
coefficient block and global connected group stay fixed when the
local embedding varies. Hodge–Tate comparison on this block gives
weights \(0,1\) with the second multiplicity one, up to the common
sign convention on weights. The Sen operator on this block therefore
gives a rank-one projector in the projected Lie algebra of
\(G^\circ\) after extending scalars to \(\mathbb C_\ell\) and
changing its sign if needed. The Lie algebra of the local monodromy
group lies in that of the global group. Serre's Theorems 2 and
\(2''\) explicitly give the Lie membership of the Hodge–Tate grading
operator used here.
[Faltings, *\(p\)-adic Hodge theory*](https://public.websites.umich.edu/~bhattb/almost_purity_2011/Faltings_-_p-adic-Hodge-theory.pdf),
[Sen, Theorem 1](https://annals.math.princeton.edu/1973/97-1/p06),
[Serre, Astérisque 65, §1.4, Theorems 2 and \(2''\)](https://www.numdam.org/article/AST_1979__65__155_0.pdf).

An irreducible Lie subalgebra \(\mathfrak l\subset\mathfrak{gl}(V)\)
containing a rank-one projector \(P\) contains \(\mathfrak{sl}(V)\).
Write \(V=L\oplus W\) for the image and kernel of \(P\). The three
eigenspaces of \(\operatorname{ad}P\) on \(\mathfrak l\) are the
block-diagonal part, a subspace of \(\operatorname{Hom}(L,W)\), and
a subspace of \(\operatorname{Hom}(W,L)\). If \(W_-\) is the span
of images of \(L\) under the second space, the bracket shows that
\(W_-\) is stable under the block-diagonal part. The three eigenspaces
all preserve \(L\oplus W_-\), so irreducibility gives \(W_-=W\).
Since \(L\) is one-dimensional the second space is all of
\(\operatorname{Hom}(L,W)\). Applying this argument to the dual
representation gives all of \(\operatorname{Hom}(W,L)\) as well.
Their brackets give the remaining traceless matrices. Applying the
lemma to each \(V_\sigma\) gives the full \(\mathrm{PGL}_h\) as the
projective image of \(D\) on that block. This conclusion allows
different blocks to see the same simple factor of \(D\).
For simple signature, the full product of projective groups over
the conjugate pairs also follows from the Mumford–Tate group
calculation and Mumford–Tate conjecture in the cited work.
[Cantoral Farfán–Li–Mantovan–Pries–Tang, Lemma 2.11 and Remark 2.13](https://arxiv.org/html/2508.11174v1#S2.SS4).

## The Frobenius power equals the field-action order

At an absolute degree-one good place \(v\) of characteristic \(p\),
away from a fixed finite set, let \(\gamma\) be the action of
arithmetic Frobenius on \(E\) and put \(t=\operatorname{ord}(\gamma)\).
Absolute degree-one places have natural density one.
For each coefficient place \(\lambda\mid p\), let
\(f_\lambda=[E_\lambda:\mathbb Q_p]\), let \(S_{\lambda,j}\) be the
sum of the \(j\) least of its \(h\) relative slopes, and let
\(B_{\lambda,j}\) be the corresponding sum for the reference polygon.
Put \(\Delta_{\lambda,j}=S_{\lambda,j}-B_{\lambda,j}\ge0\).
As usual, each relative slope occurs \(f_\lambda\) times on the
rational crystalline summand, and \(\Delta_{\lambda,0}\) and
\(\Delta_{\lambda,h}\) are zero.
[Orbit-weighted comparison](CM_NEWTON_DISCONNECTED_SPECTRUM_INEQUALITIES.md#averaged-hodge-polygons-and-compression).

**Proposition.** Under the full geometric endomorphism hypothesis
and with at least one Hodge multiplicity one in a mixed pair, the
following hold on a set of places of natural density one over \(K\):
\[
  \sum_{\lambda\mid p} f_\lambda\left\lceil t S_{\lambda,j}\right\rceil
     \le etj-1,\qquad 1\le j<h.
\]
If the smaller multiplicity of an original conjugate pair \(C\)
is \(m_C\), then
\[
  \sum_{\lambda\mid p} f_\lambda\left\lceil t\Delta_{\lambda,j}\right\rceil
     \le t\sum_C\min(j,h-j,m_C)-1.                 \tag{1}
\]

The normalizer of \(E\) contains \(G\), and the induced action on
\(E\) is constant on each component. Let \(g\) belong to a component
whose action on \(E\) has order \(t\). Then \(g^t\) preserves every
block \(V_\sigma\). Each surjection from \(D\) to an individual
\(\mathrm{PGL}(V_\sigma)\) sees exactly one simple factor of
\(D^{\mathrm{ad}}\). Every factor is seen, since a connected
semisimple factor with scalar image on every block acts trivially
on the faithful \(V_\ell(A)\).
Conjugation by \(g^t\) preserves each such surjection's kernel. On
its target it acts by the projective class of \(g^t|_{V_\sigma}\),
so on each adjoint factor it acts by an inner automorphism.
Thus its action on \(D\) is inner as well. An automorphism of the
connected semisimple group \(D\) inducing the identity on the
adjoint quotient is trivial, since its pointwise difference from
the identity takes values in the finite center.

Adjust \(g\) by an element of \(D\) to get \(g_0\) whose conjugation
preserves a pinning of \(D\). The \(t\)-th power of this pinned
automorphism is inner, and hence trivial. By Schur's lemma \(g_0^t\)
is scalar on every coefficient block. Also \(g_0\) commutes with a
nonzero cocharacter \(\nu\) obtained by summing the simple coroots
in an orbit of its pinned action.
[SGA 3, Exposé XXIV, §1](https://jcreinhold.github.io/sga/iii/24-automorphisms-of-reductive-groups.html).

Use the normalized coefficient function from the
[general orbit-weighted argument](CM_NEWTON_DISCONNECTED_SPECTRUM_INEQUALITIES.md)
on this component, taking the \(t\)-th power:
\[
 H_j(x)=
 \prod_{\sigma:E\hookrightarrow\overline{\mathbb Q}_\ell}
       \frac{a_{j,\sigma}(x^t)^h}{\delta_\sigma(x^t)^j}.
\]
Along \(x=g_0\nu(z)\) the block scalars of \(g_0^t\) cancel and the
remaining expression is that of \(\nu(z^t)\). At least one block
has unequal cocharacter weights, and for each \(1\le j<h\) the
elementary coefficient there has at least two monomials of different
weights. Evaluating the resulting Laurent polynomials at generic
complex unit \(z\) gives a strictly smaller absolute value than at
\(z=1\), by the triangle inequality in that block and weak bounds
in the others. Thus \(H_j\) is nonconstant on each component. The
coefficient norm argument applies with \(t\): the value at Frobenius
is \(N_{E/\mathbb Q}(a_j)^h/p^{ehtj}\), of absolute value at most
\(\binom hj^{2eh}\), where \(a_j\in E\) is the relative coefficient
of \(\pi_v^t\). If its numerator were divisible by its denominator,
the value would lie in a fixed finite set of integer levels. These
levels form proper closed sets on every component. The functions with
the relevant \(t\) on different components together give a regular
conjugation-invariant function, since conjugation preserves the
order of the field action. Algebraic Chebotarev therefore excludes
these levels on a density-zero set. The strict valuation bound and
\(v_\lambda(a_j)\ge \lceil tS_{\lambda,j}\rceil\) give the first
inequality. In this argument \(\pi_v^t\) centralizes \(E\), so the
relative roots have valuations \(t s_{\lambda,k}\) and the
valuation on unramified \(E_\lambda\) is integral.
[Rajan, Theorem 3](https://repository.ias.ac.in/38317/1/2-PUB.pdf).

Let \(\theta\) be a Frobenius in a normal closure of \(E\) at the
chosen reduction place. It acts on the embedding labels by
postcomposition and commutes with precomposition by \(\gamma\).
The Hodge multiplicities \(m_\sigma\), counting weights equal to
one, are invariant under
\(\sigma\mapsto\theta\circ\sigma\circ\gamma^{-1}\).
Indeed, apply the residue field Frobenius to a coefficient eigenvector
in the Hodge space of the geometric special fiber, using its
\(\gamma\)-conjugation of the specialized endomorphisms. The
characteristic-zero Hodge multiplicities specialize unchanged away
from the finite set where the action or its idempotents fail to
extend over the unramified local lattices. It follows that
\(\theta^t\) preserves \(m_\sigma\). On a \(\theta\)-orbit of length
\(f_\lambda\) the sorted Hodge vectors therefore have period dividing
\(\gcd(t,f_\lambda)\). In particular \(tB_{\lambda,j}\) is integral.
Subtract the reference sums from the first inequality. The
contribution of an original pair with smaller multiplicity \(m_C\)
to the reference sum over all embedding labels is
\(j-\min(j,h-j,m_C)\). This gives (1).

## The unique mixed pair

Let \(\tau\) denote the unique label with \(m_\tau=1\). Uniqueness
uses \(h\ge3\). Invariance of the Hodge multiplicities gives
\(\theta\circ\tau=\tau\circ\gamma\). By commutation,
\(\theta^i\circ\tau=\tau\circ\gamma^i\) for all \(i\), so the
\(\theta\)-orbit of \(\tau\) has length exactly \(t\). The same
holds for \(\bar\tau\), and its orbit can coincide with that of
\(\tau\). At a place containing one of these labels the reference
polygon is the average of definite Hodge vectors (each a constant
vector), one vector with a jump between slopes \(h-1\) and \(h\)
if \(\tau\) is present, and one with a jump between slopes \(1\)
and \(2\) if \(\bar\tau\) is present. All such jumps in the average
are \(1/t\). A place containing neither label has a straight-line
reference polygon.

Suppose the actual polygon exceeds the reference polygon. It must
exceed it at some strict kink of the reference polygon. If the
polygons agree at each kink and at the endpoints, convexity of
the actual polygon bounds it above by the reference segment between
each adjacent pair of such points. The usual lower bound by the
reference polygon then gives equality. Every strict kink here
lies in an orbit of length \(t\). A positive excess at one such
kink \(j\) gives a term of at least \(t\) on the left of (1).
The right side is \(t-1\), since there is exactly one pair with
smaller multiplicity one. This contradiction proves the theorem.
The absence of a unique oriented multiplicity-one label in rank
two prevents this argument for the orbit length in that rank.
[Rank-two semilinear orbit bounds](CM_NEWTON_SEMILINEAR_PROJECTIVE_BOUNDS.md).

For the ordinary density, if \(\gamma=1\), all oriented
multiplicities are \(\theta\)-invariant, so the reference slopes
are all zero or one. If \(\gamma\ne1\), at least one mixed orbit
has a reference slope jump \(1/t<1\), which forces a slope
strictly between zero and one in that orbit. Thus outside a set
of density zero the ordinary places are exactly those where
\(\gamma=1\). The kernel of the action on \(E\) defines \(L\),
so finite Chebotarev gives density \(1/[L:K]\).

The degree can also be read from the distinguished embedding. For
\(g\in G_K\), let \(\gamma(g)\) act on the endomorphism field. The
semilinear action of \(g\) on the characteristic-zero Hodge space
sends the eigenspace labelled \(\sigma\) to the one labelled
\(g\circ\sigma\circ\gamma(g)^{-1}\). Uniqueness of \(\tau\) gives
\(g\circ\tau=\tau\circ\gamma(g)\). Thus \(G_K\) preserves the subfield
\(\tau(E)\), and \(g\) acts trivially on the endomorphisms exactly
when it fixes that subfield pointwise. This identifies \(L\) with
\(K\tau(E)\), which is Galois over \(K\) even for non-Galois
\(E/\mathbb Q\).

### Examples with ordinary density \(1/2\)

**Corollary.** For every CM field \(E\) and integer \(h\ge3\),
there exists an absolutely simple polarized abelian variety over
a number field \(K\subset\mathbb R\) with full geometric rational
endomorphism algebra \(E\) and simple signature in relative rank
\(h\). Its minimal endomorphism field is quadratic over \(K\).
It has ordinary reduction at density \(1/2\) and
\(\mu\)-ordinary reduction at density one. On a set of places of
density \(1/2\), the mixed coefficient place has relative slopes
\(0,\tfrac12,\ldots,\tfrac12,1\).

Put \(F=E^+\), fix an embedding \(E\subset\mathbb C\), and write
\(c\) for complex conjugation on \(E\). By weak approximation at
the real places of \(F\), choose a diagonal hermitian form \(H\)
on \(E^h\) with signature \((1,h-1)\) at one real place and
definite at every other real place. On \(\Lambda=\mathcal O_E^h\)
take the alternating form
\(\psi(x,y)=\operatorname{Tr}_{E/\mathbb Q}(\delta H(x,y))\)
with \(\bar\delta=-\delta\) and \(\delta\ne0\), scaled by a
positive integer to be integral. Coordinatewise conjugation
\(c_0\) on \(\Lambda\) is \(c\)-semilinear and anti-symplectic.
On each complex coordinate of \(\Lambda\otimes\mathbb R\), choose
the complex structure \(J\) to act by \(i\) or \(-i\) according
to the sign that makes \(\psi(x,Jx)\) positive. This produces a
polarized complex torus with the prescribed simple signature,
an \(\mathcal O_E\)-action, and \(c_0J=-Jc_0\).

Choose full level \(N\ge3\) prime to the degree of the
polarization and take \(k\subset\mathbb C\) to be the compositum
of the normal closure of \(E\) and \(\mathbb Q(\zeta_N)\).
It contains the reflex field because the signature polynomial
\[
       P_a(T)=\prod_\sigma (T-\sigma(a))^{m_\sigma}
       \qquad (a\in\mathcal O_E)
\]
has its coefficients in the normal closure. Use the fine
characteristic-zero PEL moduli space over \(k\) of the
prescribed polarization type and signature. Its level map
\(\eta:\Lambda/N\Lambda\longrightarrow A[N]\) is
\(\mathcal O_E\)-linear and identifies the alternating trace
pairing with the Weil pairing via \(\zeta_N\).
[Kottwitz, §5, pp. 389–392](https://webusers.imj-prg.fr/~christophe.cornut/M2/Ref/Kottwitz_PointsShiModp.pdf).

The nontrivial automorphism of \(k/(k\cap\mathbb R)\) sends the
determinant condition \(P_a\) to \(P_{c(a)}\) and sends
\(\zeta_N\) to \(\zeta_N^{-1}\). On that conjugate moduli problem,
replace \(\iota\) by \(\iota\circ c\) and \(\eta\) by
\(\eta\circ c_0\), keeping the polarized abelian scheme.
Semilinearity of \(c_0\) gives linearity for the new
endomorphism action, and its anti-symplectic property restores
the original pairing condition. Both relabelings are involutions.
They give effective descent of the fine space and its underlying
polarized abelian scheme to \(k\cap\mathbb R\).
On the diagonal torus, \(c_0\) identifies \((J,\psi)\) with
\((-J,-\psi)\) and gives the same change of labels and level.
Thus this torus, marked by the lattice level map, is a real
point of the descended space. The geometric connected
component through it is fixed by complex conjugation and
descends over some finite extension \(k_0/(k\cap\mathbb R)\)
inside \(\mathbb R\). Indeed, its stabilizer on the finite set
of geometric components is open and contains complex
conjugation, so its fixed field lies in \(\mathbb R\).
[Blasius–Guerberoff, §4.1 and Theorem 5.2.2](https://msp.org/ant/2017/11-10/ant-v11-n10-p02-p.pdf).

The arithmetic subgroup in the complex uniformization of this
component intersects \(\operatorname{Res}_{F/\mathbb Q}\mathrm{SU}(H)\)
in an arithmetic subgroup \(\Gamma\) of that derived group.
The real group has one noncompact factor \(\mathrm{SU}(1,h-1)\)
and compact factors at the other places. The projection of
\(\Gamma\) onto the noncompact factor is a lattice. Discreteness
follows from discreteness of \(\Gamma\) and compactness of the
other factors, and finite covolume follows by pushing forward
the invariant finite measure. By Borel density the projected
group is Zariski dense in the noncompact factor. The rational
Zariski closure of \(\Gamma\) contains that factor. To see this,
the Lie algebra of the real Zariski closure surjects onto the
noncompact simple Lie algebra. The kernel of its projection
onto the compact factors must project onto the noncompact
factor, since otherwise a compact Lie algebra would have a
noncompact simple quotient. The resulting full factor over
\(\mathbb C\) has every other factor as a Galois conjugate.
Thus Betti monodromy, and hence geometric \(\ell\)-adic
algebraic monodromy by comparison, contains the full derived
group.
[Borel–Harish-Chandra, Theorem 7.8](https://www.mathi.uni-heidelberg.de/~wienhard/retreat10/references/borel_harishchandra_2.pdf),
[Borel, density theorem](https://www.jstor.org/stable/1970150),
[Noot, §1.8](https://www.numdam.org/article/CM_1995__97_1-2_161_0.pdf).

Specialize this family over \(k_0\) while preserving its
arithmetic \(\ell\)-adic image and the chosen real place. If
\(\mathcal G\) is the compact \(\ell\)-adic image, choose a
finitely generated open normal pro-\(\ell\) subgroup \(U\).
Surjection of a closed subgroup onto the finite quotient
\(\mathcal G/\Phi(U)\) implies surjection onto
\(\mathcal G\), by the defining property of the Frattini
subgroup on its intersection with \(U\). On an étale chart
over affine space, apply Hilbert irreducibility to the cover
for \(\mathcal G/\Phi(U)\) composed with the chart. After
shrinking the chart and the target so that the composed cover
is finite étale, choose a nearby real point in the chart and
use the real inverse function theorem to give a local section.
Hilbert irreducibility with weak approximation provides a
\(k_0\)-rational parameter in this real neighborhood with
connected fiber in the composed cover. The chart fiber is a
single closed point whose residue field \(K\) embeds in
\(\mathbb R\) by the local section. Connectedness of the cover
fiber gives the required specialization of the finite
quotient, hence of the full \(\ell\)-adic image.
[Noot, Proposition 1.3](https://www.numdam.org/article/CM_1995__97_1-2_161_0.pdf),
[Serre, *Topics in Galois Theory*, Theorem 3.5.3](https://indico.math.cnrs.fr/event/11410/attachments/4760/7351/SerreTopicsGaloisTheory.pdf).

Over an algebraic closure of \(\mathbb Q_\ell\), each derived
factor acts by the standard and dual representations on its
conjugate pair of \(E\)-blocks. These are inequivalent for
\(h\ge3\), and distinct factors act on distinct pairs. The
centralizer of the derived group therefore has dimension
\(2[F:\mathbb Q]\). Faltings's endomorphism theorem bounds the
dimension of the geometric endomorphism algebra by that
number, and the specialized \(E\)-action supplies equality.
Thus the algebra is \(E\) and the fiber is absolutely simple.
The action is defined over \(Kk\). This is a quadratic
extension of the real field \(K\) and its nontrivial
automorphism acts on \(E\) by \(c\), by the descent of the
labeled action. Consequently \(L=Kk\). At primes with this
Frobenius class, the unique mixed orbit consists of
\(\tau,\bar\tau\), so its reference slopes are the averages
of the sorted vectors with multiplicities \(1\) and \(h-1\).
They give the displayed slopes and the density assertions by
the theorem.

## Several mixed pairs on the same orbit

The exact power also permits
[connected projective-spectrum compression](CM_NEWTON_PROJECTIVE_SPECTRUM_FIELDS.md).
For \(x=g_0y\) with \(y\in G^\circ\), the power \(x^t\) is
\(g_0^t\) times an element of \(G^\circ\), by normality. Since
\(g_0^t\) acts by block scalars, \(x^t\) has the same projective
block spectra as that element. The projective compression then
identifies the Newton vectors, up to constant shifts, on labels
with the same restriction to the projective spectrum field \(E_0\).
The averaged Hodge vectors differ by constant shifts as well.
Equality of the Hodge and Newton endpoints identifies those shifts,
so the excesses agree on each oriented class.

In particular, suppose \(h\ge3\) and every mixed pair has smaller
multiplicity one. Then \(E_0\) is CM when a mixed class exists,
since the unequal multiplicities on a mixed conjugate pair give
distinct projective spectrum classes. Let \(B\) count the mixed
conjugate classes, and view a single \(\theta\)-orbit as an orbit
of oriented embeddings of \(E_0\). Outside the density-zero set
of the proposition, if each orbit containing a mixed label has
length at least \(tB\), the reduction is \(\mu\)-ordinary. Indeed,
divide (1) by the common number of labels over each oriented class.
The result is
\[
  \sum_{C\ {\rm oriented}}\lceil t\Delta_{C,j}\rceil\le tB-1.
\]
The excess is constant along each oriented \(\theta\)-orbit. Any
failure must be in an orbit containing a mixed label, since only
these have reference kinks. At a failed kink the orbit would
contribute at least its length to the left side. One sufficient
configuration is that all \(B\) oriented multiplicity-one labels
belong to a single orbit of length \(tB\), since its conjugate
contains all the multiplicity-\((h-1)\) labels and also has length
\(tB\). The orbit may also contain definite labels, and its
conjugate may be the same orbit.

## A local bound at an individual kink

For a coefficient place of degree \(f\), put \(q=\gcd(f,t)\).
At each integer abscissa where the actual relative polygon has a
vertex, its ordinate \(S_j\) satisfies \(f S_j,t S_j\in\mathbb Z\).
For the first integrality use the geometric isocrystal on a single
coefficient-embedding block with its \(f\)-fold semilinear Frobenius,
which has slopes scaled by \(f\). A slope break cuts out a slope
subisocrystal, whose determinant has integral valuation in the Witt
fraction field.
The second follows from the relative characteristic polynomial of
\(\pi_v^t\) over unramified \(E_\lambda\). Hence \(qS_j\) is integral
at each vertex. The period argument above also gives \(qB_j\in\mathbb Z\)
for every reference prefix. Write \(d=q(B_{j+1}-2B_j+B_{j-1})\)
at a strict reference kink \(j\), so \(d\) is a positive integer.
If the excess there is positive, then
\[
    \left\lceil t\Delta_j\right\rceil
    \ \ge\
    \left\lceil \frac{t}{q}\min(1,d/2)\right\rceil .
\]
An actual vertex at \(j\) entails \(q\Delta_j\ge1\). Otherwise
let \(a<j<b\) be the neighboring actual vertices. The actual
segment lies above the chord of the reference polygon on
\([a,b]\). The contribution of its slope jump \(d/q\) at \(j\)
alone to the chord gap at \(j\) is
\((d/q)(j-a)(b-j)/(b-a)\ge d/(2q)\), since the two distances are
positive integers. This proves the bound. It is sharp among the
local convex polygons with these vertex lattice constraints.
For \(d=1\) average the two reference slopes adjoining the kink,
removing that vertex. For \(d\ge2\) move \(1/q\) from the higher
slope to the lower one. Each change preserves order and the
vertex lattice and raises only the prefix at \(j\). At an actual
reduction the same changes may have to be made at other
coefficient places or other indices by projective equivalence and
polarization.
