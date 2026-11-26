# Pfaffian bounds for linked balanced rank-four coefficient blocks

Let \(A/K\) be a polarized abelian variety over a number field
with a Galois-stable geometric CM coefficient field \(E\).
Assume Rosati acts on \(E\) by complex conjugation \(c\), the
relative rank is four, and every conjugate Hodge pair is
definite or of signature \((2,2)\). Write \(F=E^+\).
Fix an auxiliary prime \(\ell\). Require the embedding
blocks of \(E\) in \(V_\ell(A)\) to be irreducible over
\(\overline{\mathbb Q}_\ell\) under connected monodromy
\(G^\circ\). This holds by Faltings's theorems if
\(E=\operatorname{End}^0(A_{\overline K})\).

Let \(D=(G^\circ)^{\mathrm{der}}\), let \(\Gamma\) be the
Galois image on \(E\), and let \(E_0\subset E\) be the
[connected projective spectrum field](CM_NEWTON_PROJECTIVE_SPECTRUM_FIELDS.md).
Conjugation by \(G\) preserves the projective equivalence
of embedding blocks and so \(\Gamma\) preserves \(E_0\).
Put \(F_0=E_0\cap F\). The conjugate projective classes
are the restriction fibers to \(F_0\). Here \(E_0\) is
either CM or totally real.

**Theorem.** Suppose the balanced pairs form exactly one
conjugate projective class. Suppose the kernel of
\[
       \langle\Gamma,c\rangle\longrightarrow\operatorname{Aut}(E_0)
\]
is trivial when \(E_0\) is CM, and is \(\{1,c\}\) when
\(E_0\) is totally real.

For an algebraic monodromy component, consider its
permutation of the simple adjoint factors of \(D\).
Suppose the two factors of every tensor-standard block
lie in distinct orbits of this permutation. Here a
tensor-standard block has derived projective image of
type \(A_1\times A_1\) acting as \(2\otimes2\).
Then the primes in that component, or in its conjugates,
have \(\mu\)-ordinary reduction outside a set of natural
density zero over \(K\). In particular, this conclusion
holds at density one over \(K\) if all blocks have
simple derived projective images. The restriction on
field actions is automatic when \(\Gamma\subseteq\{1,c\}\).

The individual-block hypothesis permits type \(C_2\)
and differs from requiring the full \(\mathrm{PGL}_4\).
The preceding [Pfaffian theorem](CM_NEWTON_RANK_FOUR_PFAFFIAN.md)
has one actual balanced pair and uses the full
\(\mathrm{PGL}_4\) on every embedding block. Fité's
Theorem 4 concerns positive-density ordinary reduction
for the imaginary-quadratic fourfold case. The
simple-signature hypothesis of Cantoral
Farfán–Li–Mantovan–Pries–Tang excludes \((2,2)\).
[Fité, Theorem 4](https://arxiv.org/pdf/2005.10185v3),
[Cantoral Farfán et al., Assumption 2.5 and Theorems 4.18–4.19](https://arxiv.org/html/2508.11174v1).

## Exact power on the indicated components

Write \(\gamma\in\Gamma\) for the field action of this
component and \(t=\operatorname{ord}(\gamma)\). A block
is also irreducible under \(D\), since the connected
center acts by scalars. The classification of faithful
irreducible four-dimensional semisimple representations
gives four possibilities for the derived image:
\(\mathrm{SL}_4\) standard or dual,
\(\mathrm{Sp}_4\) standard,
\(\mathrm{SL}_2\) on \(\operatorname{Sym}^3\), or
the tensor product of two \(\mathrm{SL}_2\) standards.
In this list we use the specified *images*, so any
simple factors acting trivially on a block are omitted.
[Hui, §3.1](https://arxiv.org/pdf/2208.04004).

Since \(g^t\) preserves every embedding block, its
conjugation on \(D\) preserves each block's projective
kernel and its set of supporting simple factors.
The singleton supports therefore stay fixed. For a
two-element support, the orbit hypothesis prevents
\(g^t\) from exchanging the two elements. Every simple
factor belongs to at least one support because the
representation of \(D\) is faithful. On a factor of
type \(A_3\), \(g^t\) induces inner conjugation through
the full projective image \(\mathrm{PGL}_4\).
Types \(A_1\) and \(C_2\) have no outer automorphisms.
Thus conjugation by \(g^t\) on \(D\) is inner.

Choose \(g_0\) by multiplying \(g\) by a member of \(D\)
so that its conjugation preserves a pinning. The
\(t\)-th power of the pinned automorphism is inner and
hence trivial. Consequently \(g_0^t\) is scalar on
every block. For any \(x\in G^\circ\),
\((g_0x)^t=g_0^t x'\) for some \(x'\in G^\circ\),
by normality. In particular, although a \(t\)-th power
can lie outside \(G^\circ\), its projective block
spectra satisfy all the equivalences on \(G^\circ\).
[SGA 3, Exposé XXIV, §1](https://jcreinhold.github.io/sga/iii/24-automorphisms-of-reductive-groups.html).

Take a fixed cocharacter \(\nu\) of the pinned action
by using the sum of positive coroots on every \(A_3\)
and \(C_2\) factor. For \(A_1\), use positive integer
multiples of the coroot, chosen constant within each
orbit and distinct between orbits. The weights of
the \(A_3\) and \(C_2\) blocks are \(3,1,-1,-3\).
For \(\operatorname{Sym}^3\) the weights are nonzero
multiples of those four integers. The tensor weights
have the form \(\pm a\pm b\) for two unequal positive
integers \(a,b\). In particular there are no zero
weights in any block. On \(g_0\nu(z)\) the \(t\)-th
power is \(g_0^t\nu(z^t)\). This is the cocharacter
needed both for the strict normalized coefficient norm
and for nonconstancy of the product of shifted
determinants in the Pfaffian test. In the latter
product every Laurent factor is of the form
\(1+a_\sigma z^{tw}\) for \(a_\sigma\ne0\) and
\(w\ne0\), with both positive and negative \(w\).

## Synchronization and the two field-action cases

At a good absolute degree-one place, write \(\theta\)
for Frobenius on the coefficient embedding labels and
\(\Delta_{\lambda,j}\) for the difference between the
actual and averaged Hodge prefixes at a coefficient
place. The exact power gives the norm inequality and
the projective comparison of slopes with this \(t\).
The averaged Hodge vectors for embeddings in one
oriented projective class differ by constants, since
the individual Hodge vectors differ by constants
and simultaneous \(\theta\)-translates remain in a
common class. The analogous constant for the actual
vectors is the same, by equality of their endpoints
with the respective Hodge endpoints. Thus the
excesses agree throughout each oriented class. They
also agree at \(j=2\) on conjugate classes.

The permutation \(\theta\gamma^{-1}\) preserves the
original \(0,2,4\) multiplicities. In particular
\(\theta^t\) preserves them, so \(t\) times an
averaged Hodge prefix is integral. Subtract those
prefixes from the strict norm bound at \(j=2\)
and divide by \(d=[F:F_0]\).
The single conjugate balanced class contributes
\(2td\) to the midpoint allowance before division.
For \(E_0\) CM there are \(d\) embeddings over each
oriented reduced label. For \(E_0\) real there are
\(2d\) such embeddings. Hence the respective bounds
outside a density-zero set are
\[
 \sum_{C\text{ oriented}}\lceil t\Delta_{C,2}\rceil\le 2t-1
 \quad (E_0\text{ CM}),\qquad
 2\sum_C\lceil t\Delta_{C,2}\rceil\le 2t-1
 \quad (E_0\text{ real}).
\]
Here the excess is constant on each reduced
\(\theta\)-orbit. Each mixed reference has its
single kink at the midpoint, with jump equal to
the fraction of balanced labels in its orbit.
This fraction can be computed on the reduced
orbit, because each actual orbit visits every
label on its projected orbit equally often.
[Exact-power compression](CM_NEWTON_SIMPLE_SIGNATURE_NONABELIAN.md#several-mixed-pairs-on-the-same-orbit),
[Orbit-weighted comparison](CM_NEWTON_DISCONNECTED_SPECTRUM_INEQUALITIES.md#averaged-hodge-polygons-and-compression).

**Case \(c\notin\langle\gamma\rangle\).** If \(E_0\) is
CM, restriction preserves the order of \(\gamma\) and
the absence of \(c|_{E_0}\) from its cyclic subgroup.
The two balanced oriented labels are either in
separate conjugate reduced orbits of length \(t\),
or in one reduced orbit of length \(2t\), by the
[one-pair orbit calculation](CM_NEWTON_RANK_FOUR_PFAFFIAN.md#midpoint-bounds-on-the-two-types-of-field-action)
on \(E_0\). Any failure contributes at least \(2t\).
If \(E_0\) is real, restriction again preserves the
order of \(\gamma\). The unique balanced reduced label
is fixed by \(\theta\gamma^{-1}\), so its \(\theta\)-orbit
has length \(t\). The factor two gives the same
contribution \(2t\) for a failure. Both contradict
the compressed inequality.

**Case \(c\in\langle\gamma\rangle\).** The half power
\(\pi_v^{t/2}\) exchanges every conjugate block while
commuting with the full power \(\pi_v^t\), so all
relative slope vectors are self-dual. For \(E_0\) CM
the mixed fraction in each reduced mixed orbit is
\(2/t\) by the one-pair calculation. For \(E_0\) real
the restricted field action has order \(t/2\), which
is also the orbit length of the unique mixed label.
Thus the mixed fraction is again \(2/t\). Every
mixed block has reference slopes \(1/2-1/t\) twice
and \(1/2+1/t\) twice. All other blocks are flat of
slope \(1/2\).

There is just one mixed reduced orbit or a conjugate
pair of them, and excesses agree on conjugates.
Therefore a failure anywhere gives a positive
midpoint excess on every mixed block. The argument
on the integral crystalline lattice
[giving the Pfaffian parity](CM_NEWTON_RANK_FOUR_PFAFFIAN.md#even-valuation-from-crystalline-frobenius)
applies since all one-step determinant valuations are
even. For every mixed block the valuation of
\(\det(1+\pi_v^t/p^{t/2})\) on that block is strictly
greater than \(-2\) and even, hence nonnegative.
The flat blocks already give nonnegative valuations.
The norm is then a rational integer of absolute value
at most \(2^{8[F:\mathbb Q]}\). The cocharacter from
the preceding section gives nonconstancy on the
component, so algebraic Chebotarev excludes those
integer levels outside a density-zero set.
[Rajan, Theorem 3](https://repository.ias.ac.in/38317/1/2-PUB.pdf).

## A tensor obstruction on the other components

Suppose a tensor-standard block has its two supporting
\(A_1\) factors in a single orbit of the pinned action
\(\alpha=\operatorname{Int}(g_0)|_D\).
On that block the two middle slopes of \(g^t\) coincide.
This conclusion means equality of valuations of the
two middle roots at every nonarchimedean realization
of the compatible Frobenius polynomial.

If \(\alpha^t\) interchanges the factors, the
restriction of \(g^t\) is a scalar times
\((A\otimes B)\operatorname{flip}\) for invertible
two-by-two matrices \(A,B\). Indeed the flip implements
the factor exchange, and after composing with it the
remaining automorphisms of the two adjoint groups are
inner. Schur's lemma accounts for the scalar. After
conjugation by \(1\otimes B^{-1}\) the matrix is
\((AB\otimes 1)\operatorname{flip}\) up to the same
scalar. For eigenvalues \(u,v\) of \(AB\), its
eigenvalues up to scalar are
\[
             u,\ v,\ \sqrt{uv},\ -\sqrt{uv}.
\]
This follows first for diagonalizable \(AB\) and
then for all \(AB\) by the characteristic polynomial
identity. The last two valuations are the average of
the first two.

If \(\alpha^t\) preserves the two factors individually,
their common orbit length divides \(t\).
Write \(g=g_0x\) with \(x\in G^\circ\).
The factors of \(g_0^{-t}g^t\) in the adjoint quotient
are the ordered products of the \(\alpha\)-translates
of \(x\). Under the pinned identifications of the
two \(\mathrm{PGL}_2\) factors, these products are
cyclic rotations of the same orbit product, repeated
the same number of times. They are conjugate, and
\(g_0^t\) is scalar on this block since it centralizes
its two derived factors. The resulting spectrum up
to scalar has form \(u^2,uv,uv,v^2\), again with equal
middle valuations.

Any coefficient block with this property and a
positive reference midpoint kink fails to attain
the reference. In particular, if \(E\) is defined
over \(K\) and an actual monodromy component joins the
two supports of a block with multiplicities \((2,2)\),
its conjugacy class gives failure on a positive-density
set of primes. In that case block labels are fixed
under every component, and the reference at the
selected block is ordinary. This is a conditional
arithmetic obstruction. It does not supply an
example of such a disconnected arithmetic component.
Moonen's construction realizes the connected tensor
projections in Shimura families from subsets of
cardinality two.
[Moonen, Proposition 4.4, Construction 5.6 and Theorem 5.1](https://doi.org/10.4153/S0008414X24001159).

The slope pattern of the obstruction occurs in
actual geometric reduction. The genus-four curve
\[
     y^3=(x^3+x+1)(x^3-x^2-x-1)^2
\]
has good reduction at \(7\) and signature \((2,2)\)
for the deck action, by the branch and differential
calculation for the
[previous genus-four example](CM_NEWTON_RANK_FOUR_PFAFFIAN.md#an-explicit-genus-four-jacobian-over-mathbb-q).
With the same choice of cubic character at \(7\),
the character counts in degrees one to four are
\[
     (2,3,1),\ (16,13,17),\ (135,102,102),\ (832,833,733).
\]
They give the relative polynomial
\[
     (X^2+2(\zeta_3+1)X+7\zeta_3)(X^2-7\zeta_3).
\]
The first factor has slope list \(0,1\), since its
trace is a unit and its determinant has valuation
one. The second has slopes \(1/2,1/2\), so the
midpoint excess is \(1/2\).
This also has the displayed tensor-flip spectral
form. The spectrum at this prime alone does not
identify its global monodromy component.

Run the [exact PARI/GP replay](../scripts/verify_linked_rank_four_tensor_spectrum.gp)
from the repository root. It uses the same
[character-count functions](../scripts/verify_rank_four_pfaffian_curve.gp)
as the preceding genus-four example.
