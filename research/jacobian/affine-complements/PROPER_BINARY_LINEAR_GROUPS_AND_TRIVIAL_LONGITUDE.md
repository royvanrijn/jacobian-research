# Proper binary linear groups force trivial preferred longitude

This is PF2RULED131. It applies to EVERY genuine primitive block quotient
of the monodromy of an actual polynomial plane Keller map, including
imprimitive original monodromy. It makes the common-form hypothesis in
[PF126](SEIFERT_FRAMING_CANCELLATION_AND_RESIDUAL_TWO_COVERS.md) automatic
for every PROPER binary linear group. It excludes nonidentity-longitude
profiles of those groups; it does not exclude the groups themselves.
The full finite-target, trace-energy and EMPTY ORIGINAL fibre gates remain
OPEN. The full research goal remains ACTIVE.

## 1. Actual geometric inputs and the new classification input

Retain [PF122](HEREDITARY_BLOCK_QUOTIENT_MONODROMY_REDUCTION.md) and
[PF123](INTERMEDIATE_CRITICAL_AXIS_DESCENT_AND_AFFINE_CUSP_PROFILES.md).
If a primitive quotient is binary affine, write

    Q = V : H,  V = F2^k,  b = 2^k,  H <= GL(V).

Here H is irreducible. There is ONE full quotient branch C, its meridian
has half support, and after translating the origin it is

    q(x) = x + l(x)v,   l != 0,   v != 0,   l(v)=0.

Its linear part is a transvection. Meridians normally generate Q, hence
their linear images generate H by conjugate transvections. Invisible
ORIGINAL branch components have identity image in Q and are retained in
the original normalization. Nothing here asserts that the intermediate
surface is a plane or that its involution lifts to an involution upstairs.

We now explicitly IMPORT the following narrow consequence of McLaughlin's
classification: an irreducible transvection-generated PROPER subgroup of
GL_k(F2) has even k>=4 and preserves a common nondegenerate alternating
bilinear form on its given module. GL=SL over F2.

The primary research restatement is Brumer--Kramer, *Arithmetic of division
fields*, arXiv:1102.4557v1, Section2, Proposition2.4 and Remark2.6, printed
pages2--3. The complete relevant PDF pages and the primary HTML statements
were read. Proposition2.4 concerns proper subgroups of SL(V); the explicit
symplectic-subgroup hypothesis belongs to the SEPARATE Proposition2.5.
The representation is essential: orthogonal/symplectic groups act in their
defining modules; symmetric cases use the even-sum permutation module,
modulo its constant radical when necessary. Remark2.6 specifies this
module, so this is not a deduction from abstract group names alone.

The original source is J. McLaughlin, *Some subgroups of SL_n(F2)*,
Illinois J.Math.13(1969),108--115, DOI10.1215/ijm/1256053741. Its original
PDF was NOT obtained: the publisher access returned a challenge page.
The imported result is adopted through the pinned primary research
restatement, not independently proved or independently source-reviewed.
No CFSG or finite-group census is newly performed here. We need only
the common-form consequence, not a full exact group-family catalogue.

## 2. Why the imported form cancels the longitude

Assume H is proper. Fix its common form B. For the chosen transvection,
preserving B implies l=B(v,-): for x in ker(l) and l(y)=1, comparison of
B(qx,qy) with B(x,y) gives B(x,v)=0. Nondegeneracy identifies the two
nonzero functionals over F2.

PF126 proves, in the ACTUAL Seifert framing, that the preferred longitude
fixes the moved half S={l=1} pointwise. Its entire affine permutation has
the form

    lambda(x) = x + (l(x)+1)u,  u in ker(l).

The form alone does NOT force u=0. The complementary transvection u=v
also preserves B and fixes S pointwise. The residual curve cover along
the WHOLE normalization C^nu supplies the necessary global information.

Let v_o count double values of C with odd contact. PF126 gives exactly
2*v_o finite branch parameters of the generically unramified residual
normalization R -> C^nu=A1, of total degree b/2. Each based branch
monodromy is a nonzero translation u_i of ker(l); its full affine extension
fixes S and has linear part I+u_i*l in H. The same preservation argument
therefore gives

    l = B(u_i,-) = B(v,-),  hence u_i=v for EVERY i.

Even contact contributes identity. The optional ordinary cusp contributes
an unramified section. There are no other finite residual branch points.
PF126 identifies the residual longitude with translation by sum_i u_i.
There are an EVEN number 2*v_o of these equal columns, so their sum is0.
Consequently the preferred longitude is the IDENTITY on ALL b letters.
This uses the global residual cover, not only affine centralizer algebra.

If v_o>0, its translation subgroup is <v>. Thus R is the disjoint union
of b/4 copies of the SAME double cover

    y^2 = product_(odd double values j) (t-a_j)(t-b_j).

All 2*v_o parameters are distinct. Every component has genus v_o-1 and
two points over infinity. If v_o=0, R is b/2 copies of A1. No equality of
paired columns was assumed before applying the common form.

This gives a uniform simplification of the peripheral recursion: every
critical orbit has (e,f,a)=(2,1,0), every unramified orbit has (1,1,0),
and the total number of peripheral tori is 3*b/4. These are permutation
and slope data. Actual rational linking entries and integral torsion
defects still require PF101. Identity longitude does not automatically
supply PF92's integral-null-core hypothesis or PF111's coverage condition.

## 3. The only remaining binary container for nonidentity longitude

If a genuine binary primitive quotient has nonidentity preferred longitude,
the preceding argument implies H=GL_k(F2), hence Q=AGL_k(F2). Also v_o>0,
since PF126 already proves identity longitude whenever v_o=0, without any
form hypothesis. The case k=2 has identity longitude by the retained
degree-four argument. Thus nonidentity longitude requires k>=3, full
AGL_k(F2), and at least one odd-contact double value. This is a NECESSARY
condition, not a construction or exclusion of that container.

For EVERY odd k>=3, the classification also gives H=GL_k(F2), irrespective
of the longitude. Proper linear survivors can occur only in even dimension.
For k=2 an irreducible transvection-generated subgroup is GL2(F2): two
distinct transvections generate it, whereas one fixes a line. The already
closed degree2 block quotient is not reintroduced.

For comparison, full GL_k(F2) has no common nonzero alternating form when
k>=3. Indeed all elementary transvections I+E_ij belong to it. An invariant
form B must pair e_i trivially with every e_a for a!=j. Given any a!=i,
choose j distinct from i,a; then B(e_i,e_a)=0. Alternation kills its diagonal,
so B=0. This proves the limitation of the form argument in all dimensions,
without computing a group list.

There is also an explicit algebraic control showing why one cannot simply
set every affine longitude to identity: in F2^3 take l=x_2, v=e_0 and
u=e_1. The two admissible residual columns v,u have nonzero sum v+u;
their span has order4. Their residual curve cover has genus0 and two
infinity points, and the full lambda=T_(v+u) fixes S pointwise but is
nonidentity. These permutations lie in AGL3(F2). This is ONLY a peripheral
and curve-cover control, not an actual polynomial Keller map or a global
geometric realization of all branch transports.

## 4. Checks, attribution and the remaining gate

The accompanying checker pins the primary PDF, its precise source review,
the retained geometric proofs and PF126 receipt. Cheap exact F2 controls
check invariant forms, the complementary-transvection warning, full-GL
failure, the odd/even deleted-permutation modules, all-degree cancellation
arithmetic and the full-AGL algebraic escape. The written proofs and explicit
import above establish the unbounded statements; small matrix controls
do not establish the classification. Registered predecessors are not rerun.

The user's suggested dihedral comparison is retained at
[the route map](DIHEDRAL_LINKING_RECURSION_ROUTE.md). PF92's unmodified
authors-code comparison applies on its odd-dihedral overlap. PF101 already
extends the recursion to nonidentity longitude and rational core classes.
No binary affine cover is relabeled as a dihedral cover, and no new knot
diagram expansion is needed for this reduction.

Natural alternating/symmetric quotients, full-AGL compatibility, the
identity-longitude proper affine survivors, original invisible primes,
and the actual finite-target/trace-energy/EMPTY ORIGINAL fibre problems
remain OPEN. No common affine quadratic refinement, original global
cyclicity, formal verification, independent replay, external review or
full-goal completion is asserted. Shared priorities, the CURRENT KR0/1/2
prefix and the published quartic PLANE retirement remain in force.
