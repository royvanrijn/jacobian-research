# Fibre-transitive shear recursion and the global spin boundary budget

PF2RULED104, 16 November 2026. This gives an all-degree source obstruction,
a full-shear cable theorem at arbitrary parameters and length, and an
ACTUAL Keller exclusion for complete colored iterated families. Keep the
whole finite normalization and all omitted primes throughout.

**Source theorem.** A smooth finite normalization of a nonautomorphic
plane Keller map, with EVERY omitted transverse index odd and integral
homology-sphere end, has r>0 omitted primes with

    r = 0 modulo 8,   d >= r+3 >= 11.                 (1)

More precisely, for any retained list of distinct omitted primes with
known (e_i,f_i),

    r <= d-1-sum_i(e_i*f_i-1).                       (2)

Thus the source packet is excluded whenever the interval from the number
of known primes (and at least one) to the upper bound in (2) contains no
positive multiple of eight. The r>=8 conclusion requires definiteness,
smoothness, every index odd, and H1(end)=0, not merely finite or cyclic H1.

**Cable theorem.** For a connected companion cover, let each joint
peripheral orbit have lattice (e_i,0),(a_i,f_i). Suppose

    gcd(e_i,f_i-a_i)=1 for EVERY orbit i.             (3)

Let M be a common multiple of the orders of both peripheral permutations.
For ANY coprime r,k>1 with r=k=1 modulo M, extend the cable representation
by z=mu*lambda. Its full meridian and longitude actions stay unchanged.
The exact integral transfer on each pattern component is

    m_old=r*m_new,
    L_new=r*L_old+c_i*m_new,
    c_i=a_i*(1-r^2)/e_i in Z.                        (4)

Moreover m_new,L_old are an integral H1 basis. At any finite length,
R=product of windings, the full rational linking rule is

    F_final = R^2 F_initial + diag((1-R^2)*a_i/e_i).  (5)

Closed H1 is preserved, and each normal-filled core class is multiplied
by R under that exact identification. The old zero-shear rule is the
special case a=0; dropping the diagonal term at nonzero shear is false.

**Actual applications.** Start from the natural eight-sheet or regular
twenty-four-sheet SL(2,F3) trefoil representation of
[PF101](FULL_PERIPHERAL_RATIONAL_LINKING_AND_LONGITUDE_JOINING_RECURSION.md).
Take ANY finite list of coprime ordinary cables (r_j,k_j)>1 with
r_j=k_j=1 modulo6 and the indicated extension z_j=mu_j*lambda_j.
If this is the FULL marked cover at infinity of the finite branch curve
of an actual Keller normalization, then:

* the entire natural EIGHT-sheet family is excluded when N is smooth,
  allowing ANY extra unramified omitted primes and arbitrary finite
  image singularities;
* the entire regular TWENTY-FOUR-sheet family is excluded even without
  smoothness, by the global residue budget alone.

The empty cable list is allowed. Nonempty lists are genuine iterated
knots. The parameters are independent and unbounded; for example (7,31)
and (7,43) both qualify, and k can be arbitrarily large. This is a
statement about the specified COLORED covers. Other representations of
the same uncolored knots are not excluded. No algebraic or Keller
realization of every knot in the parameter family is asserted.

## 1. The global even-lattice budget

Use the actual smooth N, U=A2, and ALL disjoint omitted embedded lines
D_1,...,D_r. [PF103](MULTIPLE_BOUNDARY_CANONICAL_CONTACT_AND_A_TYPE_END_EXCLUSION.md)
gives H2(N)=Z^r, H1(N)=H3(N)=0, negative definite compact intersection
matrix Q in the basis dual to D_i, and H1(W)=coker Q. Its canonical
characteristic evaluations are k_i=e_i-1, with Q_ii=k_i modulo2.
Every odd e_i makes Q even. If H1(W)=0, Q is unimodular.

Apply van der Blij's characteristic/signature congruence to B=-Q.
The characteristic vector zero is available for an even lattice;
its square is zero and the positive definite signature is r. Hence
r=0 modulo8. The complete three-page primary paper, including Theorem2
and its unimodular specialization (7), was rendered and read. This
classical lattice theorem is imported, not independently reproved here.
No Rokhlin invariant, Donaldson theorem, or contact classification is
needed for this step.

A nonautomorphic Keller map has a ramified omitted prime, by purity and
the finite-etale/plane-unit argument in PF103. Thus r>0. Under the
all-odd assumption some e_i>=3. The global finite-flat Euler ledger of
[PF2OMIT2](../jc2-omitted-divisors/GLOBAL_LEDGER.md) and its residue refinement
[PF2OMIT6](../jc2-omitted-divisors/RESIDUE_LEDGER.md) give

    d-1 = sum_(all D) e_D*f_D + sum_(D,q) eta_q,
    eta_q=ell_q-e_D*q_D(q)>=0.                       (6)

Here ell_q is the ambient cluster length; q_D(q) is the local degree of
the polynomial map D^nu=A1 -> C_D^nu=A1. To recall why the residue is
not inserted twice: the prime Euler excess is sum(ell_q-e_D), and
Riemann--Hurwitz gives sum(q_D(q)-1)=f_D-1. Substituting yields (6).
The underlying Euler proof uses chi(N)=1+r and finite flatness, not
generic branch smoothness. Every omitted prime, including e=f=1,
must occur in the sum. Extra finite image singularities only contribute
nonnegative eta; they cannot pay a negative budget.

Each omitted prime costs at least one. Retaining known costs e_i*f_i
in place of those ones proves (2), even if there are further unknown
ramified primes. Some odd ramified index costs at least three, so
d-1>=r+2. Together with r a positive multiple of eight this proves (1).

The necessary conditions do not assert nonexistence of even unimodular
lattices at every rank. The positive E8 lattice at r=8 is a valid
control. Formal indices (1,...,1,3) and f=1 pass the numerical budget
at d=11. This is NOT a surface or a finite Keller realization. Conversely
an odd rank-one unimodular lattice, or an even rank-one determinant-two
lattice, does not meet the hypotheses of the eight-divisibility rule.

## 2. Fibre-transitive covers of a cable space

Use the root presentation and preferred slopes already audited in PF101:

    P=<m,l,z | [m,l]=1, z^r=m^k*l^r>,
    m_new=m^u*l^v*z^(-v), r*u-k*v=1,
    l_new=m_new^(-r*k)*z^r.                          (7)

The regular Seifert fibre is h=z^r=m^k*l^r, which is central. Since
gcd(k,r)=1, the outer torus lattice has an integral basis h,b. Equivalently

    P=<b,z,h | [b,h]=[z,h]=1, z^r=h>.                (8)

Consider a CYCLIC connected N-sheet cover with h acting as an N-cycle.
Choose its character chi:P -> Z/N so chi(h) is a generator. This forces
gcd(r,N)=1, since chi(h)=r*chi(z). Choose integers A,B so that
z'=z*h^A and b'=b*h^B belong to the kernel. Put H=h^N. The kernel has
the presentation

    <b',z',H | [b',H]=[z',H]=1,
                  (z')^r=H^t>,  t=(1+r*A)/N.        (9)

Here 1+r*A is divisible by N and gcd(r,t)=1, since N*t-r*A=1.
To check generation, any word may have its character cancelled by a
central h power. Modulo h this projects onto the same base orbifold
group <b,z | z^r=1>, while the remaining kernel is <h^N>. Lifting its
single cone relation gives exactly (9), with no extra relation. This
also describes the ordinary fibrewise cover directly, without importing
a classification of closed Seifert manifolds.

Abelianizing (9) gives three generators and the primitive relation
r*z'-t*H=0. Thus

    H1(P_cover;Z)=Z^2, with NO torsion.               (10)

The pushforward to H1(P;Z)=Z^2 is injective: transfer makes it rationally
surjective, equal ranks make it rationally injective, and (10) has no
torsion. Its image is the kernel of the
surjective abelian character, hence has index N. The reverse inclusion
holds because every zero-character loop lifts closed.

Assume the old and new boundary lattices are the SAME (e,0),(a,f),
with N=e*f. The ordinary cable-space homology has basis m_new,l_old
and relations m_old=r*m_new, l_new=r*l_old. Pushforward of the actual
lifted ports is

    m_old -> e*r*m_new,      m_new -> e*m_new,
    L_old -> a*r*m_new+f*l_old,
    L_new -> a*m_new+f*r*l_old.                      (11)

The two columns for m_new,L_old have determinant e*f=N, exactly the
index of the pushforward image. Therefore they are an INTEGRAL basis
upstairs. Injectivity now proves (4), including integrality of c_i.
This is the all-parameter proof; matrix examples are only controls.

## 3. Recognize the hypotheses from full peripheral permutations

On a transitive peripheral orbit with lattice (e,0),(a,f), the element
m*l has order

    e*f/gcd(e,f-a).                                 (12)

Indeed (m*l)^q closes only when q=f*s and e divides s*(f-a).
Thus (3) says this element is a single e*f-cycle. Let M be as in the
theorem, r=k=1 modulo M, and assign z=m*l. The root relation holds;
from r*u-k*v=1 we have u-v=1 modulo M. Formula (7) then gives

    rho(m_new)=rho(m), rho(l_new)=rho(l).             (13)

All joint peripheral lattices are preserved, including shear. Each
covered pattern component has one old and one new boundary torus and
cyclic action generated by rho(h)=rho(m*l), a full cycle by (3).
Hence Section2 applies. Connectedness of the whole cover persists
because its old transitive image is retained.

The congruence restriction allows independent unbounded parameters:
for any r>1 with r=1 modulo M, k=1+M*r*T is coprime to r and is 1
modulo M for every integer T>=1. It is not restricted to a fixed
linear difference between r and k. Any additional lower bound on k
can therefore be met. No infinity/Puiseux realization is inferred.

## 4. Integral gluing, full shear, and arbitrary length

For each component use the basis in (10)--(11). Mayer--Vietoris with
the ACTUAL old cover has one pattern edge per old boundary, so introduces
no extra incidence cycle. Eliminating L_old gives exactly

    H1(E_new)=(H1(E_old) plus sum_i Z*m_i,new)
                    /<m_i,old-r*m_i,new>.            (14)

After normal filling m_new=0, (14) identifies closed H1 exactly with
old closed H1; (4) gives [K_i,new]=r*[K_i,old]. This keeps unrelated
torsion and cores of finite or infinite order. For a finite core order
q_i, its new order is q_i/gcd(q_i,r).

If all cores are rationally null, substituting L_old=F_old*m_old into
(4) yields

    F_new=r^2*F_old+diag(a_i*(1-r^2)/e_i).            (15)

In particular F_new*e=a, as required by PF101, rather than r^2*a.
Iteration telescopes to (5), and [K_i,final]=R*[K_i,initial]. There is
no upper bound on length, coefficients, degree or number of components.
The rule is conditional on (3), stable profiles and the stated coloring;
it is not claimed for arbitrary cable covers.

The torsion control (e,f,a)=(3,1,1), (r,k)=(4,7) is instructive. Its
regular fibre is NOT transitive, and the lifted cable space has
H1=Z^2 plus (Z/4)^2. Even when the displayed port differences vanish,
m_new,L_old leave those extra torsion classes. Omitting hypothesis (3)
would therefore invalidate the closed-homology transfer.

## 5. Complete eight- and twenty-four-sheet colored exclusions

The natural PF101 seed has two boundary orbits (1,2,0),(3,2,0), trivial
filled H1, and F=[[-9,3],[3,-1]]. The regular seed has four (3,2,0)
orbits, trivial filled H1, and F with diagonal -3 and off-diagonals 1.
Both satisfy (3), and M=6 works. Equations (5), (13)--(14) retain these
profiles and H1=0 at EVERY finite length and all the stated parameters.
The old PF101 Z[t] certificates are preserved as independently checked
overlap controls; their fixed subfamily is not the universal proof here.

Now suppose such a cover is the full infinity cover of the finite branch
curve of an ACTUAL plane Keller map. The branch link is the single
specified knot; its original affine curve has one place at infinity.
[PF98](PROPER_FINITE_MARKED_END_AND_COMPONENT_OPEN_AUDIT.md) identifies
the original finite radial end with its normal-filled cover. Every
ramified component lies outside U and gives one of the e>1 joint
boundary tori; its longitudinal degree is the residue degree of its
A1 normalization. A ramified affine curve has an end and cannot be
hidden at a finite value. Extra unramified omitted primes are allowed,
including primes whose images differ from the branch curve; deleting
them and retaining finiteness is NOT permitted.

For the NATURAL EIGHT-sheet family there is precisely one ramified
omitted prime R, (e,f)=(3,2). The other profile (1,2,0) cannot itself
be an omitted e1 prime: PF2BATT1 gives f=1 on every such prime. It may
be an ordinary preimage curve meeting U and is not counted as a second
omitted prime. All other omitted primes, if present, have index one.
Thus ALL omitted indices are odd. With N smooth, (1) requires r>=8,
while (2) with the actual R cost6 gives

    r <= 8-1-(6-1)=2.                               (16)

This excludes the whole family, WITHOUT assuming r=1, zero excess,
an immersive finite branch parametrization, or a finite-image degree.
Already the general bound d>=11 excludes this eight-sheet source;
(16) retains the sharper residue cost.

For the REGULAR TWENTY-FOUR-sheet family all four ramified boundary
orbits must be omitted, each costing e*f=6. Equation (6) would require

    24-1 >= 4*6=24,                                 (17)

which is impossible. This argument uses the normal affine boundary
theorem and residue budget, so does not require smooth N. Additional
unramified primes or finite exceptional excess only increase the cost.

These are ACTUAL Keller compatibility exclusions. The normal-filled
topological covers themselves exist; (17) is not a contradiction for
an arbitrary finite branched cover whose marked complement is not A2.
We do not infer that their filled 3-manifolds are S3 from H1=0.

## 6. A genuine finite polynomial countercontrol

To prevent a false extension from the marked plane to arbitrary finite
maps, consider on C2_(x,y)

    P=3x^2+12x*y^3+8y^6, Q=x*y+y^4.

The map is finite of generic degree8, and

    Jac(P,Q)=6x^2,
    P^2-64Q^3=x^3*(9x+8y^3),
    y^8-6Q*y^4+P*y^2-3Q^2=0.                       (18)

The monic equation makes y integral; the monic quadratic equation for
x from P then makes the map finite. Its degree-eight y equation is
irreducible over C(P,Q): over C(Q)[y], its coefficient of P is y^2,
coprime to y^8-6Q*y^4-3Q^2, so the polynomial primitive and linear in P
is irreducible by Gauss's lemma. Also x=(Q-y^4)/y generically.

The smooth ramification line x=0 has e=3 and residue degree2 onto the
cusp (P,Q)=(8y^6,y^4). At zero its local degree is8; the reduced
unramified companion 9x+8y^3=0 meets it there. Weighted dilation
(x,y)->(t^3*x,t*y) makes (P,Q) have weights(6,4). The proper radial
level is diffeomorphic to S3: both Euclidean and target-norm levels
meet each nonzero weighted ray once and transversely. In particular
its filled H1 is zero. No identification of its precise permutation
labeling with the retained PF101 seed is needed for this control.

Here N=C2 and N minus {x=0}=C* times C, NOT the marked Keller source
A2. The ramification line is principal; the Thom rank conclusion and
the Euler ledger with chi(U)=1 do not apply. The finite map and the
topological data are real, but it is not a Keller counterexample. This
also blocks a tempting unrestricted assertion that a smooth critical
line always has primitive image parametrization: rank-zero points can
carry residue ramification at index3.

## 7. Verification and the actual frontier

[The checker](verify_fibre_transitive_spin.py) checks integer lifted port
rows, exact H1/core torsion, both covering trees, nonzero shear, failure
when fibre transitivity is omitted, large binary exponents, mixed cable
parameters, the classical E8 countercontrol, exact budget identities,
and the finite polynomial map (18). The
[receipt](fibre_transitive_spin_receipt.json) separates the universal
written proofs and imported lattice theorem from finite controls.
The [source manifest](sources/fibre-transitive-spin/manifest.json) pins
the full primary lattice paper and retained root/peripheral sources.

The general all-degree Wright goal remains open: other colored covers,
singular eight-sheet sources, even omitted indices, nontrivial end H1,
and the ACTUAL Laurent-terminal finite-algebra implication are not
settled. Published quartic retirement remains intact. No parameter
campaign, formal verification, external review or outreach is claimed;
closure precedes external-review expense.
