# Odd dihedral covers split formally along the branch; the index can avoid it

15 November 2026. PF2RULED64 proves an unbounded reduction, not a new
quartic exclusion. Every odd cyclic etale sign cover reversed by the
sign involution splits over the **entire** ramification curve, with its
singular gluing retained, and over every finite thickening of that curve.
The degree-n reflection quotient, n odd, consequently has a polynomial
generator on an open neighbourhood of the whole affine branch. For the
actual quartic's cubic resolvent, a generator can approximate this formal
split model to any prescribed order. Its necessarily nonconstant index
can be chosen disjoint from the branch. Thus the obstruction left by
PF56 is global; it cannot be localized to a branch singularity or detected
by demanding agreement with this model to a larger finite order.

## 1. Hypotheses and odd torsor splitting

Let A=C[x,y], let G be nonconstant, reduced and irreducible, and set

    S=Spec A[z]/(z^2-G),   tau(z)=-z,
    C=V_S(z)=Spec A/(G).

Suppose T->S is a connected finite etale cyclic torsor of odd degree n>=3, with
class chi in H1_et(S,Z/n). Assume it has an involution sigma over tau
which conjugates the cyclic action by inversion. Then

    tau^*chi=-chi.

This is the dihedral hypothesis; an arbitrary cyclic cover without this
action is not being substituted. Restrict to j:C->S. Since tau*j=j
as morphisms of schemes, including at singular points of C,

    j^*chi=-j^*chi,   2*j^*chi=0,   n*j^*chi=0.

Integers 2,n are coprime, so j^*chi=0. Consequently

    T x_S C = disjoint union of n copies of C.       (1)

This is a statement about the torsor, not just its associated line
bundle on the normalization of C. In particular no node or tangency
gluing character survives. The proof does not require C to be smooth,
immersive, rational, or to have normalization A1.

For the sign surface, PF55 gives O(S)^*=C^* and identifies Kummer
torsors with Pic(S)[n]. Its norm argument shows that tau reverses
these classes. Thus the same splitting holds for every odd Kummer
torsor on this S, even when Pic(S)[n] itself is nonzero. Alternatively
the actual quartic supplies the involution and torsor directly in PF54.
It would be incorrect to infer Pic(S)[n]=0 from (1).

## 2. All infinitesimal orders, not just the reduced branch

For m>=1 put A_m=A/(G^m) and

    S_m=Spec A_m[z]/(z^2-G).

The closed immersion C->S_m has nilpotent ideal (z), since z^(2m)=0.
Etale morphisms and their maps lift uniquely across nilpotent
thickenings: Stacks Project Theorems41.15.1-2, tag06NE, or the affine
Lemma15.11.2, tag0ALI. Hence (1), with a fixed labelling of components,
lifts uniquely to

    T x_S S_m = disjoint union of n copies of S_m.   (2)

The isomorphisms for different m are compatible after fixing the
labelling on C. We use only these finite statements and their inverse
limit, not an algebraization theorem for a formal section.

On the labels i in Z/n, sigma acts by i->a-i. Since n is odd, relabel
by the unique fixed point a/2. There is one fixed label and (n-1)/2
pairs interchanged by sigma. On a fixed component sigma acts as tau;
on an interchanged pair it swaps components while applying tau to
their coefficient ring. There are no additional automorphisms over
S_m of a labelled component. The resulting reflection quotient B has

    B/G^m B = A_m x (A_m[z]/(z^2-G))^((n-1)/2).    (3)

Here and below the products are products of algebras. Formation of
invariants commutes with reduction because 2 is invertible: averaging
by (1+sigma)/2 proves surjectivity of the reduction on invariants.
The fixed component contributes A_m, since the tau-invariants of
A_m direct-sum A_m*z are A_m. Each swapped pair contributes one
copy of the quadratic algebra. This proves (3) globally along the
whole branch, not just after strictly local completion at one point.

If T is connected, its composite over Spec A has dihedral Galois
group of order2n. The reflection quotient is a normal integral
finite A-algebra of rank n. It is locally free: a normal surface is
Cohen--Macaulay, and a finite maximal Cohen--Macaulay module over the
regular surface A is locally free. This also applies to the actual
normal cubic algebra already retained in PF56. Singularities of S
are allowed; etaleness of T->S there is essential.

## 3. A global element generates near every branch point, in all odd degrees

Choose distinct constants c0,c1,...,c_h, h=(n-1)/2. In (3) use

    t_formal=(c0,c1+z,...,c_h+z).

Its formal characteristic polynomial is

    (W-c0) * product_i ((W-ci)^2-G).                 (4)

Modulo G the factors (W-c0), (W-ci)^2 are pairwise comaximal: their
resultants are nonzero constants. Each quadratic factor has generator
ci+z, so the Chinese remainder theorem proves that t_formal generates
B/GB. Equivalently, its power-basis determinant is a unit modulo G;
it is the confluent Vandermonde determinant, up to sign,

    product_(i<j) (cj-ci)^(e_i*e_j), e0=1,e_i=2.    (5)

The factors remain comaximal over A_m: a polynomial in G with
nonzero constant term is a unit in A_m. Thus t_formal generates (3)
at every m. Lift its image in B/G^m B to t_m in B. The determinant
line of B is trivial because Pic(A2)=0. Relative to a trivialization
let I(t_m) be the determinant of 1,t_m,...,t_m^(n-1). Then

    I(t_m) is a unit modulo G,
    V(I(t_m)) intersect V(G) is empty,
    B[1/I(t_m)] = A[1/I(t_m)][t_m].                 (6)

The last statement follows by the invertible power-basis matrix.
In particular I(t_m) is nonzero and t_m is a primitive integral
element of the degree-n field. This is an all-odd-degree branch
localization theorem. It is not asserted for arbitrary Sd covers
or for even dihedral degree. Nor does it say I(t_m) is a unit on
the entire plane.

## 4. Exact cubic normal form for the surviving quartic

Use the actual cubic B=C3 of PF54/PF56. At every m, (3) becomes

    B/G^m B = A_m x A_m[z]/(z^2-G).

Take the trace-zero formal element

    t_formal=(2,-1+z).

Lift it to B and subtract one third of its trace. The subtraction
changes nothing modulo G^m. Its characteristic polynomial is

    F_m(W)=W^3-3*a_m*W-2*b_m,
    a_m=1+G/3 mod G^m,   b_m=1-G mod G^m.           (7)

Here the formal polynomial is exactly

    (W-2)*((W+1)^2-G)=W^3-(3+G)*W-2+2G.            (8)

We can also normalize the global determinant of B so its trace
discriminant is 4G. PF56 gives a nonzero constant times G, and the
constant has a square root in C. The standard formal product basis
has discriminant 4G too. The ratio of these determinant
trivializations therefore has square1. For this argument use the
G-adic completion: G remains a nonzerodivisor in the completed
coefficient ring, by flatness of noetherian completion. It is
connected because its reduced closed fibre A/(G) is a domain.
Its square1 units are consequently +1 or -1. Choose the orientation
once to make this ratio +1. The power-basis determinant of (2,-1+z)
in the product basis is 9-G, by a direct 3x3 determinant calculation.

Write the global determinant as J_m=9*I_m. We have, simultaneously,

    a_m^3-b_m^2=3*G*I_m^2,
    I_m=1-G/9 mod G^m,
    (I_m,G)=A.                                     (9)

Indeed disc(F_m)=108*(a_m^3-b_m^2)=4G*J_m^2.
The first congruence in (9) follows from the chosen determinant
orientation and the matching formal element. In particular I_m=1
on the whole branch, including all its singular points. PF56 applies
to the actual A1-normalized branch and excludes a global generator,
so **each I_m is nonconstant**. Thus an unavoidable index divisor
exists but can always be placed away from C. This sharpens the
earlier statement that an arbitrary primitive element's index
support might coincide with C; that earlier statement remains valid
for arbitrary choices of primitive element.

There is a different element t_m for each prescribed finite order.
No single trace-zero element can satisfy (7) at all orders: the
intersection of (G^m) in A is zero, so it would have the reducible
characteristic polynomial (8). This contradicts generation of the
cubic field. Formal splitting and arbitrary finite approximation
do not yield a global product or a global monogenic normal algebra.

## 5. Controls and the remaining global obstruction

The splitting theorem does not imply that the Picard torsor vanishes.
PF55's smooth Oka example

    G=(xu-2)^2-u^3

has Pic(S)=Z/3 and a genuine etale triple cover. Its branch is C*,
parametrized by u=s^2, xu-2=s^3, s!=0. The Kummer function xu-2+z
restricts to s^3 on the branch, so the triple torsor really does
split there. Its nontriviality is global. Its global monogenic cubic
is also allowed, since the normalization is C*, not A1.

PF56's A1 curve G0=x^2-(x+y^2)^3 supplies the opposite control:
the displayed monogenic torus cubic has a triple-root fibre at the
inner A5 point. Its sign triple is only quasi-etale there. The
finite-etale lifting theorem cannot be applied to it. Local class
group 3-torsion or a discriminant norm equation alone is insufficient.

Even degree is excluded from the theorem: 2*chi_C=0 and n*chi_C=0
no longer force chi_C=0 when n is even. A nonzero order2 class on
C* illustrates the failed arithmetic inference; no ambient even
dihedral counterexample is asserted by this illustration.

Finally (9) has an essential square factor I_m^2. PF56's theorem for
an irreducible pure power difference cannot be applied after dropping
it. The branch-avoiding index may carry collisions of a chosen
primitive element over the target's etale locus; those collisions
are not new branch values of the normal cubic. The retained
power-difference logarithmic proof needs a morphism from the whole
affine plane to its blowup complement. Removing the index divisor
changes the source and does not preserve that hypothesis.

For a useful next exclusion, relate the branch-avoiding index and
the actual cubic order to the full quartic's shared conductor and
marked A2 source (PF53/PF57/PF58), or prove a genuinely global
Picard/linking obstruction for the branch families of PF63. Neither
larger branch jets nor an assumed injective restriction
Pic(S)[odd]->Pic(C)[odd] can finish this route: restriction is zero
for every such torsor. Full quartic and higher-degree Keller
exclusions remain open.

## Evidence and assurance

[verify_formal_branch_index.py](verify_formal_branch_index.py) checks
the exact cubic determinants/discriminants, unrestricted symbolic
error divisibilities, CRT determinant controls in several odd ranks,
the odd-degree labelling and dropped-hypothesis controls. The proofs
above, not a finite rank or jet census, establish the universal
statements. [formal_branch_index_receipt.json](formal_branch_index_receipt.json)
pins this note, checker, predecessor evidence and the primary HTML
sources in [the manifest](sources/formal-branch-index/manifest.json).
The Stacks nilpotent invariance theorem is externally imported and
was read in its primary HTML form. No full independent reproof,
formal verification, external review, novelty priority, Picard
vanishing, finite Keller construction or whole-packet exclusion
is claimed.
