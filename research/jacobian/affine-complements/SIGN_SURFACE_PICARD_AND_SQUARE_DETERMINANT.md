# Picard torsion and a square determinant obstruction for the sign surface

15 November 2026. PF2RULED55 strengthens the whole-quartic obstruction:
the knot at infinity must have determinant divisible by **9**. The proof
has no bound on source pole orders, target image degrees or cable stages.
It also closes the units alternative left explicitly open in PF2RULED54:
the actual sign surface must have nonzero Picard 3-torsion. The underlying
square-divisibility theorem concerns affine double surfaces in arbitrary
branch degree; the required triple cover is supplied here only in the
surviving quartic Keller packet. This does not exclude the whole packet.

## 1. Scope and retained inputs

Use [PF2RULED54](QUARTIC_COLOURING_AND_ETALE_TRIPLE_GATE.md), with its
original source A2, both omitted primes R,D, sole ramified index2,
residue degrees one, zero excess, immersive normalization of the full
branch curve C, and actual S4 monodromy. Its normal sign surface is

    S = Spec C[X,Y,z]/(z^2-G(X,Y)),   G irreducible.

C has normalization A1, so its sufficiently large affine link is a
knot K. Its finite singular values have exactly two smooth branches,
of contact m>=1. The corresponding sign singularity is A_(2m-1).
PF54 constructs a connected cyclic degree3 cover T->S that is etale
at **every** affine point, including those singularities. Its class
in H1_et(S,mu3) is nonzero and reversed by z->-z. These are actual
global cover inputs, not assumptions deduced from a knot colouring.

The new topological argument needs only that S is normal affine of
complex dimension2, has finitely many singularities whose resolutions
have exceptional trees of rational curves, and has a rational homology
sphere as its single outer boundary. In the displayed double-cover
setting the knot hypothesis gives that boundary condition. Thus the
topological theorem is not intrinsically restricted to field degree4.

## 2. All units are constants: the Kummer alternative closes

Ford, *The Group of Units on an Affine Variety*, Theorem2.10, applies
to C[x1,...,xm,z]/(z^p-f) for irreducible f and prime p in characteristic
zero. Its normal, potentially singular, setup is on complete PDF7;
Theorems2.3 and2.8(c), which prove2.10, are on PDF8,11-12. In particular,

    O(S)^* = C^*.

For this quadratic case there is also a short independent argument.
Write a unit as a+bz, with a,b in A=C[X,Y]. Its norm a^2-Gb^2 is a
nonzero constant, so rescale to norm1. If a is constant then b=0:
otherwise G times a square would be a nonzero constant. If any
nonconstant units exist, choose one with deg(a)>0 minimal. In the UFD A,

    (a-1)(a+1)=G b^2,    gcd(a-1,a+1)=1.

Exactly one of the two factors has odd G-valuation; every other
irreducible valuation in either factor is even. Therefore one factor
is c*u^2 and the other is d*G*v^2, c,d nonzero constants. Their
difference is 2 or -2. Absorb square roots of c,d and rescale to get
a unit of norm1 with scalar part of degree deg(a)/2>0. This contradicts
minimality. This argument uses irreducibility, not smoothness.

The Kummer sequence retained in PF54 now gives, for every n>=1,

    H1_et(S,mu_n) = Pic(S)[n],                         (1)

since every constant has an nth root. Thus **Pic(S)[3] is nonzero**.
By complex finite-cover comparison, the left side is
Hom(H1(S,Z),Z/n). This comparison applies to singular complex varieties;
it does not require that S is a manifold. The double involution acts
as -1 on Cl(S), hence also on its Cartier subgroup Pic(S): for a Weil
divisor E, pi^*pi_*E=E+tau(E), including ramified primes with their
multiplicities, and pi_*E is principal on factorial A2. Thus all of
Pic(S)[3] is anti-invariant. It is actual Cartier Picard torsion that
is required; local Weil class torsion alone does not supply an etale cover.

## 3. Resolution preserves H1 and makes H2 free

Let p:U->S be a resolution, and E its exceptional locus. Above each
singular point E is a connected tree of rational 2-spheres. It has
H1(E)=H3(E)=0 and H2(E) free. Write F=Sing(S). Collapsing each exceptional
tree to its point identifies U/E componentwise with S; equivalently
the pairs (U,E) and (S,F) have the same relative homology. This follows
also by excision using the common boundary of a resolution plumbing
and a conical singular neighbourhood. These are good triangulable
pairs, so there is no quotient-space exactness assumption missing here.

Hamm, *Zum Homotopietyp Steinscher Raeume*, Satz1 and its corollary,
printed p121, applies to the possibly singular Stein space S with
empty analytic subset A. It gives

    H3(S;Z)=0,   H2(S;Z) free abelian.                (2)

Here Z denotes integer coefficients. Finite generation also follows
from algebraic triangulation. The relative exact sequences, with
the H0(E)->H0(F) identification, give

    H1(U;integers) = H1(S;integers),
    H3(U;integers) = 0,
    0 -> H2(E) -> H2(U) -> H2(S) -> 0.              (3)

The last sequence splits as abelian groups, since H2(S) is free.
Consequently H2(U) is free. Neither U's affineness nor any Stein
structure on U is asserted or needed.

Choose a sufficiently large target ball transverse to C and beyond
all stratified critical radii. Its preimage, resolved in the interior,
is a compact connected oriented smooth four-manifold W with interior
end homotopy equivalent to U. Outside this ball the algebraic pair is
a product by stratified triviality; the finite double cover lifts the
product. Thus W is homotopy equivalent to U and

    M = boundary W = Sigma_2(K),
    |H1(M)| = D = |Delta_K(-1)| > 0.                 (4)

The determinant-cover identity is retained from Przytycki Theorem3.3
in PF54. A knot has odd, nonzero D because Delta_K(1)=1 and reduction
modulo2 identifies evaluations at1 and-1. Thus M is a rational homology
sphere and H2(M)=0. No general position at projective infinity is used.

## 4. The full square, not just one factor of three

Poincare-Lefschetz duality and the universal coefficient theorem,
using (3), give

    H1(W,M) = H^3(W) = 0.

Hence H1(M)->H1(W) is onto, and A=H1(W)=H1(S) is finite. Let
L=H2(W), a free lattice, and let Q be its integral intersection form.
Over Q the form is nonsingular since M is a rational homology sphere.
There is an exact sequence

    0 -> L -> H2(W,M) -> H1(M) -> A -> 0,
    H2(W,M) = H^2(W),
    Tor H^2(W) = Ext(A,Z),   |Tor H^2(W)|=|A|.       (5)

Modulo torsion, the first map in (5) is Q:L->L^*. Since Q is nonsingular,
the image of L meets no torsion. It follows that

    |coker(L->H2(W,M))| = |det Q| * |A|,
    D = |det Q| * |A|^2.                            (6)

The determinant of the empty rank-zero form is1. This also derives
the precise torsion correction in Owens-Strle Lemma2.1 (complete PDF2-3).
In general one cannot replace its torsion image by the whole H1(W);
the freeness proved in (3) is what permits that replacement here.
Negative definiteness, rationality of a projective compactification,
and a guessed boundary intersection matrix are unnecessary.

The required etale triple gives A->>Z/3, so 3 divides |A|. Equation(6)
therefore proves the new actual exclusion:

    9 divides D for EVERY surviving quartic packet. (7)

It excludes both v3(D)=0 and v3(D)=1. The latter is strictly stronger
than PF54's determinant test. In particular T(3,4), which had an actual
parabolic S4 lift in PF54, now fails the required global sign-surface
extension. It is a topological control, not a new polynomial candidate.

More generally (1) and finite complex-cover comparison identify the
finite torsion Picard group, abstractly, with Hom(A,Q/Z). To see that
the entire torsion group is finite, choose n divisible by the exponent
of A. The Kummer groups and their transition maps identify Pic(S)[mn]
with Pic(S)[n] for all m, so no further torsion appears. Thus

    |Pic(S)_tors| = |A|,
    |Pic(S)_tors|^2 divides D.                      (8)

In particular dim_Fell Pic(S)[ell]=r implies ell^(2r) divides D.
This is a theorem about these affine double surfaces in every branch
degree. It supplies a Keller exclusion in another field degree only
after a genuine etale cover of its actual sign surface is established;
there is no assertion that S_d with d>=5 has the S3 quotient used in PF54.

There is also a useful group refinement of the numerical gate. The
image B of Tor H^2(W) in H^2(M) is injective by (5), has order |A|,
and is isotropic for the boundary linking form. Indeed for torsion
alpha,beta in H^2(W), lift alpha by the Q/Z Bockstein to a in
H^1(W,Q/Z). The linking value of their restrictions is the evaluation
of i^*(a cup beta) on [M]; it vanishes because i_*[M]=0 in H3(W).
Therefore the linking form of Sigma_2(K) must contain a nonzero
isotropic order3 subgroup. A 3-primary group of order9 can fail this:
on (Z/3)^2 the diagonal form (1/3,1/3) has no nonzero isotropic vector,
whereas (1/3,-1/3) does. Orders alone do not settle the refined gate.

## 5. Closed propagation through all cable stages

Keep PF54's actual knot convention: C_(w,p)(J) has winding w>=2 and
coprime signed other parameter p. The proved recurrence is

    D_new = |p|  (w even),
          = w*D (w odd, p even),
          = D   (w odd, p odd).

Consequently the new gate v3(D)>=2 propagates by

    e_new = v3(p)           (w even),
          = v3(w)+e         (w odd, p even),
          = e               (both odd).

Every final even-winding stage with 9 not dividing p is excluded,
regardless of all inner stages. An odd-winding/even-p stage must have
v3(w)+e>=2. Any number of odd/odd stages preserves the condition.
Passing it proves neither a global S4 representation nor a Keller map.
The sign and Picard gates retain all finite affine relations and the
unramified omitted prime D; no factor has been deleted to build a cover.

## 6. Exact controls and the remaining question

The companion [checker](verify_sign_square.py) verifies predecessor
seals and pinned sources, the determinant recurrence against exact
satellite Alexander polynomials, finite torsion/linking-form controls,
and two algebraic warnings. It is not a formal proof of the topology.

First, irreducibility matters in the unit proof: G=x^2-1 admits the
nonconstant norm-one unit x+z. Second, Pic and Cl must be separated.
For the normal cone uv=w^3, localization at u gives the UFD
C[u,u^-1,w], and the only boundary prime E=(u,w) has div(u)=3E.
Thus Cl=Z/3. At the vertex E is not Cartier: its ideal (u,w) needs
two generators by Nakayama, and its inverse class2E cannot be Cartier
either. Hence Pic=0. The cover u=s^3,v=t^3,w=st is etale off the
vertex but is ramified at it. Local class torsion cannot replace the
everywhere-etale requirement in PF54 or the nonzero Pic[3] in (1).

The actual smooth Oka sign surface from PF54 supplies the positive
Picard control: z^2=(xu-2)^2-u^3. Put a=xu-2+z and b=xu-2-z.
Then ab=u^3 and x=(a+b+4)/(2u). On u!=0 the surface is the torus
Spec C[u^+-1,a^+-1], with Pic0. Its boundary consists of two affine
lines L_+,L_-, and div(u)=L_++L_-, div(a)=3L_+ (choose L_+ by a=0).
The localization presentation matrix is [[1,3],[1,0]], with cokernel
Z/3 and zero integral kernel. Thus O(S)^*=C^* and Pic(S)=Z/3;
the involution swaps the two lines and acts by -1. Its branch
normalization is C*, with two ends. It is not an A1-normalized
quartic branch curve, and the knot-boundary step is not applied to it.

The remaining global question is whether an A1-normalized immersive
branch curve compatible with the marked Keller source can have the
required nonzero Pic[3], with its infinity linking form and primitive
finite algebra simultaneously satisfying the retained constraints.
Square divisibility and isotropy do not resolve this for all knots.
Higher field degrees also require a new source of global cover data;
the topology here does not manufacture one. No coefficient search is
substituted for these geometric questions.

## Sources and assurance

Exact bytes, inspected complete pages and statements are pinned in
[the source manifest](sources/sign-square/manifest.json). Hamm's original
printed p121 is retained as a complete GDZ scan, not inferred from a
smooth-only statement. Ford's hypotheses and the complete unit-theorem
proof chain were inspected. Owens-Strle Lemma2.1 and its entire proof
were inspected; (5)-(6) spell out its specialization here. Kerner's
recent paper was used to locate the classical singular theorem and is
retained as background; its new relative theorem is not imported.
Kummer and knot-cover comparison facts are retained via PF54; the
singular complex finite-cover comparison is also pinned directly.

This is a written unbounded proof with exact finite controls and
explicit primary theorem imports. It is not formal verification,
an independent full reproof, external review or a novelty claim.
No finite Keller pair, full quartic exclusion, higher-degree exclusion
or completion of the research goal is asserted.
