# Trace linking, half ramification, and the missing cubic multiplication

2026-11-16. HC4MYLY0/1. ALL DEGREES. Full HC4/JC2 remains OPEN.

HC4 and the plane Jacobian conjecture remain open.
The remaining unproved implication is: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
We examine the finite discriminant module itself. Its trace linking form
is perfect globally, but a canonical half-ramification lattice makes the
odd-index part metabolic. A perfect trace form on a smaller finite lattice
therefore does not supply the full codifferent pairing. Cubic multiplication
data restore the missing condition; their regularity is NOT automatic.

## 1. LY0: the GLOBAL finite discriminant linking form

Let B=C[P,Q], F=Frac B, L/F finite separable, and S its finite normal
integral closure, locally free over B. This setup includes every ACTUAL
plane Keller normalization. Put

    Lambda={a in L:Tr(a*S) subset B}, D=Lambda/S.

The ordinary field trace gives a well-defined symmetric linking form

    beta:D times D -> F/B,
    beta([a],[b])=Tr(a*b) modulo B.                 (1)

It is PERFECT as a B-module linking form:

    D -> Hom_B(D,F/B) is an isomorphism.            (2)

Here is a global proof, including every finite point. The finite projective
trace duality gives Lambda^#=S, where M^#={a:Tr(a*M) subset B}. This proves
injectivity in (2). For surjectivity lift a homomorphism D->F/B to Lambda.
Since Lambda is B-projective, its map lifts through the surjection F->F/B
to a B-linear map Lambda->F. After extending scalars this map is Tr(c*-)
for a unique c in L. Its reduction vanishes on S exactly when c in Lambda.
Thus it comes from (1). No injectivity theorem for F/B or codimension-two
gluing assertion is assumed.

The module D is the cokernel of the finite trace Gram map S->S^vee. In a
local free basis it has a square symmetric matrix presentation. Its zeroth
Fitting ideal is the discriminant. It is Cohen--Macaulay of pure codimension
one when nonzero: the two projective lattices give projective dimension one
over the regular base surface. In particular there is NO independent
point-supported submodule. Perfect linking is not D=0: regularity of the
full pairing Tr(Lambda Lambda) is exactly beta IDENTICALLY ZERO, which by
perfection forces D=0. In the actual Keller setup this is the retained inverse
criterion. LX0 identifies the SAME D with every target-additive derivative
defect; its pairing is an additional structure, not a new vanishing theorem.

## 2. Canonical lower and upper half-ramification lattices

Let R_pi=sum_E(e_E-1)E be the finite tame ramification divisor on Z=Spec S.
Write O_Z(A) for the reflexive fractional S-module of functions with poles
bounded by an integral Weil divisor A. Relative to the pulled-back target
volume, the retained canonical-module identification is Lambda=O_Z(R_pi).
Define

    T_minus=O_Z(floor(R_pi/2)),
    T_plus =O_Z(ceil(R_pi/2)).                      (3)

They are finite S-modules and B-PROJECTIVE. To check the latter, a reflexive
module on a normal surface has depth two at closed points. A regular system
of parameters of the base is a system of parameters in each finite local
ring upstairs, and is regular on this maximal Cohen--Macaulay module. The
same depth/miracle-flatness argument used in KU1 makes it B-locally free.
This is projectivity over B, NOT invertibility or freeness over S.

Their trace duals are exactly

    T_minus^#=T_plus, T_plus^#=T_minus,
    S subset T_minus subset T_plus subset Lambda.  (4)

At a height-one prime, write a=e-1. The local lattices have exponents
u^-floor(a/2) and u^-ceil(a/2); the trace dual of u^-b S_E is
u^(-(a-b)) S_E. Normality/reflexivity gives (4) globally. The dual lattices
are reflexive too by their B-projectivity and the same surface depth
argument. No crossing or smoothness hypothesis on the whole ramification
curve, and no Galois premise, is used.

In D the submodule A=T_minus/S is isotropic and

    A^perp=T_plus/S,
    A^perp/A=T_plus/T_minus.                       (5)

The induced residual linking form on T_plus/T_minus is perfect, by the
same projective-lattice proof as in Section 1. It is supported EXACTLY on
images of EVEN-index ramification primes. At an odd index the two half
exponents coincide; at an even index e=2m their quotient is
u^-m S_E / u^(-(m-1)) S_E, a single residue-field layer. Its generic length
over the target residue field is the residue degree. As a quotient of
projective B-lattices it has no separate point-supported piece.

In a strictly split tame chart t=u^e, the even-index residual pairing is
e/t times the residue pairing on that layer. Without this split choice,
unit twists and residue-field extensions must be retained. No global
classification of those quadratic residue forms, or their Witt classes,
is inferred from this split formula. The full residual module in (5),
including special points, has been defined globally rather than reconstructed
only from generic lengths.

## 3. Odd indices: perfect trace self-duality coexists with a NONZERO defect

If EVERY ramification index is ODD, the two half lattices coincide:

    T=O_Z(R_pi/2), T^#=T,
    (T*T)^(reflexive S-hull)=Lambda.               (6)

In (5), A=A^perp. Thus beta is METABOLIC, meaning it has this Lagrangian;
its Witt class vanishes by the definition that kills metabolic forms. If
there is any ramification, D is nevertheless NONZERO. This is not a theorem
that an all-odd actual Keller counterexample exists. It proves what the
linking form would do in that case, and the explicit finite-cover controls
below realize that algebraic behavior without the whole-plane hypothesis.

Among reflexive S-STABLE fractional lattices M with S subset M subset Lambda,
trace self-duality is possible IFF every e is odd, and then M=T uniquely.
Indeed at each source DVR self-duality requires 2b=e-1 with b integral;
reflexivity determines the module from these lattices. This classification
does NOT exclude a non-S-stable B-lattice mixing different ramified factors.

For an ACTUAL plane Keller map, both half lattices lie in Lambda subset R:
they consist of source POLYNOMIALS. If all indices are odd and failure held,
T would be a finite B-projective, S-stable, trace-self-dual lattice of
polynomials containing S and 1. Its reflexive S-square would be the NONTRIVIAL
canonical coefficient module Lambda. KN0 makes that canonical class nonzero
and nontorsion. There is no contradiction: B-trace self-duality is NOT
S-duality and does not say that T or the canonical class is torsion over S.
In particular a target-unimodular trace lattice does not give KN1's torsion
canonical hypothesis. No orthogonal-bundle trivialization is imported.

These observations delimit a quadratic/Witt obstruction: vanishing Witt
class, square discriminant in the all-odd case, or existence of T alone
does not identify S with Lambda. It would still be necessary to exclude
the all-odd case using additional ACTUAL source structure, or prove a
stronger multiplication statement. No such exclusion is asserted here.

For the square-discriminant assertion, at each target DVR the inclusion
S subset T changes the Gram determinant by the SQUARE of its index
determinant, while T has unimodular Gram determinant. Thus all discriminant
valuations are even. The polynomial base is factorial and its complex
units are squares, giving the assertion without a quadratic-bundle
trivialization theorem.

## 4. LY1: CUBIC traces recover the missing algebra multiplication

The following criterion does NOT require S-stability of the supplied
intermediate lattice. Let T be ANY finite B-lattice in the actual source R,
containing S, and trace self-dual T^#=T. Then

    Tr(T*T*T) subset B
      iff T*T subset T
      iff T=S=Lambda
      iff the full plane Keller inverse exists.    (7)

The first equivalence is exactly self-duality: ab lies in T iff Tr(abT)
is regular. A finite B-module closed under multiplication is integral over
B by the determinant trick. Thus T lies in the full integral closure S;
the given containment gives T=S. Self-duality then gives S=Lambda and
the retained inverse criterion applies. The converse is immediate.

It is enough to require Tr(t^3) in B for EVERY t in T. Characteristic-zero
polarization recovers the mixed trilinear form. In a supplied local free
basis, regularity of the cubic tensor and invertibility of the quadratic
Gram matrix are precisely regularity of all multiplication structure
constants. A perfect Gram matrix alone supplies no such constants.

In the all-odd actual setup, (6) supplies T canonically. If failure held,
there would be a source polynomial t in T with

    Tr(t) in B, Tr(t^2) in B, but Tr(t^3) NOT in B. (8)

Indeed a failure of the last condition for ALL t would give (7). More
precisely, for ANY finite B-generators of T, a generic CONSTANT complex
linear combination has the bad cubic trace: the coefficients of this
cubic polynomial in the constants lie in a finite-dimensional complex
subspace of F/B. Its zero locus is a proper Zariski-closed subset, since
polarization and (7) show that not all coefficients vanish. The first two
traces are regular because 1 and every t belong to the trace-integral T.

This is an ALL-DEGREE three-moment witness, conditional on an actual
all-odd failure. It does not prove the cubic trace regular or provide a
counterexample. It strengthens the warning about replacing FULL Lambda
by a smaller good quadratic lattice; it does not replace the original
unrestricted trace-pairing goal by the odd-index subcase.

## 5. Exact constant-Gram finite-cover controls for EVERY odd index

Take e=2m+1 with m>=1 and the connected Laurent maps retained in LW:

    P=u^e, Q=v/u^(e-1), Jac_(u,v)=e,
    B=C[P,Q], S=C[u,Q], Lambda=u^(-2m)S,
    source R_Laurent=C[u,u^-1,v].

This is an ACTUAL finite free normal completion, but its open source is
Gm times A1, not the whole polynomial plane. Its half lattice is T=u^-m S.
In the B-basis tau_a=u^a, -m<=a<=m, the trace Gram matrix is the CONSTANT
anti-diagonal matrix

    Tr(tau_a*tau_b)=e if a+b=0, and 0 otherwise.    (9)

Thus it is everywhere perfect, not merely nonsingular generically. Yet
T is not a ring, T properly contains S, and the FULL codifferent pairing
is not regular. The mixed cubic witness is universal:

    Tr(tau_(-m)*tau_(-m)*tau_(-1))=e/P.            (10)

On D=Lambda/S, T/S is a NONZERO Lagrangian. In index three, D has basis
u^-1,u^-2 over the target branch residue field, with linking matrix having
zero diagonal and off-diagonal 3/P. It is a nonzero metabolic form.
For EVERY m the lattice T has constant quadratic trace data as in (9),
while its cubic tensor has the simple pole (10).

These controls are stronger than the retained KR2 rank-one Kummer IC
self-duality control: here a FINITE B-free lattice has a PERFECT ORDINARY
TRACE Gram matrix. No independent historical-priority claim is made for
this standard half-different construction. Source units, the denominator
in Q and the missing retained sheet at P0 are precisely the non-whole-plane
features; no Laurent control is an actual plane Keller or HC4 counterexample.

## 6. Literature, global proof boundary, and the open implication

Square roots of the inverse different with trace self-duality are classical.
The primary metadata/abstracts of Cassou-Nogues--Erez--Taylor (2000),
*Invariants of a quadratic form attached to a tame covering of schemes*,
and Agboola--Burns--Caputo--Kuang, arXiv:1803.09392v6, were checked for
attribution and hypotheses. Their cohomological, normal-basis and arithmetic
class formulas are NOT imported. In particular no normal-crossing, Galois
or invertible-half-different theorem is silently applied to this normal
surface. The present reflexive-module/projective-base proof was written
above; the entire external papers were not read or independently reproduced.

KU's finite duality, KN's canonical ramification vector, LX's finite
discriminant derivative defect, and KR2's distinction from holonomic
self-duality are retained with unchanged assurance limits. New algebra
controls check the half-exponent dualities, even/odd residue layers,
constant-Gram cubic witness and module reconstruction identities. No
registered mathematical checker is replayed, and no primary PDF cached.
No formal verification, independent whole-proof replay or external review
is claimed.

The linking form is not forced to vanish by its being metabolic or by a
smaller perfect trace lattice. In the all-odd case the missing data have
an exact cubic description, but their regularity is UNSUPPLIED for the
ACTUAL whole-plane source. Even-index residual forms and all mixed-index
sources remain part of the full problem; the branch-index cases are not
exhausted by a claimed inverse theorem. Full HC4/JC2 and ordinary
Tr(Lambda Lambda) regularity/full marked-energy integrality remain OPEN.
