# Flat plinth families and one-value recognition

Canonical claims: HC4MYKF0, HC4MYKF1, HC4MYKF2, HC4MYKF3.
[MATH_STATUS.json](../../MATH_STATUS.json) is the authority. Fixed
cylindrical support already makes the entire plinth scheme a finite
flat family over the value line. It does not yet prove that its length
is minimal. One value of that finite scheme now suffices to certify
the full source recovery, provided that minimal length is established.
No native action or source alignment is constructed here.

## 1. Fixed support forces saturation, finite freeness and constant nilpotent length (KF0)

Let X be a smooth irreducible factorial affine complex fourfold, R its
ring, and delta a nonzero LND with actual embedded kernel

    A=B[h], B=C[b1,b2],

the three generators algebraically independent. Put J=pl(delta).
Assume only an equality of radical ideals

    sqrt(J)=I_Z A,
    Z={b_i:1<=i<=n} a nonempty finite reduced subset of A2_b. (1)

The full ideal J is not assumed extended from B, reduced, or locally
a complete intersection. Let D=C[h], K=C(h), T=A/J.

Straighten Z to x=0,C(z)=product(z-i)=0 by the actual JV interpolation.
Some powers x^m,C^m belong to J. They supply the two local slices
used in KE. Factorial kernel closure keeps x,C coprime on the source,
so R=R_x intersect R_C. As in KD/KE, the one Cech torsor class xi on
U times A1_h, U=A2 minus Z, satisfies

    J=ann_A(xi).                                      (2)

In the two-affine cover D(x),D(C), the cohomology is computed by a
complex of C-vector spaces tensored with D. Thus

    H1(U times A1_h,O)=H1(U,O) tensor_C D.              (3)

This is a direct Cech calculation, not a proper base-change assertion.
The right side has no D-torsion. The map

    T -> H1(U,O) tensor_C D,  f mod J -> f xi           (4)

is injective by (2). Consequently

    q(h) f in J, q!=0  implies  f in J,                (5)

and T is torsion free over D. Since A is Noetherian, (1) implies
I_Z^m A subset J for some m. T is a quotient of (B/I_Z^m)[h], so
is a finite D-module. The elementary PID module theorem, explicitly
[Stacks, Lemma 15.22.11](https://stacks.math.columbia.edu/tag/0549),
now makes T finite free, hence flat. Write its rank as L.
The hypotheses of that import are exactly finiteness and torsion
freeness just proved; D=C[h] is a PID.

There is a more precise split sequence of D-modules:

    0 -> M:=I_Z A/J -> T -> (B/I_Z)[h]=D^n -> 0.       (6)

The last module is free, so the sequence splits as D-modules.
M is the nilradical of T and is finite free of rank L-n. On every
closed value h=kappa, the sequence stays exact and M_kappa is still
nilpotent with reduced quotient C^n. Hence

    dim_C T_kappa=L,
    dim_C Nil(T_kappa)=L-n, for every kappa.           (7)

No restriction to regular values is needed in (7). By the Chinese
remainder decomposition for sufficiently thick disjoint points,
T is the product of finite free D-algebras T_i supported on the
individual lines. Each has a constant positive rank L_i, with
sum L_i=L, and nilpotent rank L_i-1 in every value fibre.

In particular nilpotent thickness cannot appear only at a finite set
of critical values. Its multiplication law can still degenerate,
and its socle dimension or number of local generators can jump.
This is a statement about the finite plinth scheme, not about the
plinth of a derivation obtained after specializing the source.

Saturation also gives the exact recovery formula

    J=(J tensor_D K) intersect A.                     (8)

Therefore generic reducedness cannot conceal an exceptional nilpotent
correction on the value line under hypothesis (1).

## 2. The action is Nambu-determined and its support counts generic fibre cycles (KF1)

First take the actual source X=A4, with ordered coordinates X1,...,X4.
Define the polynomial Jacobian derivation

    N(f)=Jac_(X1,...,X4)(b1,b2,h,f).                  (9)

Then for some c in C*,

    delta=c N.                                      (10)

Indeed on the complement of E=V(I_Z R), the quotient is covered by
actual local-slice charts R_g=A_g[t], delta(t)=1, for g in J.
There both delta and N are nonzero: b1,b2,h,t are local coordinates
and the source volume differs from db1 db2 dh dt by a unit. The
set E has source codimension at least two, because it is contained
in V(x,C) and these invariants remain coprime in R. Thus the four
coefficients of each derivation have no common prime divisor.
Both delta and N are primitive on the entire source.

Over Frac(R), the derivations annihilating b1,b2,h have dimension
one, so N=r delta for r in Frac(R)*. Write r=u/v in coprime
polynomials. Polynomiality and primitivity of delta force v a unit;
primitivity of N then forces u a unit. Source units are C*, proving
(10). Local nilpotence of N is a consequence of (10), not an assumption
about arbitrary Jacobian derivations. This is the elementary primitive
comparison from JT, applied to the full embedded quotient. Two chosen
kernel generators and h therefore determine the candidate action up
to a constant once these quotient hypotheses hold.

We next identify the generic whole fibre without assuming J descends.
Localize the source and kernel at D minus 0:

    R_K=R tensor_D K, ker(delta_K)=B_K,
    pl(delta_K)=J_K.                                 (11)

These equalities use localization by invariant elements. They are not
claimed for specialization at a closed h-value. R_K is factorial,
being a localization of R, and has dimension three over K. The
support of J_K consists of the same n K-rational points b_i.

Apply KE0's local negative-Laurent pairing over this field. The
annihilator of xi_i is the corresponding primary ideal of J_K;
its finite multiplication pairing has one-dimensional socle. The
two-variable Koszul calculation gives two local generators F_i,G_i
over B_K. Choose local slices with derivatives G_i,-F_i and numerator
a_i in B_K. The colon identity J_K:a_i=J_K makes its residue at b_i
a nonzero element of K. The local whole source is consequently

    B_(K,W_i)[P,Q]/(F_i P+G_i Q-a_i).                (12)

There are no points in (12) over b_i, since F_i=G_i=0 there and
a_i(b_i)!=0. For completeness, (12) follows by the same two-chart
intersection: its removed set is empty, so the displayed hypersurface
is smooth and normal; it agrees with the normal factorial source
on both charts. No three-variable polynomial-kernel or principal-
plinth theorem has been applied to R_K.

Thus no component of E dominates the value line. There are finitely
many components, so after removing finitely many values to obtain
S subset A1_h, the entire family is the actual Ga torsor

    X_S -> U times S.                                (13)

The topological fibre C is contractible, so for every value in S

    H3(X_s,Z)=H3(U,Z)=Z^n.                            (14)

The point links e_i are a constant integral basis on S. The assertion
is about the whole family, not just a selected open part of a fibre.

There is a useful native consequence even before reducedness of J
is known. Suppose h has N finite Morse points whose vanishing cycles
form an integral basis of the generic H3, as proved in KA for the
remaining frontiers. Then n=N by (14). Every critical point lies in
E, since h is a local coordinate off E. A small Morse vanishing sphere
near a point mapping to b_i projects into a small punctured ball about
b_i, so its global class projects to an integer d_j e_i. Transport
through (13) preserves these labels. The matrix of all N finite
cycles in the e_i basis has just one possibly nonzero entry per column.
It is unimodular by the native basis theorem. Therefore each row has
exactly one column, all d_j=+-1, and the critical points map bijectively
to the n support lines, including collided critical values.

This primitive degree is derived from the full integral cycle basis.
It does not assert that D_b(b1,b2) has rank two at the Morse point,
that a normal-bundle map is transverse, or that a generic Artinian
plinth component has length one. Those are different claims.

## 3. One value of the finite scheme certifies the entire cotangent inverse (KF2)

Under (1), the following are equivalent for this same action:

- J=I_Z A;
- L=n;
- J_K is reduced;
- for some kappa in C, T_kappa is reduced;
- for some kappa in C, dim_C T_kappa=n.

They also imply that every value fibre T_kappa is reduced. To prove
the equivalences, use (6)--(7): the nilradical is finite free of
rank L-n, on the whole line, at the generic point and at every closed
point. It is zero in any one of these places exactly when L=n.
This is stronger than constancy of total length alone and requires
the fixed support in (1). Equivalently it is enough, for each i,
to find a value kappa_i with dim_C (T_i)_(kappa_i)=1; the values for
different components may be different.

On actual A4, any of these conditions invokes KD to recover the
full simple cotangent frame and polynomial inverse. No separate
source codimension, quotient-image, slice or fixed-plane certificate
is needed. The new test uses the fibre of the actual ideal J:

    T_kappa=B/J_kappa,
    J_kappa=(J+(h-kappa))/(h-kappa).                  (15)

It does not substitute pl(delta specialized at kappa), or the
annihilator of a specialized torsor class, for J_kappa.

For the native frontiers, (14) makes n the known critical count N.
Consequently pure source alignment with the explicit R_(N,c) is
equivalent to existence of one action with the actual embedded
kernel B[h], fixed finite cylindrical radical support (1), and

    dim_C A/(J,h-kappa)=N for one kappa.              (16)

For S_rho the number is three, for K2/h1 it is two, including all
critical-value collisions. Necessity comes from transporting the
explicit normal-form translation. Sufficiency is the preceding
minimal-length argument, KD and the native common critical Hessian
condition in KB. Under these native hypotheses, exact finite-ideal
descent from KE is equivalent to the same action satisfying (16):
KE forces a descended ideal reduced, while (16) forces J=I_Z A.
Unit-Hessian source existence still equals E_N AND this pure action
criterion. The global plane Keller gate has not been absorbed.

The value kappa in (16) may itself be critical, since (7) holds at
every value. However, no such value fibre or native action has been
constructed here. Fixed support alone is not proved to imply L=N.
The remaining fixed-support problem is now precisely whether the
native geometry forces minimal length, or permits a family with
positive nilpotent rank at every value. Both endpoint inversions,
E2/E3, source alignment and unrestricted HC4/JC2 remain open.

## 4. Flatness does not make annihilator specialization exact

The following is an explicit local-cohomology control. Discard all
Laurent monomials with a nonnegative u or v exponent, and put

    zeta=1/(u^2 v)+h/(u v^3).

Its exact annihilator over C[h,u,v] is

    J_zeta=(u^2,uv,v^2-hu).                           (17)

The quotient is finite free over C[h] with basis 1,u,v and rules
u^2=uv=0, v^2=hu. The three images of this basis multiplied by zeta
are independent over C[h], which proves (17). Its radical is the
constant cylinder (u,v). Generically it is the curvilinear algebra
K[v]/(v^3), with u=v^2/h, so is a complete intersection with socle
dimension one. At h=0 the actual ideal fibre is

    (J_zeta)_0=(u,v)^2,

of length three with socle dimension two. But zeta_0=1/(u^2v) has
annihilator (u^2,v), of length two. Thus

    (ann zeta)_0 != ann(zeta_0).                      (18)

All total and nilpotent fibre lengths of the finite scheme (17)
stay constant, as KF0 requires. Its Gorenstein property and local
number of generators jump, and the specialized inverse system loses
length. Generic residue duality alone cannot justify applying the
two-generator argument at every closed value when the ideal has not
descended from B. This is a local Cech example, not a claimed native
four-space action. It neither disproves nor proves automatic native
minimality in (16).

An even simpler false-positive control is

    theta=h/(u^2v)+(h-1)/(uv).

Its full annihilator is (u^2,v), of constant length two. At h=0,
theta_0=-1/(uv) has reduced annihilator (u,v), of length one, whereas
the actual ideal fibre is still (u^2,v). Thus substituting the
annihilator after specializing the class could falsely pass the
minimal-length test. Even a residue coefficient linear in h does
not justify that substitution. This remains a local-cohomology
control with no native realization asserted.

## 5. A line-bundle lift cannot directly principalize this plinth (KF3)

A focused comparison was made with the research manuscript
[cancellation construction, Section 4](https://github.com/openai/math/blob/main/preprints/An-explicit-failure-of-complex-affine-space-cancellation-September-23-2026/build/sections/04-bundle.tex)
and its [rigidity section](https://github.com/openai/math/blob/main/preprints/An-explicit-failure-of-complex-affine-space-cancellation-September-23-2026/build/sections/06-rigidity.tex).
The actual coefficient inclusion there is

    R0=C[x,y,z,u]/(xy-z(z+1)),
    s_R=y-x(x-u^3), I0=(x-s_R^2-u^3,z-s_R x),
    R0 subset G=R0[tau,I0/tau^2] subset R0[tau,tau^-1].

Its line-bundle torsor has coordinate ring

    Rtilde=C[a,d,b,c,u]/(ac-bd-1),
    x=ab,y=dc,z=db,
    v=a^3b-a^2u^3-d^2.

The paper supplies actual Laurent charts and identities
I0 Rtilde=(v), with v nonzero and nonunit. It tracks fibre weights
(a,d,b,c,u)=(1,1,-1,-1,0) and valuation degrees tau=-1,V=2 in
Gtilde=Rtilde[tau,V]/(tau^2 V-v). The additive lift is over Spec G;
it need not preserve R0 or Rtilde. Its separate degree action fixes
those coefficient rings. The rigidity step also needs a specific
leading derivation E of nonpositive fibre shift with E^2(v)=0.
None of these ring inclusions or leading-degree conditions has
been constructed for the native S_rho action gate.

There is a direct structural mismatch even before those latter steps.
A line-bundle torsor is locally a Laurent polynomial extension and
is faithfully flat and surjective, with one-dimensional fibres.
The inverse image of the nonempty codimension-two V(J) in Spec A
therefore still has codimension two. A nonzero nonunit principal
ideal in this smooth integral total space has a codimension-one
zero set. Thus no line-bundle torsor over Spec A can make J into
a nonunit principal ideal. The same argument applies to a nonempty
source pullback of codimension at least two. Equivalently, invertibility
of an ideal descends faithfully flat locally; a non-Cartier point-
cylinder ideal cannot become Cartier through this pullback.

The manuscript's principalization is consistent with this fact: its
I0 pulls back to the divisor v=0, so it is in a codimension-one
Cartier situation, rather than the present codimension-two plinth.
On the actual polynomial kernel A3 and source A4, line bundles are
already trivial (invertible ideals in these UFDs are principal),
so their complements merely adjoin a Laurent variable. A lift of
an additive action on a trivial line bundle, normalized at the
identity and commuting with fibre scaling, fixes that fibre variable:
its multiplier is a unit in R[t], hence is independent of the action
parameter and equals one at zero. The lift retains the original
action; it supplies no new principalization or contradiction.

A blow-up can principalize a codimension-two ideal, but it is a
different, nonflat construction and changes the source/quotient.
Alternatively an established filtration might produce a different
invertible ideal in a different coefficient ring. Either could be
useful only after proving the new ring inclusion, action lift,
leading derivation, weights and descent. The present mismatch does
not rule out every use of filtrations or line bundles.

The companion [stable-coordinate manuscript introduction](https://github.com/openai/math/blob/main/preprints/A-stable-coordinate-that-is-not-a-coordinate-in-four-variables-October-5-2026/build/source/sections/01-introduction.tex)
reports an example separating affine-space fibres and stabilization
from ambient rectifiability. Only that introduction was inspected
here; its main theorem is not imported as a proved counterexample.
Our recognition uses exact ideal and coordinate-ring generation,
so it does not infer ambient coordinates from fibre geometry alone.
No manuscript lifting or rigidity theorem is imported here. The
local codimension argument is proved above, and the displayed
principalization/weight identities are checked as source controls.

## Sources and verification

The sole new proved algebra import is the finite torsion-free/PID
conclusion of [Stacks 15.22.11](https://stacks.math.columbia.edu/tag/0549),
read with its proof. All source/intersection, residue duality, generic
local-presentation, primitive-comparison and native-cycle inputs
are retained from unchanged KD, KE, JT, KA and KB evidence.
The three manuscript TeX sources are retained with their Git blob
identifiers as a focused applicability audit, separately from that
proved import. Their cancellation and stable-coordinate claims are
not used as premises for any HC4 conclusion.

[The checker](verify_flat_plinth_families.py) checks the torsion-free
Laurent module model, full annihilator and specialization discrepancy,
generic/special socles, primitive Jacobian derivations, integral cycle
matrix consequences and one-value gate logic. The all-degree and
all-value proofs above do not rest on a bounded computation.
[The context](sources/flat-plinth-families-context.json) pins prior
evidence and the exact new import; [the receipt](receipt-flat-plinth-families.json)
records the controls. Independent replay, formal verification and
external review are not claimed. The original goal is not complete.
