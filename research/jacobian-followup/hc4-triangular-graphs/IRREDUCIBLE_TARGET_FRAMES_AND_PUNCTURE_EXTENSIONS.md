# Actual irreducible target frames and the puncture extension

2026-11-16. HC4MYLR0/1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
External review/outreach is DEFERRED. LQ0/1 made PROGRESS by using the
actual Keller mate in discriminant and canonical-boundary constraints.
Here a published reduction supplies its all-irreducible frame, and the
remaining cohomological extension is made explicit in every degree.

## 1. LR0: the all-irreducible frame is actually supplied

We explicitly adopt Kaliman, *On the Jacobian Conjecture*, Proc. Amer.
Math. Soc.117(1993),45-51, DOI10.1090/S0002-9939-1993-1106179-7:
for EVERY plane Keller map F there is a polynomial TARGET automorphism
a such that P=first(a composed with F) has EVERY fibre irreducible.
Its other component Q remains an actual polynomial Keller mate.

The complete primary proof, Lemmas1-8 and final argument, was read in
the [public primary-paper mirror](https://www.scribd.com/document/1040189362/On-the-Jacobian-Conjecture-1993).
It combines the bad-line locus with a target shear and a generic final
linear combination. This is an actual construction on the given map,
not only an implication between conjectures. The publisher PDF returned
403; the cached mirror is OCR with transcription defects. We import
the published statement, not an exact OCR shear formula. No original
PDF visual inspection, independent reproduction or novelty is claimed.

Write R=C[x,y], B=C[F1,F2], L=Frac R, K=Frac B. A target automorphism
gives EXACTLY

    C[P,Q]=B,   C(P,Q)=K.                           (1)

Consequently the embedded full finite normalization S of B in L and
its trace dual Lambda are IDENTICAL before and after this change.
The ordinary trace-pairing gate is literally unchanged. The source
plane is unchanged as well. Scale Q by a nonzero constant if needed
to normalize its Jacobian with P to one; this also preserves B.

The generic genus, punctures and fibre monodromy of the ORIGINAL first
component need NOT be preserved. An original puncture hypothesis cannot
be transferred through this reduction. Apply LQ1 to the NEW P instead:
its geometric compact generic curve admits NO nonzero REGULAR positive
pluricanonical section whose zeros are entirely at its punctures, under
failure. In particular the specified hyperelliptic Weierstrass-point
or complete hyperelliptic-pair puncture configurations are excluded in
this SUPPLIED frame. Effectivity is still not automatic.

## 2. Fixed cohomology and residues: an elementary source-plane proof

First allow ANY primitive polynomial submersion P on the whole plane.
Choose a nonempty finite set T of b removed values containing all
atypical values and any further values needed for a smooth projective
compactification with finite etale puncture scheme. Put U=C minus T.
The family X_U=P^{-1}(U) is a locally trivial smooth connected-fibre
bundle. Write C for a fibre, Cbar its compactification, S_infty its
punctures, s their number, and h their number of monodromy orbits.
Use COMPLEX cohomology throughout. Let G=pi1(U).

Let r_c be the number of irreducible components of P=c and
delta=sum_c(r_c-1). The removed divisor D=P^{-1}(T) is a disjoint union
of N=b+delta connected smooth curves. The Thom-Gysin sequence for
D inside C^2, using H^1(C^2)=H^2(C^2)=0, gives

    H^1(X_U) = C^N.                                (2)

The meridians, or equivalently logarithmic differentials of the
irreducible factors of the removed fibres, give a basis. The fibration
low-degree cohomology sequence is

    0 -> H^1(U) -> H^1(X_U) -> H^1(C)^G
      -> H^2(U)=0.                                 (3)

Thus dim H^1(C)^G=N-b=delta. The b pullbacks dlog(P-c)
are exactly the one product relation for each removed fibre. This
recovers, with classical credit, Friedland1999 Theorem2.7; its full
statement/proof pp229-230 and Corollary5.2 p242 were also reread.
Equations(2)-(3) supply the proof used here independently of that import.

Set V=H^1(Cbar), H=H^1(C), W=ker(C^{S_infty}->C), where the last map
sums coordinates. Residue gives the equivariant exact sequence

    0 -> V -> H -> W -> 0.                         (4)

The Tate twist has trivial topological action here. W is a finite
permutation quotient, and dim W^G=h-1. The fixed classes in H are the
logarithms of generic-fibre units from the removed reducible factors.
Their residues are their integer boundary divisors. By LQ1 the unit
lattice modulo base constants has rank delta. Its divisor map into
Div(S_infty) is injective: a function with no zero or pole on the
compact curve is constant. Integer divisor rank is unchanged over C.
Hence H^G -> W^G is INJECTIVE, and V^G=0.

The connecting map of group cohomology therefore has EXACT image rank

    rank( W^G -> H^1(G,V) ) = h-1-delta.            (5)

In particular delta<=h-1. Since G is free of rank b and V^G=0,
its standard length-one cochain complex also gives

    dim H^1(G,V)=(b-1)*2g,
    0 <= h-1-delta <= (b-1)*2g.                    (6)

These statements concern affine FIBRE cohomology, not the finite
monodromy of the plane map. They hold after a further shrink of U;
extra removed typical values add one component and one base generator,
so do not change delta or the fixed-cohomology calculation.

## 3. LR1: the actual mate forces two kinds of puncture

Now P,Q is an actual plane Keller pair. On C, Q is nonconstant and
dQ is nowhere zero. Extend Q to a meromorphic map Cbar -> P^1 of
degree d_Q>=1. It has a nonempty pole set, contained in S_infty.

If EVERY puncture were a pole of Q, its differential would have no
zeros anywhere and its divisor degree would be

    2g-2 = -d_Q-s.

Thus d_Q+s=2-2g, forcing g=0,d_Q=s=1. Then C is A^1 and Q has
degree one on it. Such an actual typical closed fibre is A^1; retained
GR0 therefore makes the whole Keller pair invertible. Consequently
under failure at least one puncture is a pole and at least one has a
finite Q value. After shrinking U the meromorphic map extends over the
projective family and pole degrees are constant. These two nonempty
sets are invariant under monodromy, so h>=2.

In the ACTUAL Kaliman-supplied frame, delta=0. Equations(4)-(5) give

    H^G=0,  V^G=0,  dim W^G=h-1>=1,
    W^G -> H^1(G,V) injective.                     (7)

The extension(4) is NONSPLIT. H is NOT a semisimple G representation;
its image is INFINITE, because finite complex group representations
are semisimple. Also g cannot be zero: then V=0 and H=W, contradicting
(7). LQ1 excludes genus one in this all-irreducible frame, so g>=2.
This genus restriction is CLASSICAL, already Friedland Corollary5.2;
no new genus bound or exclusion of all higher genera is asserted.

By retained LG1 the point-supported torsion of the raw module
M=R/(P_y partial_x-P_x partial_y)R has total multiplicity delta.
Thus M is TORSION-FREE in this supplied frame. Under failure its
generic rank is 2g+s-1>=5. Its connection is bijective and its marked
class [1] vanishes because of the actual mate (LP0). These facts do
NOT force M=0: that remaining pointed implication is still JC2.

## 4. Why compact semisimplicity does not close the gap

Even finite compact monodromy and finite puncture monodromy permit
a nonsplit affine extension. For a free group on two generators take
V=C^2, W=C, with matrices on H=V plus W

    A = [ -I_2  e_1 ],     B = [ -I_2  0 ].
        [   0    1 ]          [   0   1 ]

Both diagonal representations have finite image and V is symplectic.
But H^G=0 whereas W^G=C, and AB is a nontrivial unipotent matrix of
infinite order. This is an ABSTRACT representation control, not a
geometric Keller counterexample. It shows precisely why a theorem
about the compact quotient alone cannot split(4).

No affine splitting, canonical-boundary effectivity, finite invariant
source lattice or ordinary trace-pairing landing is proved automatic.
A proof forcing a splitting in SOME actual supplied all-irreducible
frame would finish, but adding it as a hypothesis does not advance the
last implication. The next useful question is whether the commuting
polynomial Keller lifts constrain the NONZERO connecting classes in
(7), or directly force the unchanged finite trace-dual pairing to land.

## 5. Evidence and boundaries

The new checker only tests new shear/target-ring identities, local
differential divisor bookkeeping and the explicit nonsplit block
representation, plus retained byte pins and scope flags. The written
all-degree proof, not those finite controls, establishes(1)-(7).
No registered mathematical checker is replayed. The cached Kaliman
mirror is a published proved import with the reading limitations
stated above. All assurance flags remain false. Full HC4/JC2 is OPEN;
full quartic plane realization queues remain retired by the retained
published result, and external review remains DEFERRED.
