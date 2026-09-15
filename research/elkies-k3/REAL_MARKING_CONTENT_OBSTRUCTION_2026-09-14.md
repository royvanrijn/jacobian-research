# A real-place obstruction to a full rational K3 marking

Let X be a projective K3 surface over a subfield of R, with every geometric
Néron–Severi class Galois invariant. Let T be its **actual integral**
transcendental lattice and let m be the gcd of the entries of any Gram
matrix of T. Then **m divides 2**.

This is a necessary condition, not a sufficiency theorem. In particular,
content one or two supplies neither a rational marking nor an MW17
fibration. The argument does not use a norm-one subgroup, a choice of
Shimura descent, or a rational-point census.

## Proof

Write L=H²(X(C),Z), with its unimodular intersection pairing, S=NS(X(C)),
and T=S-perp. These are the usual primitive orthogonal complements. Let c
be the pullback of complex conjugation on integral Betti cohomology, and
put s=-c. The minus sign is essential: the cycle-class map is Galois
equivariant in H²(1), so c acts as -1 on classes of invariant divisors.
Thus s fixes S pointwise and preserves T. It is an integral involutive
isometry.

Unimodularity identifies A_S and A_T by the gluing anti-isometry. Because
s extends to L and fixes S, its action on A_T=T*/T is the identity.
Equivalently, (s-I)T* is contained in T. As all pairings in T are divisible
by m, (1/m)T is contained in T*. Consequently (s-I)T is contained in mT.

The involution s on T is not the identity. Complex conjugation exchanges
H^(2,0) and H^(0,2); on their real period plane its eigenvalues are +1
and -1. Multiplication by -1 leaves one eigenvalue of each sign. This
plane lies in T tensor R. Since s is rational, its nonzero -1 eigenspace
contains a primitive integral vector v in T. Now (s-I)v=-2v belongs to
mT. Primitivity of v implies m divides 2, as required.

This uses only the standard integral cohomology, cycle-class, Hodge, and
primitive-gluing facts for K3 surfaces. The involution argument is written
out; no torsion-freeness theorem or modular-curve identification is imported.
It applies in every Picard rank, provided T is the actual complement and
the **full** NS is invariant. It does not exclude a surface having an
unsaturated displayed subgroup with a larger-content apparent complement,
or realizations over fields with no real embedding.

## Retained-catalogue application

The [audit](../artifacts/generated-results/elkies-k3-real-marking-content-v1.json)
reads the literal Gram matrices of the retained 827-row T catalogue.
Seventy have content greater than two. Three already have recorded
arithmetic exclusions; **67 additional rows** in the retained 820-row
queue are excluded by this real-place obstruction. The determinant-948
positive control passes the necessary condition.

In particular, the two formerly leading rows are

| Row | Literal T | Content | Outcome |
|---|---|---:|---|
| K3-be060aee9b10a819 | U(4) + <16> | 4 | no full real marking |
| K3-0e522728bed79087 | U(4) + <32> | 4 | no full real marking |

Removing these 67 rows leaves 753 unresolved rows in that historical
queue. Of its 21 coarse-genus-at-most-two diagnostics, only determinants
800, 736, and 480 pass this gate. These are priorities for further
arithmetic analysis, **not positive handoffs**. The retained earlier planner
and its certificate remain unchanged; this audit supplies a sourced filter
on that snapshot, not a reconstruction of its lattice searches.

The input snapshot contains only its original determinant-948 positive
control. The separately recovered determinant-1092 parent must also remain
part of the comparison baseline; this audit does not establish independence
from that mechanism or exhaust all rank-19 K3 lattices.

## Reproduction and boundary

From the repository root:

```sh
python3 research/elkies-k3/scripts/audit_real_marking_content.py --check
```

The checker recomputes exact integer gcds and determinants, reconciles
surface identifiers, and checks the necessary-condition filter against the
retained positive and negative controls. It does not prove the Hodge/gluing
argument by machine or certify any new positive surface. No independent
implementation, formal verification, or external review is claimed.
