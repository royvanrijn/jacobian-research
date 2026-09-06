# An independent principal relation, but no additional strict class yet

Follow-up: the [bounded near-root circuit test](MATCHED103B2_RELATION_CIRCUIT_CAPACITY.md)
adds one independent principal relation but still has zero strict image.
It derives a capacity bound separating unit dependencies from possible
non-unit contributions.

The fixed103b2 pair now has a bounded test of a genuinely independent
class-arithmetic source. Exact binary cubic norm coordinates greatly reduce
coefficient size. The frozen box produces no fully400000-smooth norms,
but one already retained small cofactor can be resolved and its complete
principal ideal verified. It gives a nonzero seven-prime relation in the
localized ideal-class presentation of the high fibre.

It does **not** yield a strict Selmer class. Its odd valuations outside S
cannot be removed by multiplying generic classes. Thus this single element,
even together with all17 generic sections, adds no strict direction.
The additional class basis and quotient CT data remain UNKNOWN.

This extends the [nine-direction necessity](MATCHED103B2_JUMP_REQUIRES_NINE_STRICT_DIRECTIONS.md):
the successful fibre needs at least nine strict rational classes, but an
independent construction must pass a relation-dependency gate before CT
can even be evaluated on them.

## Fixed experiment

The [protocol](MATCHED103B2_NORM_RELATION_PROTOCOL.json) keeps the two masked
fibres, **3726/881** and **−1049/2296**, and tests the same5039 primitive pairs

\[
 -64\le m\le64,\qquad1\le n\le64,\qquad\gcd(m,n)=1.
\]

It compares their already reduced monic cubics with binary cubics derived
from the maximal-order multiplication table. The existing norm-form helper
is imported unchanged; no active class-relation campaign is modified.
Equations, maximal orders and bad primes supply every arithmetic input.
Exceptional points and their class representatives are excluded.

Normalize an integral basis to1,w,z so that

\[
 w^2=-ac+bw-az,\quad wz=-ad,\quad z^2=-bd+dw-cz.
\]

The associated binary cubic is F(u,v)=au³+bu²v+cuv²+dv³, with discriminant
Disc(K). Its positive Hessian supplies an exact determinant-one change
(u,v)=M(m,n). The complete norm identity is

\[
 N_{K/\mathbf Q}(a u+v w)=a^2 F(M(m,n)).
\]

The fixed square factor a² is retained in the subsequent ideal audit.
A square rational factor in the norm cannot simply be discarded from a
principal-ideal calculation.

| Fibre | Arm | Maximum coefficient bits | Tested pairs | Fully smooth at400000 | Smallest retained cofactor bits |
|---|---|---:|---:|---:|---:|
| high | reduced monic | 256 | 5039 | 0 | 174 |
| high | maximal-order binary | 122 | 5039 | 0 | 20 |
| low | reduced monic | 253 | 5039 | 0 | 182 |
| low | maximal-order binary | 175 | 5039 | 0 | 91 |

The [norm certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_matched103b2_norm_relations_v1.json)
binds the forms, transports, complete-enumeration digests and smallest
remainders. Each remainder is obtained by exact repeated gcd with the
primorial through400000. No incomplete remainder is called prime or smooth.

These are different explicitly defined sets of principal elements. The
results measure the computational cost of this class-arithmetic method,
not an incidence or rational-solubility difference between the curves.
In particular, better norm smoothness is not proposed as a rank predictor.
Nor does this experiment measure elliptic point-search visibility.

## Audit of the already retained small cofactor

A separately frozen [audit protocol](MATCHED103B2_RETAINED_NORM_AUDIT_PROTOCOL.json)
permits only saved binary-arm remainders of at most32 bits. There is exactly
one: the high fibre's pair **(m,n)=(57,8)** leaves the prime **735311** after
the400000-smooth part is removed. The low fibre has no qualifying saved
remainder. No larger box or new norm element is evaluated.

The [exact audit](../../artifacts/generated-results/elliptic-curves/rank_jump_matched103b2_retained_norm_audit_v1.json)
constructs β=au+vw, verifies its norm including a², and verifies its entire
principal ideal against the product of prime-ideal powers. After removing
S-supported factors and even exponents, the remaining relation is

\[
 [\mathfrak p_{491}]+[\mathfrak p_{569}]+[\mathfrak p_{22307}]
 +[\mathfrak p_{30557}]+[\mathfrak p_{37831}]
 +[\mathfrak p_{276587}]+[\mathfrak p_{735311}]=0
 \quad\text{in }\mathrm{Cl}(\mathcal O_{K,S_K})/2.
\]

Each displayed prime ideal has residue degree one. The certificate retains
its PARI prime index in the pinned cubic model, so the rational-prime
subscript is not being used as a complete ideal identifier.

This is a verified **relation among ideal classes**, not seven independent
classes or a square-norm Selmer representative. The outside-S parity vector
has seven nonzero coordinates. Every generic point class has even valuations
outside S. Consequently any product β^e times generic classes can be strict
only if e is even; in that case β^e is a square and adds no squareclass.
Since the generic strict kernel is zero, the generated strict image is zero.
This closes the incidence gate for this particular retained element.

## Generation and relations are separate certificates

The cutoff400000 was a smoothness bound, never a generating-set theorem.
For comparison, interval arithmetic gives conservative Bach cutoffs

\[
 B_{\rm high}=843814,\qquad B_{\rm low}=1080103.
\]

Under the relevant GRH assumption, prime ideals of norm below12 log²|Disc(K)|
generate the ordinary class group and hence its S-localized quotient.
This is the classical bound recalled in
[Belabas–Diaz y Diaz–Friedman, introduction](https://www.math.u-bordeaux.fr/~kbelabas/research/OnBach.pdf).
All prime-ideal factors of the audited relation lie within the high cutoff.
The ideal identity itself is unconditional; the generating interpretation
of the cutoff is conditional.

No exhaustive factor base at either cutoff was constructed. One supported
relation does not provide a class-group dimension, a useful upper envelope,
or the desired c_S≤9 certificate. The exact-rank consequence of such an
upper certificate remains hypothetical. An unconditional upper bound also
needs an unconditional generation certificate, not an unlabelled use of GRH.

## Lessons and next gate

1. **Useful independent construction step:** reduced maximal-order binary
   forms produce much smaller norm values on the fixed pair. One retained
   value now gives a fully verified relation without exceptional inputs.
2. **Failed incidence gate:** no strict class is generated by the retained
   relation together with G. Norm smoothness and a principal relation are
   insufficient; dependencies of outside-S ideal parity are required.
3. **Missing upper-bound gate:** a proven generating set and enough supported
   relations. A generating cutoff alone supplies neither a small class rank
   nor its character basis. Any GRH-based result must remain conditional.
4. **Highest-value next test:** a bounded, equation-defined set of relations
   targeting dependencies of outside-S parity on this fixed pair, with local
   squareness and nonsquareness checked on any resulting product. Only a
   certified nonzero strict product reaches the requested additional CT
   computation. No candidate-selection change for Agent1 follows yet.

The two norm boxes and the retained-cofactor audit are complete. The audit
replays all20156 exact values/remainders with independent polynomial
arithmetic, checks the full norm identities and the retained principal
ideal, and validates the outward-rounded conditional cutoffs. A failed
initial audit caused by Sage/Python Fraction interoperability was corrected
before the certificate was produced.

Replay:

```sh
timeout 45 sage -python elliptic-curves/rank-jump/verify_matched103b2_norm_relations.py check
```

The scripts checkpoint under rank-jump-specific local paths. They import the
existing norm-form algebra without editing it and do not read the active
small-conductor campaign's point or relation outputs. No active-search file,
worker limit or mathematical status entry changes.
