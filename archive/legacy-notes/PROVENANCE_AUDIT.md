# Provenance audit

Audited 20--21 July 2026. This note separates the earliest public announcement
located so far from exact mathematical verification and from later structural
exposition.

## Public trail located

The earliest primary public item located is Levent Alpöge's X post, status
[`2079028340955197566`](https://x.com/__alpoge__/status/2079028340955197566).
The X page was not readable through the audit environment, so its text and
timestamp could not be captured directly.

Two same-day sources link that exact post and independently transcribe the same
map:

- [Zihan Zhang's source audit](https://zzhang-iu.github.io/papers/direct-consequences-jacobian/index.html)
  dates the post to 20 July 2026 UTC and says it credits Akhil with posing the
  question and Fable with producing the example.
- [Alexis Gallagher's explainer](https://jacobianfun.org/jacobian-explained)
  describes Alpöge's announcement as 19 July 2026 and credits “Levent Alpöge +
  Fable.” The one-day discrepancy is consistent with timezone presentation but
  cannot be resolved without a direct timestamp or archive of the post.

Both sources give, in the same coordinate order,

\[
\begin{aligned}
F_1&=(1+xy)^3z+y^2(1+xy)(4+3xy),\\
F_2&=y+3x(1+xy)^2z+3xy^2(4+3xy),\\
F_3&=2x-3x^2y-x^3z,
\end{aligned}
\]

the determinant `-2`, and the same three-point collision stored here. This
agrees with the repository's independent finite certificates.

## Lean formalization

[Dean Cureton's `deancureton/jacobian`](https://github.com/deancureton/jacobian)
provides a separately authored Lean 4 formalization.  Commit
[`0d4a9212d874226ad81ce5a926becddfa94e6a88`](https://github.com/deancureton/jacobian/commit/0d4a9212d874226ad81ce5a926becddfa94e6a88),
authored 20 July 2026, proves the determinant and three rational evaluations
for the displayed map, then proves a determinant-one rescaling noninjective
over every field and specializes the result to `C`.  Its README credits
“Levent Alpöge/Fable 5's counterexample”; this repository credits Cureton for
the formalization and does not infer discovery priority from the later proof.

The audited commit has no license file.  The integration therefore records
authorship and an immutable source link and supplies a target that builds the
upstream checkout, rather than copying or adapting its Lean source.  See
[Lean foundational-map audit](../../verified/LEAN_FOUNDATIONAL_MAP.md) for the theorem-by-theorem scope and reproduction
instructions.

## Later exposition, not discovery provenance

Gallagher's page explicitly distinguishes the externally announced map from
its own weighted-lift interpretation and higher-degree seed family. It names
its `RESEARCH.md` and verification scripts as the source of truth for those
follow-on claims, and says it does not establish whether the discoverer's
derivation used that mechanism.

Zhang's page is a same-day consequence audit, not the discovery record. Its
value for provenance is its direct link to the announcement and exact formula
transcription.

## Marked-projective-root interpretation

Qiaochu Yuan's same-day MathOverflow update links a geometric interpretation
to [Andy Jiang (`@davikrehalt`), status
`2079175065695035442`](https://x.com/davikrehalt/status/2079175065695035442).
The status ID encodes a posting time of 20 July 2026 at 12:02:19 UTC. The link
and attribution were located for this audit at 20 July 2026 22:43 UTC
(21 July 00:43 CEST).

The repository adopts from that post the organization of the source as binary
cubics with a marked simple projective root. The affine cubic, its derivative
reconstruction, generic degree, discriminant, root at infinity, and `3/1/0`
fiber results substantially overlap the pre-existing cubic marked-root and image results material in
`IMAGE_AND_NONPROPERNESS.md` and its exact scripts. The new
`MARKED_ROOT_MODEL.md` supplies the two-chart global isomorphism proof and
credits the public post for the projective marked-root formulation. This
audit has not established priority for that formulation, and the repository
makes no such claim.

## Juntang Zhuang's quartic compilation

On 21 July 2026, Juntang Zhuang published
[`jzkay12/jacobian_conjecture`](https://github.com/jzkay12/jacobian_conjecture)
with a PDF titled *Explicit Polynomial Maps with Constant Jacobian and
Verified Collisions* and parallel Python and Wolfram Language checkers.  This
audit pins commit
[`1ff68e870f66afec8c6611f910fcc8f5522fdbce`](https://github.com/jzkay12/jacobian_conjecture/commit/1ff68e870f66afec8c6611f910fcc8f5522fdbce).
The source names three quartic examples `F4a`, `F4b`, and `F4c`, labels them
Islands A, B, and C, displays their expanded rational-coefficient maps, and
checks their constant Jacobians and listed rational collisions.  The PDF is
authored by Juntang Zhuang and says it was inspired by Levent Alpöge's public
post.

This repository credits Zhuang for that public compilation, notation, and
the external finite certificates.  No license file was present at the pinned
commit, so no upstream source code or prose is copied here.  This audit does
not infer discovery priority for
the individual quartic formulas.  The local integration does not copy the
upstream checker: it reconstructs the maps from compact formulas and derives
their weighted-seed resolvents and canonical boundaries independently.  The
result is the [quartic-islands audit](../../extended-geometry/EXTERNAL_QUARTIC_ISLANDS.md):
all three have degree four and `S_4` monodromy; Island A is the canonical
triple-zero seed, while B and C are split seeds with extra roots `3` and
`-1/2`; none is a cancellation `(m,r,h)` normal form.

## Not yet located

- the original Fable conversation or full prompt;
- discovery/search code, model and sampling settings, random seeds, or
  intermediate candidates;
- a repository or immutable archive owned by the discoverer;
- a paper or preprint giving the construction and attribution history; or
- a direct archived copy of the announcement with an unambiguous UTC timestamp.

Accordingly, the formula and public attribution trail are corroborated, but
the construction provenance is still incomplete.  Mathematical validity rests
on the executable core certificate together with the written proofs of the
subsequent structural theorems, not on attribution or source authority.
