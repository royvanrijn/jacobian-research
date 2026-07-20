# Provenance audit

Audited 20 July 2026. This note separates the earliest public announcement
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

## Later exposition, not discovery provenance

Gallagher's page explicitly distinguishes the externally announced map from
its own weighted-lift interpretation and higher-degree seed family. It names
its `RESEARCH.md` and verification scripts as the source of truth for those
follow-on claims, and says it does not establish whether the discoverer's
derivation used that mechanism.

Zhang's page is a same-day consequence audit, not the discovery record. Its
value for provenance is its direct link to the announcement and exact formula
transcription.

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
