# Erdős #942 and coprime occupants of square intervals

Write \(h(q)=\#\{t:q^2\le t<(q+1)^2,\ t\text{ powerful}\}\), as in Erdős #942. The problem asks for its size, in particular matching logarithmic-power upper and infinitely-often lower behavior. The left endpoint \(q^2\) is always powerful, so the count in the open interval is exactly \(h(q)-1\).

The [short-power-window theorem, CPW1](ABC_COPRIME_SHORT_POWER_WINDOWS.md#a-positive-density-primitive-construction), specialized to \(k=2\) and centre \(q+1\), already proves: for each **fixed** \(m\), a positive proportion of sufficiently large \(q\) have at least \(m\) pairwise coprime, non-square powerful integers strictly in \((q^2,(q+1)^2)\). Its proof uses fixed primes and a fixed-dimensional coprimality sieve; the density constant and threshold may depend on \(m\).

[Narumi–Tachiya, Theorem 1 and Lemmas 1–2](https://arxiv.org/html/2611.20656v1#S2) give growing ordinary occupancy in successive power intervals, without pairwise coprimality. The [repository's existing comparison](ABC_COPRIME_SHORT_POWER_WINDOWS.md) records this distinction. An exact novelty comparison with other coprime constructions is still needed before presenting CPW1 as a new #942 contribution.

The fixed-\(m\) theorem already implies **qualitative unbounded coprime occupancy**: choose recursively an increasing centre \(q_m\) from the positive-density set for each \(m\). Thus \(q_m\to\infty\) and the number of pairwise coprime non-square occupants at \(q_m\) is at least \(m\).

The open refinement is a **quantitative** growing multiplicity, for example at least \((\log q)^\eta\) for some fixed \(\eta>0\) along an unbounded sequence. The fixed-\(m\) theorem gives no rate: its constants are not uniform as \(m\) grows. A proof needs quantitative simultaneous approximation and a coprimality sieve uniform in the number and size of chosen classes. Even this refinement would leave Erdős's logarithmic-power upper bound and matching lower exponent open.
