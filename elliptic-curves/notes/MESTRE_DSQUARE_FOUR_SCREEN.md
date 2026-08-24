# Four split-infinity Mestre families and the `u=197` rank-17 fiber

Status: exact family identities, bounded search, exact conductor computation,
and an unconditional Mordell--Weil rank lower bound.  The search is not an
upper bound and does not meet the rank-21 or rank-30 target.

## The four families

For each ordered root tuple, let `Q_T(X)` be the primitive quartic remainder
in the square approximation to

\[
\prod_i(X-T-r_i)\prod_i(X+T-r_i).
\]

Write its coefficients in ascending order as

\[
Q_T(X)=e(T)+d(T)X+c(T)X^2+b(T)X^3+a(T)X^4.
\]

The exact reconstructed coefficients are as follows.

For `(0,7,225,232,235,265)`:

\[
\begin{aligned}
e={}&9T^6-811682T^4+14880036425T^2+159815739330625,\\
d={}&2(2892T^4-86747347T^2-2590866182800),\\
c={}&-18T^4+1338504T^2+49472753311,\\
b={}&-6(964T^2+31474383),\\
a={}&3(3T^2+84878)=9(T^2+84878/3).
\end{aligned}
\]

For `(0,9,213,247,256,291)`:

\[
\begin{aligned}
e={}&9T^6-873878T^4+22064800401T^2+355052317780224,\\
d={}&4(1524T^4-66633679T^2-2753992953792),\\
c={}&-18T^4+1686984T^2+108001331119,\\
b={}&-12(508T^2+35471675),\\
a={}&3(3T^2+196250)=9(T^2+196250/3).
\end{aligned}
\]

For `(0,25,95,143,168,205)`:

\[
\begin{aligned}
e={}&T^6-36418T^4+371623025T^2+3140781450625,\\
d={}&2(212T^4-4178297T^2-76432519800),\\
c={}&-2T^4+82536T^2+2523380759,\\
b={}&-2(212T^2+8415079),\\
a={}&T^2+39146.
\end{aligned}
\]

For `(0,43,128,197,231,289)`:

\[
\begin{aligned}
e={}&9T^6-650822T^4+10500535681T^2+160595196283456,\\
d={}&12(444T^4-13952701T^2-455859342752),\\
c={}&-3(6T^4-443672T^2-21219509653),\\
b={}&-36(148T^2+8394945),\\
a={}&9(T^2+55950).
\end{aligned}
\]

The fixed square contents before primitive normalization are respectively
`1800^2`, `1872^2`, `6600^2`, and `5712^2`.  In every case the Mestre
obstruction is zero and the primitive discriminant is even of degree 20 in
`T`.  The leading coefficients split after the base changes

\[
T=\frac{C-u^2}{2u},\qquad
\sqrt{a(T)}=m\frac{C+u^2}{2u},
\]

with `(C,m)` equal to `(84878/3,3)`, `(196250/3,3)`, `(39146,1)`, and
`(55950,3)`.  At `u=1`, the twelve visible affine sections and split infinity
have exact mod-3 dimension 12 in each family.  For the third family a separate
affine-section elimination finds six further generic companion sections.  At
`u=197`, one companion raises the visible-plus-infinity dimension from 12 to
13, proving generic rank at least 13.  That theorem is owned by
[the affine-section moduli note](MESTRE_AFFINE_SECTION_MODULI.md) and
[`verify_mestre_02595143168205_rank13_section.py`](../cas/verify_mestre_02595143168205_rank13_section.py),
not by the bounded specialization search below.

## Frozen bounded screen

The screened population was all positive reduced `u=a/b` with `a<=512` and
`b<=32`: 10,170 parameters per family, or 40,680 family-parameter pairs.  The
C++ scorer used every usable prime from 11 through 251.  For each family the
top 48 local-score rows and 16 smallest-`T`-height anchors were selected;
overlap left 255 conductor calls.

PARI/GP completed 153 conductor computations and hit the explicit 30-second
cap on 102.  Exactly 61 completed fibers had `ln(N)<182.72`; all 61 received
the three frozen ratpoints charts `x=z`, `x=T+z`, and `x=-T+z` with
`H=2,000,000` and denominator at most 13,000.  There were no point-search
timeouts.  The exact mod-3 lower-bound distribution was

| lower bound | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| fibers | 2 | 2 | 12 | 12 | 14 | 14 | 4 | 1 |

The 102 conductor timeouts remain unclassified; they are not negative
conductor results.  Likewise, no absence of points outside the ratpoints box
is asserted.

Within the 61 exactly conductor-qualified and point-searched fibers, the
global rank/conductor Pareto frontier is:

| family | `u` | `T` | exact rank lower bound | exact conductor | PARI `ln(N)` |
| ---: | ---: | ---: | ---: | --- | ---: |
| 2 | 74 | 455/2 | 12 | 5989694742451216979801855127121230 | 77.7753485184... |
| 2 | 148 | 233/4 | 13 | 2195708281586324417451772556822101309014 | 90.5873233008... |
| 0 | 186 | -305/18 | 15 | 107596381215984028381677371798362506621890021610 | 108.2947162001... |
| 0 | 333/2 | 185/108 | 16 | 6702328565434011423906254252272908920085073667739673087338570 | 140.0575605927... |
| 2 | 197 | 337/394 | 17 | 2462086522751621334987931952469307556796057284118717977320345864383117775914 | 173.5948911450... |

The familywise maximum lower bounds are 16, 16, 17, and 16, attained at
`u=333/2`, `350`, `197`, and `234`, respectively.  Their corresponding
`ln(N)` values are 140.0575605927..., 153.3577925264...,
173.5948911450..., and 160.9683340795....

## Exact `u=197` certificate

For roots `(0,25,95,143,168,205)` and `u=197`, the base parameter is

\[
T=\frac{39146-197^2}{2\cdot197}=\frac{337}{394}.
\]

PARI/GP independently reconstructs the global minimal model

```text
[1,1,1,
 -1163348683373499147707371416562962,
 15227131493689013260364706485730874765958430844575]
```

with root number `-1` and exact conductor

```text
2462086522751621334987931952469307556796057284118717977320345864383117775914.
```

The discovery pool has 28 distinct Jacobian representatives modulo inversion:
13 displayed points and 15 novel abscissae from the three charts.  The pinned
artifact retains the 17-column independent subset.  Its images in
`E(F_p)/3E(F_p)` at

```text
37, 41, 61, 67, 79, 83, 101, 103,
137, 139, 149, 163, 167, 173, 181, 193
```

have full column rank 17 over `F_3`.  The good reduction at `p=13` has group
order 20 and therefore excludes rational 3-torsion.  Infinite descent proves
the 17 points Z-independent, hence `rank E(Q) >= 17` unconditionally.

The strict cutoff needs no floating-point logarithm.  The 76-digit conductor
is less than `10^76`; the degree-seven positive exponential partial sum at
`231/100` is `80381233705038797/8000000000000000 > 10`.  Thus
`ln(10)<231/100` and

\[
\log N < 76\frac{231}{100}=\frac{4389}{25}
<\frac{4568}{25}=182.72.
\]

This is four independent points short of the operational rank-21 target.  No
saturation or rank upper bound is claimed.

## Reproduction

The self-contained certificate is
[the generated rank-17 artifact](../../artifacts/generated-results/elliptic-curves/elliptic_mestre_dsquare_four_u197_rank17.json).
It contains the specialized coefficients, minimal model and conductor, all 17
point representatives, every finite-reduction row, exact cutoff proof, and
the hashes of the untracked discovery outputs.

From the repository root:

```bash
PYTHONPATH=elliptic-curves/cas python3 \
  elliptic-curves/cas/verify_mestre_dsquare_four_u197.py

python3 -m unittest elliptic-curves/tests/test_mestre_dsquare_four.py
```

The first command replays the pinned certificate and independently calls
PARI/GP.  If the local discovery outputs are present, their complete 28-point
replay is added with:

```bash
PYTHONPATH=elliptic-curves/cas python3 \
  elliptic-curves/cas/verify_mestre_dsquare_four_u197.py \
  --discovery-root artifacts/local/elliptic-curves/mestre-dsquare-four-v1
```

The complete bounded discovery is reproduced by:

```bash
PYTHONPATH=elliptic-curves/cas python3 \
  elliptic-curves/cas/search_mestre_dsquare_four.py --workers 8
```

That archival discovery replay expects the locally installed `ratpoints`
bundle at `tmp/ratpoints/root/usr/bin/ratpoints`; the tracked theorem
certificate and its default verifier do not depend on that untracked bundle.
