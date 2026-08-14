# Elkies rank-18 source-recovery audit

Audit date: 2026-08-14.

## Outcome

The published existence claim for Elkies's generic-rank-18 curve was
recovered, but its explicit rank-17 K3 model, 17 sections, quadratic base
change, and eighteenth section were not recovered from any publicly readable
primary source inspected here. The official arXiv bundles contain no ancillary
coefficient data. The most likely historical source is Elkies's September 2007
Keio lecture; its archived page retains only three `.ram` pointers to an
unavailable RTSP host.

The best newly recovered explicit geometric-rank-18 alternative is Kumar and
Kuwata's ancillary basis for Example 9.1. Its constant-field Galois audit is
replayable, but it is not competitive for the present arithmetic search: the
computed fixed rank is 5 over `Q(t)`, and the largest quadratic-character
eigenspace has rank 3. No specialization search was run.

## Target source audit

| Source | Pinned object | What was recovered | Missing data |
|---|---|---|---|
| [arXiv:0709.2908](https://arxiv.org/abs/0709.2908) | official [source bundle](https://export.arxiv.org/e-print/0709.2908), SHA-256 `daa77bdb1da8ea01ebdd7e7add7ed37af5f7057c90e7221b7ddac39fd65f785c` | Existence/construction discussion, Shimura curve `u^2=16t^6-19t^4+88t^2-48`, and the non-CM parameters `t=+/-14/13` | The K3 coefficients, sections, base-change map, and rank-18 basis |
| [arXiv:0802.1301](https://arxiv.org/abs/0802.1301) | official [source bundle](https://export.arxiv.org/e-print/0802.1301), SHA-256 `b8a9cab6c5c723c12d56b11a2060a0e9462bd1c4b5e8742f2cfcc2a7f77fc994` | Reference back to the `X(6,79)` construction and a statement that full results could be made available online | No ancillary files and no link to the model data |
| [BIRS 2018 Elkies slides](https://www.birs.ca/workshops/2018/18w5190/files/Noam_D_Elkies.pdf) | SHA-256 `adc982e6e39cc9c3de49de5f225f6fd28fd4ff19542ce4c72edec027a121a44b` | A rank-17 K3 with `(rho,disc)=(19,948)` is identified as the source of a rank-28 specialization | No Weierstrass model or sections |
| [Elkies NMBRTHRY post](https://listserv.nodak.edu/cgi-bin/wa.exe?A2=NMBRTHRY%3Bb9d018b1.2409&S=b) | SHA-256 `a2234a502ad21de063a1a5b5267f169c0444696131a7b5f93d43fff3870c5ec5` | The rank-29 search used a rank-17 fibration on the same K3; the post says a write-up of the K3 computation is intended | The promised computation/model is not attached |

The two official arXiv inventories and their semantic markers are independently
replayed by
`elliptic-curves/cas/audit_elkies_rank18_sources.py`; its pinned artifact is
`artifacts/generated-results/elliptic_elkies_rank18_source_audit.json`.

## Keio lecture archive

Shioda's [2004 lecture notes](https://www.ms.u-tokyo.ac.jp/publication/docs/lecturenotes04-shioda.pdf),
SHA-256 `9a06ade57ff4f3641daa18368a377a711e0635c2e718a9af1a9f7f81099ce09e`,
cite the relevant source as “N.D. Elkies, Elliptic surfaces and curves of high
rank, Transparency texts and diagrams, Seminar at Keio Univ. Sept. (2007).”

The current Keio [archive index](https://www.math.keio.ac.jp/coe/index.htm) and
[lecture flyer](https://www.math.keio.ac.jp/coe/2007/pathways/elkies.pdf), flyer
SHA-256 `ed474bfa3578bff432586df01da48114dd30ee5c1ddf4a7ba1589cea52d14c28`,
do not contain the transparencies. The NDL WARP
[archived lecture page](https://warp.ndl.go.jp/20110625/20110620071032/http://www.math.hc.keio.ac.jp/coe/videos/elkies2007/)
has normalized SHA-256
`d5099df2cec1fb1e8ed57422837fa93100da5e55c7617c0f3f0d49db2f8ed740`.
It exposes these three pointer files:

- `elkies001.ram`, SHA-256 `93a7ac2846c2e6389c4948d9ed971dbe173ecb42e303dbfa7767b743644d2d2c`;
- `elkies002.ram`, SHA-256 `1b1dbf35f6e2c16f9b0d80400d9ad7cc3737e4554f32b849d8c8c97db222425d`;
- `elkies003.ram`, SHA-256 `7dd42bd6c245302de63425a2ca89d94dc65132a75357a40b326a0d76e3dd9127`.

They point respectively to `elkies001.rm`, `elkies002.rm`, and `elkies003.rm`
on `rtsp://ega.math.hc.keio.ac.jp/pathwaylec/elkies2007/`. During this audit,
the host timed out and the WARP time maps exposed no archived memento for the
stream targets. This is a source blocker, not proof that no copy exists in any
archive.

## Nagao Section 7 and the proposed hidden section

Nagao's primary paper is available from the
[Kobe repository](https://da.lib.kobe-u.ac.jp/da/kernel/E0003610/E0003610.pdf).
The server changes the second PDF trailer ID between requests. After replacing
only that 32-hex-digit ID by zeroes, the stable normalized SHA-256 is
`383ebb58ecd9df35367170fa74b959844beb11e3fe023a5225065cf7493cd2af`; the
stable first ID is `414d0bcd51ca2943bc89211e1b6fc091`.

Section 7 records only that the specialization for
`A=(346,260,255,146,55,0)` and `t=5081/94` has rank at least 20. It does not
state generic rank 13 and does not print a thirteenth section for this tuple.
The paper's explicit generic-rank-13 section belongs to the different tuple
`(148,116,104,57,25,0)`. Exact searches for the Section 7 tuple and parameter
found no further primary formula.

Nagao's 1996 thesis record is at
[Kyushu University](https://catalog.lib.kyushu-u.ac.jp/opac_detail_md/?amode=MD823&bibid=453565)
and [NDL](https://ndlsearch.ndl.go.jp/books/R100000002-I000000305729), DOI
`10.11501/3120538`. The catalog exposes metadata only; NDL marks the text for
on-site/personal transmission rather than public download. It was therefore
not inspected. The matching Artin--Tate square class `-21` in two reductions
is compatible with a section over `Q(sqrt(-21))(t)` and does not by itself
establish a hidden `Q(t)`-section.

## Explicit alternatives recovered

Kuwata's 2000 primary paper is available at
[Rikkyo](https://rikkyo.repo.nii.ac.jp/record/9860/files/AA00610867_49-01_08.pdf),
DOI [`10.14992/00009828`](https://doi.org/10.14992/00009828), SHA-256
`116ff83dc1ef2cc19f7016de8c969daee274f55f3c1d88efa14d290fd27c5cc1`.
It prints four rank-18 elliptic K3 models over `Q`:

```text
Y^2 = X^3 - 33*t^4*X + t*(8*t^10 + 1)
Y^2 = X^3 - 33*t^4*X + 8*t^12 + 1
Y^2 = X^3 + 432*t*(16*t^10 + 44*t^5 - 1)
Y^2 = X^3 + 432*(16*t^12 + 44*t^6 - 1)
```

Here rank 18 is geometric rank over `Qbar(t)`, not arithmetic rank over
`Q(t)`, and that paper does not print an arithmetic 18-section basis.

Kumar and Kuwata later supplied an exact saturated basis in the official
[arXiv:1409.2931 source bundle](https://export.arxiv.org/e-print/1409.2931),
SHA-256 `d3903ffa610826528ee87fa1872a46b7fccb1bcc2b0b2b6240bcdd14d71a73dd`.
The relevant `auxfiles/Example9.1.txt` has SHA-256
`53848ce6b2205353a03b134bf4086bf2ab4d6671d6134a6858800262b4444b2f`.
It gives a geometric rank-18 basis for

```text
Y^2 = X^3 - 33*X + t^6 + 8/t^6,
```

over `Q(i,sqrt(2),3^(1/4))`. The exact finite-fibre Galois/lattice audit in
`elliptic-curves/cas/audit_kumar_kuwata_f6_galois.py` finds:

- fixed rank 5 for the original curve over `Q(t)`;
- quadratic-character ranks `3,2,2,1,1,0,0`, specifically
  `d=2:3`, `d=-1:2`, `d=-3:2`, `d=3:1`, `d=-6:1`, `d=-2:0`, and
  `d=6:0` for the twists by squareclass `d`;
- maximum quadratic-twist rank 3, attained only by squareclass `d=2`.

The action is identified by exact degree-16 number-field arithmetic at
`t=2,3,5,7`, then checked against the Galois group relations and both published
height lattices. This is an exact finite-fibre computational audit, not a
formal symbolic proof of all section identities over the function field. Its
scope and limitation are recorded in the generated artifact.

For comparison, Kloosterman's [arXiv:math/0502439 source](https://export.arxiv.org/e-print/math/0502439),
SHA-256 `34412c4ed858af0c085994be491d4d50e78630d0828ee895a896bcd644c30dfd`,
gives the explicit geometric-rank-15 family

```text
y^2 = x^3 + 2*(t^8 + 14*t^4 + 1)*x + 4*t^2*(t^8 + 6*t^4 + 1).
```

It likewise supplies no rank-15 guarantee over `Q(t)`.

## Replay

```sh
.venv/bin/python elliptic-curves/cas/audit_kumar_kuwata_f6_galois.py \
  --output artifacts/generated-results/elliptic_kumar_kuwata_f6_galois.json
.venv/bin/python -m unittest \
  elliptic-curves/tests/test_kumar_kuwata_f6_galois_audit.py
```

The first command requires `gp` on `PATH` for exact minimal-vector
enumeration. It performs no specialization search and makes no external
state changes beyond writing the named artifact.
