# Jacobian Research

This is the public landing page for the Jacobian research repository. The
proofs, programmes, certificates, and reproducible calculations are gathered
under [research/](research/README.md), so this root remains easy to scan on
GitHub.

- [Research index](research/README.md)
- [Mathematical status](research/MATH_STATUS.json) — authoritative claim ledger
- [Reproduction guide](research/REPRODUCE.md)
- [Elliptic K3 / high-rank programme](research/elkies-k3/README.md)
- [Elliptic-curve programme](research/elliptic-curves/README.md)
- [Elliptic-curve inventory](research/elliptic-curves/INVENTORY.md)

For local work, enter `research/`, create its `.venv`, and install
`requirements.txt`. `make` continues to work from the repository root and
forwards its targets to that directory.

<!-- BEGIN GENERATED ELLIPTIC CURVE TABLE -->
## Elliptic curve inventory

**201 research curves · 129 exact conductors · 72 unresolved**.
Includes seven ICARM matches; the rank-28 row is a public-point reproduction. Rank values are proved lower bounds.

Columns and height conventions follow [ICARM’s table](https://elliptic-rank.icarm.cloud/curves). Logs are natural and shown to two decimals. A dash means the exact conductor is unknown; bounds and partial primes are available on the linked curve page. Coefficients are clipped here; each page contains the complete equation and data.

[Download JSON](research/elliptic-curves/data/research_curves/database.json) · [Download CSV](research/elliptic-curves/data/research_curves/database.csv) · [Arithmetic and replay notes](research/elliptic-curves/notes/INVENTORY201_TABLE_AND_CONDUCTORS_2026-09-07.md)

<details>
<summary>Show all 201 curves</summary>

| Curve | a-invariants | Rank | log N | Naive height | Faltings height | log abs(Δ) |
|---|---|---:|---:|---:|---:|---:|
| [new-20260906-188](research/elliptic-curves/data/research_curves/new-20260906-188.md) [#619](https://elliptic-rank.icarm.cloud/curve/619) | `[1, 0, 0, -1382175065636…, 20242130248207…]` | ≥ 28 | 318.98 | 420.19 | 32.92 | 409.64 |
| [new-20260906-40](https://elliptic-rank.icarm.cloud/curve/630) | `[0, 1, 0, -4780847590649…, 11914273121243…]` | ≥ 27 | 275.84 | 340.97 | 26.37 | 331.42 |
| [new-20260906-41](research/elliptic-curves/data/research_curves/new-20260906-41.md) | `[1, 0, 0, -1534141051320…, 72409877873967…]` | ≥ 27 | 282.27 | 344.47 | 26.57 | 333.09 |
| [new-20260906-186](https://elliptic-rank.icarm.cloud/curve/628) | `[1, -1, 1, -1286302193798…, 17802205930101…]` | ≥ 27 | 282.75 | 350.85 | 27.05 | 338.12 |
| [new-20260906-71](https://elliptic-rank.icarm.cloud/curve/629) | `[1, 0, 0, -2061464727961…, 11053114193386…]` | ≥ 27 | 283.62 | 359.17 | 27.85 | 348.88 |
| [new-20260906-72](research/elliptic-curves/data/research_curves/new-20260906-72.md) | `[1, 0, 0, -6373259334285…, 79686118348726…]` | ≥ 27 | 288.14 | 335.43 | 25.96 | 327.05 |
| [new-20260906-48](research/elliptic-curves/data/research_curves/new-20260906-48.md) | `[1, 0, 0, -3285612539947…, 50166999670818…]` | ≥ 27 | 295.12 | 353.66 | 27.50 | 345.56 |
| [new-20260906-90](https://elliptic-rank.icarm.cloud/curve/627) | `[1, 0, 0, -8810888018746…, 31850830142264…]` | ≥ 27 | 295.67 | 356.62 | 27.47 | 342.37 |
| [new-20260906-92](research/elliptic-curves/data/research_curves/new-20260906-92.md) | `[1, -1, 1, -9220735229473…, 10935366812635…]` | ≥ 26 | 240.05 | 336.06 | 25.89 | 325.06 |
| [new-20260906-74](research/elliptic-curves/data/research_curves/new-20260906-74.md) | `[1, 0, 0, -7397203003011…, 23402190894252…]` | ≥ 26 | 249.51 | 328.47 | 25.31 | 318.57 |
| [new-20260905-37](research/elliptic-curves/data/research_curves/new-20260905-37.md) | `[1, -1, 1, -2712997815437…, 18284322660546…]` | ≥ 26 | 255.18 | 332.49 | 25.65 | 322.87 |
| [new-20260906-104](research/elliptic-curves/data/research_curves/new-20260906-104.md) | `[0, 1, 0, -3846365223002…, 92240176516892…]` | ≥ 26 | 255.34 | 326.51 | 25.04 | 314.36 |
| [new-20260906-49](research/elliptic-curves/data/research_curves/new-20260906-49.md) | `[1, 0, 1, -1969463772621…, 33282711533330…]` | ≥ 26 | 255.51 | 324.50 | 24.91 | 313.19 |
| [new-20260906-192](research/elliptic-curves/data/research_curves/new-20260906-192.md) | `[1, 0, 0, -2868763178245…, 18701058442112…]` | ≥ 26 | 274.78 | 346.35 | 26.54 | 329.78 |
| [new-20260906-42](research/elliptic-curves/data/research_curves/new-20260906-42.md) | `[1, 0, 0, -3302148024364…, 71176199439053…]` | ≥ 26 | 276.89 | 353.68 | 27.38 | 343.23 |
| [new-20260906-63](research/elliptic-curves/data/research_curves/new-20260906-63.md) | `[1, 0, 1, -2445843071846…, 14040344901689…]` | ≥ 26 | 277.01 | 332.05 | 25.61 | 322.20 |
| [new-20260906-75](research/elliptic-curves/data/research_curves/new-20260906-75.md) | `[1, 0, 0, -5433430543708…, 15466891962023…]` | ≥ 26 | 278.76 | 341.36 | 26.27 | 328.89 |
| [new-20260906-73](research/elliptic-curves/data/research_curves/new-20260906-73.md) | `[1, -1, 1, -1135454119938…, 14657867194438…]` | ≥ 26 | 279.02 | 350.47 | 27.04 | 338.34 |
| [new-20260906-103](research/elliptic-curves/data/research_curves/new-20260906-103.md) | `[1, 0, 0, -4735718613418…, 12205008133792…]` | ≥ 26 | 281.16 | 354.76 | 27.47 | 344.37 |
| [new-20260906-106](research/elliptic-curves/data/research_curves/new-20260906-106.md) | `[1, -1, 1, -2225634284449…, 70081981858019…]` | ≥ 26 | 282.76 | 353.59 | 27.51 | 345.74 |
| [new-20260906-50](research/elliptic-curves/data/research_curves/new-20260906-50.md) | `[0, 1, 0, -2429928713309…, 14456788823722…]` | ≥ 26 | — | 332.03 | 25.53 | 320.49 |
| [new-20260906-105](research/elliptic-curves/data/research_curves/new-20260906-105.md) | `[1, 0, 0, -4738382172613…, 11718902875043…]` | ≥ 26 | — | 327.13 | 25.21 | 317.62 |
| [new-20260906-102](research/elliptic-curves/data/research_curves/new-20260906-102.md) | `[0, 1, 0, -2497472982592…, 11650272964469…]` | ≥ 26 | — | 359.75 | 28.00 | 351.41 |
| [new-20260906-91](research/elliptic-curves/data/research_curves/new-20260906-91.md) | `[1, 0, 0, -5452412384522…, 50997797154510…]` | ≥ 26 | — | 375.98 | 29.26 | 365.96 |
| [new-20260906-99](research/elliptic-curves/data/research_curves/new-20260906-99.md) | `[1, 0, 0, -6703231592513…, 21158562400025…]` | ≥ 26 | — | 397.25 | 30.90 | 384.07 |
| [new-20260906-189](research/elliptic-curves/data/research_curves/new-20260906-189.md) | `[1, -1, 1, -1663545350154…, 25999962336300…]` | ≥ 26 | — | 393.07 | 30.59 | 380.88 |
| [new-20260906-54](research/elliptic-curves/data/research_curves/new-20260906-54.md) | `[1, 0, 0, -2903206315434…, 18059949234337…]` | ≥ 25 | 228.23 | 318.75 | 24.50 | 309.00 |
| [new-20260906-43](research/elliptic-curves/data/research_curves/new-20260906-43.md) | `[1, 0, 0, -2430269761279…, 14025401250208…]` | ≥ 25 | 228.74 | 290.59 | 22.14 | 280.54 |
| [new-20260906-79](research/elliptic-curves/data/research_curves/new-20260906-79.md) | `[1, 0, 1, -6937735796735…, 77561917052089…]` | ≥ 25 | 230.23 | 307.75 | 23.61 | 298.56 |
| [new-20260906-82](research/elliptic-curves/data/research_curves/new-20260906-82.md) | `[1, 0, 0, -1993742509891…, 44432685174151…]` | ≥ 25 | 231.36 | 311.24 | 23.95 | 302.88 |
| [new-20260906-77](research/elliptic-curves/data/research_curves/new-20260906-77.md) | `[1, 0, 0, -1685796931974…, 29884403100220…]` | ≥ 25 | 234.06 | 310.44 | 23.85 | 301.41 |
| [new-20260905-38](research/elliptic-curves/data/research_curves/new-20260905-38.md) | `[1, -1, 1, -1934034388435…, 33485883901975…]` | ≥ 25 | 239.82 | 310.67 | 23.79 | 300.10 |
| [new-20260906-66](research/elliptic-curves/data/research_curves/new-20260906-66.md) | `[0, 1, 0, -2089947409724…, 11550142079794…]` | ≥ 25 | 241.16 | 303.95 | 23.18 | 292.20 |
| [new-20260906-110](research/elliptic-curves/data/research_curves/new-20260906-110.md) | `[0, 1, 0, -4320720252161…, 12251842242517…]` | ≥ 25 | 241.92 | 299.45 | 22.93 | 290.41 |
| [new-20260906-107](research/elliptic-curves/data/research_curves/new-20260906-107.md) | `[1, 0, 0, -3015986594896…, 63902234024899…]` | ≥ 25 | 244.96 | 325.78 | 24.95 | 312.96 |
| [new-20260906-121](research/elliptic-curves/data/research_curves/new-20260906-121.md) | `[1, 0, 0, -1873386508355…, 98832501227118…]` | ≥ 25 | 247.81 | 331.26 | 25.39 | 317.93 |
| [new-20260905-01](research/elliptic-curves/data/research_curves/new-20260905-01.md) | `[1, -1, 1, -1059556401049…, 14325193788023…]` | ≥ 25 | 248.48 | 322.79 | 24.86 | 313.38 |
| [new-20260906-81](research/elliptic-curves/data/research_curves/new-20260906-81.md) | `[0, 1, 0, -7231579511386…, 23288656695013…]` | ≥ 25 | 251.03 | 314.58 | 24.10 | 303.68 |
| [new-20260906-114](research/elliptic-curves/data/research_curves/new-20260906-114.md) | `[1, 0, 1, -1973799652407…, 10332174018545…]` | ≥ 25 | 252.53 | 345.23 | 26.69 | 335.01 |
| [new-20260905-03](research/elliptic-curves/data/research_curves/new-20260905-03.md) | `[0, 1, 0, -1785836566136…, 91071136753903…]` | ≥ 25 | 255.46 | 317.29 | 24.30 | 305.77 |
| [new-20260906-53](research/elliptic-curves/data/research_curves/new-20260906-53.md) | `[1, 1, 1, -1073411973939…, 42737476273169…]` | ≥ 25 | 258.90 | 329.58 | 25.26 | 316.37 |
| [new-20260906-120](research/elliptic-curves/data/research_curves/new-20260906-120.md) | `[1, 0, 0, -3499002624654…, 26167673898815…]` | ≥ 25 | 270.20 | 333.20 | 25.69 | 323.13 |
| [new-20260906-118](research/elliptic-curves/data/research_curves/new-20260906-118.md) | `[1, 0, 0, -2757470848258…, 55282180676235…]` | ≥ 25 | 270.38 | 353.14 | 27.28 | 341.55 |
| [new-20260906-44](research/elliptic-curves/data/research_curves/new-20260906-44.md) | `[1, 0, 0, -4261786634678…, 34777610170353…]` | ≥ 25 | 279.56 | 347.59 | 26.87 | 337.17 |
| [new-20260906-200](research/elliptic-curves/data/research_curves/new-20260906-200.md) | `[1, 0, 0, -1600850384897…, 77921536229961…]` | ≥ 25 | 281.00 | 386.04 | 29.92 | 371.68 |
| [new-20260906-78](research/elliptic-curves/data/research_curves/new-20260906-78.md) | `[1, 0, 0, -1108875710517…, 13910170292293…]` | ≥ 25 | 292.47 | 364.22 | 28.25 | 353.60 |
| [new-20260906-113](research/elliptic-curves/data/research_curves/new-20260906-113.md) | `[1, 0, 0, -1346744695971…, 19422896536707…]` | ≥ 25 | 293.96 | 378.66 | 29.45 | 368.00 |
| [new-20260906-116](research/elliptic-curves/data/research_curves/new-20260906-116.md) | `[1, 1, 1, -1917418746764…, 10068368819109…]` | ≥ 25 | 295.69 | 358.95 | 27.80 | 347.97 |
| [new-20260906-117](research/elliptic-curves/data/research_curves/new-20260906-117.md) | `[0, 1, 0, -2007390083834…, 10726523259192…]` | ≥ 25 | 299.12 | 359.09 | 27.82 | 348.42 |
| [new-20260906-52](research/elliptic-curves/data/research_curves/new-20260906-52.md) | `[1, 0, 0, -2367619162358…, 47019730697497…]` | ≥ 25 | 305.99 | 380.43 | 29.65 | 370.77 |
| [new-20260905-02](research/elliptic-curves/data/research_curves/new-20260905-02.md) | `[0, 1, 0, -2294922577037…, 43151908014471…]` | ≥ 25 | — | 297.36 | 22.67 | 286.65 |
| [new-20260906-64](research/elliptic-curves/data/research_curves/new-20260906-64.md) | `[1, 0, 0, -2294485253390…, 13137493111195…]` | ≥ 25 | — | 318.05 | 24.40 | 307.26 |
| [new-20260906-76](research/elliptic-curves/data/research_curves/new-20260906-76.md) | `[1, 0, 1, -1253339005394…, 57810766480572…]` | ≥ 25 | — | 316.37 | 24.32 | 306.85 |
| [new-20260906-101](research/elliptic-curves/data/research_curves/new-20260906-101.md) | `[1, 0, 0, -1786694850006…, 28897710108108…]` | ≥ 25 | — | 324.20 | 24.86 | 312.30 |
| [new-20260906-65](research/elliptic-curves/data/research_curves/new-20260906-65.md) | `[1, 0, 0, -3350508572775…, 77321628884879…]` | ≥ 25 | — | 326.16 | 25.10 | 316.02 |
| [new-20260906-80](research/elliptic-curves/data/research_curves/new-20260906-80.md) | `[1, -1, 1, -5795582957697…, 53759930078101…]` | ≥ 25 | — | 334.64 | 25.66 | 321.04 |
| [new-20260906-111](research/elliptic-curves/data/research_curves/new-20260906-111.md) | `[1, 0, 0, -1597682356873…, 24607383817479…]` | ≥ 25 | — | 337.69 | 25.92 | 324.12 |
| [new-20260906-51](research/elliptic-curves/data/research_curves/new-20260906-51.md) | `[0, 1, 0, -1913863892853…, 97801937897196…]` | ≥ 25 | — | 331.32 | 25.54 | 321.32 |
| [new-20260906-119](research/elliptic-curves/data/research_curves/new-20260906-119.md) | `[1, 0, 0, -3049366837259…, 54419841957103…]` | ≥ 25 | — | 339.62 | 26.30 | 330.95 |
| [new-20260905-39](research/elliptic-curves/data/research_curves/new-20260905-39.md) | `[1, 0, 0, -1384488735192…, 14791922357121…]` | ≥ 25 | — | 351.07 | 27.28 | 342.80 |
| [new-20260906-115](research/elliptic-curves/data/research_curves/new-20260906-115.md) | `[1, -1, 1, -4616790007514…, 12051978149259…]` | ≥ 25 | — | 354.68 | 27.35 | 341.62 |
| [new-20260906-109](research/elliptic-curves/data/research_curves/new-20260906-109.md) | `[1, 0, 0, -1413662557339…, 66332414636833…]` | ≥ 25 | — | 358.09 | 27.75 | 347.61 |
| [new-20260906-112](research/elliptic-curves/data/research_curves/new-20260906-112.md) | `[1, 0, 1, -9445216260008…, 36055694157722…]` | ≥ 25 | — | 329.24 | 25.33 | 318.56 |
| [new-20260906-108](research/elliptic-curves/data/research_curves/new-20260906-108.md) | `[0, 1, 0, -7414389144151…, 24140726469159…]` | ≥ 25 | — | 356.10 | 27.57 | 345.29 |
| [new-20260906-129](research/elliptic-curves/data/research_curves/new-20260906-129.md) | `[1, 0, 0, -1117426831395…, 43673028426873…]` | ≥ 24 | 218.32 | 288.26 | 21.95 | 278.24 |
| [new-20260906-127](research/elliptic-curves/data/research_curves/new-20260906-127.md) | `[1, 0, 0, -3428548612984…, 36377001247496…]` | ≥ 24 | 218.91 | 292.42 | 22.40 | 284.36 |
| [new-20260906-56](research/elliptic-curves/data/research_curves/new-20260906-56.md) | `[1, 0, 0, -2574004602740…, 17125620800021…]` | ≥ 24 | 222.54 | 290.91 | 22.20 | 281.48 |
| [new-20260906-94](research/elliptic-curves/data/research_curves/new-20260906-94.md) | `[1, 0, 0, -7739802787653…, 79702914556097…]` | ≥ 24 | 226.98 | 321.69 | 24.73 | 311.65 |
| [new-20260906-126](research/elliptic-curves/data/research_curves/new-20260906-126.md) | `[0, 1, 0, -1995056018165…, 10414654163957…]` | ≥ 24 | 230.99 | 276.18 | 20.94 | 266.17 |
| [new-20260906-85](research/elliptic-curves/data/research_curves/new-20260906-85.md) | `[1, 0, 0, -2729141972253…, 54281574759682…]` | ≥ 24 | 232.00 | 311.66 | 23.84 | 300.37 |
| [new-20260905-07](research/elliptic-curves/data/research_curves/new-20260905-07.md) | `[1, 0, 0, -3947112552094…, 52338831877962…]` | ≥ 24 | 232.60 | 306.96 | 23.62 | 299.10 |
| [new-20260905-06](research/elliptic-curves/data/research_curves/new-20260905-06.md) | `[1, 0, 0, -4267429672548…, 99032113848645…]` | ≥ 24 | 234.53 | 299.18 | 22.89 | 289.82 |
| [new-20260905-33](research/elliptic-curves/data/research_curves/new-20260905-33.md) | `[1, 0, 0, -1307563683235…, 16719966175686…]` | ≥ 24 | 242.08 | 323.27 | 24.90 | 313.95 |
| [new-20260905-09](research/elliptic-curves/data/research_curves/new-20260905-09.md) | `[1, -1, 1, -3500362298043…, 24696101281366…]` | ≥ 24 | 248.49 | 319.31 | 24.51 | 308.64 |
| [new-20260906-133](research/elliptic-curves/data/research_curves/new-20260906-133.md) | `[1, 1, 1, -8866059595848…, 10448976726496…]` | ≥ 24 | 248.53 | 308.34 | 23.61 | 297.97 |
| [new-20260906-83](research/elliptic-curves/data/research_curves/new-20260906-83.md) | `[1, 1, 1, -2661829654680…, 66308254642014…]` | ≥ 24 | 250.33 | 325.85 | 25.16 | 317.39 |
| [new-20260906-95](research/elliptic-curves/data/research_curves/new-20260906-95.md) | `[1, 0, 0, -3651786279336…, 11839717098456…]` | ≥ 24 | 252.70 | 313.20 | 24.12 | 305.02 |
| [new-20260906-143](research/elliptic-curves/data/research_curves/new-20260906-143.md) | `[0, 1, 0, -2980460120883…, 19594862083066…]` | ≥ 24 | 253.90 | 305.02 | 23.29 | 293.70 |
| [new-20260906-137](research/elliptic-curves/data/research_curves/new-20260906-137.md) | `[1, 0, 0, -1677449248023…, 25756310241953…]` | ≥ 24 | 254.00 | 337.83 | 26.06 | 327.41 |
| [new-20260906-87](research/elliptic-curves/data/research_curves/new-20260906-87.md) | `[0, 1, 0, -2164196205172…, 46300576290568…]` | ≥ 24 | 258.12 | 325.13 | 25.09 | 316.47 |
| [new-20260906-197](research/elliptic-curves/data/research_curves/new-20260906-197.md) | `[1, 0, 0, -3400460644190…, 24469282084097…]` | ≥ 24 | 263.06 | 360.70 | 27.94 | 349.64 |
| [new-20260906-130](research/elliptic-curves/data/research_curves/new-20260906-130.md) | `[1, 0, 0, -3232044040627…, 68516021242878…]` | ≥ 24 | 263.26 | 325.98 | 25.08 | 315.74 |
| [new-20260906-135](research/elliptic-curves/data/research_curves/new-20260906-135.md) | `[1, 0, 1, -2743349956370…, 17605033030357…]` | ≥ 24 | 264.43 | 332.41 | 25.55 | 320.62 |
| [new-20260905-08](research/elliptic-curves/data/research_curves/new-20260905-08.md) | `[1, 0, 0, -2697433537331…, 17159737861992…]` | ≥ 24 | 264.60 | 346.18 | 26.69 | 334.34 |
| [new-20260906-132](research/elliptic-curves/data/research_curves/new-20260906-132.md) | `[1, 0, 0, -1824691664341…, 28842900192061…]` | ≥ 24 | 268.40 | 338.08 | 26.10 | 328.05 |
| [new-20260906-128](research/elliptic-curves/data/research_curves/new-20260906-128.md) | `[0, 1, 0, -2430184479771…, 20218699907363…]` | ≥ 24 | 270.26 | 332.69 | 25.75 | 324.50 |
| [new-20260906-45](research/elliptic-curves/data/research_curves/new-20260906-45.md) | `[0, 1, 0, -5588454304361…, 49148140969994…]` | ≥ 24 | 270.71 | 320.72 | 24.65 | 310.54 |
| [new-20260906-193](research/elliptic-curves/data/research_curves/new-20260906-193.md) | `[0, 1, 0, -2067788020150…, 36475186922808…]` | ≥ 24 | 278.18 | 352.29 | 27.21 | 340.67 |
| [new-20260906-86](research/elliptic-curves/data/research_curves/new-20260906-86.md) | `[1, 0, 0, -1812698883857…, 27934229722153…]` | ≥ 24 | 280.05 | 338.06 | 26.12 | 328.45 |
| [new-20260906-55](research/elliptic-curves/data/research_curves/new-20260906-55.md) | `[1, 0, 0, -2173884683749…, 13454337848007…]` | ≥ 24 | 281.20 | 359.50 | 27.92 | 350.21 |
| [new-20260906-139](research/elliptic-curves/data/research_curves/new-20260906-139.md) | `[1, 0, 0, -9847579800901…, 10946309691352…]` | ≥ 24 | 283.76 | 336.23 | 25.98 | 326.90 |
| [new-20260906-134](research/elliptic-curves/data/research_curves/new-20260906-134.md) | `[0, 1, 0, -5288141449191…, 15667139535776…]` | ≥ 24 | 286.29 | 369.02 | 28.70 | 359.33 |
| [new-20260906-97](research/elliptic-curves/data/research_curves/new-20260906-97.md) | `[1, 0, 0, -1313271189490…, 59892209008991…]` | ≥ 24 | 289.36 | 344.07 | 26.59 | 333.88 |
| [new-20260906-196](research/elliptic-curves/data/research_curves/new-20260906-196.md) | `[1, 0, 0, -1240823078810…, 53158040941219…]` | ≥ 24 | 293.64 | 371.46 | 28.72 | 357.56 |
| [new-20260906-131](research/elliptic-curves/data/research_curves/new-20260906-131.md) | `[1, -1, 1, -8630149593675…, 97860332866838…]` | ≥ 24 | 294.53 | 363.47 | 28.10 | 350.84 |
| [new-20260905-10](research/elliptic-curves/data/research_curves/new-20260905-10.md) | `[1, 0, 0, -9941757705488…, 12011979816231…]` | ≥ 24 | — | 308.63 | 23.55 | 296.45 |
| [new-20260906-123](research/elliptic-curves/data/research_curves/new-20260906-123.md) | `[1, 0, 0, -2198240392712…, 39077091419249…]` | ≥ 24 | — | 324.83 | 24.95 | 313.85 |
| [new-20260906-140](research/elliptic-curves/data/research_curves/new-20260906-140.md) | `[1, 0, 0, -2619197841993…, 16231107460154…]` | ≥ 24 | — | 304.63 | 23.22 | 292.60 |
| [new-20260905-04](research/elliptic-curves/data/research_curves/new-20260905-04.md) | `[1, 0, 0, -4069458158284…, 99419607730772…]` | ≥ 24 | — | 312.86 | 23.91 | 300.80 |
| [new-20260906-96](research/elliptic-curves/data/research_curves/new-20260906-96.md) | `[1, -1, 1, -3421843534327…, 24831278244948…]` | ≥ 24 | — | 305.47 | 23.35 | 294.73 |
| [new-20260905-05](research/elliptic-curves/data/research_curves/new-20260905-05.md) | `[0, 0, 0, -1325752240777…, 58208200456254…]` | ≥ 24 | — | 302.59 | 23.08 | 291.14 |
| [new-20260906-125](research/elliptic-curves/data/research_curves/new-20260906-125.md) | `[0, 1, 0, -1315092262900…, 58082690700794…]` | ≥ 24 | — | 316.38 | 24.12 | 302.21 |
| [new-20260906-136](research/elliptic-curves/data/research_curves/new-20260906-136.md) | `[1, 0, 0, -4043691855637…, 98134444950998…]` | ≥ 24 | — | 326.65 | 25.08 | 315.12 |
| [new-20260906-57](research/elliptic-curves/data/research_curves/new-20260906-57.md) | `[1, 1, 1, -6436606488599…, 18805343670831…]` | ≥ 24 | — | 328.05 | 25.28 | 318.34 |
| [new-20260906-122](research/elliptic-curves/data/research_curves/new-20260906-122.md) | `[1, 0, 0, -1302443483471…, 18653957646221…]` | ≥ 24 | — | 323.32 | 24.86 | 313.04 |
| [new-20260906-84](research/elliptic-curves/data/research_curves/new-20260906-84.md) | `[1, 0, 0, -3155622638734…, 21559871690156…]` | ≥ 24 | — | 319.00 | 24.35 | 305.06 |
| [new-20260906-67](research/elliptic-curves/data/research_curves/new-20260906-67.md) | `[1, 0, 0, -6426079574194…, 21005251302244…]` | ≥ 24 | — | 341.97 | 26.44 | 332.30 |
| [new-20260906-93](research/elliptic-curves/data/research_curves/new-20260906-93.md) | `[1, -1, 1, -2136819806863…, 40938394829016…]` | ≥ 24 | — | 338.70 | 26.18 | 329.27 |
| [new-20260906-194](research/elliptic-curves/data/research_curves/new-20260906-194.md) | `[1, 0, 0, -1795680405567…, 28583980700528…]` | ≥ 24 | — | 351.85 | 27.23 | 341.35 |
| [new-20260906-124](research/elliptic-curves/data/research_curves/new-20260906-124.md) | `[1, 0, 0, -1719616684361…, 29269725196439…]` | ≥ 24 | — | 338.03 | 26.12 | 328.46 |
| [new-20260906-142](research/elliptic-curves/data/research_curves/new-20260906-142.md) | `[1, 0, 0, -5649562455320…, 17061296868837…]` | ≥ 24 | — | 355.37 | 27.54 | 345.42 |
| [new-20260906-138](research/elliptic-curves/data/research_curves/new-20260906-138.md) | `[0, 1, 0, -1140167018600…, 50664698433370…]` | ≥ 24 | — | 357.55 | 27.75 | 348.16 |
| [new-20260906-46](research/elliptic-curves/data/research_curves/new-20260906-46.md) | `[0, 1, 0, -2451091942964…, 14729131656777…]` | ≥ 24 | — | 345.88 | 26.64 | 333.23 |
| [new-20260906-141](research/elliptic-curves/data/research_curves/new-20260906-141.md) | `[0, 1, 0, -2568174338905…, 15911676041496…]` | ≥ 24 | — | 359.84 | 27.82 | 347.66 |
| [new-20260906-162](research/elliptic-curves/data/research_curves/new-20260906-162.md) | `[1, 0, 0, -3230152894772…, 27145586613507…]` | ≥ 23 | 198.13 | 250.38 | 18.87 | 241.80 |
| [new-20260905-12](research/elliptic-curves/data/research_curves/new-20260905-12.md) [#600](https://elliptic-rank.icarm.cloud/curve/600) | `[1, 0, 0, -2673631332732…, 16765556663649…]` | ≥ 23 | 205.75 | 277.06 | 20.91 | 264.68 |
| [new-20260906-150](research/elliptic-curves/data/research_curves/new-20260906-150.md) | `[1, 0, 0, -1139693092949…, 14458937284597…]` | ≥ 23 | 211.71 | 295.22 | 22.51 | 284.71 |
| [new-20260906-161](research/elliptic-curves/data/research_curves/new-20260906-161.md) | `[0, -1, 0, -7495108039460…, 78286697905930…]` | ≥ 23 | 214.82 | 252.52 | 18.90 | 241.02 |
| [new-20260906-156](research/elliptic-curves/data/research_curves/new-20260906-156.md) | `[1, 0, 0, -2827155291397…, 66384479645955…]` | ≥ 23 | 222.75 | 298.22 | 22.84 | 289.34 |
| [new-20260905-34](research/elliptic-curves/data/research_curves/new-20260905-34.md) | `[1, 0, 0, -4164297933286…, 34099946132101…]` | ≥ 23 | 224.39 | 319.92 | 24.59 | 309.94 |
| [new-20260906-68](research/elliptic-curves/data/research_curves/new-20260906-68.md) | `[1, 0, 0, -2377369611505…, 13611534170192…]` | ≥ 23 | 226.20 | 304.34 | 23.28 | 294.21 |
| [new-20260906-195](research/elliptic-curves/data/research_curves/new-20260906-195.md) | `[1, 0, 0, -7820082157402…, 41015532490086…]` | ≥ 23 | 230.39 | 315.68 | 24.34 | 307.68 |
| [new-20260905-11](research/elliptic-curves/data/research_curves/new-20260905-11.md) | `[1, -1, 1, -2252486937486…, 41143497194498…]` | ≥ 23 | 231.40 | 297.27 | 22.47 | 281.22 |
| [new-20260905-20](research/elliptic-curves/data/research_curves/new-20260905-20.md) | `[1, -1, 1, -6742501320324…, 21424925831398…]` | ≥ 23 | 236.27 | 286.75 | 21.73 | 274.76 |
| [new-20260906-157](research/elliptic-curves/data/research_curves/new-20260906-157.md) | `[1, -1, 1, -1408079364950…, 68041512951665…]` | ≥ 23 | 246.72 | 330.51 | 25.49 | 320.82 |
| [new-20260905-13](research/elliptic-curves/data/research_curves/new-20260905-13.md) | `[1, -1, 1, -5853400087974…, 52978911118699…]` | ≥ 23 | 249.59 | 320.86 | 24.65 | 310.51 |
| [new-20260906-146](research/elliptic-curves/data/research_curves/new-20260906-146.md) | `[1, 0, 0, -1610440628057…, 24558192145702…]` | ≥ 23 | 266.96 | 337.71 | 26.02 | 326.58 |
| [new-20260906-191](research/elliptic-curves/data/research_curves/new-20260906-191.md) | `[1, -1, 1, -6477959409604…, 20115266830084…]` | ≥ 23 | 269.90 | 328.07 | 25.14 | 315.26 |
| [new-20260906-158](research/elliptic-curves/data/research_curves/new-20260906-158.md) | `[1, 0, 0, -1012016381028…, 39839083266819…]` | ≥ 23 | 271.73 | 343.25 | 26.49 | 332.37 |
| [new-20260906-98](research/elliptic-curves/data/research_curves/new-20260906-98.md) | `[1, -1, 1, -4181985850957…, 10617824252209…]` | ≥ 23 | 276.37 | 326.79 | 25.13 | 316.09 |
| [new-20260906-159](research/elliptic-curves/data/research_curves/new-20260906-159.md) | `[1, 0, 0, -1720067963124…, 27818990195436…]` | ≥ 23 | 276.85 | 337.93 | 26.04 | 326.82 |
| [new-20260906-163](research/elliptic-curves/data/research_curves/new-20260906-163.md) | `[1, 0, 0, -1882488345207…, 31047247388211…]` | ≥ 23 | 283.66 | 338.18 | 26.06 | 327.02 |
| [new-20260906-153](research/elliptic-curves/data/research_curves/new-20260906-153.md) | `[0, 0, 0, -1221508605764…, 16795651673703…]` | ≥ 23 | 284.37 | 350.74 | 27.13 | 340.13 |
| [new-20260906-58](research/elliptic-curves/data/research_curves/new-20260906-58.md) | `[1, 0, 0, -1112832408013…, 14096997571593…]` | ≥ 23 | 285.07 | 336.60 | 25.93 | 325.52 |
| [new-20260906-155](research/elliptic-curves/data/research_curves/new-20260906-155.md) | `[1, -1, 1, -4838793547941…, 12772587366610…]` | ≥ 23 | 285.12 | 341.01 | 26.30 | 329.98 |
| [new-20260906-149](research/elliptic-curves/data/research_curves/new-20260906-149.md) | `[1, 0, 0, -1349138305607…, 58257951475835…]` | ≥ 23 | 287.33 | 357.90 | 27.75 | 347.74 |
| [new-20260906-154](research/elliptic-curves/data/research_curves/new-20260906-154.md) | `[1, 0, 0, -7718651134788…, 25289202550782…]` | ≥ 23 | 288.22 | 356.22 | 27.60 | 345.98 |
| [new-20260906-164](research/elliptic-curves/data/research_curves/new-20260906-164.md) | `[1, 0, 0, -7221629950470…, 23572775966072…]` | ≥ 23 | 294.24 | 356.02 | 27.47 | 343.07 |
| [new-20260906-60](research/elliptic-curves/data/research_curves/new-20260906-60.md) | `[0, 1, 0, -1719129188397…, 26744869524665…]` | ≥ 23 | 294.51 | 351.72 | 27.22 | 341.26 |
| [new-20260906-199](research/elliptic-curves/data/research_curves/new-20260906-199.md) | `[0, -1, 0, -2914075924876…, 17035551226407…]` | ≥ 23 | 307.35 | 360.21 | 28.00 | 351.19 |
| [new-20260906-187](research/elliptic-curves/data/research_curves/new-20260906-187.md) | `[0, 1, 0, -8385913570454…, 29012491279858…]` | ≥ 23 | 311.87 | 370.29 | 28.75 | 359.53 |
| [new-20260906-89](research/elliptic-curves/data/research_curves/new-20260906-89.md) | `[1, 0, 0, -6267191125441…, 61123619929522…]` | ≥ 23 | — | 307.27 | 23.48 | 296.08 |
| [new-20260905-18](research/elliptic-curves/data/research_curves/new-20260905-18.md) | `[0, 1, 0, -5076486514628…, 44423038931901…]` | ≥ 23 | — | 306.63 | 23.41 | 295.15 |
| [new-20260906-151](research/elliptic-curves/data/research_curves/new-20260906-151.md) | `[1, 0, 0, -9667688991118…, 11022338501743…]` | ≥ 23 | — | 308.55 | 23.65 | 298.71 |
| [new-20260906-59](research/elliptic-curves/data/research_curves/new-20260906-59.md) | `[1, 0, 0, -3162549426908…, 21169341254272…]` | ≥ 23 | — | 319.01 | 24.48 | 308.42 |
| [new-20260906-147](research/elliptic-curves/data/research_curves/new-20260906-147.md) | `[0, 0, 0, -5122075337152…, 14101619749366…]` | ≥ 23 | — | 313.55 | 23.88 | 299.32 |
| [new-20260906-70](research/elliptic-curves/data/research_curves/new-20260906-70.md) | `[1, -1, 0, -1148664320555…, 48293606439039…]` | ≥ 23 | — | 316.01 | 24.23 | 305.26 |
| [new-20260905-15](research/elliptic-curves/data/research_curves/new-20260905-15.md) | `[1, 0, 0, -2412778483332…, 44611178931588…]` | ≥ 23 | — | 352.74 | 27.30 | 342.15 |
| [new-20260905-14](research/elliptic-curves/data/research_curves/new-20260905-14.md) | `[1, 0, 0, -3379115227832…, 34155989418695…]` | ≥ 23 | — | 319.92 | 24.68 | 311.79 |
| [new-20260906-88](research/elliptic-curves/data/research_curves/new-20260906-88.md) | `[1, -1, 1, -2067795634560…, 11255749467158…]` | ≥ 23 | — | 317.73 | 24.37 | 306.86 |
| [new-20260906-145](research/elliptic-curves/data/research_curves/new-20260906-145.md) | `[1, 0, 1, -3814722202196…, 30886229735545…]` | ≥ 23 | — | 333.54 | 25.75 | 324.10 |
| [new-20260906-160](research/elliptic-curves/data/research_curves/new-20260906-160.md) | `[1, -1, 1, -4964097511666…, 13375026222408…]` | ≥ 23 | — | 354.90 | 27.42 | 343.09 |
| [new-20260906-152](research/elliptic-curves/data/research_curves/new-20260906-152.md) | `[0, 1, 0, -5910250997331…, 55857924302997…]` | ≥ 23 | — | 334.72 | 25.76 | 323.34 |
| [new-20260906-144](research/elliptic-curves/data/research_curves/new-20260906-144.md) | `[1, 0, 0, -5070870113187…, 14972565780909…]` | ≥ 23 | — | 341.30 | 26.40 | 331.86 |
| [new-20260905-19](research/elliptic-curves/data/research_curves/new-20260905-19.md) | `[1, 0, 0, -3693473800452…, 27354798579804…]` | ≥ 23 | — | 347.11 | 26.71 | 333.64 |
| [new-20260906-190](research/elliptic-curves/data/research_curves/new-20260906-190.md) | `[1, 0, 0, -1510167912489…, 71494663866974…]` | ≥ 23 | — | 372.06 | 28.77 | 358.27 |
| [new-20260907-201](research/elliptic-curves/data/research_curves/new-20260907-201.md) | `[1, 0, 0, -5779692239484…, 53344140755227…]` | ≥ 23 | — | 389.90 | 30.30 | 377.17 |
| [new-20260905-17](research/elliptic-curves/data/research_curves/new-20260905-17.md) | `[0, 1, 0, -1542808627338…, 23479098197119…]` | ≥ 23 | — | 365.22 | 28.28 | 353.43 |
| [new-20260905-16](research/elliptic-curves/data/research_curves/new-20260905-16.md) | `[1, -1, 1, -5693101643316…, 18063073394371…]` | ≥ 23 | — | 369.30 | 28.74 | 360.03 |
| [new-20260906-148](research/elliptic-curves/data/research_curves/new-20260906-148.md) | `[0, 1, 0, -4002927137327…, 97079439842254…]` | ≥ 23 | — | 368.07 | 28.50 | 355.81 |
| [new-20260906-69](research/elliptic-curves/data/research_curves/new-20260906-69.md) | `[0, 1, 0, -9270078556395…, 10822713038341…]` | ≥ 23 | — | 363.68 | 28.13 | 351.34 |
| [new-20260905-36](https://elliptic-rank.icarm.cloud/curve/626) | `[1, 0, 0, -1824519766025…, 71000315025379…]` | ≥ 22 | 174.43 | 234.47 | 17.56 | 226.19 |
| [new-20260906-181](research/elliptic-curves/data/research_curves/new-20260906-181.md) | `[1, -1, 1, -5860321714788…, 17761183148463…]` | ≥ 22 | 179.38 | 244.93 | 18.32 | 234.57 |
| [new-20260905-32](research/elliptic-curves/data/research_curves/new-20260905-32.md) | `[1, 0, 0, -3618975169760…, 26426786986322…]` | ≥ 22 | 186.42 | 250.34 | 18.67 | 237.67 |
| [new-20260905-25](research/elliptic-curves/data/research_curves/new-20260905-25.md) | `[1, 0, 0, -3319768398911…, 95350749122919…]` | ≥ 22 | 189.82 | 257.50 | 19.47 | 249.14 |
| [new-20260905-29](research/elliptic-curves/data/research_curves/new-20260905-29.md) | `[1, 0, 0, -2883710112976…, 59903201646862…]` | ≥ 22 | 205.59 | 284.20 | 21.52 | 272.14 |
| [new-20260905-27](research/elliptic-curves/data/research_curves/new-20260905-27.md) | `[1, 0, 0, -4149207237066…, 31844205037279…]` | ≥ 22 | 231.00 | 292.19 | 22.25 | 281.56 |
| [new-20260905-26](research/elliptic-curves/data/research_curves/new-20260905-26.md) | `[0, 0, 0, -1160513473040…, 47080765238776…]` | ≥ 22 | 239.01 | 288.37 | 21.93 | 277.76 |
| [new-20260906-167](research/elliptic-curves/data/research_curves/new-20260906-167.md) | `[1, 0, 0, -1328080246714…, 58602998026106…]` | ≥ 22 | 246.94 | 316.41 | 24.20 | 304.38 |
| [new-20260905-24](research/elliptic-curves/data/research_curves/new-20260905-24.md) | `[0, 1, 0, -8722313749028…, 31145701431169…]` | ≥ 22 | 247.61 | 301.33 | 22.96 | 289.55 |
| [new-20260905-22](research/elliptic-curves/data/research_curves/new-20260905-22.md) | `[1, 0, 0, -8465384516501…, 30037954732579…]` | ≥ 22 | 250.54 | 315.06 | 24.05 | 302.06 |
| [new-20260905-28](research/elliptic-curves/data/research_curves/new-20260905-28.md) | `[1, 0, 1, -1154218044306…, 48102665846221…]` | ≥ 22 | 250.61 | 302.19 | 23.04 | 290.56 |
| [new-20260905-21](research/elliptic-curves/data/research_curves/new-20260905-21.md) | `[1, 0, 1, -2793535264635…, 18068004952949…]` | ≥ 22 | 262.09 | 304.83 | 23.24 | 292.84 |
| [new-20260906-61](research/elliptic-curves/data/research_curves/new-20260906-61.md) | `[1, 0, 0, -4344853162488…, 34825353767971…]` | ≥ 22 | 264.00 | 333.78 | 25.59 | 320.06 |
| [new-20260906-182](research/elliptic-curves/data/research_curves/new-20260906-182.md) | `[0, 1, 0, -3715952530259…, 88275691895514…]` | ≥ 22 | 276.10 | 326.43 | 25.08 | 315.26 |
| [new-20260906-173](research/elliptic-curves/data/research_curves/new-20260906-173.md) | `[1, -1, 0, -1307663998272…, 18479916622205…]` | ≥ 22 | 276.33 | 323.30 | 24.82 | 312.34 |
| [new-20260906-184](research/elliptic-curves/data/research_curves/new-20260906-184.md) | `[1, 0, 0, -4070386412748…, 97600281585843…]` | ≥ 22 | 281.78 | 326.67 | 25.13 | 316.15 |
| [new-20260906-180](research/elliptic-curves/data/research_curves/new-20260906-180.md) | `[1, -1, 1, -5021161073238…, 23950656040417…]` | ≥ 22 | 289.19 | 348.03 | 27.05 | 340.21 |
| [new-20260906-174](research/elliptic-curves/data/research_curves/new-20260906-174.md) | `[0, 1, 0, -4863170889922…, 41593109262160…]` | ≥ 22 | 290.44 | 334.13 | 25.70 | 322.48 |
| [new-20260905-31](research/elliptic-curves/data/research_curves/new-20260905-31.md) | `[1, 0, 0, -2671116209953…, 16819388980554…]` | ≥ 22 | 294.91 | 373.77 | 28.92 | 360.07 |
| [new-20260905-30](research/elliptic-curves/data/research_curves/new-20260905-30.md) | `[1, 0, 0, -5984090776937…, 17434344010559…]` | ≥ 22 | 298.06 | 369.28 | 28.67 | 358.66 |
| [new-20260906-179](research/elliptic-curves/data/research_curves/new-20260906-179.md) | `[0, 1, 0, -6474392930107…, 11481778223250…]` | ≥ 22 | 298.86 | 355.70 | 27.69 | 347.85 |
| [new-20260906-172](research/elliptic-curves/data/research_curves/new-20260906-172.md) | `[0, 0, 0, -9839149617583…, 12153084261771…]` | ≥ 22 | 302.05 | 363.91 | 28.23 | 353.34 |
| [new-20260906-198](research/elliptic-curves/data/research_curves/new-20260906-198.md) | `[0, 1, 0, -1189330734292…, 14210689411316…]` | ≥ 22 | 303.63 | 401.06 | 31.47 | 393.47 |
| [new-20260906-47](research/elliptic-curves/data/research_curves/new-20260906-47.md) | `[1, 0, 0, -5353402805519…, 16299172758828…]` | ≥ 22 | 305.98 | 369.10 | 28.72 | 359.71 |
| [new-20260906-165](research/elliptic-curves/data/research_curves/new-20260906-165.md) | `[1, -1, 1, -7431792442473…, 24626838413328…]` | ≥ 22 | 307.19 | 369.93 | 28.61 | 356.54 |
| [new-20260906-177](research/elliptic-curves/data/research_curves/new-20260906-177.md) | `[0, 1, 0, -3945704642437…, 28923117563071…]` | ≥ 22 | 322.20 | 374.93 | 29.17 | 364.96 |
| [new-20260906-100](research/elliptic-curves/data/research_curves/new-20260906-100.md) | `[0, 1, 0, -1092859705443…, 11547168992023…]` | ≥ 22 | 329.29 | 391.81 | 30.65 | 383.18 |
| [new-20260905-23](research/elliptic-curves/data/research_curves/new-20260905-23.md) | `[1, 0, 0, -4315415756119…, 34825741607359…]` | ≥ 22 | — | 306.14 | 23.37 | 294.69 |
| [new-20260906-166](research/elliptic-curves/data/research_curves/new-20260906-166.md) | `[1, 0, 0, -1572876882936…, 76398008240300…]` | ≥ 22 | — | 316.93 | 24.26 | 305.07 |
| [new-20260905-35](research/elliptic-curves/data/research_curves/new-20260905-35.md) | `[0, 1, 0, -5614225203918…, 15025074802140…]` | ≥ 22 | — | 341.45 | 26.41 | 332.03 |
| [new-20260906-168](research/elliptic-curves/data/research_curves/new-20260906-168.md) | `[1, 0, 0, -5121214450763…, 43480846311916…]` | ≥ 22 | — | 320.46 | 24.61 | 310.00 |
| [new-20260906-169](research/elliptic-curves/data/research_curves/new-20260906-169.md) | `[1, 0, 0, -1623411907128…, 80083724511219…]` | ≥ 22 | — | 344.65 | 26.56 | 332.75 |
| [new-20260906-171](research/elliptic-curves/data/research_curves/new-20260906-171.md) | `[0, 1, 0, -6141588560520…, 62031149318571…]` | ≥ 22 | — | 348.75 | 27.01 | 339.07 |
| [new-20260906-62](research/elliptic-curves/data/research_curves/new-20260906-62.md) | `[1, 0, 0, -2588484100489…, 43231998909016…]` | ≥ 22 | — | 352.95 | 27.41 | 344.19 |
| [new-20260906-178](research/elliptic-curves/data/research_curves/new-20260906-178.md) | `[1, 0, 0, -1221330541006…, 16231053681057…]` | ≥ 22 | — | 350.69 | 27.10 | 339.50 |
| [new-20260906-185](research/elliptic-curves/data/research_curves/new-20260906-185.md) | `[1, 0, 0, -6889419529793…, 21968763016388…]` | ≥ 22 | — | 369.70 | 28.61 | 356.66 |
| [new-20260906-170](research/elliptic-curves/data/research_curves/new-20260906-170.md) | `[0, 1, 1, -9435946073212…, 33990140963693…]` | ≥ 22 | — | 343.01 | 26.51 | 332.92 |
| [new-20260906-183](research/elliptic-curves/data/research_curves/new-20260906-183.md) | `[1, 0, 0, -2410326942603…, 10732128134291…]` | ≥ 22 | — | 345.83 | 26.84 | 337.56 |
| [new-20260906-175](research/elliptic-curves/data/research_curves/new-20260906-175.md) | `[1, -1, 1, -4178960739123…, 31658080454283…]` | ≥ 22 | — | 375.11 | 29.18 | 365.04 |
| [new-20260906-176](research/elliptic-curves/data/research_curves/new-20260906-176.md) | `[1, 0, 1, -6924559956107…, 21911629220112…]` | ≥ 22 | — | 369.71 | 28.68 | 358.53 |

</details>

<!-- END GENERATED ELLIPTIC CURVE TABLE -->
