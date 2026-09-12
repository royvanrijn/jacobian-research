# Elliptic K3 script map

Start with the [K3 programme](../README.md), [claim catalogue](../../index/elkies-k3.md)
and [algorithmic memory](../../knowledge/ALGORITHMS.md). The selected claim's
canonical source, exact inputs and software lock determine the replay.
Script names, directory placement and old `ACTIVE_SEARCH` labels do not
establish that a worker is running or that another campaign is needed.

| Task | Current source or implementation |
|---|---|
| Load a completed rational MW17 equation | [Published/q12 endpoint](../PRIOR_WORK_SHORTCUT_AND_ENDPOINT_CERTIFICATION_2026-08-25.md) · [direct alternate11952](../R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md) |
| Reuse the Curve302 parent | [Full recovered equation and basis](../../elliptic-curves/notes/CURVE302_RECOVERED_MW17_PARENT_2026-09-07.md) |
| Reuse Curve398 recovery and presentation deduplication | [Canonical recovery](../../elliptic-curves/notes/ICARM_CURVE398_RANK30_AND_CONSTRUCTION.md) |
| Plan a marked fibration hop | [Rank/lift theorem contracts](../RANK_MUTATION_AND_LIFT_THEOREMS.md) · [marked-U planner](../MARKED_U_REALIZATION_PLANNER_2026-09-03.md) |
| Understand route corrections | [Process atlas](../ELKIES_K3_PROCESS_ATLAS.md) · [historical H3 prefix](success-path/README.md) |
| Compile the noncyclic4A1 bridge | [Exact equation and certificate](../R17_NONCYCLIC_4A1_DIRECT_FIBRATION_2026-09-04.md) |
| Screen a different NS lattice | [Full marked-curve gate](../DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md) · [classifier](../RANK19_ARITHMETIC_MARKING_CLASSIFIER_2026-09-04.md) |
| Reuse finite-field and formal foundry experiments | [Retained outcomes and exact NS0031 precursor](../LATTICE_FOUNDRY_EQUATION_FIRST_SHORTLIST_2026-09-02.md) |
| Decide the determinant1236 next step | [Known candidate cover; open CM branch-orbit identification](../DET1236_MARKED_SHIMURA_CURVE_2026-09-04.md) |
| Use rootless genus filters | [Genus theory and bounds](../ROOTLESS_GENUS_THEORY_2026-09-03.md) |
| Interpret bisection/character searches | [Scoped failures](../../knowledge/FAILED_ROUTES.md) · [elliptic rank-jump programme](../../elliptic-curves/README.md) |
| Reproduce orbit1222 | [Its separate characteristic-zero pipeline](q80-orbit1222-char0/README.md) |

<a id="fixed-corridor-reverse-lift-from-the-q12o5867-endpoint"></a>
The [fixed-corridor reverse lift](../FIXED_CORRIDOR_REVERSE_LIFT_2026-08-26.md)
is retained with its exact route and marking hypotheses. A fallback or historical
intermediate equation is not a missing requirement for the completed endpoint.

## Find a script without repeating its discovery

From the repository root:

```sh
python3 research/scripts/research.py search "q8 denominator"
python3 research/scripts/research.py search "target coset"
python3 research/scripts/research.py routes --area elkies-k3
python3 research/scripts/research.py show EC-K3-H3-Q12O5867-ENDPOINT-QQ
```

Run legacy relative-path CAS commands from `research/`. The local Sage10.9
environment, if installed, is `/home/royvanrijn/.local/share/jacobian-sage-10.9`.
Some files are Python with Sage imports; others require Sage preparsing.
Use the exact invocation in the canonical note and [replay guide](../../REPRODUCE.md).

`--check` is not a cost guarantee: some programs reconstruct class fields,
enumerate full coset populations or spawn hundreds of lifts before comparing
outputs. Inspect the mode and reuse retained evidence first. Exactly verified
independent points establish lower bounds; incomplete descent cannot veto a
separately scoped point search or establish an upper bound.

The [original 4,394-line map](../../archive/repository-cleanup-2026-09-12/research__elkies-k3__scripts__README.md.txt)
preserves every previous script description and historical command. It contains
superseded frontiers. [Older scripts](archive/README.md) and exact proof inputs
remain available; consult current status before using them.

<!-- status-consumer: EC-K3-DET1236-GENUS2-RATIONAL-POINTS 5a3c84eb9f7f0604 -->
<!-- status-consumer: EC-K3-DET1236-MARKED-SHIMURA-CURVE e482668e1208f764 -->
<!-- status-consumer: EC-K3-DET1236-V4-LOCAL-CONSISTENCY deb09ed7326145f9 -->
<!-- status-consumer: EC-K3-DET1236-CANDIDATE-DOUBLE-COVER fecec75b4ff1e8e1 -->
<!-- status-consumer: EC-K3-DET1236-RATIONAL-CM-LOCUS bd6ab0e86ca70ab2 -->
<!-- status-consumer: EC-K3-R17-NORM12-HIGHEST-RANK-TRANSPORTS c4c8a81fc735fea2 -->
<!-- status-consumer: EC-K3-R17-NORM12-PROSPECTIVE-FAMILY-HOLDOUT 8fb7417663ea1d98 -->
<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-MOD131-HORIZONTAL 2249c509c1217d7c -->
<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-QQ-HORIZONTAL 688f0a5f6d989e9c -->
<!-- status-consumer: EC-K3-ELKIES-2026-R17 9208e67f51fc8c97 -->
<!-- status-consumer: EC-K3-ELKIES-2026-R18-COVER 6b4ee5bbc1afc01e -->
<!-- status-consumer: EC-K3-ELKIES-2026-R19-PAIRED f1e135d2ba803e80 -->
<!-- status-consumer: EC-K3-BISECT-BIQUADRATIC-R19 707bffd8b85f8f3e -->
<!-- status-consumer: EC-K3-H3-Q12O5867-POINT-FACTORY 9399c93ee42ee2a4 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-TWO-PRIMARY-BOUNDARY 783482d8f700105d -->
<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-LINEAR-CONDUCTOR 957479f39bedd57b -->
<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-SPECIALIZED-QUARTICS 725664f9e36ae8a7 -->
<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-GENERIC-QUARTICS aa704dc4685e4c9b -->
<!-- status-consumer: EC-K3-H3-ROOTLESS-J2-COMPLETE c6f054948b04b507 -->
<!-- status-consumer: EC-K3-H3-ROOTLESS-J1-UNIFORM-BOUND b71330a75ad2c9ad -->
<!-- status-consumer: EC-K3-CUSTOM-NS-HALF-LATTICE-SWEEP 9dc0e4d23f677392 -->
<!-- status-consumer: EC-K3-NS0024-DIRECT-QQ-INOSE-OBSTRUCTION e87afc1b3529a07f -->
<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: EC-K3-NS0031-MARKED-SOURCE-PRECURSOR 2e115b35c30a8cea -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
<!-- status-consumer: EC-K3-NS0031-MARKED-RATIONAL-PARAMETER-SCAN ca678e520745dd3c -->
<!-- status-consumer: EC-K3-NS0031-QQ-MARKING-OBSTRUCTION 8e2dc35cdf9b6bc3 -->
<!-- status-consumer: EC-K3-GOLAY-DET720-QQ-MARKING-OBSTRUCTION 972f591d2885f9ba -->
<!-- status-consumer: EC-K3-DIFFERENT-NS-ARITHMETIC-GATE-RERANK fd9549d1fcb2e9e7 -->
<!-- status-consumer: EC-K3-RANK19-ARITHMETIC-MARKING-CLASSIFIER eec5710ee1b498ab -->
<!-- status-consumer: EC-K3-DET378-QQ-MARKING-OBSTRUCTION 1e910f72f54ac228 -->
<!-- status-consumer: EC-K3-ARITHMETIC-FIRST-MARKED-T-FOUNDRY 2f7b65d586e96394 -->
<!-- status-consumer: EC-K3-DET500-DET750-QQ-MARKING-OBSTRUCTIONS 14498ad134ffa60e -->
<!-- status-consumer: EC-K3-LATTICE-FOUNDRY-PRESCRIBED-ROOT-MW1-CENSUS 01298fec30fa94a3 -->
<!-- status-consumer: EC-K3-UNIVERSAL-DEGREE2-FIBRATION-COMPILER fd4b5d71c9497eaf -->
<!-- status-consumer: EC-K3-H3-Q4O208-Q4O1599-QQ-A3-2A2 0937fa8496fb99a6 -->
<!-- status-consumer: EC-K3-H3-Q4O208-Q4O323-QQ-A3-2A2 d16db02fc62ee239 -->
<!-- status-consumer: EC-K3-H3-Q4O164-Q8O376-QQ-SMOOTH-RR c0bc752848961743 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-DEGREE1-SECTION-QQ 28056772646e9fc7 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-QQ-ROOTLESS 6a3dbee5942ddd0a -->
<!-- status-consumer: EC-K3-H3-Q12O5867-QQ-R17-BASIS a2097150acf00645 -->
<!-- status-consumer: EC-K3-H3-Q12O5867-ENDPOINT-QQ a83b08acd921c32b -->
<!-- status-consumer: EC-K3-H3-FIXED-REVERSE-A1-2A1-QQ 26bce707a77972c4 -->
<!-- status-consumer: EC-K3-H3-FIXED-REVERSE-4A1-QQ 8f3863b630d27e16 -->
<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-EQUATION 827d75cb8d14d7f4 -->
<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-ARITHMETIC-RANK2 387d6237125637a3 -->
<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT103-SPECIALIZATION-RANK7 bf1d025228805b31 -->
<!-- status-consumer: EC-K3-E6A1-RHO19-ORBIT96-A7D7-GALOIS ba008502f0e5533f -->
<!-- status-consumer: EC-K3-RES-QBC-E6-II-RANK3-RHO19 5b10608e230145e9 -->
<!-- status-consumer: EC-K3-RES-QBC-E6-II-Q2-MW4 3aa5084463780acc -->
<!-- status-consumer: EC-K3-RES-D6-RATIONALIZED-SECTION-CHART a94042dd2d76797c -->
<!-- status-consumer: EC-K3-RES-QBC-E6-RANK4-LINEAR-CHORD 3bcfe3534656b26f -->
<!-- status-consumer: EC-K3-E6-RANK4-ROOTLESS-Q2Q4-CENSUS 2351738f44774cfe -->
<!-- status-consumer: EC-K3-E6-RANK4-DET78-GLOBAL-ROOTFUL 648ec884ce7152bb -->
<!-- status-consumer: EC-K3-ROOTLESS-GENUS-MASS e7589727ca8f7e50 -->
<!-- status-consumer: EC-K3-ARITHMETIC-RANK-TRANSFER 3031dd2365a29cd5 -->
<!-- status-consumer: EC-K3-R17-ALTERNATE-Q80-ARITHMETIC-RANK17 a304934727bb3f87 -->
<!-- status-consumer: EC-K3-RELATIVE-U-BRIDGE-LIFTING 800e22abf69b91aa -->
<!-- status-consumer: EC-K3-LOCAL-BRIDGE-MUTATION-H1C 7421328afadcf61f -->
<!-- status-consumer: EC-K3-PRIME-LOCAL-BRIDGE-MUTATION-H1D d4c6c84967a8fbc5 -->
<!-- status-consumer: EC-K3-NS0024-RELATIVE-U-FIRST-EDGE-OBSTRUCTION d57544697149506f -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-CALCULUS 7eeeeaa80d9b2bf3 -->
<!-- status-consumer: EC-K3-INTEGRAL-CHARACTER-GLUE 0b76d65366279037 -->
<!-- status-consumer: EC-K3-R17-NORM12-103B2-INTEGRAL-GLUE 52de13c8443f2b7d -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-BRIDGE-PREDICTOR-BENCHMARK 3127e24cc505f646 -->
<!-- status-consumer: EC-K3-E6-DET78-PROSPECTIVE-BRIDGE-NEGATIVE d23a0abd146c2ed9 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-THETA-CONVOLUTION 5ebbd3d242fdb3db -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-CORE-GENERATION d0d78c49b44f55ac -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-REVERSE-THETA eee16ce986ec0a1f -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-WEIL-COMPRESSION 34d2abea91a265f4 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MODULAR-DIMENSION-SIEVE 9622c6eb4d8522bd -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MASKED-CORE-GENERATION 9a7a1e01cb22f62e -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MASKED-CORE-CONTROLS 3cbde45fb2cb0f17 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-DIRECTED-Q80 80de8b6727cd3409 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-BIRTH-DEATH a755a3956c4c97cb -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-ROOT-SYSTEM-SIGNATURE d32b35b66a35627c -->
<!-- status-consumer: EC-K3-NS0024-INVERSE-ADE-MUTATION 5c56f07d14129837 -->
<!-- status-consumer: EC-K3-INVERSE-ADE-PROJECTIVE-BIRTH-STRATA b4a7edb452e6dcc7 -->
<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-GRAPH-REACHABILITY e02f950eba79b32a -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-SINGLETON-PO0-TOP200 80ab545a98b4e2d7 -->
<!-- status-consumer: EC-K3-RES-D5-TWO-MARKED-TWO-TWIST-POLYNOMIAL ea0496c9566cfdc3 -->
<!-- status-consumer: EC-K3-RES-D5-TWO-MARKED-LOW-SLICE-ELIMINANTS 43d297285eb3655b -->
<!-- status-consumer: EC-K3-RES-A4-TWO-POINT-TATE-SLICE-OBSTRUCTION b9729a0a8f2f17be -->
<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-BISECTION-PILOT 80fa6e59107cc9e6 -->
<!-- status-consumer: EC-K3-R17-MULTISECTION-VISIBILITY-FILTRATION 5569f42aac3ab952 -->
<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-SIMULTANEOUS-SPLITTING-H10000 40fb0bc465e3e95c -->
<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-MIXED-TRACE-SPLITTING-H10000 c7aa09836b842b60 -->
<!-- status-consumer: EC-K3-R17-GENUS1-HIGH-THROUGHPUT-SPLITTING cad3d98ce58c89e7 -->
<!-- status-consumer: EC-K3-R17-NORM12-103B2-MW-LATTICE-SIEVE aa0d8718eb57de6f -->
<!-- status-consumer: EC-K3-R17-NORM12-103B2-HARD-FIBRE-PRODUCT-H300000 b4fef7ab54b922e0 -->
<!-- status-consumer: EC-K3-R17-NORM12-103B2-ISOTROPIC-FRAME 47f3a0eb7ee50bcb -->
<!-- status-consumer: EC-K3-R17-POINTED-COVER-JACOBIAN-CONTROL-H10000 4bb087b3a1ebc684 -->
<!-- status-consumer: EC-K3-H3-ROOTLESS-J2-MINIMAL-ACCESSIBILITY 631f50389e0a3283 -->
<!-- status-consumer: EC-K3-R17-NORM12-ICARM-573-REFRESH a93ce35de34fde21 -->
<!-- status-consumer: EC-K3-R17-074D9-RIGID-CROSS-FIBRE-TRANSFER abbedd192865f172 -->
<!-- status-consumer: EC-K3-R17-074D9-NORM8-CROSS-FIBRE-TRANSFER-16 262e405b0adbbb73 -->
<!-- status-consumer: EC-K3-R17-074D9-LATE-POINT-HOLDOUT 284e0f92def23419 -->
<!-- status-consumer: EC-K3-R17-074D9-RECORD-TWIST-MW-OBSTRUCTION c794f827e9a8ac36 -->
<!-- status-consumer: EC-K3-R17-074D9-LOCAL-KUMMER-SEPARATION 375cb897b59e077f -->
<!-- status-consumer: EC-K3-R17-074D9-QUOTIENT-ARITHMETIC-BLOCK-OBSTRUCTION af45468d1b7d831a -->
<!-- status-consumer: EC-K3-R17-074D9-PROSPECTIVE-CRT-LOCAL-STABILITY 0edaaa6f05041634 -->
<!-- status-consumer: EC-K3-R17-074D9-PROSPECTIVE-CRT-ESCAPE-EXPERIMENT 021a952efb9ea0f4 -->
<!-- status-consumer: EC-K3-R17-074D9-HALF-LATTICE-PROMOTION-GATE 9a1f080523c9ecae -->
<!-- status-consumer: EC-K3-R17-074D9-QUOTIENT-RANK-ESCAPE-DETECTOR-V2 eda7a0053b31b7c9 -->
<!-- status-consumer: EC-K3-R17-KUMMER-CLASSGROUP-PRESSURE-COMPARISON 74b1dae24470b531 -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-DIRECT-Q80-EQUATION 077c6409d76cbe63 -->
<!-- status-consumer: EC-K3-R17-NONCYCLIC-4A1-DIRECT-EQUATION f657620e07f8f3f0 -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-CONTROL-J-PREIMAGES 1ef38474a0d7f629 -->
<!-- status-consumer: EC-K3-CURVE398-A1-MW16-RECOVERY 75978a18cc26690f -->
<!-- status-consumer: EC-K3-CURVE398-TWO-PARENT-COLLISION 626a440519ff77f3 -->
<!-- status-consumer: EC-K3-ICARM-A1-MW16-ATLAS df43ab6cf2f46d4e -->
<!-- status-consumer: EC-K3-R17-EXTREME-ANCHORED-MW18-COVERS 4763babcbf5d923c -->
<!-- status-consumer: EC-K3-R17-PRODUCT-19BAD-083AD-ARITHMETIC-RANK-ZERO fe572bd5979b5d2c -->
<!-- status-consumer: EC-K3-R17-PRODUCT-REGULATOR-OBSTRUCTION-SWEEP f86dead53d55babe -->
