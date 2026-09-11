# Accessibility closure and the independent 11952 rank-bound target

Current research direction and the first bounded results are in the
[2026-09-12 ancestry / arithmetic / strict-class comparison](RANK_TRIANGLE_ANCESTRY_AND_DESCENT_2026-09-12.md).
The closure theorem and rank-25 criterion below remain valid. Its earlier
single-route priority has been replaced by that three-lane investigation.

## Conclusion

Close the original signed-unit/pair experiment as **sparse one-step
accessibility improvements with diminishing returns**, not an avalanche or
an intrinsic Curve302 rank-amplification mechanism. The main arithmetic target
is now an independent class-group 2-rank upper bound for the **11952 fibre
at `t=921/653`**, not the different 11952 rank-25 fibre at `102/1525`.

For the selected fibre, the new equation-only calculation establishes

\[
\operatorname{rank}E(\mathbf Q)\le\dim\operatorname{Sel}_2(E)
\le g+9,\qquad g=\dim_{\mathbf F_2}\operatorname{Cl}(K)[2].
\]

Consequently **an independently certified `g <= 16` would prove rank 25**
when combined with the existing 25-point independence certificate.
No such upper bound has been obtained. The initial 512 MB BNF probe failed
before class-group invariants. The separately authorized8 GiB attempt cleared
that overflow but hit its600-second time limit during relation collection,
without a provisional BNF or cyclic invariants. Quotient certification was
not admitted. The present rank claim remains **rank at least25**, with no
new GRH-conditional or unconditional upper bound.

The lower bound belongs to case `b-7bb187bc9254e81c6283`, batch002:
the [25 exact points](../../artifacts/local/elliptic-curves/broad-rank-v1/runtime/research/broad-cases/b-7bb187bc9254e81c6283/batch-002/search-00/terminal.json)
are bound by the [independent replay](../../artifacts/local/elliptic-curves/broad-rank-v1/runtime/research/broad-cases/b-7bb187bc9254e81c6283/batch-002/search-00/verified.json)
with status `PASS_INDEPENDENT_PRODUCTIVE_V3_REPLAY` and terminal SHA256
`3a50e6d7c92ebc6ab6c0d74bf08079c8dc9dd93906ed71cb62760cfaa3fab676`.

The [new verifier and certificate](../../artifacts/generated-results/elliptic-curves/accessibility_rank_pivot_v1/verified.json)
check the finite rule reduction and the integer arithmetic of this bound
criterion. The original exact slope/group-law replay remains its provenance,
and PARI remains the implementation used for maximal-order and local reduction
certification. This is not an independently implemented full descent.

## The exact finite-atlas conclusion

The [sealed all-subset experiment](../../artifacts/generated-results/elliptic-curves/rank_accessibility_subsets_v1/REPORT.md)
contains these five rules on 302:

\[
E_1\longrightarrow E_{10},E_{14},\qquad
E_2\longrightarrow E_7,E_{10},\qquad
E_6\longrightarrow E_{12}.
\]

The source set `{1,2,6}` and target set `{7,10,12,14}` are disjoint.
There is no directed path of length two. The single original-basis rule on
11952 is `E4 -> E2`, likewise with no two-step path. These are rules for
improved chart costs, not certificates that an algorithm acquires a point.

For each hidden target `j`, let `H0_j` be its smallest generic chart parameter
height and `H_ij` its smallest height among charts with exceptional support
exactly `{i}`. The new checker verifies that every support of size two is
dominated by one of its singleton supports or by the empty support. It also
checks the following equality against all **114,688 + 1,024** stored cells:

\[
\kappa_T(P_j)=\min\bigl(k_j^0,\min_{i\in T}k_{ij}\bigr),
\quad
f_j(T):=k_j^0-\kappa_T(P_j)=\max\bigl(\{0\}\cup\{w_{ij}:i\in T\}\bigr).
\]

Here `k=log2(H)` and `w_ij=max(0,k0_j-k_ij)`. All comparisons are made with
exact integer heights, not rounded logarithms. For `A` contained in `B`,

\[
f_j(A\cup\{i\})-f_j(A)=\max(0,w_{ij}-f_j(A))
\ge\max(0,w_{ij}-f_j(B)).
\]

Thus each improvement function is monotone submodular: diminishing returns,
not cooperative enhancement. Positive subset cells repeat these few rules.
The five-rule system does not explain the historical 13/14 recovery.
No conclusion is asserted for other mixed targets, atlases, or bases.

Retain the [four exact event decompositions](../../artifacts/generated-results/elliptic-curves/rank_accessibility_subsets_v1/events.json)
unchanged. Their residual-height, finite-cancellation and coordinate-distortion
contributions differ; there is no common fitted explanation. The largest
parameter-height ratio is approximately 229,097, not a measured runtime gain.
The [common integral trajectory core](CURVE302_EXACT_SHARED_CORE_2026-09-11.md)
also remains established. This experiment does not explain that core.

## Audit of existing descent evidence

| Existing evidence | What it establishes | What it does not establish |
|---|---|---|
| This fibre's 25-point replay, and displayed `M17` inside `D25` | Rank at least 25; primitive inclusion inside the displayed lattice; known quotient rank 8 | Full-group saturation, full exceptional quotient, or rank upper bound |
| [302 descent anatomy](CURVE302_DESCENT_ANATOMY_2026-09-10.md) | `dim Sel2 = 21+cS`; known strict classes give `cS >= 10`; associated ordinary class lower information | An upper bound for `cS`, or an independent measurement distinguishing available ranks |
| [Small-conductor descent diagnostics](SMALL_CONDUCTOR_DESCENT_SHORTCUT_2026-09-06.md) | Certified order/local calculations and retained failed global computations for MW16 at `3/17` | A certificate for the current 11952 fibre |
| Bounded point-search continuations | Results within their stated chart budgets | Nonexistence of further points or upper rank bounds |

The measured class-group object below is constructed from the equation alone.
No exceptional point, Kummer image, trajectory, or known-rank generator is
passed to either arithmetic worker. The source point artifact is used by the
preparer only to extract and hash-bind the curve equation. In particular,
`g >= 16`, which follows from the known rank and `rank <= g+9`, would be a
**point-derived lower bound**, not an independent class-group measurement.

## Equation-derived criterion for 11952 at 921/653

The original short equation is `y^2=x^3+A*x+B`, with

```text
A = -84196621147240320617047979399750704971769668763392
B = 206937637967498573533732205060490773473330035432339348502618816917963546624
```

Its global minimal model is

```text
y^2 = x^3 + x^2 + a4*x + a6
a4 = -4060408041437129659386958883089829522172534180
a6 = 69302996254333102097575929282499740768182420645699427883330979661600
x_short = 144*x_min + 48; y_short = 1728*y_min
```

Let `K=Q(theta)`, where `theta^3+theta^2+a4*theta+a6=0`. This cubic has no
root modulo 23, proving irreducibility and the absence of rational 2-torsion.
The equation-only field diagnostic certifies the maximal order using
`nfcertify(nf)=[]`, after proving primality of the supplied factor hints.

```text
signature(K) = (3,0)
[O_K : Z[theta]] = 2112935718757500
disc(K) = 30932024439246801322787810077880389223317162788400703970784702350579246728265890412316825281477068142364105
```

The local terms in the Brumer–Kramer bound are:

| Term | Exact input | Contribution |
|---|---|---:|
| `u` | Positive discriminant | 2 |
| Multiplicative primes with even minimal discriminant valuation | `3,7,13,19,83` | 5 |
| Additive primes | Only `2`, with three primes of `K` above it | `3-1=2` |
| Total `u+n` | `2+5+2` | **9** |

The minimal discriminant valuations at these multiplicative primes are
`8,4,4,2,2`. The additive reduction at 2 has conductor exponent 3 and
minimal discriminant valuation 8. The full factorization, all local splitting
types, exact basis and timings are in
[field.log](../../artifacts/generated-results/elliptic-curves/11952_rank_bound_audit_v1/field.log).
The root number is `-1`; no algebraic-rank parity assumption is used.

The bound is [Proposition 3.1 of Klagsbrun–Sherman–Weigandt](https://arxiv.org/html/1606.07178#S3.SS1),
attributed there to Brumer–Kramer. Their work supplies a relevant precedent
for proving class-group 2-rank bounds rather than computing everything about
the class group. A mod-2 relation matrix gives a global upper bound only
after factor-base generation is justified; their stated generation bounds
use GRH. No GRH-dependent output is to be relabelled unconditional.

## Fixed feasibility probe and its boundary

The [frozen protocol](../../artifacts/generated-results/elliptic-curves/11952_rank_bound_audit_v1/protocol.json)
allows one 30-second field diagnostic and one 60-second class-group probe,
one worker, 1.5 GiB RSS and 512,000,000 bytes of PARI stack. No point search
or automatic retry is allowed. PARI version: `2.17.3`.

The field diagnostic completed in about 0.13 seconds, including maximal-order
certification. The class-group probe failed in about 0.39 seconds with a
**PARI stack overflow in `bnfinit`**, before producing invariants. Later GP
lines continued after that error; the zero process exit and printed `DONE`
marker are therefore explicitly rejected by the result validator. The
supervisor's generic stack field has its default 256 MB value; the frozen GP
argv and the error log establish the actual 512 MB stack used.

The retained [class-group log](../../artifacts/generated-results/elliptic-curves/11952_rank_bound_audit_v1/class_group.log)
contains no certified class group, no conditional class-group upper bound,
and no rank upper bound. `bnfcertify` was never successfully reached with a
BNF object. Its intended flag-1 success would certify that the true class
group is a quotient of the computed group, sufficient for an unconditional
2-rank upper bound; see the [PARI documentation](https://pari.math.u-bordeaux.fr/dochtml/html/General_number_fields.html#bnfcertify).

The next arithmetic gate is the separately budgeted quotient-only lane below.
A larger/unresolved rank bound would not establish that additional rational
points exist.

## Larger-memory quotient-only lane

The user explicitly prioritized this lane and paused the height-ball
benchmark. The [new frozen protocol](../../artifacts/generated-results/elliptic-curves/11952_class_quotient_v1/protocol.json)
uses the accessible host `Roy5080` (31 GiB total RAM, approximately28 GiB
available at preflight), one worker, an **8 GiB PARI stack**, and a **12 GiB
RSS cap**. Provisional `bnfinit(nf,0)` gets at most600 seconds. Only a completed
provisional 2-rank at most16 admits one additional `bnfcertify(bnf,1)` stage,
capped at300 seconds. No automatic retry, time extension, full-unit
certification or point search is authorized by this protocol. The source
maximal order is hash-bound and its polynomial, discriminant and maximality
are checked again before BNF. The earlier failed attempt is not overwritten.

The quotient argument is sufficient: if `C_comp` surjects onto `Cl(K)`, then
`C_comp/2*C_comp` surjects onto `Cl(K)/2*Cl(K)`. For a finite abelian group,
the latter dimension equals the dimension of its 2-torsion. Hence counting
the even cyclic factors of `C_comp` gives the required upper bound after
quotient certification. No equality of class numbers, regulator certificate
or certified fundamental units is needed. This saves work at the
**certification** stage; it does not bypass provisional BNF relation collection.
The original512 MB probe already intended to use flag1 but never reached it.

| Completed output | Admissible conclusion |
|---|---|
| Provisional BNF with `g_comp=16`, no successful quotient certificate | GRH-conditional exact rank25 only |
| Same, and `bnfcertify(bnf,1)=1` | Unconditional `g<=16`, rank at most25, and exact rank25 after combining the independent lower bound |
| Provisional `g_comp>16` | `NOT_SUFFICIENT`; stop before certification; no unconditional inference that true `g>16` |
| Provisional `g_comp<16` | Inconsistent with the certified lower bound and local criterion; fail closed for investigation |
| Incomplete BNF | No provisional class2-rank and no new conditional or unconditional bound |

The [launch worker](../cas/run_11952_class_quotient.py) freezes inputs and
captures process-group resource limits. A [separate output validator](../cas/validate_11952_class_quotient.py)
retains the raw result while allowing only two known benign PARI diagnostic
forms that the launcher's conservative `***` filter rejects: the exact8 GiB
stack-size announcement and the numeric Bach-constant diagnostic. All other
warnings/errors, missing markers, changed invariant factors, missing BNF
checkpoints and non-completed supervision results remain failures. This
validator does not change the running GP program, resource limits or
arithmetic acceptance criteria. Its ten synthetic tests cover the success,
conditional-only, insufficient, inconsistent, missing-checkpoint, duplicate
marker, backend-error and resource-limit cases.

Completed execution evidence is in
[provisional.log](../../artifacts/generated-results/elliptic-curves/11952_class_quotient_v1/provisional.log)
and the
[supervisor checkpoint](../../artifacts/generated-results/elliptic-curves/11952_class_quotient_v1/provisional.supervisor.json).
The8 GiB stack passed the previous initialization failure and reached relation
collection over9,170 factor-base ideals. The supervisor stopped it at the
600-second limit (601.01 seconds including termination grace), with peak
observed RSS1,448,108,032 bytes, approximately1.35 GiB. It produced **no
provisional BNF checkpoint, no cyclic invariants and no provisional2-rank**.
The [validated result](../../artifacts/generated-results/elliptic-curves/11952_class_quotient_v1/provisional.validation.json)
is `RESOURCE_LIMIT`; both conditional and unconditional rank-upper fields
are null. No quotient-certification call was launched. No relation count is
itself a class-group2-rank or a proof of a global class bound.

This identifies relation collection, not the observed initial memory failure,
as the unfinished stage within this budget. It does not predict how long
completion would take or rule out later memory demands. The next resource
decision is a separately bounded longer provisional-BNF run or an explicitly
scoped relation/refinement approach, not merely more stack allocation. The
present attempt is stopped, not scheduled to resume; its incomplete internal
BNF state is not a resumable arithmetic checkpoint. All logs and protocol
evidence remain available, and the unrelated broad search is untouched.

One-shot execution in a clean output location:

```bash
python3 research/elliptic-curves/cas/run_11952_class_quotient.py prepare
python3 research/elliptic-curves/cas/run_11952_class_quotient.py provisional
python3 research/elliptic-curves/cas/validate_11952_class_quotient.py
# Only after status PROVISIONAL_GRH_SUFFICIENT:
python3 research/elliptic-curves/cas/validate_11952_class_quotient.py --run-quotient
```

Cheap gate regression and retained-result replay:

```bash
python3 -m unittest discover -s research/elliptic-curves/tests -p test_11952_class_quotient.py
python3 research/elliptic-curves/cas/validate_11952_class_quotient.py --check
```

## Height-bounded benchmark: paused by the user

The [small fixed benchmark gate](../../artifacts/generated-results/elliptic-curves/accessibility_rank_pivot_v1/height_benchmark_gate.json)
specifies `A_H(S)`, both signs and nonprimitive multiples, fixed anchor-only
normalization, absolute endpoint versus chart count, exact metric transport,
and a 120-second/one-worker cap. It is now **paused**, not queued behind the
arithmetic run. It was not run before the pause: the atlas contains
numerical Gram entries, not rigorously bounded canonical-height enclosures.
Decimal formatting and a matrix symmetry check supply no error certificate.
More displayed digits or re-enumerating a frozen rounded form does not repair this gap.

A rigorous lower bounding form and boundary membership enclosures are
required before claiming complete canonical-height-ball coverage. That
optional benchmark stays gated; it must not displace the equation-derived
rank work. No rounded-metric substitute is represented as its completion.

## Reproduction

The original one-shot arithmetic attempt is retained and cannot be silently
overwritten. In a clean output location, its preparation and execution are:

```bash
python3 research/elliptic-curves/cas/audit_11952_rank_bound.py prepare
python3 research/elliptic-curves/cas/audit_11952_rank_bound.py field
python3 research/elliptic-curves/cas/audit_11952_rank_bound.py class_group
```

The cheap certificate check is safe on the current repository:

```bash
python3 research/elliptic-curves/cas/verify_accessibility_rank_pivot.py --check
```

This check replays exact model changes, discriminants, the mod-23 test, bound
arithmetic, source hashes, singleton domination and every original subset
cell. It does not repeat the preceding 149,916 rational slope checks or
pretend that a failed class computation is a certificate.
