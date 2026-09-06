# Two individually soluble native covers cannot split together over Q

**Solubility:** the native R17 parameter carrier A,D is everywhere locally
soluble but has **no rational point**. Its particular torsor class has a
nonzero Cassels–Tate pairing. This closes the class-identification question
left open by the [preceding carrier comparison](GLOBAL_CARRIER_SOLUBILITY_AND_SPECIALIZATION.md).
It gives a global obstruction to assembling two specified directions,
without exceptional points or a search for new original-family parameters.
The choice of these cover labels remains retrospective.

Write A=`orbit-030cb`, D=`orbit-11278`. Their primitive equations are
\[
\begin{aligned}
u^2=f_A(t)&=4865126421024+2514185838528t+320914613929t^2,\\
v^2=f_D(t)&=-144492039-201200094t+18383017t^2.
\end{aligned}
\]
Each conic has a retained rational point separately:
\[
(t,u)=(-5877/1690,6956427/130),\qquad
(t,v)=(-70665/36751,647424714/36751).
\]
Nevertheless, no rational t makes both values squares. More strongly, the
smooth projective genus-one carrier \(C_{AD}\) has no rational point,
including over infinity. Its local nonemptiness at every place was already
proved with complete bad-prime support in the preceding certificate.

## Identifying the obstruction, rather than just finding Sha nearby

The carrier's minimal Jacobian has Weierstrass coefficients
```
[1, 0, 0,
 -22346724249689819463954606277705347361,
 40659077139706222675001173966969442586598667246064627785]
```
Its certified auxiliary rank is 2, rational 2-torsion dimension is 1,
full 2-Selmer dimension is 5, and Sha[2] dimension is 2. Those dimensions
alone did not identify the carrier class. Earlier bounded attempts returned
too few rational generators to settle membership in the rational Kummer
subgroup; their frozen UNKNOWN outputs are retained.

The new calculation observes the **labelled** full Selmer basis used inside
the pinned PARI 2.17.3 descent. In that basis the alternating Cassels matrix,
with entries in F2, is
\[
M=\begin{pmatrix}
0&1&1&0&0\\
1&0&1&0&0\\
1&1&0&0&0\\
0&0&0&0&0\\
0&0&0&0&0
\end{pmatrix}.
\]
The earlier exact degree-four Jacobian map and cubic square identity supply
the carrier's descent representative beta. Transporting it to this model
and evaluating finite split-root characters gives
\[
[\beta]=e_3,\qquad e_3^TM=(1,1,0,0,0).
\]
The character map has rank five on the complete Selmer basis, so the
coordinates are unique: this is not an inference from an incomplete point
subgroup. Beta belongs to Selmer by the verified cover map and complete
local solubility. A nonzero pairing annihilates the possibility that beta
comes from a rational point. Thus its image in Sha is nonzero and
\(C_{AD}(\mathbf Q)=\varnothing\).

As a consistency check, the matrix radical has dimension three and is
spanned by \(e_1+e_2+e_3,e_4,e_5\). Given the preceding exact rank and
torsion dimensions, it equals the rational Kummer subgroup. The negative
conclusion itself only needs the nonzero pairing, not this identification
of the entire radical or recovery of the missing rational generator.

## An explicit reciprocity certificate

The pairing entry (3,1) was separately expanded using
[Fisher's binary-quartic Cassels–Tate formula, Theorem 3.1](https://www.dpmms.cam.ac.uk/~taf1000/papers/bq-ctp.pdf).
The three quartics represent \(e_3,e_1,e_3+e_1\), independently checked
against the same injective character map. Their cubic-algebra invariants
have a verified product square. Its square root determines an auxiliary
quadratic gamma, and the pairing is a product of local Hilbert symbols
\((a,\gamma(z_v))_v\) at points of the first quartic.

The exact retained values are
\[
\begin{aligned}
a&=112055369960752,\\
\gamma(z)&=8500097789885614889679
 +659438809664462289110z-220534595387434756313z^2.
\end{aligned}
\]
The product over the complete finite support of 22 primes is -1; the real
symbol is +1. **Only the contribution at 2 is negative.** At its local
witness z=4,
\[
\gamma(4)=7609299502344507945111,
\quad v_2(a)=4,
\quad a/16\equiv3\pmod8,
\quad\gamma(4)\equiv7\pmod8.
\]
Consequently the dyadic symbol is \((3,7)_2=-1\). This z is a coordinate
on the normalized quartic, not the original parameter t.

The carrier does have a Q2-point. The obstruction is global reciprocity:
a rational point would give Hilbert-symbol product +1, whereas the
Cassels–Tate evaluation gives -1. A list of local points therefore cannot
be assembled into a rational point here.

## What this changes in the rank-jump model

Together with the preceding positive controls, the comparison is now:

| Parameter carrier | Complete local solubility | Exact auxiliary rank | Global carrier solubility | Frozen cohort observation |
|---|---|---:|---|---|
| 1795d, 11278 | Yes | 3 | Yes, Sha[2]=0 | No simultaneous split in 32 fibres |
| 1795d, 0911e | Yes | 2 | Yes, Sha[2]=0 | Simultaneous split on observed +8 fibre |
| 030cb, 11278 | Yes | 2 | **No, labelled nonzero Sha class** | No simultaneous split in 32 fibres |

All three carriers arise from the existing generic bisection atlas, not
point-fitted constructions. The third row proves a stronger fact than a
bounded miss: this pair is forbidden at every rational specialization.
The first row also proves that absence from the sample need not mean a
global obstruction. Auxiliary rank two occurs in both the soluble and
insoluble rows; auxiliary rank alone does not decide the issue.

Ranked implications for continuing the investigation:

1. **Solubility, strongest established gate:** a class-specific reciprocity
   obstruction can prevent simultaneous rational lifting even when every
   local condition passes and each individual cover is rationally soluble.
   A positive mechanism must trivialize the relevant torsor classes, not
   merely exhibit a large Selmer space or common descent geometry.
2. **Solubility, still missing:** after a parameter carrier is proved
   globally soluble, the prescribed parameter must belong to its rational
   image. The first two rows show why global nonemptiness alone does not
   distinguish the sampled high-gain fibre. Conditions synchronizing several
   such images remain to be found.
3. **Incidence, still missing:** the singleton constructions yield two
   generic directions after base change, but this does not certify extra
   product-character sections or enough surviving specialized quotient
   directions to explain an extreme jump. These independence implications
   must accompany any proposed solubility mechanism.
4. **Weak explanations:** independent local pass/fail rates, auxiliary
   Jacobian rank, and nontrivial Sha dimension without identifying the
   actual class. None separates all three rows.
5. **Visibility:** no chart or point-height conclusion follows. This is a
   proof about existence, not exposure by the current point search.

The small falsifiable experiment fixed A,D, required an injective full-Selmer
character map, and tested whether its labelled pairing row vanished. The
observed nonzero row rejects global solubility for this carrier. It does
not yet test a prospective high-rank selector. Agent 1 could eventually use
certified forbidden cover combinations to eliminate impossible *construction
routes* in a predeclared atlas. They cannot exclude high rank of an original
fibre: other covers and directions remain available. No selection policy
was changed.

## Reproducibility and proof boundary

```sh
python3 elliptic-curves/rank-jump/verify_labelled_carrier_ct.py check
sage -python elliptic-curves/rank-jump/carrier_ct_local_witness.py check
```

The first command uses standard Python to replay model transport, class
coordinates, norm-square identities, the cubic square root and quadratic
formula, exact local square evaluations, Hilbert symbols, and recursive
primality certificates. The second cross-checks the local arithmetic with
Sage/PARI. Both read frozen inputs and run no descent or point search.

The full Selmer basis and completeness of the Fisher place support inherit
the pinned PARI computation. The verifier does not claim to reimplement
full descent. Source URL and hashes, compiler invocation, linked runtime,
observational patches and transcripts are retained. The separate executable
adds observations to the upstream algorithm and namespaces its exported
functions; it changes neither arithmetic nor the installed runtime.

- [Labelled-pairing protocol](LABELLED_CARRIER_CT_PROTOCOL.json), [local-witness protocol](CARRIER_CT_LOCAL_WITNESS_PROTOCOL.json)
- [Labelled basis and matrix inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_labelled_carrier_ct_inputs_v2.json)
- [Carrier class certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_labelled_carrier_ct_v2.json)
- [Local-formula inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_carrier_ct_local_inputs_v1.json), [local Hilbert certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_carrier_ct_local_witness_v1.json)
- [Primality certificates](../../artifacts/generated-results/elliptic-curves/rank_jump_carrier_ct_primality_v1.json), [portable verification](../../artifacts/generated-results/elliptic-curves/rank_jump_labelled_carrier_ct_verification_v1.json)

The initial capture and failed coefficient-conversion attempt remain in
`rank_jump_labelled_carrier_ct_inputs_v1.json` and
`rank_jump_labelled_carrier_ct_failure_v1.json`; the failure artifact retains
that producer's source. Corrected v2 artifacts were created separately.
All bounds above concern fixed auxiliary Jacobians. Original-fibre full
ranks remain UNKNOWN; their certified witness subgroups are unchanged.
