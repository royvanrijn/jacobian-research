# Seven low-degree targets for the global cyclic A14 gate

The [field-valued strategy](FIELD_COVER_AND_GENERIC_COMPONENT_STRATEGY.md)
also allows later prefixes, arbitrary radical exponents and additional
generically safe pure components. This exact length51 gate is an optional
stronger route, not a requirement for the full field-valued Image theorem.

17 November 2026. The user supplied the complete seven-target certificate,
including the global selected-locus proof, exact normalized local coordinate
adapter, rational Bezout data, both gates and deterministic replay scripts.
The archive and all supplied files are retained unchanged.

**SIC2B33SEVENGATE20261117.** In the full cyclic model normalized by
\(D=1/8\), the following seven memberships are equivalent to global length51:

\[
 \boxed{K_1^3,\ \Psi,\ s_1^5,\ U_{000}^3,\ K_5^2,\ K_6^3,\ K_7^2
             \ \in I_{14}.}
\]

Their raw polynomial degrees are \((3,3,5,3,2,3,4)\). These are now the
preferred global targets. None of their unrestricted \(I_{14}\) memberships
has been proved. Witness degrees, global finiteness and length51 remain open.
Both \(\mu_{13}\) and \(\mu_{14}\) remain in the ideal.

## Targets and local identities

Use

\[
\begin{aligned}
K_1&=4p+5s_0,&K_2&=s_0,&K_3&=s_1,\\
K_4&=U_{000},&K_5&=U_{002}-6U_{011},\\
K_6&=4U_{022}-13U_{112},&
K_7&=U_{222}+128s_2+\frac{896}{75}p^2.
\end{aligned}
\]

The exact least nilpotence exponents in the actual normalized local
length17 algebra are \((3,7,5,3,2,3,2)\), with preceding powers nonzero.
The supplied adapter checks the entire cyclic frame and all ten symmetric
cubic coordinates, including the amplitude normalization. It does not use
the weighted initial multiplication table.
In its JSON, the key `U222` denotes the nilpotent target
\(T=U_{222}+128s_2\); raw \(U_{222}\) is recovered by subtracting \(128s_2\).

The integer cubic is

\[
\begin{aligned}
\Psi={}&3814472897542065s_0^3+9204109983515700s_1^2\\
&-9625759731192528s_2K_1
-10061220068131439K_1s_0^2-5712042093615K_6.
\end{aligned}
\]

The five displayed local vectors have rank4 and this primitive integer
dependence vanishes exactly. The monic rational version in the supplied
proof differs by the nonzero scalar3814472897542065, a unit over \(\mathbb Q\).

The original seven-power gate
\(K_1^3,s_0^7,s_1^5,U_{000}^3,K_5^2,K_6^3,K_7^2\in I_{14}\)
is equivalent as well. The cubic version removes its degree-seven target.
This says nothing about the necessary degrees of witness multipliers.

## The new global support and scheme identity

Let \(P\) be the ideal of the three reduced normalized Rodrigues points,
and put

\[
 \mathcal Z=(p,s_0,s_1,U_{000},U_{002},U_{022},U_{222}+128s_2).
\]

The new global conclusion is the augmented identity

\[
 \boxed{I_9+(\mathcal Z)=I_9+(K_1,\ldots,K_7)=P.}
\]

Already the first six moments have exactly those three points as support
**after adjoining the seven selected equations**. This is not a first-six
or first-nine exhaustion statement for the unrestricted cyclic quotient.

On the selected locus the structural relations force
\(s_2^3=1/64\), \(U_{011}=U_{112}=0\) and \(U_{111}=-2U_{012}\).
After transporting the entire oriented frame on the branch \(s_2=1/4\),
the full remaining phase is

\[
 F=F_0+Rw^2+Atw^2+\beta(t^3-3ztw)+Cz^2t.
\]

Here \(R=4r\), \(A=3U_{001}/2\), \(\beta=-U_{012}/4\),
\(C=3U_{122}/32\), and the coefficient map has determinant \(-9/64\).
No value of \(r\) or remaining cubic coefficient is discarded.

The cases \(\beta=0\) or \(C=0\) reduce to the origin. On the complement,
the first moments solve for parameters using units
\(\beta,C,2\beta C+3\). A further zero-denominator case is excluded by
an exact rational Bezout identity. The remaining fifth/sixth equations
give coprime degree8/degree10 polynomials in \(\beta C\), again with an
explicit rational Bezout identity. This proves the first-six support result.
The supplied replay separately compares these moment polynomials with
direct spatial integration, so no modular global-length inference enters.

The Jacobian of the primitive integer-scaled odd moments3,5,7,9 has
determinant49411883664. This is the normalization used in the archive.
For the raw normalized spherical moments the determinant is
\(2048/676350675\), also nonzero. The repository audit checks the row
scalings explicitly. Thus the first-nine selected intersection is reduced;
the three simple roots of \(s_2^3-1/64\) give exactly \(P\).

The equality \((K_1,\ldots,K_7)=(\mathcal Z)\) is scheme-theoretic in the
ambient structural ring. Modulo \((p,s_0,s_1)\), the Gram relation makes
\(s_2\) a unit, and the traces give \(U_{011}=-2U_{002}\),
\(U_{112}=-2U_{022}\). Then \(K_5=13U_{002}\), \(K_6=30U_{022}\),
\(K_7=U_{222}+128s_2\). The checker additionally verifies equality of the
two ambient ideals by exact rational Groebner reduction of the structural
relations and seven generators. This small auxiliary calculation does not
compute the unrestricted moment ideal.

## Why the seven tests are equivalent to length51

Assume all seven cubic-version memberships. At any complex closed point
of \(A_{14}\), the six power tests force
\(K_1=s_1=U_{000}=K_5=K_6=K_7=0\). There \(\Psi\) reduces to the nonzero
constant3814472897542065 times \(s_0^3\); its vanishing forces \(s_0=0\).
The augmented support theorem leaves precisely the three Rodrigues points.
Finite support makes the quotient Artinian. Its actual local lengths17
then give51.

Conversely, length51 is exhausted by those three known local algebras.
All seven targets vanish there by actual multiplication. The cyclic
grading below transports their homogeneous identities among the three
points. Consequently all seven polynomials belong to \(I_{14}\).

## The graded module exporter

The full normalized ring has the rational \(\mathbb Z/12\mathbb Z\) grading

\[
\deg p=2,\quad \deg r=3,\quad \deg s_j=j+2,\quad
\deg U_{ijk}=-i-j-k-2\pmod{12}.
\]

It is induced by \(F(x)\mapsto\lambda F(\lambda^6x)\),
\(\lambda^{12}=1\). The oriented determinant stays fixed; moment \(\mu_n\)
has weight \(n\). The Gram equation has weight0 and the traces have
weights0,11,10. The audit checks the grading gate-for-gate in all1341228
nodes of the retained first14 moment circuit.

The seven target weights are \((6,6,3,6,4,6,8)\), and the monic Gram
basis \((1,s_2,s_2^2)\) has weights \((0,4,8)\). No extra scalar or
orthogonal normalization is introduced.

After removing the Gram equation, the full \(I_{14}\) has16 source
generators: three traces and thirteen moments. Their three Gram-basis
multiples give48 **source-generator recipes**. This count is independent
of the number of targets. Actual finite ansatz columns must still be
recomputed for the selected degree bound and target; no universal
factor-of-twelve reduction is asserted.

For a generator of weight \(v\), its \(s_2^k\) source column has weight
\(v+4k\). A target of weight \(w\) therefore needs multiplier monomials
only of weight \(w-v-4k\pmod{12}\). Its three remainder components have
weights \(w,w-4,w-8\). Projection of a witness identity onto its target
weight proves this filtering exact over \(\mathbb Q\).

The [machine-readable module program](cyclic_prefix14_seven_target_module_program.json)
exports all seven polynomials, exact remainder vectors, component weights
and per-source-column multiplier filters. It specifies no witness-degree
bound and records zero verified memberships. The repository polynomial
symbols are lowercase `uijk`; the supplied mathematical notation is
\(U_{ijk}\).

## Retained evidence and replay

- [Supplied full proof](rodrigues_global_seven_target_certificate/PROOF.md).
- [Supplied support and denominator-case proof](rodrigues_global_seven_target_certificate/odd_slice_support_certificate/ODD_SLICE_PROOF.md).
- [Supplied target definitions](rodrigues_global_seven_target_certificate/global_gate_targets.json).
- [Supplied adapter and normalized nilpotence convention](rodrigues_global_seven_target_certificate/GLOBAL_TARGET_NILPOTENCE_ADAPTER.md).
- [Supplied manifest](rodrigues_global_seven_target_certificate/manifest.json).
- [Repository model](seven_target_gate.py), [checker](verify_seven_target_gate.py)
  and [verification receipt](seven_target_gate_verification.json).
- [Replay log](rodrigues_seven_targets_replay.log) and
  [bounded run record](rodrigues_seven_targets_replay_run.json).

The preserved archive SHA256 is
`5ac529c2338c5404bace25091c2623e9a75710a5660455e1d58dc52b0361b020`.
All input hashes match, the included original local certificate is
unchanged, and all four regenerated JSON outputs match byte-for-byte.
The full replay took about13 seconds under180-second wall/CPU and3GiB
address-space bounds. The supplied packet is preserved while replay
writes to a separate workspace directory.

```sh
python verify_seven_target_gate.py
python verify_seven_target_gate.py --replay
```

The first performs the additional finite audits and exporter checks; the
second also reruns the supplied packet. These calculations have no msolve
dependency. The older separate global \(A_{12}\) lower bound retains its
own CAS dependency. No unrestricted seven-target membership, global
finiteness, length51, witness-degree bound, extra pure component or full
(3,3) Image theorem is proved here. Independent specialist and formal
verification are not claimed.
