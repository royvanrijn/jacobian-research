# From Euler's restricted-prime sums to Eisenstein norms

**Unpublished · draft.**

The paper follows a mathematical exploration: the released Gaussian theorem, its Eisenstein variation, the mixed case forced by the prime3, and the primitive-norm and higher-sumset consequences. The complete statements and proofs remain precise; the surrounding prose explains why each extension is worth trying.

For `A = {a >= 1: every prime divisor of a is 1 mod 3}`, including `1`, the paper proves both exact binary identities: `A+A` is the positive `2 mod 6` class, and `A+(3·A)` is the positive `4 mod 6` class. Dilation `3·A` is distinct from the threefold sumset `3A`.

Primitive Eisenstein norm values are `N = A ⊔ (3·A)`. Their binary sumset consists of the positive even integers outside `0,12 mod 18`; for `k >= 3`, `kN = {n >= k: n = k mod 2}`. Four primitive norm values suffice for every positive integer, and four is sharp. For sums from `A`, the corresponding sharp bound is seven.

[PDF](paper.pdf) · [LaTeX](main.tex) · [Metadata](metadata.json) · [Source bundle](source.zip) · [Verification and review](verification.json)

The proof counts quaternion norm representations, retains only odd component norms, and uses an all-ideal Hecke estimate for odd inert primes. The mixed case adds an exact six-point calculation at the ramified prime `3`; its formula requires the reduced-target hypothesis. The total spoiled proportion is strictly less than `591/625` for reduced targets at least `10^8`. The standalone program printed in the PDF checks every remaining target in both families.

The abstract orbit estimate specifies a subgroup orbit `y₀K`, positive degree and the exceptional-space orthogonality. Its arithmetic criterion requires uniform local estimates, summable main terms and a strict exceptional-prime margin. It gives a parameterized split and an explicit cutoff when the error has an effective convergence modulus. A CRT proposition proves that no fixed finite small-summand menu works globally. It does not automatically prove a theorem for another norm form.

Build from this directory:

```sh
latexmk -pdf -jobname=paper -interaction=nonstopmode -halt-on-error main.tex
```

The canonical printed finite verifier needs only Python 3:

```sh
mkdir -p ../../.scratch/eisenstein-replay
python3 verify_binary_certificates.py > ../../.scratch/eisenstein-replay/binary_certificates_receipt.json
```

It covers all `16,666,667` unrestricted and `16,666,666` mixed targets below `10^8`, with no exceptions. The least small summands are at most `211` and `199`, respectively. The explicit failure checks remain active under `python3 -O`; direct-divisor and corrupted-certificate controls are recorded in the verification record.

The separate NumPy implementations generate and recheck all witnesses:

```sh
python3 compute_finite_certificate.py --limit 99999999 --output ../../.scratch/eisenstein-replay/unrestricted
python3 compute_mixed_certificate.py --limit 99999999 --output ../../.scratch/eisenstein-replay/mixed
python3 verify_mixed_constants.py > ../../.scratch/eisenstein-replay/mixed_constants_receipt.json
python3 check_ramified_counts.py --limit 6000 --output ../../.scratch/eisenstein-replay/ramified_counts_receipt.json
```

The ramified-count program is a bounded control of the written all-index proof. It checks the four count identities on 444 reduced even targets through 6000, the complete modulo-3 quotient and a wrong-sign control. These finite checks do not establish an unrestricted identity by themselves. The mixed constants program uses exact fractions. Original supplemental theta and finite-group eigenvalue diagnostics retain their original scope; numerical diagnostics do not prove the spectral theorem.

The top-level [receipts](receipts/) and [research notes](research_notes/) are preserved from the supplied source package. Fresh receipts are in [local replay](receipts/local-replay/), with source hashes and controls in the verification record. The first re-upload duplicated the earlier package. The subsequent extended package is preserved exactly as [supplied_extension.zip](supplied_extension.zip), with its Java/Python verifiers and supporting calculations integrated here. The supplied mixed Python program is in `supplied/compute_mixed_certificate.py`; the canonical verifier and robust NumPy witness checker remain separate. All supplied mixed witness indices agree exactly with the fresh calculation. The archived mixed Java receipt and source hash also agree, but no Java runtime was available for a fresh Java run. Current imported receipts are in [supplied extension receipts](receipts/supplied-extension/).

Additional supplied verification sources can be replayed without overwriting their imported receipts:

```sh
python3 supplied/compute_mixed_certificate.py --limit 99999999 --output ../../.scratch/eisenstein-replay/supplied-mixed
python3 compute_mixed_ramified_checks.py --limit 3000 --output ../../.scratch/eisenstein-replay/supplied-ramified.json
cp compute_mixed_odd_checks.py compute_mixed_wheel_checks.py mixed_spectral_check.py ../../.scratch/eisenstein-replay/
python3 ../../.scratch/eisenstein-replay/compute_mixed_odd_checks.py
python3 ../../.scratch/eisenstein-replay/compute_mixed_wheel_checks.py
python3 ../../.scratch/eisenstein-replay/mixed_spectral_check.py
java -Xmx512m VerifyMixedEisenstein.java 100000000 ../../.scratch/eisenstein-replay/mixed-java.json
```

The last command requires an existing Java 21 runtime and was not run here. The Python bound is inclusive; the Java bound is exclusive. The supplied checker sources have some `assert` checks and were run normally. The canonical printed verifier uses explicit exceptions and has passed invalid-input controls under `-O`. Exact wheel checks establish the all-real discrepancy bounds; exact modulo-3 checks support the six-point proof. They do not replace the all-index proof or certify imported spectral theorems.

The current source is split by topic. `main.tex` holds the preamble, abstract and section order. `sections/mixed-local.tex` contains the odd filter and six-point proof; `sections/primitive-sumsets.tex` contains the algebraic classification and exact sumsets; `sections/orbit-criterion.tex` contains the reusable theorem. Other sections retain the common spectral input, elementary estimates and finite completion. The appendix imports the canonical verifier directly so its printed and distributed code stay identical.

Extend the relevant section, keep theorem scope and metadata aligned, and update the existing verification record rather than creating another current edition. To refresh the source archive:

```sh
python3 -m zipfile -c source.zip main.tex references.tex sections \
  verify_binary_certificates.py verify_certificate.py VerifyEisenstein.java VerifyMixedEisenstein.java \
  compute_finite_certificate.py compute_mixed_certificate.py \
  compute_elementary_constants.py verify_mixed_constants.py check_ramified_counts.py \
  compute_theta_checks.py spectral_check.py compute_mixed_odd_checks.py \
  compute_mixed_ramified_checks.py compute_mixed_wheel_checks.py mixed_spectral_check.py \
  README.md metadata.json verification.json receipts research_notes supplied supplied_extension.zip
```

Keep caches, build auxiliaries and regenerated witness arrays out of the archive; include any new essential proof inputs explicitly. `verify_certificate.py` is the retained earlier unrestricted verifier.

Next: external specialist review of the spectral and ramified calculations, a complete novelty audit for primitive norm sumsets and basis bounds, external replay, and a persistent source-archive identifier before submission. Internal source and PDF-only review are recorded separately from external review. No public deposit or formal verification is claimed. The current PDF contains the finite balanced experiment. New canonical research
now proves [eventual balanced representations and fixed-interval angular distribution](../../research/arithmetic/EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md)
for both binary families, using this manuscript's odd-prime Hecke lemma and
fixed-level harmonic theta series. That qualitative extension has received
internal AI proof/source checks and is not yet incorporated into the PDF. It
gives no effective onset, so identifying the complete finite balanced exception
sets remains open, as do external review and historical priority.
