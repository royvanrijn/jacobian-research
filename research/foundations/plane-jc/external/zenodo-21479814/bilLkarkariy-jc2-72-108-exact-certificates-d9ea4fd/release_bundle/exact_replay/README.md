# Exact replay workspace

The missing characteristic-zero certificates were independently reconstructed
on 2026-07-21 from the archived exact equation generators and the deterministic
fixed minor selected modulo 71.

The eight requested replay artifacts are now present:

```text
verify_all.sh
verify_serialized_certificates.py
verify_hne0_branch_symmetry.py
case2_exact_certificate.pkl
h0_branch1_exact_certificate.pkl
h0_exact_certificate.pkl
hard/h_certificate_exact.txt
hard/verify_certificate_gmpy2.py
```

Run the complete local replay with:

```bash
./verify_all.sh
```

The expected terminal marker is:

```text
JC2_72_108_EXACT_REPLAY_PASS
```

`RECONSTRUCTED_CERTIFICATES.sha256` records the hashes of these newly generated
artifacts. They are mathematically equivalent reconstructions, not byte-for-byte
copies of the historical files whose expected hashes appear in the old report;
their hashes therefore differ.

The main reconstruction programs are:

```text
reconstruct_case2_certificate.py
reconstruct_h0_certificates.py
hard/reconstruct_h_certificate.py
```

The replay proves the saved algebraic identities for the transcribed systems.
The claimed degree-125 consequence remains conditional on the correctness and
exhaustiveness of the published Laurent/Newton reduction and its transcription.
