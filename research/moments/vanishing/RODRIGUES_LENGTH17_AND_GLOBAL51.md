# Actual Rodrigues local algebra and the global length51 gate

17 November 2026. The user supplied
`rodrigues_local_length17_certificate.zip`, including the proof, exact
rational Taylor scripts, finite relation matrices, actual multiplication
matrices and verification log. The archive is retained unchanged with SHA256

```text
920fc19dc1a43f89b7de809d0b3a3c1098e3472b973649e1d392561c8b291f87
```

**SIC2B33LOCALTHICK20261117.** In the full harmonic coefficient space,
transverse to the scalar/rotation Rodrigues orbit,

\[
 d_{12}=22,\qquad d_{13}=18,\qquad d_N=d_\infty=17\quad(N\ge14).
\]

The all-moment ordinary local ideal has minimal consecutive generating
prefix14. The full transverse algebra has nilpotence index7, ordinary
Hilbert function \((1,2,3,4,3,3,1)\), and the supplied actual17-by17
multiplication matrices. There is no msolve dependency in this local result.

The normalized global cyclic target is now exactly

\[
 \boxed{\dim_{\mathbb C}A_{14}\stackrel{?}=51.}
\]

Global finiteness and length remain uncomputed. The local thickness,
actual multiplication and infinite-tail questions are closed by this
certificate; the next computation concerns global support and finiteness.

The subsequent [seven-target certificate](SEVEN_LOW_DEGREE_GLOBAL_GATE.md)
reduces the preferred global test to seven memberships of degree at most5,
using an exact augmented first-nine ideal identity. The earlier fourteen
seventh-power formulation below stays valid and is retained for provenance.

## Retained certificate and replay

- [Unchanged supplied proof](rodrigues_local_length17_certificate/PROOF.md).
- [Unchanged supplied scripts and conventions](rodrigues_local_length17_certificate/README.md).
- [Exact rational data and multiplication matrices](rodrigues_local_length17_certificate/certificate_data.json).
- [Original manifest](rodrigues_local_length17_certificate/manifest.json).
- [Repository audit and replay entry point](verify_rodrigues_length17.py).
- [Replay log](rodrigues_length17_replay.log) and [bounded run record](rodrigues_length17_replay_run.json).
- [Additional exact row witnesses](rodrigues_length17_row_witnesses.json).
- [Verification receipt](rodrigues_length17_verification.json).

All manifest hashes match. Replaying the supplied verifier in a separate
workspace directory regenerated `certificate_data.json` byte-for-byte:
SHA256 `bc17438be6954d88dcf0569f94cdf7fbe5cb7afa616cce16b78f27cfbe709194`.
The replay took about13 seconds with180-second wall/CPU and2GiB address-space
bounds. The original uploaded files and original log remain unchanged.

The separate repository audit reconstructs both rational relation matrices
from residual coefficient dictionaries and checks the exported actual
multiplication operators. It also exports explicit32-coordinate row
witnesses for the six generators of \(K\) and for \(R_{15},\ldots,R_{24}\).
These are direct identities against the supplied32-by36 matrix, not just
an assertion that an augmented rank was unchanged.

## Why the finite calculation is the full local algebra

Use the certificate's explicit transverse slice and its parameter conventions.
Its first-ten Jacobian is invertible, and its twelve slice directions plus
four orbit directions have determinant \(-3/50\) in the full harmonic basis.
Solving the first ten moments therefore gives actual formal residuals
\(R_n(a,b)\). The two parameters are normalized by

\[
 Q_a=w^2-\tfrac1{10}(t^3-3zwt)+\tfrac{17}{70}t,
 \qquad Q_b=wt+\tfrac25z(t^2-zw/2)-\tfrac15z.
\]

Use these coordinates for the new matrices; no multiplication in the
earlier weighted initial quotient is substituted for them.

The finite ambient vector space is
\(V=\mathbb Q[a,b]/(a,b)^8\), of dimension36. The exact prefix12,13,14
relation ranks are16,18,19, giving truncation lengths20,18,17.
The prefix12 value20 is a truncation, not the retained full length22.
For prefix13 the weighted initial upper bound18 matches the truncation
lower bound18, proving \(d_{13}=18\).

At prefix14 let

\[
 K=(a^5,a^4b,a^3b^2,a^2b^4,ab^5,b^7).
\]

This is exactly the ideal of monomials of weight at least13 for
\(\operatorname{wt}(a,b)=(3,2)\). The six exact row memberships give
\(KL_{14}\subseteq\mathfrak m^8L_{14}\), while
\(\mathfrak m^7\subseteq K\) gives
\(\mathfrak m^8L_{14}\subseteq\mathfrak m KL_{14}\).
The finitely generated ideal \(KL_{14}\) is therefore equal to its
maximal-ideal multiple. [Nakayama's lemma](https://stacks.math.columbia.edu/tag/07RC)
gives \(KL_{14}=0\). Thus the computed17-dimensional truncation is the
actual completed local algebra, with both \(F^{13}L_{14}=0\) and
\(\mathfrak m^7L_{14}=0\). The checked nonzero class \(b^6\) makes the
nilpotence index exactly7.

The actual multiplication matrices commute, have a cyclic unit and the
given17-element basis, and have maximal-ideal power dimensions

\[
 (17,16,14,11,7,4,1,0).
\]

This proves the stated ordinary Hilbert function, without identifying
actual multiplication with a weighted initial multiplication table.

For the infinite tail, the inverse Rodrigues shear sends \(F_0\) to
\(w\). The audit independently checks its maximum torus weights:
\(4,3\) for \(Q_a,Q_b\), and \(9,8,7\) for the three actual quadratic
implicit-correction markers. It verifies the quadratic correction in
the first ten moment equations. Every later correction marker still
lies in the original degree-at-most-three space, whose pulled-back
maximum weight is at most9. Consequently

\[
 [a^ib^j]R_n\ne0\ \Longrightarrow\ 5i+4j\ge n.
\]

Since \(5i+4j\le2(3i+2j)\) and \(F^{13}L_{14}=0\), moments with
\(n>24\) vanish automatically in \(L_{14}\). Exact row witnesses verify
the remaining \(R_{15},\ldots,R_{24}\); the supplied replay also checks
through28. Hence \(J_\infty=J_{14}\) in the completion. Equality descends
to the ordinary local ring by
[faithful flatness of completion](https://stacks.math.columbia.edu/tag/00MC).
The strict drop \(18\to17\) proves minimality of the consecutive prefix14.

## The exact global certificate to pursue

At each of the three normalized Rodrigues points the local contribution
is now17. The scalar direction is transverse to \(D=1/8\), and the cyclic
frame removes rotations without a stabilizer. A normalization root near
the identity identifies the supplied transverse slice with this fixed-D
quotient locally. The relative sign acts as an automorphism there; it is
not an additional quotient in \(A_N\). There is no division by two.

Thus global length51 is equivalent to \(A_{14}\) having precisely these
three points as its support, and supplies the desired cyclic-open pure
exhaustion certificate. It is not yet known whether other finite-prefix
points or positive-dimensional pieces survive.

The settled nilpotence index sharpens the earlier fourteen-membership
formulation. Keep exactly the same linear forms

\[
 \mathcal L=\{p,r,s_0,s_1\}
 \cup\{U_{ijk}:(i,j,k)\ne(2,2,2)\}
 \cup\{U_{222}+128s_2\}.
\]

Let \(I_{14}\) be the four structural relations and
\(\mu_2,\ldots,\mu_{14}\) in the complete normalized quotient. Then

\[
 \boxed{\dim_{\mathbb C}A_{14}=51
 \quad\Longleftrightarrow\quad
 \ell^7\in I_{14}\ \hbox{for every }\ell\in\mathcal L.}
\]

The reverse direction follows because setting \(\mathcal L=0\) leaves
exactly \(\mathbb C[s_2]/(s_2^3-1/64)\), so nilpotence forces finite support
at the three known points. Their actual local lengths then sum to51.
The forward direction uses the actual maximal-ideal power7 at all three
locals. This corollary asserts equivalence, not any successful membership.

The monic Gram relation still gives the rank-three module formulation.
After removing it there are16 generators: three traces and thirteen
moments. Multiplying each by \(1,s_2,s_2^2\) gives48 module columns.
The [exact module program](cyclic_prefix14_nilpotence_module_program.json)
uses seventh-power targets and the retained full moment circuit through14.
Both \(\mu_{13}\) and \(\mu_{14}\) remain in the global input.

The known anti-halfturn first12 scheme is already removed by its inherited
\(\mu_{14}\)-unit theorem. Its global lower bound13818 for \(A_{12}\)
retains its existing msolve dependency; this local length17 calculation
adds none. No new all-order pure component, full global finiteness result
or full(3,3) Image theorem is claimed.

## Reproduction

```sh
python verify_rodrigues_length17.py
python verify_rodrigues_length17.py --replay
```

The first checks retained provenance, exact finite arithmetic, the
shear bounds and the new module program. The second first reruns the
supplied Taylor scripts in a separate workspace directory, preserving
the uploaded bundle. No new local jet campaign is needed.
