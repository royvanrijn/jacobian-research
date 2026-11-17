# Endpoint obstruction for a nonzero-target collision family

For the foundational map \(F\) in [KFRAME1's proof](../sparse-nine/PROOF.md), put

\[
a(k)=((1-k)/2,-2,-4(1+3k)/k^2),\qquad
b(k)=((1+k)/2,-2,-4(1-3k)/k^2),\quad k\ne0.
\]

Both points map to \((0,-2,(1-k^2)/(2k^2))\). They are distinct for
\(k\ne0\). Let \(A_k=JF(a(k))\), \(L_k=[b(k)-a(k),e_y,e_z]\), and

\[
H_k(X)=(A_kL_k)^{-1}\bigl(F(a(k)+L_kX)-F(a(k))\bigr).
\]

Then \(\det A_k=-2\), \(\det L_k=k\), \(JH_k(0)=I\),
\(\det JH_k=1\), and \(H_k(0)=H_k(e_1)=0\). At \(k=1\) this is the
49-term frame from which the 48-term witness is obtained by a transverse
shear. Generic support in this natural frame is 105.

**Claim KFRAME3.** No invertible meromorphic matrix \(P(k)\) satisfying
\(P(k)e_1=e_1\) can make
\(P(k)^{-1}H_k(P(k)X)\) have a finite coefficientwise polynomial limit at
either \(k=0\) or \(k=\infty\). The statement also holds after finite
ramification of the local parameter. It concerns this specified pair
family; it does not exclude other pair families or finite-parameter
support drops.

## Two polynomial invariants prove the obstruction

For a polynomial map \(H\), define the scalar polynomial functions of its
coefficients

\[
I(H)=\operatorname{tr}JH(2e_1),\qquad
T(H)=\operatorname{tr}JH(H(2e_1)).
\]

For \(G=P^{-1}H(PX)\) and \(Pe_1=e_1\), the chain rule gives

\[
JG(2e_1)=P^{-1}JH(2e_1)P,\quad
JG(G(2e_1))=P^{-1}JH(H(2e_1))P.
\]

Thus both traces are invariant. In the fixed degree-at-most-seven
coefficient space, a finite coefficientwise limit makes these polynomial
functions finite as well.

Exact computation gives

\[
I(H_k)=864k^2-2736k+723.
\]

Its pole at infinity excludes that endpoint. At zero the second invariant
has Laurent expansion

\[
T(H_k)=-49460507836416\,k^{-3}+O(k^{-2}).
\]

This nonzero pole excludes zero. A finite ramification multiplies the pole
order by a positive integer and cannot remove it. This proves the claim.

The full Laurent polynomial for \(T\), rather than only its leading term,
is retained in the receipts below. No specialization at \(k=0\) of the
original collision formulas is used.

## Reproducible verification

- [Primary exact construction](nonzero_target_probe.py) checks the collision,
  tangent normalization, all 105 generic nonlinear coefficients and both traces.
- [Primary data](nonzero-target-probe.json) retains the full formulas and source hash.
- [Independent checker](verify_nonzero_target.py) uses only standard-library
  rational Laurent arithmetic, the original factored map and explicitly
  differentiated formulas. It computes the source point
  \(a+A^{-1}(F(a+2(b-a))-F(a))\) directly, evaluates
  \(\operatorname{tr}(A^{-1}JF)\) there, and checks the entire Laurent
  polynomial against the primary result.
- [Independent receipt](nonzero-target-validation.json) records the comparison
  and all Laurent coefficients.

From the repository root:

~~~sh
.venv/bin/python research/jacobian-followup/collision-frame-search/nonzero_target_probe.py
.venv/bin/python research/jacobian-followup/collision-frame-search/verify_nonzero_target.py
~~~

## Earlier restricted tests, retained

The primary data also retains diagonal-rescaling obstructions in the natural
frame and in a canonical quartic-axis frame. For the latter, if \(C\) is the
matrix of the three axis coefficient vectors for \(x^2,x^3,x^4\), then

\[
\det C=248832k^2,\quad
D=\begin{pmatrix}-1&0&0\\1&-1&0\\0&1&-1\end{pmatrix},\quad M=CD^{-1}.
\]

The matrix \(M\) fixes \(e_1\), and \(M^{-1}H_k(MX)\) has axis
\((x-x^2,x^2-x^3,x^3-x^4)\) and 312 generic nonlinear occurrences. Exact
valuation certificates exclude diagonal limits there; its leading
\(k^{-4}\) output tensor has rank one. These limited tests motivated the
invariant calculation. The trace obstruction above now excludes every
meromorphic collision frame at the two endpoints, so revisiting diagonal
scalings in this family cannot advance the search.

The support interval remains **[10,48]** in the normalized three-variable,
raw-degree-at-most-seven domain. This obstruction does not change that bound.

