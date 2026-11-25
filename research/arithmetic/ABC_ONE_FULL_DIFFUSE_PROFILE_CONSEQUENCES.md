# Low degrees and narrow residue bands at the one-full \(abc\) endpoint

At a limiting box with normalized count exponent \(d=1/3\),
suppose one letter has weak scaled-degree measure supported on
\([3,\infty)\). Its measure is \(\frac13\delta_3\). For each
of the other two letters let \(\nu\) be its diffuse measure.
The one-full density estimate gives absolute constants \(R,C_0>0\)
such that \(\nu\) has mass \(1/3\), support in \([0,R]\),
density at most \(C_0\), and first moment at most one.
[The one-full density estimate](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md#uniformly-bounded-diffuse-density-at-the-one-full-endpoint)

## Corollary. Quantitative low-degree mass on both diffuse letters

For every fixed \(s>0\), each of the two diffuse letters satisfies
\[
 \nu([0,s))\ \ge\
 \max\left(0,\ \frac13-C_0\bigl(\sqrt{s^2+2/C_0}-s\bigr)\right).
\]
In particular the lower bound is positive for
\(s>3-1/(6C_0)\), including \(s=3\).

Let \(x=\nu([s,\infty))\). At distance \(u\ge0\) to the right
of \(s\), at least \((x-C_0u)_+\) of this mass remains, by the
density bound. Integrating in \(u\) shows
\(\int_{[s,\infty)} r\,d\nu(r)\ge sx+x^2/(2C_0)\).
The first moment is at most one. Solving this quadratic upper
bound on \(x\) proves the assertion.

## Corollary. A saving with a near-full letter and narrow residue bands

Fix a positive integer \(K\) and constants \(\kappa,b>0\).
There are positive \(\eta,\tau,\delta,\lambda_*\) with the
following property. For fixed \(0<\lambda\le\lambda_*\), count
ordered positive primitive triples of height at most \(X\) with
\(\operatorname{rad}(abc)<c^\lambda\). Require a prescribed
letter to be \(p\)-full with \(p\lambda\ge3-\tau\). On a different
prescribed letter, fix a positive integer \(q\) with
\(q\lambda\ge\kappa\) and at most \(K\) residues modulo \(q\).
Capture primes whose positive valuations \(e\) satisfy
\[
 \operatorname{dist}\left(\lambda e,
                 \bigcup_{r\in\mathcal R}\lambda(r+q\mathbb Z)\right)
 \le\eta .
\]
If the product of captured primes has logarithm at least
\(b\lambda\log c\), the count is
\(O_{\lambda,p,q,\mathcal R,K,\kappa,b,\varepsilon}
 (X^{(1/3-\delta)\lambda+\varepsilon})\).
The width \(\eta\) is in scaled degree and is independent of the
height and small budget.

In a limiting endpoint as above, every union of \(\eta\)-neighborhoods
of \(K\) progressions of spacing at least \(\kappa\) has
\(\nu\)-mass at most
\[
 2 C_0\eta K\left(1+\frac{R+2\eta}{\kappa}\right).
 \tag{NR}
\]
Only the centers in \([-\eta,R+\eta]\) can meet the support.
The same bound covers a union with some of the progressions
replaced by singleton centers.

Choose \(0<\eta<1\) so the right side of (NR) is strictly below
\(b\). If no \(\tau,\delta,\lambda_*>0\) works, choose failing
fixed-parameter counts with \(\lambda_j\to0\), \(p_j\lambda_j\ge3-o(1)\),
and exponent at least \((1/3-o(1))\lambda_j\). The large-step
box argument gives \(L_j/\lambda_j\to1\) and a \(d=1/3\) profile.
The prescribed full letter's measure is supported on \([3,\infty)\).
Choose the anatomic compression index with \(k_jL_j\to\infty\)
and \(k_j\ge p_j\). The captured submeasure on exact, nonpacked
coordinates has mass at least \(b-o(1)\) by the same packed
root and truncation estimates. [The large-step box argument](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md#proof-of-the-endpoint-transfer)

In scaled degree \(eL_j\) the captured coordinates lie within
\(\eta L_j/\lambda_j\) of the tested progressions, whose spacing
is \(q_jL_j\ge\kappa-o(1)\). Along a further subsequence the spacing
has a finite limit at least \(\kappa\), or tends to infinity.
In the finite case, pass to limits of the finitely many phase
offsets. In the infinite case, retain at most one limiting center
per phase in any fixed bounded interval. The captured weak
submeasure is supported on the resulting closed bands where the
full weak measure lives. This support assertion follows on each
compact interval by convergence of the relevant centers, and the
first-moment bound gives tightness. The captured weak submeasure
therefore has mass at least \(b\), contradicting (NR).
