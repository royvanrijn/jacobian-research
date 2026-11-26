# Fourteen-row coherent detection at integral places

For the conductor-18 formal response with the two opposite fourth-order
roots, \(\deg B\leq5\), phase \(D=-z^{-1}+c+ez\), and correction degree at
most two, shifts \(0,\ldots,13\) detect every nonzero coherent input if
the algebraic coefficients \(c,e\) are integral at a common place above
one of the primes
\[
\mathcal P=\{331,337,347,349,353,359,367,401,65521\}.
\]
The coefficients of \(B\) may be arbitrary algebraic numbers.

Use the [three-pole reduction](CM_NEWTON_FULL_WIDTH_STRIP_MOMENTS.md#three-extra-pole-inputs-for-odd-degree-five-and-even-degree-one).
The branch with \(b=0\) is detected in shifts \(0,1\). If \(b\ne0\),
failure through shift 13 gives a common zero of the seven primitive
residual polynomials \(S_7,\ldots,S_{13}\) in the normalized coordinates
\((r,\ell,h,v,w)=(b^2,b_1/b,b_2/b,c,e)\). We use the same rows as in the
[seventeen-row proof](CM_NEWTON_SEVENTEEN_ROW_LOCUS.md).

For each \(p\in\mathcal P\), two exact Gröbner basis computations
over \(\mathbb F_p\) verify the following.

* The ideal of \(S_7,\ldots,S_{13}\) is the unit ideal in
  \(\mathbb F_p[r,\ell,h,v,w]\).
* With weights \((2,1,2)\) on \((r,\ell,h)\), the coefficient-free
  leading forms of \(S_7,S_8,S_9\) have zero-dimensional homogeneous
  quotient of length 84. At \(p=359\) use the first four forms
  \(S_7,S_8,S_9,S_{10}\), whose quotient has length 68. The leading
  forms are obtained by setting \(v=w=0\) and retaining weight \(k-1\)
  in \(S_k\). Since each such ideal is homogeneous for positive
  weights and has finite-dimensional quotient, its only geometric
  zero is the origin.

Suppose an algebraic common zero has \(v,w\) integral at a place over
\(p\in\mathcal P\). Let \(\nu\) be the corresponding additive valuation
and take the minimum of \(\nu(r)/2,\nu(\ell),\nu(h)/2\), with
\(\nu(0)=+\infty\). If this minimum is negative, in a finite ramified
extension choose a scaling parameter with valuation equal to the
minimum. Divide \(r,h\) by its square and \(\ell\) by the parameter, and
divide each \(S_k\) by the parameter to the power \(k-1\). Each term of
smaller weight in \(r,\ell,h\) then vanishes on reduction because \(v,w\)
are integral. The scaled triple reduces to a nonzero common zero of the
leading forms, a contradiction. Thus all five coordinates are integral
and their residues would be a common zero of the seven rows. The first
computation excludes that zero over every residue field extension.

This affine computation is also used at \(331\) in the
[universal fourteen-row proof](CM_NEWTON_FOURTEEN_ROW_LOCUS.md).

## Replay

    sage -python research/foundations/extended-geometry/verify_cm_newton_fourteen_row_integral_places.py
