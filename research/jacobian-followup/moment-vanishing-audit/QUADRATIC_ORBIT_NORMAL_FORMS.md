# Five complete quadratic normal forms for the full(3,3) field problem

17 November 2026. **SIC2B33QUADFORMS20261117.** After the first moment
removes the constant, every phase in the full field-valued problem can
be carried, by an oriented linear isometry and a nonzero phase amplitude
over an algebraic field extension, into one of the following five families.
In **every** row retain an arbitrary linear harmonic and all seven cubic
harmonics:

\[
 F=L+H_2+\sum_{j=-3}^{3}\gamma_jC_j,\qquad
 L=\ell_z z+\ell_t t+\ell_w w,
\]

with no coefficient restriction on \(L\) or the cubic.

| Quadratic operator type | Complete normalized \(H_2\) |
| --- | --- |
| Zero | \(0\) |
| Nilpotent rank one | \(w^2\) |
| Nilpotent rank two | \(zt\) |
| Non-nilpotent, nonsemisimple | \(t^2-zw+w^2\) |
| Nonzero semisimple | \(z^2+w^2+\eta(t^2-zw),\quad\eta\in\overline K\) |

**No value of \(\eta\) is omitted.** In particular the semisimple row
retains the repeated-eigenvalue values \(\eta=\pm2/3\), and rank-two
values \(\eta=0,\pm2\). The largest family has eleven phase parameters:
\(\eta\), three linear coefficients and seven cubic coefficients.

This is exhaustive **field-valued orbit coverage**, not a purity
classification, generic safety theorem, ordinary quotient identification
or global length certificate. It retains every \(D=0\) boundary point.

## Completeness proof

Work first over an algebraic closure of a characteristic-zero field.
For a self-adjoint operator, distinct generalized eigenspaces are
orthogonal: Bezout for their relatively prime primary polynomials moves
the annihilating factors between the two arguments of the pairing.
Each primary summand is nondegenerate, since their orthogonal direct sum
is the whole nondegenerate space.

If all eigenvalues coincide, trace zero makes them zero. The operator is
nilpotent and has rank zero, one or two. Rank one gives a square of an
isotropic linear form, normalized to \(w^2\). Rank two has a length-three
Jordan chain. The exact complete Jordan-chain Gram normalization in the
[nilpotent boundary proof](NILPOTENT_QUADRATIC_ZERO_LINEAR_BOUNDARY.md)
puts it in the form \(zt\). These arguments concern the quadratic
operator alone and impose no condition on the other harmonics.

Suppose the operator is semisimple and nonzero. Its eigenvalues are not
all equal. Choose two distinct ones \(\alpha,\gamma\), and set
\(r=(\alpha-\gamma)/2\ne0\), \(c=(\alpha+\gamma)/2\).
An orthogonal eigenbasis, followed by an isotropic basis on the plane of
those two eigenlines, puts its matrix in the form

\[
 \begin{pmatrix}c&0&r\\0&-2c&0\\r&0&c\end{pmatrix}.
\]

Scale the entire phase by \(1/r\) and put \(\eta=-2c/r\).
This gives

\[
 B_S(\eta)=
 \begin{pmatrix}-\eta/2&0&1\\0&\eta&0\\1&0&-\eta/2\end{pmatrix},
\]

whose quadratic is the last row of the table. The eigenvalues are
\(1-\eta/2,\eta,-1-\eta/2\); the eigenvector matrix is invertible for
every \(\eta\). Thus this normalization keeps the repeated-eigenvalue
cases as well as the distinct ones. It divides only by a chosen nonzero
eigenvalue gap.

Finally suppose the operator is non-nilpotent and nonsemisimple. Its
spectrum is \((\lambda,\lambda,-2\lambda)\), \(\lambda\ne0\), with a
length-two Jordan block at the repeated eigenvalue. On that nondegenerate
two-plane the nilpotent part is self-adjoint, rank one and square-zero,
so it is an isotropic square. Normalize the semisimple amplitude to
\(t^2-zw\); the remaining orthogonal plane torus normalizes the nonzero
isotropic-square coefficient to one. This gives

\[
 B_J=
 \begin{pmatrix}-1/2&0&1\\0&1&0\\0&0&-1/2\end{pmatrix},
\]

and \(H_2=t^2-zw+w^2\). A central sign corrects orientation where
necessary. These cases exhaust all three-dimensional self-adjoint Jordan
types. Phase scaling changes the other coefficients, which were arbitrary
throughout; it never discards them.

## Exact cyclic determinant on every family

Use \(q=2zw+t^2\), coordinates \((z,t,w)\), and
\(v=(\ell_w,\ell_t,\ell_z)^T\). Direct determinant identities give

| \(H_2\) | \(D=\det(v,Bv,B^2v)\) |
| --- | --- |
| \(0\) | \(0\) |
| \(w^2\) | \(0\) |
| \(zt\) | \(\ell_w^3/8\) |
| \(t^2-zw+w^2\) | \(9\ell_z^2\ell_t/4\) |
| \(z^2+w^2+\eta(t^2-zw)\) | \(\ell_t(\ell_z-\ell_w)(\ell_z+\ell_w)(9\eta^2-4)/4\) |

Hence the complete \(D=0\) equations in these normalized inputs are
explicit products of linear factors and the repeated-eigenvalue factor.
No moving-frame inverse or nonvanishing assumption on \(v\) is needed.

For the semisimple row,

\[
 p=\tfrac12\operatorname{tr}(B^2)=1+3\eta^2/4,\qquad
 r=\det B=\eta(\eta^2/4-1),\qquad
 4p^3-27r^2=(9\eta^2-4)^2/4.
\]

These formulas distinguish the complete scalar/eigenvalue boundaries
without a coefficient-support threshold. For the regular nilpotent row,
the cyclic open is simply \(\ell_w\ne0\); its complement retains every
cubic coupling.

## What a proof on these families would establish

Pure moments and all-original-multiplier safety are preserved by the
linear isometries and nonzero phase amplitudes. Under rotation, spherical
markers stay polynomial with the same degree. Hopf coefficient balancing
therefore preserves the full marked Image conclusion. Vanishing over an
algebraic extension descends to the original field by injectivity.

It is consequently enough to establish purity-implies-safety on all five
**full** coefficient families. The first is already closed by the
[odd theorem](ODD_HOPF_CUBIC_IMAGE_TRANSFER.md). The others remain open
with arbitrary nonzero linear harmonic. At zero linear harmonic the
nilpotent rows are closed, and the
[axial theorem](AXIAL_QUADRATIC_ZERO_LINEAR_EXCLUSION.md) excludes the
semisimple repeated-eigenvalue values. The non-nilpotent Jordan row and
other semisimple values remain open even on that zero-linear boundary.

These coefficient families are not asserted to be pure components.
Generic safety descent still requires an actual pure component or a
certified finite-prefix radical. Local Rodrigues multiplicities and
ordinary coefficient-ring extensions keep their separate obligations.
In particular **do not divide51 or17 by a residual normalization symmetry**:
these field-cover normalizations do not identify the fixed-D cyclic
ordinary algebra or its scheme lengths.

## Replay

Run `python verify_quadratic_orbit_normal_forms.py`. The checker verifies
each self-adjoint trace-free matrix, harmonic quadratic, nilpotent index
or nonsemisimple minimal polynomial, the semisimple eigenvectors,
characteristic invariants and every determinant factor above. The
[machine-readable full phase inputs](quadratic_orbit_normal_form_input.json)
retain all ten nonquadratic coefficients on every row and every
\(\eta\) value. The
[receipt](quadratic_orbit_normal_forms_verification.json) records hashes.

The normalization and determinant identities have no Singular or msolve
dependency. Assurance of imported safety subcases remains as stated in
their own proofs. There is no independent specialist, external or formal
verification. The full(3,3) Image theorem remains **OPEN**.
