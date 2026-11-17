# The complete nilpotent cyclic slice and a length-one global gate

Current refinement: the [contraction-zero theorem](NILPOTENT_CONTRACTION_ZERO_AND_THREE_GLOBAL_TARGETS.md)
now proves J12+(gminus2,gminus1,g0)=P0. For N>=12, three linear exact
or radical memberships suffice and are equivalent to branch length1.
All three remain open. The five-target result below remains valid and
retains its earlier prefix9 augmented identity.

17 November 2026. **SIC2B33NILPCYCLIC20261117.** On the entire regular
nilpotent quadratic stratum with \(D\ne0\), field-valued phase scaling
and oriented isometry reduce the full phase to

\[
 F=w+az+bt+zt+\sum_{j=-3}^{3}\gamma_j C_j.
\]

All nine coefficients \(a,b,\gamma_{-3},\ldots,\gamma_3\) are free.
There is no extra cubic restriction. Let \(J_N=(\mu_2,\ldots,\mu_N)\)
in this nine-variable polynomial ring and write \(B_N=\mathbb C[\cdot]/J_N\).
The Rodrigues point has \(a=b=0\), \(\gamma_3=-1/2\), all other
\(\gamma_j=0\); let its maximal ideal be \(P_0\).

The proved conclusions are:

- The normalized slice is polynomially isomorphic to the complete cyclic
  structural branch \(p=r=0,s_2=1/4\), including every cubic coupling.
- The Rodrigues local ring on this slice is **reduced of length one**.
  The first ten moments are the minimal consecutive ordinary local
  generating prefix. Independently, restricting the supplied actual
  seventeen-dimensional algebra by \((p,r)\) gives length one.
- There is a global augmented ideal identity
  \[
  J_9+(a,b,\gamma_{-3},\gamma_{-1},\gamma_1)=P_0.
  \]
- Consequently, for every \(N\ge10\), global length one is equivalent
  to **five exact linear memberships**, or five separate radical tests,
  on these five displayed functions. Those unrestricted memberships
  remain **OPEN**. Global length one is not established.

The full(3,3) goal remains open. This is an exact global-input and
certificate reduction on a complete intrinsic open stratum, not an
all-pure-point exhaustion theorem.

## 1. A complete intrinsic slice

Normalize a nonzero rank-two nilpotent self-adjoint quadratic operator
to \(H_2=zt\), retaining every linear and cubic coefficient. In coordinates
\((z,t,w)\),

\[
 B=\begin{pmatrix}0&0&0\\1/2&0&0\\0&1/2&0\end{pmatrix},\qquad
 L=az+bt+cw,\qquad v=(c,b,a)^T,
\]

and \(D=c^3/8\). Thus the cyclic open is exactly \(c\ne0\).
The torus \((z,t,w)\mapsto(\rho z,t,\rho^{-1}w)\), combined with phase
amplitude \(\rho^{-1}\), fixes \(zt\) and sends \(c\mapsto c\rho^{-2}\).
Choose \(\rho^2=c\). This gives the displayed slice over an algebraic
field extension. Every remaining coefficient stays arbitrary. Marker
degrees are preserved under the isometry, and a nonzero phase amplitude
preserves vanishing. Resulting polynomial vanishing descends to the
original characteristic-zero field.

For \(c=1\), the cyclic frame and its dual coordinates are

\[
 E=(v,Bv,B^2v)=
 \begin{pmatrix}1&0&0\\b&1/2&0\\a&b/2&1/4\end{pmatrix},\qquad
 y_0=w+bt+az,\quad y_1=(t+bz)/2,\quad y_2=z/4.
\]

Its determinant is the constant \(1/8\), and its polynomial inverse is

\[
 z=4y_2,\qquad t=2y_1-4by_2,\qquad
 w=y_0-2by_1+4(b^2-a)y_2.
\]

The Gram coordinates are exactly

\[
 p=r=0,\qquad s_0=2a+b^2,\quad s_1=b,\quad s_2=1/4.
\]

The inverse on linear coefficients is \(b=s_1\),
\(a=(s_0-s_1^2)/2\). Substituting the polynomial inverse into the full
cubic gives its ten symmetric \(U_{ijk}\) coefficients. The three trace
relations solve \(U_{011},U_{111},U_{112}\) with constant pivots; every
other coefficient is retained. Conversely their trace-zero cubic pulls
back to an arbitrary harmonic cubic. The checker verifies both directions,
so this is a polynomial isomorphism, not a moving-frame localization with
unexamined denominator cases.

On the full normalized cyclic restriction \(p=r=0,D=1/8\), the Gram
determinant gives \(s_2^3=1/64\). Its three distinct branches are
isomorphic over \(\mathbb C\) to the slice above. In the fixed-\(B\)
coordinates they have \(c^3=1,s_2=c^2/4\); the torus/amplitude action
relates the three branch algebras. Therefore

\[
 \dim\bigl(A_N/(p,r)\bigr)=3\dim B_N
\]

when finite. The optional nilpotent cyclic exhaustion benchmark is thus
**length three**, rather than51. This changes the ideal by adjoining
\(p,r\); it does not divide out a relative sign or identify the ordinary
full cyclic algebra with a smaller quotient by symmetry.

## 2. The Rodrigues point is reduced on this whole slice

In the variable order \((a,b,\gamma_{-3},\ldots,\gamma_3)\), direct
spatial integration gives

\[
 \det\left(\frac{\partial(\mu_2,\ldots,\mu_{10})}
 {\partial(a,b,\gamma_{-3},\ldots,\gamma_3)}\right)_{F_0}
 =-\frac{134217728}{94732977184601460705}\ne0.
\]

Hence these nine equations cut out a reduced isolated local point.
Earlier consecutive prefixes contain at most eight equations in nine
variables and cannot isolate that point. The first eight Jacobian rows
have rank eight. The retained all-order Rodrigues purity makes every
later moment zero in this local ring. Thus the ordinary local prefix ten
is exact, on this **restricted slice**.

There is a separate finite-algebra check using the supplied actual
seventeen-dimensional multiplication matrices. The image matrix
\((M_p\;M_r)\) has rank16, with zero constant row. A retained sixteen
by sixteen minor has determinant

\[
 16553/4043047157760\ne0.
\]

Thus \(L_{17}/(p,r)\cong\mathbb C\). This restriction accounts for the
loss of local thickness; it is not a replacement of the original algebra
by its radical. The full local length remains17, and the full cyclic
length51 remains an independent open benchmark.

## 3. Five linear targets: a global augmented ideal identity

Adjoin
\(a=b=\gamma_{-3}=\gamma_{-1}=\gamma_1=0\). The complete remaining
family is

\[
 F=w+zt+A w^2t+B(t^3-3ztw)+Cz^2t+D z^3.
\]

Its primitive second and third moments are

\[
 M_2=4AC+15B^2,\qquad M_3=-B.
\]

The fourth moment reduces modulo \((B,AC)\) to
\(3432(D+1/2)\). Thus \(D+1/2\) already belongs to the selected
first-four ideal, with constant pivots throughout.

After setting \(D=-1/2\), the primitive fifth residual is

\[
 R_5=B(24A^2-252AC-153B^2)-26(3A-11C).
\]

Together with \(M_2,M_3\), this gives the exact ideal
\((B,3A-11C,C^2)\). Every residual through8 belongs to it. The ninth
residual has remainder \(19228C\), so the first-nine ideal is exactly
\((A,B,C,D+1/2)\). The
[certificate](nilpotent_cyclic_slice_five_target_certificate.json) exports
four rational polynomial multiplier rows proving membership of each
generator. Every row is checked by direct multiplication over
\(\mathbb Q\). No Groebner unit oracle is used for these identities.

Taking the preimage under the five-coordinate quotient proves the full
augmented identity \(J_9+(a,b,\gamma_{-3},\gamma_{-1},\gamma_1)=P_0\).
It does **not** prove any of those five functions already belongs to
\(J_9\) or a later unrestricted moment ideal.

## 4. The precise optional global gate

For \(N\ge10\) the following are equivalent:

1. \(\dim_{\mathbb C}B_N=1\).
2. Each of the five linear targets belongs to \(J_N\).
3. Each belongs to \(\sqrt{J_N}\).
4. The reduced zero set of the first \(N\) moments is the Rodrigues point.

The augmented identity proves the support implication. Finite support
and the known reduced local point then force length one. Conversely a
length-one quotient containing the rational Rodrigues point is its residue
field, so every target vanishes exactly. Separate radical tests may use
\(1\in J_N+(1-u f)\); they must not be combined into a single inverse
test. Later moment prefixes are allowed. Witness degrees are unknown.

The [exported program](nilpotent_cyclic_slice_five_target_program.json)
uses the exact first14 input and five degree-one targets. It also records
the mod2 symmetry: variable characters are \((0,1,0,1,0,1,0,1,0)\), and
\(\mu_n\) has character \(n\). Each inverse variable gets its target's
character. It records **zero verified unrestricted memberships**.

If this sufficient gate is proved at any finite prefix, every pure phase
on the entire regular-nilpotent cyclic stratum is a Rodrigues copy and
therefore safe. Additional pure components could instead be acceptable
if their marked safety is proved; the full Image goal does not require
Rodrigues-only exhaustion.

## Evidence and computational boundary

Run `python verify_nilpotent_cyclic_slice.py`. The checker regenerates the
[exact sparse moment input through14](nilpotent_cyclic_slice_prefix14.json.gz),
independently checks spatial moments through6 and every Rodrigues
derivative through10, verifies the polynomial frame adapter, restricts the
actual17 algebra, and multiplies all four selected-family witness rows.
The [receipt](nilpotent_cyclic_slice_verification.json) records hashes.
All proved results have no Singular or msolve dependency.

Full modular std and F4 scouts stopped at45 seconds without a certified
basis. An ell_z elimination-block F4 scout also timed out. A modular
F4SAT test on ell_z was first rejected at the small prime, then terminated
with signal11 at a supported prime under the3GiB resource cap. All
inputs, logs and receipts are retained. None proves membership,
nonmembership, global finiteness or a characteristic-zero outsider.
Do not blindly restart them.

The \(D=0\) complement, other quadratic types and the full(3,3) Image
theorem remain open. Independent specialist, external and formal
verification are not claimed.

The [ker(B) linear-vector boundary checkpoint](NILPOTENT_KERNEL_LINEAR_REDUCTION.md)
retains a separate complete noncyclic chart reduction, with two exact
rational exclusions and one unresolved five-variable chart. The
[highest-cubic global chart](HIGHEST_CUBIC_GLOBALIZATION_AND_GRADED_GATE.md)
includes both cyclic and boundary cases in a positively graded finite-module
framework; its dimension theorem does not settle either purity classification.
