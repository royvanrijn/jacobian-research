# Full anti-halfturn scheme length and the next A14 gate

16 November 2026. The full-length lifting argument and the finite-cover
degree calculation were supplied by the user. The symbolic verification
below checks the oriented coordinate maps, including the invertible
cubic change, and retains their nilpotent multiplicities.

**SIC2B33CYCLICLENGTH20261116.** Taking the previously retained exact
msolve computation of modular chart length9168 as input,

\[
 \boxed{\dim_{\mathbb C}A_{12}\ge13818.}
\]

The normalized non-null anti-halfturn first12 scheme has special-fibre
length13752 and characteristic-zero vector-space dimension at least13752.
Its characteristic-zero scheme is removed entirely by adding \(\mu_{14}\),
using the inherited anti-Weyl unit theorem. The next global comparison is
\(\dim_{\mathbb C}A_{14}=51\), with both \(\mu_{13}\) and
\(\mu_{14}\) retained: the [17 November local certificate](RODRIGUES_LENGTH17_AND_GLOBAL51.md)
now proves \(d_{14}=d_\infty=17\). Global finiteness and exact lengths
remain uncomputed.

The claim inherits the exact msolve dependency of
[the first12 outsider certificate](CYCLIC_NILPOTENCE_GATE_AND_PREFIX12_OUTSIDERS.md).
No new Groebner run, independent Buchberger verification or all-order
pure component is asserted here.

## 1. Global torsion-freeness and lifting all9168 classes

Let \(\ell=1073741827\), \(R=\mathbb Z_{(\ell)}\), and let \(B\) be the
seven-equation algebra of the preceding certificate over \(R\). Its
variables are \((\eta,a_z,a_w,e_z,e_w,b_z,b_w)\); its equations are the
six primitive even moments through12 and

\[
 \eta(b_z-b_w)(b_z+b_w)-1=0.
\]

Every moment scaling factor and input denominator is a unit in \(R\).
The retained exact modular computation proves

\[
 \dim_{\mathbb F_\ell}B/\ell B=9168<\infty.
\]

At every special closed point, localizing the polynomial ambient ring
gives an eight-dimensional regular local ring. The seven equations
followed by \(\ell\) are a system of parameters, since the special fibre
is zero-dimensional. They form a regular sequence, so \(\ell\) is a
nonzerodivisor after imposing the equations. This is the dimension
criterion in [Stacks02JN](https://stacks.math.columbia.edu/tag/02JN).

This conclusion is global in \(B\). The finite \(B\)-module
\(\operatorname{Ann}_B(\ell)\) is supported on the special fibre, a finite
set of closed points, and all its localizations there vanish. It is
therefore zero. Hence \(B\) is torsion-free, and thus flat, over the DVR.

Choose lifts \(b_1,\ldots,b_{9168}\) of a special-fibre basis. A nonzero
\(\mathbb Q\)-linear relation among their images in \(B_\mathbb Q\) can
be cleared of denominators. Torsion-freeness puts that relation in \(B\).
Cancel the smallest common power of \(\ell\) among its coefficients.
At least one remaining coefficient is a unit, so reduction modulo
\(\ell\) contradicts the chosen basis. Therefore

\[
 \dim_{\mathbb Q}B_\mathbb Q\ge9168.                    \tag{1}
\]

Neither reducedness, finite generation as an \(R\)-module nor generic
finiteness is used. Other affine generic components may still avoid the
special fibre, so this argument supplies a lower bound, not an equality.

## 2. The entire normalized anti-halfturn locus

Fix \(D=D_0=1/8\) in the complete oriented cyclic-frame quotient, and
write \(a=s_0\), \(c=s_2\). On \(a\ne0\), impose

\[
 r=s_1=U_{001}=U_{012}=U_{111}=U_{122}=0.               \tag{2}
\]

The Gram matrix and companion are compatible with the half-turn
\(J=\operatorname{diag}(-1,1,-1)\):
\(J^TGJ=G\), \(JCJ=-C\) and \(\det J=1\).
Its complete anti-invariant cubic eigenspace is exactly the remaining
six symmetric tensor coordinates. Thus (2) describes the full normalized
anti-halfturn locus with non-null linear harmonic, scheme-theoretically.
In particular all odd moments vanish on it as polynomial identities.

The determinant relation becomes

\[
 c(c^2-apc)=D_0^2.
\]

It forces \(c\) to be a unit, with inverse
\((c^2-apc)/D_0^2\); no additional open condition has been imposed. Since
\(a\) is already a unit, solve

\[
 p=\frac{c-D_0^2/c^2}{a}.
\]

Two harmonic traces solve

\[
\begin{aligned}
 U_{011}&=-\frac acU_{000}-2U_{002}-pU_{022},\\
 U_{112}&=-\frac acU_{002}-2U_{022}-pU_{222};
\end{aligned}
\]

the third trace vanishes identically. Consequently the ambient normalized
locus is the Laurent polynomial algebra

\[
 C_0=R[a^{\pm1},c^{\pm1},U_{000},U_{002},U_{022},U_{222}].
\]

Let \(C_{12}\) denote its quotient by the even moments through12. It is
also the quotient of \(A_{12,R}[a^{-1}]\) by (2). The moment conventions
are the ordered symmetric tensor sum of the preceding cyclic-frame note.

## 3. The rank-six and rank-four finite etale covers

Adjoin to the original chart the amplitude satisfying

\[
 E=B[\lambda]/(\lambda^6D-D_0),\qquad
 D=\frac{(b_z-b_w)(b_z+b_w)^2}{8}.
\]

Both \(D\) and \(D_0\) are units. This is a monic degree-six root equation
after division by \(D\); \(\lambda\) is a unit and its derivative
\(6\lambda^5\) is a unit. Thus \(E\) is finite etale and free of rank6
over \(B\). This uses that2 and3 are units in \(R\). The elementary monic
criterion is recorded in [Stacks03PA](https://stacks.math.columbia.edu/tag/03PA).

For the same cover over \(C_{12}\), put

\[
 k=\frac{b_z+b_w}{2},\qquad d=b_z-b_w.
\]

Scaling the phase by \(\lambda\) gives

\[
 a=2\lambda^2,\qquad c=\lambda^4k^2,
 \qquad D_0=\lambda^6dk^2/2.
\]

The inverse formulas are

\[
 \lambda^2=a/2,\qquad k^2=4c/a^2,\qquad
 d=\frac{4D_0}{ac},
\]

\[
 b_z=k+d/2,\qquad b_w=k-d/2,\qquad \eta=(2kd)^{-1}.
\]

In particular the sign of the oriented determinant is retained in the
formula for \(d\). Every denominator is an already declared unit. The
two monic quadratic adjunctions give

\[
 E\cong C_{12}[\lambda,k]/(\lambda^2-a/2,\ k^2-4c/a^2),             \tag{3}
\]

free of rank4 with basis \(1,\lambda,k,\lambda k\). Their derivatives
are units, so this cover is finite etale and faithfully flat as well.

To verify that (3) recovers the full cubic rather than just the base,
the normalized dual cyclic coordinates are

\[
 y_0=\lambda(z+w),\quad y_1=\lambda^2kt,\quad
 y_2=\frac{\lambda^3k}{2}\big((k+d/2)z+(k-d/2)w\big).
\]

The original cubic is
\(\lambda^{-1}\sum_{i,j,h}U_{ijh}y_i y_j y_h\). After the two trace
substitutions it is a linear function of the four free \(U\)'s, with
coefficients \((a_z,a_w,e_z,e_w)\) in the original harmonic basis.
The exact determinant of that coefficient change is

\[
 \det T=\frac9{16}d^6k^6\lambda^{20},
\]

a unit. The inverse is \(\operatorname{adj}(T)/\det T\), so this is an
invertible coordinate change over arbitrary coefficient algebras here.
The exact phase identity is \(F_{\rm normalized}=\lambda F_{\rm chart}\).
Every moment therefore pulls back by the unit factor \(\lambda^n\).
This proves the covers are compatible with the moment ideals, including
their nilpotent structure, at every order.

The [model](anti_involution_cover.py),
[checker](verify_anti_involution_cover.py) and
[receipt](anti_involution_cover_verification.json) verify the determinant,
trace, inverse and full cubic identities. The receipt retains the full
four-by-four cubic matrix and the input length's provenance hashes.

## 4. Scheme-length transfer and the global lower bound

Reducing the two finite free covers modulo \(\ell\) gives

\[
 \dim_{\mathbb F_\ell}E/\ell E=6\cdot9168=55008
   =4\dim_{\mathbb F_\ell}C_{12}/\ell C_{12}.
\]

Thus

\[
 \dim_{\mathbb F_\ell}C_{12}/\ell C_{12}=13752.           \tag{4}
\]

Since \(E\) is free over the torsion-free \(B\), it is torsion-free over
\(R\). The second finite free cover is faithfully flat, so flatness
descends to \(C_{12}\); one can also use the injection
\(C_{12}\hookrightarrow E\) to check torsion-freeness directly.
Repeating the basis-lifting argument from Section1 gives

\[
 \dim_{\mathbb Q}(C_{12})_\mathbb Q\ge13752.
\]

These are scheme lengths, with no replacement by reduced point counts.

If \(A_{12}\) is finite-dimensional, its \(a\ne0\) Artinian factor
surjects onto \((C_{12})_\mathbb C\), so that factor has length at least13752.
The three normalized Rodrigues local algebras have length22 each and
\(a=s_0\) nilpotent; they belong to the disjoint \(a=0\) factor. Therefore

\[
 \dim_{\mathbb C}A_{12}\ge13752+3\cdot22=13818.
\]

If \(A_{12}\) is infinite-dimensional, the same numerical lower bound
holds automatically. Finiteness and exact global length are not decided.

## 5. Removing the whole chart and preparing A14

The inherited theorem SIC2B33AW14 excludes common zeros of the first12
moments and \(\mu_{14}\) on the entire characteristic-zero non-null
anti-halfturn chart. Finite type and the Nullstellensatz imply

\[
 (\mu_{14})=(1)\quad\hbox{in }(C_{12})_\mathbb Q.
\]

Thus \(\mu_{14}\) is a unit in that algebra and adding it kills the full
scheme thickening. This conclusion does not require reducedness or an
explicit inverse. The historical moment14 artifact is not reconstructed
or replayed here; its theorem is an inherited dependency. No analogous
unit statement is asserted over the entire DVR or its special fibre.

For the full quotient the next exhaustion certificate is

\[
 \dim_{\mathbb C}A_{14}=51<\infty.
\]

The [actual local calculation](RODRIGUES_LENGTH17_AND_GLOBAL51.md) now
establishes \(d_{13}=18\), \(d_{14}=d_\infty=17\), minimal prefix14 and
nilpotence index7, including actual filtered multiplication. The earlier
qualification about constructing that multiplication is resolved.

Keep \(\mu_{13}\) in this global calculation. Its vanishing on the
anti-halfturn chart does not make it redundant on the full quotient.
The [A14 builder/checker](verify_cyclic_prefix14.py) now retains
[all fourteen exact moment circuits](cyclic_frame_moments_14.circuit.json.gz)
and a [verification receipt](cyclic_prefix14_verification.json).
It checks the entire old first12 prefix gate-for-gate, the Rodrigues
point through14, and direct integration of the non-symmetric control
\(F=tz+t+w^3+w+z\). Both \(\mu_{13}\) and \(\mu_{14}\) are nonzero
on that control. The full A14 input has15 variables and17 equations:
four structural equations and \(\mu_2,\ldots,\mu_{14}\).

The circuit has1341228 nodes and5036929 compressed bytes. Building this
input performs no global Groebner calculation. In particular the square
DVR argument is not automatically applicable to this overdetermined
global A14 presentation.

## Reproduction and assurance

From this directory:

```sh
python verify_anti_involution_cover.py
python verify_cyclic_prefix14.py
```

The cover check reads and validates the prior msolve receipt and primary
hashes, then recounts the retained special-fibre standard monomials.
The existing optional exact CAS replay remains available through
`verify_anti_involution_prefix.py --replay`; it is unnecessary to rerun
that completed computation for the cover argument.

No exact global A12/A14 length, new all-order
pure component or full(3,3) Image theorem is claimed. Independent
Buchberger, external specialist and formal verification remain absent.
