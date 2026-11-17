# Pure critical-boundary multipliers force the matched chart

Powers of the critical source coordinate do not provide a surviving
unmatched quadratic input. The obstruction is independent of the
power: a polynomial weighted boundary matrix cannot be nilpotent unless
its offending linear form vanishes.

The result closes the pure-u scalar class across all four arbitrary
matrix families, not the whole unmatched architecture.
Mixed multipliers retain an additional logarithmic rank-one term,
so they require a separate analysis.

## 1. What a boundary factor forces, including mixed scalars

Work in rank4 over algebraically closed characteristic-zero constants.
Let V be homogeneous quadratic, detDV=delta*u^4 nonzero, and suppose
a homogeneous scalar F gives a normalized input.
The determinant identity gives f=-lambda*u*F.
Let u^e be the exact power of u in F and write
\[
                    F=u^eG,\quad u\nmid G,\quad e>0.
\]

**KMFI146.** A critical-boundary factor forces
\[
                         V_u=uL,\qquad L\text{ linear}.   \tag{1}
\]
For the four forms V=A*K of KMFI143, this means the first output row is
\[
 A_{02}=A_{03}=0\quad(K_A,K_B),\qquad
 A_{03}=0\quad(K_C,K_D).                                  \tag{2}
\]
These linear restrictions do not yet make the two charts match.

Indeed eta=dlogF has residue e*du/u, while the other scalar factors
are regular at a generic point of u0. The trace identity
trDV+eta*V=4lambda*u has no pole, so e*V_u/u is regular and u divides V_u.
One may equivalently use the separate boundary-factor residue from KMFI144.
The quadratic terms in each first output row give precisely(2).

There is also an exact boundary matrix which retains mixed factors.
Put L=a*u+ell*x, Q=V_high|_(u=0), and G0=G(0,x).
In the boundary function field set xi=dlogG0.
For
\[
 N=DV+V\,d\log F-\lambda uI,\qquad T=\operatorname{diag}(u,1,1,1),
\]
the rational matrix T^-1*N*T is regular at a generic point of u0.
Its limit is
\[
 C_e+
 \begin{pmatrix}\ell x\\Q\end{pmatrix}(0,\xi),\qquad
 C_e=\begin{pmatrix}(e+1)\ell x&\ell\\eQ&DQ\end{pmatrix}.   \tag{3}
\]
Normalization makes this limit nilpotent.
The extra term in(3) is generally nonzero; dropping it would incorrectly
extend the pure-power obstruction to mixed multipliers.

If the target boundary alpha of KMFI142 differs from u, then Q has
linear image span of dimension at most2. Both u and alpha annihilate
the full boundary vector V|_(u=0), and their independence leaves a
nonzero constant row on the three high components.

## 2. The weighted plane-image matrix cannot be nilpotent

**KMFI147.** Let Q:k^3->k^3 be homogeneous quadratic with image in
a linear plane, and ell a constant row. For a constant weight
s not in{-1,-2}, the polynomial matrix
\[
                   C_s=\begin{pmatrix}
                    (s+1)\ell x&\ell\\sQ&DQ
                    \end{pmatrix}                        \tag{4}
\]
cannot be nilpotent identically unless ell0.
In particular this holds for every nonnegative integer s.

If ell annihilates the image plane, choose that plane as z0.
The last output row is zero, and the upper left three-dimensional
block has (s+1)ell*x as an eigenvalue. Nilpotency gives ell0.
If ell is nonzero on the plane, choose simultaneous source/target
coordinates x=ell*x, y and z with that plane z0.
The last row is again zero. Restrict to z0.
The first three rows and columns give the nilpotent matrix
\[
 C=\begin{pmatrix}
       (s+1)x&1&0\\
       sR_1&(R_1)_x&(R_1)_y\\
       sR_2&(R_2)_x&(R_2)_y
    \end{pmatrix},                                       \tag{5}
\]
where R1,R2 are binary homogeneous quadratics.
Euler's identity gives exactly
\[
               \det C=\frac{s+2}{2}\,x\,\det D(R_1,R_2).
\]
Since s!=-2, nilpotency makes the binary Jacobian zero.
This makes the two quadratics proportional: writing
Ri=x^2*ri(y/x), their Jacobian is
2x^2*(r1*r2'-r1'*r2), and a rational function with derivative zero
is constant in characteristic zero.

If R1=0, matrix(5) has (s+1)x as an eigenvalue, impossible for s!=-1.
Otherwise replace y by y-k*x to make the second quadratic zero,
preserving the first coordinate x.
The remaining two-dimensional block has trace and determinant
\[
                     R_x+(s+1)x=0,\qquad
                     (s+1)xR_x-sR=0.
\]
Homogeneity and the trace give
\[
                        R=-\frac{s+1}{2}x^2+b y^2.
\]
The determinant's x^2 coefficient is
\[
                            -\frac{(s+1)(s+2)}2,
\]
nonzero under the hypothesis. This is the contradiction.

The two negative weights are real exceptions to the matrix lemma:
ell=(1,0,0), Q=0 works for s=-1;
ell=(1,0,0), Q=(x^2/2,0,0) works for s=-2.
These rational-weight controls are not polynomial scalar multipliers.

## 3. All pure powers of u are excluded to primitive full rank

**KMFI148.** Let V be primitive homogeneous quadratic in rank4 with
detDV=delta*u^4, delta nonzero.
No scalar F=constant*u^s, for any integer s>=0, produces a normalized
input of generic rankDX4.

More precisely, for s>=1 normalization alone forces the source and
target boundary hyperplanes to match. Full rank is needed to finish
the exclusion of the resulting primitive matched direction.

For s0 the reconstructed X has degree<=2 and KMFI117 excludes full rank.
For s>=1, KMFI146 gives V_u=uL.
Write L=a*u+ell*x and put Q=V_high|_(u=0).
If the target boundary alpha is already proportional to u, the matched
class is closed by KMFI141.
Otherwise Q has image in a plane by KMFI146.

For a pure power there is no extra term in(3).
In more explicit coordinates,
\[
 N=DV+(s/u)V\,du-\lambda uI
\]
is nilpotent, and conjugation by T=diag(u,1,1,1) gives a polynomial
matrix. Its u0 limit is exactly C_s.
Nilpotency extends to that limit. KMFI147 forces ell0.
Thus V_u=a*u^2, and a is nonzero because detDV is nonzero.
After scalar rescaling the component is matched, contradicting either
the assumed unmatched target boundary or, under full rank, KMFI141.
This proves the all-s obstruction.

The full-rank qualification is sharp. For every positive integer s,
set lambda=s+2 and
\[
 V=(u^2,\lambda ux,\lambda uy+x^2,\lambda uz+xy),\quad F=u^s.
\]
Then f=-lambda*u^(s+1), g=F(V-lambda*u*p) give a normalized polynomial
input with whole sweep Jacobian gamma^4.
Its matrix P is strictly lower triangular and generically regular
nilpotent, but X is independent of z and has rank3.
This is a boundary control, not a new Keller action.

## 4. Evidence and the remaining scalar arrangements

[Helper](boundary_multiplier.py) ·
[checker](verify_boundary_multiplier.py) ·
[model](boundary_multiplier_model.json) ·
[receipt](boundary_multiplier_receipt.json).

The exact controls include generic boundary conjugation, mixed-factor
correction, both binary identities and exceptional weights, all four
output-row restrictions, and the normalized rank3 family at several s.
The proof is unbounded in the scalar exponent and arbitrary output matrix.

Every remaining primitive quadratic full-rank candidate therefore needs
at least one non-u scalar factor, of degree1,2 or3 by KMFI145.
It cannot be the unmatched target hyperplane alpha by KMFI144.
A power of u can still accompany other factors, subject to(2)--(3);
such mixed factors have not been excluded.
The fixed-direction feedback, integer residues, coframe, exact generic
action and ordinary-affine completion gates remain necessary.

The absolute exceptional11 action and full faithful-transitive Keller
monodromy spectrum remain open. No new occurrence group is claimed.
Formal verification, independent whole-proof replay and external review
are not claimed; the goal remains active.
