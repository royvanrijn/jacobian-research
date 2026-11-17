# Nilpotent B with zero cubic contraction: a complete safe container and three global targets

17 November2026. Claim `SIC2B33THREEGATE20261117`.
**The conditional classification below is proved. All three unrestricted
global moment radical memberships remain OPEN.**

Remove the constant by the first moment and write \(F=H_1+H_2+H_3\).
Let \(B\) be the trace-free self-adjoint operator for \(H_2\), put
\(p=\operatorname{tr}(B^2)/2\), \(r=\det B\), and define the linear
harmonic covariant

\[
 \mathcal C(F)=\Delta^2(H_2H_3),\qquad
 \Delta=2\partial_z\partial_w+\partial_t^2.
\]

This is a nonzero scalar multiple of the spin-one projection of
\(H_2H_3\), equivalently the contraction of the cubic tensor with \(B\).
It is bilinear in the quadratic and cubic coefficient vectors, and
commutes with every orthogonal change of coordinates.

**Theorem.** Every all-order pure field point of the closed container

\[
 p=r=0,\qquad \mathcal C(F)=0                         \tag{1}
\]

is one-sided or an orthogonal/scalar Rodrigues copy. It is safe for
every original multiplier, with uniform cutoff \(m>3\deg_{Z,Y}Q\),
sharp on Rodrigues. All three linear coefficients are unrestricted
before purity. Every cyclic and noncyclic nilpotent boundary is included.

For **nonzero** \(B\), twelve initial moments already force this
classification under (1). For \(B=0\), we use the existing all-order
odd theorem; no assertion that twelve odd moments suffice is made.
The container is not entirely pure, and we have not proved that every
pure phase satisfies (1). No general coefficient-ring safety is inferred.

The [checker](verify_nilpotent_contraction_zero.py) and
[receipt](nilpotent_contraction_zero_verification.json) support the
exact algebra below. The largest four-variable exclusion retains an
explicit exact QQ Singular dependency. Smaller identities are directly
multiplied. There is no msolve dependency or independent specialist,
external or formal verification.

## 1. The exhaustive nilpotent quadratic split

The characteristic polynomial is \(T^3-pT-r\), so \(p=r=0\) is exactly
nilpotence. The [complete quadratic normal forms](QUADRATIC_ORBIT_NORMAL_FORMS.md)
give all possibilities: \(B=0\), rank-one \(H_2=w^2\), or regular
nilpotent \(H_2=zt\), after field extension and allowed normalizations.
All linear and cubic coefficients stay in the input until the stated
covariant equations are imposed.

At \(B=0\), the retained odd classification applies to
\(H_3+qH_1\) by Gaussian polar decomposition. Its sphere pure points
are exactly one-sided, with cutoff coefficient one.

## 2. Rank one, with arbitrary linear harmonic

In the circular cubic basis, direct differentiation gives

\[
 \Delta^2(w^2H_3)=-8\gamma_1w+16\gamma_2t+48\gamma_3z.
\]

Thus the contraction condition kills exactly the three positive cubic
coefficients. Write \(H_1=l_zz+l_tt+l_ww\); the other four cubic
coefficients remain free. The third moment is then
\(\mu_3=2l_z^2/5\). At a field point this forces \(l_z=0\).
All remaining noncentral terms have negative circular charge. Every
pure moment therefore equals the moment of

\[
 l_tt+\gamma_0(t^3-3ztw)
 =l_tt+\gamma_0(5t^3-3t)/2
\]

on the unit sphere. If \(\gamma_0=0\), its second moment forces
\(l_t=0\). Otherwise put \(u=l_t/\gamma_0\). The primitive second
and fourth moments are

\[
 7u^2+3,\qquad
 3003u^4+3432u^3+6578u^2+1872u+723.
\]

They are coprime. An [explicit rational Bezout identity](nilpotent_contraction_zero_rank1_central_unit.json)
is checked by direct multiplication. Consequently \(\gamma_0=l_t=0\).
The entire remaining phase has strictly negative circular charge:

\[
 F=l_ww+w^2+\gamma_{-3}w^3+\gamma_{-2}w^2t
       +\gamma_{-1}w(t^2-zw/2).
\]

It is pure and full-Image safe by the existing weight argument. This
proof retains arbitrary \(H_1\), rather than assuming it zero.

## 3. Regular nilpotent, with arbitrary linear harmonic

Now

\[
 \Delta^2(ztH_3)=-24\gamma_0z+16\gamma_{-1}t+16\gamma_{-2}w.
\]

All phases satisfying the contraction equation have exactly the form

\[
 F=az+bt+cw+zt+dw^3+ez(t^2-zw/2)+fz^2t+gz^3.          \tag{2}
\]

The [exact moment adapter](nilpotent_contraction_zero.py) restores the
uniquely recoverable power of the linear-\(w\) coefficient in the
previous unrestricted cyclic generator. It does not incorrectly read
that coefficient from an already specialized polynomial. Its first
equations are

\[
 \mu_2=(2ac+b^2)/3+4dg/35,\qquad
 \mu_3=2bc/5+4df/105.                                \tag{3}
\]

Direct spatial integrations through six independently check the adapter.

### The case d=0

If \(c=0\), (3) forces \(b=0\). The phase is exactly

\[
 az+zt+ez(t^2-zw/2)+fz^2t+gz^3,
\]

which is strictly positive in circular charge and hence safe.

If \(c\ne0\), (3) forces \(b=a=0\). The torus followed by phase
amplitude sends

\[
 (a,b,c,d,e,f,g)\mapsto
 (a,b/\rho,c/\rho^2,d/\rho^4,e,\rho f,\rho^2g),
\]

while fixing the coefficient of \(zt\). Choose \(\rho^2=c\) to
normalize \(c=1\). The primitive residuals through eight are

\[
\begin{aligned}
 R_4&=e^2+6g+3,& R_5&=f,\\
 R_6&=-e^3-39eg+13e+13f^2,& R_7&=ef,\\
 R_8&=2e^4+34e^2g+17e^2+68g^2+68g+17.
\end{aligned}
\]

Their ordinary ideal is exactly \((e,f,2g+1)\). The
[three explicit multiplier columns](nilpotent_contraction_zero_rodrigues_witness.log)
are checked by exact polynomial multiplication. For example, after
eliminating \(f,g\), the last constraints are proportional to
\(e(11e^2+65)\) and \(e^4\), whose gcd is \(e\).
Thus this case is Rodrigues. Undoing the normalization gives

\[
 F=cw+zt-z^3/(2c),\qquad c\ne0.
\]

### The case d is nonzero

Choose \(\rho^4=d\). Equations (3) have the constant-pivot solutions

\[
 d=1,\qquad f=-21bc/2,\qquad g=-35(2ac+b^2)/12.
\]

The remaining four variables \((a,b,c,e)\) are completely unrestricted,
including \(c=0\). Exact primitive moment residuals4..12 are retained
in [the chart input](nilpotent_contraction_zero_charts.json).
Their exact QQ standard basis is \((1)\): see the
[complete input](nilpotent_contraction_zero_lowest_minus3_prefix12.sing)
and [output](nilpotent_contraction_zero_lowest_minus3_prefix12.log).
The checker regenerates every coefficient and replays this calculation
in about15 seconds. It excludes this entire case in characteristic zero.

An optional independent expanded unit export stopped at40 seconds and
is preserved. The exact QQ standard-basis dependency remains explicit;
the timeout is not presented as an exported unit certificate. No finite-field
unit is used for this exclusion.

This completes all nilpotent quadratic ranks and all linear/cubic cases
under (1). The existing full original-multiplier one-sided and Rodrigues
theorems give the asserted cutoff. It is a field-valued classification;
the field step \(l_z^2=0\Rightarrow l_z=0\), for example, does not
claim an ordinary coefficient-ring identity.

## 4. Three targets on the complete global chart

Use the twelve-variable normalized chart of the
[global finite-module theorem](HIGHEST_CUBIC_GLOBALIZATION_AND_GRADED_GATE.md).
Its phase has \(\gamma_3=1,\gamma_2=0\), and \(\mu_2\) eliminates
\(\gamma_{-3}\) polynomially. Define

\[
\begin{aligned}
 P={}&4p=3b_0^2+2b_1b_{-1}+4b_2b_{-2},\\
 R={}&4r=b_0^3+b_0b_1b_{-1}-4b_0b_2b_{-2}
                    +b_1^2b_{-2}+b_2b_{-1}^2,\\
 T={}&-\tfrac18[z]\mathcal C
       =b_2c_{-1}+3b_1c_0-3b_0c_1-6b_{-2}.
\end{aligned}                                                    \tag{4}
\]

Their ordinary coefficient degrees are \((2,3,2)\) and their positive
weights are \((6,9,5)\). The first and third already occur in the earlier
31-target ideal. These are the new preferred global targets:

\[
 \boxed{P,R,T\in\sqrt{(R_3,\ldots,R_N)}}              \tag{5}
\]

for any workable finite \(N\ge12\). All three memberships remain open.
Each is a separate radical test; combining inverse conditions would
test the wrong open set.

The extra point requiring proof is that the single coefficient \(T\)
controls the **whole** covariant. Here is a finite argument. Before
normalization its general coefficient is

\[
 -\tfrac18[z]\mathcal C
 =b_2\gamma_{-1}+3b_1\gamma_0-3b_0\gamma_1
       -2b_{-1}\gamma_2-6b_{-2}\gamma_3.
\]

It is fixed by the lowering rotation \(\exp(uM)\), since that
rotation fixes the \(z\) coefficient of any linear harmonic.
Phase amplitude \(\lambda\) multiplies it by \(\lambda^2\).

Take the **nine** fixed rotations: first \(\exp(sM)\), then the
oriented swap, for \(s=0,1,\ldots,8\). For a nonzero cubic harmonic
its resulting \(z^3\) coefficient is a nonzero polynomial of degree at
most six in \(s\), so at least three of these rotations enter the
normalization chart. Membership of \(T\) in (5) implies that the
rotated covariant's \(z\) coefficient vanishes at each such value.
That coefficient is a polynomial of degree at most two in \(s\): for
\(\mathcal C=v_zz+v_tt+v_ww\) it is
\(v_w+sv_t-s^2v_z/2\). Three distinct zeros force all three vector
coefficients to vanish. The checker verifies the covariant and
normalization identities exactly. Using only the earlier seven rotations
would not guarantee three good values and is not the argument used here.

The scalar conditions \(P=R=0\) transfer directly through a good
normalization, because they are rotational invariants of amplitude
degrees two and three. Thus (5) forces every pure cubic phase into the
closed safe container (1). The lower \(H_3=0\) complement is already
covered by the full(2,2) theorem. This proves the sufficiency of (5) for
the entire field-valued(3,3) Image theorem, including all \(D=0\).

This gate allows first-\(N\) odd false positives on \(B=0\): the
all-order odd theorem handles the pure points there. We do not need to
make a small prefix classify the whole odd locus first. At the level of
the **full pure locus**, the three radical conditions are equivalent to
exhaustion by the known nullcone/Rodrigues families. That classification
is stronger than Image if other safe pure components exist.

For a fixed power \(f^e\) of weight \(ed\), only moments
\(n\le\lfloor ed/3\rfloor\) and cofactor monomials of weight
\(ed-3n\) are needed. The known necessary exponent bounds are
\((7,5,7)\): the first two follow from actual multiplication in the
supplied17-dimensional Rodrigues algebra; amplitude changes multiply
\(p,r\) by units and preserve their nilpotence indices. The third
comes from the retained exact18-dimensional power obstruction.
No assertion that those powers suffice globally is made.

The [three-target program](highest_cubic_three_target_gate.json) records
all polynomials, weights, lower bounds, nine rotations and **zero proved
global memberships**.

## 5. The cyclic nilpotent gate also shrinks to three linear functions

On the previously certified complete cyclic nilpotent slice
\(F=w+az+bt+zt+H_3\), the covariant is zero precisely when
\(\gamma_{-2}=\gamma_{-1}=\gamma_0=0\). The theorem above forces
every first-twelve point of this intersection to be Rodrigues. The
already certified first-ten Jacobian makes its actual local ring reduced
of length one. Finite support and that local ring prove the **ordinary**
global augmented identity

\[
 J_{12}+(\gamma_{-2},\gamma_{-1},\gamma_0)=P_0.
\]

This is a proof by global support plus actual local length, not an
export of full nine-variable multiplier rows. Consequently for every
\(N\ge12\), the following are equivalent:

* the branch algebra \(B_N\) has length one;
* the three displayed linear functions belong to \(J_N\);
* all three belong to \(\sqrt{J_N}\).

The [cyclic three-target program](nilpotent_cyclic_three_target_program.json)
keeps the full prefix14 source and permits later prefixes. Its mod2
characters are \((1,0,1)\). All three unrestricted memberships and the
normalized three-branch length3 remain open. No ordinary witness-degree
bound is known in that inhomogeneous cyclic presentation. There is no
division of local or global length by a sign symmetry.

## Evidence and remaining global obligation

Run `python verify_nilpotent_contraction_zero.py`. It checks the two
rank covariant formulas; the rank-one central Bezout identity; the full
seven-parameter moment adapter and independent direct moments; every
coefficient of both regular-nilpotent branches; the Rodrigues multiplier
columns; the four-variable exact QQ unit; the determinant/invariant and
spin-one targets; the finite rotation adapters; and the actual local
nilpotence bounds. The dependent closed-family classifications and the
cyclic local Jacobian are imported at their retained assurance scopes.

An earlier weighted-homogeneous modular reconstruction attempt for the
larger **unrestricted** ker\(B\) boundary stopped at90 seconds. That
source, log and receipt remain preserved. It proves no exclusion and
does not close that unrestricted boundary. The new theorem closes it
only after the declared contraction condition is imposed.

The remaining task is to prove that all pure phases have nilpotent
quadratic part and zero contraction, or find and control components
outside this container. Equivalently, pursue the three global radical
tests (5), or a broader safe-component argument if they fail. The new
conditional theorem and gate do not complete full Image. Goal ACTIVE.
