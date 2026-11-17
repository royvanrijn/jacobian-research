# Centered odd concurrent cubics: class36 and complete odd prefix100

15 November 2026. Internal written refinement of the [conic96/odd104 cover](CONIC_ROW_FACTORS_96_ODD104.md), retaining the original Gaussian Wick functional and full coefficient/frame cover.

Every complex odd ternary polynomial of total degree at most three has all positive pure Gaussian moments zero if and only if its first **100** pure moments vanish. The least odd prefix is therefore in **[20,100]**. The unrestricted interval remains **[20,130]** and the homogeneous interval remains **[14,96]**. Conic/full-W-linear96, smooth-leading[20,94], homogeneous-smooth[14,56] and both exact first-nonzero20/14 witnesses remain unchanged. All final upper bounds are even; parity gives no99 reduction from100.

The new restricted class bound is **36** for an odd cubic whose leading form is three distinct concurrent lines. The generic concurrent104 bound with arbitrary lower terms is retained. This refinement uses oddness to remove all quadratic polynomial terms; it does not replace arbitrary lower terms by centered data without justification.

The [checker](verify_odd_concurrent_prefix.py), [certificate](odd_concurrent_prefix_certificate.json), [receipt](odd_concurrent_prefix_verification.json) and [source audit](odd_concurrent_prefix_source_audit.json) support exact auxiliary algebra. Universal all-degree reduction, complete vertex/frame coverage, the retained nonisotropic no-constant-period argument and formal Airy--Riccati obstruction remain written inputs. There is no independent replay, external review, formal verification or expanded global radical certificate.

## 1. Oddness forces constant H on the nonisotropic vertex branch

Choose the unique concurrent vertex as the T-axis. All three leading linear factors vanish on that axis, so the cubic leading part depends only on u,v. A binary linear change, with no translations, normalizes the squarefree binary cubic to u^3+v^3. Since P is odd, its quadratic and constant parts vanish. Thus

    P=h*T+u^3+v^3+e*u+f*v,
    Q=m*T^2/2+T*J+q2, J=ju*u+jv*v,
    q2=(a*u^2+2c*u*v+b*v^2)/2.

The Gaussian stays centered; every cross coefficient is retained. If m!=0, exact completion gives

    tP-Q=-m*(T-(t*h-J)/m)^2/2
                       +t*(u^3+v^3+e*u+f*v)-q2+(t*h-J)^2/(2m).

The Schur quadratic q=q2-J^2/(2m) is nondegenerate, since det Hess Q=m*det Hess q!=0. The normalized eliminated Gaussian factor is independent of t. Rescale u,v by1/t and retain the binary volume factor t^-2. The normalized period phase numerator is

    W=u^3+v^3-q+t^2*L+c0*t^4,
    L=e*u+f*v-h*J/m, c0=h^2/(2m), delta=t^2.

No affine mean is introduced by the leading-form normalization. The t-dependent constant from completing T is the displayed c0*t^4 term and is retained.

## 2. Keep the full ordinary-volume radial Euler source

The original rescaled connection has source S0=t*W_t-2W-2t^2 over h0=t^3. Add the exact ordinary twisted derivatives with witnesses(2/3)u*f,(2/3)v*f. The same quotient connection is

    nabla(f)=f_t+(2/(3t))*(u*f_u+v*f_v)+S*f/t^3,
    S=t*W_t-2W+(2/3)*(u*W_u+v*W_v)-(2/3)*t^2
      =(2/3)*q+(2/3)*t^2*L+2*c0*t^4-(2/3)*t^2.

The scalar divergence from two ordinary variables is4t^2/3 before combining it with the old -2t^2. It is not a logarithmic-measure divergence. This cancels the top cubic source while preserving the same original Wick seed.

Give t,u,v weight1. The source has joint degree at most4, or at most3 if h=0. Its amplitude degree is at most2. The full actual numerator recurrence is

    R0=1,
    R_(j+1)=t^3*(R_j)_t-3j*t^2*R_j
             +(2/3)*t^2*(u*(R_j)_u+v*(R_j)_v)+S*R_j.

Every derivative term adds at most2 to total degree. Hence joint degree R_j<=c*j with c=4, or c=3 on h=0, and amplitude degree<=2j. No derivative of the scalar denominator is omitted.

Use the full rank-four basis1,u,v,u*v and principal gradients3u^2,3v^2. Dividing a monomial outside that basis produces a witness of amplitude degree input-2. The lower derivatives -q_u+t^2*L_u and -q_v+t^2*L_v have joint degree at most2 and amplitude degree at most1; the delta-divergence correction has smaller joint degree. Thus each branch preserves joint degree and strictly lowers ordinary amplitude degree. Depth2j-1 suffices for every seed. All leading divisors are nonzero constants in time. Completed normal rows obey

    deg_t N_(e,j)<=c*j-|e|, H_j=t^(3j).

Every actual seed rank and selected row set containing1 gives zeroth-cofactor degree at most

    c*s*(s+1)/2-sum_selected |e|.

At full rank4 the row degree sum is4, giving4*10-4=36, or3*10-4=26 on h=0. The same retained [concurrent nonisotropic no-constant-period exclusion](CONCURRENT_LINE_EXCLUSION.md) applies to the original connection and seed, so the zeroth scalar coefficient is nonzero. Higher scalar coefficients retain H_j of valuation3j, and the initial-order contradiction adds no derivative-order loss. Thus36 excludes the entire m!=0 odd concurrent complement.

## 3. Treat the isotropic vertex rather than discarding it

If m=0, the vertex is isotropic. A centered Gaussian Witt frame puts the exponent into Z*Wc+T^2/2, with E(Wc^i*Z^j)=i!*delta_ij. Because the leading cubic is independent of the vertex direction, oddness gives

    P=a*Wc+G3(Z,T)+r*Z+s*T,
    G3=g30*Z^3+g21*Z^2*T+g12*Z*T^2+g03*T^3.

There is no quadratic W coefficient on this branch. The retained formal Gaussian--Lagrange identity has the exact constant fixed point Z=a*t and determinant1, giving

    I=E_T exp(t*G3(a*t,T)+r*a*t^2+s*t*T).

This identity is formal and coefficientwise; no singular Gaussian integration is used. Rescale T=y/t, retaining volume t^-1. Its phase numerator is

    V=g30*a^3*t^6+g21*a^2*t^4*y+g12*a*t^2*y^2
                         +g03*y^3+r*a*t^4+s*t^2*y-y^2/2.

If g03!=0, retain h0=t^3 and add the exact ordinary Euler witness(2/3)y*f. The new source is

    S=t*V_t-2V+(2/3)*y*V_y-t^2/3.

It has amplitude degree<=2 and weighted degree<=6 under weights(t,y)=(1,2). In the complete basis1,y, principal division by3g03*y^2 preserves this weight: both lower derivative terms have weight<=4, while the t^2-divergence term is smaller. The actual recurrence retains the t^3 derivative, -3j*t^2 and(2/3)*t^2*y*partial_y terms. Completed normal rows satisfy6j-2k with H_j=t^(3j). Every rank/row choice gives gate6*3-2=16.

The zeroth coefficient is nonzero by the retained [formal cubic Gaussian Airy--Riccati obstruction](CUBIC_HEIGHT_AIRY_OBSTRUCTION.md). Its statement applies to any nonzero solution of the cubic Gaussian system. Here its cubic source A=t*g03 has positive order, B=-1/2+a*g12*t^2 has unit constant term and C=a^2*g21*t^3+s*t vanishes at zero. A constant normalized scalar period would make the cubic Gaussian solution have rational logarithmic derivative after removing exp(g30*a^3*t^4+r*a*t^2), contradicting that obstruction. This also covers a=0.

If g03=0, the remaining Gaussian is exact. Put

    q=1-2*g12*a*t^2, E=g30*a^3*t^4+r*a*t^2,
    F=(g21*a^2*t^3+s*t)^2,
    I=q^(-1/2)*exp(E+F/(2q)).

Since q(0)=1, its logarithmic derivative is rational:

    I_prime/I=N/q^2,
    N=-q_prime*q/2+E_prime*q^2+(F_prime*q-F*q_prime)/2,
    deg N<=7.

If the first8 pure moments vanish, I-1 has order at least9 and I_prime/I has order at least8. Its denominator is a unit and its numerator degree at most7, so N=0 identically. Then I=I(0)=1. This proves gate8 on every flat case, including q=1 and all pure-null members. Combining16 and8 gives an isotropic-vertex gate16, below the nonisotropic36 bound. The squarefree leading restriction is retained in the full concurrent claim; the constant-W frame calculations themselves allow any G3.

## 4. Complete odd100 cover and remaining bounds

The complete restricted cover is now

| Leading type | Odd sufficient prefix |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible node | 100 |
| Irreducible cusp | 94 |
| Transverse line-smooth-conic | 98 |
| Tangent line-smooth-conic | 94 |
| Three independent lines | 96 |
| Three distinct concurrent lines | **36** |
| Repeated leading line(s) | 98 |

Its maximum100 proves the odd-only sufficient prefix. The odd nodal complement now controls this estimate. No Gaussian-frame condition is added to the original polynomial. This is a full restricted coefficient/class cover, not a parity deduction from even global130. Homogeneous96 and conic/full-W-linear96 stay unchanged. The unrestricted cover keeps its generic concurrent104 branch and maximum130.

## Evidence and assurance

The accepted checker verifies the centered completion, Schur determinant and rescaling;15 binary source-gauge and62 principal identities;9 cubic-frame source-gauge and11 principal identities;16 binary and12 cubic-frame actual seed identities in eight models; the generic flat Gaussian logarithmic derivative and32 direct positive moments in four models; and every cofactor rank/row calculation. Binary models retain nonzero full metric determinant and Gaussian cross coefficients, with h=0/nonzero and m=1 and2. Frame models include a=0 and g03=1 or2; flat models include a=0 and q=1. The first run passed, in about3 seconds; no failed input is hidden or promoted.

These exact identities support auxiliary algebra. Universal all-degree reduction, full vertex/frame choices, minimal scalar seed/period comparison, the nonisotropic analytic exclusion and algebraic Airy--Riccati theorem remain written, independently unreviewed inputs. No fresh literature-coverage claim, external theorem or producer package is imported; the selected primary-source refresh is retained. Earlier canonical proof bodies and receipts remain byte-identical.

Least global[20,130],odd[20,100],homogeneous[14,96],20/14 sufficiency, expanded I130 radical cofactors, independent analytic review, exceptional schemes, global term minima and lower-power GVC remain open. Generic cusp130/node129, odd node100 and homogeneous conic96 are the next controlling bounds. Qualitative cubic coverage and both lower witnesses stay closed internally.
