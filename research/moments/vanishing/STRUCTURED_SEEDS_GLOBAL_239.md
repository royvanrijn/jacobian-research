# Structured seed modules: global239, homogeneous216 and mixed-parent210

15 November 2026. Internal written refinement of the
[shared-divisor258 proof](SHARED_DIVISOR_GLOBAL_258.md), retaining its complete
class/rank cover, original Wick functional and analytic exclusions.

| Scope | Sufficient pure gate/exclusion | Least-prefix interval |
|---|---:|---:|
| Every complex ternary polynomial of total degree at most three | **239** | **[20,239]** |
| Odd ternary cubics | **238** | **[20,238]** |
| All homogeneous ternary cubics | **216** | **[14,216]** |
| Full sixteen-coefficient W-linear family | **216** | Not computed |
| Nondegenerate quadratic W coefficient, cubic G | **216**, an exclusion | Not computed |
| Fourteen-coefficient mixed parent | **210** | Not computed |

The global test requires only the initial pure moments. No extra Gaussian
frame hypothesis is imposed. Pure nullity then gives the retained W-linear
frame and sharp mixed cutoff m>3*deg Q. Smooth-leading[20,94] and
homogeneous-smooth[14,56] are unchanged. Both isolated algebraic lower
witnesses with first nonzero moments20 and14 remain complete. No numerical
leastness or20/14 sufficiency is proved.

The [checker](verify_structured_seed_modules.py),
[certificate](structured_seed_modules_certificate.json),
[receipt](structured_seed_modules_verification.json) and
[source audit](structured_seed_modules_source_audit.json) bind auxiliary exact
algebra. Complete coefficient/rank coverage, all-degree reduction inductions
and analytic nonzero zeroth coefficients remain written mathematics.
Independent replay, external review, formal verification and expanded global
radical cofactors remain outstanding.

## 1. The conic phase has matching powers at every negative endpoint

Use the retained balanced logarithmic normalization

    I=constant*t^-1*integral exp(W/t^2) dy/y.

Both rescaled conic coordinates have the form

    constant*y+(c0+c1*t)+constant*K(t)/y,
    deg_t K<=2.

Consequently, for cubic G including arbitrary lower terms, the entire phase
has the stronger structure

    W=A0(t)+sum_(k=1)^3 (A_k(t)*y^k+K(t)^k*B_k(t)*y^(-k)),
    deg A0<=3, deg A_k,deg B_k<=3-k.                 (1)

The Gaussian quadratic satisfies stronger bounds. To prove(1), expand a
coordinate monomial: an exponent -k requires at least k negative factors,
so K^k divides its coefficient. After removing those factors, the remaining
t-degree is at most3-k. This includes both f!=0 and f=0,e!=0 charts from
the preceding conic proof. The trapped formal component retains negative
order0, and identically zero K belongs to a lower actual endpoint case.

If either actual endpoint order is below3, the retained all-rank bounds
apply: p<3,q=3 has maximum201 and q<3 has maximum145. It remains to treat
p=q=3, where A3 and B3 are nonzero constants. No generic seed rank or
coefficient specialization is assumed.

Retain the previous exact logarithmic Euler representative, with

    B=t*K_prime-2*K/3, h=t^3*K,
    S_K=K*(t*W_t-2W-t^2)+B*y*W_y.

The actual seed recurrence is

    R_0=1,
    R_(j+1)=h*R_j_prime-j*h_prime*R_j
                     +t^2*B*y*(R_j)_y+S_K*R_j,     (2)

with support[-2j,3j]. Both the Euler amplitude derivative and the unchanged
logarithmic volume normalization are required.

## 2. Low reductions divide factors already in the amplitudes

Consider the C[t]-module of Laurent amplitudes whose coefficient of y^k,
for k<0, is divisible by K^(-k). Equation(1) and the factor K in h show
that the numerator operator in(2) preserves this module. Differentiating
a coefficient can lower its K divisibility by one, which multiplication
by K restores. The other derivative and source terms preserve the required
negative powers. Thus every R_j lies in the module.

Use the unchanged full twisted basis y^-3,...,y^2 and the exact operator

    t^2*Dlog(f)=t^2*y*f_y+y*W_y*f.

A high reduction divides by3*A3, a nonzero constant. It preserves the
module: if a resulting exponent is negative, its phase coefficient
contains at least the required power of K. A low reduction of c(t)*y^n,
n<-3, divides c by -3*B3*K^3. Since c is already divisible by K^(-n),
the division is polynomial. Its new coefficient retains K^(-n-3), and
all phase or Euler terms again have the required powers at their resulting
negative exponents. This proves by induction that the complete normal
forms of R_j have polynomial t coefficients.

This is a statement about actual seed amplitudes. Arbitrary unscaled
negative monomial basis columns need not have that polynomial property.
No assertion about a denominator-free full connection matrix is made.
The quotient, its full rank and the original formal functional are unchanged;
no singular t=0 integral or fibrewise basis is substituted.

In particular the seed denominator is only H_j=h^j, without padding it
by powers of A_-3. Its valuation is at least3j even when K(0)=0.
The sufficient high and low depths remain3j-2 and max(0,2j-3).

## 3. Stratify by the actual degree of K

Let ell=deg_t K, which is0,1 or2. Use the integer weights and whole-step
costs below. A high step has cost phase_weight-3*weight_y. A low step
has cost phase_weight+3*weight_y-3*ell*weight_t, because its polynomial
division removes K^3. The Euler term has no larger weight.

| ell | weight(t),weight(y) | Phase weight at most | Seed R_j weight at most | High cost | Low cost |
|---|---|---:|---:|---:|---:|
| 0 | 1,1 | 3 | 3j | 0 | 6 |
| 1 | 2,1 | 6 | 8j | 3 | 3 |
| 2 | 1,1 | 3 | 5j | 0 | 0 |

For a chosen row exponent k, the reduced coefficient degree satisfies
weight_t*deg_t N_(k,j)+weight_y*k<=b_j, where

    ell=0: b_j=3j+6*max(0,2j-3);
    ell=1: b_j=8j+3*(3j-2+max(0,2j-3));
    ell=2: b_j=5j.

At actual seed rank s select rows containing0 and use the retained
zeroth-cofactor estimate

    deg a0 <= floor((sum_(j=1)^s b_j-weight_y*sum_selected k)/weight_t).

All smaller ranks and row choices are included. The largest values occur
at full rank6 and are216,199,108 for ell=0,1,2. For ell=1 the unfloored
bound is199.5; polynomial degree gives199. The selected row sum at full
rank is-3. These three bounds apply to p=q=3 only; the other actual
endpoint cases retain201/145 as stated above.

The analytic nonzero zeroth coefficient is unchanged because the connection
and seed are unchanged. Higher scalar coefficients contain h^j and have
valuation at least3j, so no extra derivative-order loss is needed.
The resulting uniform cubic conic exclusion is max(216,201,145)=216.
This number is not asserted for arbitrary-degree G.

## 4. Mixed-parent degrees5/6 have monomial leading coefficients

Keep the normalized mixed-parent interior phase

    F=z/t-a-z^2, Phi=t*G(z,F)-F^2/2,
    I=constant*t^-1*integral exp(Phi) dz,
    V=t^2*Phi.

The actual degree d is4,5 or6. On d=5 or6, the leading coefficient of
V is a nonzero constant times t^3, by the retained triangular degree split.
The old source is S=t*V_t-2V-t^2. Add the exact ordinary-volume Euler
operator -Dop(z*f)/(d*t), where Dop(f)=f_z+V_z*f/t^2. The same quotient
connection has representative

    nabla_E(f)=f_t-z*f_z/(d*t)+(S_E/t^3)*f,
    S_E=S-t^2/d-z*V_z/d.                           (3)

The term -t^2/d is its volume divergence and cannot be omitted. It cancels
the highest source amplitude together with -z*V_z/d, giving source
amplitude degree<=d-1 and coefficient t-degree<=3.
Actual seeds have denominator t^(3j) and recurrence

    R_(j+1)=t^3*R_j_prime-3j*t^2*R_j
                            -t^2*z*(R_j)_z/d+S_E*R_j.

They have amplitude degree<=(d-1)j and coefficient t-degree<=3j.
Leading reduction uses depth(d-1)*(j-1)+1 and a degree-three principal
coefficient. Thus the numerator bounds are3j+3*((d-1)*(j-1)+1), giving
114 for degree5 and210 for degree6 at their full actual seed ranks.
There is no basis-grade subtraction from this parameter-only estimate.
The degree4 chart retains72, while c=0 and d=0 faces retain35/92.
Hence the full fourteen-coefficient mixed family has sufficient gate210.

With rank-one140 retained, the complete W-linear gate is now216. Combining
homogeneous smooth56,cusp10,reducible10,nonisotropic-node101 and conic216
at an isotropic node gives all-homogeneous216. Parity does not lower even216.

## 5. A further exact shear lowers the tangent source weight

For the tangent D!=0 complement use the retained weighted Euler connection,
with delta=t^2D, source S of joint degree6 and weighted degree7, and
beta_x=delta_prime/4, beta_y=tD-delta_prime/8. Write C(t)=[x^2*y]S,
so deg_t C<=3. The only monomial of weight7 under weights(2,3) is x^2*y.

Add -C*G_y(x*f)/(2*D*delta^2). The vertical field x*partial_y has zero
ordinary volume divergence. Clearing h=delta^2*D gives

    S_sh=D*S-C*x*W_y/2,
    nabla_sh(f)=f_t+(beta_x/delta)*x*f_x+(beta_y/delta)*y*f_y
                       -C*x*f_y/(2*D*delta)+(S_sh/h)*f. (4)

Since W_y has leading term2D*x*y, the weight7 source cancels. Thus
S_sh has joint degree<=7 and weighted degree<=6. No new parameter open
is introduced: D is already nonzero on this complement. C may vanish
identically, and D(0)=0 is allowed.

Actual numerators over h^j satisfy

    R_(j+1)=h*R_j_prime-j*h_prime*R_j+S_sh*R_j
         +delta*D*(beta_x*x*(R_j)_x+beta_y*y*(R_j)_y)
         -delta*C*x*(R_j)_y/2.

The difference from the retained weighted connection on the same amplitude
is the exact image G_y(-C*x*R_j/2). The shear derivative lowers weight by1;
all other amplitude derivatives preserve it. Hence R_j has joint degree
at most7j and weighted degree at most6j. The retained principal divisor
D^2 and depth6j-4 give row degree19j-8-d_i. Higher scalar valuations
are at least4j, so the tangent complement gate is

    19*15-8*5-6=239.

The D=0 branches now use mixed-parent210 or rank-one140, so the full tangent
gate is239. The cusp complement remains204, with full cusp max(204,210,140)=210.

## 6. The global239 cover and assurance

| Leading type, arbitrary lower terms | Full-class sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible nodal cubic | max(227,216)=227 |
| Irreducible cuspidal cubic | 210 |
| Transverse line-smooth-conic | max(220,216)=220 |
| Tangent line-smooth-conic | 239 |
| Three independent lines | max(204,216)=216 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The maximum239 proves the pure-only global test. For odd polynomials
moment239 vanishes automatically, giving the odd gate238. In the rational
twenty-coefficient ring, radical(I239)=radical(I_infinity), without ordinary
ideal equality or expanded global cofactors. At most ceil(log2(239-19))=8
exact monotone search rounds would find the least K in[20,239]; none has run.

The checker verifies12 coordinate-support cases, all three p=q=3 K-degree
cofactor strata and retained other endpoint bounds. It checks polynomial
coefficients and K divisibility after every Laurent cancellation in nine
models, with27 direct seeds. Ten symbolic mixed-parent identities and three
polynomial models give9 seed certificates. Fifteen symbolic tangent shear
identities and four full-Gaussian models give20 connection and12 seed
certificates, including D(0)=0 and C identically zero. The initial successful
checker and its outputs are preserved; the strengthened receipt clarifies
the p=q=3 scope and checks parent denominator valuations explicitly.

No new external theorem or producer code is imported. The retained period
comparison, logarithmic and shrinking/escape exclusions, full rank and case
cover remain independently unreviewed written inputs. The historical proof
bodies stay hash-identical. The next largest gates are tangent239,
ordinary-node227 and transverse220. Exact leastness,20/14 sufficiency,
expanded radical certificates, independent review, exceptional schemes,
term-count minima and lower-power GVC remain open.
