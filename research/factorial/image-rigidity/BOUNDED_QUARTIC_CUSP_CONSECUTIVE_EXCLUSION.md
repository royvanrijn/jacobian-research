# Actual consecutive-moment exclusion in a bounded quartic-cusp window

The theorem in this note concerns each fixed integer cusp and every fixed bounded scaled field window. It supplies no global root localization, no effective order cutoff, and no all-order Strong Factorial theorem for a cusp support.

## 1. Statement and actual moment normalization

Fix an integer s>=1, put a=s^2 and d=2s+1, and define

    A_m(x)=sum_(k=0)^m binom(m,k)*(a*m+1)_(d*k)*x^(m-k),

where the subscript denotes a rising factorial. These are the actual ordinary-binomial moment polynomials: for B!=0 and Lambda(X^q)=q!,

    Lambda([X^a(A+B X^d)]^m)=B^m*(a*m)!*A_m(A/B).

Let R=(d+1)/(d-1), delta=2d/(d-1),

    z_c=(R^(d+1)+1)/delta,
    w_c=a^d*R^(d+1)>0,
    beta=-R^(-d-1/2),  K_w=a^d*delta/beta.

**Theorem BoundedQuarticCuspConsecutiveExclusion.** For each finite L>=0 there is N_0(s,L) such that for every integer N>=N_0 the actual pair A_N,A_(N-1) has no common complex root x satisfying

    |x/N^d-w_c| <= L*N^(-3/4).                 (B1)

All constants in this result may depend on s and L. The domain is a complex analytic window, not an algebraic coefficient restriction over an unspecified field.

For real m>0 set

    I_m(w)=integral_0^infinity [S^a(w+S^d)]^m*exp(-m*S)dS.

On a sufficiently small complex w-disc about w_c, use the branch that is positive for real w. This integral is holomorphic in w, and locally holomorphic in complex m with Re(m)>0. The claim follows by domination near S=0 by S^(a*Re(m)) and at infinity by a polynomial times exp(-Re(m)*S); a small neighbourhood of each fixed m has uniform dominating bounds. For integer m,

    I_m(w)=(a*m)!/m^((a+d)*m+1)*A_m(m^d*w).    (B2)

Thus a common actual x gives I_N(w_N)=I_n(w_n)=0 with n=N-1 and

    w_N=x/N^d,  w_n=(N/n)^d*w_N.              (B3)

## 2. Exact cusp jets and a linear scaling

Use the actual normalized phase

    p=1+delta*v/d,
    F(v,z)=d^2/delta^2*(log(p)-delta*v/d)
                 +log(z+(p^d-1)/delta),
    S=a+sqrt(a)*v, w=a^d*(-1+delta*z).

The original I_m(w) is the v-integral multiplied by the nonzero scalar

    sqrt(a)*[a^(a+d)*delta*exp(-a)]^m.

Indeed factor S^a(w+S^d)=a^(a+d)*delta*p^a*(z+h(v)) and exp(-m*S)=exp(-a*m)*exp(-m*sqrt(a)*v).

At v=1,z=z_c its first three v derivatives vanish, and

    F_4=-2/R^2,
    F_5=24/[(d-1)*R^3],
    F_vz=-R^(-d-1), F_vvz=0.                  (B4)

The exact proof and Fraction replay are in QUARTIC_CUSP_UNFOLDING_AND_ACTUAL_ORDER_DRIFT.md. Put A=sqrt(R), epsilon=m^(-1/4), and scale linearly,

    v=1+A*epsilon*u,
    z=z_c+epsilon^3*h/beta,
    beta=A*F_vz.

After removing the nonzero scalar exp(m*F(1,z))*A*epsilon from the v-integral, denote the resulting period by J_epsilon(h). The linear coordinate avoids a spurious real-temperature assumption: F_vvz=0 holds in these raw coordinates, even though the full nonlinear quartic normal form has nonzero mixed temperature derivative.

**Lemma ActualQuarticTwoTermExpansion.** Uniformly for h in each fixed compact complex disc, as real epsilon decreases to zero,

    J_epsilon(h)=Z(h)+epsilon*gamma_d*Z^(5)(h)+O_(d,H)(epsilon^2),
    gamma_d=1/[5*sqrt(d^2-1)],
    Z(h)=integral_R exp(-u^4/12+h*u)du.         (B5)

The function extends C^2 in real epsilon to epsilon=0, locally uniformly holomorphically in h, with the same assertion for each fixed h derivative. In particular the O(epsilon^2) remainder is uniform on a slightly larger h-disc. No analyticity at epsilon=0 is asserted.

### Local domination and differentiated expansion

Fix a sufficiently small real v-neighbourhood |v-1|<r. At the cusp its real phase is F_c-c(v-1)^4+O((v-1)^5), with c>0. Shrink r so the quartic coefficient retains a negative bound. By subtracting F(1,z), write the scaled exponent exactly as

    E_epsilon(u,h)=epsilon^(-4)*[F(1+A*epsilon*u,z)-F(1,z)].

The base phase is x^4*g(x), with analytic g and g(0)=F_4/24. The parameter-dependent difference after subtraction is x*(z-z_c)*C(x,z-z_c), with C analytic and C(0,0)=F_vz. Therefore, on the local integration interval,

    Re(E_epsilon) <= -c_1*u^4+C_H*|u|,
    |partial_epsilon E_epsilon|
          <= C_H*(|u|^5+|u|^2+|u|),
    |partial_epsilon^2 E_epsilon|
          <= C_H*(|u|^6+|u|^3+|u|^2+|u|).    (B6)

The constants apply for all small positive epsilon. Derivatives of exp(E) are consequently bounded by an integrable polynomial times exp(-c_1*u^4+C_H*|u|). Endpoints u=+-r/(A*epsilon) contribute exponentially small boundary terms when differentiating the local integral. Taylor expansion and dominated differentiation give

    E_0=-u^4/12+h*u,
    partial_epsilon E_0=gamma_d*u^5.           (B7)

Indeed the fifth base term has coefficient F_5*A^5/120=gamma_d, whereas the possible order-epsilon parameter term is F_vvz*A^2*h/(2*beta)=0. This proves the local part of (B5), including bounded second epsilon derivatives and h derivatives.

### Complete original-cycle tails

At z=z_c the positive real cycle has a unique global quartic maximum v=1. For completeness, with D=a+d and S_c=s(s+1), the stationary map

    W(S)=S^d*(D-S)/(S-a)

on a<S<D has derivative

    W'(S)=-d*S^(d-1)*(S-S_c)^2/(S-a)^2.

It decreases from infinity to zero; the phase derivative is positive before S_c and negative after it, also outside (a,D). Away from the fixed local interval there is therefore a strict phase gap g_0>0. At the left endpoint the phase tends to minus infinity, and at infinity it is a logarithm minus S, so

    integral exp(F(v,z_c)-F_c)*(1+|F(v,z_c)-F_c|^2)dv < infinity.

On the whole original real cycle z_c+(p^d-1)/delta >= z_c-1/delta>0. Hence changing z by epsilon^3*h/beta changes F and F(1,z) uniformly by O_H(epsilon^3). The nonlocal normalized integral is bounded by a power of epsilon^(-1) times exp(-g_0/(2*epsilon^4)). Differentiating it once or twice in epsilon only introduces further powers and integrable powers of the base phase. Such terms are flat at epsilon=0. This proves the complete-period assertion (B5), not just a truncated local approximation.

The apparent moving endpoints in u are merely endpoints of the local decomposition. The original v endpoints remain fixed; the original integral can be differentiated there before the local change of variables.


## 3. The model zeros and their actual finite-window branches

Newman's primary Theorem2, Proceedings AMS61(1976),245-251, printed246-248, applies to the empty-product quartic density in his formula(1.8). Taking alpha=1/12,beta=-tau/2-1,b=1 in that theorem gives imaginary h zeros for

    Z(h,tau)=integral_R exp(-u^4/12+tau*u^2/2+h*u)du

for every real tau. The finite measure premise holds for arbitrary real tau because its quartic coefficient is negative. The actual primary PDF was read at quartic_lee_yang_sources/newman_fourier_transforms_1976.pdf; the primary PDF is hosted at https://sites.math.northwestern.edu/~auffing/papers/Newman.pdf . No Curie-Weiss limiting theorem is needed for this import.

These model zeros are simple. Integration by parts gives Z_hhh=3*tau*Z_h+3*h*Z and Z_tau=Z_hh/2. Set f(x,tau)=Z(ix,tau). A zero of multiplicity at least three would give zero initial data for f_xxx+3*tau*f_x-3*x*f=0, impossible since f(0,tau)>0. At a hypothetical double zero (x_0,tau_0), let c=f_xx(x_0,tau_0)/2!=0. Then

    eta^(-2)*f(x_0+eta*q,tau_0-eta^2)
        tends to c*(q^2+1)

uniformly on compact q sets. Rouche gives nonreal zeros for real tau_0-eta^2, contradicting Newman. In particular the zeros of Z(h)=Z(h,0) are simple. Since Z is a nonzero entire function, only finitely many lie in a fixed compact h set.

Choose H_1>H+2 with no model zero on |h|=H_1. Take disjoint small discs around the finitely many model zeros h_j inside it. Uniform convergence in (B5) and Rouche give exactly one actual simple zero h_j(epsilon) in each disc, and no other actual zeros in |h|<=H+1, for sufficiently small epsilon. The one-sided C^2 function at epsilon=0 may be extended to negative epsilon by its quadratic Taylor jet, whose coefficients are holomorphic in h. This is a C^2 extension, so the ordinary real-parameter implicit-function theorem gives C^2 branches through epsilon=0, with bounded derivatives on a common small interval. Explicitly,

    h_j(epsilon)=h_j-epsilon*gamma_d*
                   Z^(5)(h_j)/Z'(h_j)+O_(d,H)(epsilon^2).  (B8)

The scalar ratio is finite because Z'(h_j)!=0. This construction includes every actual zero in the chosen bounded h window; no assertion is made about indices escaping that window.

For consistency with the nonlinear canonical dictionary, at tau=0 the exact ODE implies

    Z^(5)=6*Z'+3*h*Z''.

Thus the coefficient in (B8) is exactly -b*[1+h_j*Z''(h_j)/(2*Z'(h_j))], b=6*gamma_d. The complex-temperature and amplitude corrections of the canonical normal form are both present in this raw quintic term; neither has been discarded.

## 4. Adjacent-order separation

For m sufficiently large, the actual zero branch corresponding to h_j is

    w_j(m)=w_c+K_w*m^(-3/4)*h_j(m^(-1/4)).     (B9)

Since the finitely many h_j and their epsilon derivatives are bounded,

    dw_j/dm=-(K_w/4)*[3*m^(-7/4)*h_j
                              +m^(-2)*h_j']
             =O_(d,H)(m^(-7/4)).              (B10)

Here h_j' denotes the epsilon derivative. This differentiates the proved C^2 branch, not a bare asymptotic remainder. The integral itself is locally holomorphic for positive real m, but no holomorphic extension through epsilon=0 is required.

Suppose now that x is a common actual root satisfying (B1), and put h_N=K_w^(-1)*N^(3/4)*(w_N-w_c). It is bounded. For n=N-1, kappa=N/n, the exact lower scaled field is

    h_n=kappa^(d-3/4)*h_N
                    +(w_c/K_w)*n^(3/4)*(kappa^d-1).        (B11)

Consequently h_n-h_N=O_(d,H)(N^(-1/4)), and h_n lies in |h|<=H+1 after choosing H to bound the upper fields. As the finitely many distinct model zeros are separated, the upper and lower zeros must belong to the same branch j for all sufficiently large N. Formally, choose the branch discs smaller than one quarter the minimum model-zero separation, and then make both their perturbations and the difference (B11) smaller than that separation margin.

The common-root equation would then be

    kappa^d*w_j(N)=w_j(n).

But (B9)-(B10) give

    kappa^d*w_j(N)-w_j(n)
       =(kappa^d-1)*w_c
          +(kappa^d-1)*(w_j(N)-w_c)
          +w_j(N)-w_j(n)
       =d*w_c/N+O_(d,H)(N^(-7/4)).             (B12)

The leading coefficient d*w_c is positive and nonzero, so this cannot vanish for all sufficiently large N. This proves the theorem.

Equivalently, the explicit two-term expansion is

    w_j(m)=w_c+K_w*h_j*m^(-3/4)
        -K_w*gamma_d*Z^(5)(h_j)/Z'(h_j)*m^(-1)
        +O_(d,H)(m^(-5/4)).                   (B13)

Even without (B10), comparing (B13) at N and N-1 makes the root difference o(1/N), while the exact radial transport in (B3) has size d*w_c/N. Thus no derivative estimate for an unspecified O-term is being assumed.

## 5. Scope

This proves actual consecutive-root exclusion on each fixed bounded quartic scale for each fixed integer cusp. It closes the actual-cycle and finite-index perturbation issues in that range using the original real integral. It does not provide an all-index relative quartic remainder, a uniform estimate when d or the scaled window grows, an effective N_0, or localization of every possible common root. In particular it does not prove the full SFC for the support (s^2,2s+1).

The finite exact checker supports the universal written jet, scaling and differential algebra; it is not a formal verification of dominated convergence, Newman, or the analytic root argument.

