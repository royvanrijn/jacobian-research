# Joint Gaussian-safe algebras, finite inverse gates and nonclosure

16 November 2026. Internal written deductions. Independent replay, external review and formal verification are pending. No novelty, new moment-preserving map or smaller GMC endpoint is claimed. The proof uses the retained weighted inverse and nilpotent-section results, not a new external analytic theorem.

## 1. A tuple is an algebra question, not a list of scalar checks

Fix a characteristic-zero field K and a normalized Keller map F in n variables, with degree at most D>=1 and JF(x)Wx=WF(x), where W is diagonal with nonzero K-weights. Write g for its formal inverse and Y=JF(x)^(-T)p. The split Wick functional L has pairings L(x_i p_j)=W_i^(-1)delta_ij. Let A be ANY commutative K-algebra, and let V be an n by r matrix over A, r>=1. Define

    q_j=(W V_j)^T Y, q(t)=sum_j t_j q_j, b(t)=g(Vt).

Here V_j is column j, and t_1,...,t_r are independent formal parameters. Rank loss is allowed; no inverse of a minor of V is used. Let A[q]^+ denote the image of polynomials H(T_1,...,T_r) with H(0)=0. Call the tuple jointly safe when every such H(q) is safe against every polynomial source marker. All elements of this augmented algebra are pure-null, regardless of joint safety: the known cotangent map preserves moments and the target linear momentum forms have all joint positive moments zero.

The retained completed formal Gaussian-Lagrange identity gives, coefficientwise,

    L[Q(x) exp(q(t))]=Q(b(t)).                              (1)

The determinant cancellation remains identically1 after replacing the single linear direction by W Vt. In particular this is an identity over A, with nilpotents retained.

The following are equivalent:

1. The tuple is jointly safe.
2. Every constant K-linear combination sum_j c_j q_j is safe over A.
3. b(t) is a polynomial tuple in A[t].
4. The finite inverse-section ideal J_B(V) of Section2 has nilpotent image in A, where B=D^(n-1).

We first prove3 implies1, including arbitrary momentum markers. Introduce source rho for Wp. The preceding scalar formal proof now gives

    L[Q(x) exp(q(t)+rho^T Wp)]=Q(b(t,rho))/j_t(b(t,rho)),
    b=rho+W^(-1)JF(b)^(-1)W Vt,
    j_t(x)=det(I-J_x(W^(-1)JF(x)^(-1)W Vt)).                (2)

At rho=0, b=b(t) and j_t(b)=1. If b(t) is polynomial, every finite rho derivative is polynomial in t: implicit differentiation uses the polynomial adjugate of a determinant-one matrix. This proves that the full source-marker generating polynomial

    M_Q(t)=L[Q(x,p) exp(q(t))]

is polynomial for every Q. If e>=1 bounds deg_t b and H0=(n-1)(D-1), the same TOTAL-degree estimate as before gives

    deg_t M_Q <= C(e)*deg Q,
    C(e)=max(e,(2n-1)*(1+e*max(H0-1,0))).                   (3)

For a polynomial H with H(0)=0 and origin order nu>=1,

    L[Q H(q)^m]=(H(partial_t)^m M_Q)(0)=0
    whenever nu*m>C(e)*deg Q.                              (4)

This includes every nonlinear polynomial composition in the tuple and all source markers. Zero polynomials are harmless. It is a sufficient bound, not a universal sharpness assertion. Implication1 to2 is immediate. Section2 proves2 to4 and4 to3 without assuming a cutoff uniform in c beforehand.

This classifies JOINT SAFETY of the whole augmented algebra in this cotangent setting. It does not classify which individual nonlinear H(q) may be safe when the joint algebra fails; Section3 gives positive and negative examples in one failing algebra.

## 2. Finite ordinary ideals, nilpotent lifting and a finite line test

Work universally with V having independent matrix entries and coefficient ring R=K[V_ij]. Let b_N(t) be the total-degree-N truncation of g(Vt), and put

    J_N=(all t-coefficients of F(b_N(t))-Vt) in R, N>=1.

The residual has total degree at most DN. The quotient/formal-uniqueness proof of MVINVTAILIDEAL20261116 works word for word with several formal parameters, and yields literal ordinary equality

    J_N=(all coefficients of b_i(t) in total degrees>N).     (5)

Thus J_N vanishes over ANY coefficient algebra iff the inverse restriction has degree at most N. No radical, reduced-point substitution or finite numerical tail is hidden in (5).

If b(t) is polynomial over a field, restrict it to the generic line t=sc with algebraically independent c. The retained affine-line inverse-section degree bound gives deg_s b(sc)<=B=D^(n-1). A nonzero homogeneous component in t remains nonzero on that generic line, so deg_t b<=B. Consequently J_B=0 is the exact FIELD joint-safe locus. It is a finite Zariski-closed condition on frames V.

Over an arbitrary A, if all constant combinations in1.2 are safe, reduce modulo any prime of A and pass to its fraction field. For every c in K^r the inverse line b(sc) is polynomial by GMCALLRINGINV20261116, hence of degree at most B. A polynomial in c over a field extension of the infinite field K that vanishes at every c in K^r is zero. Every homogeneous component of b of degree>B therefore vanishes at that prime. By (5), each of the finitely many generators of J_B(V) is nilpotent. Their ideal I is nilpotent even if A is not Noetherian or its whole nilradical is not nilpotent.

Conversely I^L=0 permits finite Newton correction of b_B, now with total t-degree. With

    H0=(n-1)(D-1), M=H0+D,
    E_L=B*M^ceil(log_2 L),

the retained polynomial inverse-Jacobian iteration gives an exact b(t) of degree at most E_L. Its error ideal squares at each step. Thus1.4 implies1.3, and (3)--(4) give an explicit common all-marker bound. The joint-safe tuple functor on discrete A is the formal-completion functor along J_B, and

    J_(E_L) contained in J_B^L                              (6)

as ordinary coefficient ideals. No general holonomic field-stratification or symbolic-power theorem is needed here.

Even the infinite list of constant combinations has a finite replacement. Choose any set G of DB+1 distinct elements of K. Then

    joint safety over A iff every q(c), c in G^r, is safe.   (7)

For the nontrivial direction, safe line phases make every coefficient of F(b_B(sc))-sVc nilpotent, by the scalar theorem. Taking all c in the finite grid, their coefficient ideals have nilpotent sum I_grid. In particular every entry of the polynomial residual R(t)=F(b_B(t))-Vt takes values in I_grid on G^r. Its degree in each t_j is at most DB. Tensor-product Vandermonde interpolation over K, whose denominators are units in A, puts every coefficient of R in I_grid. Thus J_B has nilpotent image. If each scalar line residual ideal has nilpotence exponent at most L0, their sum has exponent at most |G|^r*(L0-1)+1, supplying a concrete thickness bound. This is a theoretical finite grid; no large grid is enumerated.

For a change of free frame C in GL_r(A), b_N(VC;t)=b_N(V;Ct), and coefficient ideals of the finite residual agree under the invertible parameter substitution. Hence J_N(VC) and J_N(V) generate the same ideal. On full-rank frame charts the field condition is independent of the basis of the represented subspace; it gives closed equations for joint-safe subspaces on Grassmannian charts. The frame formula suffices over arbitrary rings; no unsupported assertion about globally free projective modules is used.

## 3. Two individually safe generators with an unsafe sum and product

Use the retained normalized foundational map f and W=diag(1,-1,-2). Put

    S=Y1, T=-Y2, Y=Jf^(-T)p.

These are explicit polynomials specified by the already retained adjugate. Their total degrees are12 and10. All joint moments L[S^i T^j] vanish for i+j>0, by full Gaussian-state preservation. Nevertheless the safe subset of this jointly pure-null algebra is neither additive nor multiplicatively closed.

Give source positions (x,y,z) weights (1,-1,-2) and momenta (p1,p2,p3) weights(-1,1,2). The two phases S and T have weights-1 and+1 respectively. A source marker of degree delta has all weight components between-2delta and2delta. The Wick pairing has total weight0, so

    L[Q S^m]=L[Q T^m]=0 for m>2deg Q.

This proof works over every coefficient algebra. The coefficient2 is sharp for both phases, since direct exact contractions give

    L[p3 S^2]=-15, L[z T^2]=-8.                             (8)

Every zero-constant polynomial in S alone, or T alone, is safe by (4), or directly by this same weight separation and its origin order.

The joint inverse slice is explicit. Write target coordinates (u,v,0), and let

    U=(1-2uv)^(-1/2),
    x=(U-1)/v, y=v,
    z=-v^2*(1+3U)/U^2.                                    (9)

All expressions are regular formal series at(0,0); U-1 is divisible by uv. The algebraic relation (1-2uv)U^2=1 and substitution into the three f components give exactly (u,v,0). This extends the earlier diagonal branch, without a new map construction.

For S+T, the earlier fixed linear marker x has nonzero odd-index moments indefinitely. Thus the sum is unsafe. For the product, (1) and x*y=U-1 give

    L[xy S^m T^m]=(m!)^2*[u^m v^m](U-1)
                   =(2m)!/2^m !=0 for every m>=1.           (10)

Hence ST, of degree22, is also pure-null and unsafe, with the fixed quadratic marker xy. The corresponding target linear phases p1,-p2 and target quadratic phase -p1p2 are individually safe; the known noninvertible cotangent pullback creates these failures. The target quadratic phase has cutoff m>deg Q/2 by its momentum origin order. Fixed complex linear covariance conjugation gives the same assertions for standard independent Gaussian inputs, preserving all marker degrees.

There is also an exact two-parameter safety classification. For lambda,mu in any coefficient algebra A, the phase lambda*S+mu*T is safe iff lambda*mu is nilpotent. Over fields, (9) on(u,v)=(lambda*s,mu*s) is polynomial exactly on the two coordinate axes: if lambda*mu!=0, its x branch has infinitely many nonzero coefficients. Over A, polynomiality descends to every prime field, so safety implies lambda*mu is nilpotent. Conversely suppose (lambda*mu)^L=0. In an expansion term S^iT^j with a degree-delta marker, a nonzero moment requires |i-j|<=2delta. If min(i,j)>=L its coefficient lambda^i mu^j is zero. Therefore

    L[Q(lambda*S+mu*T)^m]=0 for m>2deg Q+2L-2,

and, using pure-nullness for constant markers, m>2L*deg Q is a common sufficient linear bound. This is a direct charge proof; it avoids expanding the large field gate. The exact safe-parameter functor is the formal completion of the coordinate cross lambda*mu=0. Individually safe axes do not certify their two-dimensional span, and (1,1) already detects its failure.

Joint failure does not make every mixed nonlinear phase unsafe. For example S*(1+ST) has weight-1 and is safe with coefficient2, despite using both generators. The inverse-subspace criterion describes the safety of the ENTIRE augmented algebra, not membership of each single H(S,T) in its surviving safe subset.

## 4. Exact ordinary cutoff ideals for the product failure line

Put R0=ST and consider e*R0 over K[e]. The phase has weight0. Every source variable has nonzero weight, so all linear-marker moments vanish identically; constants vanish by purity. Thus any monomial with a nonzero mixed moment has degree at least2. For an integer N>=0, m>N*deg Q therefore implies m>=2N+1 whenever the moment is nonzero. All generators of the full ordinary cutoff ideal are divisible by e^(2N+1). At m=2N+1, the fixed marker xy of degree2 and (10) give a nonzero K-unit times that exact power. Hence

    I_N^product=(e^(2N+1)) for every integer N>=0.            (11)

Over any A, eta*R0 is safe exactly when eta is nilpotent. If eta has exact nilpotence order L>=1, its least nonnegative REAL linear cutoff coefficient is (L-1)/2, and its least nonnegative INTEGER coefficient is floor(L/2). Indeed all moments with m>=L vanish, linear-marker moments vanish for every m, and (10) at m=L-1 with marker degree2 prevents any smaller real coefficient when L>=2. For eta=0 the least coefficient is0. The discrete-point completion remains K[[e]], while the individual ordinary filtration differs from that of the earlier linear-witness phase, whose ideals are(e^(N+1)). Cofinality does not identify their levels or sharp slopes.

## 5. Verification and remaining questions

The [checker](verify_joint_safe_algebras.py) verifies the two-parameter algebraic inverse relation, Euler charges, independent sparse mixed contractions, the sharp scalar witnesses, small product moments and small multivariate ordinary gates retaining nilpotents. It implements guarded multivariate composition-ideal generators; no large grid or degree49 branch is computed. The [receipt](joint_safe_algebras_verification.json) and [source audit](joint_safe_algebras_source_audit.json) bind this evidence. Finite controls check transcription; the written proofs establish all-degree, all-index and arbitrary-ring statements.

Individual nonlinear safety in a jointly failing algebra,full moment-preservation necessity without inverse,actual differential modules and OP-GVC-FORMAL-LINEAR remain open. All previous scalar inverse-section criteria,negative closures,exact scaling ideals,pure prefix94 and protected smaller endpoints remain valid. Existing map,cotangent,weighted moment-map and classical Gaussian-Lagrange mechanisms retain their prior source credits. New assurance is internal only.
