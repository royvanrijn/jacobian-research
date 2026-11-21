# Primitive nonpower Image witnesses in every balanced degree for at least three pairs

Work over every characteristic-zero field K, and over every coefficient-field extension. For n>=3 use the original contraction pairs

    (tau,t), (w,z), (v,y), (xi_i,z_i) for 4<=i<=n.

Put S_n=K[tau,w,v,xi_4,...,xi_n,t,z,y,z_4,...,z_n] and

    R_n=tau*t+w*z+v*y+sum_(i=4)^n xi_i z_i.

Balanced degree D means homogeneous degree D in the first member of every pair together and independently degree D in the second members together; the ordinary total degree is 2D. Define the original polynomial-valued contraction by

    E_n(xi^alpha z^beta)=partial_z^alpha z^beta,
    M_n=sum_i (partial_(z_i)-xi_i) S_n.

Here xi_1=tau,xi_2=w,xi_3=v and z_1=t,z_2=z,z_3=y. A radial-primitive polynomial is one not divisible by R_n in the FULL original ring. An absolutely nonproper-power polynomial cannot equal a nonzero scalar times an r-th power, r>=2, even after a coefficient-field extension.

**Theorem AllBalancedDegreesThreeOrMorePairs.** For EVERY n>=3 and EVERY D>=2 there is an explicit balanced-degree-D, radial-primitive, absolutely nonproper-power F_(n,D), together with the SAME balanced linear marker Q=tau*y for all n,D, such that for EVERY m>=1

    E_n(F_(n,D)^m)=0,
    E_n(Q F_(n,D)^m)
       =(-1)^(m-1)(m+1)!m! (n+2m+1)_(m(D-2)) !=0.       (N1)

Thus every positive power lies in the original Image M_n, while this one fixed marker fails at every positive order. No general polynomial automorphism classification, all-automorphism trace theorem, asymptotic result, finite-prefix certificate or probability integral is a logical premise of the construction.

## 1. The retained degree-two seed, with its exact all-order proof

Take

    F_2=tau*(t-y)*(w*z+v*t),     Q=tau*y.                 (N2)

This is the degree-two seed of the retained preprint `results/three-pair-image-counterexample/main.tex`, Theorem `thm:main`, lines134--145, with its kernel lemma at92--120. That source uses marker y and records the coefficient of t in E_3(yF_2^m). Multiplying its marker by tau differentiates in t and gives the balanced constant used here. The seed and its moment cancellation are credited antecedents. For completeness the needed identities are proved directly next, so no opaque source theorem is imported.

In expanding (t-y)^m let s count the chosen factors -y, and in (wz+vt)^m let j count the chosen factors vt. A term of F_2^m has exponents

    tau^m w^(m-j) v^j t^(m-s+j) z^(m-j) y^s

and coefficient (-1)^s binom(m,s)binom(m,j). Balanced contraction is nonzero exactly when j=s. Its contribution is

    (-1)^s binom(m,s)^2 m!s!(m-s)!
      =(-1)^s binom(m,s)(m!)^2.

Consequently E_n(F_2^m)=(m!)^2 sum_(s=0)^m (-1)^s binom(m,s)=0. The extra n-3 pairs do not appear, so their zero exponents introduce no factor.

For QF_2^m, the tau exponent is m+1 and the y exponent is s+1. The matching condition is j=s+1, where0<=s<=m-1. Its contribution is

    (-1)^s binom(m,s)binom(m,s+1)
                   (m+1)!(s+1)!(m-s-1)!
      =(-1)^s binom(m,s)(m+1)!m!.

The truncated binomial sum is (-1)^(m-1), proving

    E_n(QF_2^m)=(-1)^(m-1)(m+1)!m!.                       (N3)

These calculations are literal in the original 2n-variable polynomial ring, at every order.

## 2. The normalized trace and a direct replica invariance

Let C_n be the balanced subring, with homogeneous piece C_(n,D), and put A_n=C_n/(R_n-1). For H in C_(n,D), direct monomial contraction gives

    E_n(R_n H)=(n+D)E_n(H).                              (N4)

For a matching monomial xi^alpha z^alpha the sum is sum_i(alpha_i+1)alpha!=(D+n)alpha!; for a nonmatching balanced monomial both sides vanish. Thus

    Lambda_n(H)=E_n(H)/(n)_D                             (N5)

extends linearly and compatibly to A_n, with Lambda_n(1)=1. This is a linear functional, not a ring homomorphism. On balanced homogeneous polynomials the contraction is a scalar.

Define the original-ring linear derivation and its invariant coefficient

    L=w partial_v-y partial_z,    h=w*y.                 (N6)

It preserves R_n and balanced degree, kills h and Q, and is locally nilpotent: it only differentiates v,z while its coefficients w,y are fixed. Direct factorial cancellation gives

    Lambda_n(La)=0  for every a in A_n.                  (N7)

Indeed on a monomial the two possible nonzero contractions require the same exponent condition beta=alpha+e_w-e_v. The first contribution is alpha_v beta!, the second beta_z alpha!, and these are equal because beta_z=alpha_w+1 and beta_y=alpha_v-1. All other terms contract to zero. This establishes (N7) before and after passing to R_n=1.

For every integer k>=0, Lh=0 gives

    Lambda_n(h^k La)=Lambda_n(L(h^k a))=0.                (N8)

The derivation h^kL is locally nilpotent and preserves R_n-1. Its finite exponential sigma_k=exp(h^kL) is an automorphism of A_n. Every positive derivative term in that exponential has zero trace by(N8), hence

    Lambda_n(sigma_k a)=Lambda_n(a),
    Lambda_n(Q*(sigma_k a)^m)=Lambda_n(Q*a^m).           (N9)

The second identity uses that sigma_k fixes Q and is an algebra automorphism. This is the special replica argument needed here. It does not depend on the parallel, stronger all-Aut(A_n) trace theorem.

## 3. The explicit degree-D lift

For D=2 use F_(n,2)=F_2. For D>=3 set k=D-2>=1 and define

    F_(n,D)=tau*(t-y)*[
         R_n^k*(w*z+v*t)+(w*y)^k*w*(t-y)].               (N10)

Each summand has balanced degree k+2=D. The exact original derivative identities are

    L F_2=tau*w*(t-y)^2,    L^2 F_2=0.

Therefore the restriction of (N10) to R_n=1 is exactly

    f_D=f_2+h^k Lf_2=exp(h^kL)f_2.                     (N11)

Equations(N3),(N5),(N9) now give

    E_n(F_(n,D)^m)=0,
    E_n(QF_(n,D)^m)
       =[(n)_(mD+1)/(n)_(2m+1)] (-1)^(m-1)(m+1)!m!.

The quotient of rising factorials is (n+2m+1)_(m(D-2)), giving precisely(N1). It is a positive integer before mapping into K and is nonzero in every characteristic-zero extension. The marker Q is unchanged for every n,D and for every order m; it is not selected separately after observing an order.

## 4. Radial primitivity and exact least balanced degree

The homogeneous quadratic R_n is absolutely irreducible: its rank is2n>=6, while a product of two linear forms has rank at most2. It is prime in the full polynomial ring over any coefficient extension.

For D>=3, reduction of(N10) modulo R_n leaves

    tau*w^(k+1)*y^k*(t-y)^2.                            (N12)

Every displayed linear factor has a nonzero class modulo R_n. Since S_n/(R_n) is a domain, their product is nonzero. Thus R_n does not divide F_(n,D). For D=2, R_n cannot divide any linear factor of F_2 and is not proportional to w*z+v*t, so it does not divide the seed either.

This also gives exact least degree D of f_D in the matrix-entry/balanced filtration of A_n. If f_D had a representative of degree at mostD-1, homogenizing that representative to degreeD with powers of R_n would give a balanced degreeD polynomial divisible by R_n. Its difference from F_(n,D) would lie in C_(n,D) intersect (R_n-1)C_n. That intersection is zero: a nonzero product (R_n-1)H has distinct lowest and highest graded components, since C_n is a domain. The resulting R_n divisibility contradicts(N12).

## 5. Absolute nonproper-power status in the full original ring

The factor tau occurs exactly once. For D>=3, the bracket in(N10), specialized at tau=w=0, is

    (v*y+sum_(i=4)^n xi_i z_i)^k *v*t,

which is a nonzero polynomial. Hence the bracket is not divisible by tau. The other factor t-y is also not divisible by tau. The same assertion is immediate for F_2. Thus the tau-adic valuation of every F_(n,D) is exactly1 over every coefficient-field extension.

If F_(n,D)=c G^r with c nonzero and r>=2 in that full original ring, this valuation would equal r times an integer, an impossibility. This proves absolute scalar nonproper-power status directly. No descent from a quotient factorization, generic point or algebraically closed base field is used.

## 6. Literal original Image membership and scope

For all n, E_n(partial_(z_i)H)=partial_(z_i)E_n(H)=E_n(xi_i H), so M_n is contained in ker E_n. Conversely, reducing a monomial's dual degree with xi_i H congruent to partial_(z_i)H modulo the VECTOR subspace M_n gives p-E_n(p) in M_n by induction. E_n restricts to the identity on K[z], so M_n=ker E_n. This argument does not treat M_n as an ideal.

Thus(N1) supplies literal original Image-Mathieu failures, with all the claimed primitivity and scalar-power qualifications, in every balanced degree D>=2 for every n>=3. The statement is about these explicit balanced phases, not all polynomials of those degrees. It gives no degree-one claim, no resolution of the two-pair degrees5/7 frontier, no Jacobian consequence, and no equivalence between finite checks and all-order cancellation.

The retained degree-two source is a deposited/preprint account, not cited here as a journal publication. Its seed cancellation is credited and reproved. All higher-degree propagation and all normalization factors needed above are displayed explicitly. Finite exact original contractions may support these formulas but are not the proof of their universal quantifiers.
