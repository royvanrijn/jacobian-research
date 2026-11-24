# Uniform-domain normalizers and sharp Laurent-marker bounds

A uniform polynomial functional on a bounded open domain gives an exact multivariable endpoint cancellation. It yields a pure Laurent phase with one fixed unit marker at every order and sharp bounds for all Laurent-monomial markers. The uniform moment law and actual regular coordinates are hypotheses; no homogeneous-model application is inferred from this criterion.

## 1. The conditional theorem

Let Ω⊂R^r, r>=1, be ANY bounded nonempty open set. No boundary regularity, semialgebraicity, connectedness, or star-shapedness is assumed. Write V=vol(Ω)>0. Let f∈Q[x_1,...,x_r] vanish on the whole topological boundary and satisfy f(0)=1. In particular0∈Ω is assumed. Let A_Q be a commutative Q-algebra with normalized functional λ and actual regular coordinates x=(x_1,...,x_r) such that for every polynomial P,

    λ(P(x))=V^(-1)∫_Ω P(u)du.

All polynomial moments on the right are therefore rational. Assume λ and the coordinates admit the stated scalar extensions when an all-field conclusion is needed. For absolute nonpower also assume X=Spec A_Q is smooth geometrically integral and the tuple x has a common geometric zero. These geometric hypotheses are independent of an abstract functional identity.

Theorem 1. Under the stated functional hypotheses, put T=1+t and

    p=T t^(-1) f(Tx)^2,
    F=p^r(1+p), H=t^r.

These are actual rational regular functions on X×splitGm, with H a unit. Then for every m>=1,

    (λ⊗CT)(F^m)=0,
    (λ⊗CT)(HF^m)=Σ_(h=0)^m binom(m,h) μ_(rm+h)>0,
    μ_k=V^(-1)∫_Ω f(u)^(2k)du=λ(f(x)^(2k))∈Q_(>0).

Positivity certifies nonzero rational values before any characteristic-zero base extension; no ordering or probability interpretation of an arbitrary field is asserted.

## 2. Fubini zero-flux and the exact endpoint order

We first record an elementary zero-flux lemma. If W=(W_1,...,W_r) is a polynomial vector field vanishing on the whole topological boundary of a bounded open set Ω, then

    integral_Ω div W = 0.                              (Z)

Indeed fix a coordinate i and all the other coordinates y. The line fiber Ω_y is an open subset ofR, hence a countable disjoint union of bounded open intervals (a_j,b_j). Each endpoint (y,a_j),(y,b_j) lies in the topological boundary of Ω: it is approached by points of Ω along the line and cannot be in Ω, because Ω is open. Thus W_i vanishes at both endpoints. The ordinary fundamental theorem of calculus gives

    integral_(a_j,b_j) partial_i W_i(y,s)ds=0.

The derivative is bounded on a box containing Ω. Consequently its integral over every fiber and over Ω is absolutely integrable, so summing the countable intervals and applying Fubini yields integral_Ω partial_i W_i=0. Summing i proves(Z). For r=1 this is the same interval argument without the other coordinates. This proof requires no perimeter, rectifiability or boundary-area assumption.

Define the polynomial

    J_k(T)=V^(-1)integral_(TΩ) f(w)^(2k)dw
          =T^r V^(-1)integral_Ω f(Tu)^(2k)du,

initially for positive real T and then by its polynomial expression. Differentiation of this bounded polynomial integral gives the EXACT algebraic identity

    J_k'(T)=T^(r-1)V^(-1)integral_Ω div_u[u f(Tu)^(2k)]du.   (E)

This is obtained by differentiating T^r and f(Tu)^(2k); it is not an application of a boundary integral formula.

Write the finite polynomial expansion

    u f((1+t)u)^(2k)=sum_(j>=0)t^j W_j(u).

For j<2k, Leibniz's rule leaves at least2k-j factors of f(u) undifferentiated. Equivalently every component of W_j(u) is divisible, in R[u], by f(u)^(2k-j). Since f vanishes on the whole boundary, W_j vanishes there. Lemma(Z) shows integral_Ω div W_j=0 for every j<2k. In(E), multiplication by T^(r-1)=(1+t)^(r-1) cannot produce terms of smaller degree. Therefore

    J_k'(1+t) is divisible by t^(2k),
    J_k(1+t)-J_k(1) is divisible by t^(2k+1).           (R)

This is the exact multivariable endpoint cancellation. It uses neither a nesting relation between TΩ and Ω nor any shape assumption. No Gauss-Green or finite-perimeter theorem is needed.

## 3. The mixed coefficient identity and its domain

Coefficientwise functional substitution gives for every k>=1 and integer ℓ

    (λ⊗CT)(t^ℓ p^k)=[t^(k−ℓ)]T^(k−r)J_k(T), T=1+t.

If ℓ>=−k, the requested degree is at most2k. The endpoint remainder can then be discarded even if k−r is negative, since T is a unit in R[[t]]. Hence the exact identity is

    (λ⊗CT)(t^ℓ p^k)=μ_k binom(k−r,k−ℓ),  ℓ>=−k.

Here binom(n,j)=n(n-1)...(n-j+1)/j! for integers j>=0, including negative n, while binom(n,j)=0 for j<0.

For k>=r the binomial at ℓ0 is zero and the one at ℓr is1. Expanding F^m proves the claims, with a strictly positive sum of nonzero rational μ-values. The unbounded proof is the endpoint identity, not a finite moment prefix.

## 4. Sharp powers and monomial markers

For positive integer c, p^c is pure IF AND ONLY IF c>=r. If c<r its first moment is μ_c(−1)^c binom(r−1,c), nonzero. Given c>=r, a fixed monomial t^ℓ is unsafe IF AND ONLY IF ℓ>=r, and survives EVERY positive order of p^c IF AND ONLY IF r<=ℓ<=c. These are the ordinary binomial coefficient bounds once k=cm is large; the negative-ℓ domain ℓ>=−k eventually holds. For F the unsafe fixed monomials are exactly ℓ>=r; among them survival at EVERY order holds exactly for r<=ℓ<=r+1. At m1 only k=r,r+1 occur; for generalm the range k∈[rm,(r+1)m] contains a k>=ℓ whenever r<=ℓ<=r+1, and its binomial and μ_k are positive. Everyℓ<r is eventually safe, so cannot survive every order.

## 5. Actual simple-divisor nonpower and Laurent support

At an ACTUAL geometric common zero x0, p(x0,t)=T/t. At t=−1/2, p=−1 and ∂_t(1+p)=−4. Smoothness gives a prime divisor of1+p of multiplicity1 through that point, while p is a unit there. Thus F is geometrically not a unit times any proper power. At the same point

    F(x0,t)=(1+t^(−1))^r(2+t^(−1)).

The coefficients at0 and−1 are2 and2r+1, so actual Laurent support differences have gcd1. Both conclusions persist after every characteristic-zero field extension under the stated model hypotheses.

## 6. A conditional rational-coordinate recipe

Let D be a compact character image with nonempty interior Ω=int(D), and suppose its boundary has Lebesgue measure zero. If an actual Q-map has uniform law on that interior Ω, choose a rational q in Ω away from a rational boundary polynomial d=0. Translate coordinates y=x−q and take f(y)=d(q+y)/d(q). The law remains uniform, f(0)=1, and the whole boundary vanishes. If q is an actual real image point, its real fiber supplies the geometric zero. The abstract theorem does not itself establish this image-point condition, the rational character map, or the domain law. A Coxeter representative can potentially specify q explicitly, but is not necessary for existence of a rational interior parameter.

## 7. Method context

The construction uses the normalizer mechanism of C. D. Long, A Classification of Mathieu's Second Conjecture for Reductive Spherical Pairs, versioned source [revision e9edac2](https://github.com/long-mathematics/mathieu-second-conjecture-spherical-pairs/blob/e9edac262b4b64eb68feb0e6e80496cefe4a1dee/mathieu_second_conjecture_spherical_pairs.tex), Section5, Proposition5.1, equations5.1-5.3. The polynomial square-root-density endpoint circuit and its precise coefficient domain are also recorded in [the classical compact-product proof](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Section4. These are method context. The multivariable Fubini cancellation, exact binomial identity, both sharp directions and actual geometric conclusions are proved here; neither an existing product theorem nor a finite moment prefix is a premise. The criterion does not assert a model classification or priority for the normalizer circuit.
