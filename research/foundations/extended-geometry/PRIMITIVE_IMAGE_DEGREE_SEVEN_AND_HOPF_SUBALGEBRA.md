# A primitive degree-seven Image counterexample and a complete Hopf subalgebra

## 1. Explicit counterexample

Work over any field K of characteristic zero in K[xi_1,xi_2,z_1,z_2]. Define
\[
 R=\xi_1z_1+\xi_2z_2,\quad Z=\xi_1z_2,\quad
 W=2\xi_2z_1,\quad T=\xi_1z_1-\xi_2z_2,
 \qquad T^2+2ZW=R^2.
 \tag{1}
\]
Each bilinear has bidegree (1,1). Put
\[
 L=R+Z,\quad
 D_3=2WR^2-(2R+Z)T^2,\quad E_2=2WR-T^2.
 \tag{2}
\]

**Theorem 1 (primitive septic with a linear marker).** The integer polynomial
\[
                         F_7=L^2D_3E_2                              \tag{3}
\]
has balanced degree seven, is not divisible by R, and is not a proper power even over an algebraic closure. For every integer m>=1,
\[
 \mathcal E_2(F_7^m)=0,\qquad
 \mathcal E_2(ZF_7^m)=(7m+2)!\frac{4^m(m!)^2}{(2m+1)!}\ne0.             \tag{4}
\]
Here
\[
 \mathcal E_2(\xi^\alpha z^\beta)=\partial_z^\alpha z^\beta.
\]
Thus F_7 is a radially primitive, nonproper-power two-pair Image counterexample with the same fixed marker Z of balanced degree one at every positive order.

The central construction is a cancellation inside an algebra of jointly pure functions. A quartic seed e and a degree-two function w satisfy a joint moment identity. Their combination e(e+w^2) loses one degree because e+w^2 is cubic on the quadric. This gives degree seven without multiplying a known seed by an affine function or applying an automorphism.

## 2. A joint endpoint identity

Use the quadric and normalized functional
\[
 A_K=K[x,y,t]/(t^2+2xy-1),\qquad
 \Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
 f\left(x,\frac{1-t^2}{2x},t\right)dt.
 \tag{5}
\]
Integration is the rational polynomial functional. Set U=1+x and
\[
 e=\frac Ux(1-t^2U^2)=U\bigl(2y-(2+x)t^2\bigr),\quad
 w=tU,\quad
 g=e+w^2=U(2y-t^2).
 \tag{6}
\]
In particular e,w,g are regular elements of A_K, of degrees four, two and three.

**Lemma 2 (joint pure and marked laws).** For every r>=1 and s>=0,
\[
 \Lambda(e^rw^s)=0,\qquad
 \Lambda(xe^rw^s)=\frac12\int_{-1}^1z^s(1-z^2)^r\,dz.                  \tag{7}
\]
For s>=0 one also has
\[
                     \Lambda(w^s)=\frac12\int_{-1}^1z^s\,dz.           \tag{8}
\]

Indeed let
\[
 J_{r,s}(V)=\frac12\int_{-V}^{V}z^s(1-z^2)^r\,dz,\quad
 c_{r,s}=J_{r,s}(1).
\]
The derivative of J_{r,s} is divisible by (V-1)^r, since both endpoints have that contact order. Therefore
\[
                  (V-1)^{r+1}\mid J_{r,s}(V)-c_{r,s}.                \tag{9}
\]
The change of variable z=tU gives
\[
 \Lambda(e^rw^s)=[x^r]U^{r-1}J_{r,s}(U),\quad
 \Lambda(xe^rw^s)=[x^{r-1}]U^{r-1}J_{r,s}(U).
\]
The endpoint error in (9) contributes to neither coefficient. The polynomial U^{r-1} has coefficient zero at order r and one at order r-1, proving (7). Formula (8) follows directly from w^s=t^sU^s.

Consequently every element of e K[e,w] is pure: each positive power expands into monomials e^r w^s with r>=1. More generally, for every polynomial P in two variables and every m>=1,
\[
 \Lambda\bigl((eP(g,w))^m\bigr)=0,
\]
\[
 \Lambda\bigl(x(eP(g,w))^m\bigr)
    =\frac12\int_{-1}^1\bigl[(1-z^2)P(1,z)\bigr]^m dz.                \tag{10}
\]
This is a joint moment statement; it does not infer purity merely from purity of individual summands.

For P(g,w)=g, the right side is
\[
 c_m=\frac12\int_{-1}^1(1-z^2)^m dz
     =\frac{4^m(m!)^2}{(2m+1)!}.
 \tag{11}
\]
Thus f_7=eg is pure with x-marked value c_m at every order.

## 3. Homogeneous lift, radial cancellation and absolute nonpower

Let
\[
 F_4=LD_3,\qquad B_2=TL,\qquad G_3=LE_2.
\]
They restrict at R=1 to e,w,g. The exact homogeneous identity is
\[
                    F_4+B_2^2=R G_3.                               \tag{12}
\]
Multiplying by F_4 gives F_4^2+F_4B_2^2=R F_7. The homogeneous degree-eight expression loses exactly one radial factor, with F_7=F_4G_3 as in (3).

For any balanced homogeneous polynomial H of degree d,
\[
             \mathcal E_2(H)=(d+1)!\Lambda(H|_{R=1}).                 \tag{13}
\]
This can be checked on original monomials. Phase averaging kills unequal xi and z multi-indices; the equal-index monomial with exponents a,b has normalized beta integral a!b!/(a+b+1)!, while its contraction is a!b!. The identity is rational and therefore holds over every characteristic-zero field. Applying (13) to F_7^m and ZF_7^m proves (4).

At (xi_1,xi_2,z_1,z_2)=(1,1,1,-1),
\[
 R=0,\ Z=-1,\ W=2,\ T=2,\ L=-1,\ D_3=4,\ E_2=-4,
\]
so F_7=-16. Hence R does not divide F_7. Specializing instead xi_1=z_1=1, z_2=0 and xi_2=u gives
\[
                         F_7\longmapsto(4u-2)(4u-1).                \tag{14}
\]
It has simple roots. A nonzero scalar times a polynomial proper power would specialize to a scalar times a proper power or to zero. Thus (14) excludes every proper power, over every extension field.

For completeness the Image interpretation is exact. With Theta_i=xi_i-partial_{z_i}, one has E_2 Theta_i=0. Induction on xi-degree gives h-E_2(h) in the sum of the images of the Theta_i, by replacing one xi_i with partial_{z_i} at each step. Hence
\[
                   \ker\mathcal E_2=\sum_i\operatorname{Im}\Theta_i.
\]
The pure and marked values in (4) therefore prove failure of the Mathieu property for this Image space.

## 4. Every balanced degree at least six

**Theorem 3 (one integer family, one linear marker).** For nonnegative integers a,b define
\[
 A_2=R^2+TL,\qquad
 F_{a,b}=F_4 A_2^aG_3^b,\qquad d=4+2a+3b.
 \tag{15}
\]
Each F_{a,b} has balanced degree d, is R-primitive and is not a proper power over any extension field. For every m>=1,
\[
 \mathcal E_2(F_{a,b}^m)=0,
\]
\[
 \boxed{\ 
 \mathcal E_2(ZF_{a,b}^m)
   =(dm+2)!
       \frac{2^{(a+2)m}((a+1)m)!\,m!}{((a+2)m+1)!}\ne0.\ }           \tag{16}
\]
Thus every balanced degree d>=6, as well as d=4, has such an integer two-pair Image witness with the fixed marker Z.

At R=1 the phase is e(1+w)^a g^b. Formula (10) gives the marked integral
\[
 \frac12\int_{-1}^1(1+z)^{am}(1-z^2)^m dz
   =2^{(a+2)m}\int_0^1u^{(a+1)m}(1-u)^m du
   =\frac{2^{(a+2)m}((a+1)m)!\,m!}{((a+2)m+1)!}.
\]
All these rational values are nonzero; no ordering on K is required.

The R=0 point used above gives
\[
                       F_{a,b}=-4(-2)^a4^b\ne0.
\]
The specialization used in (14) gives
\[
                 F_{a,b}\longmapsto2^a(4u-2)(4u-1)^b,               \tag{17}
\]
with a simple root at u=1/2, proving absolute nonpower. The displayed nonzero homogeneous polynomial has balanced degree d and no radial factor, so this is its primitive balanced degree. For even d>=6 choose a=(d-4)/2,b=0; for odd d>=7 choose a=(d-7)/2,b=1. Degree four is a=b=0.

In particular degree seven is realized, and all prime degrees at least seven are realized. No degree-five existence or exclusion is implied.

## 5. Exact classification inside the joint subalgebra

This section uses one additional classical input: the non-eventual vanishing of interval power moments of any nonzero polynomial, [Müger–Tuset, arXiv:1910.12148v1, Theorem 2.1](https://arxiv.org/abs/1910.12148v1). Over an arbitrary characteristic-zero field this consequence follows by embedding the finitely generated coefficient field into C. The construction and degree-completion Theorems 1 and 3 above are independent of this input.

The highest symbols of w and g in the associated graded ring of A_K are xt and -xt^2. For distinct monomials g^i w^j, the pairs of exponents (i+j,2i+j) are distinct. The subring K[x,t] injects into K[x,y,t]/(t^2+2xy), so no cancellation of these symbols is possible. Therefore g,w are algebraically independent, and
\[
 \mathcal B=K[e,w]=K[g,w]
\]
is a polynomial algebra. A nonzero P(g,w) has quadric degree equal to its weighted degree with weights deg g=3, deg w=2.

**Theorem 4 (purity, ambient safety and degree spectrum).** In the subalgebra B:

1. A phase is pure, or even eventually pure, exactly when it belongs to e B.
2. A phase is safe against every fixed marker from the entire ambient quadric algebra A_K exactly when it belongs to e(g-1) B.
3. For 0!=f=eP(g,w), its least quadric degree is
\[
                         4+\deg_{(g:3,w:2)}P.                       \tag{18}
\]
Consequently the nonzero pure degree spectrum in B is exactly {4,6,7,8,...}; in particular there is no pure degree-five phase in this specified subalgebra.

For any Q(e,w), (7)--(8) give
\[
                  \Lambda(Q(e,w))=\frac12\int_{-1}^1Q(0,z)dz.
\]
The same identity for Q^m and the classical interval theorem prove statement 1 in both directions. If f=eP(g,w) is safe, (10) and that theorem force (1-z^2)P(1,z)=0, hence P(1,z)=0. Polynomial division in K[g,w] gives (g-1)|P.

Conversely suppose f=e(g-1)Q(g,w). The exact identity
\[
                              g-1=e/U                              \tag{19}
\]
gives f^m=e^{2m}U^{-m}Q(e+w^2,w)^m. Expand the last factor into finitely many e^j w^s. Each term has the form U^{-m}e^r w^s with r=2m+j>=2m. Fix an ambient monomial marker x^A y^B t^C. After the substitution z=tU its marked moment is
\[
 2^{-B}[x^{r+B-A}]U^{r-m-C-2B-1}J(U),
\]
where
\[
 J(V)=\frac12\int_{-V}^{V}
       z^{s+C}(1-z^2)^r(V^2-z^2)^B dz.
\]
Replacing endpoints by -1,1 changes J by a multiple of (V-1)^{r+B+1}: at an endpoint the integrand contains v^r(u-v)^B, and integration raises the total order by one. The fixed-endpoint polynomial has degree at most 2B in V. If m>=2D+1, where A+B+C<=D, then the outer U exponent is nonnegative, the requested coefficient lies below the endpoint error, and the fixed-endpoint term has degree at most r-m-C-1, strictly below r+B-A. The marked moment is zero. This bound is uniform in j,s and Q, so it applies to the entire finite expansion and every fixed ambient polynomial marker. This proves statement 2.

Finally the highest symbol of e is -x^2t^2. Multiplying it by the independent weighted top symbols of P(g,w) proves (18). The weights 2 and 3 generate {0,2,3,4,...}, giving statement 3 and its degree-five exclusion. The shifted monomials in Theorem 3 realize every listed degree with unsafe phases. These are classifications inside B; they do not classify the entire ambient polynomial algebra or exclude arbitrary primitive quintics.

## 6. Classical relations and proof scope

The quartic seed and endpoint mechanism are the existing Long/Hopf construction. [Long, arXiv:2609.16178v1, Theorem 3.1](https://arxiv.org/abs/2609.16178v1), gives its coefficient identity; the earlier [arXiv:2607.19012v1](https://arxiv.org/abs/2607.19012v1) gives the related interval–Laurent cancellation. The local two-pair notes retain the original contraction conventions and quartic/sextic witnesses. The joint height-mixed identity (7), radial cancellation (12) and all-degree family are proved here explicitly; a bounded literature search does not establish worldwide priority.

The non-eventual interval theorem is used only for the exact subalgebra converses in Section 5. No support-converse conjecture, low-degree search or finite moment prefix is a premise of Theorems 1 or 3. Finite original-variable and Laurent checks corroborate the displayed identities; the unrestricted statements follow from the written endpoint proofs and the explicit radial/nonpower certificates. Formal verification, independent replay assurance and external human review are unclaimed.

## 7. Reproducible finite checks

The [standard-library checker](verify_primitive_septic_hopf_degree_family.py) expands the original four-variable forms, checks the radial identity, the R=0 values and simple-root specializations, and computes the original contraction directly on selected powers. A separate literal Laurent calculation checks the septic angular moments. Both the signed auxiliary family and the shifted common-Z family are included, together with a corrupted candidate whose first pure contraction is nonzero. The [receipt](primitive_septic_hopf_degree_family_verification.json) records the complete finite scope.

These checks corroborate explicit identities. The all-order statements are the endpoint and beta-integral proofs above, and the degree/primitivity/nonpower proofs are exact polynomial arguments.
