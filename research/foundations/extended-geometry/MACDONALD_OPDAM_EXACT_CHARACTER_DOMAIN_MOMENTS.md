# Half-integer torus moments and the Coxeter discriminant

The compact torus integral evaluated by Opdam gives exact rational moments for the established uniform-character-domain normalizer. The theorem below is conditional on an actual trace and discriminant-coordinate identity; the subsequent application uses the specified split compact rational models. The half-integer evaluation uses the published complex-parameter theorem. The Coxeter normalization is a separate root-theoretic calculation.

## 1. The precise torus integral and the primary theorem

Let R be a REDUCED crystallographic root system of rankr, choose positive roots, and let T be its compact weight-lattice torus with probability Haar measure da. Put

D_R(a)=product_(alpha>0)(1-alpha(a))(1-alpha(a)^(-1))
      =product_(alpha>0)|1-alpha(a)|^2,
I_R(s)=integral_T D_R(a)^s da.                        (1.1)

Off the root hypertori, D_R is positive; use the real logarithm to define complex powers there. Root hypertori have Haar measure zero. For real s>=0 the integral is ordinary nonnegative integration. The same normalized integral is obtained on any isogenous compact root torus: a surjective finite isogeny pushes probability Haar to probability Haar and preserves all root factors. There is no extra lattice-index or Weyl-order factor.

E. M. Opdam, Some applications of hypergeometric shift operators, Inventiones Mathematicae98(1989),1-18, DOI10.1007/BF01388841, Theorem4.1, printed8, evaluates the compact torus integral for a Weyl-invariant COMPLEX multiplicity function k_alpha with Re(k_alpha)>=0, even for a possibly nonreduced root system. The actual integral uses product_(alpha>0)|a^(alpha/2)-a^(-alpha/2)|^(2k_alpha), with probability Haar explicitly normalized in Section2, printed4. The absolute-value expression is well-defined even when the individual half-root characters require a finite cover.

For reducedR and UNIFORM k_alpha=s, write h_alpha=<rho,alpha^vee>=ht(alpha^vee). The theorem specializes directly to

I_R(s)=product_(alpha>0)
  Gamma(1+s(h_alpha+1))Gamma(1+s(h_alpha-1))
  /Gamma(1+s h_alpha)^2,  Re(s)>=0.                  (1.2)

In particular s=1/2 and s=2k+1/2 for every integerk>=0 are inside the published theorem's domain. No Carlson argument, analytic uniqueness assumption, or integer interpolation is added by this note.

In Opdam's proof, shift-operator/Jacobi norm relations give the multiplicity-shift quotient (4.2); reflection of gamma factors produces a periodic holomorphic quotient. Stirling bounds and the torus integral estimate make that periodic quotient constant by Liouville, and its small-parameter limit fixes normalization. Thus this primary theorem already proves noninteger parameters. The shift-operator existence and Jacobi-norm theorems are inputs to that published proof. The Gaussian Macdonald--Mehta identity is a DIFFERENT result in Section6 and is not being substituted for the torus theorem.

## 2. Reduction to invariant degrees, including both root lengths

Let d_1,...,d_r be the invariant degrees of the Weyl group. The classical Macdonald--Bott--Solomon height/Poincare identity, explicitly displayed in Opdam's proof at printed10, is

product_(alpha>0)(1-z^(ht(alpha^vee)+1))/(1-z^ht(alpha^vee))
  =product_i(1-z^d_i)/(1-z).                          (2.1)

The coroot system has the SAME Weyl group and invariant degrees. Thus no simply-laced identification of roots with coroots is required.

Write H_j for the number of positive coroots of heightj. Comparison of the finite products in(2.1) gives

H_j=#{i:d_i-1>=j}.

One can justify the comparison without treating the factors1-z^j as irreducible: any finite identity product_j(1-z^j)^b_j=1 forces all integerb_j zero, by examining the first nonzero coefficient at the least indexj withb_j!=0, and induction. Applying this to(2.1) gives H_1=r and H_j=H_(j-1)-#{i:d_i=j}.

Now put G_j=Gamma(1+s j), G_0=1. The height product in(1.2) becomes

product_j(G_(j+1)G_(j-1)/G_j^2)^H_j
 =product_i product_(j=1)^(d_i-1)(G_(j+1)G_(j-1)/G_j^2)
 =product_i G_(d_i)/(G_1 G_(d_i-1)).

Therefore

I_R(s)=product_i Gamma(1+s d_i)
                  /[Gamma(1+s)Gamma(1+s(d_i-1))],
                  Re(s)>=0.                         (2.2)

This applies to every reduced crystallographic type A,B,C,D,E,F,G, and to products using the union of degree lists. Uniform multiplicity means the SAMEs on long and short roots, so the B/C/F/G cases really are included. For nonreducedBC systems Opdam's general formula still applies, but the additional k_(alpha/2)/2 terms must be retained; the degree-only formula(2.2) is not asserted there. Nor is a compact root-torus assertion made for noncrystallographic reflection groups.

## 3. Exact positive rational seed moments

**Theorem 1 (conditional moment formula).** Under the actual trace and coordinate hypotheses below, every moment mu_k has the positive rational formula (3.2). When the uniform-domain normalizer hypotheses also hold, its every-order marked trace has formula (3.3). For the established split compact character-domain models these formulas apply over every characteristic-zero field. They evaluate the existing pair; they do not construct its model or prove its geometric hypotheses.

Assume an actual normalized trace on a stated rational model has the compact polar expression

lambda(P(x))=I_R(1/2)^(-1) integral_T D_R(a)^(1/2) P(Phi(a)) da,

together with an actual rational polynomial d satisfying

d(Phi(a))=D_R(a).

These are APPLICATION hypotheses supplied by the character-domain proof, not inferred from Opdam. Let q be a rational regular interior parameter with d(q)!=0, and let f(u)=d(q+u)/d(q) after the coordinate translation x-q. Then the normalizer seed moment is

mu_k=lambda(f(x-q)^(2k))
    =I_R(2k+1/2)/[I_R(1/2)d(q)^(2k)].                (3.1)

Substituting(2.2) and using Gamma(a+n)/Gamma(a)=(a)_n yields the fully explicit formula

mu_k=d(q)^(-2k) product_i
 (1+d_i/2)_(2k d_i)
 /[(3/2)_(2k)(1+(d_i-1)/2)_(2k(d_i-1))].             (3.2)

All Pochhammer lengths are integers and all starting parameters are positive rationals. Thus(3.2) is a positive rational for EVERYintegerk>=0, withmu_0=1. Possible pi factors in I_R(1/2) cancel in the ratio; no unknown volume constant remains. This identity is unbounded in k because it comes from Opdam's parameter theorem, not a finite list of computed moments.

For rankr and the actual uniform-domain construction F=p^r(1+p),H=t^r, its exact all-order marked value is therefore the finite CLOSED positive rational sum

Lambda(H F^m)=sum_(h=0)^m binom(m,h) mu_(rm+h),        (3.3)

with every mu given by(3.2). The criterion's purity, sharp monomial bounds, actual common-zero condition and simple-divisor certificate remain separate proved hypotheses/results; Opdam does not supply their model geometry.

## 4. A simply connected Coxeter-lift calculation

Let G be simply connected semisimple with irreducible reduced root systemR. Work over an algebraic closure. Choose simple roots alpha_1,...,alpha_r, corresponding reflectionss_i, and w=s_1...s_r of orderh. Define

beta_i^vee=s_1...s_(i-1)(alpha_i^vee).

Each beta_i^vee equals alpha_i^vee plus an integer combination of earlier simple coroots. Hence the beta_i^vee form a triangular integral basis of the cocharacter lattice, because G is SIMPLY CONNECTED. For t inT, the root-reflection formula telescopes to

w(t)t^-1=product_i beta_i^vee(alpha_i(t)^(-1)).

As these cocharacters form a basis, w(t)=t if and only if alpha_i(t)=1 for everyi. This proves the full group-scheme identity

T^w=Z(G).                                           (4.1)

For ANY lift n in N_G(T) ofw, n^h belongs toT and is fixed byw because it commutes withn. Thus n^h lies inZ(G) and Ad(n)^h=1. No assertion that the lift itself has orderh, and no separate Tits-lift theorem, is needed. For an arbitrary isogenousG one may lift to its simply connected cover; the adjoint calculation and discriminant are unchanged.

The same telescoping formula on cocharacters gives, in the simple-coroot basis,

1-w=B C,

where C is the Cartan matrix with C_ij=<alpha_i,alpha_j^vee> and B has columns beta_i^vee. Since detB=1,

det(1-w)=detC.                                       (4.2)

In particularw has no eigenvalue1. This does not require an exponents table.

## 5. Free root cycles and the discriminant

A Coxeter element is regular: it has an eigenvector outside all root hyperplanes with eigenvalue of orderh. This is the classical regular-Coxeter fact recalled explicitly in Springer, Regular Elements of Finite Reflection Groups, Inventiones Mathematicae25(1974),159-198, Section5; Section5.4 gives the exceptional degree data. The exact root-orbit consequence and its ENTIRE short proof are Proposition4.10(i), printed173: if w^j alpha=alpha and v is such a regular eigenvector with eigenvaluezeta, then (v,alpha)=zeta^j(v,alpha). Regularity makes this pairing nonzero, so h dividesj. Thus every root orbit has lengthh. There are exactlyr orbits, either by the usual identity |R|=rh or directly from the remaining parts of the same proposition. Its part(ii) gives an ordering in which the word length is the orbit count; since rank(1-w)=r and every reflection changes rank by at most1, that count is at leastr. Its part(iii), together with the chosen expression w=s_1...s_r, makes the count at mostr. This also proves |R|=rh at the needed scope.

On the direct sum of the root lines in the ORIGINAL Cartan decomposition, Ad(n) is a weighted permutation with those h-cycles. Its h-th power is1 by Section4; consequently the product of the weights around each cycle is1. The characteristic polynomial of every cycle is z^h-1, so this root-line part has characteristic polynomial(z^h-1)^r. On the original Cartan subalgebra, Ad(n) acts byw, with no fixed vectors by(4.2).

Thus the eigenvalue1 of Ad(n) has multiplicity exactlyr, n is regular semisimple, and removing those r invariant directions gives

det'(1-Ad(n))=h^r det(1-w)=h^r detC.                 (5.1)

For a compact real representative, conjugate n into a compact torus. The remaining adjoint eigenvalues are precisely its root characters, and therefore

d_T(n)=product_(alpha>0)(1-alpha(n))(1-alpha(n)^(-1))
      =|Delta(n)|^2=h^r detC.                        (5.2)

This uses the POSITIVE compact discriminant. If Delta=e^rho product_(alpha>0)(1-e^-alpha), then d_T=(-1)^|R+| Delta^2. In particular literal Delta^2 has a sign change forE7; the positive d_T in(1.1) has none. Even-power normalization in(3.2) is unaffected by that sign convention.

**Corollary 2 (Coxeter normalization).** For an IRREDUCIBLE root system, at the actual Coxeter parameter of the split compact model, d(qstar)=h^r detC. The exact exceptional data are

E6: r6, h12, degrees2,5,6,8,9,12, detC3, d_T(n)=3*12^6;
E7: r7, h18, degrees2,6,8,10,12,14,18, detC2, d_T(n)=2*18^7;
E8: r8, h30, degrees2,8,12,14,18,20,24,30, detC1, d_T(n)=30^8.

Springer Section5.4 tables, printed177, record these degree/order data. The Cartan determinants follow directly from the standard simply-laced E diagrams: the three arms from the branching vertex have lengths (1,2,2), (1,2,3) and (1,2,4). Eliminating each A-arm gives determinant product_i(a_i+1) times [2-sum_i a_i/(a_i+1)], namely 3,2 and1. This is an exact determinant calculation, not an assumed finite-control certificate. For semisimple products, use one Coxeter element in each irreducible factor: d(qstar)=product_nu h_nu^(r_nu)detC_nu. Formula(3.2) uses the union of the factors' degree lists. No single common Coxeter number is substituted into h^r for a product.

Formula(5.2) computes the discriminant at the Coxeter class. It does NOT by itself prove that a chosen real character-coordinate tuple qstar is rational or that it has an actual fiber in a specified quotient. Those are distinct application-model facts. WHEN that actual identification qstar=Phi(n) has been supplied, one can substitute the displayed positive integer d(qstar) into(3.2)-(3.3).

## 6. Application and references

For the ACTUAL split compact character-domain Q-models, the normalized trace, compact polar law and adjoint discriminant identity are proved in [the character-domain theorem](SPLIT_COMPACT_SYMMETRIC_CHARACTER_DOMAIN_PRODUCT_COUNTEREXAMPLES.md), Sections2--5. They give the trace hypothesis of Theorem1. Its specified Coxeter parameter is an actual rational interior value, by Section7 of that theorem. Consequently (3.2) and (3.3), with the value from Corollary2, give fully explicit positive rational marked traces over every characteristic-zero field. The existing central-quotient phase has the same values by injective pullback and normalized trace restriction, as established there. No individual character coordinate is presumed to descend to that quotient.

The normalizer coefficient formula used in (3.3) is [the bounded-domain criterion](UNIFORM_BOUNDED_DOMAIN_NORMALIZER_CRITERION.md), Sections2--3. Purity, sharp thresholds, actual common-zero fiber, nonpower and support remain the separate results of those proofs; the classical integral evaluation does not replace them.

[O] E. M. Opdam, *Some applications of hypergeometric shift operators*, Inventiones Mathematicae98(1989),1--18, DOI10.1007/BF01388841. [Retained primary](sources/opdam-hypergeometric-shift-operators-1989.pdf). The mathematical input is Theorem4.1 with its complex multiplicity domain and normalized Section2 torus convention, together with the height identity displayed in its proof. This is the torus integral theorem; the Gaussian identity in Section6 is separate. Nonreduced systems require the extra half-root terms and are not covered by the degree-only formula here.

[S] T. A. Springer, *Regular Elements of Finite Reflection Groups*, Inventiones Mathematicae25(1974),159--198, DOI10.1007/BF01390173. [Retained primary](sources/springer-regular-elements-1974.pdf). The mathematical inputs are classical regularity of Coxeter elements in Section5, Proposition4.10 and its proof, and the exceptional degree/order data of Section5.4. Lift centrality and det(1-w)=detC are proved in Section4 above.

The formulas are classical specializations and elementary consequences of these inputs. Reduced crystallographic systems of both root lengths are included, with uniform multiplicity. There is no integer-only interpolation, nonreduced degree-only shortcut, arbitrary twisted-form assertion, maximum-discriminant assertion or finite-prefix inference. Rational qstar and an actual fiber are application hypotheses, not consequences of the Coxeter scalar value alone.
