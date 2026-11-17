// Polynomial composition rigidity via critical values
#set page(paper: "a4", margin: (x: 1.1in, y: 1.1in), numbering: "1")
#set text(size: 11pt)
#set par(justify: true, leading: 0.62em, first-line-indent: 1.2em)
#set heading(numbering: "1.")
#show heading: it => { v(0.6em); it; v(0.3em) }
#set math.equation(numbering: "(1)")
#show ref: it => {
  let el = it.element
  if el != none and el.func() == math.equation {
    link(el.location(), numbering(el.numbering, ..counter(math.equation).at(el.location())))
  } else { it }
}

#let deg = math.op("deg")
#let CV = math.op("CV")
#let ord = math.op("ord")
#let un(body) = [#set math.equation(numbering: none); #body]

#let thm(kind, num, name: none, body) = block(width: 100%, above: 1.1em, below: 1.1em)[
  *#kind #num*#if name != none [ (#name)]. #emph(body)
]
#let rem(kind, num, name: none, body) = block(width: 100%, above: 1.1em, below: 1.1em)[
  *#kind #num*#if name != none [ (#name)]. #body
]
#let proof(body) = block(width: 100%, above: 0.8em, below: 1.1em)[_Proof._ #body #h(1fr) $square$]

#align(center)[
  #text(size: 17pt, weight: "bold")[Polynomial composition rigidity via critical values]
  #v(0.6em)
  #text(size: 12pt)[liqsweep]
  #v(0.2em)
  #text(size: 10pt)[X: \@liqsweep]
  #v(0.2em)
  #text(size: 10pt)[September 2026]
]
#v(0.6em)

#block(inset: (x: 2.2em))[
  #set par(first-line-indent: 0em)
  *Abstract.* Let $a, b in bb(C)[z]$ be polynomials with $a compose b(0) = 0$ and $(a compose b)'(0) = 1$. We prove that either $a compose b = z$ or
  #un[$ ord_0 (a compose b - z) <= deg a + deg b - 1 . $]
  The proof is dynamical: scalar multiplication splits a degenerate fixed point into attracting fixed points without changing the number of distinct critical values, attracting basins are separated by critical values, and critical values are subadditive under composition. The bound is sharp for every pair of degrees. It proves the rigidity conjecture $R(m,n)$ of Furter for all $m, n >= 1$, previously known only for $min(m,n) <= 2$. Consequently the composition-coefficient map $bb(A)^(m+n) -> bb(A)^(m+n)$ is finite and flat of degree $binom(m+n, m)$, the closure of every length-two multidegree stratum of $"Aut"(bb(A)^2_bb(C))$ is the one predicted by Furter, and the Strong Factorial Conjecture of Edo and van den Essen holds for all polynomials $X_1 dots.c X_m (mu_1 X_1 + dots.c + mu_m X_m)$.
]
#v(0.4em)

= Introduction

For a nonzero convergent power series $u$ at the origin, $ord_0 u$ denotes the exponent of its lowest-order term. For a nonconstant polynomial $f in bb(C)[z]$ let
#un[$ CV(f) = { f(zeta) : zeta in bb(C), f'(zeta) = 0 } $]
be the set of distinct finite critical values of $f$; multiplicities are not counted, so $\# CV(f) <= deg f - 1$.

In @Fur15 Furter introduced the following statement about polynomial composition.

#thm("Conjecture", "1.1", name: [$R(m,n)$, Furter])[
  Let $m, n >= 1$ and let
  #un[$ a(z) = z + sum_(i=1)^m a_i z^(i+1), quad b(z) = z + sum_(j=1)^n b_j z^(j+1) $]
  be polynomials with complex coefficients. If $a compose b(z) = z + O(z^(m+n+2))$, then $a = b = z$.
]

Furter proved $R(m,n)$ when $m <= 2$ or $n <= 2$ @Fur15[Theorem A], verified it by Gröbner basis computations for $m <= 8$ and $n <= 5$, and showed that $R(m,n)$ for $d = (m+1, n+1)$ determines the Zariski closure of the set of plane polynomial automorphisms of multidegree $d$ @Fur15[Theorem B]. Edo and van den Essen @EvdE14 showed that the family $R(m)$ of inverse-coefficient statements equivalent to the $R(m,n)$ is a special case of their Strong Factorial Conjecture. Lewis, Perry and Straub @LPS19 proposed an algorithmic approach to the multidegree closure problem and recovered the known cases. The case $R(3)$, equivalently $R(3,n)$ for all $n$, remained open.

The main result of this paper is a contact bound for compositions from which $R(m,n)$ follows for all $m, n$.

#thm("Theorem", "1.2", name: "Composition contact bound")[
  Let $f_1, dots, f_s in bb(C)[z]$ be nonconstant polynomials and let $F = f_s compose dots.c compose f_1$. Suppose $F(0) = 0$, $F'(0) = 1$ and $F != z$. Then
  $ ord_0 (F - z) <= 1 + sum_(i=1)^s (deg f_i - 1). $ <eq:contact>
  In particular, for two factors,
  $ ord_0 (a compose b - z) <= deg a + deg b - 1 . $ <eq:two>
]

Only the normalization of $F$ is used; the factors need not fix the origin. If $deg F = 1$ the normalization forces $F = z$, so the hypothesis $F != z$ implies $deg F >= 2$.

#thm("Corollary", "1.3")[
  $R(m,n)$ holds for all $m, n >= 1$.
]

#proof[
  Let $a, b$ be as in Conjecture 1.1 with $a compose b(z) = z + O(z^(m+n+2))$. If $a compose b != z$, then @eq:two gives
  #un[$ m + n + 2 <= ord_0 (a compose b - z) <= deg a + deg b - 1 <= (m+1) + (n+1) - 1 = m + n + 1, $]
  a contradiction. Hence $a compose b = z$. Degrees multiply under composition, so $deg a = deg b = 1$; the normalizations then give $a = b = z$.
]

The bound @eq:two is attained for every pair of degrees (Corollary 6.2). Theorem 1.2 holds over every field of characteristic zero and fails in every positive characteristic (Remarks 4.2 and 4.3).

The proof of Theorem 1.2 occupies Sections 2 to 4. Let $F(z) = z + c z^(k+1) + O(z^(k+2))$ with $c != 0$. Multiplying $F$ by a scalar $1 + delta^k$ with $delta > 0$ small splits the degenerate fixed point at the origin into $k$ distinct attracting fixed points and does not change the number of distinct critical values (Lemma 3.1). The basin of each attracting fixed point contains a critical point, and distinct basins are disjoint and fully invariant, so distinct attracting fixed points require distinct critical values (Lemma 2.1). Therefore $k <= \# CV(F)$. Finally $\# CV$ is subadditive under composition, which gives $k <= sum (deg f_i - 1)$.

Sections 5 to 8 record the consequences.

- Corollary 5.1: the formal inverse of a polynomial $a = z + O(z^2)$ of degree at most $m+1$ with $a != z$ has no block of $m$ consecutive vanishing coefficients. This is Furter's statement $R(m)$.
- Corollary 6.1: the polynomial map $(a_1, dots, a_m, b_1, dots, b_n) |-> (c_1, dots, c_(m+n))$ sending the coefficients of $a, b$ to those of $a compose b$ is finite, flat and surjective of degree $binom(m+n, m)$.
- Corollary 7.2: for all $d_1, d_2 >= 2$ the closure of the set of plane polynomial automorphisms of multidegree $(d_1, d_2)$ is
  #un[$ overline(cal(G)_((d_1,d_2))) = cal(G)_nothing union union.big_(2 <= s <= d_1+d_2-1) cal(G)_((s)) union union.big_(2 <= u <= d_1, 2 <= v <= d_2) cal(G)_((u,v)) . $]
  This is the Polydegree Conjecture of @LPS19 in length two, in the form conjectured by Furter.
- Corollary 8.2: the Strong Factorial Conjecture holds for every polynomial $X_1 dots.c X_m (mu_1 X_1 + dots.c + mu_m X_m)$.

= Attracting fixed points require distinct critical values

The following is the polynomial case of a theorem of Fatou; see @Mil06[Section 8] for rational maps. We include the short proof.

#thm("Lemma", "2.1")[
  Let $P in bb(C)[z]$ have degree at least two. The number of distinct attracting fixed points of $P$ is at most $\# CV(P)$.
]

#proof[
  Let $p$ be a fixed point with multiplier $lambda = P'(p)$, $abs(lambda) < 1$. Its basin
  #un[$ B_p = { z in bb(C) : P^(compose n)(z) -> p } $]
  is open, since it is the union of the preimages of a small disk around $p$ on which $P$ is a contraction, and bounded, since for an escape radius $R$ with $abs(P(z)) > 2 abs(z)$ whenever $abs(z) > R$ every point outside the disk of radius $R$ escapes to infinity.

  We claim that $B_p$ contains a critical point of $P$. If $lambda = 0$, then $p$ itself is critical. Assume $0 < abs(lambda) < 1$ and suppose that no critical point lies in $B_p$.

  Choose a disk $D = D(p, rho) subset B_p$. Every point of $P^(-n)(D)$ lies in $B_p$. The derivative of $P^(compose n)$ does not vanish on $P^(-n)(D)$: a zero would mean that some intermediate iterate $P^(compose i)(z)$, $0 <= i < n$, is a critical point of $P$ whose orbit later enters $D$, hence lies in $B_p$. Since a polynomial is proper, the restriction
  #un[$ P^(compose n) : P^(-n)(D) -> D $]
  is a proper unramified holomorphic map. On each connected component of $P^(-n)(D)$ its image is open and, by properness, closed in $D$, hence equal to $D$; a proper local biholomorphism onto the simply connected disk $D$ is a covering, hence a biholomorphism on each component. The component containing $p$ supplies an inverse branch
  #un[$ u_n : D -> B_p, quad u_n (p) = p, quad P^(compose n) compose u_n = "id"_D, quad (u_n)'(p) = lambda^(-n) . $]
  All $u_n$ take values in $B_p subset { abs(z) <= R }$. Cauchy's estimate on the circle of radius $rho \/ 2$ about $p$ gives $abs((u_n)'(p)) <= 2 R \/ rho$ for every $n$, contradicting $abs(lambda)^(-n) -> oo$.

  Choose a critical point $zeta_p in B_p$. Its image $P(zeta_p)$ is a finite critical value lying in $B_p$, because $B_p$ is invariant under $P$. The basins of distinct fixed points are disjoint, since an orbit cannot converge to two points. Thus distinct attracting fixed points give distinct critical values.
]

= Splitting a degenerate fixed point

#thm("Lemma", "3.1", name: "Contact versus critical values")[
  Let $f in bb(C)[z]$ satisfy
  #un[$ f(z) = z + c z^(k+1) + z^(k+2) r(z), quad c != 0, quad k >= 1, $]
  with $r$ a polynomial. Then
  $ k <= \# CV(f) . $ <eq:kcv>
]

#proof[
  For a small real parameter $delta > 0$ put
  #un[$ g_delta (z) = (1 + delta^k) f(z) . $]
  Scalar multiplication does not change the critical points and multiplies every critical value by the same nonzero scalar, so
  $ CV(g_delta) = (1 + delta^k) CV(f), quad \# CV(g_delta) = \# CV(f) . $ <eq:cvscale>

  We look for nonzero fixed points of $g_delta$ of the form $z = delta w$. The fixed-point equation $g_delta (z) = z$ reads $delta^k z + (1 + delta^k)(c z^(k+1) + z^(k+2) r(z)) = 0$; substituting $z = delta w$ and dividing by $delta^(k+1) w$ gives
  $ H(delta, w) = 1 + (1 + delta^k) c w^k + delta (1 + delta^k) w^(k+1) r(delta w) = 0 . $ <eq:H>
  The function $H$ is holomorphic in $(delta, w)$, including at $delta = 0$, where $H(0, w) = 1 + c w^k$ has $k$ distinct nonzero roots $alpha_1, dots, alpha_k$. Since $partial_w H(0, alpha_j) = k c alpha_j^(k-1) != 0$, the holomorphic implicit function theorem gives distinct branches
  #un[$ w_j (delta) = alpha_j + O(delta), quad z_j (delta) = delta w_j (delta), quad j = 1, dots, k, $]
  and for small $delta > 0$ these are $k$ distinct nonzero fixed points of $g_delta$.

  These fixed points are attracting. For $z != 0$ one has the exact identity
  $ g'_delta (z) = 1 - k delta^k + (1 + delta^k) z^(k+1) (r(z) + z r'(z)) + (k+1) (g_delta (z) - z) / z , $ <eq:mult>
  which follows by expanding $(k+1)(g_delta (z) - z)\/z = (k+1) delta^k + (1 + delta^k)(k+1)(c z^k + z^(k+1) r(z))$ and comparing with $g'_delta (z) = (1 + delta^k)(1 + (k+1) c z^k + (k+2) z^(k+1) r(z) + z^(k+2) r'(z))$. At $z = z_j (delta)$ the last term of @eq:mult vanishes, and since $z_j (delta) = O(delta)$,
  $ g'_delta (z_j (delta)) = 1 - k delta^k + O(delta^(k+1)), $ <eq:multfp>
  uniformly in $j$. Choose $C$ with $abs(g'_delta (z_j (delta)) - 1 + k delta^k) <= C delta^(k+1)$ for all $j$ and small $delta$, and then $delta$ so small that $k delta^k < 1$ and $C delta < k \/ 2$. Then
  #un[$ abs(g'_delta (z_j (delta))) <= 1 - k delta^k + C delta^(k+1) < 1 - k/2 delta^k < 1 . $]
  Thus $g_delta$ has at least $k$ distinct attracting fixed points, and $deg g_delta = deg f >= k + 1 >= 2$. Lemma 2.1 and @eq:cvscale give $k <= \# CV(g_delta) = \# CV(f)$.
]

#rem("Remark", "3.2")[
  The origin itself becomes a repelling fixed point of $g_delta$, with multiplier $1 + delta^k$; the $k$ nonzero branches supply the attracting fixed points. The first-order term $-k delta^k$ in @eq:multfp, with its sign, is what makes the argument work; nothing is inferred from a count of fixed points with multiplicity.
]

#rem("Remark", "3.3")[
  Lemma 3.1 can also be read off from the theory of parabolic fixed points: the fixed point $0$ of $f$ has exactly $k$ attracting petals, each parabolic basin contains a critical point @Mil06[Section 10], and since the multiplier is exactly $1$ the $k$ basins are disjoint and $f$-invariant, so their critical values are distinct. The perturbation argument above uses only Lemma 2.1.
]

= Proof of the composition bound

For nonconstant $a, b in bb(C)[z]$ the chain rule gives
$ CV(a compose b) subset.eq CV(a) union a(CV(b)) . $ <eq:chain>
Indeed, if $(a compose b)'(z) = a'(b(z)) b'(z) = 0$, then either $b'(z) = 0$, in which case $a(b(z)) in a(CV(b))$, or $a'(b(z)) = 0$, in which case $a(b(z)) in CV(a)$. Hence
$ \# CV(a compose b) <= \# CV(a) + \# CV(b) <= (deg a - 1) + (deg b - 1) . $ <eq:subadd>

#proof[
  (of Theorem 1.2.) Write $F(z) = z + c z^(k+1) + O(z^(k+2))$ with $c != 0$; here $k >= 1$ because $F(0) = 0$, $F'(0) = 1$ and $F != z$. Applying @eq:chain repeatedly to $F = f_s compose (f_(s-1) compose dots.c compose f_1)$ gives $\# CV(F) <= sum_i \# CV(f_i)$. Lemma 3.1 then gives
  #un[$ ord_0 (F - z) - 1 = k <= \# CV(F) <= sum_(i=1)^s \# CV(f_i) <= sum_(i=1)^s (deg f_i - 1) . $]
  Adding one proves @eq:contact, and @eq:two is the case $s = 2$.
]

#rem("Remark", "4.1", name: "Sharpness in equal degrees")[
  For $d >= 2$ let $a(z) = z + z^d$ and $b(z) = z - z^d$. Then
  #un[$ a compose b(z) - z = -d z^(2d-1) + O(z^(3d-2)), $]
  so equality holds in @eq:two. Corollary 6.2 gives sharpness for every pair of degrees.
]

#rem("Remark", "4.2", name: "Fields of characteristic zero")[
  Theorem 1.2 and Corollary 1.3 hold over every field $K$ of characteristic zero. The finitely many coefficients of a putative counterexample generate a finitely generated extension of $bb(Q)$, which embeds into $bb(C)$, and the contact order and the degrees are preserved by the embedding.
]

#rem("Remark", "4.3", name: "Positive characteristic")[
  In characteristic $p > 0$ the bound fails: for $a = z + z^p$ and $b = z - z^p$ one has $a compose b = z - z^(p^2)$, of contact order $p^2 > 2p - 1$. Accordingly $R(m,n)$ is a statement about characteristic zero.
]

= Formal inverses

Furter's original formulation is in terms of the formal inverse of a single polynomial @Fur15[Section 1]; see also @EvdE14[Section 2.2].

#thm("Corollary", "5.1", name: [$R(m)$])[
  Let $a(z) = z + sum_(i=1)^m a_i z^(i+1) in bb(C)[z]$ and let $a^(-1)(z) = z + sum_(j >= 1) alpha_j z^(j+1)$ be its formal compositional inverse. If $alpha_(n+1) = dots.c = alpha_(n+m) = 0$ for some $n >= 0$, then $a = z$.
]

#proof[
  Let $b(z) = z + sum_(j=1)^n alpha_j z^(j+1)$. The vanishing block gives $a^(-1)(z) - b(z) = O(z^(m+n+2))$. Substitution into the polynomial $a$ preserves this order, so $a compose b(z) = a compose a^(-1)(z) + O(z^(m+n+2)) = z + O(z^(m+n+2))$. If $n >= 1$, Corollary 1.3 gives $a = z$. If $n = 0$, then $b = z$ and $a(z) = z + O(z^(m+2))$; as $deg a <= m+1$, again $a = z$.
]

Conversely $R(m)$ implies $R(m,n)$ for every $n$ @Fur15[Lemma 2], so Corollaries 1.3 and 5.1 are equivalent statements.

= The composition-coefficient map

Put $M = m + n$. For $a, b$ normalized as in Conjecture 1.1 write
#un[$ a compose b(z) = z + sum_(j >= 1) c_j z^(j+1) . $]
Let $S = bb(C)[a_1, dots, a_m, b_1, dots, b_n]$, graded by $"wt"(a_i) = i$ and $"wt"(b_j) = j$. Then $c_j in S$ is homogeneous of weight $j$, and $c_j = 0$ for $j > (m+1)(n+1) - 1$. Consider the polynomial map
#un[$ C_(m,n) : bb(A)^M_bb(C) -> bb(A)^M_bb(C), quad (a_1, dots, a_m, b_1, dots, b_n) |-> (c_1, dots, c_M) . $]
Furter observed that $R(m,n)$ is equivalent to $c_1, dots, c_M$ being a homogeneous system of parameters of $S$, and to $C_(m,n)$ being quasi-finite @Fur15[Introduction and Lemma 11]. Corollary 1.3 therefore yields the following.

#thm("Corollary", "6.1")[
  The map $C_(m,n)$ is finite, flat and surjective, of degree $binom(m+n, m)$. Every fiber has this scheme-theoretic length, and a generic fiber consists of $binom(m+n, m)$ distinct points.
]

#proof[
  By Corollary 1.3 the only common zero of $c_1, dots, c_M$ in $bb(C)^M$ is the origin, so
  #un[$ B = S \/ (c_1, dots, c_M) $]
  is a finite-dimensional graded $bb(C)$-algebra and $c_1, dots, c_M$ is a homogeneous system of parameters of $S$ @Sta96[Chapter I, Section 5]. Let $T = bb(C)[c_1, dots, c_M] subset S$. Lifting a homogeneous $bb(C)$-basis of $B$ to $S$ and inducting on weighted degree shows that these lifts generate $S$ as a $T$-module. Thus $T subset S$ is finite, $dim T = dim S = M$, and since $T$ is generated by $M$ elements they are algebraically independent: $T$ is a polynomial ring with generators of weights $1, dots, M$.

  The polynomial ring $S$ is Cohen–Macaulay, so its homogeneous system of parameters $c_1, dots, c_M$ is a regular sequence @Sta96[Chapter I, Section 5] @BH93[Section 2.1]. Hence the Hilbert series of $B$ is
  $ H_B (t) = (product_(k=1)^M (1 - t^k)) / (product_(i=1)^m (1 - t^i) product_(j=1)^n (1 - t^j)) . $ <eq:hilb>
  The lifted basis of $B$ defines a surjection from a graded free $T$-module $Phi$ onto $S$, where $Phi$ has one generator in each degree of the basis of $B$. By @eq:hilb the Hilbert series of $Phi$ equals $H_B (t) \/ product_(k=1)^M (1 - t^k) = H_S (t)$, so the surjection is an isomorphism and $S$ is a finite free $T$-module. Its rank is
  #un[$ dim_bb(C) B = lim_(t -> 1) H_B (t) = (M!) / (m! thin n!) = binom(m+n, m) . $]
  Freeness gives flatness and the constancy of the fiber length. A finite injective ring map is surjective on spectra, so $C_(m,n)$ is surjective. Since the characteristic is zero the finite extension of function fields is separable, so $C_(m,n)$ is generically étale and a generic fiber is reduced, hence consists of $binom(m+n,m)$ distinct points.
]

#thm("Corollary", "6.2", name: "Sharpness for all degree pairs")[
  For every $m, n >= 1$ there exist $a, b in bb(C)[z]$ with $deg a = m+1$, $deg b = n+1$, $a(0) = b(0) = 0$, $a'(0) = b'(0) = 1$ and $ord_0 (a compose b - z) = m + n + 1$.
]

#proof[
  By surjectivity in Corollary 6.1 there is a point $(a, b)$ with $(c_1, dots, c_M) = (0, dots, 0, 1)$, that is, $a compose b(z) = z + z^(M+1) + O(z^(M+2))$, of contact order $M + 1 = m + n + 1$. Theorem 1.2 forces $deg a + deg b >= M + 2 = (m+1) + (n+1)$. Since $deg a <= m+1$ and $deg b <= n+1$, both bounds are attained.
]

= Closures of multidegree strata

Let $cal(G) = "Aut"(bb(A)^2_bb(C))$ be the group of polynomial automorphisms of the complex affine plane, let $cal(A) subset cal(G)$ be the affine subgroup and $cal(B) subset cal(G)$ the triangular subgroup of automorphisms $(x, y) |-> (alpha x + p(y), beta y + gamma)$ with $alpha beta != 0$. By the theorem of Jung and van der Kulk @Jun42 @vdK53, $cal(G)$ is the amalgamated product of $cal(A)$ and $cal(B)$ over $cal(A) inter cal(B)$. Every $g in cal(G) without cal(A)$ can be written as $g = a_0 compose b_1 compose a_1 compose dots.c compose b_l compose a_l$ with $a_i in cal(A)$ and $b_j in cal(B) without cal(A)$, and the sequence $(deg b_1, dots, deg b_l)$ of integers $>= 2$ depends only on $g$. It is the _multidegree_ of $g$ @Fur97 @Fur15, also called the _polydegree_ @LPS19; affine automorphisms have multidegree $nothing$. For a multidegree $d$ let $cal(G)_d subset cal(G)$ be the set of automorphisms of multidegree $d$. The group $cal(G)$ is an ind-variety, the union of the algebraic varieties $cal(G)_(<= N)$ of automorphisms of degree at most $N$, and closures are taken in this topology; since $cal(G)_d$ consists of automorphisms of degree $d_1 dots.c d_l$, its closure is its closure in $cal(G)_(<= d_1 dots.c d_l)$.

Furter defined the partial order $prec.eq$ on multidegrees generated by the relations

+ $nothing prec.eq d$ for every $d$;
+ $(d_1, dots, d_k) prec.eq (e_1, dots, e_k)$ whenever $d_j <= e_j$ for all $j$;
+ $(d_1, dots, d_(j-1), d_j + d_(j+1) - 1, d_(j+2), dots, d_k) prec.eq (d_1, dots, d_k)$ for $1 <= j <= k-1$,

and proved the following.

#thm("Theorem", "7.1", name: [Furter, @Fur15[Theorem B]])[
  Let $m, n >= 1$ and $d = (m+1, n+1)$. If $R(m,n)$ holds, then
  #un[$ overline(cal(G)_d) = union.big_(e prec.eq d) cal(G)_e . $]
]

#thm("Corollary", "7.2")[
  For all integers $d_1, d_2 >= 2$,
  $ overline(cal(G)_((d_1,d_2))) = cal(G)_nothing union union.big_(s=2)^(d_1+d_2-1) cal(G)_((s)) union union.big_(2 <= u <= d_1, thin 2 <= v <= d_2) cal(G)_((u,v)) . $ <eq:closure>
]

#proof[
  Corollary 1.3 supplies $R(d_1 - 1, d_2 - 1)$, so Theorem 7.1 applies with $d = (d_1, d_2)$. The multidegrees $e prec.eq (d_1, d_2)$ are exactly: $nothing$; the pairs $(u, v)$ with $2 <= u <= d_1$ and $2 <= v <= d_2$, by relation 2; and the singletons $(s)$ with $2 <= s <= d_1 + d_2 - 1$, obtained from such a pair by relation 3 followed by relation 2. No other multidegree is reachable: none of the relations increases the length, a multidegree of length two below $(d_1, d_2)$ arises only through relation 2, and relations 2 and 3 do not increase the sum of the entries minus the length, which for a singleton $(s)$ gives $s - 1 <= d_1 + d_2 - 2$. This gives @eq:closure.
]

Corollary 7.2 is the length-two case of the equality $overline(cal(G)_d) = union.big_(e prec.eq d) cal(G)_e$, known in length one by @Fur97 and previously in length two when $d_1 <= 3$ or $d_2 <= 3$ @Fur15[Theorem C]. In length three the equality is false @EF04, so Corollary 7.2 completes the picture for the lengths in which the closure is described by $prec.eq$. As Furter remarked, the length-two equality implies by induction on the length that $union.big_(e prec.eq d) cal(G)_e subset.eq overline(cal(G)_d)$ for every multidegree $d$. Corollary 7.2 also settles the Polydegree Conjecture of @LPS19, which is @eq:closure in their notation.

= The Strong Factorial Conjecture

Let $bb(C)^([m]) = bb(C)[X_1, dots, X_m]$ and let $L : bb(C)^([m]) -> bb(C)$ be the linear map with $L(X_1^(l_1) dots.c X_m^(l_m)) = l_1 ! dots.c l_m !$. For a nonzero polynomial $f$ let $N(f)$ be its number of monomials. The Factorial Conjecture asserts that $L(f^k) = 0$ for all $k >= 1$ forces $f = 0$. Edo and van den Essen @EvdE14 formulated the stronger statement below and connected it to Furter's conjecture.

#thm("Conjecture", "8.1", name: [Strong Factorial Conjecture, @EvdE14])[
  Let $f in bb(C)^([m])$ be nonzero. For every $n >= 1$ there is $k$ with $n <= k <= n + N(f) - 1$ and $L(f^k) != 0$.
]

Let $E^([m]) = { X_1 dots.c X_m (mu_1 X_1 + dots.c + mu_m X_m) : mu_1, dots, mu_m in bb(C) }$. Theorem 2.25(5) of @EvdE14 states that Conjecture 8.1 holds for every $f in E^([m])$ if and only if $R(m')$ holds for all $m' <= m$. Corollary 5.1 therefore gives the following.

#thm("Corollary", "8.2")[
  For every $m >= 1$ the Strong Factorial Conjecture holds for every polynomial in $E^([m])$: if $f = X_1 dots.c X_m (mu_1 X_1 + dots.c + mu_m X_m)$ is nonzero, then for every $n >= 1$ there is $k$ with $n <= k <= n + N(f) - 1$ and $L(f^k) != 0$.
]

Previously the windows $n = 1, 2, 3$ were known, corresponding to $R(1)$ and $R(2)$, and with them the Factorial Conjecture on $E^([m])$ @EvdE14[Corollary 2.28]. Corollary 8.2 gives every window.

#v(1em)
*Acknowledgments.* The exposition was prepared with the assistance of GPT-6 Astra. The small cases of Corollaries 1.3 and 6.1 with $m + n <= 6$ were verified by Gröbner basis computations in SymPy, in agreement with Furter's computations.

#v(0.5em)
#bibliography("refs.yml", title: "References", style: "ieee")
