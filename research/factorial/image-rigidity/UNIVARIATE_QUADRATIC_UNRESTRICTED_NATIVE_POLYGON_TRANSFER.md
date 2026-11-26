# Unrestricted prime-unit native polygon transfer for quadratic factorial moments

The proof uses monotone native coefficient heights and two elementary vertex-unit classes. No restriction on the high digits of the moment index is needed.

Let m>=1 be an integer, let C be an integer, and let p be an odd prime with p not dividing C. Define

    F_m(A;C)=Lambda((A x^2+x+C)^m),  Lambda(x^a)=a!.

For0<=i<=m put

    B_i=binom(m,i)(2i)! C^(m-i),  b_i=v_p(B_i),
    H_i=sum_(l=0)^(m-i) binom(m-i,l)(2i+1)_l C^(-l),
    (2i+1)_l=(2i+1)...(2i+l),  (2i+1)_0=1.             (1)

Directly selecting the quadratic and then linear factors gives

    [A^i]F_m=B_i H_i.                                  (2)

The native lower Newton polygon is the lower convex hull of all (i,b_i). Collinear interior points need not be listed as vertices. Set

    r=min(m mod p,(p-1)/2).                            (3)

**Theorem 1 (unrestricted prime-unit transfer).** The actual Newton polygon of F_m agrees with the native polygon at every positive-slope face. Its initial part is confined to indices0 through r and has only nonpositive slopes, allowing rootsA=0 as part of that initial block. If p divides m, then r=0, the constant coefficient is a p-adic unit, and the complete actual and native Newton polygons are equal. These assertions require neither no-carry digits nor a height bound.

This theorem gives an exact native polygon, not a blanket degree-divisibility byp^v_p(m). Carrying native faces can have other slope denominators; their complete degree capacities must be retained.

## 1. The native heights are nondecreasing

For every i<m,

    B_(i+1)/B_i=2(m-i)(2i+1)/C,
    b_(i+1)-b_i=v_p(m-i)+v_p(2i+1)>=0.                (4)

All native slopes are therefore nonnegative. The first positive height increment occurs at i=r, unless r=m: for0<=i<r neither factor in (4) is divisible byp; at i=r one of them is, if r<m. Hence

    b_0=...=b_r=0;  b_i>0 for i>r.                    (5)

The native polygon has initial horizontal edge from0 to r when r>0, followed by strictly positive-slope faces. If r=m there are no such later faces.

Let i<m be any native vertex whose outgoing face has positive slope. Necessarily

    b_(i+1)>b_i.                                       (6)

Otherwise the native point(i+1,b_i) would lie strictly below the outgoing supporting line, contradicting that it supports the lower hull. By(4), every such vertex satisfies at least one of

    m-i=0 modp,   or   2i+1=0 modp.                   (7)

In particular every native vertex after the initial horizontal block lies in these two automatic unit classes. The endpoint r, when r<m, lies in the classes as well. The final endpoint m will be handled directly.

This is stronger than the false assertion that all positive vertices satisfy i=m modp. For example, at m=15,p=3 the native polygon has vertices

    (0,0),(3,2),(12,10),(13,11),(15,14).

Vertex13 is not congruent to m modulo3. It satisfies instead2i+1=27=0 mod3.

## 2. Both automatic vertex classes are protected

Every H_i in(1) is p-integral because C is a p-adic unit. Thus actual coefficient heights are never below their corresponding native points, including when an actual coefficient is zero.

If p divides m-i, then every term with1<=l<p contains p in binom(m-i,l), unless its value is zero: its numerator contains m-i and l! is a unit. Every term with l>=p has a multiple ofp among the consecutive factors(2i+1)_l. Thus

    H_i=1 modp if p divides m-i.                       (8)

If p divides2i+1, every positive-l rising factorial already contains p in its first factor. Hence also

    H_i=1 modp if p divides2i+1.                       (9)

At i=m there are no positive-l summands and H_m=1 exactly. Therefore every vertex bounding a positive native face remains at exactly its native height. All other actual coefficient heights are at least the native heights and hence above every native supporting line. It follows that every positive native face is present, with the same width and slope, in the actual lower hull.

When r<m, the protected point(r,0) is retained and all actual heights are nonnegative. Consequently no positive-slope actual face can end at or before r. The actual lower hull through indices0..r may have negative or horizontal slopes, but joins the unchanged later positive faces at(r,0). If r=m, all actual points still have nonnegative heights and the final endpoint(m,0) is retained, so every actual face is nonpositive. This proves the initial-block assertion including that boundary case.

If p divides m, then r=0 and(8) protects the constant point: H_0=1 modp, so the actual constant is a unit. All native vertices, including both endpoints, are now protected, proving complete polygon equality. The leading coefficient is always(2m)!, nonzero, so full degree is retained.

## 3. Exact factor-degree capacities

List the positive-slope native edges as(width w_j,slope sigma_j), and let d_j be the reduced denominator of sigma_j. The native endpoint coordinates are integers, so d_j divides w_j. Any rational factor of F_m has degree in

    { b+sum_j d_j e_j:
          0<=b<=r, 0<=e_j<=w_j/d_j,
          b,e_j integers }.                           (10)

Here factors and capacities include multiplicity. No claim that these necessary capacities are sufficient is made.

For completeness, the local argument handles a vanishing constant coefficient. Let z be the least nonzero coefficient index of the actual polynomial. It has exactly z rootsA=0 with multiplicity. Since the coefficient at r is nonzero, z<=r. The nonpositive-slope initial faces have combined width r-z, and their nonzero roots have nonnegative valuation. Including the z zero roots, the initial block accounts for exactly r roots. It can contribute an arbitrary number b between0 and r to a factor; no stronger claim about it is needed.

A positive face of width w_j accounts for exactly w_j nonzero roots of valuation -sigma_j, by multiplication of lowest-weight initial forms. Every monic irreducible factor over Q_p has conjugates of a common valuation. Its nonzero constant coefficient has integer valuation equal to its degree times that root valuation, so that degree is divisible by d_j. Its contribution from this face is at most w_j. Summing local factors gives(10). A zero-root factor is part of b. No identification of local factors across different primes is assumed.

When p divides m, r=0 and there is no initial exception. Formula(10) uses the entire exact native polygon without a deformed expansion. It need not reduce to a single p-power degree divisor. At m=6,p=3,C=1, for example, the exact polygon is

    (0,0),(3,2),(6,5),

with face denominators3 and1. The denominator1 face must not be discarded.

## 4. Scope for resonant consecutive moments

At an odd resonance q>=3, take C=-q and m=q-2 orq-1. Any odd prime p dividing either m is prime to C. The p|m part of Theorem1 gives the full native polygon and exact degree capacities for that moment, with every carry allowed. At a generic odd prime p prime toq but not necessarily dividingm, the initial-block form remains valid, with r given by(3). These are additional necessary constraints on a common rational factor of the two moments.

**Corollary 2 (necessary resonant degree intersection).** Suppose a dense quadratic f=ax^2+bx+c over a characteristic-zero field has three consecutive zero factorial moments starting at n>=1. Put q=-c/b=n+2. The separately proved exact resonance, even-constant and binary-cluster results force q odd with at least four set bits in q-1. Put D=2 to the power of the fourth-lowest set-bit position of q-1. There is a positive integer d<=q-2 such that D divides d and

    d belongs to C_(q-2,p) intersect C_(q-1,p)
    for every odd prime p not dividing q,              (11)

where C_(m,p) is the finite degree-capacity set in(10), formed from the native polygon of the mth moment and the initial width r_(m,p). Thus an empty intersection for any selected finite set of such primes proves coprimeness at that fixed resonance. This is a necessary-condition criterion, not an assertion that all intersections are empty.

Indeed the nonzero rational moment polynomials F_(q-2)(A;-q) and F_(q-1)(A;-q) make a common value-zero algebraic over Q. Its minimal polynomial gives a nonconstant rational common factor G of degree d<=q-2. The binary common-factor lemma gives D|d. At each selected oddp prime toq, C=-q is a p-adic unit, so Theorem1 and(10) apply to BOTH moments. This proves(11), including generic p not dividing their moment indices and possible A=0 contributions to the local initial capacity. The dense common root itself is nonzero, but permitting a larger local capacity only weakens the necessary test.

The exact local resonance/even/binary statements used only for this corollary are proved in [Univariate quadratic Strong Factorial resonance](UNIVARIATE_QUADRATIC_STRONG_FACTORIAL_RESONANCE.md), Sections1--5. They are not premises of Theorem1 or(10). No universal coprimeness or full dense quadratic closure is asserted.

The proof uses only the native coefficient ratio, p-integrality and the two protected vertex classes. There is no irreducibility theorem, no no-carry assumption and no height assumption. No full dense quadratic Strong Factorial theorem, full binary quadratic theorem, SFC(2), SFC(3), formal verification or novelty claim is made.

## 5. Supplementary exact controls

The [standard-library checker](verify_quadratic_unrestricted_native_polygon_transfer.py) and [retained finite receipt](quadratic_unrestricted_native_polygon_transfer_verification.json) reconstruct native and deformed coefficients on explicitly finite ranges, compare complete p-divisor polygons and generic positive faces, and check the necessary qualifications. They do not infer the theorem from a prefix. From the directory containing this note run:

```text
python verify_quadratic_unrestricted_native_polygon_transfer.py --output fresh_unrestricted_native_verification.json
```

Optimized Python is refused and an existing output file is never overwritten. The universal result is the written ratio, unit-class and initial-block proof above. No large finite scan of resonances is a premise or conclusion.
