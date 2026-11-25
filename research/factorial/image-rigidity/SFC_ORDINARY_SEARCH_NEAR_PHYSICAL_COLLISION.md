# Ordinary three-monomial search near the physical row-merger collision

For a nonnegative integer row matrix $A$ with three distinct columns, let $U,V,W$ be the associated monomials. Write $M_m(x,y,t)=\Lambda((xU+yV+tW)^m)$ for the ordinary factorial moment. The main search targets a common projective zero of $M_2,M_3,M_4$, allowing arbitrary complex coefficients. None of the ordinary supports checked here has such a zero. The scope is the row families below.

The reference support and squared correlation $q$ come from the [physical row-merger collision](SFC_PHYSICAL_ROW_MERGE_INTERIOR_COLLISION.md). Keep the base rows $H=((5,4,4),(4,5,4),(5,5,4))$. The split block is $P=((5,2,3),(2,5,3))$ and its merged row is $J=(7,7,6)$. The ordinary endpoints are $P,H$ and $J,H$.

The window length is the number of nonzero monomials in the original conjecture. Edo and van den Essen’s Theorem 3.2 treats two distinct monomials in two variables. Rocks’s Theorems 3.3 and 3.10 treat powers of linear forms and powers of sums of own-variable prime powers with a common prime. The mixed rows here fall outside those hypotheses. [Edo and van den Essen, Definition 2.7 and Theorem 3.2](https://arxiv.org/pdf/1304.3956v2), [Rocks, Theorems 3.3 and 3.10](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf).

## Arithmetic of the physical collision

The exact [elimination checker](verify_sfc_physical_collision_obstructions.py) computes the degree-130 $t$-resultant of the cubic and quartic in the original certificate. It repeatedly removes factors shared with $V_2=\mathbb E_q(U-V)^2$, which is strictly positive on the physical interval. Taking the primitive part leaves $F(q)$ of degree 102. Modulo 389 the degree is preserved and repeated Frobenius powering certifies irreducibility. In detail, the reduction divides $q^{389^{102}}-q$ and has constant gcd with $q^{389^j}-q$ for $j\in\{6,34,51\}$. Hence $F$ is irreducible over $\mathbb Q$. The Bernstein coefficients on $[21/25,17/20]$ have one sign variation. Substituting $q=(21/25+(17/20)u)/(1+u)$ and clearing the denominator converts that variation to the power coefficient variation, so Descartes' rule gives at most one interior real root. The original certificate supplies existence and the unique cubic parameter in its field interval.

Two exact obstructions follow.

- The ratio $\mathbb E_q U^2\,\mathbb E_q V^2/(\mathbb E_q UV)^2$ is invariant under independent column rescalings and a common normalization of the second moments. Its numerator and denominator are polynomials in $q$ of degrees 8 and 14. A rational value would give a nonzero polynomial of degree at most 14 vanishing at the degree-102 correlation. Every ordinary integer support has a rational value of the ratio. Thus even the second-moment tensor cannot be reproduced at the physical collision by ordinary rows and these rescalings.
- In the original coordinates the projective point $[1+is:1-is:t]$ first lies on a nonzero rational homogeneous form in degree 19. To see this, put $K=\mathbb Q[q]/(F)$ and reduce in the rank-102 algebra $\mathbb Z_{(389)}[q]/(F)$, using the unit leading coefficient of $F$. Its reduction is a field. An element with nonzero reduction is a unit, by the determinant of multiplication on the free module. Over the residue field the remainder of $Q_4$ by $Q_3$ has degree 2 and the next remainder has degree 1. Both are the greatest possible remainder degrees, so these degrees and the unit leading coefficients lift. The common root $t$ is therefore the root of the linear remainder, and $s^2=M_2/V_2$ also lies in the localized algebra. The 100 powers $t^k(s^2)^j$ with $k+2j\le18$ have rank 100 after reduction. They are independent over $\mathbb Q$, since any rational relation can be cleared to primitive integer coefficients. Dehomogenize a rational form by setting $r=(x+y)/2=1$ and use $z=(x-y)/2=is$. The even powers of $z$ evaluate to a real linear combination of the tested powers. After dividing by $is\ne0$, the odd powers give a second real linear combination. In degree at most 18 the independence forces both coefficient lists to vanish. In degree 19, the 110 even powers lie in a space of dimension 102, so a nonzero rational relation exists and can be homogenized. In particular an ordinary quadratic moment cannot vanish at the exact physical coefficient direction.

These obstructions concern exact transfer of the exhibited correlation or direction. An ordinary counterexample at a different direction remains possible.

## Exact ordinary row scans

For a row $a$ the factor at the monomial indexed by $k=(i,j,l)$ is $(a\cdot k)!$. Multiply the row factors and the multinomial coefficient to obtain each integer coefficient. A prime excludes a support if the reduced homogeneous forms have no common projective zero over the algebraic closure. Indeed a characteristic-zero projective zero can be taken over a number field. At a place above the prime, scale its coordinates to make their minimum valuation zero and reduce.

### Entry perturbations and arbitrary splits

In the first two domains choose each entry independently within the displayed radius from its endpoint. In the third domain choose the three entries of a total row within radius one of $J$, then enumerate every unordered decomposition of that total into two nonnegative rows. Independently perturb each entry of $H$ by at most one. Tuples with repeated columns are omitted. The domains can overlap.

| Domain | Row tuples | Repeated-column tuples | Survived 65521 | Survived both primes |
|---|---:|---:|---:|---:|
| Split endpoint, radius one | 14,348,907 | 11,664 | 296 | 0 |
| Merged endpoint, radius two | 244,140,625 | 649,488 | 3,738 | 0 |
| Arbitrary splits of 27 nearby totals, with base perturbations | 119,082,150 | 199,132 | 2,028 | 0 |

The [sieve](verify_sfc_collision_endpoint_boxes.cpp) normalizes the quadratic as $x^2+B x+C$, where $B,C$ are homogeneous of degrees one and two in $(y,t)$. Its leading coefficient is a product of factorials with arguments below either prime. The cubic and quartic remainders are $Dx+E$ and $Gx+L$. Over a fixed projective pair $[y:t]$, if at least one remainder has nonzero $x$-coefficient, a shared zero must satisfy $DL-EG=0$ and the norms $E^2-BDE+CD^2=L^2-BGL+CG^2=0$. The same equations force $E=L=0$ if $D=G=0$. If $y=t=0$ the pure $x^2$ coefficient excludes the point. The sieve tests the three binary polynomials by univariate gcd and the point at infinity. An [independent implementation](crosscheck_sfc_collision_endpoint_boxes.py) recomputes factorial coefficients and tests degree-seven Macaulay rank for each of the 6,062 first-prime survivor tuples and 245 controls at both primes.

### Swap-symmetric row multiplicities

Use $u$ copies of $P$, $v$ copies of $J$ and base rows $(a,b,c),(b,a,c),(d,d,e)$. Put $B=U-V$ and $R=r(U+V)+tW$. In the coordinates $(r,z,t)$ for the phase $R+zB$, odd powers of $B$ disappear from the moment expansions. The second moment gives $z^2=-\mathbb E R^2/V_2$ when $V_2\ne0$, leaving homogeneous binary forms of degrees three and four as in the physical certificate. If $r=t=0$, the second moment excludes $z\ne0$. The [sieve](verify_sfc_collision_symmetric_rows.cpp) conservatively retains the primes with $V_2=0$ and otherwise checks homogeneous binary coprimality, including $r=0$.

- For $0\le u,v\le24$, $3\le a\le9$, $1\le b<a$, $1\le c\le8$ and $1\le d,e\le9$, the 14,175,000 parameter tuples leave 199 survivors at 65521 and zero after 65519.
- For the exact base $H$ and **all** $u,v\ge0$, use primes 37, 73 and 109. Each split and merged row factor is a nonzero residue, so the powers are periodic in each multiplicity modulo $p-1$. Start with the $36^2$ residue classes modulo 36. For each prime retain a class only if some lift modulo $p-1$ is a survivor. The respective counts are 38, 4 and 0. The four classes after 73 are $(11,16),(23,12),(35,17),(35,32)$. A [standard-library replay](verify_sfc_seed_multiplicities.py) uses homogeneous Sylvester matrices in place of the C++ gcd.
- Over the same set of 22,680 base choices in the bounded box, the C++ sieve applies this projection with 16 primes $p\equiv1\pmod{36}$ and excludes all multiplicities for 22,010 bases. The remaining 670 bases are unresolved for unbounded $u,v$ by this projection. The bounded test covers all base choices for $u,v\le24$.

### Supports ranked by proximity to the physical moments

The [ranking search](search_sfc_collision_moment_signatures.py) retains $H$ and replaces the correlated pair by up to two reflected pairs $(a,b,c),(b,a,c)$, $0\le b<a\le18$, $0\le c\le18$, and an optional row $(d,d,e)$, $0\le d,e\le18$. It scores 1,907,117,875 blocks, taking the pairs in unordered order. The score is Euclidean distance between normalized log-moment vectors. For moments $m_{ijl}$ of total order $m=i+j+l$, the entry is
\[
\log m_{ijl}-\log m_{m00}+\tfrac l2\bigl(\log m_{200}-\log m_{002}\bigr).
\]
This normalization fixes the pure $U$ moments and the pure $W$ second moment at one. Log signatures add for independent row blocks. The numerical target has $q\simeq0.8408936845642393$ and normalized field $t\simeq-2.059696792934789$.

The 2,048 selected tuples all have constant exact cubic–quartic gcd and no zero in the $r=0$ chart. The tuple with smallest absolute normalized fourth moment on the tested real cubic-root branch adds the rows
\[
(11,10,11),\ (10,11,11),\ (12,11,11),\ (11,12,11),\ (18,18,17)
\]
to $H$. It has normalized $t\simeq-2.06493514$ and fourth moment about $-1.81\times10^{-4}$ after imposing the second and third moments. The 2,048 selected tuples also pass exact tests at starting orders 1 and 3 through 12. In those windows the verification reduces the even-in-$z$ equations in the $r=1$ chart to polynomials in $z^2,t$ and takes two eliminants in $t$. A constant gcd excludes the chart. The $r=0,t\ne0$ chart uses the gcd in $z^2$, and a nonzero pure difference coefficient excludes $r=t=0$. Of the 22,528 additional tuple–window checks, 22,526 pass at 65521 and the other two at 65519. The ranking and the residuals are numerical. The projective exclusions apply to the selected tuples.

## Reproduction

Commands are relative to the repository root. The internal runs used Python 3.12.3, python-flint 0.9.0, NumPy 1.26.4, mpmath 1.3.0 and GCC 13.3.0. The two C++ sieves use the standard library. Set the Python executable to one with the indicated packages.

```sh
source_dir=research/factorial/image-rigidity
run_dir=$(mktemp -d)
python "$source_dir/verify_physical_row_merge_interior_collision.py"
python "$source_dir/verify_sfc_physical_collision_obstructions.py"
python "$source_dir/verify_sfc_seed_multiplicities.py"
g++ -std=c++20 -O3 "$source_dir/verify_sfc_collision_endpoint_boxes.cpp" -o "$run_dir/first"
g++ -std=c++20 -O3 -DPRIME=65519 "$source_dir/verify_sfc_collision_endpoint_boxes.cpp" -o "$run_dir/second"
"$run_dir/first" 0 1 > "$run_dir/source"
"$run_dir/first" 1 2 > "$run_dir/merged"
"$run_dir/first" split 1 > "$run_dir/splits"
cat "$run_dir/source" "$run_dir/merged" "$run_dir/splits" | "$run_dir/second" filter
python "$source_dir/crosscheck_sfc_collision_endpoint_boxes.py" \
  "$run_dir/first" "$run_dir/second" "$run_dir/source" "$run_dir/merged" "$run_dir/splits"
g++ -std=c++20 -O3 "$source_dir/verify_sfc_collision_symmetric_rows.cpp" -o "$run_dir/symmetric"
"$run_dir/symmetric" 24 9 9
"$run_dir/symmetric" --infinite
"$run_dir/symmetric" --infinite-box
OPENBLAS_NUM_THREADS=1 python "$source_dir/search_sfc_collision_moment_signatures.py" \
  --bound 18 --keep 2048 --max-start 12
```
