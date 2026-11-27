# Quadratic torsion moments and generalized dihedral field counts

The signed quadratic moments give discriminant counts for extensions
with generalized dihedral Galois group. Relative ramification can be
summed uniformly. Away from primes dividing the odd kernel order,
a rational prime in the relative conductor accounts for total residue
degree two in the quadratic subfield.

## The transfer for an odd abelian kernel

Fix an abelian group \(A\) of odd order \(n>1\). Let \(G=A\rtimes C_2\)
with inversion action, and write \(N_A(Y)\) for the number of
\(\mathbb Q\)-isomorphism classes of degree-\(n\) fields of absolute
discriminant at most \(Y\) with normal closure group \(G\).
Write \(\widetilde N_A(Y)\) for the corresponding count of the degree-\(2n\)
Galois fields, with their own absolute discriminants at most \(Y\).
For a quadratic field \(K\) set
\(W_A(K)=|\operatorname{Hom}(\mathrm{Cl}(K),A)|\).

**Theorem (generalized dihedral moment transfer).** If \(\Theta\ge1\) and
\[
 \sum_{K:\,|D_K|\le X} W_A(K)\ll_{A,\eta} X^{\Theta+\eta}
 \quad\text{for every }\eta>0,
\]
then, for every \(\varepsilon>0\),
\[
             N_A(Y)\ll_{A,\varepsilon}Y^{2\Theta/(n-1)+\varepsilon},
             \qquad
             \widetilde N_A(Y)\ll_{A,\varepsilon}Y^{\Theta/n+\varepsilon}.
                                                                    \tag{1}
\]
The implication also holds with \(K\) restricted to either sign. The
prime cyclic kernel case is the counting architecture of
[Klüners, equation (1) and proof of Theorem 2.5](https://jtnb.centre-mersenne.org/item/10.5802/jtnb.561.pdf).
Here the fixed kernel can have arbitrary odd abelian structure.

Let \(L/\mathbb Q\) be a normal closure with this group, and
\(K=L^A\). It is the unique quadratic subfield since any index-two
subgroup contains \(A\). All order-two subgroups are conjugate: a
reflection acts on the set \(A\) by \(x\mapsto a-x\), and the map
\(x\mapsto 2x\) on \(A\) is invertible. Thus \(L\) supplies one
isomorphism class of degree-\(n\) fields \(F\) by taking a fixed field of
a reflection subgroup. Its normal closure is \(L\), since no such subgroup
is normal for \(n>1\). In permutation characters one has
\[
 \chi_{\mathrm{reg}}=2\chi_{G/C_2}+\chi_{G/A}-2\mathbf1 .
\]
At a reflection the two nonconstant characters on the right have values
\(1\) and \(0\), at a nontrivial element of \(A\) they have values \(0\)
and \(2\), and the dimensions agree. By the Artin conductor–discriminant
formula, applied prime by prime including wild inertia,
\[
 |D_L|=|D_F|^2|D_K|,\qquad
 |D_L|=|D_K|^n N_{K/\mathbb Q}(\mathfrak d_{L/K}).                 \tag{2}
\]
The second identity is the discriminant tower formula.
[Jones–Roberts, §3C, equation (3-5)](https://msp.org/ant/2014/8-3/ant-v8-n3-p04-p.pdf)
states the linear discriminant formula for permutation characters.
If \(K\) is real, \(L/K\) cannot ramify at real places since its degree
is odd, so \(F\) is totally real. If \(K\) is imaginary, complex
conjugation acts as a reflection on the \(n\) embeddings of \(F\);
its signature is \((1,(n-1)/2)\).

Let \(b\) be the squarefree product of rational primes \(p\nmid n\)
lying under ramified primes of \(L/K\). No such \(p\) can be ramified
in \(K\). In fact if \(p\) ramifies in \(K\), its unique prime has
residue field \(\mathbb F_p\), on which quadratic conjugation is
trivial. Since \(p\nmid n\), the image of the local units in
\(\operatorname{Gal}(L/K)\) factors through the residue field units.
The local Artin map is equivariant under conjugation, so that image
would be fixed by inversion, and is therefore trivial in the odd group.
For each \(p\mid b\), either both primes above \(p\) ramify in \(L/K\),
or \(p\) is inert in \(K\) and its prime has residue degree two.
This follows from conjugation symmetry and unramifiedness in \(K\).

Let \(p_0\) be the smallest prime divisor of \(n\).
For a prime \(\mathfrak p\) of \(K\) above \(p\mid b\), the tame inertia
of \(L/K\) has order \(e\ge p_0\). The relative discriminant exponent is
\(n(1-1/e)\), since each local different has tame exponent \(e-1\).
Total residue degree two above \(p\) gives
\[
 N_{K/\mathbb Q}(\mathfrak d_{L/K})\ge
 b^{2n(1-1/p_0)}.                                                \tag{3}
\]
Thus the degree-\(n\) discriminant bound uses
\(|D_K|^{(n-1)/2}b^{n(1-1/p_0)}\le Y\); for the regular degree
use \(|D_K|^n b^{2n(1-1/p_0)}\le Y\).

For fixed \(K,b\), global class field theory counts these \(L/K\)
among \(A\)-quotients of a single ray class group. At each
\(\mathfrak p\mid b\) use conductor exponent \(1\).
At each \(\mathfrak p\) above \(p\mid n\) use exponent
\(1+4v_p(n)\), independent of \(K\). This suffices because every
local unit congruent to \(1\) to that depth has an \(n\)th root
by Hensel at \(1\), with valuation \(>2v_\mathfrak p(n)\).
There are no real-place conditions for an odd target group.
In the ray-to-class exact sequence, every kernel factor at a tame
prime is a quotient of a cyclic residue unit group, with at most
\(n\) homomorphisms to \(A\). The factors at \(p\mid n\) have total
bounded size depending only on \(A\), since the field is quadratic and
the exponents are bounded. The restriction map on homomorphisms has
kernel consisting of homomorphisms from the ordinary class group.
Thus the number of \(L/K\) is
\[
                         \ll_A n^{2\omega(b)} W_A(K).             \tag{4}
\]
Counting all maps covers the possible surjections and uses no lower
bound on the number of maps per field.
[Stevenhagen, §2, equations (2.1)–(2.3)](https://pmb.centre-mersenne.org/item/10.5802/pmb.a-56.pdf)

Combining (3)–(4) with the moment assumption, the sum over \(b\) in the
degree-\(n\) count is bounded by a constant times
\[
 Y^{2(\Theta+\eta)/(n-1)}
       \sum_{b\ge1} n^{2\omega(b)}
          b^{-\frac{2n}{n-1}(1-1/p_0)(\Theta+\eta)} .
\]
For the regular count the outer exponent is
\((\Theta+\eta)/n\) and the exponent in the series is
\(2(1-1/p_0)(\Theta+\eta)\). Both series converge because
\(p_0\ge3\), \(\Theta\ge1\), and \(n^{2\omega(b)}\ll_{A,\eta} b^\eta\)
for every \(\eta>0\). This is the uniform step that permits variation
in the relative conductor.

## Consequences at orders 75, 125, 625 and 3125

If \(A\) has \(p\)-part of order \(p^r\), its contribution to \(W_A(K)\)
is at most \(h_p(K)^r\). Indeed the cyclic factors \(C_{p^e}\)
give \(|\mathrm{Cl}(K)[p^e]|\le h_p(K)^e\), by iterating the
multiplication-by-\(p\) map. This also handles non-elementary kernels of
the same order. Substituting the [signed fifth-torsion moments](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md#signed-tails-and-moments)
gives the following exponents, each with an arbitrarily small addition.

| Order of the kernel \(A\) | Moment exponent \(\Theta\) | \(N_A(Y)\) | \(\widetilde N_A(Y)\) |
| ---: | ---: | ---: | ---: |
| \(125\) | \(7/4\) | \(7/248\) | \(7/500\) |
| \(625\) | \(53/26\) | \(53/8112\) | \(53/16250\) |
| \(3125\) | \(5/2\) | \(5/3124\) | \(1/1250\) |

## A mixed moment from the product torsion exponent

**Theorem (mixed moment from the class-group product).** Fix a prime
\(\ell\ge5\), put \(e=3\ell\), and write
\[
                  \Theta_\ell=\frac{2e\ell+e-1}{e(\ell+2)+1}.
\]
Then for either sign the signed quadratic mixed moment is
\[
     \sum_{K:\,|D_K|\le X} h_3(K)h_\ell(K)^2
                 \ll_{\ell,\varepsilon} X^{\Theta_\ell+\varepsilon}. \tag{5}
\]
In particular for a fixed abelian kernel of order \(3\ell^2\) the
transfer (1) applies with this \(\Theta_\ell\).

The product \(h_3 h_\ell\) is exactly \(h_e\) in one class group, since
\(3\) and \(\ell\) are coprime. All three torsion quantities in (5) are
at most the class number, \(O_\varepsilon(X^{1/2+\varepsilon})\).
Use dyadic boxes in \(h_\ell,h_e\), with lower endpoints \(X^b,X^s\),
respectively, inside a dyadic field window. Each summand in a box
has size \(O(X^{s+b})\). The first moment of \(3\)-torsion bounds the
number of fields therein by \(O(X^{1-s+b+\varepsilon})\), since
\(h_3=h_e/h_\ell\). [Davenport–Heilbronn, as stated in
Bhargava–Shankar–Tsimerman, Theorem 2](sources/bst-2013-davenport-heilbronn.pdf)

For \(k=\ell,e\), apply the pair arm of the
[signed higher odd tail](ABC_COHERENT_PRIME_PATH_HIGHER_ODD_TORSION.md#coherent-witnesses-with-a-bounded-unit-logarithm)
at threshold \(H=X^u\). The constant exceptional term is absorbed by
\(X/H\) in this class-number range. If \(u\ge(k-2)/(2k)\), the
pair tail is bounded by \(O(X^{1-u+\varepsilon})\), which gives at most
\(X^{3/2+O(\varepsilon)}\) for an entire box in (5) by using either
\(u=b\) or \(u=s\). Below the respective thresholds the pair arm has
exponent \(k/2-(k+1)u\). The first moment and the two pair bounds
give the following three exponents for the box, up to subpower factors:
\[
                   1+2b,\qquad
                   s+\ell/2-\ell b,\qquad
                   b+e/2-e s .
\]
The minimum is at most the weighted average with weights
\[
           \frac{e\ell-1}{e(\ell+2)+1},\qquad
           \frac{2e}{e(\ell+2)+1},\qquad
           \frac{2}{e(\ell+2)+1}.
\]
This cancels \(s,b\) and gives \(\Theta_\ell\). Also \(\Theta_\ell>3/2\)
since \(e(\ell-4)>5\). The logarithmic number of boxes and the dyadic
field windows are absorbed in the allowance. The argument uses
the prime-pair tail at the composite exponent \(e\); that tail has the
same proof as for the prime exponents.

For \(n=75\), \(\Theta_5=82/53\), and the field-count exponents are
\(82/1961\) in degree \(75\) and \(82/3975\) in regular degree \(150\).
These counts and the table allow arbitrary relative ramification.

At \(\ell=5\), the mixed bound reaches the limit of just the three tails
used in this proof. In an abstract window put \(\asymp X^{91/106}\)
objects at \(h_3\asymp X^{15/106}\), \(h_5\asymp X^{29/106}\),
and \(h_{15}=h_3 h_5\asymp X^{44/106}\), with all other objects of
constant weight. All three displayed box bounds allow these objects at
their upper thresholds and hence at smaller thresholds. The
rank-one and rank-two terms in (1) for \(5,15\), and the
[signed quadratic \(3\)-torsion tail](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md)
at \(15/106\), also allow them. Their contribution is
\(\asymp X^{82/53}\). A further saving requires a restriction on
the actual simultaneous fields and witnesses near these thresholds.

## Moments on the combined \(15\)-torsion group

**Corollary (balanced \(3\)- and \(5\)-torsion moments).** For either sign
and every fixed real \(1\le k\le6\),
\[
            \sum_{K:\,|D_K|\le X} (h_3(K)h_5(K))^k
               \ll_{k,\varepsilon} X^{(91+44k)/106+\varepsilon}.
                                                                    \tag{6}
\]
In particular for any fixed abelian kernel of order \(15^j\) with
integer \(1\le j\le6\), the two generalized dihedral count exponents
are \((91+44j)/(53(15^j-1))\) and \((91+44j)/(106\cdot15^j)\).

Keep the same boxes in \(h_5,h_{15}\), now multiplying a tail by
\(X^{ks}\). In the region \(b\le3/10\) and \(s\le13/30\) the exponents
from the first moment of \(3\)-torsion and the two pair bounds are
\[
            1+(k-1)s+b,\qquad
            5/2+ks-6b,\qquad 15/2+(k-16)s .
\]
Their minimum is bounded by the weighted average with weights
\(6(16-k)/106,\ (16-k)/106,\ (7k-6)/106\), respectively. The weights
are nonnegative and sum to one for \(1\le k\le6\). Both \(s,b\) cancel
and the constant is \((91+44k)/106\).

Outside \(O(X^\varepsilon)\) exceptional fields the small-split-prime
cap for exponent \(15\) gives \(s\le7/15+O(\varepsilon)\).
For \(s\ge13/30\), its pair tail \(X/H\) gives exponent
\(1+7(k-1)/15\). For \(b\ge3/10\) and \(s\le13/30\), the same tail for
\(5\)-torsion gives exponent \(7/10+13k/30\). The exceptional fields
contribute at most \(X^{k/2+O(\varepsilon)}\), by the class-number bound
on \(h_{15}\). All three exponents are at most \((91+44k)/106\) in
the stated range. The abstract window from (5) attains (6) as well.
Hölder is unnecessary here: the argument keeps the product threshold
of the two torsion factors on every field.
