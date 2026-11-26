# Three distinct prime tests beyond the absolute completion bound

The primitive false incidences in the balanced small-budget box persist beyond
total modulus exponent \(3/2\) when the test modulus has three appropriately
sized prime factors. On the dual side a Fourier coefficient of the Fermat cone
is a centered count on a moving line. Shifted projections with disjoint branch
sets give the four-shift correlation that is needed to apply the arithmetic
van der Corput method twice. The resulting signed lower count concerns the
modular noise in a box with no genuine integer equality.

## The three-prime false-incidence count

Fix an integer \(d\ge3\) and \(3/2\le\beta<38/25\). Define
\[
 h=\beta-1,\qquad \sigma=\frac{3\beta-4}{14},\qquad
 a_1=h-2\sigma,\qquad a_2=h-4\sigma,\qquad a_3=\beta-a_1-a_2.
\]
Independently choose primes \(p_i\in[N^{a_i},2N^{a_i}]\), for
\(1\le i\le3\), and put \(Q=p_1p_2p_3\). The three intervals are eventually
disjoint. Use the integer dictionaries
\[
 A=x^d,\quad B=y^d,\quad C=z^d,\qquad
 x,y\in[N,(1+s)N],\quad z\in[(1+2s)N,(1+3s)N],
 \qquad s=(100d)^{-1}.
\]
For all but a proportion tending to zero of the prime triples there are
\[
                    \gg_{d,\beta} N^{3-\beta}                 \tag{FI3}
\]
distinct three-radical records of pairwise coprime \(A,B,C\) in these
dictionaries satisfying \(A+B\equiv C\bmod Q\).
The three representatives are the same at every prime. The roots in
these records are squarefree, so their radicals are exactly \(x,y,z\).
They obey \(\operatorname{rad}(ABC)=xyz<z^3=C^{3/d}\), and the
one-target minimal repair at fixed order \(M\ge d\) is \(z^M\).
The choice of \(s\) gives \(A+B>C>\max(A,B)\) throughout the
dictionaries. In particular the result holds at arbitrarily small
fixed positive budgets \(3/d\). All parameters are held fixed before
\(N\) tends to infinity.
[Previous absolute-completion range](ABC_MODULAR_FALSE_INCIDENCES.md)

## Local four-shift correlation for the dual cone

Write \(e_m(t)=\exp(2\pi i t/m)\), and let \(V_p\) be the set of
unit-coordinate triples on the Fermat cone \(X^d+Y^d=Z^d\) modulo \(p\).
Its projectivization \(C_p^\circ\) is the smooth projective Fermat curve
with its coordinate-zero points removed. For \(\xi\ne0\) define
\[
 f_p(\xi)=p^{-1}\sum_{w\in V_p}e_p(\xi\cdot w)
  =\#(C_p^\circ\cap\{\xi\cdot w=0\})(\mathbb F_p)-c_p,\qquad
 c_p=|C_p^\circ(\mathbb F_p)|/p=1+O_d(p^{-1/2}).
\]
The identity follows by summing on each line of the cone, and the
bound on \(c_p\) follows from the Weil bound for smooth projective
curves. Thus \(f_p\) is real, invariant under scalar multiplication
by units, and bounded in magnitude by a constant depending only
on \(d\) away from \(\xi=0\). At zero use the actual value
\(p^{-1}|V_p|=O_d(p)\) unless explicitly replacing it by zero.

Fix \(u,v\in\mathbb F_p\) away from the exceptional directions
\[
          u\zeta+v=0 \quad(\zeta^d=1\text{ over }\overline{\mathbb F}_p).
\]
This excludes \((u,v)=(0,0)\). Set \(f(t)=f_p(t,u,v)\).
On \(X=1\), \(f(t)\) is the number of points of
\(1+y^d=z^d\), with \(yz\ne0\), over the level \(t\) of the map
\(\pi(y,z)=-uy-vz\), minus \(c_p\). The excluded directions
are precisely those for which the projective pencil has a base
point on the curve. Otherwise \(\pi\) has degree \(d\),
is separable for \(p>d!\), and is unramified over infinity because
\(X=0\) meets the Fermat curve transversely. Let \(B\) be the finite
geometric branch set of \(\pi\), of size \(O_d(1)\).
If \(t+s_1,\ldots,t+s_j\) are shifts with pairwise distinct
entries and \(s_i-s_{i'}\notin B-B\) for \(i\ne i'\), for
\(1\le j\le4\), then
\[
 \sum_{t\in\mathbb F_p}e_p(bt)\prod_{i=1}^{j} f(t+s_i)
                         =O_d(\sqrt p)\qquad(b\in\mathbb F_p).
                                                               \tag{FI4}
\]
To see this, first replace \(c_p\) by \(1\), at total cost
\(O_d(\sqrt p)\) because the local functions are bounded.
The shifted covers have pairwise disjoint finite branch sets.
Their geometric Galois closures are tame of degree at most \(d!\)
and are linearly disjoint. Indeed, a given closure and the
compositum of the others can intersect only in a tame cover of
\(\mathbb P^1\) unramified at every finite point. Such a cover is
trivial by Riemann--Hurwitz, even if one allows ramification at
infinity. Induction on the number of closures gives the assertion.
It follows that the normalized fiber product of each nonempty
subset of the shifted covers is geometrically integral. Its
genus is bounded in terms of \(d\), by Riemann--Hurwitz using
the bounded degrees and branch sets.

Expand the product of the centered counts. The term for each
nonempty subset is the point count of the corresponding fiber
product, with finitely many fibers removed or corrected. These
exceptions include the deleted coordinate-zero points and the
branch sets and cost \(O_d(1)\) per term. For \(b=0\), the Weil
bound gives \(p+O_d(\sqrt p)\) for each subset, also allowing
the empty term \(p\). The binomial sum of the main terms is zero.
For \(b\ne0\), apply the Weil bound for additive sums on a smooth
curve to the same normalization and the function \(bt\).
Its poles are simple because the map to the common line is
unramified at infinity, and thus it is not of
Artin--Schreier form even over the algebraic closure. The number
of poles and the genus are bounded in terms of \(d\).
The empty term in this case vanishes. This gives (FI4).
Completion by additive Fourier inversion bounds the same correlation
on any interval of length at most \(p\) by \(O_d(\sqrt p\log(2p))\).
In particular this applies to the shifts \(0,a,b',a+b'\) when
the increments \(a,b',a+b',a-b'\) all avoid the symmetric finite
set \(\{0\}\cup(B-B)\).
[Bombieri's curve estimate, Theorem 5 in Section VI](https://doi.org/10.2307/2373048)

## Two arithmetic differencing steps

Write \(X_0=Q/N\asymp N^h\). Fix unit progression steps \(d_i\)
modulo \(Q\), with \(1\le d_i\le N^{2\rho}\), where \(\rho>0\)
will be sufficiently small. In the Fourier sum change variables from
\(\xi_i\) to \(t_i=d_i\xi_i\bmod Q\), using signed least representatives.
For fixed \(t_2,t_3\) in a generic direction modulo \(p_3\), view the
product of local normalized cone Fourier coefficients as a function of
\(t=t_1\). Put its local factors at zero equal to zero in this
paragraph. Each factor is bounded by \(O_d(1)\) and has period its
own prime. The \(p_3\) factor is of the form in (FI4), with \(t\)
multiplied by a unit modulo \(p_3\).
Suppose \(J\) is an interval of length at most \(H\), where
\(X_0\le H\ll X_0N^{2\rho+\epsilon}\), with a sufficiently small
fixed \(\epsilon>0\). Uniformly for any linear exponential twist
\(e(\theta t)\) on the product, the sum on \(J\) satisfies
\[
  \ll_d H\left\{
       (p_1/H)^{1/2}+(p_2/H)^{1/4}
                  +\left(\sqrt{p_3}\log(2p_3)/H\right)^{1/4}
            \right\}.                                       \tag{FI5}
\]
For the first step average shifts \(a=p_1 j\) with
\(1\le j\le \lfloor H/p_1\rfloor\) and use Cauchy--Schwarz
to remove the factor of period \(p_1\). The diagonal cost in
the squared normalized bound is \(O_d(p_1/H)\).
The other terms average normalized correlation sums in the
two remaining factors. The original linear twist disappears
from each correlation apart from a scalar of modulus one.
Among the nonzero increments \(a\) only \(O_d(1)\) are
bad for \(p_3\), since the increment parameters lie in an
interval shorter than \(p_3\) and \(p_1\) is a unit.
Use the trivial \(O_d(H)\) bound for those correlations.

On a good first increment apply the same step with shifts
\(b'=p_2 j'\), now removing the first-differenced factor
of period \(p_2\). The second diagonal cost is \(O_d(p_2/H)\).
For fixed good \(a\), only \(O_d(1)\) second increments
have \(b'\) or \(a\pm b'\) in the bad set, because the
second increment parameters are also shorter than \(p_3\).
All other sums in the twice-differenced \(p_3\) factor
are bounded by \(O_d(\sqrt{p_3}\log(2p_3))\) from (FI4) and
completion. Intersections of shifted copies of \(J\) are
intervals, so no loss arises from boundaries in this bound.
Taking the two square roots proves (FI5). For completeness,
a differencing step follows by writing a sum \(S\) multiplied
by \(T\) as the sum of its \(T\) translates with the
interval indicator, pulling the periodic factor outside
the translate average, and applying Cauchy--Schwarz over
a union of intervals of size \(O(H)\). Diagonal contribution
\(O(H^2/T)\) and grouping pairs by their increment give
the bounds just used.
[Irving, Lemmas 2.1--2.2 for the arithmetic differencing steps](https://arxiv.org/pdf/1503.07156#page=5)

The parameters give
\[
   (h-a_1)/2=(h-a_2)/4=(h-a_3/2)/4=\sigma,\qquad
                       \sigma>2\beta-3 .
\]
Indeed the last inequality is equivalent to \(25\beta<38\).
On the stated range \(a_1,a_2>0\), both are below \(h\), and
\[
                  a_3>h,\qquad a_1+a_2>h.
\]
Choose \(\rho,\epsilon\) small against these strict margins.
Then \(H\ll p_3\), the increment conditions used above hold,
and (FI5) is \(O(HN^{-\sigma+o(1)})\) uniformly in \(J\).

## Weighted progression counts

Take fixed smooth nonnegative functions \(\phi_i\) in the interiors of
the three root intervals scaled by \(N\), each bounded by one and with
positive integral. For a residue class \(x_i=r_i\bmod d_i\) with
\(0\le r_i<d_i\), denote the mass of \(\phi_i(x_i/N)\) in the class
by \(W_i(0)\). We show that the weighted count of unit-residue cone solutions modulo \(Q\), before requiring squarefree
or pairwise coprime roots, is
\[
 \frac{1}{Q}
   \left(\prod_{p\mid Q}\frac{|V_p|}{p^2}\right)\prod_i W_i(0)
           +O_{d,\boldsymbol\phi,\beta,\epsilon}
                     \left(QN^{-\sigma+2\epsilon}\right).      \tag{FI6}
\]
Here unit residues refers to the three coordinates separately at each
of the test primes. The formula is uniform in all the progression steps
and residues above.

For the changed frequency \(t_i\), the Fourier transform of the
progression weight is a linear phase in \(t_i\) times
\[
    \sum_\ell \phi_i((r_i+d_i\ell)/N)\,e_Q(-t_i\ell).
\]
Write \(X_i=Qd_i/N\). The latter has magnitude
\[
 \ll_{\phi_i,A} (N/d_i)(1+|t_i|/X_i)^{-A}
                  \qquad (|t_i|\le Q/2),
\]
and its derivative with respect to real \(t_i\) has the same bound
divided by \(X_i\). These follow from discrete summation by parts
for the function \(\phi_i((r_i+d_i\ell)/N)\), and for the derivative
use the same calculation with the additional factor \(\ell/Q\).
The derivative bound applies to the part of the transform with the
linear phase removed. Hence each Fourier transform has total
absolute mass \(O(Q)\), and the phase-free part has total variation
\(O(N/d_i)\). Truncate all three frequencies to
\(|t_i|\le X_i N^\epsilon\) at a negligible cost, taking the decay
order large in terms of \(\epsilon\). This bound is uniform because
\(d_i\le N^{2\rho}\) and the three smooth profiles are fixed.

By the Chinese remainder theorem the normalized dual coefficient
away from the origin is \(Q^{-2}\) times the product of the three
\(f_p\) factors on the relevant frequencies, including their actual
values at zeros. The common local scalar factors in the CRT do not
matter by projective scalar invariance. If a frequency vector is
zero at one of \(p_1,p_2\), the sum of the absolute Fourier weights
over multiples of that prime in each coordinate is \(O(Q/p_i)\),
since \(X_i\gg p_1,p_2\). With the \(O_d(p_i)\) local zero value
these terms cost \(O_d(Q/p_i^2)\). No nonzero frequency vector
of the truncated region can be zero at \(p_3\) or at both
\(p_1,p_2\), by the margins \(a_3>h,\ a_1+a_2>h\) and the
choice of \(\rho,\epsilon\). The zero frequency contributes the
main term in (FI6).

In the remaining sum replace the local zero values by zero.
For any fixed \(t_2\), there are \(O_d(1)\) choices of \(t_3\)
in the truncated interval that give an exceptional direction at
\(p_3\): there are \(O_d(1)\) offending residues and that interval
is shorter than \(p_3\). Their contribution after Fourier
normalization is \(O_d(Q/X_3)\), using total absolute mass
\(O(Q)\) for each of the first two weights and maximum
\(O(N/d_3)\) for the third. For the other directions apply
(FI5) to the \(t_1\) sum, including its linear phase.
Partial summation and the variation bound give an upper
bound \(O_d(QN^{-\sigma+\epsilon+o(1)})\) for this inner
sum, before the other two Fourier weights and the normalization.
Indeed one may use \(H=O(X_1N^\epsilon)\) and the saving
factor in (FI5) is \(O(N^{-\sigma+o(1)})\) for every
\(H\ge X_0\). Multiply by the other two absolute masses
\(O(Q^2)\) and by \(Q^{-2}\). Together with the zero and
exceptional direction costs, this proves (FI6).

## Removing square factors and common divisors

First choose a large fixed \(P_0\). For each prime \(p\le P_0\),
require that all roots are nonzero modulo \(p^2\) and at most one
of them is zero modulo \(p\). Each condition has local density
\(1-O(p^{-2})>0\), so the product of the local densities is
bounded below uniformly in \(P_0\). Apply (FI6) separately to
the allowed CRT progressions. Their weighted main term is
\(\gg_{d,\boldsymbol\phi}N^3/Q\), and their errors are
\(o(N^3/Q)\) because \(P_0\) is fixed and
\(\sigma>2\beta-3\).

Take \(U=N^\rho\). For any prime \(P_0<p\le U\), apply (FI6)
to each of the three single square conditions and three common
prime conditions, bounding from above without the small-prime
filters. The relevant progression steps are at most \(U^2\)
and are units mod \(Q\). Each main term is \(O(N^3/(Qp^2))\)
and the total errors over these primes are
\(O(UQN^{-\sigma+2\epsilon})=o(N^3/Q)\) on choosing
\(\rho,2\epsilon\) smaller than the strict saving margin.
Hence these exclusions remove at most
\(O(N^3/(QP_0))+o(N^3/Q)\).

The unconditional number of root triples violating one of the
squarefree or pairwise coprime conditions at some \(p>U\)
is \(O(N^3/U+N^{5/2})\). For each such fixed root triple, the
nonzero integer \(x^d+y^d-z^d\) has size \(O_d(N^d)\).
It has only \(O_{d,\beta}(1)\) prime factors at each test scale,
and therefore passes only \(O_{d,\beta}(1)\) ordered test triples.
The number of prime triples is \(\asymp_\beta N^\beta/(\log N)^3\)
by the prime number theorem. Thus the average loss to this tail
is \(o(N^{3-\beta})\). Choose \(P_0\) sufficiently large and use
Markov's inequality to retain all but a proportion tending to zero of the tests.
The nonnegative weights are at most one and the surviving roots
are squarefree, so the weighted lower bound implies the number
of distinct radical records in (FI3).
[Squarefree and pairwise coprime sieve with the same prime averaging](ABC_MODULAR_FALSE_INCIDENCES.md#preserving-squarefreeness-and-primitivity)

This gives signed modular incidence beyond the boundary of the
absolute Fourier estimate, for a specified factorization of the
test modulus. The total exponent stays below the scale greater
than two at which average tests recover true occurrence with a
power-small error. The construction supplies an obstruction in
that smaller-modulus range, with no genuine exceptional triples.
