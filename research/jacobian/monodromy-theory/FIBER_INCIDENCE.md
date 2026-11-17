# Fiber-incidence constraints on an affine completion

Work over an algebraically closed characteristic-zero field k.
This note gives further necessary conditions for an ordinary-affine-space
Keller realization of an **independently stabilized polynomial curve cover**.
It does not assume that the curve coordinate is a polynomial output.

**KMFI23 (connected fiber-incidence theorem).** Suppose a polynomial Keller
map F:A^m->A^m realizes, as a field pair,

\[
 K=k(f(r),t_1,\ldots,t_{m-1})\subset
 L=k(r,t_1,\ldots,t_{m-1}),\qquad \deg f=n>1,
\]

with independent t_i and f a polynomial. Write
f(r)-a=c product_i(r-c_i)^{e_i} for a finite value a.
Then the source fibers r=c_i and the prime components of the
target pencil fiber f(r)=a form a connected incidence graph.
Edges are source prime divisors, and an edge on row i and
target component q of multiplicity N_q has source-pencil
multiplicity m satisfying

\[
 e_i m=N_q. \tag{1}
\]

If this branch fiber has at least two unramified points, then
the gcd of its target component multiplicities is one.
In particular, if its passport is e^s 1^t with e>=2, s>=1, t>=2,
some **unramified** source-pencil fiber r=c_i is reducible,
with one prime component of multiplicity divisible by e and another
whose multiplicity is not divisible by e.

For the retained degree-eleven PSL2(11) Shabat pair, this forces:

* over branch value 0, one of the two unramified tail fibers
  has both a component of multiplicity divisible by three
  and a component of multiplicity not divisible by three;
* over branch value 1, one of the three unramified tail fibers
  has both an even-multiplicity and an odd-multiplicity component;
* the two target branch fibers are both reducible and nonreduced;
* their source total transforms have at least six and eight distinct
  prime components, respectively.

Thus any such affine completion requires two distinct reducible,
nonreduced source-pencil fibers with these mixed multiplicities.
Reduced source pencils, source pencils with irreducible fibers at
all the relevant tails, and target pencils with an irreducible branch
fiber are all excluded.

The theorem concerns the independent-stabilization field class.
It does not exclude the abstract PSL2(11) actions or a different
higher-dimensional cover. The full objective remains open.

## 1. Polynomial representatives after denominator descent

Write r=U/V with relatively prime source polynomials U,V.
Write t=f(r)=a_0/b_0 using relatively prime target polynomials.
Affine Keller pullback preserves coprimeness and pulls every prime
divisor back to a nonempty squarefree divisor (KMFI1).

The reduced source denominator of f(U/V) is V^n, up to a constant.
Consequently b_0(F)=c V^n. The divisor-power injection of KMFI1
makes b_0 a constant times an n-th power b^n.
Absorbing constants, we can arrange

\[
 V=b(F),\quad t=a/b^n,\quad
 a(F)=f_h(U,V),\quad \gcd(a,b)=\gcd(U,V)=1. \tag{2}
\]

Here f_h is the degree-n homogeneous numerator of f.
The scaling of U,V is adjusted together, so r remains unchanged.
No polynomiality of r itself is asserted.
This is the denominator step already used in KMFI9.

For a finite value a_* put A=a-a_*b^n.
Then

\[
 A(F)=c\prod_i H_i^{e_i},\qquad H_i=U-c_iV. \tag{3}
\]

The H_i are pairwise coprime. Each is nonconstant:
if H_i were a nonzero constant, (2) would give
U=c_iV+H_i in K, hence r in K, contradicting n>1.
The identically zero case contradicts nonconstant r.

We also have two field facts:

\[
 k(r)\cap K=k(f(r)),\qquad
 k(r)\text{ is relatively algebraically closed in }L. \tag{4}
\]

Both follow from the displayed independent transcendence bases.
For the first, K is a purely transcendental, hence regular,
extension of k(f(r)); intersect with the finite extension k(r).
For the second use L=k(r)(t_1,...,t_{m-1}).

## 2. The incidence graph and its connectedness

Factor A=c product_q q^{N_q} into target primes.
Make one left vertex for each c_i, one right vertex for each q.
A prime factor h of H_i is an edge joining i to the unique q
whose pullback q(F) contains h. These are all the source primes
in (3), because target pullbacks have disjoint supports.
Every edge has multiplicity one in q(F).
Comparing its valuation in (3) proves (1).

There are no isolated vertices. Suppose the graph had a component
with left set I and right set J, with I a nonempty proper subset.
The weighted product of its target primes would satisfy

\[
 \prod_{q\in J}q(F)^{N_q}
   =c_I\prod_{i\in I}H_i^{e_i}.
\]

Divide by V^{E_I}, E_I=sum_{i in I}e_i, which belongs to K by (2).
It follows that product_{i in I}(r-c_i)^{e_i} lies in
K intersect k(r)=k(f(r)).

This is impossible. A nonconstant polynomial in r belonging to
k(f(r)) is a polynomial in f(r): a finite pole of a rational
function in f would pull back to a finite pole in r.
Its zero divisor must therefore contain a complete fiber of f
whenever it contains any point of that fiber.
The displayed polynomial contains only a proper part of f^{-1}(a_*).
This proves connectedness.

If d is the number of distinct c_i and z the number of target
prime components, connectedness gives at least d+z-1 edges.
Equivalently,

\[
 \sum_i(\#\{\text{prime factors of }H_i\}-1)\ge z-1. \tag{5}
\]

This is a necessary divisorial incidence condition, not a
sufficiency assertion.

## 3. Two tails force a primitive target fiber

Suppose two distinct unramified points c_i,c_j occur, so e_i=e_j=1.
Let g=gcd_q N_q. Equation (1) shows that all prime multiplicities
in both H_i and H_j are divisible by g.
If g>1, unique factorization and algebraically closed constants give
H_i=h_i^g and H_j=h_j^g after scaling. Therefore

\[
 (r-c_i)/(r-c_j)=c(h_i/h_j)^g.
\]

The root on the right is algebraic over k(r) and lies in L,
so (4) puts it in k(r). But a Mobius function on P^1 has a simple
zero and a simple pole and cannot be a nontrivial g-th power.
Thus g=1.

Now assume the passport e^s1^t, e>=2, s>=1, t>=2.
Every target vertex adjacent to a ramified left vertex has
e|N_q, by (1). Since gcd_q N_q=1, some target vertex has
e not dividing N_q. In a connected graph, take a path from
such a vertex to a ramified left vertex. At the first transition
from a target multiplicity not divisible by e to one divisible
by e, the intervening left vertex must be unramified.
Its fiber has two different prime components whose multiplicities
are N_q, one divisible by e and the other not.

This proves the mixed-multiplicity assertion on a single tail
fiber. It also proves z>=2, and some N_q>=e, so the target
branch fiber is both reducible and nonreduced.

## 4. The exceptional PSL2(11) completion gate

The retained Shabat polynomials have passports

\[
 3^3 1^2\mid 2^4 1^3\mid 11.
\]

Their explicit factors and exact identities are in the
[canonical benchmark](../../foundations/extended-geometry/PSL2_11_KELLER_ACTION_SPECTRUM.md).
KMFI9--10 force this polynomial passport for surviving
independently stabilized exceptional degree-eleven curve inputs.

Apply Section 3 twice. There are d=5 left vertices over 0
and d=7 over 1. Both target fibers have z>=2.
Equation (5) gives at least six and eight distinct source primes.
The two sets of primes are disjoint, since the H_i for distinct
r-values are pairwise coprime. There are therefore at least
fourteen distinct source prime components in these two fibers.

The reducible mixed tail above 0 is one of the two roots of p_3;
the one above 1 is one of the three roots of q_3.
These are different fibers. A source pencil with only reduced
fibers, or with every p_3/q_3 tail fiber irreducible, cannot
give an affine completion, even if all local inertia tests pass.
A failed criterion excludes the completion architecture.
Passing it does not construct a map.

## 5. An actual S4 control with the required mixed tail

The conditions allow known realizations. Use KMFI14 with p=2,d=3.
Let s^2=-2, zeta=-1+s, delta=4+s, and put

\[
 g(L)=1+\frac{4+s}{3}L+\frac{-10+2s}{3}L^2,\quad
 \rho=g(L)+L^2v,\quad Z=\rho^3-\delta L,\quad Q=(Z-1)/3,
\]
\[
 P=\frac6{\delta^2}
       \left((\rho^4-1)/4-Z(\rho-1)\right).
\]

Then A=Q/L^2 and B=P/L^2 are polynomials vanishing at (1,0).
For L=1+x^2y,v=x^2z, the map

\[
 (A/x^2,B/x^2,xL)
\]

is polynomial, has Jacobian -2, generic degree four, and group S4.
The fixed-curve coordinate is r=w=rho^3/Z, with

\[
 \tau=w(w/4-1)^3,\qquad
 64(\tau+27/64)=(w-1)^2(w^2-10w+27).
\]

The latter fiber has passport 2 1^2.
Its tails are c=5+s and 5-s.
For c=5+s, put H=rho^3-cZ and substitute L,v as above.
The checker verifies

\[
 H=x^2 H_{\rm res},\qquad x\nmid H_{\rm res},
 \qquad H_{\rm res}|_{x=0}=\alpha y+\beta z,\quad \beta\ne0.
\]

Moreover H is squarefree in k[L,v].
Indeed H_v=3(1-c)rho^2L^2, while H is coprime to both rho and L.
Away from x=0, the change from y,z to L,v is invertible.
Thus H_res is squarefree globally and nonconstant.
The tail fiber has one component of multiplicity two and residual
components of multiplicity one, exactly as the theorem requires.

This control also exhibits the full incidence tree.
In target coordinates T1,T2,C set

\[
 Z_t=1+3C^2T_1,\quad
 R_t=\delta^2 C^2T_2/6-3C^2T_1-3/4.
\]

Then 64R_t^3+27Z_t^4=C^2 D, with
D|_{C=0}=18delta^2 T2.
The target primes are C of multiplicity two and D of multiplicity one.
Indeed away from C=0 the variable change to R_t,Z_t is invertible
and the binomial 64R_t^3+27Z_t^4 is irreducible (3 and 4 are coprime);
D is not divisible by C.

Put H_minus=rho^3-(5-s)Z after source substitution.
The exact pullback is D(F)=delta^2 H_res H_minus.
Both residual factors are irreducible.
To see this first in k[L,v,L^{-1}], use the invertible change to
rho,L: each H becomes a nonzero constant times rho^3+kappa L.
Its curve there is a copy of G_m, and it has no L factor.
This proves geometric irreducibility in k[L,v], and after the
invertible x-localized change to y,z proves irreducibility of
H_res and H_minus. Neither has an x factor.

Thus the graph has four prime edges:

| source-pencil point | source prime | source multiplicity | target prime | target multiplicity |
|---|---|---:|---|---:|
| w=1, index two | L | 1 | C | 2 |
| w=5+s, index one | x | 2 | C | 2 |
| w=5+s, index one | H_res | 1 | D | 1 |
| w=5-s, index one | H_minus | 1 | D | 1 |

It is connected and attains d+z-1=3+2-1=4.
The criterion records an actual feature of an affine Keller map.

The [checker](verify_fiber_incidence.py), [model](fiber_incidence_model.json)
and [receipt](fiber_incidence_receipt.json) retain these exact controls.
The universal graph and field proofs are written above; finite controls
are not formal verification or independent whole-proof review.
