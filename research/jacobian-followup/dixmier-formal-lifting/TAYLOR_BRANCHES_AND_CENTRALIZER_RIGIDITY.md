# Taylor branches and centralizer rigidity without finite generation

15 November 2026. **QTB1.** Formal Taylor maps remove the finite-generation
hypothesis from QET1's birational and Galois conclusions. An actual scalar
centralizer admitting individual conjugates has no nontrivial finite
Galois factor over its position field. If its fraction field is birational,
the centralizer is polynomial and finite full completion follows.

There is a finite obstruction certificate: an actual ratio of commuting
operators cannot supply a new radical of the position field. This tests
all possible individual conjugates without constructing the whole
centralizer or proving it finitely generated. No such certificate or
birational field proof is supplied for a new Long pair.

## The differential algebra hypotheses

Let L have characteristic zero, let P=L[t1,...,tm], and let B be a
commutative domain containing P, with commuting derivations nabla_i
extending all coordinate partials. Put Q=Frac P and H=Frac B.
Assume H/Q is algebraic, not necessarily finite, and
\[
 \bigcap_i\ker(\nabla_i:H\to H)=L,\qquad
                              P\cap B^\times=L^\times.   \tag{1}
\]
No noetherian or finite-generation hypothesis is imposed on B.

These hypotheses hold for the actual generic scalar centralizer of m
commuting images in the rank-m Weyl algebra once QFC1's polynomial
classical-centralizer hypothesis and actual individual mates hold.
Commutativity is QSC1/QFC1. Algebraicity follows by applying
[Bavula's commutative-subalgebra bound, Corollary1.6](https://arxiv.org/html/math/0401263)
to each finitely generated P[b], not to all of B at once. Scalar extension
to an algebraic closure preserves that bound. Weyl units are scalar by
additivity of Bernstein degree.

For Long, let L=K((h)). QFC1 gives a unique Laurent formal expansion
sum h^j*P_j(F) for each actual centralizer element. Hence H embeds in the
commutative h-valued position-field completion K(t1,...,tm)((h)).
The nabla_i act there coefficientwise as partial_(ti), whose joint
constants are K((h))=L. This proves the first part of (1) for H itself,
without first proving B normal or finitely generated. QET1 removes
actual parameter poles in the mates, so generic individual mates suffice.

For an arbitrary projective endomorphism algebra, algebraicity and the
unit condition need their own proof. No scalar unit theorem is being
transferred automatically.

## Formal Taylor maps over arbitrary residue fields

As proved in QET1, B is differentially simple. If a nonzero b lies in
a differential ideal I, clearing denominators in an algebraic equation
for b with nonzero constant term puts a nonzero element of P in I.
Coordinate partials then extract a nonzero scalar, so I=B.

Choose any maximal ideal n of B, put kappa=B/n, and write ci for the
image of ti. Its residue field may have arbitrary transcendence degree
over L. Commutativity and the Leibniz rule give the differential ring map
\[
 \tau_n(b)=\sum_{\alpha\ge0}
       {\overline{\nabla^\alpha b}\over\alpha!}s^\alpha
       \ \in R=\kappa[[s_1,\ldots,s_m]],\qquad
                         \tau_n(t_i)=c_i+s_i.             \tag{2}
\]
Every coefficient uses finitely many operations. The kernel is a
differential ideal, and tau(1)=1, so differential simplicity makes tau
injective. It extends to an embedding of H in Frac R.

The restricted coordinate map
\[
                         P\longrightarrow R,\qquad
                                      t_i\mapsto c_i+s_i  \tag{3}
\]
is flat. Factor it through the field extension kappa[t], translation
by c, localization at (t-c), and completion. The final flatness input is
[Stacks00MB](https://stacks.math.columbia.edu/tag/00MB). This argument
uses no finite presentation of B and no assumption that n contracts
to a closed point of Spec P.

## Birational differential overrings are localizations

Here is a useful independent lemma. If a differential subring A of Q
contains P and is stable under all partial_i, then
\[
                    A=S^{-1}P,\qquad S=P\cap A^\times.    \tag{4}
\]
Finite generation of A is unnecessary.

Write b=p/q in A with coprime p,q in P. If q were not a unit of A,
choose a maximal ideal containing q and apply the Taylor ring map (2).
Its restriction (3) is flat even without differential simplicity.
Coprimeness makes multiplication by p injective on P/(q), and flatness
preserves this injection on R/(q). But p=q*tau(b) in R, so multiplication
by p there is zero. Since q maps into the maximal ideal of R,
R/(q) is nonzero, a contradiction. Thus q is a unit of A. This proves
both containments in (4).

In particular H=Q and (1) imply B=P, with no appeal to an etale
finite-type model.

## Every base divisor has an unramified formal branch

Let pi be any irreducible polynomial in P. By (1), it is not a unit
of B. Choose a maximal ideal n containing pi and use (2). The nonzero
series pi(c+s) lies in the maximal ideal of R. The ring R is regular
local and factorial; see
[Stacks0AG0](https://stacks.math.columbia.edu/tag/0AG0).
Choose an irreducible factor rho of pi(c+s).

For any a in P not divisible by pi, the sequence pi,a is regular in
P, and flatness of (3) preserves injectivity of multiplication by a
on R/(pi). In a factorial ring this says that rho cannot divide a.
Moreover, choose i with partial_i pi nonzero and not divisible by pi,
possible in characteristic zero. If rho occurred with multiplicity
at least two in pi(c+s), it would divide every s derivative of that
series, contradicting the same regular-sequence argument for
a=partial_i pi. Consequently
\[
         v_\rho(\pi)=1,\qquad v_\rho(a)=0
                                  \quad(a\in P\setminus(\pi)). \tag{5}
\]
The discrete valuation of Frac R at rho therefore restricts to exactly
the pi-adic valuation on Q.

For any finite intermediate field Q subset J subset H, restrict the
same valuation to J. Its value group is contained in Z and contains
v(pi)=1, hence equals Z. Its ramification index over Q is one:
\[
 \text{for every base prime divisor and every finite }J/Q,
              \text{ some extension valuation has }e=1. \tag{6}
\]
The residue extensions are separable in characteristic zero, so the
corresponding prime in the finite normalization is unramified.
This is the divisor-retention statement needed by QET1, proved directly
without constructing Spec B as a finite-type etale scheme.

The point c may lie in the singular locus of pi=0. The argument uses
a prime factor in the formal completion, not total vanishing order
at that point. Equation (5) still holds for each such factor.

## No finite Galois factor

Suppose J/Q in (6) is finite Galois. Galois transitivity gives the
same ramification index at all primes over each base divisor.
Thus every such prime is unramified. Normalize P in J; that
normalization is finite and normal. Branch purity over the regular
base makes it finite etale, hence a constant cover of affine space.
The explicit external inputs are
[Stacks09E3](https://stacks.math.columbia.edu/tag/09E3),
[Stacks0BMB](https://stacks.math.columbia.edu/tag/0BMB), and
[Milne's characteristic-zero covering theorem](https://www.jmilne.org/math/CourseNotes/LECc.pdf#page=27).

Every resulting constant element is algebraic over L, so all nabla_i
kill it. Equation (1), now applied in H, excludes every new constant.
Therefore
\[
                  Q\subset J\subset H,\quad J/Q
                    \text{ finite Galois}\quad\Longrightarrow\quad J=Q.
                                                               \tag{7}
\]
Only the normalization of the specified finite field J is used.
There is no normalization-finiteness assumption on all of B or H.

If H/Q is Galois, even an infinite algebraic Galois extension, it is
the union of its finite Galois subextensions, so H=Q. Equation (4)
and (1) then give B=P. In the quantum setting,
\[
  \operatorname{Cent}_{E}(F)=K[[h]][F],\qquad
                    \text{finite full completion follows}.   \tag{8}
\]
Use QET1 for actual parameter-pole removal and QFC1 for integral
polynomial curvature homotopy. A prescribed classical flat tuple can
be retained. For a fixed pair, a nontrivial finite Galois subfield
of its actual centralizer field instead forbids simultaneous
individual mates in all orders and degrees.

## A stronger monodromy condition

For a finite intermediate J/Q let M be its normal closure,
G=Gal(M/Q), and S=Gal(M/J), the point stabilizer in the transitive
action. The largest Galois subfield of J is fixed by the normal
closure of S in G. Equation (7) therefore requires
\[
                           \langle gSg^{-1}:g\in G\rangle=G. \tag{9}
\]
Every codimension-one inertia generator must also have a fixed
sheet, by (6). Neither condition supplies a common sheet fixed by G.

The natural degree-four D4 action fails (9): a vertex stabilizer
is a reflection, and its normal closure has order four, index two.
Thus a D4 quartic centralizer subfield is excluded, although that
quartic field itself is non-Galois. In contrast the degree-three S3
stabilizer normally generates S3. QDI1/QNI1's retained cubic therefore
passes this condition. It remains a sharp degree-three boundary
in rank at least three.

For an explicit non-Galois field with a forbidden Galois factor,
the generic equation
\[
                  z^4+u z^2+v=0,\qquad
                      (z^2)^2+u z^2+v=0                 \tag{10}
\]
contains the quadratic field generated by z^2. The discriminant
u^2-4v is not a square. The quartic has degree four by the generic
polynomial map v=-z^4-u*z^2 over L(u). Its other squared root is
-u-z^2, not a square in L(u,z), so the quartic field itself is
non-Galois. If such a subfield occurred in an actual
centralizer under these hypotheses, individual mates would be
impossible. This is a field certificate, not a constructed Long
centralizer.

## Radical certificates in the scalar Weyl algebra

In the scalar Weyl application there are no new radicals:
\[
             \theta\in H,\quad \theta^k\in Q,\quad k\ge2
                                      \quad\Longrightarrow\quad \theta\in Q.
                                                               \tag{11}
\]
First suppose L contains the kth roots of unity. Q(theta)/Q is
Galois: all roots of the relevant binomial lie in Q(theta), even
if that binomial is reducible. Equation (7) makes this extension
trivial.

For general K, adjoin the kth roots of unity to the constant field.
The classical centralizer identity survives field extension, being
the kernel of finitely many K-linear Poisson maps. The actual
commuting operators and mates extend to the scalar Weyl algebra;
units remain scalar and the same argument applies there. Thus theta
lies in the constant extension Q'=L(mu_k)(t). Since Q'/Q is a
finite abelian Galois constant extension, its intersection with H
is a finite Galois subextension of H/Q. Equation (7) makes that
intersection Q.
This proves (11) without assuming K algebraically closed.

Accordingly a finite, checkable obstruction certificate consists
of commuting b,c in the actual joint centralizer, c nonzero, and
\[
                 q(F)b^k=p(F)c^k,\qquad
                     p/q\notin Q^{\,k}.                  \tag{12}
\]
Commutativity makes theta=b/c well-defined in H. Its kth power is
p/q, so (11) contradicts the nonscalar radical. To certify the last
condition it suffices to exhibit an irreducible polynomial pi
with valuation of p/q not divisible by k. The certificate needs
neither a finite presentation of B nor a description of all H.
A rational commuting operator outside Frac B cannot be substituted
without its own proof of membership.

For the older scalar boundary F1=x+h^2*x^2,F2=y, take
\[
                    b=1+2h^2x,\quad c=1,\qquad
                            b^2=1+4h^2F_1.              \tag{13}
\]
The right side has a simple spectral divisor. This finite square-root
certificate excludes every actual individual mate, including ones
with arbitrary parameter poles. Its rational mate remains outside
the actual generic Weyl algebra, as recorded in QET1.

## Controls for the nonfinite arguments

A multivariate denominator can be forced invertible even when
ordinary polynomial Bezout does not apply. For
\[
 b={x\over x^2+y},\qquad
               (\partial_x-2x\partial_y)b={1\over x^2+y}, \tag{14}
\]
the numerator and denominator both vanish at the origin. The Taylor
flatness proof of (4) handles this case uniformly in all dimensions.

At a singular divisor take
\[
 y^2-x^2-x^3=
          (y-x\sqrt{1+x})(y+x\sqrt{1+x})\quad
                                       \text{in }K[[x,y]]. \tag{15}
\]
The two factors are distinct formal primes, each with multiplicity
one, although the total vanishing order at the origin is two.
This is the distinction required in (5). Exact series prefixes
control the factorization; multiplicity and the all-order square
root identity use the written formal argument.

There are truly non-finitely-generated differential localizations:
\[
 K[[h]][x,y,(1-nh x)^{-1}:n=1,2,3,\ldots].                \tag{16}
\]
They have polynomial h-adic completion and correct classical
reduction, with commuting coordinate derivations and flat Ore mates.
Distinct pole divisors prove that no finite set of algebra generators
can invert them all. They violate (1)'s unit condition and are not
scalar finite Weyl algebras. Thus removing finite generation does
not remove the unit hypothesis or equate completed polynomiality
with finite polynomial realization.

## Evidence and remaining construction gate

The standalone controls check exact Taylor product and differential
identities at a declared point, denominator extraction, regular-sequence
gcd witnesses, singular-divisor factor prefixes, quadratic and radical
certificates, the quartic quadratic subfield, S3/D4 stabilizer normal
closures and selected localization derivatives. Retained QET1/QFC1
hashes are validated. Universal Taylor injectivity, differential
localization, all-divisor valuations, all finite Galois factors,
monodromy, radical descent and finite quantum completion are written
proofs with the stated flatness, factoriality, purity and covering
inputs. No finite prefix proves those universal statements.

For a new actual Long pair, finite generation is no longer required
for the birational, Galois-factor or radical tests. Seek an actual
birational centralizer field and individual mates for positive
completion, or a finite Galois-factor/radical certificate for an
all-partner obstruction. A non-Galois field of degree at least three
can still survive. No new Long field certificate, finite full tuple,
finite carrier, all-representative nonlifting, DC2/JC2 resolution,
independent theorem replay, formal verification, external review or
novelty priority is asserted.

Run the declared controls from the repository root:

```sh
.venv/bin/python research/jacobian-followup/dixmier-formal-lifting/check_taylor_centralizer.py
```
