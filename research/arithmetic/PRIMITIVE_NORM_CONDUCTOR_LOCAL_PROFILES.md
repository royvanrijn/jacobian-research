# Primitive norm images and finite profiles at conductor primes

Arithmetic support for [primitive norm sums in all orders](PRIMITIVE_NORM_ALL_ORDERS.md). All orders and local fields are fixed.

## Scope and primary literature inspected first

Let K be imaginary quadratic with fundamental discriminant d<0 and canonical maximal
basis 1,omega (trace1 for odd d and trace0 for even d). Write T=Tr(omega), U=N(omega), so T²-4U=d. For a fixed
order O_f=Z+f O_K and c=v_p(f), changing y by the p-adic unit f/p^c
identifies its local norm with

    Q_c(x,y)=x²+p^c Txy+p^(2c)U y².                         (C1)

Primitivity here is (x,y) not in pZ_p² in THIS order basis. It is not
primitivity in a maximal-order basis. All images below are full Z_p images,
including the value zero when the local algebra splits.

Closest inspected primary sources:

* Earnest--Gunawardana, On locally primitively universal quadratic forms,
  arXiv:2005.11268v1, actual PDF https://arxiv.org/pdf/2005.11268:
  Section2 Jordan splitting; Examples3.1--3.3 and Proposition3.8 discuss
  primitive binary local images and bounded valuation in anisotropic
  lattices. They do not print the conductor tables below. Their local
  square argument is the same elementary tool used in this reconstruction.
* Bitan--Schein, On the flat cohomology of binary norm forms, actual author
  PDF https://u.cs.biu.ac.il/~scheinm/classnumber69.pdf:
  Definition4.4 and Section4.2 explicitly give quadratic order bases and
  their norm forms. Corollary5.23 discusses squares/principal genera in
  their torsor setting. The latter is NOT imported here as an unqualified
  identification of every ring-class/genus interface, especially its
  stated qualifications at d=-3 and non-fundamental discriminants.
* The retained Milne Class Field Theory norm and Hilbert-symbol statements
  remain the inputs for the MAXIMAL local layers c=0. Conductor layers
  below are derived directly, not by assuming maximal-order norm images.

No literature-priority claim is made. The exact conductor computations,
finite-profile reduction and their limits are written out below.

## 1 Odd conductor primes: the complete table

Suppose p is odd and c>=1. Set eta=(d/p) in {-1,0,1}. For a p-adic unit
u, let [nu] mean its quadratic residue class modulo p. For n!=0 put
e=v_p(n), u=n/p^e. The primitive image of(C1) is exactly:

1. If eta=-1 (inert):
   * e=0,2,...,2c-2 and u a square modulo p; or
   * e=2c and u arbitrary.
2. If eta=0 (ramified):
   * e=0,2,...,2c and u a square modulo p; or
   * e=2c+1 and [u]=[-d/p].
3. If eta=1 (split):
   * e=0,2,...,2c-2 and u a square modulo p;
   * e=2c and u arbitrary if p>=5;
   * e=2c and u a NONSQUARE modulo3 if p=3;
   * every e>=2c+1, with arbitrary unit u.
   Zero is also in the split primitive image and in neither nonsplit image.

In particular nonsplit valuation is bounded by2c or2c+1, whereas split
valuation is unbounded. The split p=3 clause cannot be dropped.

Proof. Complete the square by X=x+p^c Ty/2, a GL2(Z_p) change, obtaining
X²+p^a epsilon y², where a=2c+v_p(d) and epsilon=-d/(4p^v_p(d)) is a unit.
If X is a unit the norm has square unit residue, and every such unit is
represented with y=0. Otherwise primitivity forces y to be a unit.
For k=v_p(X) and2k<a, the norm has valuation2k and square unit quotient.
Every such quotient lifts with y=1 by taking the unit square root of
n/p^(2k)-p^(a-2k)epsilon. For2k>a it has valuation a and epsilon-square
unit quotient; take X=0 for existence.

Only2k=a remains. Then a=2c and n=p^(2c)(z²+epsilon y²), with y a unit.
If the residual form is inert, each nonzero target over F_p has p+1
solutions, at most two having y=0, so every unit has a y-unit root;
the root is nonsingular and lifts. There are no zero roots.
If it splits, a nonzero target has p-1 roots, again at most two with
y=0. For p>=5 this gives every unit. For p=3 direct enumeration with
y!=0 gives exactly the nonsquare units. For every target divisible by
p, write z²-s²y²=AB with A=1, B=target; their difference is a unit,
so y=(B-A)/(2s) is a unit. This gives all higher valuations and zero.
All cases in the table and their converse are now accounted for.

## 2 A complete dyadic diagonal template

Write U_r={u in Z_2^*:u=r mod8}; write U_[r,4] for u=r mod4.
For a>=2 and odd epsilon, put F=X²+2^a epsilon Y². Its primitive image
is the union of the following disjoint valuation families:

(D0) Unit norms: U_1 union U_5 if a=2; U_1 if a>=3.
(D1) For every integer k>=1 with2k<a:

    2^(2k) U_(1+2^(a-2k)epsilon mod8).

(D2) At valuation a:

    2^a (U_epsilon union U_(epsilon+4))   if a is even,
    2^a (U_epsilon union U_(epsilon+2))   if a is odd.

(D3) If a is even, also include2^a B_epsilon, where

    B_epsilon = 2 U_[(1+epsilon)/2,4]      if epsilon=1 or5 mod8,
                4 Z_2^*                 if epsilon=3 mod8,
                8 Z_2                   if epsilon=7 mod8.

All residue subscripts are reduced modulo8 or4 as indicated. There are
no other values. In particular the only unbounded case is even a and
epsilon=7 mod8, the split binary space.

Proof. If X is odd, reduce modulo8; Y even gives residue1 and Y odd gives
1+2^a epsilon. Every listed unit is represented by fixing Y=0 or1 and
using the unit-square criterion. If X is even, Y must be odd. For
k=v_2(X), the2k<a and2k>a cases give(D1) and(D2) by reduction modulo8.
Every listed coset is represented by fixing Y=1 and solving for a square,
or fixing X=0 or2^(floor(a/2)+1) and solving for the unit square Y².
Rational unit scaling gives the whole listed square coset.

If2k=a, then a is even and the remaining factor is z²+epsilon Y² with
both z,Y odd. For epsilon=1 or5 mod8 it has valuation one and quotient
u=(1+epsilon)/2 mod4. Conversely z²=2u-epsilon with Y=1 has an odd square
root exactly under that condition. For epsilon=3 mod8, put z=Y+2t.
The quotient by4 is t²+tY+((1+epsilon)/4)Y². With Y=1 its derivative
in t is odd, and its reduction is every odd target; Hensel gives all
unit quotients. For epsilon=7 mod8, choose an odd s with s²=-epsilon.
The factors z-sY and z+sY have valuations1 and at least2 in either order.
Conversely for a target w in8Z_2 take the factors2,w/2, giving odd z,Y;
this includes w=0. This proves(D3) and necessity.

## 3 Applying the dyadic template to every quadratic order

Let c=v_2(f)>=1.

* If d is odd, d=1 or5 mod8, use the integral unimodular change

      X=x+2^(c-1)y,
      Q_c=X²-2^(2c-2)d y².

  For c>=2 apply Section2 with a=2c-2 and epsilon=-d. For c=1 the
  exponent is zero, and the primitive image is

      Z_2^* union8Z_2         if d=1 mod8 (split),
      Z_2^* union4Z_2^*       if d=5 mod8 (inert).          (C2)

  Indeed exactly one odd diagonal coordinate gives all four odd square
  classes; both odd give B_epsilon with epsilon=7 or3.
* If v_2(d)=2, put epsilon=-d/4 (epsilon=1 or5 mod8). Apply Section2
  with a=2c and this epsilon. All norm valuations are at most2c+1.
* If v_2(d)=3, put epsilon=-d/8. Apply Section2 with a=2c+1. Its
  especially short form is:

      e=0,2,...,2c-2, u=1 mod8;
      e=2c,          u=1+2epsilon mod8;
      e=2c+1,        u=epsilon or epsilon+2 mod8.           (C3)

  There are no other values, so the maximal valuation is2c+1.

For odd d and c>=2, the template says the last fixed-valuation layer is
2c-2 with unit residues3,7 mod8, followed either by the inert layer
valuation2c with ALL unit quotients or by the split tail8*2^(2c-2)Z_2.
Lower even layers carry the singleton square classes prescribed by(D1).
This spells out unbounded split conductor valuations, rather than
silently bounding them by a ramified exponent.

## 4 The actual order-unit norm groups

For c>=1, O_f,p is local with residue field F_p under x+p^c yomega ->x.
An element is a unit exactly when x is a p-adic unit. Its norm is then a
unit, and every order element of unit norm is a unit. Sections1--3 give

    N(O_f,p^*) = {u in Z_p^*: (u/p)=1}       (p odd),     (C4)

and at two

    d odd, c=1:                 all Z_2^*;
    d odd, c=2:                 1+4Z_2;
    d odd, c>=3:                1+8Z_2;
    v_2(d)=2, c=1:              1+4Z_2;
    all other even d,c>=1:      1+8Z_2.                  (C5)

These really are the entire groups. Necessity follows modulo p or8,
and sufficiency follows by multiplying one residue representative by a
rational unit square. Unit multiplication on O_f,p is GL2(Z_p) on the
order basis, so it preserves coordinate primitivity. All rational unit
squares belong to(C4)--(C5).

This differs from the maximal local field norm group. For example an odd
inert conductor prime can occur at valuation2c with EITHER unit class,
although its ORDER-UNIT norm group has only square unit class. A valuation
layer must therefore sometimes be refined into several unit-norm cosets.
Maximal Hilbert symbols alone miss the extra conductor restrictions.

## 5 Finite square-class states on every fixed valuation stratum

Let Bad be the primes dividing2fd. Put L=8 times the product of their
odd members. For any fixed positive prefix r supported on Bad, define

    A_r={a mod L: (a,L)=1 and r*a is a primitive O_f,p norm
                         for every p in Bad}.            (C6)

This depends only on the indicated residue a mod L, despite arbitrarily
large conductor exponent in f or fixed valuations in r. The reason is
not a guessed conductor bound: every local primitive image is stable
under multiplication by a rational unit square. At an odd prime unit
square classes are determined modulo p; at two they are determined
modulo8. Thus every valuation stratum is a UNION of these finite unit
square classes. The tables make these unions explicit. After CRT, A_r
is a finite union of states of the full vector of quadratic characters
mod L. A state can be empty; no nonemptiness is presumed.

Let chi=chi_d modulo L. One has A_r subset ker chi. Indeed local norms
r*a give H_p(a)=H_p(r) at every bad prime; reciprocity gives
product_bad H_p(r)=1 because r is supported there, and the omitted
unramified unit symbols are one. The product of unit Hilbert characters
is chi(a). This compatibility does not eliminate the extra conductor
characters in(C4)--(C5).

## 6 Uniform finite-prefix coverage for additive targets

Here is a general lemma, valid even without the explicit conductor
tables. Let Q be any nonsingular integral binary form, A its integral
Gram matrix, and S_p=Q(Z_p²\pZ_p²). Put

    t=v_p(det A), M=2t+1.

For every primitive z, adj(A) grad Q(z)=det(A)z, so at least one
partial derivative has valuation r<=t. If n-Q(z) is divisible by
p^M, the one-variable Hensel theorem applies to that coordinate with
error valuation at least2r+1. The new coordinate differs by a multiple
of p^(r+1), so the pair remains primitive. Consequently

    S_p is a union of complete cosets modulo p^M.          (C7)

In particular, if its zero residue is allowed, all of p^M Z_p is in S_p.
This allows a bounded-valuation replacement in EVERY additive local
solution a+b=N. Keep any nonzero residue summand modulo p^M, whose
valuation is already <M. If a residue is zero, the full p^M ball is
allowed, and replace as follows while preserving the original residues:

* If exactly one residue is zero, choose its value p^M, or p^(M+1);
  choose the other value N minus it. The other residue stays unchanged
  and has valuation<M. The zero-residue value has valuation M or M+1.
* If both residues are zero and v_p(N)>M (including N=0), choose
  a=p^M,b=N-p^M; both have valuation M.
* If both residues are zero and v_p(N)=M, choose
  a=p^(M+1),b=N-p^(M+1); their valuations are M+1 and M.

Thus every pair-locally-admissible target has a local representation
with NONZERO summand values and

    v_p(a),v_p(b)<=M+1=2v_p(det A)+2.                    (C8)

No assertion that ALL individual primitive norms obey this bound is
made. Split conductor places have the unbounded tails already exhibited.
For the order form det A=-f²d, so(C8) is an explicit fixed bound.
At good primes no additional prefix is required. Taking the product of
these finite exponent ranges over Bad gives finitely many prefixes r,s.
The local solutions give a,b in A_r,A_s for

    (r/g)*a+(s/g)*b=N/g mod L,  g=(r,s), g|N.             (C9)

The all-state quadratic splitting filter applies to each
finite pair of square-class states in A_r,A_s and the coprime bad
coefficients r/g,s/g. Hence(C6)--(C9) are sufficient local profiles for
a positive lower genus-shadow convolution. All target valuations are
retained; no growing-conductor or growing-state assertion is made.

A sharper conductor-specific cap also follows from the tables. At a
nonsplit conductor prime every primitive norm already has valuation
<=2c+1. At a split conductor prime the full p^(2c+1) ball is allowed,
and the low layers are stable under additions from that ball. The same
replacement argument gives the sufficient bound2c+2 there. The larger
uniform bound(C8) avoids case bookkeeping and is enough for the theorem.

## 7 Proper ideal prefixes and consistency of class states

This section records the interface, not an unproved replacement of
ring-class genus theory by maximal-order Hilbert symbols.

For fixed bounded norm valuation e, primitive nonzero local elements
have finitely many orbits under O_f,p^*. In a nonsplit local field,
the element valuation is determined by the norm valuation. In the split
algebra, its two nonnegative component valuations sum to e. There are
finitely many possibilities, and in every case the remaining unit
quotient is finite because O_f,p^* has finite index in O_K,p^*:
1+p^c O_K,p is already contained in it. This proves the finite-orbit claim.

Choose local primitive representatives alpha_{j,p} for a selected finite
prefix profile and form the integral order ideal I_j with these local
components alpha_{j,p} O_f,p and components O_f,p elsewhere. This ideal
is locally principal and hence invertible. Its norm r_j is the product
of the p-adic norm valuations, not a claim that r_j has a principal
generator. Put

    beta_{j,p}=N(alpha_{j,p})/r_j in Z_p^*/N(O_f,p^*).       (C10)

The coset is independent of the chosen generator of the same local
ideal. For a good norm m, principalization of I_j times an oriented
prime-to-Bad ideal of norm m produces a primitive principal norm r_j*m
only with the corresponding unit condition m in beta_{j,p}N(O_f,p^*).
Conversely a principal ideal with this local type has a generator whose
local primitive coordinates are preserved by order-unit multiplication.
The local layers in(C6) can be refined into these finitely many cosets
when explicit ring-class labels are required.

One must NOT define a general ring-class character by the maximal
formula H_p(N I) without checking independence at the conductor. The
extra order-unit conditions are visible in(C4)--(C5), and arbitrary
fractional principal ideals need not have order-unit generators at the
conductor. A correct proof can use coprime representatives or the
idelic/local-unit-norm quotient and establish that its genus kernel is
2Pic(O_f). This is proved by OrderGenusKernel and OrderPrimitiveShadow in the all-orders note.

For the class-error UPPER bound no finite valuation cap is necessary.
After a primitive local ideal type for a value has been fixed, good
prime-ideal orientation flips preserve every bad local ideal and its
primitivity, even when a split conductor valuation is arbitrarily large.
If the independently proved ring-class genus criterion supplies one
primitive invertible ideal in2Pic(O_f), the generic signed saturation
argument then applies unchanged using primes away from2fd. Every hidden
failure is deficient in one exact-one good class-prime family. Dropping
ALL bad-prime conditions enlarges the upper set, and the random-marked
prime-square sieve uses only inert primes away from Bad. The omitted
finite prime set changes fixed constants and the definition of J only;
it causes no unbounded target-valuation loss.

## 8 Exact finite challenges and concrete warning examples

The portable checker verify_primitive_norm_order_profiles.py
compares the tables with complete finite Taylor-lifted primitive coordinate
images in99 cases, involving282260 target residues. It includes every
c=1,...,5 at two for d=-3,-4,-7,-8,-11,-15,-20,-23,-24,-40,-56,
and the odd cases(p,c)=(3,1),(3,2),(5,1),(7,1) for those fields.
All comparisons pass. These finite checks challenge, but do not replace,
the all-depth proofs above.

Concrete exact examples:

* d=-8,f=3,Q=x²+18y² at3:9 is NOT a primitive local norm;
  18=Q(0,1) and27=Q(3,1) are primitive. The missing square-unit
  threshold at split3 is real.
* d=-4,f=2,Q=x²+4y² at2: primitive norms have valuations0,2,3,
  with unit conditions from Section2; no primitive norm has valuation4.
  In particular16 is not a primitive local norm although4=Q(0,1)
  and8=Q(2,1) are.
* d=-7,f=2: Q=x²+2xy+8y² is equivalent to X²+7y². Every valuation
  at least3 occurs primitively at2, while valuations1,2 do not.
  Thus excluding all high conductor valuations would be incorrect.
* d=-3,f=2: Q=x²+2xy+4y² is equivalent to X²+3y².
  Its primitive2-adic image is all odd units together with4Z_2^*;
  valuation3 and above are absent. Split/inert behavior cannot be merged.

The same checker, with receipt primitive_norm_order_profile_checks.json,
checks26 explicit orders, including Gaussian conductors5 and25, with52 actual
local images,1850 complete profile families and4098200 target residues across
the complete stated periods. Every full profile mask equals an independent
primitive-coordinate local selfsum mask. It checks uniform unit-square-state
periods, the selected caps against the coarse Hensel cap, and52
full-mask corruptions. These are finite validations, not a theorem over all
discriminants.

There is a particularly important nonmaximal warning at discriminant-100.
The primitive forms x²+25y² and2x²+2xy+13y² both represent25, via(0,1) and
(2,1). Their nonzero unit images modulo5 are respectively{1,4} and{2,3},
so they belong to different local genera. Thus a numerical conductor norm
need not determine a unique genus. The correct class-error input is existence
of ONE primitive invertible ideal in the principal genus, never a claim that
all ideals with that norm have the same genus.

No conductor-order global exact-sumset classification, leading constant,
effective cutoff or historical novelty assertion follows merely from
these local calculations. The finite profiles and local-unit groups are
ready for the separate genus, analytic convolution and class-error modules.
