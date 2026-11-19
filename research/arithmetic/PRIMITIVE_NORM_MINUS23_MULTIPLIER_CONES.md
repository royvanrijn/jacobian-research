# Universal multiplier cones for primitive discriminant -23 sums

Let (Q(x,y)=x^2+xy+6y^2), and let (P) be its positive
primitive numerical values. Write (chi(p)=(-23/p)), and let


    S={positive integers supported on primes p with chi(p)=1},
    C=S union 23S.

Thus q belongs to C exactly when it has no inert prime factor and
v23(q)<=1. The multiplier q need not itself be a primitive or ordinary
principal norm.

**Theorem CN3MultiplierCones.** For every q in C,

    16q, 19q, 39q belong to P+P.

In each case one of the following three fixed seed pairs works after
multiplication by q:

| Seed | Pair A | Pair B | Pair C |
|---|---|---|---|
| 16 | (8,8) | (3,13) | (4,12) |
| 19 | (1,18) | (3,16) | (6,13) |
| 39 | (3,36) | (12,27) | (13,26) |

For example, the theorem settles every target whose inert part is
exactly19 and whose23-adic valuation is at most one: all other split
prime classes and exponents are arbitrary. It also settles all
inert-free targets divisible by16, and all inert-free targets divisible
by39, with the same ramified restriction. In particular every power
of two at least16 is represented.

## The complete finite class reduction

Use the retained [C3 primitive criterion](PRIMITIVE_NORM_CLASS_NUMBER_THREE.md),
Lemma CN3Individual. For a positive n with no inert prime factor and
v23(n)<=1, let r(n) count DISTINCT nonprincipal split primes whose
whole exponent is nonzero modulo three. Then

    n belongs to P if and only if r(n)!=1.                 (MC1)

Both r=0 and r>=2 are allowed. The distinction concerns whole
prime-power orientations; an ordinary nonprimitive representation
does not replace (MC1).

The primes2,3,13 are split and nonprincipal. The nonprincipal form
2x²+xy+3y² represents them at(1,0),(0,1),(2,1), whereas the
principal form has no representation of any of them. The exact checker
verifies these small facts with complete coordinate bounds. It also
checks that29 is split and nonprincipal: the nonprincipal form represents
it at(-2,3), and a complete principal-form check finds no representation.

Every summand in the seed table has factors only in this bank. Set

    e=(v2(q),v3(q),v13(q)) modulo3,
    R=number of active nonprincipal split primes of q outside the bank.

Principal split primes of q contribute no class term. If a seed summand
u has bank-exponent vector a(u), then

    r(qu)=R+#{i:(e_i+a_i(u)) mod3 !=0}.                   (MC2)

This is an equality: multiplying by u changes no outside-bank
exponent. Since every seed summand is prime to23 and has no inert
factor, qu satisfies all the other hypotheses of (MC1), including
when23 divides q once.

If R>=2 every pair in the table works, because both final active
counts are at least two. It remains to check exactly R=0 and R=1
for the27 bank phases. The following table records a valid pair label
for every case, using (MC2). It is a complete unbounded reduction,
not a scan through a finite collection of integers q.

**Lemma CN3DisplayedPairCover.** For any specified finite list of
positive bank-smooth pairs with sumM, the complete R=0,1 phase cover
is equivalent to the assertion that, for every q in C, one of THOSE
displayed pairs (u,v) has both qu and qv in P.

The forward direction is (MC2), with R>=2 automatic. For the reverse,
realize any bank phase e with q0=2^e2*3^e3*13^e13, using
e_i in{0,1,2}. This q0 belongs to C and has outside active count
zero. The multiplier29q0 belongs to C with the same bank phase
and outside active count exactly one. The asserted displayed-pair
property and (MC1) then give the required finite witnesses for both
cases. This proves an equivalence for the specified pair list; it does
not assert that every universal P+P membership theorem must use
those pairs.

| e2 | e3 | e13 | 16:R0 | 16:R1 | 19:R0 | 19:R1 | 39:R0 | 39:R1 |
|---:|---:|---:|:---:|:---:|:---:|:---:|:---:|:---:|
| 0 | 0 | 0 | A | B | A | B | B | A |
| 0 | 0 | 1 | C | A | B | A | A | A |
| 0 | 0 | 2 | B | A | B | A | A | A |
| 0 | 1 | 0 | C | A | C | A | C | A |
| 0 | 1 | 1 | A | A | A | A | A | A |
| 0 | 1 | 2 | A | A | A | A | A | A |
| 0 | 2 | 0 | B | A | B | A | A | B |
| 0 | 2 | 1 | A | A | A | A | B | A |
| 0 | 2 | 2 | A | A | A | A | B | A |
| 1 | 0 | 0 | B | A | C | A | C | A |
| 1 | 0 | 1 | A | A | A | A | A | A |
| 1 | 0 | 2 | A | A | A | A | A | A |
| 1 | 1 | 0 | A | A | B | A | A | B |
| 1 | 1 | 1 | A | A | A | A | B | A |
| 1 | 1 | 2 | A | A | A | A | B | A |
| 1 | 2 | 0 | A | A | A | A | B | A |
| 1 | 2 | 1 | A | A | A | A | A | A |
| 1 | 2 | 2 | A | A | A | A | A | A |
| 2 | 0 | 0 | B | A | B | A | A | A |
| 2 | 0 | 1 | A | A | A | A | A | A |
| 2 | 0 | 2 | A | A | A | A | A | A |
| 2 | 1 | 0 | A | A | A | B | B | A |
| 2 | 1 | 1 | A | A | B | A | A | A |
| 2 | 1 | 2 | A | A | B | A | A | A |
| 2 | 2 | 0 | A | A | C | A | C | A |
| 2 | 2 | 1 | A | A | A | A | A | A |
| 2 | 2 | 2 | A | A | A | A | A | A |


For each cell, the corresponding two counts in(MC2) are different
from one; (MC1) therefore supplies positive primitive principal norm
representations of BOTH qu and qv. Since u+v is the seed, their sum
is16q,19q or39q. This proves the theorem, including all split
prime powers and the optional single ramified factor23.

## Direct partition for the seed16

The first cone also has a short partition without the phase table.
If r(q)!=1, multiplication by8=2³ preserves every active exponent
residue, so8q+8q works. Otherwise q has exactly one active
nonprincipal split prime p.

If p is different from3 and13, multiplying by3 or13 adds a second
active prime, so3q+13q works. The same pair works if p is3 or13
with exponent two modulo three: one summand cancels that active
prime to r=0, and the other has two active primes. In the remaining
cases p is3 or13 with exponent one modulo three. Use4q+12q.
For p=3 its two active sets are{2,3} and{2,3}; for p=13 they
are{2,13} and{2,3,13}. This exhausts all possibilities.

## Scope and exact certificate

The [stdlib checker](verify_primitive_norm_minus23_multiplier_cones.py)
and [receipt](primitive_norm_minus23_multiplier_cone_checks.json)
record all81 phase cases per seed:27 phases times R=0,1,>=2.
They verify the bank factorizations, sums, active counts and the small
nonprincipal-prime controls. Deleting any one of a seed's three
listed pairs produces an explicit uncovered phase, retained as a
negative control. This is necessity within the displayed pair list;
no minimality among other possible constructions is claimed.

This proof uses only the exact individual C3 criterion and finite
arithmetic. It uses no Hecke estimate, analytic asymptotic or numerical
large-target cutoff. A multiplier containing an inert prime is outside
the theorem's hypotheses. The result does not settle the complete
40-exception sumset classification or the remaining middle interval.
Internal proof and finite replay do not establish external human
review, formal verification or historical priority.
