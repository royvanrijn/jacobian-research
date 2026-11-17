# Oriented pattern linking excludes the complementary source in the all-odd class

15 November 2026. PF2RULED66 strengthens PF65 for the ENTIRE actual
primitive-binomial two-stage class with coprime odd first pair p,q>=3.
Every remaining quartic end forces n_R=2 mod3, where -n_R is the
ramified diagonal of the full source lattice. The complementary
source configuration is therefore impossible at ALL source orders
and unramified pole heights. The only remaining source configuration
in this target class has the COMMON intrinsic coordinate pencil,
with the SAME source boundary value on R and D, and n_R=2 mod3.

The argument adds the orientation of the torsion linking form to
PF65's homology and zero ramification mark. Its general substitution
lemma applies to homology-solid-torus pieces at any cover degree.
No coefficient search, contact-filling classification, negative
plumbing claim about the full satellite, or assumed lens end is used.
The full quartic and higher-degree Keller problems remain open.

## 1. Integral homology torus replacement preserves oriented linking

Let E be a compact oriented integral homology solid torus, with
boundary basis (m,l) such that m generates H1(E)=Z and l bounds an
integral surface in E. Glue E to an oriented compact 3-manifold P
along one of P's boundary tori. Replace E by a solid torus with the
same boundary identification and the same killed longitude l.
If the resulting closed manifold M is a rational homology sphere,
the induced homology isomorphism preserves its torsion linking form.

Here is a chain proof. The boundary meridian generates H1(E), so
Mayer--Vietoris represents every H1(M) class by a cycle in P. For a
multiple of a torsion cycle, a bounding 2-chain in P can have extra
boundary only in the kernel Z*l of H1(boundary E)->H1(E); cap these
longitudes by integral surfaces in E. The second cycle can stay in
the interior of P. All intersections used to compute linking are
then in P, and changing the caps inside E does not change them.
Exactly the same construction works in the replacing solid torus.
It also preserves the class of any marked knot lying in P. Apply
this one piece at a time if there are several such attachments.

PF65 gives precisely this peripheral basis for a cyclic a-sheeted
companion exterior: m=M^a, l=L. When H1(Sigma_a(companion))=0 it
is an integral homology solid torus. Consequently in ALL finite
cover degrees, every such cyclic companion replacement preserves
the oriented end linking form, as well as homology and outer
markings. Nontrivial companion homology cannot simply be deleted.
This lemma asserts no homeomorphism or contact equivalence between
the two ends.

## 2. An explicit complex model fixes the pattern orientation

Let s,t>=1 satisfy gcd(3s,4t)=1. Consider the normal surface germ

    X_(s,t): w^4+u^t*w+v^s=0,
    pi(u,v,w)=(u,v).                                         (1)

It is smooth away from the origin, and at s=1 or t=1 the possible
smoothness at the origin causes no problem. The derivative equations
show that any singular point is the origin; a reducible surface germ
would have a positive-dimensional intersection of components and
hence a singular curve. Thus the germ is integral and normal. The
map pi is finite of degree4, with discriminant

    256*v^(3s)-27*u^(4t).                                   (2)

Its oriented branch knot is the POSITIVE torus knot T(3s,4t).
A small base sphere pulls back to the positively oriented link of
X_(s,t): the positive weighted radial field is transverse to these
level sets, so this choice of link agrees with the complex boundary
orientation. It is not the orientation-reversed link.

Along u=0 the four roots are permuted in a 4-cycle, since s is odd;
along v=0 there is one fixed root and a 3-cycle, since 3 does not
divide t. Local branching is simple. The meridian has image a
transposition. Taking the coprime generator powers in a torus
meridian reduces to a 3-cycle and a 4-cycle whose product is a
transposition; S4 has one simultaneous-conjugacy class of these
pairs. Thus this is the unique pattern class up to relabelling.
This identifies (1) with the unknotted pattern cover in PF65 whenever
the parameters match, including its orientation.

The weighted circle action has primitive weights

    wt(u)=3s, wt(v)=4t, wt(w)=st, degree=4st.                  (3)

There are three exceptional orbit types: on v=0,w^3+u^t=0 the
stabilizer order is s; on v=w=0 it is 3s; on u=0 it is t. An
order1 type is omitted as a regular orbit. The base is S2. For
example the degree4 map of base surfaces has branch partitions
(3,1),(4),(2,1,1), with total ramification6; Riemann--Hurwitz gives
2g-2=4*(-2)+6=-2. Equivalently, in the applications PF65 already
gives finite H1, which forces genus zero for a Seifert base.

We retain the following explicit orientation convention. In a
negative star plumbing, the central weight is -b and a leg has
negative continued fraction -alpha/beta, with 0<beta<alpha.
The rational Euler number is

    e=-b+sum beta_i/alpha_i.                                 (4)

Neumann--Raymond, complete PDF3-4/printed165-166, defines the Euler
number and proves its naturality under oriented fibre-preserving
maps, including branched maps. Their nonnormalized leg parameter
is minus the positive beta used here. Complete PDF11-13/printed173-175
gives the star plumbing, continued fractions and intersection form.
We convert these signed leg conventions explicitly; a group
presentation without this conversion would not fix linking signs.

The weighted Hopf fibration of the base S3 has Euler number
-1/((3s)*(4t)). Indeed the power map from the ordinary Hopf S3 to
this weighted S3 has base degree(3s)*(4t) and fibre degree1; the
ordinary positive Hopf circle bundle has Euler number -1. The
degree4 projection (1) has base degree4 and fibre degree1. Euler
naturality therefore gives

    e(X_(s,t))=-1/(3st).                                    (5)

For the order-s orbit, the transverse complex coordinate is v,
with stabilizer character 4t modulo s. Its positive normalized
negative-plumbing parameter satisfies

    beta=-(4t)^(-1) mod s.                                  (6)

The signed Neumann--Raymond parameter is the inverse of the normal
slice rotation. Passing from their signed leg continued fraction
to the negative continued fraction -s/beta negates that parameter,
giving (6). This convention is also tested against the A-chain
orientation below; a signless stabilizer order would not suffice.

The other normalized parameters can be taken as

    beta_2=-(4t)^(-1) mod(3s),
    beta_3=-(3s)^(-1) mod t,
    b=beta/s+beta_2/(3s)+beta_3/t+1/(3st).                    (7)

For t=1 omit the third leg. These expressions have integral b:
beta_2=beta mod s makes 3*beta+beta_2 equal to 4*beta_2 mod(3s),
so the numerator in (7) vanishes mod(3s); it also vanishes mod t.
The two moduli are coprime.
The normalized inverses exist by gcd(3s,4t)=1. The plumbing has
negative definite intersection form because its legs are negative
and its central Schur complement is the negative number (5).
We use this plumbing for the PATTERN link, not for the full source
or the full satellite after a topological substitution.

## 3. The pattern linking form is <-t/s>

Use the outward-boundary convention of PF61: a negative plumbing
matrix Q gives linking -Q^(-1) modulo Z. For the first vertex of
the order-s leg, continued-fraction elimination gives

    (Q^(-1))_11=-beta/s+(beta/s)^2/e.

Combining with (5), its boundary self-linking is

    beta*(1+3t*beta)/s.                                      (8)

By (6), multiplying this generator by the unit 4t changes its
self-linking numerator modulo s to

    16t^2*beta*(1+3t*beta) = -t mod s.                       (9)

In particular this element has order s and generates the cyclic
H1 of order s. Thus the oriented pattern linking form is

    lambda_pattern is isometric to <-t/s>.                   (10)

Here <a/s> denotes the nonsingular cyclic linking form on Z/s
whose generator has self-linking a/s modulo Z. Its coefficient
is meaningful up to multiplication by a unit square.

There is an independent orientation control: at t=1, the polynomial
coordinate change U=u+w^3, W=w, V=v makes (1) exactly UW+V^s=0.
This is A_(s-1), with negative
(-2)-chain resolution. Its outward linking is <-1/s>, agreeing
with (10). Reversing orientation would give <+1/s>, a genuinely
different form on its 3-primary part. Thus this sign is material.

## 4. Transfer to the surviving all-odd satellite class

Retain PF63/PF65: p,q>=3 coprime odd, r>=2, c=pq-p-q, k=rc+1,
and an ACTUAL primitive-binomial two-stage key. No general lower
weight key or longer sequence is assumed to have this form.
The only remaining peripheral alternatives are

* C4 companion: r=9 mod36. Set D=s=r/3, t=k/4.
* C3 companion: 4|r,9|k and 3 does not divide pq. Set D=s=k/3,
  t=r/4, interchanging the two torus parameters in the pattern.

The C3 cases with 3|pq were excluded in PF65 and are not reinstated.
The remaining cyclic companion pieces have ZERO cyclic branched
homology: Sigma4(T(p,q)) in C4, Sigma3(T(p,q)) plus the degree1
piece in C3. Section1 therefore transfers (10) to the FULL satellite
end, retaining the zero ramification-knot class proved in PF65.

In C4, t=k/4 equals 1/4 modulo D, since k=rc+1 and D=r/3.
This is a unit square, so its full end linking is <-1/D>.
In C3, (10) gives <-(r/4)/D>; equivalently it is <c/D>, since
r*c=-1 mod D and 4 is a square. The complete forms are

| Remaining alternative | D=order H1 | Oriented linking form | Ramification class |
|---|---|---|---|
| C4 | r/3 | <-1/D> | 0 |
| C3, 3 does not divide pq | k/3 | <-(r/4)/D>, equivalently <c/D> | 0 |

Both have D odd and divisible by3. In C4, t=1 mod3 because r=0
mod3 and 4t=k=rc+1. In remaining C3, k=0 mod3 forces c!=0 mod3.
For p,q both nonzero mod3, c is either0 (both p,q=2 mod3) or2;
the first case is impossible. Thus c=2 mod3 and r=1 mod3, so
t=r/4=1 mod3 as well. In either remaining alternative the coefficient
of (10) is a NONSQUARE modulo3.

## 5. All complementary source orders are excluded

The actual source matrix and canonical covector are

    Q=[[-n_R,c_src],[c_src,-n_D]],  k_source=(1,0).

PF65 gives [K_R]=0. Hence e_D generates H1(end)=Z/D and PF61
gives its self-linking as n_R/D. Matching (10) requires

    n_R=-t*u^2 mod D for some unit u modulo D.                (11)

Restrict (11) modulo3. Since t=1 mod3 and every unit square is1
there, the actual ramified dual diagonal MUST satisfy

    n_R=2 mod3.                                             (12)

This is stronger than PF65's gcd(n_R,D)=1. It is independent of
every source principal-part coefficient, pole cutoff and choice
of primitive inverse equation.

PF50's complete coordinate dichotomy says the source intrinsic
pencils are common or complementary. In the complementary case:
source order r_s=1 has q_s=2 and c_src=2*p_D, which is even,
contradicting the odd D determinant as already noted in PF65.
For EVERY r_s>=2 it has q_s=r_s-1 and

    n_R=r_s^2-r_s+1.

This is0 or1 modulo3 for every integer r_s, never2. Equation(12)
therefore excludes the ENTIRE complementary source configuration
for this ENTIRE all-odd target class, with arbitrary unramified
orders/poles and all source automorphisms retained. We do not
infer that the source end is a lens space from its homology,
its ruling, or its number of omitted affine lines.

## 6. The surviving common pencil has coincident boundary values

Here is a general useful orthogonality fact for source valuations.
Let two omitted valuations have the same regular source-coordinate
pencil f, but distinct finite constants f=a and f=b. Use a common
completion of P1_x times P1_f with x a complementary source
coordinate. The extraction centres over x=infinity,f=a and over
x=infinity,f=b are disjoint. If the respective x pole orders are
q,p, their integral duals push forward to q*F and p*F, where F is
the fibre class of f. To check this, a generic x fibre pairs with
the dual by q or p (the negative valuation of x-beta); a generic
f fibre pairs by0. Since F^2=0 and the exceptional clusters are
disjoint, their full dual cross pairing is ZERO. Further boundary
blowups preserve it. No ramification-index assumption is needed
for this calculation.

Apply this to R,D in the common quartic pencil. Distinct boundary
constants would give c_src=0. With zero ramification class the
diagonal presentation would then force n_R=1: equivalently
D=n_R*n_D and D|n_D. But the actual index2 whole-chart formula
gives n_R>=3. Indeed r_s=1 has n_R=2q_s+1>=5, and r_s>=2 has
S>= (r_s-1)q_s and q_s>=r_s-1, so

    n_R=S+2q_s-r_s+2 >= r_s^2-r_s+1 >=3.

Thus the common pencil must have the SAME boundary value on R,D.
This is a value of an intrinsic SOURCE coordinate, not an equality
of the distinct target image curves C,H. In all degrees with a
sole fold prime and unramified primes in one common pencil, the
same orthogonality shows that a zero canonical mark requires at
least one unramified prime over the ramified prime's source value;
otherwise its lattice block is isolated and its mark is nonzero.
The smoothness and integral-lattice hypotheses of PF61 are retained.

## Evidence, abandoned routes and the remaining obligation

[verify_oriented_pattern_linking.py](verify_oriented_pattern_linking.py)
checks the unrestricted Schur-complement identity, modular generator
rescaling, exact continued-fraction plumbing matrices, cyclic linking
and canonical-mark controls, the A-chain orientation control, both
remaining parameter alternatives, and the source congruence and
different-value orthogonality controls. The universal arguments are
written above; finite samples are not an interpolation proof.
[oriented_pattern_linking_receipt.json](oriented_pattern_linking_receipt.json)
seals the inputs. [The primary manifest](sources/oriented-pattern-linking/manifest.json)
pins Neumann--Raymond, with complete PDF3-4 and11-13 visually read;
complete PDF2,5,6 were also read for conventions/context and the
signed inverse slice parameter. No broad
symplectic/contact filling theorem is imported.

Exploratory ruling and contact approaches supplied no further
theorem: an A1 ruling does not by itself make this end a lens space,
and homology replacement does not identify contact structures.
Those inferences are not adopted. The exact linking comparison is
what proves the exclusion above.

For the remaining all-odd branch class, the source now has a common
pencil, coincident boundary value and n_R=2 mod3, as well as the full
square-class constraint(11), zero canonical boundary mark and PF61's
integral determinant conditions. These do not supply an affine cover.
The common-pencilled affine controls still have no actual finite
Keller algebra. Even first pairs, non-binomial keys, general longer
sequences and higher-degree exclusions remain open. No full goal
completion, independent full topological reproof, formal verification,
external review, novelty priority or finite Keller construction is
claimed.
