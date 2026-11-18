# Zero-branch ramification obstructs the degree-256 Hasse construction

The economical element cover H1,H3,W5 does not give an everywhere locally
soluble family on any full arithmetic progression of degree-ten parameters.
In fact its locally soluble integer parameters have density zero in every
arithmetic progression. The obstruction comes from a density-7/16 set of
primes at the single rational branch value s=0.

This does not rule out a thinner infinite family. It identifies the
additional arithmetic restriction that such a family must satisfy and
shows why the rank-301 construction cannot be reduced to rank 256 merely
by imposing finitely many congruences.

Use the [regular degree-ten pencil](RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md)

\[
 P_s(z)=9z^{10}-864z^8+28161z^6-350150z^4
                  +1127500z^2-sz+250000.
\]

For a full-A10 specialization with splitting field M_s/Q, let

\[
 A_s^{256}=M_s^{H_1}\times M_s^{H_3}\times M_s^{W_5}.
\tag{1}
\]

The subgroups and the local fixed-coset criterion are defined in
[the subgroup-cover note](ALTERNATING_SUBGROUP_COVERS_AND_HASSE.md),
Theorems 2 and 5. The three field degrees are 10,120,126. Write

\[
\begin{split}
R(u)&=9u^3-414u^2+1836u+400,\\
g(z)&=R(z^2)=9z^6-414z^4+1836z^2+400,\\
B&=\{2,3,5,13,17\}.
\end{split}
\tag{2}
\]

All finite algebra and group counts below have an
[exact checker](verify_alternating_ramification_density.py), an
[actual replay receipt](alternating_ramification_density_verification.json),
and a [source ledger](alternating_ramification_density_source_ledger.json).
The density theorem and integral Hilbert theorem are explicit imports.

## Theorem 1. The residual sextic has full signed-pair group

The splitting field of g over Q has group

\[
                            C_2\wr S_3
\]

of order 48 in its natural six-letter action. The set E of primes p outside
B for which every irreducible factor of g modulo p has even degree has
natural density

\[
                          \frac{21}{48}=\frac7{16}.
\tag{3}
\]

The possible factor-degree patterns in E, and their densities among all
primes, are

| Pattern | Group elements | Prime density |
| --- | ---: | ---: |
| (2,2,2) | 7 | 7/48 |
| (2,4) | 6 | 1/8 |
| (6) | 8 | 1/6 |

Proof. Exact discriminant calculations give

\[
\begin{split}
 \operatorname{Disc}R&=2^6 3^6 5^5 13^2 17,\\
 \operatorname{Disc}g&=-2^{22}3^{14}5^{12}13^4 17^2.
\end{split}
\tag{4}
\]

The monic reduction of R at 7 is u^3+3u^2+u+4, which has no root in F7.
It is irreducible. The cubic discriminant has nonsquare class 85, so its
splitting field L/Q has group S3. Its discriminant is positive; its three
roots u1,u2,u3 are real and distinct, and L is contained in R.

Let M be the splitting field of g. Its roots come in the three pairs
{sqrt(ui),-sqrt(ui)}. Thus Gal(M/Q) embeds into C2 wr S3 and maps onto
Gal(L/Q)=S3. The kernel N is an S3-invariant subgroup of the sign space
F2^3.

The cubic roots have product -400/9 and sum 46. Since they are all real
and nonzero, an odd number is negative. All three cannot be negative
because their sum is positive. Exactly one is negative and the other
two are positive. Complex conjugation fixes L and flips precisely the
square-root pair belonging to the negative root. It is a weight-one
element of N. Conjugating by the surjective S3 action gives all three
basis flips. They generate F2^3, so N=F2^3 and the Galois group is the
full wreath group of order 48.

As a separate exact control, at 29 one has the squarefree factorization

\[
                9^{-1}g(z)=(z^2+3)(z^4+9z^2+3)
                         \quad\text{in }\mathbb F_{29}[z],
\tag{5}
\]

with both factors irreducible. It exhibits one of the obstructing residual
patterns but is not needed to prove completeness of the Galois group.

Enumerate its 6*8 signed permutations of the three pairs. There are
exactly 7 elements of type (2,2,2), 6 of type (2,4), and 8 of type (6).
Those are exactly its elements with all cycles even. The property is
conjugacy-invariant. Apply the inspected primary statement of Chebotarev:
Lenstra, *The Chebotarev Density Theorem*,
[Fact 2.1 and Theorem 3.1, page 3](sources/lenstra-chebotarev-density-notes.pdf).
It identifies factor degrees with Frobenius cycle lengths and gives each
conjugacy class its size divided by the group order as natural density.
This proves (3) and the table. For the nonmonic integral g one can use
the monic polynomial 9^5*g(T/9); outside 3 the variable scaling does not
change the factor degrees or splitting field.

The modular irreducibility controls use exact modular powers and gcds.
The real-root and complex-conjugation argument above supplies the full group; no inference of
that group from a single factorization is made.

## Theorem 2. Exact local criterion above the zero branch

Let s be a nonzero integer, p outside B, and v_p(s)>0. Assume the global
specialization has full A10, so (1) is defined. Then

\[
 \operatorname{Spec}(A_s^{256})(\mathbb Q_p)=\varnothing
 \quad\Longleftrightarrow\quad
       p\in E\ \text{and}\ v_p(s)\ \text{is odd}.
\tag{6}
\]

In the obstructed case the decomposition-group orbits on the ten roots
have sizes (2,2) together with the factor degrees of g modulo p. Inertia
is generated by a double transposition. The possible groups are C2 x C2,
C2 x C4, and C2 x C6 for the three residual patterns in Theorem 1.

Proof. At s=0 one has

\[
 P_0(z)=(z-5)^2(z+5)^2g(z),\qquad
 g(5)=g(-5)=-5^2 13^2 17.
\tag{7}
\]

The three monic factors (z-5)^2, (z+5)^2, and g/9 are pairwise coprime
modulo p, and the last factor is squarefree. Hensel factorization of
P_s/9 gives monic factors Q_+,Q_-,G_s of degrees 2,2,6 reducing to those
three factors. The splitting field of G_s/Q_p is unramified; the orbit
lengths on its six roots are the factor degrees of g modulo p.

Put v=v_p(s). Uniqueness of the coprime Hensel factorization modulo p^v
gives Q_+(z)=(z-5)^2 modulo p^v and
Q_-(z)=(z+5)^2 modulo p^v. For r=5 or -5, the other two factors are
units at r, while P_s(r)=-sr has valuation v. Thus the translated
quadratic is

\[
                  Q_r(r+Y)=Y^2+b_rY+c_r,\qquad
                  v_p(b_r)\geq v,\quad v_p(c_r)=v.
\]

Since p is odd, its discriminant b_r^2-4c_r has valuation exactly v.
When v is even, both quadratic splitting extensions are unramified
(possibly trivial). The whole local splitting field is then unramified,
its decomposition group is cyclic, and the element cover supplies a
point on one factor of (1).

Suppose v is odd. Both quadratic factors are irreducible and ramified.
Their splitting extensions are tame quadratic extensions. The compositum
with the unramified residual splitting field has cyclic inertia of
order two, acting nontrivially on each quadratic and trivially on the six
residual roots. Its permutation is a double transposition d. Since each
quadratic factor is defined over Q_p, its pair is invariant under the
whole decomposition group D. Those are two separate D-orbits of size two.

If a residual factor has odd degree, a partition of six containing an
odd part necessarily contains a part of size one or three. That D-orbit
gives a fixed sheet or an invariant three-subset, so D lies in a
conjugate of H1 or H3.

Otherwise all D-orbits have even size, and D preserves no subset of size
one, three, or five. It also preserves no unordered partition into two
five-element blocks. Indeed d fixes the six residual letters, so it
must preserve each block individually; its two transposition pairs
must each be monochromatic. Every element of D preserves each of those
pairs as a set. It therefore cannot exchange the two blocks, since it
preserves a nonempty monochromatic pair. Consequently an invariant
unordered partition would give an invariant individual five-subset,
which has already been excluded. D is contained in none of the
conjugates of H1,H3,W5, proving the local failure.

For completeness, write D=<d,h>, where h lifts residue Frobenius, and
label the two ramified pairs separately. On the six residual letters h
has the indicated even-cycle pattern of order m=2,4,6. On each ramified
pair it is the identity or a flip. Since D lies in A10, the parity of
those two flips is fixed by the residual sign; changing h by d switches
both flips. Thus D is conjugate to the explicit C2 x C_m models in the
checker. These have orders 4,8,12, respectively. This also verifies
directly the claimed group list.

The added H2 factor in the rank-301 construction supplies a local point
in every obstructed case, precisely through either ramified pair.

## Theorem 3. Density zero in every arithmetic progression

Fix integers a,M with M nonzero. Among the parameters s=a+Mn for which
P_s has full A10 splitting-field group, those satisfying

\[
                \operatorname{Spec}(A_s^{256})(\mathbb Q_v)
                  \ne\varnothing\quad\text{for every place }v
\]

have natural density zero as a set of integer indices n. Every such
progression nevertheless contains infinitely many full-A10 parameters
failing this local condition. Their splitting fields can be chosen
linearly disjoint.

Proof of density zero. For any p in E not dividing M, the condition
v_p(a+Mn)=1 occupies exactly p-1 of the p^2 residue classes of n modulo
p^2. By Theorem 2 these classes must be avoided. For any finite subset
E0 of these primes, the Chinese remainder theorem bounds the upper
density of everywhere locally soluble indices by

\[
                    \prod_{p\in E_0}
                         \left(1-\frac1p+\frac1{p^2}\right).
\tag{8}
\]

The positive prime density in Theorem 1 implies that the sum of 1/p
over E diverges. Therefore the products (8) tend to zero as E0 exhausts
E away from the finitely many divisors of M. This proves density zero.
Only the necessary zero-branch conditions were used; ramification at
the four other finite branch values may impose further restrictions.

For infinitely many failures, choose p in E outside the divisors of M.
Choose n0 modulo p^2 so that a+Mn0=p modulo p^2. Every parameter

\[
                    s=a+M(n_0+p^2t),\qquad t\in\mathbb Z,
\tag{9}
\]

has valuation one at p and hence the local failure. This affine
reparametrization preserves the regular generic A10 group. The retained
integral Hilbert theorem, Zywina
[Theorem 1.2, page 2](sources/zywina-quantitative-hit-2010.pdf), supplies
infinitely many full-group integer t. To obtain linear disjointness,
apply the same full arithmetic-group theorem to the compositum of the
regular generic cover with the finite Galois compositum of previously
selected fields, exactly as in the subgroup-cover note, Theorem 5.
Its generic group is the direct product, and t remains an integer
parameter in (9).

Thus finitely many congruences cannot make a whole progression work for
the rank-256 construction. This is compatible with the possible
existence of a thinner infinite admissible set. It neither proves nor
disproves that existence.

## Proposition 4. A rational base change cannot erase all five finite inertias

For any nonconstant rational map phi:P1 -> P1 in characteristic zero,
at most three target values can have every point in their fibre of
even ramification index. Consequently the pullback of this degree-ten
cover by any rational phi leaves nontrivial double-transposition inertia
above at least two of its five original finite branch values.

Proof. Put d=deg(phi), and suppose q distinct target values have only
even indices in their fibres. Each such fibre has at most d/2 points
and contributes at least d/2 to the ramification sum. Riemann--Hurwitz
therefore gives

\[
                         qd/2\leq2d-2.
\]

Hence q<4. In the local tame model of a double-transposition branch,
ramification index e in the base change replaces the order-two inertia
generator by its e-th power. It is killed precisely when e is even.
This proves the assertion. For odd d, no fibre can have all indices
even, so none of the five branch values can be erased completely.

The obvious quadratic base curve killing all finite order-two inertia
is w^2=sW(s), where W is the branch quartic from the pencil note. Its
squarefree degree-five right side gives genus two. Thus this construction
does not provide a rational parameter line. Other ways of controlling
decomposition groups, rather than killing all inertia, are not excluded.

## Scope and verification

The theorem concerns the precise subgroup product (1). It is not a
universal obstruction to all A10 Hasse constructions, nor to individual
rank-256 examples. The rank-301 infinite family remains valid; its H2
factor was designed to handle the groups that occur here. The fixed-map
compiler questions are separate from these local arithmetic conditions.

The checker verifies the two discriminants, the branch factorization,
the exact exceptional primes, both modular irreducibility certificates,
the sign-space generation, all 48 signed permutations and their 21
all-even elements, and all three local-group models. A finite
base-change-degree sanity table supplements the written
Riemann--Hurwitz proof. The infinite density and progression claims use
the written Hensel/CRT argument and the two named primary theorem
imports. These are internal exact controls, not external review or
formal verification.
