# A quartic Pell curve with two independent class directions

The replacement problem in
[LY9](LEVIN_YAN_WILJANEN_DEGREE_TEN_MIRROR_SLICE.md)
asks for a Pell polynomial with several norm-cube sections and one
congruence slice that preserves their independence and supplies an
integral unit noncubic at three. The construction below does this for
two sections. The exponent and rank lower bound concern this particular
Pell family. For unconditional quadratic \(3\)-rank counts without this
Pell condition, see [Levin–Yan–Wiljanen, Theorem 1 and Corollary 10](sources/levin-yan-wiljanen-large-three-rank.pdf).
Their Theorem 3 has degree divisible by six, and their Theorem 2 assumes
a rational Weierstrass point. The even octic constructed here has
neither property. The argument uses direct arithmetic on the sections.

## Two sections on a common integral slice

**Theorem CP1.** There is an integral quartic polynomial
\(\mathcal P\) for which the smooth genus-three curve
\(R^2=\mathcal P^2+1\) has rational Jacobian \(3\)-torsion \(C_3^2\).
Its conjugate infinity difference has order four. For every sufficiently
large positive integral specialization, the real quadratic field
\(K=\mathbb Q(\sqrt{\mathcal P^2+1})\) is ramified at two and split at
three. The norm-\(-1\) integral unit \(U=\mathcal P+R\) is noncubic at
both places over three. Two explicit virtual cube classes are locally
cubic at both places over three and give independent classes in
\(\operatorname{Cl}(K)[3]\). The specialization produces
\(\gg X^{1/8}/\log X\) distinct fields of discriminant at most \(X\).

Start with rational polynomials in \(t\). For \(w\ne0\) put
\[
 v=\frac1{24(1+2w^2)},\qquad z=wv,\qquad s=t-\frac3{8v},
 \qquad Z=z(s^2+3),\qquad H=\frac{s+1}{2}+v(s^2+3).
\]
Define \(P\) by
\[
 16zP=s+3+6v(s+1)^2+12v^2(s+1)(s^2+3)
                  +8v^3(s^2+3)^2-8z^2(s^2+3).
\]
The factorization \(((s+1)/2)^3+1=(s+3)(s^2+3)/8\)
gives \(2PZ+Z^2-1=H^3\). The shift in \(s\) removes the cubic
coefficient of \(P\). At the new origin the derivative of the displayed
right-hand side is \((24v^2+48z^2-v)/(8v)=0\).
Thus \(P(-t)=P(t)\), and the reflected pair supplies a second
section over the same quartic. We take \(w=1/48\). Explicitly,
\[
 \begin{split}
 P(t)&=\frac{55296}{1153^2}t^4-\frac{2670347}{4\cdot1153^2}t^2
                    +\frac{1037316177403}{2^{17}\cdot1153^2},\\
 Z(t)&=\frac{t^2}{1153}-\frac t{64}+\frac{1378561}{2^{14}\cdot1153},\\
 H(t)&=\frac{48t^2}{1153}-\frac t4-\frac{591617}{2^{10}\cdot1153}.
 \end{split}
\]
Use \(Z_\pm(t)=Z(\pm t)\), \(H_\pm(t)=H(\pm t)\),
\(F=P^2+1\), and \(\alpha_\pm=1-Z_\pm U\) where \(U=P+R\).
Then \(N U=-1\) and \(N\alpha_\pm=-H_\pm^3\).
The polynomials \(F,P\) have degrees eight and four, and \(F\) is
squarefree.

Take \(t=N/128\) with \(N\) in the following common congruence class.

| Modulus | Residue |
| ---: | ---: |
| \(2^{12}\) | 1105 |
| \(3^5\) | 1 |
| 5 | 2 |
| 13 | 6 |
| 17 | 1 |
| 19 | 2 |
| 29 | 2 |
| 31 | 1 |
| \(1153^2\) | 657 |
| 2677 | 1 |

Write \(N=N_0+L j\) for any representative \(N_0\) of these
congruences and their product \(L\). The coefficients of
\(\mathcal P(j)=P((N_0+L j)/128)\) are integers. In particular,
the five coefficients of \(P((1105+2^{12}j)/128)\) have
dyadic valuations \(0,1,9,21,31\) in increasing order.
At \(1153\) the three polynomials \(P,Z_\pm,H_\pm\) have integral
values on the selected residue class, and both \(H_\pm\) are units.
For \(P\), substitution \(N=657+1153^2j\) also gives integral
coefficients there. All other odd primes are away from the
denominators. The radicand for large \(j\) is strictly between
consecutive squares, and \(U\) and its inverse are integral.
At two, \(\mathcal P\) is odd, \(F\equiv2\bmod 8\), and
\(K_2\) is ramified. There is only one prime over two, so the
cube norm makes the valuation of each \(\alpha_\pm\) a multiple of
three there even when \(Z_\pm\) has a denominator.

At an odd split prime with integral \(Z_\pm\), both values of an
\(\alpha_\pm\) can vanish modulo the prime only if
\(Z_\pm^2+1\equiv H_\pm\equiv0\). Indeed the product of the two
reductions of \(U\) is \(-1\), and the product of the two reductions
of \(\alpha_\pm\) is \(-H_\pm^3\). Away from the coefficient
denominators the set of possible common roots is detected by
\[
 \operatorname{Res}_t(H,Z^2+1)
       =\frac{5\cdot13\cdot17\cdot29\cdot2677}{2^4\cdot1153^4}.
\]
The reflected pair gives the same set of primes. At \(5,17,29,2677\)
both \(H_\pm\) are units on the chosen slice. At \(13\), the plus
section has \(H_+\equiv0\) and \(Z_+^2+1\equiv11\); the minus
section has \(H_-\equiv5\). At \(1153\) use the integral unit values
of \(H_\pm\) given above. At a nonsplit prime the cube norm
already forces virtuality since the residue degree divides two.
At the remaining split primes at most one of the two valuations
is nonzero, so the cube norm again suffices.

At three the branch with \(R\equiv1\bmod3\) satisfies
\[
 P\equiv87,\quad R\equiv100,\quad
 \alpha_+\equiv223,\quad \alpha_-\equiv27\pmod{243}.
\]
Here \(U\equiv7\bmod9\) and \(\alpha_-/27\equiv1\bmod9\).
A unit in \(\mathbb Q_3\) is a cube exactly when it is congruent
to \(1\) or \(-1\) modulo nine. Thus
\(\beta_+=\alpha_+U^2\) and \(\beta_-=\alpha_-\) are both
locally cubic, and \(U\) is noncubic. The other completion has the
same assertions because the rational norms of all three elements
are cubes. More generally, for a split prime at three the two local
classes of a virtual element with cube norm are opposite unit
classes. A single global unit noncubic there supplies a simultaneous
adjustment at both places.

For independence, at each of the following split primes choose the
displayed root \(R\) and primitive cube root \(\omega\) in the residue
field. The last three columns specify \(\log_\omega(a^{(q-1)/3})\)
for \(a=U,\alpha_+,\alpha_-\). The three elements are units at
the indicated place.

| \(q\) | \(R\bmod q\) | \(\omega\) | \(U\) | \(\alpha_+\) | \(\alpha_-\) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 13 | 11 | 3 | 2 | 1 | 0 |
| 19 | 5 | 7 | 2 | 2 | 0 |
| 31 | 6 | 25 | 0 | 1 | 1 |

The matrix has determinant \(2\bmod3\). Hence \(U,\beta_+,\beta_-\)
are independent modulo cubes in every field of this slice.
Since the real quadratic unit group modulo cubes is a line,
the images of \(\beta_+,\beta_-\) in the ordinary class group are
independent. Equivalently, the local cubic condition separates them
from the noncubic unit line. A single bad prime for a section can
therefore be useful for the independence test. At 13 one branch
of the plus section is nonunit, and the chosen branch still
detects its cubic class.

For the Jacobian statement, each function \(\alpha_\pm\) has divisor
divisible by three. The two finite numerator factors \(R-P-Z_\pm\)
and \(R+P+Z_\pm\) have product \(-H_\pm^3\). A common zero would force
a multiple root of \(F\) by differentiating \(F=(P+Z_\pm)^2-H_\pm^3\).
The denominator \(R-P\) is a unit on the affine curve, and the two
infinity orders of \(\alpha_\pm\) are \(-6,0\).
If the corresponding divided divisors were dependent, a nontrivial
product of the two functions would be a cube times a rational
constant. Taking norms forces this constant to be a rational cube.
Specialization would contradict the residue test.
At the good prime 11, the respective point counts in degrees one
through three are \(12,140,1332\). Newton's identities and the
reciprocal functional equation give
\[
 X^6+9X^4+99X^2+1331.
\]
The Jacobian order \(1440\) is divisible by exactly \(3^2\), giving
the upper bound by prime-to-11 torsion specialization. The function
\(U\) gives the infinity difference of exact order four by the
constant-norm degree argument of LY9.

As \(N\) ranges up to \(T\) in the fixed progression there are
\(\gg T\) parameters and \(D_K\ll T^8\). The positive units \(U\)
are at most \(O(T^4)\). In any fixed real quadratic field the
ratio of successive positive norm-\(-1\) units is at least
\((3+\sqrt5)/2\), since it is a norm-one unit greater than one
with integral trace at least three. Each unit determines \(P\)
via \(2P=U-U^{-1}\), with at most four preimages. This bounds
the multiplicity by \(O(\log T)\) uniformly in the field and proves
the discriminant count.

## Exact mixed mirror discriminants

**Corollary CP2.** For every field in CP1 let \(D=\operatorname{Disc}(K)\)
and \(M=\mathbb Q(\sqrt{-3D})\), of discriminant \(-3D\) since
three splits in \(K\). The three cube classes give one
extension over \(M\) with group \(C_3^3\), unramified away from
three. Exactly one dimension of its local quotient over three
is nonzero, and the corresponding local order-three characters
have conductor exponent two over \(M_3\). Its maximal quotient
split at three has group \(C_3^2\).
For each \(k=1,2,3\) and each \(k\)-dimensional subspace \(B\)
of the three classes, the associated degree-\(3^k\) scalar
affine field has signature \((1,(3^k-1)/2)\) and
\[
 |\operatorname{Disc}(L_B)|
 = (3D)^{(3^k-1)/2}
   \begin{cases}
      1 & B\subseteq\langle \beta_+,\beta_-\rangle,\\
      3^{2\cdot3^{k-1}} & \text{otherwise}.
   \end{cases}
\]
At \(k=1\) the respective multiplicities are four and nine, and at
\(k=2\) they are one and twelve. The unique \(k=3\) subspace has
the larger discriminant.

Over \(E=K(\zeta_3)\), of degree two over \(K\), the three classes
stay independent and yield the Kummer extension of degree \(27\).
Conjugation on \(K\) and conjugation on \(\zeta_3\) each act by
\(-1\) on its translation group, the first by the cube norms
and the second by its action on the Kummer pairing. Their product
has a central order-two lift, whose fixed field is the asserted
extension of \(M\). Every subspace \(B\) gives the scalar affine
group \(C_3^k\rtimes C_2\) over \(\mathbb Q\) and the degree-\(3^k\)
field fixed by a reflection. The reflection has one fixed point
and \((3^k-1)/2\) pairs. This gives the signature at infinity.
Away from three, Kummer theory for virtual units shows
unramifiedness over \(E\), and hence over \(M\) by the coprime
degree of \(E/M\). At each such prime only the quadratic
inertia contributes, with discriminant exponent that of \(M\)
times \((3^k-1)/2\), including wild quadratic inertia at two.

At three, \(E/M\) splits because \(K_3\) splits.
The \(\beta_\pm\) are cubes and the only potentially nonzero
local translation line comes from \(U\). Its cubic field
over \(\mathbb Q_3\) is generated by a root of
\(X^3+3X-2P\), the sum of cube roots chosen with product \(-1\).
This polynomial is Eisenstein, and its discriminant has
valuation three because \(P^2+1\) is a unit. A root is a
uniformizer, so this is the field discriminant exponent.
Its normal closure has local group \(C_3\rtimes C_2\) and
discriminant exponent seven. Indeed its regular representation
is the sum of the trivial and quadratic characters and two
copies of the standard representation. Their conductor exponents
are \(0,1,3\), the last read from the cubic permutation
representation. Alternatively the quadratic base has discriminant
valuation one, so the relative cyclic order-three characters
have conductor exponent \((7-3)/2=2\) by the
conductor-discriminant formula.

For a subspace with only cubic local classes, the quadratic
involution gives exponent \((3^k-1)/2\) at three. Otherwise
the local decomposition group translates one coordinate and
negates all coordinates. On the quotient by the translation
line, there is one zero and \((3^{k-1}-1)/2\) opposite pairs.
These give one orbit of length three and the remaining orbits
of length six. The exponent at three is therefore
\(3+7(3^{k-1}-1)/2\), larger by \(2\cdot3^{k-1}\).
Different subspaces have different normal closures in the
common scalar extension. The multiplicities follow by counting
the lines and planes of a three-dimensional space over
\(\mathbb F_3\).

## Factorization of the whole reflected family

**Theorem CP3.** For every \(w\in\mathbb Q^\times\), the quartic
\(P_w\) defined in CP1 has \(P_w^2+1\) irreducible over \(\mathbb Q\).
In particular, changing the rational symmetry parameter of these two
norm sections cannot replace the octic by smaller rational factors.

Put \(P_w=a t^4+b t^2+c\) and take \(i^2=-1\). The leading coefficient
is \(a=v^2/(2w)\ne0\). For \(P_w-i\) use the monic constant and
discriminant in the quadratic variable \(t^2\)
\[
 C=(c-i)/a,\qquad D=(b^2-4a(c-i))/a^2.
\]
The defining formula gives
\[
\begin{split}
 C&=-19440(w^2+iw-1/5)(w^2-iw/3+1/3)^3,\\
 D&=147456(w-3i/4)(w+i/4)^3(w^2+1/2)^2.
\end{split}
\]
Their rational norm squareclasses are respectively
\[
  (25w^4+15w^2+1)(9w^4+7w^2+1),\qquad
  (16w^2+1)(16w^2+9).
\]
Both are nonsquares. Here are the two descents that enforce this.

For the second norm, put \(X=4w\). If
\((X^2+1)(X^2+9)=Y^2\), then \((X^2,XY)\) lies on
\(E:y^2=x(x^2+10x+9)\). Under the two-isogeny rank formula,
the \(x\)-squareclass images on \(E\) and its isogenous curve
\(E':y^2=x(x-4)(x-16)\) have product of sizes \(2^{2+r}\).
The primitive quartic criterion for a class \(d\) on
\(y^2=x(x^2+A x+B)\) is
\(n^2=d u^4+A u^2v^2+(B/d)v^4\) for coprime integers \(u,v\).
The possible classes on \(E\) are \(1,-1,3,-3\). On \(E'\)
only \(1,2\) are possible because \(x\ge0\). The equation for
\(2\) is \(n^2=2u^4-20u^2v^2+32v^4\). If \(u\) is odd, the
right side is \(2\) or \(6\bmod8\). Otherwise \(v\) is odd. For
\(4\mid u\) its dyadic valuation is five, and for \(u=2u_0\) with
\(u_0\) odd its quotient by \(16\) is \(7\bmod8\). Thus \(E\) has
rank zero. Good reduction at five and seven bounds the torsion
order by eight since each count is eight. The eight displayed
coordinates are the point at infinity, \(x=0,-1,-9\) with
\(y=0\), and \((3,\pm12),(-3,\pm6)\). The only positive \(x\)
is \(3\), which cannot equal \(X^2\).

For the first norm write \(w=A/B\) in lowest terms. The integers
\[
 I=25A^4+15A^2B^2+B^4,\qquad J=9A^4+7A^2B^2+B^4
\]
are positive, odd and coprime. A common odd prime would divide
\(I-J=8A^2(2A^2+B^2)\). The case of a prime dividing \(A\) is
excluded by coprimality of \(A,B\). The other case leaves
\(I\equiv-A^4\) at that prime. If the norm were square, \(I\)
would be square. Writing \(Y^2=I/B^4\) gives a point
\((25w^2,25wY)\) on \(E_1:y^2=x(x^2+15x+25)\).
The isogenous curve is \(y^2=x(x-5)(x-25)\), whose squareclass
image has size two with possible classes \(1,5\).
On \(E_1\) the candidates are \(1,-1,5,-5\). The class \(5\)
is represented by \((5,25)\). The class \(-1\) is impossible
because \(-u^4+15u^2v^2-25v^4\) is \(3,5\) or \(7\bmod8\)
for coprime \(u,v\). This also excludes \(-5\), so \(E_1\) has
rank zero. The reductions at three and seven both have order four.
The nonzero positive \(x\) of a rational point is therefore \(5\),
from the four known points \(O,(0,0),(5,\pm25)\). It cannot
equal \(25w^2\). The formula and squareclass criterion in these
descents are those in [Dixit–Murty–Pathak, Section 2](sources/dixit-murty-pathak-mordell-weil-rank.pdf).

Consequently \(C,D\) are both nonsquares in \(\mathbb Q(i)\).
The nonsquare \(D\) excludes a root of \(P_w-i\) and a
factorization into even quadratics. In every factorization into
monic quadratics, evenness gives
\((t^2+\lambda t+\mu)(t^2-\lambda t+\nu)\) and
\(\lambda(\nu-\mu)=0\). If \(\lambda\ne0\), then \(C=\mu^2\).
Thus \(P_w-i\) is irreducible over \(\mathbb Q(i)\). Every root
of \(P_w^2+1\) generates a field containing \(i\), so its
rational degree is eight.

## The fixed slice and its squarefree sieve

**Corollary CP4.** At \(w=1/48\), the degree-eight splitting
group is the wreath product \(D_4\wr C_2\) in its action on two
blocks of four, where \(D_4\) has order eight. Among unramified
primes where the model retains degree, the numbers \(0,2,4,6,8\)
of roots have respective natural densities
\[
              (89,20,14,4,1)/128.
\]
On refining CP1 to \(N\equiv2\bmod25\), the integral octic
\(f(j)=P((N_0'+L'j)/128)^2+1\) has no fixed square divisor.
Its squarefree sieve has a positive product of local densities.

For the group calculation, use the classes \(C,D,\bar C,\bar D\)
in \(\mathbb Q(i)^\times/\mathbb Q(i)^{\times2}\). At the primes
over \(101\) with \(i\equiv\pm10\) the only odd valuations occur
for \(C\) and \(\bar C\), on opposite branches. At the primes over
\(1297\) with \(i\equiv\pm36\), the only odd valuations occur for
\(D\) and \(\bar D\), on opposite branches. Thus the four classes
are independent. Each even quartic splitting field is contained
in a \(D_4\) action on its two opposite pairs and contains the
quadratic roots of its \(C,D\). The independent two quadratic
classes force its group to be \(D_4\). Any nontrivial intersection
of the two splitting fields would give a common quadratic
subfield, since a nontrivial quotient of a finite \(2\)-group
has such a quotient. Each \(D_4\) has three quadratic subfields
and they are precisely the ones from the two independent
classes and their product. Hence the two splitting fields are
disjoint over \(\mathbb Q(i)\). Conjugation exchanges the blocks,
giving the full order-\(128\) wreath product over \(\mathbb Q\).
In a block preserved by Frobenius, the numbers \(0,2,4\) of
fixed points have multiplicities \(5,2,1\) in \(D_4\). Block
exchange precludes fixed points. Convolving the two block
counts and adjoining the exchange coset gives the densities
by [Chebotarev](sources/milne-cft403.pdf).

On the refined slice \(f\equiv2\bmod8\), \(f\equiv15\bmod25\),
and \(f\) is a unit at all other primes dividing \(L'\). At
\(p\nmid L'\), the leading coefficient is a unit. If
\(p\equiv3\bmod4\), \(f\) has no roots modulo \(p\). Every
remaining \(p\nmid L'\) is at least \(13\), so a degree-eight
polynomial cannot vanish identically on \(\mathbb F_p\).
Writing \(\rho(p^2)=\#\{j\bmod p^2:f(j)\equiv0\bmod p^2\}\),
we get \(\rho(p^2)<p^2\) everywhere, and \(\rho(p^2)\le8\)
at all but finitely many primes by simple-root lifting.
This proves positivity of \(\delta=\prod_p(1-\rho(p^2)/p^2)\).

**Conditional consequence.** Assume the \(abc\) conjecture.
[Granville's Theorem 1](sources/granville-count-squarefrees-1998-excerpt.pdf)
applies to \(f\), which is separable, integral and has fixed
divisor squarefree. The number of squarefree values for
\(1\le j\le T\) is asymptotic to \(\delta T\).
If \(A_f\) is the leading coefficient, the number of
distinct quadratic fields from these squarefree values of
discriminant at most \(X\) is asymptotic to
\(\delta(X/(4A_f))^{1/8}\). Their discriminant is \(D=4f(j)\).
The corresponding Mordell models of
[PM2](PELL_NORM_SECTIONS_MORDELL_RANK.md) have conductor
\(144 f(j)^2=9D^2\), since \(v_2(f)=1\) and \(v_3(f)=0\).
Thus the same subfamily gives asymptotically
\(\delta(X/(144 A_f^2))^{1/16}\) distinct conductors up to \(X\).
Both counts use eventual strict monotonicity of \(f\).
The unconditional argument gives all fixed finite sieves.
For squarefree values themselves it leaves the large-prime-square
tail of Granville's sieve.

## Scope and verification

The [exact checker](verify_quartic_pell_two_class_directions.py)
verifies the polynomial identities, denominators and local data,
the residue-class independence and finite-field Frobenius calculation,
the symbolic squareclasses, two local descents and the fixed octic sieve data.
Run it with Sage's Python. The prime-to-three ramification and the
class-group and field-count conclusions use the written arguments.
The five-direction replacement in LY9 remains open. In particular
the exact \(3\)-torsion calculation prevents obtaining a third
Jacobian direction by looking for more sections on this fixed
quartic.
