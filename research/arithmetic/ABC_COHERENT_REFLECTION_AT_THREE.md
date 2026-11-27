# The local condition at 3 for the reflected triangle

For the [reflected cubic of a coherent norm triangle](ABC_COHERENT_NORM_CUBIC_REFLECTION.md#the-exact-reflected-kummer-class),
the class sum of the three vertices specifies the field. Its image
in a small group of 3-adic units specifies the discriminant and
splitting at \(3\). This lets the prime-path count select fundamental
reflected discriminants while retaining its established number
of triangles per field.

## The residue test

**Proposition RT1.** Use the notation of RF1 with all three
vertex primes different from \(3\), and suppose the class sum is
nonzero. Write \(K=\mathbb Q(\sqrt{-d})\) and put
\[
  \gamma=\overline\alpha_{12}\alpha_{23}
     =\frac{t+s\sqrt{-d}}2,\qquad
  t=\frac{u_{12}u_{23}+d v_{12}v_{23}}2,\quad
  s=\frac{u_{12}v_{23}-u_{23}v_{12}}2.
\]
Let \(M\) be the reflected cubic field,
\(L=\mathbb Q(\sqrt{\operatorname{sf}(3d)})\), and let
\(\Delta_K,\Delta_L\) denote the fundamental discriminants.
The following conditions decide the local cubic algebra at \(3\).

| Condition | Test | Factorization at \(3\) |
| --- | --- | --- |
| \(3\nmid d\) | \(\gamma^8\equiv1\bmod9\mathcal O_K\) | One rational factor and a ramified quadratic factor |
| \(3\nmid d\) | Test fails | Totally ramified cubic |
| \(3\mid d,\ d/3\equiv2\bmod3\) | \(3\mid s\) | Unramified factors of degrees one and two |
| \(3\mid d,\ d/3\equiv1\bmod3\) | \(9\mid s\) | Completely split |
| \(3\mid d,\ d/3\equiv1\bmod3\) | \(3\mid s,\ 9\nmid s\) | Unramified cubic |
| \(3\mid d\) | \(3\nmid s\) | Totally ramified cubic |

The global discriminant is \(\Delta_L\) in the first, third,
fourth and fifth rows. In the two totally ramified cases it is
\(27|\Delta_K|\).

### The two local thresholds

Let \(T\) be a completion of \(K(\zeta_3)\) at \(3\) and put
\(\lambda=\zeta_3-1\). Its ramification index over \(\mathbb Q_3\)
is two: the field is a compositum of at most two quadratic
extensions, the tame inertia is cyclic of exponent at most two,
and \(T\) contains \(\zeta_3\). Thus \(\lambda\) is a uniformizer.
For a unit \(\eta\in T\), the Kummer extension given by \(\eta\)
is unramified, including the trivial extension, precisely when
\(\eta\) is a cube modulo \(\lambda^3\). It is trivial precisely
when \(\eta\) is a cube modulo \(\lambda^4\).

In fact for \(\eta=1+\lambda^3 c\), substitution of \(1+\lambda X\)
for a cube root leads to
\[
               X^3+(3/\lambda)X^2+(3/\lambda^2)X-c=0.
\]
The reduction is \(X^3-X-\bar c\) with simple roots, so all roots
lift over an unramified extension. Choosing \(\bar c\) outside the
image of \(X^3-X\) gives the unique unramified cubic extension.
Kummer theory then places every unramified class on the same
line, with all representatives on that line congruent to cubes
modulo \(\lambda^3\). If \(c\equiv0\bmod\lambda\), the reduced
polynomial has a root in the residue field, so the class is
trivial. Conversely an actual cube is a cube modulo \(\lambda^4\).
This proves both thresholds.

At a completion \(K_v\) below \(T\), the congruences descend by
the norm. In fact \(N_{T/K_v}(1+\lambda^j\mathcal O_T)\) lies
in \(1+\mathfrak p_v^{\lceil j/e(T/K_v)\rceil}\), and the
degree one or two is invertible on the cube quotient. The converse
is immediate by extending the congruence to \(T\).
Thus both thresholds become \(9\mathcal O_{K_v}\) when
\(K_v/\mathbb Q_3\) is unramified. If \(K_v\) is ramified they
become \(\mathfrak p_v^3\) and \(\mathfrak p_v^4=9\mathcal O_{K_v}\),
respectively. Ramification of index three in \(N/L\) survives
the quadratic base extension to \(T\).

### The computation in the quadratic field

The hypotheses on the vertices make \(\gamma\) a unit at all
places above \(3\). Its rational norm is a cube of a 3-adic
unit. If \(3\nmid d\), the reduction ring modulo \(3\) is
\(\mathbb F_9\) or \(\mathbb F_3\times\mathbb F_3\).
Each residue unit has a unique multiplicative lift of order
dividing eight modulo \(9\). The kernel of the reduction on
units modulo \(9\) has exponent three. Cubeness is therefore
equivalent to \(\gamma^8=1\) modulo \(9\), as claimed.

If \(3\mid d\), write \(r=d/3\) and \(\theta=\sqrt{-d}\).
The reduction modulo \(\mathfrak p_v^3\) has cube units
represented by \(+1,-1\). The rational part of
\(\gamma-\varepsilon\) has even valuation at \(v\), and the
\(\theta\)-part has odd valuation, so they cannot cancel.
Consequently \(\gamma\equiv\varepsilon\bmod\mathfrak p_v^3\)
requires \(3\mid s\) and \(t\equiv2\varepsilon\bmod9\).
The latter follows from the norm condition as soon as
\(3\mid s\): writing \((t^2+d s^2)/4=n^3\) gives
\(n\equiv1\bmod3\) and hence \(n^3\equiv1\bmod9\).
To test cubes modulo \(9\), let \(a\) be a rational 3-adic unit
and \(b\in\mathbb Z_3\). Then
\[
         (a+b\theta)^3\equiv
                a^3+3b(a^2-r b^2)\theta\pmod9 .
\]
The possible values of the coefficient of \(3\theta\) modulo
\(3\) are the values of \((1-r)b\). If \(r\equiv2\bmod3\),
every unramified class is therefore a cube locally after the
test \(3\mid s\). The resolvent \(L\) is inert at \(3\), so the
unramified Frobenius is a transposition. If \(r\equiv1\bmod3\),
the resolvent is split at \(3\), and the class is a local cube
exactly when \(9\mid s\). This yields the split and inert cubic
cases. For \(3\nmid d\) in the unramified row, \(L\) is ramified
and local inertia is a transposition as in RF3.

### The ramified discriminant

In either ramified case the inertia in \(N/\mathbb Q\) above
\(3\) contains a three-cycle, and its action on the three roots
is transitive. Hence \(M\) is totally ramified at \(3\).
If \(3\mid d\), the completions of \(L\) at \(3\) are
unramified over \(\mathbb Q_3\). Every unit congruent to \(1\)
modulo \(9\) there is a local cube, by the logarithm on principal
units. The cubic character of a ramified extension is nontrivial
on \(1+3\mathcal O_{L_v}\) since the residue units have order
prime to three. Its conductor exponent is thus two, giving
relative conductor \((9)\) at \(3\).

If \(3\nmid d\), choose a cube root of the residue of \(\gamma\)
in \(K_v\) and divide its cube out. Failure of the test modulo
\(9\) yields a unit \(\beta\) with \(v_T(\beta-1)=2\).
Write \(\xi=\sqrt[3]\beta\). The equation for
\(X=(\xi-1)/\lambda\) above, with constant term of valuation
\(-1\), has one Newton slope \(1/3\). Thus
\(\pi=\lambda/(\xi-1)\) is a uniformizer in the degree-three
extension of \(T\). For a nontrivial generator
\(\sigma(\xi)=\zeta_3\xi\), normalize valuation in that extension
and compute
\[
          v\bigl(\sigma(\pi)-\pi\bigr)
                       =v(\lambda^2)-2v(\xi-1)=6-4=2 .
\]
The ramification break is one, so the character conductor
exponent is two over \(T\). Here \(T/L_v\) is unramified or
trivial, so the conductor exponent over \(L_v\) is also two.
Since \(L_v\) itself is ramified quadratic over \(\mathbb Q_3\),
the relative conductor is \((3)\) at \(3\).

The relative cubic character has no other finite ramification
by RF1, and no infinite ramification because it has odd order.
The cubic discriminant formula is \(\Delta_L f^2\) for relative
conductor \(f\mathcal O_L\). Thus the ramified cases give
\((|\Delta_K|/3)9^2\) and \((3|\Delta_K|)3^2\), respectively.
This agrees with the virtual-unit discriminant alternatives in
[Cohen–Thorne, Section 3 and Proposition 4.1](https://arxiv.org/pdf/1301.3563).

## Large coherent parts retain the fundamental sector

**Corollary RT2.** In a coherent part with \(m\) distinct class
labels, an absolute positive fraction of its ordered triples
with distinct vertices and different coefficient colours on their
consecutive edges give nontrivial reflected cubic fields of
fundamental discriminant, for \(m\) sufficiently large.
If \(3\mid d\) and \(d/3\equiv1\bmod3\), one can additionally
require that \(3\) split completely. The class labels are
automatically distinct at the prime-path scales \(4\max p_i^2<d\).

In \(K\), the virtual-unit quotient is \(\mathrm{Cl}(K)[3]\)
as in RF1. Its elements have representatives that are units at
all places over \(3\): adjust the valuations by a cube using
weak approximation. Such changes keep the induced local classes
well-defined. If \(3\nmid d\), the residue map
\[
                \eta\longmapsto (\eta^8-1)/3\bmod3
\]
is additive on these unit classes and has trace zero, since
\(N(\eta)\) is a cube. The trace-zero subspace has dimension one
over \(\mathbb F_3\) in both reduction rings. If \(3\mid d\),
write \(\eta=(t+s\theta)/2\); the map \(s/t\bmod3\) is additive
and kills cubes. On its kernel, when \(d/3\equiv1\bmod3\),
the map \(s/(3t)\bmod3\) is again additive and kills cubes by
the explicit cubing formula. Thus the subspace of classes
with fundamental reflected discriminant has codimension at most
one. The subspace with the extra complete-splitting condition
in the stated case has codimension at most two. These subspaces
include the zero class.

Group the \(m\) vertices by their images in the appropriate
quotient of dimension at most two. A triple from one fibre
has class sum in the subspace, since the three images agree.
In each sufficiently large fibre the original monochromatic
\(K_4\) obstruction gives a positive fraction of paths with
different colours on their two adjacent edges. Discarding the
at most nine fibres of bounded size and applying Hölder gives
\(\gg m^3\) such ordered paths. The triples with zero class sum
number only \(O(m^2)\), because once two labels are specified
there is at most one third vertex with the required label. This proves
the fraction assertion. For the last claim, equality of the
classes of different primes would give an integral principal
ideal of norm \(p_i p_j<d/4\) by multiplying their oriented
prime ideals with one conjugated. Such an ideal cannot have
a rational generator by its one-sided supports, while the
norm of a nonrational integral generator is at least \(d/4\).

The triangle count in the
[prime-path tail proof](ABC_QUADRATIC_THREE_TORSION_FERMAT_TAIL.md#witness-selection)
can therefore use fundamental reflected cubic fields at the
same scale after discarding parts below a sufficiently large
fixed size. This restriction still leaves the global problem
of summing over coefficients and marked prime lines belonging
to the same cubic field.

## Exact controls

Here \((e,f)\) are the ramification index and residue degree
of the factors above \(3\). The traces and \(E\) are recovered
from the listed primes and normalized coefficients using the
three-coefficient reconstruction. The last two columns can be
replayed in PARI/GP from the binary form (1) of RF1 using
`polredabs`, `nfinit` and `idealprimedec` at \(3\).

| Primes | \((x,y,z)\) | \(d\) | \(s\) | Reflected discriminant | \((e,f)\) at \(3\) |
| --- | --- | ---: | ---: | ---: | --- |
| \(7,11,17\) | \(1,2,1\) | \(605107\) | \(-3541\) | \(1815321\) | \((1,1),(2,1)\) |
| \(7,11,19\) | \(1,4,1\) | \(114127\) | \(-5484\) | \(3081429\) | \((3,1)\) |
| \(7,13,19\) | \(1,3,2\) | \(578847\) | \(-3772\) | \(15628869\) | \((3,1)\) |
| \(7,13,31\) | \(1,1,4\) | \(245202\) | \(-14724\) | \(326936\) | \((1,1),(1,2)\) |
| \(7,19,31\) | \(1,5,2\) | \(563979\) | \(-42432\) | \(187993\) | \((1,3)\) |
| \(7,37,67\) | \(2,7,5\) | \(1785027\) | \(-711747\) | \(595009\) | \((1,1),(1,1),(1,1)\) |
