# Prime balanced Image failures are automatically absolutely nonproper powers

Let K be a field of characteristic zero, n>=1, and

    S_n=K[xi_1,...,xi_n,z_1,...,z_n].

Use the polynomial-valued contraction

    E_n(xi^alpha z^beta)=partial_z^alpha z^beta,
    M_n=sum_i(partial_(z_i)-xi_i)S_n.

A polynomial F is pure if E_n(F^m)=0 for EVERY positive integer m.
It is safe if for EVERY fixed polynomial marker Q in S_n there exists
a cutoff, allowed to depend on Q, beyond which E_n(QF^m)=0. An Image
counterexample is pure and unsafe: one fixed Q fails at arbitrarily
large orders. No finite nonzero prefix alone asserts unsafety.

Balanced degree D means homogeneous of degree D in the xi variables
and separately of degree D in the z variables. All assertions below
remain valid after EVERY extension of the coefficient field.

**Theorem LinearSampledVanishingNilpotenceAndSafety.** For
H=xi^T Bz, B in Mat_n(K), and each FIXED positive integer r, the
following are equivalent:

1. B is nilpotent.
2. E_n(H^k)=0 for every positive integer k.
3. E_n(H^(rm))=0 for every positive integer m.
4. E_n(H^(rm))=0 for all sufficiently large integers m.
5. H is safe for every original polynomial marker.
6. H^r is safe for every original polynomial marker.

In particular a pure balanced-degree-one phase is always safe.

**Theorem PrimeBalancedCounterexampleAutomaticNonpower.** Let p be
a prime integer. Every balanced-degree-p original Image counterexample
is absolutely scalar-nonproper-power in the FULL polynomial ring.
Explicitly, for every field extension L/K there are no alpha in L*,
integer r>=2 and G in L[xi,z] such that F=alpha G^r. Radial primitivity
is NOT a hypothesis: the conclusion applies to primitive and nonprimitive
counterexamples alike.

These are conditional statements about actual counterexamples. Neither
theorem supplies a starting counterexample in an unproved degree. They
do not assert irreducibility, squarefreeness, or a two-pair degree5/7
existence or exclusion theorem.

## 1. The full Image kernel and extension boundary

For every monomial a, direct differentiation gives
E_n((partial_(z_i)-xi_i)a)=0. Conversely, if alpha_i>0, the literal identity

    xi^alpha z^beta
      =partial_(z_i)(xi^(alpha-e_i)z^beta)
          -(partial_(z_i)-xi_i)(xi^(alpha-e_i)z^beta)

reduces the number of xi factors by one modulo M_n. Repeating this
linearly gives P-E_n(P) in M_n. Since E_n is the identity on K[z],

    M_n=ker E_n.                                      (P1)

This is a VECTOR-space reduction; M_n is not being treated as an ideal.
Thus the definitions in the introduction are exactly the original
Image membership and Mathieu-failure conditions. The proof of(P1)
works over every extension L/K. A nonzero polynomial E_n(QF^m) over K
remains nonzero over L, so a given counterexample and its fixed marker
remain a counterexample after every coefficient-field extension.

## 2. MacMahon's determinant identity, proved here directly

For H=xi^T Bz put ell_i(z)=sum_j B_ij z_j. The multinomial expansion gives

    H^k/k!=sum_(|alpha|=k) xi^alpha ell(z)^alpha/alpha!.

Each ell(z)^alpha has total degree k. After contraction, the corresponding
scalar is [z^alpha]ell(z)^alpha. Therefore

    E_n(H^k)/k!
      =sum_(|alpha|=k)[z^alpha](Bz)^alpha
      =tr(Sym^k B).                                  (P2)

The last equality is the trace of the substitution operator
Q(z)->Q(Bz) on degree-k homogeneous polynomials in the monomial basis.
Over an algebraic closure, triangularize B. The induced operator has
diagonal entries lambda_1^alpha_1 ... lambda_n^alpha_n, |alpha|=k.
Trace is basis independent, so summing its generating series gives

    sum_(k>=0) E_n(H^k)t^k/k!
       =product_i(1-lambda_i t)^(-1)
       =1/det(I-tB).                                 (P3)

Both sides are over K, so the identity descends. This is a formal
power-series identity; convergence, a positive measure, and diagonalizability
of B are unnecessary. It is a scalar specialization of the classical
MacMahon master theorem, not a new determinant formula claim.

## 3. Multiples of one order detect nilpotence

Fix r>=1, let D(t)=det(I-tB), and pass to a constant extension containing
a primitive rth root of unity zeta. Averaging(P3) gives the rational identity

    (1/r)sum_(j=0)^(r-1) 1/D(zeta^j t)
       =sum_(m>=0) E_n(H^(rm))t^(rm)/(rm)!.           (P4)

This follows coefficientwise from the finite geometric sum of roots of
unity. If the sampled contractions vanish at every positive order, the
right side is1. If they vanish only eventually, it is a polynomial with
constant coefficient1.

If D is nonconstant, each rational function 1/D(zeta^j t) is regular
and vanishes at infinity. Their average does too. A nonzero polynomial
cannot vanish at infinity, giving a contradiction in either case.
This is just the valuation at infinity of rational functions, not an
analytic limiting assertion. Thus D=1 because D(0)=1.

The characteristic polynomial of B is s^n D(1/s)=s^n; Cayley--Hamilton
gives B^n=0. Conversely a nilpotent B has D=1 and(P3) gives every positive
contraction zero. This proves the equivalence of1--4.

The multiples condition is important. For B=diag(1,-1), the determinant
is1-t^2 and all ODD contractions vanish, while B is not nilpotent.
One must not replace the multiples rm by an arbitrary shifted progression.

## 4. Full polynomial-valued equivariance

For an invertible S in GL_n(K), define

    T_S Q(xi,z)=Q(S^(-T)xi,Sz).

The chain rule, applied to normally ordered monomials, gives the FULL
polynomial identity

    E_n(T_S Q)(z)=E_n(Q)(Sz).                          (P5)

For instance this also follows by taking coefficients in the formal
identity E_n(exp(xi^T u+z^T v))=exp(u^T v+z^T v).
Under T_S the matrix of H becomes S^(-1)BS. Thus vanishing or eventual
vanishing of original contractions is preserved by this conjugation,
including for arbitrary unbalanced polynomial markers.

E_n is NOT scalar-invariant on arbitrary markers. It is only a scalar
on balanced inputs. Likewise it does not kill every nonzero weight;
for example E_2(z_1 xi_1 z_2)=z_2 is nonzero. The next argument keeps
track of the bounded output degree rather than using either false assertion.

## 5. Nilpotent linear phases have safe powers for all original markers

A nilpotent matrix is similar over K itself to a strictly upper triangular
matrix (a nilpotent Jordan-chain basis suffices). By(P5), assume B is
strictly upper triangular. Assign integer weights

    wt(z_i)=i,   wt(xi_i)=-i.

Each monomial xi_i z_j of H then has j>i, so its weight is at least1.
The contraction E_n preserves weights on every surviving monomial.

First take a marker Q homogeneous of bidegree(a,b). If ell=b-a<0,
E_n(QH^k)=0 for every k by total derivative order. If ell>=0, every
nonzero output monomial would have total z degree ell and therefore
weight at most n ell. But every input monomial has weight at least
k+w_min(Q), where w_min(Q) is the least weight of a nonzero monomial
of Q. Consequently

    E_n(QH^k)=0 whenever k>n ell-w_min(Q).             (P6)

For a general fixed marker, apply this to its finitely many bihomogeneous
components and take the largest cutoff. Thus H is safe. Its power H^r
is also safe by restricting to k=rm. Conversely safety of H or H^r,
using marker1, gives eventual sampled vanishing and hence nilpotence
by Section3. This proves the equivalence of1--6 completely.

The bound is concrete and can be sharp: with H=xi_1 z_2 and Q=z_1^M,

    E_2(QH^k)=M!/(M-k)! z_1^(M-k)z_2^k  for0<=k<=M,

and it is zero for k>M. In particular the boundary k=M need not vanish.

## 6. Both grading extremes force a linear root in prime degree

Suppose F is nonzero of bidegree(p,p), p prime, and after an arbitrary
extension L/K one has F=alpha G^r, alpha!=0,r>=2. Decompose G by total
xi degree. Its nonzero lowest and highest homogeneous components have
degrees a and b. The corresponding extreme components of G^r have
degrees ra and rb and cannot vanish, since the polynomial ring is a
domain. Homogeneity of F forces a=b=p/r. Thus r divides p, r=p and
G is homogeneous of xi degree1.

Independently repeat the argument with total z degree. It forces G to be
homogeneous of z degree1 as well. Hence

    G=xi^T Bz, B in Mat_n(L).                          (P7)

Both gradings are essential; this is a full polynomial-ring argument,
not a conclusion merely in the quotient R=1. No radial primitivity is
used at any point.

## 7. Prime-degree counterexamples cannot be proper powers

Purity of F=alpha G^p gives

    E_n(G^(pm))=alpha^(-m)E_n(F^m)=0, m>=1.

Sections2--3 make B nilpotent. Section5 then makes G^p, and hence F,
safe over L for EVERY original polynomial marker. This contradicts the
unchanged fixed original counterexample marker after extension, as explained
in Section1. Thus no such scalar power factorization exists over any L.
This proves absolute scalar nonproper-power, for primitive and nonprimitive
F alike.

Equivalently, for each fixed n and prime p, existence of a balanced
degree-p Image counterexample is equivalent to existence of one that is
absolutely scalar-nonproper-power. The same equivalence holds with radial
primitivity imposed on BOTH sides. The reverse direction merely forgets
the nonpower property. These are universal equivalence PROPOSITIONS;
they do not unconditionally assert existence in any degree.

## 8. Necessary scope controls

Purity alone does not force nonpower, even with a primitive lift and prime
degree: for n>=2, (xi_1 z_2)^p is pure, radial-primitive and a pth power,
but is safe.
The UNSAFETY assumption in the main theorem is indispensable.

Primality also matters. In three pairs write pairs as(tau,t),(w,z),(v,y)
and let F_2=tau(t-y)(wz+vt), Q=tau y. Expanding the two binomials gives

    E_n(F_2^m)=0,
    E_n(QF_2^m)=(-1)^(m-1)(m+1)!m!, m>=1.

Indeed pure matching indices are j=s, giving the full alternating
binomial sum; the marker changes matching to j=s+1, giving the alternating
sum with its last term removed. Hence F_2^2 is an Image counterexample
of composite balanced degree4 and a genuine proper square, with the SAME
Q surviving at every order. This is a control of the prime restriction,
not a new seed construction claim. The seed and its binomial proof are
already retained in the repository's all-balanced-degree construction.

Finally characteristic zero is essential for the displayed sampled test:
in characteristic2 the nonnilpotent scalar B=[1] has E_1((xi z)^(2m))
equal to the zero scalar (2m)! for every m>=1. No positive-characteristic
extension is asserted.

## 9. Classical and repository antecedents

The determinant identity is a specialization of MacMahon's master theorem.
An actual primary proof inspected here is Doron Zeilberger, *A combinatorial
approach to matrix algebra*, Discrete Mathematics56(1985),61--72,
Section2, printed63 equation(2):
https://sites.math.rutgers.edu/~zeilberg/mamarimY/DM85.pdf .
It states the coefficient/determinant identity for commuting matrix entries
and presents Foata's classical combinatorial proof. The source is credit;
Section2 above proves the needed contraction identity directly, without
an external mathematical premise.

The retained TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md, even-parity stratum
equations(E1)--(E2), already uses a quadratic Gaussian determinant kernel,
nilpotence and strict weights to prove safety in its specific ternary
quadratic setting. The all-pair normalized contraction/kernel and original
grading arguments are also retained. The present note gives the sampled
linear test and full-marker prime-degree deduction literally, rather than
importing a coarse earlier principal. No historical-priority claim is made
from a repository search or selected primary reading.

The separate exact checker authenticates this proof and typed scopes,
replays literal contractions/determinants and the countermodels, and keeps
finite controls separate from the written all-field/all-order arguments.
