# Polynomial degree bounds for orders four and six

**Theorem PolynomialFourSixInfinityDegreeGate.** Let \(k\) have
characteristic zero and let \(P,Q\in k[x][D]\), with \([D,x]=1\),
have orders four and six and nonzero scalar commutator. Suppose
their leading coefficients are \(\alpha h^2,\beta h^3\), where
\(d=\deg h\). A polynomial shear of \(D\) fixes these leading
coefficients and gives the estimates
\[
                  \deg p_i\le 2d-4+i,\qquad
                  \deg q_i\le 3d-6+i.
\]
Consequently \(d\ge14\).

The proportional square and cube at the leading differential orders
are forced by the order-nine commutator identity. The
[constant-leading exclusion](MONIC_POLYNOMIAL_FOUR_SIX_STRING_EXCLUSION.md)
and the [degree-one and degree-two exclusions](LOW_DEGREE_LEADING_FACTOR_POLYNOMIAL_FOUR_SIX_EXCLUSION.md)
already give \(d\ge3\), which we use below. The theorem concerns
ordinary polynomial coefficients. The Laurent substitution in its
proof is used only to derive estimates back in the original
polynomial ring.

## Normalize the cubic coefficient at infinity

Apply \(D\mapsto D-A(x)\), where \(A(x)\) is the polynomial part of
\(p_3/(4p_4)\) at \(x=\infty\). In the new pair the ratio \(p_3/p_4\)
has degree at infinity at most \(-1\). Work over an algebraic closure
of \(k\) and put \(x=t^e\), with \(e=1\) if \(d\) is even and \(e=2\)
if \(d\) is odd. Set \(u=e(d-2)/2>0\). In \(k((t^{-1}))\) choose
\(\lambda=\sqrt h\) and use \(E\) proportional to \(\lambda D\),
scaled by a nonzero constant so that
\[
                        [E,t]=t^{1+u}(1+O(t^{-1})).
\]
After constant scaling of the operators they are monic in \(E\).
The coefficient of \(E^3\) in the monic quartic has \(t\)-degree
at most \(u\): for \(E=\lambda D\) before the extra constant scaling,
the coefficient is \(p_3/(\alpha\lambda^3)-6\,d\lambda/dx\).
Both summands have degree at infinity in \(x\) at most \(d/2-1\),
using the bound on \(p_3/p_4\). Depress this monic quartic by
a derivative shear and keep the resulting monic sextic as its
scalar-commutator mate.

## Restrict the growth of the depressed coefficients

Suppose some lower coefficient of the depressed quartic has
\(t\)-degree larger than \(u\) times its order deficit. Let \(r>u\)
be the maximum ratio of these degrees to their deficits.
Give \(E\) weight \(r\) and \(t\) weight one. Contractions
lower weight by \(r-u>0\). The fractional-root recursion and
leading-coefficient descent show that the depressed quartic has
weight \(4r\) and its sextic mate weight \(6r\). They use that
\([P,(P^{j/4})_+]\) has order at most two. Vanishing at
commutator order \(j+3\) makes the leading coefficient of an
order-\(j\) residual constant, so each stage can be subtracted.
The weight bound for every root coefficient is preserved by
the products and contractions in the defining recursion.

The leading Poisson bracket uses \(\{\xi,t\}=t^{1+u}\)
and would have weight \(9r+u>0\), so it is zero. The weighted
Euler comparison gives symbols \(H^2,H^3\), where
\(H=\xi^2+\sigma t^n\), \(n=2r\), and \(\sigma\ne0\).
As in the [higher-root descent](POLYNOMIAL_FOUR_SIX_HIGHER_ROOT_LOCAL_GATE.md#the-hamiltonian-equation-on-a-generic-level),
start with the square of the mate minus the cube of the quartic.
Its degree is less than \(12r\) and its commutator with the
quartic has weight \(6r\). At weight larger than \(3r-u\)
its leading symbol commutes with \(H\), and is a constant
multiple of \(H^j\) for \(2\le j\le5\). The corresponding
constant multiples of \(P,Q,P^2,PQ\) lower its weight without
changing the leading commutator. At weight \(3r-u\) the
leading symbol \(F\) satisfies \(\{H,F\}=\kappa_*H^2\).
Here \(\kappa_*\) is the nonzero scalar commutator after scaling.
Restricting to \(t=w^2,\xi=\pm\sqrt{-\sigma}\,w^n\) twice,
first for \(F\) and then for \(F/H\), forces \(H^2\mid F\).
Both weights are positive, so constants on the parametrizations
must be zero. Thus \(\{H,G\}=\kappa_*\) for a Laurent polynomial
\(G\).

Writing \(G=A(t,H)+\xi B(t,H)\) on the indeterminate level \(T=H\)
gives the zero-\(\xi\) equation
\[
 2T B_t-\sigma t^{n-1}(2tB_t+nB)=\kappa_*t^{-u-1}.
\]
If the least exponent of \(B\) is less than \(-u\), its
derivative gives a nonzero term strictly below the right
side and below all terms from the \(\sigma\) part.
Thus every exponent is at least \(-u\).
For the largest exponent \(j\), we have
\(2j+n\ge-2u+2r>0\), and its \(\sigma\) term has
exponent \(j+n-1>-u-1\). That highest term cannot cancel.
This rules out \(r>u\).

## Return to polynomial coefficients

Use the weight \(u\) on \(E\) and weight one on \(t\).
Contractions do not increase weight. The depressed quartic
has weight \(4u\) by the established coefficient bounds.
The same root recursion now gives \(\operatorname{wt}(L^j)_+
\le ju\), and the same leading-coefficient descent gives
the sextic weight at most \(6u\). The coefficient used in
cubic depression has degree at most \(u\), so the reverse
shear preserves these bounds.

The coefficient of \(D^i\) in \(E^j\) has \(t\)-degree at
most \(ei+ju\). To see this, multiply on the left by \(E\).
Differentiating in \(x\) lowers \(t\)-degree by at least \(e\),
and multiplication by \(\lambda\) increases degree by \(e+u\).
If the new \(D\) is retained instead of differentiating
the coefficient, its extra index accounts for the remaining
\(e\). This proves the bound inductively. The coefficient
of \(D^i\) in the quartic therefore has \(t\)-degree at most
\((4-j)u+ei+ju=4u+ei\) in each contributing summand,
and analogously the sextic bound is \(6u+ei\).
The claimed polynomial degree estimates follow by dividing by \(e\).

The total degrees are exactly \(2d+4,3d+6\), because the
top differential coefficients have degrees \(2d,3d\),
and all lower differential terms have strictly smaller
total degree. This pair cannot be the images of the Weyl
generators under an automorphism. Such an image of a generator
has locally nilpotent adjoint action, by conjugating the
locally nilpotent adjoint action of \(x\). Here the unique
Bernstein leading monomial of \(P\) is a nonzero multiple
of \(x^{2d}\xi^4\). Its iterated Hamiltonian bracket on
\(x\) is nonzero at every step: at step \(j\) the monomial
has exponents \((1+j(2d-1),3j)\), and the new coefficient
multiplier, apart from the leading scalar of \(P\), is
\(4+j(2d-4)>0\). Each nonzero bracket is the leading
symbol of the corresponding iterated commutator. Hence
the adjoint action of \(P\) is not locally nilpotent.

Normalize the mate by the nonzero scalar to give the defining
relation of the first Weyl algebra. The transformed pair
is a hypothetical nonautomorphic pair of total-degree gcd
\(d+2\). By Guccione--Guccione--Valqui
[Corollary 7.4, printed page 42](sources/guccione-guccione-valqui-1111.6100v3.pdf),
this gcd exceeds fifteen. Therefore \(d\ge14\).
