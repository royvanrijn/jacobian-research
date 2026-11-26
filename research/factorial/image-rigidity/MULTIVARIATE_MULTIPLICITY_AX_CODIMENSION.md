# A codimension bound for several varying tail multiplicities

Let \(L_{i\ell}(u)\), \(1\le i\le s,\ 1\le\ell\le k\), be analytic
logarithms of nonzero algebraic functions of a nonconstant complex
variable \(u\), and assume their \(sk\) logarithms are linearly
independent over \(\mathbb Q\) modulo constants. Work on a
nonvertical analytic curve, allowing ramified local parameters.
Let \(r_1,\ldots,r_k\) be meromorphic functions there, and set
\[
           G_i=\sum_{\ell=1}^k r_\ell L_{i\ell}(u),
           \qquad Z_i=\exp(G_i).
\]
It suffices to use a germ where all these expressions are regular.

**Theorem.** If \(Z=(Z_1,\ldots,Z_s)\) belongs to an algebraic
subvariety of the coefficient torus over \(\overline{\mathbb C(u)}\)
of codimension strictly greater than \(k\), then all \(r_\ell\) are rational constants.
For ordinary binomial moment polynomials this bounds the
degrees of possible common factors under simultaneous
variation of several tail multiplicities. At \(k=1\) the bound
is strong enough to control any common root.

## Ax and the rank of the multiplicity relations

We use the one-derivation exponential theorem of
[Ax, On Schanuel's conjectures](https://annals.math.princeton.edu/1971/93-2/p03).
For analytic germs \(f_1,\ldots,f_b\) which are linearly
independent over \(\mathbb Q\) modulo constants, it gives
\(\operatorname{trdeg}_{\mathbb C}\mathbb C(f_i,e^{f_i})\ge b+1\)
when \(b>0\). Here the differential field has constants
\(\mathbb C\), and \(D e^{f_i}/e^{f_i}=Df_i\). Some derivative
is nonzero by the independence assumption. For an arbitrary
list one applies the theorem to a basis modulo constants,
giving the same lower bound with \(b\) its rational rank.

Applied first to the \(L_{i\ell}\) this shows that they are
algebraically independent, even over \(\mathbb C(u)\).
Indeed all their exponentials are algebraic over \(\mathbb C(u)\),
so adding \(u\) and the logarithms to the field gives
transcendence degree at least \(sk+1\).

Let \(\lambda\) be the dimension over \(\mathbb Q\) of the
space of relations modulo constants among the \(L_{i\ell}\)
and \(G_i\). For a relation write
\[
                       \alpha\cdot L+\beta\cdot G=c .
 \tag{1}
\]
The projection to its rational coefficient row \(\beta\) is
injective, since the \(L\) have no rational relation modulo
constants. For \(\lambda\) independent relations the
\(\beta\)'s therefore form independent rows of a rational
matrix \(B\). Substitution of \(G=Lr\) gives a linear system
in the \(r_\ell\) with coefficient matrix \(BL\).
Its rank over \(\mathbb C(L)\) is \(\min(\lambda,k)\):
the entries of \(BL\) are independent indeterminates as
polynomials in \(L\), by extending the row change of
variables with a rational invertible matrix.

If \(\lambda\le k\), at most \(k-\lambda\) entries of \(r\)
are needed to generate the others over \(\mathbb C(L)\).
Writing \(d\) for the codimension of the algebraic target
for \(Z\), an upper bound on the transcendence degree of
\(\mathbb C(L,G,e^L,Z)\) is therefore
\[
                  sk+1+(k-\lambda)+(s-d).
\]
Here the \(1\) accounts for \(u\), over which \(e^L\) and
the coefficients of the target equations are algebraic. Ax gives the lower bound \(sk+s-\lambda+1\).
The inequalities are incompatible with \(d>k\).

If \(\lambda>k\), take \(k+1\) independent relations. They give
\(Mr=y\), where \(M\) is a generic \((k+1)\)-by-\(k\) matrix,
and \(y\) is affine-linear in the entries of \(M\) and the
remaining independent variables. The last statement follows
from the invertible rational change of variables on each
column of \(L\). We claim \(y=Mv\) for a rational constant
vector \(v\). In \(\det[M\mid y]=0\) the coefficients of
each remaining variable, and the term of degree one in
each column of \(M\) with no remaining variables, show
that the corresponding constant vectors of \(y\) vanish:
a vector in the span of every full-rank \(M\) is zero.
For the part depending on column \(M_j\), write it \(T_j M_j\).
The separate column gradings distinguish these terms, so
\(\det[M\mid T_j M_j]=0\) separately for each \(j\).
Fixing \(M_j=w\ne0\) and varying all other columns forces
\(T_j w\) into the intersection of all such column spaces,
which is \(\mathbb Cw\). Thus \(T_j\) is a scalar matrix.
Its entries are rational from (1), so the scalar is rational.
The resulting identity \(y=Mv\) gives \(r=v\) by the rank
of \(M\). This finishes the proof.

The codimension hypothesis is sharp for this abstract
logarithmic argument. For example, take \(s>k\) and
\(L_{i\ell}(u)=\log(1+c_{i\ell}u)\) with distinct nonzero
rational \(c_{i\ell}\), chosen so that their first
\(k\)-by-\(k\) coefficient block is invertible. The divisors
give the required independence. Solve the first \(k\)
equations \(Lr=v_0\) for a nonzero constant column \(v_0\).
Then \(r\) has a pole at zero and \(ur\) is bounded there.
The \(Z\)'s belong to the torus defined by fixing their
first \(k\) coordinates to \(\exp(v_0)\), of codimension \(k\).

## A uniform height consequence for ordinary coefficients

Fix degrees \(n,N=n+1\) and a core. Normalize each parent
to constant coefficient one. Let its \(s=n+N\) other
coefficients be nonzero algebraic functions of \(u\) times
\[
                     \prod_{\ell=1}^k R_{i\ell}(u)^{r_\ell}.
 \tag{2}
\]
The algebraic prefactors must be analytic units at zero.
Here \(u=1/M\), each \(R_{i\ell}\) is rational and equals
\(1+O(u^{a_\ell})\) with a fixed positive integer \(a_\ell\),
and the \(r_\ell\) are nonnegative
integer multiplicities. For example (2) holds for affine
height blocks with bounded gaps and zero signed gap in
each block. Assume the \(sk\) rational functions have
independent divisor columns over \(\mathbb Q\), so
\(L_{i\ell}=\log R_{i\ell}\) satisfy the theorem.
Assume also that \(u^{a_\ell}r_\ell\) are bounded in the
profiles in question, as holds in the bounded-curvature regime
for balanced blocks at scales \(M^{a_\ell}\) with positive
curvature increments.

The locus of coefficient pairs with a common divisor of
degree at least \(d\) is algebraic in this torus of exact
degrees. It has codimension at least \(d\): factor out
a common degree-\(d\) factor with constant coefficient one,
and parametrize the factor and both remaining factors by
\(d+(n-d)+(N-d)=s-d\) coefficients. Multiplication by the nonzero algebraic prefactors
preserves this codimension over \(\overline{\mathbb C(u)}\).

If \(d>k\), any sequence \(M\to\infty\) of profiles with such
a divisor has a subsequence with one fixed integer vector of
multiplicities. For that vector the common divisor condition
holds identically as a function of the height. To see this,
write \(\tau_\ell=u^{a_\ell}r_\ell\), so (2) becomes
\(\exp(\sum_\ell \tau_\ell L_{i\ell}(u)/u^{a_\ell})\), analytic also
at \(u=0\). Pass to a limit of the bounded \(\tau\)'s.
The pulled-back algebraic conditions form an analytic set
in \((u,\tau)\) with finitely many irreducible components
locally. A component with \(u\) varying and some
\(\tau_\ell/u^{a_\ell}\) varying would contain a nonvertical analytic
curve with both varying: choose a smooth point with these
two differentials nonzero and a local tangent avoiding both
kernels. The theorem forbids it. All nonvertical components
therefore lie on curves \(\tau_\ell=c_\ell u^{a_\ell}\) with rational
constant \(c_\ell\). A component with \(u\) constant cannot
supply the sequence. Consequently a subsequence lies on a
fixed such curve with all \(c_\ell\) integers. The factor
condition on that curve is an analytic identity near zero,
by its infinitely many zeros.

In particular, if a separate fixed-multiplicity argument
excludes such a persistent divisor for every multiplicity
vector, there is a uniform height cutoff above which the gcd has degree at most \(k\). No effective cutoff in
the several-multiplicity setting is asserted here.
For two multiplicities this gives the bound two whenever
fixed-multiplicity persistent divisors of degree at least
three have been excluded. At adjacent ordinary orders a
nonconstant common factor has even degree, because the
even parent has no real zeros.

When the heights are all affine functions of \(M\), the
[all-directions affine-ray theorem](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md)
supplies the fixed-multiplicity exclusion, provided each admissible
profile has distinct exponent columns and at least one moving
unequal row. Rational affine functions can be restricted to an
arithmetic progression of integral heights. If a fixed-multiplicity
identity held, it would persist along that progression, in
contradiction with the affine-ray theorem. If the normalized
divisor columns are independent at every adjacent pair below the
absolute cutoff, there is one uniform height beyond which the
gcd degree is at most \(k\) in every degree.
