# A uniform real obstruction in reciprocal degrees 4 modulo 8

**Theorem.** Let N in R[u] be monic of degree 2k, and put
f(z)=N(z^2)/z. Assume its critical polynomial R_N=2uN'-N has 2k
distinct nonzero roots, and the values N(u)^2/u at these roots comprise
k distinct nonzero values, each taken twice. Suppose their equal-value
pairs are the fixed-point-free reciprocal pairs u <-> lambda/u for one
lambda != 0. If every critical point of f is real, then k is even.

Consequently there is no all-real good reciprocal model of degree
n=4k congruent to 4 modulo 8, over R or any subfield of R. The conclusion
includes transcendental real coefficients. It does not require rational
coefficients, elimination, or a good-reduction hypothesis.

The [degree-twelve note](RECIPROCAL_DEGREE_TWELVE_OBSTRUCTION.md) supplies
the even critical-degree normalization and weighted-Rolle bounds.
The proof below includes those bounds and the additional endpoint sign
that makes them contradictory in every odd degree k. Its degree-twelve
rational obstruction and retained eliminant remain separate results.

## Normalization and the negative first critical value

All roots u of R_N are positive, because the critical points of f have
the form z=+-sqrt(u) and are real. Thus lambda>0. Choose rho=sqrt(lambda)>0,
and write

\[
 F(v)=\rho^{-2k}N(\rho v),\qquad
 L=4k-1,\qquad
 2vF'(v)-F(v)=L\Sigma(v).
\tag{1}
\]

The polynomial Sigma is monic of degree 2k, has distinct positive roots,
and is reciprocal with constant coefficient one. Evaluating (1) at zero
gives F(0)=-L. Let alpha be its smallest root. Reciprocity and the absence
of fixed roots give 0<alpha<1, with largest root alpha^{-1}.

Put h(v)=F(v)/sqrt(v) for v>0. Equation (1) gives

\[
 h'(v)=\frac{L\Sigma(v)}{2v^{3/2}},
\]

and hence

\[
 h(v)=-\frac{L}{\sqrt v}
       +\frac L2\int_0^v\frac{\Sigma(t)-1}{t^{3/2}}\,dt.
\tag{2}
\]

The integral converges because Sigma(t)-1=O(t); the integration constant
is zero since (F(v)+L)/sqrt(v) tends to zero as v tends to zero.

Writing the positive roots as r_1,...,r_(2k), the constant coefficient
one gives

\[
 \Sigma(t)=\prod_{j=1}^{2k}(1-t/r_j).
\]

For 0<t<alpha, every factor is strictly between zero and one. Thus
0<Sigma(t)<1 there. Formula (2) proves the strict endpoint sign

\[
                         h(\alpha)<0.                       \tag{3}
\]

This endpoint argument uses positivity of all the critical roots. It
does not assume the prescribed equalities of critical values.

## The forced order of the two signs

For t>=0 define

\[
 A(t)=h(e^t)+h(e^{-t}),\qquad
 B(t)=h(e^t)-h(e^{-t}).
\tag{4}
\]

Write Sigma(v)=v^k s(v+v^{-1}), with s monic of degree k.
Its roots are distinct and exceed 2. If 0<t_1<...<t_k are the positive
logarithms of the roots of Sigma greater than one, direct differentiation
gives

\[
 \begin{split}
 A'(t)&=L\sinh((k-\tfrac12)t)s(2\cosh t),\\
 B'(t)&=L\cosh((k-\tfrac12)t)s(2\cosh t).
 \end{split}
\tag{5}
\]

On t>0, both derivative weights are strictly positive. Thus their only
zeros are t_1,...,t_k.

At each reciprocal pair, equality of the nonzero squared critical values
means exactly one of A(t_j)=0 or B(t_j)=0. Call these signs P and Q,
respectively, and let their counts be p,q, with p+q=k. This is the same
sign convention as the degree-twelve note. The nonzero-value hypothesis
ensures that a point cannot have both signs.

Between consecutive P zeros, Rolle's theorem requires a Q point.
Therefore q>=p-1. Since B(0)=0, Rolle applied before the first Q zero
and between consecutive Q zeros requires p>=q.

Suppose now that k is odd. These inequalities force p=q+1. There are
p-1=q gaps between consecutive P points, and each must contain a Q
point. They exhaust all Q points. Consequently the entire ordered
sign pattern is

\[
                         P,Q,P,Q,\ldots,Q,P.                 \tag{6}
\]

In particular t_k=log(alpha^{-1}) has sign P. If q>=1, the point
t_(k-1) has sign Q, so B(t_(k-1))=0. The monic polynomial
s(2cosh t) is strictly negative between its last two roots. Equation (5)
therefore gives

\[
                         B(t_k)<0.                          \tag{7}
\]

For k=1, the same conclusion follows from B(0)=0 and the negative sign
of s(2cosh t) on 0<t<t_1.

On the other hand, the last P condition says

\[
                  h(\alpha^{-1})+h(\alpha)=0.
\]

Using (3), it follows that

\[
             B(t_k)=h(\alpha^{-1})-h(\alpha)
                         =-2h(\alpha)>0,
\tag{8}
\]

contradicting (7). This proves the theorem.

For even k the Rolle bounds give p=q, and the last sign is Q.
The argument does not exclude those degrees, and asserts no existence
there. It also does not exclude nonreciprocal all-real paired pencils.

## Exact controls and scope

The [checker](verify_odd_k_real_reciprocal_obstruction.py) checks the
coefficient reconstruction, both derivative identities, and the
integral-antiderivative identity symbolically in bounded degrees. It
checks (3) exactly on positive reciprocal root sets and enumerates the
finite sign-order controls. These controls supplement the uniform proof;
they are not its quantifiers. Its
[receipt](odd_k_real_reciprocal_obstruction_verification.json) and
[source ledger](odd_k_real_reciprocal_obstruction_source_ledger.json)
retain the actual source pins and replay scope.

No novelty comparison, formal verification, independent external replay
or external human review is asserted. The mathematical argument has
received internal critical review.
