# Translated two-divisor chart

Consider

\[
u=y+1/x=(xy+1)/x,\qquad v=x/u=x^2/(xy+1).
\]

Then

\[
x=uv,\qquad y=u-(uv)^{-1},\qquad J\Phi=-1/u.
\]

Thus the outer bracket target is `[A,B]=-u`. The planted inverse
`(A,B)=(uv,u-(uv)^-1)` confirms reciprocal Jacobians and exact cancellation.

## Exact symmetric-box theorem

Let `L_n` be the Laurent space supported on the square
`-n<=i,j<=n`. Its polynomialization kernel under this chart is exactly the
pullback of the ordinary polynomial space `C[x,y]_{<=n}`.

Indeed,

\[
x^a y^b=(uv)^a\left(u-(uv)^{-1}\right)^b
=\sum_{k=0}^b(-1)^k\binom bk
u^{a+b-2k}v^{a-k}.
\]

If `a+b<=n`, every displayed exponent lies in the square. Conversely, for a
nonzero polynomial of degree `d`, the degree-`d` homogeneous terms contribute
distinct Laurent monomials `u^d v^a` at `k=0`; they cannot cancel each other.
Therefore square support forces `d<=n`.

The exact kernel computation verifies this through `n=4`, with dimensions
`3,6,10,15`. Consequently, symmetric boxes in this chart merely re-express the
ordinary dense bounded-degree search. They do not exploit the boundary geometry
enough to make high-degree searches smaller.

## Executed nonsymmetric families

- 13-term diamond: kernel `span(1,x,y)`.
- 25-term square: complete degree-two space.
- 30-term asymmetric box: degree-two space plus
  `xy^2+2y` and `x^2y+x`.
- 42-term shifted box: nine polynomial directions.
- 49-term square: complete degree-three space.

All five collision-normalized systems have exact reduced Gröbner basis `[1]`
over `Q` and are empty at the three tested primes.

The useful next supports are therefore not larger squares. They should be thin
cones or edge bands selected from admissible Newton/Puiseux valuations, where a
large ordinary degree is represented by relatively few Laurent directions.

