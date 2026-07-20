# A one-parameter adversarial continuation benchmark

Restrict the target of the degree-`(7,6,4)` map to

\[
(a,b,c)=(-1/4+s,0,0),\qquad s<1/4.
\]

Every such real target has three distinct real affine preimages.  Nevertheless,
as `s` tends to `1/4`, two paths escape to infinity, the third stays finite,
and the affine Jacobian determinant remains the constant `-2`.  This separates
path divergence from the usual diagnostics for a singular affine endpoint.

Constant determinant does not imply good numerical conditioning.  Along an
unbounded branch the Jacobian and its inverse are exactly

\[
DF=\begin{pmatrix}
\frac92t^3&\frac34t&-\frac18\\
\frac32t^2&\frac{25}4&-\frac3{8t}\\
-\frac{17}2&-\frac3{4t^2}&\frac1{8t^3}
\end{pmatrix},\quad
(DF)^{-1}=\begin{pmatrix}
-\frac1{4t^3}&0&-\frac14\\
-\frac3{2t}&\frac14&-\frac34t^2\\
-26&\frac32t&-\frac{27}2t^3
\end{pmatrix}.
\]

Hence normwise condition numbers can grow on the order of `|t|^-6 =
delta^-3`, despite `det(DF)=-2`.  This is a conditioning test as well as a
path-at-infinity test.

## Exact solution paths

For `x != 0`, the primitive-element model reduces to

\[
P(T)=-2T^2+\tfrac12-2s,
\qquad t^2=\delta:=\tfrac14-s.
\]

The reconstruction formulas therefore give

\[
(x,y,z)=\left(-\frac1{2t},3t,26t^2\right),
\qquad t=\pm\sqrt\delta.
\]

The `x=0` chart supplies the third path

\[
(x,y,z)=(0,0,-\delta).
\]

Thus the two unbounded paths have the exact growth laws

\[
|x|=\frac12\delta^{-1/2},\quad |y|=3\delta^{1/2},
\quad z=26\delta,
\]

while the bounded path converges to `(0,0,0)`.  At `s=0` these specialize to
the three collision points already recorded in `FACTS.md`.

## Affine and projective endpoint counts

At `s=1/4`, the target is `(0,0,0)`.  Its affine fiber contains only the
bounded endpoint `(0,0,0)`.  The other two paths have not collided at a
singular affine solution: `det(DF)=-2` everywhere.

In homogeneous source coordinates `[X:Y:Z:W]`, multiply
`[x:y:z:1]` by `2t` to obtain the regular parametrization

\[
[X:Y:Z:W]=[-1:6t^2:52t^3:2t].
\]

Both signs of `t` approach the same point `[-1:0:0:0]` on the hyperplane at
infinity.  They are two local branches because the target parameter satisfies
`s=1/4-t^2`.  Consequently the nearby affine count is three, the endpoint
affine count is one, and the missing multiplicity two is supported at one
projective point at infinity.

This is a useful warning for projective output: deduplicating endpoint
coordinates alone merges two paths.  Path labels, local parameters, winding,
or endgame multiplicities must be retained.

## Benchmark protocol

`scripts/nonproper_fiber_homotopy.jl` starts from the exact three-point fiber
at `s=0`, attempts each labeled path through `delta=10^-1,...,10^-10`, and
emits CSV rows comparing observed and predicted norms, coordinates, and
normalized projective endpoints.  A lost affine path is recorded as `LOST`
rather than aborting the benchmark.  A successful adapter should report:

- three nonsingular affine solutions at every positive `delta`;
- two norms asymptotic to `1/(2*sqrt(delta))` and one norm equal to `delta`;
- two projective residuals tending to zero at `[-1:0:0:0]`;
- Jacobian condition-number growth of order `delta^-3` on those paths;
- only one finite affine solution at the exact endpoint;
- projective multiplicity/branch count two at the common infinite endpoint.

The numerical script intentionally tracks the original high-degree equations,
not the reduced quadratic.  The quadratic is used only as an exact oracle.
This exposes path-loss thresholds, overflow or scaling problems, incorrect
finite/infinite classification, and accidental merging of projective branches.

Run the exact oracle first:

```bash
.venv/bin/python scripts/nonproper_fiber_benchmark.py
```

Then run the numerical adapter in an environment with Julia and
HomotopyContinuation.jl:

```bash
julia --project=. scripts/nonproper_fiber_homotopy.jl
```

The parameter-homotopy call follows the current HomotopyContinuation.jl
parameter-system interface. Numerical results are environment-dependent and
are not an exact certificate.

## Verified reference run

On 20 July 2026 with Julia 1.12.6 and HomotopyContinuation.jl 2.21.0, the
double-precision affine tracker followed all three paths through `delta=10^-8`.
At that point the two unbounded paths had norm about `5000`, relative coordinate
error `2.6e-10`, projective residual `5.3e-14`, and Jacobian condition number
`3.125e22`.  Both were reported `LOST` at `delta=10^-9`; the bounded path
continued through `delta=10^-10`.  This is the intended adversarial behavior:
the loss threshold reflects extreme affine conditioning, while the exact
certificate proves that neither path terminates at a singular affine point.
