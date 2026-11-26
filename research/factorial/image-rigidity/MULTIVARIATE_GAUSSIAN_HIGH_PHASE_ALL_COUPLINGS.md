# Strict adjacent Gaussian phase order at every coupling for \(\phi\ge\pi/4\)

Use the [exact Gaussian saddle phase](MULTIVARIATE_GAUSSIAN_EXACT_SADDLE_PHASE.md)
for the common normalized time \(t>0\), with raw times \(t/n\)
and \(t/(n+1)\). At every angle with saddle phase
\(\phi_t(\theta)\ge\pi/4\), the exact phase increment between
consecutive degrees is strictly between zero and \(\pi/2\).
Thus the same angular interlacing and quantitative phase-gap
conclusions as in the high-phase section there hold for every
positive coupling.

Only the radial convexity step needs to be replaced.
For the right thimble, write \(z=iA+w\) with \(A=\pi-\theta\)
and \(w=R e^{-i\alpha}\). Its phase equation gives
\(0<\alpha<\pi/2-\phi\le\pi/4\).
Put \(K_k=4\pi^2k^2\).
For fixed \(\alpha\), differentiating the hyperbolic-sine
product formula for the imaginary phase \(I\) gives
\[
 \frac{I_R}{R}
     =-\frac{A\cos\alpha}{2tR}
       +\frac{\sin(2\alpha)}{2t}
       -\sum_{k\ge1}
           \frac{2K_k\sin(2\alpha)}
                {K_k^2+2K_k R^2\cos(2\alpha)+R^4}.
\]
This is strictly increasing in \(R>0\). The first term
strictly increases, the middle term is constant and each
negative summand increases because \(\cos(2\alpha)\ge0\).
The series and differentiation are locally uniform.
Consequently a radial stationary point is a unique strict
minimum. The imaginary phase itself need not be convex.

At any ray different from the saddle ray, the point where
the ray meets the height of the saddle has phase \(<\phi\),
by strict concavity in the horizontal direction.
Hence such a ray cannot have a radial stationary point on
the selected level \(\phi\). On the saddle ray the saddle
is its unique point at level \(\phi\).
The polar angle of the graph therefore has no stationary
point away from the saddle. Its decrease follows from the
endpoint limits as in the
[original high-phase argument](MULTIVARIATE_GAUSSIAN_EXACT_SADDLE_PHASE.md).
Thus \(0<-y'<\cot\phi\).
Both amplitude corrections belong to \((\phi-\pi/2,0)\),
so the phase increment is in \((2\phi-\pi/2,\pi/2)\),
as claimed.
