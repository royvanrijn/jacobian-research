# Balanced four-symbol Gaussian quantization

Let \(X\sim\gamma_n\) be standard Gaussian, and let a measurable encoder partition its space into four cells \(A_i\) of probability \(1/4\). The reconstruction associated with symbol \(i\) is \(r_i\in\mathbb R^n\). Write \(m_i=\int_{A_i} x\,d\gamma_n\) and \(W_1(A)=\sum_i\|m_i\|^2\).

## Exact equivalence with the first Hermite objective

Expanding the expected squared error gives the identity
\[
 \sum_i\int_{A_i}\|x-r_i\|^2\,d\gamma_n
  =n-4W_1(A)+\tfrac14\sum_i\|r_i-4m_i\|^2 .
\]
For each partition the unique optimal reconstruction is \(r_i=4m_i\). Hence among balanced encoders, a smaller minimal distortion is equivalent to a larger \(W_1\), and equality cases correspond with those reconstructions. The regular tetrahedral quantizer is interpreted in \(n\ge3\) by adding unused coordinates. By the [centered affine reduction](BALANCED_AFFINE_CENTROID_OPTIMIZATION.md#affine-reduction), a first-Hermite improvement in any dimension would already give an improvement in dimension three. Conversely the dimension-three problem embeds isometrically for \(W_1\) in every higher dimension. The same reduction handles the upper-bound interpretation in dimensions one and two by embedding their partitions into dimension three.

In transport terms, the minimum for a fixed reconstruction tuple is the quadratic transport cost from \(\gamma_n\) to \(\frac14\sum_i\delta_{r_i}\). The affine reduction proves attainment by a balanced assignment even when some sites coincide. The global distortion inequality is therefore a reformulation of the existing first-Hermite gap. The unconstrained propeller optimum of [Heilman–Jagannath–Naor, Theorem 1.1](https://web.math.princeton.edu/~naor/homepage%20files/PropR3.pdf) permits empty cells, and its equality partition does not meet the balance constraint.

## Sharp bound for antipodally paired reconstructions

Suppose the balanced cells have first moments \(u,-u,v,-v\) up to relabeling. This includes partitions with cells paired by \(x\mapsto-x\). For their centroid sites, the winning score bounds the linear assignment, so
\[
\begin{aligned}
W_1(A)
 &\le {\mathbb E}\max(|u\cdot X|,|v\cdot X|)\\
 &= \frac{\|u+v\|+\|u-v\|}{\sqrt{2\pi}}
 \le \sqrt{\frac{W_1(A)}{\pi}}.
\end{aligned}
\]
The equality on the second line follows from the identity \(2\max(|s|,|t|)=|s+t|+|s-t|\). The last step is Cauchy–Schwarz and the parallelogram identity. Thus \(W_1\le1/\pi\) and the distortion is at least \(n-4/\pi\). For \(n\ge2\) equality is attained by the four quadrants in two orthogonal Gaussian coordinates. Their first moments are the four sign choices of \((1,1)/(2\sqrt{2\pi})\) in those coordinates.

For comparison, for four independent standard normals \(G_i\), the tetrahedral value is \(({\mathbb E}\max_i G_i)^2/3\). Indeed the unit tetrahedral directions have the same scores in law as \(\sqrt{4/3}(G_i-\bar G)\), and by symmetry their cone first moments are all the same common multiple of their unit site vectors. The winning-score expectation is the sum of the four first moments dotted with those vectors. Writing \(\phi,\Phi\) for the standard normal density and distribution function, integration by parts gives \({\mathbb E}\max_i G_i=12\int\phi^2 \Phi^2\). The integral is \(1/(2\sqrt\pi)\) times \({\mathbb E}\Phi(Y)^2\) with \(Y\sim N(0,1/2)\). This expectation is the probability that \(Y-G_1,Y-G_2\) are both positive, a centered Gaussian orthant of correlation \(1/3\). By angular measure in the standard Gaussian plane this probability is \(1/4+\arcsin(1/3)/(2\pi)\). Consequently
\[
 W_1(T)=\frac{3}{4\pi}\left(1+\frac{2}{\pi}\arcsin(1/3)\right)^2
 >\frac1\pi ,
\]
where \(\arcsin(1/3)>1/3\) and \(\pi<4\) imply strictness. Thus a better balanced quantizer would have to break antipodal pairing of the optimal reconstruction tuple. General rank-two tuples and nonlocal full-rank configurations are still subject to the compact affine optimization problem.
