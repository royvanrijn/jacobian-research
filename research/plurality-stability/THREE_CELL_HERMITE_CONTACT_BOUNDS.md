# Cumulative Hermite weights at a simplex tie

Let $A$ be a balanced three-cell Gaussian partition,
$T$ the regular simplex, and $W_j(A)$ the total Hermite
weight at degree $j$. These estimates use the exact
deficit $D_\rho=S_\rho(T)-S_\rho(A)$ of a fixed partition.

## A first-moment bound at a tie with nonnegative slope

If $D_r=0$ and $D_r'\ge0$ at $0<r<1$, then
\[
 \begin{split}
 W_1(T)-W_1(A)
    &\le S'_r(T)-2\frac{S_r(T)-1/3}{r}+W_1(T)\\
    &=\sum_{m\ge3}(m-2)r^{m-1}W_m(T).                      \tag{HC1}
 \end{split}
\]
Indeed $S'_r(A)\le S'_r(T)$ and, at the tie, regrouping
the derivative by degree gives
\[
 rS'_r(A)=2(S_r(T)-1/3)-rW_1(A)
                 +\sum_{m\ge3}(m-2)r^m W_m(A).
\]
Discard the nonnegative sum. The bound applies to a
tangential tie and to a right endpoint approached by
simplex-optimal correlations from above. The right side
is explicit from the [simplex value (D0)](THREE_CELL_DEFICIT_AND_PERIMETER.md#exact-correlation-derivative-for-the-three-cell-simplex).

## Eventual cumulative domination for a fixed partition

If $A$ has finite Gaussian perimeter, write $P_\gamma(A)$
for its Gaussian interface area, counting each interface once.
Then its Hermite tail $V_N(A)=\sum_{m>N}W_m(A)$ satisfies
\[
       V_N(A)\sim \frac{2P_\gamma(A)}{\pi\sqrt N}
                                             \qquad(N\to\infty).       \tag{HC2}
\]
Use the finite-perimeter Gaussian heat-content limit of
[Ambrosio, Miranda and Pallara, *Some fine properties of BV functions on Wiener spaces*](https://cvgmt.sns.it/media/doc/paper/2601/finepropertiesfin.pdf),
(3.1) and Theorem 3.2, with time $-\log\rho$. For one
cell the exiting probability divided by the square root
of this time tends to its Gaussian perimeter divided by
$\sqrt\pi$. Summing over cells counts every interface
twice, so $(1-S_\rho(A))/\sqrt{1-\rho}\to
C=2P_\gamma(A)/\sqrt\pi$. Parseval then gives
$\sum_{N\ge0}\rho^N V_N(A)\sim C(1-\rho)^{-1/2}$.

For completeness the monotone Tauberian step here has
a direct proof. Form the measures
$\mu_t=\sqrt t\sum_{N\ge0}V_N(A)\delta_{tN}$ as $t\downarrow0$.
For fixed $s>0$ their Laplace transforms tend to
$C/\sqrt{s}$, the transform of
$C x^{-1/2}\,dx/\sqrt\pi$ on $x>0$.
Tilt by $e^{-x}$ and push forward under $x\mapsto e^{-x}$
to finite measures on $[0,1]$. Their polynomial moments
converge to those of the tilted limit, including the
total mass, so uniform polynomial approximation gives
weak convergence. Undo the tilt on $[0,1]$ in the
original coordinate, whose endpoints have zero limit
mass. With $t=1/N$ this yields
$\sum_{j\le N}V_j(A)\sim 2C\sqrt N/\sqrt\pi$.
Since $V_j(A)$ is nonincreasing, bound $V_N(A)$
between the adjacent mean differences over
$[\lfloor(1-\epsilon)N\rfloor,N]$ and
$[N,\lfloor(1+\epsilon)N\rfloor]$. Take $N\to\infty$
and then $\epsilon\downarrow0$ to obtain (HC2).

For $A$ outside the simplex orbit the
[Gaussian multi-bubble perimeter inequality and uniqueness](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#external-endpoint-and-application-audit)
give $P_\gamma(A)>P_\gamma(T)$. Thus the simplex
strictly dominates the cumulative weights of this
fixed competitor for all sufficiently large $N$:
\[
 \sum_{m\le N}(W_m(T)-W_m(A))
                  =V_N(A)-V_N(T)>0.
\]
This applies to all actual maximizers, by their finite
Gaussian perimeter. The threshold depends on $A$.
Also
\[
 \frac{D_\rho}{\rho(1-\rho)}
   =\sum_{N\ge1}\rho^{N-1}\sum_{m\le N}(W_m(T)-W_m(A)).
\]
Therefore a fixed finite-perimeter partition strictly
improving the simplex would have to reverse the cumulative
comparison at some intermediate degree, despite the
strict comparisons in degree one and at all high degrees.

## Finite cumulative domination fails arbitrarily close to the simplex

For every even cutoff $N\ge2$ and every $N\equiv3\pmod4$,
there are balanced partitions arbitrarily close to the
simplex in Gaussian measure whose cumulative weight through
$N$ is strictly larger. One balanced partition can also
have arbitrarily long alternating sequences of signs
among its odd-cutoff cumulative differences. All these
examples have strictly lower full noise stability at
every positive correlation.

To see this, use the pair coordinates $(a,b)$ of
[the simplex score (CT2)](THREE_CELL_CORRELATION_CONTACT.md#pointwise-score-growth)
with cell $i$ on $a>0$ and cell $j$ its reflection.
Write $D(A,B)$ for the pair gap with conditional mean
coordinates $(A,B)$ and unit covariance. The coefficient
of $\rho^m$ in the unrotated pair score is its degree-$m$
Hermite projection, a polynomial in $(a,b)$. Its highest
homogeneous term is the coefficient of $u^m$ in
$D(ua,ub)$ with $u=\rho/\sqrt{1-\rho^2}$, since all
lower powers of $u$ contribute only lower degrees.
For $a=1,b=0$ use
\[
 D_A(A,0)=\phi(0)\big(e^{-A^2/2}+\tfrac12 e^{-A^2/8}\big).
\]
This follows from the derivative in (CT3) by the Gaussian
integral of $\phi(|V|/\sqrt3)$ for $V\sim N(A,1)$.
Thus the highest term in degree $2k+1$ at $(1,0)$
is nonzero with sign $(-1)^k$. For the even degree
$2k+2$, the coefficient at $b=0$ is zero, but its
derivative in $b$ is nonzero. In fact completing the
square in (CT2) gives
\[
 D_B(A,0)=\tfrac{\sqrt3}{2}\phi(0)e^{-A^2/8}
                      (2\Phi(\tfrac{\sqrt3}{2}A)-1).
\]
Its coefficient of $A^{2k+1}$ has sign $(-1)^k$,
as the two Taylor factors alternate with strictly
positive absolute coefficients. Take $a=1$ and $b$
small on the side making the degree-$2k+2$ homogeneous
term negative. Both sides are inside cell $i$ for
small $|b|$.

Choose a sufficiently large multiple of the chosen
point. When the top homogeneous term is negative,
the sum of the score projections through degree $N$
is negative in a neighborhood of that point.
Choose a small positive-measure patch $E$ there and
exchange its label with that of its reflected copy.
This preserves balance. If $D_{\le N}$ denotes the
sum of the score projections through $N$, the
cumulative gain in the partition objective is
$-4\int_E D_{\le N}\,d\gamma$ plus the nonnegative
cumulative interaction of the transfer with itself.
Consequently it is positive. For any fixed finite list
of odd cutoffs, choose one sufficiently distant point
with $b=0$, where the top-degree signs give the
desired alternation. On a fixed neighborhood of that
point the cumulative interaction is $O(\gamma(E)^2)$
for each finite polynomial projection. Taking $E$
small enough makes the linear term decide all these
signs simultaneously.

The swap leaves the union of the two affected cells
fixed. On the half-plane the label difference of this
union is odd under reflection, and the interaction
uses $M(p,q)-M(p,Rq)>0$. A proper swap reverses its
sign on a positive-measure subset of that union while
preserving the sign elsewhere. Every cross interaction
of the two parts strictly decreases in the odd energy.
The other energy terms are fixed, by the decomposition
(RF2). This proves the full-noise loss and shows why
finite-degree comparisons alone cannot certify the
all-correlation inequality.
