# Entrywise powers of Gaussian monomial-fibre moments

Use the notation of the
[original Wick coefficient calculation](GAUSSIAN_NONLINEAR_MONOMIAL_FIBRE_FULL_IMAGE.md).
Let \(m\ge2\), \(k_i\ge1\), and let \(\ell\) be the smallest positive
multiple of \(m\) such that all \(k_i\ell\) are even. Write
\[
 d=\ell/m,\quad e=\ell-d,\quad r_i=k_i\ell/2,\quad \rho=\sum_i r_i.
\]
The nonzero coefficient subsequence of the unit for
\(P=W(1+Z^m)\prod_iY_i^{k_i}\), and its angular normalization, are
\[
 a_q=\binom{\ell q}{dq}\prod_i\frac{(2r_iq)!}{2^{r_iq}(r_iq)!},
 \qquad b_q=\binom{\ell q}{dq}\qquad(q\ge0).
\]

For every real \(\theta>0\), both \((a_q^\theta)\) and \((b_q^\theta)\)
are Stieltjes moment sequences with strictly positive Hankel minors of
every size on increasing nonnegative integer rows and columns.
The powered raw sequence is Stieltjes determinate exactly when
\(\theta\rho\le2\). The powered angular sequence is determinate with a
compact representing probability measure. There is also an exact
positive leading coefficient at \(\theta=0\) for every nonempty minor,
given in (6) below.

This extends the \(\theta=1\) results
\(\texttt{GAUSSIAN-NONLINEAR-MONOMIAL-FIBRE-STIELTJES}\) and
\(\texttt{GAUSSIAN-NONLINEAR-FIBRE-HANKEL-STRICT-POSITIVITY}\).
The parameter \(\theta\) acts on the nonzero coefficient subsequences
above. It does not alter the phase or the normalization of the
original Gaussian expectation.

## 1. The exact factorial and beta splitting

Write \((u)_q=\Gamma(u+q)/\Gamma(u)\). Let \(\beta_1\le\cdots\le\beta_\ell\)
be the ordering with multiplicity of
\(\{j/d:1\le j\le d\}\cup\{j/e:1\le j\le e\}\), and put
\(\alpha_j=j/\ell\). For all \(0\le x\le1\), the number of
\(\beta_j\) at most \(x\) is \(\lfloor dx\rfloor+\lfloor ex\rfloor
\le\lfloor\ell x\rfloor\), so \(\alpha_j\le\beta_j\) for each \(j\).
Consequently, with \(L=\ell^\ell/(d^d e^e)\),
\[
 b_q=L^q\prod_{j=1}^{\ell}\frac{(\alpha_j)_q}{(\beta_j)_q}.          \tag{1}
\]
This follows by grouping the factorial factors by residue class,
equivalently using the gamma multiplication formula.
A strict pair \(\alpha_j<\beta_j\) gives the moments of
\(\mathrm{Beta}(\alpha_j,\beta_j-\alpha_j)\). Equal pairs give the
constant one. There is at least one strict pair since
\(1/\ell<\min(1/d,1/e)\).

For a Gaussian factor, direct grouping gives
\[
 \frac{(2rq)!}{2^{rq}(rq)!}
     =(2r)^{rq}\prod_{j=0}^{r-1}
                    \left(\frac{j+1/2}{r}\right)_q.
\]
Each \(u=(j+1/2)/r\) lies in \((0,1)\), and
\((u)_q=q!\,(\mathbb E B^q)\) for \(B\sim\mathrm{Beta}(u,1-u)\).
Combining with (1) gives a finite list of nondegenerate beta
variables \(B_j\) and a constant \(C=L\prod_i(2r_i)^{r_i}\) such that
\[
             a_q=C^q(q!)^\rho\prod_j\mathbb E B_j^q.               \tag{2}
\]

## 2. The beta powers and the determinacy threshold

For a beta law with parameters \(u,v>0\), differentiation of the
gamma ratio, followed by integration from \(s=0\), gives
\[
 \log\mathbb E B^s
   =\int_0^\infty(e^{-s x}-1)
     \frac{e^{-u x}-e^{-(u+v)x}}{x(1-e^{-x})}\,dx,\qquad
        \Re s\ge0.                                              \tag{3}
\]
For the derivative, use
\(\psi(u+v+s)-\psi(u+s)=
 \int_0^\infty e^{-(u+s)x}(1-e^{-v x})/(1-e^{-x})\,dx\).
The integral in (3) converges at zero by its \(O(x)\) cancellation
against an \(O(1/x)\) density, and converges exponentially at infinity.
It is the negative-jump Lévy exponent of \(\log B\), with no drift
in the uncompensated convention. Scaling the density by \(\theta\)
therefore gives \(B_\theta=\exp(-S_\theta)\in(0,1]\), where
\(S_\theta\) is a nonnegative infinitely divisible variable, with
\(\mathbb E B_\theta^s=(\mathbb E B^s)^\theta\).
One construction is a Poisson sum with the displayed positive
intensity, which has finite first moment on jumps smaller than \(1\)
and finite mass above \(1\). Its logarithmic characteristic function
is the exponential of a finite integral for every real frequency,
and hence is everywhere nonzero.

The proof step being transferred is the conversion of a gamma
Mellin quotient to a positive logarithmic Lévy density.
The external proof uses the compensated gamma integral and compares
the positive densities under variation of a stable parameter.
Here the beta difference has finite first moment on small jumps,
so the drift cancels in (3), and scaling the density produces the
required powers.
[Pinned Simon stable-infimum proof, §§2–3, equations (11)–(19)](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/simon-stable-infimum/PROOF.md)

Let \(U_c>0\) have the Urbanik law with moments \((q!)^c\), for
\(c=\rho\theta\). Equations (1)–(3) give the canonical raw measure
as the law of
\[
 C^\theta U_c\prod_j B_{j,\theta}                                \tag{4}
\]
with independent factors. The compact angular measure is the analogous
scaled product from (1) without the Urbanik factor. Its beta factors
have infinite support. For example, the Lévy density in (3) is
positive on every positive interval, so there are arbitrarily many
jumps in a fixed bounded interval with positive probability, which
forces \(S_\theta\) to have unbounded support. A finite product of the
powered beta factors therefore has positive values arbitrarily near
zero and cannot have finite support. The canonical raw law also has
infinite support, since the Urbanik factor has an unbounded support
and the independent beta product is strictly positive.

At large \(q\), Stirling's formula directly in the Wick expression
gives \(\log a_q=\rho q\log q+O(q)\). Hence
\(\sum_{q\ge1}(a_q^\theta)^{-1/(2q)}=\infty\) when \(\theta\rho\le2\),
which proves Stieltjes determinacy by Carleman's criterion. For
\(\theta\rho>2\), the factorial-power sequence has two distinct
Stieltjes probability laws by the Berg–López result verified below.
Multiplicative convolution with the beta product preserves this
distinction. If the two laws have different masses at zero, the
masses remain different after multiplication by a strictly positive
factor. If their masses at zero agree, restrict both laws to
\((0,\infty)\) and push these finite measures forward by \(\log\).
Their transforms differ somewhere by Fourier uniqueness. Multiplication
by the powered beta product convolves these log measures with a law
having nowhere-zero characteristic function. The difference therefore
survives. Both product laws have the same moments by independence.
The angular beta-product measures are compactly supported and hence
determinate.

## 3. Primary verification of the factorial-power input

Berg and López give the Urbanik law via Mellin inversion with
characteristic transform \(\Gamma(1-i\xi)^c\), using the branch
normalized at \(\xi=0\). Their introduction derives the Lévy exponent
from Malmsten's integral and records the moments in (4) of that paper.
Theorem 2.1 and Corollary 2.3 give its large-argument density and
Stieltjes indeterminacy for \(c>2\).
[Berg–López, *Asymptotic behaviour of the Urbanik semigroup*, pp. 1–4 and 6–8](https://arxiv.org/pdf/1310.0252)

The sign and decay needed for this import can be checked directly
in their proof. Move the Mellin line from \(\Re z=1\) to
\(\Re z=T=t^{1/c}\) in
\((2\pi i)^{-1}\int t^{-z}\Gamma(z)^c\,dz\). The intervening half-plane
has an analytic log-Gamma branch, and the horizontal sides vanish by
vertical Stirling for each fixed \(t\). Put \(z=T(1-iu)\) and use uniform
Stirling in the right half-plane. The positive prefactor for the real
integral is
\((2\pi)^{c/2-1} t^{1/c-1/2}e^{-cT}\), and the integral is
\[
 \int_{\mathbb R}e^{cT f(u)}(1-iu)^{-c/2}
                 (1+O_c(T^{-1}))\,du,\qquad
 f(u)=iu+(1-iu)\log(1-iu).
\]
The error is uniform in \(u\). Here
\(\Re f(u)=\tfrac12\log(1+u^2)-u\arctan u\) has derivative
\(-\arctan u\). On \(|u|\le1\) it is at most \(-u^2/4\), and beyond
that interval it decreases linearly in \(|u|\) at a rate at least
\(\pi/4\). On rescaling \(u=w/\sqrt{cT}\), dominated convergence
therefore gives integral asymptotic \(\sqrt{2\pi/(cT)}\), since
\(f(u)=-u^2/2+O(u^3)\). Thus the density is eventually positive and
\(\log e_c(t)=-ct^{1/c}+O_c(\log t)\).
The tail Stieltjes–Krein test now converges when \(c>2\), because
\(t^{1/c}/(\sqrt t(1+t))\) is integrable at infinity. This uses the
tail version, independent of any behavior near zero.
[Hörfelt with Pedersen's appendix, Theorems 4–5](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/7DACBC49C9AF30F7CF3B12C60214EA61/S0021900200000826a.pdf/the-moment-problem-for-some-wiener-functionals-corrections-to-previous-proofs-with-an-appendix-by-h-l-pedersen.pdf)

## 4. Strict minors and their exact small-power coefficient

Let \(s_q\) denote either \(a_q\) or \(b_q\). Andréief applied to the
representing measure of \(s_q^\theta\) gives positivity of every
Hankel minor with increasing nonnegative integer row and column
indices at \(\theta>0\). The two generalized Vandermonde determinants
have positive sign on ordered distinct positive nodes, and every
fixed number of such nodes has positive product measure by infinite
support.

For a more precise certificate, write the index sets as
\(I=(i_0<\dots<i_{r-1})\), \(J=(j_0<\dots<j_{r-1})\). Set
\(L_q=\log s_q\) and
\[
 K_{ij}=L_{i+j}-L_{i+j_0}-L_{i_0+j}+L_{i_0+j_0}
       \qquad(i>i_0,\ j>j_0).                                    \tag{5}
\]
Factor \(e^{\theta L_{i+j_0}}\) from row \(i\) and
\(e^{\theta(L_{i_0+j}-L_{i_0+j_0})}\) from column \(j\). The
remaining matrix has entries \(e^{\theta K_{ij}}\), with ones in
its first row and column. Subtract the first row from each later row
and expand along the first column to obtain
\[
        \det(s_{i+j}^\theta)_{i\in I,j\in J}
                =\theta^{r-1}\det K+O(\theta^r).                 \tag{6}
\]
For \(r=1\) the empty determinant is one.

To see strict positivity for \(r>1\), use the Lévy exponent of the
unpowered sequence. Equations (1)–(3) express \(L_q\), up to an affine
term and the linear compensation needed for the Gamma factor, as an
integral with negative jumps and a positive measure \(\nu\) on
\((0,\infty)\). The Gamma part has density
\(\rho e^{-x}/[x(1-e^{-x})]\). The beta part is a sum of the densities
in (3), and at least one is strictly positive on \((0,\infty)\).
The difference in (5) annihilates all affine and compensating terms
and gives
\[
 K_{ij}=\int_0^\infty e^{-(i_0+j_0)x}
                  (1-e^{-(i-i_0)x})(1-e^{-(j-j_0)x})\,\nu(dx).
\]
The paired factors vanish to second order at \(x=0\), where
\(\nu\) has density \(O(x^{-2})\). Under \(z=e^{-x}\), for strictly
ordered \(z_1,\ldots,z_{r-1}\) in \((0,1)\), the determinant of the
functions \(1-z_b^{i-i_0}\) is positive. Adjoin the exponent \(0\)
and last node \(1\), subtract the last column from each previous
column, and expand along the exponent-\(0\) row. The resulting signs
cancel and the augmented determinant is a positive generalized
Vandermonde. The column-index determinant has the same property.
Andréief gives \(\det K>0\), using truncated measures first if
necessary. Each product of these determinants is nonnegative on
ordered nodes and its integral converges by the entrywise integrability
just checked. This proves the exact order \(r-1\) in (6).
