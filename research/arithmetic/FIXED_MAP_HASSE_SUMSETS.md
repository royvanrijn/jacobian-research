# Binary sums of Hasse-failing fibre parameters

Fix the Jacobian-one map \(\Phi:\mathbb A^3_{\mathbb Q}\to\mathbb A^3_{\mathbb Q}\)
and the targets \(y_a=(-1,32a/9,(8a+1)/3)\) from the fixed-map manuscript.
For \(a\in\mathbb Q\setminus\{0,1\}\), the complete fibre is the root
scheme of \((X^3-a)(X^2+X+1)\). Call the fibre Hasse-failing if it has a
point over every completion of \(\mathbb Q\) and no point over \(\mathbb Q\).
[Fixed map and fibre identification](../../unpublished/fixed-map-hasse-failures/sections/02-fixed-map.tex)

Put
\[
 C_{18}=\{n\in\mathbb Z_{\ge1}:p\mid n\Longrightarrow p\equiv1\pmod{18}\},
 \qquad \mathcal H_{18}=C_{18}\setminus\{m^3:m\in\mathbb Z_{\ge1}\}.
\]
Every \(a\in\mathcal H_{18}\) is a Hasse-failing parameter. Indeed,
\(v_p(a)=0\) for every \(p\equiv2\pmod3\), and \(a\equiv1\pmod9\)
makes \(a\) a \(3\)-adic unit cube. A positive integer is a rational
cube exactly when it is an integer cube, since a reduced rational cube
with denominator greater than one is nonintegral. In particular the
excluded values \(0,1\) never enter \(\mathcal H_{18}\).
[Rational classification](../../unpublished/fixed-map-hasse-failures/sections/03-local-global.tex)

## Binary count

Define the ordered count and its arithmetic factors by
\[
 \begin{aligned}
 R_{\mathcal H_{18}}(N)&=\#\{(a,b)\in\mathcal H_{18}^2:a+b=N\},\\
 \mathcal E_{18}&=\{p:p\nmid18,\ p\not\equiv1\pmod{18}\},\\
 J_{18}(N)&=\prod_{\substack{p\mid N\\p\in\mathcal E_{18}}}\frac{p-1}{p-2},\\
 A_{18}&=18c_{18}^{\,2}\prod_{p\in\mathcal E_{18}}
                            \left(1-\frac{1}{(p-1)^2}\right).
 \end{aligned}
\]
Here
\[
 c_{18}=\frac{1}{\Gamma(1/6)}
       \lim_{s\downarrow1}(s-1)^{1/6}
              \prod_{p\equiv1\ (18)}(1-p^{-s})^{-1}.
\]
Equivalently, with \(\chi_0\) the principal character modulo \(18\)
and \(r(p)\) the order of \(p\) in \((\mathbb Z/18\mathbb Z)^\times\),
\[
 c_{18}=
 \frac{3^{-1/6}}{\Gamma(1/6)}
 \left(\prod_{\substack{\chi\bmod18\\\chi\ne\chi_0}}L(1,\chi)\right)^{1/6}
 \prod_{p\in\mathcal E_{18}}(1-p^{-r(p)})^{1/r(p)}.
\]
The sixth root is taken by continuation from real \(s>1\). The Euler
product in the latter expression and the product in \(A_{18}\) converge
absolutely because their factors are \(1+O(p^{-2})\). The marginal
constant \(c_{18}\) is positive by the subgroup theorem, and every factor
in \(A_{18}\) is positive since \(p\ge5\) on \(\mathcal E_{18}\).

**Corollary.** Assume the fixed Dirichlet zero-free half-plane and
weighted-dilation inputs of Section 3 of the binary splitting-semigroup
paper. Uniformly as \(N\to\infty\) through \(N\equiv2\pmod{18}\),
\[
 R_{\mathcal H_{18}}(N)=A_{18}J_{18}(N)
   \frac{N}{(\log N)^{5/3}}
   \left(1+O\bigl((\log N)^{-1/49152}\bigr)\right).
\]
On every other residue class the count is identically zero. In particular,
every sufficiently large \(N\equiv2\pmod{18}\) is a sum of two parameters
of Hasse-failing fibres of the same fixed map \(\Phi\).
[Analytic inputs](../../results/pointwise-binary-splitting-semigroups/inputs.tex)

**Proof.** Apply the general binary theorem with \(L_0=18\),
\(H_1=H_2=U=V=\{1\}\), \(A=B=1\) and \(I=[0,1]\).
Both cores are \(C_{18}\), with \(\rho_1=\rho_2=1/6\) since
\(\varphi(18)=6\). Thus \(2-\rho_1-\rho_2=5/3\).
There is exactly one allowed residue pair \((1,1)\) when \(N\equiv2\pmod{18}\),
so the finite factor \(\Gamma_N\) is \(18\) on this class and zero elsewhere.
The good primes forbidden in both slots are precisely \(\mathcal E_{18}\).
The theorem therefore gives the displayed main term for the ordered
\(C_{18}\)-pairs, with relative error \(O((\log N)^{-1/49152})\)
by its choice \(\delta=1/(8192\varphi(18))\).
The two formulae for \(c_{18}\) are the corresponding specializations
of the marginal formulae. [Binary splitting-semigroup theorem, Theorem 11.1](../../results/pointwise-binary-splitting-semigroups/theorems.tex) · [Logarithmic saving, Section 8](../../results/pointwise-binary-splitting-semigroups/assembly.tex)

For each positive cube below \(N\), choosing it in either slot determines
the other summand. Hence the removal of pairs containing at least one
cube costs at most \(2\lfloor (N-1)^{1/3}\rfloor\). Since \(J_{18}(N)\ge1\),
this is \(o(NJ_{18}(N)/(\log N)^{5/3+1/49152})\), so the main term
and the relative error are preserved. The main term tends to infinity
uniformly on the admissible class, giving eventual existence.
Every member of \(C_{18}\) is \(1\pmod{18}\), which proves the exact
zero outside that class. \(\square\)

The formula counts pairs from the prescribed subclass \(\mathcal H_{18}\)
of the positive integer Hasse-failing parameters. For example,
\(38=19+19\) uses the prime \(19\equiv1\pmod{18}\), which is a noncube.
The fibre at \(19\) is Hasse-failing by the local classification alone.
