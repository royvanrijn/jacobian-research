# Hasse–Witt block ranks and marked sphere moments

The source and target dimensions of a character block give lower
bounds for the one-step rank at primes where the Jacobian has
\(p\)-rank zero. In a different prime class, four marked sphere
moments recover a full two-dimensional target block.

## Rank from the coherent dimensions

Let \(W^n=z^{n-1}P(z)\) define the connected tame cyclic cover and
let \(d_i=\dim H^1(\mathcal O)_i\), with characters indexed by
powers of \(W\). Write \(r_i=d_i+d_{-i}\) for the de Rham dimension.

**Proposition (rank from the two coherent dimensions).**
The block \(\Phi_i\) of coherent Frobenius from character \(i\)
to \(pi\) satisfies
\[
        \operatorname{rank}\Phi_i
             \ge \max(0,d_i+d_{pi}-r_i).
\]
The bound is attained if \(d_i=r_i\) or if \(d_{pi}=r_i\).

After the Frobenius twist, the map factors through an injection
of the \(d_i\)-dimensional coherent source into the de Rham space
at \(pi\), followed by the projection onto the \(d_{pi}\)-dimensional
coherent target. The de Rham dimensions match since multiplication
by \(p\) preserves which branch exponents pair nontrivially with
the character. The kernel of the projection has dimension
\(r_i-d_{pi}\), proving the bound by intersection of subspaces.
If either space has dimension \(r_i\), the rank equals the
dimension of the other space. This is the factorization used in
the [two-input strip proof](CM_NEWTON_TWO_INPUT_STRIP_MOMENTS.md#surjectivity-in-every-good-reduction).

For the Jacobian \(J\), its \(a\)-number is \(g-\operatorname{rank}\Phi\).
Indeed, on the contravariant de Rham Dieudonné module \(M\),
\(\ker F=\operatorname{im}V\) and \(\ker V=\operatorname{im}F\).
The coherent space is \(M/\ker F\), and the cokernel of its
induced Frobenius has dimension
\(\dim M/(\operatorname{im}F+\operatorname{im}V)\).
The \(p\)-rank uses the stable rank of coherent Frobenius.
A full Ekedahl–Oort type uses the extra map \(\Psi\) of the
Hasse–Witt triple.
[Moonen, §§2.1, 2.5 and 3.1](https://www.math.ru.nl/~bmoonen/Papers/CDIPart1.pdf).

In the [unequal-multiplicity family](CM_NEWTON_TWO_INPUT_STRIP_MOMENTS.md#unbounded-conductor-and-jacobian-p-rank-zero)
with \(n=6N\) and \(P=(z-1)^N(z-\lambda)^{N+1}\), the eigensheaf
formula gives, for \(1\le i<n\),
\[
       d_i=\left\lceil i/3\right\rceil-\left\lfloor i/6\right\rfloor
                          -\left\lfloor i(N+1)/n\right\rfloor .
\]
The eigensheaf generator at zero removes \(z^{i-1}\),
at the finite roots it removes the two factors with floor
exponents, and the generator at infinity accounts for the
ceiling term. For \(N=9,27\) the gcds of the branch exponents
with \(n\) are \(1,N,2,2N\), so Riemann–Hurwitz gives
\(g=(9N-1)/2\).

**Corollary (one-step ranks at 54 and 162).**
In the two whole-Jacobian \(p\)-rank-zero classes below,
the indicated bounds hold at every good reduction.

| \(n\) | \(g\) | Prime class | \(\operatorname{rank}\Phi\) | \(a(J)\) |
|---:|---:|:---|---:|:---|
| 54 | 40 | \(p\equiv11\bmod54\) | \(\ge8\) | \(\le32\) |
| 162 | 121 | \(p\equiv65\bmod162\) | \(\ge30\) | \(\le91\) |

Summing the block bound for \(1\le i<n\) gives the table.
At 54 the blocks from dimension one to two and two to one
contribute six. Two further dimension-one to dimension-one
blocks have \(r_i=1\). At 162 the corresponding contributions
are 22 and eight. Neither class has a two-to-two coherent block.
The [checker](verify_cm_newton_block_rank_moments.py) evaluates
these finite sums and checks the genera from the signatures.

These bounds are universal. The generic ranks at each fixed good
prime in the respective classes are exactly \(30,90\), and the
generic \(a\)-numbers are \(10,31\). The sum of the per-block
upper bounds \(\min(d_i,d_{pi})\) gives these maximal ranks.
For a block with both coherent dimensions one, let \(b_\nu\)
be the four branch exponents and put
\(e_\nu=\lfloor p\{ib_\nu/n\}\rfloor\).
The sums of the source and target fractional parts are both two,
so \(\sum_\nu e_\nu=2(p-1)\). The block entry is the coefficient
at \(z^{p-1}\) of \(z^{e_0}(z-1)^{e_1}(z-\lambda)^{e_2}\).
This follows by raising the finite normalized eigensheaf
generator divided by \(z\) to the \(p\)-th power, then comparing
with the normalized target generator. Each \(e_\nu\) is between
zero and \(p-1\), so \(0\le p-1-e_0\le e_1+e_2\).
There is at least one binomial term at that index. Every such
term has a distinct exponent of \(\lambda\) and a nonzero
binomial factor in the residue field. Hence the block entry
is a nonzero polynomial in the branch ratio. The one-to-two
and two-to-one blocks already have rank one by the dimension
bound, and the two-to-two case does not occur. Removing the
finitely many roots of the scalar block polynomials therefore
gives the stated generic maxima.

## A two-dimensional target recovered at 229

**Proposition (four marked moments recover a rank-two block).**
Let \(n=162\), \(P=(z-1)^{27}(z+1)^{28}\) and \(p=229\).
Take any phase with that top polynomial on the full strip with
coefficients integral at a place above \(229\). Put
\[
       H_i=t^{10}z^{i-1}(z^2-1),\qquad i=1,2 .
\]
The four residues \(229\mathcal A(H_i f^k)\) at \(k=7,236\)
recover all entries of the coherent Frobenius block from
character \(5\) to character \(11\). Its determinant is
\(137\bmod229\), so the block has rank two.

Here \(5p=162\cdot7+11\) and both eigenspaces have dimension
two. For a character \(j\), write
\[
  R_j=(z-1)^{\lfloor27j/162\rfloor}(z+1)^{\lfloor28j/162\rfloor},
  \qquad E_j=\frac{W^j}{z^{j-1}R_j}.
\]
The two coherent bases are \([E_j/z^i]\), \(i=1,2\), at
\(j=5,11\). We have \(R_5=1,\ R_{11}=z^2-1\).
Raising \(E_5/z^j\), now with \(j=1,2\), to the \(p\)-th power
and dividing by \(E_{11}\) gives
\(z^{p-8-pj}P^7R_{11}\). Thus the matrix entries are
\[
       L_{ij}=[z^{8-i+p(j-1)}]P^7R_{11}.
\]
The polynomial is \((z-1)^{190}(z+1)^{197}\). The four coefficient
extractions give
\[
          L=\begin{pmatrix}148&141\\16&27\end{pmatrix}\pmod{229}.
\]
The same calculation works with \(z+1\) replaced by \(z-\lambda\)
and \(R_{11}\) replaced by \((z-1)(z-\lambda)\).
The determinant is a polynomial in \(\lambda\bmod229\) that is
nonzero at \(-1\), so the fixed-prime result holds off its zero
locus for good branch ratios.

The markers are regular since their angular degrees are at most
three and their corresponding literal heights are nonnegative.
For \(m=7\), the top height in \(H_i f^m\), in the \((z,t)\)
chart, is \(5p-1\). The possible deficits from the top are at
most \(2m=14<p\). Thus the only interval pole modulo \(p\)
comes from the top part of \(H_i f^m\), for both moments.
In the second moment, use \(f^p\) modulo \(p\). Its angles
are multiples of \(p\) and the only angular pairs are
\((0,0),(p,-p)\). For the latter pair the two terms in the
negative lift have even height. The odd lower terms do not
enter either residue by parity. With \(P_j=[z^j]P\) and
\(D_0=[z^0]D\), the formulas are
\[
 \begin{aligned}
    p\mathcal A(H_i f^7)&\equiv L_{i1}/5,\\
    p\mathcal A(H_i f^{236})&\equiv
         L_{i1}\left(\frac{P_1^p}{167}+\frac{D_0^p}{165}\right)
        +L_{i2}P_0^p\left(\frac1{167}-\frac1{165}\right).
 \end{aligned}
\]
All contributing interval denominators are less than \(p^2\).
Since \(P_0\) is a unit, the second-column pivot is nonzero.
Thus the ordinary marked moments recover a block with two
independent coherent output entries. This uses the prime
class \(67\bmod162\) and prescribed markers.

### Unmarked evaluations across marker deformations

The four marked residues can also be obtained from unmarked
moments on the deformed phases \(g_\eta=f+\eta H_i\).
For each marker and each \(k=0,1\), evaluate the normalized
moments of order \(8+kp\) at \(\eta=0,1,\ldots,5\) modulo
\(p=229\). In the low factor \(g_\eta^8\), a term selecting
\(j\) markers has \(t\)-degree between
\(1296-152j-2(8-j)\) and \(1296-152j\). Only \(j=1,4\)
can hit a degree congruent to \(-1\) modulo \(p\).
The high factor \(g_\eta^p\), when present, contributes
either no \(\eta\) or \(\eta^p=\eta\) at these evaluation
points. All interval denominators are less than \(p^2\).
Thus the normalized residues on those points interpolate
a polynomial supported in degrees \(1,2,4,5\).
Its linear coefficient is \(8p\mathcal A(H_i f^{7+kp})\)
reduced modulo \(p\). These are unmarked measurements across
the indicated deformations of \(f\). The four marked
measurements in the proposition use a single phase.

### Rigidity in the two marker directions

**Corollary (two-marker rigidity at 229).**
Under the hypotheses of the rank-two proposition, set
\(g_{\boldsymbol\eta}=f+\eta_1H_1+\eta_2H_2\) and
\[
 M(\boldsymbol\eta)
   =\begin{pmatrix}229\mathcal A(g_{\boldsymbol\eta}^{8})\\
                    229\mathcal A(g_{\boldsymbol\eta}^{237})
     \end{pmatrix}.
\]
Then
\[
 \det DM(0)
   \equiv 8^2\cdot137\cdot\frac{2}{5\cdot167\cdot165}
   \equiv26\pmod{229}.
\]
In particular \(M\) has rank two at the base phase in characteristic
zero. If both moments vanish at the base phase, the common zero is
isolated and nonsingular in this deformation plane.

Differentiation multiplies the two marked residue rows by \(8\)
and \(237\equiv8\), respectively. The first row of the marked
matrix uses the first column of \(L\) divided by 5. The second
row adds a multiple of that column to the second column
multiplied by \(2/(167\cdot165)\), since \(P_0=-1\).
Taking determinants proves the congruence. The Jacobian
criterion at a zero of the two polynomial functions gives the
characteristic-zero assertion.

For a single marker \(H_i\), the eighth-moment residue is
\(a_i\eta+d_i\eta^4\), with
\(a_i=8L_{i1}/5\ne0\). This follows from the same height count used
for interpolation. Its values at \(1\) and \(-1\) differ by
\(2a_i\), so at least one of \(f+H_i,f-H_i\) has a nonzero eighth
moment.
