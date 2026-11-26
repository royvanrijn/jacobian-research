# Hasse–Witt witnesses for factorial moments at native orders

Let \(f\) be a nonzero polynomial in \(m\) variables, and
let \(d\) and \(D\) be its least and greatest total
degrees. Write \(f_d\) for its lowest homogeneous part
and \(M_n=\Lambda_m(f^n)\) for its ordinary factorial
moments. This note gives a finite residue-field test for
infinitely many nonzero \(M_n\), and therefore for a
finite nonvanishing window at every starting order.

This proves a finite all-start window over \(\mathbb C\)
when the lowest part is a smooth plane cubic in three
variables or a smooth quartic surface in four variables.
The higher-degree terms are arbitrary. The resulting window can depend on the
phase and can exceed its number of monomials.

## Exact factorial normalization

Work over a finite extension of \(\mathbb Q_p\).
Assume the coefficients of \(f\) are integral and
\(d\ge m\). Write \(k\) for the residue field and
use bars for coefficient reduction. Fix integers
\(h\ge1\) and \(s\ge1\), with \(p>\max(hd,D)\), and
put \(n=p^s-h\). Choose any \(\sigma\in\mathbb Z_{>0}^m\)
with \(|\sigma|=d\). This vector serves as a factorial
normalization and need not be an exponent in the support.
Set
\[
 \Omega_{s,h}=\prod_{i=1}^m(n\sigma_i)!,
 \qquad R_{s,h}=M_n/\Omega_{s,h},
\]
and define the index sets
\[
  T=\{t\in\mathbb Z_{>0}^m:|t|=d\},\qquad
  U_h=\{z\in\mathbb Z_{>0}^m:|z|\le hd\}.
\]
For a positive vector \(v\) put \((v-\mathbf1)!=\prod_i(v_i-1)!\).
Define a row \(\ell_t=(t-\mathbf1)!\), a column
\(\rho_h(z)=(-1)^{hd-|z|}/(z-\mathbf1)!\), and
\[
 C_{\sigma,h}=\prod_i\frac{(h\sigma_i-1)!}{(\sigma_i-1)!}.
\]
All three can be taken in \(k\), and \(C_{\sigma,h}\)
is a unit.

**Theorem.** The quotient \(R_{s,h}\) is integral.
Let \(W_{s,h}\) have rows indexed by \(T\), columns
indexed by \(U_h\), and entries
\[
           W_{s,h}(t,z)=[X^{p^st-z}]\bar f^{\,p^s-h}.
\]
Then its exact residue is
\[
                  \overline{R_{s,h}}
                       =C_{\sigma,h}\,\ell W_{s,h}\rho_h.      \tag{1}
\]

To prove the formula, consider an exponent vector \(x\)
contributing to \(f^n\). Its total is at least \(dn\)
and at most \(Dn<p^{s+1}\). At a floor level
\(p^j\), \(1\le j\le s\), the baseline exponents
\(n\sigma_i\) have residues \(p^j-h\sigma_i\)
because \(h\sigma_i<p\). Thus
\[
 \begin{aligned}
 p^j\!\left(\sum_i\left\lfloor x_i/p^j\right\rfloor
        -\sum_i\left\lfloor n\sigma_i/p^j\right\rfloor\right)
 &= |x|-dn+mp^j-hd-\sum_i(x_i\bmod p^j)\\
 &\ge -(hd-m)>-p^j.
 \end{aligned}
\]
The floor difference is a nonnegative integer.
At all higher levels both floor sums are zero.
Legendre's formula therefore makes each factorial ratio
\(\prod_i x_i!/\Omega_{s,h}\) integral, and it is a unit
exactly when all the floor differences vanish.

In that equality case, every digit of each \(x_i\) at positions
\(1,\ldots,s-1\) is \(p-1\). Indeed successive equal floor
sums give equal total digits, and the baseline digit at
each such position is \(p-1\) in every coordinate.
The top floor sum is \(d-m\), so write
\[
 x_i=p^s t_i-z_i,\qquad
 t_i=1+\lfloor x_i/p^s\rfloor,\quad z_i=p-(x_i\bmod p).
\]
Then \(t\in T\) and \(|z|=p^sd-|x|\le hd\).
All entries of \(z\) are positive. The description also
works when \(s=1\), with no intermediate digit positions.
Conversely these forms with \(t\in T,z\in U_h\)
have the same floor sums as the baseline, by their
top, intermediate and low digits.

The unit part of \(x!\) is
\((-1)^{v_p(x!)}\) times the product of its digit
factorials modulo \(p\), by Wilson's theorem applied
after each base-\(p\) division. In a surviving factorial
ratio the total valuation signs cancel. The intermediate
digit factorials cancel as well, and the top and bottom
factors give, for each coordinate,
\[
 \frac{(t_i-1)!}{(\sigma_i-1)!}
    \frac{(p-z_i)!}{(p-h\sigma_i)!}
 \equiv
 (-1)^{h\sigma_i-z_i}
 \frac{(t_i-1)!}{(\sigma_i-1)!}
 \frac{(h\sigma_i-1)!}{(z_i-1)!}\pmod p .
\]
The last identity follows from
\((p-r)!\equiv(-1)^r/(r-1)!\) for \(1\le r<p\).
This gives the precise relative factorial weights in (1).

## The finite Frobenius test

Let \(A\) be the square matrix on \(T\) with entries
\[
            A(t,z)=[X^{pt-z}]\bar f_d^{\,p-1},
\]
and let \(B_h=W_{1,h}\). For a matrix \(J\) over \(k\),
write \(J^{(p^e)}\) for entrywise \(p^e\)-th power.
Then
\[
                   W_{s,h}
                = A^{(p^{s-1})}\cdots A^{(p)} B_h.             \tag{2}
\]
The empty product at \(s=1\) is the identity. For the
proof, at \(s>1\) split the power in characteristic \(p\)
as \(\bar f^{p-h}(\bar f^{p^{s-1}-1})^p\).
An exponent \(x'\) chosen from the first factor must
satisfy \(x'+z=py\) for an integral vector \(y\).
Since \(x'\ge0,z>0\), every \(y_i\ge1\).
Also \(|y|\ge((p-h)d+|z|)/p>d-1\).
The second factor requires \(|y|\le d\), since its
coefficient at \(p^{s-1}t-y\) has to have total degree
at least \((p^{s-1}-1)d\).
Hence \(y\in T\), and the second factor uses only the
lowest homogeneous part. Repeating the same argument
with \(h=1\) there proves (2).

The matrix \(A\), up to transposing the index convention,
is the standard coefficient matrix for top coherent
Frobenius on the homogeneous hypersurface. Higher
coefficient matrices \( [X^{p^sv-u}]f^{p^s-1}\) and
their Frobenius products are part of
[Vlasenko's higher Hasse–Witt construction, equations (1)–(3)](https://arxiv.org/pdf/1605.06440).
The projective coefficient formula follows from the
[Čech computation in Vlasenko, equations (2)–(3)](https://www.matrix-inst.org.au/wp_Matrix2016/wp-content/uploads/2017/09/Vlasenko.pdf).
Equation (1) is a factorial-weighted pairing of these
matrices. The use of the initial rectangular block
\(B_h\) accounts for nonhomogeneity and the order offset.

Here is a finite decision for infinitely many nonzero
residues in (1). Write \(|k|=p^e\) and set
\[
 L_r=A^{(p^{r-1})}\cdots A^{(p)}
   \quad(1\le r\le e),\qquad
 C=A^{(p^e)}\cdots A^{(p)}.
\]
At \(s=ej+r\), the matrix in (2) is \(L_r C^jB_h\).
Let \(q=|T|\). The images and kernels of \(C^j\)
stabilize by \(j=q\), with \(C\) invertible on its
stable image. Thus for each \(r\) the pairing is
periodic for \(j\ge q\). It is identically zero there
exactly when
\[
       \ell L_r C^{q+i}B_h\rho_h=0
                      \quad\text{for }0\le i<q.              \tag{3}
\]
The sufficiency follows by applying the characteristic
polynomial of \(C\) to the scalar sequence. Every
nonzero tested entry repeats infinitely often, since
an invertible matrix over a finite field has finite
order. Therefore a failure of (3) at any \(r\) gives
infinitely many exact nonzero moments.

An especially simple case has \(h=1\) and a gap
greater than \(d-m\) from degree \(d\) to every higher
degree in the support. Then the columns of \(B_1\)
outside \(T\) vanish. Those in \(T\) form \(A\), so
(1) uses only \(f_d\). If \(A\) is invertible, take
\(s\) a multiple of \(e\) times the order of
\(A^{(p^{e-1})}\cdots A\). Then \(W_{s,1}\) restricted
to these columns is the identity. Its pairing is
\[
                    \ell\rho_1|_T=|T|=\binom{d-1}{m-1},
\]
a unit because \(p>d\). This gives an infinite
nonzero subsequence for arbitrary homogeneous support
with an invertible \(A\) at one prime, and also for the
higher-degree extensions just specified.

The degree-gap condition is relevant to this last
invertibility argument. For instance, take
\(f=XY(X+Y)-2XY^3\) at \(p>4\) and \(h=1\).
On \(T=\{(1,2),(2,1)\}\), the matrix \(A\) is the
identity. The extra column \(z=(1,1)\) of \(B_1\)
has entries \(2,0\), from using the term \(-2XY^3\)
once. Its weight in \(\rho_1\) is \(-1\).
Consequently the residue (1) vanishes at every
\(p^s-1\) in this example. Its even real moments are
strictly positive.

## From an infinite witness to a fixed window

The factorial moments of any fixed polynomial \(f\)
over a characteristic-zero field form a polynomial
recurrence sequence. Consider the two formal series
\[
       (1-Tf(X))^{-1},\qquad
       \frac1{1-T}\prod_{i=1}^m \sum_{k\ge0}k!X_i^k .
\]
The first is rational. For \(U(X)=\sum_{k\ge0}k!X^k\)
one has \(X^2U'+(X-1)U=-1\), which yields a homogeneous
differential equation by differentiating. Thus the second
series is also D-finite. Multiply the two series in
disjoint variable sets, which preserves D-finiteness
since each variable differentiates only one factor.
Take a primary diagonal in every pair of corresponding
variables to obtain their full Hadamard product.
D-finiteness is preserved at each diagonal step by the
[primary diagonal theorem of Chen, Chyzak, Ma and Zhu, Theorem 4.2](https://arxiv.org/pdf/2506.22011).

The product has coefficients \([X^\alpha]f^n\,\alpha!\)
at \(T^nX^\alpha\), so each \(T\)-coefficient is a
polynomial in \(X\). A nonzero differential equation in
\(T\) for this product, with polynomial coefficients,
specializes to a nonzero one at \(X_i=1\). At each
specialization first remove the largest common power
of \(X_i-1\) from the equation coefficients.
This is a valid cancellation before evaluation because
the series lies in the integral domain \(K[X][[T]]\).
The specialized series is \(\sum_{n\ge0}M_nT^n\).
Comparing coefficients now gives a polynomial
recurrence. Indeed a term \(T^k\partial_T^j\) produces
\((n+j-k)^{\underline j}M_{n+j-k}\) at large \(n\).
For a fixed shift \(j-k\), the polynomials for distinct
\(j\) have distinct degrees, so the resulting recurrence
is nontrivial. This proof applies after any actual
complex coefficient specialization.

If \(M_n\) has infinitely many nonzero entries, choose
a homogeneous recurrence \(\sum_{j=0}^w p_j(n)M_{n+j}=0\)
with \(p_w\ne0\), and an index \(N_0\) after which the
recurrence holds and \(p_w\) does not vanish at integers.
The order \(w\) must be positive, since an order-zero
recurrence would force eventual vanishing. A block of
\(w\) zero moments starting at \(n\ge N_0\) would force
the entire subsequent sequence to vanish by forward
propagation. There are only finitely many earlier
starts. For example, \(N_0+w\) is a sufficient window
at every positive start. The result is finite-width
nonvanishing for the particular phase.

## Smooth ternary cubic lowest part

Suppose \(m=d=3\) and the lowest part defines a
smooth projective plane cubic. First take all coefficients
in a number field. Over a further number field choose a
point on this cubic, identify it with an elliptic curve,
and include its full \(3\)-torsion and a primitive cube
root of unity. Choose a rational prime \(p>\max(D,4)\)
splitting completely in a Galois field containing these
data, excluding the finitely many primes of bad reduction
or denominators. Such primes exist. One elementary
argument uses a monic minimal polynomial of a primitive
integral generator of the Galois field. A nonconstant
integer polynomial has infinitely many prime divisors of
its values: if all belonged to a finite set with product
\(P_0\), choose \(f(t_0)\ne0\) and use the nonconstant
integer polynomial \(f(t_0+f(t_0)P_0u)/f(t_0)\), which
is one modulo every prime in that set. Away from the
discriminant and index, a root of the generator polynomial
modulo \(p\) gives a residue-degree-one prime. Normality
then gives complete splitting.

The reduced curve has its full \(3\)-torsion over
\(\mathbb F_p\), by injectivity of prime-to-\(p\) torsion
under good reduction.
[Silverman, Proposition VII.3.1](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf) Also \(p\equiv1\bmod3\).
Thus \(a_p=p+1-\#E(\mathbb F_p)\) is \(2\bmod3\).
The [Hasse bound, Silverman, *The Arithmetic of Elliptic Curves*, V.1.1](https://www.math.ens.psl.eu/~obenoist/refs/Silverman.pdf)
gives \(|a_p|\le2\sqrt p<p\). In particular \(a_p\)
is nonzero modulo \(p\).

For this degree the set \(T\) contains only
\((1,1,1)\). The scalar \(A\) is \(a_p\bmod p\)
directly by point counting. In \(\mathbb F_p\),
the number of affine zeros is
\(-\sum_{x\in\mathbb F_p^3}\bar f_3(x)^{p-1}\).
Only the exponent \((p-1,p-1,p-1)\) survives in
this sum, since each coordinate exponent must be a
positive multiple of \(p-1\) and the degree is \(3(p-1)\).
Each coordinate sum contributes \(-1\). On the other hand,
the affine zero count is one plus \(p-1\) times the
projective count, hence congruent to \(a_p\).
It follows that \(A\) is invertible. Since \(d-m=0\),
the higher terms in \(f\) do not affect the \(h=1\)
witness. This supplies infinitely many nonzero moments
over the number field.

Finally take arbitrary complex coefficients with a smooth
cubic lowest part and suppose the moments vanish past some
index. The moments for this fixed finite support are
integer polynomials in the coefficient variables.
Their common equations past this index, together with
smoothness of the cubic part, have a point over
\(\mathbb C\). Smoothness is enforced by adjoining an
inverse for the cubic discriminant. Hilbert's
Nullstellensatz gives a point over \(\overline{\mathbb Q}\)
satisfying the same conditions. Equivalently the infinite
moment ideal is finitely generated by Noetherianity, and
the resulting finite algebraic system has an algebraic
point. Its phase contradicts the number-field case.
There are therefore infinitely many nonzero moments
for the original complex phase, and the recurrence
argument gives its finite window.

## Smooth quartic lowest part in four variables

Suppose \(m=d=4\) and the lowest part defines a smooth
quartic in \(\mathbb P^3\). By adjunction and the
cohomology sequence of a projective hypersurface it is
a K3 surface. First assume number-field coefficients.
After a finite extension it has ordinary good reduction
at a density-one set of places. We can choose one with
\(p>D\) at which the quartic model is smooth and the
coefficients of \(f\) are integral.
[Bogomolov and Zarhin, Theorem 0.1](https://arxiv.org/pdf/0902.1548)

Write \(q=p^e\) for that residue field size and \(X\) for
the reduced quartic. Here ordinary means height one.
If \(a_1\) is the linear coefficient of
\(\det(1-t\operatorname{Fr}\mid H^2_{\mathrm{et}})\), then
\(p\nmid a_1\) by the ordinary-prime criterion.
The point-count trace formula gives
\(\#X(\mathbb F_q)=1+q^2-a_1\).
The odd cohomology vanishes for a smooth quartic by
weak Lefschetz, duality and smooth proper base change.
[Bogomolov and Zarhin, §1.2 and Lemma 1.4](https://arxiv.org/pdf/0902.1548)

The direct point count used for cubics applies in four
variables as well. The affine count is congruent to
minus the coefficient of \((X_1X_2X_3X_4)^{q-1}\) in
\(\bar f_4^{q-1}\), since the product of the four
coordinate sums is now \(+1\). Thus this coefficient is
\(\#X(\mathbb F_q)-1\equiv-a_1\bmod p\), which is nonzero.
Formula (2), with its one-element set \(T\), implies
\(A\ne0\). This produces the \(h=1\) witness for the
full phase because \(d-m=0\).

The algebraization of a putative eventually-zero complex
phase proceeds as in the cubic case, using the nonvanishing
quartic discriminant. The recurrence argument therefore
gives a finite all-start window for every such complex phase.

## Exact replay

The [standard-library checker](verify_factorial_hasse_witt_native_orders.py)
compares complete factorial moments to the normalized
pairing, including the relative factorial units for
individual exponent vectors. It checks the rectangular
Frobenius product, different positive normalization
vectors, native offsets larger than one, and the direct
cubic and quartic point-count identities.

    python3 research/factorial/image-rigidity/verify_factorial_hasse_witt_native_orders.py
