# Short resultant bounds for projective binomial primes

Suppose two consecutive ordinary factorial moment polynomials have
a monic irreducible common factor \(H(z)\in\mathbb Q[z]\) of degree
\(d\). Even positivity makes \(d\) even. Its roots \(\alpha_i\) are
nonzero. Set
\[
 c=\prod_{i=1}^d\alpha_i,\qquad
 \beta_i=\alpha_i^d/c,\qquad
 B_H(Y)=\prod_{i=1}^d(Y-\beta_i).
\]
These are the invariants of the
[projective factor sieve](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md).
Let \(L_H\) be the least positive integer such that \(L_HB_H\) has
integer coefficients.

Fix \(g\ge1\) and a set \(\mathcal T_g\) of primes where the
projective roots are all units and their reductions satisfy
\(\overline{\beta_i}^{\,2g}=1\). The
[uniform shared slope with binomial residue](MULTIVARIATE_DIFFUSE_PROJECTIVE_FACTOR_SIEVE.md#projective-residues-from-a-uniform-shared-slope)
is one way to supply such primes directly from the two moment
arrays. The primes may depend on the actual profile. Suppose
\(|\beta_i|\le R\) in the complex absolute value for some \(R\ge1\).
At positive diffuse contacts of a fixed core one can use a
uniform \(R\), by the compact annulus bound in the projective sieve.

**Theorem.** Let \(e\) be the number of distinct projective roots.
Either every projective root is a \(2g\)-th root of unity, or
\[
           \left(\prod_{p\in\mathcal T_g}p\right)^e
              \le (1+R^{2g})^e(2R)^{e(e-1)}
                       L_H^{(2g+2e-2)e/d}.                    \tag{1}
\]
At primes where the distinct-root projective polynomial has
squarefree reduction, the bound for their product improves to
\[
           \left(\prod_{p\in\mathcal T_g^{\rm sf}}p\right)^e
                    \le (1+R^{2g})^e L_H^{2ge/d}.             \tag{2}
\]
In particular (1) replaces the exponent \(2gd-d+1\) in the
projective sieve's general monic remainder bound by
\((2g+2e-2)/d\) when estimating \(\prod p\) for an irreducible
common factor. It allows inseparable stepped residues and all
prime characteristics. The conclusion gives a torsion condition
on any sequence where the product of certified moving primes
outgrows this power of the projective denominator.

## Distinct conjugates and denominators

All \(\beta_i\) are conjugates over \(\mathbb Q\). Indeed the
Galois group of a splitting field acts transitively on the
\(\alpha_i\), and \(c\) is rational and fixed. The map
\(\alpha_i\mapsto\beta_i\) is equivariant, so each distinct
\(\beta_i\) has the same fiber size \(r=d/e\).
Consequently \(B_H=B_0^r\), where \(B_0\) is the monic irreducible
polynomial of a projective root. Write its primitive integer
form as \(L_0B_0\), with \(L_0\) the positive leading coefficient.
Gauss's lemma makes \((L_0B_0)^r\) primitive, so \(L_H=L_0^r\).
All primes in \(\mathcal T_g\) are outside \(L_0\), since the
projective roots are units at those primes and each coefficient
of the monic \(B_0\) is a sum of products of these roots.

Put \(F(Y)=Y^{2g}-1\). If \(B_0\) has any common root with \(F\)
then \(B_0\mid F\) by irreducibility, giving the torsion alternative.
Otherwise both rational numbers
\[
               U=\operatorname{Res}(B_0,F),\qquad
                            V=\operatorname{disc}(B_0)
\]
are nonzero. The resultant is homogeneous of degree \(2g\)
in the coefficients of \(B_0\). Thus \(L_0^{2g}U\) is an
integer from the primitive integer form. The discriminant is
a homogeneous integer polynomial of degree \(2e-2\) in all
the coefficients of a degree-\(e\) polynomial, including the
leading coefficient. Applying it to the primitive integer form
gives \(L_0^{2e-2}V\in\mathbb Z\). Equivalently this follows from
the Sylvester formula \(\operatorname{disc}(A)=(-1)^{e(e-1)/2}
\operatorname{Res}(A,A')/a_e\). The formal Sylvester resultant
has total coefficient degree \(2e-1\). Setting \(a_e=0\)
zeros its leading-degree column, because both \(A\) and \(A'\)
lose their declared leading coefficients. Hence the integer
polynomial is divisible by the formal variable \(a_e\).
Over \(\mathbb C\) the root product formulas give
\[
                       |U|\le(1+R^{2g})^e,\qquad
                       |V|\le(2R)^{e(e-1)}.                  \tag{3}
\]
For degree one take \(V=1\).

## The local budget from two Sylvester determinants

Fix \(p\in\mathcal T_g\), and let \(s_p\) be the number of
distinct roots of \(\overline{B_0}\) over \(\overline{\mathbb F}_p\).
Each lies among the roots of \(\bar F\). Hence
\(\deg\gcd(\overline{B_0},\bar F)\ge s_p\).
Also \(\deg\gcd(\overline{B_0},\overline{B_0'})\ge e-s_p\).
The latter follows by factoring over the algebraic closure
of the residue field. A root of multiplicity \(m\) has
multiplicity at least \(m-1\) in the derivative, including
when \(p\mid m\). If the derivative reduces identically to
zero, the gcd has degree \(e\).

For any such polynomial pair, its fixed-degree Sylvester
matrix modulo \(p\) has corank at least the degree of the
residual gcd. Indeed its image consists of polynomials
divisible by that gcd in the space of polynomials of degree
below the sum of the two declared degrees. The divisible
subspace has the indicated codimension. This argument also
works when reduction drops a declared degree. A matrix of
integer \(p\)-adic entries and mod-\(p\) corank \(s\) has
determinant divisible by \(p^s\), by unit row and column
operations followed by factoring \(p\) from each of the
last \(s\) rows. Apply this to the two resultants, using
the monic discriminant identity, to obtain
\[
           v_p(U)\ge s_p,\qquad v_p(V)\ge e-s_p,\qquad
                            v_p(UV)\ge e.                    \tag{4}
\]
Both denominators are units at \(p\). Thus the product of
the \(e\)-th powers of all the certified primes divides the
nonzero integer \(L_0^{2g+2e-2}UV\). Formula (3) proves (1)
since \(L_0=L_H^{e/d}\). When the reduction is squarefree
\(s_p=e\), so \(U\) alone yields (2).

If the binomial condition supplies the finer information
\(\overline{\beta_i}^{\,g}=\sigma\) for one fixed
\(\sigma\in\{1,-1\}\) on a subset of these primes, use
\(F(Y)=Y^g-\sigma\) on that subset. Under the corresponding
torsion alternative the same bounds hold with \(2g\) replaced
by \(g\). In particular no separability hypothesis on the
binomial or projective reductions is needed in (1).

## A prime-power bound at a deep uniform face

Fix an integer \(\ell\ge1\). Suppose at a prime \(p\), all roots
of \(H\) have the same valuation \(s\), and a parent polynomial
has a [binomial face at \(s\)](MULTIVARIATE_DEEP_BINOMIAL_PROJECTIVE_QUOTIENTS.md#face-depth-and-the-all-degree-inequality)
of width \(w\mid\ell\) and depth \(\delta_p>0\). For every root
\(\alpha_i\), the face equation gives
\(\alpha_i^w/a=1+O_p(p^{\delta_p})\) for the same nonzero
face constant \(a\), where the error denotes valuation at
least \(\delta_p\). Write \(A_i=\alpha_i^w/a\). The closure
of \(1+\{x:v_p(x)\ge\delta_p\}\) under multiplication and
division gives
\[
    \beta_i^w=\frac{A_i^d}{\prod_{h=1}^d A_h}
                                  =1+O_p(p^{\delta_p}).
\]
Then \(\beta_i^\ell-1\) has valuation at least \(\delta_p\)
by the integral binomial expansion. This uses the same face
constant across all roots at \(p\).

Let \(U_\ell=\operatorname{Res}(B_0,Y^\ell-1)\) and let
\(\mathcal P_\ell\) be any set of eligible primes with these
depths. Either \(B_0\mid Y^\ell-1\) or the rational nonzero
\(U_\ell\) satisfies
\[
            v_p(U_\ell)\ge\lceil e\delta_p\rceil,
                \qquad p\in\mathcal P_\ell .
\]
Since \(L_0\) is a unit at these primes and \(L_0^\ell U_\ell\)
is an integer, the bound
\[
       \sum_{p\in\mathcal P_\ell}\frac{\lceil e\delta_p\rceil}{e}
                    \log p
             \le \frac{\ell}{d}\log L_H+\log(1+R^\ell)           \tag{5}
\]
follows from the product over the \(e\) distinct projective
roots. The discriminant is unnecessary for this prime-power
bound, including in inseparable residue characteristic.
The deeper quotient test above also applies with several
different common-factor valuations at the same prime.

Near a fixed positive diffuse reference contact, if the
\(d\) projective roots counted with multiplicity include
\(d_*\) limiting \(\ell\)-torsion values, the left side of (5)
can be increased by
\(\frac{d_*}{d}\gamma\log(1/\eta)\), with an adjusted bounded
constant on the right. Here \(\gamma\) is the projective
sieve's effective root-matching exponent. To see this in
the non-torsion alternative, the distinct actual roots
occur with the common multiplicity \(r=d/e\).
For sufficiently small \(\eta\), equal actual roots must
be matched within the same reference cluster. Exactly
\(d_*/r\) of the distinct roots therefore contribute factors
\(\beta_i^\ell-1=O(\eta^\gamma)\) to \(U_\ell\).
The remaining factors stay bounded. This gives the
stated strengthened Archimedean estimate.
