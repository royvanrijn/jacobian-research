# Gaussian factor counts for square summands and a power target

When both primitive summands are squares, an exact power on the
target makes their Gaussian integer a power before repair. The real
and imaginary coordinates then give a family of disjoint
cyclotomic radical records. Counting every nonlinear record and
retaining the two linear records for a smaller primitive triple
gives a saving in the minimally repaired count. An odd target
power does not require an opposite pair of exact powers of that
order.

## Original counts and the repair transfer

Fix an integer \(q\ge2\) and nonnegative radical caps
\(\ell_{ab},\ell_c\). Count positive primitive triples
\[
       a=x^2,\qquad b=y^2,\qquad c=R^q,\qquad a+b=c\le X
\]
with \(\operatorname{rad}(ab)\le c^{\ell_{ab}\lambda}\) and
\(\operatorname{rad}(c)\le c^{\ell_c\lambda}\).
For every \(\eta>0\) and all sufficiently small fixed
\(\lambda>0\), the count is
\[
       O_{q,\ell_{ab},\ell_c,\eta,\lambda,\varepsilon}
                   \bigl(X^{(\gamma_q+\eta)\lambda+\varepsilon}\bigr).
                                                               \tag{GS1}
\]
For odd \(q\) one may take
\[
      \gamma_q=\min\left\{\frac{\ell_{ab}}{2\tau(q)},
                       \frac{2\ell_c+3\ell_{ab}}{6\tau(q)+2}\right\}.
                                                               \tag{GS2}
\]
For even \(q\) put
\[
 J_q=\sum_{\substack{e\mid 2q\\ e\ge3,\ e\ne4}}
          \bigl(1+\mathbf 1_{4\mid e}\bigr).
\]
Then one may take
\[
      \gamma_q=\min\left\{\frac{\ell_{ab}}{J_q+4},
                       \frac{4\ell_c+3\ell_{ab}}{16+3J_q}\right\}.
                                                               \tag{GS3}
\]
The bounds may be combined with the original cyclotomic
opposite-power estimates if there is also an exact
power on a summand.
[Cyclotomic counts for opposite powers](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md)

In the central minimal-repair band
\(\mathcal S_{\lambda,\xi}\), use
\(\ell_c=1/3+\xi\) and \(\ell_{ab}=2/3+\xi\).
If \(H_M(c)\le Y\), then
\(c\le Y^{1/(M(1/3-\xi)\lambda)}\).
Applying (GS1) gives at most
\[
              O_{q,M,\lambda,\xi,\eta,\varepsilon}
               \left(
                Y^{(\gamma_q+\eta)/(M(1/3-\xi))+\varepsilon}
               \right)
                                                               \tag{GS4}
\]
primitive quotients. At vanishing balance tolerances the
coefficient relative to \(1/M\) tends to \(1/2\) for a cube
target and \(1/3\) for a fourth-power target. The corresponding
normalized original coefficients are \(1/6\) and \(1/9\).
Each counted primitive quotient has a unique mandatory repair.
All square and power hypotheses are on the original values.
[Least-target-repair equivalence](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## Gaussian splitting on the primitive values

By primitivity \((x,y)=1\). The norm \(x^2+y^2=R^q\) is odd.
Indeed if \(x,y\) were both odd its \(2\)-adic order would
be one, which is impossible for a \(q\)-th power. The
Gaussian integers \(x+iy\) and \(x-iy\) are coprime. Any
common Gaussian prime away from \(2\) divides both \(x\)
and \(y\), and the prime over \(2\) is excluded by the odd
norm. Unique factorization therefore gives
\[
                         x+iy=\epsilon (r+is)^q,
                \qquad r^2+s^2=R
\]
for a Gaussian unit \(\epsilon\) and coprime integers \(r,s\).
A rational prime dividing both \(r,s\) would divide
both \(x,y\). Given \(r,s\) there are at most four choices
of the original triple, coming from the unit. Signs and the
interchange of \(x,y\) leave their joint radical cap unchanged.

Put \(U=r+is\), \(V=r-is\), and use homogeneous
cyclotomic polynomials. Since \(\Phi_e\) is reciprocal for
\(e>1\), the binary forms
\[
                    G_e(r,s)=\Phi_e(U,V)\qquad(e>1)
\]
have integer coefficients. The linear factors from \(e=1,2\)
are \(2is\) and \(2r\), respectively. Taking the product of
all factors of \(U^{2q}-V^{2q}\) gives
\[
                 \left|rs\prod_{\substack{e\mid2q\\e\ge3}}
                      G_e(r,s)\right|=xy .
\]
In particular none of these values is zero. For even \(q\),
the \(e=4\) factor is \(2(r-s)(r+s)\).

For clarity the rational factor counts follow from the
projective roots of \(G_e\). For \(\zeta\) a primitive
\(e\)-th root of unity, its slope is
\(\rho=i(\zeta+1)/(\zeta-1)\), a real number. The
fields \(\mathbb Q(i,\rho)\) and
\(\mathbb Q(i,\zeta)\) are equal, and adjoining \(i\)
to \(\mathbb Q(\rho)\) has degree two. The composite
cyclotomic field is
\(\mathbb Q(\zeta_{\operatorname{lcm}(4,e)})\).
Thus each slope has degree \(\varphi(e)\) over \(\mathbb Q\)
if \(4\nmid e\), and degree \(\varphi(e)/2\) if \(4\mid e\).
All slopes are distinct. It follows that \(G_e\) is
irreducible over \(\mathbb Q\) in the first case, and
has exactly two irreducible rational factors in the second.
These are nonlinear for \(e\ge8\). In counting values, take
primitive integer representatives of the rational factors.

The radical records from distinct factors have disjoint
support outside the primes of \(2q\). Indeed for
\(p\nmid 2q\) work in a residue field of \(\mathbb Z[i]\)
above \(p\). The substitution giving \(U,V\) is invertible
there, and \(r,s\) cannot both vanish by primitivity.
If either \(U,V\) vanishes with the other nonzero, none of
the homogeneous cyclotomic factors vanishes. Otherwise
\(U/V\) is a unit and the separability of \(T^{2q}-1\)
prevents roots of distinct cyclotomic factors from
coinciding. The same separability after the invertible
substitution prevents two rational factors of one \(G_e\)
from having a common root. Their integral scalar contents
can only involve the fixed bad primes: away from these
primes the transformed cyclotomic polynomial is nonzero
and squarefree over the residue field. This argument
includes the linear records \(r,s\) and, for even \(q\),
\(r-s,r+s\). Consequently the sum of the logarithmic
radical masses of all these records differs from
\(\log\operatorname{rad}(xy)\) by \(O_q(1)\).

Each fixed nonzero value of a nonlinear irreducible factor
gives only \(X^{o_X(1)}\) choices of \(r,s\), by the
fixed-value norm argument for a binary form. The coefficients
and the number field are fixed in terms of \(q\).
A specification of an absolute factor value allows two
signs. Enumerating the factor's radical and then its full
value uses only the corresponding radical cap and the
uniform fixed-support enumeration.
[Fixed-value norm recovery for binary factors](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md#primitive-cyclotomic-factors-and-their-fixed-value-counts)

## Counting in one radical box

Work in a dyadic original-height shell with endpoint \(X\).
Take upper cap bins of arbitrarily small fixed normalized
mesh for all the linear and nonlinear factor radicals, in
units of \(\lambda\log X\). Denote the normalized count
exponent in a bin by \(d\). For a nonlinear factor of
cap \(\beta\), the preceding fixed-value count gives
\(d\le\beta\) up to an arbitrarily small allowance.
All bin and small-budget allowances below can be taken
smaller than an arbitrarily prescribed \(\eta\).

For odd \(q\) there are \(J=2(\tau(q)-1)\) nonlinear
factors. Let \(\beta_0\) be the sum of the two linear
caps for \(r,s\). Since the sum of all factor caps is at
most \(\ell_{ab}\) up to the bin errors,
\(\beta_0\le\ell_{ab}-Jd\) with the same convention.
Use the derived primitive sum
\(r^2+s^2=R\). Its total radical cap is at most
\(\ell_c+\beta_0\), and one of its letters has radical
cap at most \(\beta_0/2\). This derived triple determines
the original \(r,s,x,y,R\) up to bounded signs and units.
Its height is comparable to \(X^{1/q}\). The two-letter
count (GD3) on the two linear records gives \(2d\le\beta_0\),
so \(d\le\ell_{ab}/(J+2)\). For a smaller target cap the
one-letter quarter inequality can improve this. It gives
\[
          4d\le\ell_c+\tfrac32\beta_0
                  \le\ell_c+\tfrac32(\ell_{ab}-Jd)
\]
up to the allowances. This proves (GS2).

For even \(q\) use the two pairs of linear records
\((r,s)\) and \((r-s,r+s)\). Write \(\beta_0,\beta_1\)
for their cap sums. The other \(J_q\) rational factors
are nonlinear, so
\(\min(\beta_0,\beta_1)\le(\ell_{ab}-J_qd)/2\)
up to the allowances. For the first pair keep
\(r^2+s^2=R\). For the second use
\[
                     (r-s)^2+(r+s)^2=2R
\]
and divide by its common gcd, which is at most four
because \((r,s)=1\). No linear record vanishes, as its
vanishing would put \((r+is)^q\) on a coordinate axis.
In either case the derived height is comparable to
\(X^{1/q}\). The derived triple determines the original
parameters up to bounded choices, since the gcd and
the signs of the two linear factors suffice to recover
\(r,s\). Its total radical cap is at most
\(\ell_c+\min(\beta_0,\beta_1)\), and its smaller linear
cap is at most half the selected sum. The two-letter count
gives \(2d\le\min(\beta_0,\beta_1)\) and hence
\(d\le\ell_{ab}/(J_q+4)\). The quarter inequality also gives
\[
        4d\le\ell_c+\tfrac32\min(\beta_0,\beta_1)
                 \le\ell_c+\tfrac34(\ell_{ab}-J_qd)
\]
up to the allowances, proving (GS3).
[Two-letter bound (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md#group-dictionary-count)
· [One-letter quarter inequality (LQ7)](ABC_WIDE_COFACTOR_CONGRUENCES.md#small-radical-on-any-one-letter)

All derived budgets are fixed before the height limit.
For a cap \(t\lambda\) in the original shell and
numerical height comparable to \(X^{1/q}\), use a relative
budget \(q(t+\delta)\lambda\) at a fixed positive allowance
\(\delta\), and similarly enlarge the selected letter cap.
The factor \(q\) then cancels on returning to the base
\(X\). There are finitely many possible mesh profiles
at fixed \(q\) and fixed mesh. Taking \(\lambda\)
sufficiently small for the finitely many applications
and summing the shells proves (GS1).

The Gaussian splitting uses the exact target order before repair.
A common multiplier does not preserve the coprime Gaussian ideals
or the displayed power factorization of the primitive values. The
count therefore supplies a bound for this exact-power sector and
does not require a power decomposition of an arbitrary repaired
representative.
