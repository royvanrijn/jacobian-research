# Simultaneous quadratic torsion from a small shared norm

For a fixed multiplicative twist, both imaginary quadratic fields
can have primitive small norm generators witnessing nontrivial
torsion. A norm over the fixed real quadratic twist field gives
these simultaneous witnesses as binary forms. The factor
degrees allow a squarefree sieve and a quantitative lower count.

**Theorem SN.** Fix a positive squarefree integer \(t>1\) and
\(\ell\in\{3,5\}\). For all large \(X\), there are
\(\gg_{t,\ell}X^{1/\ell}\) negative squarefree \(d\) with \(1\le|d|\le X\) for which both
\[
        K=\mathbf Q(\sqrt d),\qquad K'=\mathbf Q(\sqrt{td})
\]
have ideal classes of exact order \(\ell\). More precisely, both
classes have representative ideals of the same norm
\(m\asymp_{t,\ell}|d|^{1/\ell}\), whose \(\ell\)-th powers have
integral generators with bounded, nonzero quadratic coefficients.
The generators are primitive at every rational prime dividing \(m\).

[Iizuka–Konomi–Nakano, Theorem 1](https://www.jstage.jst.go.jp/article/jmath/68/2/68_899/_pdf)
prove infinitude at every nonzero integer twist for \(\ell=3,5,7\)
using elliptic curves. [Aoki–Kishi](https://doi.org/10.1016/j.jnt.2016.12.007)
construct an infinite imaginary family at the specific fifth twist
with Fibonacci parameters. The statement here gives a power-count
lower bound in absolute discriminant for each fixed positive twist
at \(\ell=3,5\), with the displayed small-norm certificate.

## The two norms from one binary form

Take integers \(a,b\) with \(m=a^2-tb^2>0\), and write
\[
                A+B\sqrt t=(a+b\sqrt t)^\ell .
\]
Thus \(A^2-tB^2=m^\ell\) and \(A,B\) are homogeneous binary forms
of degree \(\ell\). Put
\[
                  D=B\bigl(2A+(t+1)B\bigr).
 \tag{1}
\]
Multiplying \(A+B\sqrt t\) by \(1+\sqrt t\) gives the
rational part \(A+tB\) and the quadratic coefficient \(A+B\).
The resulting norm identity yields
\[
                  (A+B)^2-D=(A+tB)^2-tD=m^\ell.                 \tag{2}
\]
In a rectangle \(Y<a\le2Y,\ -2\delta Y\le b<-\delta Y\)
with sufficiently small fixed \(\delta>0\), both embeddings of
\(a+b\sqrt t\) are positive, \(A\gg_{t,\ell}Y^\ell\) and
\(B<0\). By first taking \(\delta\) small, the positive term
\(2A\) dominates \(|(t+1)B|\), and then
\(|B|\asymp_{t,\ell}Y^\ell\) for this fixed \(\delta\).
Hence \(D<0\), \(|D|\asymp_{t,\ell}Y^{2\ell}\) and
\(m\asymp_{t,\ell}Y^2\).

## Square factors and primitivity

The two factors \(B\) and \(2A+(t+1)B\) of (1) each have
degree \(\ell\). Over \(\mathbf C\), use the invertible linear
coordinates \(a+b\sqrt t,a-b\sqrt t\). In these coordinates each
factor is a linear combination of the \(\ell\)-th powers
with both coefficients nonzero. For the second factor the
coefficients are \(1\pm(t+1)/(2\sqrt t)\), both nonzero for
\(t>1\). So each factor has simple projective roots. They
share no projective root, since \(A=B=0\) would make the
two new coordinates zero. In particular (1) has no repeated
factor and all irreducible factors have degree at most five.

Let \(S\) contain the primes dividing \(2t(t-1)\), all primes
at most \(2\ell\), and the prime factors of the content and
discriminant of the integer form (1). Choose once a residue
base point \((1,b_0)\) with \(b_0\) a multiple of all primes
of \(S\) and \(D(1,b_0)\ne0\). Impose its congruence class modulo
\[
             M=\prod_{p\in S}p^{1+v_p(D(1,b_0))}.
\]
This fixes the valuation of \(D\) at every prime in \(S\), and
\(m\) is a unit at all those primes. Outside \(S\), no prime
square divides every allowed value of \(D\). In fact the form
has nonzero reduction of total degree less than \(p\), so some
pair makes it nonzero modulo \(p\), with the prescribed
congruence modulo \(M\) independent by CRT.

Apply [Manduchi, Proposition 1.1 and Remark 1.1](https://www.numdam.org/article/CM_1995__99_1_33_0.pdf),
the squarefree binary-form sieve with Greaves's degree-six input,
to require \(D\) squarefree outside \(S\). Apply it to the form
\(D(a,-y)\) in positive \(a,y\) and use differences of rectangles
to get the stated rectangle. The local factors are positive
by the preceding residue check, and outside \(S\) their defect
is \(O_{t,\ell}(1/p^2)\), from simple projective roots and
the pairs with both entries divisible by \(p\). This supplies
\(\gg_{t,\ell}Y^2\) parameter pairs. The squarefree reductions
\(D=d s^2\), \(tD=d' s'^2\) have fixed positive integers
\(s,s'\) on this set, since all removable squares are at \(S\)
with fixed valuations. In particular \(d'=\operatorname{sf}(td)\)
and \(|d|\asymp_{t,\ell}|d'|\asymp_{t,\ell}Y^{2\ell}\).
The sieve also implies \(\gcd(a,b)=1\): a common prime
outside \(S\) would divide \(D\) to at least its homogeneous
degree, and \(a\) is a unit at primes of \(S\).

At a prime \(p\mid m\), one has \(p\notin S\). The two
coordinates of \(a+b\sqrt t\) in the reduction modulo \(p\)
are nonzero in precisely one of the two components: \(m=0\),
the pair \((a,b)\) is primitive, and \(\mathbf F_p[\sqrt t]\)
is split because \(a^2=tb^2\) with both \(a,b\) units. Thus
\(B\) is a unit and \(A/B\equiv\pm\sqrt t\bmod p\). Neither
\(A+B\) nor \(A+tB\) can vanish since \(p\nmid t(t-1)\).
By (2) the numbers \(D,tD\) are nonzero squares modulo \(p\).
In each imaginary quadratic field the relevant generator of
(2) is a unit at one prime above \(p\), and contains all
\(\ell v_p(m)\) of its norm valuation at the other. For example
the two residues of \(A+B+s\sqrt d\) sum to \(2(A+B)\), a
unit modulo \(p\). This gives integral ideals
\[
             (A+B+s\sqrt d)=\mathfrak a^\ell,\qquad
             (A+tB+s'\sqrt{d'})=\mathfrak a'^\ell,
             \qquad N\mathfrak a=N\mathfrak a'=m.                 \tag{3}
\]

If either ideal in (3) were principal, its generator would
have norm \(m\). The ideal is not fixed under conjugation,
since \(m>1\) and it uses only one prime above each \(p\mid m\).
A generator therefore has nonzero quadratic coefficient
and gives \(4m\ge |d|\) or \(4m\ge |d'|\) in the respective
field. Both contradict the parameter scales for large \(Y\).
As \(\ell\) is prime, the two classes have exact order \(\ell\).

## Distinct radicands and the fifth twist

By [Stewart–Xiao, Lemma 2.4](https://arxiv.org/pdf/1605.03427),
the parameter pairs with \(0<|D|\le Z\) for which \(D\) has
two representations unrelated by a rational automorphism of
the binary form are \(O_{t,\ell,\epsilon}(Z^{\beta+\epsilon})\).
The form has nonzero discriminant and degree \(2\ell\), and
has the real factor \(b\). The corresponding \(\beta\) is
\(3/(4\sqrt6+3)<1/3\) in degree six and \(1/9<1/5\)
in degree ten. Since \(Z\asymp_{t,\ell}Y^{2\ell}\), these
exceptional pairs are \(o(Y^2)\). The remaining values have
at most the fixed number of rational automorphisms of the
form as their number of representations. The group is finite,
because it permutes the \(2\ell\) distinct projective roots
and any scalar preserving the form has finite order. Since
\(s\) is fixed, the surviving \(\gg_{t,\ell}Y^2\) pairs give
as many distinct \(d\), up to a constant. Choosing \(Y\) a
suitable constant times \(X^{1/(2\ell)}\) proves the lower count.
Both discriminants are bounded by a constant times \(|d|\)
for fixed \(t\), so the same exponent holds if their maximum
is used to order the pairs.

At \(t=\ell=5\), one may relabel the radicands so the first is
prime to 5. This gives the same count in the family for the
[cyclotomic moving-kernel reduction](ABC_FIFTH_REFLECTION_MOVING_KERNELS.md).
Each pair has \(h_5(K)h_5(K')\ge25\), so (1) in that reduction
gives \(h_5(\mathbf Q(\zeta_5,\sqrt d))\ge25\).
The shared norms in (3) may use many prime factors. This
construction does not give two-prime witnesses at the small
split-prime cutoff of the twisted moment problem.

## Upper bound at the same shared norm scale

**Proposition SN2.** Fix squarefree \(t>1\), \(\ell\ge1\), and
constants \(C,L>0\). For \(1\le |d|\le X\), the total number of
records of the form
\[
 \begin{split}
 x^2-dy^2 &=m^\ell,\qquad z^2-tdw^2=m^\ell,\qquad
 1\le m\le C X^{1/\ell},\\
 0<|y|,|w|&\le L,\qquad \gcd(m,xz)=1
 \end{split}                                                     \tag{4}
\]
with integer entries is \(O_{t,\ell,C,L}(X^{1/\ell})\).
One may include signed nonsquarefree \(d\). In the family
of Theorem SN the lower bound matches this exponent with both norm
records in (4) having nonprincipal ideal roots.

Fix \(y,w\) and use \(F=\mathbf Q(\sqrt t)\). Set
\[
       \theta=zy+xw\sqrt t,\qquad
       N_F(\theta)=(y^2-tw^2)m^\ell .                             \tag{5}
\]
The norm factor \(y^2-tw^2\) is a fixed nonzero integer.
Let \(S\) contain the primes dividing \(2ty w(y^2-tw^2)\).
For \(p\mid m\) outside \(S\), the reduction of \(\theta\)
has both its rational and irrational coefficients nonzero
modulo \(p\), while its norm vanishes. Thus \(p\) is split
in \(F\) and precisely one prime above \(p\) divides \(\theta\).
Its order there is \(\ell v_p(m)\). Outside \(S\) all principal
ideal valuations of \(\theta\) are multiples of \(\ell\).

Write \((\theta)=\mathfrak c\mathfrak a^\ell\) by extracting
the nonnegative remainders modulo \(\ell\) at primes above \(S\).
There are finitely many \(\mathfrak c\), and both ideals are
integral. For each class of \(\mathfrak a\) fix one integral
representative \(\mathfrak j\), so
\(\mathfrak a=(\lambda)\mathfrak j\) and
\(\lambda\in\mathfrak j^{-1}\). The ideal
\(\mathfrak c\mathfrak j^\ell\) is principal whenever this
case occurs. Fix a generator \(\theta_0\). Then
\(\theta/(\theta_0\lambda^\ell)\) is a unit. Choose one of
finitely many representatives modulo \(\ell\)-th powers and
absorb the remaining unit power into \(\lambda^\ell\). This gives
\(\theta=\gamma\lambda^\ell\) with finitely many fixed
\(\gamma\in F^\times\), and \(\lambda\) in one of the fixed
fractional ideal lattices.

Equation (4) gives \(|x|,|z|\ll_{t,\ell,C,L}\sqrt X\) for
both signs. By (5) both embeddings of \(\theta\) have the same
bound. Hence each embedding of \(\lambda\) is
\(O_{t,\ell,C,L}(X^{1/(2\ell)})\), giving \(O(X^{1/\ell})\)
points in the fixed rank-two lattice. For the fixed nonzero
\(y,w\), the rational and irrational coefficients of \(\theta\)
recover \(z,x\). Its norm recovers the positive \(m\) and
then (4) recovers \(d\). Summation over the finitely many
coefficients and lattices proves the bound.

For (3), \(y=w=s\) in (4) when the second radical is written
\(\sqrt{td}\), and the proof gives the required coprimality.
The construction fixes \(s\) and the necessary bound \(C\).
The two-sided exponent is thus at the scale where the norm of
a nonrational generator is at least a constant times \(|d|\)
in an imaginary quadratic field.

## Exact control

For an exact fifth-power control, use
\((t,a,b)=(5,16,-1)\). Then
\[
 \begin{gathered}
 m=251,\quad A=1255376,\quad B=-340505,\\
 D=-159261679610,\quad
 -D=2\cdot5\cdot11\cdot41\cdot151\cdot233861,\qquad
 d'=-31852335922.
 \end{gathered}
\]
The two generators are \(914871+\sqrt D\) and
\(-447149+5\sqrt{d'}\). The factorization, equal norms \(251^5\),
coprimality with 251, and the inequalities
\(\min(|D|,|d'|)>4\cdot251\) verify both exact order-five
classes by (3), without computing their entire class groups.
