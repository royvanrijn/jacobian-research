# Powerful codegrees with moving targets

An omitted-power graph for two summands can allow the target
to move. Fix a divisor of the two actual targets and use it
to recover the reduced ratio of the omitted roots. After that
recovery, a weighted ternary equation completes the targets
and the common kernel on the other side. The coefficient
record uses the defect of the target *ratio*, because any
extra prime present on both targets divides a known integer
from the first recovery.

## Actual representations and kernel dictionaries

Fix an integer \(e\ge2\), a sign \(\sigma\in\{1,-1\}\), and
a height parameter \(X\ge2\). Consider actual primitive
positive rows
\[
           P_i y_i^e+\sigma Q x_i^e=T_i,\qquad
           \gcd(P_i y_i^e,Qx_i^e)=1,\qquad i=1,2,            \tag{MP1}
\]
where the complete values and targets are at most \(X^B\)
for fixed \(B>0\). The two positive left labels \(P_1,P_2\)
are prescribed and distinct modulo rational \(e\)-th powers.
The common positive right label \(Q\) belongs to a dictionary
whose labels are distinct modulo rational \(e\)-th powers.
One may take the \(e\)-free kernels, or the canonical full
kernels supported on valuation orders \(e+1,\ldots,2e-1\).
In the full case these positive orders have distinct nonzero
residues modulo \(e\), so the kernel is determined by the
power class. Both options allow different targets for the
two rows.

Prescribe \(D\mid\gcd(T_1,T_2)\), with \(D=1\) allowed.
Let \(Y,Z\ge1\) bound the \(y_i,x_i\), respectively. We
may take both at most \(X^{B/e}\). For a prescribed prime
pool \(\mathcal P\) whose prime product is at most \(X^{B'}\),
define the relative target defect by
\[
 H_e(T_1,T_2;\mathcal P)
   =\prod_{\substack{p\notin\mathcal P\\
             v_p(T_1)-v_p(T_2)\not\equiv0\bmod e}}p .
\]
The pool can depend on the prescribed labels and modulus.
No assumption on the size of \(\gcd(T_1,T_2)\) beyond the
common divisor condition is needed.

## A bound for a fixed pair of left labels

**MP2.** For \(e=3\), the number of actual ordered pairs
in (MP1) with \(H_e(T_1,T_2;\mathcal P)\le H\) is
\[
         O_{B,B',\varepsilon}\!\left(
               H X^\varepsilon\left(1+\frac{(YZ)^2}{D}\right)
                                \right).                  \tag{MP2}
\]
For each fixed \(e\ge36\) the same bound holds with an
additional factor \(R^{3/\sqrt e}\), where
\(R=2\max\{2,Y,ZX^{B/e}\}\). With the stated height caps
this extra term is \(O_e(X^{6B/e^{3/2}})\). The constants
are uniform in the prescribed labels, the pool and \(D\).
Further restrictions on original primitive rows only reduce
the count.

Write
\[
     y_i=g u_i,\quad x_i=h v_i,\qquad
     g=\gcd(y_1,y_2),\quad h=\gcd(x_1,x_2).
\]
Both \(g\) and \(h\) are units modulo \(D\), as are all
the roots and kernel factors. The ratio
\[
                            \frac{u_1v_2}{u_2v_1}
\]
is reduced. A common prime across numerator and denominator
would divide both \(u_i\), both \(v_i\), or a same-edge
pair \(u_i,v_i\). Primitivity rules out these possibilities.
Taking the two row congruences modulo \(D\) gives
\[
              (u_1v_2/u_2v_1)^e\equiv P_2/P_1\pmod D .
\]
At an odd prime power the unit group is cyclic, and at a
power of two it is a product of at most two cyclic groups.
Thus there are at most \((2e)^{\omega(D)}=X^{o(1)}\) unit
root residues, including primes dividing \(e\), for fixed
\(e\). At \(D=1\) there is one unrestricted class.

In each class, the determinant of two distinct reduced
fractions is a nonzero multiple of \(D\). Each numerator
and denominator is at most \(YZ\). Scaling to the unit
square yields angular spacing at least \(D/(2(YZ)^2)\),
so there are \(O(1+(YZ)^2/D)\) fractions in a class.
There are divisor-many allocations into the four factors
\(u_i,v_i\) for each fraction, uniformly at polynomial
height. Thus all four reduced roots are now fixed at the
stated cost.
[Omitted-coordinate recovery for prescribed targets](ERDOS_943_ALL_K_FULL_COHERENT_CYCLE_SAVINGS.md#2-fixed-pairs-omitted-coordinate-codegrees-and-all-heights)

Set
\[
                     E=P_1(u_1v_2)^e-P_2(u_2v_1)^e.
\]
The power-class distinction for \(P_1,P_2\) gives \(E\ne0\).
It has polynomial height. Elimination of the term \(Q h^e\)
from the actual equations gives
\[
                     E g^e=v_2^eT_1-v_1^eT_2 .            \tag{MP3}
\]
In particular \(\gcd(T_1,T_2)\mid E\), since \(g\) is
coprime to both targets. Write \(T_i=A_i z_i^e\) with
\(A_i\) canonically \(e\)-free. Any prime of \(A_1A_2\)
not in the relative target defect belongs to both targets,
and hence divides \(E\). To choose \(A_1,A_2\), first
choose the relative defect radical outside \(\mathcal P\)
at cost \(O(H)\), then enumerate both kernels on the
resulting support together with the fixed primes of
\(\mathcal P E\). Their heights are polynomial, so this
last step costs \(X^{o(1)}\). The amount of mass shared
by the targets does not increase that price.
[Uniform polynomial-height support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

At fixed kernels, (MP3) is the actual weighted diagonal
equation in the roots \(g,v_2z_1,v_1z_2\):
\[
                   E g^e-A_1(v_2z_1)^e
                           +A_2(v_1z_2)^e=0 .             \tag{MP4}
\]
The three roots have collective gcd one. In fact \(g\)
is coprime to each \(v_i\) by same-edge primitivity and
to each \(z_i\) by its coprimality with the targets.
All three coefficients are nonzero. At \(e=3\),
Heath-Brown's coefficient-uniform theorem counts the root
triples at subpower cost. At fixed \(e\ge36\), use the
coefficient-uniform generalized Fermat bound of
Browning–Verzobio in the equal box of side \(R\). Its
three terms are bounded by \(O(R^{3/\sqrt e})\), before
the usual arbitrarily small height allowance.
[Heath-Brown, Theorem 1.3](https://arxiv.org/html/2601.07817v1#S1)
· [Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

Since the \(v_i\) are known, each root tuple recovers both
targets and both \(y_i\). The first original row
determines \(Qh^e\). Reject nonpositive or nonintegral
candidates. The positive \(Q\) is then unique in its
prescribed dictionary by the distinct power-class condition,
and the positive \(h\) is unique. This recovers both rows
and proves the bound. If the two left kernels have equal
power class, \(E\) may vanish and this completion argument
does not apply. A common target is an allowed special
case, with relative defect one. At fixed target the direct
integer recovery from (MP3) does not require the ternary
theorem.

## Consequence for an actual selected-edge graph

Let \(A,B_1\) be the sizes of the left and right kernel
dictionaries for one such graph and let \(M\) be its edge
count. Select one actual primitive row per represented
kernel pair. Every selected target must be divisible by
the prescribed \(D\). Suppose at least a fraction \(\rho>0\)
of the ordered distinct left pairs incident to common right
vertices have relative target defect at most \(H\).
Write \(\Delta\) for the right side of (MP2), including
the larger-order factor when used. Summing over pairs gives
at most \(A(A-1)\Delta\) eligible ordered pairs at the
right vertices. On the other hand their number is at least
\(\rho(M^2/B_1-M)\) by Cauchy–Schwarz. Consequently
\[
                       M\ll B_1+A\sqrt{B_1\Delta/\rho}.    \tag{MP5}
\]
There is a symmetric count when the counterpart eligibility
condition holds on the other side. This statement does not
bound the number of unselected representations at a kernel
pair. A theorem on the original row count also needs to
control that multiplicity, or use an independently imposed
one-witness condition.

For example take \(e=3\), \(n=X^L\), \(A,B_1\le n^{1/6}\)
and \(D\ge(YZ)^2\). If \(H\le n^\beta\) and
\(\rho\ge n^{-\theta}\), the selected-edge exponent is at
most
\[
           \max\{1/6,\, 1/4+(\beta+\theta)/2\}.
\]
It is strictly below \(1/3\) if \(\beta+\theta<1/6\)
with a fixed margin. The inner exponent allowance is taken
after fixing the budget \(L\). One sufficient in-column
hypothesis is that the complete targets have \(e\)-free
kernels of the form \(S_j C_{xj}\) in column \(j\), for a
factor \(S_j\) common to that column and with
\(\operatorname{rad}_{\overline{\mathcal P}}(C_{xj})\le
n^\gamma\) on its selected rows. The relative defect cap
is \(H\le n^{2\gamma}\) independently of the size or number
of common factors \(S_j\). For \(\theta=0\), \(\gamma<1/12\)
suffices here. The common-factor budgets need not give an
individual target class dictionary of comparable size.

These conditions isolate a moving-target form of the
powerful graph count. They retain actual targets and their
shared divisor. A small relative defect, affordable kernel
parts and control of unselected multiplicities remain
arithmetic inputs when applying the graph estimate at the
unrestricted balanced endpoint.
