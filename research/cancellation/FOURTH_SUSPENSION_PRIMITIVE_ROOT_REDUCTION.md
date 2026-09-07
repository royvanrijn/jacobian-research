# Primitive-root reduction for the fourth suspension ansatz

## 0. Result and scope

Continue with

\[
T=\frac{U_0(S)+PU_1(S)}{V_0(S)+PV_1(S)},\qquad
\deg U_i,\deg V_i\le2.                                 \tag{0.1}
\]

The old-boundary fan leaves six effective-degree pairs.  This note imposes
the next proposed condition: \(T\) itself is a primitive root variable over
\(k(P)\).  It proves:

> **Primitive-root reduction theorem.**  In the degree-two ansatz, every
> primitive \(T\) reduces over \(k(P)\) to a Möbius function of \(S\).
> The \(m=3\) row is nonbirational, while the \((2,0)\) and \((2,1)\) rows
> either are nonbirational or cancel to lower rows.  Among the remaining
> rows:
>
> 1. a \(1/P\) root translation makes the controlled-divisor coefficient
>    nonpolynomial;
> 2. the only new reciprocal scaling normalizes to \(T=S/P\), but its
>    forced quadratic incidence has an unavoidable \(Q S^2/P\) pole;
> 3. a unimodular genuinely nonmonomial Möbius denominator requires a
>    third-order pole in the controlled-divisor shear.

Thus no fourth suspension occurs with one primitive root variable in these
three reduced strata.  What remains is precisely:

* a nonunimodular Möbius matrix whose determinant primes are cancelled by a
  complete multi-prime ledger; or
* one additional primitive reconstruction variable.

The theorem is scoped to the reciprocal controlled divisor

\[
\mathcal D=1-SQ+PS^2                                  \tag{0.2}
\]

and the affine-in-conjugate incidence used in the existing suspension
search.  It is not a classification of arbitrary non-fibre-preserving
threefold birational charts.

## 1. Birationality collapses the six-row fan

Let

\[
K=k(P).
\]

A rational function \(T=N(S)/D(S)\in K(S)\) generates \(K(S)\) over \(K\)
if and only if its reduced degree is one:

\[
[K(S):K(T)]=\max\{\deg_S N_{\rm red},\deg_S D_{\rm red}\}=1. \tag{1.1}
\]

Apply this to the six old-boundary survivors.

\[
\begin{array}{c|c|c}
(\lambda_U,\lambda_V)&m&\text{primitive-root consequence}\\ \hline
(0,-1)&1&\text{already affine in }S,\\
(1,0)&1&\text{Möbius after deleting quadratic numerator terms},\\
(2,1)&1&\text{degree two or a common-factor collapse to a lower row},\\
(1,-1)&2&\text{affine only when the }PU_1\text{ quadratic term vanishes},\\
(2,0)&2&\text{degree two or a common-factor collapse},\\
(2,-1)&3&\text{quadratic polynomial in }S,\text{ hence nonprimitive}.
\end{array}                                             \tag{1.2}
\]

Consequently every genuinely reduced primitive case is represented by a
Möbius matrix

\[
T=\frac{A(P)S+B(P)}{C(P)S+D(P)}.                       \tag{1.3}
\]

Common-factor strata must be reconstructed before orbit removal, but they
do not define additional primitive-root rows.

## 2. Affine \(1/P\) translations are not controlled

The reduced \((0,-1)\) row has the form, after constant rescaling,

\[
T=aS+b(P)+\frac{c}{P},\qquad a\ne0.                    \tag{2.1}
\]

Since \(T_S=a\), a determinant-one cotangent lift is

\[
\widetilde Q=Q/a+H(P,T).
\]

Substituting

\[
S=\frac{T-b(P)-c/P}{a},\qquad
Q=a(\widetilde Q-H)
\]

in (0.2), the coefficient of \(\widetilde Q\) is

\[
-(T-b(P)-c/P).                                         \tag{2.2}
\]

The shear \(H\) changes only the constant term of \(\mathcal D\), not this
coefficient.  Polynomiality of the transformed controlled divisor therefore
forces

\[
c=0.                                                   \tag{2.3}
\]

The remaining chart is an ordinary polynomial affine root change and lies
in the known orbit.  Thus the apparent one-prime \(P\)-dependent escape is
empty.

## 3. The reciprocal scaling reaches the incidence and then fails

The reduced \((1,-1)\) row has

\[
T=\frac{A(P)S+B(P)}{P}.
\]

The cotangent coefficient contains \(P/A(P)\).  Any zero of \(A(P)\) away
from \(P=0\) is a new Wronskian prime: there \(T_S=0\), while the \(Q\)-linear
residue is nonzero and cannot be cancelled by a shear in \(k[P,T]\).
Therefore the one-prime stratum forces \(A(P)\in k^\times\).

After constant rescaling and a polynomial translation in \(T\), controlled
divisor polynomiality reduces the new case to

\[
\boxed{T=S/P,\qquad \widetilde Q=PQ+H(P,T).}            \tag{3.1}
\]

Indeed

\[
\mathcal D
=1-T\widetilde Q+T H+P^3T^2.                          \tag{3.2}
\]

Its conjugate coefficient is \(\kappa=T\).  The incidence determinant
criterion therefore forces, up to a function of \(P\),

\[
X_T=\lambda T,\qquad
X=\frac{\lambda}{2}T^2+f(P),\qquad\lambda\ne0.          \tag{3.3}
\]

Write the incidence coordinates as

\[
\widetilde B=\widetilde Q+\beta(P,T),\qquad
\widetilde C=Y(P,T)-\widetilde B X(P,T).
\]

The coefficient of the old \(Q\) in \(\widetilde C\) is \(-PX\).  Pulling
back \(T=S/P\) gives

\[
-PXQ=
-\frac{\lambda}{2}\frac{QS^2}{P}-P f(P)Q.             \tag{3.4}
\]

At the component \(q=0\) of \(P=tq=0\), the functions \(t,S,Q\) are
generically units.  The first term of (3.4) has exact order \(-1\).
Neither \(Y\) nor \(\beta\) contains \(Q\), so its residue cannot cancel.

Hence \(T=S/P\) is a valid polynomial controlled-divisor rechart but not a
polynomial incidence suspension.  This distinction is why the incidence
condition must remain in the search after determinant cancellation.

## 4. Unimodular nonmonomial denominators fail one step earlier

Suppose (1.3) has

\[
\Delta=A(P)D(P)-B(P)C(P)\in k^\times.                  \tag{4.1}
\]

Put

\[
L=D(P)T-B(P),\qquad K=A(P)-C(P)T.
\]

Then

\[
S=L/K,\qquad
\partial_ST=K^2/\Delta.
\]

For the determinant-one cotangent lift

\[
\widetilde Q=Q\Delta/K^2+H(P,T),
\]

substitution in (0.2) gives

\[
\mathcal D=
1-\frac{LK}{\Delta}\widetilde Q
+\frac{LK}{\Delta}H
+\frac{PL^2}{K^2}.                                    \tag{4.2}
\]

If \(C\ne0\), then \(K\) is a genuine linear polynomial in \(T\).
Since \(\Delta\) is a unit, \(L\) and \(K\) are coprime.  The last term of
(4.2) has an order-two pole at \(K=0\), while the coefficient multiplying
\(H\) has a simple zero there.  Cancelling the pole requires

\[
H=-\frac{\Delta P L}{K^3}+\text{a term regular at }K,  \tag{4.3}
\]

which is not polynomial in \(P,T\).  Therefore controlled-divisor
polynomiality forces \(C=0\).

When \(C=0\), (4.1) says \(AD\in k^\times\), so \(A,D\) are constants and
the chart is affine in \(S\).  Section 2 then removes its \(1/P\) translation
and leaves only the known affine orbit.

This excludes the genuinely nonmonomial denominator stratum for unimodular
Möbius matrices without any Gröbner basis.

## 5. The reduced continuation

The proposed search order now becomes:

1. **One denominator prime:** closed in Sections 2--3.
2. **Two primes, one primitive root:** only nonunimodular Möbius matrices
   remain.  Compile the primes of
   \[
   \Delta(P),\quad C(P)S+D(P),
   \]
   together with their intersection valuations.  At a prime of \(\Delta\)
   not supported on the denominator, the \(Q\)-linear cotangent residue is
   immediately fatal.
3. **One additional primitive:** this is now the smallest genuinely open
   stratum.  Its valuation equation is
   \[
   i-mj+k\rho=-\delta.
   \]
   The new primitive must either cancel the \(QS^2/P\) residue in (3.4) or
   supply the third-order \(K^{-3}\) principal part in (4.3).
4. **Genuinely nonmonomial denominator:** closed for \(\Delta\in k^\times\);
   only determinant-supported multi-prime ledgers remain.

This is a useful stopping boundary.  Further progress should enlarge the
reconstruction algebra by exactly one primitive variable, rather than
increase the degrees of \(U_i,V_i\) or launch an undifferentiated Gröbner
search.
