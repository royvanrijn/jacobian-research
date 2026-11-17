# Universal quintic calculator: a height-21 witness card

For

\[
t=1+xy,\qquad q=t^2z-\frac45y^2(1+3t),
\]

fix the single polynomial map

\[
F(x,y,z)=\left(
tq,\;
y-\frac{15}{4}xq+\frac54t^2x^3q^5,\;
x(5-3t)+\frac54x^3z-\frac34(xq)^5
\right).
\]

It has \(\det DF=-2\).  At target \((\Pi,B,C)\), its inverse polynomial is

\[
E_{\Pi,B,C}(S)=
\Pi^5S^5-5\Pi S^3-2BS^2+4S-2C.
\]

Write a rational target as a primitive projective quadruple
\([W:P:B_0:C_0]\), meaning
\((\Pi,B,C)=(P/W,B_0/W,C_0/W)\), and put
\(H=\max(|W|,|P|,|B_0|,|C_0|)\).
In the table, \(e(S)\) is the positive-leading primitive integral multiple
of \(E(S)\), and \(h(e)\) is its coefficient height.  A notation such as
\(5:(5)\) is the squarefree factor-degree partition modulo \(5\).

| group | target \([W:P:B_0:C_0]\) | \(H\) | primitive inverse polynomial \(e(S)\) | \(h(e)\) | at most three exact checks | resolvent degree |
|---|---:|---:|---|---:|---|---:|
| \(C_5\) | \([10:10:0:-7]\) | 10 | \(5S^5-25S^3+20S+7\) | 25 | \(2:(5)\); explicit \(\sigma\) below has order \(5\) | none |
| \(D_5\) | \([10:4:-21:20]\) | 21 | \(32S^5-6250S^3+13125S^2+12500S-12500\) | 13125 | \(\Delta=(2040064/78125)^2\); \(\mathcal D=\mathcal D_1\mathcal D_5\), with \(\mathcal D_5\) irreducible by \(3:(5)\); \(11:(2,2,1)\) | 6 |
| \(F_{20}\) | \([10:5:15:4]\) | 15 | \(5S^5-400S^3-480S^2+640S-128\) | 640 | \(29:(5)\); \(\mathcal D(-13/2)=0\); \(\Delta=5(2127/320)^2\) | 6 |
| \(A_5\) | \([5:5:0:-2]\) | 5 | \(5S^5-25S^3+20S+4\) | 25 | \(3:(5)\); \(\Delta=232^2\); \(23:(3,1,1)\) | none |
| \(S_5\) | \([1:-1:-1:-1]\) | 1 | \(S^5-5S^3-2S^2-4S-2\) | 5 | \(5:(5)\); \(43:(2,1,1,1)\); \(\Delta=-16\cdot3\cdot61813\) | none |

These checks use only the five transitive subgroups of \(S_5\).  A
transposition singles out \(S_5\).  A square discriminant and a \(3\)-cycle
single out \(A_5\).  For the dihedral row, Dummit's sextic has orbit
pattern \(1+5\): the rational root puts the group in \(F_{20}\), while
irreducibility of the quintic cofactor forces an element of order five and
hence transitivity.  The square discriminant leaves \(C_5\) or \(D_5\),
and a \((2,2,1)\)-element selects \(D_5\).  For the Frobenius row, an
irreducible quintic with a rational root of Dummit's sextic has group in
\(F_{20}\); its nonsquare discriminant excludes \(C_5,D_5\).  Finally, an
irreducible degree-five algebra with a nontrivial automorphism of order five
is cyclic.

For the cyclic row, in the quotient by
\(5S^5-25S^3+20S+7\), take

\[
\sigma(S)=
\frac{45S^4-20S^3-240S^2+78S+174}{43}.
\]

Direct remainders give \(e(\sigma(S))=0\),
\(\sigma^5(S)=S\), and no smaller positive iterate fixes \(S\).

For the dihedral row, the centered monic polynomial is

\[
f(T)=T^5-5T^3+\frac{42}{25}T^2+\frac{32}{125}T-\frac{128}{3125}.
\]

Its Dummit sextic has primitive integral form

\[
\begin{aligned}
\mathcal D_{\rm int}(X)
&=(125X+516)(
30517578125X^5-63476562500X^4\\
&\quad-1586171875000X^3+2798725000000X^2\\
&\quad+15539492530000X-28777893593024).
\end{aligned}
\]

The quintic cofactor is irreducible modulo \(3\), and the modular pattern
at \(11\) is \((2,2,1)\).  Thus the two solvability certificates in the
five-row card have total resolvent degree \(6+6=12\).

For the Frobenius row, the centered monic polynomial is

\[
f(T)=T^5-5T^3-\frac32T^2+\frac12T-\frac1{40}.
\]

Dummit's sextic from formula (2) in
[Solving Solvable Quintics](https://site.uvm.edu/ddummit/files/2021/04/Solving_Solvable_Quintics__Math_Comp_57_no195_1991__pp_387_401.pdf)
has primitive integral form

\[
\begin{aligned}
\mathcal D_{\rm int}(X)
&=(2X+13)(
20480X^5-51200X^4-1497600X^3\\
&\qquad\qquad
+1947840X^2+24055920X+24084531).
\end{aligned}
\]

Thus \(-13/2\) is the promised rational root.

The exact, oracle-free certificate is:

```bash
.venv/bin/python scripts/verify_universal_quintic_calculator.py
```

The separate exhaustive discovery audit enumerates primitive projective
targets through height \(21\), modulo
\((B_0,C_0)\mapsto(-B_0,-C_0)\), prefilters by exact discriminants and
Frobenius patterns, and uses PARI/GP only to classify the survivors:

```bash
.venv/bin/python scripts/search_universal_quintic_calculator.py --bound 21
```

It finds first heights
\[
H(S_5),H(A_5),H(C_5),H(F_{20}),H(D_5)=(1,5,10,15,21).
\]
Hence \(21\) is the bounded computational optimum for the common height.
The verifier proves the displayed rows independently; the minimality
statement retains the status of an exhaustive bounded PARI computation.
