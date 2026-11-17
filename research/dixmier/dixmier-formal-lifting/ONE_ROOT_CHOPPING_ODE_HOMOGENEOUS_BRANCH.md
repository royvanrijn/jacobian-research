# The homogeneous branch in the one-root chopping ODE

This is a local countercontrol to deriving a pure-power polynomial from
the one-root ODE alone. It is **not** a DC-pair, a simultaneous Weyl
operator construction, or a refutation of Zheglov's main theorem.
The additional global/normal-form conditions of an actual pair remain
to be checked before concluding that this branch occurs there.

The primary source is Zheglov,
[arXiv:2410.06959v5](https://arxiv.org/html/2410.06959v5), dated
19 January 2026. Its exact [cached PDF](sources/zheglov-2410.06959v5.pdf)
has 78 pages and SHA256
`897913ef95e011332e86e348803c6c9da5dbd196db17a04f8f951110e9b095e3`.
The relevant printed locations are Theorem 4.2 and its prelude on page
44, equation (4.13) on page 45, Case 2b on page 46, equation (4.14)
on page 47, and the corresponding induction case b on page 49.
The prelude has \(1<l<d\); the numerical control below respects this.
The overall adoption boundary is recorded in the
[primary-source audit](ZHEGLOV_DC1_PRIMARY_SOURCE_AUDIT.md).

## The complete local polynomial solution

Let \(k\) have characteristic zero and let \(l,d,m,n\) be integers
with \(1<l<d\). Set
\(\beta=mn-m-n\geq0\), \(v=l\beta+1\). Put
\[
 u=b+\alpha y,\quad \alpha\ne0,\quad g=u^d,\qquad
 \mathcal L_v(H)=-l gH'+v g'H.
\]
Equation (4.13), after the source's scalar normalization, is
\[
                         \mathcal L_v(H)=g^{\beta+1}.     \tag{1}
\]
Its particular polynomial solution is
\[
       a u^r,\qquad a=\frac1{\alpha(d-l)},\quad r=d\beta+1.
\]
The homogeneous equation has polynomial kernel
\[
 \ker\mathcal L_v\big|_{k[y]}=
 \begin{cases}
    k\,u^s,& l\mid d,\quad s=d\beta+d/l,\\
    0,&l\nmid d.
 \end{cases}                                                   \tag{2}
\]
Indeed \(\mathcal L_v(u^j)=\alpha(dv-lj)u^{j+d-1}\).
The monomials in \(u\) are a basis, so this diagonal calculation
proves the entire kernel statement and the general solution
\[
                         H=a u^r+C u^s\quad(l\mid d).     \tag{3}
\]
For \(C\ne0\), its degree is \(s>r\). No coefficient condition
setting \(C=0\) is stated in the displayed Case 2b ODE inference.
For equation (4.14), the same coefficient calculation applies with
its actual coefficient \(v\) and scalar right side; polynomial
exponents must be checked rather than assuming \(M-n\) integral.

The extra term does not automatically defeat subrectangularity.
If the ordinary order is \(s\), the corner
\((v,s)\) is proportional to \((l,d)\), since
\(s=dv/l\). Thus a repair can retain the higher-degree branch
in a form resembling Case 2a. The pure-power conclusion requires
an additional argument.

## The next equation does not force this coefficient to vanish

Assume \(l\mid d\) and put \(w=s-r=d/l-1>0\). For an integer
\(e\geq2\), define
\[
                   H_1=H^e-C^e u^{es}.                   \tag{4}
\]
The second term is homogeneous for \(\mathcal L_{ev}\), so the
product rule and (1) prove exactly
\[
                  \mathcal L_{ev}(H_1)
                           =eH^{e-1}g^{\beta+1}.         \tag{5}
\]
All particular exponents in (4) are
\(er+jw\), \(0\leq j\leq e-1\), strictly below the homogeneous
exponent \(es\). Consequently there is no resonance obstruction
to (5), and, for \(C\ne0\),
\[
                 \deg H_1=es-w=es-d/l+1.                 \tag{6}
\]
This is exactly the terminal degree pattern with \(M_1=d/l\)
and \(M_2=1\). A terminal degree ceiling excludes a *new*
homogeneous term proportional to \(u^{es}\); it does not exclude
the earlier parameter \(C\).

Moreover, (4) is the required local top-\(x\) chopping recursion,
not merely an independently selected solution of (5). If the next
chop is \(Q_2=Q_1^e-C^eP^{m_1}\), where \(m_1=es/(dn)\) is integral,
the two top-\(x\) terms have the same degree, and their coefficient
polynomial is exactly (4). Derivative contractions lower the
\(x\)-degree, so multiplication of these top faces is commutative.
This statement concerns those faces; it constructs no full operators.

## Exact shifted numerical control

Take
\[
 l=2,\ d=4,\ n=2,\ m=3,\ \beta=1,\quad u=1+y,\quad e=4.
\]
Then \(d_2=d/\gcd(l,d)=2\), \(v=3\), and the polynomial family is
\[
 \begin{split}
 g&=u^4,\\
 H&=\tfrac12u^5+Cu^6,\\
 H_1&=\tfrac1{16}u^{20}+\tfrac C2u^{21}
                      +\tfrac{3C^2}{2}u^{22}+2C^3u^{23}.
 \end{split}                                                   \tag{7}
\]
The following identities hold over \(\mathbb Q[C,y]\):
\[
 \begin{split}
 -2H'g+3Hg'&=g^2,\\
 H_1&=H^4-C^4g^6,\\
 -2H_1'g+12H_1g'&=4H^3g^2.                               \tag{8}
 \end{split}
\]
For \(C\ne0\), the orders and chopping data are
\[
 p=8,\ q=12,\ q_1=6,\ q_2=23,\quad n_0=2,\ n_1=4,
 \quad m_1=3,\quad M_1=2,\ M_2=1.
\]
They obey
\[
 2\cdot12-6=18,\quad4\cdot6-23=1,\quad18+1=p+q-1=19,
\]
and sequential gcds are \(4,2,1\). In particular \(q_1\) lies
in the divisible-order branch allowed by Lemma 4.2, while the next
order is terminal. No stated numerical degree ceiling in these local
alternatives kills \(C\).

Both individual axis-face patterns can also coexist. Put \(a=1+x\)
and consider the ordinary commuting-variable polynomial controls
\[
 P_*=a^4u^8,\quad Q_*=a^6u^{12},\quad
 U_*=a^3H,\quad V_*=a^{12}H_1.
\]
Their top-\(x\) faces have the coefficients \(g^2,g^3,H,H_1\);
their top-\(y\) faces have coefficients
\(a^4,a^6,Ca^3,2C^3a^{12}\). Thus both corner descriptions and
the ordinary versus \(\max(j-i)\) order equalities are compatible
with the local family. Also
\(V_*=U_*^4-C^4P_*^3\). These controls do **not** give the first
chop: \(Q_*^2-P_*^3=0\ne U_*\). They are therefore not admitted
simultaneous symbols of the original Weyl pair and its full chopping
sequence. Their purpose is to show that the separately displayed
axis-face and order conditions do not by themselves remove the term.

## The precise centering obstruction in Lemma 6.1

Lemma 6.1 on printed page 64 uses the pure-power terminal statement
to put both edge polynomials into the same centered form by one
translation. The family (7) shows why the extra branch matters in
the exceptional case. Centering \(g=u^4\) puts \(u=y\), after which
\[
 g=y^4,\qquad H_1=2C^3y^{23}+\tfrac32C^2y^{22}+\text{lower terms}.
\]
For \(C\ne0\), the coefficient of \(y^{22}\) is nonzero.
Under a further common translation \(y\mapsto y+t\), the two
normalized subleading coefficients become
\[
 \frac{[y^3]g(y+t)}{4[y^4]g(y+t)}=t,\qquad
 \frac{[y^{22}]H_1(y+t)}{23[y^{23}]H_1(y+t)}
                         =t+\frac{3}{92C}.
\]
Their difference cannot be removed by translation. Thus the local
family also fails the simultaneous centering claimed from the
pure-power inference, without failing either ODE or the next chop.
Here \(d_2=2,l_2=1\), and the prior sequential gcd is \(2\), matching
the source's exceptional arithmetic. The orders nevertheless differ:
the centroid mismatch becomes a coefficient at head minus \(N\),
whereas the [conditional Laurent lemma](CONDITIONAL_LAURENT_COEFFICIENT_INVARIANCE.md)
uses head minus one. A monomial-face transport with gap \(N>1\)
already supplies the latter zero, even with unequal centroids; see the
[exact PBW transport calculation](PBW_TRANSPORT_AND_FINAL_COEFFICIENT_COMPARISON.md).
Thus this control does not obstruct that local invariance statement.
It remains a countercontrol to the displayed simultaneous-centering
inference, without realizing the actual operators of Lemma 6.1.

## Concrete remaining obligation

The retained controls show that the one-root ODE inference and the
next ODE/order/face constraints alone do not imply \(C=0\).
The specific next verification obligation is to prove that an
additional condition of the **actual simultaneous** Weyl pair and its
normal forms forces \(C=0\), or to repair Theorem 4.2's induction and
terminal-edge conclusion while retaining \(C\). In particular the
pure-power terminal statement, later used in Lemma 6.1, cannot be
obtained from (1) and (5) alone.

[The checker](check_one_root_chopping_ode.py) uses independent exact
Fraction polynomial arithmetic, with \(C\) left indeterminate. The
[receipt](one_root_chopping_ode_controls.json) pins the primary PDF,
this note and checker, equations (8), the numerical data, and the
individual face controls. The evidence is a scoped local countercontrol,
not a counterexample to Theorem 4.2 under all its global hypotheses,
a refutation of DC1, formal verification, or external review.
