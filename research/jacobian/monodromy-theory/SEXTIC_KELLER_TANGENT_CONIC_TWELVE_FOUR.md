# The free-degree pair (12,4) for the sextic tangent-conic scalar

**Theorem SexticTangentConicTwelveFourExclusion.** Let a
characteristic-zero three-variable Keller map of ordinary total
degree at most six have a quartic component
\(R=g(x)+x^2(xz+y^2)\). In the localized coordinates
\(z=v-y^2/x\), the other two components cannot have free
\(y\)-degrees \((12,4)\) and respective leading coefficients
proportional to \(x^{-6},x^{-2}\).

Use the [high-mate restriction](SEXTIC_KELLER_TANGENT_CONIC_HIGH_MATES.md)
and its setup. Work over \(\mathbb C\) by embedding the finitely
generated coefficient field. Subtract \(g(0)\) from \(R\), retain the
notation \(g\) for the resulting polynomial, and set \(g_1=g'(0)\ne0\).
In the localized coordinates \(r=g(x)+x^3v\), the bracket of the
mates is \([P,Q]=D(P)Q_y-P_yD(Q)=\kappa\ne0\), with
\(D=-x^3\partial_x|_r\). The coefficient ring is
\(\mathbb C[r,x,x^{-1}]=\mathbb C[x,x^{-1},v]\).

## Cancellation and the nonvertical component

Scale \(P,Q\) to the indicated leading coefficients and put
\(P_1=P-Q^3\). For any positive free degree \(j\) of \(P_1\)
with leading coefficient \(a\), the leading bracket equation
makes \(a^4/x^{-2j}\) invariant under \(D\). Its \(x\)-valuation at
fixed \(r\) gives \({\rm ord}_x(a)=-j/2\), so \(j\) is even.
In particular the \(y^{11}\)-coefficient vanishes. For the
original \(P\), that coefficient is a scalar times \(x^{-5}\):
in the support grouping of the high-mate restriction the only
total \(y,z\) weight that contributes is \(d=11\), and the
Laurent exponent \(N\) has \(-5\le N\le 6-11=-5\).
Consequently \([y^3]Q\) is a scalar multiple of \(x^{-1}\) by
comparison with \([y^{11}]Q^3\). Center \(Q\) by a translation
of \(y\) by a scalar multiple of \(x\), at fixed \(v\).
This is an affine change in the original coordinates, since
it translates \(z\) by a linear form in \(x,y\).

If \(\deg_y P_1=10\), its leading coefficient is
\(F(r)x^{-5}\) with \(F\in\mathbb C[r]\), by the same leading
equation and the coefficient ring. We claim \(F(0)=0\).
For this claim use the support grouping with no ordinary
degree bound. For an original polynomial of free degree at most
\(2b\), a term \(y^a z^j\) on \(x=0\) lies in the group
\(N=-j,\ d=a+2j\). The polynomial in \(xz/y^2\) for that
group is divisible by \(t^j(1+t)^{\max(0,\lceil(d-2b)/2\rceil)}\)
and has degree at most \(\lfloor d/2\rfloor\). Thus \(j\le b\).
When the coefficient of \(x^{-b}y^{2b}v^0\) is nonzero, the group
\(N=-b,d=2b\) consists of a nonzero multiple of \(z^b\),
so the original \(z\)-degree on \(x=0\) is exactly \(b\).
For the term in question this coefficient is
\([x^0]F(g(x))=F(0)\). If nonzero, \(P_1|_{x=0}\) has
\(z\)-degree five and \(Q|_{x=0}\) has \(z\)-degree two.
They are still a plane automorphism because \(P-Q^3\) is a
polynomial target shear on that plane. A generic line \(y=y_0\)
would give a polynomial embedding of the line with coordinate
degrees five and two, contrary to
[Abhyankar–Moh degree divisibility, García Barroso–Płoski Theorem 6.9](sources/garcia-barroso-ploski-1208.0913v1.pdf).

On the nonvertical component \(H=R/x=0\), the coordinate \(x\)
is a unit since \(H=g_1+xT\), where
\[
               T=y^2+xz+\frac{g(x)-g_1x}{x^2}.
\]
The localized equation is \(r=0\), with free coordinates
\(x,x^{-1},y\). Here \(Q\) retains its centered free degree
four and \(P_1\) has free degree at most nine. Over \(\mathbb C(x)\),
the leading bracket equation again forces its positive degree
to be even. Successively cancel degree eight or four by a
constant multiple of \(Q^2\) or \(Q\), respectively. The
multipliers are constant because \(D\) on this component has
constant field \(\mathbb C\). A degree-two mate would have
leading \(x\)-valuation \(-1\), contrary to the
[quadratic-mate square relation](KELLER_REPEATED_PARALLEL_SCALAR_QUADRATIC_MATE.md#parity-after-centering).
Degree zero cannot give a nonzero constant bracket with a
quartic. The only remaining pair on this component is
\((6,4)\) with leading coefficients proportional to
\(x^{-3},x^{-2}\).

## The centered quartic on the nonvertical component

Set \(x=t^2\), \(w=t^{-1}\) and \(Z=wy\). Rescale the degree-six
mate to be monic, so the [centered sextic-quartic equations
(1)–(3)](SEXTIC_KELLER_TRIPLE_LINE_SEXTIC_QUARTIC_MATES.md#leading-powers-and-centered-equations)
apply with \(H=1\) in their notation and
\(D=-t^5\partial_t/2\). Those coefficient identities hold for any
derivation on the coefficient field. The involution \(t\mapsto-t\)
sends \(Z\) to \(-Z\). The same equations give
\(a=e_0=g_0=0\), since these are anti-invariant constants.
Writing \(Q=Z^4+cZ^2+dZ+h\), the nonsquare reduction gives
\[
        h=c^2/4-2f_0/3,\qquad cd^2\in\mathbb C,\qquad
        D(d^3/4)=\kappa_1/w,\quad \kappa_1\ne0 .
\]
In fact \(d(h-c^2/4+2f_0/3)\) is an anti-invariant constant and
thus zero. Equation (3) forces \(d\ne0\), giving the displayed
relations by (2)–(3).

Integration gives \(d^3=\delta^3 t^{-3}\) for some \(\delta\ne0\);
its additive integration constant vanishes by anti-invariance.
Thus \(d=\delta/t\) and \(c=2\alpha x\) for a constant \(\alpha\).
Substitution gives the following exact restriction of the
original degree-six mate:
\[
     Q|_H=x^{-2}y^4+2\alpha y^2+\delta x^{-1}y+\alpha^2x^2+\beta ,
     \qquad \beta\in\mathbb C .
\]
In the original coordinate ring modulo \(H\), replace
\(x^{-1}\) by \(-T/g_1\). This gives a polynomial representative
of \(Q|_H\) of degree eight whose highest part modulo \(x\) is
\(y^8/g_1^2\), since \(T_2|_{x=0}=y^2\). Such a term cannot
be removed by a multiple of \(H\), whose highest homogeneous
part is divisible by \(x\). The polynomial \(H\) generates
the restriction kernel: \(x\) is already a unit modulo \(H\)
and the localized quotient is the Laurent polynomial ring
\(\mathbb C[x,x^{-1},y]\). This contradicts the original
ordinary degree bound on \(Q\).
