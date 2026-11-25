#!/usr/bin/env python3
"""Exact endpoint checks for the Padé lemma and the large-order appendix.

All arithmetic uses rational intervals. The proof reduces each unbounded
comparison to U=15 by monotonicity, and the parameter comparisons in that
appendix to their stated endpoints. No moment polynomials are sampled.
"""
from fractions import Fraction as Q

SCALE = 10**24


class Interval:
    def __init__(self, lo, hi=None):
        if isinstance(lo, Interval):
            self.lo, self.hi = lo.lo, lo.hi
            return
        a = Q(lo) * SCALE
        b = Q(lo if hi is None else hi) * SCALE
        self.lo = Q(a.numerator // a.denominator, SCALE)
        self.hi = Q(-((-b.numerator) // b.denominator), SCALE)

    def __add__(self, other):
        other = Interval(other)
        return Interval(self.lo + other.lo, self.hi + other.hi)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -Interval(other)

    def __rsub__(self, other):
        return Interval(other) + -self

    def __mul__(self, other):
        other = Interval(other)
        values = [a*b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(values), max(values))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Interval(other)
        assert other.lo * other.hi > 0
        return self * Interval(1/other.hi, 1/other.lo)

    def __rtruediv__(self, other):
        return Interval(other)/self

    def __pow__(self, exponent):
        if isinstance(exponent, int) and exponent >= 0:
            result = Interval(1)
            for _ in range(exponent):
                result = result * self
            return result
        return exp(log(self) * Q(exponent))


def log_reduced(x):
    assert 1 <= x <= 2
    z = (x - 1)/(x + 1)
    total = 2 * sum(z**(2*j+1)/Q(2*j+1) for j in range(20))
    return Interval(total, total + 2*z**41/(41*(1-z*z)))


LOG2 = log_reduced(Q(2))


def log_point(x):
    assert x > 0
    k = 0
    while x > 2:
        x /= 2
        k += 1
    while x < 1:
        x *= 2
        k -= 1
    return log_reduced(x) + k*LOG2


def log(x):
    x = Interval(x)
    return Interval(log_point(x.lo).lo, log_point(x.hi).hi)


def exp_point(x):
    if x < 0:
        return 1 / exp_point(-x)
    k = 0
    while x > 1:
        x /= 2
        k += 1
    total, term = Q(1), Q(1)
    for j in range(1, 31):
        term = term*x/j
        total += term
    result = Interval(total, total + term*x/31/(1-x/32))
    for _ in range(k):
        result = result*result
    return result


def exp(x):
    x = Interval(x)
    return Interval(exp_point(x.lo).lo, exp_point(x.hi).hi)


def atan_inverse(n):
    t = Q(1, n)
    total = sum((-1)**j*t**(2*j+1)/(2*j+1) for j in range(24))
    return Interval(total, total + t**49/49)


PI = 16*atan_inverse(5)-4*atan_inverse(239)
C = Q("1.2762")
U = Q(15)
N = exp(U)
SIGMA = 2*U+2


def above(actual, bound, label="lower bound"):
    value = Interval(actual)
    assert value.lo > Q(bound), (label, float(value.lo), float(value.hi), bound)


def below(actual, bound, label="upper bound"):
    value = Interval(actual)
    assert value.hi < Q(bound), (label, float(value.lo), float(value.hi), bound)


def report_lower(actual, bound, label):
    above(actual, bound, label)
    lower = Interval(actual).lo
    scaled = lower * 10**8
    integer = scaled.numerator // scaled.denominator
    whole, decimals = divmod(integer, 10**8)
    print(f"{label} >= {whole}.{decimals:08d} (required > {bound})")


def maximum(a, b):
    a, b = Interval(a), Interval(b)
    return Interval(max(a.lo,b.lo), max(a.hi,b.hi))


def minimum(a, b):
    a, b = Interval(a), Interval(b)
    return Interval(min(a.lo,b.lo), min(a.hi,b.hi))


def check_core():
    for alpha, bound in [(".03", ".128"), (".1", ".134"),
                         (Q(1,3), ".160"), (".5", ".176"), (".75", Q(7,32))]:
        alpha = Q(alpha)
        kp, km = Q(1,2)+alpha/14, Q(1,2)-alpha/14
        plus = (-log(1-alpha*kp)-alpha*kp)/alpha**2
        minus = (alpha*km-(1-alpha)*log(1+alpha*km/(1-alpha)))/alpha**2
        below(plus, bound, "positive loss")
        below(minus, bound, "negative loss")
    above(2+log(2*PI/199)/2, 0, "Stirling")
    alpha = Q(".03")
    x = Interval(alpha)**Q(1,2)
    # Integrate the mean-value bound on (1-y*y)^(-1/2).
    inv_sine_upper = x + x**3/(6*Interval(1-alpha)**Q(3,2))
    below(2*inv_sine_upper/PI, ".111", "discarded angles")
    below(4*LOG2*x/PI, ".153", "angle loss")
    sin_argument = PI/18
    above((sin_argument-sin_argument**3/6)**2, alpha, "small cosets")
    above((1-Q(".111")-Q(1,2))*LOG2-Q(".153"), 0, "weaker angle bound")

    primes = []
    excesses = []
    prime_sum = Interval(0)
    for n in range(2, 3001):
        if any(n % p == 0 for p in primes if p*p <= n):
            if n == 350:
                below(prime_sum-log(n), "-.499", "prime-sum lower endpoint")
            continue
        primes.append(n)
        log_n = log(n)
        prime_sum += log_n/(n-1)
        excesses.append(prime_sum-log_n)
        below(excesses[-1], ".144", "prime sums")
        if n >= 350:
            below(excesses[-1], "-.499", "prime sums from 350")
    for n, value in zip(primes, excesses):
        if n != 3:
            above(excesses[1]-value, 0, "maximum prime-sum excess")
    above(N, 3269017)
    below(N, 3269018)
    above(exp(Q(".24")*U)-SIGMA-1/N**3, 0, "upper bound on V")
    above(Q(".24")*exp(Q(".24")*U)-2, 0, "upper bound on V derivative")


def check_parameter_reduction():
    above(log(1 - 1/N), "-.001", "log order loss")
    above(-log(Q(".03")), "3.5", "small-angle t threshold")
    below(log(4/Q(".97")), "1.417", "weak angular coefficient")
    below(Q(".001") + Q(2,3)*(Q("1.01") + Q(".128")*Q(".27"))
          - log(Q(".97"))/3, ".71", "left-half-plane loss")
    below(Q("1.01")/25 + 2*Q(".128")*Q(".03"), ".049", "proper-coset baseline per rotation")
    below(Q(".111")*Q(".256"), ".0285", "proper-coset angular slope loss")
    below(Q(".001") + Q(".96")*Q("1.417") + Q(".049") + Q(".111")*Q("1.266"),
          "1.608", "proper-coset angular constant")
    below(Q(".211")*Q(".256"), ".055", "small b angular slope loss")
    below(Q(".001") + Q(".9")*Q("1.417") + Q(".211")*Q("1.266"),
          "1.70", "small b angular constant")
    below(Q(".001") + Q(".98")*Q("1.417") + Q(".025") + Q(".111")*Q("2.459"),
          "1.72", "small B angular constant")
    below(Q(".001") + Q(".98")*Q("1.417") + Q(".025") + Q(".19"),
          "1.64", "small delta angular constant")
    below(Q(".001") + Q("1.39") + Q(".009") + Q(".19"),
          "1.60", "final angular constant")
    below(1/(exp(8)-1), ".001", "large-prime sum remainder")
    below(Q("3.001")*(1+C/log(8)), "4.848", "proper-coset digit ratio")
    below(Q("4.848")/2, "2.43")
    below(Q("4.848")*(50/N)**Q(1,3), ".121")
    for ell, kappa in [(".03",".134"), (".1",".160"), (Q(1,3),Q(7,32))]:
        ell, kappa = Q(ell), Q(kappa)
        surplus = (U+log(4*ell))/2-Q("1.011")-kappa*SIGMA+Q(".49") \
                  -Q("3.6")/(4*ell*N)**Q(1,2)
        above(surplus, ".12", "proper coset at moderate alpha")
    smallest_x = (3*N/4)**Q(1,2)
    below(Q("3.001")*(1+C/log(smallest_x)), "3.6")
    X = exp(2*U/3)/Interval(4)**Q(1,3)
    below(Q("4.001")*(1+C/log(X)), "4.6")
    above(Q(".744")*U-Q("1.267")-log(X)+Q(".54")-Q("4.6")*X/N, ".85")
    upper = maximum(8, (2*U+log(50))/3)-Q(".54")+Q(".121")
    above(U-Q("1.22")-upper, 0)
    t = Q("3.5")
    above(Q(".9715")*U+Q(".96")*t-Q("1.608")
          -maximum(8,(U+2*t+log(SIGMA))/3)-Q("1.89"), 0)

    below(SIGMA/Q(".03")**2/N, Q(1,4))
    below(Q("2.001")*(1+C/(U/2)), "2.35")
    below(Q("2.001")*Q("1.03")*(1+C/(U/2)), "2.42")
    below(log(Q("1.03"))-Q(".54")+Q("2.42"), "1.91", "small b norm charge")
    below(1/(2*N*log(8)), ".001", "V rounding in prime sum ratios")
    above((3*N)**Q(1,2), exp(8).hi, "log square-root cutoff for b >= 3")
    above(Q("10.08")*(8*U/3+Q("1.62"))*exp((U+Q("3.5"))/3),
          exp(8).hi, "small B norm split threshold")
    above(Interval(102400)**Q(1,4)*(2*U+log(51))*exp(U/2),
          exp(8).hi, "small d norm split threshold")
    for ell,u,kappa,lower,upper in [
            (".03",".1",".134","13.07","12.45"),
            (".1",".5",".176","11.54","11.08"),
            (".5",".75",Q(7,32),"10.09","9.41"),
    ]:
        ell,u,kappa = Q(ell),Q(u),Q(kappa)
        b = SIGMA/ell**2
        right = (U+log(b))/2-Q(".54")+Q("2.35")*(b/N)**Q(1,2)
        if u == Q(".5"):
            left = U-Q("1.011")-Q(13,30)*kappa*SIGMA
        else:
            left = U-Q(".001")-Q(".55")*(Q("1.01")+kappa*SIGMA) \
                +Q(".45")*log((1-u)/u)
        above(left, lower, "full-coset lower")
        below(right, upper, "full-coset upper")
    above(U-Q("1.011")-Q(7,32)*Q("6.75")
          -(U+log(300))/2+Q(".54")-Q("2.35")*(300/N)**Q(1,2), 0)
    above(U-Q("1.011")-Q(".176")/2
          -(U+LOG2)/2+Q(".49")-Q("2.35")*(2/N)**Q(1,2), 0)
    above(Q(".395")*U-Q(".45")*log(SIGMA)-Q("3.61"), 0)
    above(U/3-log(9)/2-Q("2.62"), 0)
    above(U/2-LOG2/2-Q("1.195")-Q("1.91"), 0)
    below(Q(".001")+log(2)+Q(".5"), "1.195", "adjacent-power gain")

    below(Interval(1024)**Q(1,3), "10.08")
    below(log(128)/3+Q(".001"), "1.62")
    below(Q("1.01")+Q(".128")*8*Interval(2)**Q(1,2), "2.459")
    below(Q("1.01")/50+Q(".128")*Q(".03"), ".025")
    vv = 8*U/3+Q("1.62")
    z = U-log(16)
    below(vv*(1+C/z)/(16*z), ".236")
    below(vv**2/(U*z), "9.46")
    for t, bound in [(Q("3.5"),"3.59"), (U,"1.60")]:
        gap = 2*U/3+(Q(".98")-Q(1,3))*t-log(Q("10.08")*vv) \
              -Q("2.417")-Q("95.3568")*U*exp((t-2*U)/3)
        above(gap, bound, "small B")
    n_term = Interval(102400)**Q(1,4)
    below(Q(".236")+Q("9.46")*U*n_term*exp(-U/2), "1.65")
    above(4*U/3-log(50)/3-Q(".71")
          -U/2-log(102400)/4-log(2*U+log(51))-Q(".461")-Q("1.65"), 0)
    above(log(Q(".97")/Q("1.03"))+Q(".54"), ".47")
    below(Q("2.06")*(1+C/U), "2.24")


def check_degree_and_lattice():
    W = exp(U/2)/Interval(SIGMA)**Q(1,2)
    above(W, 308)
    above(log(Q(307,309))-1-LOG2, "-1.70")
    above(log(Q(307,309))-1, "-1.01")
    above(log(Q(307,4))-(1+log(309))/2, ".97")
    delta_at_top = log((U+Q(".47"))/Q("2.24"))+log(1-1/W) \
                   -2*(U+Q(".47"))/Q("2.24")*(1+log(W+1))/W
    above(delta_at_top, "1.6")
    above(LOG2+log(2048)/4, "2.59")
    above(4*Q("2.59")-Q(".001"), 4*C, "prime-count slack for V")
    above((U+Q(".47"))/Q("2.24"), 2, "spacing of progression valuation jumps")
    above(Interval(2048)**Q(1,4)*exp(U/2)*2*U,
          exp(8).hi, "norm split threshold after B lower bound")
    below(6*(Q("3.24")*U+1)/(U+Q(".95")-log(Q("2.5"))-10-Q("3.24")*U/8192), 61)
    above((U+Q(".47"))/Q("2.24"), "6.9")
    below(1-log(6)+log(1+Q("1.5")/Q("6.9")), "-.595")
    t_min = (U-log(SIGMA))/2
    x = exp(-t_min/2)
    below((2/PI)*x/(1-x*x)**Q(1,2)*(Q("1.266")+Q(".256")*U), ".19")
    j0 = Q("6.65")+6*Q("3.24")*U/16384
    delta = Q(".045")
    above((Q(".98")-9*delta)*Q(15,32)-1/U-delta*6*Q("3.24")/16384,
          0, "small delta surplus derivative")
    above((Q(".98")-9*delta)*t_min-log(delta)-log(Q("3.24")*U)
          -Q("2.101")-j0*delta, ".10")
    below(3000/N*(9*U+j0), ".14")
    above(U+(U-log(50))/3-Q(".71")-log(3000*Q("3.24")*U)-Q(".601"), 0)
    below((U+Q("1.266")+Q(".256")*U)*60/(Q(".045")*N), ".009")
    above(log(1-(SIGMA/N)**Q(1,2))-log(4), "-1.39")

    eta = Interval(Q(9,8))
    below((N*(N-1))**Q(-1,6), ".007")
    below((2/(N-1))**Q(1,3), ".01")
    below(Q(".264")/(1-eta**Q(-1,3)), 7)
    below(log(eta)/(1-eta**Q(-1,3)), "3.06")
    coefficient = Q(".132")*eta**2*Interval(Q(2,3))**Q(1,2) \
                  *Interval(Q(32,7))**Q(2,9) \
                  *(1/(1-eta**Q(-2,3))+1/(1-eta**Q(-7,3)))
    below(coefficient, "3.4")
    above(log(eta)+log(Q(32,7)*N*(N-1))/3, Q(3,7),
          "decrease of lattice log-weighted summand")
    below(log(eta)+log(Q(32,7))/3
          +log(eta)*eta**Q(-7,3)/(1-eta**Q(-7,3)), 1)


def wt(x, w):
    return exp(-U/3)*(14*(U+log(Q("2.5")*x)+Q("3.06"))/((5*w)**Q(1,3)*x)
        +Q("6.8")*(U+log(x/(2*w))/3+1)/((2*w)**Q(1,18)*x**Q(4,9)))


def v0(u, w):
    return 3*U+log(u/(2*w))+Q(".001")


def norm_lower(ell, u, w):
    b = Interval(u/(2*w))
    zz = sum((b**Q(1,j))*exp((Q(3,j)-1)*U)/ell for j in range(4,9))
    zz += 3*U/(ell*LOG2)*b**Q(1,9)*exp(-2*U/3)
    coefficient = Q("1.017")*(5+minimum(Q("2.5"),Interval(2*w*ell**2)**Q(-1,3))+zz)
    return U-Q("1.01")-log(Q("2.5"))-coefficient-v0(u,w)/8192


def check_last_comparison():
    log_p = U+log(Q(".09"))
    numerator = U-Q("1.01")-log(Q("2.5")) \
        -Q("2.5")*Q("3.24")*U/log_p*(1+C/log_p)-Q("3.24")*U/8192
    above(numerator/(Q("3.24")*U), ".0502", "first multiplier bound")
    w0 = Q(6,5)
    ell,u = Q(".21"),Q("1.4")
    report_lower(norm_lower(ell,u,w0)-wt(Interval(ell),Interval(w0)),
                 ".074", "delta cutoff surplus")
    above(norm_lower(Q(".045"),Q(".21"),w0), 0, "positive outer norm lower bound")
    below(Q(2,25)*v0(Q(".21"),w0)/norm_lower(Q(".045"),Q(".21"),w0), ".675")
    x_min, w_max = Q(".045"), Q(9)
    above(U+log(Q("2.5")*x_min)+Q("3.06"), 14, "first lattice numerator")
    above(U+log(x_min/(2*w_max))/3+1, 14, "second lattice numerator")
    above(U+log(Q("2.5")*x_min), 8, "logarithm in prime-count derivative")
    above(15*x_min*N, exp(8).hi, "large-prime split starts above e^8")
    above(2*x_min*N, 7481, "prime-sum lower bound applies to P")
    below(1+log(Q(5,6)+1/(6*N))-LOG2, ".127", "equal-height loss")
    below(Q(3,40), ".08", "first norm logarithmic derivative")
    below((1+Q(3,8192))/(U+Q(".127")), ".08", "second norm logarithmic derivative")
    above((C-Q(".02"))*8-2*Q(".02")*C, 0, "decrease of rational prime-count factor")
    below(Q(".21")*(2+Q(3,4096)), ".43", "digit loss derivative")
    for ell,u,w,J,test,margin in [
        (".045",".10",None,15,None,".36"),
        (".10",".14","9","1.26",".109","1.63"),
        (".14",".18","3.6",".648",".117",".83"),
        (".18",".21","1.8",".378",".093",".30"),
    ]:
        ell,u,J = Q(ell),Q(u),Q(J)
        if w is not None:
            w = Q(w)
            report_lower(norm_lower(ell,u,w)-wt(Interval(ell),Interval(w)),
                         test, f"v cutoff surplus with delta in [{ell}, {u}]")
        av = v0(u,w0)
        above(av, 40, "V0 bound on final row")
        threshold = minimum(av,U+Q(".127")+wt(Interval(ell),Interval(w0))+av/8192)
        z = U+log(Q("2.5")*ell)
        extra = maximum(0,Q("1.017")*Interval(u/(2*w0))**Q(1,3)-Q("1.96")*ell)
        loss = u*(Q("2.5")*av/z*(1+C/z)+2*(U+Q(".127"))+av/4096) \
            +2*u*wt(Interval(u),Interval(w0))+extra
        report_lower(U-log(J*threshold)-Q("2.061")-loss, margin,
                     f"final Padé surplus with delta in [{ell}, {u}]")
        below(log(Q("2.5")*ell)-(log(u/Q("2.4"))+Q(".001"))/3, ".02")


if __name__ == "__main__":
    if not __debug__:
        raise SystemExit("Verification requires assertions. Run Python without -O.")
    check_core()
    check_parameter_reduction()
    check_degree_and_lattice()
    check_last_comparison()
    print("Large-order Padé endpoint comparisons verified by exact rational bounds.")

