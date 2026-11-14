"""Solver-free exact replay primitives for polynomial trace certificates."""
from flint import fmpq, fmpq_poly, fmpq_mpoly_ctx

CTX=fmpq_mpoly_ctx.get(('U','V','a'))
ZERO=fmpq_poly([])
ONE=fmpq_poly([1])


def exact_div(a,b):
    quotient,remainder=divmod(a,b)
    assert not remainder
    return quotient


def lcm(a,b):
    result=exact_div(a,a.gcd(b))*b
    return result/result.leading_coefficient()


def embed(poly):
    return CTX.from_dict({(0,0,i):co for i,co in enumerate(poly.coeffs()) if co})


def encode_expression(expr,U,V,parameter):
    """Used by discovery only; the replay decodes rational sparse rows."""
    import sympy as s
    domain=s.QQ.frac_field(parameter)
    data=s.Poly(expr,U,V,domain=domain).rep.to_dict()

    def univariate(poly):
        if not poly:
            return ZERO
        return fmpq_poly([fmpq(str(poly.get((i,),0))) for i in range(poly.degree()+1)])

    coefficients={m:(univariate(co.numer),univariate(co.denom)) for m,co in data.items()}
    denominator=ONE
    for num,den in coefficients.values():
        denominator=lcm(denominator,den)
    numerator={}
    for m,(num,den) in coefficients.items():
        for exponent,coefficient in enumerate((num*exact_div(denominator,den)).coeffs()):
            if coefficient:
                numerator[m+(exponent,)]=coefficient
    return dict(numerator=[[list(m),str(co)] for m,co in sorted(numerator.items())],
                denominator=[str(co) for co in denominator.coeffs()])


def decode(rows):
    monomials=[tuple(m) for m,co in rows['numerator']]
    assert len(monomials)==len(set(monomials))
    assert all(len(m)==3 and all(type(e) is int and e>=0 for e in m) for m in monomials)
    numerator=CTX.from_dict({tuple(m):fmpq(co) for m,co in rows['numerator']})
    denominator=fmpq_poly([fmpq(co) for co in rows['denominator']])
    assert denominator
    return numerator,denominator


def constant(value):
    return CTX.constant(value),ONE


def monomial(i,j):
    return CTX.from_dict({(i,j,0):fmpq(1)}),ONE


def product(terms):
    numerator=CTX.constant(1)
    denominator=ONE
    for num,den in terms:
        numerator*=num
        denominator*=den
    return numerator,denominator


def sum_products_zero(terms):
    """Check sum sign*product(factors)=0 after exact denominator clearing."""
    products=[(sign,product(factors)) for sign,factors in terms]
    denominator=ONE
    for sign,(num,den) in products:
        denominator=lcm(denominator,den)
    total=CTX.constant(0)
    for sign,(num,den) in products:
        total+=sign*num*embed(exact_div(denominator,den))
    return not total


def leading_UV(pair):
    numerator,denominator=pair
    exponents=[(i,j) for i,j,k in numerator.to_dict()]
    # Singular dp is degree reverse lexicographic, U before V.
    return max(exponents,key=lambda m:(sum(m),-m[1]))


def univariate(pair):
    numerator,denominator=pair
    values=numerator.to_dict()
    assert all(i==j==0 for i,j,k in values)
    degree=max((k for i,j,k in values),default=-1)
    return fmpq_poly([values.get((0,0,k),fmpq(0)) for k in range(degree+1)]),denominator


def trace_powers(entries,size):
    """Return denominator-cleared traces of powers1,2,3; no CAS solver."""
    denominator=ONE
    for pair in entries.values():
        numerator,den=univariate(pair)
        denominator=lcm(denominator,den)
    matrix=[[ZERO for j in range(size)] for i in range(size)]
    for (i,j),pair in entries.items():
        num,den=univariate(pair)
        matrix[i][j]=num*exact_div(denominator,den)
    tr1=sum((matrix[i][i] for i in range(size)),ZERO)
    square=[[sum((matrix[i][k]*matrix[k][j] for k in range(size)
                   if matrix[i][k] and matrix[k][j]),ZERO)
             for j in range(size)] for i in range(size)]
    tr2=sum((square[i][i] for i in range(size)),ZERO)
    tr3=sum((square[i][j]*matrix[j][i] for i in range(size) for j in range(size)
             if square[i][j] and matrix[j][i]),ZERO)
    return [(tr1,denominator),(tr2,denominator**2),(tr3,denominator**3)]


def same_univariate(left,right):
    return left[0]*right[1]==right[0]*left[1]
