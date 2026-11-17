"""Exact Gaussian square-product/Kloosterman and symmetric-power constructors.

Galois groups use the explicitly cited classical theorem, not finite jet tests.
The full marked-image and symmetric-power faithfulness proofs are in the note.
"""
import math
import sympy as s

t, theta = s.symbols('t theta')


def positive_integer(k):
    if type(k) is not int or k < 1: raise ValueError('positive integer required')


def family(k, scale=1, parameterized=False):
    positive_integer(k); scale = s.Rational(scale)
    if not scale: raise ValueError('nonzero square multiplier required')
    # A positive rational scale is a rational square only when its root is
    # rational, but all nonzero rational moments remain units for any scale.
    n = k+3; terms = []
    for spatial,epower,c in [((0,1,0),0,1),((1,1,0),1,1),
                             ((0,0,2),1,-1),((1,0,2),2,s.Rational(-3,2)),((2,0,2),3,s.Rational(-1,2))]:
        terms.append([list(spatial)+(k*[2])+[epower if parameterized else 0],str(scale*c)])
    H = s.diag(s.Matrix([[0,1],[1,0]]),s.eye(k+1))
    return dict(dimension=n,parameter_count=1,phase_terms=terms,
        parameter_ideal=[] if parameterized else [[[[1],'1']]],
        precision_matrix=[[str(v) for v in row] for row in H.tolist()])


def independent_sum_family(k,copies,scale=1,parameterized=False):
    positive_integer(copies); data = family(k,scale,parameterized); n = data['dimension']
    terms = []
    for j in range(copies):
        for e,c in data['phase_terms']:
            spatial = [0]*(n*copies);spatial[j*n:(j+1)*n] = e[:n]
            terms.append([spatial+[e[-1]],c])
    H = s.diag(*[s.Matrix(data['precision_matrix']) for _ in range(copies)])
    return dict(data,dimension=n*copies,phase_terms=terms,
                precision_matrix=[[str(v) for v in row] for row in H.tolist()])


def double_factorial_even_moment(m):
    if type(m) is not int or m < 0: raise ValueError('nonnegative moment index required')
    return math.factorial(2*m)//(2**m*math.factorial(m))


def coefficients(k,N,scale=1):
    positive_integer(k)
    if type(N) is not int or N < 0: raise ValueError('nonnegative jet order required')
    scale = s.Rational(scale)
    return [s.Integer(double_factorial_even_moment(m))**k*scale**m for m in range(N+1)]


def quotient_operator(k,scale=1):
    positive_integer(k); scale = s.Rational(scale)
    if not scale: raise ValueError('nonzero pole scale required')
    return s.expand(1-scale*t*(2*theta+1)**k)


def connection(k,scale=1):
    """Rows give theta*basis_i in (1,h,theta*h,...,theta^(k-1)*h)."""
    positive_integer(k); scale = s.Rational(scale)
    if not scale: raise ValueError('nonzero pole scale required')
    M = s.zeros(k+1)
    for i in range(1,k): M[i,i+1] = 1
    M[k,0] = -1/(scale*2**k*t)
    for j in range(k): M[k,j+1] = -s.Rational(math.comb(k,j)*2**j,2**k)
    M[k,1] += 1/(scale*2**k*t)
    return M


def base_marker(a,b,c):
    """Actual B period of Z^a W^b T^(2c), via its full momentum source."""
    if any(type(j) is not int or j < 0 for j in [a,b,c]): raise ValueError('nonnegative marker powers required')
    beta = s.Symbol('beta')
    binomial = sum(s.binomial(-s.Rational(1,2)-c,j)*(t*(3-t)*beta+t*beta**2)**j for j in range(b+1))
    co = s.expand((t+beta)**a*binomial).coeff(beta,b)
    return s.cancel(math.factorial(b)*double_factorial_even_moment(c)*(1-t)**(2*c-a)*co)


def marker_operator(a,b,c,spectator_powers,scale=1):
    """Actual marked period = polynomial(t)+operator(theta)*h_scale.

    T has exponent 2c; spectator_powers are FULL original exponents.
    Each rational source pole is expanded with exact Euler principal parts.
    """
    k = len(spectator_powers); positive_integer(k); scale = s.Rational(scale)
    if any(type(j) is not int or j < 0 for j in spectator_powers): raise ValueError('invalid spectator marker')
    if any(j%2 for j in spectator_powers): return s.S.Zero,s.S.Zero
    rs = [j//2 for j in spectator_powers]
    moment_shift = s.sympify(s.prod(2*theta+2*j+1 for r in rs for j in range(r)))
    f = base_marker(a,b,c); num,den = s.fraction(f)
    polynomial,proper = s.div(num,den,t); proper = s.cancel(proper/den)
    degree = s.degree(den,t); operator = s.S.Zero
    if proper:
        q = s.Symbol('q'); expansion = s.Poly(s.cancel(proper*(1-t)**degree).subs(t,1-q),q)
        for (j,),co in expansion.terms():
            pole = degree-j
            if pole < 1: raise ValueError('unexpected source pole polynomial')
            operator += co*s.prod(theta+i for i in range(1,pole))/math.factorial(pole-1)
    polynomial_out = s.S.Zero
    for (m,),co in s.Poly(polynomial,t).terms():
        polynomial_out += co*(scale*t)**m*double_factorial_even_moment(m)**k*moment_shift.subs(theta,m)
    return s.expand(polynomial_out),s.expand(operator*moment_shift)


def compositions(n,m):
    if n == 1: yield (m,); return
    for j in range(m+1):
        for e in compositions(n-1,m-j): yield (j,)+e


def symmetric_connection(k,m,scale=1):
    positive_integer(k); positive_integer(m)
    M = connection(k,scale); basis = list(compositions(k+1,m)); positions = {e:i for i,e in enumerate(basis)}
    result = s.zeros(len(basis))
    for row,e in enumerate(basis):
        for i,a in enumerate(e):
            if not a: continue
            for j in range(k+1):
                if not M[i,j]: continue
                z = list(e);z[i] -= 1;z[j] += 1
                result[row,positions[tuple(z)]] += a*M[i,j]
    return basis,result


def galois_model(k,m=1):
    positive_integer(k); positive_integer(m)
    if k < 2: raise ValueError('classical non-solvable statement requires k>=2')
    H = 'Sp_'+str(k) if k%2 == 0 else 'SL_'+str(k)+' times central mu_2'
    return dict(rank=math.comb(k+m,m),base_rank=k+1,symmetric_power=m,
        quotient_rank=k,quotient_full_group=H,
        connected_reductive_group='Sp_'+str(k) if k%2 == 0 else 'SL_'+str(k),
        full_group='Ga^'+str(k)+' semidirect ('+H+')',unipotent_radical_dimension=k,
        faithful_symmetric_realization=True,identity_component_solvable=False,
        galois_basis='Katz theorem as explicitly restated by Fresan-Sabbah-Yu1810.06454v5,Proposition2.7 proof;own Kummer twist and nonsplit-extension argument.',
        finite_jets_prove_galois_group=False,independent_replay=False,formal_verification=False,external_review=False)
