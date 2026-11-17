"""Exact finite certificates for the homogeneous radial-pencil theorem.

Mathematical scopes and proofs are in SIMPLEX_PENCILS_AND_AFFINE_LIMITS.tex.
The executable input representation uses exact algebraic coefficients and
positive integer Dirichlet shapes. Outside-class results make no FC decision.
"""
from collections import Counter
import sympy as s


def exact_nodes(coefficients, shapes=None):
    coefficients = [s.sympify(c) for c in coefficients]
    if not coefficients:
        raise ValueError('at least one node required')
    if any(c.has(s.Float) or c.is_algebraic is not True for c in coefficients):
        raise ValueError('this certificate requires exact algebraic nodes')
    shapes = [1]*len(coefficients) if shapes is None else list(shapes)
    if len(shapes) != len(coefficients) or any(type(a) is not int or a <= 0 for a in shapes):
        raise ValueError('positive integer shapes required')
    field = s.QQ if all(c.is_Rational for c in coefficients) else s.QQ.algebraic_field(*coefficients)
    grouped = Counter()
    for c,a in zip(coefficients,shapes):
        grouped[field.from_sympy(c)] += a
    return [(field.to_sympy(c),m) for c,m in grouped.items()], sum(shapes)


def chain_representation(coefficients, shapes=None, normalized=True):
    """Return polynomial path weights, or the separately typed point case.

    Every path starts at zero; a chain entry means integral_0^node w(z)F(z)dz.
    Normalized=True gives Dirichlet probability expectation. False gives the
    coordinate-volume functional on the expanded uniform simplex.
    """
    grouped,total = exact_nodes(coefficients,shapes)
    scale = s.factorial(total-1) if normalized else s.Integer(1)
    if len(grouped)==1:
        return dict(kind='point',node=grouped[0][0],mass=scale/s.factorial(total-1),shape_sum=total)
    z = s.Symbol('chain_z')
    entries = []
    for node,multiplicity in grouped:
        regular = 1/s.prod((z-other)**m for other,m in grouped if other != node)
        weight = 0
        for j in range(1,multiplicity+1):
            order = multiplicity-j
            partial = s.diff(regular,z,order).subs(z,node)/s.factorial(order)
            weight += partial*(node-z)**(total-j-1)/(s.factorial(j-1)*s.factorial(total-j-1))
        entries.append((node,s.expand(scale*weight)))
    return dict(kind='chains',variable=z,entries=entries,shape_sum=total)


def integrate_polynomial(rep, polynomial, variable):
    if rep['kind']=='point':
        return s.expand(rep['mass']*s.sympify(polynomial).subs(variable,rep['node']))
    z = rep['variable']
    value = 0
    for node,weight in rep['entries']:
        primitive = s.Poly(s.expand(weight*s.sympify(polynomial).subs(variable,z)),z,extension=True).integrate()
        value += primitive.eval(node)-primitive.eval(0)
    return s.simplify(value)


def exponential_certificate(coefficients, shapes=None):
    """Group an affine exponential period by algebraic exponents.

    Nonvanishing follows from the written pole-budget proof and
    Lindemann--Weierstrass, not numerical evaluation of the exponential sum.
    Individual grouped coefficients are permitted to be zero.
    """
    rep = chain_representation(coefficients,shapes)
    if rep['kind']=='point':
        groups = {rep['node']:rep['mass']}
    else:
        groups = {}
        z = rep['variable']
        for node,weight in rep['entries']:
            p = s.Poly(weight,z,extension=True)
            primitive = s.expand(sum((-1)**j*s.diff(weight,z,j) for j in range(p.degree()+1)))
            if s.simplify(s.diff(primitive,z)+primitive-weight) != 0:
                raise ArithmeticError('bad exponential primitive')
            groups[node] = groups.get(node,0)+primitive.subs(z,node)
            groups[s.Integer(0)] = groups.get(s.Integer(0),0)-primitive.subs(z,0)
    groups = {node:s.simplify(c) for node,c in groups.items() if s.simplify(c) != 0}
    if not groups:
        raise ArithmeticError('contradicts the all-degree pole-budget theorem')
    return dict(status='nonzero_by_written_theorem',groups=groups,shape_sum=rep['shape_sum'],
                theorem='FIRAFFINESIMPLEX20261116')


def homogeneous_pencil_certificate(polynomial, variables, coefficients, shapes=None):
    """Recognize F(Y)=H(sum Y,sum c_iY_i), with homogeneous H.

    The returned finite identity can be rechecked independently. A successful
    nonzero certificate invokes the written all-degree theorem; it does not
    claim an explicit first surviving exponent or a Strong Factorial bound.
    """
    variables = tuple(variables)
    if not variables or len(coefficients) != len(variables):
        raise ValueError('one node per variable required')
    exact_nodes(coefficients,shapes)
    coefficients = [s.simplify(s.sympify(c)) for c in coefficients]
    polynomial = s.sympify(polynomial)
    if polynomial.has(s.Float) or polynomial.free_symbols-set(variables):
        raise ValueError('exact algebraic phase coefficients required')
    p = s.Poly(polynomial,*variables,extension=True)
    if any(c.is_algebraic is not True for _,c in p.terms()):
        raise ValueError('exact algebraic phase coefficients required')
    if p.is_zero:
        return dict(status='zero',degree=0,binary_terms=[],phase_terms=[],theorem='FIRPENCILGAMMA20261116')
    degrees = {sum(es) for es,_ in p.terms()}
    if len(degrees) != 1:
        return dict(status='outside_class',reason='phase is not homogeneous')
    degree = degrees.pop()
    S,L = s.Dummy('pencil_S'),s.Dummy('pencil_L')
    pivot = next((j for j in range(len(variables)-1)
                  if s.simplify(coefficients[j]-coefficients[-1]) != 0),None)
    if pivot is None:
        c = p.as_expr().subs({v:0 for v in variables[:-1]}).subs(variables[-1],1)
        binary = c*S**degree
        if s.expand(p.as_expr()-c*sum(variables)**degree) != 0:
            return dict(status='outside_class',reason='coincident nodes span only the radial sum')
    else:
        remaining = [j for j in range(len(variables)-1) if j != pivot]
        pivot_value = (L-coefficients[-1]*S-sum((coefficients[j]-coefficients[-1])*variables[j] for j in remaining))/(coefficients[pivot]-coefficients[-1])
        replacements = {variables[pivot]:pivot_value,
                        variables[-1]:S-pivot_value-sum(variables[j] for j in remaining)}
        binary = s.expand(p.as_expr().subs(replacements,simultaneous=True))
        binary = s.Poly(binary,S,L,*[variables[j] for j in remaining],extension=True).as_expr()
        if binary.free_symbols & {variables[j] for j in remaining}:
            return dict(status='outside_class',reason='phase needs additional independent coordinates')
    terms = s.Poly(binary,S,L,extension=True).terms()
    if any(sum(es) != degree for es,_ in terms):
        raise ArithmeticError('coordinate conversion lost homogeneity')
    reconstructed = binary.subs({S:sum(variables),L:sum(c*y for c,y in zip(coefficients,variables))},simultaneous=True)
    if not s.Poly(s.expand(reconstructed-p.as_expr()),*variables,extension=True).is_zero:
        raise ArithmeticError('pencil membership identity failed')
    return dict(status='nonzero_infinitely_many_pure_moments',degree=degree,
                coefficients=[str(c) for c in coefficients],
                shapes=[1]*len(variables) if shapes is None else list(shapes),
                phase_terms=[[list(es),str(c)] for es,c in p.terms()],
                binary_terms=[[list(es),str(c)] for es,c in terms],
                theorem='FIRPENCILGAMMA20261116',
                assures_strong_factorial=False,assures_first_nonzero_index=False)
