"""Exact sharp Gaussian real regions over presented QQ coefficient rings.

Full primary localization is trusted CAS. Sparse Weyl identities are checked
against the input-derived PBW-induced source. Uncapped enumeration terminates
by holonomicity; a finite cap returns UNKNOWN, never a mathematical verdict.
"""
from pathlib import Path
from itertools import count
from fractions import Fraction
import copy, json, math, re, time
import sympy as s
import gaussian_real_cutoff_producer as G

O, W, S = G.O, G.W, G.S


class FieldWeyl(W.Weyl):
    """The sealed PBW multiplication with exact QQ(U) coefficients."""
    def __init__(self, n, parameters):
        super().__init__(n)
        self.parameters = tuple(s.Symbol(x) for x in parameters)
        self.K = s.QQ.frac_field(*self.parameters)
        self.locals = {str(x): x for x in self.parameters}

    def coefficient(self, c):
        if isinstance(c, (float, bool)):
            raise ValueError('exact coefficient required')
        return self.K.from_sympy(s.sympify(str(c), locals=self.locals))

    def scalar(self, c):
        c = self.coefficient(c)
        return {self.zero: c} if c else {}

    def var(self, i):
        e = list(self.zero); e[i] = 1
        return {tuple(e): self.K.one}

    def parse(self, terms):
        out = {}
        for e, c in terms:
            if len(e) != 2*self.n or any(type(j) is not int or j < 0 for j in e):
                raise ValueError('invalid Weyl exponent vector')
            e = tuple(e); out[e] = out.get(e, self.K.zero)+self.coefficient(c)
        return {e: c for e, c in out.items() if c}

    def serialize(self, p):
        return [[list(e), str(self.K.to_sympy(c))] for e, c in sorted(p.items()) if c]


def source(data):
    allowed = {'dimension', 'field_parameters', 'coefficient_variables',
               'coefficient_relations', 'phase_terms', 'precision_matrix'}
    if set(data)-allowed:
        raise ValueError('unknown induced-source input')
    n = data['dimension']; names = data['field_parameters']+data['coefficient_variables']
    if type(n) is not int or n < 1 or len(names) != len(set(names)) or any(
            not isinstance(x, str) or not re.fullmatch(r'c[0-9]+', x) for x in names):
        raise ValueError('invalid source dimensions or coefficient names')
    q = len(data['coefficient_variables']); b = n+q; A = FieldWeyl(b+1, data['field_parameters'])
    vs = [A.var(i) for i in range(2*A.n)]
    u, y, v = vs[:n], vs[n:b], vs[b]
    du, dy, dv = vs[A.n:A.n+n], vs[A.n+n:A.n+b], vs[-1]
    H = s.Matrix(data['precision_matrix'])
    if H.shape != (n,n) or H != H.T or not H.det() or any(not c.is_Rational for c in H):
        raise ValueError('fixed nonsingular rational precision required')
    phase = data['phase_terms']; f = A.substitute(phase, du+y); gens = []
    for i in range(n):
        derivative = []
        for e, c in phase:
            if len(e) != b or any(type(j) is not int or j < 0 for j in e):
                raise ValueError('invalid phase exponent')
            if e[i]:
                z = list(e); z[i] -= 1
                derivative.append([z, str(A.K.to_sympy(A.coefficient(c)*e[i]))])
        gens.append(A.add(A.scale(u[i], -1),
            *[A.scale(du[j], A.coefficient(H[i,j])) for j in range(n)],
            A.mul(A.substitute(derivative, du+y), dv)))
    gens.extend(A.substitute(r, y) for r in data['coefficient_relations'])
    gens.append(A.sub(v, f))
    ES = A.scale(A.mul(dv,v), -1)
    E = A.add(*[A.mul(x,d) for x,d in zip(u,du)])
    return A, gens, u+y+du+dy, ES, E, f


def relation(columns, A, budget):
    """Exact dependence with a nonzero scalar-polynomial coordinate."""
    basis = {}
    for index, (tag, column) in enumerate(columns):
        budget.check(); row = dict(column); comb = {index: A.K.one}
        while row:
            budget.check(); lead = max(row)
            if lead not in basis: break
            b, c = basis[lead]; scale = row[lead]
            for target, values in [(row,b),(comb,c)]:
                for k,z in values.items():
                    target[k] = target.get(k,A.K.zero)-scale*z
                    if not target[k]: del target[k]
        if not row:
            if any(columns[j][0][0] == 'p' and c for j,c in comb.items()): return comb
            continue
        z = row[max(row)]
        basis[max(row)] = ({k:v/z for k,v in row.items()}, {k:v/z for k,v in comb.items()})
    return None


def search_identity(data, kind, budget):
    if kind not in ('ordinary','Euler'): raise ValueError('unknown certificate kind')
    A, gens, base, ES, E, f = source(data); n = data['dimension']; b = len(base)//2
    for box in count(1):
        columns = []
        for exponent in G.monomials(2*A.n,box):
            budget.check(); term = {exponent:A.K.one}
            for i,g in enumerate(gens):
                columns.append((('w',i,exponent), A.scale(A.mul(term,g),-1)))
        for degree in range(box+1):
            for exponent in G.compositions(2*b,degree):
                budget.check()
                if kind == 'Euler' and sum(exponent[:n])-sum(exponent[b:b+n]) < 1: continue
                term = A.substitute([[exponent,'1']],base)
                for i in range(box-degree+1):
                    value = A.mul(term,A.power(ES,i))
                    if kind == 'ordinary': value = A.mul(value,f)
                    columns.append((('r',exponent,i), A.scale(value,-1)))
        for i in range(box+1):
            for j in range(box-i+1 if kind == 'Euler' else 1):
                columns.append((('p',i,j), A.mul(A.power(ES,i),A.power(E,j))))
        rel = relation(columns,A,budget)
        if rel is None: continue
        scalar, operator = {}, {}; witnesses = [{} for _ in gens]
        for index,c in rel.items():
            tag = columns[index][0]
            if tag[0] == 'p': scalar[(tag[1],tag[2])] = c
            elif tag[0] == 'r': operator[tag[1]+(tag[2],)] = c
            else: witnesses[tag[1]] = A.add(witnesses[tag[1]],{tag[2]:c})
        encode = lambda d:[[list(e),str(A.K.to_sympy(c))] for e,c in sorted(d.items()) if c]
        return dict(kind=kind, method='exhaustive QQ(U) graph search', box=box,
            scalar_terms=encode(scalar), operator_terms=encode(operator),
            left_witnesses=[A.serialize(w) for w in witnesses])


def verify(data, certificate):
    A,gens,base,ES,E,f = source(data); n = data['dimension']; b = len(base)//2
    for kind in ('ordinary','Euler'):
        c = certificate[kind]; p = A.substitute(c['scalar_terms'],[ES,E])
        R = A.substitute(c['operator_terms'],base+[ES])
        if not p: raise ValueError('zero scalar polynomial')
        if kind == 'ordinary' and any(e[1] for e,co in c['scalar_terms'] if A.coefficient(co)):
            raise ValueError('ordinary scalar depends on Euler')
        if kind == 'Euler' and any(sum(e[:n])-sum(e[b:b+n]) < 1
                for e,co in c['operator_terms'] if A.coefficient(co)):
            raise ValueError('Euler remainder does not strictly raise u-order')
        if len(c['left_witnesses']) != len(gens): raise ValueError('wrong witness length')
        total = A.add(*[A.mul(A.parse(w),g) for w,g in zip(c['left_witnesses'],gens)])
        if total != A.sub(p,A.mul(R,f) if kind == 'ordinary' else R):
            raise ValueError('induced-source graph identity failed')
    if certificate.get('graph_generators') != [A.serialize(g) for g in gens]:
        raise ValueError('changed original graph generators')
    if certificate.get('specialization_start') != 0: raise ValueError('specialization must start at zero')
    return True


def source_input(data, algebra):
    n,r,phase,ideal,H = O.F.normalize(data); xs = s.symbols('x:'+str(n))
    ys = algebra.dependent
    K = s.QQ.frac_field(*algebra.independent)
    P = sum(c*s.prod(z**i for z,i in zip(xs+algebra.cs,e)) for e,c in phase.items())
    encode = lambda f,vs:[[list(e),str(K.to_sympy(c))] for e,c in s.Poly(f,*vs,domain=K).terms() if c]
    # When there are no dependent variables, a scalar ordinary relation is zero.
    rels = [encode(O.polynomial(g,algebra.cs),ys) for g in algebra.component['primary_ideal']] if ys else []
    return dict(dimension=n,field_parameters=[str(c) for c in algebra.independent],
        coefficient_variables=[str(c) for c in ys],coefficient_relations=rels,
        phase_terms=encode(P,xs+ys),precision_matrix=[[str(c) for c in row] for row in H.tolist()])


def membership(target, relations, variables, A, budget):
    """Return exact polynomial LEFT coefficients in a commutative ideal.

    This search is used only after membership is certified by the finite table.
    It retains a witness in the original, localized primary relations.
    """
    if target == 0: return [s.S.Zero for _ in relations]
    def poly(f):
        if not variables: return {():A.coefficient(f)} if f else {}
        return {e:A.coefficient(c) for e,c in s.Poly(f,*variables,domain=A.K).terms() if c}
    target_col = poly(target)
    for box in count(0):
        cols = []; tags = []
        for exp in G.monomials(len(variables),box):
            m = s.prod(y**j for y,j in zip(variables,exp))
            for j,r in enumerate(relations):
                cols.append((('w',j,exp),poly(m*r))); tags.append((j,m))
        cols.append((('p',0),{e:-c for e,c in target_col.items()}))
        rel = relation(cols,A,budget)
        if rel is None: continue
        z = rel[len(cols)-1]; out = [s.S.Zero for _ in relations]
        for i,c in rel.items():
            if i < len(tags):
                j,m = tags[i]; out[j] += A.K.to_sympy(c/z)*m
        out = [s.expand(v) for v in out]
        if s.expand(sum(v*r for v,r in zip(out,relations))-target) != 0:
            raise ValueError('ordinary ideal membership witness failed')
        return out


def linear_nilvariance(data, algebra, budget):
    """Exact all-degree linear-phase shortcut; None delegates to general search."""
    n = data['dimension']; xs = s.symbols('x:'+str(n)); ys = algebra.dependent
    A,gens,base,ES,E,f = source(data); b = n+len(ys); N = A.n
    P = sum(A.K.to_sympy(A.coefficient(c))*s.prod(z**j for z,j in zip(xs+ys,e)) for e,c in data['phase_terms'])
    poly = s.Poly(P,*xs)
    if any(sum(e) != 1 for e,c in poly.terms() if c): return None
    l = [s.diff(P,x) for x in xs]; inverse = None
    for i,c in enumerate(l):
        budget.check(); vector = algebra.vector(c)
        matrix = s.Matrix.hstack(*[s.Matrix([algebra.K.to_sympy(v) for v in algebra.mul(vector,algebra.basis_vector(j))]) for j in range(algebra.rank)])
        if matrix.det() == 0: continue
        v = matrix.inv()*s.Matrix([algebra.K.to_sympy(z) for z in algebra.unit])
        inverse = s.expand(sum(z*s.prod(y**j for y,j in zip(ys,e)) for z,e in zip(v,algebra.basis)))
        chosen = i; break
    if inverse is None: return None
    C = s.Matrix(data['precision_matrix']).inv()
    variance = s.expand((s.Matrix(l).T*C*s.Matrix(l))[0])
    ell = next((k for k in range(1,algebra.rank+1) if not any(algebra.vector(variance**k))), None)
    if ell is None: return None
    relations = [sum(A.K.to_sympy(A.coefficient(c))*s.prod(y**j for y,j in zip(ys,e)) for e,c in r) for r in data['coefficient_relations']]
    inv_w = membership(s.expand(inverse*l[chosen]-1),relations,ys,A,budget)
    var_w = membership(s.expand(variance**ell),relations,ys,A,budget)
    def yoperator(expr):
        if not ys: return A.scalar(expr)
        terms = [[list(e),str(c)] for e,c in s.Poly(expr,*ys,domain=A.K).terms() if c]
        return A.substitute(terms,base[n:b])
    vs = [A.var(i) for i in range(2*N)];u = vs[:n]; du = vs[N:N+n];dv = vs[-1]
    H = C.inv(); T = A.add(A.scale(u[chosen],-1),*[A.scale(du[j],A.coefficient(H[chosen,j])) for j in range(n)])
    Q = A.mul(yoperator(inverse),T)
    ordinary_w = [{} for _ in gens]
    ordinary_w[chosen] = A.scale(A.mul(yoperator(inverse),f),-1)
    ordinary_w[-1] = A.scale(dv,-1)
    fd_minus_one = A.sub(A.mul(f,dv),A.scalar(1))
    for j,z in enumerate(inv_w): ordinary_w[n+j] = A.mul(yoperator(z),fd_minus_one)
    # Aop+Vop=sum_i((C*u)_i+(C*l)_i*Dv)G_i + Dv*G0.
    quadratic = A.add(*[A.scale(A.mul(u[i],u[j]),A.coefficient(C[i,j])) for i in range(n) for j in range(n)])
    Aop = A.sub(A.sub(E,quadratic),ES); Vop = A.mul(yoperator(variance),A.power(dv,2))
    wb = [A.add(*[A.scale(u[j],A.coefficient(C[i,j])) for j in range(n)],
              A.mul(yoperator(sum(C[i,j]*l[j] for j in range(n))),dv)) for i in range(n)]
    wb += [{} for _ in relations]+[dv]
    witnesses = [{} for _ in gens]; scalar = s.S.One
    # Keep the central S exponent separate from the full graph Weyl pair.
    epoly, npoly = {0:A.scalar(1)}, {0:A.scalar(1)}
    def multiply_linear(poly, constant):
        result = {}
        for j,op in poly.items():
            result[j] = A.add(result.get(j,{}),A.mul(constant,op))
            result[j+1] = A.sub(result.get(j+1,{}),op)
        return result
    ss,zz = s.symbols('ss zz')
    for k in range(ell):
        budget.check(); factor = A.add(Aop,A.scalar(2*k))
        witnesses = [A.add(A.mul(factor,w),A.scale(A.mul(A.power(Vop,k),w0),(-1)**k)) for w,w0 in zip(witnesses,wb)]
        scalar *= zz-ss+2*k
        epoly = multiply_linear(epoly,A.add(E,A.scalar(2*k)))
        npoly = multiply_linear(npoly,A.add(A.sub(E,quadratic),A.scalar(2*k)))
    for j,z in enumerate(var_w):
        witnesses[n+j] = A.add(witnesses[n+j],A.scale(A.mul(yoperator(z),A.power(dv,2*ell)),(-1)**ell))
    def base_terms(operator, central_degree=0):
        terms = []
        for e,c in sorted(operator.items()):
            if e[b] or e[-1]: raise ValueError('unexpected graph variables in base operator')
            terms.append([list(e[:b]+e[N:N+b]) + [central_degree],str(A.K.to_sympy(c))])
        return terms
    return dict(ordinary=dict(kind='ordinary',method='linear unit coefficient',scalar_terms=[[[0,0],'1'],[[1,0],'1']],
                operator_terms=base_terms(Q),left_witnesses=[A.serialize(w) for w in ordinary_w]),
        Euler=dict(kind='Euler',method='linear nilpotent variance',variance_nilpotence_index=ell,
                scalar_terms=[[list(e),str(c)] for e,c in s.Poly(scalar,ss,zz).terms()],
                operator_terms=[term for j in sorted(epoly) for term in base_terms(A.sub(epoly[j],npoly[j]),j)],
                left_witnesses=[A.serialize(w) for w in witnesses]),
        shortcut_witnesses=dict(unit_coordinate=chosen,inverse=str(inverse),inverse_relation_coefficients=[str(z) for z in inv_w],
            variance=str(variance),variance_power_relation_coefficients=[str(z) for z in var_w]))


def generate(data, algebra, budget):
    A,gens,*_ = source(data)
    cert = linear_nilvariance(data,algebra,budget)
    if cert is None:
        cert = {kind:search_identity(data,kind,budget) for kind in ('ordinary','Euler')}
    cert.update(graph_generators=[A.serialize(g) for g in gens],specialization_start=0,
        torsion='Retained in exact unlocalized scalar/LEFT-witness identities; no division by scalar factors.')
    verify(data,cert); return cert


def rational_component(terms, parameters):
    """Clear QQ(U) denominators, retain one nonzero QQ component and its audit."""
    U = tuple(s.Symbol(x) for x in parameters); X,Y = S.X,S.Y
    local = {str(x):x for x in U}; expr = sum(s.sympify(c,locals=local)*X**e[0]*Y**e[1] for e,c in terms)
    numerator,denominator = s.fraction(s.cancel(expr))
    if U:
        components = s.Poly(numerator,*U).terms()
        exponent,selected = next((e,c) for e,c in components if c)
    else: exponent,selected = (),numerator
    p = s.Poly(selected,X,Y,domain=s.QQ)
    if p.is_zero: raise ValueError('zero rational polynomial component')
    return p,dict(cleared_denominator=str(denominator),cleared_numerator=str(numerator),
        selected_parameter_monomial=list(exponent),rational_polynomial_terms=S.terms(p))


def order_law(data, certificate, order):
    verify(data,certificate); n = data['dimension']; b = n+len(data['coefficient_variables'])
    D = max(1,max(sum(e[:n]) for e,c in data['phase_terms']))
    p,pa = rational_component(certificate['Euler']['scalar_terms'],data['field_parameters'])
    B,ba = rational_component(certificate['ordinary']['scalar_terms'],data['field_parameters'])
    X,Y = S.X,S.Y; B = s.Poly(B.as_expr(),X,domain=s.QQ)
    shifted = s.Poly(s.expand(p.as_expr().subs(Y,Y-D*X)),X,Y,domain=s.QQ)
    J = max((sum(e[b:b+n]) for e,c in certificate['ordinary']['operator_terms'] if c != '0'),default=0)
    start = max([0]+S.integer_roots(B))+1
    plan = S.make_plan(dict(polynomial_terms=S.terms(shifted),validity_start=start,step_bound=D+J))
    K = plan['query_operator_power']; d = order(K)
    if d is None: raise ValueError('nonnull primary phase became nilpotent')
    strips = [z for z in plan['strips'] if abs(z['b']*(d+D*K)-z['a']*K) < z['normal_radius']]
    if len(strips) != 1: raise ValueError('order failed to select a unique rational strip')
    z = strips[0]; rho = s.Rational(z['a'],z['b'])-D
    if rho < 0: raise ValueError('negative marker slope')
    offsets = []; threshold = K
    for c in range(-z['normal_radius']+1,z['normal_radius']):
        line = s.Poly(s.expand(shifted.as_expr().subs(Y,(z['a']*X+c)/z['b'])),X,domain=s.QQ)
        if line.is_zero: offsets.append(str(s.Rational(c,z['b'])))
        else: threshold = max(threshold,int(s.floor(S.root_bound(line)))+1)
    return dict(marker_slope=str(rho),safe=bool(rho>0),asymptotic_coefficient=str(1/rho) if rho else None,
        exact_line_threshold=threshold,late_offsets=offsets,shift_D=D,ordinary_upper_step=J,
        ordinary_step_start=start,query_power=K,original_order_at_query=d,shifted_strip=z,
        ordinary_rational_projection=ba,Euler_rational_projection=pa,
        arithmetic_scope='Shifted INTEGER detecting orders; no old formal-operator safety labels or periodicity imported.')


def line_envelope(law, order, C, E):
    """A sufficient positive vertex proof using EVERY possible late line.

    A failed envelope is inconclusive; exact rational zero separation follows.
    No recurrence or periodicity of the actual line selection is assumed.
    """
    C,E = Fraction(C),Fraction(E); rho = Fraction(law['marker_slope'])
    if rho <= 0 or C*rho < 1: raise ValueError('outside positive-slope recession cone')
    K = law['exact_line_threshold']; early = []
    for m in range(1,K):
        q = order(m); early.append(dict(power=m,least_marker_degree=q))
        if q is not None and m > C*q+E: return None
    late = []; a,b = rho.numerator,rho.denominator
    for offset in law['late_offsets']:
        c = Fraction(offset); bc = c*b
        if bc.denominator != 1: raise ValueError('incompatible rational line offset')
        lower = max(1,K,math.ceil(-c/rho))
        residue = (-int(bc)*pow(a,-1,b))%b if b > 1 else 0
        m = lower+(residue-lower)%b; q = rho*m+c
        if q.denominator != 1 or q < 0: raise ValueError('invalid first admissible line point')
        late.append(dict(offset=str(c),first_admissible_power=m,least_marker_degree=int(q),deficit=str(m-C*q)))
        if m-C*q > E: return None
    return dict(status='DECIDED',zero=True,method='certified finite parallel-line envelope',C=str(C),E=str(E),
        early_orders=early,late_line_bounds=late,exact_line_threshold=K,
        justification='Deficit is nonincreasing along each candidate line,including lines never actually selected.')


def solve(data, out, timeout=30):
    """Full primary decomposition, faithful component certificates and regions."""
    out = Path(out); out.mkdir(parents=True,exist_ok=True)
    n,r,phase,ideal,H = O.F.normalize(data); budget = O.Budget(timeout)
    assurance = dict(independent_replay=False,formal_verification=False,external_review=False)
    cs = s.symbols('c:'+str(r)); gb = s.groebner([O.polynomial(g,cs) for g in data.get('parameter_ideal',[])],*cs,domain=s.QQ)
    results = []; components = []
    if gb.reduce(s.S.One)[1] == 0:
        return dict(status='DECIDED',input=data,method='zero ring',region=dict(asymptotic_coefficient='0',pure_linear_minimum='0',vertices=[['0','0']],supporting_original_points=[],recession_generators=[[1,0],[0,1]]),**assurance)
    try:
        components = O.primary_components(data.get('parameter_ideal',[]),r,out/'primary',budget)
        for j,component in enumerate(components):
            budget.check(); A = O.FiniteAlgebra(r,component,n); moments = O.Moments(data,A)
            family = dict(data,parameter_ideal=copy.deepcopy(component['primary_ideal']))
            order_cache = {}
            def order(m):
                budget.check()
                if m not in order_cache: order_cache[m] = G.heat_order(family,m)
                return order_cache[m]
            item = dict(component_index=j,finite_algebra=A.manifest())
            if not moments.power(A.rank):
                orders = [order(m) for m in range(1,A.rank+1)]
                points = [(m,q) for m,q in enumerate(orders,1) if q is not None]
                region = dict(status='DECIDED',asymptotic_coefficient='0',vertices=[[str(C),str(E)] for C,E in G.vertices(0,points)],
                    supporting_original_points=[dict(power=m,least_marker_degree=q) for m,q in points],recession_generators=[[1,0],[0,1]])
                item.update(method='finite nilpotence table',nilpotence_index=next(m for m,q in enumerate(orders,1) if q is None),
                    original_orders=orders,region=region)
            else:
                inp = source_input(data,A); cert = generate(inp,A,budget)
                law = order_law(inp,cert,order); item.update(source_input=inp,certificate=cert,order_law=law)
                if not law['safe']:
                    results.append(item)
                    return dict(status='UNSAFE',input=data,primary_components=components,component_results=results,
                        unsafe_component=j,region_empty=True,**assurance)
                calls = [0]
                def oracle(C,E):
                    proof = line_envelope(law,order,C,E)
                    if proof is not None: return proof
                    calls[0] += 1
                    remaining = budget.remaining()
                    return G.rational_zero(family,C,E,out/('component_'+str(j))/('vertex_'+str(calls[0])),
                        timeout=max(1,math.ceil(remaining)) if remaining is not None else 0)
                region = G.refine_region(Fraction(law['asymptotic_coefficient']),order,oracle)
                item['region'] = region
                if region['status'] != 'DECIDED':
                    results.append(item)
                    return dict(status='UNKNOWN',input=data,primary_components=components,component_results=results,**assurance)
            results.append(item)
        lower = max([Fraction(0)]+[Fraction(z['region']['asymptotic_coefficient']) for z in results])
        powers = sorted({z['power'] for item in results for z in item['region']['supporting_original_points']})
        points = [(m,G.heat_order(data,m)) for m in powers]
        if any(q is None for m,q in points): raise ValueError('faithful primary intersection lost a nonzero moment')
        vs = G.vertices(lower,points)
        pure = [C for C,E in vs if E == 0]
        return dict(status='DECIDED',input=data,primary_components=components,component_results=results,
            region=dict(asymptotic_coefficient=str(lower),vertices=[[str(C),str(E)] for C,E in vs],
                supporting_original_points=[dict(power=m,least_marker_degree=q) for m,q in points],recession_generators=[[1,0],[0,1]],
                pure_linear_minimum=str(min(pure)) if pure else None,
                equality='Intersection of ALL faithful primary component regions; support orders recomputed in original ring.'),
            general_presented_ring_backend_implemented=True,primary_decomposition_trusted_CAS=True,
            actual_image_kernel_computed=False,**assurance)
    except (O.BudgetExpired,TimeoutError):
        return dict(status='UNKNOWN',input=data,primary_components=components,component_results=results,
            reason='finite primary/graph/vertex cap; retained partial results are not a decision',**assurance)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path); parser.add_argument('output',type=Path)
    parser.add_argument('--timeout',type=int,default=30)
    args = parser.parse_args(); result = solve(json.loads(args.input.read_text()),args.output.with_suffix('.work'),args.timeout)
    args.output.write_text(json.dumps(result,indent=2)+'\n')


if __name__ == '__main__': main()
