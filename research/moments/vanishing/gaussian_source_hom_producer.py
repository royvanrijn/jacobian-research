"""Whole unlocalized Gaussian integration and polynomial Hom, using Singular.

Exact QQ / simple number-field inputs. Singular's integration implementation is
instrumented to retain its pre-pruning presentation and add recurrence/IBP
certificate data, preserving the original integration computations.
No time localization or actual-image kernel computation occurs. The source
presentation relies on CAS integration; Hom bounds and LEFT witnesses are
checked separately with exact arithmetic.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
import subprocess

import sympy as sp

LIB = Path('/usr/share/singular/LIB/dmodapp.lib')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


class Field:
    def __init__(self, spec=None):
        self.a = sp.Symbol('a')
        self.spec = spec or {'kind': 'QQ'}
        if self.spec == {'kind': 'QQ'}:
            self.K = sp.QQ
            self.root = None
            self.singular = '0'
            self.minpoly = ''
        else:
            if set(self.spec) != {'kind', 'minimal_polynomial'} or self.spec['kind'] != 'number_field':
                raise ValueError('field must be QQ or a simple number_field')
            f = sp.Poly(self.expression(self.spec['minimal_polynomial']), self.a, domain=sp.QQ)
            if f.degree() < 2 or not f.is_irreducible:
                raise ValueError('minimal polynomial must be irreducible of degree >=2')
            self.root = sp.CRootOf(f.as_expr(), 0)
            self.K = sp.QQ.algebraic_field(self.root)
            self.singular = '(0,a)'
            self.minpoly = 'minpoly=' + str(f.monic().as_expr()).replace('**', '^') + ';'

    def expression(self, value):
        value = str(value)
        if not re.fullmatch(r'[0-9a +*/^().\-]+', value):
            raise ValueError('coefficient syntax: rational expressions in a only')
        e = sp.sympify(value.replace('^', '**'), locals={'a': self.a})
        if e.has(sp.Float) or e.free_symbols - {self.a}:
            raise ValueError('exact coefficients required')
        return e

    def parse(self, value):
        e = self.expression(value)
        if self.root is None:
            if self.a in e.free_symbols:
                raise ValueError('a requires a number field')
            return self.K.from_sympy(e)
        return self.K.from_sympy(e.subs(self.a, self.root))

    def encode(self, value):
        if self.root is None:
            return str(value)
        coeffs = value.to_list()
        expr = sum(sp.Rational(c.numerator, c.denominator) * self.a**(len(coeffs)-i-1)
                   for i, c in enumerate(coeffs))
        return str(expr).replace('**', '^')


def normalize(data):
    if set(data) - {'dimension', 'phase_terms', 'field'}:
        raise ValueError('unknown input keys')
    n = data['dimension']
    if type(n) is not int or n < 1:
        raise ValueError('dimension must be positive')
    f = Field(data.get('field'))
    terms = {}
    for powers, coefficient in data['phase_terms']:
        if len(powers) != n or any(type(k) is not int or k < 0 for k in powers):
            raise ValueError('invalid phase exponents')
        p = tuple(powers)
        terms[p] = terms.get(p, f.K.zero) + f.parse(coefficient)
    terms = {p: c for p, c in terms.items() if c}
    return f, n, terms


def instrument_library(destination):
    """Preserve integration computations; retain maps and add certificate data."""
    original = LIB.read_text()
    replacements = [
        ('    vector v = pIntersect(s,inG);\n    list L = bFactor(vec2poly(v));',
         '    vector v = pIntersect(s,inG);\n'
         '    vector gaussianIntegrationBVector=v; export(gaussianIntegrationBVector);\n'
         '    list L = bFactor(vec2poly(v));'),
        ('    kill inG,v,s;',
         '    ideal gaussianInitial=inForm(G,ww);\n'
         '    poly gaussianEulerB=0;\n'
         '    for(i=1;i<=nrows(v);i++){gaussianEulerB=gaussianEulerB+v[i]*s^(i-1);}\n'
         '    matrix gaussianLift=lift(gaussianInitial,ideal(gaussianEulerB));\n'
         '    poly gaussianRecurrenceOperator=0;poly gaussianCoefficient,gaussianPart;\n'
         '    for(i=1;i<=size(G);i++){\n'
         '      gaussianCoefficient=gaussianLift[i,1];gaussianPart=0;\n'
         '      while(gaussianCoefficient<>0){\n'
         '        if(deg(lead(gaussianCoefficient),ww)==-deg(G[i],ww)){\n'
         '          gaussianPart=gaussianPart+lead(gaussianCoefficient);\n'
         '        }\n'
         '        gaussianCoefficient=gaussianCoefficient-lead(gaussianCoefficient);\n'
         '      }\n'
         '      gaussianRecurrenceOperator=gaussianRecurrenceOperator+gaussianPart*G[i];\n'
         '    }\n'
         '    export(gaussianRecurrenceOperator);\n'
         '    kill inG,v,s;'),
        ('  ideal B,M;\n  poly p;',
         '  ideal B,M;\n  ideal gaussianUnprojectedRelations;\n  poly p,gaussianUnprojected;'),
        ('        p = B[k]*G[i];\n        p = f(p);\n        M[size(M)+1] = p;',
         '        p = B[k]*G[i];gaussianUnprojected=p;\n        p = f(p);\n'
         '        gaussianUnprojectedRelations[size(M)+1]=gaussianUnprojected;\n'
         '        M[size(M)+1] = p;'),
        ('  return(list(Bl0,M));',
         '  export(gaussianUnprojectedRelations);\n  return(list(Bl0,M));'),
        ('  if (l0 < 0)                  // maximal integral root is < 0',
         '  poly gaussianIntegrationCutoff=l0; export(gaussianIntegrationCutoff);\n'
         '  ideal gaussianAdaptedBasis=G; export(gaussianAdaptedBasis);\n'
         '  if (l0 < 0)                  // maximal integral root is < 0'),
        ('  module M = coeffs(N,B,prodvar);',
         '  ideal gaussianAmplitudeBasis=B; export(gaussianAmplitudeBasis);\n'
         '  module M = coeffs(N,B,prodvar);'),
        ('  M = engine(M,eng);\n  M = prune(M);',
         '  module gaussianWholeRelations=M; export(gaussianWholeRelations);\n'
         '  M = engine(M,eng);\n  M = prune(M);')]
    adapted = original
    for old, new in replacements:
        if adapted.count(old) != 1:
            raise ValueError('unsupported dmodapp.lib layout')
        adapted = adapted.replace(old, new)
    Path(destination).write_text(adapted)
    return dict(original_library=str(LIB), original_sha256=sha(LIB),
                instrumented_sha256=sha(destination), original_integration_computation_changes=False,
                additional_certificate_computations=True,
                added_exports=['gaussianIntegrationCutoff', 'gaussianAdaptedBasis',
                               'gaussianAmplitudeBasis', 'gaussianWholeRelations',
                               'gaussianIntegrationBVector', 'gaussianRecurrenceOperator',
                               'gaussianUnprojectedRelations'])


EMITTER = r'''
proc emitPolynomial(poly p,string tag)
{
  print("POLY|"+tag);
  while(p<>0){
    print("TERM|"+string(leadcoef(p))+"|"+string(leadexp(p)));
    p=p-lead(p);
  }
  print("ENDPOLY");
}
'''


def source_script(f, n, terms, library):
    xs = ['x'+str(i+1) for i in range(n)]
    ds = ['Dx'+str(i+1) for i in range(n)]
    phase = '+'.join('('+f.encode(c)+')*' + '*'.join(
        [x+'^'+str(e) for x, e in zip(xs, p) if e] or ['1']) for p, c in terms.items()) or '0'
    coordinates = ','.join(xs+['t']+ds+['Dt'])
    generators = [d+'+'+x+'-t*diff(P,'+x+')' for d, x in zip(ds, xs)] + ['Dt-P']
    return f'''LIB "dmodapp.lib";
LIB "{library}";
short=0;
{EMITTER}
ring R={f.singular},({coordinates}),dp;
{f.minpoly}
def W=Weyl();setring W;
poly P={phase};
ideal I={','.join(generators)};
intvec w={','.join(['1']*n+['0'])};
def integrated=integralModule(I,w);
setring W;
print("CUTOFF|"+string(gaussianIntegrationCutoff));
int i,j,k;
for(i=1;i<=nrows(gaussianIntegrationBVector);i++){{emitPolynomial(gaussianIntegrationBVector[i],"IB|"+string(i));}}
for(i=1;i<=size(gaussianAmplitudeBasis);i++){{
  emitPolynomial(gaussianAmplitudeBasis[i],"B_RAW|"+string(i));
  emitPolynomial(gaussianAmplitudeBasis[i]/leadcoef(gaussianAmplitudeBasis[i]),"B|"+string(i));
}}
for(i=1;i<=size(gaussianAdaptedBasis);i++){{emitPolynomial(gaussianAdaptedBasis[i],"G|"+string(i));}}
intvec integratedDirections=1..{n};
ideal RC=inverseFourier(ideal(gaussianRecurrenceOperator),integratedDirections);
ideal UP=inverseFourier(gaussianUnprojectedRelations,integratedDirections);
for(i=1;i<=ncols(UP);i++){{emitPolynomial(UP[i],"UNPROJECTED_RELATION|"+string(i));}}
emitPolynomial(RC[1],"RECURRENCE");
matrix recurrenceLift=lift(I,RC);
for(i=1;i<=size(I);i++){{emitPolynomial(recurrenceLift[i,1],"RECURRENCE_WITNESS|"+string(i));}}
setring integrated;
module M=gaussianWholeRelations;
int rr=nrows(M);int cc=ncols(M);
print("SHAPE|"+string(rr)+"|"+string(cc));
for(i=1;i<=rr;i++){{for(j=1;j<=cc;j++){{emitPolynomial(M[i,j],"M|"+string(i)+"|"+string(j));}}}}
module E,S;poly A;int selected;
for(i=1;i<=rr;i++){{
  E=M;E[ncols(E)+1]=gen(i);S=syz(E);selected=0;
  for(j=1;j<=ncols(S);j++){{if(S[ncols(E),j]<>0){{selected=j;break;}}}}
  if(selected==0){{ERROR("no nonzero coordinate annihilator");}}
  A=S[ncols(E),selected];emitPolynomial(A,"A|"+string(i));
  for(j=1;j<=cc;j++){{
    emitPolynomial(-S[j,selected],"W|"+string(i)+"|"+string(j));
  }}
}}
print("GAUSSIAN_SOURCE_COMPLETE");quit;
'''


def parse_output(output, f):
    if 'GAUSSIAN_SOURCE_COMPLETE' not in output or re.search(r'^\s*\?', output, re.M):
        raise ValueError('Singular did not complete cleanly')
    polys = {}
    current = None
    shape = cutoff = None
    for line in output.splitlines():
        if line.startswith('POLY|'):
            current = tuple(line.split('|')[1:]);polys[current] = []
        elif line.startswith('TERM|'):
            _, c, es = line.split('|')
            powers = [int(v.strip()) for v in es.split(',')]
            polys[current].append([powers, f.encode(f.parse(c))])
        elif line == 'ENDPOLY':
            current = None
        elif line.startswith('SHAPE|'):
            shape = tuple(map(int, line.split('|')[1:]))
        elif line.startswith('CUTOFF|'):
            cutoff = int(line.split('|')[1])
    if shape is None or cutoff is None:
        raise ValueError('incomplete source output')
    rank, columns = shape
    basis = [polys['B', str(i)] for i in range(1, rank+1)]
    if any(len(ts[0][0]) % 2 for ts in basis if ts):
        raise ValueError('invalid amplitude basis')
    return dict(rank=rank, relation_columns=columns, integration_cutoff=cutoff,
                amplitude_basis=basis,
                inverse_fourier_raw_basis=[polys['B_RAW', str(i)] for i in range(1, rank+1)],
                relations=[[polys['M', str(i), str(j)] for j in range(1, columns+1)]
                           for i in range(1, rank+1)],
                adapted_fourier_basis=[v for k, v in polys.items() if k[0] == 'G'],
                integration_b_coefficients=[v for k, v in polys.items() if k[0] == 'IB'],
                original_recurrence_operator=polys['RECURRENCE',],
                original_recurrence_left_witnesses=[v for k, v in polys.items() if k[0] == 'RECURRENCE_WITNESS'],
                unprojected_relation_operators=[v for k, v in polys.items() if k[0] == 'UNPROJECTED_RELATION'],
                coordinate_annihilators=[polys['A', str(i)] for i in range(1, rank+1)],
                coordinate_annihilator_left_witnesses=[
                    [polys['W', str(i), str(j)] for j in range(1, columns+1)]
                    for i in range(1, rank+1)])


def nullspace(rows, columns, K):
    """Exact field row reduction, including the empty-system case."""
    mat = [list(row) for row in rows if any(row)]
    pivots = []
    row = 0
    for col in range(columns):
        found = next((k for k in range(row, len(mat)) if mat[k][col]), None)
        if found is None:
            continue
        mat[row], mat[found] = mat[found], mat[row]
        lead = mat[row][col]
        mat[row] = [v / lead for v in mat[row]]
        for k in range(len(mat)):
            if k != row and mat[k][col]:
                scale = mat[k][col]
                mat[k] = [x-scale*y for x, y in zip(mat[k], mat[row])]
        pivots.append(col);row += 1
        if row == len(mat):
            break
    result = []
    for free in sorted(set(range(columns))-set(pivots)):
        v = [K.zero]*columns;v[free] = K.one
        for k, p in enumerate(pivots):
            v[p] = -mat[k][free]
        result.append(v)
    return result


def operator_action(terms, degree, f):
    out = {}
    for (power, derivative), coefficient in terms:
        if degree >= derivative:
            e = power+degree-derivative
            c = f.parse(coefficient)*f.K.convert(sp.ff(degree, derivative))
            out[e] = out.get(e, f.K.zero)+c
    return {e: c for e, c in out.items() if c}


def weyl_product(left, right, f):
    result = {}
    for (a, b), x in left:
        for (c, d), y in right:
            for j in range(min(b, c)+1):
                e = (a+c-j, b+d-j)
                factor = math.comb(b, j)*math.factorial(c)//math.factorial(c-j)
                result[e] = result.get(e, f.K.zero)+f.parse(x)*f.parse(y)*f.K.convert(factor)
    return {e: c for e, c in result.items() if c}


def certify_coordinate_bounds(source, f):
    """Check LEFT identities; bound degrees by the leading power at infinity.

    For A=sum c_ab t^a d^b and w=max(a-b), A(t^m) has highest
    coefficient c(m)=sum_(a-b=w)c_ab(m)_b. Every polynomial solution's
    degree is an integer root of c. A rational norm polynomial provides a
    computable superset of these roots over a simple number field.
    """
    bounds, certificates = [], []
    z = sp.Symbol('z')
    for i, operator in enumerate(source['coordinate_annihilators']):
        if not operator:
            raise ValueError('zero coordinate annihilator')
        for k in range(source['rank']):
            lhs = {}
            for j, witness in enumerate(source['coordinate_annihilator_left_witnesses'][i]):
                for e, c in weyl_product(witness, source['relations'][k][j], f).items():
                    lhs[e] = lhs.get(e, f.K.zero)+c
            lhs = {e: c for e, c in lhs.items() if c}
            rhs = {tuple(es): f.parse(c) for es, c in operator} if k == i else {}
            if lhs != rhs:
                raise ValueError('invalid LEFT coordinate annihilator witness')
        weight = max(a-b for (a, b), c in operator if f.parse(c))
        leading = sp.expand(sum(f.expression(c)*sp.ff(z, b) for (a, b), c in operator if a-b == weight))
        if leading == 0:
            raise ValueError('zero infinity indicial polynomial')
        if f.root is None:
            norm = sp.Poly(leading, z, domain=sp.QQ)
        else:
            minimal = f.expression(f.spec['minimal_polynomial'])
            norm = sp.Poly(sp.resultant(minimal, leading, f.a), z, domain=sp.QQ)
        if norm.is_zero:
            raise ValueError('zero indicial norm')
        roots = sorted(int(r) for r in norm.ground_roots() if r.is_Integer and r >= 0)
        bound = max(roots, default=-1)
        bounds.append(bound)
        certificates.append(dict(weight=weight, leading_indicial=str(leading).replace('**', '^'),
                                 rational_norm=str(norm.as_expr()).replace('**', '^'),
                                 nonnegative_integer_norm_roots=roots, degree_bound=bound,
                                 left_witness_verified=True))
    source['polynomial_degree_bounds'] = bounds
    source['coordinate_degree_certificates'] = certificates


def check_integration_basis(source, n, f):
    cutoff = source['integration_cutoff']
    if cutoff < 0 or source['rank'] != math.comb(n+cutoff, n):
        raise ValueError('incomplete integration amplitude basis')
    seen = set()
    for terms, raw in zip(source['amplitude_basis'], source['inverse_fourier_raw_basis']):
        if len(terms) != 1:
            raise ValueError('integration basis must be signed monomials')
        powers, c = terms[0]
        alpha = tuple(powers[:n])
        if len(powers) != 2*n+2 or any(powers[n:]) or sum(alpha) > cutoff or alpha in seen:
            raise ValueError('invalid integration amplitude basis')
        if f.parse(c) != f.K.one or raw != [[powers, str((-1)**sum(alpha))]]:
            raise ValueError('wrong coefficient-module amplitude normalization')
        seen.add(alpha)
    z = sp.Symbol('z')
    coefficients = []
    for terms in source['integration_b_coefficients']:
        if any(any(powers) for powers, c in terms):
            raise ValueError('nonconstant integration b coefficient')
        coefficients.append(sum((f.expression(c) for powers, c in terms), sp.Integer(0)))
    b = sum(c*z**i for i, c in enumerate(coefficients))
    norm = b if f.root is None else sp.resultant(f.expression(f.spec['minimal_polynomial']), b, f.a)
    norm = sp.Poly(norm, z, domain=sp.QQ)
    if norm.is_zero:
        raise ValueError('zero integration b polynomial')
    roots = sorted(int(r) for r in norm.ground_roots() if r.is_Integer)
    if not roots or max(roots) != cutoff:
        raise ValueError('integration cutoff does not equal maximal integer root')
    source['integration_b_rational_norm'] = str(norm.as_expr()).replace('**', '^')
    source['integration_integer_roots'] = roots
    source['all_cutoff_amplitudes_retained'] = True


def sparse_product(left, right, f):
    out = {}
    for ex, c in left.items():
        m = len(ex)//2;a, b = ex[:m], ex[m:]
        for ey, d in right.items():
            g, h = ey[:m], ey[m:]
            for ell in itertools.product(*(range(min(x, y)+1) for x, y in zip(b, g))):
                scalar = c*d
                for bi, gi, li in zip(b, g, ell):
                    scalar *= f.K.convert(math.comb(bi, li)*math.factorial(gi)//math.factorial(gi-li))
                powers = tuple(x+y-l for x, y, l in zip(a, g, ell))+tuple(x+y-l for x, y, l in zip(b, h, ell))
                out[powers] = out.get(powers, f.K.zero)+scalar
    return {e: c for e, c in out.items() if c}


def sparse_add(target, other, f, scale=None):
    scale = f.K.one if scale is None else scale
    for e, c in other.items():
        target[e] = target.get(e, f.K.zero)+scale*c
        if not target[e]:
            del target[e]


def commutative_product(left, right, f):
    out = {}
    for a, c in left.items():
        for b, d in right.items():
            es = tuple(x+y for x, y in zip(a, b))
            out[es] = out.get(es, f.K.zero)+c*d
    return {e: c for e, c in out.items() if c}


def derivative(poly, i, f):
    out = {}
    for es, c in poly.items():
        if es[i]:
            powers = list(es);powers[i] -= 1
            out[tuple(powers)] = c*f.K.convert(es[i])
    return out


def connections(n, phase, f):
    m = n+1
    P = {tuple(es)+(0,): c for es, c in phase.items()}
    hs = []
    for i in range(n):
        powers = [0]*m;powers[i] = 1
        h = {tuple(powers): -f.K.one}
        for es, c in derivative(P, i, f).items():
            powers = list(es);powers[n] += 1
            sparse_add(h, {tuple(powers): c}, f)
        hs.append(h)
    return hs+[P]


def connection_apply(poly, i, hs, f):
    out = derivative(poly, i, f)
    sparse_add(out, commutative_product(hs[i], poly, f), f)
    return out


def exponential_action(operator, amplitude, hs, f):
    out = {};m = len(hs)
    cache = {}
    for es, coefficient in operator.items():
        a, b = es[:m], es[m:]
        if b not in cache:
            value = dict(amplitude)
            for i, count in enumerate(b):
                for _ in range(count):
                    value = connection_apply(value, i, hs, f)
            cache[b] = value
        value = {tuple(x+y for x, y in zip(a, e)): coefficient*c for e, c in cache[b].items()}
        sparse_add(out, value, f)
    return out


def lift_time_operator(terms, n, f):
    return {(0,)*n+(a,)+(0,)*n+(b,): f.parse(c) for (a, b), c in terms}


def divide_spatial_derivatives_left(operator, n, f):
    """R=sum D_i Q_i+remainder; use X^a D_i=D_i X^a-a_i X^(a-e_i)."""
    queue = dict(operator);quotients = [{} for _ in range(n)];remainder = {}
    m = n+1
    while queue:
        es, c = queue.popitem()
        i = next((i for i in range(n) if es[m+i]), None)
        if i is None:
            sparse_add(remainder, {es: c}, f)
            continue
        powers = list(es);powers[m+i] -= 1
        sparse_add(quotients[i], {tuple(powers): c}, f)
        if powers[i]:
            factor = powers[i];powers[i] -= 1
            sparse_add(queue, {tuple(powers): -c*f.K.convert(factor)}, f)
    return quotients, remainder


def certify_original_ibp(source, n, phase, f):
    """Build and check polynomial IBP witnesses for relations and scalar ODEs."""
    m = n+1;hs = connections(n, phase, f);one = {(0,)*m: f.K.one}
    encode = lambda op: [[list(e), f.encode(c)] for e, c in sorted(op.items())]
    basis = [{tuple(es[:m]): f.parse(c) for es, c in b} for b in source['amplitude_basis']]
    columns = source['relation_columns']
    if len(source['unprojected_relation_operators']) < columns:
        raise ValueError('missing unprojected integration relations')
    relation_primitives = []
    for j in range(columns):
        op = {}
        for i in range(source['rank']):
            sparse_add(op, sparse_product(lift_time_operator(source['relations'][i][j], n, f),
                                         {e+(0,)*m: c for e, c in basis[i].items()}, f), f)
        full = {tuple(e): f.parse(c) for e, c in source['unprojected_relation_operators'][j]}
        sparse_add(full, op, f, -f.K.one)
        quotients, remainder = divide_spatial_derivatives_left(full, n, f)
        if remainder:
            raise ValueError('projection did not remove only left spatial derivatives')
        primitives = [exponential_action(q, one, hs, f) for q in quotients]
        primitives = [{e: -c for e, c in q.items()} for q in primitives]
        rhs = {}
        for i, q in enumerate(primitives):
            sparse_add(rhs, connection_apply(q, i, hs, f), f)
        if rhs != exponential_action(op, one, hs, f):
            raise ValueError('invalid original Gaussian relation IBP witness')
        relation_primitives.append(primitives)
    scalar_primitives = []
    for i, ann in enumerate(source['coordinate_annihilators']):
        primitives = [{} for _ in range(n)]
        for j, witness in enumerate(source['coordinate_annihilator_left_witnesses'][i]):
            W = lift_time_operator(witness, n, f)
            for k in range(n):
                sparse_add(primitives[k], exponential_action(W, relation_primitives[j][k], hs, f), f)
        rhs = {}
        for k, q in enumerate(primitives):
            sparse_add(rhs, connection_apply(q, k, hs, f), f)
        if rhs != exponential_action(lift_time_operator(ann, n, f), basis[i], hs, f):
            raise ValueError('invalid original Gaussian scalar-annihilator IBP witness')
        scalar_primitives.append(primitives)
    source['original_relation_IBP_primitives'] = [[encode(q) for q in cs] for cs in relation_primitives]
    source['original_coordinate_IBP_primitives'] = [[encode(q) for q in cs] for cs in scalar_primitives]
    source['all_original_IBP_witnesses_verified'] = True


def selected_scalar_periods(source, n, phase, f):
    """Decide ACTUAL base-period polynomiality by finite original origin jets."""
    z = sp.Symbol('z');periods = [];certificates = []
    P = {tuple(es): c for es, c in phase.items()}
    for i, ann in enumerate(source['coordinate_annihilators']):
        weight = min(a-b for (a, b), c in ann)
        indicial = sp.expand(sum(f.expression(c)*sp.ff(z, b) for (a, b), c in ann if a-b == weight))
        norm = indicial if f.root is None else sp.resultant(f.expression(f.spec['minimal_polynomial']), indicial, f.a)
        norm = sp.Poly(norm, z, domain=sp.QQ)
        if norm.is_zero:
            raise ValueError('zero origin indicial norm')
        roots = sorted(int(r) for r in norm.ground_roots() if r.is_Integer and r >= 0)
        jet_bound = max(roots, default=-1)
        amplitude = {tuple(es[:n]): f.parse(c) for es, c in source['amplitude_basis'][i]}
        powers = {(0,)*n: f.K.one};jets = []
        for k in range(jet_bound+1):
            product = commutative_product(amplitude, powers, f)
            expectation = sum((c*f.K.convert(wick(es)) for es, c in product.items()), f.K.zero)
            jets.append(expectation/f.K.convert(math.factorial(k)))
            powers = commutative_product(powers, P, f)
        bound = source['polynomial_degree_bounds'][i]
        images = [operator_action(ann, k, f) for k in range(bound+1)]
        rows = [[v.get(e, f.K.zero) for v in images]+[f.K.zero]
                for e in sorted(set().union(*(set(v) for v in images)))]
        for k, jet in enumerate(jets):
            rows.append([f.K.convert(int(k == j)) for j in range(bound+1)]+[-jet])
        kernel = nullspace(rows, bound+2, f.K)
        chosen = next((v for v in kernel if v[-1]), None)
        period = None if chosen is None else [[[k], f.encode(c/chosen[-1])]
                                             for k, c in enumerate(chosen[:-1]) if c]
        mismatch = None
        if chosen is None:
            coefficients = [row[:-1] for row in rows]
            target = [-row[-1] for row in rows]
            transposed = [list(col) for col in zip(*coefficients)]
            left_kernel = nullspace(transposed, len(rows), f.K)
            separating = next((w for w in left_kernel if sum((x*y for x, y in zip(w, target)), f.K.zero)), None)
            if separating is None:
                raise ValueError('inconsistent scalar system lacks a separating functional')
            mismatch = dict(left_null_weights=[f.encode(c) for c in separating],
                            value_on_original_target=f.encode(sum((x*y for x, y in zip(separating, target)), f.K.zero)))
        periods.append(period)
        certificates.append(dict(origin_weight=weight, origin_indicial=str(indicial).replace('**', '^'),
                                 origin_indicial_rational_norm=str(norm.as_expr()).replace('**', '^'),
                                 nonnegative_integer_origin_roots=roots, original_jet_bound=jet_bound,
                                 original_wick_jets=[f.encode(c) for c in jets], polynomial=period is not None,
                                 no_polynomial_match_witness=mismatch))
    return dict(all_original_base_periods_polynomial=all(v is not None for v in periods),
                selected_original_generator_periods=periods, certificates=certificates,
                source_relation_completeness_required=False,
                proof='verified scalar IBP, infinity degree bounds, origin-jet injectivity, and verified recurrence coverage')


def certify_recurrence(source, f, n, phase):
    m = n+1;zero = (0,)*(2*m)
    mon = lambda index: tuple(int(j == index) for j in range(2*m))
    original = []
    for i in range(n):
        g = {mon(m+i): f.K.one, mon(i): f.K.one}
        for powers, c in phase.items():
            if powers[i]:
                es = list(powers)+[1]+[0]*m;es[i] -= 1
                sparse_add(g, {tuple(es): -c*f.K.convert(powers[i])}, f)
        original.append(g)
    g = {mon(2*m-1): f.K.one}
    for powers, c in phase.items():
        sparse_add(g, {tuple(powers)+(0,)+(0,)*m: -c}, f)
    original.append(g)
    if len(source['original_recurrence_left_witnesses']) != len(original):
        raise ValueError('incomplete integration recurrence witness')
    total = {}
    for witness, g in zip(source['original_recurrence_left_witnesses'], original):
        sparse_add(total, sparse_product({tuple(e): f.parse(c) for e, c in witness}, g, f), f)
    recurrence = {tuple(e): f.parse(c) for e, c in source['original_recurrence_operator']}
    if total != recurrence:
        raise ValueError('invalid original-exponential recurrence identity')
    theta = {zero: f.K.convert(-n)}
    for i in range(n):
        es = list(zero);es[i] = es[m+i] = 1;theta[tuple(es)] = -f.K.one
    euler_b = {};power = {zero: f.K.one}
    for coefficients in source['integration_b_coefficients']:
        c = sum((f.parse(c) for e, c in coefficients), f.K.zero)
        sparse_add(euler_b, power, f, c)
        power = sparse_product(power, theta, f)
    remainder = dict(recurrence);sparse_add(remainder, euler_b, f, -f.K.one)
    slope = sp.Rational(0)
    for es in remainder:
        drop = sum(es[m:m+n])-sum(es[:n])
        if drop < 1:
            raise ValueError('integration recurrence does not lower original marker degree')
        slope = max(slope, sp.Rational(es[n]-es[-1], drop))
    source['recurrence_degree_slope'] = str(slope)
    source['original_recurrence_identity_verified'] = True
    source['original_degree_lowering_remainder'] = [[list(e), f.encode(c)] for e, c in sorted(remainder.items())]


def all_marker_cutoff(source, normalized):
    slope = sp.Rational(source['recurrence_degree_slope'])
    offset = sp.Rational(0)
    pure_slope = slope
    unit = None
    for basis, period in zip(source['amplitude_basis'], normalized):
        degree = sum(basis[0][0])
        time_degree = max((es[0] for es, c in period), default=-1)
        offset = max(offset, time_degree-slope*degree)
        if degree:
            pure_slope = max(pure_slope, sp.Rational(time_degree, degree))
        else:
            unit = period
    pure = unit == [([0], '1')]
    return dict(slope=str(slope), intercept=str(offset),
                convention='For every original Q, E(Q P^m)=0 for every integer m > slope*deg(Q)+intercept.',
                all_positive_pure=pure,
                pure_linear_slope=str(pure_slope) if pure else None,
                proof='checked degree-lowering recurrence plus normalized periods of ALL cutoff monomials')


def polynomial_hom(source, f):
    """Impose every module relation on all degree-bounded generator values."""
    slots = [(i, k) for i, b in enumerate(source['polynomial_degree_bounds']) for k in range(b+1)]
    rows = []
    for j in range(source['relation_columns']):
        column = [operator_action(source['relations'][i][j], k, f) for i, k in slots]
        for e in sorted(set().union(*(set(v) for v in column))):
            rows.append([v.get(e, f.K.zero) for v in column])
    kernel = nullspace(rows, len(slots), f.K)
    maps = []
    for vector in kernel:
        values = [[] for _ in range(source['rank'])]
        for c, (i, k) in zip(vector, slots):
            if c:
                values[i].append([[k], f.encode(c)])
        maps.append(values)
    return dict(dimension=len(kernel), generator_images=maps,
                scalar_unknowns=len(slots), scalar_equations=len(rows),
                degree_bound_method='leading infinity indicial polynomial of a certified nonzero coordinate annihilator; rational norm for number fields')


def wick(powers):
    if any(k % 2 for k in powers):
        return 0
    out = 1
    for k in powers:
        out *= int(sp.factorial2(k-1)) if k else 1
    return out


def decide(data, workdir, timeout=30):
    if type(timeout) is not int or timeout < 0:
        raise ValueError('timeout must be a nonnegative integer; zero means uncapped')
    f, n, phase = normalize(data)
    out = Path(workdir).resolve();out.mkdir(parents=True, exist_ok=True)
    lib = out/'gaussian_dmodapp.lib'
    provenance = instrument_library(lib)
    script = out/'source.sing';script.write_text(source_script(f, n, phase, lib))
    try:
        run = subprocess.run(['Singular', '-q', '--random=1', str(script)],
                             capture_output=True, text=True, timeout=timeout or None)
    except subprocess.TimeoutExpired as exc:
        (out/'source.stdout').write_bytes((exc.stdout or b'')+(exc.stderr or b''))
        return dict(status='UNKNOWN', reason='integration/Hom CAS timeout', input=data)
    output = run.stdout+run.stderr;(out/'source.stdout').write_text(output)
    if run.returncode:
        raise RuntimeError('Singular failed; retained stdout')
    source = parse_output(output, f)
    check_integration_basis(source, n, f)
    certify_recurrence(source, f, n, phase)
    certify_coordinate_bounds(source, f)
    certify_original_ibp(source, n, phase, f)
    hom = polynomial_hom(source, f)
    selected = selected_scalar_periods(source, n, phase, f)
    if bool(hom['dimension']) != selected['all_original_base_periods_polynomial']:
        raise ValueError('whole-source Hom disagrees with certified actual-period selection')
    if hom['dimension'] > 1:
        raise ValueError('whole Gaussian source failed its formal-Hom dimension control')
    normalized = None
    if hom['dimension']:
        values = hom['generator_images'][0]
        a = None
        for b, image in zip(source['amplitude_basis'], values):
            if b == [[[0]*(2*n+2), '1']]:
                a = sum((f.parse(c) for (k,), c in image if k == 0), f.K.zero)
        if not a:
            raise ValueError('missing unit amplitude or nonzero Gaussian normalization')
        normalized = [[(powers, f.encode(f.parse(c)/a)) for powers, c in image] for image in values]
        for b, image in zip(source['amplitude_basis'], normalized):
            expected = f.K.zero
            for powers, c in b:
                if any(powers[n:]):
                    raise ValueError('amplitude basis must consist of original X monomials')
                expected += f.parse(c)*f.K.convert(wick(powers[:n]))
            actual = sum((f.parse(c) for (k,), c in image if k == 0), f.K.zero)
            if expected != actual:
                raise ValueError('normalized polynomial Hom disagrees with Wick values at zero')
        if json.loads(json.dumps(normalized)) != selected['selected_original_generator_periods']:
            raise ValueError('normalized Hom disagrees with certified original origin jets')
    return dict(status='DECIDED', input=data, coefficient_field=f.spec,
                all_original_marker_eventual_safety=bool(hom['dimension']),
                original_bad_marker_terms=next(([[es[:n], c] for es, c in source['amplitude_basis'][i]]
                                                for i, value in enumerate(selected['selected_original_generator_periods']) if value is None), None),
                exact_unlocalized_source=source, polynomial_hom=hom,
                actual_period_selection=selected,
                normalized_generator_periods=normalized,
                all_marker_cutoff=all_marker_cutoff(source, normalized) if normalized is not None else None,
                source_contract='whole Gaussian integration quotient; all polynomial amplitudes retained',
                actual_image_kernel_computed=False, time_localization_used=False,
                library_instrumentation=provenance,
                cas_files_sha256={p.name: sha(p) for p in (script, out/'source.stdout', lib)},
                verification_boundary='CAS produces the full integration presentation; its completeness is NOT needed for the decision certificate: original Gaussian IBP, scalar indicial bounds, actual Wick jets and all-marker recurrence coverage replay independently')


def verify_algebraic_parts(result):
    """Replay finite identities/Hom; does not certify CAS presentation completeness."""
    import copy
    if result['status'] != 'DECIDED':
        raise ValueError('only decided results have algebraic certificates')
    f, n, phase = normalize(result['input'])
    source = copy.deepcopy(result['exact_unlocalized_source'])
    before = copy.deepcopy(source)
    check_integration_basis(source, n, f)
    certify_recurrence(source, f, n, phase)
    certify_coordinate_bounds(source, f)
    certify_original_ibp(source, n, phase, f)
    selected = selected_scalar_periods(source, n, phase, f)
    if selected != result['actual_period_selection']:
        raise ValueError('wrong original-period selection certificate')
    if source != before or polynomial_hom(source, f) != result['polynomial_hom']:
        raise ValueError('source algebraic certificate or Hom mismatch')
    safe = result['polynomial_hom']['dimension'] == 1
    if safe != result['all_original_marker_eventual_safety']:
        raise ValueError('decision does not match whole-source Hom')
    if safe != selected['all_original_base_periods_polynomial']:
        raise ValueError('decision does not match actual original periods')
    bad_marker = next(([[es[:n], c] for es, c in source['amplitude_basis'][i]]
                       for i, value in enumerate(selected['selected_original_generator_periods']) if value is None), None)
    if bad_marker != result['original_bad_marker_terms']:
        raise ValueError('wrong original bad marker')
    periods = result['normalized_generator_periods']
    if safe:
        values = result['polynomial_hom']['generator_images'][0]
        unit = next(i for i, b in enumerate(source['amplitude_basis']) if sum(b[0][0]) == 0)
        scalar = sum((f.parse(c) for (k,), c in values[unit] if k == 0), f.K.zero)
        expected = [[[[*powers], f.encode(f.parse(c)/scalar)] for powers, c in image] for image in values]
        if json.loads(json.dumps(periods)) != expected:
            raise ValueError('wrong Gaussian normalization')
        if json.loads(json.dumps(periods)) != selected['selected_original_generator_periods']:
            raise ValueError('normalized periods do not match actual origin jets')
        for b, image in zip(source['amplitude_basis'], periods):
            expected_value = sum((f.parse(c)*f.K.convert(wick(powers[:n])) for powers, c in b), f.K.zero)
            value = sum((f.parse(c) for (k,), c in image if k == 0), f.K.zero)
            if expected_value != value:
                raise ValueError('wrong original Wick value')
        # JSON round trips replace tuple pairs by lists; normalize for comparison.
        normalized_pairs = [[(powers, c) for powers, c in image] for image in periods]
        if all_marker_cutoff(source, normalized_pairs) != result['all_marker_cutoff']:
            raise ValueError('wrong all-original-marker cutoff')
    elif periods is not None or result['all_marker_cutoff'] is not None:
        raise ValueError('unsafe result contains safe periods or cutoff')
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--workdir', type=Path)
    parser.add_argument('--timeout', type=int, default=30)
    args = parser.parse_args()
    result = decide(json.loads(args.input.read_text()), args.workdir or args.output.with_suffix('.cas'), args.timeout)
    args.output.write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
