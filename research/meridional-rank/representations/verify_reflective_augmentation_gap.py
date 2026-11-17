#!/usr/bin/env python3
"""Exact reflective cord witnesses and universal2x2 identity controls.

Standard library only. It checks rank3 cyclotomic witnesses on all8 knots and4 independent
quadratic-field control witnesses. The uniform upper bound and deck
character argument are written proofs, not finite-search conclusions.
"""
from copy import deepcopy
from fractions import Fraction as F
from itertools import combinations,permutations
from pathlib import Path
import hashlib,json,math

HERE=Path(__file__).resolve().parent
ZERO=(F(0),F(0));ONE=(F(1),F(0))


def add(*values):return (sum(x[0] for x in values),sum(x[1] for x in values))
def neg(x):return (-x[0],-x[1])
def multiply(x,y):
    a,b=x;c,d=y
    return (a*c+b*d,a*d+b*c+b*d)
def inverse(x):
    a,b=x;n=a*a+a*b-b*b
    assert n
    return ((a+b)/n,-b/n)
def scalar(x,v):return (x[0]*v,x[1]*v)


def rank(matrix):
    a=[[tuple(x) for x in row] for row in matrix]
    position=0
    for column in range(len(a[0])):
        pivot=next((i for i in range(position,len(a)) if a[i][column]!=ZERO),None)
        if pivot is None:continue
        a[position],a[pivot]=a[pivot],a[position]
        scale=inverse(a[position][column])
        a[position]=[multiply(scale,x) for x in a[position]]
        for i in range(position+1,len(a)):
            factor=a[i][column]
            if factor!=ZERO:
                a[i]=[add(x,neg(multiply(factor,y))) for x,y in zip(a[i],a[position])]
        position+=1
    return position


def determinant3(matrix):
    result=ZERO
    for p in permutations(range(3)):
        sign=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
        term=(F(sign),F(0))
        for i,j in enumerate(p):term=multiply(term,matrix[i][j])
        result=add(result,term)
    return result


def crossings(case):
    components=[{i} for i in range(1,27)]
    for row in case['pd']:
        a=next(c for c in components if row[1] in c)
        b=next(c for c in components if row[3] in c)
        if a is not b:
            components.remove(a);components.remove(b);components.append(a|b)
    components=sorted(sorted(c) for c in components)
    arc={v:i for i,c in enumerate(components) for v in c}
    return [(arc[row[1]],arc[row[0]],arc[row[2]]) for row in case['pd']]


def verify_cord_matrix(C,cross):
    assert len(C)==13 and all(len(row)==13 for row in C)
    assert all(C[i][i]==(F(2),F(0)) for i in range(13))
    assert all(C[i][j]==C[j][i] for i in range(13) for j in range(13))
    for o,a,b in cross:
        for l in range(13):
            assert add(C[l][a],C[l][b])==multiply(C[l][o],C[o][a]), 'cord relation fails'


def ring_add(a,b,p):
    out=dict(a)
    for e,c in b.items():
        out[e]=out.get(e,0)+c
        if not out[e]:del out[e]
    return out

def ring_multiply(a,b,p):
    out={}
    for e,c in a.items():
        for f,d in b.items():
            index=(e+f)%p
            out[index]=out.get(index,0)+c*d
            if not out[index]:del out[index]
    return out


class CyclotomicPrime:
    """Z[zeta] for prime p, using1+zeta+...+zeta^(p-1)=0.

    Multiplication first reduces exponents modulo p and then subtracts
    the coefficient of zeta^(p-1) times the cyclotomic polynomial.
    No floating-point roots or algebraic-number package is used.
    """
    def __init__(self,p):
        assert p>2 and all(p%d for d in range(2,math.isqrt(p)+1))
        self.p=p
        self.zero=(0,)*(p-1)
        self.one=(1,)+(0,)*(p-2)
    def const(self,a):return (a,)+(0,)*(self.p-2)
    def monomial(self,e):
        e%=self.p
        return tuple(int(i==e) for i in range(self.p-1)) if e<self.p-1 else (-1,)*(self.p-1)
    def add(self,*values):return tuple(sum(v[i] for v in values) for i in range(self.p-1))
    def neg(self,a):return tuple(-x for x in a)
    def scale(self,a,n):return tuple(n*x for x in a)
    def multiply(self,a,b):
        coeff=[0]*self.p
        for i,c in enumerate(a):
            if not c:continue
            for j,d in enumerate(b):
                if d:coeff[(i+j)%self.p]+=c*d
        last=coeff[-1]
        return tuple(c-last for c in coeff[:-1])
    def det3(self,A):
        out=self.zero
        for perm in permutations(range(3)):
            sign=(-1)**sum(perm[i]>perm[j] for i in range(3) for j in range(i+1,3))
            term=self.const(sign)
            for i,j in enumerate(perm):term=self.multiply(term,A[i][j])
            out=self.add(out,term)
        return out


def cyclotomic_witness(record,diagram):
    p=record['prime'];field=CyclotomicPrime(p)
    zero,one=field.zero,field.one
    x=field.add(field.monomial(1),field.monomial(-1))
    gram=([[field.const(2),one,x],[one,field.const(2),zero],[x,zero,field.const(2)]]
          if record['Gram_type']=='first' else
          [[field.const(2),one,zero],[one,field.const(2),field.neg(x)],
           [zero,field.neg(x),field.const(2)]])
    def dot(u,v):
        return field.add(*(field.multiply(field.multiply(u[i],gram[i][j]),v[j])
                           for i in range(3) for j in range(3)))
    e0,e1,V=[[one if i==j else zero for i in range(3)] for j in range(3)]
    U=[field.add(a,field.neg(b)) for a,b in zip(e0,e1)]
    r=[U,V]
    for _ in range(1,p+1):
        r.append([field.add(field.multiply(x,a),field.neg(b)) for a,b in zip(r[-1],r[-2])])
    assert r[p]==U and r[p+1]==V
    colours={int(i):int(v) for i,v in record['core_colours'].items()}
    assert colours[0]==0 and colours[6]==1
    def fold(i):return 0 if i==4 else 6 if i==7 else i
    roots=[None]*13
    for i in range(13):
        if i not in (1,2,3):roots[i]=r[colours[fold(i)]]
    roots[1]=e0;roots[3]=e1
    roots[2]=e1 if record['Gram_type']=='first' else e0
    signs=record['arc_root_signs']
    assert len(signs)==13 and all(s in (1,-1) for s in signs)
    roots=[[field.scale(v,s) for v in root] for root,s in zip(roots,signs)]
    assert all(dot(v,v)==field.const(2) for v in roots)
    C=[[dot(u,v) for v in roots] for u in roots]
    cross=crossings(diagram)
    def key(step):
        o,a,b=step
        return (o,)+tuple(sorted((a,b)))
    local=([(0,1,2),(1,3,4),(2,6,7),(3,1,0),(7,2,3)]
           if record['Gram_type']=='first' else
           [(3,1,0),(1,3,4),(0,3,2),(2,6,7),(7,1,2)])
    assert {key(s) for s in cross if set(s)&{1,2,3}}==set(map(key,local))
    for o,a,b in cross:
        if not {o,a,b}&{1,2,3}:
            assert (2*colours[fold(o)]-colours[fold(a)]-colours[fold(b)])%p==0
        for l in range(13):
            assert field.add(C[l][a],C[l][b])==field.multiply(C[l][o],C[o][a])
    indices=(1,2,6) if record['Gram_type']=='first' else (1,3,6)
    minor=field.det3([[C[i][j] for j in indices] for i in indices])
    expected=field.scale(field.add(field.const(3),field.neg(field.multiply(x,x))),2)
    assert minor==expected and minor!=zero
    # Factorization C=R G R^T gives rank<=3; this minor gives rank>=3.
    # Independent finite control of the abstract dihedral conjugation law.
    reflections=[tuple((2*a-v)%p for v in range(p)) for a in range(p)]
    for i in range(p):
        for j in range(p):
            product=tuple(reflections[i][reflections[j][reflections[i][v]]] for v in range(p))
            assert product==reflections[(2*i-j)%p]
    return {'name':record['name'],'cyclotomic_prime':p,'cord_matrix_rank':3,
            'principal_minor_indices':indices,'principal_minor_coefficients':minor,
            'crossing_row_equations':169,'triangle_Coxeter_type':[2,3,p],
            'dihedral_conjugation_pairs_checked':p*p}


# Sparse integer polynomials for generic2x2 controls, independent of the
# quadratic-field arithmetic used for the retained numerical witnesses.
def pc(n,value):return {(0,)*n:value} if value else {}
def pv(n,i):
    e=[0]*n;e[i]=1
    return {tuple(e):1}
def pa(*ps):
    out={}
    for p in ps:
        for e,c in p.items():
            out[e]=out.get(e,0)+c
            if not out[e]:del out[e]
    return out
def pn(p):return {e:-c for e,c in p.items()}
def pm(a,b):
    out={}
    for e,c in a.items():
        for f,d in b.items():
            k=tuple(x+y for x,y in zip(e,f));out[k]=out.get(k,0)+c*d
            if not out[k]:del out[k]
    return out

def mm(a,b):return [[pa(*(pm(a[i][k],b[k][j]) for k in range(2))) for j in range(2)] for i in range(2)]
def ma(a,b):return [[pa(a[i][j],b[i][j]) for j in range(2)] for i in range(2)]
def ms(p,a):return [[pm(p,a[i][j]) for j in range(2)] for i in range(2)]
def adj(a):return [[a[1][1],pn(a[0][1])],[pn(a[1][0]),a[0][0]]]
def tr(a):return pa(a[0][0],a[1][1])
def det(a):return pa(pm(a[0][0],a[1][1]),pn(pm(a[0][1],a[1][0])))


def matrix_identity_controls():
    a,b,c,d,e,f,g,h=[pv(8,i) for i in range(8)]
    X=[[a,b],[c,d]];Y=[[e,f],[g,h]]
    pairing=tr(mm(X,adj(Y)))
    assert pairing==tr(mm(Y,adj(X)))
    assert tr(mm(X,adj(X)))==pa(det(X),det(X))
    lhs=mm(mm(X,adj(Y)),X)
    rhs=ma(ms(pairing,X),ms(pn(det(X)),Y))
    assert lhs==rhs
    a,b,c,d,e,f=[pv(6,i) for i in range(6)]
    X=[[a,b],[c,pn(a)]];Y=[[d,e],[f,pn(d)]]
    anticommutator=ma(mm(X,Y),mm(Y,X))
    t=tr(mm(X,Y));zero={}
    assert anticommutator==[[t,zero],[zero,t]]
    # Gram matrices on Mat2 and sl2 have determinants1 and-2.
    mat2=[[0,0,0,1],[0,0,-1,0],[0,-1,0,0],[1,0,0,0]]
    sl2=[[2,0,0],[0,0,1],[0,1,0]]
    def dint(A):
        n=len(A);out=0
        for p in permutations(range(n)):
            sign=(-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
            term=sign
            for i,j in enumerate(p):term*=A[i][j]
            out+=term
        return out
    assert dint(mat2)==1 and dint(sl2)==-2
    return ['B(X,Y)=B(Y,X)','B(X,X)=2detX',
            'X adjY X=B(X,Y)X-detX Y','XY+YX=tr(XY)I for tracelessX,Y',
            'Mat2 and sl2 trace-form Gram determinants1,-2']


def main():
    data=json.loads((HERE/'reflective_augmentation_inputs.json').read_text())
    diagrams={c['name']:c for c in json.loads((HERE/'seven_sheet_cover_inputs.json').read_text())['cases']}
    positive=[];matrices={}
    for case in data['cases']:
        roots=[[(F(a),F(b)) for a,b in row] for row in case['normal_vectors']]
        assert len(roots)==13 and all(len(r)==4 for r in roots)
        C=[[scalar(add(*(multiply(a,b) for a,b in zip(u,v))),F(1,2))
            for v in roots] for u in roots]
        verify_cord_matrix(C,crossings(diagrams[case['name']]))
        assert rank(C)==3
        minor=None
        for indices in combinations(range(13),3):
            value=determinant3([[C[i][j] for j in indices] for i in indices])
            if value!=ZERO:
                minor={'indices':indices,'determinant':[str(x) for x in value]};break
        assert minor
        positive.append({'name':case['name'],'field':'Q(phi)',
                         'cord_matrix_rank':3,'nonzero_principal_minor':minor,
                         'crossing_row_equations':169})
        matrices[case['name']]=C
    assert [p['name'] for p in positive]==['13a_373','13a_392','13a_3264','13a_3276']
    cyclotomic=[cyclotomic_witness(c,diagrams[c['name']])
                for c in data['cyclotomic_rank3_witnesses']]
    assert len(cyclotomic)==8 and len({c['name'] for c in cyclotomic})==8
    controls=matrix_identity_controls()
    bad=deepcopy(matrices['13a_373'])
    bad[0][1]=bad[1][0]=add(bad[0][1],ONE)
    try:verify_cord_matrix(bad,crossings(diagrams['13a_373']))
    except AssertionError as error:assert 'cord relation fails' in str(error)
    else:raise AssertionError('corrupted cord matrix accepted')
    files=['verify_reflective_augmentation_gap.py','reflective_augmentation_inputs.json',
           'seven_sheet_cover_inputs.json','REFLECTIVE_AUGMENTATION_RANK_GAP.tex',
           'UNIVERSAL_REFLECTION_OBSTRUCTION.tex',
           'sources/cornwell-character-cord-1509.04962v1.pdf']
    receipt={'status':'PASS','quadratic_rank3_control_certificates':positive,
             'cyclotomic_rank3_positive_certificates':cyclotomic,'generic_matrix_identity_controls':controls,
             'corrupt_symmetric_cord_matrix_rejected':True,
             'scope':'Exact rank3 reflective cord matrices on all8 named knots in prime cyclotomic fields, with4 independent quadratic-field controls and universal2x2 polynomial identities. Uniform rank<=3 uses the written Gram-radical quotient and retained reflection theorem. Deck invariance uses its written proof. Explicit dihedral core colours and outer A2/commutation relations also prove onto triangle Coxeter quotients of rank3. No unrestricted augmentation-rank bound follows.',
             'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'reflective_augmentation_gap_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','exact_rank3_cases':len(cyclotomic),
                      'independent_quadratic_controls':len(positive),
                      'cord_equations_checked':12*169,'generic_matrix_identities':len(controls),
                      'corrupt_matrix_rejected':True},indent=2))


if __name__=='__main__':main()
