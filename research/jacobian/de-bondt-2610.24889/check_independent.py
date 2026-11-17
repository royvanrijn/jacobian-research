"""Standard-library-only replay, separate from SymPy's construction.

Reuses the archived Fraction sparse arithmetic engine, but transcribes the
new map separately. Checks a full 12x12 polynomial determinant, the 18-row
cube certificate, the integer 90x90 matrix and its exact rational collision.
"""
from __future__ import annotations
import ast
from fractions import Fraction as F
import importlib.util
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
ENGINE=ROOT.parents[1]/'foundations/scripts/audit_macfarlane_f12_independent.py'
spec=importlib.util.spec_from_file_location('sparse_independent',ENGINE)
p=importlib.util.module_from_spec(spec);spec.loader.exec_module(p)

FORMULAS=[
    'x1-3*x1**2*x2+x6*x7+x1**2*x7+x1*x3*x6',
    'x2+3*x1*x3+36*x1*x2**2+12*x1*x8+12*x4*x10+36*x1*x2*x10+12*x1*x4*x8-12*x7*x12-12*x1*x3*x12-12*x1*x2*x7-12*x1*x7*x10',
    'x3+8*x2**2+6*x1*x2*x3+8*x1*x5-12*x4*x8-36*x1*x2*x8-12*x2**2*x4+8*x4*x9+24*x1*x2*x9+8*x1*x4*x5+12*x7*x10+12*x1*x3*x10+12*x1*x7*x8-4*x7*x11-4*x1*x3*x11-8*x1*x7*x9',
    'x4+x1**2*x3+3*x1*x2','x5+x2**3','x6+x1**2',
    'x7+x1*x3','x8+x2**2','x9+x1*x5','x10+x1*x8',
    'x11+2*x1*x9','x12+x1*x2+x1*x10',
]


def power(a,n):
    result=p.const(1)
    for _ in range(n):result=p.multiply(result,a)
    return result


def parse(text):
    def visit(node):
        if isinstance(node,ast.Expression):return visit(node.body)
        if isinstance(node,ast.Constant):return p.const(node.value)
        if isinstance(node,ast.Name):return p.var(int(node.id[1:])-1)
        if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return p.scale(visit(node.operand),-1)
        if isinstance(node,ast.BinOp):
            if isinstance(node.op,ast.Add):return p.add(visit(node.left),visit(node.right))
            if isinstance(node.op,ast.Sub):return p.add(visit(node.left),p.scale(visit(node.right),-1))
            if isinstance(node.op,ast.Mult):return p.multiply(visit(node.left),visit(node.right))
            if isinstance(node.op,ast.Pow):return power(visit(node.left),ast.literal_eval(node.right))
        raise ValueError(ast.dump(node))
    return visit(ast.parse(text,mode='eval'))


def homogeneous(poly,d):return {m:c for m,c in poly.items() if sum(m)==d}
def extend(poly):return {m+(0,)*6:c for m,c in poly.items()}
def matvec(M,v):return [sum(F(c)*a for c,a in zip(row,v)) for row in M]
def cube(v):return [a**3 for a in v]
def addvec(a,b):return [x+y for x,y in zip(a,b)]


def rank(matrix):
    a=[[F(c) for c in row] for row in matrix];row=0
    for col in range(len(a[0])):
        pivot=next((i for i in range(row,len(a)) if a[i][col]),None)
        if pivot is None:continue
        a[row],a[pivot]=a[pivot],a[row]
        d=a[row][col];a[row]=[v/d for v in a[row]]
        for i in range(row+1,len(a)):
            if a[i][col]:
                c=a[i][col];a[i]=[v-c*w for v,w in zip(a[i],a[row])]
        row+=1
        if row==len(a):break
    return row


def main():
    k=[parse(f) for f in FORMULAS]
    dk=[[p.derivative(f,j) for j in range(12)] for f in k]
    determinant,states=p.determinant(dk)
    assert determinant==p.const(1)
    assert all(sum(m)<=3 for f in k for m in f)
    c=[extend(homogeneous(f,3)) for f in k[:5]]
    q=[extend(homogeneous(f,2)) for f in k]
    p.N=18;p.ZERO=(0,)*18
    x=[p.var(i) for i in range(18)]
    h=[p.multiply(x[17],f) for f in q]
    for i in range(5):h[i]=p.add(h[i],p.scale(p.multiply(power(x[17],2),x[12+i]),-1))
    h+=c+[{}]
    data=json.loads((ROOT/'polarization.json').read_text())
    L=data['forms'];B=[[F(a) for a in row] for row in data['B']]
    powers=[power(p.add(*(p.scale(v,co) for v,co in zip(x,row))),3) for row in L]
    target=h[:-1]+[power(x[17],3)]
    for row,want in zip(B,target):
        assert p.add(*(p.scale(f,co) for co,f in zip(row,powers)))==want
    assert rank(B)==18
    assert all(any(row[j] for row in B) for j in range(90))
    # Independent coefficient rank of all selected cubic powers.
    mons=sorted(set(m for f in powers for m in f))
    assert rank([[f.get(m,0) for m in mons] for f in powers])==90
    matrix=json.loads((ROOT/'druzkowski90.json').read_text())
    A=matrix['A'];B6=matrix['projection_B'];pad=matrix['padding_index_1_based']-1
    for i,row in enumerate(B):
        assert B6[i]==[co*(6 if i<17 else 1) for co in row]
    L0=[row if i!=pad else [0]*18 for i,row in enumerate(L)]
    shear=list(map(F,matrix['target_shear_vector']))
    assert shear==[row[pad]-(1 if i==17 else 0) for i,row in enumerate(B6)]
    assert shear[17]==0
    t3=power(x[17],3)
    # Verify the actual sheared companion map coefficient by coefficient.
    for i,row in enumerate(B6):
        got=p.add(*(p.scale(f,co) for j,(co,f) in enumerate(zip(row,powers)) if j!=pad))
        assert got==p.add(p.scale(h[i],6),p.scale(t3,-shear[i]))
    for a,row in zip(A,L0):
        assert a==[sum(row[i]*B6[i][j] for i in range(18)) for j in range(90)]
    u=list(map(F,matrix['u']));v=list(map(F,matrix['v']))
    fu=addvec(u,cube(matvec(A,u)));fv=addvec(v,cube(matvec(A,v)))
    assert u!=v and fu==fv==list(map(F,matrix['image']))
    # Symbolic cube factorization + direct determinant above certify det DF=1
    # by the companion Schur complement and Sylvester identities in README.md.
    powerA=A;power_ranks=[]
    for degree in range(1,5):
        power_ranks.append(rank(powerA))
        if degree<4:
            powerA=[[sum(row[k]*A[k][j] for k in range(90)) for j in range(90)] for row in powerA]
    assert power_ranks==[18,2,1,0]
    result={'independent_engine':'stdlib Fraction sparse polynomials; no SymPy',
            'determinant12':1,'determinant_memo_states':states,
            'cube_identity_rows':18,'cube_power_rank':90,'B_rank':18,
            'integer_matrix90':True,'rational_collision90':True,
            'projection_to_sheared_companion':True,
            'matrix_power_ranks':power_ranks,
            'assurance':'exact identity replay plus written determinant proof; no formal verification or external review'}
    (ROOT/'independent.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
