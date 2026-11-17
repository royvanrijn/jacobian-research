#!/usr/bin/env python3
"""Independent Laurent-rational replay of the nonzero-target endpoint obstructions."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations
import ast,hashlib,json
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class L:
    def __init__(self,value=0):
        self.p={e:Q(c) for e,c in value.items() if c} if isinstance(value,dict) else ({0:Q(value)} if value else {})
    def __add__(self,other):
        other=asL(other);out=dict(self.p)
        for e,c in other.p.items():out[e]=out.get(e,Q(0))+c
        return L(out)
    __radd__=__add__
    def __neg__(self):return L({e:-c for e,c in self.p.items()})
    def __sub__(self,other):return self+-asL(other)
    def __rsub__(self,other):return asL(other)+-self
    def __mul__(self,other):
        other=asL(other);out={}
        for e,c in self.p.items():
            for f,d in other.p.items():out[e+f]=out.get(e+f,Q(0))+c*d
        return L(out)
    __rmul__=__mul__
    def __truediv__(self,other):
        other=asL(other)
        assert len(other.p)==1,"Only monomial denominators are admitted"
        f,d=next(iter(other.p.items()))
        return L({e-f:c/d for e,c in self.p.items()})
    def __pow__(self,n):
        assert type(n) is int
        if n<0:return L(1)/(self**(-n))
        r=L(1)
        for _ in range(n):r=r*self
        return r
    def __eq__(self,other):return self.p==asL(other).p
def asL(value):return value if isinstance(value,L) else L(value)
def parse(text):
    def visit(n):
        if isinstance(n,ast.Constant) and type(n.value) is int:return L(n.value)
        if isinstance(n,ast.Name) and n.id=="k":return L({1:1})
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.USub,ast.UAdd)):
            return -visit(n.operand) if isinstance(n.op,ast.USub) else visit(n.operand)
        if isinstance(n,ast.BinOp):
            a=visit(n.left)
            if isinstance(n.op,ast.Pow):
                b=visit(n.right);assert set(b.p)<={0} and b.p.get(0,0).denominator==1
                return a**int(b.p.get(0,0))
            b=visit(n.right)
            if isinstance(n.op,ast.Add):return a+b
            if isinstance(n.op,ast.Sub):return a-b
            if isinstance(n.op,ast.Mult):return a*b
            if isinstance(n.op,ast.Div):return a/b
        raise ValueError("Unsupported exact expression")
    return visit(ast.parse(text,mode="eval").body)
def determinant(A):
    value=L()
    for p in permutations(range(3)):
        term=L((-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3)))
        for i in range(3):term*=A[i][p[i]]
        value+=term
    return value
def inverse(A):
    d=determinant(A);assert d==-2
    out=[]
    for i in range(3):
        row=[]
        for j in range(3):
            rr=[q for q in range(3) if q!=j];cc=[q for q in range(3) if q!=i]
            row.append((-1)**(i+j)*(A[rr[0]][cc[0]]*A[rr[1]][cc[1]]-A[rr[0]][cc[1]]*A[rr[1]][cc[0]])/d)
        out.append(row)
    assert matmul(A,out)==[[L(int(i==j)) for j in range(3)] for i in range(3)]
    return out
def matmul(A,B):
    return [[sum((A[i][q]*B[q][j] for q in range(3)),L()) for j in range(len(B[0]))] for i in range(3)]
def F(point):
    x,y,z=point;u=1+x*y
    return [u**3*z+y**2*u*(4+3*x*y),
            y+3*x*u**2*z+3*x*y**2*(4+3*x*y),
            2*x-3*x**2*y-x**3*z]
def J(point):
    # Explicit differentiated formulas; no symbolic engine or primary coefficients.
    x,y,z=point;u=1+x*y;v=4+3*x*y
    return [
      [3*u**2*y*z+y**3*(v+3*u),3*u**2*x*z+2*y*u*v+x*y**2*(v+3*u),u**3],
      [3*u**2*z+6*x*y*u*z+3*y**2*(v+3*x*y),1+6*x**2*u*z+6*x*y*v+9*x**2*y**2,3*x*u**2],
      [2-6*x*y-3*x**2*z,-3*x**2,-x**3]]
def value_rows(p):return [dict(exponent=e,coefficient=str(c)) for e,c in sorted(p.p.items())]
def main():
    probe=json.loads((HERE/"nonzero-target-probe.json").read_text())
    assert probe["source_sha256"]==sha(HERE/"nonzero_target_probe.py")
    k=L({1:1})
    a=[(1-k)/2,L(-2),-4*(1+3*k)/k**2]
    b=[(1+k)/2,L(-2),-4*(1-3*k)/k**2]
    target=[L(),L(-2),(1-k**2)/(2*k**2)]
    assert F(a)==F(b)==target
    A=J(a);Ai=inverse(A)
    c=[a[i]+2*(b[i]-a[i]) for i in range(3)]
    axis_trace=sum((matmul(Ai,J(c))[i][i] for i in range(3)),L())
    assert axis_trace==864*k**2-2736*k+723
    delta=[F(c)[i]-target[i] for i in range(3)]
    physical_image=[a[i]+sum((Ai[i][j]*delta[j] for j in range(3)),L()) for i in range(3)]
    trace=sum((matmul(Ai,J(physical_image))[i][i] for i in range(3)),L())
    assert min(trace.p)==-3 and trace.p[-3]==-49460507836416
    invariants=probe["frame_invariants"]
    assert parse(invariants["trace_at_axis2"])==axis_trace
    assert parse(invariants["trace_at_image_axis2"])==trace
    receipt=dict(status="pass",source_sha256=sha(Path(__file__)),
       primary_source_sha256=probe["source_sha256"],probe_sha256=sha(HERE/"nonzero-target-probe.json"),
       method="Independent standard-library Fraction Laurent arithmetic; direct foundational map and explicit derivatives",
       collision_family=True,normalizing_det=-2,
       trace_at_axis2=value_rows(axis_trace),
       trace_at_image_axis2=value_rows(trace),
       zero_valuation=-3,zero_leading_coefficient="-49460507836416",
       infinity_axis_trace_degree=2,infinity_axis_trace_leading_coefficient="864",
       scope="No finite coefficientwise limit at either endpoint after any meromorphic frame fixing e1 for this specified collision pair family; finite-parameter support and other pair families remain open.")
    (HERE/"nonzero-target-validation.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+"\n")
    print(json.dumps({k:v for k,v in receipt.items() if k not in ("trace_at_axis2","trace_at_image_axis2")}))
if __name__=="__main__":main()

