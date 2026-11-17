"""Small exact polynomial utilities used only by independent weighted-map audits."""
from fractions import Fraction as Q
from collections import Counter
import ast
from verify import add,mul
from verify_power_single_compress import decode,scale

def power(p,k,n):
 out={(0,)*n:Q(1)}
 for i in range(k):out=mul(out,p)
 return out

def substitution(p,values,n):
 out={};cache={}
 for e,c in p.items():
  term={(0,)*n:c}
  for i,k in enumerate(e):
   if not k:continue
   if (i,k) not in cache:cache[i,k]=power(values[i],k,n)
   term=mul(term,cache[i,k])
  out=add(out,term)
 return out

def expression(text,names):
 n=len(names);one={(0,)*n:Q(1)};lookup={name:{tuple(int(i==j) for i in range(n)):Q(1)} for j,name in enumerate(names)}
 def visit(node):
  if isinstance(node,ast.Constant):assert type(node.value) is int;return scale(one,Q(node.value))
  if isinstance(node,ast.Name):return lookup[node.id]
  if isinstance(node,ast.UnaryOp):return scale(visit(node.operand),Q(-1) if isinstance(node.op,ast.USub) else Q(1))
  assert isinstance(node,ast.BinOp)
  if isinstance(node.op,ast.Add):return add(visit(node.left),visit(node.right))
  if isinstance(node.op,ast.Sub):return add(visit(node.left),visit(node.right),Q(-1))
  if isinstance(node.op,ast.Mult):return mul(visit(node.left),visit(node.right))
  if isinstance(node.op,ast.Div):
   denominator=visit(node.right);assert len(denominator)==1 and (0,)*n in denominator
   return scale(visit(node.left),1/denominator[(0,)*n])
  assert isinstance(node.op,ast.Pow) and isinstance(node.right,ast.Constant) and type(node.right.value) is int and node.right.value>=0
  return power(visit(node.left),node.right.value,n)
 return visit(ast.parse(text.replace('^','**'),mode='eval').body)

def normal(p):
 if not p:return ()
 return tuple(sorted(scale(p,1/p[max(p)]).items()))

def equal_generators(a,b):
 assert {normal(p) for p in a if p}=={normal(p) for p in b if p}

def rank(rows,width):
 A=[[Q(c) for c in row] for row in rows];assert all(len(row)==width for row in A);r=0
 for j in range(width):
  k=next((k for k in range(r,len(A)) if A[k][j]),None)
  if k is None:continue
  A[r],A[k]=A[k],A[r];c=A[r][j];A[r]=[v/c for v in A[r]]
  for k in range(r+1,len(A)):
   if A[k][j]:
    c=A[k][j];A[k]=[v-c*w for v,w in zip(A[k],A[r])]
  r+=1
  if r==len(A):break
 return r

def coefficient_buckets(poly,source_vars):
 out={}
 for e,c in poly.items():out.setdefault(e[:source_vars],{})[e[source_vars:]]=c
 return {e:p for e,p in out.items() if p}
