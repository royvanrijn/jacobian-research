"""Alexander-quandle colorings of the 293 open 13-crossing knots over F_p[t]/(g).
Nullity of the n x n coloring matrix at an irreducible g = number of invariant
factors (Smith form over F_p[t]) divisible by g.  Nullity >= 4 gives meridional
rank >= 4, hence bridge index 4."""
import sys, json, ast
sys.path.insert(0,'.')
from wirt import load
from quandle import crossings
from flint import nmod_poly
def matrix(n,cr,p):
    t=nmod_poly([0,1],p); one=nmod_poly([1],p); z=nmod_poly([],p)
    M=[]
    for o,u,v,e in cr:
        r=[z]*n
        def add(i,c):
            r[i]=r[i]+c
        if e==1:  # v = t u + (1-t) o
            add(v,one); add(u,-t); add(o,t-one)
        else:     # u = t v + (1-t) o
            add(u,one); add(v,-t); add(o,t-one)
        M.append(r)
    return M
def snf_diag(M,p):
    M=[row[:] for row in M]; m=len(M); n=len(M[0]); d=[]
    r=0
    for c in range(n):
        if r>=m: break
        while True:
            piv=None
            for i in range(r,m):
                for j in range(c,n):
                    if M[i][j]!=0 and (piv is None or M[i][j].degree()<M[piv[0]][piv[1]].degree()):
                        piv=(i,j)
            if piv is None: return d
            i,j=piv; M[r],M[i]=M[i],M[r]
            for row in M: row[c],row[j]=row[j],row[c]
            a=M[r][c]; done=True
            for i in range(r+1,m):
                if M[i][c]!=0:
                    q=M[i][c]//a
                    M[i]=[x-q*y for x,y in zip(M[i],M[r])]
                    if M[i][c]!=0: done=False
            for j in range(c+1,n):
                if M[r][j]!=0:
                    q=M[r][j]//a
                    for row in M: row[j]=row[j]-q*row[c]
                    if M[r][j]!=0: done=False
            if done:
                # divisibility: ensure a divides rest
                bad=None
                for i in range(r+1,m):
                    for j in range(c+1,n):
                        if M[i][j]%a!=0: bad=(i,j);break
                    if bad: break
                if bad:
                    M[r]=[x+y for x,y in zip(M[r],M[bad[0]])]; continue
                d.append(a); r+=1; break
    return d
if __name__=='__main__':
    todo=json.load(open('../data/still_open.json')); ps=[int(x) for x in sys.argv[1].split(',')]
    rows,I=load(); S=set(todo); sel=[x for x in rows if x[I['name']] in S]
    res={}
    for x in sel:
        name=x[I['name']]; pd=ast.literal_eval(x[I['pd_notation']]); n,cr=crossings(pd)
        for p in ps:
            d=snf_diag(matrix(n,cr,p),p)
            nonunit=[f for f in d if f.degree()>0]
            nz=len(d)  # rank over F_p(t)
            # nullity at g = (n-nz) + #{f in d: g|f}; need >=4
            k=n-nz+len(nonunit)
            if k>=4:
                # max multiplicity count of a common factor
                g=nonunit[-3] if len(nonunit)>=3 and n-nz==1 else None
                res.setdefault(name,[]).append({'p':p,'invariants':[str(f) for f in nonunit],'corank':n-nz})
        print(name, res.get(name,''), flush=True)
    json.dump(res,open('alexq_%s.json'%sys.argv[1].replace(',','_'),'w'),indent=1)
    print('DONE',len(sel),len(res))
