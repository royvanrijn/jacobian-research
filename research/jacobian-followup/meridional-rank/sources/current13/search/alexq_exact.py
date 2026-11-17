import sys,ast,json
sys.path.insert(0,'.')
from wirt import load
from quandle import crossings
from alexq import matrix, snf_diag
def nullmax(n,cr,p):
    d=snf_diag(matrix(n,cr,p),p); free=n-len(d); best=free
    nu=[f for f in d if f.degree()>0]
    if nu:
        for g,_ in nu[-1].factor()[1]:
            best=max(best,free+sum(1 for f in nu if f%g==0))
    return best
if __name__=='__main__':
    names=json.load(open(sys.argv[1])); ps=[p for p in range(2,int(sys.argv[2])) if all(p%q for q in range(2,p))]
    rows,I=load(); S=set(names); hits={}
    for x in rows:
        if x[I['name']] not in S: continue
        pd=ast.literal_eval(x[I['pd_notation']]); n,cr=crossings(pd)
        for p in ps:
            if nullmax(n,cr,p)>=4: hits.setdefault(x[I['name']],[]).append(p)
    print(len(hits),hits)
