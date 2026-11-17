"""Signature-selected discovery of an arithmetic membership DAG.

Only a verified DAG ending in one is a certificate; a capped attempt
saves summary statistics rather than expanded failed proof expressions.
"""
from pathlib import Path
import inspect,json,signal,time,sys
import sympy as S
from sympy.polys.rings import ring
import sympy.polys.groebnertools as G

HERE=Path(__file__).resolve().parent
class Found(Exception):pass
class Capped(Exception):pass
def encode(p):return [[list(m),str(c)] for m,c in sorted(p.items())]


def main():
    radical='--radical' in sys.argv
    source=json.loads((HERE/('universal_r6_input.json' if radical else 'universal_r6_closed_input.json')).read_text())
    R,zw,zB,zc,zd,za=ring(','.join(source['variables']),S.QQ,order='grevlex')
    fs=[R.from_expr(S.sympify(f)) for f in source['generators']]
    nodes=[dict(input=i,polynomial=f) for i,f in enumerate(fs)]
    node_of={f:i for i,f in enumerate(fs)};solution=None
    targets=[f.monic() for f in [zc,6*zB+27*zd-10*za,za*za,zd*zd,90*zw*zw-75*zw+4]] if radical else []
    found={}
    start=time.monotonic()
    def record(p,parents):
        nonlocal solution
        if not p:return p
        if p not in node_of:
            node_of[p]=len(nodes)
            nodes.append(dict(polynomial=p,parents=[(node_of[f],g) for f,g in parents if g]))
            if len(nodes)%100==0:print('nodes',len(nodes),'seconds',round(time.monotonic()-start,2),flush=True)
            if len(nodes)>10000:raise Capped
        if p==R.one:
            solution=node_of[p];raise Found
        if p in targets:
            found[targets.index(p)]=node_of[p]
            print('target',targets.index(p),'node',node_of[p],flush=True)
            if len(found)==len(targets):raise Found
        return p
    def rem(p,divisors):
        if not p:return p
        quotients,r=p.div(divisors)
        return record(r,[(p,R.one)]+[(f,-g) for f,g in zip(divisors,quotients) if g])
    def monic(p):
        if not p:return p
        multiplier=R.ground_new(R.domain.one/p.LC)
        return record(p.monic(),[(p,multiplier)])
    def mul(f,term):
        result=G.lbp_mul_term(f,term)
        record(G.Polyn(result),[(G.Polyn(f),R.from_dict({term[0]:term[1]}))])
        return result
    def sub(f,g):
        result=G.lbp_sub(f,g)
        record(G.Polyn(result),[(G.Polyn(f),R.one),(G.Polyn(g),-R.one)])
        return result
    ns=dict(G.__dict__)
    ns.update(lbp_mul_term=mul,lbp_sub=sub,recorded_remainder=rem,recorded_monic=monic)
    ns['red_groebner']=lambda polynomials,ring:polynomials
    for name in ['s_poly','f5_reduce','_f5b']:
        text=inspect.getsource(getattr(G,name))
        if name=='_f5b':
            text=text.replace('r = p.rem(F[:i])','r = recorded_remainder(p, F[:i])')
            text=text.replace('Polyn(p).monic()','recorded_monic(Polyn(p))')
            text=text.replace('Polyn(g).monic()','recorded_monic(Polyn(g))')
        exec(compile(text,'<recorded installed SymPy F5B>','exec'),ns)
    def cap(*_):raise Capped
    budget=90 if '--extended' in sys.argv else 30
    signal.signal(signal.SIGALRM,cap);signal.alarm(budget)
    status='unfinished'
    try:
        basis=ns['_f5b'](fs,R)
        if radical:
            for i,target in enumerate(targets):
                if i in found:continue
                qs,rr=target.div(basis)
                if not rr:record(target,[(f,q) for f,q in zip(basis,qs) if q])
        status='no-unit'
    except Found:status='radical-targets' if radical else 'unit'
    except Capped:status='capped'
    finally:signal.alarm(0)
    output=dict(status=status,variables=source['variables'],input=source['generators'],
                generated_nodes=len(nodes),seconds=round(time.monotonic()-start,3))
    if solution is not None or found:
        needed=set()
        def visit(i):
            if i in needed:return
            for j,_ in nodes[i].get('parents',[]):visit(j)
            needed.add(i)
        roots=([solution] if solution is not None else [])+list(found.values())
        for i in roots:visit(i)
        order=sorted(needed);renumber={old:new for new,old in enumerate(order)}
        encoded=[]
        for i in order:
            node=nodes[i];new=dict(polynomial=encode(node['polynomial']))
            if 'input' in node:new['input']=node['input']
            else:new['parents']=[dict(node=renumber[j],multiplier=encode(q)) for j,q in node['parents']]
            encoded.append(new)
        output.update(solution=renumber[solution] if solution is not None else None,nodes=encoded,
                      targets=[dict(index=i,node=renumber[n],polynomial=encode(targets[i])) for i,n in sorted(found.items())])
    p=HERE/('universal_r6_radical_trace.json' if radical else 'universal_r6_f5_trace.json');p.write_text(json.dumps(output,separators=(',',':'))+'\n')
    print('status',status,'nodes',len(nodes),'bytes',p.stat().st_size,'seconds',output['seconds'],flush=True)


if __name__=='__main__':main()
