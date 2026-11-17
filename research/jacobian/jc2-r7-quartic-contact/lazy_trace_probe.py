"""Generate a local-arithmetic ideal-membership DAG without expanded cofactors.

Discovery invokes the installed SymPy Buchberger routine with its arithmetic
steps recorded. Correctness of pair selection is not a proof dependency:
the final certificate is checked by local polynomial identities alone.
"""
from pathlib import Path
import inspect,json,signal,time
import sympy as S
from sympy.polys.rings import ring
import sympy.polys.groebnertools as G

HERE=Path(__file__).resolve().parent


def encode(poly):return [[list(m),str(c)] for m,c in sorted(poly.items())]


class Found(Exception):pass
class Capped(Exception):pass


def main():
    source=json.loads((HERE/'universal_r6_closed_input.json').read_text())
    R,*_=ring(','.join(source['variables']),S.QQ,order='grevlex')
    fs=[R.from_expr(S.sympify(x)) for x in source['generators']]
    nodes=[dict(input=i,polynomial=encode(f)) for i,f in enumerate(fs)]
    node_of={f:i for i,f in enumerate(fs)};solution=None
    start=time.monotonic()

    def record(poly,parents):
        nonlocal solution
        if not poly:return poly
        if poly not in node_of:
            node_of[poly]=len(nodes)
            nodes.append(dict(polynomial=encode(poly),parents=[dict(node=node_of[p],multiplier=encode(q)) for p,q in parents if q]))
            if len(nodes)%25==0: print('nodes',len(nodes),'seconds',round(time.monotonic()-start,2),flush=True)
        if poly==R.one:
            solution=node_of[poly];raise Found
        return poly

    def remainder(poly,divisors):
        if not poly:return poly
        quotients,reduced=poly.div(divisors)
        if not reduced:return reduced
        scale=R.domain.one/reduced.LC
        return record(reduced*scale,[(poly,R.ground_new(scale))]+[(g,-q*scale) for g,q in zip(divisors,quotients) if q])

    def spoly(f,g,ring):
        lcm=R.monomial_lcm(f.LM,g.LM)
        first=R.from_dict({R.monomial_div(lcm,f.LM):R.domain.one})
        second=R.from_dict({R.monomial_div(lcm,g.LM):-R.domain.one})
        return record(first*f+second*g,[(f,first),(g,second)])

    text=inspect.getsource(G._buchberger)
    text=text.replace('h = g.rem([ f[j] for j in J ])','h = recorded_remainder(g, [ f[j] for j in J ])')
    text=text.replace('r = p.rem(f[:i])','r = recorded_remainder(p, f[:i])')
    namespace=dict(spoly=spoly,recorded_remainder=remainder)
    exec(compile(text,'<recorded installed SymPy Buchberger>','exec'),namespace)
    def cap(*_):raise Capped
    signal.signal(signal.SIGALRM,cap);signal.alarm(30)
    status='unfinished'
    try:
        namespace['_buchberger'](fs,R);status='no-unit'
    except Found:status='unit'
    except Capped:status='capped'
    finally:signal.alarm(0)
    if solution is not None:
        needed=set()
        def visit(n):
            if n in needed:return
            for parent in nodes[n].get('parents',[]):visit(parent['node'])
            needed.add(n)
        visit(solution)
        ordered=sorted(needed);renumber={old:new for new,old in enumerate(ordered)}
        selected=[]
        for old in ordered:
            node=nodes[old]
            if 'parents' in node:
                node=dict(node,parents=[dict(parent,node=renumber[parent['node']]) for parent in node['parents']])
            selected.append(node)
        solution=renumber[solution];nodes=selected
    output=dict(status=status,variables=source['variables'],nodes=nodes,solution=solution,
                input=source['generators'],seconds=round(time.monotonic()-start,3))
    (HERE/'universal_r6_trace.json').write_text(json.dumps(output,separators=(',',':'))+'\n')
    print('status',status,'retained nodes',len(nodes),'bytes',(HERE/'universal_r6_trace.json').stat().st_size,'seconds',output['seconds'],flush=True)


if __name__=='__main__':main()
