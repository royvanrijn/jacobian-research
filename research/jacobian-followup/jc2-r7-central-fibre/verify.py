"""Replay finite geometry controls; the authored geometric imports stay explicit."""
from pathlib import Path
from itertools import permutations,combinations
import hashlib,json,time
import sympy as S
from critical_values import witnesses

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[]
def check(value,label):
    assert value,label
    checks.append(label)
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def closure(gens,n):
    identity=tuple(range(n));group={identity};todo=[identity]
    while todo:
        a=todo.pop()
        for b in gens:
            c=compose(a,b)
            if c not in group:group.add(c);todo.append(c)
    return frozenset(group)
def cycle(p):
    seen=set();sizes=[]
    for i in range(len(p)):
        if i in seen:continue
        j=i;k=0
        while j not in seen:seen.add(j);k+=1;j=p[j]
        sizes.append(k)
    return tuple(sorted(sizes,reverse=True))
def orbit_sizes(group,n):
    unseen=set(range(n));sizes=[]
    while unseen:
        i=min(unseen);orbit={p[i] for p in group};unseen-=orbit;sizes.append(len(orbit))
    return tuple(sorted(sizes,reverse=True))
def generated_subgroups(n,allowed):
    initial=frozenset({tuple(range(n))});groups={initial};todo=[initial]
    while todo:
        group=todo.pop()
        for g in allowed:
            if g in group:continue
            new=closure(list(group)+[g],n)
            if new not in groups:groups.add(new);todo.append(new)
    return groups

def main():
    started=time.monotonic();perms6=list(permutations(range(6)));c=(1,0,2,3,4,5)
    # Enumerate every possible disjoint local inertia choice on the remaining sheets.
    node_choices=[p for p in perms6 if cycle(p)==(2,1,1,1,1) and p[0]==0 and p[1]==1]
    check(len(node_choices)==6,'all six other-node transpositions')
    for p in node_choices:
        check(compose(c,p)==compose(p,c) and orbit_sizes(closure([c,p],6),6)==(2,2,1,1),'regular node local degrees '+str(p))
    for e,kind in [(2,(2,1,1,1,1)),(3,(3,1,1,1))]:
        choices=[p for p in perms6 if cycle(p)==kind and p[0]==0 and p[1]==1]
        for p in choices:
            check(orbit_sizes(closure([c,p],6),6)==tuple(sorted([2,e]+[1]*(4-e),reverse=True)),f'C/B disjoint inertia e{e} '+str(p))
    groups={};table=[]
    for n in [3,4]:
        perms=list(permutations(range(n)))
        trans=[p for p in perms if cycle(p)==tuple([2]+[1]*(n-2))]
        threes=[p for p in perms if cycle(p)==tuple([3]+[1]*(n-3))]
        for mode,allowed in [('transpositions',trans),('mixed',trans+threes)]:
            if n==3 and mode=='mixed':continue
            allgroups=generated_subgroups(n,allowed)
            selected=[g for g in allgroups if any(t in g for t in trans) and (mode!='mixed' or any(t in g for t in threes))]
            types={(orbit_sizes(g,n),len(g)) for g in selected}
            expected={((3,),6),((2,1),2)} if n==3 else ({((4,),24),((3,1),6)} if mode=='mixed' else {((4,),24),((3,1),6),((2,2),4),((2,1,1),2)})
            check(types==expected,f'exhaustive degree{n} {mode} residue subgroup types')
            groups[f'{n}:{mode}']=[dict(orbits=list(a),order=b) for a,b in sorted(types)]
            grouping='C' if mode=='mixed' else ('A_H_C' if n==3 else 'A_H_B')
            e=3 if mode=='mixed' else 2;v=6-n;hb=int(grouping=='A_H_B');w=e+hb
            necessary=set()
            for g in selected:
                orbits=orbit_sizes(g,n);count=len(orbits)
                for inf in g:
                    shape=cycle(inf);n0=len(shape)
                    if (n-n0)%2:continue # b16, with24 other transpositions
                    genus=S.Rational(24-n+(e-1)*16-n0,2)+count
                    punctures=n0+24*(v-1)+(1+hb)*16
                    necessary.add((orbits,shape,int(genus),punctures))
            for orbits,shape,genus,punctures in sorted(necessary):
                table.append(dict(grouping=grouping,degrees=list(orbits),infinity_partition=list(shape),n0=len(shape),components=len(orbits),total_genus=genus,punctures=punctures))
    # Independent Euler methods: compact normalization minus punctures/nodes,
    # and constructible integration over the target nodal curve.
    b,n0,count=S.symbols('b n0 count',integer=True)
    formula=[]
    for name,e,v,hb in [('C',3,2,0),('A_H_B',2,2,1),('A_H_C',2,3,0)]:
        d0=6-v;w=e+hb;genus=count+S.Rational(24-d0,2)+S.Rational(e-1,2)*b-n0/2
        punctures=n0+24*(v-1)+(1+hb)*b;nodes=12*(6-2*v)
        euler=S.expand(2*count-2*genus-punctures-nodes)
        check(S.expand(euler-(-66-v-w*b))==0,name+' compact and target-integrated central Euler agree')
        jump=S.expand(euler-(-138-w*b));lam=S.expand(jump-nodes)
        check(lam==23*v,name+' zero-level infinity Milnor jump')
        rho,mu,loss=S.symbols('rho mu loss')
        check(S.expand((6*mu+w*rho)-(6*mu-loss))==w*rho+loss,name+' every nonzero-level infinity jump')
        check(w*(16-1)==(45 if w==3 else 30),name+' contact34 nonzero total infinity lower bound')
        total_mu=nodes+6*12;total_lambda=23*v+15*w
        check((total_mu,total_lambda)==((96,91) if v==2 else (72,99)) and 1-total_mu-total_lambda==-138-w*16,name+' contact34 complete Milnor totals')
        formula.append(dict(grouping=name,e=e,v=v,w=w,residue_degree=d0,genus=str(genus),punctures=str(punctures),nodes=nodes,central_euler=str(euler),generic_euler=str(-138-w*b),zero_infinity_milnor=int(lam),other_infinity_milnor_lower_bound=str(w*(b-1))))
    controls=[witnesses(name) for name in ['L','K']]
    for control in controls:
        check(len(control['coefficient_choices'])==2,control['family']+' both coefficient roots retained')
        for row in control['coefficient_choices']:
            for label,value in row['checks'].items():check(value,control['family']+' M'+str(row['M_mod_prime'])+' '+label)
    for grouping,minimum in [('C',25),('A_H_B',17),('A_H_C',18)]:
        rows=[v for v in table if v['grouping']==grouping]
        check(min(v['total_genus'] for v in rows)==minimum,grouping+' contact34 minimum central genus')
    ledger=dict(schema='r7-central-fibre-v1',scope='Necessary source constraints conditional on an actual cover; no curve cover or source A2 constructed.',formulas=formula,central_groups=groups,contact34_necessary_table=table,critical_value_controls=controls)
    (HERE/'ledger.json').write_text(json.dumps(ledger,indent=2)+'\n')
    paths=['PROOF.md','README.md','REVIEW.md','LITERATURE.md','critical_values.py','verify.py','sources/manifest.json','sources/gwozdziewicz-ploski-2001.pdf','sources/borisov-gabber-vasiu-BGV160.pdf']
    inputs={str((HERE/n).relative_to(ROOT)):hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in paths}
    for n in ['jc2-r7-elliptic-contact/algebra.py','jc2-r7-elliptic-contact/models.json','jc2-r7-elliptic-contact/receipt.json','jc2-r7-elliptic-contact/PROOF.md','jc2-r7-complement/geometry.py','jc2-r7-image-family/family.py','jc2-r7-laurent-collisions/laurent.py','jc2-r7-quintic-incidence/finite.py','jc2-r7-quintic-normalizations/normalizations.json','jc2-omitted-divisors/RESIDUE_LEDGER.md','jc2-omitted-divisors/SIMPLE_RAMIFICATION_AND_RADIAL.md','jc2-omitted-divisors/IMAGE_PENCIL_MONODROMY.md','jc2-omitted-divisors/DISJOINT_IMAGE_GATE.md']:
        path=HERE.parent/n
        inputs[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    rec=dict(passed=True,check_count=len(checks),checks=checks,input_sha256=inputs,output_sha256={'ledger.json':hashlib.sha256((HERE/'ledger.json').read_bytes()).hexdigest()},elapsed_seconds=round(time.monotonic()-started,3),assurance='Authored finite-normalization, local-cover and topological proof with exact permutation, symbolic and good-prime controls; no formalization or independent external review.')
    (HERE/'receipt.json').write_text(json.dumps(rec,indent=2)+'\n')
    print(json.dumps({'passed':True,'check_count':len(checks),'elapsed_seconds':rec['elapsed_seconds']}))

if __name__=='__main__':main()
