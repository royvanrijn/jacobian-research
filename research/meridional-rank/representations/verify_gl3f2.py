#!/usr/bin/env python3
"""Independent projective group, orbit marking, diagram closure and modules.

Does not import the producer's group, orbit, core-extension, bit-kernel,
Steinberg, Fox-prefix or reflection constructors.
"""
from collections import Counter
from itertools import product
import hashlib
import json
from pathlib import Path
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent


def perm_mul(a,b):return tuple(a[b[i]] for i in range(len(a)))
def perm_inv(a):return tuple(a.index(i) for i in range(len(a)))
def mobius(a):
    x,y,z,w=a;out=[]
    for t in range(8):
        num,den=(x,z) if t==7 else ((x*t+y)%7,(z*t+w)%7)
        out.append(num*pow(den,-1,7)%7 if den else 7)
    return tuple(out)


def projective():
    matrices={}
    for a in product(range(7),repeat=4):
        if (a[0]*a[3]-a[1]*a[2])%7:matrices.setdefault(mobius(a),a)
    full=sorted(matrices);unit=tuple(range(8));full.remove(unit);full.insert(0,unit)
    elements=[a for a in full if sum(a[i]>a[j] for i in range(8) for j in range(i+1,8))%2==0]
    assert len(full)==336 and len(elements)==168
    index={a:i for i,a in enumerate(elements)}
    table=[[index[perm_mul(a,b)] for b in elements] for a in elements]
    inv=[index[perm_inv(a)] for a in elements]
    orders=[]
    for i in range(168):
        g=i;n=1
        while g:g=table[g][i];n+=1
        orders.append(n)
    autos=[tuple(index[perm_mul(perm_mul(a,b),perm_inv(a))] for b in elements) for a in full]
    assert len(set(autos))==336
    return elements,table,inv,orders,autos,matrices


def generate(values,table):
    found={0};queue=[0]
    for a in queue:
        for b in values:
            z=table[a][b]
            if z not in found:found.add(z);queue.append(z)
    return frozenset(found)


def marked_orbits(arity,table,autos):
    count=168**arity;marked=bytearray(count);result=[]
    for code in range(count):
        if marked[code]:continue
        x=code;values=[]
        for _ in range(arity):values.append(x%168);x//=168
        values=tuple(reversed(values))
        for a in autos:
            c=0
            for v in values:c=168*c+a[v]
            marked[c]=1
        if len(generate(values,table))==168:result.append(values)
    assert all(marked)
    return result


def matmul(a,b,p):
    columns=list(zip(*b));return [[sum(x*y for x,y in zip(row,col))%p for col in columns] for row in a]
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def reduce_rows(rows,columns,p):
    a=[list(r) for r in rows];k=0;piv=[]
    for j in range(columns):
        i=next((i for i in range(k,len(a)) if a[i][j]%p),None)
        if i is None:continue
        a[k],a[i]=a[i],a[k];c=pow(a[k][j],-1,p);a[k]=[x*c%p for x in a[k]]
        for i in range(len(a)):
            if i!=k:
                c=a[i][j];a[i]=[(x-c*y)%p for x,y in zip(a[i],a[k])]
        piv.append(j);k+=1
        if k==len(a):break
    return a[:k],piv
def rank(rows,columns,p):return len(reduce_rows(rows,columns,p)[1])
def kernel(rows,columns,p):
    reduced,piv=reduce_rows(rows,columns,p);out=[]
    for j in range(columns):
        if j in piv:continue
        v=[0]*columns;v[j]=1
        for i,k in enumerate(piv):v[k]=-reduced[i][j]%p
        out.append(v)
    return out
def invmat(a,p):
    n=len(a);a,piv=reduce_rows([list(r)+e for r,e in zip(a,eye(n))],2*n,p)
    assert piv==list(range(n));return [r[n:] for r in a]
def unpack(a,n):return [[(r>>j)&1 for j in range(n)] for r in a]


def identify(catalogue,table,orders):
    gl=[unpack(a,3) for a in catalogue['matrices']]
    # Independently enumerate every dense invertible binary matrix.
    independently=[]
    for values in product(range(2),repeat=9):
        a=[list(values[3*i:3*i+3]) for i in range(3)]
        if rank(a,3,2)==3:independently.append(tuple(sum(x<<j for j,x in enumerate(r)) for r in a))
    independently.sort();independently.remove((1,2,4));independently.insert(0,(1,2,4))
    assert independently==[tuple(a) for a in catalogue['matrices']]
    idx={tuple(x for r in a for x in r):i for i,a in enumerate(gl)}
    oldtable=[[idx[tuple(x for r in matmul(a,b,2) for x in r)] for b in gl] for a in gl]
    a,b=catalogue['std_pair'];isom=None
    for u in range(168):
        if orders[u]!=2:continue
        for v in range(168):
            if orders[v]!=3 or orders[table[u][v]]!=7:continue
            mapping={0:0};queue=[0];failed=False
            for g in queue:
                for x,y in [(a,u),(b,v)]:
                    target,image=oldtable[g][x],table[mapping[g]][y]
                    if target in mapping:
                        if mapping[target]!=image:failed=True;break
                    else:mapping[target]=image;queue.append(target)
                if failed:break
            if not failed and len(mapping)==len(set(mapping.values()))==168:isom=[mapping[i] for i in range(168)];break
        if isom is not None:break
    assert isom is not None
    assert all(isom[oldtable[i][j]]==table[isom[i]][isom[j]] for i in range(168) for j in range(168))
    return gl,oldtable,isom


def adjoint8(gl):
    off=[(i,j) for i in range(3) for j in range(3) if i!=j]
    basis=[]
    for i,j in off:
        a=[[0]*3 for _ in range(3)];a[i][j]=1;basis.append(a)
    basis.extend([[[int(i==j and i in [0,1]) for j in range(3)] for i in range(3)],
                  [[int(i==j and i in [1,2]) for j in range(3)] for i in range(3)]])
    result=[]
    for a in gl:
        inv=invmat(a,2);columns=[]
        for b in basis:
            x=matmul(matmul(a,b,2),inv,2);columns.append([x[i][j] for i,j in off]+[x[0][0],x[2][2]])
        result.append(list(map(list,zip(*columns))))
    return result


def module_isomorphism(ad,building,std):
    rows=[]
    for g in std:
        for i in range(8):
            for j in range(8):
                row=[0]*64
                for k in range(8):row[8*k+j]^=ad[g][i][k];row[8*i+k]^=building[g][k][j]
                rows.append(row)
    basis=kernel(rows,64,2);assert len(basis)==1
    x=[basis[0][8*i:8*i+8] for i in range(8)];assert rank(x,8,2)==8
    assert all(matmul(a,x,2)==matmul(x,b,2) for a,b in zip(ad,building))
    return x


def affine_rows(images,inverses,words,p):
    d=len(images[0]);out=[]
    for word in words:
        block=[[0]*(d*len(images)) for _ in range(d)]
        for variable in range(d*len(images)):
            point=[0]*d
            for letter in reversed(word):
                i=abs(letter)-1;a=images[i] if letter>0 else inverses[i]
                point=[sum(x*y for x,y in zip(r,point))%p for r in a]
                if i==variable//d:
                    point=[(x+(int(j==variable%d) if letter>0 else -inverses[i][j][variable%d]))%p for j,x in enumerate(point)]
            for j,x in enumerate(point):block[j][variable]=x
        out.extend(block)
    return out


def check_cocycles(test,images,inverses,words):
    d=len(images[0]);columns=d*len(images);rows=affine_rows(images,inverses,words,2)
    assert rows==[[(r>>j)&1 for j in range(columns)] for r in test['rows']]
    basis=[[(v>>j)&1 for j in range(columns)] for v in test['kernel_basis']]
    assert all(not any(sum(x*y for x,y in zip(r,v))%2 for r in rows) for v in basis)
    z=columns-rank(rows,columns,2);assert rank(basis,columns,2)==z==test['Z1_dimension']
    principal=[[a[i][j]^int(i==j) for j in range(d)] for a in images for i in range(d)]
    invariants=d-rank(principal,d,2)
    assert invariants==test['invariants']==0 and z-d==test['H1_dimension']
    return z


def check_binary_roots(test,images,crossings):
    roots=[[(r>>j)&1 for j in range(3)] for r in test['roots']]
    rows=[]
    for a,r in zip(images,roots):
        assert matmul(a,a,2)==eye(3)
        assert rank([[x^int(i==j) for j,x in enumerate(rr)] for i,rr in enumerate(a)],3,2)==1
        assert [sum(x*y for x,y in zip(rr,r))%2 for rr in a]==r
    for o,a,b,sign in crossings:
        word=[-(o+1),a+1,o+1,-(b+1)] if sign else [o+1,a+1,-(o+1),-(b+1)]
        block=[[0]*13 for _ in range(3)]
        for variable in range(13):
            point=[0]*3
            for l in reversed(word):
                j=abs(l)-1;point=[sum(x*y for x,y in zip(rr,point))%2 for rr in images[j]]
                if j==variable:point=[x^y for x,y in zip(point,roots[j])]
            for j,x in enumerate(point):block[j][variable]=x
        rows.extend(block)
    assert rows==[[(r>>j)&1 for j in range(13)] for r in test['rows']]
    basis=[[(v>>j)&1 for j in range(13)] for v in test['kernel_basis']]
    assert rank(rows,13,2)==10 and rank(basis,13,2)==test['dimension']==3
    assert all(sum(x*y for x,y in zip(r,v))%2==0 for r in rows for v in basis)


def algebra_rank(generators):
    d=len(generators[0]);basis={};queue=[eye(d)]
    def insert(a):
        v=[x for row in a for x in row]
        for j,row in sorted(basis.items()):
            if v[j]:v=[x^y for x,y in zip(v,row)]
        j=next((j for j,x in enumerate(v) if x),None)
        if j is None:return False
        basis[j]=v;return True
    insert(queue[0])
    for a in queue:
        for b in generators:
            c=matmul(a,b,2)
            if insert(c):queue.append(c)
    return len(basis)


def torus_control(control):
    x,y=[unpack(control[k],4) for k in ['homogeneous_X','homogeneous_Y']]
    def power(a,n):
        out=eye(4)
        for _ in range(n):out=matmul(out,a,2)
        return out
    assert power(x,4)==power(y,7)==eye(4)
    m=matmul(invmat(x,2),power(y,2),2)
    assert m==unpack(control['homogeneous_meridian'],4)
    assert rank([[v^int(i==j) for j,v in enumerate(r)] for i,r in enumerate(m)],4,2)==1
    assert rank([[v^int(i==j) for j,v in enumerate(r)] for a in [x,y] for i,r in enumerate(a)],4,2)==4
    def action(a):
        out=[]
        for v in range(8):
            point=[(v>>j)&1 for j in range(3)]+[1]
            image=[sum(c*t for c,t in zip(r,point))%2 for r in a]
            assert image[3]==1;out.append(sum(c<<j for j,c in enumerate(image[:3])))
        return tuple(out)
    generators=[action(a) for a in [x,y]];unit=tuple(range(8));seen={unit};queue=[unit]
    for a in queue:
        for b in generators:
            c=perm_mul(a,b)
            if c not in seen:seen.add(c);queue.append(c)
    assert len(seen)==control['affine_group_order']==1344
    principal=[[(v>>j)&1 for j in range(6)] for v in control['principal_basis']]
    basis=[[(v>>j)&1 for j in range(6)] for v in control['kernel_basis']]
    constraints=[[(v>>j)&1 for j in range(6)] for v in control['admissibility_rows']]
    assert rank(principal,6,2)==3 and rank(basis,6,2)==4 and rank(constraints,6,2)==2
    assert all(sum(a*b for a,b in zip(r,v))%2==0 for r in constraints for v in basis)
    # Every homogeneous meridian image has direction in the first3 coordinates.
    # Select four conjugate meridians with independent fixed-hyperplane equations.
    selected=[];functionals=[];directions=[]
    for h in queue:
        image=perm_mul(perm_mul(h,action(m)),perm_inv(h))
        # Recover its affine matrix from images of0 and the three basis points.
        translation=image[0]
        linear_columns=[image[1<<j]^translation for j in range(3)]
        a=[[(linear_columns[j]>>i)&1 for j in range(3)]+[(translation>>i)&1] for i in range(3)]+[[0,0,0,1]]
        difference=[[v^int(i==j) for j,v in enumerate(r)] for i,r in enumerate(a)]
        f=next(r for r in difference if any(r))
        if rank(functionals+[f],4,2)>len(functionals):
            functionals.append(f);pivot=next(j for j,v in enumerate(f) if v)
            directions.append([difference[i][pivot] for i in range(4)])
            selected.append(a)
            if len(selected)==4:break
    assert len(selected)==4 and all(u[3]==0 for u in directions)
    gram=[[sum(a*b for a,b in zip(f,u))%2 for u in directions] for f in functionals]
    assert rank(gram,4,2)<=3
    return {'homogeneous_X':x,'homogeneous_Y':y,'meridian':m,'affine_group_order':1344,
            'root_lift_dimension':4,'global_fixed_codimension':4,'selected_conjugate_meridians':selected,
            'directions':directions,'functionals':functionals,'seed_pairing':gram,'seed_pairing_rank':rank(gram,4,2)}


def extend(values,seeds,crossings,core):
    labels={s:v for s,v in zip(seeds,(0,)+tuple(values))}
    changed=True
    while changed:
        changed=False
        for o,a,b,_ in crossings:
            if o not in labels:continue
            if a in labels:
                value=core[labels[o]][labels[a]]
                if b in labels:
                    if labels[b]!=value:return None
                else:labels[b]=value;changed=True
            elif b in labels:
                labels[a]=core[labels[o]][labels[b]];changed=True
    assert len(labels)==13
    return [labels[i] for i in range(13)]


def check_reflections(report,panel):
    # Compute Ad from its action on three basis matrices, without the producer formula.
    elements=report['projective_elements'];matrices=[]
    basis=[[[1,0],[0,6]],[[0,1],[0,0]],[[0,0],[1,0]]]
    for a in elements:
        dense=[a[:2],a[2:]];inv=invmat(dense,7)
        columns=[]
        for b in basis:
            x=matmul(matmul(dense,b,7),inv,7);columns.append([-x[0][0]%7,-x[0][1]%7,-x[1][0]%7])
        matrices.append(list(map(list,zip(*columns))))
    results=[]
    for row in report['results']:
        images=[matrices[i] for i in row['projective_meridians']]
        assert images==[[m[3*i:3*i+3] for i in range(3)] for m in row['meridian_matrices']]
        roots=row['roots'];crossings=oriented_diagram(panel[row['name']]['pd'])
        for a,r in zip(images,roots):
            assert matmul(a,a,7)==eye(3)
            assert rank([[(x-int(i==j))%7 for j,x in enumerate(rr)] for i,rr in enumerate(a)],3,7)==1
            assert [sum(x*y for x,y in zip(rr,r))%7 for rr in a]==[-x%7 for x in r]
        rows=[]
        for o,a,b,sign in crossings:
            word=[-(o+1),a+1,o+1,-(b+1)] if sign else [o+1,a+1,-(o+1),-(b+1)]
            product_matrix=eye(3)
            for l in word:product_matrix=matmul(product_matrix,images[abs(l)-1],7)
            assert product_matrix==eye(3)
            block=[[0]*13 for _ in range(3)]
            for variable in range(13):
                point=[0]*3
                for l in reversed(word):
                    j=abs(l)-1;point=[sum(x*y for x,y in zip(rr,point))%7 for rr in images[j]]
                    if j==variable:point=[(x+y)%7 for x,y in zip(point,roots[j])]
                for j,x in enumerate(point):block[j][variable]=x
            rows.extend(block)
        assert rows==row['matrix']
        basis=row['kernel_basis'];assert rank(rows,13,7)==10 and rank(basis,13,7)==3
        assert all(sum(x*y for x,y in zip(r,v))%7==0 for r in rows for v in basis)
        results.append({'name':row['name'],'orbit_index':row['orbit_index'],'dimension':3})
    return results


def main():
    catalogue=json.loads((HERE/'gl3f2_catalogue.json').read_text())
    source=json.loads((HERE/'gl3f2_full_panel_results.json').read_text())
    reflections=json.loads((HERE/'reflection7_full_panel_results.json').read_text())
    assert source['status']==reflections['status']=='COMPLETE'
    for report in [source,reflections]:
        for p,h in report['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    panel=json.loads((HERE/'joint_a5_panel.json').read_text())['cases']
    psl,table,inverses,orders,autos,pgl_matrices=projective()
    pgl=sorted(pgl_matrices);pgl.remove(tuple(range(8)));pgl.insert(0,tuple(range(8)))
    gl,oldtable,isom=identify(catalogue,table,orders);back={v:i for i,v in enumerate(isom)}
    transported=[tuple(back[a[isom[i]]] for i in range(168)) for a in autos]
    assert set(transported)=={tuple(a) for a in catalogue['automorphisms']}
    pairs=marked_orbits(2,table,autos);triples=marked_orbits(3,table,autos)
    assert len(pairs)==catalogue['pair_capacity']==57 and len(triples)==13368
    unseen=set(range(168));regular_classes=[]
    while unseen:
        g=min(unseen)
        conjugates={table[table[h][g]][inverses[h]] for h in range(168)}
        unseen-=conjugates
        if orders[g]%2:regular_classes.append({'order':orders[g],'size':len(conjugates)})
    assert sorted((c['order'],c['size']) for c in regular_classes)==[(1,1),(3,56),(7,24),(7,24)]
    orbit_values=[]
    for values in triples:
        old=tuple(back[i] for i in values)
        action=min(transported,key=lambda a:tuple(a[i] for i in old))
        orbit_values.append((values,tuple(action[i] for i in old),action))
    assert {v for _,v,_ in orbit_values}=={tuple(v) for v in catalogue['generating_triple_orbits']}
    print('PASS independent PSL2(F7) group and marked tuple-orbit catalogue',flush=True)
    core=[[table[table[b][inverses[a]]][b] for a in range(168)] for b in range(168)]
    results=[]
    for case,row in zip(panel,source['results']):
        assert case['name']==row['name'];crossings=oriented_diagram(case['pd']);found={}
        for values,key,action in orbit_values:
            labels=extend(values,row['seeds'],crossings,core)
            if labels is not None:found[key]=[action[back[i]] for i in labels]
        expected={tuple(m['seed_images']):m['arc_images'] for m in row['maps']}
        assert found==expected and len(found)==row['GL3_quotient_kernels']
        results.append({'name':row['name'],'kernels':len(found)})
        if len(results)%50==0:print('PASS independent diagram closure',len(results),flush=True)
    first=json.loads((HERE/'gl3f2_results.json').read_text())
    building=[unpack(a,8) for a in first['Steinberg_matrices']]
    ad=adjoint8(gl);intertwiner=module_isomorphism(ad,building,catalogue['std_pair'])
    assert algebra_rank([gl[i] for i in catalogue['std_pair']])==9
    assert algebra_rank([building[i] for i in catalogue['std_pair']])==64
    natural=gl;dual=[list(map(list,zip(*invmat(a,2)))) for a in gl]
    modules=[];maps_done=0;root_checks=0;deck_types=Counter()
    for row,case in zip(source['results'],panel):
        crossings=oriented_diagram(case['pd']);words=[[o+1,-(a+1),o+1,-(b+1)] for o,a,b,_ in crossings]+[[row['seeds'][0]+1]]
        trivial=[[sum(l==j+1 or l==-(j+1) for l in w)%2 for j in range(13)] for w in words]
        assert rank(trivial,13,2)==13
        for mapping in row['maps']:
            labels=mapping['arc_images'];values=[]
            psl_labels=[isom[i] for i in labels]
            deck=[j for j,a in enumerate(autos) if all(a[i]==inverses[i] for i in psl_labels)]
            assert len(deck)==1
            dtype='inner' if pgl[deck[0]] in psl else 'outer'
            assert mapping['deck_type']==dtype;deck_types[dtype]+=1
            if dtype=='inner':
                c=mapping['deck_element'];assert isom[c]==psl.index(pgl[deck[0]])
                for test,model in zip(mapping['root_lifts'],[natural,dual]):
                    meridians=[oldtable[i][c] for i in labels]
                    if test['module']=='dual':
                        # The stored indices already incorporate inverse-transpose.
                        images=[dual[i] for i in meridians]
                    else:images=[natural[i] for i in meridians]
                    check_binary_roots(test,images,crossings);root_checks+=1
            else:assert not mapping['root_lifts']
            for test,model in zip(mapping['ordinary_cocycles'],[natural,dual,building]):
                images=[model[i] for i in labels];inv=[invmat(a,2) for a in images]
                values.append(check_cocycles(test,images,inv,words))
            assert values[0]==values[1] and values[0] in [4,5] and values[2] in [8,10]
            modules.append({'name':row['name'],'orbit_index':mapping['orbit_index'],'Z1_dimensions':values})
            maps_done+=1
            if maps_done%50==0:print('PASS independent affine module matrices',maps_done,flush=True)
    reflection_results=check_reflections(reflections,{r['name']:r for r in panel})
    print('PASS all191 field7 reflection root matrices',flush=True)
    positive_control=torus_control(source['T47_control'])
    paths=['verify_gl3f2.py','gl3f2_catalogue.json','gl3f2_full_panel_results.json','reflection7_full_panel_results.json',
           'gl3f2_results.json','joint_a5_panel.json','parabolic_lifts.py',
           'sources/current13/search/glsearch.py','sources/current13/search/w3.py']
    receipt={'status':'PASS','targets':293,'PSL_to_GL_identification':{str(k):v for k,v in back.items()},
             'pair_capacity':57,'triple_capacity':13368,'results':results,'module_checks':modules,
             'characteristic2_regular_classes':regular_classes,'absolute_simple_matrix_algebra_dimensions':{'natural':9,'dual':9,'Steinberg':64},
             'Steinberg_to_adjoint_intertwiner':intertwiner,'reflection_checks':reflection_results,
             'binary_root_checks':root_checks,'deck_type_histogram':dict(deck_types),'T47_control':positive_control,
             'kernel_histogram':dict(Counter(r['kernels'] for r in results)),
             'module_Z1_histograms':{n:dict(Counter(r['Z1_dimensions'][j] for r in modules)) for j,n in enumerate(['natural','dual','Steinberg'])},
             'new_sharp_cases':0,'scope':'Independent projective group, full tuple orbit marking, iterative crossing closure, adjoint8 transport and reverse affine module/root matrices. Same source PD parser and input coefficient conventions; no actual rank2 or all-field exclusion.',
             'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'gl3f2_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','targets':293,'quotient_maps':len(modules),'reflection_bases':len(reflection_results),'histogram':receipt['kernel_histogram']},indent=2))


if __name__=='__main__':main()
