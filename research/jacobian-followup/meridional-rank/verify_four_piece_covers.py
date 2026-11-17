#!/usr/bin/env python3
"""Independent amalgam normal forms, ports, unimodular slopes and upper4 replay.

Finite inputs are checked here. Seifert/JSJ realization and the published
rank2 obstruction are explicit uniform imports in the companion proof.
"""
from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
from verify_double_transpositions import geometric_crossings,check_trace

HERE=Path(__file__).resolve().parent


def inverse(w):return [-x for x in w[::-1]]


def stack(w):
    result=[]
    for x in w:
        if result[-1:]==[-x]:result.pop()
        else:result.append(x)
    return result


def subst(w,images):return stack([x for letter in w for x in (images[abs(letter)] if letter>0 else inverse(images[-letter]))])


def cyclic_key(w):
    w=stack(w)
    while w and w[0]==-w[-1]:w=w[1:-1]
    if not w:return ()
    return min(tuple(v[i:]+v[:i]) for v in (w,inverse(w)) for i in range(len(w)))


def power(w,n):return (w if n>=0 else inverse(w))*abs(n)


def H_normal(word):
    """Normal form for H=<a,t | a^3=t^2=h>, central power plus C3*C2 syllables."""
    central=0;syllables=[]
    for letter in word:
        g=abs(letter);assert g in (1,2)
        value=1 if letter>0 else -1
        if syllables and syllables[-1][0]==g:value+=syllables.pop()[1]
        q,r=divmod(value,3 if g==1 else 2);central+=q
        if r:syllables.append((g,r))
    return central,tuple(syllables)


def root_identity(word,d):
    """Identity in H*_(d=s^2)<s>, via cyclic edge-membership and normal form."""
    assert len(word)<=10000
    abel_d=sum((2 if abs(x)==1 else 3)*(1 if x>0 else -1) for x in d)
    assert abel_d in (7,-1)
    def member(w):
        value=sum((2 if abs(x)==1 else 3)*(1 if x>0 else -1) for x in w)
        if value%abel_d:return None
        k=value//abel_d
        assert len(d)*abs(k)<=10000
        return k if H_normal(w)==H_normal(power(d,k)) else None
    pieces=[]
    for x in word:
        kind='S' if abs(x)==3 else 'H';value=(1 if x>0 else -1) if kind=='S' else [x]
        if pieces and pieces[-1][0]==kind:
            pieces[-1][1]=pieces[-1][1]+value
        else:pieces.append([kind,value])
    while True:
        merged=[]
        for kind,value in pieces:
            if (kind=='S' and value==0) or (kind=='H' and H_normal(value)==(0,())):continue
            if merged and merged[-1][0]==kind:merged[-1][1]+=value
            else:merged.append([kind,value])
        pieces=merged;changed=False
        for i,(kind,value) in enumerate(pieces):
            if kind=='S' and value%2==0:
                pieces[i]=['H',power(d,value//2)];changed=True;break
        if changed:continue
        for i in range(1,len(pieces)-1):
            kind,value=pieces[i]
            if kind=='H' and pieces[i-1][0]==pieces[i+1][0]=='S':
                k=member(value)
                if k is not None:pieces[i]=['S',2*k];changed=True;break
        if not changed:break
    if not pieces:return True
    if len(pieces)==1 and pieces[0][0]=='H':return H_normal(pieces[0][1])==(0,())
    return False


def equals(a,b,d):return root_identity(a+inverse(b),d)


def commutator(a,b):return a+b+inverse(a)+inverse(b)


def det(M):return M[0][0]*M[1][1]-M[0][1]*M[1][0]


def check_slope(M,fiber_source,fiber_target):
    assert len(M)==2 and all(len(r)==2 and all(isinstance(x,int) for x in r) for r in M)
    assert abs(det(M))==1
    image=[sum(M[i][j]*fiber_source[j] for j in range(2)) for i in range(2)]
    intersection=abs(image[0]*fiber_target[1]-image[1]*fiber_target[0])
    assert intersection>0
    return intersection


def check_models(result,assembly):
    P,Q=result['models'];assert P['name']=='P' and Q['name']=='Q'
    expectedP=assembly[0]['pieces'][0]['relators'][0]
    expectedQ=[x-2 if x>0 else x+2 for x in assembly[0]['pieces'][1]['relators'][0]]
    assert P['source_relator']==expectedP and Q['source_relator']==expectedQ
    forward={int(k):v for k,v in P['free_basis_forward'].items()}
    back={int(k):v for k,v in P['free_basis_inverse'].items()}
    assert forward=={1:[1],2:[-2,1,1]} and back=={1:[1],2:[1,1,-2]}
    assert all(subst(forward[g],back)==[g] and subst(back[g],forward)==[g] for g in (1,2))
    tP=[-1,-1,2,2];tQ=[2,2,1,1]
    assert cyclic_key(subst(expectedP,forward))==cyclic_key(tP+tP+[-1]*3)==cyclic_key(P['power_relator_two_generator'])
    assert cyclic_key(expectedQ)==cyclic_key(tQ+tQ+[-1]*3)==cyclic_key(Q['power_relator_two_generator'])
    assert P['root_word_d']==[1,1,2] and Q['root_word_d']==[2,-1,-1]
    h=[1]*3;eP=[-1,2];eQ=[2,1]
    assert H_normal(P['root_word_d'])==H_normal(eP+h)
    assert H_normal(Q['root_word_d'])==H_normal(eQ+inverse(h))
    for model in (P,Q):
        d=model['root_word_d'];assert H_normal(commutator(h,d))==(0,())
        assert root_identity(commutator(h,[3,3]),d)
        assert not root_identity(commutator(h,[3]),d)
    return P,Q


def check_root(row,assembly,P,Q):
    name=row['name'];epsilon=row['epsilon'];assert epsilon=={'13a_129':1,'13a_135':-1}[name]
    assert row['orbifold_Euler_characteristics']==['-1/6','-1/2','-1/2','-1/6']
    assert row['JSJ_bases_proposed']==['D(2,3)','A(2)','A(2)','D(2,3)']
    assert all(Fraction(x)<0 for x in row['orbifold_Euler_characteristics'])
    Pimages={1:[1],2:[-3,1,1]};Qimages={3:[1],4:[3]}
    h=[1]*3;f=[3,3];aP=[3]+h;aQ=[3]+power(h,epsilon)
    targetP=[inverse(f),aP,stack(aP+inverse(f)),[]]
    targetQ=[[],f,stack(inverse(aQ)+f),aQ]
    assert row['P_port_expected']==targetP and row['Q_port_expected']==targetQ
    assert row['P_outer_basis_a_f']==[aP,f] and row['Q_outer_basis_a_f']==[aQ,f]
    for i,port in enumerate(assembly['ports']):
        oldP=subst(port['inner_word'],Pimages);oldQ=subst(port['outer_word'],Qimages)
        assert row['P_ports_in_model'][i]=={'edge':port['edge'],'word':oldP}
        assert row['Q_ports_in_model'][i]=={'edge':port['edge'],'word':oldQ}
        assert equals(oldP,targetP[i],P['root_word_d'])
        assert equals(oldQ,targetQ[i],Q['root_word_d'])
    assert root_identity(commutator(aP,f),P['root_word_d'])
    assert root_identity(commutator(aQ,f),Q['root_word_d'])
    matrices=[[[0,1],[1,1]],[[0,1],[1,1]],[[0,1],[epsilon,-1]]]
    assert row['gluing_matrices']==matrices
    intersections=[check_slope(M,[0,1],[0,1]) for M in matrices]
    assert intersections==row['fiber_intersections']==[1,1,1]
    abel=[[14,-(1-12*epsilon)],[19,-(1+12*epsilon)]]
    assert row['abelian_gluing_matrix']==abel
    order=abs(det(abel));assert order==row['expected_H1_order']=={'13a_129':391,'13a_135':401}[name]
    # Primitive H peripheral pairs: d=e*h forP and d=e*h^-1 forQ.
    # Their integer exponents also reproduce the source homology without Smith code.
    return {'name':name,'epsilon':epsilon,'edge_maps':matrices,'fiber_intersections':intersections,
        'abelian_gluing_determinant':det(abel),'source_H1_order':order,'proposed_JSJ_pieces':4,
        'algebra_and_slope_checks':'PASS'}


def main():
    result=json.loads((HERE/'four_piece_cover_results.json').read_text())
    plan=json.loads((HERE/'four_piece_cover_plan.json').read_text())
    assembly=json.loads((HERE/'port_assembly_results.json').read_text())['targets']
    assert result['status']=='COMPLETE_PERIPHERAL_MODEL_INPUTS'
    for path,digest in result['input_sha256'].items():assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest,path
    assert len(result['roots'])==len(assembly)==2
    P,Q=check_models(result,assembly)
    roots=[check_root(row,data,P,Q) for row,data in zip(result['roots'],assembly)]
    assert [r['name'] for r in roots]==plan['roots']
    diagrams=json.loads((HERE/'mutation_cover_results.json').read_text())['diagrams']
    sources={r['name']:r for r in json.loads((HERE/'joint_a5_results.json').read_text())['results']}
    upper=[]
    for names in plan['blocks'].values():
        for name in names:
            row=sources[name];triples=geometric_crossings(diagrams[name]['pd'])
            check_trace(triples,row['seeds'],row['seed_trace'])
            upper.append({'name':name,'seeds':row['seeds'],'propagation_steps':len(row['seed_trace']),'upper4_input_checked':True})
    assert len(upper)==8 and len({r['name'] for r in upper})==8
    rejected=[]
    def reject(label,fn):
        try:fn()
        except (AssertionError,IndexError,KeyError,ValueError):rejected.append(label)
        else:raise AssertionError('corruption accepted: '+label)
    reject('fiber matching cannot certify4JSJpieces',lambda:check_slope([[1,0],[0,1]],[0,1],[0,1]))
    reject('determinant2 is not a torus isomorphism',lambda:check_slope([[0,1],[2,1]],[0,1],[0,1]))
    bad=deepcopy(result['roots'][0]);bad['P_port_expected'][0][0]*=-1
    reject('changed root port word',lambda:check_root(bad,assembly[0],P,Q))
    bad=deepcopy(result);bad['models'][0]['free_basis_inverse']['2']=[1,1,2]
    reject('wrong free basis inverse',lambda:check_models(bad,assembly))
    bad=deepcopy(result['roots'][0]);bad['gluing_matrices'][2][1][0]*=-1
    reject('wrong endpoint fiber sign',lambda:check_root(bad,assembly[0],P,Q))
    # Positive2piece controls have exact2generator presentations and embed H nonabelian.
    # They show that a pieces list alone, without4 canonical pieces, is not the rank gate.
    controls={'two_piece_P_Q_presentations_generators':[2,2],
        'H_relation_positive':H_normal([2,2,-1,-1,-1])==(0,()),
        'H_noncommuting_positive':H_normal(commutator([1],[2]))!=(0,()),
        'root_relation_P_positive':root_identity([3,3]+inverse(P['root_word_d']),P['root_word_d']),
        'root_relation_Q_positive':root_identity([3,3]+inverse(Q['root_word_d']),Q['root_word_d'])}
    assert all(v for k,v in controls.items() if k!='two_piece_P_Q_presentations_generators')
    files=['verify_four_piece_covers.py','four_piece_cover_plan.json','four_piece_cover_results.json',
        'port_assembly_results.json','port_assembly_verification.json','compact_cover_verification.json',
        'mutation_cover_results.json','mutation_cover_verification.json','joint_a5_results.json','joint_a5_verification.json',
        'verify_double_transpositions.py','parabolic_lifts.py']
    receipt={'status':'PASS independent four-piece peripheral and upper4 inputs','roots':roots,'upper4':upper,
        'controls':controls,'corruptions_rejected':rejected,
        'assurance':'Exact free-basis/cyclic relator identities, complete cyclic-amalgam identity reductions, port/commutation checks and unimodular nonparallel slopes. Existing geometric arc parser used only for legacy labels and upper4 replay, disclosed. Four canonical JSJ pieces, rank>=3 via classical Boileau-Weidmann, and knot equalities require the written manifold realization/JSJ/classification proof; finite checks alone do not establish these uniform imports.',
        'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'four_piece_cover_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':receipt['status'],'roots':len(roots),'upper4':len(upper),'corruptions':len(rejected)}))


if __name__=='__main__':main()
