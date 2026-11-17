"""Small exact presentation operations; producers retain every accepted move."""
from collections import Counter
from itertools import product
import time


def inv(w):return [-x for x in w[::-1]]


def free(w):
    out=[]
    for x in w:
        if out and x==-out[-1]:out.pop()
        else:out.append(x)
    return out


def cyclic(w):
    w=free(w)
    while len(w)>1 and w[0]==-w[-1]:w=w[1:-1]
    if not w:return []
    return list(min(tuple(v[i:]+v[:i]) for v in (w,inv(w)) for i in range(len(w))))


def normalize(rels):return [list(w) for w in sorted({tuple(cyclic(w)) for w in rels}-{()})]


def substitute(w,images):return free([x for g in w for x in (images[abs(g)] if g>0 else inv(images[abs(g)]))])


def eliminate(rels,g,index):
    row=rels[index];positions=[i for i,x in enumerate(row) if abs(x)==g];assert len(positions)==1
    k=positions[0];tail=row[k+1:]+row[:k];value=inv(tail) if row[k]>0 else tail
    images={abs(x):[abs(x)] for r in rels for x in r};images[g]=value
    result=normalize([substitute(r,images) for i,r in enumerate(rels) if i!=index])
    return value,result


def autos(gens):
    alphabet=[s*g for g in gens for s in (1,-1)]
    for a in alphabet:
        other=[x for x in alphabet if abs(x)!=abs(a)]
        for flags in product((0,1),repeat=len(other)):
            chosen={a}|{x for x,f in zip(other,flags) if f};images={};back={}
            for g in gens:
                if g==abs(a):images[g]=[g];back[g]=[g]
                else:
                    images[g]=([-a] if -g in chosen else [])+[g]+([a] if g in chosen else [])
                    back[g]=([a] if -g in chosen else [])+[g]+([-a] if g in chosen else [])
            assert all(substitute(images[g],back)==[g] and substitute(back[g],images)==[g] for g in gens)
            yield {'images':images,'inverse_images':back}


def dehn_move(rels,deadline):
    best=None;best_cost=sum(map(len,rels))
    for i,target in enumerate(rels):
        assert time.monotonic()<deadline
        for j,source in enumerate(rels):
            if i==j:continue
            needed=len(source)//2+1
            if needed>len(target):continue
            for sign in (1,-1):
                s=source if sign==1 else inv(source)
                for sr in range(len(s)):
                    S=s[sr:]+s[:sr];prefix=S[:needed]
                    for tr in range(len(target)):
                        T=target[tr:]+target[:tr]
                        if T[:needed]!=prefix:continue
                        matched=needed
                        while matched<min(len(S),len(T)) and S[matched]==T[matched]:matched+=1
                        new=normalize([r if k!=i else inv(S[matched:])+T[matched:] for k,r in enumerate(rels)])
                        cost=sum(map(len,new))
                        if cost<best_cost:
                            best_cost=cost;best=({'kind':'relator_dehn','target':i,'source':j,'source_sign':sign,
                                'source_rotation':sr,'target_rotation':tr,'matched':matched},new)
    return best


def simplify(gens,rels,limits):
    started=time.monotonic();deadline=started+limits['seconds_per_target']
    gens=list(gens);rels=normalize(rels);initial={'generators':gens[:],'relators':rels[:]}
    trace=[];status='LOCAL_MINIMUM';original_words={g:[g] for g in gens};generator_words={g:[g] for g in gens}
    while len(trace)<limits['maximum_steps']:
        if time.monotonic()>=deadline:status='TIME_LIMIT';break
        old_cost=sum(map(len,rels));assert old_cost<=limits['maximum_intermediate_total_length']
        choices=[]
        for i,row in enumerate(rels):
            for g,count in Counter(map(abs,row)).items():
                if count==1:
                    value,new=eliminate(rels,g,i)
                    choices.append((sum(map(len,new)),g,i,value,new))
        if choices:
            _,g,i,value,new=min(choices)
            trace.append({'kind':'eliminate','generator':g,'relator':i,'solved_word':value})
            images={h:[h] for h in gens};images[g]=value
            original_words={h:substitute(w,images) for h,w in original_words.items()}
            gens.remove(g);generator_words.pop(g);rels=new;continue
        try:move=dehn_move(rels,deadline)
        except AssertionError:
            assert time.monotonic()>=deadline;status='TIME_LIMIT';break
        if move:
            step,rels=move;trace.append(step);continue
        if len(gens)>limits['maximum_automorphism_rank']:status='AUTOMORPHISM_RANK_LIMIT';break
        best=None;best_cost=old_cost
        for aut in autos(gens):
            if time.monotonic()>=deadline:status='TIME_LIMIT';break
            new=normalize([substitute(w,aut['images']) for w in rels]);cost=sum(map(len,new))
            if cost<best_cost:best_cost=cost;best=(aut,new)
        if status=='TIME_LIMIT':break
        if best is None:break
        aut,rels=best;trace.append({'kind':'automorphism',**aut})
        original_words={h:substitute(w,aut['images']) for h,w in original_words.items()}
        generator_words={g:substitute(aut['inverse_images'][g],generator_words) for g in gens}
    else:status='STEP_LIMIT'
    return {'initial':initial,'trace':trace,'final':{'generators':gens,'relators':rels},
        'original_generators_in_final':original_words,'final_generators_in_original':generator_words,
        'status':status,'seconds':round(time.monotonic()-started,6)}


def graph_core(graph):
    edges={i:list(e) for i,e in enumerate(graph['edges'])};vertices=set(range(graph['vertices']));trace=[]
    while True:
        groups={}
        for i,e in edges.items():groups.setdefault(tuple(sorted(e)),[]).append(i)
        pair=next(((ends,ids) for ends,ids in sorted(groups.items()) if len(ids)>1),None)
        if pair:
            ends,ids=pair;kind='parallel';step={'kind':kind,'edges':ids,'ends':list(ends)}
        else:
            candidate=None
            for v in sorted(vertices):
                ids=[i for i,e in edges.items() if v in e]
                if len(ids)==2:
                    other=[next(x for x in edges[i] if x!=v) for i in ids]
                    if other[0]!=other[1]:candidate=(v,ids,other);break
            if candidate is None:break
            v,ids,ends=candidate;step={'kind':'series','vertex':v,'edges':ids,'ends':ends};vertices.remove(v)
        new=max(edges)+1;step['new_edge']=new
        for i in ids:edges.pop(i)
        edges[new]=list(ends);trace.append(step)
    assert len(vertices)==4 and len(edges)==6
    assert {tuple(sorted(e)) for e in edges.values()}=={(a,b) for a in vertices for b in vertices if a<b}
    return {'trace':trace,'vertices':sorted(vertices),'edges':edges,'conclusion':'Abstract K4 series/parallel core only; no manifold/tangle decomposition is inferred.'}
