#!/usr/bin/env sage-python
"""All153 old-section pairs, no new lattice or point enumeration.25s cap."""
import hashlib,json,signal,itertools
from pathlib import Path
from sage.all import QQ,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_brauer_section_gate_v3'
def read(p):return json.loads(p.read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def retain(p,d):
    s=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==s
    else:p.write_text(s)
def main():
    paths=[ART/'curve302_recovered_mw17_parent_v1.json',
      ART/'det1092_brauer_section_gate_v1/frame.json',ART/'det1092_brauer_section_gate_v2/frame.json',
      ART/'det1092_seed_local_code_v5/protocol.json',Path(__file__)]
    protocol={'classification':'generic-only rational-section intersection graph',
      'rule':'Use O and all17 original generic sections. Check all153 unordered distinct pairs from the exact Gram. Keep intersection degree one and use deterministic breadth-first traversal from O. No section or point is selected from a specialized outcome.',
      'retained_input_failures':'v1 raw4096-row prefix and v2 norm6-filtered old bisection table both expose no norm6 words; the saved table lists degree-two bisection classes, not the omitted short section classes. Neither is an arithmetic exclusion. This version uses only the original18 curves.',
      'limits':{'wall_seconds':25,'curves':18,'pairs':153,'new_lattice_nodes':0,'new_addresses':0,
        'exceptional_coordinate_inputs':0,'point_searches':0,'Brauer_group_computation':0,
        'Selmer_groups':0,'class_groups':0,'V3_inputs':0,'pilot_changes':0},
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths}}
    retain(OUT/'protocol.json',protocol)
    G=matrix(QQ,read(paths[0])['generic_height_gram']);words=[vector(QQ,17)]+[vector(QQ,[int(i==j) for i in range(17)]) for j in range(17)]
    labels=['O']+[f'S{i}' for i in range(17)];adj=[[] for _ in words];pairs=[]
    for i,j in itertools.combinations(range(18),2):
        w=words[i]-words[j];h=w*G*w;intersection=h/2-2
        assert intersection>=0 and intersection.denominator()==1
        row={'i':i,'j':j,'labels':[labels[i],labels[j]],'difference_height':str(h),'intersection_degree':int(intersection)}
        pairs.append(row)
        if intersection==1:adj[i].append(j);adj[j].append(i)
    queue=[0];seen={0};tree=[]
    for i in queue:
        for j in adj[i]:
            if j in seen:continue
            seen.add(j);queue.append(j);tree.append({'from':labels[i],'to':labels[j],'i':i,'j':j})
    W=matrix(QQ,[words[r['j']]-words[r['i']] for r in tree])
    result={'classification':'exact rational-intersection connectivity certificate',
      'status':'CONNECTED_DEGREE_ONE_SECTION_GRAPH' if len(seen)==18 else 'NOT_CONNECTED',
      'labels':labels,'pairs':pairs,'tree':tree,'visited':queue,
      'norm6_tree_word_determinant':str(W.det()) if len(tree)==17 else None,
      'argument':'S_i.S_j=height(P_i-P_j)/2-2 on the proved24I1 K3. An effective rational intersection divisor of degree one is a rational point. Restrictions of a Brauer class to the rational curves S_i and S_j are constants and agree there. Connectivity makes all old sections have the same constant as O.',
      'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [paths[0],OUT/'protocol.json']}}
    retain(OUT/'frame.json',result)
    print(result['status'],'pairs',len(pairs),'tree',len(tree),'det',result['norm6_tree_word_determinant'],flush=True)
if __name__=='__main__':
    signal.alarm(25);OUT.mkdir(parents=True,exist_ok=True);main()
