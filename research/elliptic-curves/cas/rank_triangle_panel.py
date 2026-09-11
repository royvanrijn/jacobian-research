#!/usr/bin/env python3
"""Frozen equation-only panel from the live wide-search lower-bound snapshot."""
import argparse
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import time

from research_runtime.store import checkpoint
from research_runtime.supervisor import Limits, run

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
OUT = ART / 'rank_triangle_v1'
BROAD = ROOT / 'artifacts/local/elliptic-curves/broad-rank-v1'
RT = BROAD / 'runtime/research'
GP = Path('/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/gp')


def read(p):
    return json.loads(p.read_text())


def digest(p):
    return sha256(p.read_bytes()).hexdigest()


def save(p, data):
    if p.exists():
        assert read(p) == data, p
    else:
        checkpoint(p, data)


def freeze():
    assert not OUT.exists(), 'preserve the frozen panel'
    rows = read(BROAD / 'queue.json')['rows']
    census = []
    for row in rows:
        states = []
        for p in sorted((RT / 'broad-cases' / row['id']).glob('batch-*/broad-state.json')):
            s = read(p)
            if s.get('packet') and s.get('packet_sha256') and isinstance(s.get('rank'), int):
                states.append((s['rank'], s.get('round', 0), p, s))
        if not states:
            continue
        _, _, path, state = max(states, key=lambda r: (r[0], r[1], str(r[2])))
        census.append({'id': row['id'], 'family': row['family'], 'arm': row['arm'],
                       'parameter': row['parameter'], 'rank_lower_bound': state['rank'],
                       'jump_lower_bound': state['rank']-17, 'calls': state['calls'],
                       'round': state.get('round'), 'status': state['status'],
                       'packet': state['packet'], 'packet_sha256': state['packet_sha256'],
                       'state_path': str(path.relative_to(ROOT)), 'state_sha256': digest(path),
                       'backend': row['backend'], 'model': row['model']})
    selected = []
    for jump in range(9):
        candidates = sorted((r for r in census if r['jump_lower_bound'] == jump),
                            key=lambda r: sha256(('triangle-v1:'+r['id']).encode()).hexdigest())
        chosen = []
        for r in candidates:
            if not chosen or r['family'] != chosen[0]['family']:
                chosen.append(r)
            if len(chosen) == 2:
                break
        if len(chosen) < 2:
            chosen += [r for r in candidates if r not in chosen][:2-len(chosen)]
        selected += chosen
    target = 'b-7bb187bc9254e81c6283'
    if not any(r['id'] == target for r in selected):
        selected.append(next(r for r in census if r['id'] == target))
    bindings = {}
    for r in selected:
        packet_path = RT / r['packet']
        packet = read(packet_path)
        replay_path = packet_path.parent / ('packet-verified.json' if r['backend'] == 'native' else 'verified.json')
        replay = read(replay_path)
        assert digest(packet_path) == r['packet_sha256'] == replay['packet_sha256']
        assert replay['status'] == 'PASS_TWO_FINITE_IMPLEMENTATIONS'
        assert packet['curve'] == r['model'] and packet['rank_lower_bound'] == r['rank_lower_bound']
        bindings[str(packet_path)] = digest(packet_path)
        bindings[str(replay_path)] = digest(replay_path)
    import icarm_curve302 as curve
    selected.append({'id': '302', 'family': 'x1092-original', 'arm': 'reference', 'parameter': '0',
                     'model': list(map(str, curve.GENERAL_WEIERSTRASS_COEFFICIENTS)),
                     'rank_lower_bound': 31, 'jump_lower_bound': 14, 'calls': None, 'round': None})
    OUT.mkdir(parents=True)
    save(OUT / 'population.json', {'captured_unix_time': time.time(), 'selected_slots': len(rows),
                                  'queue_sha256': digest(BROAD / 'queue.json'), 'census': census})
    save(OUT / 'panel.json', {'rows': selected, 'point_replay_bindings': bindings})
    save(OUT / 'protocol.json', {
        'schema': 'rank-triangle.v1', 'producer_sha256': digest(Path(__file__)),
        'population_sha256': digest(OUT / 'population.json'), 'panel_sha256': digest(OUT / 'panel.json'),
        'selection': 'At most two per known-jump stratum0..8, SHA256(triangle-v1:id) order, different families where possible; append named11952 control if absent and302 reference. Never refill after arithmetic failures.',
        'limits': {'field_seconds_per_case': 15, 'bnf_seconds_per_complete_field': 3,
                   'pari_stack_bytes': 2*2**30, 'rss_bytes': 3*2**30, 'workers': 1,
                   'geometry_seconds': 180, 'geometry_traces_per_target': 34,
                   'strict_seconds': 60, 'new_point_searches': 0},
        'boundary': 'Snapshot lower-bound jumps, not exact ranks; adaptive search exposure and parent differ. Equation workers receive no points or ranks. Without g, g+local is symbolic, not a numerical Selmer capacity. Partial factorizations and failed BNF remain UNKNOWN.'})
    print('FROZEN', len(selected), 'rows; population strata', dict(Counter(r['jump_lower_bound'] for r in census)), flush=True)


def program(row, folder):
    curve = ','.join(row['model'])
    return f'''
default(realprecision,80);setrand(1);
main()={{
 my(E,M,v,D,f,nf,S,n,u,L,Q,p,c,fac);
 E=ellinit([{curve}]);v=0;M=ellminimalmodel(E,&v);
 print("MINIMAL_MODEL|",vector(5,i,M[i]));print("TRANSFORM|",v);
 print("DISCRIMINANT|",M.disc);print("C4|",M.c4);
 D=factor(abs(M.disc));for(i=1,matsize(D)[1],if(!isprime(D[i,1]),error("unproved factor")));
 S=vector(matsize(D)[1],i,D[i,1]);addprimes(S);
 for(i=1,#S,print("FACTOR|",[D[i,1],D[i,2]]));
 f=x^3+M.b2*x^2+8*M.b4*x+16*M.b6;
 print("CUBIC|",vector(4,i,polcoef(f,i-1)));
 if(!polisirreducible(f),error("rational2torsion"));
 nf=nfinit([f,concat([2],S)]);if(nfcertify(nf)!=[],error("uncertified order"));
 print("FIELD_DISC|",nf.disc);print("FIELD_INDEX|",nf.index);print("SIGNATURE|",nf.sign);
 n=0;for(i=1,#S,p=S[i];L=elllocalred(M,p);Q=idealprimedec(nf,p);c=if(L[1]==1,D[i,2]%2==0,#Q-1);n+=c;print("LOCAL|",[p,D[i,2],L[1],L[2],#Q,c]);print("SPLIT|",p,"|",vector(#Q,j,[Q[j].e,Q[j].f])));
 Q=idealprimedec(nf,2);print("LOCAL2_DIM|",#Q);
 print("REAL_LOCAL_DIM|",if(nf.sign[1]==3,1,0));
 print("BK_OFFSET|",n+if(M.disc>0,2,1));print("ROOT_NUMBER|",ellrootno(M));
 print("CONDUCTOR|",ellglobalred(M)[1]);
 writebin("{folder / 'nf.bin'}",nf);print("DONE_FIELD|1");
}};
iferr(main(),e,print("FAIL|",e);quit(1));quit(0)
'''


def execute():
    protocol = read(OUT / 'protocol.json')
    assert protocol['producer_sha256'] == digest(Path(__file__))
    assert protocol['panel_sha256'] == digest(OUT / 'panel.json')
    limits = protocol['limits']
    results = []
    for row in read(OUT / 'panel.json')['rows']:
        folder = OUT / 'arithmetic' / row['id']
        folder.mkdir(parents=True, exist_ok=True)
        result_path = folder / 'result.json'
        if result_path.exists():
            results.append(read(result_path));continue
        assert not (folder / 'field.supervisor.json').exists(), 'do not repeat an interrupted attempt'
        gp = folder / 'field.gp'
        gp.write_text(program(row, folder))
        sup = run([str(GP), '-q', '-f', '-s', str(limits['pari_stack_bytes']), str(gp)],
                  limits=Limits(limits['field_seconds_per_case'], limits['rss_bytes'], pari_stack_bytes=limits['pari_stack_bytes']),
                  log_path=folder / 'field.log', checkpoint_path=folder / 'field.supervisor.json', cwd=ROOT)
        log = (folder / 'field.log').read_text()
        record = {'id': row['id'], 'status': 'UNKNOWN_FIELD', 'field_supervision': sup,
                  'field_program_sha256': digest(gp), 'field_log_sha256': digest(folder / 'field.log'),
                  'g_upper': None, 'grh_g': None, 'rank_upper': None}
        if sup['outcome'] == 'completed' and 'DONE_FIELD|1' in log and 'FAIL|' not in log:
            single = dict(line.split('|',1) for line in log.splitlines() if '|' in line)
            local = [json.loads(line.split('|')[1]) for line in log.splitlines() if line.startswith('LOCAL|')]
            record.update(status='PASS_EQUATION_LOCAL_ARITHMETIC', bk_offset=int(single['BK_OFFSET']),
                          field_discriminant=single['FIELD_DISC'], field_index=single['FIELD_INDEX'],
                          signature=json.loads(single['SIGNATURE']), root_number=int(single['ROOT_NUMBER']),
                          conductor=single['CONDUCTOR'], minimal_model=json.loads(single['MINIMAL_MODEL']),
                          cubic=json.loads(single['CUBIC']), local=local,
                          local_2_descent_dimensions={str(r[0]): (int(single['LOCAL2_DIM']) if r[0]==2 else r[4]-1) for r in local},
                          dyadic_local_dimension=int(single['LOCAL2_DIM']), real_local_dimension=int(single['REAL_LOCAL_DIM']),
                          split=[line.split('|')[1:] for line in log.splitlines() if line.startswith('SPLIT|')])
            # A separate short provisional probe; incomplete relation data never
            # supply g. Only successful bnfinit invariants are GRH-conditional.
            bgp = folder / 'bnf.gp'
            bgp.write_text(f'default(realprecision,80);setrand(1);\nmain()={{my(nf,b);nf=read("{folder / "nf.bin"}");b=bnfinit(nf,0);print("CYC|",b.cyc);print("DONE_BNF|1");}};\niferr(main(),e,print("FAIL|",e);quit(1));quit(0)\n')
            bsup = run([str(GP),'-q','-f','-s',str(limits['pari_stack_bytes']),str(bgp)],
                       limits=Limits(limits['bnf_seconds_per_complete_field'],limits['rss_bytes'],pari_stack_bytes=limits['pari_stack_bytes']),
                       log_path=folder/'bnf.log',checkpoint_path=folder/'bnf.supervisor.json',cwd=ROOT)
            blog=(folder/'bnf.log').read_text();record['bnf_supervision']=bsup
            record['bnf_log_sha256']=digest(folder/'bnf.log');record['bnf_program_sha256']=digest(bgp)
            if bsup['outcome']=='completed' and 'DONE_BNF|1' in blog and 'FAIL|' not in blog:
                cyc=json.loads(next(line.split('|',1)[1] for line in blog.splitlines() if line.startswith('CYC|')))
                record['grh_g']=sum(v%2==0 for v in cyc);record['provisional_cyclic_factors']=list(map(str,cyc))
                record['grh_rank_upper']=record['grh_g']+record['bk_offset']
        save(result_path,record);results.append(record)
        print(row['id'],row['jump_lower_bound'],record['status'],record.get('bk_offset'),flush=True)
    save(OUT/'arithmetic.json',{'protocol_sha256':digest(OUT/'protocol.json'),'results':results})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','run']);a=p.parse_args()
    freeze() if a.mode=='freeze' else execute()
