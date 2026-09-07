"""Exact MW-lattice obstruction to a positive-rank quadratic source.

One worker, 90 seconds. PARI's full lattice automorphism calculation is
independently bounded by automorphisms of the exact spanning norm4 graph.
The resulting theorem concerns this pointed elliptic fibration only.
"""
import argparse,json,runpy,signal
from hashlib import sha256
from pathlib import Path
import numpy as np
from sage.all import ZZ,matrix,pari,Graph

ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'artifacts/generated-results/elliptic-curves'
PARENT=BASE/'curve302_recovered_mw17_parent_v1.json'
PROOF=BASE/'curve302_recovered_mw17_parent_proof_v1.json'
RECOVERY=BASE/'curve302_det1092_recovery_v1.json'
ENUM=ROOT/'elkies-k3/scripts/certify_curve302_anchor6_triangle_gate.sage'
OUT=BASE/'curve302_parent_quadratic_descent_gate_v1.json'


def build():
    parent=json.loads(PARENT.read_text());proof=json.loads(PROOF.read_text())
    assert proof['status']=='PASS_FULL_ARITHMETIC_MW17_PARENT'
    for path in [PARENT,RECOVERY]:
        assert proof['input_sha256'][str(path.relative_to(ROOT))]==sha256(path.read_bytes()).hexdigest()
    G=matrix(ZZ,json.loads(RECOVERY.read_text())['recovery']['gram'])
    W=matrix(ZZ,parent['basis_words_in_recovered_core'])
    assert abs(W.det())==1 and W*G*W.transpose()==matrix(ZZ,parent['generic_height_gram'])
    assert G.is_positive_definite() and G.det()==1092
    aut=pari(G).qfauto();assert aut[0]==2
    generators=[matrix(ZZ,x) for x in aut[1]]
    assert generators and all(M.transpose()*G*M==G for M in generators)
    assert all(M in [matrix.identity(ZZ,17),-matrix.identity(ZZ,17)] for M in generators)
    shell,nodes=runpy.run_path(str(ENUM))['exact_shell'](G,4)
    assert len(shell)==2436 and matrix(ZZ,shell).rank()==17
    # Any lattice isometry acts faithfully on this spanning shell and
    # preserves the graph relation inner product2. Thus graph order2 is an
    # independent upper bound; both +/-identity are actual isometries.
    v=np.array(shell,dtype=np.int64);g=np.array(G.rows(),dtype=np.int64)
    assert 17**2*int(abs(v).max())**2*int(abs(g).max())<2**63
    pairs=v@g@v.T;i,j=np.where(np.triu(pairs,1)==2)
    graph=Graph(len(shell));graph.add_edges(zip(map(int,i),map(int,j)))
    assert graph.size()==158496 and graph.automorphism_group().order()==2
    paths=[Path(__file__),PARENT,PROOF,RECOVERY,ENUM]
    return {'schema':'curve302.quadratic-descent-gate.v1','status':'PASS',
            'sources':{str(path.relative_to(ROOT)):sha256(path.read_bytes()).hexdigest() for path in paths},
            'rank':17,'height_determinant':1092,'lattice_automorphism_order':2,
            'lattice_automorphisms':['I','-I'],'independent_graph_vertices':2436,
            'independent_graph_edges':158496,'independent_graph_automorphism_order':2,
            'spanning_shell_rank':17,'exact_shell_nodes':nodes,
            'possible_quadratic_descent_arithmetic_ranks':[0,17],
            'excluded_quadratic_descent_arithmetic_ranks':list(range(1,17)),
            'scope':'For a quadratic Q(u)/Q(s) descent of this pointed elliptic fibration, the deck involution is an integral height isometry. It is +/-I, so its invariant MW space has rank17 or0. This excludes every source with arithmetic generic rank1..16, in particular an MW8 rational elliptic source.',
            'boundary':'Does not exclude a rank-zero source, other fibrations of this K3, higher-degree covers, multiparameter constructions or identify original discoverer provenance.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--build',action='store_true');args=parser.parse_args();signal.alarm(90)
    result=build()
    if args.build:
        assert not OUT.exists();OUT.write_text(json.dumps(result,indent=2)+'\n')
    assert result==json.loads(OUT.read_text())
    print('PASS Aut(MW17)={+I,-I}; independent graph replay; positive-rank quadratic source1..16 excluded')
