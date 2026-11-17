#!/usr/bin/env python3
"""Bounded auxiliary recognition, deliberately not a target rank certificate."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent


def worker(index):
    import regina
    cues=json.loads((HERE/'compact_cover_geometry_cues.json').read_text())['results']
    cue=cues[index];T=regina.SnapPeaTriangulation(cue['triangulation'])
    assert T.countFilledCusps()==1
    closed=T.filledAll();before=closed.size();closed.simplify()
    standard=regina.StandardTriangulation.recognise(closed)
    manifold=standard.manifold() if standard is not None else None
    result={'name':cue['name'],'status':'CUE_ONLY','filled_tetrahedra_before':before,'filled_tetrahedra_after':closed.size(),
        'closed_isosig':closed.isoSig(),'homology':str(closed.homology()),'recognition':str(standard) if standard else None,
        'manifold':str(manifold) if manifold else None,'regina_wheel':importlib.metadata.version('regina')}
    if index:
        row=next(r for r in json.loads((HERE/'compact_cover_results.json').read_text())['targets'] if r['name']==cue['name'])
        gens=row['final']['generators'];words=[''.join(chr(97+gens.index(abs(x))) if x>0 else chr(65+gens.index(abs(x))) for x in w) for w in row['final']['relators']]
        G=regina.GroupPresentation(len(gens),words);mapping=G.simplify()
        result['compact_algebra_simplification']={'generators':G.countGenerators(),'presentation':str(G),'mapping':str(mapping) if mapping else None}
    print(json.dumps(result))


def main():
    plan=json.loads((HERE/'compact_cover_recognition_plan.json').read_text())
    files=['compact_cover_recognition_plan.json','probe_compact_cover_recognition.py']+plan['inputs']
    result={'status':'RUNNING','results':[],'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    output=HERE/'compact_cover_recognition_cues.json'
    for i in range(3):
        try:
            process=subprocess.run([sys.executable,str(Path(__file__).resolve()),'--case',str(i)],capture_output=True,text=True,timeout=plan['limits']['seconds_per_case'])
            assert process.returncode==0,process.stderr
            row=json.loads(process.stdout)
        except subprocess.TimeoutExpired:row={'case':i,'status':'CAPPED_NO_CONCLUSION'}
        result['results'].append(row);output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(row),flush=True)
    result['status']='COMPLETE_CUES';result['assurance']='Exploratory only; no target geometry, rank or knot equality certificate.'
    output.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':
    if len(sys.argv)>1:worker(int(sys.argv[2]))
    else:main()
