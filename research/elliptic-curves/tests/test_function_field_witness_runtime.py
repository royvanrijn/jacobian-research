"""Retained finite function-field group laws reject algebraic tampering."""
import copy
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from research_runtime.function_field_witness import GroupLawDiscovery,replay_group_law,RationalPair
from sage.all import GF,PolynomialRing,EllipticCurve

class GroupLawWitnessTests(unittest.TestCase):
    def test_group_dag_replays_without_group_addition(self):
        R=PolynomialRing(GF(17),'t');K=R.fraction_field();E=EllipticCurve(K,[-7,10])
        basis=[E(1,2),E(2,2)];graph=GroupLawDiscovery(E,basis)
        words=[(1,1),(2,1),(-3,2),(3,-5),(12,0)]
        indices=[graph.trace(w) for w in words]
        with patch.object(GroupLawDiscovery,'add',side_effect=AssertionError('addition rediscovery forbidden')):
            points,labels=replay_group_law(R,R(-7),R(10),basis,graph.nodes)
        self.assertEqual([labels[i] for i in indices],words)
        resumed=GroupLawDiscovery(E,basis,retained=graph.nodes)
        self.assertEqual(resumed.nodes,graph.nodes)
        bad=copy.deepcopy(graph.nodes)
        j=next(i for i,r in enumerate(bad) if r['kind']=='sum' and r['point'] is not None)
        coefficients=bad[j]['point']['jacobian'][0].split(',');coefficients[0]=str(int(coefficients[0])+1)
        bad[j]['point']['jacobian'][0]=','.join(coefficients)
        with self.assertRaises(ArithmeticError):replay_group_law(R,R(-7),R(10),basis,bad)
        bad=copy.deepcopy(graph.nodes);bad[j]['word'][0]+=1
        with self.assertRaises(ArithmeticError):replay_group_law(R,R(-7),R(10),basis,bad)

    def test_reciprocal_rational_identity(self):
        R=PolynomialRing(GF(101),'t');t=R.gen();K=R.fraction_field()
        for n,d,weight in [(1+t,1+t*t,4),(t**8+3*t,1+t,4),(R(0),t+1,6)]:
            n,d=R(n),R(d);result=RationalPair(n,d).inverted_parameter(weight)
            expected=t**weight*K(n/d)(1/t)
            self.assertEqual(K(result.n/result.d),expected)

class CensusWitnessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        import gzip,json,runpy,zipfile
        cls.root=Path(__file__).resolve().parents[2]
        cls.api=runpy.run_path(str(cls.root/'elkies-k3/scripts/retain_r17_mw18_census_witnesses.sage'))
        with zipfile.ZipFile(cls.api['BUNDLE']) as bundle:
            cls.manifest=json.loads(bundle.read('manifest.json'))
            cls.entry=cls.manifest['entries'][0]
            cls.row=json.loads(gzip.decompress(bundle.read(cls.entry['member'])))

    def check_batch(self,row):
        import contextlib,io
        from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
        entry=self.entry
        with contextlib.redirect_stdout(io.StringIO()),patch.object(EllipticCurvePoint_field,'_add_',side_effect=AssertionError('elliptic addition forbidden during census replay')),patch.object(GroupLawDiscovery,'__init__',side_effect=AssertionError('discovery forbidden during replay')):
            return self.api['replay'](entry['chart'],entry['prime'],row,batch_index=entry['batch_index'])

    def test_real_census_batch_replays_without_elliptic_addition(self):
        result=self.check_batch(self.row)
        self.assertEqual(result['verified_roots'],len(self.row['roots']))
        self.assertGreater(result['verified_nonresidues'],0)

    def test_census_chord_scope_and_missing_root_tampering(self):
        bad=copy.deepcopy(self.row);bad['roots'][0]['chord']['q']='0'
        with self.assertRaises(ArithmeticError):self.check_batch(bad)
        bad=copy.deepcopy(self.row);bad['roots'][0]['obstructed_curves'].append(-1)
        with self.assertRaises(ArithmeticError):self.check_batch(bad)
        bad=copy.deepcopy(self.row);bad['roots'].pop()
        with self.assertRaises(ArithmeticError):self.check_batch(bad)

    def test_aggregate_rejects_missing_batch(self):
        import json,tempfile,zipfile
        manifest=copy.deepcopy(self.manifest);manifest['entries'].pop()
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'incomplete.zip'
            with zipfile.ZipFile(path,'w') as bundle:bundle.writestr('manifest.json',json.dumps(manifest))
            with self.assertRaisesRegex(ArithmeticError,'missing or duplicate'):
                self.api['replay_bundle'](path)

if __name__=='__main__':unittest.main()
