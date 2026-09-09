"""Editorial selection must retain strongest proofs and distinguish Q-twists."""
import hashlib
import json
from pathlib import Path
import sys
import unittest
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from refresh_foundry_curve_ledger import select,validate_record


def row(rank,model,parameter):
    return dict(rank_lower_bound=rank,packet=dict(curve=model),family='test',parameter=parameter)


class Selection(unittest.TestCase):
    @patch('refresh_foundry_curve_ledger.sha',return_value='testhash')
    @patch('refresh_foundry_curve_ledger.validate_record')
    def test_strongest_and_twists(self,*mocks):
        model=['0','0','0','-1','1'];twist=['0','0','0','-4','8']
        candidates=[('a',row(22,model,'1')),('b',row(23,model,'1')),('c',row(22,twist,'2'))]
        selected=select(candidates,[])
        self.assertEqual(len(selected),2)
        self.assertEqual(sorted(r['result']['rank_lower_bound'] for r in selected.values()),[22,23])

    @patch('refresh_foundry_curve_ledger.sha',return_value='testhash')
    @patch('refresh_foundry_curve_ledger.validate_record')
    def test_existing_only_strengthens(self,*mocks):
        model=['0','0','0','-1','1']
        baseline=[dict(id='existing',ainvs=model,rank_lower_bound=23)]
        self.assertEqual(select([('a',row(22,model,'1'))],baseline),{})
        selected=select([('b',row(25,model,'1'))],baseline)
        self.assertEqual(list(selected),['existing'])
        self.assertEqual(selected['existing']['previous_rank_lower_bound'],23)

    def test_packet_binding(self):
        packet=dict(points=[[1,2]],rank_lower_bound=1)
        digest=hashlib.sha256((json.dumps(packet,indent=2,sort_keys=True)+'\n').encode()).hexdigest()
        result=dict(status='PASS_CERTIFIED_SEARCH',packet=packet,rank_lower_bound=1,packet_sha256=digest,
                    certificate_replay=dict(status='PASS_TWO_FINITE_IMPLEMENTATIONS',packet_sha256=digest,rank_lower_bound=1))
        validate_record(result)
        result['packet']['points'][0][0]=3
        with self.assertRaisesRegex(ValueError,'hash differs'):validate_record(result)


if __name__=='__main__':unittest.main()
