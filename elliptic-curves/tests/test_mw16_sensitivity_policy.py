import copy
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from check_mw16_sensitivity_policy import check_policy


class PolicyBindingTests(unittest.TestCase):
    def fixture(self):
        chart={'mask':1,'specification':'metric:16','height_bound':100000,
               'centre_construction':{'centre':'generic'}}
        return {'declared_budget':{'mode':'initial','parent_ids':['blind'],
            'centres':['generic'],'specifications':['metric:16'],'heights':[100000]},
            'results':[{'parent_id':'blind','deepest_masks':[1],
                'settings':[{'centre':'generic','specification':'metric:16','height':100000,'charts':[chart]}]}]}

    def test_bound_and_coordinate_mislabelling_are_rejected(self):
        source=self.fixture();check_policy(source)
        for key,value in (('height_bound',10000),('specification','gauss')):
            changed=copy.deepcopy(source);changed['results'][0]['settings'][0]['charts'][0][key]=value
            with self.assertRaises(ArithmeticError):check_policy(changed)

    def test_missing_and_duplicate_classes_are_rejected(self):
        for charts in ([],[1,1]):
            source=self.fixture();setting=source['results'][0]['settings'][0]
            setting['charts']=[setting['charts'][0] for _ in charts]
            with self.assertRaises(ArithmeticError):check_policy(source)


if __name__=='__main__':unittest.main()
