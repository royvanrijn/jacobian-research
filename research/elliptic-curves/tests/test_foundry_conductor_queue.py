"""The successor queue must not start while its predecessor is incomplete."""
import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cas'))
from run_foundry_conductor_queue_v1 import ready_to_start


class QueueGate(unittest.TestCase):
    def test_requires_complete_sealed_and_idle(self):
        self.assertTrue(ready_to_start(dict(status='COMPLETE',active=[]),99,99))
        self.assertFalse(ready_to_start(dict(status='RUNNING',active=[]),99,99))
        self.assertFalse(ready_to_start(dict(status='COMPLETE',active=['last']),99,99))
        self.assertFalse(ready_to_start(dict(status='COMPLETE',active=[]),98,99))
        self.assertFalse(ready_to_start(dict(status='STOPPED_REVIEW_REQUIRED',active=[]),99,99))


if __name__=='__main__':unittest.main()
