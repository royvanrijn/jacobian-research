#!/usr/bin/env python3
"""Detached autorun for transfer wrapper v3."""
from importlib.machinery import SourceFileLoader
from pathlib import Path

SELF = Path(__file__).resolve()
CAS = SELF.parent
base = SourceFileLoader('v3_transfer_autorun_v1_preserved', str(CAS/'run_v3_transfer_autorun.py')).load_module()
base.TRANSFER = base.LOCAL/'v3-transfer-11952-v3'
base.AUTO = base.LOCAL/'v3-transfer-autorun-v3'
base.STATE = base.AUTO/'state.json'
base.LOG = base.AUTO/'autorun.log'
base.CAMPAIGN = CAS/'v3_transfer_campaign_v3.sage'
base.__file__ = str(SELF)

if __name__ == '__main__':
    base.main()
