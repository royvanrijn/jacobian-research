#!/usr/bin/env python3
"""Detached autorun for the corrected v2 11952 transfer wrapper."""
from importlib.machinery import SourceFileLoader
from pathlib import Path

SELF = Path(__file__).resolve()
CAS = SELF.parent
base = SourceFileLoader('v3_transfer_autorun_v1_preserved', str(CAS/'run_v3_transfer_autorun.py')).load_module()
base.TRANSFER = base.LOCAL/'v3-transfer-11952-v2'
base.AUTO = base.LOCAL/'v3-transfer-autorun-v2'
base.STATE = base.AUTO/'state.json'
base.LOG = base.AUTO/'autorun.log'
base.CAMPAIGN = CAS/'v3_transfer_campaign_v2.sage'
base.__file__ = str(SELF)  # detached worker re-enters this corrected wrapper

if __name__ == '__main__':
    base.main()
