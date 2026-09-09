"""Presentation-only minima for the rounded numeric curve-table columns."""
from decimal import Decimal
import re

LEGEND = ('Bold marks the lowest displayed value in each metric within each '
          'certified-rank group in this table (rounded ties included; missing '
          'values ignored). This is not a conductor-record or exact-rank claim.')


def highlight_rank_minima(section):
    lines=section.splitlines();rows=[];minima={}
    for index,line in enumerate(lines):
        cells=line.split('|')
        if len(cells)!=9 or not re.fullmatch(r'≥\s*\d+',cells[3].strip()):continue
        rank=int(cells[3].strip()[1:].strip());values={}
        for column in range(4,8):
            value=cells[column].strip().removeprefix('**').removesuffix('**')
            cells[column]=' '+value+' '
            if value=='—':continue
            if not re.fullmatch(r'-?\d+(?:\.\d+)?',value):
                raise ValueError('unexpected curve-table metric: '+value)
            values[column]=Decimal(value)
            key=(rank,column)
            minima[key]=min(minima.get(key,values[column]),values[column])
        rows.append((index,cells,rank,values))
    for index,cells,rank,values in rows:
        for column,value in values.items():
            if value==minima[(rank,column)]:cells[column]=' **'+cells[column].strip()+'** '
        lines[index]='|'.join(cells)
    return '\n'.join(lines)
