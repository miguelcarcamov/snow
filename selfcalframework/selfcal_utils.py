import numpy as np
from __casac__.table import table

tb = table()
def getTableRows(mstable=""):
    tb.open(tablename=mstable)
    rows = len(tb.rownumbers())
    tb.close()
    return rows
