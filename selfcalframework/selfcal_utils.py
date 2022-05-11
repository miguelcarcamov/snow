import numpy as np
from casatools import table

tb = table()

def getTableRows(mstable=""):
    tb.open(tablename=mstable)
    rows = tb.nrows()
    tb.close()
    return rows
