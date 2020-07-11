import numpy as np
import tb as tb


def getTableRows(mstable=""):
    tb.open(mstable)
    rows = len(tb.rownumbers())
    tb.close()
