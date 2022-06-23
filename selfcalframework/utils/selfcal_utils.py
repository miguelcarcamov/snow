import numpy as np
from casatools import table

tb = table()


def getTableRows(mstable=""):
    tb.open(tablename=mstable)
    rows = tb.nrows()
    tb.close()
    return rows


def calculate_number_antennas(ms_name: str = ""):
    if ms_name != "":
        if os.path.exists(ms_name):
            tb.open(tablename=ms_name + "/ANTENNA")
            query_table = tb.taql("select NAME from " + ms_name + "/ANTENNA" + " where !FLAG_ROW")
            nrows = len(query_table.getcol("NAME"))
            tb.close()
            return nrows
        else:
            raise FileNotFoundError("The Measurement Set File does not exist")
    else:
        raise ValueError("Measurement Set File cannot be empty")
