import os

from casatools import ms, table

tb = table()


def is_column_in_ms(ms_name: str = "", column_name: str = "") -> bool:
    """
        Function that returns True if column_name is present in Measurement Set file

        Parameters
        ----------
        ms_name :
            Measurement set name
        column_name:
            Column name

        Returns
        -------
        True if column_name exists False otherwise
    """
    if ms_name != "":
        if os.path.exists(ms_name):
            # Check if data_column is present in measurement set file
            tb.open(tablename=ms_name)
            col_names = tb.colnames()
            if column_name in col_names:
                return True
            else:
                return False
        else:
            raise FileNotFoundError("The Measurement Set File does not exist")
    else:
        raise ValueError("Measurement Set File cannot be empty")


def get_table_rows(ms_table: str = "") -> int:
    """
    Function that returns the number of rows of a measurement set table

    Parameters
    ----------
    ms_table :
        Measurement set table name

    Returns
    -------
    rows :
        The number of rows of the measurement set table

    """
    tb.open(tablename=ms_table)
    rows = tb.nrows()
    tb.close()
    return rows


def calculate_number_antennas(ms_name: str = "") -> int:
    """
    Function that calculates the number of non-flagged antennas of a measurement set file

    Parameters
    ----------
    ms_name :
        Absolute file name to the measurement set file

    Returns
    -------
    nrows :
        Number of non-flagged antennas
    """
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


def count_spectral_windows(ms_name: str = "") -> int:
    """
    Function that returns the number of spectral windows of a measurement set

    Parameters
    ----------
    ms_name :
        Measurement set name

    Returns
    -------
    The number of rows of the SPECTRAL_WINDOW table
    """
    return get_table_rows(ms_name + "/SPECTRAL_WINDOW")


def get_spw_ids(ms_name: str = "", spw: str = "") -> list:
    """
    Function that returns the ids of the spectral windows selected by a CASA spw selection

    Parameters
    ----------
    ms_name :
        Measurement set name
    spw :
        CASA spectral window selection, e.g. "0~3" or "1,3:10~20". Empty selects every spectral
        window

    Returns
    -------
    Sorted list of selected spectral window ids
    """
    if spw in ("", "*"):
        return list(range(count_spectral_windows(ms_name)))
    ms_tool = ms()
    try:
        selection = ms_tool.msseltoindex(vis=ms_name, spw=spw)
    finally:
        ms_tool.done()
    return sorted(int(spw_id) for spw_id in selection["spw"])
