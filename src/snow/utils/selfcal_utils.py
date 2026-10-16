import os
from contextlib import contextmanager

from casatools import ms, table


@contextmanager
def open_table(table_name: str):
    """
    Context manager that opens a CASA table read-only and always closes it, also on errors, so
    no table lock is left behind

    Parameters
    ----------
    table_name :
        Path to the CASA table, e.g. a measurement set or one of its subtables

    Yields
    ------
    The opened casatools table tool
    """
    table_tool = table()
    table_tool.open(tablename=table_name)
    try:
        yield table_tool
    finally:
        table_tool.close()
        table_tool.done()


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
            with open_table(ms_name) as ms_table:
                return column_name in ms_table.colnames()
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
    with open_table(ms_table) as table_tool:
        return table_tool.nrows()


def calculate_number_antennas(ms_name: str = "") -> int:
    """
    Function that counts the antennas of a measurement set that are not flagged in its ANTENNA
    table (FLAG_ROW). Antennas whose data are all flagged in the main table are still counted.

    Parameters
    ----------
    ms_name :
        Absolute file name to the measurement set file

    Returns
    -------
    nrows :
        Number of antennas not flagged in the ANTENNA table
    """
    if ms_name != "":
        if os.path.exists(ms_name):
            with open_table(ms_name + "/ANTENNA") as antenna_table:
                unflagged = antenna_table.query(query="!FLAG_ROW")
                try:
                    return unflagged.nrows()
                finally:
                    unflagged.close()
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
