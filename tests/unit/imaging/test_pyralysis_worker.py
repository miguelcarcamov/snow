"""The Pyralysis worker script (#24), without Pyralysis: the parts that only inspect its data."""
from types import SimpleNamespace

import pytest

from snow.imaging.backends import _pyralysis_worker as worker


def dataset(*columns_per_partition):
    return SimpleNamespace(
        ms_list=[
            SimpleNamespace(visibilities=SimpleNamespace(dataset=dict.fromkeys(columns)))
            for columns in columns_per_partition
        ]
    )


def test_loaded_corrected_data_passes():
    worker.check_corrected(dataset(("DATA", "CORRECTED_DATA"), ("DATA", "CORRECTED_DATA")))


def test_a_partition_without_corrected_data_raises_instead_of_fitting_data():
    with pytest.raises(ValueError, match=r"partitions \[1\].*DATA"):
        worker.check_corrected(dataset(("DATA", "CORRECTED_DATA"), ("DATA", )))
