import pytest

from snow.selfcalibration.caltable import CalTable, build_gaintable_arguments, derive_spwmap


class TestDeriveSpwmap:

    @pytest.mark.parametrize("combine", ["", "scan", "field,scan"])
    def test_without_spw_combination_no_map_is_needed(self, combine):
        assert derive_spwmap(combine, [0, 1, 2, 3], 4) == []

    @pytest.mark.parametrize("combine", ["spw", "scan,spw", "scan, spw"])
    def test_combined_solutions_map_to_lowest_selected_spw(self, combine):
        assert derive_spwmap(combine, [3, 1, 2], 5) == [1, 1, 1, 1, 1]

    def test_combining_without_selected_spws_raises(self):
        with pytest.raises(ValueError, match="No spectral windows"):
            derive_spwmap("spw", [], 4)


class TestBuildGaintableArguments:

    def test_one_entry_per_table_in_application_order(self):
        tables = [CalTable("pcal", "p", (0, 0, 0, 0)), CalTable("apcal", "ap")]

        assert build_gaintable_arguments(tables, "linear") == {
            "gaintable": ["pcal", "apcal"],
            "spwmap": [[0, 0, 0, 0], []],
            "interp": ["linear", "linear"],
        }

    def test_no_tables_gives_empty_lists_not_empty_names(self):
        assert build_gaintable_arguments([], "linear") == {
            "gaintable": [],
            "spwmap": [],
            "interp": [],
        }


def test_calmode_properties():
    assert CalTable("t", "p").solves_phase and not CalTable("t", "p").solves_amplitude
    assert CalTable("t", "a").solves_amplitude and not CalTable("t", "a").solves_phase
    assert CalTable("t", "ap").solves_phase and CalTable("t", "ap").solves_amplitude
