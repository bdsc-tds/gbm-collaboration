"""Test definitions for the gbmcollaboration package."""

from gbmcollaboration.definitions import MOSAIC_DATA_DICT_FILENAME


def test_mosaic_data_dict_filename():
    """Test the MOSAIC data dictionary filename."""
    assert MOSAIC_DATA_DICT_FILENAME.exists(), f"{MOSAIC_DATA_DICT_FILENAME} not found."
