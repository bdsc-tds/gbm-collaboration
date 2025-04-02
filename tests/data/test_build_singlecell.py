"""Test module for the MOSAIC spatial data loader."""

import os

import pytest

from gbmcollaboration import CHUVDataCenter, UKERDataCenter


IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


@pytest.mark.skipif(
    IN_GITHUB_ACTIONS,
    reason="Skipping test that requires access to CHUV data files on abstra.",
)
@pytest.mark.parametrize("center", [CHUVDataCenter, UKERDataCenter])
def test_load_singlecell_rna(
    center,
):
    """Test for CHUV and UKER single cell data loader.

    Test that the single cell data loader correctly works on the CHUV and UKER cohort.

    Parameters
    ----------
    center : class
        Data center class to load the data from.
    """
    singlecell_object = center.load_singlecell()
    if center == CHUVDataCenter:
        assert len(singlecell_object.obs["orig.ident"].unique()) == 63
    else:
        assert len(singlecell_object.obs["orig.ident"].unique()) == 64
    # in both cases check that the adata.X is not empty
    assert singlecell_object.X.shape[0] > 0
