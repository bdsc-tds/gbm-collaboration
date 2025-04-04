"""Test module for the MOSAIC spatial data loader."""

import os

import pytest

from gbmcollaboration import CHUVDataCenter, UKERDataCenter


IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


@pytest.mark.skipif(
    IN_GITHUB_ACTIONS,
    reason="Skipping test that requires access to CHUV data files on abstra.",
)
@pytest.mark.parametrize(
    "normalization", ["log_norm_count", "raw_count", "norm_count", "SCTransform"]
)
@pytest.mark.parametrize("radius", [1, 2])
@pytest.mark.parametrize("center", [CHUVDataCenter, UKERDataCenter])
@pytest.mark.parametrize("exclude_low_quality_visium", [True, False])
@pytest.mark.parametrize("exclude_idh1_mut", [True, False])
def test_load_count(
    normalization, radius, center, exclude_low_quality_visium, exclude_idh1_mut
):
    """Test for CHUV and UKER spatial data loader.

    Test that the spatial data loader correctly works on the CHUV and UKER cohort.

    Parameters
    ----------
    normalization : str
        Normalization method.
    radius : int
        Radius for the spatial data to define neighboring spots for
        downstream analyses.
    center : class
        Data center class to load the data from.
    """
    spatial_object = center.load_spatial(
        normalization=normalization,
        radius=radius,
        exclude_low_quality_visium=exclude_low_quality_visium,
        exclude_idh1_mut=exclude_idh1_mut,
        deconvolution_slot_name="deconv_GMBatlas_2steps",
    )
    dict_counts = spatial_object.load_counts_cohort()
    if center == CHUVDataCenter:
        if exclude_low_quality_visium:  # no need to test IDH1 as there is none
            assert len(dict_counts.keys()) == 61
        else:
            assert len(dict_counts.keys()) == 63  # as of 2024-15-11
    elif center == UKERDataCenter:
        if exclude_low_quality_visium and not exclude_idh1_mut:
            assert len(dict_counts.keys()) == 43
        elif exclude_low_quality_visium and exclude_idh1_mut:
            assert len(dict_counts.keys()) == 42
        elif not exclude_low_quality_visium and not exclude_idh1_mut:
            assert len(dict_counts.keys()) == 68  # as of 2024-15-11
    # read the df and check that it is not empty
    assert not spatial_object.genecounts_df.empty
    assert not spatial_object.celltype_df.empty
