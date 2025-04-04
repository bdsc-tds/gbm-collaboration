"""Test module for the MOSAIC tabular data loader."""

import os

import pytest

from gbmcollaboration import PATIENTS_IDH1_MUT, CHUVDataCenter, UKERDataCenter


IN_GITHUB_ACTIONS = os.getenv("GITHUB_ACTIONS") == "true"


@pytest.mark.skipif(
    IN_GITHUB_ACTIONS,
    reason="Skipping test that requires access to CHUV data files on abstra.",
)
@pytest.mark.parametrize("center", [CHUVDataCenter, UKERDataCenter])
@pytest.mark.parametrize("exclude_idh1_mut", [True, False])
def test_loader_tabular(center, exclude_idh1_mut):
    """Test for CHUV/UKER tabular data loader and cleaner.

    Test that the tabular data loader and cleaner correctly work on the both
    Mosaic cohorts (CHUV, UKER) for the available sources.

    Parameters
    ----------
    center : class
        Data center class to load the data from.
    """
    data_dict = center.load_tabular(exclude_idh1_mut=exclude_idh1_mut)

    for source_name in (
        "clinical",
        "bulk_rna",
        "wes",
    ):
        source_dict = data_dict[source_name]
        assert len(source_dict) > 0, "No dataframes loaded for source {source_name}"
        for file_name, source_df in source_dict.items():
            print("file_name", file_name)
            print("source_df", source_df)
            assert source_df is not None, (
                f"Dataframe for source {source_name}, {file_name} is None"
            )
            assert source_df.index.names in (["Subject Id"], ["EnsemblID"]), (
                f"Index for source {source_name} should be 'Subject Id',"
                f" or 'EnsemblID', got {source_df.index.names}"
            )
            if exclude_idh1_mut:
                # assert that none of the IDH1 mutation samples are found in the columns
                assert not any(idh in source_df.columns for idh in PATIENTS_IDH1_MUT), (
                    f"IDH1 mutation samples found in the {file_name} columns."
                )
                # assert that none of the IDH1 mutation samples are found in the index
                assert not any(idh in source_df.index for idh in PATIENTS_IDH1_MUT), (
                    f"IDH1 mutation samples found in the {file_name} index."
                )
