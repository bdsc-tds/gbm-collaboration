"""Top level package for the :mod:`gbmcollaboration`."""

import os
#import git
from gbmcollaboration.data import CHUVDataCenter, UKERDataCenter
from gbmcollaboration.definitions import (
    BAD_QUALITY_CHUV_VISIUM,
    BAD_QUALITY_UKER_VISIUM,
    PATIENTS_IDH1_MUT,
    PROBE_ISSUE_CHUV,
    PROBE_ISSUE_UKER,
    UKER_MAX_ID,
    CHUV_MAX_ID
)

#git_repo = git.Repo(__file__, search_parent_directories=True)
#git_root = os.path.abspath(git_repo.git.rev_parse("--show-toplevel"))
