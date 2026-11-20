"""Binned clinical table, one row per included tumour block (= patient; 273 patients).

Columns: sample_id, patient_id, centre (CH-1/CR-2/ER-3), os_months, os_event,
survival_eligible, mgmt, p53, til_group.

Dates, free text, sex and age (exact or binned) are not released, for privacy reasons:
together with centre and survival they would raise re-identification risk. OS is set to NA for the 13 Charité
patients whose follow-up was pending (SURVIVAL_EXCLUSION), matching the paper's 260
survival-evaluable patients.
"""

import numpy as np
import pandas as pd

import config as C
from common import block_of, centre_code, definitions_list, load_clinical, out, patient_of, til_labels

AVG_DAYS_PER_MONTH = 365.25 / 12


def main():
    clin = load_clinical()
    clin = clin.drop(index=[b for b in C.DUPLICATE_BLOCK_DROP if b in clin.index])
    assert patient_of(clin.index).is_unique
    surv_excl = set(definitions_list("SURVIVAL_EXCLUSION"))

    df = pd.DataFrame(index=block_of(clin.index))
    df.index.name = "sample_id"
    df["patient_id"] = patient_of(df.index)
    df["centre"] = centre_code(df.index)
    df["survival_eligible"] = ~df.index.isin(surv_excl)
    os_m = clin["os_years"].to_numpy() * 12
    df["os_months"] = np.where(df.survival_eligible, np.round(os_m, 1), np.nan)
    df["os_event"] = np.where(df.survival_eligible, clin["os_event"].to_numpy(), np.nan)
    df["mgmt"] = clin["MGMT promoter methylation"].to_numpy()
    df["p53"] = clin["p53"].to_numpy()

    til = til_labels()
    # TIL-high = top quartile of TIL+ fraction over the 228 sections; TIL-low = the other 171,
    # including the 23 without any TIL+ spot. NA = no Visium.
    df["til_group"] = til["til_group"].reindex(df.index)

    df.to_csv(out("clinical", "clinical_binned.csv"))
    print(df.shape)
    print(df.drop(columns=["patient_id", "os_months"]).apply(lambda s: s.value_counts(dropna=False).to_dict()))


if __name__ == "__main__":
    main()
