import numpy as np
import pandas as pd
from util import docking_score_frame, load_zinc_scores


def compute_num_mols_better(basename, ext, reference_value, mode = "htvs"):
    frame = docking_score_frame(basename, ext, mode)
    values = frame[f"{mode}_ds"]
    return np.size(np.where(values < reference_value))


if __name__ == '__main__':
    # we need the SI table for known binding energies
    KEYS = ["", "_filter", "_sa_norm", "_filter_sa_norm"]

    si_table = pd.read_csv('data/known_binders.csv')
    id_bnd_2cht = 0  # transition state analog binder
    id_bnd_2rh1 = 4  # most potent binder
    id_bnd_3zos = 9  # ponatinib has one charge
    id_bnd_bcd = 12  # most potent binder

    erg_htvs_2cht = si_table["htvs_ds"][id_bnd_2cht]
    erg_sp_2cht = si_table["sp_ds"][id_bnd_2cht]
    erg_htvs_2rh1 = si_table["htvs_ds"][id_bnd_2rh1]
    erg_htvs_3zos = si_table["htvs_ds"][id_bnd_3zos]
    erg_htvs_bcd = si_table["htvs_ds"][id_bnd_bcd]

    # ZINC Data
    zinc_2cht = load_zinc_scores("data/ZINC_2CHT.csv")
    zinc_2cht_data = np.array(zinc_2cht['data'])

    zinc_2rh1 = load_zinc_scores("data/ZINC_2RH1.csv")
    zinc_2rh1_data = np.array(zinc_2rh1['data'])

    zinc_3zos = load_zinc_scores("data/ZINC_3ZOS.csv")
    zinc_3zos_data = np.array(zinc_3zos['data'])

    zinc_bcd = load_zinc_scores("data/ZINC_BCD.csv")
    zinc_bcd_data = np.array(zinc_bcd['data'])

    # Data processing and columns
    cnt_2cht = [compute_num_mols_better("2CHT", s, erg_htvs_2cht) for s in KEYS]
    cnt_2cht += [compute_num_mols_better("2CHT", "_filter_sa_norm", -9),
                 compute_num_mols_better("2CHT", "_filter_sa_norm", -10)]
    cnt_2cht += [np.size(np.where(zinc_2cht_data < erg_htvs_2cht)),
                 np.size(np.where(zinc_2cht_data < -9)),
                 np.size(np.where(zinc_2cht_data < -10))]

    cnt_2cht_sp = [compute_num_mols_better("2CHT", s, erg_sp_2cht, mode="sp") for s in [""]] + [""]*8

    cnt_2rh1 = [compute_num_mols_better("2RH1", s, erg_htvs_2rh1) for s in KEYS]
    cnt_2rh1 += [compute_num_mols_better("2RH1", "_filter_sa_norm", -9),
                 compute_num_mols_better("2RH1", "_filter_sa_norm", -10)]
    cnt_2rh1 += [np.size(np.where(zinc_2rh1_data < erg_htvs_2rh1)),
                 np.size(np.where(zinc_2rh1_data < -9)),
                 np.size(np.where(zinc_2rh1_data < -10))]

    cnt_3zos = [compute_num_mols_better("3ZOS", s, erg_htvs_3zos) for s in KEYS]
    cnt_3zos += [compute_num_mols_better("3ZOS", "_filter_sa_norm", -9),
                 compute_num_mols_better("3ZOS", "_filter_sa_norm", -10)]
    cnt_3zos += [np.size(np.where(zinc_3zos_data < erg_htvs_3zos)),
                 np.size(np.where(zinc_3zos_data < -9)),
                 np.size(np.where(zinc_3zos_data < -10))]

    cnt_bcd = [compute_num_mols_better("BCD", s, erg_htvs_bcd) for s in KEYS]
    cnt_bcd += [compute_num_mols_better("BCD", "_filter_sa_norm", -9),
                compute_num_mols_better("BCD", "_filter_sa_norm", -10)]
    cnt_bcd += [np.size(np.where(zinc_bcd_data < erg_htvs_bcd)),
                 np.size(np.where(zinc_bcd_data < -9)),
                 np.size(np.where(zinc_bcd_data < -10))]

    frame = pd.DataFrame(data=(cnt_2cht, cnt_2cht_sp, cnt_2rh1, cnt_3zos, cnt_bcd),
                 index=["2CHT", "2CHT*", "2RH1", "3ZOS", "BCD"],
                 columns=("HTVS", "+Filter", "+Sa'", "+Filter+Sa'", "< -9.0", "< -10.0", "ZINC", "< -9.0", "< -10.0"))

    print(frame)
