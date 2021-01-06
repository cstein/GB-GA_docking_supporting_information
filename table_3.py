import numpy as np
import pandas as pd
from util import docking_score_frame, load_zinc_scores


def population_means(basename, ext, mode = "htvs"):
    frame = docking_score_frame(basename, ext, mode)
    scores = np.array(frame[f"{mode}_ds"])
    indices = np.where(scores < 0)
    return np.mean(scores[indices]), np.std(scores[indices])


if __name__ == '__main__':
    # we need the SI table for known binding energies
    KEYS = ["", "_filter", "_sa_norm", "_filter_sa_norm"]

    ms_2cht = [population_means("2CHT", s) for s in KEYS]
    ms_2cht_sp = [population_means("2CHT", s, mode="sp") for s in [""]]
    ms_2rh1 = [population_means("2RH1", s) for s in KEYS]
    ms_3zos = [population_means("3ZOS", s) for s in KEYS]
    ms_bcd = [population_means("BCD", s) for s in KEYS]

    # ZINC Data
    zinc_2cht = load_zinc_scores("data/ZINC_2CHT.csv")
    zinc_2cht_data = np.array(sorted(zinc_2cht['data']))[:8000]

    zinc_2rh1 = load_zinc_scores("data/ZINC_2RH1.csv")
    zinc_2rh1_data = np.array(sorted(zinc_2rh1['data']))[:8000]

    zinc_3zos = load_zinc_scores("data/ZINC_3ZOS.csv")
    zinc_3zos_data = np.array(sorted(zinc_3zos['data']))[:8000]

    zinc_bcd = load_zinc_scores("data/ZINC_BCD.csv")
    zinc_bcd_data = np.array(sorted(zinc_bcd['data']))[:8000]

    ms_2cht.append((zinc_2cht_data.mean(), zinc_2cht_data.std()))
    ms_2rh1.append((zinc_2rh1_data.mean(), zinc_2rh1_data.std()))
    ms_3zos.append((zinc_3zos_data.mean(), zinc_3zos_data.std()))
    ms_bcd.append((zinc_bcd_data.mean(), zinc_bcd_data.std()))

    table = pd.DataFrame(data=(ms_2cht, ms_2cht_sp, ms_2rh1, ms_3zos, ms_bcd),
                             index=["2CHT", "2CHT*", "2RH1", "3ZOS", "BCD"],
                             columns=("HTVS", "+Filter", "+SA", "+Filter+SA", "ZINC"))

    print(table.to_string(formatters={"HTVS": "{0[0]:.1f} +/- {0[1]:.1f}".format,
                                          "+SA": "{0[0]:.1f} +/- {0[1]:.1f}".format,
                                          "+Filter": "{0[0]:.1f} +/- {0[1]:.1f}".format,
                                          "+Filter+SA": "{0[0]:.1f} +/- {0[1]:.1f}".format,
                                          "ZINC": "{0[0]:.1f} +/- {0[1]:.1f}".format
                                          }))
