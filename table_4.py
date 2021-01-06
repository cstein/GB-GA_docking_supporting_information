import json

import numpy as np
import pandas as pd


def fractions_from_reactions(reactions):
    num_reactions = len(reactions)
    if num_reactions == 0:
        return 0
    num_possible_reactions = len([1 for reaction in reactions if reaction['result'] < 10])
    return num_possible_reactions / num_reactions


def from_dict(d):
    score_fractions_table = {}
    for key in mone_res:
        score_fractions_table[key] = {}
        for method in mone_res[key]:
            score_fractions_table[key][method] = fractions_from_reactions(mone_res[key][method])
    return score_fractions_table


def score_from_filename(filename):
    """ hardcoded """
    data = pd.read_csv(filename)
    return np.size(np.array(np.where(data['r_moleculeone_score'] < 10))) / 100


if __name__ == '__main__':
    with open("data/molecule_one.json", 'r') as dat_file:
        # GA results
        mone_res = json.load(dat_file)
        score_fractions_table = from_dict(mone_res)
        frame = pd.DataFrame(score_fractions_table).transpose()

        # load ZINC results
        zinc_scores = [score_from_filename(f) for f in
                       ["data/ZINC_2CHT_MONE.csv",
                        "data/ZINC_2RH1_MONE.csv",
                        "data/ZINC_3ZOS_MONE.csv",
                        "data/ZINC_BCD_MONE.csv"]
                       ]

        frame2 = pd.DataFrame(frame[["htvs", "htvs_filter", "htvs_sa_norm", "htvs_filter_sa_norm"]])
        frame3 = frame2.drop(index="6LU7")

        frame3["ZINC"] = zinc_scores

        frame4 = frame3.rename(
            columns={"htvs": "GA", "htvs_filter": "+Filter", "htvs_sa_norm": "+SA", "htvs_filter_sa_norm": "+Filter+SA"}
        )
        print(frame4)
