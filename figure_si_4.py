import json
import os

import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

import pandas as pd

import neutralize
from util import top_n_score_smiles_as_frame, smiles_to_molecule


def get_ECFP4(mol):
    return AllChem.GetMorganFingerprint(mol, 2)


def rediscovery(mol,args):
    target = args[0]
    fp_mol = get_ECFP4(mol)
    fp_target = get_ECFP4(target)
    score = DataStructs.TanimotoSimilarity(fp_mol, fp_target)
    return score
    # except:
    #     print('Failed ', Chem.MolToSmiles(mol))
    #    return None


def molecule_one_results(filename):
    """ Computes the molecule.one score fraction

        note: should probably store results as .json to not overload server
    """
    if os.path.exists(filename):
        with open(filename, 'r') as dat_file:
            return json.load(dat_file)
    else:
        raise FileNotFoundError(f"Could not find the file {filename}")


if __name__ == '__main__':
    mone_res = molecule_one_results("data/molecule_one_results.json")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(18, 6))

    basenames = ["2CHT", "2RH1", "3ZOS", "BCD"]
    labels = ["CM", r"$\beta_2$AR", "DDR1", "BCD"]
    for index, basename in enumerate(basenames):
        frame = pd.read_csv(f"data/ZINC_{basename}_MONE.csv")
        docking_scores = list(frame['r_i_docking_score'])
        molecule_scores = list(frame['r_moleculeone_score'])
        ax.scatter(docking_scores, molecule_scores, label=labels[index], alpha=0.9, s=48)

        smiles = list(frame['s_smiles'])
        mols = [Chem.MolFromSmiles(s) for s in smiles]
        tanimoto_scores = [rediscovery(m, mols) for m in mols]

        bx.scatter(docking_scores, tanimoto_scores)

    ax.legend(loc='upper left', fontsize=14)
    ax.set_ylim((0, 10.2))
    ax.set_xlim((-13, -6))
    ax.set_xlabel("Docking Score [kcal/mol]", fontsize=14)
    ax.set_ylabel("Molecule.one Score", fontsize=14)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

    ax.text(0.95, 0.92, "a", transform=ax.transAxes, fontsize=20)

    # now bx
    bx.set_xlim((-13, -6))
    bx.set_xlabel("Docking Score [kcal/mol]", fontsize=14)
    bx.set_ylabel("Tanimoto Similarity", fontsize=14)
    for tick in bx.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in bx.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

    bx.text(0.95, 0.92, "b", transform=bx.transAxes, fontsize=20)

    fig.subplots_adjust(left=0.06, right=0.96)
    fig.savefig('figure_si_4.png', dpi=72)
