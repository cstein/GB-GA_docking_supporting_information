import json
import os

import matplotlib.pyplot as plt

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem

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

    targets = ["2CHT", "2RH1", "3ZOS", "BCD"]
    labels = ["CM", r"$\beta_2$AR", "DDR1", "BCD"]
    for index, target in enumerate(targets):
        top_n = top_n_score_smiles_as_frame(target, '_filter_sa_norm')
        top_m = mone_res[target]['htvs_filter_sa_norm']

        # left panel has docking scores matched with molecule.one scores
        docking_scores = []
        molecule_scores = []
        for i in range(len(top_n)):
            smi_n = top_n["SMILES"][i]
            neu_smi_n = neutralize.neutralize_smiles([smi_n])
            mol_n = smiles_to_molecule(neu_smi_n[0])
            sco_n = float(top_n["htvs_ds"][i])
            can_smi_n = Chem.MolToSmiles(mol_n, canonical=True)

            for j in range(100):
                smi_m = top_m[j]['targetSmiles']
                sco_m = float(top_m[j]['result'])
                mol_m = smiles_to_molecule(smi_m)
                can_smi_m = Chem.MolToSmiles(mol_m, canonical=True)
                if can_smi_m == can_smi_n:
                    docking_scores.append(sco_n)
                    molecule_scores.append(sco_m)
                    break

        ax.scatter(docking_scores, molecule_scores, label=labels[index], alpha=0.9, s=48)
        # left panel done

        # right panel has has docking scores matched with tanimoto similarity to best scoring
        target_top_smiles = neutralize.neutralize_smiles(list(top_n["SMILES"]))
        target_top_scores = list(top_n["htvs_ds"])

        target_top_mols = [smiles_to_molecule(s) for s in target_top_smiles]

        tanimoto_scores = [rediscovery(m, target_top_mols) for m in target_top_mols]
        bx.scatter(target_top_scores, tanimoto_scores)
        # right panel done

    # fix axes. ax first
    ax.legend(loc='upper left', fontsize=14)
    ax.set_ylim((0, 10.2))
    ax.set_xlim((-13, -6))
    ax.set_xlabel("Docking Score", fontsize=14)
    ax.set_ylabel("Molecule.one Score", fontsize=14)
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

    ax.text(0.95, 0.92, "a", transform=ax.transAxes, fontsize=20)

    # now bx
    bx.set_xlim((-13, -6))
    bx.set_xlabel("Docking Score", fontsize=14)
    bx.set_ylabel("Tanimoto Similarity", fontsize=14)
    for tick in bx.xaxis.get_major_ticks():
        tick.label.set_fontsize(14)
    for tick in bx.yaxis.get_major_ticks():
        tick.label.set_fontsize(14)

    bx.text(0.95, 0.92, "b", transform=bx.transAxes, fontsize=20)

    fig.savefig('figure_2.png', dpi=72)
