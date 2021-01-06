import json
import os
import math

import matplotlib.pyplot as plt
import cairosvg
import pandas as pd

from rdkit.Chem.Draw import rdMolDraw2D
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


def show_mols(mols, filename, mols_per_row = 5, size=200, min_font_size=12, legends=[]):
  rows = math.ceil(len(mols)/mols_per_row)
  d2d = rdMolDraw2D.MolDraw2DSVG(mols_per_row*size,rows*size,size,size)
  d2d.drawOptions().minFontSize = min_font_size
  if legends:
    d2d.DrawMolecules(mols, legends=legends)
  else:
    d2d.DrawMolecules(mols)
  d2d.FinishDrawing()

  with open(filename, 'w') as f:
    f.write(d2d.GetDrawingText())


if __name__ == '__main__':
    mone_res = molecule_one_results("data/molecule_one.json")
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(18, 6))

    molecs = []
    t_scores = []
    d_scores = []
    m_scores = []

    targets = ["2CHT", "2RH1", "3ZOS", "BCD"]
    labels = ["CM", r"$\beta_2$AR", "DDR1", "BCD"]
    for index, target in enumerate(targets):
        top_n = top_n_score_smiles_as_frame(target, '_filter_sa_norm')
        top_m = mone_res[target]['htvs_filter_sa_norm']

        # left panel has docking scores matched with molecule.one scores
        docking_scores = []
        molecule_scores = []
        mols = []
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
                    mols.append(Chem.MolFromSmiles(can_smi_m))
                    docking_scores.append(sco_n)
                    molecule_scores.append(sco_m)
                    break

        target_top_smiles = neutralize.neutralize_smiles(list(top_n["SMILES"]))
        target_top_scores = list(top_n["htvs_ds"])

        target_top_mols = [smiles_to_molecule(s) for s in target_top_smiles]

        tanimoto_scores = [rediscovery(m, target_top_mols) for m in target_top_mols]

        molecs.append(mols)
        d_scores.append(docking_scores)
        m_scores.append(molecule_scores)
        t_scores.append(tanimoto_scores)

    # prepare to make figure with RDKit
    top_molecules = []
    legends = []
    num_molecules = 10
    for target, d_score, t_score, m_score in zip(molecs, d_scores, t_scores, m_scores):
        top_molecules += target[:num_molecules]
        legends += ["{0:.1f} ({1:.2f}) [{2:.1f}]".format(x, y, z) for x, y, z in
                    zip(d_score[:num_molecules], t_score[:num_molecules], m_score[:num_molecules])]

    show_mols(top_molecules, "figure_si_2.svg", legends=legends, min_font_size=12)
    cairosvg.svg2png(url="figure_si_2.svg", write_to="figure_si_2.png")