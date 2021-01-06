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
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(18, 6))

    molecs = []
    t_scores = []
    d_scores = []
    m_scores = []

    basenames = ["2CHT", "2RH1", "3ZOS", "BCD"]
    labels = ["CM", r"$\beta_2$AR", "DDR1", "BCD"]
    for index, basename in enumerate(basenames):
        frame = pd.read_csv(f"data/ZINC_{basename}_MONE.csv")
        smiles = list(frame['s_smiles'])

        docking_scores = list(frame['r_i_docking_score'])
        molecule_scores = list(frame['r_moleculeone_score'])
        mols = [Chem.MolFromSmiles(s) for s in smiles]
        tanimoto_scores = [rediscovery(m, mols) for m in mols]

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

    show_mols(top_molecules, "figure_si_3.svg", legends=legends, min_font_size=12)
    cairosvg.svg2png(url="figure_si_3.svg", write_to="figure_si_3.png")