import math
import cairosvg

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from rdkit.Chem.Draw import rdMolDraw2D

import pandas as pd

from util import docking_score_frame, smiles_to_molecule, show_mols


def format_legend(letter, method, score):
    return f"{letter}. {method} ({score:.1f})"


def get_ECFP4(mol):
    return AllChem.GetMorganFingerprint(mol, 2)


if __name__ == '__main__':
    df = pd.read_csv('data/covid_moonshot_submission.csv')
    df_ours = df.loc[
        df['Submission Creator'] == 'Casper Steinmann, Department of Chemistry and Bioscience, Aalborg University'
        ]

    df_ordered = df.loc[df['Ordered'] == True]
    mols = [Chem.MolFromSmiles(smiles) for smiles in df_ordered['SMILES']]
    fps = [get_ECFP4(mol) for mol in mols]

    from string import ascii_lowercase

    all_sims = []
    mol_pairs = []
    legends = []
    scores = ['Round 2 (-8.0)', 'Round 2 (-7.8)', 'Round 2 (-7.8)', 'Round 2 (-7.6)', 'Round 2 (-7.6)',
              ' Round 2 (-7.6)',
              'Round 3 (-7.7)', 'Round 3 (-7.6)', 'Round 3 (-7.3)', 'Round 3 (-7.1)', '', '', '', '']
    scores = ['{0}. {1}'.format(x, y) for x, y in zip(ascii_lowercase, scores)]

    our_mols = [Chem.MolFromSmiles(smiles) for smiles in df_ours['SMILES']]
    for mol, score in zip(our_mols, scores):
        fp_target = get_ECFP4(mol)
        sim = [(mol, TanimotoSimilarity(fp, fp_target)) for mol, fp in zip(mols, fps)]
        sim.sort(key=lambda x: x[1], reverse=True)
        all_sims.append(sim)
        mol_pairs.append(mol)
        mol_pairs.append(sim[0][0])
        legends.append(score)
        legends.append(str(round(sim[0][1], 2)))

    show_mols(mol_pairs[:20], "figure_3.svg", legends=legends[:20], mols_per_row=4)
    cairosvg.svg2png(url="figure_3.svg", write_to="figure_3.png")
