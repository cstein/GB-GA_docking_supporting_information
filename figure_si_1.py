import math
import cairosvg

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from rdkit.Chem.Draw import rdMolDraw2D

from util import docking_score_frame, smiles_to_molecule, show_mols


def get_ECFP4(mol):
    return AllChem.GetMorganFingerprint(mol, 2)


import pandas as pd

if __name__ == '__main__':
    from string import ascii_lowercase
    df_exp = pd.read_excel("data/mpro_xchem_screen_ver-2020-06-12.xlsx")
    df_nc = df_exp.loc[df_exp['Site'] == 'A - active']

    df = pd.read_csv('data/covid_moonshot_submission.csv')
    df_ours = df.loc[
        df['Submission Creator'] == 'Casper Steinmann, Department of Chemistry and Bioscience, Aalborg University'
        ]
    scores = ['Round 2 (-8.0)', 'Round 2 (-7.8)', 'Round 2 (-7.8)', 'Round 2 (-7.6)', 'Round 2 (-7.6)',
              ' Round 2 (-7.6)',
              'Round 3 (-7.7)', 'Round 3 (-7.6)', 'Round 3 (-7.3)', 'Round 3 (-7.1)', '', '', '', '']
    scores = ['{0}. {1}'.format(x, y) for x, y in zip(ascii_lowercase, scores)]
    our_mols = [Chem.MolFromSmiles(smiles) for smiles in df_ours['SMILES']]

    fragments = [Chem.MolFromSmiles(x) for x in df_nc['Compound SMILES']]
    fps = [get_ECFP4(x) for x in fragments]

    mol_pairs = []
    legends = []
    show_mols(fragments, "figure_si_1.svg", legends=ascii_lowercase[:len(fragments)])
    cairosvg.svg2png(url="figure_si_1.svg", write_to="figure_si_1.png")
