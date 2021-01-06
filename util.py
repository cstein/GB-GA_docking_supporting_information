import csv
import zipfile

import pandas as pd
import numpy as np

from rdkit import Chem


def smiles_to_molecule(s):
    """ Converts a smiles string to RDKit Molecule """
    try:
        mol = Chem.MolFromSmiles(s)
    except TypeError:
        mol = None
    finally:
        return mol


def load_frame_from_zip(basename, ext, mode):
    """ Loads a pandas.DataFrame from .csv data in a zip-file """
    with zipfile.ZipFile(f"data/{basename}.zip") as z_file:
        return pd.read_csv(z_file.open(f"{basename}_{mode}{ext}.csv"))


def docking_score_frame(basename, ext, mode="htvs"):
    """ SA Score corrected docking score frame """
    raw_docking_score_frame = load_frame_from_zip(basename, ext, mode)
    frame = raw_docking_score_frame[['i', 'SMILES']].copy()

    docking_scores = np.array(raw_docking_score_frame[f"{mode}_ds"])

    # Load the SA data
    sa_frame = pd.read_pickle(f"data/{basename}.pkl")
    sa_score_coeff = np.array(sa_frame[f"{ext}"])

    # compute the updated docking score in the new frame
    frame[f"{mode}_ds"] = docking_scores  # * sa_score_coeff[:len(docking_scores)]

    # drop duplicates and sort by docking score values
    frame.drop_duplicates(subset="SMILES", keep=False, inplace=True)
    frame.sort_values(by=f"{mode}_ds", ascending=True, inplace=True, ignore_index=True)
    return frame


def top_n_score_smiles_as_frame(basename, ext, n=100, mode='htvs'):
    frame = docking_score_frame(basename, ext, mode)
    return frame[:n]


def load_csv_iter(filename):
    with open(filename, 'r') as csv_file:
        csv_reader = csv.reader(csv_file)
        for i, row in enumerate(csv_reader):
            yield i, row


def load_zinc_scores(filename):
    data = []
    imin = 0
    vmin = 0.0
    n_failed = 0
    n_noposes = 0
    for i, row in load_csv_iter(filename):
        status = row[3]
        if "Done" in status:
            value = float(row[5])
            if value < 0.0:
                if value < vmin:
                    vmin = value
                    imin = i
                data.append(value)
        if "GRID-ENERGY MIN FAILED" in status:
            n_failed += 1
        if "NO GOOD POSES FOUND" in status:
            n_noposes += 1
    result = {'data': data, 'idx_min': imin, 'val_min':  vmin, 'failed': n_failed, 'noposes': n_noposes}
    return result

