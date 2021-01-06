import cairosvg

from util import docking_score_frame, smiles_to_molecule, show_mols

def format_legend(letter, method, score):
    return f"{letter}. {method} ({score:.1f})"

if __name__ == '__main__':
    KEYS = ["", "_filter", "_sa_norm", "_filter_sa_norm"]
    LABELS = ["HTVS", "+Filter", "+SA", "+Filter+SA"]
    SYSTEMS = ["2CHT", "2RH1", "3ZOS", "BCD"]
    LETTERS = " abcdefghijklmnop"
    i = 0

    molecules = []
    legends = []
    for system in SYSTEMS:
        for key, label in zip(KEYS, LABELS):
            i += 1
            frame = docking_score_frame(system, key)
            score = frame['htvs_ds'][0]
            smiles = frame['SMILES'][0]
            molecules.append(smiles_to_molecule(smiles))
            legends.append(format_legend(LETTERS[i], label, score))

    show_mols(molecules, "figure_1.svg", legends=legends, mols_per_row=4)
    cairosvg.svg2png(url="figure_1.svg", write_to="figure_1.png")
    