import numpy as np

from fingerprints.apps import FingerprintsConfig


def prediction(one_smile_string: str) -> float:
    model = FingerprintsConfig.ml_model
    return model.predict(one_smile_string)


def generating(n: int, energy: float) -> np.array:
    model = FingerprintsConfig.ml_model
    return model.generate(n, energy)


def predict_array_energy(smiles_string: str):
    smiles = set(smiles_string.splitlines())
    output_energy = {}
    for one_smile_string in smiles:
        try:
            output_energy[one_smile_string] = prediction(one_smile_string)
        except RuntimeError as e:
            print(e)
            output_energy[one_smile_string] = "SMILES was written incorrectly! Check the input."
    print(output_energy)
    return output_energy
