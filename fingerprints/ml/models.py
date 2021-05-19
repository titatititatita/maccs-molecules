# import tensorflow as tf
import numpy as np
import os
# from typing import Tuple
# from rdkit import Chem
# from rdkit.Chem import MACCSkeys

DEFAULT_DECODER_MODEL_PATH = "fingerprints/ml/decoder.h5"
DEFAULT_REGRESSION_MODEL_PATH = "fingerprints/ml/HIV_keras_energy_18_5.h5"
LATENT_LAYER_SIZE = 10


class TF_Models:
    "Container class for tensorflow models"

    def __init__(
        self, decoder_model_path=None, energy_regression_model_path=None
    ):

        self.decoder_model_path = DEFAULT_DECODER_MODEL_PATH
        self.energy_regression_model_path = DEFAULT_REGRESSION_MODEL_PATH

        if decoder_model_path is not None:
            self.decoder_model_path = decoder_model_path

        if energy_regression_model_path is not None:
            self.energy_regression_model_path = energy_regression_model_path

        self.decoder_model, self.energy_regression_model = self._get_models(
            self.decoder_model_path, self.energy_regression_model_path
        )

    @staticmethod
    def _get_models(
        decoder_model_path: str, regression_model_path: str
    ): #-> Tuple[tf.keras.Model, tf.keras.Model]:
        """
        Function for loading a keras model from a .h5 file

        Parameters
        ----------
        decoder_model_path: str
             string path to a file with a saved decoder model
        regression_model_path: str
            string path to a file with a saved energy energy_regression model
        Returns
        -------
            tf.keras model
        """
        return None, None
        if os.path.exists(decoder_model_path) and os.path.exists(
            regression_model_path
        ):
            return tf.keras.models.load_model(
                decoder_model_path
            ), tf.keras.models.load_model(regression_model_path)
        else:
            raise OSError(
                "Path {} or path {} doesn't exist".format(
                    decoder_model_path, regression_model_path
                )
            )

    def predict(self, smi):
        return 100

        m = Chem.MolFromSmiles(smi)
        if m is None:
            raise RuntimeError("Invalid smiles string provided")
        fp = MACCSkeys.GenMACCSKeys(m)
        fp_array = []
        for byte in fp:
            fp_array.append(byte)
        fp_array = np.array(fp_array[1:]).reshape((1, -1))

        return self.energy_regression_model(fp_array).numpy()[0][0]

    def generate(self, n: int, sample_energy: float) -> np.array:
        """Generate molecular descriptors given energy and a number of descriptors
        Parameters
        ----------
        n: int
            Number of descriptor vectors to generate
        sample_energy: float
            Energy threshold to generate compound descriptors with

        Returns
        -------
        np.array
            Descriptors generated

        """
        return np.array([100])
        z = tf.random.normal([n, LATENT_LAYER_SIZE], mean=0.0, stddev=1.0)
        energy = np.zeros(shape=(n, 1))
        energy.fill(sample_energy)
        decoder_input = tf.concat([z, energy], axis=1)
        return self.decoder_model(decoder_input, training=False).numpy()
