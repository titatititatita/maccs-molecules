from django.apps import AppConfig

from fingerprints.ml import models


class FingerprintsConfig(AppConfig):
    name = 'fingerprints'
    ml_model = models.TF_Models()
