#!/usr/bin/env sh

python manage.py makemigrations fingerprints
python manage.py migrate
