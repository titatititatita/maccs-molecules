release: python manage.py makemigrations fingerprints && python manage.py migrate --fake-initial
web: gunicorn maccs.wsgi --log-file -
