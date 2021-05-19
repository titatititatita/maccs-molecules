import csv
import json

from django.core.exceptions import ViewDoesNotExist
from django.core.paginator import Paginator
from django.http import StreamingHttpResponse
from django.shortcuts import render, redirect, get_object_or_404

from . import utils
from .models import ContactRequest, NewsEntry

INPUT_TEXT = 'text'
INPUT_FILE = 'file'
INPUT_DRAW = 'draw'
allowed_types = {INPUT_TEXT, INPUT_FILE, INPUT_DRAW}
DEFAULT_TYPE = INPUT_TEXT


def resolve_type(input_type):
    if input_type in allowed_types:
        return input_type
    return DEFAULT_TYPE


def predict(request, input_type=DEFAULT_TYPE):
    if request.method == 'POST':
        return make_prediction(request, input_type)
    return render(request, 'fingerprints/prediction.html', {'type': resolve_type(input_type)})


def make_prediction(request, input_type):
    if input_type == INPUT_FILE:
        smiles = request.FILES['smiles'].read().decode("utf-8")
    else:
        smiles = request.POST['smiles']
    result = utils.predict_array_energy(smiles)
    context = {
        'type': resolve_type(input_type),
        'fingerprints': result
    }
    print(context)
    return render(request, 'fingerprints/prediction.html', context)


def generate(request):
    return render(request, 'fingerprints/generate.html')


def index(request):
    context = {
        'news_entries': NewsEntry.objects.all().order_by('-created_at')[:2]
    }
    return render(request, 'fingerprints/index.html', context)


def news(request):
    news_list = NewsEntry.objects.all().order_by('-created_at')
    paginator = Paginator(news_list, 6)

    page_number = request.GET.get('page')
    page_obj = paginator.get_page(page_number)
    return render(request, 'fingerprints/news.html', {'page_obj': page_obj})


def news_details(request, news_id):
    news_entry = get_object_or_404(NewsEntry, pk=news_id)

    return render(request, 'fingerprints/news_entry.html', {'news_entry': news_entry})


def contacts(request):
    return render(request, 'fingerprints/contacts.html')


class Echo:
    """An object that implements just the write method of the file-like
    interface.
    """

    def write(self, value):
        """Write the value by returning it, instead of storing in a buffer."""
        return value


def download_predict(request):
    check_POST(request)

    fingerprints = request.POST['fingerprints'].replace('"', '\\"')
    fingerprints = json.loads(fingerprints.replace("'", '"'))
    pseudo_buffer = Echo()
    writer = csv.writer(pseudo_buffer)
    response = StreamingHttpResponse((writer.writerow(row) for row in fingerprints.items()),
                                     content_type="text/csv")
    response['Content-Disposition'] = 'attachment; filename="result.csv"'
    return response


def generate_download(request):
    check_POST(request)

    n_fingerprints = int(request.POST['n-fingerprints'])
    free_energy = float(request.POST['free-energy'])
    generated = utils.generating(n_fingerprints, free_energy)
    pseudo_buffer = Echo()
    writer = csv.writer(pseudo_buffer)
    response = StreamingHttpResponse((writer.writerow(row) for row in generated),
                                     content_type="text/csv")
    response['Content-Disposition'] = 'attachment; filename="result.csv"'
    return response


def contact(request):
    check_POST(request)

    first_name = request.POST['first_name']
    second_name = request.POST['second_name']
    email = request.POST['email']
    message = request.POST['message']

    contact_request = ContactRequest(first_name=first_name, second_name=second_name, email=email, message=message)
    contact_request.save()

    return redirect('fingerprints:index')


def check_POST(request):
    if request.method != 'POST':
        raise ViewDoesNotExist("Only POST method is supported.")
