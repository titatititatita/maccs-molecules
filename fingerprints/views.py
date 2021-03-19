from django.shortcuts import render
from . import utils


def predict(request):
    if request.method == 'POST':
        return make_prediction(request)
    return render(request, 'fingerprints/prediction.html')


def make_prediction(request):
    f = request.POST['text']
    result = utils.make_fingerprints(f)
    context = {
        'fingerprints': result
    }
    return render(request, 'fingerprints/prediction.html', context)


def index(request):
    return render(request, 'fingerprints/index.html')


def news(request):
    return render(request, 'fingerprints/news.html')


def feedback(request):
    return render(request, 'fingerprints/feedback.html')
