from django.urls import path

from . import views

app_name = 'fingerprints'
urlpatterns = [
    path('', views.index, name='index'),
    path('predict/', views.predict, name='predict'),
    # path('make_prediction/', views.make_prediction, name='make_prediction'),
    path('news/', views.news, name='news'),
    path('feedback/', views.feedback, name='feedback')
]
