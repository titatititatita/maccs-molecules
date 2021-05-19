from django.urls import path

from . import views

app_name = 'fingerprints'
urlpatterns = [
    path('', views.index, name='index'),
    path('predict/', views.predict, name='predict'),
    path('generate/', views.generate, name='generate'),
    path('generate_download/', views.generate_download, name='generate_download'),
    path('predict/<str:input_type>/', views.predict, name='predict'),
    path('download_predict/', views.download_predict, name='download_predict'),
    path('news/', views.news, name='news'),
    path('news/<int:news_id>', views.news_details, name='news_details'),
    path('contacts/', views.contacts, name='contacts'),
    path('contact/', views.contact, name='contact')
]

