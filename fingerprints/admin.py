from django.contrib import admin
from .models import ContactRequest, NewsEntry

admin.site.register(ContactRequest)
admin.site.register(NewsEntry)
