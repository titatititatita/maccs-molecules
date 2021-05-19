from ckeditor_uploader.fields import RichTextUploadingField
from django.db import models


class ContactRequest(models.Model):
    first_name = models.CharField(null=False, max_length=64, blank=False)
    second_name = models.CharField(null=True, max_length=64, blank=True)
    email = models.EmailField(null=False, blank=False)
    message = models.TextField(null=False, blank=False)


class NewsEntry(models.Model):
    title = models.CharField(null=False, blank=False, max_length=64)
    preview_image = models.ImageField()
    content = RichTextUploadingField()
    author = models.CharField(null=False, blank=False, max_length=128)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
