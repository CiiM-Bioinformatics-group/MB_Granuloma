from django.contrib import admin



# contact/admin.py
from django.contrib import admin
from .models import Contact

@admin.register(Contact)
class ContactAdmin(admin.ModelAdmin):
    list_display = ('full_name', 'company_name', 'email_address')
    search_fields = ('full_name', 'company_name', 'email_address')
