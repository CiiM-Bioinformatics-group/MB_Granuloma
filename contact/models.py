from django.db import models

# Create your models here.


class Contact(models.Model):
    company_name = models.CharField(max_length=255)
    full_name = models.CharField(max_length=255)
    email_address = models.EmailField()

    def __str__(self):
        return f'{self.full_name} ({self.company_name})'


