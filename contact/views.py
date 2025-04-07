from django.shortcuts import render

# Create your views here.




def home(request):
    return render(request, 'home.html')

def about(request):
    return render(request, 'about.html')

def contact(request):
    return render(request, 'contact.html')

def datasets(request):
    return render(request, 'datasets.html')

def lungs(request):
    return render(request, 'lungs.html')



# contact/views.py

from rest_framework.views import APIView
from rest_framework.response import Response
from .models import Contact
from .serializers import ContactSerializer

class ContactCreate(APIView):
    def post(self, request, format=None):
        print("Received Data:", request.data)  # Debugging
        serializer = ContactSerializer(data=request.data)
        if serializer.is_valid():
            serializer.save()
            return Response(serializer.data, status=201)
        print("Errors:", serializer.errors)  # Debugging
        return Response(serializer.errors, status=400)

    



