from django.shortcuts import render
from django.template import loader
from django.http import HttpResponse

# Create your views here.
def index(request):
    return render(request,"polls/index.html",
        {
        'title':"Penelitian",
        }
    )

def convert(arg):
    return