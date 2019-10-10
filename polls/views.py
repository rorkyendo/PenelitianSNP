from __future__ import division
from django.shortcuts import render
from django.template import loader
from django.http import HttpResponse, HttpResponseRedirect
from django.urls import reverse
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.contrib import messages

import pandas as pd
import csv
import numpy as np




# Create your views here.
def index(request):
    return render(request,"polls/index.html",
        {
        'title':"Penelitian",
        }
    )

def convert(request):
    if request.method == 'POST':
            mafthreshold = request.POST.get('mafthreshold')
            mrthreshold = request.POST.get('mrthreshold')
            hwethreshold = request.POST.get('hwethreshold')
            file = request.FILES['csv_file']
            fs = FileSystemStorage(location='static/datasnp')
            filename = fs.save(file.name, file)
            file_url = fs.url(filename)

            messages.success(request, 'Data berhasil diupload '+file_url)

            return HttpResponseRedirect(reverse('index'))
