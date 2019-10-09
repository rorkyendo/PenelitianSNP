from django.shortcuts import render
from django.template import loader
from django.http import HttpResponse, HttpResponseRedirect
from django.urls import reverse
from __future__ import division

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
    mafthreshold = request.POST.get('mafthreshold')
    mrthreshold = request.POST.get('mrthreshold')
    hwethreshold = request.POST.get('hwethreshold')

    return HttpResponse(mafthreshold)
