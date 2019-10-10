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
            file = request.FILES['csv_file']
            fs = FileSystemStorage(location='static/datasnp')
            filename = fs.save(file.name, file)
            file_url = fs.url(filename)

            df=pd.read_csv('.'+settings.STATIC_URL+'datasnp/'+filename) #import datasnp.csv
            df.columns=["ID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10","V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21"]
            mrthreshold = float(request.POST.get('mrthreshold')) #nilai mrthreshold diinput user
            mr = pd.DataFrame(columns=['SNP', 'MissingRate'])
            for i in range(0, int(len(df))):  #counter dari baris (ID-SNP) 1 s.d. banyak data SNP (1546)
                missingrt = (df.iloc[i] == '--').sum()/len(df.iloc[i]) #menghitung missing rate setiap SNP
                mr = mr.append({'SNP': df.ID[i], 'MissingRate': missingrt}, ignore_index=True) #memasukkan nilai missing rate ke list mr

            filteredSNP = mr['MissingRate'] < mrthreshold  # Create dataframe with TRUE if MissingRate is less than 0.05
            mrSNP=df[filteredSNP]

            mrSNP.to_csv('static/datasnp/drop_mr.csv', index=False)

            df=pd.read_csv('.'+settings.STATIC_URL+'datasnp/drop_mr.csv') #import datasnp.csv
            df.columns=["ID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10","V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21"]

            mafSNP = pd.DataFrame(columns=['SNP', 'MAF'])
            mafthreshold = float(request.POST.get('mafthreshold'))
            for i in range(0, int(len(df))):  #counter dari baris (ID-SNP) 1 s.d. banyak data SNP (1546)
                a=df.iloc[i][:].str.count("A").sum()
                t=df.iloc[i][:].str.count("T").sum()
                g=df.iloc[i][:].str.count("G").sum()
                c=df.iloc[i][:].str.count("C").sum()
                allelfreq=[a,t,g,c]
                allelfreq.sort(reverse=True)
                maf=allelfreq[1]/sum(allelfreq)
                mafSNP = mafSNP.append({'SNP': df.ID[i], 'MAF': maf}, ignore_index=True)

            filtermaf = mafSNP['MAF'] > mafthreshold  # Create variable with TRUE if MAF is less than 0.05
            mafSNP=df[filtermaf]
            mafSNP.to_csv('static/datasnp/drop_maf.csv', index=False)

            dfr=pd.read_csv('.'+settings.STATIC_URL+'datasnp/drop_maf.csv') #import datasnp.csv
            dfr.columns=["ID", "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10","V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21"]
            df=dfr
            df = df.drop('ID', 1)
            listhwe = []
            hwethreshold = float(request.POST.get('hwethreshold')) #input nili hwe dari user
            for i in range(0, int(len(df))):  #counter dari baris (ID-SNP) 1 s.d. banyak data SNP (1546)
                a=df.iloc[i][:].str.count("A").sum() #menghitung jumlah allel A yang muncul di data SNP ke-i
                t=df.iloc[i][:].str.count("T").sum() #menghitung jumlah allel T yang muncul di data SNP ke-i
                g=df.iloc[i][:].str.count("G").sum() #menghitung jumlah allel G yang muncul di data SNP ke-i
                c=df.iloc[i][:].str.count("C").sum() #menghitung jumlah allel C yang muncul di data SNP ke-i
                allelfreq=[a,t,g,c]
                allelfreq.sort(reverse=True) #mengurutkan
                s=sum(allelfreq)
                allel=pd.DataFrame(columns=['p','q'])
                allel = allel.append({'p': (allelfreq[0])/s, 'q': (allelfreq[1])/s}, ignore_index=True)
                AAe = allel.p**2
                BBe = allel.q**2
                ABe = 2 * allel.p * allel.q

                gen_count=df.iloc[i].value_counts() #menghitung frekuensi kemunculan allel.
                a = gen_count.rename_axis('genotype').reset_index(name='counts')
                a = a[a.genotype != "--"].reset_index(drop=True)
                if len(a) > 2:
                    if a.genotype[0][0] == a.genotype[0][1]: #apakah allel frekuensi terbanyak pertama homozigot dominan?
                        AAo = a.counts[0]/a.counts.sum(axis=0) # homozygot dominan observed
                        if a.genotype[1][0] == a.genotype[1][1]: # apakah di urutan kedua homozigot?
                            BBo = a.counts[1]/a.counts.sum(axis=0) #homozygot resesif observed
                            ABo = a.counts[2]/a.counts.sum(axis=0) #heterozygot observed
                        else:
                            ABo = a.counts[1]/a.counts.sum(axis=0) #heterozygot observed
                            BBo = a.counts[2]/a.counts.sum(axis=0) #homozygot resesif observed
                    else:
                        ABo = a.counts[0]/a.counts.sum(axis=0)
                        AAo = a.counts[1]/a.counts.sum(axis=0)
                        BBo = a.counts[2]/a.counts.sum(axis=0)

                if len(a) == 2:
                    if a.genotype[0][0] == a.genotype[0][1]: #apakah allel frekuensi terbanyak pertama homozigot dominan?
                        AAo = a.counts[0]/a.counts.sum(axis=0) # homozygot dominan observed
                        if a.genotype[1][0] == a.genotype[1][1]: # apakah di urutan kedua homozigot?
                            BBo = a.counts[1]/a.counts.sum(axis=0) #homozygot resesif observed
                        else:
                            ABo = a.counts[1]/a.counts.sum(axis=0) #heterozygot observed
                    else:
                        ABo = a.counts[0]/a.counts.sum(axis=0)
                        AAo = a.counts[1]/a.counts.sum(axis=0)
                hwe=((AAo-AAe)**2/AAe)+((ABo-ABe)**2/ABe)+((BBo-BBe)**2/BBe)
                listhwe.append(hwe)

            listhwe = pd.DataFrame(listhwe)
            listhwe.columns=["HWE"]
            filterhwe = listhwe['HWE'] > hwethreshold
            hweSNP=dfr[filterhwe]
            hweSNP.to_csv('static/datasnp/drop_hwe.csv', index=False)

            df=pd.read_csv('.'+settings.STATIC_URL+'datasnp/drop_hwe.csv')
            df.columns=["ID","V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10","V11", "V12", "V13", "V14", "V15", "V16", "V17", "V18", "V19", "V20", "V21"]

            for i in range(0, int(len(df))):  #counter dari baris (ID-SNP) 1 s.d. panjang SNP
                allel=df.iloc[i].drop('ID').value_counts()
                a = allel.rename_axis('allel').reset_index(name='counts')
                a = a[a.allel != "--"].reset_index(drop=True)
                if len(a) > 2:
                    if a.allel[0][0] == a.allel[0][1]: #apakah allel frekuensi terbanyak pertama homozigot?
                        df.iloc[i]=df.iloc[i].replace([a.allel[0][0]+a.allel[0][1]], [0]) # urutan pertama homozgt => 0
                        if a.allel[1][0] == a.allel[1][1]: # apakah di urutan kedua homozigot?
                            df.iloc[i]=df.iloc[i].replace([a.allel[1][0]+a.allel[1][1]], [2]) #jika ya, ganti urutan kedua jadi 2
                            df.iloc[i]=df.iloc[i].replace([a.allel[2][0]+a.allel[2][1]], [1]) #urutan 3 jadi 1
                        else:
                            df.iloc[i]=df.iloc[i].replace([a.allel[1][0]+a.allel[1][1]], [1]) #jika tidak ganti 1
                            df.iloc[i]=df.iloc[i].replace([a.allel[2][0]+a.allel[2][1]], [2])
                    else:
                        df.iloc[i]=df.iloc[i].replace([a.allel[0][0]+a.allel[0][1]], [1])
                        df.iloc[i]=df.iloc[i].replace([a.allel[1][0]+a.allel[1][1]], [0])
                        df.iloc[i]=df.iloc[i].replace([a.allel[2][0]+a.allel[2][1]], [2])
                if len(a) == 2:
                    if a.allel[0][0] == a.allel[0][1]: #apakah allel frekuensi terbanyak pertama homozigot?
                        df.iloc[i]=df.iloc[i].replace([a.allel[0][0]+a.allel[0][1]], [0]) # urutan pertama homozgt => 0
                        if a.allel[1][0] == a.allel[1][1]: # apakah di urutan kedua homozigot?
                            df.iloc[i]=df.iloc[i].replace([a.allel[1][0]+a.allel[1][1]], [2]) #jika ya, ganti jadi 2
                        else:
                            df.iloc[i]=df.iloc[i].replace([a.allel[1][0]+a.allel[1][1]], [1]) #jika tidak ganti 1
                    else:
                        df.iloc[i]=df.iloc[i].replace([a.allel[0][0]+a.allel[0][1]], [1])
                        df.iloc[i]=df.iloc[i].replace([a.allel[1][0]+a.allel[1][1]], [0])
                if len(a) == 1:
                    if a.allel[0][0] == a.allel[0][1]:
                        df.iloc[i]=df.iloc[i].replace([a.allel[0][0]+a.allel[0][1]], [0])
                    else:
                        df.iloc[i]=df.iloc[i].replace([a.allel[0][0]+a.allel[0][1]], [1])
            df.to_csv('static/datasnp/snp_encoded.csv', index=False)

            messages.success(request, 'Silahkan Download data.')

            return HttpResponseRedirect(reverse('index'))
