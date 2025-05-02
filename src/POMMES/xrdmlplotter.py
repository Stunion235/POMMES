import xml.dom.minidom as md
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
import csv

#PARAMETERS
#logarithmic y-axis?
logy=True
#plot title; empty string for no title
title="Intensity vs 2Theta"
#x and y axis limits; empty tuple to pick default limits
xlim=()
ylim=()
#list of xrdml files to plot, and matplotlib.pyplot.plot kwargs to customize each one
#for best results, highest offsetY to lowest
files=[{'path':'../../data/EuTaO3Sample.xrdml',
        
        'offsetY': 5,
        'kw':{'color':'blue','label':"Example of offset"}},
       {'path':'../../data/EuTaO3Sample.xrdml',
        
        'offsetY': 0,
        'kw':{'color':'green','label':"Sample"}}]
# files=[{'path':'/Users/oiseau/Documents/Research/quantifyingorder/TS_ETO_0031_XRD.csv',
#         'offsetY': 0,
#         'kw':{'color':'orange','label':"TS_ETO_0031_csv"}}]
#PEAK FITTING SETTINGS
#whether or not to find peaks and print their 2Theta and intensity
printpeaks=True
#whether or not to plot peaks on the graph; looks best with only 1 data file
plotpeaks=True
#find_peaks args, None=default
prominence=6
threshold=None
distance=None
#data smoothing strength: number of points to "merge"
smooth_strength=15

ax = plt.subplots()[1]
ax.set_xlabel("2Theta (°)")
ax.set_ylabel("Intensity")
offset=False
for f in files:
    p=f.get("path")
    try:
        root=md.parse(p)
        counts=root.getElementsByTagName('counts')[0]
        y=np.array(list(map(int,counts.childNodes[0].data.split(" "))))
        twotheta=list(filter(lambda n : (n.getAttribute("axis")=="2Theta"),root.getElementsByTagName('positions')))[0]
        twothetamin=float(twotheta.getElementsByTagName("startPosition")[0].childNodes[0].data)
        twothetamax=float(twotheta.getElementsByTagName("endPosition")[0].childNodes[0].data)
        x=np.linspace(twothetamin, twothetamax, num=len(y))
    except:
        print(p,"could not be parsed as xrdml. Trying csv.")
        def ask_col(prompt):
            print(prompt)
            while True:
                try:
                    c=int(input())
                    if c<0:
                        print("Column index must be positive.")
                    else:
                        return c
                except:
                    print("Column index must be an integer.")
        xcol=ask_col("Which column of the file is 2Theta, where 0 is the leftmost?")
        ycol=ask_col("Which column of the file is Intensity, where 0 is the leftmost?")
        with open(p,'r',newline='') as file:
            read=csv.reader(file)
            x=[]
            y=[]
            for r in read:
                try:
                    r=list(map(float,r))
                    x.append(r[xcol])
                    y.append(r[ycol])
                except:
                    pass
            x=np.array(x)
            y=np.array(y)
        
    try:
        ax.set_xlim(xlim)
    except:
        pass
    try:
        ax.set_ylim(ylim)
    except:
        pass
    offsetY=f.get("offsetY")
    if offsetY!=0:
        offset=True
    if logy:
        ax.semilogy(x,(y+offsetY)*(10**offsetY),**f.get("kw"))
        pass
    else:
        ax.plot(x,y+offsetY,**f.get("kw"))
        pass
    if printpeaks or plotpeaks:
        #smooth to avoid reading background noise as peaks
        ymodified=np.copy(y)
        for i in range(smooth_strength,len(y)-smooth_strength):
            ymodified[i]=max(y[i-smooth_strength:i+smooth_strength])
        peaks=sig.find_peaks(ymodified,prominence=prominence,width=2*smooth_strength,threshold=threshold,distance=distance)[0]
        if printpeaks:
            print("~"*45)
            print("\t\t",f.get("kw").get("label") or "File "+str(files.index(f)+1))
            print("Peak\t\t|\t2Theta(°)\t|\tIntensity")
            print("-"*39)
        for p in range(len(peaks)):
            if plotpeaks:
                plt.vlines(x[peaks[p]], (offsetY)*(10**offsetY) if logy else offsetY, (y[peaks[p]]+offsetY)*(10**offsetY) if logy else y[peaks[p]]+offsetY,'red')
                plt.hlines((y[peaks[p]]+offsetY)*(10**offsetY) if logy else y[peaks[p]]+offsetY, x[max(peaks[p]-4*smooth_strength,0)], x[min(peaks[p]+4*smooth_strength,len(x)-1)],'red')
                plt.text(x[peaks[p]], (y[peaks[p]]+offsetY)*(10**offsetY) if logy else y[peaks[p]]+offsetY, str(p),color='red',ha='center',va='bottom')
            if printpeaks:
                print("%4d\t\t|\t%8.4f\t\t|\t%9d" % (p, x[peaks[p]], y[peaks[p]]))
plt.title(title)
plt.legend()
if offset:
    ax.set_yticks([])