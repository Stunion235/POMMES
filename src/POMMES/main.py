import readCIF
import numpy as np
import XRD
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
#PARAMETERS
#CIF file of disordered structure
disordered_file='../../data/EuTaO3disordered.cif'
#CIF file of ordered structure
ordered_file='../../data/EuTaO3ordered.cif'
#The amounts of order you want to plot (0=disordered, 1=ordered)
orderlevels=np.linspace(0,1,num=21)
#X-ray wavelength in your measurement
wavelength_in_Å=1.5406
#hkl is an array of arrays of 2 [h,k,l] arrays. Each element of the outer array is a ratio
#with its first element is the numerator of a peak ratio and the second is the denominator.
#e.g. [[[0,0,1],[0,0,2]],[[0,0,1],[0,0,4]]] will plot the ratios 001/002 and 001/004.
# hkl=[[[0,0,1],[0,0,2]],[[0,0,1],[0,0,4]],[[0,0,1],[0,0,6]]]
hkl=[[[0,0,1],[0,0,2]],[[0,0,1],[0,0,4]]]
#Array of colors in which to plot the ratios in the same order. If there are more ratios
#than colors, they will cycle.
colors=['magenta','purple','green']
#y-axis limits on plot
ymin,ymax=(-.1,2)
#Your experimentally measured peak ratio values in the same order as hkl
experimental_ratios=[63/2913,63/1873]
#If true, load matrices from `preload.py` instead of reading the above CIF files.
#To save time reading the CIF files every time the code must be run, see `readCIF.py`
preload_matrices=False
#If true, interpolate experimental ratios above to estimate order; else just plot ratios
interpolate_experimental_values=True
if preload_matrices:
    from preload import sites, disordered, ordered
else:
    sites, disordered, ordered = readCIF.make_cif_matrices(disordered_file,ordered_file)
#Get lattice contstants of both structures
disorderedABC = np.array(readCIF.get_lattice_constants(disordered_file))
orderedABC = np.array(readCIF.get_lattice_constants(ordered_file))
abc_matrix=np.array([orderedABC,disorderedABC])
fig,ax1=plt.subplots()
estimates=[]
for i in range(len(hkl)):
    e=[experimental_ratios[i]]
    p=hkl[i][0]
    p2=hkl[i][1]
    structurefactors=[]
    intensities=[]
    structurefactors2=[]
    intensities2=[]
    ratios=[]
    for o in orderlevels:
        #Get structure factors of both structures, compute intensity, and find ratio
        d=XRD.plane_spacing(p,[o,1-o]@abc_matrix)
        structurefactors.append(abs(XRD.structure_factor(XRD.form_factors(sites[:,1], d), disordered, ordered, p, o)))
        intensities.append(XRD.mult(p)*XRD.intensity(structurefactors[-1],d,wavelength_in_Å))
        d=XRD.plane_spacing(p2,[o,1-o]@abc_matrix)
        structurefactors2.append(abs(XRD.structure_factor(XRD.form_factors(sites[:,1], d), disordered, ordered, p2, o)))
        intensities2.append(XRD.mult(p2)*XRD.intensity(structurefactors2[-1],d,wavelength_in_Å))
        ratios.append(intensities[-1]/intensities2[-1])
    # Unused code to plot the structure factors themselves
    # structurefactors=list(map(lambda o:abs(XRD.structure_factor(XRD.form_factors(sites[:,1], p, XRD.plane_spacing(p,[o,1-o]@abc_matrix)), disordered, ordered, p, o)),orderlevels))
    # ax1.plot(orderlevels*100,structurefactors,label=str(p[0])+str(p[1])+str(p[2]),color=colors[(i+1)%len(colors)])
    # ax1.plot(orderlevels*100,structurefactors2,label=str(p2[0])+str(p2[1])+str(p2[2]),color=colors[(i+2)%len(colors)])
    ax1.plot(orderlevels*100,ratios,label=str(p[0])+str(p[1])+str(p[2])+"/"+str(p2[0])+str(p2[1])+str(p2[2]),color=colors[i%len(colors)])
    if interpolate_experimental_values:
        inter_func1 = interp1d(ratios,orderlevels*100, kind='nearest')
        #Interpolate
        interp_order=inter_func1(e)
        estimates.append(interp_order)
        ax1.scatter(interp_order, e, color=colors[i%len(colors)])
        print('interpolated order=' + str(interp_order/100) + ', ratio=' + str(e))
#Make plot
if interpolate_experimental_values:
    avg=sum(estimates)/len(estimates)
    plt.vlines(avg,ymin,ymax,linestyles='dashed',colors='black')
    print('Average order estimate:',avg/100)
ax1.xaxis.set_ticks_position('both')
ax1.yaxis.set_ticks_position('both')
ax1.tick_params(axis='both', direction='in')
ax1.set_ylabel('Intensity Ratio')
ax1.set_xlabel('% order')
ax1.set_xlim(-2,102)
ax1.set_ylim(ymin,ymax)
plt.title('Intensity Ratio vs. Order')
plt.legend()
plt.show()
