import readCIF
import numpy as np
import XRD
import tkinter as tk
from tkinter import filedialog as fd
from tkinter import messagebox as mb
from tkinter import ttk
has_plt=True
try:
    import matplotlib.pyplot as plt
except:
    has_plt=False
from scipy.interpolate import interp1d
#PARAMETERS
# #CIF file of disordered structure
# disordered_file='/Users/oiseau/Documents/Research/quantifyingorder/OneDrive_1_9-5-2024/Ba2ScNbO6/Ba2ScNbO6-occ5.cif'
# #CIF file of ordered structure
# ordered_file='/Users/oiseau/Documents/Research/quantifyingorder/OneDrive_1_9-5-2024/Ba2ScNbO6/Ba2ScNbO6-occ10.cif'
#The amounts of order you want to plot (0=disordered, 1=ordered)
orderlevels=np.linspace(0,1,num=21)
#X-ray wavelength in your measurement
# wavelength_in_Å=1.5406
#hkl is an array of arrays of 2 [h,k,l] arrays. Each element of the outer array is a ratio
#with its first element is the numerator of a peak ratio and the second is the denominator.
#e.g. [[[0,0,1],[0,0,2]],[[0,0,1],[0,0,4]]] will plot the ratios 001/002 and 001/004.
# hkl=[[[0,0,1],[0,0,2]],[[0,0,1],[0,0,4]],[[0,0,1],[0,0,6]]]
# hkl=[[[1,1,1],[2,2,2]],[[1,1,1],[4,4,4]],[[3,3,3],[2,2,2]],[[3,3,3],[4,4,4]]]
#Array of colors in which to plot the ratios in the same order. If there are more ratios
#than colors, they will cycle.
# colors=['red','blue','green','violet']
#y-axis limits on plot
# ymin,ymax=(-.1,2)
#Your experimentally measured peak ratio values in the same order as hkl
# experimental_ratios=[21/159862,21/61577,9/159862,9/61577]
# We'd never be doing this
# preload_matrices=False
# We'd always be doing this
# interpolate_experimental_values=True
# When you think about it,neither of the below should happen. We don't preload and we don't read the given filepaths.
# if preload_matrices:
#     from preload import sites, disordered, ordered
# else:
#     sites, disordered, ordered = readCIF.make_cif_matrices(disordered_file,ordered_file)

class PommesGUI(tk.Frame):
    def __init__(self, root=None):
        self.root=root
        tk.Frame.__init__(self, root)
        self.grid()
        self.createWidgets()
        self.disorderedfile=None
        self.orderedfile=None
        self.orderlevels=np.linspace(0,1,num=21)
        
    def createWidgets(self):
        # self.medialLabel=tk.Label(self.root,text="Hello World")
        # self.medialLabel.config(bg="#0f0")
        # self.medialLabel.grid()
        # self.quitButton=tk.Button(self.root,text="Quit", command=self.quit)
        # self.quitButton.grid()
        #my tests
        # self.test=tk.Entry(self.root,bg="#000",fg="#00f")
        # self.test.grid()
        #what I think is ACTUALLY needed
        #TODO make sure things that aren't needed as variables, like submit button, aren't stored to save space
        self.filelabel=tk.Label(self.root,text="Parameters",font=("",16))
        self.filelabel.grid(row=0,column=0,columnspan=2)
        self.disorderedfilebutton=tk.Button(self.root,command=self.ask_disordered_file,text="Select disordered CIF")
        self.disorderedfilebutton.grid(row=1,column=0)
        self.disorderedfilename=tk.Label(self.root,text="No file selected")
        self.disorderedfilename.grid(row=2,column=0)
        self.orderedfilebutton=tk.Button(self.root,command=self.ask_ordered_file,text="Select ordered CIF")
        self.orderedfilebutton.grid(row=3,column=0)
        self.orderedfilename=tk.Label(self.root,text="No file selected")
        self.orderedfilename.grid(row=4,column=0)
        self.wllabel=tk.Label(self.root,text="Wavelength (Å)")
        self.wllabel.grid(row=1,column=1)
        self.wl=tk.Entry(self.root)
        self.wl.grid(row=2,column=1)
        self.wl.insert(0,'1.5406')
        if has_plt:
            self.collabel=tk.Label(self.root,text='Matplotlib colors, comma-separated; will cycle')
            self.collabel.grid(row=3,column=1)
            self.col=tk.Entry(self.root)
            self.col.grid(row=4,column=1)
            self.col.insert(0,'red,blue,green')
        self.line1=ttk.Separator(self.root)
        self.line1.grid(row=5,column=0,columnspan=2,sticky="we")
        self.title1=tk.Label(self.root,text='Single Ratio Mode',font=("",16))
        self.title1.grid(row=6,column=0,columnspan=2)
        self.peaklabel=tk.Label(self.root,text='Miller indices of peaks to compare. 3 indices separated by spaces, e.g. "1 1 1"',font=("",12))
        self.peaklabel.grid(row=7,column=0,columnspan=2)
        self.m1l=tk.Label(self.root,text="Numerator")
        self.m1l.grid(row=8,column=0)
        self.m1l=tk.Label(self.root,text="Denominator")
        self.m1l.grid(row=8,column=1)
        self.miller1=tk.Entry(self.root)
        self.miller1.grid(row=9,column=0)
        self.miller1.insert(0,'0 0 1')
        self.miller2=tk.Entry(self.root)
        self.miller2.grid(row=9,column=1)
        self.miller2.insert(0,'0 0 2')
        self.intlabel=tk.Label(self.root,text='Intensities of peaks to compare',font=("",16))
        self.intlabel.grid(row=10,column=0,columnspan=2)
        self.i1l=tk.Label(self.root,text="Numerator")
        self.i1l.grid(row=11,column=0)
        self.i1l=tk.Label(self.root,text="Denominator")
        self.i1l.grid(row=11,column=1)
        self.int1=tk.Entry(self.root)
        self.int1.grid(row=12,column=0)
        self.int1.insert(0,63)
        self.int2=tk.Entry(self.root)
        self.int2.grid(row=12,column=1)
        self.int2.insert(0,2913)
        self.run=tk.Button(self.root,text="Calculate",command=lambda: self.validate(False) and self.calc(False))
        self.run.grid(row=13,column=0,columnspan=2)
        self.line2=ttk.Separator(self.root)
        self.line2.grid(row=14,column=0,columnspan=2,sticky="we")
        self.title2=tk.Label(self.root,text='Multi-Ratio Mode',font=("",16))
        self.title2.grid(row=15,column=0,columnspan=2)
        self.peak2label=tk.Label(self.root,text='Numerator peak Miller indices,Denominator peak Miller indices,Numerator peak intensity,Denominator peak intensity',font=("",11),padx=15)
        self.peak2label.grid(row=16,column=0,columnspan=2)
        self.millers=tk.Text()
        self.millers.insert("1.0","0 0 1,0 0 2,63,2913\n0 0 1,0 0 4,63,1873")
        self.millers.grid(row=17,column=0,columnspan=2)
        self.run2=tk.Button(self.root,text="Calculate",command=lambda: self.validate(True) and self.calc(True))
        self.run2.grid(row=18,column=0,columnspan=2)
    
    def ask_disordered_file(self):
        self.disorderedfile=fd.askopenfile()
        self.disorderedfilename.config(text=self.disorderedfile.name if self.disorderedfile != None else "No file selected")
        self.disorderedlines=self.disorderedfile.readlines()
        return self.disorderedfile
    
    def ask_ordered_file(self):
        self.orderedfile=fd.askopenfile()
        self.orderedfilename.config(text=self.orderedfile.name if self.orderedfile != None else "No file selected")
        self.orderedlines=self.orderedfile.readlines()
        return self.orderedfile
    
    def validate(self,multi):
        errors=""
        if self.disorderedfile==None:
            self.disorderedfilename.config(text="Need a file!")
            errors+="No disordered file.\n"
        if self.orderedfile==None:
            self.orderedfilename.config(text="Need a file!")
            errors+="No ordered file.\n"
        try:
            float(self.wl.get())
        except:
            errors+="Wavelength must be a number.\n"
        if has_plt:
            try:
                cs=self.col.get().split(',')
                f,a=plt.subplots()
                for c in cs:
                    plt.plot(0,0,color=c)
                plt.clf()
            except ValueError as v:
                errors+="Color unrecognized by matplotlib: "+str(v).split("'")[1]+"\n"
        if multi:
            try:
                rows=self.millers.get("1.0","end-1c").splitlines()
                map(lambda p:list(map(lambda h:list(map(int,h.split())),p)),map(lambda a:a.split(',')[0:2],rows))
                map(lambda p:float(p[0])/float(p[1]),map(lambda l:l.split(',')[2:4],rows))
            except:
                errors+="Multi-ratio entry is not formatted correctly.\n"
        else:
            try:
                float(self.int1.get())
            except:
                errors+="Numerator intensity must be a number.\n"
            try:
                float(self.int2.get())
            except:
                errors+="Denominator intensity must be a number.\n"
        if errors!="":
            print("Errors:",errors,sep='\n')
            mb.showerror("Errors",'Errors:\n'+errors)
        return errors==""
    
    def calc(self,multi):
        disorderedABC = np.array(readCIF.get_lattice_constants(self.disorderedfile.name,True))
        orderedABC = np.array(readCIF.get_lattice_constants(self.orderedfile.name,True))
        abc_matrix=np.array([orderedABC,disorderedABC])
        sites, disordered, ordered = readCIF.make_cif_matrices(self.disorderedfile.name,self.orderedfile.name,True)
        colors=self.col.get().split(',') if self.col.get().strip() != '' else ['red','blue','green']
        wavelength_in_Å=float(self.wl.get())
        if has_plt: 
            fig,ax1=plt.subplots()
        estimates=[]
        experimental_ratios=[]
        hkl=[]
        if multi:
            rows=self.millers.get("1.0","end-1c").splitlines()
            hkl.extend(map(lambda p:list(map(lambda h:list(map(int,h.split())),p)),map(lambda a:a.split(',')[0:2],rows)))
            experimental_ratios.extend(map(lambda p:float(p[0])/float(p[1]),map(lambda l:l.split(',')[2:4],rows)))
            # raise Warning(hkl)
        else:
            hkl.append([list(map(int,self.miller1.get().split())),list(map(int,self.miller2.get().split()))])
            experimental_ratios.append(float(self.int1.get())/float(self.int2.get()))
        for i in range(len(hkl)):
            e=[experimental_ratios[i]]
            p=hkl[i][0]
            p2=hkl[i][1]
            structurefactors=[]
            intensities=[]
            structurefactors2=[]
            intensities2=[]
            ratios=[]
            for o in self.orderlevels:
                #Get structure factors of both structures, compute intensity, and find ratio
                d=XRD.plane_spacing(p,[o,1-o]@abc_matrix)
                structurefactors.append(abs(XRD.structure_factor(XRD.form_factors(sites[:,1], d), disordered, ordered, p, o)))
                intensities.append(XRD.mult(p)*XRD.intensity(structurefactors[-1],d,wavelength_in_Å))
                d=XRD.plane_spacing(p2,[o,1-o]@abc_matrix)
                structurefactors2.append(abs(XRD.structure_factor(XRD.form_factors(sites[:,1], d), disordered, ordered, p2, o)))
                intensities2.append(XRD.mult(p2)*XRD.intensity(structurefactors2[-1],d,wavelength_in_Å))
                ratios.append(intensities[-1]/intensities2[-1])
            has_plt and ax1.plot(self.orderlevels*100,ratios,label=str(p[0])+str(p[1])+str(p[2])+"/"+str(p2[0])+str(p2[1])+str(p2[2]),color=colors[i%len(colors)])
            inter_func1=interp1d(ratios,self.orderlevels*100,kind='nearest')
            #Interpolate
            interp_order=inter_func1(e)
            estimates.append(interp_order)
            has_plt and ax1.scatter(interp_order, e, color=colors[i%len(colors)])
        avg=sum(estimates)/len(estimates)
        mb.showinfo('Result','interpolated order=' + str(avg/100))
        if has_plt:
            ax1.xaxis.set_ticks_position('both')
            ax1.yaxis.set_ticks_position('both')
            ax1.tick_params(axis='both', direction='in')
            ax1.set_ylabel('Intensity Ratio')
            ax1.set_xlabel('% order')
            ax1.set_xlim(-2,102)
            yl=ax1.get_ylim()
            plt.vlines(avg,yl[0],yl[1],linestyles='dashed',colors='black')
            plt.title('Intensity Ratio vs. Order')
            plt.legend()
            plt.show()
        raise Warning("unfinished")
            
def onclose():
    if mb.askyesno("Quit","Sure? Unsaved data will be lost."):
        root.destroy()
root = tk.Tk()
root.protocol("WM_DELETE_WINDOW", onclose)
root.title("POMMES")
app = PommesGUI(root)
app.mainloop()