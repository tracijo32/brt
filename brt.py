'''
brt.py
Created by Traci Johnson
tljohn@umich.edu
version 1.1 -- 7.15.2017

Buffalo Ray Trace (BRT) is an interactive GUI for plotting image 
predictions for a lens model. BRT is written in Python, utilizing the
tkinter GUI library, the matplotlib plotting library, the astropy
library of tools for astrophysical data analysis. All are available
through Anaconda.

The only required inputs for BRT are the x and y deflection files (FITS),
in units of arcseconds, and a PNG color image or FITS image of the field of view.
These two sets of inputs need to have the same field of view. The program provides
helper functions to create these files.

VERSION HISTORY:
1.1 -- 7.15.2017: Fixed minor bugs. Fixed bug with computing dls/ds using proper
                  cosmology. Added postage stamp feature. Added feature that inserts
                  the redshift of the selected arcs from the arc list into the boxes
                  for ray tracing and plotting the critical curve.

'''

import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import os
import sys
if sys.version_info[0] < 3:
    from Tkinter import *
    import tkMessageBox as tkMB
else:
    from tkinter import *
    from tkinter import messagebox as tkMB
import pickle
from astropy.io import fits
from astropy.wcs import WCS
from astropy.cosmology import FlatLambdaCDM
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from itertools import cycle
import warnings
import time
import datetime
import platform
from PIL import Image

def dlsds(zl,zs,Om0=0.3,H0=70):
    cosmo = FlatLambdaCDM(Om0=Om0,H0=H0)

    ratio = np.zeros_like(zs)
    for i in range(len(zs)):
        dls = cosmo.angular_diameter_distance_z1z2(zl,zs[i]).value
        ds = cosmo.angular_diameter_distance(zs[i]).value
        ratio[i] = dls/ds

    return ratio

def predict_images(xim,yim,deflectx,deflecty,dlsds=1,maxdist=0.5):
    dims = deflectx.shape
    source_x = np.zeros_like(deflectx)
    source_y = np.zeros_like(deflecty)
    if dims[0] == dims[1]:
        for i in range(dims[0]):
            source_x[:,i] = i + 1 - deflectx[:,i]*dlsds
            source_y[i,:] = i + 1 - deflecty[i,:]*dlsds
    else:
        for j in range(dims[0]): source_x[:,j] = j + 1 - deflectx[:,j]*dlsds
        for k in range(dims[1]): source_y[k,:] = k + 1 - deflecty[k,:]*dlsds
    
    xs = source_x[int(np.round(yim))-1,int(np.round(xim))-1]
    ys = source_y[int(np.round(yim))-1,int(np.round(xim))-1]
    
    d = np.sqrt((source_x-xs)**2+(source_y-ys)**2)
    indices = np.where(d<maxdist)
    
    ximp = []
    yimp = []
    for i,j in zip(indices[1],indices[0]): ximp.append(i+1),yimp.append(j+1)
    ximp = np.array(ximp)
    yimp = np.array(yimp)
    
    return ximp, yimp

def update(f):
    data = fits.getdata(f)
    h = fits.getheader(f)
    
    h['CRPIX1'] += 0.5
    h['CRPIX2'] += 0.5
    
    fits.writeto(f,data,header=h,clobber=True)

class StartWindow:
    def __init__(self):
        self.root = Tk()
        self.root.wm_title("Start up")
        self.root.geometry("380x380")

        titleFrame = Frame(self.root)
        titleFrame.pack()
        title = Label(titleFrame,text="Buffalo Ray Trace",fg='blue')
        title.config(font=("Helvetica", 24))
        title.pack()

        Label(titleFrame,text='Enter model parameters',fg='red').pack()

        entryFrame = Frame(self.root)
        entryFrame.pack()
        
        Label(entryFrame, text = "Cluster redshift: ").grid(row=0, column=0,sticky=E)
        self.entry_zl = Entry(entryFrame,width=15)
        self.entry_zl.grid(row=0, column=1)
        
        Label(entryFrame, text = "Image file: ").grid(row=1, column=0,sticky=E)
        self.entry_imagefile = Entry(entryFrame,width=15)
        self.entry_imagefile.grid(row=1, column=1)
        Button(entryFrame, text='Create',command=self.tiff_window,padx=5,pady=5).grid(row=1,column=2)
        
        Label(entryFrame, text = "X deflection file: ").grid(row=2, column=0,sticky=E)
        self.entry_dxfile = Entry(entryFrame,width=15)
        self.entry_dxfile.grid(row=2, column=1)
        
        Label(entryFrame, text = "Y deflection file: ").grid(row=3, column=0,sticky=E)
        self.entry_dyfile = Entry(entryFrame,width=15)
        self.entry_dyfile.grid(row=3, column=1)
        
        Label(entryFrame, text = "Deflection file redshift: ").grid(row=4, column=0,sticky=E)
        self.entry_dz = Entry(entryFrame,width=15)
        self.entry_dz.grid(row=4, column=1)
        
        self.check = IntVar()
        self.check.set(0)
        Checkbutton(entryFrame, variable=self.check, onvalue=1, offvalue=0, text='Deflection files are dls/ds = 1').grid(row=5,columnspan=2)
        
        Label(entryFrame,text='Enter cosmological parameters',fg='red').grid(row=6,columnspan=2)
        
        Label(entryFrame, text = "Omega M (1 - Omega L): ").grid(row=7, column=0,sticky=E)
        self.entry_Om0 = Entry(entryFrame,width=15)
        self.entry_Om0.grid(row=7, column=1)
        self.entry_Om0.insert(0,'0.3')
        
        Label(entryFrame, text = "Hubble constant (km/s/Mpc): ").grid(row=8, column=0,sticky=E)
        self.entry_H0 = Entry(entryFrame,width=15)
        self.entry_H0.grid(row=8, column=1)
        self.entry_H0.insert(0,'70.0')

        submitFrame = Frame(self.root)
        submitFrame.pack()

        Button(submitFrame, text = "Enter", command = self.getParams,padx=5,pady=5).pack()

        Label(submitFrame, text='Or').pack()

        Button(submitFrame,text="Load previous model",command=self.loadPrevious,padx=5,pady=5).pack()

        self.root.mainloop()

    def tiff_window(self):
        self.toplevel = Toplevel(self.root)
        Label(self.toplevel, text='Open a fits file in ds9 (or three if RGB) and\nscale to the desired output.',fg='blue').grid(row=0,columnspan=2)
        Label(self.toplevel, text='TIFF file name: ').grid(row=1,column=0,sticky=E)
        self.file_entry = Entry(self.toplevel,width=10)
        self.file_entry.grid(row=1,column=1)
        Button(self.toplevel,text='Write TIFF',command=self.createTiff,padx=5,pady=5).grid(row=2,columnspan=2)

    def createTiff(self):
        tiffname = self.file_entry.get()
        os.system('xpaset -p ds9 export tiff '+os.getcwd()+'/'+tiffname+' none')
        fitsname = os.popen('xpaget ds9 file').readlines()[0].rsplit()[0]
    
        htiff = fits.getheader(fitsname)
        wtiff = WCS(htiff)
        pickle.dump(wtiff,open(tiffname+'.wcs','wb'))
        self.entry_imagefile.delete(0,END)
        self.entry_imagefile.insert(0,tiffname)
        self.toplevel.destroy()


    def getParams(self):
        self.zl = self.entry_zl.get()
        self.imagefile = self.entry_imagefile.get()
        self.dxfile = self.entry_dxfile.get()
        self.dyfile = self.entry_dyfile.get()
        self.dz = self.entry_dz.get()
        self.isInf = self.check.get()
        self.Om0 = self.entry_Om0.get()
        self.H0 = self.entry_H0.get()
    
        errors = []
        try:
            self.zl = float(self.zl)
        except ValueError or self.zl < 0:
            errors.append('Cluster redshift must be a number > 0.')
    
        for file in [self.imagefile, self.dxfile, self.dyfile]:
            if not os.path.isfile(file):
                errors.append('File "'+file+ '" does not exist.')

        if self.isInf == 0:
            try:
                self.dz = float(self.dz)
            except ValueError or self.dz < self.zl:
                    errors.append('Deflect file redshift must be a number > cluster redshift.')
    
        try:
            self.Om0 = float(self.Om0)
        except ValueError or Om0 < 0 or Om0>1:
            errors.append('Omega M must be a number between 0 and 1')
    
        try:
            self.H0 = float(self.H0)
        except ValueError or self.H0 < 0 or self.H0 > 100:
            errors.append('H0 must be a number between 0 and 100.')
    
        if len(errors) > 0:
            tkMB.showinfo('Error','\n\n'.join(errors))
        else:
            pickle.dump((self.zl, self.imagefile, self.dxfile, self.dyfile, self.dz, self.isInf,self.Om0,self.H0),open('last.brt','wb'))
            self.startUp()
            self.root.destroy()

    def loadPrevious(self):
        if os.path.isfile('last.brt'):
            self.zl, self.imagefile, self.dxfile, self.dyfile, self.dz, self.isInf, self.Om0, self.H0 = pickle.load(open('last.brt','rb'))
            self.entry_zl.delete(0,END)
            self.entry_zl.insert(0,str(self.zl))
            self.entry_imagefile.delete(0,END)
            self.entry_imagefile.insert(0,self.imagefile)
            self.entry_dxfile.delete(0,END)
            self.entry_dxfile.insert(0,self.dxfile)
            self.entry_dyfile.delete(0,END)
            self.entry_dyfile.insert(0,self.dyfile)
            self.entry_dz.delete(0,END)
            self.entry_dz.insert(0,str(self.dz))
            self.check.set(self.isInf)
            self.entry_Om0.delete(0,END)
            self.entry_Om0.insert(0,str(self.Om0))
            self.entry_H0.delete(0,END)
            self.entry_H0.insert(0,str(self.H0))
        else:
            tkMB.showinfo('Error','Could not locate previous model. Enter new parameters.')


    def startUp(self):
        global zl, image, deflectx, deflecty, dDXdx, dDXdy, dDYdx, dDYdy, Om0, H0, wcs, scalefactor, xoff, yoff
        zl = self.zl
        Om0 = self.Om0
        H0 = self.H0
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            im = Image.open(self.imagefile)
            wtiff = pickle.load(open(self.imagefile+'.wcs','rb'))

        deflectx = fits.getdata(self.dxfile)
        deflecty = fits.getdata(self.dyfile)
        h = fits.getheader(self.dxfile)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            wcs = WCS(h)  

        ra,dec = wcs.wcs_pix2world([1,h['NAXIS1']+1],[1,h['NAXIS2']+1],1)
        x,y = wtiff.wcs_world2pix(ra,dec,1)

        xoff = 0.5-(x[0]-int(x[0])+x[1]-int(x[1]))/2
        yoff = 0.5-(y[0]-int(y[0])+y[1]-int(y[1]))/2

        image = im.crop((int(x[0]),im.height-int(y[1])+1,int(x[1])+1,im.height-int(y[0])))

        scalefactor = image.width/deflectx.shape[0]
                 
    
        ps = h['CDELT2']*3600
        deflectx /= ps
        deflecty /= ps


        if self.isInf == 0:
            ratio = dlsds(zl,[self.dz],Om0=Om0,H0=H0)
            deflectx /= ratio[0]
            deflecty /= ratio[0]
    
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            dDXdy,dDXdx = np.gradient(deflectx)
            dDYdy,dDYdx = np.gradient(deflecty)

class MainWindow:
    def __init__(self):
        self.root = Tk()
        self.root.wm_title("BRT")
        self.root.geometry("1000x750")
        self.points = []
        self.labels = []
        self.curves = []
        self.arcs = []
        self.arc_annotate = []
        self.arc_labels = np.array([])
        self.arc_x = np.array([])
        self.arc_y = np.array([])
        self.stamps = []
        self.cid = None
        self.X,self.Y = np.meshgrid(np.arange(deflectx.shape[1])*scalefactor,np.arange(deflecty.shape[0])*scalefactor)

        plotFrame = Frame(self.root,width=800,height=750)
        plotFrame.grid(rowspan=10,column=0)

        self.fig = Figure(figsize=(7, 7), dpi=100)
        self.fig.subplots_adjust(left=0.05,right=0.85,top=0.95,bottom=0.05)
        self.ax = self.fig.add_subplot(111)
        self.ax.imshow(image,extent=(-xoff,image.width-xoff,-yoff,image.height-yoff),interpolation='none')
        self.ax.set_xlim(-xoff,image.width-xoff)
        self.ax.set_ylim(-yoff,image.height-yoff)
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)
        self.ax.set_title('cluster z = '+str(zl),size=20)

        canvas = FigureCanvasTkAgg(self.fig, master=plotFrame)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP,fill=BOTH,expand=1)

        self.toolbar = NavigationToolbar2TkAgg(canvas, plotFrame)
        self.toolbar.update()
        canvas._tkcanvas.pack(side=TOP,fill=BOTH,expand=1)

        dropDownFrame = Frame(self.root,padx=10,pady=15)
        dropDownFrame.grid(row=0,column=1,sticky=S)
        self.plot_var = StringVar(dropDownFrame)
        self.plot_var.set("Ray Trace")
        self.connect()
        Label(dropDownFrame,text='Click Mode: ').grid(row=0,column=0,sticky=E)
        option = OptionMenu(dropDownFrame, self.plot_var, "Ray Trace", "Postage Stamp",
            command= lambda x: self.connect() if self.plot_var.get() == 'Ray Trace' else self.connect_postage_stamp())
        option.grid(row=0,column=1)
        Button(dropDownFrame,text='Kill all windows',padx=5,pady=5,command=self.kill_all).grid(row=1,columnspan=2)
        Button(dropDownFrame,text='Clear all points',padx=5,pady=5,command=self.clear_all).grid(row=2,columnspan=2)

        clickFrame = Frame(self.root)
        clickFrame.grid(row=1,column=1)
        
        Label(clickFrame,text='Click to add images',fg='blue').grid(row=0,columnspan=2)
        Label(clickFrame,text = 'source z: ').grid(row=1,column=0,sticky=E)
        self.zsource_entry = Entry(clickFrame,width=10)
        self.zsource_entry.grid(row=1,column=1)
        self.zsource_entry.insert(0,'2.0')
        
        Label(clickFrame,text = 'z step: ').grid(row=2,column=0,sticky=E)
        self.zstep_entry = Entry(clickFrame,width=10)
        self.zstep_entry.grid(row=2,column=1)
        self.zstep_entry.insert(0,'0.5')
        
        Label(clickFrame,text = '# iterations: ').grid(row=3,column=0,sticky=E)
        self.nz_entry = Entry(clickFrame,width=10)
        self.nz_entry.grid(row=3,column=1)
        self.nz_entry.insert(0,'6')

        Button(clickFrame,text="Clear points",padx=5,pady=5,command=self.clear_points).grid(row=4,columnspan=2)
        
        critFrame = Frame(self.root,padx=10)
        critFrame.grid(row=2,column=1)
        
        Label(critFrame,text="Plot critical curve",fg='blue').grid(row=0,columnspan=4)
        Label(critFrame,text='z crit: ').grid(row=1,column=0,sticky=E)
        self.z_crit_entry = Entry(critFrame,width=10)
        self.z_crit_entry.insert(0,'2.0')
        self.z_crit_entry.grid(row=1,column=1)
        
        Button(critFrame,text="Plot",padx=5,pady=5,command=self.plot_crit).grid(row=1,column=2)
        Button(critFrame,text="Clear",padx=5,pady=5,command=self.clear_curve).grid(row=1,column=3)
        
        saveFrame = Frame(self.root,padx=10)
        saveFrame.grid(row=3,column=1)
        Label(saveFrame,text='Save points to region file',fg='blue').grid(row=0,columnspan=3)
        Label(saveFrame,text='file name: ').grid(row=1,column=0,sticky=E)
        self.reg_entry = Entry(saveFrame,width=10)
        self.reg_entry.grid(row=1,column=1)
        self.reg_entry.insert(0,'points.reg')
        
        self.status_var = StringVar(self.root)
        self.status_var.set(' ')
        Label(saveFrame,textvariable=self.status_var).grid(row=2,columnspan=3)
        Button(saveFrame,text='Save',padx=5,pady=5,command=self.plot_ds9).grid(row=1,column=2)

        loadFrame = Frame(self.root,padx=10)
        loadFrame.grid(row=4,column=1)
        Label(loadFrame,text='Load arcs from file',fg='blue').grid(row=0,columnspan=3)
        Label(loadFrame,text='file name: ').grid(row=1,column=0,sticky=E)
        self.file_entry = Entry(loadFrame,width=10)
        self.file_entry.grid(row=1,column=1)
        self.file_entry.insert(0,'arcs.dat')

        Button(loadFrame,text='Load list',padx=5,pady=5,command=self.load_arcs).grid(row=1,column=2)
        Button(loadFrame,text='Clear list',command=self.clear_list,padx=5,pady=5).grid(row=1,column=3)

        scrollFrame = Frame(self.root,padx=10)
        scrollFrame.grid(row=5,column=1)
        Label(scrollFrame,text='arc list').pack()
        scrollbar = Scrollbar(scrollFrame)
        scrollbar.pack(side=RIGHT, fill=Y)

        self.listbox = Listbox(scrollFrame, width=10, height=5, selectmode='multiple',yscrollcommand=scrollbar.set)        
        self.listbox.pack(side=LEFT, fill=BOTH)

        scrollbar.config(command=self.listbox.yview)

        arcFrame = Frame(self.root,padx=10)
        arcFrame.grid(row=6,column=1)
        Button(arcFrame,text='Use selected arc redshift',command=self.insert_redshift,padx=5,pady=5).grid(row=0,columnspan=2)
        Button(arcFrame,text='Plot',command=self.plot_arcs,padx=5,pady=5).grid(row=1,column=0)
        Button(arcFrame,text='Clear',command=self.clear_arcs,padx=5,pady=5).grid(row=1,column=1)
        Button(arcFrame,text='Select all',command=self.select_all,padx=5,pady=5).grid(row=2,column=0)
        Button(arcFrame,text='Deselect all',command=self.deselect_all,padx=5,pady=5).grid(row=2,column=1)        

        self.root.mainloop()

    def connect(self):
        if self.cid is not None:
            self.disconnect()
        self.cid = self.fig.canvas.callbacks.connect('button_press_event', self.on_click)

    def on_click(self,event):
        if event.xdata is None or event.ydata is None: return
        zsource = self.zsource_entry.get()
        zstep = self.zstep_entry.get()
        nz = self.nz_entry.get()
        if self.toolbar._active is not None: return
        for l in self.labels:
            l.remove()
        self.labels = []

        errors = []
        try:
            zsource = float(zsource)
        except ValueError or zsource < zl:
            errors.append('Please enter a number greater than z='+str(zl)+' for source redshift.')
        try:
            zstep = float(zstep)
        except ValueError:
            errors.append('Please enter a number z step.')
        try:
            nz = int(nz)
        except ValueError:
            errors.append('Please enter a number greater than 0 for # iterations.')

        if len(errors) > 0:
            tkMB.showinfo('Error','\n\n'.join(errors))
        else:
            zvalues = zsource + np.arange(nz)*zstep
            ratios = dlsds(zl,zvalues,Om0=Om0,H0=H0)

            colors = ['red','green','magenta','cyan','yellow','blue']
            colorcycle = cycle(colors)

            for i in range(len(ratios)):
                c = colorcycle.next()
                ximp,yimp = predict_images(event.xdata/scalefactor,event.ydata/scalefactor,deflectx,deflecty,dlsds=ratios[i],maxdist=1)
                p, = self.ax.plot(ximp*scalefactor,yimp*scalefactor,'o',color=c,markeredgewidth=0,markersize=3)
                self.points.append(p)
                l = self.fig.text(0.92,0.8-i*0.05,'z='+str(zvalues[i]),color=c,size=14,ha='center',va='center')
                self.labels.append(l)

                self.fig.canvas.draw()

    def connect_postage_stamp(self):
        if self.cid is not None:
            self.disconnect()
        self.cid = self.fig.canvas.callbacks.connect('button_press_event', self.postage_stamp)
        
    def postage_stamp(self,event):
        xcent = event.xdata
        ycent = event.ydata
        if xcent is None or ycent is None: return
        
        stamp = Toplevel(self.root)
        self.stamps.append(stamp)
        fig = Figure(figsize=(2, 2), dpi=100)
        fig.subplots_adjust(left=0.05,right=0.95,top=0.95,bottom=0.05)
        ax = fig.add_subplot(111)
        ax.imshow(image,extent=(-xoff,image.width-xoff,-yoff,image.height-yoff),interpolation='none')
        ax.set_xlim(xcent-200,xcent+200)
        ax.set_ylim(ycent-200,ycent+200)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
    
        canvas = FigureCanvasTkAgg(fig, master=stamp)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP,fill=BOTH,expand=1)
     
        stamp.attributes("-topmost", True)
        
    def kill_all(self):
        for s in self.stamps:
            s.destroy()
        self.stamps = []

    def disconnect(self):
        self.fig.canvas.callbacks.disconnect(self.cid)
        self.cid = None

    def clear_points(self):
        for p in self.points:
            p.remove()
        self.points = []
        for l in self.labels:
            l.remove()
        self.labels = []
        self.fig.canvas.draw()

    def plot_crit(self):
        if self.cid is None: return

        z = self.z_crit_entry.get()
        try:
            z = float(z)
        except ValueError or z < zl:
            errors.append('Please enter a number greater than z='+str(zl)+' for source redshift.')
            return
        
        ratio = dlsds(zl,[z])

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            A = (1-dDXdx*ratio)*(1-dDYdy*ratio)-(dDXdy*ratio)*(dDYdx*ratio)
            mag = 1.0/np.abs(A)
        m = self.ax.contour(self.X,self.Y,mag,levels=[100],colors='white')
        self.curves.append(m)
        self.fig.canvas.draw()

    def clear_curve(self):
        for m in self.curves:
            for coll in m.collections:
                coll.remove()
        self.curves = []
        self.fig.canvas.draw()

    def plot_ds9(self):
        if self.cid is None: return
        f = open(self.reg_entry.get(),'w')
        for p in self.points:
            xy = p.get_xydata()
            
            ra,dec = wcs.wcs_pix2world(xy[:,0]/scalefactor,xy[:,1]/scalefactor,1)
            c = p.get_markerfacecolor()
            for i in range(xy.shape[0]):
                f.write('fk5;circle({0:0.6f},{1:0.6f},0.1") # color={2} width=2\n'.format(ra[i],dec[i],c))
        f.close()
        ts = time.time()
        self.status_var.set('file saved at: '+datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d %H:%M:%S'))
        self.root.update_idletasks()

    def load_arcs(self):
        arcFile = self.file_entry.get()
        if os.path.isfile(arcFile): 
            self.arc_labels,ra,dec,self.arc_z = np.loadtxt(arcFile,usecols=(0,1,2,6),dtype='S8,<f8,<f8,S8',unpack=True)
            self.arc_x, self.arc_y = wcs.wcs_world2pix(ra,dec,1)
            self.arc_x *= scalefactor
            self.arc_y *= scalefactor
            for a,z in zip(self.arc_labels,self.arc_z):
                if float(z) == 0:
                    s = ''
                else:
                    s = ' (z=%s)' % z
                self.listbox.insert(END,a+s)
        else:
            tkMB.showinfo('Error','File "'+arcFile+ '" does not exist.')

    def plot_arcs(self):
        if self.cid is None: return
        for l in self.labels:
            l.remove()
        self.labels = []
        select = self.listbox.curselection()
        if len(select) == 0: return

        errors = []
        a, = self.ax.plot(self.arc_x[[select]],self.arc_y[[select]],'sw')
        self.arcs.append(a)
        for i in range(len(select)):
            x =  self.arc_x[select[i]]
            y = self.arc_y[select[i]]
            t = self.arc_labels[select[i]]
            l = self.ax.text(x,y+50,t,color='white',size=12,ha='center',va='center')
            if x < 0 or y < 0 or x > image.width or y > image.height:
                errors.append(t+' is out of range of deflection matrices.')
            self.arc_annotate.append(l)
        self.fig.canvas.draw()
        if len(errors) > 0:
            tkMB.showinfo('Error','\n\n'.join(errors))

    def clear_arcs(self):
        for a in self.arcs:
            a.remove()
        self.arcs = []
        for l in self.arc_annotate:
            l.remove()
        self.arc_annotate = []
        self.fig.canvas.draw()

    def clear_all(self):
        self.clear_points()
        self.clear_curve()
        self.clear_arcs()

    def select_all(self):
        self.listbox.select_set(0, END)

    def deselect_all(self):
        self.listbox.selection_clear(0, END)

    def clear_list(self):        
        for i in range(self.arc_labels.size):
            self.listbox.delete(0)
        self.arc_labels = np.array([])
        self.arc_x = np.array([])
        self.arc_y = np.array([])

    def insert_redshift(self):
        select = self.listbox.curselection()
        if len(select) == 0: return
        self.zsource_entry.delete(0,END)
        self.zsource_entry.insert(0,self.arc_z[select[0]])
        self.nz_entry.delete(0,END)
        self.nz_entry.insert(0,'1')
        self.z_crit_entry.delete(0,END)
        self.z_crit_entry.insert(0,self.arc_z[select[0]])

if __name__ == '__main__':
    StartWindow()
    MainWindow()

