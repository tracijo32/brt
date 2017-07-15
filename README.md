# brt
Buffalo Ray Trace (BRT) Documentation

About: Buffalo Ray Trace (BRT) is an interactive GUI for plotting image predictions for lens model. BRT is written in Python, utilizing the tkinter GUI library, the matplotlib plotting library, the astropy library of tools for astrophysical data analysis. All are available through Anaconda.

How to get started: The program is run as a python script. Simply type "python brt.py" into the terminal while in the same directory as the file.

---- The start window ---
The start window is GUI used to load the files needed for BRT to run. Enter the following into the boxes:

    Cluster redshift
    Image file: a PNG, TIFF, or JPG of the cluster
    X deflection file: FITS file with the same field of view as the image file
    Y deflection file: FITS file with the same field of view as the image file
    Deflection file redshift (if this is infinity, leave blank and check the 
        box that says this file is set to dls/ds=1)
    Omega M (for a flat universe, Omega M = 1 - Omega L)
    Hubble constant
    
Note that the deflection matrices and the image file do not have to be the same pixel
scale, they only need to be the same field of view. I have a python script helper_functions.py that shows an example of how these files are made.

Click enter to submit the parameter values.

You only need to enter these values the first time you start an application with a given model. The next time you open it, you may click "Load previous model" to load the previous values. (The file is saved as last.brt.)

--- The main window ---
After submitting the parameter values, a new window with your image on the left hand side will appear. This is a matplotlib plot embedded into a canvas in tkinter. The navigation toolbar down at the bottom works just like it would any normal matplotlib plot. Hover your cursor over these so know what they do.

The right panels are where all the action takes place.

Switch the value in the drop down menu under plotting to toggle off all plotting. This is useful if you plan on doing a lot of panning and zooming without plotting.

The clear all button clears everything plotted.

Click to add images
- Indicate the source redshift in the first box.
- The second box is for the step in redshift if you choose to plot more than one redshift.
- The third box is the number of redshifts you want to plot.
- If you chose 2.0,0.5, and 6 for these three boxes, respectively, you would plot 2.0,2.5,3.0,2.5,3.0,3.5.
- Click anywhere on the plot to start plotting values.
- Clear all added points by clicking the clear points button.

Plot critical curve
- Indicate the redshift you want to plot the critical curve for and click the button to plot.
- Clear will only remove the critical curve from the plot.

Save points to region file
- This function saves all plotted points to a region file that can be displayed in ds9. This is helpful in case you want to identify an arc in ds9 without having to flip back and forth between BRT and ds9.
- Clicking the save button will create a file of the name you indicate. A status line should appear showing that a file was created.
- Note: you may need to zoom in quite a bit on ds9 to see these points.

Load arcs from file
- If you have an arcs file that you use with lenstool, you may load those arcs into ds9. This would be helpful to see which arcs you've already identified and how well your model predicts them.
- After loading, a list of your arcs will appear in the box below. Highlight arcs to plot them. There are also convenience buttons for selecting and deselecting all the arcs.
- You can clear all the plotted arcs with clear. This will not clear the predicted images or the critical curve if they are already plotted.

Known bugs or issues:
- If on a Mac, clicking on the icon in the dock will cause the program to suddenly crash. Instead, try clicking the window if you need to revisit it. I'm still trying to figure out why this is seg faulting, and it probably has to do with tkinter rather than something in the code I can fix...Sorry, I know this one is really annoying.
- While a button in the toolbar is toggled on (i.e., while zooming and panning), the window will flash as you move your cursor in and out of the window. It's annoying, but know that the flashing only occurs when you enter and exit the plot frame. Also, it will stop if you toggle off anything in the toolbar.
- All has been tested using Python 2.7. I believe I have covered my bases in case you want to use Python 3+. Similarly, I'm not sure whether or not 
- The GUI can be run over ssh with X11 forwarding on.

FYI what's the deal with the name?
It's a play on Buffalo Trace Kentucky Straight Bourbon Whiskey. It's corny, but delicious. The taste is rich and complex, with hints of vanilla, toffee and candied fruit. The smooth finish lingers on the palate.
