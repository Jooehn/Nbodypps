import matplotlib.pyplot as plt
import numpy as np
import datetime

##### plot setting for python ####
fs0 = 10
fs1 = 13
fs2 = 15
fs3 = 17
lw0 = 0.5
lw1 = 0.8
lw2 = 1.2
lw3 = 2.5
lw4 = 4.2
lw5 = 5

plt.rcParams['font.size']= 16
plt.rcParams['xtick.minor.visible'], plt.rcParams['xtick.top'] = True,True
plt.rcParams['ytick.minor.visible'], plt.rcParams['ytick.right'] = True,True
plt.rcParams['xtick.direction'], plt.rcParams['ytick.direction'] = 'in','in'
plt.rcParams['xtick.labelsize'] = plt.rcParams['ytick.labelsize'] = 14
plt.rcParams['axes.labelsize'] = 18
plt.rcParams['mathtext.fontset'] = 'cm'

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

def add_date(fig,xcoord=0.88,ycoord=0.945):
    """Adds a box with the current date in the upper right corner of
    the figure"""
    date = datetime.datetime.now()
    
    datestr = '${0}$-${1}$-${2}$'.format(date.day,date.month,date.year)
    
    fig.text(xcoord,ycoord,datestr,bbox=dict(facecolor='None'),fontsize=14)