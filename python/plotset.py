import matplotlib.pyplot as plt
import numpy as np


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



def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

