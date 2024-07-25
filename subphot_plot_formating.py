import matplotlib as mpl

#mpl.use('pgf')

'''
    def figsize(scale):
    fig_width_pt = 504.0/2                        # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size
    '''
pgf_with_latex = {                      # setup matplotlib to use latex for output
    #"pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
    #"font.family": "sans-serif",
    #"font.family": "serif",
    #"text.usetex": True,
    #"font.serif": ['Computer Modern'],                   # blank entries should cause plots to inherit fonts from the document
    #"font.sans-serif": ['Helvetica'],
    #"font.monospace": [],
    "axes.labelsize": 16,               # LaTeX default is 10pt font.
    "font.size": 16,
    "legend.fontsize": 13,               # Make the legend/label fonts a little smaller
    "xtick.labelsize": 15,
    "ytick.labelsize": 15,
    "xtick.top" :True,
    "ytick.right" :True,
    "ytick.direction" : "in",
    "xtick.direction" : "in",
    "ytick.minor.visible" : True,
    "xtick.minor.visible" : True,
    "ytick.major.size": 6,
    "xtick.major.size": 6,
    "ytick.minor.size": 3,
    "xtick.minor.size": 3,
    "ytick.major.width"   : 1,
    "xtick.major.width"   : 1,
    
    #"figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
    #"pgf.preamble": [
    #                 r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
    #                 r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
    #                 ]
}

mpl.rcParams.update(pgf_with_latex)
import matplotlib.pyplot as plt
