# Figure formatting function using seaborn

def formatfigures():
    import seaborn as sns

    sns.set(style='whitegrid', font_scale=1.5, rc={'legend.frameon': True})

    sns.set_style('ticks')

    sns.set_style({"xtick.direction": "in", "ytick.direction": "in"})

    sns.set_context(rc={'lines.markeredgewidth': 0.1})

    mplcolors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728','#9467bd', 
                '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

    sns.set_palette(mplcolors)

    colors = sns.color_palette()

    
    params = {
        # FONT SIZES
        'axes.labelsize': 30,  # Axis Labels
        'axes.titlesize': 30,  # Title
        'font.size': 28,  # Textbox
        'xtick.labelsize': 22,  # Axis tick labels
        'ytick.labelsize': 22,  # Axis tick labels
        'legend.fontsize': 24,  # Legend font size
        'font.family': 'serif',
        'font.fantasy': 'xkcd',
        'font.serif': 'Computer Modern',
        'font.sans-serif': 'Computer Modern Sans serif',
        'font.monospace': 'Courier',

        # AXIS PROPERTIES
        'axes.titlepad': 2*6.0,  # title spacing from axis
        'axes.grid': True,  # grid on plot
        'figure.figsize': (8, 8),  # square plots
        'savefig.bbox': 'tight',  # reduce whitespace in saved figures

        # LEGEND PROPERTIES
        'legend.framealpha': 0.5,
        'legend.fancybox': True,
        'legend.frameon': True,
        'legend.numpoints': 1,
        'legend.scatterpoints': 1,
        'legend.borderpad': 0.1,
        'legend.borderaxespad': 0.1,
        'legend.handletextpad': 0.2,
        'legend.handlelength': 1.0,
        'legend.labelspacing': 0,
    }

    import matplotlib
    matplotlib.rcParams.update(params)
    # print('Function has run')
    from matplotlib import rc
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)


#    print('-------------------------------------------------------------------------')
#    print('\n\n------------ Figure format has been set by formatfigures.py -------------\n\n')
#    print('-------------------------------------------------------------------------')


    return (colors, params)
    

    


