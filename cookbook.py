"""
Recipes or proto-tools

"""
import os


def view_function():
    import visuals

    folder = '/home/ipora/Documents/bin'
    visuals.griewank_2d_plots(folder='/home/ipora/PycharmProjects/evolve3/docs', filename='grie', show=False)
    #visuals.himmelblaus_2d_plots(folder='/home/ipora/PycharmProjects/evolve3/docs', filename='himm', show=False)
    #visuals.rastrigin_2d_plots(folder='/home/ipora/PycharmProjects/evolve3/docs', filename='rastr_2d', show=False)
    #visuals.paraboloid_2d_plots(folder='/home/ipora/PycharmProjects/evolve3/docs', filename='parab', show=False)

def plot_2d_evolution():


def evolution_2d_recipe():
    """
    Recipe for 2d evolution
    :return:
    """
    from evolution import evolve_2d_function, evolve
    from backend import create_rundir, status
    from visuals import convergence, pannel_2d_generation
    from out import export_gif
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import numpy as np
    import pandas as pd
    # load parameters
    folder = '/home/ipora/Documents/bin'
    #folder = 'C:/bin'
    kinds = ['paraboloid', 'rastrigin', 'himmelblaus', 'griewank']
    kind = kinds[2]
    wkpl = True
    label = kind
    # folder setup
    if wkpl:  # if the passed folder is a workplace, create a sub folder
        if label != '':
            label = label + '_'
        folder = create_rundir(label=label + 'EC', wkplc=folder)
    # parameters setup
    generations = 20
    lo_x = -5
    hi_x = 5
    lo_y = -5
    hi_y = 5
    mid = 0
    popsize = 100
    r_std = 0.8
    trace = True
    # definir df
    ranges_df = pd.DataFrame({'Lo': [lo_x, lo_y],
                              'Hi': [hi_x, hi_y],
                              'Labels': ['X', 'Y']})
    # run evolution
    out = evolve(df_genes=ranges_df,
                 n_generations=generations,
                 n_popsize=popsize,
                 b_trace=trace,
                 r_std=r_std,
                 kind = kind,
                 b_coarse=True,
                 b_recomb=False,
                 b_explore=True,
                 upper=80,
                 lower=40)
    # retrieve curve dataframe
    curve_df = out['Curve']
    # plot convergence
    convergence(curve_df=curve_df, folder=folder, show=False)
    # export curve csv:
    fpath = folder + '/convergence.txt'
    curve_df.to_csv(fpath, sep=';', index=False)
    # retrieve traced dataframe
    trace_df = out['Traced']
    # export trace csv:
    fpath = folder + '/evolution.txt'
    trace_df.to_csv(fpath, sep=';', index=False)

    # plotting
    if trace:
        lo = lo_x
        hi = hi_x
        if kind == 'rastrigin':
            from benchmark import rastrigin_2d
        elif kind == 'paraboloid':
            from benchmark import paraboloid_2d
        elif kind == 'himmelblaus':
            from benchmark import himmelblaus
        elif kind == 'griewank':
            from benchmark import griewank_2d
        else:
            from benchmark import paraboloid_2d
        #
        # define density
        density = 500
        mid = (hi + lo) / 2

        # generate background image
        xs = np.linspace(lo, hi, density)
        ys = np.linspace(lo, hi, density)
        zs = np.zeros(shape=(len(ys), len(xs)))
        for i in range(len(ys)):
            lcl_y = ys[i]
            for j in range(len(xs)):
                lcl_x = xs[j]
                # -- compute local score
                if kind == 'rastrigin':
                    zs[i][j] = rastrigin_2d(x=lcl_x, y=lcl_y, x0=mid, y0=mid)
                elif kind == 'paraboloid':
                    zs[i][j] = paraboloid_2d(x=lcl_x, y=lcl_y, x0=mid, y0=mid)
                elif kind == 'himmelblaus':
                    zs[i][j] = himmelblaus(x=lcl_x, y=lcl_y, x0=mid, y0=mid, level=100)
                elif kind == 'griewank':
                    zs[i][j] = griewank_2d(x=lcl_x, y=lcl_y, x0=mid, y0=mid, level=100)
                else:
                    zs[i][j] = paraboloid_2d(x=lcl_x, y=lcl_y, x0=mid, y0=mid)

        # plot generations
        for g in range(generations):
            # plot images
            status('plot {}'.format(g))
            lcl_df = trace_df.query('Gen == {}'.format(g))
            # plot
            pannel_2d_generation(lcl_df, xs, ys, zs, g, hi, lo, lo_x, hi_x, popsize, folder=folder, show=False)
        status('plotting gif')
        export_gif(dir_output=folder, dir_images=folder, nm_gif='animation', kind='png', suf='G')
        # auto open folder in Windows
        #os.startfile(folder)

evolution_2d_recipe()

#view_function()