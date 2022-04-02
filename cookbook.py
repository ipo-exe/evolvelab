"""
Recipes or proto-tools

"""

def evolution_2d_recipe():
    from evolution import evolve_2d_function, evolve
    from backend import create_rundir, status
    from visuals import convergence, pannel_2d_generation
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import numpy as np
    import pandas as pd
    import imageio
    import os

    # load parameters
    folder = '/home/ipora/Documents/bin'
    kinds = ['paraboloid', 'rastrigin', 'himmelblaus']
    kind = kinds[1]
    wkpl = True
    label = kind
    # folder setup
    if wkpl:  # if the passed folder is a workplace, create a sub folder
        if label != '':
            label = label + '_'
        folder = create_rundir(label=label + 'EC', wkplc=folder)
    generations = 20
    lo_x = -5
    hi_x = 5
    lo_y = -5
    hi_y = 5
    mid = 0
    popsize = 20
    trace = True
    ranges_df = pd.DataFrame({'Lo': [lo_x, lo_y],
                              'Hi': [hi_x, hi_y],
                              'Labels': ['X', 'Y']})
    out = evolve(genes_df=ranges_df,
                 generations=generations,
                 popsize=popsize,
                 trace=trace,
                 mutt=5,
                 coarse=True)

    curve_df = out['Curve']
    convergence(curve_df=curve_df, folder=folder, show=False)

    if trace:
        lo = lo_x
        hi = hi_x
        if kind == 'rastrigin':
            from benchmark import rastrigin_2d
        elif kind == 'paraboloid':
            from benchmark import paraboloid_2d
        elif kind == 'himmelblaus':
            from benchmark import himmelblaus
        else:
            from benchmark import paraboloid_2d
        trace_df = out['Traced']
        density = 500
        mid = (hi + lo) / 2
        # get background image
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
                    zs[i][j] = himmelblaus(x=lcl_x, y=lcl_y, x0=mid, y0=mid)
                else:
                    zs[i][j] = paraboloid_2d(x=lcl_x, y=lcl_y, x0=mid, y0=mid)

        for g in range(generations):
            # plot images
            status('plot {}'.format(g))
            lcl_df = trace_df.query('Gen == {}'.format(g))
            # plot
            pannel_2d_generation(lcl_df, xs, ys, zs, g, hi, lo, lo_x, hi_x, popsize, folder=folder, show=False)
        status('plotting gif')
        png_dir = folder
        gifname = png_dir + '/animation.gif'
        images = []
        for file_name in sorted(os.listdir(png_dir)):
            if file_name.endswith('.png') and file_name.startswith('G'):
                file_path = os.path.join(png_dir, file_name)
                images.append(imageio.imread(file_path))
        imageio.mimsave(gifname, images)


evolution_2d_recipe()