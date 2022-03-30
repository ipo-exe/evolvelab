"""
Recipes or proto-tools

"""

def evolution_2d_recipe():
    from evolution import evolve_2d_function
    from backend import create_rundir, status
    from visuals import convergence, pannel_generation
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import numpy as np
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
    lo_x = 0
    hi_x = 4
    lo_y = 0
    hi_y = 10
    mid = 5.0
    popsize = 100
    genesize = 10
    mutrate = 0.5
    elite = True
    trace = True

    out = evolve_2d_function(lo_x=lo_x,
                             hi_x=hi_x,
                             mid_x=mid,
                             lo_y=lo_y,
                             hi_y=hi_y,
                             mid_y=mid,
                             generations=generations,
                             popsize=popsize,
                             genesize=genesize,
                             mutrate=mutrate,
                             elitism=elite,
                             trace=trace,
                             kind=kind,
                             tui=True)

    curve_df = out['Curve']
    convergence(curve_df=curve_df, folder=folder, show=False)

    if trace:
        lo = 0
        hi = 10
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
        mid = (hi - lo) / 2
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
            # plot
            pannel_generation(trace_df, xs, ys, zs, g, hi, lo, lo_x, hi_x, popsize, folder=folder, show=False)
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