"""
Recipes or proto-tools

"""

def evolution_2d_recipe():
    from evolution import evolve_2d_function
    from backend import create_rundir, status
    import matplotlib.pyplot as plt
    import numpy as np
    import imageio
    import os

    def id_label(id):
        if id < 10:
            return '000' + str(id)
        elif id >= 10 and id < 100:
            return '00' + str(id)
        elif id >= 100 and id < 1000:
            return '0' + str(id)
        elif id >= 1000 and id < 10000:
            return str(id)

    # load parameters
    folder = '/home/ipora/Documents/bin'
    kind = 'himmelblaus' #'paraboloid' #'rastrigin' #'himmelblaus'
    wkpl = True
    label = kind
    # folder setup
    if wkpl:  # if the passed folder is a workplace, create a sub folder
        if label != '':
            label = label + '_'
        folder = create_rundir(label=label + 'EC', wkplc=folder)
    generations = 10
    lo_x = 3
    hi_x = 7
    lo_y = 0
    hi_y = 10
    mid = 5.0
    popsize = 500
    genesize = 12
    mutrate = 0.5
    elite = False
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
    fig = plt.figure(figsize=(7, 5), )  # Width, Height
    plt.plot(curve_df['Gen'], curve_df['Best_S'])
    plt.ylabel('Best Score')
    plt.xlabel('Generation')
    expfile = folder + '/' + 'convergence.png'
    plt.savefig(expfile)
    plt.close(fig)

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
            fig = plt.figure(figsize=(6, 6), )  # Width, Height
            plt.title('Generation: {}'.format(g))
            im = plt.imshow(zs, cmap='Spectral', origin='lower')
            plt.colorbar(im, shrink=0.4)
            ax_ticks = np.arange(start=0, stop=density, step=density / 10)
            ax_labels = np.arange(start=lo, stop=hi, step=hi / 10)
            plt.xticks(ticks=ax_ticks, labels=ax_labels)
            plt.yticks(ticks=ax_ticks, labels=ax_labels)
            plt.ylim(0, len(ys))
            plt.xlim(0, len(xs))
            #plt.plot(mid * density / (hi - lo), mid * density / (hi - lo), 'o', color='magenta', zorder=1)
            plt.scatter(x=trace_df['G{}_x'.format(g)] * density / (hi - lo),
                        y=trace_df['G{}_y'.format(g)] * density / (hi - lo),
                        marker='+', c='k', zorder=2)
            expfile = folder + '/' + 'G' + id_label(g) + '.png'
            #plt.show()
            plt.savefig(expfile)
            plt.close(fig)
        status('plotting gif')
        png_dir = folder
        gifname = png_dir + '/spoiler.gif'
        images = []
        for file_name in sorted(os.listdir(png_dir)):
            if file_name.endswith('.png') and file_name.startswith('G'):
                file_path = os.path.join(png_dir, file_name)
                images.append(imageio.imread(file_path))
        imageio.mimsave(gifname, images)


evolution_2d_recipe()