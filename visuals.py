import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def paraboloid_2d_plots(folder='C:/bin', filename='view_para', suff='', show=True, dark=True):
    from benchmark import paraboloid_2d

    z_cmap = 'Spectral'
    if dark:
        plt.style.use('dark_background')
        marker_color = 'white'
        z_cmap = 'inferno'

    density = 500
    lo = 0
    hi = 10
    mid = (hi + lo) / 2
    xs = np.linspace(lo, hi, density)
    ys = np.linspace(lo, hi, density)
    zs = np.zeros(shape=(len(ys), len(xs)))
    for i in range(len(ys)):
        y = ys[i]
        for j in range(len(xs)):
            x = xs[j]
            zs[i][j] = paraboloid_2d(x=x, y=y, x0=mid, y0=mid)
    fig = plt.figure(figsize=(10, 5), )  # Width, Height
    gs = mpl.gridspec.GridSpec(1, 2, wspace=0.3, hspace=0.45, left=0.05, bottom=0.05, top=0.95, right=0.95)
    ax = fig.add_subplot(gs[0, 0])
    plt.imshow(zs, cmap=z_cmap, origin='lower')
    ax_ticks = np.arange(start=0, stop=density, step=density / 10)
    ax_labels = np.arange(start=lo, stop=hi, step=(hi - lo) / 10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.title('a. upper view')
    # 3D Projection
    x_2d, y_2d = np.meshgrid(xs, ys)
    ax = fig.add_subplot(gs[0, 1], projection='3d')
    ax.plot_surface(x_2d, y_2d, zs, cmap=z_cmap)
    plt.title('b. 3-D view')
    if show:
        plt.show()
        plt.close(fig)
    else:
        # export file
        if suff == '':
            filepath = folder + '/' + filename + '.png'
        else:
            filepath = folder + '/' + filename + '_' + suff + '.png'
        plt.savefig(filepath)
        plt.close(fig)
        return filepath


def rastrigin_2d_plots(folder='C:/bin', filename='view_rast', suff='', show=True, dark=True):
    from benchmark import rastrigin_2d

    z_cmap = 'Spectral'
    if dark:
        plt.style.use('dark_background')
        marker_color = 'white'
        z_cmap = 'inferno'

    density = 500
    lo = 0
    hi = 10
    mid = (hi + lo) / 2
    xs = np.linspace(lo, hi, density)
    ys = np.linspace(lo, hi, density)
    zs = np.zeros(shape=(len(ys), len(xs)))
    for i in range(len(ys)):
        y = ys[i]
        for j in range(len(xs)):
            x = xs[j]
            zs[i][j] = rastrigin_2d(x=x, y=y, x0=mid, y0=mid)
    fig = plt.figure(figsize=(10, 5), )  # Width, Height
    gs = mpl.gridspec.GridSpec(1, 2, wspace=0.3, hspace=0.45, left=0.05, bottom=0.05, top=0.95, right=0.95)
    ax = fig.add_subplot(gs[0, 0])
    plt.imshow(zs, cmap=z_cmap, origin='lower')
    ax_ticks = np.arange(start=0, stop=density, step=density / 10)
    ax_labels = np.arange(start=lo, stop=hi, step=(hi - lo) / 10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.title('a. upper view')
    # 3D Projection
    x_2d, y_2d = np.meshgrid(xs, ys)
    ax = fig.add_subplot(gs[0, 1], projection='3d')
    ax.plot_surface(x_2d, y_2d, zs, cmap=z_cmap)
    plt.title('b. 3-D view')
    if show:
        plt.show()
        plt.close(fig)
    else:
        # export file
        if suff == '':
            filepath = folder + '/' + filename + '.png'
        else:
            filepath = folder + '/' + filename + '_' + suff + '.png'
        plt.savefig(filepath)
        plt.close(fig)
        return filepath


def himmelblaus_2d_plots(folder='C:/bin', filename='view_himm', suff='', show=True, dark=True):
    from benchmark import himmelblaus

    z_cmap = 'Spectral'
    if dark:
        plt.style.use('dark_background')
        marker_color = 'white'
        z_cmap = 'inferno'

    density = 500
    lo = 0
    hi = 10
    mid = (hi + lo) / 2
    xs = np.linspace(lo, hi, density)
    ys = np.linspace(lo, hi, density)
    zs = np.zeros(shape=(len(ys), len(xs)))
    for i in range(len(ys)):
        y = ys[i]
        for j in range(len(xs)):
            x = xs[j]
            zs[i][j] = himmelblaus(x=x, y=y, x0=mid, y0=mid)
    fig = plt.figure(figsize=(10, 5), )  # Width, Height
    gs = mpl.gridspec.GridSpec(1, 2, wspace=0.3, hspace=0.45, left=0.05, bottom=0.05, top=0.95, right=0.95)
    ax = fig.add_subplot(gs[0, 0])
    plt.imshow(zs, cmap=z_cmap, origin='lower')
    ax_ticks = np.arange(start=0, stop=density, step=density/10)
    ax_labels = np.arange(start=lo, stop=hi, step=(hi - lo)/10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.title('a. upper view')
    # 3D Projection
    x_2d, y_2d = np.meshgrid(xs, ys)
    ax = fig.add_subplot(gs[0, 1], projection='3d')
    ax.plot_surface(x_2d, y_2d, zs, cmap=z_cmap)
    plt.title('b. 3-D view')
    if show:
        plt.show()
        plt.close(fig)
    else:
        # export file
        if suff == '':
            filepath = folder + '/' + filename + '.png'
        else:
            filepath = folder + '/' + filename + '_' + suff + '.png'
        plt.savefig(filepath)
        plt.close(fig)
        return filepath


def griewank_2d_plots(folder='C:/bin', filename='view_himm', suff='', show=True, dark=True):
    from benchmark import griewank_2d

    z_cmap = 'Spectral'
    if dark:
        plt.style.use('dark_background')
        marker_color = 'white'
        z_cmap = 'inferno'

    density = 500
    lo = 0
    hi = 10
    mid = (hi + lo) / 2
    xs = np.linspace(lo, hi, density)
    ys = np.linspace(lo, hi, density)
    zs = np.zeros(shape=(len(ys), len(xs)))
    for i in range(len(ys)):
        y = ys[i]
        for j in range(len(xs)):
            x = xs[j]
            zs[i][j] = griewank_2d(x=x, y=y, x0=mid, y0=mid)
    fig = plt.figure(figsize=(10, 5), )  # Width, Height
    gs = mpl.gridspec.GridSpec(1, 2, wspace=0.3, hspace=0.45, left=0.05, bottom=0.05, top=0.95, right=0.95)
    ax = fig.add_subplot(gs[0, 0])
    plt.imshow(zs, cmap=z_cmap, origin='lower')
    ax_ticks = np.arange(start=0, stop=density, step=density/10)
    ax_labels = np.arange(start=lo, stop=hi, step=(hi - lo)/10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.title('a. upper view')
    # 3D Projection
    x_2d, y_2d = np.meshgrid(xs, ys)
    ax = fig.add_subplot(gs[0, 1], projection='3d')
    ax.plot_surface(x_2d, y_2d, zs, cmap=z_cmap)
    plt.title('b. 3-D view')
    if show:
        plt.show()
        plt.close(fig)
    else:
        # export file
        if suff == '':
            filepath = folder + '/' + filename + '.png'
        else:
            filepath = folder + '/' + filename + '_' + suff + '.png'
        plt.savefig(filepath)
        plt.close(fig)
        return filepath


def pannel_nd_scattergram(df_dvars, df_trace, s_ttl,
                          folder='C:/bin',
                          filename='scattergram',
                          suff='',
                          show=True,
                          dark=True):

    marker_color = 'k'
    z_cmap = 'Spectral'
    if dark:
        plt.style.use('dark_background')
        marker_color = 'white'
        z_cmap = 'inferno'

    fig = plt.figure(figsize=(14, 6), )  # Width, Height
    fig.suptitle(s_ttl)

    # smart gridding
    n_cols = len(df_dvars)
    n_rows = 1
    if len(df_dvars) > 4:
        n_cols = 4
        n_div = len(df_dvars) // 4
        n_rem = len(df_dvars) % 4
        if n_rem > 0: ## odd number
            n_rows = n_div + 1
        else:
            n_rows = n_div
    gs = mpl.gridspec.GridSpec(n_rows, n_cols, wspace=0.2, hspace=0.25, left=0.05, right=0.95, bottom=0.1, top=0.9)

    counter = 0
    for i in range(n_rows):
        for j in range(n_cols):
            if counter == len(df_dvars):
                break
            lcl_var = df_dvars['Labels'].values[counter]
            ax = fig.add_subplot(gs[i, j])
            plt.plot(df_trace[lcl_var], df_trace['Score'], '.', c=marker_color, alpha=0.2)
            plt.xlabel(lcl_var)
            plt.xlim((df_dvars['Lo'].values[counter], df_dvars['Hi'].values[counter]))
            counter = counter + 1

    if show:
        plt.show()
        plt.close(fig)
    else:
        # export file
        if suff == '':
            filepath = folder + '/' + filename + '.png'
        else:
            filepath = folder + '/' + filename + '_' + suff + '.png'
        plt.savefig(filepath)
        plt.close(fig)
        return filepath



def pannel_2d_generation(df_trace, xs, ys, zs, g, hi, lo, lo_x, hi_x, popsize,
                         x_lbl='X', y_lbl='Y', folder='C:/bin', suff='', show=True, dark=True):

    marker_color = 'k'
    z_cmap = 'Spectral'
    if dark:
        plt.style.use('dark_background')
        marker_color = 'white'
        z_cmap = 'inferno'

    if df_trace['Score'].min() > 0:
        zmin = 0
    else:
        zmin = df_trace['Score'].min()

    zmax = np.max(zs)

    fig = plt.figure(figsize=(10, 5), )  # Width, Height
    gs = mpl.gridspec.GridSpec(2, 3, wspace=0.3, hspace=0.45, left=0.05, bottom=0.1, top=0.9, right=0.95)
    # plot 1
    ax = fig.add_subplot(gs[:2, :2])
    fig.suptitle('Generation: {}'.format(g))
    im = plt.imshow(zs,
                    cmap=z_cmap,
                    origin='lower',
                    vmin=np.percentile(zs, 5))
    plt.colorbar(im, shrink=0.4)
    ax_ticks = np.arange(start=0, stop=len(xs), step=len(xs) / 10)
    ax_labels = np.arange(start=lo, stop=hi, step=(hi - lo) / 10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.ylim(0, len(ys))
    plt.xlim(0, len(xs))
    plt.xlabel('x')
    plt.ylabel('y')
    aux_a = (len(zs[0]) / (hi - lo))
    aux_b = -aux_a * lo
    plt.scatter(x=(aux_a * df_trace[x_lbl]) + aux_b,
                y=(aux_a * df_trace[y_lbl]) + aux_b,
                marker='+', c='darkgrey', zorder=2)
    # plot 2
    ax = fig.add_subplot(gs[:2, 2])
    plt.scatter(x=df_trace[x_lbl],
                y=df_trace['Score'],
                color=marker_color,
                alpha=0.3,
                edgecolors='none')
    plt.xlim(np.min(xs), np.max(xs))
    plt.ylim(zmin, zmax * 1.1)
    plt.xlabel('x')
    plt.ylabel('Score')
    filename = 'G' + str(g).zfill(6)
    if show:
        plt.show()
        plt.close(fig)
    else:
        # export file
        if suff == '':
            filepath = folder + '/' + filename + '.png'
        else:
            filepath = folder + '/' + filename + '_' + suff + '.png'
        plt.savefig(filepath)
        plt.close(fig)
        return filepath


def convergence(curve_df, folder='C:/bin', filename='convergence', suff='', show=True, dark=True):

    marker_color = 'black'
    color_median = 'tab:red'
    color_range = 'lightsteelblue'
    if dark:
        plt.style.use('dark_background')
        marker_color = 'white'
        color_median = 'cyan'
        color_range = 'silver'

    fig = plt.figure(figsize=(7, 5), )  # Width, Height
    plt.plot(curve_df['Gen'], curve_df['Best_S'], '.', c=marker_color, label='best')
    plt.plot(curve_df['Gen'], curve_df['p50'], c=color_median, label='median')
    plt.fill_between(x=curve_df['Gen'],
                     y1=curve_df['p05'],
                     y2=curve_df['p95'],
                     color=color_range,
                     label='90% range')
    plt.legend(loc='lower right')
    plt.ylabel('Score')
    plt.xlabel('Generation')
    if show:
        plt.show()
        plt.close(fig)
    else:
        # export file
        if suff == '':
            filepath = folder + '/' + filename + '.png'
        else:
            filepath = folder + '/' + filename + '_' + suff + '.png'
        plt.savefig(filepath)
        plt.close(fig)
        return filepath