import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


def paraboloid_2d_plots(folder='C:/bin', filename='convergence', suff='', show=True):
    from benchmark import paraboloid_2d
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
    plt.imshow(zs, cmap='Spectral', origin='lower')
    ax_ticks = np.arange(start=0, stop=density, step=density / 10)
    ax_labels = np.arange(start=lo, stop=hi, step=(hi - lo) / 10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.title('a. upper view')
    # 3D Projection
    x_2d, y_2d = np.meshgrid(xs, ys)
    ax = fig.add_subplot(gs[0, 1], projection='3d')
    ax.plot_surface(x_2d, y_2d, zs, cmap='Spectral')
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


def rastrigin_2d_plots(folder='C:/bin', filename='convergence', suff='', show=True):
    from benchmark import rastrigin_2d
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
    plt.imshow(zs, cmap='Spectral', origin='lower')
    ax_ticks = np.arange(start=0, stop=density, step=density / 10)
    ax_labels = np.arange(start=lo, stop=hi, step=(hi - lo) / 10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.title('a. upper view')
    # 3D Projection
    x_2d, y_2d = np.meshgrid(xs, ys)
    ax = fig.add_subplot(gs[0, 1], projection='3d')
    ax.plot_surface(x_2d, y_2d, zs, cmap='Spectral')
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


def himmelblaus_2d_plots(folder='C:/bin', filename='convergence', suff='', show=True):
    from benchmark import himmelblaus
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
    plt.imshow(zs, cmap='Spectral', origin='lower')
    ax_ticks = np.arange(start=0, stop=density, step=density/10)
    ax_labels = np.arange(start=lo, stop=hi, step=(hi - lo)/10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.title('a. upper view')
    # 3D Projection
    x_2d, y_2d = np.meshgrid(xs, ys)
    ax = fig.add_subplot(gs[0, 1], projection='3d')
    ax.plot_surface(x_2d, y_2d, zs, cmap='Spectral')
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


def pannel_2d_generation(trace_df, xs, ys, zs, g, hi, lo, lo_x, hi_x, popsize,
                         x_lbl='X', y_lbl='Y', folder='C:/bin', suff='', show=True):

    def id_label(id):
        if id < 10:
            return '000' + str(id)
        elif id >= 10 and id < 100:
            return '00' + str(id)
        elif id >= 100 and id < 1000:
            return '0' + str(id)
        elif id >= 1000 and id < 10000:
            return str(id)

    zmin = np.min(zs)
    zmax = np.max(zs)

    fig = plt.figure(figsize=(10, 5), )  # Width, Height
    gs = mpl.gridspec.GridSpec(2, 3, wspace=0.3, hspace=0.45, left=0.05, bottom=0.1, top=0.9, right=0.95)
    # plot 1
    ax = fig.add_subplot(gs[:2, :2])
    fig.suptitle('Generation: {}'.format(g))
    im = plt.imshow(zs, cmap='Spectral', origin='lower')
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
    plt.scatter(x=(aux_a * trace_df[x_lbl]) + aux_b,
                y=(aux_a * trace_df[y_lbl]) + aux_b,
                marker='+', c='k', zorder=2)
    # plot 2
    ax = fig.add_subplot(gs[:2, 2])
    lcl_ax = np.linspace(lo_x, hi_x, popsize)  # np.random.randint(low=-100, high=100, size=popsize)
    plt.scatter(x=trace_df[x_lbl],
                y=trace_df['Score'],
                color='black',
                alpha=0.3,
                edgecolors='none')
    plt.xlim(np.min(xs), np.max(xs))
    plt.ylim(0, zmax * 1.1)
    plt.xlabel('x')
    plt.ylabel('Score')
    filename = 'G' + id_label(g) + '.png'
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


def convergence(curve_df, folder='C:/bin', filename='convergence', suff='', show=True):
    fig = plt.figure(figsize=(7, 5), )  # Width, Height
    plt.plot(curve_df['Gen'], curve_df['Best_S'], 'k.', label='best')
    plt.plot(curve_df['Gen'], curve_df['p50'], 'tab:red', label='median')
    plt.fill_between(x=curve_df['Gen'],
                     y1=curve_df['p05'],
                     y2=curve_df['p95'],
                     color='lightsteelblue',
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