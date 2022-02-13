import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def griewank_2d_plots():
    from benchmark import griewank_2d
    density = 500
    lo = 0
    hi = 100
    mid = (hi - lo) / 2
    xs = np.linspace(lo, hi, density)
    ys = np.linspace(lo, hi, density)
    zs = np.zeros(shape=(len(ys), len(xs)))
    for i in range(len(ys)):
        y = ys[i]
        for j in range(len(xs)):
            x = xs[j]
            zs[i][j] = griewank_2d(x=x, y=y, x0=mid, y0=mid)
    plt.imshow(zs, cmap='jet', origin='lower')
    ax_ticks = np.arange(start=0, stop=density, step=density/10)
    ax_labels = np.arange(start=lo, stop=hi, step=hi/10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.show()
    # 3D Projection
    x_2d, y_2d = np.meshgrid(xs, ys)
    fig = plt.figure(figsize=(6,5))
    ax = fig.add_subplot(111, projection='3d')
    # Surface Plot
    ax.plot_surface(x_2d, y_2d, zs, cmap='jet')
    plt.show()


def paraboloid_2d_plots():
    from benchmark import paraboloid_2d
    density = 500
    lo = 0
    hi = 10
    mid = (hi - lo) / 2
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
    ax_labels = np.arange(start=lo, stop=hi, step=hi / 10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.title('a. upper view')
    # 3D Projection
    x_2d, y_2d = np.meshgrid(xs, ys)
    ax = fig.add_subplot(gs[0, 1], projection='3d')
    ax.plot_surface(x_2d, y_2d, zs, cmap='Spectral')
    plt.title('b. 3-D view')
    plt.show()


def rastrigin_2d_plots():
    from benchmark import rastrigin_2d
    density = 500
    lo = 0
    hi = 10
    mid = (hi - lo) / 2
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
    ax_labels = np.arange(start=lo, stop=hi, step=hi / 10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.title('a. upper view')
    # 3D Projection
    x_2d, y_2d = np.meshgrid(xs, ys)
    ax = fig.add_subplot(gs[0, 1], projection='3d')
    ax.plot_surface(x_2d, y_2d, zs, cmap='Spectral')
    plt.title('b. 3-D view')
    plt.show()


def himmelblaus_2d_plots():
    from benchmark import himmelblaus
    density = 500
    lo = 0
    hi = 10
    mid = (hi - lo) / 2
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
    ax_labels = np.arange(start=lo, stop=hi, step=hi/10)
    plt.xticks(ticks=ax_ticks, labels=ax_labels)
    plt.yticks(ticks=ax_ticks, labels=ax_labels)
    plt.title('a. upper view')
    # 3D Projection
    x_2d, y_2d = np.meshgrid(xs, ys)
    ax = fig.add_subplot(gs[0, 1], projection='3d')
    ax.plot_surface(x_2d, y_2d, zs, cmap='Spectral')
    plt.title('b. 3-D view')
    plt.show()

paraboloid_2d_plots()