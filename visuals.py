import numpy as np
import matplotlib.pyplot as plt


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

def sphere_2d_plots():
    from benchmark import sphere_2d
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
            zs[i][j] = sphere_2d(x=x, y=y, x0=mid, y0=mid)
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
