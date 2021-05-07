import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import numpy as np
from matplotlib import cm


def surface_plot3d(X, Y, eta, file_name, t=None):
    fig = plt.figure(figsize=(11, 7))
    ax = Axes3D(fig)
    surf = ax.plot_surface(X, Y, eta, rstride=1, cstride=1,
                           cmap=cm.coolwarm, linewidth=0, antialiased=True)

    ax.set_title("Mô tả bề mặt tại $t={:.2f}$ seconds".format(t), fontname="serif", fontsize=17)
    ax.set_xlabel("x [m]", fontname="serif", fontsize=16)
    ax.set_ylabel("y [m]", fontname="serif", fontsize=16)
    ax.set_zlabel("h [m]", fontname="serif", fontsize=16)
    plt.savefig(file_name)

    # plt.show()
    plt.close(fig)


def eta_animation3d(X, Y, eta_list, frame_interval, filename):
    fig = plt.figure(figsize=(8, 8), facecolor="white")
    ax = fig.add_subplot(111, projection='3d')

    surf = ax.plot_surface(X, Y, eta_list[0], cmap=plt.cm.RdBu_r)

    def update_surf(num):
        ax.clear()
        surf = ax.plot_surface(X / 1000, Y / 1000, eta_list[num], cmap=plt.cm.RdBu_r)
        ax.set_title("Mô tả bề mặt H(x,y,t) sau $t={:.2f}$ seconds".format(
            num * frame_interval), fontname="serif", fontsize=19, y=1.04)
        ax.set_xlabel("x [m]", fontname="serif", fontsize=14)
        ax.set_ylabel("y [m]", fontname="serif", fontsize=14)
        ax.set_zlabel("h [m]", fontname="serif", fontsize=16)
        ax.set_xlim(X.min() / 1000, X.max() / 1000)
        ax.set_ylim(Y.min() / 1000, Y.max() / 1000)
        ax.set_zlim(-0.3, 0.7)
        plt.tight_layout()
        return surf,

    anim = animation.FuncAnimation(fig, update_surf,
                                   frames=len(eta_list), interval=10, blit=False)
    mpeg_writer = animation.FFMpegWriter(fps=24, bitrate=10000,
                                         codec="libx264", extra_args=["-pix_fmt", "yuv420p"])
    anim.save("{}.mp4".format(filename), writer=mpeg_writer)
    return anim


if __name__ == '__main__':
    eta_list = []
    x = np.asarray([i for i in range(100)])
    y = np.asarray([i for i in range(100)])

    x, y = np.meshgrid(x, y)
    for i in range(0, 10001, 10):
        eta = []
        file = open(f"result/outputH_{i}.txt")
        lines = file.readlines()[1:]
        for line in lines:
            tmp = [float(x) for x in line.split()]
            eta.append(tmp)

        eta = np.asarray(eta) - 1
        eta_list.append(eta)

        if i in [0, 10, 100, 1000, 10000]:
            surface_plot3d(x, y, eta, f"img/H{i}.png", i)
    eta_animation3d(x, y, eta_list, 10, "h")


    for i in range(0, 10001, 10):
        eta = []
        file = open(f"/Users/macbookpro/Workspace/TTHNC/ShallowWater/result_mpi/outputH_{i}.txt")
        lines = file.readlines()[1:]
        for line in lines:
            tmp = [float(x) for x in line.split()]
            eta.append(tmp)

        eta = np.asarray(eta) - 1
        eta_list.append(eta)

        if i in [0, 10, 100, 1000, 10000]:
            surface_plot3d(x, y, eta, f"img/H_mpi_{i}.png", i)
    eta_animation3d(x, y, eta_list, 10, "h_mpi")