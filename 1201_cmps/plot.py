# PYTHON_MATPLOTLIB_3D_PLOT_02

# 3次元散布図
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D

for i in glob.glob("*.prof"):
    with open(i,"r") as f:
        # Figureを追加
        print(i)
        fig = plt.figure(figsize = (8, 8))

        # データ
        data = np.loadtxt(f, skiprows=1).T
        x = data[1]
        y = data[2]
        z = data[3]

        if (y is not None) and np.ptp(y)>0 and (z is not None) and np.ptp(z)>0 :
            # z,y=y,z
            ax = fig.add_subplot(111, projection='3d')
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            ax.set_zlabel("z")
            ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
            # ax.view_init(elev=0, azim=-90)
        else:
            ax = fig.add_subplot(111)
            ax.set_xlabel("x")
            ax.set_ylabel("y")
            if y is None or np.ptp(y)==0:
                y,z=z,None
                ax.set_ylabel("z")
            else:
                z=None
            ax.set_box_aspect(np.ptp(y)/np.ptp(x))
        # Axesのタイトルを設定
        # ax.set_title("", size = 20)

        # 軸目盛を設定
        # ax.set_xticks([-5.0, -2.5, 0.0, 2.5, 5.0])
        # ax.set_yticks([-5.0, -2.5, 0.0, 2.5, 5.0])

        # 曲線を描画
        sc = ax.scatter(*([x, y] if z is None else [x, y, z]), s = 3**2, c = (data[0] if len(data)<9 else data[8]))
        fig.colorbar(sc, aspect=40, pad=0.08, orientation='vertical')

        if 0:
            plt.show()
        else:
            plt.savefig(i+".png")
        plt.close()