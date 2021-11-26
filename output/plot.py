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
        fig = plt.figure(figsize = (8, 8))

        # 3DAxesを追加
        ax = fig.add_subplot(111, projection='3d')

        # Axesのタイトルを設定
        ax.set_title("", size = 20)

        # 軸ラベルを設定
        ax.set_xlabel("x", size = 14, color = "r")
        ax.set_ylabel("y", size = 14, color = "r")
        ax.set_zlabel("z", size = 14, color = "r")

        # 軸目盛を設定
        # ax.set_xticks([-5.0, -2.5, 0.0, 2.5, 5.0])
        # ax.set_yticks([-5.0, -2.5, 0.0, 2.5, 5.0])

        # データ
        data = np.loadtxt(f,skiprows=1).T
        x = data[2]
        y = data[3]
        z = data[4]
        ax.set_box_aspect((np.ptp(x), np.ptp(y), np.ptp(z)))
        # 曲線を描画
        ax.scatter(x, y, z, s = 3**2, c = data[1], cmap = ListedColormap(['b', 'b', '#00ff0011']))

        # plt.show()
        plt.savefig(i+".png")
        plt.close()