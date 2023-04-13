import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.switch_backend('agg')


h=[8,16,32,64,128,256,512]

for fun in range(1,4):
    for H in h:

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        X = np.arange(0, 1, 1.0/H)
        Y = np.arange(0, 1, 1.0/H)
        X, Y = np.meshgrid(X, Y)

# 打开文件
        with open('../data/fun'+str(fun)+str(H)) as f:
    # 使用loadtxt()函数读取数据
            Z = np.loadtxt(f)
        ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
        ax.contourf(X,Y,Z,zdir='z',offset=-2)
        ax.set_zlim(-2,np.max(Z))
        plt.savefig('../pic/fun'+str(fun)+str(H)+'.eps')
