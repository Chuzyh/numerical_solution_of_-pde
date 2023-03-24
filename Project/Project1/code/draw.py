import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.switch_backend('agg')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
h=16
X = np.arange(0, 1, 1.0/h)
Y = np.arange(0, 1, 1.0/h)
X, Y = np.meshgrid(X, Y)

# 打开文件
with open('result') as f:
    # 使用loadtxt()函数读取数据
    Z = np.loadtxt(f)
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
ax.contourf(X,Y,Z,zdir='z',offset=-2)
ax.set_zlim(-2,2)
plt.savefig('../pic/A.eps')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X = np.arange(0, 1, 1.0/h)
Y = np.arange(0, 1, 1.0/h)
X, Y = np.meshgrid(X, Y)
Z=np.exp(X+np.sin(Y))
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=plt.get_cmap('rainbow'))
ax.contourf(X,Y,Z,zdir='z',offset=-2)
ax.set_zlim(-2,2)
plt.savefig('../pic/B.eps')