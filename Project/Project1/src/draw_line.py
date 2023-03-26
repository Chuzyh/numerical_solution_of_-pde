import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.switch_backend('agg')


h=[8,16,32,64,128]

for fun in range(1,4):
    with open('../data/fun'+str(fun)+'_regu_points') as f:
    # 使用loadtxt()函数读取数据
        Z = np.loadtxt(f)
    radio=np.log2(Z[0:4,:]/Z[1:5,:])
    #print(radio)
    plt.figure()
    x=np.arange(1,5,1)
    for i in range(0,3):
        yi=radio[:,i]
        plt.plot(x,yi)

    l=["Dirichlet","Neumann","mix"]
    plt.legend(l)
    plt.savefig('../pic/fun'+str(fun)+'_regu_points.eps')

for fun in range(1,4):
    with open('../data/fun'+str(fun)+'_irregu_points') as f:
    # 使用loadtxt()函数读取数据
        Z = np.loadtxt(f)
    radio=np.log2(Z[0:4,:]/Z[1:5,:])
    #print(radio)
    plt.figure()
    x=np.arange(1,5,1)
    for i in range(0,2):
        yi=radio[:,i]
        plt.plot(x,yi)

    l=["Dirichlet","mix"]
    plt.legend(l)
    plt.savefig('../pic/fun'+str(fun)+'_irregu_points.eps')


for fun in range(1,4):
    with open('../data/fun'+str(fun)+'_regu_error_norm') as f:
    # 使用loadtxt()函数读取数据
        Z = np.loadtxt(f)
    radio=np.log2(Z[0:4,:]/Z[1:5,:])
    #print(radio)
    plt.figure()
    x=np.arange(1,5,1)
    for i in range(0,3):
        yi=radio[:,i]
        plt.plot(x,yi)

    l=["1","2","inf"]
    plt.legend(l)
    plt.savefig('../pic/fun'+str(fun)+'_regu_error_norm.eps')

for fun in range(1,4):
    with open('../data/fun'+str(fun)+'_irregu_error_norm') as f:
    # 使用loadtxt()函数读取数据
        Z = np.loadtxt(f)
    radio=np.log2(Z[0:4,:]/Z[1:5,:])
    #print(radio)
    plt.figure()
    x=np.arange(1,5,1)
    for i in range(0,3):
        yi=radio[:,i]
        plt.plot(x,yi)

    l=["1","2","inf"]
    plt.legend(l)
    plt.savefig('../pic/fun'+str(fun)+'_irregu_error_norm.eps')