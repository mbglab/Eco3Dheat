'''
Functions used to analyze E. coli 3C data under heat stress
Author: Xu-Ting Wang, revised by Bin-Guang Ma
'''

from pylab import *
from scipy import stats,ndimage
import math
from scipy.stats.stats import spearmanr
import pandas as pd
from scipy.stats import pearsonr
import os
import pymol
from scipy.stats import gaussian_kde

import scipy.stats as st
from sklearn import metrics

# bin divide scheme a
def div_a(fil_nam):
    A = loadtxt(fil_nam)
    B = multiply(ones(A.shape),A)
    B[0,1:] += B[-1,1:]
    B[:,0] +=B[:,-1]
    C = B[:-1,:-1]
    return(C)

# bin divide scheme b
def div_b(fil_nam):
    A = loadtxt(fil_nam)
    B = multiply(ones(A.shape),A)
    B[-2,:-2] += B[-1,:-2]
    B[:,-2] += B[:,-1]
    B[-2,-2] += B[-1,-1]
    C = B[:-1,:-1]
    return(C)

# Convert matrix to list
def mat2lis(mat_nam):
    A = loadtxt(mat_nam)
    n = A.shape[0]
    Lis = []
    for i in range(n):
        for j in range(i,n):
            Lis.append(A[i,j])
    return Lis

# Normalization method scm
def scm(raw_mat):
    A = loadtxt(raw_mat)
    n = A.shape[0]
    B = multiply(ones((n,n)),sum(A,axis=0))
    C = (multiply(ones((n,n)),sum(A,axis=1))).T
    D = divide(A,multiply(B,C))
    D[isnan(D)]=0
    return D

# Normalization method scn
def scn(fil_nam):
    A = loadtxt(fil_nam)
    n = A.shape[0]
    B = multiply(ones((n,n)),linalg.norm(A,axis=0))
    C = (multiply(ones((n,n)),linalg.norm(A,axis=1))).T
    D = divide(A,multiply(B,C))
    D[isnan(D)]=0
    return D

# Interactive frequency calculation
def dto(fil_nam):
    A = loadtxt(fil_nam)
    n = A.shape[0]
    to = 0
    for i in range(n):
        for j in range(i,n):
            to += A[i,j]
    D = divide(A,to)
    D[isnan(D)]=0
    return D

# Short distance interaction frequency calculation
def SIpro(A,nw):
    #A: dto_matrix
    n = A.shape[0]
    sipro = zeros((1,n))
    for i in range(n):
        ipro = 0.0
        for j in range(i-nw,i+nw+1):
            ipro += A[i,j%n]
        sipro[0,i] = ipro/sum(A[i,])
    return sipro
def SIfrq(A,nw):
    n = A.shape[0]
    si = zeros((1,n))
    for i in range(n):
        jl = (i-nw)%n
        jr = (i+nw)%n
        si[0,i] = (A[i,jl]+A[i,jr])/sum(A[i,])
    return si
    
# Ratio matrix of two matrices, n/m
def rationm(n,m):
    ratio = zeros_like(m)
    ratio[logical_and(m == 0, n != 0)]=2
    ratio[logical_and(m == 0, n == 0)]=1
    ratio[m != 0]=n[m!=0]/m[m!=0]
    ratio[ratio > 2] = 2
    ratio[ratio < 0.5] = 0.5
    return ratio

# Rotating ratio matrix
def rota(A):
    A2 = zeros((200,A.shape[0]))
    for i in range(200):
        for j in range(A.shape[0]):
            A2[i,j] = A[j,(j+i)%A.shape[0]]
    return A2

# Transpose the interaction matrix by 45 degree
def mat_tri(raw_mat):
    A = scn(raw_mat)
    B = ndimage.rotate(A,45)
    C = delete(B,arange(int(len(B)/2)+1,len(B),1),axis=0)
    D = where(C,C,nan)
    return(D)    

# Calculate Directionality Index
def cacudi(A, nw): #A: matrix normalized; nw: caculate range of upstream or downstream 
    n = A.shape[0]
    di = np.zeros((n,1))
    for i in range(n):
        L = []
        R = []
        for j in range(i-1,i-nw-1,-1):
            j = j%n
            if A[i,j] <= 0:
                L.append(0)
                continue
            L.append(log(A[i,j]))
        for j in range(i+1,i+nw+1):
            j = j%n
            if A[i,j] <= 0:
                R.append(0)
                continue
            R.append(log(A[i,j]))
        a = sum(L)
        b = sum(R)
        if a != 0 and b != 0:
            di[i] = stats.ttest_rel(R,L)[0]
        else:
            di[i] = 0
    return di

# Find the boundaries of Directionality Index
def dibds(X,x):
    #X: DI list
    n = len(X)
    di_front=[]
    di_bds = []
    for i in range(n):
        try:
            d = X.index(x,i)
            if X[i] == -x and X[i+1:d].count(-x) == 0:
                di_front.append([i,d])
                for j in range(d,i,-1):
                    if X[j]>=0 and X[j-1]<=0:
                        di_bds.append(j)
                        break
        except ValueError:
            pass
    return(di_bds)

# Convert PDB format to XYZ coordinate format
def ReadPDB(PDBfile):
    f_o = open(PDBfile, "r")
    f_r = f_o.read().split("\n")
    f_o.close()   
    XYZ = []
    X = []
    Y = []
    Z = []
    for line in f_r:
        if len(line) < 10:
            continue
        if "ATOM" not in line:
            continue
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        XYZ.append([x, y, z])
    XYZ = array(XYZ)
    return XYZ

# Convert XYZ corrdinate format to pdb format
def WritePDB(positions, pdb_file, ctype = "0"):
    o_file = open(pdb_file, "w")
    o_file.write("\n")
    col1 = "ATOM"
    col3 = "CA MET"
    col8 = "0.20 10.00"
    bin_num = len(positions)
    for i in range(1, bin_num+1):
        col2 = str(i)
        col4 = "B"+col2
        col5 = "%.3f" % positions[i-1][0]
        col6 = "%.3f" % positions[i-1][1]
        col7 = "%.3f" % positions[i-1][2]
        col2 = " "*(5 - len(col2)) + col2
        col4 = col4 + " " * (6 - len(col4))
        col5 = " " * (8 - len(col5)) + col5
        col6 = " " * (8 - len(col6)) + col6
        col7 = " " * (8 - len(col7)) + col7
        col = (col1, col2, col3, col4, col5, col6, col7,col8)
        line = "%s  %s   %s %s   %s%s%s  %s\n" % col
        o_file.write(line)
    col1 = "CONECT"
    for i in range(1, bin_num+1):
        col2 = str(i)
        j = i + 1
        if j > bin_num:
            if ctype == "1":
                continue
            j = 1
        col3 = str(j)
        col2 = " " * (5 - len(col2)) + col2
        col3 = " " * (5 - len(col3)) + col3
        line = "%s%s%s\n" % (col1, col2, col3)
        o_file.write(line)
    o_file.write("END")
    o_file.close()

def StructNor(PDBfile):
    f_o = open(PDBfile)
    f_r = f_o.read().split("\n")
    f_o.close()
    x = []
    y = []
    z = []
    c_r = {}
    for e in f_r:
        if "ATOM" in e:
            x_p = float(e[30:38])
            y_p = float(e[38:46])
            z_p = float(e[46:54])
            x.append(x_p)
            y.append(y_p)
            z.append(z_p)
        if "CONECT" in e:
            c1 = int(e[6:11]) - 1
            c2 = int(e[11:16]) - 1
            c_r[c1] = c2
    chr_len = 0
    x_c = 0
    y_c = 0
    z_c = 0
    for i in range(len(x)):
        j = (i + 1) % len(x)
        chr_len += math.sqrt((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2 + (z[j] - z[i]) ** 2)
        x_c += x[i]
        y_c += y[i]
        z_c += z[i]
    x_c = x_c / len(x)
    y_c = y_c / len(x)
    z_c = z_c / len(x)
    for i in range(len(x)):
        x[i] -= x_c
        y[i] -= y_c
        z[i] -= z_c
    for i in range(len(x)):
        x[i] = x[i] / chr_len * 500
        y[i] = y[i] / chr_len * 500
        z[i] = z[i] / chr_len * 500
    positions = []
    for i in range(len(x)):
        positions.append([x[i], y[i], z[i]])
    WritePDB(positions, PDBfile[:-4]+ "_s.pdb")

# Generate the distance matrix of the model
def GetDistances(XYZ):
    bin_num = XYZ.shape[0]
    D = np.zeros((bin_num, bin_num))
    for i in range(bin_num):
        for j in range(bin_num):
            if i== j:
                continue            
            dis = math.sqrt((XYZ[i, 0] - XYZ[j, 0]) ** 2 + (XYZ[i, 1] - XYZ[j, 1]) ** 2 + (XYZ[i, 2] - XYZ[j, 2]) ** 2)
            D[i, j] = dis
            D[j, i] = dis
    return D

# Generate the distance matrix of the reconstruction circle
def Circdistances(D):
    bin_num = D.shape[0]
    Len_Chr = 0
    for i in range(bin_num):
        j = (i+1)%bin_num
        Len_Chr += D[i,j]
    print('Len of the chromosome is '+str(Len_Chr))
    r = Len_Chr/(2*math.pi)
    XY = []
    bow = 0
    bin_index = 0
    while bow < Len_Chr:
        A = (bow/Len_Chr)*360
        x = r*math.cos(math.radians(A))
        y = r*math.sin(math.radians(A))
        XY.append([x,y])
        bow += D[bin_index,(bin_index+1)%bin_num]
        bin_index += 1
    d = np.zeros((bin_num,bin_num))
    for i in range(bin_num):
        for j in range(bin_num):
            if i == j:
                continue
            dis = math.sqrt((XY[i][0]-XY[j][0])**2+(XY[i][1]-XY[j][1])**2)
            d[i,j]=dis
            d[j,i]=dis
    return d

# Calculate the Global Compactness of chromosome
def CompactGlobal(D, d):
    bin_num = D.shape[0]
    D_all = []
    d_all = []
    for i in range(bin_num):
        for j in range(i+1,bin_num):
            D_all.append(D[i,j])
            d_all.append(d[i,j])
    d_sum = sum(d_all)
    D_sum = sum(D_all)
    CG = -log2(D_sum/d_sum)
    return CG

# Caculate the Local Compactness of chromosome
def cacuCL(D, bin_num, ii, isw):
    D_local = 0
    d_local = 0
    icoor = []
    x = 0
    for j in range(ii,isw):
        icoor.append(x)
        x += D[j%bin_num,(j+1)%bin_num]
        for n in range(j+1,isw):
            D_local += D[j%bin_num,n%bin_num]
    for p in range(len(icoor)):
        for q in range(p+1,len(icoor)):
            d_local += abs(icoor[q]-icoor[p])
    cl = -log2(D_local/d_local)
    return cl

# Caculate the Local Compactness for all bins
def CompactLocal(D, thr):
    bin_num = len(D)
    CL = zeros((int(thr),bin_num))
    for sw in range(1,int(thr)+1):
        C_local = zeros((1,bin_num))
        for i in range(bin_num):
            # D_local = 0
            # d_local = 0
            # icoor = []
            # x = 0
            # for j in range(i-sw,i+sw+1):
                # icoor.append(x)
                # x += D[j%bin_num,(j+1)%bin_num]
                # for n in range(j+1,i+sw+1):
                    # D_local += D[j%bin_num,n%bin_num]
            # for p in range(len(icoor)):
                # for q in range(p+1,len(icoor)):
                    # d_local += abs(icoor[q]-icoor[p])
            # C_local[0,i] = D_local/d_local
            C_local[0,i] = cacuCL(D,bin_num,i-sw,i+sw+1)
        CL[sw-1,] = C_local
    return CL

# Calculate the Local compactness index for all bins
def CompactIndex(D, thr):
    bin_num = D.shape[0]
    CI = zeros((int(thr),bin_num))
    for sw in range(2,int(thr)+1):
        C_index = zeros((1,bin_num))
        for i in range(bin_num):
            il = cacuCL(D,bin_num,i-sw,i+1)
            ir = cacuCL(D,bin_num,i,i+sw+1)
            C_index[0,i] = (ir-il)/(ir+il)
        CI[sw-1,] = C_index
    return CI
    
def smoothListGaussian(list,strippedXs=False,degree=5):
    window=degree*2-1
    weight=np.array([1.0]*window)
    weightGauss=[]
    for i in range(window):
        i=i-degree+1
        frac=i/float(window)
        gauss=1/(np.exp((4*(frac))**2))
        weightGauss.append(gauss)
    weight=np.array(weightGauss)*weight
    smoothed=[0.0]*(len(list)-window)
    for i in range(len(smoothed)):
        smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)
    return smoothed

class DelongTest():
    def __init__(self, preds1, preds2, label, threshold=0.05):
        self._preds1=preds1
        self._preds2=preds2
        self._label=label
        self.threshold=threshold
        self._show_result()
    def _auc(self,X, Y)->float:
        return 1/(len(X)*len(Y)) * sum([self._kernel(x, y) for x in X for y in Y])
    def _kernel(self,X, Y)->float:
        return .5 if Y==X else int(Y < X)
    def _structural_components(self,X, Y)->list:
        V10 = [1/len(Y) * sum([self._kernel(x, y) for y in Y]) for x in X]
        V01 = [1/len(X) * sum([self._kernel(x, y) for x in X]) for y in Y]
        return V10, V01
    def _get_S_entry(self,V_A, V_B, auc_A, auc_B)->float:
        return 1/(len(V_A)-1) * sum([(a-auc_A)*(b-auc_B) for a,b in zip(V_A, V_B)])    
    def _z_score(self,var_A, var_B, covar_AB, auc_A, auc_B):
        return (auc_A - auc_B)/((var_A + var_B - 2*covar_AB )**(.5)+ 1e-8)
    def _group_preds_by_label(self,preds, actual)->list:
        X = [p for (p, a) in zip(preds, actual) if a]
        Y = [p for (p, a) in zip(preds, actual) if not a]
        return X, Y
    def _compute_z_p(self):
        X_A, Y_A = self._group_preds_by_label(self._preds1, self._label)
        X_B, Y_B = self._group_preds_by_label(self._preds2, self._label)
        V_A10, V_A01 = self._structural_components(X_A, Y_A)
        V_B10, V_B01 = self._structural_components(X_B, Y_B)
        auc_A = self._auc(X_A, Y_A)
        auc_B = self._auc(X_B, Y_B)
        var_A = (self._get_S_entry(V_A10, V_A10, auc_A, auc_A) * 1/len(V_A10)+ self._get_S_entry(V_A01, V_A01, auc_A, auc_A) * 1/len(V_A01))
        var_B = (self._get_S_entry(V_B10, V_B10, auc_B, auc_B) * 1/len(V_B10)+ self._get_S_entry(V_B01, V_B01, auc_B, auc_B) * 1/len(V_B01))
        covar_AB = (self._get_S_entry(V_A10, V_B10, auc_A, auc_B) * 1/len(V_A10)+ self._get_S_entry(V_A01, V_B01, auc_A, auc_B) * 1/len(V_A01))
        z = self._z_score(var_A, var_B, covar_AB, auc_A, auc_B)
        p = st.norm.sf(abs(z))*2
        return z,p
    def _show_result(self):
        z,p=self._compute_z_p()
        print(f"z score = {z:.5f};\np value = {p:.5f};")
        if p < self.threshold :print("There is a significant difference")
        else:
            print("There is NO significant difference")