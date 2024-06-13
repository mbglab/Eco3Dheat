'''
Matrix generate
python Matrix_Generate.py <filtered reads> <bin_size>
Author: Xu-Ting Wang
'''

from pylab import *
import sys

BIN = int(sys.argv[2])

mat={}
maxb = 0

i = 0
with open(sys.argv[1]) as f:
    for line in f:
        i += 1
        if i%1000000 ==0:
            print(i)
        chr1,locus1,sens1,indice1, chr2,locus2,sens2,indice2 = line.split()
        locus1=int(locus1); sens1=int(sens1)
        locus2=int(locus2); sens2=int(sens2)
        bin1 = int(locus1/BIN)
        bin2 = int(locus2/BIN)
        key1=(bin1,bin2)
        key2=(bin2,bin1)
        
        if key1 in mat:
            mat[key1] += 1
        else:
            mat[key1] = 1
        if key2 in mat:
            mat[key2] += 1
        else:
            mat[key2] = 1 
        if bin1 > maxb:
            maxb = bin1
        if bin2 > maxb:
            maxb = bin2

N_BINS = maxb+1
print(N_BINS)

matr = zeros((N_BINS,N_BINS))
for i in range(0,N_BINS):
    for j in range(0,N_BINS):
        key2 = (i,j)
        if key2 in mat:
            matr[i,j] = mat[key2]
            matr[j,i] = mat[key2]
savetxt(sys.argv[1].split('/')[1].split('.')[0]+'_'+str(int(int(sys.argv[2])/1000))+'kb.matr',matr)
print(sys.argv[1]+' matrix generated!')