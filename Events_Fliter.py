'''
3C library events analysis
python Events_Fliter.py <input_file> 
Author: Xu-ting Wang
'''

import sys
from pylab import *

i = 0
fout = open(sys.argv[1]+'.filtered','w')
nb_weirds = 0
nb_loops = 0
nb_uncuts = 0
lrange_intra = 0
thr_loop = 1
thr_uncut = 1
with open(sys.argv[1]) as f:#input_file: output_alignment_idpt.dat.indices from 3C_toturiol
    for line in f:
        i = i+1
        if i%1000000 ==0:
            print(i)
        chr1, locus1, sens1, indice1, chr2, locus2, sens2, indice2 = line.split()
        locus1=int(locus1); sens1=int(sens1); indice1=int(indice1)
        locus2=int(locus2); sens2=int(sens2); indice2=int(indice2)
        if locus2 < locus1:
            sens=sens1;sens1=sens2;sens2=sens;
            locus=locus1;locus1=locus2;locus2=locus;
            indice=indice1;indice1=indice2;indice2=indice
        nsites = indice2-indice1
        if indice1 == indice1 and ((sens1==0 and sens2==0) or (sens1==16 and sens2==16)):
            nb_weirds+=1
        elif nsites <= thr_loop and (sens1==16 and sens2==0):
            nb_loops+=1
        elif nsites <= thr_uncut and (sens1==0 and sens2==16):
            nb_uncuts+=1
        else:
            lrange_intra+=1
            fout.write(str(chr1)+"\t"+str(locus1)+"\t"+str(sens1)+"\t"+str(indice1)+"\t"+str(chr2)+"\t"+str(locus2)+"\t"+str(sens2)+"\t"+str(indice2)+"\n")
fout.close()
labels = ['weirds','loops','uncuts','valids']
fracs = [nb_weirds,nb_loops,nb_uncuts,lrange_intra]
colors = ['gold', 'lightskyblue', 'lightcoral', 'darkblue']
pie(fracs , labels=labels, colors=colors, autopct='%1.1f%%', shadow=True, startangle=90)
title('Distribution of different events in the 3C library')
plt.text(-1.5, -1.2, "Total number of reads ="+str(i))
plt.text(-1.5, -1.3, "Ratio uncuts ="+str((nb_uncuts/i)*100)+"%")
plt.text(-1.5, -1.4, "Ratio loops ="+str((nb_loops/i)*100)+"%")
plt.text(-1.5, -1.5, "Ratio valids ="+str((lrange_intra/i)*100)+"%")
savefig(sys.argv[1]+'.events.png')
fe = open(sys.argv[1]+'.events','w')
print('Total number of reads: ',str(i),file=fe)
print('Uncuts number: ',str(nb_uncuts),file=fe)
print('Loops number: ',str(nb_loops),file=fe)
print('Weirds number: ',str(nb_weirds),file=fe)
print('Valids number: ',str(lrange_intra),file=fe)
fe.close()
print('done!')        