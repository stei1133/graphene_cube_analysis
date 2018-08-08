import numpy as np
import matplotlib.pyplot as plt


def readFile(filename):
   data=[]
   headercount=0
   
   global NX,NY,NZ
   with open(filename) as f:
        for line in f:
             line = line.split() # to deal with blank 
             if line:            # lines (ie skip them)
                 try:
                       line = [float(i) for i in line]
                       if len(line)==6:
                            for i in line:
                                data.append(i)
                       elif len(line)==4:
                            if headercount==0:
                                 numofatoms=int(line[0])
                            elif headercount==1:
                                 NX=int(line[0])
                            elif headercount==2:
                                 NY=int(line[0])
                            elif headercount==3:
                                 NZ=int(line[0])
                            headercount+=1
                            
                 except ValueError:
                       total=0
   return data

def xplanaraverage(imported_data):
     averagedpotential=np.zeros((NY,NZ))
     total=0
     for y in range(NY):
          for z in range(NZ):
               for i in np.linspace(y*NZ+z,y*NZ+z+NX*NY*NZ,num=NX,endpoint=False):
                   total+=imported_data[int(i)]
               avg=total/NX
               averagedpotential[y][z]=avg
               total=0      
     return averagedpotential
def makeLineProfile(data,z_index):
     linedata=change_in_potential[0:NY-1][int(NZ/2+.5)]
     plt.plot(linedata)
     plt.xlabel('X-position (Angstroms)')
     plt.ylabel('Electrostatic Potential (eV)')
     tmparrayx=[0, NZ/4,NZ/2,3*NZ/4, NZ]
     labelx=[round(tmparrayx[i]*zstep,1) for i in range(len(tmparrayx))]
     plt.xticks(tmparrayx,labelx)
     plt.show()

xdimension=2.48140559433563
ydimension=35
zdimension=21.3000040156398
count=0
imported_overall_data=[]
imported_overall_data=readFile('2Ang3up3down_potpp.cube')
overallpotential=xplanaraverage(imported_overall_data)
#subtract HF
imported_HF_data=readFile('2Ang3up3down_HF.cube')
HF_potential=xplanaraverage(imported_HF_data)
graph_sub_HFpotential=overallpotential-HF_potential

#subtract graphene
imported_graphene_data=readFile('2Ang3up3down_Cell.cube')
Cell_potential=xplanaraverage(imported_graphene_data)
vacuumdiff=graph_sub_HFpotential[0][0]-Cell_potential[0][0]
change_in_potential=(graph_sub_HFpotential-Cell_potential)*13.6056980659

zstep=zdimension/NZ
ystep=ydimension/NY
makeLineProfile(change_in_potential,NZ/2+.5)

levels=np.linspace(-1,1,num=21)
CS=plt.contourf(change_in_potential,levels,cmap='bwr',origin='lower',extend='both')
plt.xlabel('X-position (Angstroms)')
plt.ylabel('Z-position (Angstroms from the surface)')
tmparray=[0, NY/4,NY/2.0,3*NY/4,NY]
label=[-1*round(NY/2*ystep,2), -1*round(NY/4*ystep,2), 0,
       round(NY/4*ystep,2), round(NY/2*ystep,2)]

plt.yticks(tmparray, label)
tmparrayx=[0, NZ/4,NZ/2,3*NZ/4, NZ]
labelx=[round(tmparrayx[i]*zstep,1) for i in range(len(tmparrayx))]
plt.xticks(tmparrayx,labelx)
cbar=plt.colorbar(CS)
cbar.set_label('Electrostatic Potential (eV)')
plt.show()



