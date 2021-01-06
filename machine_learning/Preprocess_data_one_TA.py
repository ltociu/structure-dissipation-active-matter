import numpy as np
import math
from tqdm import tqdm
import sys

# The parameter i selects which Pe to pre-process the coordinates for, from the list of desired Pe's.
i = int(sys.argv[1])

L = 25
num_snapshots = 20
nparticles = 312
ndistances = 25
tau = 1

TA = ['0', '3', '6', '9', '12', '15', '18', '21', '24', '27', '30']

# Pe = \sqrt(T_A/tau)

Pe = ['0.00','1.73','2.45', '3.00', '3.46', '3.87', '4.24', '4.58', '4.90', '5.20', '5.48']

def get_ndistances_in_polar_coord(coord, nparticles, ndistances):
    # This function returns, for a given set of particles, the distances from each particles to its neighbors and the angles \theta between them.
    # In other words, for each particle, the vectors pointing from that particle to all its neighbors is transformed to polar coordinates.
    # The output has dimension [2*Nparticles, Nparticles] and lists the distances first, and angles \theta second.
    
    output = np.zeros((nparticles*2, ndistances))
    for i in range(len(coord)):
        output_local = np.zeros((2, nparticles))
        index = 0
        for j in range(len(coord)):
            if j != i:
                if coord[i][0]-coord[j][0] < -L/2:
                    coord_x = coord[i][0]-coord[j][0] + L
                elif coord[i][0]-coord[j][0] > L/2:
                    coord_x = coord[i][0]-coord[j][0] - L
                else:
                    coord_x = coord[i][0]-coord[j][0]

                if coord[i][1]-coord[j][1] < -L/2:
                    coord_y = coord[i][1]-coord[j][1] + L
                elif coord[i][1]-coord[j][1] > L/2:
                    coord_y = coord[i][1]-coord[j][1] - L
                else:
                    coord_y = coord[i][1]-coord[j][1]

                output_local[0][index] = math.sqrt(coord_y**2 + coord_x**2) # These should be the distances
                output_local[1][index] = math.atan2(coord_y, coord_x) # These should be the angles \theta
                index += 1

        output_local_2 = [(output_local[0][i], output_local[1][i]) for i in np.argsort(output_local[0][:])]
        output_local = np.transpose(output_local_2)

        output[i][:] = output_local[0][:ndistances]
        output[nparticles+i][:] = output_local[1][:ndistances]
    return output


data_x = []    
data_y = []
fpy = open('dw.AOU.density0.50.Tau'+str(tau)+'.00.Pe'+Pe[i]+'.Gamma100.00', 'r')
fpcoord = open('positions.AOU.density0.50.Tau'+str(tau)+'.00.Pe'+Pe[i]+'.Gamma100.00', 'r')
linesy = fpy.readlines()
linescoord = fpcoord.readlines()
nsnapshot = 0
for j in range(len(linescoord)):
    if linescoord[j].split('\t')==['Snapshot',str(nsnapshot+1)+'\n']:
        nsnapshot += 1
        index = -1
        coordinates = np.zeros((nparticles,2))
    else:
        index += 1
        coord = linescoord[j].split('\t')
        coordinates[index][0] = float(coord[0])
        coordinates[index][1] = float(coord[1])
        if index == nparticles-1:
            data_x.append(get_ndistances_in_polar_coord(coordinates,nparticles,ndistances))
            data_y.append(float(linesy[0].split('\t')[1]))

np.save('X_TA_' + TA[i]+'.npy', data_x)
np.save('Y_TA_' + TA[i] +'.npy', data_y)
                                  
