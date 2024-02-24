#select landmark coordinates in python itself instead of separate image viewing software
import SimpleITK as sitk
import os
import matplotlib.pyplot as plt
import argparse
import numpy as np
from tkinter import *

# Flag subject code
parser = argparse.ArgumentParser(description='input subject stl files')
parser.add_argument('-s', '--subject', help='Input subject name (code)', required=True)
args = parser.parse_args()
subject = args.subject

# load in image w/ SimpleITK
path = '/home/msoo935/Code/MRI-to-model_platform_v1/subjects/'
image = sitk.ReadImage(os.path.join(path, subject, 'metaimage', subject + '_lung_raw.mha'))
size = image.GetSize()
print(size)
X = sitk.GetArrayFromImage(image)

coords = [] # initialise array to store clicked coords

def array2exnode(output_path, nodes):
    with open(output_path + '.exnode', 'w') as file:
        file.write(" Group name: \n")
        file.write(" !#nodeset nodes\n")
        file.write(' #Fields=1\n')
        file.write(" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
        file.write("  x.  Value index= 1, #Derivatives=0, #Versions=1\n")
        file.write("  y.  Value index= 2, #Derivatives=0, #Versions=1\n")
        file.write("  z.  Value index= 3, #Derivatives=0, #Versions=1\n")

        # nodes taken in as list. Convert to numpy array
        # nodes = np.asarray(nodes)

        for i in range(len(nodes)):
            file.write(f" Node:            {i+1}\n")
            # line = str(nodes[i]).strip("[]")
            file.write(f"  {nodes[i][0]}\n")
            file.write(f"  {nodes[i][1]}\n")
            file.write(f"  {nodes[i][2]}\n")


def array2ipnode(output_path, nodes):
    with open(output_path + '.ipnode', 'w') as file:
        file.write(' CMISS Version 2.0  ipnode File Version 2\n')
        file.write(' Heading: \n')
        file.write('\n')
        file.write(f' The number of nodes is [   {len(nodes)}]:    {len(nodes)}\n')
        file.write(' Number of coordinates [ 3]: 3\n')
        file.write(' Do you want prompting for different versions of nj=1 [N]? Y\n')
        file.write(' Do you want prompting for different versions of nj=2 [N]? Y\n')
        file.write(' Do you want prompting for different versions of nj=3 [N]? Y\n')
        file.write(' The number of derivatives for coordinate 1 is [0]: 0\n')
        file.write(' The number of derivatives for coordinate 2 is [0]: 0\n')
        file.write(' The number of derivatives for coordinate 3 is [0]: 0\n')
        file.write('\n')
        for i in range(len(nodes)):
            file.write(f' Node number [    {i+1}]:     {i+1}\n')
            file.write(f' The Xj(1) coordinate is [ 0.00000E+00]:    {nodes[i][0]}\n')
            file.write(f' The Xj(2) coordinate is [ 0.00000E+00]:    {nodes[i][1]}\n')
            file.write(f' The Xj(3) coordinate is [ 0.00000E+00]:    {nodes[i][2]}\n')
            file.write('\n')


## Function to collect coords on cursor click.
def onclick(event):
    global ix, iy
    ix, iy = event.xdata, event.ydata
    z = tracker.ind
    # print(f'x = {ix}, y = {iy}, z = {z}') # for axial view
    print(f'x = {ix}, y = {z}, z = {iy}') # for coronal view

    global coords
    # coords.append((ix, iy))
    # coords.append((ix, iy, z)) # for axial view
    coords.append((ix, z, iy)) # for coronal view


    # Change values of selected pixel and its neighbours
    x = int(ix) #+ x_offset
    y = int(iy) #+ y_offset
    z = int(z)
    # X[z][(y-2):(y+2),(x-2):(x+2)] = 1300 # for axial view
    X[y][(z - 2):(z + 2), (x - 2):(x + 2)] = 1300 # for coronal view

    return coords



# Scrollable plot of 2D image slices w/ mpldatacursor
class IndexTracker(object):
    def __init__(self, ax, X):
        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')

        self.X = X
        ## When changing views, rmb to change here and in update()
        rows, self.slices, cols, = X.shape # for coronal view
        # self.slices, rows, cols = X.shape # for axial view
        self.ind = self.slices//2

        # self.im = ax.imshow(self.X[:, :, self.ind]) # for sagittal view
        self.im = ax.imshow(self.X[:, self.ind,:], origin='lower') # for coronal view
        # self.im = ax.imshow(self.X[self.ind, :, :]) # for axial view
        self.update()

    def onscroll(self, event):
        # print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()

    def update(self):
        # self.im.set_data(self.X[:, :, self.ind]) # for sagittal view
        self.im.set_data(self.X[:, self.ind, :]) # for coronal view
        # self.im.set_data(self.X[self.ind, :, :])  # for axial view
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()



# open a figure window
fig, ax = plt.subplots(1, 1)


tracker = IndexTracker(ax, X) # Show plot of 2D slice


fig.canvas.mpl_connect('scroll_event', tracker.onscroll) # connect plot to scrollable slice viewer
fig.canvas.mpl_connect('button_press_event', onclick) # connect plot to coord collector
plt.show()


print("Done")


## Need to convert coords to global positions first
# multiply by -1 bc flipped from RAI to LPI
spacing = image.GetSpacing()
origin = image.GetOrigin()
z_spacing = spacing[2] #1.5
z_offset = origin[2] #-41.221
y_spacing = spacing[1] * -1 #-1.4063
y_offset = origin[1] * -1 #136.797
x_spacing = spacing[0] * -1 #-1.4063
x_offset = origin[0] * -1#176.364

coords = np.asarray(coords) #convert list to np array

for i in range(len(coords)):
    # translate z coord
    coords[i][2] = coords[i][2] * z_spacing + z_offset
    # translate y coord
    coords[i][1] = coords[i][1] * y_spacing + y_offset
    # translate x coord
    coords[i][0] = coords[i][0] * x_spacing + x_offset


## Arrange selected landmarks in order of: trachea nodes, left lung nodes, right lung nodes.
landmarks = np.zeros((11,3))

# Sort along z axis in descending order. Airway entry is the 1st node; as node number increases, z value decreases (growing downwards into the lungs)
coords = coords[coords[:, 2].argsort()] # sorts in ascending order
coords = np.flip(coords, axis=0) # flip array vertically

# 3 trachea landmarks stored. del stored landmarks
landmarks[:3] = coords[:3]
coords = np.delete(coords, slice(0,3), axis=0)


# In Tairawhiti data, global orientation is LPI
# Isolate right lung nodes based on trachea bifurcation node.
right_lung_nodes = []
left_lung_nodes = []
for i in range(len(coords)):
    if coords[i][0] > landmarks[2][0]:
        right_lung_nodes.append(coords[i])
    else:
        left_lung_nodes.append(coords[i])

## convert to np array
right_lung_nodes = np.asarray(right_lung_nodes)
left_lung_nodes = np.asarray(left_lung_nodes)

# Sort along z axis in descending order.
right_lung_nodes = right_lung_nodes[right_lung_nodes[:, 2].argsort()] # sorts in ascending order
right_lung_nodes = np.flip(right_lung_nodes, axis=0) # flip array vertically

# Arrange Nodes 16 (RUL), 11 (right lung hilum), & 17 (lower right lung bifurcation) by z value.
n = 6
for i in range(3):
    landmarks[n] = right_lung_nodes[i]
    n = n+1

right_lung_nodes = np.delete(right_lung_nodes, slice(0,3), axis=0)

## Arrange Node 24 (RLL), smallest y value in right_lung_nodes.
smallest_y = np.min(right_lung_nodes[:,1])
idx = np.where(right_lung_nodes[:,1] == smallest_y)
landmarks[9] = right_lung_nodes[idx]
right_lung_nodes = np.delete(right_lung_nodes, idx, axis=0)

## Assign Node 25 (RML)
landmarks[10] = right_lung_nodes[0]

# Sort along z axis in descending order.
left_lung_nodes = left_lung_nodes[left_lung_nodes[:, 2].argsort()] # sorts in ascending order
left_lung_nodes = np.flip(left_lung_nodes, axis=0) # flip array vertically

# Arrange Nodes 12 (LUL), 8 (left lung hilum), & 13 (LLL) by z value.
n = 3
for i in range(3):
    landmarks[n] = left_lung_nodes[i]
    n = n+1

################################################################################################################
final_centreline = np.zeros((25,3))

final_centreline[0] = landmarks[0] # top trachea
final_centreline[2] = landmarks[1] # middle trachea
final_centreline[4] = landmarks[2] # trachea bifurcation
final_centreline[11] = landmarks[3] # LUL
final_centreline[7] = landmarks[4] # left lung hilum
final_centreline[12] = landmarks[5] # LLL
final_centreline[15] = landmarks[6] # RUL
final_centreline[10] = landmarks[7] # right lung hilum
final_centreline[16] = landmarks[8] # lower right lung bifurcation
final_centreline[23] = landmarks[9] # RLL
final_centreline[24] = landmarks[10] # RML

################################################################################################################
################################################################################################################


## Generate trachea nodes
# find Node 2, the midpoint of Node 1 & Node 3
final_centreline[1] = [(final_centreline[0][0]+final_centreline[2][0])/2.0, (final_centreline[0][1]+final_centreline[2][1])/2.0, (final_centreline[0][2]+final_centreline[2][2])/2.0]

# find Node 4, the midpoint of Node 3 & Node 5
final_centreline[3] = [(final_centreline[2][0]+final_centreline[4][0])/2.0, (final_centreline[2][1]+final_centreline[4][1])/2.0, (final_centreline[2][2]+final_centreline[4][2])/2.0]

print("Done assigning trachea nodes.")



# Generate left lung nodes
## Generate Node 14
x_12_14 = 4.5
y_12_14 = 1.5
z_12_14 = -10
final_centreline[13] = [final_centreline[11][0]-x_12_14, final_centreline[11][1]-y_12_14, final_centreline[11][2]-z_12_14]

## Generate Node 15
x_12_15 = 5.2
y_12_15 = -2.2
z_12_15 = 1.6
final_centreline[14] = [final_centreline[11][0]-x_12_15, final_centreline[11][1]-y_12_15, final_centreline[11][2]-z_12_15]

# Generate Nodes 18 & 19
x_14_18 = 6.4
y_14_18 = 1.5
z_14_18 = -7
final_centreline[17] = [final_centreline[13][0]-x_14_18, final_centreline[13][1]-y_14_18, final_centreline[13][2]-z_14_18]

x_14_19 = -4.7
y_14_19 = -4.2
z_14_19 = -8.4
final_centreline[18] = [final_centreline[13][0]-x_14_19, final_centreline[13][1]-y_14_19, final_centreline[13][2]-z_14_19]



## Generate Nodes 6 & 7
x_interval = (final_centreline[4][0] - final_centreline[7][0])/3.0
y_interval = (final_centreline[4][1] - final_centreline[7][1])/3.0
z_interval = (final_centreline[4][2] - final_centreline[7][2])/3.0

final_centreline[5] = [final_centreline[4][0]-x_interval, final_centreline[4][1]-y_interval, final_centreline[4][2]-z_interval]
final_centreline[6] = [final_centreline[5][0]-x_interval, final_centreline[5][1]-y_interval, final_centreline[5][2]-z_interval]

#Generate Node 20
x_15_20 = 3.6
y_15_20 = -4.0
z_15_20 = 1.7
final_centreline[19] = [final_centreline[14][0]-x_15_20, final_centreline[14][1]-y_15_20, final_centreline[14][2]-z_15_20]

#Generate Node 21
x_13_21 = 2.5
y_13_21 = 3.3
z_13_21 = 7.4
final_centreline[20] = [final_centreline[12][0]-x_13_21, final_centreline[12][1]-y_13_21, final_centreline[12][2]-z_13_21]


print("Done assigning left lung nodes.")



## Generate right lung nodes
## Generate Nodes 9 & 10 in right lung 1st generation
x_interval = (final_centreline[4][0] - final_centreline[10][0])/3.0
y_interval = (final_centreline[4][1] - final_centreline[10][1])/3.0
z_interval = (final_centreline[4][2] - final_centreline[10][2])/3.0

final_centreline[8] = [final_centreline[4][0]-x_interval, final_centreline[4][1]-y_interval, final_centreline[4][2]-z_interval]
final_centreline[9] = [final_centreline[8][0]-x_interval, final_centreline[8][1]-y_interval, final_centreline[8][2]-z_interval]

# Generate Nodes 22 & 23
x_16_22 = -8.5
y_16_22 = -10.3
z_16_22 = -0.45
final_centreline[21] = [final_centreline[15][0]-x_16_22, final_centreline[15][1]-y_16_22, final_centreline[15][2]-z_16_22]

x_16_23 = -1.3
y_16_23 = 1.2
z_16_23 = -7.5
final_centreline[22] = [final_centreline[15][0]-x_16_23, final_centreline[15][1]-y_16_23, final_centreline[15][2]-z_16_23]

print("Done assigning right lung nodes.")

################################################################################################################

## Write exnode file
output_ex_path = os.path.join('/home/msoo935/Code/MRI-to-model_platform_v1/subjects', subject, 'exfile')
if not os.path.exists(output_ex_path):
    os.makedirs(output_ex_path)
array2exnode(os.path.join(output_ex_path, subject+'_upper_airway'), final_centreline)

output_ip_path = os.path.join('/home/msoo935/Code/MRI-to-model_platform_v1/subjects', subject, 'ipfile')
if not os.path.exists(output_ip_path):
    os.makedirs(output_ip_path)
array2ipnode(os.path.join(output_ip_path, subject+'_upper_airway'), final_centreline)