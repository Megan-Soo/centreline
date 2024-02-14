#select landmark coordinates in python itself instead of separate image viewing software
import SimpleITK as sitk
import os
import matplotlib.pyplot as plt
import argparse
import numpy as np

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
    print(f'x = {ix}, y = {iy}, z = {z}')

    global coords
    # coords.append((ix, iy))
    coords.append((ix, iy, z))

    # if len(coords) == 2:
    #     fig.canvas.mpl_disconnect(cid)

    # rect = (ix,iy,2,2)
    # fig.canvas.drawRectangle(rect)

    # plt.gca().add_patch(Rectangle((ix-1,iy-1), 2, 2, edgecolor='red', facecolor='none', lw=4))
    x = int(ix) #+ x_offset
    y = int(iy) #+ y_offset
    X[z][(y-2):(y+2),(x-2):(x+2)] = 1300 # max value in the image stack is 1301

    return coords

# cid = fig.canvas.mpl_connect('button_press_event', onclick)



# Scrollable plot of 2D image slices w/ mpldatacursor
class IndexTracker(object):
    def __init__(self, ax, X):
        self.ax = ax
        ax.set_title('use scroll wheel to navigate images')

        self.X = X
        # rows, cols, self.slices= X.shape
        self.slices, rows, cols = X.shape
        self.ind = self.slices//2

        # self.im = ax.imshow(self.X[:, :, self.ind]) # for sagittal view
        self.im = ax.imshow(self.X[self.ind, :, :]) # for axial view
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
        self.im.set_data(self.X[self.ind, :, :])  # for axial view
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()



# open a figure window
fig, ax = plt.subplots(1, 1)


tracker = IndexTracker(ax, X) # Show plot of 2D slice


fig.canvas.mpl_connect('scroll_event', tracker.onscroll) # connect plot to scrollable slice viewer
fig.canvas.mpl_connect('button_press_event', onclick) # connect plot to coord collector
plt.show()

print("Done")


# remove coords w/ z values beyond slice range (happens due to scrolling past end-slices)
# or just don't scroll beyond slice range



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


## Sort coords according to template node numbering
landmarks = np.zeros((6,3))

# Sort along z axis in descending order. Airway entry is the 1st node; as node number increases, z value decreases (growing downwards into the lungs)
coords = coords[coords[:, 2].argsort()] # sorts in ascending order
coords = np.flip(coords, axis=0) # flip array vertically

# 3 trachea landmarks stored. del stored landmarks
landmarks[:3] = coords[:3]
coords = np.delete(coords, slice(0,3), axis=0)


# In Tairawhiti data, global orientation is LPI
## Arrange Node 11 ie the only right lung node.
largest_x = np.max(coords[:,0]) # find point w/ largest x val after bifurcation node (Node 5)
idx = np.where(coords[:,0]==largest_x) # idx is a tuple: idx = ([1],)
landmarks[5] = coords[idx]
coords = np.delete(coords, idx, axis=0)

# Arrange Nodes 8 & 12 in the left lung
smallest_x = np.min(coords[:,0])
idx = np.where(coords[:,0]==smallest_x)
landmarks[4] = coords[idx]
coords = np.delete(coords, idx, axis=0)
landmarks[3] = coords

################################################################################################################
################################################################################################################
################################################################################################################
final_centreline = np.zeros((25,3))

## Create a node between selected trachea points
# find Node 2, the midpoint of Node 1 & Node 3
final_centreline[1] = [(landmarks[0][0]+landmarks[1][0])/2.0, (landmarks[0][1]+landmarks[1][1])/2.0, (landmarks[0][2]+landmarks[1][2])/2.0]

# find Node 4, the midpoint of Node 3 & Node 5
final_centreline[3] = [(landmarks[1][0]+landmarks[2][0])/2.0, (landmarks[1][1]+landmarks[2][1])/2.0, (landmarks[1][2]+landmarks[2][2])/2.0]

final_centreline[0] = landmarks[0]
final_centreline[2] = landmarks[1]
final_centreline[4] = landmarks[2]
print("Done assigning trachea nodes.")

## Relative distance of Node 13 from Node 8 (Left Lung)
x_8_13 = 6.0 # node 8 coords minus node 13 coords
y_8_13 = 0.0
z_8_13 = 25.0
final_centreline[12] = [landmarks[3][0]-x_8_13, landmarks[3][1]-y_8_13, landmarks[3][2]-z_8_13]

## Relative distance of Node 14 from Node 12 (Left Lung)
x_12_14 = 4.5
y_12_14 = 1.5
z_12_14 = -10
final_centreline[13] = [landmarks[4][0]-x_12_14, landmarks[4][1]-y_12_14, landmarks[4][2]-z_12_14]

# Assign nodes 18 & 19
x_14_18 = -4.618196097785565e+01 - (-5.258204112862408e+01)
y_14_18 = -5.329925357443111e+01 - (-5.483413465178744e+01)
z_14_18 = 1.102828459712099e+02 - (1.172760683630398e+02)
final_centreline[17] = [final_centreline[13][0]-x_14_18, final_centreline[13][1]-y_14_18, final_centreline[13][2]-z_14_18]

x_14_19 = -4.618196097785565e+01 - (-4.150666826874128e+01)
y_14_19 = -5.329925357443111e+01 - (-4.913028954352593e+01)
z_14_19 = 1.102828459712099e+02 - (1.186481923270559e+02)
final_centreline[18] = [final_centreline[13][0]-x_14_19, final_centreline[13][1]-y_14_19, final_centreline[13][2]-z_14_19]

## Relative distance of Node 15 from Node 12 (Left Lung)
x_12_15 = -4.598829210581377e+01 - (-5.113998437694327e+01)
y_12_15 = -5.122125033167512e+01 - (-4.903304273674225e+01)
z_12_15 = -6.408140371586414e+00 - (-8.019036835019175e+00)
final_centreline[14] = [landmarks[4][0]-x_12_15, landmarks[4][1]-y_12_15, landmarks[4][2]-z_12_15]

final_centreline[7] = landmarks[3]
final_centreline[11] = landmarks[4]

## Assign Nodes 6 & 7 in 1st generation
x_interval = (final_centreline[4][0] - final_centreline[7][0])/3.0
y_interval = (final_centreline[4][1] - final_centreline[7][1])/3.0
z_interval = (final_centreline[4][2] - final_centreline[7][2])/3.0

final_centreline[5] = [final_centreline[4][0]-x_interval, final_centreline[4][1]-y_interval, final_centreline[4][2]-z_interval]
final_centreline[6] = [final_centreline[5][0]-x_interval, final_centreline[5][1]-y_interval, final_centreline[5][2]-z_interval]

#Assign Node 20
x_15_20 = -5.613998437694327e+01 - (-5.976084630276873e+01)
y_15_20 = -4.903304273674225e+01 - (-4.505820059438779e+01)
z_15_20 = 9.198096316498082e+01 - (9.023212758569248e+01)
final_centreline[19] = [final_centreline[14][0]-x_15_20, final_centreline[14][1]-y_15_20, final_centreline[14][2]-z_15_20]

#Assign Node 21
x_13_21 = -4.513724479993320e+01 - (-4.760721332561324e+01)
y_13_21 = -6.446167035559129e+01 - (-6.773530801912059e+01)
z_13_21 = 8.290107704304990e+01 - (7.554822188005851e+01)
final_centreline[20] = [final_centreline[12][0]-x_13_21, final_centreline[12][1]-y_13_21, final_centreline[12][2]-z_13_21]


print("Done assigning left lung nodes.")

## Relative distance of Node 16 from Node 11 (Right Lung)
x_11_16 = 2.495796935409081e+01 - (4.550830645417334e+01)
y_11_16 = -5.696376219466988e+01 - (-5.281966131862813e+01)
z_11_16 = -10
final_centreline[15] = [landmarks[5][0]-x_11_16, landmarks[5][1]-y_11_16, landmarks[5][2]-z_11_16]

## Relative distance of Node 17 from Node 11 (Right Lung)
x_11_17 = 2.495796935409081e+01 - (3.926763733386907e+01)
y_11_17 = -5.696376219466988e+01 - (-5.370586795628331e+01)
z_11_17 = 1.225307303958735e+01 - (0.1e+01) #(-2.054374806372181e+01)
final_centreline[16] = [landmarks[5][0]-x_11_17, landmarks[5][1]-y_11_17, landmarks[5][2]-z_11_17]

final_centreline[10] = landmarks[5]

## Assign Nodes 9 & 10 in right lung 1st generation
x_interval = (final_centreline[4][0] - final_centreline[10][0])/3.0
y_interval = (final_centreline[4][1] - final_centreline[10][1])/3.0
z_interval = (final_centreline[4][2] - final_centreline[10][2])/3.0

final_centreline[8] = [final_centreline[4][0]-x_interval, final_centreline[4][1]-y_interval, final_centreline[4][2]-z_interval]
final_centreline[9] = [final_centreline[8][0]-x_interval, final_centreline[8][1]-y_interval, final_centreline[8][2]-z_interval]

# Assign Nodes 22 & 23
x_16_22 = 3.534039360879932e+01 - (4.380813391713374e+01)
y_16_22 = -5.510146657023744e+01 - (-4.479868485609614e+01)
z_16_22 = 1.127523540266833e+02 - (1.132094549030188e+02)
final_centreline[21] = [final_centreline[15][0]-x_16_22, final_centreline[15][1]-y_16_22, final_centreline[15][2]-z_16_22]

x_16_23 = 3.534039360879932e+01 - (3.666848246997668e+01)
y_16_23 = -5.510146657023744e+01 - (-5.634325194272495e+01)
z_16_23 = 1.127523540266833e+02 - (1.202465628182297e+02)
final_centreline[22] = [final_centreline[15][0]-x_16_23, final_centreline[15][1]-y_16_23, final_centreline[15][2]-z_16_23]

# Assign Nodes 24 & 25
x_17_24 = 3.053089452003725e+01 - (3.797513892108935e+01)
y_17_24 = -5.411148920595319e+01 - (-6.235725610016527e+01)
z_17_24 = 9.140545396347598e+01 - (7.238945904848460e+01)
final_centreline[23] = [final_centreline[16][0]-x_17_24, final_centreline[16][1]-y_17_24, final_centreline[16][2]-z_17_24]

x_17_25 = 3.053089452003725e+01 - (4.738006890717470e+01)
y_17_25 = -5.411148920595319e+01 - (-3.504017256739004e+01)
z_17_25 = 12
final_centreline[24] = [final_centreline[16][0]-x_17_25, final_centreline[16][1]-y_17_25, final_centreline[16][2]-z_17_25]

print("Done assigning right lung nodes.")

################################################################################################################

## Write exnode file
output_ex_path = os.path.join('/home/msoo935/Code/MRI-to-model_platform_v1/subjects', subject, 'exfile',subject+'_upper_airway')
array2exnode(output_ex_path, final_centreline)
output_ip_path = os.path.join('/home/msoo935/Code/MRI-to-model_platform_v1/subjects', subject, 'ipfile',subject+'_upper_airway')
array2ipnode(output_ip_path, final_centreline)