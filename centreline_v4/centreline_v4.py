# centreline v4 takes 12 nodes selected from raw airway image, assembles a 14-node centreline template of the upper airway, and writes out the exnode ipnode files.
import SimpleITK as sitk
import os
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import argparse
import numpy as np
from matplotlib.widgets import Slider


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


# Scrollable plot of 2D image slices w/ mpldatacursor
class IndexTracker(object):
    def __init__(self, ax, X):
        self.ax = ax
        ax.set_title('Use scroll wheel to navigate images', loc='center', wrap=True)

        self.X = X
        if view == 's':
            rows, cols, self.slices = X.shape  # for sagittal view
        elif view == 'c':
            rows, self.slices, cols, = X.shape  # for coronal view
        elif view == 'a':
            self.slices, rows, cols = X.shape  # for axial view
        else:  # default to coronal view
            rows, self.slices, cols, = X.shape
        self.ind = self.slices // 2

        if view == 's':
            self.im = ax.imshow(self.X[:, :, self.ind])  # for sagittal view
        elif view == 'c':
            self.im = ax.imshow(self.X[:, self.ind, :])  # for coronal view
        elif view == 'a':
            self.im = ax.imshow(self.X[self.ind, :, :])  # for axial view
        else:  # default to coronal view
            self.im = ax.imshow(self.X[:, self.ind, :])
            # Tairawhiti MRI orientations incorrect. Manually flip for now. With imgs w/ correct header oreintations, use direction to decide if flip (in next version).
        self.update()

    def onscroll(self, event):
        # print("%s %s" % (event.button, event.step))
        if event.button == 'up':
            self.ind = (self.ind + 1) % self.slices
        else:
            self.ind = (self.ind - 1) % self.slices
        self.update()
        amp_slider.set_val(self.ind)

    def update(self):
        if view == 's':
            self.im.set_data(self.X[:, :, self.ind])  # for sagittal view
        elif view == 'c':
            self.im.set_data(self.X[:, self.ind, :])  # for coronal view
        elif view == 'a':
            self.im.set_data(self.X[self.ind, :, :])  # for axial view
        else:  # default to coronal view
            self.im.set_data(self.X[:, self.ind, :])
        self.ax.set_ylabel('slice %s' % self.ind)
        self.im.axes.figure.canvas.draw()


    ## Function to collect coords on cursor click.
    def onclick(self, event):
        global ix, iy
        ix, iy = event.xdata, event.ydata
        z = tracker.ind

        xmin = 0
        ymin = 0
        xmax = ax.dataLim.bounds[2]
        ymax = ax.dataLim.bounds[3]

        try:
            # Change values of selected pixel and its neighbours
            x = int(ix)  # + x_offset
            y = int(iy)  # + y_offset
            z = int(z)
            if x<=xmin or x>=xmax or y<=ymin or y>=ymax:  # if click outside image slice, don't store coords
                return
        except TypeError:
            return

        try:
            if view == 's':
                X[x][(y - 2):(y + 2), (x - 2):(x + 2)] = 1300  # for sagittal view
            elif view == 'c':
                X[y][(z - 2):(z + 2), (x - 2):(x + 2)] = 1300  # for coronal view
            elif view == 'a':
                X[z][(y - 2):(y + 2), (x - 2):(x + 2)] = 1300  # for axial view
            else:  # default to coronal view
                X[y][(z - 2):(z + 2), (x - 2):(x + 2)] = 1300

            global count
            count = count+1
            print(f'{count}) ', end='')
            if view == 's':
                print(f'x = {z}, y = {ix:.2f}, z = {iy:.2f}') # for sagittal view
            elif view == 'c':
                print(f'x = {ix:.2f}, y = {z}, z = {iy:.2f}')  # for coronal view
            elif view == 'a':
                print(f'x = {ix:.2f}, y = {iy:.2f}, z = {z}') # for axial view
            else: # default to coronal view
                print(f'x = {ix:.2f}, y = {z}, z = {iy:.2f}')

            global coords
            if view == 's':
                coords.append((z, ix, iy)) # for sagittal view
            elif view == 'c':
                coords.append((ix, z, iy))  # for coronal view
            elif view == 'a':
                coords.append((ix, iy, z)) # for axial view
            else: # default to coronal view
                coords.append((ix, z, iy))

        except IndexError: # if clicked out of bounds of img slice, no coords to collect, return
            return

        return coords

# The function to be called anytime a slider's value changes
def update(val):
    tracker.ind = int(amp_slider.val)
    fig.canvas.mpl_connect('scroll_event', tracker.update())  # REconnect plot to scrollable slice viewer

def global_coords(coords):
    ## Need to convert coords to global positions first (LPI)
    direction = image.GetDirection()
    # Tairawhitit - multiply by -1 bc flipped from RAI to LPI
    spacing = image.GetSpacing() # keep spacing as abs values
    origin = image.GetOrigin()

    z_spacing = spacing[2] * -1
    z_offset = origin[2]
    y_spacing = spacing[1] * -1
    y_offset = origin[1] #* -1
    x_spacing = spacing[0] * -1
    x_offset = origin[0] #* -1

    coords = np.asarray(coords) #convert list to np array

    for i in range(len(coords)):
        # translate z coord
        coords[i][2] = coords[i][2] * z_spacing + z_offset
        # translate y coord
        coords[i][1] = coords[i][1] * y_spacing + y_offset
        # translate x coord
        coords[i][0] = coords[i][0] * x_spacing + x_offset

    return coords

def sort_centreline_nodes(coords):
    try:
        ## Arrange selected landmarks in order of: trachea nodes, left lung nodes, right lung nodes.
        centreline = np.zeros((14,3))

        # Sort along z axis in descending order. Airway entry is the 1st node; as node number increases, z value decreases (growing downwards into the lungs)
        coords = coords[coords[:, 2].argsort()] # sorts in ascending order
        coords = np.flip(coords, axis=0) # flip array vertically

        # Assign selected trachea nodes
        ## in CMISS template, node 1 is airway opening, node 3 is trachea midpoint, node 5 is trachea bifurcation/carina
        centreline[0] = coords[0]
        centreline[2] = coords[1]
        centreline[4] = coords[2]

        # Generate trachea nodes
        ## find Node 2, the midpoint of Node 1 & Node 3
        centreline[1] = [(centreline[0][0]+centreline[2][0])/2.0, (centreline[0][1]+centreline[2][1])/2.0, (centreline[0][2]+centreline[2][2])/2.0]
        # find Node 4, the midpoint of Node 3 & Node 5
        centreline[3] = [(centreline[2][0]+centreline[4][0])/2.0, (centreline[2][1]+centreline[4][1])/2.0, (centreline[2][2]+centreline[4][2])/2.0]
        print("Trachea nodes assigned.\n")

        # delete processed trachea landmarks from coords
        coords = np.delete(coords, slice(0,3), axis=0)


        # global orientation is LPI
        # Separate right & left lung nodes based on trachea bifurcation node.
        right_lung_nodes = []
        left_lung_nodes = []
        for i in range(len(coords)):
            if coords[i][0] > centreline[4][0]:
                right_lung_nodes.append(coords[i])
            else:
                left_lung_nodes.append(coords[i])

        ## convert to np array
        right_lung_nodes = np.asarray(right_lung_nodes)
        left_lung_nodes = np.asarray(left_lung_nodes)

        #--------------------- Start right lung. -----------------

        # Sort along z axis in descending order.
        right_lung_nodes = right_lung_nodes[right_lung_nodes[:, 2].argsort()] # sorts in ascending order
        right_lung_nodes = np.flip(right_lung_nodes, axis=0) # flip array vertically

        # Top 2 nodes are right bronchi bifurcation & RUL node
        ## in CMISS template, node 8 is right bronchi bif, node 11 is RUL node.
        ## sort by x value - x_RUL > x_RBronchi
        if right_lung_nodes[0][0] > right_lung_nodes[1][0]:
            centreline[10] = right_lung_nodes[0]
            centreline[7] = right_lung_nodes[1]
        elif right_lung_nodes[0][0] < right_lung_nodes[1][0]:
            centreline[10] = right_lung_nodes[1]
            centreline[7] = right_lung_nodes[0]
        print('Right bronchus bifurcation and RUL node assigned.')

        # RML-RLL bifurcation node is 3rd node
        centreline[11] = right_lung_nodes[2]
        print('RML-RLL bifurcation assigned.')

        # Bottom 2 nodes are RML & RLL
        ## in CMISS template, node 13 is RML, node 14 is RLL
        ## sort by y value = y_RML > y_RLL
        if right_lung_nodes[3][1] > right_lung_nodes[4][1]:
            centreline[12] = right_lung_nodes[3]
            centreline[13] = right_lung_nodes[4]
        elif right_lung_nodes[3][1] < right_lung_nodes[4][1]:
            centreline[12] = right_lung_nodes[4]
            centreline[13] = right_lung_nodes[3]
        print('RML & RLL nodes assigned.')
        print('All right lung nodes sorted.\n')

        #--------------------- Right lung done. -----------------
        #--------------------- Start left lung. -----------------

        # Sort along z axis in descending order.
        left_lung_nodes = left_lung_nodes[left_lung_nodes[:, 2].argsort()] # sorts in ascending order
        left_lung_nodes = np.flip(left_lung_nodes, axis=0) # flip array vertically

        # Top node is left bronchus node.
        ## in CMISS template, node 6 is left bronchi node.
        centreline[5] = left_lung_nodes[0]
        print('Left bronchus node assigned.')
        left_lung_nodes = np.delete(left_lung_nodes, 0, axis=0)

        ## in CMISS template, node 7 is left bronchi bif, node 9 is LUL, node 10 is LLL
        ## sort by x value - x_LBifurcation > x_LUL & x_LLL
        idx = np.argmax(left_lung_nodes, axis=0)[0]
        centreline[6] = left_lung_nodes[idx] # argmin finds the idx of min vals in each col i.e., min x, y,z values. 0th col is x col, so take [0]
        print('Left bronchus bifurcation assigned.')

        # remove processed nodes
        left_lung_nodes = np.delete(left_lung_nodes, idx, axis=0)

        ## sort LUL & LLL nodes by z values - z_LUL > z_LLL
        if left_lung_nodes[0][2] > left_lung_nodes[1][2]:
            centreline[8] = left_lung_nodes[0]
            centreline[9] = left_lung_nodes[1]
        elif left_lung_nodes[0][2] < left_lung_nodes[1][2]:
            centreline[8] = left_lung_nodes[1]
            centreline[9] = left_lung_nodes[0]
        print('LUL & LLL nodes assigned.')
        print('All left lung nodes sorted.\n')

        #--------------------- Left lung done. -----------------

    except IndexError:  # if wrong number of nodes, exit() to cancel operation
        print('Cancelled. IndexError.')
        exit()

    return centreline


def get_node_translation(moving, fixed, node):
    translation = np.zeros(3)
    translation[0] = fixed[node - 1][0] - moving[node - 1][0]
    translation[1] = fixed[node - 1][1] - moving[node - 1][1]
    translation[2] = fixed[node - 1][2] - moving[node - 1][2]

    return translation

################################################################################################################
if __name__ == '__main__':
    # Flag subject code
    parser = argparse.ArgumentParser(description='input subject stl files')
    parser.add_argument('-s', '--subject', help='Input subject name (code)', required=True)
    parser.add_argument('-view', '--view', help = 'Indicate coronal (c), axial (a), or sagittal (s) view')
    args = parser.parse_args()
    subject = args.subject
    view = args.view

    # load in image w/ SimpleITK
    path_hpc = '/run/user/57900/gvfs/sftp:host=hpc5.bioeng.auckland.ac.nz,user=msoo935/hpc/msoo935/'
    path_local = '/home/msoo935/Code/MRI-to-model_platform_v1/subjects/'
    # proj_folder = 'MRI-to-model_platform_v1/subjects'
    # image = sitk.ReadImage(os.path.join(path_hpc, proj_folder, subject, 'metaimage', subject + '_lung_raw.mha'))
    # proj_folder = 'Human_Lung_Atlas_EI'
    # image = sitk.ReadImage(os.path.join(path_hpc, proj_folder, subject + '_raw.nii'))
    proj_folder = 'Human_Lung_Atlas_EE'
    image = sitk.ReadImage(os.path.join(path_hpc, proj_folder, subject, 'metaimage', subject + '_raw.mha'))
    image = sitk.DICOMOrient(image, 'LPI') # SITK reorient to LPI
    # NB: in SITK, 'LPI' = 'RAS'. In global_coords(), spacing * -1 thus 'RAS' flipped back to 'LPI' to match sort_centreline_coords() assumption (as in ITK-Snap).
    size = image.GetSize()
    print(size)
    if view=='a':
        slices = size[2]
    elif view=='c':
        slices = size[0]
    elif view=='s':
        slices = size[1]
    else:
        slices = size[0]
    X = sitk.GetArrayFromImage(image)


    coords = []  # initialise array to store clicked coords
    count=0 # initialise counter for number of coords collected

    # open a figure window
    fig, ax = plt.subplots()
    # adjust the main plot to make room for the sliders
    fig.subplots_adjust(left=0.25)
    # Make a vertically oriented slider to control the amplitude
    axamp = fig.add_axes([0.1, 0.25, 0.0225, 0.63])
    amp_slider = Slider(
        ax=axamp,
        label="Slice",
        valmin=0,
        valmax=slices-1,
        valinit=0,
        orientation="vertical",
        initcolor='None',
        valfmt='%0.0f' # set to discrete integer values
    )

    # register the update function with each slider
    amp_slider.on_changed(update)

    ## Call scrollable slice viewer
    tracker = IndexTracker(ax, X)
    fig.canvas.mpl_connect('scroll_event', tracker.onscroll) # connect plot to scrollable slice viewer
    fig.canvas.mpl_connect('button_press_event', tracker.onclick) # connect plot to coord collector
    # plt.axis([max(x), min(x), max(y), min(y)])
    plt.show()

    if len(coords) != 12:
        print('Incomplete set of nodes.\nCancelled.')
        exit()

    print("Collected slice coordinates.")


    coords = global_coords(coords) # apply origin offset and spacing
    print('Derived spatial coordinates.')
    centreline = sort_centreline_nodes(coords)
    print('Derived upper airway centreline template.')

    ## Write exnode file
    output_ex_path = os.path.join(path_hpc, proj_folder, subject, 'centreline_v4')
    if not os.path.exists(output_ex_path):
        os.makedirs(output_ex_path)
    array2exnode(os.path.join(output_ex_path, subject+'_centrelinev4'), centreline)

    output_ip_path = os.path.join(path_hpc, proj_folder, subject, 'centreline_v4')
    if not os.path.exists(output_ip_path):
        os.makedirs(output_ip_path)
    array2ipnode(os.path.join(output_ip_path, subject+'centrelinev4'), centreline)
