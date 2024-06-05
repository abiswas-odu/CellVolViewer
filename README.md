## CellVolViewer for 3D Epithilial Cell Layer Reconstruction (Latest Release 3.11)

Documentation on running the the 3D cell layer reconstruction tool, CellVolViewer. 

### Download the Code

1.	Obrain the latest release from `https://github.com/abiswas-odu/CellVolViewer/releases`. 

2. Unzip the file downloaded file and extract the release. Here we assume you are building in new subdirectory in your home directory:

```shell
   mkdir ~/CellVolViewer 
   cd ~/CellVolViewer 
   tar -xvf ../CellVolViewer.tar.gz
```

Note the CellVolViewer tar file will be have a different name based on the version number.

2.  Let us assume that the unzipped path on the filesystem is PATH. This is the install path. 

### Data Processing

1. Convert the z-stack of the image from the microscope into individual images in TIF or PNG format in Fiji. 

2. Each of those images should be converted into 8-bit greyscale images and saved in a folder with the z-stack index as the filename. For example, 1.tif, 2.tif ... .

3. Each of the 8-bit greyscale images need to be segmented using CellViewer(see above) to generate the masks. This will produce a directory structure with a folder for each image and the handCorrected.tif mask inside it. 

4. Please see examples of the expected folder structure in the install path. For example "test_2".

5. Now that the folder structure is setup. Let us assume it is called DATA. Create a text file inside DATA named config.txt. The config.txt contains various input parameters reqired to run the CellVolViewer commands. The expected parameters are described below. Please see an example of a config.txt inside the test folders. 

```
// Image Pixel Dimensions
// x-axis in pixels
{x_px} 
1024
// y-axis in pixels
{y_px}
1024
// number of z frames
{z_stack}
15
// number of interpolated frames in the 3D volume 
{inprep_frames}
10

// Real Dimensions of the image in microns 
{x_real}
176.68
{y_real}
176.68
{z_real}
15

// Placode area box. Dimemsion of a box that captures the placode cells only
{x_placode_min}
20
{x_placode_max}
1000
{y_placode_min}
20
{y_placode_max}
1000
```

### How to Run 

Prerequisites: MATLAB > 2019b

1.  On MAC OSX you need to grant permission for the app to run. To grant permissions launch the ```Terminal``` app from Launchpad and execute the following commands:
    ```
    cd PATH
	chmod a+x gentracemap.exe
    ```

2.  Open the install folder PATH in MATLAB. In the MATLAB commandline run ```gen_3d_vol DATA```. 

3.  Run the command ```gen_slice_view DATA``` to look at the sliceView of the cells. The cells are tracked across the frames and the same cell across frames id colored the same.

4.  The cell tracking can be corrected by editing the ```DATA\output\cell_tracks.csv```. The file has the same number of columns as the frames in the z-stack. The values are the cell IDs corresponding to the labels in the slice view.

5.  The masks can also be corrected to perform better 2D segmentation. Run ```gen_slice_view DATA``` again to show and label the new cells and assign them to tracks in cell_tracks.csv. 

6.  After the cell tracking is complete, run ```gen_cell_metrics DATA``` to generate the final metics.
    a.  It generates a file DATA\output\cell_areas.csv which has the cross section area of the cells in each frame and volumes. 
    b.  It generates a plot of the cross sectional area of each cell across each frame. 
7.  To generate the 3D volume image of the cells, run ```show_volume DATA```. It generates a 3D volume re-construction of the cells. 
