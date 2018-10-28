# MM-ellipsometry
Small library, which allow fast visualization and further decomposition of measured Mueller Matrices.
Written on Python3.

Takes *.dat files created and formated by Woolam CompleteEASE ellipsometer containing the Mueller Matrix measurement data, 
with one or several angles of incidence and one or complete azimuthal rotation of the smaple with the step of N degrees. 

One can perform Cloude, Differential and Polar decomposition of the matrices and visualize the result in fast and convinient manner.

For example one can visualize the data at specific azimuthal or incidence angle of specific elements, or plot the complete Mueller matrix.
![alt text](https://github.com/Slimco86/MM-ellipsometry/blob/master/Depol.png)

Complete azimuthal rotation can also be easily visualized:

![alt text](https://github.com/Slimco86/MM-ellipsometry/blob/master/Fig2.png)




Basic GUI is avaliable if you run the GUI version.
User can choose between options:

1) Plot Mueller Matrices
plots 16 elements of the mueller matrix as polar plots with complete azimuthal roatation

2)Cloude decomposition:

Returnes wavelength dependent graph of the Mueller Matrix eigenvalues

3)Differential decomposition:

returns two matrices L_m and L_u. L_m is the mean value of the Mueller element and L_u is the uncertainty of this value.

4)Polar decomposition:
test
Returns Diatennuation matrix, retardance matrix, and depolarization matrix.
