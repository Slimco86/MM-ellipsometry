# MM-ellipsometry
Small library, which allows fast visualization and further decomposition of measured Mueller Matrices.
Written on Python3.

Implements a basic class of Mueller matrix MM, which serves as a container for Mueller matrix data. The data can be aploaded from different sources such as *.dat files created and formated by Woolam CompleteEASE ellipsometer containing the Mueller Matrix measurement data, or *.txt files with simulations, with one or several angles of incidence and one or complete azimuthal rotation of the smaple with the step of N degrees. 

Allows to easily perform Cloude, Differential and Polar decomposition of the matrices and visualize the result in fast and convinient manner.

For example one can visualize the data at specific azimuthal or incidence angle of specific elements, or plot the complete Mueller matrix.
![alt text](https://github.com/Slimco86/MM-ellipsometry/blob/master/Depol.png)

Complete azimuthal rotation can also be easily visualized:

![alt text](https://github.com/Slimco86/MM-ellipsometry/blob/master/Fig2.png)

The depolarization index, or degree of polarization can be easily calculated with dedicated fuctions. The outcomming Stokes vector can be calculated given an incoming Stokes vector. One can also check the physicality of the Mueller matrix and find its Hermitian representation.






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
