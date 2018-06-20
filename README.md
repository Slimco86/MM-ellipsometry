# MM-ellipsometry
My first project, which allow fast visualization and further decomposition of measured Mueller Matrices.
Written on Python3.

Takes *.dat files created by Woolam CompleteEASE ellipsometer containing the Mueller Matrix measurement data, 
with one or several angles of incidence and complete azimuthal rotation of the smaple with the step of 5 degrees. 

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
