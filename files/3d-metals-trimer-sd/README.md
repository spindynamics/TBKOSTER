# TBKOSTER
Tests made for the triangular trimers of Cr, Fe and Co and a comparisson with UppASD code.

Data found in this folder for a simulation with and without damping. Base job file can be found in job.sh in this folder.

Data for TBKOSTER can be read as
#timestep #mx #my #mz

Data for UppASD can be read as
#stepnumber #atomnumber #mx #my #mz #mnorm

Note that the timstep of TBKOSTER is synchronized with the stepnumber of UppASD.
Note that #mx, #my and #mz of UppASD is normalized, whereas it is not on TBKOSTER.
