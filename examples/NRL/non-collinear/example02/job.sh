#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "SCF NON-collinear spin calculation of a 5-atom Fe wire with magnetic penalization on atom 1"
$ECHO "The penalization on atome 1 is theta,phi=(30,0) and the other atoms are free "
$ECHO "during scf the spin tend do align and the magnetization of atom 1 is the same as the other atoms"

# set the needed environment variables
. ../../../environment_variables

mkdir scf
rm -f out*

cat > in_master.txt<<EOF
&units
energy = 'ev'
length = 'ang'
time = 'fs'
mass='hau'
/
&calculation
processing = 'scf'
post_processing='txt2xyz'
post_processing_dir='scf'
/
&element
ne = 1
symbol(1) = 'Fe'
q(1)   = 8.0
q_d(1) = 7.0
u_lcn(1) = 20.0
i_stoner_d(1) = 0.95
xi_so_d(1) = 0.06
/
&element_tb
filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
/
&lattice
v_factor = 2.27
v(1,:) = 10.0,  0.0,  0.0
v(2,:) =  0.0, 10.0,  0.0
v(3,:) =  0.0,  0.0,  5.0
/
&atom
ns = 4
na = 5
ntag = 2
stag(1) = 1
stag(2)= 4
tag(1) = 'Fe_pen'
tag(2)='Fe'
pbc = 0, 0, 1
r(1,:) = 0.0, 0.0, 0.0
r(2,:) = 0.0, 0.0, 0.2
r(3,:) = 0.0, 0.0, 0.4
r(4,:) = 0.0, 0.0, 0.6
r(5,:) = 0.0, 0.0, 0.8
m_listing = 'by_tag'
m_coord = 'spherical'
m(1,:) = 3.0, 30.0, 0.0
m(2,:) = 3.0, 0.0, 0.0
lambda_pen_listing = 'by_tag'
lambda_pen(1) = 5.0
/
&mesh
type = 'mp'
gx = 1, 1, 100
dx = 0, 0, 0
/
&hamiltonian_tb
m_penalization = 'theta,phi'
/
&energy
smearing = 'mv'
degauss = 0.2
/
&mixing
alpha = 0.1
/
&scf
delta_en = 0.00001
delta_q  = 0.00001
verbose = .true.
ni_max=200
/
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 
