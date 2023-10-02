#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "Fe monatomic wire: calculation of spin spiral by two different ways:"
$ECHO "non-collinear calculation with 18 atoms per unit cell and theta_pen(i)=20*(i-1) i=1,18 "
$ECHO "spin spiral calculation with one atom per unit-cell and  k_spiral=(0,0,1/2) (xyz)=(-1/2,1/2,1/2) (direct)"

# set the needed environment variables
. ../../environment_variables

rm -f out*
rm -rf results
mkdir results
rm -rf scf
mkdir scf

$ECHO "super-cell (18 atoms) calculation"

cat > in_master.txt<<EOF
&units
energy = 'ev'
length = 'ang'
time = 'fs'
mass='hau'
/
&calculation
processing = 'scf'
/
&element
ne = 1
symbol(1) = 'Fe'
q(1)   = 8.0
q_d(1) = 7.0
u_lcn(1) = 20.0
i_stoner_d(1) = 0.95
xi_so_d(1) = 0.0
/
&element_tb
filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
/
&lattice
v_factor = 2.20
v(1,:) = 10.0,  0.0,  0.0
v(2,:) =  0.0, 10.0,  0.0
v(3,:) =  0.0,  0.0,  18
/
&atom
ns = 4
na = 18
ntag = 1
stag(1) = 18
tag(1) = 'Fe'
pbc = 0, 0, 5
r(1,:)  = 0.0, 0.0, 0.0
r(2,:)  = 0.0, 0.0, $(echo "1.0/18" |bc -l)
r(3,:)  = 0.0, 0.0, $(echo "2.0/18" |bc -l)
r(4,:)  = 0.0, 0.0, $(echo "3.0/18" |bc -l)
r(5,:)  = 0.0, 0.0, $(echo "4.0/18" |bc -l)
r(6,:)  = 0.0, 0.0, $(echo "5.0/18" |bc -l)
r(7,:)  = 0.0, 0.0, $(echo "6.0/18" |bc -l)
r(8,:)  = 0.0, 0.0, $(echo "7.0/18" |bc -l)
r(9,:)  = 0.0, 0.0, $(echo "8.0/18" |bc -l)
r(10,:) = 0.0, 0.0, $(echo "9.0/18" |bc -l)
r(11,:) = 0.0, 0.0, $(echo "10.0/18" |bc -l)
r(12,:) = 0.0, 0.0, $(echo "11.0/18" |bc -l)
r(13,:) = 0.0, 0.0, $(echo "12.0/18" |bc -l)
r(14,:) = 0.0, 0.0, $(echo "13.0/18" |bc -l)
r(15,:) = 0.0, 0.0, $(echo "14.0/18" |bc -l)
r(16,:) = 0.0, 0.0, $(echo "15.0/18" |bc -l)
r(17,:) = 0.0, 0.0, $(echo "16.0/18" |bc -l)
r(18,:) = 0.0, 0.0, $(echo "17.0/18" |bc -l)
m_listing = 'by_atom'
m_coord = 'spherical'
m(1,:) =  3.0, 90.0, 0.0
m(2,:) =  3.0, 90.0, $(echo "1.0*20" |bc -l)  
m(3,:) =  3.0, 90.0, $(echo "2.0*20" |bc -l)
m(4,:) =  3.0, 90.0, $(echo "3.0*20" |bc -l)
m(5,:)  = 3.0, 90.0, $(echo "4.0*20" |bc -l)
m(6,:)  = 3.0, 90.0, $(echo "5.0*20" |bc -l)
m(7,:)  = 3.0, 90.0, $(echo "6.0*20" |bc -l)
m(8,:)  = 3.0, 90.0, $(echo "7.0*20" |bc -l)
m(9,:)  = 3.0, 90.0, $(echo "8.0*20" |bc -l)
m(10,:) = 3.0, 90.0, $(echo "9.0*20" |bc -l)
m(11,:) = 3.0, 90.0, $(echo "10.0*20" |bc -l)
m(12,:) = 3.0, 90.0, $(echo "11.0*20" |bc -l)
m(13,:) = 3.0, 90.0, $(echo "12.0*20" |bc -l)
m(14,:) = 3.0, 90.0, $(echo "13.0*20" |bc -l)
m(15,:) = 3.0, 90.0, $(echo "14.0*20" |bc -l)
m(16,:) = 3.0, 90.0, $(echo "15.0*20" |bc -l)
m(17,:) = 3.0, 90.0, $(echo "16.0*20" |bc -l)
m(18,:) = 3.0, 90.0, $(echo "17.0*20" |bc -l)
/
&mesh
type = 'mp'
gx = 1, 1, 20
dx = 0, 0, 0
/
&hamiltonian_tb
m_penalization = 'none'
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
ni_max=500
/
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

cp -f out_log.txt results/out_log_super_cell.txt


$ECHO "spin-spiral calculation"

cat > in_master.txt<<EOF
&units
energy = 'ev'
length = 'ang'
time = 'fs'
mass='hau'
/
&calculation
processing = 'scf'
/
&element
ne = 1
symbol(1) = 'Fe'
q(1)   = 8.0
q_d(1) = 7.0
u_lcn(1) = 20.0
i_stoner_d(1) = 0.95
xi_so_d(1) = 0.0
/
&element_tb
filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
/
&lattice
v_factor = 2.20
v(1,:) = 10.0,  0.0,  0.0
v(2,:) =  0.0, 10.0,  0.0
v(3,:) =  0.0,  0.0,  1.0
/
&atom
ns = 4
na = 1
ntag = 1
stag(1) = 1
tag(1) = 'Fe'
pbc = 0, 0, 5
k_spiral = 0, 0, 0.055
r(1,:) = 0.0, 0.0, 0.0
m_listing = 'by_tag'
m_coord = 'spherical'
m(1,:) = 3.0, 90.0, 0.0
/
&mesh
type = 'mp'
gx = 1, 1, 360
dx = 0, 0, 0
/
&hamiltonian_tb
m_penalization = 'none'
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
ni_max=500
/
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

cp -f out_log.txt results/out_log_spin_spiral.txt


