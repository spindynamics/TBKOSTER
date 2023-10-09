#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate band structure of a Fe wire with SOC"

# set the needed environment variables
. ../../environment_variables

rm -rf results
mkdir results
rm -f out*

IStoner=0.95


for e_e_interaction in stoner  ; do

$ECHO 'electronic interaction' $e_e_interaction
cat > results/Etot_vs_theta.$e_e_interaction.dat << EOF
@# theta  Etot(eV)
EOF


for theta in 0 10 20 30 40 50 60 70 80 90 ; do

$ECHO 'theta=' $theta
cat > in_master.txt<<EOF
&units
energy = 'ev'
length = 'ang'
time = 'fs'
mass='hau'
/
&calculation
processing = 'scf'
post_processing='band'
/
&element
ne = 1
symbol(1) = 'Fe'
q(1)   = 8.0
q_d(1) = 7.0
u_lcn(1) = 20.0
i_stoner_d(1) = $IStoner
xi_so_d(1) = 0.06
/
&element_tb
filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
/
&lattice
v_factor = 2.25
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
r(1,:) = 0.0, 0.0, 0.0
m_listing = 'by_tag'
m_coord = 'spherical'
m(1,:) = 3.0, $theta, 0.0
lambda_pen_listing = 'by_tag'
lambda_pen(1) = 5.0
/
&mesh
type = 'mp'
gx = 1, 1, 500
dx = 0, 0, 0
/
&hamiltonian_tb
m_penalization = 'theta,phi'
e_e_interaction = '$e_e_interaction'
/
&energy
smearing = 'mv'
degauss = 0.05
/
&mixing
alpha = 0.1
/
&scf
delta_en = 0.00001
delta_q  = 0.00001
verbose = .false.
ni_max=500
/
EOF


cat > band/in_mesh.txt<<EOF
&mesh
 type = 'path'
 gxs = 100
 nxs = 2
 xs_label(1) = 'G'
 xs_label(2) = 'X'
 xs(1,:) =  0  ,  0.0 ,  0
 xs(2,:) =  0.0,  0.0  ,  0.5
 /
EOF

cat > band/in_energy.txt<<EOF
&energy
 smearing = 'mv'
 degauss = 0.2
 en_min = -10.0
 en_max =  10.0
 /
EOF

cat > band/in_dos.txt<<EOF
&dos
 nen=100
 na_dos=1
 ia= 1
en_min=-10
 en_max=10
 /
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

cp band/band.dat results/band$theta.dat

done


done
