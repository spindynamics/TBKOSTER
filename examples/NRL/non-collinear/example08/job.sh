#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "SCF NON-collinear spin calculation with SOC of a Fe monolayer"
$ECHO " The E(theta) curve is evaluated for various elecronic interaction"
$ECHO " Stoner and UJB "

# set the needed environment variables
. ../../../environment_variables

rm -rf results
mkdir results
rm -f out*

IStoner=0.95
U=$(echo "5.0/7.0*$IStoner" |bc -l)
J=$(echo "$U" |bc -l)
B=$(echo "0.14*$J" |bc -l)

for e_e_interaction in stoner   ; do

$ECHO 'electronic interaction' $e_e_interaction
cat > results/Etot_vs_theta.$e_e_interaction.dat << EOF
@# theta  Etot(eV)
EOF


rm -f tempo tempo2

#for theta in 0 10 20 30 40 50 60 70 80 90 ; do
for theta in 0  ; do
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
/
&element
ne = 1
symbol(1) = 'Fe'
q(1)   = 8.0
q_d(1) = 7.0
u_lcn(1) = 20.0
i_stoner_d(1) = $IStoner
u_dd(1)= $U
j_dd(1)=$J
b(1) = $B
xi_so_d(1) = 0.06
/
&element_tb
filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
/
&lattice
v_factor = 2.25
v(1,:) =  1.0,  0.0,  0.0
v(2,:) =  0.0,  1.0,  0.0
v(3,:) =  0.0,  0.0,  10
/
&atom
ns = 4
na = 1
ntag = 1
stag(1) = 1
tag(1) = 'Fe'
pbc = 5, 5, 0
r(1,:) = 0.0, 0.0, 0.0
m_listing = 'by_tag'
m_coord = 'spherical'
m(1,:) = 3.0, $theta, 0.0
lambda_pen_listing = 'by_tag'
lambda_pen(1) = 5.0
/
&mesh
type = 'mp'
gx = 20, 20, 1
dx = 0, 0, 0
/
&hamiltonian_tb
m_penalization = 'theta'
e_e_interaction = '$e_e_interaction'
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

cat > tempo << EOF
theta= $theta
EOF

cat tempo out_energy.txt>>tempo2
done

grep -e 'theta=' -e 'en =' tempo2 | awk '/theta/{theta = $(NF)}/en/{print theta, $(NF)}' >> results/Etot_vs_theta.$e_e_interaction.dat

rm -f tempo tempo2 

done
