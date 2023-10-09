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
. ../../environment_variables

mkdir scf
rm -f out*
rm -f  results/Etot_vs_theta.dat

#a=$(echo "0.5291771*2.77*e(l(8.0*3.14159/3.0)/3.0)" |bc -l)
a=3.1

for a in 2.8 2.9  3.0  3.1  ; do
rm -f in_charge.txt  tempo tempo2
echo a=$(echo "$a" | bc -l)

cat >> results/Etot_vs_theta.dat << EOF
@#  theta   Etot  a=$a
EOF

for theta in  0 5 10 20 30 40 50 60 70 80 90 ; do 
#for theta in  0  ; do 

echo "theta=$theta"


cat > in_master.txt<<EOF
&units
energy = 'eV'
length = 'ang'
time = 'fs'
mass='hau'
/
&calculation
processing = 'scf'
/
&element
 ne = 2
 symbol(1) = 'Fe'
 symbol(2) = 'Rh'
 u_lcn(1) = 20.0
 i_stoner_d(1) = 0.95
 u_lcn(2) = 20.0
 i_stoner_d(2) = 0.85
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
 filename(2) = '$TBPARAM_DIR/rh_par_fcc_bcc_sc_gga_fl'
 /
&lattice
v_factor = $a
 v(1,:) = 0.0 1.0 1.0
 v(2,:) = 1.0 0.0 1.0
 v(3,:) = 1.0 1.0 0.0
/
&atom
 ns = 4
 na = 4
 ntag = 3
 stag(1) = 1
 stag(2) = 1
 stag(3) = 2
 tag(1) = 'Fe1'
 tag(2) = 'Fe2' 
 tag(3) = 'Rh'
 r_coord='cartesian'
 r(1,:) = 0.0, 0.0, 0.0
 r(2,:) = $a,  $a, $a
 r(3,:) = $(echo "$a/2" |bc -l), $(echo "$a/2" |bc -l), $(echo "$a/2" |bc -l) 
 r(4,:) = $(echo "3*$a/2" |bc -l), $(echo "3*$a/2" |bc -l),$(echo "3*$a/2" |bc -l)  
 m_listing = 'by_tag'
 m_coord = 'spherical'
 m(1,:) = 3.0, $(echo "180-$theta" |bc -l), 0.0
 m(2,:) = 3.0,  $theta,  0.0
 m(3,:) = 0.0,  90.0, 0.0
 lambda_pen_listing = 'by_tag'
 lambda_pen(1) = 10.0
 lambda_pen(2) = 10.0
 lambda_pen(3) = 0.0 
 /
&mesh
type = 'mp'
gx = 10, 10, 10
dx = 0, 0, 0
/
&hamiltonian_tb
m_penalization = 'theta,phi'
/
&energy
smearing = 'mv'
degauss = 0.05
/
&mixing
alpha = 0.2
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

cp -f out_charge.txt in_charge.txt

cat > tempo << EOF
theta= $theta
EOF

cat tempo out_energy.txt>>tempo2


done
grep -e 'theta=' -e 'en =' tempo2 | awk '/theta/{theta = $(NF)}/en/{print theta, $(NF)}' >> results/Etot_vs_theta.dat

rm -f tempo tempo2

done

