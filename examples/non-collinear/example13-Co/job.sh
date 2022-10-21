#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "SCF NON-collinear spin calculation of a 4-atom Fe wire"
$ECHO "The initial magnetization is a spin spiral of period 5a "
$ECHO "Due to periodic boundary conditions the spin spiral configuration is kept during scf"

# set the needed environment variables
. ../../environment_variables

rm -f out*
rm -rf results
mkdir results

rm -f tempo tempo2 

$ECHO "spin-spiral calculation"

cat > results/Etot_vs_kspiral.dat << EOF
@# kspiral  Etot
EOF

for k_spiral in 0.0 0.025 0.05 0.075 0.1 0.125 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 ; do

$ECHO "k_spiral=" $k_spiral
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
symbol(1) = 'Co'
q(1)   = 9.0
q_d(1) = 8.0
u_lcn(1) = 20.0
i_stoner_d(1) = 1.10
xi_so_d(1) = 0.0
/
&element_tb
filename(1) = '$TBPARAM_DIR/co_par_fcc_bcc_sc_gga_fl'
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
tag(1) = 'Co'
pbc = 0, 0, 5
k_spiral = 0, 0, $k_spiral
r(1,:) = 0.0, 0.0, 0.0
m_listing = 'by_tag'
m_coord = 'spherical'
m(1,:) = 3.0, 90.0, 0.0
/
&mesh
type = 'mp'
gx = 1, 1, 400
dx = 0, 0, 0
/
&hamiltonian_tb
m_penalization = 'none'
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

# Set DyNaMol root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run DyNaMol
$BIN_DIR/DyNaMol.x 

cat > tempo << EOF
k_spiral= $k_spiral
EOF

cat tempo out_energy.txt>>tempo2


done

grep -e 'k_spiral=' -e 'en =' tempo2 | awk '/k/{k = $(NF)}/en/{print k, $(NF)}' >> results/Etot_vs_kspiral.dat

rm -f tempo tempo2

