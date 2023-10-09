#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "Dampled spin dynamics of a 5 atoms Fe chain at several lattice parameters"

# set the needed environment variables
. ../../../environment_variables

for a in 2.05 2.10 2.15 2.20 2.25 2.30 ; do

mkdir a$a

$ECHO "a= $a"

cat > in_master.txt<<EOF
&units
 energy = 'eV'
 length = 'ang'
 time = 'fs'
 mass='hau'
 /
&calculation
 processing = 'sd'
 post_processing = 'txt2xyz'
 post_processing_dir = 'sd'
 /
&element
 ne = 1
 symbol(1) = 'Fe'
 q(1) = 8.0
 q_d(1) = 7.0
 u_lcn(1) = 20.0000000000000000
 i_stoner_d(1) = 0.95
 xi_so_d(1) = 5.9999999999999998E-002
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = 10.000000000000000, 0.0000000000000000, 0.0000000000000000
 v(2,:) = 0.0000000000000000, 10.000000000000000, 0.0000000000000000
 v(3,:) = 0.0000000000000000, 0.0000000000000000, 5.0000000000000000
 /
&atom
 ns = 4
 na = 5
 ntag = 1
 stag(1) = 5
 tag(1) = 'Fe_wire'
 pbc = 0, 0, 1
 r(1,:) = 0.0000000000000000, 0.0000000000000000, 0.0
 r(2,:) = 0.0000000000000000, 0.0000000000000000, 0.2
 r(3,:) = 0.0000000000000000, 0.0000000000000000, 0.4
 r(4,:) = 0.0000000000000000, 0.0000000000000000, 0.6
 r(5,:) = 0.0000000000000000, 0.0000000000000000, 0.8
 m_listing = 'by_atom'
 m(1,:) = 3.1800000000000000, 1.0000000000000000, 0.0000000000000000
 m(2,:) = 3.1800000000000000, 0.0000000000000000, 0.0000000000000000
 m(3,:) = 3.1800000000000000, 0.0000000000000000, 0.0000000000000000
 m(4,:) = 3.1800000000000000, 0.0000000000000000, 0.0000000000000000
 m(5,:) = 3.1800000000000000, 0.0000000000000000, 0.0000000000000000
 lambda_pen_listing= 'by_tag'
 lambda_pen(1) = 30.0
 m_coord = 'spherical'
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
 verbose = .false.
 /
&sd
 integrator = 'st_1'
 t_i = 0
 t_f = 30
 dt = 0.05
 fixed_time_step = .true.
 alpha = 1.0
 temp = 0.0
 verbose = .true.
 /
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

cp -rf sd a$a

done

