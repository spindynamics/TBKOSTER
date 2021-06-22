#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use DyNaMol.x to calculate a Fe cluster"

# set the needed environment variables
. ../../environment_variables

rm -rf scf/*.*
rm -f out*
mkdir scf
rm -rf xyz
mkdir xyz

cat > in_master.txt<<EOF
&calculation
 processing = 'scf'
 post_processing = 'txt2xyz'
 post_processing_dir = 'scf'
 /
&units
 energy = 'eV'
 length = 'ang'
 time = 'fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Cr'
 q(1) = 6.0
 q_d(1) = 5.0
 u_lcn(1) = 20.0000000000000000
 i_stoner_d(1) = 0.82
 xi_so_d(1) = 0.0
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/cr_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = 2.
 v(1,:) = 1.000000000000000, 0.0000000000000000, 0.0000000000000000
 v(2,:) = 0.0000000000000000, 1.000000000000000, 0.0000000000000000
 v(3,:) = 0.0000000000000000, 0.0000000000000000, 1.0000000000000000
 /
&atom
 ns = 4
 na = 3
 ntag = 1
 stag(1) = 3
 tag(1) = 'Cr_trimer'
 pbc = 0, 0, 0
 r(1,:) = 0.0000000000000000, 0.0000000000000000, 0.0
 r(2,:) = 0.86602540378443864676,-0.5, 0.0
 r(3,:) = 0.86602540378443864676, 0.5, 0.0
 m_listing = 'by_atom'
 m_coord = 'spherical'
 m(1,:) = 3.1800000000000000, 30.0000000000000000,    180.0000000000000000
 m(2,:) = 3.1800000000000000, 30.0000000000000000,   -60.0000000000000000
 m(3,:) = 3.1800000000000000, 30.0000000000000000,    60.0000000000000000
 lambda_pen_listing= 'by_tag'
 lambda_pen(1) = 0.0
 m_coord = 'spherical'
 /
&mesh
 type = 'mp'
 gx = 1, 1, 1
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
 verbose = .true.
 ni_max = 200
 /
EOF

# Set DyNaMol root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run DyNaMol
$BIN_DIR/DyNaMol.x 
