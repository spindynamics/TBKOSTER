#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use DyNaMol.x to model Crbcc AF by two approaches"
$ECHO "super-cell spin collinear (up down) or non-collinear spin spiral with k_spiral in zone border"

# set the needed environment variables
. ../../environment_variables

rm -f out*
rm -rf results
mkdir results
mkdir scf

cat > in_master.txt<<EOF
&calculation
 processing='scf'
 /
&units
 energy = 'ev'
 length = 'ang'
 time = 'fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Cr'
 q(1)   = 6.0
 q_d(1) = 5.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 0.82
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/cr_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = 2.87
 v(1,:) = 1.0	0.0  0.0
 v(2,:) = 0.0	1.0  0.0
 v(3,:) = 0.0	0.0  1.0
 /
&atom
 ns = 2
 na = 2
 ntag = 2
 stag(1) = 1
 stag(2) = 1
 tag(1) = 'Cr_up'
 tag(2) = 'Cr_dn'
 r(1,:) = 0.0, 0.0, 0.0
 r(1,:) = 0.5, 0.5, 0.5
 m_listing = 'by_tag'
 m_coord = 'spherical'
m(1,:) =  1.0, 0.0, 0.0
m(2,:) = -1.0, 0.0, 0.0
 /
&mesh
 type = 'mp'
 gx = 15, 15, 15
 dx = 0, 0, 0
 /
&hamiltonian_tb
 /
&energy
 smearing = 'mv'
 degauss = 0.2
 /
&mixing
 alpha = 0.1
 /
&scf
 delta_en = 0.0001
 delta_q  = 0.0001
 verbose = .true.
 /
EOF

# Set DyNaMol root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run DyNaMol
$BIN_DIR/DyNaMol.x 

cp out_log.txt results/out_log_super_cell.txt


cat > in_master.txt<<EOF
&calculation
 processing='scf'
 /
&units
 energy = 'ev'
 length = 'ang'
 time = 'fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Cr'
 q(1)   = 6.0
 q_d(1) = 5.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 0.82
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/cr_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = 2.87
 v(1,:) = -0.5,  0.5,  0.5
 v(2,:) =  0.5, -0.5,  0.5
 v(3,:) =  0.5,  0.5, -0.5
 /
&atom
 ns = 4
 na = 1
 ntag = 1
 stag(1) = 1
 tag(1) = 'Cr'
 k_spiral= -0.5, 0.5, 0.5
 r(1,:) = 0.0, 0.0, 0.0
 m_listing = 'by_tag'
 m_coord = 'spherical'
m(1,:) =  1.0, 90.0, 0.0
 /
&mesh
 type = 'mp'
 gx = 30, 30, 30
 dx = 0, 0, 0
 /
&hamiltonian_tb
 /
&energy
 smearing = 'mv'
 degauss = 0.2
 /
&mixing
 alpha = 0.1
 /
&scf
 delta_en = 0.0001
 delta_q  = 0.0001
 verbose = .true.
 /
EOF

# Set DyNaMol root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run DyNaMol
$BIN_DIR/DyNaMol.x 

cp out_log.txt results/out_log_bcc.txt
