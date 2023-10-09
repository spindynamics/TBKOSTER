#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate a Fe cluster"

# set the needed environment variables
. ../../../environment_variables

rm -f out*

cat > in_master.txt<<EOF
&calculation
 processing = 'scf'
 post_processing= 'txt2xyz'
 post_processing_dir= 'scf'
 /
&units
 length='ang'
 energy='ev'
 time='fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Fe'
 q(1)   = 8.0
 q_d(1) = 7.0
 u_lcn(1)=20
 i_stoner_d(1) = 0.95
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/fe_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor =1.0
 v(1,:) = 1 0 0
 v(2,:) = 0 1 0
 v(3,:) = 0 0 1 
 /
&atom
 ns = 2
 na = 13
 ntag = 1
 tag(1) = 'Fe'
 stag(1)=13
 pbc = 0, 0, 0
 r_coord='direct'
 r(1,:) =      0.0000000     0.0000000     0.0000000
 r(2,:) =      0.0000000     1.7575089     1.7575089
 r(3,:) =      0.0000000    -1.7575089     1.7575089
 r(4,:) =      0.0000000     1.7575089    -1.7575089
 r(5,:) =      0.0000000    -1.7575089    -1.7575089
 r(6,:) =      1.7575089     0.0000000     1.7575089
 r(7,:) =     -1.7575089     0.0000000     1.7575089
 r(8,:) =      1.7575089     0.0000000    -1.7575089
 r(9,:) =     -1.7575089     0.0000000    -1.7575089
 r(10,:) =     1.7575089     1.7575089     0.0000000
 r(11,:) =    -1.7575089     1.7575089     0.0000000
 r(12,:) =     1.7575089    -1.7575089     0.0000000
 r(13,:) =    -1.7575089    -1.7575089     0.0000000
 m_listing = 'by_tag'
 m(1,:) = 1.0, 0.0, 0.0
 lambda_pen_listing = 'by_atom'
 lambda_pen(1)=10
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
 delta_en=0.0001
 delta_q=0.0001
 verbose=.true.
 ni_max=200
 /
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 
