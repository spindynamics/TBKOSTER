#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate band structure with SOC of Cufcc"

# set the needed environment variables
. ../../../environment_variables

rm -f tempo tempo2 tempo3 tempo4


a=3.61

$ECHO "a= $a"


cat > in_master.txt<<EOF
&calculation
 processing='scf'
 post_processing='band'
 /
&units
 energy = 'ev'
 length = 'ang'
 time = 'fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Cu'
 q(1)   = 11.0
 q_d(1) = 10.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 0.60
 xi_so_d(1) = 0.11
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/cu_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = 0.0 0.5 0.5
 v(2,:) = 0.5 0.0 0.5
 v(3,:) = 0.5 0.5 0.0
 /
&atom
 ns = 4
 na = 1
 ntag = 1
 stag(1) = 1
 tag(1) = 'Cu'
 r(1,:) = 0.0, 0.0, 0.0
 m_listing = 'by_tag'
 m_coord = 'spherical'
m(1,:) =  1.0, 0.0, 0.0
 /
&mesh
 type = 'mp'
 gx = 10, 10, 10
 dx = 0, 0, 0
 /
&hamiltonian_tb
 /
&energy
 smearing = 'mv'
 degauss = 0.1
 /
&mixing
 alpha = 0.1
 /
&scf
 delta_en = 0.0001
 delta_q  = 0.0001
 verbose = .false.
 ni_max = 100
 /
EOF

cat > band/in_mesh.txt<<EOF
&mesh
 type = 'path'
 gxs = 100
 nxs = 12
 xs_label(1) = 'G'
 xs_label(2) = 'X'
 xs_label(3) = 'W'
 xs_label(4) = 'K'
 xs_label(5) = 'G'
 xs_label(6) = 'L'
 xs_label(7) = 'U'
 xs_label(8) = 'W'
 xs_label(9) = 'L'
 xs_label(10) = 'K'
 xs_label(11) = 'U'
 xs_label(12) = 'X' 
 xs(1,:) =  0  ,  0.0 ,  0
 xs(2,:) =  0.5,  0.0  ,  0.5
 xs(3,:) =  0.5,  0.25 ,  0.75
 xs(4,:) =  0.365, 0.375  , 0.75
 xs(5,:) =  0  ,  0    ,  0
 xs(6,:) =  0.5, 0.5,  0.5
 xs(7,:) =  0.625  ,  0.25 ,  0.625
 xs(8,:) =  0.5,  0.25  ,  0.75
 xs(9,:) =  0.5,  0.5 ,  0.5
 xs(10,:) =  0.365, 0.375  , 0.75
 xs(11,:) =  0.625  ,  0.25 ,  0.625
 xs(12,:) =  0.5,  0.0  ,  0.5
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


