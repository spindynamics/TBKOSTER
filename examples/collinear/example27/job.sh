#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate the total energy vs a of Ptfcc"

# set the needed environment variables
. ../../environment_variables


for a in  3.94 ; do

  echo "a= $a"

cat > in_master.txt<<EOF
&calculation
 processing='scf'
 post_processing='band'
 /
&units
 length='ang'
 energy='ev'
 time='fs'
 mass='hau'
 /
&element
 ne = 2
 symbol(1) = 'Ti'
 symbol(2) = 'H'
 no(1) = 3 
 q(1) = 3.0
 o(1,1:3) =  5, 6, 7
 no(2)= 1
 o(2,1) =  5
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/dd_t2g_par'
 filename(2) = '$TBPARAM_DIR/dd_t2g_par'
 /
&lattice
 v_factor = $a
 v(1,:) = 1.0 0.0 0.0
 v(2,:) = 0.0 1.0 0.0
 v(3,:) = 0.0 0.0 1.0
 /
&atom
 ns = 1
 na = 2
 ntag = 2
 tag(1) = 'Ti'
 tag(2) = 'H' 
 stag(1)=1
 stag(2)=1
 pbc = 1, 1, 0
 r_coord='direct'
 r(1,:) = 0.0, 0.0, 0.0
 r(2,:) = 0.0, 0.0, 1.0 
 /
&mesh
 type = 'mp'
 gx = 10, 10, 1
 dx = 0, 0, 0
 /
&hamiltonian_tb
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

cat > band/in_energy.txt<<EOF
&energy
 smearing = 'mv'
 degauss = 0.2
 en_min = -10.0
 en_max =  10.0
 /
EOF

cat > band/in_mesh.txt<<EOF
&mesh
 type = 'path'
 nxs = 3
 gxs =100
 xs_label(1) = 'X'
 xs_label(2) = 'G'
 xs_label(3) = 'M'
 xs(1,:) = 0.5  , 0 , 0
 xs(2,:) = 0.0 , 0 , 0
 xs(3,:) = 0.5 , 0.5 , 0
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
sed "s|TBKOSTER_ROOT_DIR|$TBKOSTER_ROOT_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt
# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

done

