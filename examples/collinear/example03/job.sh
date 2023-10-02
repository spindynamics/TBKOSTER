#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate the band structure of Pt(111)"

# set the needed environment variables
. ../../environment_variables


a=2.77185858225127
rm -f out* band/out*

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
 ne = 1
 symbol(1) = 'Pt'
 q(1)   = 10.0
 q_d(1) = 8.8
 u_lcn(1)=20
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/pt_par_fcc_bcc_sc_lda_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = 1.0 0.0 0.0
 v(2,:) =-0.5 0.866025403784438 0
 v(3,:) = 0.0 0.0 10.0 
 /
&atom
 ns = 1
 na = 13
 ntag = 1
 tag(1) = 'Pt'
 stag(1)=13
 pbc = 5, 5, 0
 r_coord='cartesian'
 r(1,:) =     0.0000000     0.0000000	  0.0000000
 r(2,:) =     1.3859293    -0.8001667	 -2.2632132
 r(3,:) =    -1.3859293     0.8001667	 -4.5264263
 r(4,:) =     0.0000000     0.0000000	 -6.7896395
 r(5,:) =     1.3859293    -0.8001667	 -9.0528526
 r(6,:) =    -1.3859293     0.8001667	-11.3160658
 r(7,:) =     0.0000000     0.0000000	-13.5792789
 r(8,:) =     1.3859293    -0.8001667	-15.8424921
 r(9,:) =    -1.3859293     0.8001667	-18.1057053
 r(10,:) =    0.0000000     0.0000000	-20.3689175
 r(11,:) =    1.3859293    -0.8001667	-22.6321316
 r(12,:) =   -1.3859293     0.8001667	-24.8953457
 r(13,:) =    0.0000000     0.0000000	-27.1585579
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
 verbose=.false.
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
 nxs = 4
 gxs =100
 xs_label(1) = 'G'
 xs_label(2) = 'X'
 xs_label(3) = 'M'
 xs_label(4) = 'G'
 xs(1,:) = 0  , 0 , 0
 xs(2,:) = 0.5 , 0 , 0
 xs(3,:) = 0.666666666 , 0.33333333 , 0
 xs(4,:) = 0 ,  0 , 0
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


