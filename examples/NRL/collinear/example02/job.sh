#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate PDOS of a Pt(111) surface"

# set the needed environment variables
. ../../environment_variables

a=2.77185858225127
rm -fr xyz dos scf *.txt
mkdir xyz dos scf

cat > in_master.txt<<EOF
&calculation
 pre_processing='txt2xyz'
 pre_processing_dir='xyz'
 processing = 'scf'
 post_processing='dos'
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
 gx = 15, 15, 1
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

cat > dos/in_energy.txt<<EOF
&energy
 smearing = 'g'
 degauss = 0.1
 en_min = -10.0
 en_max =  10.0
 /
EOF

cat > dos/in_dos.txt<<EOF
&dos
 nen=500
 na_dos=2
 ia(1)= 1
 ia(2)=7
 en_min=-10
 en_max=10
 /
EOF

cat > dos/in_mesh.txt<<EOF
&mesh
 type = 'mp'
 gx = 15, 15, 1
 dx = 0, 0, 0
 /
EOF

# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 
