#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate the PDOS of Ni monolayer of Au(111)"

# set the needed environment variables
. ../../../environment_variables

a=2.88499566724111
rm -f out*
mkdir scf

cat > in_master.txt<<EOF
&calculation
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
 ne = 2
 symbol(1) = 'Ni'
 symbol(2) = 'Au'
 q(1)  = 10.0
 q_d(1) = 9.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 1.1
 q(2)  = 11.0
 q_d(2) = 10.0
 u_lcn(2)= 20.0
 i_stoner_d(2) = 0.8 
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/ni_par_fcc_bcc_sc_gga_fl'
 filename(2) = '$TBPARAM_DIR/au_par_fcc_bcc_sc_lda_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = 1.0 0.0 0.0
 v(2,:) =-0.5 0.866025403784438 0
 v(3,:) = 0.0 0.0 0.1061445555206044D+02 
 /
&atom
 ns = 2
 na = 15
 ntag = 3
 tag(1) = 'Au' 
 tag(2) = 'Ni'
 tag(3) = 'Au_surf'
 stag(1)=1
 stag(2)=13 
 stag(3)=1
 pbc = 5, 5, 0
 r_coord='cartesian'
 r(1,:) =      0.0000000     0.0000000     0.0000000
 r(2,:) =      1.4424978    -0.8328265    -2.3555892
 r(3,:) =     -1.4424978     0.8328265    -4.7111783
 r(4,:) =      0.0000000     0.0000000    -7.0667677
 r(5,:) =      1.4424978    -0.8328265    -9.4223566
 r(6,:) =     -1.4424978     0.8328265   -11.7779455
 r(7,:) =      0.0000000     0.0000000   -14.1335354
 r(8,:) =      1.4424978    -0.8328265   -16.4891243
 r(9,:) =     -1.4424978     0.8328265   -18.8447132
 r(10,:) =    -0.0000000     0.0000000   -21.2003021
 r(11,:) =     1.4424978    -0.8328265   -23.5558910
 r(12,:) =    -1.4424978     0.8328265   -25.9114799
 r(13,:) =     0.0000000     0.0000000   -28.2670708
 r(14,:) =     1.4424978    -0.8328265   -30.6226578
 r(15,:) =    -1.4424978     0.8328265   -32.9782486
 m(1,:)  =     1.0           0            0
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

cat > dos/in_dos.txt<<EOF
&dos
 nen=100
 na_dos=1
 ia= 1
 en_min=-10
 en_max=10
 /
EOF

cat > dos/in_energy.txt<<EOF
&energy
 smearing = 'mv'
 degauss = 0.2
 en_min = -10.0
 en_max =  10.0
 /
EOF

cat > dos/in_mesh.txt<<EOF
&mesh
 type = 'mp'
 gx = 15, 15, 1
 dx = 0, 0, 0
 /
EOF


# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 
