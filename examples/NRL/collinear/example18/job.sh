#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate the PDOS of Co-Pt L10"

# set the needed environment variables
. ../../../environment_variables

a=2.68378136158666993601
covera=1.35340237919105196170

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
 symbol(1) = 'Co'
 symbol(2) = 'Pt'
 q(1)  = 9.0
 q_d(1) = 8.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 1.1
 q(2)  = 10.0
 q_d(2) = 9.0
 u_lcn(2)= 20.0
 i_stoner_d(2) = 0.6
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/co_par_fcc_bcc_sc_gga_fl'
 filename(2) = '$TBPARAM_DIR/pt_par_fcc_bcc_sc_lda_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = 1.0 0.0 0.0
 v(2,:) = 0.0 1.0 0.0
 v(3,:) = 0.0 0.0 $covera 
 /
&atom
 ns = 2
 na = 2
 ntag = 2
 tag(1) = 'Co' 
 tag(2) = 'Pt'
 stag(1)=1
 stag(2)=1 
 pbc = 5, 5, 5
 r(1,:) =      0.0000000     0.0000000     0.0000000
 r(2,:) =      0.5000000     0.5000000     0.5000000
 m(1,:) =      1.0000000     0.0000000     0.0000000
 m(2,:) =      0.5000000     0.0000000     0.0000000
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
 nen=200
 na_dos=2
 ia= 1 , 2
 en_min=-10
 en_max=10
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

cat > dos/in_mesh.txt<<EOF
&mesh
 type = 'mp'
 gx = 15, 15, 15
 dx = 0, 0, 0
 /
EOF


# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 
