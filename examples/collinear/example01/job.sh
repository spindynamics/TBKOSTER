#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use DyNaMol.x to calculate the total energy vs a of Ptfcc"

# set the needed environment variables
. ../../environment_variables

rm -f tempo tempo2


cat > Etot_vs_a.dat << EOF
@#  a   Etot
EOF

for a in 3.80 3.82 3.84 3.86 3.88 3.90 3.92 3.94 3.96 3.98 4.00; do

  echo "a= $a"

cat > in_master.txt<<EOF
&calculation
 processing = 'scf'
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
 u_lcn(1)=20
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/pt_par_fcc_bcc_sc_lda_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = 0.0 0.5 0.5
 v(2,:) = 0.5 0.0 0.5
 v(3,:) = 0.5 0.5 0.0
 /
&atom
 ns = 1
 na = 1
 ntag = 1
 tag(1) = 'Pt'
 stag(1)=1
 pbc = 5, 5, 5
 r_coord='direct'
 r(1,:) = 0.0, 0.0, 0.0
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
 verbose=.false.
 ni_max=200
 /
EOF

# Set DyNaMol root directory in in_master.txt
sed "s|DYNAMOL_ROOT_DIR|$DYNAMOL_ROOT_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run DyNaMol
$BIN_DIR/DyNaMol.x

cat > tempo << EOF
a= $a
EOF

cat tempo out_energy.txt>>tempo2


done
grep -e 'a=' -e 'en =' tempo2 | awk '/a/{a = $(NF)}/en/{print a, $(NF)}' >> Etot_vs_a.dat

rm -f tempo tempo2
