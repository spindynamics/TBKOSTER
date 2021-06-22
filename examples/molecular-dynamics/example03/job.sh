#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use DyNaMol.x to compute the forces on a trimer"

# checking scf folder existence
if [ ! -d "./scf" ]; then 
mkdir -p ./scf
fi

# removing things before computation
rm -f out_* in_master* results.dat

# set the needed environment variables
. ../../environment_variables

y=$(echo "sqrt(3.0)" | bc -l)

cat > in_master2.txt<<EOF
&calculation
 processing = 'scf'
 post_processing= 'forces'
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
 symbol(1) = 'Pt'
 q(1)   = 10.0
 q_d(1) = 8.8
 u_lcn(1)=20
 /
&element_tb
 filename(1) = 'TBPARAM_DIR/pt_par_fcc_bcc_sc_lda_fl'
 /
&lattice
 v_factor = 1.0
 v(1,:) = 1 0 0
 v(2,:) = 0 1 0
 v(3,:) = 0 0 1 
 /
&atom
 ns = 1
 na = 3
 ntag = 1
 tag(1) = 'Pt'
 stag(1)= 3
 pbc = 0, 0, 0
 r_coord='cartesian'
 r(1,:) = 0.0 0.0 0.0
 r(2,:) = 0.0 ${y} 1.0
 r(3,:) = 0.0 0.0 2.0
 /
&mesh
 type = 'mp'
 gx = 1, 1, 1
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
&forces
 computed=.true.
 /
EOF

# Set DyNaMol root directory in in_master.txt
sed "s|TBPARAM_DIR|$TBPARAM_DIR|g" in_master2.txt >in_master.txt
rm -f in_master2.txt

# Run DyNaMol
$BIN_DIR/DyNaMol.x

en=$(grep "en =" ./out_energy.txt | cut -d"=" -f2)

Fx1=$(grep "Atom 1" ./scf/out_log.txt | cut -d"=" -f2 | cut -d"F" -f1)
Fy1=$(grep "Atom 1" ./scf/out_log.txt | cut -d"=" -f3 | cut -d"F" -f1)
Fz1=$(grep "Atom 1" ./scf/out_log.txt | cut -d"=" -f4)
F1=$(echo "scale=4; sqrt(${Fx1}*${Fx1}+${Fy1}*${Fy1}+${Fz1}*${Fz1})" | bc -l)
Fx2=$(grep "Atom 2" ./scf/out_log.txt | cut -d"=" -f2 | cut -d"F" -f1)
Fy2=$(grep "Atom 2" ./scf/out_log.txt | cut -d"=" -f3 | cut -d"F" -f1)
Fz2=$(grep "Atom 2" ./scf/out_log.txt | cut -d"=" -f4)
F2=$(echo "scale=4; sqrt(${Fx2}*${Fx2}+${Fy2}*${Fy2}+${Fz2}*${Fz2})" | bc -l)
Fx3=$(grep "Atom 3" ./scf/out_log.txt | cut -d"=" -f2 | cut -d"F" -f1)
Fy3=$(grep "Atom 3" ./scf/out_log.txt | cut -d"=" -f3 | cut -d"F" -f1)
Fz3=$(grep "Atom 3" ./scf/out_log.txt | cut -d"=" -f4)
F2=$(echo "scale=4; sqrt(${Fx2}*${Fx2}+${Fy2}*${Fy2}+${Fz2}*${Fz2})" | bc -l)

echo ${z} ${en} ${Fx1} ${Fy1} ${Fz1} ${Fx2} ${Fy2} ${Fz2} ${Fx3} ${Fy3} ${Fz3} >> results.dat
