#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to compute the forces"

# set the needed environment variables
. ../../../environment_variables

# checking scf folder existence
if [ ! -d "./scf" ]; then 
    mkdir -p ./scf
fi

rm -f results.dat

print_forces(){
    Fx1=$(grep "Atom 1" ./scf/out_log.txt | cut -d"=" -f2 | cut -d"F" -f1)
    Fy1=$(grep "Atom 1" ./scf/out_log.txt | cut -d"=" -f3 | cut -d"F" -f1)
    Fz1=$(grep "Atom 1" ./scf/out_log.txt | cut -d"=" -f4)
    Fx2=$(grep "Atom 2" ./scf/out_log.txt | cut -d"=" -f2 | cut -d"F" -f1)
    Fy2=$(grep "Atom 2" ./scf/out_log.txt | cut -d"=" -f3 | cut -d"F" -f1)
    Fz2=$(grep "Atom 2" ./scf/out_log.txt | cut -d"=" -f4)
    en=$(grep "en =" ./out_energy.txt | cut -d"=" -f2)
    x=$(echo "sqrt(3)*${z}*${v_factor}" | bc -l)
    echo ${k} ${x} ${en} ${Fx1} ${Fy1} ${Fz1} ${Fx2} ${Fy2} ${Fz2} >> results.dat
}

v_factor=6.0
epsilon=0.001
evaluated_z=0.45

for k in 10; do

for i in 0 1 2 3 4 5 6 7 8 9 10; do

z=$(echo "${evaluated_z}*(1.0+${epsilon}*(${i}-5))" | bc -l)

echo "Compute Energy for z=${z}"

cat > in_master.txt<<EOF
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
 v_factor = ${v_factor}
 v(1,:) = 1 0 0
 v(2,:) = 0 1 0
 v(3,:) = 0 0 1 
 /
&atom
 ns = 1
 na = 2
 ntag = 1
 tag(1) = 'Pt'
 stag(1)= 2
 pbc = 5, 5, 5
 r(1,:) =      0.00     0.00     0.00
 r(2,:) =      0.5     0.5     ${z}
 /
&mesh
 type = 'mp'
 gx = ${k}, ${k}, ${k}
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

# Set TBKOSTER root directory in in_master.txt
sed "s|TBPARAM_DIR|$TBPARAM_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt

# Run TBKOSTER
$BIN_DIR/TBKOSTER.x

print_forces

done
done
