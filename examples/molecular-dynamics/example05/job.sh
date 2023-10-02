#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate the forces for 4 slabs of Au(111)"


# set the needed environment variables
. ../../environment_variables

a=2.88499566724111
rm -f out*
mkdir scf
rm -f out* ./scf/out*


cat > Etot_vs_forces.dat << EOF
@# a f Etot
EOF


for dx in -0.1 -0.09 -0.08 -0.07 -0.06 -0.05 -0.04 -0.03 -0.0288 -0.02304 -0.01728 -0.01152 -0.00576 0 0.00576 0.01152 0.01728 0.02304 0.0288 0.03 0.04 0.05; do

  echo "dx= ${dx}"

cat > in_master.txt<<EOF
&calculation
 processing = 'scf'
 post_processing='forces'
 post_processing_dir = 'scf'
 /
&units
 length='ang'
 energy='ev'
 time='fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Au'
 q(1)   = 11.0
 q_d(1) = 10
 u_lcn(1)=20
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/au_par_fcc_bcc_sc_lda_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = 1.0 0.0 0.0
 v(2,:) =-0.5 0.866025403784438 0
 v(3,:) = 0.0 0.0 0.1061445555206044D+02 
 /
&atom
 ns = 1
 na = 4 
 ntag = 1
 tag(1) = 'Au'
 stag(1)=4 
 pbc = 5, 5, 0 
 r_coord='cartesian'
 r(1,:) =      0.0000000     0.0000000     $dx
 r(2,:) =      1.4424978    -0.8328265    -2.3555892
 r(3,:) =     -1.4424978     0.8328265    -4.7111783
 r(4,:) =      0.0000000     0.0000000    -7.0667677
 /
&mesh
 type = 'mp'
 gx = 10 , 10 , 1
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
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x


cat > alat << EOF
a= $(echo "$dx" | bc -l)
EOF

grep "Atom 1" ./scf/out_log.txt | awk '{print $8}' | head -1 > fz

grep -e 'a=' alat | awk '{print $2}' > a

grep -e 'en =' out_energy.txt | awk '{print $3}' > en

paste a fz en >> Etot_vs_forces.dat

cp ./scf/out_log.txt ./scf/out_log$(echo "$dx" | bc -l).txt
done

rm -f alat a en fz
 
