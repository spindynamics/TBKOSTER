#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to compute the forces on a Pt dimer"

# set the needed environment variables
. ../../environment_variables

rm -f tempo tempo2
rm -f Etot_vs_force*

cat > Etot_vs_forces.dat << EOF
@# a f Etot
EOF

for a in 1.90 1.92 1.94 1.96 1.98 2.00 2.02 2.04 2.06 2.08 2.10 ; do


  echo "a= ${a}"


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
 v_factor =${a}
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
 pbc = 0, 0, 0
 r_coord='direct'
 r(1,:) =      0.0000000     0.0000000     0.0000000
 r(2,:) =      0.0000000     0.0000000     1.0000000
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

# Set TBKOSTER root directory in in_master.txt
sed "s|TBPARAM_DIR|$TBPARAM_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt

# Run TBKOSTER
$BIN_DIR/TBKOSTER.x

cat > alat << EOF
a= $a
EOF

grep Fz ./scf/out_log.txt | awk '{print $8}' | head -1 > fz

grep -e 'a=' alat | awk '{print $2}' > a

grep -e 'en =' out_energy.txt | awk '{print $3}' > en

paste a fz en >> Etot_vs_forces.dat

done

rm -f fz en a

