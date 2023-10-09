#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate the total energy and magnetization vs a of Crbcc AF"
$ECHO "penalization technique"

# set the needed environment variables
. ../../environment_variables

rm -f tempo tempo2 tempo3 tempo4
rm -f out*
cat > Etot_vs_a.dat << EOF
@#  a(A)  Etot(eV)
EOF

cat > M_vs_a.dat << EOF
@#  a(A)  Mtot(muB)
EOF

#for a in  2.71  2.72  2.73  2.74  2.75  2.76  2.77  2.78 2.79 2.80 2.81 2.82 2.83  2.84 2.85 2.87 2.88 2.89 2.90 2.91 2.92 2.93 2.94 2.95 3.0 3.10 3.20 3.30 3.40; do


for a in  2.91  ; do


$ECHO "a= $a"



cat > in_master.txt<<EOF
&calculation
 processing='scf'
 /
&units
 energy = 'ev'
 length = 'ang'
 time = 'fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Cr'
 q(1)   = 6.0
 q_d(1) = 5.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 0.82
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/cr_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = 1.0	0.0  0.0
 v(2,:) = 0.0	1.0  0.0
 v(3,:) = 0.0	0.0  1.0
 /
&atom
 ns = 2
 na = 2
 ntag = 2
 stag(1) = 1
 stag(2) = 1
 tag(1) = 'Cr_up'
 tag(2) = 'Cr_dn'
 r(1,:) = 0.0, 0.0, 0.0
 r(2,:) = 0.5, 0.5, 0.5
 m_listing = 'by_tag'
 m_coord = 'spherical'
m(1,:) =  1.0, 0.0, 0.0
m(2,:) = -1.0, 0.0, 0.0
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
 degauss = 0.2
 /
&mixing
 alpha = 0.1
 /
&scf
 delta_en = 0.0001
 delta_q  = 0.0001
 verbose = .true.
 /
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

cat > tempo << EOF
a= $a
EOF

grep m_r_tot out_log.txt | tail -1 >tempo3

cat tempo out_energy.txt>>tempo2
cat tempo tempo3>>tempo4



done
grep -e 'a=' -e 'en =' tempo2 | awk '/a/{a = $(NF)}/en/{print a, $(NF)}' >> Etot_vs_a.dat
grep -e 'a=' -e 'm_r_tot =' tempo4 | awk '/a/{a = $(NF)}/m_r_tot/{print a, $(NF-1)}' >> M_vs_a.dat

rm -f tempo tempo2 tempo3 tempo4
