#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use DyNaMol.x to calculate the total energy and magnetization vs a of Nifcc"
$ECHO "the calculation are made at various Stoner parameters"


# set the needed environment variables
. ../../environment_variables

rm -f out*

mkdir results

for IStoner in 0.90 0.95 1.00 1.05 1.10 1.15 1.20 ; do
echo "IStoner= $IStoner"

rm -f tempo tempo2 tempo3 tempo4

cat > results/Etot_vs_a.I$IStoner.dat << EOF
@#  a(A)  Etot(eV)
EOF

cat > results/M_vs_a.I$IStoner.dat << EOF
@#  a(A)  Mtot(muB)
EOF

for a in 3.20  3.30 3.40  3.45 3.50 3.55 3.60 3.65 3.70 3.75 3.80 3.85 3.90 3.95 4.0 4.50 5.00; do

echo "a= $a"


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
 symbol(1) = 'Ni'
 q(1)   = 10.0
 q_d(1) = 9.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = $IStoner
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/ni_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = 0.0 0.5 0.5
 v(2,:) = 0.5 0.0 0.5
 v(3,:) = 0.5 0.5 0.0
 /
&atom
 ns = 2
 na = 1
 ntag = 1
 stag(1) = 1
 tag(1) = 'Ni'
 r(1,:) = 0.0, 0.0, 0.0
 m_listing = 'by_tag'
 m_coord = 'spherical'
m(1,:) =  1.0, 0.0, 0.0
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
 ni_max = 100
 /
EOF

# Set DyNaMol root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run DyNaMol
$BIN_DIR/DyNaMol.x 

cat > tempo << EOF
a= $a
EOF

grep m_r_tot out_log.txt | tail -1 >tempo3

cat tempo out_energy.txt>>tempo2
cat tempo tempo3>>tempo4


done
grep -e 'a=' -e 'en =' tempo2 | awk '/a/{a = $(NF)}/en/{print a, $(NF)}' >> results/Etot_vs_a.I$IStoner.dat
grep -e 'a=' -e 'm_r_tot =' tempo4 | awk '/a/{a = $(NF)}/m_r_tot/{print a, $(NF-1)}' >> results/M_vs_a.I$IStoner.dat

rm -f tempo tempo2 tempo3 tempo4
done
