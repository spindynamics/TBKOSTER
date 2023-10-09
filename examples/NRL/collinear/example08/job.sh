#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi


# set the needed environment variables
. ../../../environment_variables


$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate the magnetization of Rhfcc as a function of lattice constant"
$ECHO "magnetization appears for lattice parameters above 4.3.."


rm -f tempo tempo2 tempo3 tempo4 
rm -f in_charge.txt
rm -f out*
rm -f scf
mkdir scf
rm -f results
mkdir results

cat > results/Etot_vs_a_fcc.dat << EOF
@#  a(A)  Etot(eV)
EOF

cat > results/M_vs_a_fcc.dat << EOF
@#  a(A)  Mtot(muB)
EOF

for i in {36..70..1} ; do


a=$(echo "$i/10.0" | bc -l)
echo a=$(echo "$a" | bc -l)


cat > in_master.txt<<EOF
&calculation
 processing='scf'
 /
&units
 energy = 'ev'
 length = 'ang'
 time = 'fs'
 mass = 'g/mol'
 /
&element
 ne = 1
 symbol(1) = 'Rh'
 q(1)   = 9.0
 q_d(1) = 8.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 0.85
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/rh_par_fcc_bcc_sc_gga_fl'
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
 tag(1) = 'Rh'
 r(1,:) = 0.0, 0.0, 0.0
 m_listing = 'by_tag'
 m_coord = 'spherical'
m(1,:) =  2.0, 0.0, 0.0
 /
&mesh
 type = 'mp'
 gx = 20, 20 , 20
 dx = 0, 0, 0
 /
&hamiltonian_tb
 /
&energy
 smearing = 'mv'
 degauss = 0.1
 /
&mixing
 alpha = 0.05
  n_init = 1
  n_hist = 50
 /
&scf
 delta_en = 0.0001
 delta_q  = 0.0001
 verbose = .true.
 ni_max = 200
 /
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

#cp -f out_charge.txt  in_charge.txt

cat > tempo << EOF
a= $a
EOF

grep m_r_tot out_log.txt | tail -1 >tempo3

cat tempo out_energy.txt>>tempo2
cat tempo tempo3>>tempo4



done
grep -e 'a=' -e 'en =' tempo2 | awk '/a/{a = $(NF)}/en/{print a, $(NF)}' >> results/Etot_vs_a_fcc.dat
grep -e 'a=' -e 'm_r_tot =' tempo4 | awk '/a/{a = $(NF)}/m_r_tot/{print a, $(NF-1)}' >> results/M_vs_a_fcc.dat

rm -f tempo tempo2 tempo3 tempo4


$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate the magnetization of Rhbcc as a function of lattice constant"
$ECHO "magnetization appears for lattice parameters above 4.3.."


cat > results/Etot_vs_a_bcc.dat << EOF
@#  a(A)  Etot(eV)
EOF

cat > results/M_vs_a_bcc.dat << EOF
@#  a(A)  Mtot(muB)
EOF

for i in {29..70..1} ; do


a=$(echo "$i/10.0" | bc -l)
echo a=$(echo "$a" | bc -l)


cat > in_master.txt<<EOF
&calculation
 processing='scf'
 /
&units
 energy = 'ev'
 length = 'ang'
 time = 'fs'
 mass = 'g/mol'
 /
&element
 ne = 1
 symbol(1) = 'Rh'
 q(1)   = 9.0
 q_d(1) = 8.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 0.85
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/rh_par_fcc_bcc_sc_gga_fl'
 /
&lattice
 v_factor = $a
 v(1,:) = -0.5,  0.5,  0.5
 v(2,:) =  0.5, -0.5,  0.5
 v(3,:) =  0.5,  0.5, -0.5
 /
&atom
 ns = 2
 na = 1
 ntag = 1
 stag(1) = 1
 tag(1) = 'Rh'
 r(1,:) = 0.0, 0.0, 0.0
 m_listing = 'by_tag'
 m_coord = 'spherical'
m(1,:) =  2.0, 0.0, 0.0
 /
&mesh
 type = 'mp'
 gx = 20, 20 , 20
 dx = 0, 0, 0
 /
&hamiltonian_tb
 /
&energy
 smearing = 'mv'
 degauss = 0.1
 /
&mixing
 alpha = 0.05
  n_init = 1
  n_hist = 50
 /
&scf
 delta_en = 0.0001
 delta_q  = 0.0001
 verbose = .true.
 ni_max = 200
 /
EOF

# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 

#cp -f out_charge.txt  in_charge.txt

cat > tempo << EOF
a= $a
EOF

grep m_r_tot out_log.txt | tail -1 >tempo3

cat tempo out_energy.txt>>tempo2
cat tempo tempo3>>tempo4



done
grep -e 'a=' -e 'en =' tempo2 | awk '/a/{a = $(NF)}/en/{print a, $(NF)}' >> results/Etot_vs_a_bcc.dat
grep -e 'a=' -e 'm_r_tot =' tempo4 | awk '/a/{a = $(NF)}/m_r_tot/{print a, $(NF-1)}' >> results/M_vs_a_bcc.dat

rm -f tempo tempo2 tempo3 tempo4





$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate band structure of Rhfcc"


a=3.80

$ECHO "a= $a"


cat > in_master.txt<<EOF
&calculation
 processing='scf'
 post_processing='band'
 /
&units
 energy = 'ev'
 length = 'ang'
 time = 'fs'
 mass='hau'
 /
&element
 ne = 1
 symbol(1) = 'Rh'
 q(1)   = 9.0
 q_d(1) = 8.0
 u_lcn(1) = 20.0
 i_stoner_d(1) = 0.85
 xi_so_d(1) = 0.0
 /
&element_tb
 filename(1) = '$TBPARAM_DIR/rh_par_fcc_bcc_sc_gga_fl'
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
 tag(1) = 'Rh'
 r(1,:) = 0.0, 0.0, 0.0
 m_listing = 'by_tag'
 m_coord = 'spherical'
m(1,:) =  1.0, 0.0, 0.0
 /
&mesh
 type = 'mp'
 gx = 15, 15, 15
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
 delta_en = 0.0001
 delta_q  = 0.0001
 verbose = .false.
 ni_max = 500
 /
EOF

cat > band/in_mesh.txt<<EOF
&mesh
 type = 'path'
 gxs = 100
 nxs = 12
 xs_label(1) = 'G'
 xs_label(2) = 'X'
 xs_label(3) = 'W'
 xs_label(4) = 'K'
 xs_label(5) = 'G'
 xs_label(6) = 'L'
 xs_label(7) = 'U'
 xs_label(8) = 'W'
 xs_label(9) = 'L'
 xs_label(10) = 'K'
 xs_label(11) = 'U'
 xs_label(12) = 'X' 
 xs(1,:) =  0  ,  0.0 ,  0
 xs(2,:) =  0.5,  0.0  ,  0.5
 xs(3,:) =  0.5,  0.25 ,  0.75
 xs(4,:) =  0.365, 0.375  , 0.75
 xs(5,:) =  0  ,  0    ,  0
 xs(6,:) =  0.5, 0.5,  0.5
 xs(7,:) =  0.625  ,  0.25 ,  0.625
 xs(8,:) =  0.5,  0.25  ,  0.75
 xs(9,:) =  0.5,  0.5 ,  0.5
 xs(10,:) =  0.365, 0.375  , 0.75
 xs(11,:) =  0.625  ,  0.25 ,  0.625
 xs(12,:) =  0.5,  0.0  ,  0.5
 /
EOF

cat > band/in_energy.txt<<EOF
&energy
 smearing = 'mv'
 degauss = 0.2
 en_min = -10.0
 en_max =  10.0
 /
EOF

cat > band/in_dos.txt<<EOF
&dos
 nen=100
 na_dos=1
 ia= 1
en_min=-10
 en_max=10
 /
EOF




# Set TBKOSTER root directory in in_master.txt
sed "s|BIN_DIR|$BIN_DIR|g" in_master.txt >in_master2.txt
mv -f in_master2.txt in_master.txt


# Run TBKOSTER
$BIN_DIR/TBKOSTER.x 


