#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate band structure of simple s-band cubic model"

# set the needed environment variables
. ../../../environment_variables

rm -fr band dos scf *.txt *.dat
mkdir dos scf band

cat > mod.dat<<EOF
6
-1  0  0   1 1 1 1   -1.0
0  -1  0   1 1 1 1   -1.0
0   0 -1   1 1 1 1   -1.0
0   0  1   1 1 1 1   -1.0
0   1  0   1 1 1 1   -1.0
1   0  0   1 1 1 1   -1.0
EOF

cat > in_master.txt<<EOF
&calculation
processing = 'scf'
post_processing = 'band'
/
&units
 energy = 'ev'
 length = 'ang'
 time = 'fs'
 mass = 'hau'
 /
&element
 ne = 1
 symbol(1) = 'J1'
 no(1) = 1
 o(1,1:1) = 1
 q(1) = 1.0000000000000000
 q_s(1) = 1.0
 q_p(1) = 0.0
 q_d(1) = 0.0
 /
&element_tb
  tb_type='mod'
 /
&lattice
 v_factor=1
 v(1,:) = 1.0, 0.0, 0.0
 v(2,:) = 0.0, 1.0, 0.0
 v(3,:) = 0.0, 0.0, 1.0
 /
&atom
 ns = 1
 na = 1
 ntag = 1
 stag(1) = 1
 tag(1) = 'J'
 pbc = 1, 1, 1
 r_coord = 'cartesian'
 r(1,:) = 0.00000000000000, 0.00000000000000, 0.00000000000000
 /
&mesh
 type = 'mp'
 gx = 20, 20, 20
 dx = 0, 0, 0
 /
&hamiltonian_tb
 e_e_interaction = 'stoner'
 m_penalization = 'none'
 /
&energy
 smearing = 'mv'
 degauss = 0.05
 fixed_fermi_level = .false.
 /
&mixing
 type = 'broyden'
 alpha = 0.10000000000000001
 n_init = 1
 n_hist = 50
 /
&scf
 delta_en = 0.0001
 delta_q = 0.0001
 verbose = .true.
 ni_max = 200
 /
EOF

cat > band/in_band.txt<<EOF
&band
 na_band = 1
 ia_band(1) = 1
 /
EOF

cat > band/in_mesh.txt<<EOF
&mesh
 type = 'path'
 nxs = 4
 gxs =100
 xs_label(1) = 'G'
 xs_label(2) = 'X'
 xs_label(3) = 'M'
 xs_label(4) = 'G'
 xs(1,:) = 0  , 0 , 0
 xs(2,:) = 0.5 , 0 , 0
 xs(3,:) = 0.5 , 0.5 , 0
 xs(4,:) = 0 ,  0 , 0
/
EOF

# Run TBKOSTER
$BIN_DIR/TBKOSTER.x

# Run bands.x
$BIN_DIR/bands.x

# Plot the results
cat > band/band_weight.gnuplot<<EOF
set term png enh size 700,500
set out 'band/projbands.png'
set xrange [*:*] ; set yrange [*:*]
set grid xtics
stats 'band/band_weight.dat' u 1:2 nooutput
set xra [STATS_min_x:STATS_max_x]
set yra [STATS_min_y:STATS_max_y]
set xlabel "k"
set ylabel "E - E_F (eV)"
set xzeroaxis
set key opaque box width 1.0
set style fill solid noborder
radius(proj)=proj/200.
plot 'band/band_weight.dat' u 1:2 lc rgb "grey" ,'band/band_weight.dat' u 1:2:(radius(\$3)) w circles lc rgb "red" t "{total}"-
EOF

gnuplot band/band_weight.gnuplot
