#!/usr/bin/env bash
# run from directory where this script is
cd `echo $0 | sed 's/\(.*\)\/.*/\1/'` # extract pathname
EXAMPLE_DIR=`pwd`

# check whether echo has the -e option
if test "`echo -e`" = "-e" ; then ECHO=echo ; else ECHO="echo -e" ; fi

$ECHO
$ECHO "$EXAMPLE_DIR : starting"
$ECHO
$ECHO "This example shows how to use TBKOSTER.x to calculate the total energy vs a of Ptfcc"

# set the needed environment variables
. ../../../environment_variables

# remove existing stuffs
rm -f *.dat *.txt *.gnuplot *.png

cat > Etot_vs_a.dat <<EOF
#  a   Etot
EOF

for a in 3.80 3.82 3.84 3.86 3.88 3.90 3.92 3.94 3.96 3.98 4.00; do
    rm -f tempo tempo2
    cat > in_master.txt<<EOF
    &calculation
    processing = 'scf'
    /
    &units
    length = 'ang'
    energy = 'ev'
    time = 'fs'
    mass = 'hau'
    /
    &element
    ne = 1
    symbol(1) = 'Pt'
    u_lcn(1) = 20
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
    stag(1) = 1
    pbc = 5, 5, 5
    r_coord = 'direct'
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
    delta_en = 0.0001
    delta_q = 0.0001
    verbose = .false.
    ni_max = 200
    /
EOF

    # Run TBKOSTER
    $BIN_DIR/TBKOSTER.x

    ENERGY=$(grep "en =" out_energy.txt | cut -d "=" -f2)
    echo $a $ENERGY >> Etot_vs_a.dat
done

# Display the results
if ! command -v gnuplot &> /dev/null
then 
    $ECHO "The gnuplot command cound not be found. Please install gnuplot."
    exit 1
else
    cat > Etot_vs_a.gnuplot<<EOF
	set encoding utf8
	set xlabel "Lattice Parameter (ang)"
	set ylabel "Energy (eV)"
	set terminal png
	set output "Etot_vs_a.png"
	plot 'Etot_vs_a.dat' w lines lw 2
EOF
    gnuplot Etot_vs_a.gnuplot
fi
