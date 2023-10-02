# Collinear spin examples (ns=1 and 2)

Two types of constraints:

m_penalization='none'   no constraint
m_penalization='r'      constraint on the magnetization vector (M)


* **example01**

Calculation of energy vs lattice parameter curve of fcc Pt

* **example02** 

Self consistent calculation of the (111) surface of Pt.
projected DOS calculation for the surface atom

* **example03**

Self consistent calculation of the (111) surface of Pt.
band structure calculation

* **example04**

Self consistent calculation of magnetic bcc Fe at several lattice parameters.
then as an output we can plot Etot_vs_a.dat M_vs_a.dat

* **example05**

Self consistent calculation of bcc Fe for fixed-spin-moment (fixed_spin_moment=T) at various values of M to get E(M)

* **example06**

Comparison with
Self consistent calculation of Febcc with penalization on the total moment of the type lambda(M-Mo)^2 
Here the bcc is described by a cubic lattice with 2 atoms per unit cell and the penalization is only on atom 1.
Atom is free to adapt..

* **example07**

Self consistent calculation of iron cluster (cuboctahedron)
verbose=.true. and the relaxation process is saved in scf->vizualization with ovito

* **example08**

Magnetization of Rh as a function of the lattice parameter..
and  band structure of Rh.

* **example09**

Study of the influence of the number of k points on the total energy
AND
Study of the influence of the smearing technique on the total energy

* **example10**

Self consistent calculation of a 5 atom Fe wire with penalization on atom 1 m1=3 m2,m3,m4,m5 free...

* **example11**

Chromium bcc: magnetism as a function of lattice parameter.

* **example12**

Self consistent calculation of Cofcc at several lattice parameters

* **example13**

Apparition of magnetic moment with expansion of the lattice parameter for Ptfcc

* example14

Au(111) fcc surface scf calculation.
projected DOS calculation for the surface atom

* **example15**

Au(111) fcc surface band structure calculation. To accelerate the scf calculation one can cp the out_charge.txt of example14 into in_charge.txt of example15 so that the scf calculation will restart from a converged charge. The band structure calculation is performed. Remakable agreement with pwscf calculation
(Surf. Sci. 602, (2008) 893, Fig. 2)

* **example16**

Ni fcc magnetic with respect to lattice parameter: study of the influence of the Stoner parameter->comparison with _ab-initio_

* **example17**

Ni monolayer on Au fcc(111) interface, magnetic calculation

* **example18**

Co-Pt L10   scf +DOS

* **example19**

Fe-Pt L10
Etot(c/a) for FM and AF

* **example25**
This example shows how to use TBKOSTER.x to calculate the total energy of zirconium fcc, bcc and hcp(c/a)

* **example26**
This example shows how to use TBKOSTER.x to explore the magnetic configuration FM and AFM (+ continuous transition between the two) of B2 FeRh. For lattice parameters below a=3A the system is AFM while it becomes FM above a=3A. This is in perfect agreement with DFT calculations.

# Non collinear spin examples (ns=4)

Different types of constraints:

magnetic_penalization='none'        no constraint
magnetic_penalization='r'           constraint on the magnetization component  (r,  .  , . )
magnetic_penalization='r,theta'     constraint on the magnetization components (r,theta, . )
magnetic_penalization='r,theta,phi' constraint on the magnetization vector     (r,theta,phi)
magnetic_penalization='theta'       constraint on the magnetization component  (.,theta, . )
magnetic_penalization='theta,phi'   constraint on the magnetization components (.,theta,phi)
magnetic_penalization='phi'         constraint on the magnetization component  (.,  .  ,phi)

* **example01**

Band structure of Cu fcc with SOC

* **example02**

SCF non-collinear spin calculation of a 5-atom Fe wire with magnetic penalization on atom 1. The penalization on atom 1 is theta,phi=(3,30,0) and the other atoms are free. During scf the spin tend do align.


* **example03**

Same system as example01 but with no penalization. However the initial magnetization of atom 1 is opposite to the other atoms. For symmetry reasons the magnetization remains opposite which simulate the simplest magnetic excitation.

* **example04**

SCF NON-collinear spin calculation of a 4-atom Fe wire
The initial magnetization is a spin spiral of period 4a 
Due to periodic boundary conditions the spin spiral configuration is kept during scf

* **example05**

Magnetic anisotropy, spin and orbital moment of Fe wire for several TB models (Stoner UJ, UJB).

* **example06**

band-structure of Fe wire with SOC. Magnetization along the wire and perpendicular to the wire.

* **example07 **

Same as example 13 but with SOC. Apparition of magnetic moment with expansion of the lattice parameter for Ptfcc with and without spin-orbit coupling. Strong influence of SOC.

* **example08**

Anisotropy of of Fe monolayer with several magnetic interaction: Stoner & UJB

* **example09**

same as example15 but with SOC
Au(111) scf calculation with SOC (xi_p=1.0eV xi_d=0.65eV)
 +
Au(111) band structure calculation (Rashba splitting of shockley surface state).

* **example10**

same as example19 but with SOC
Fe-Pt L10 system: magnetic anisotropy of FM solution versus c/a at fixed volume.

* **example11**

Fe monatomic wire: calculation of spin-spiral configuration by two different ways:
- standard calculation with 4 atoms per unit cell
- spin spiral calculation with one atom per unit-cell and  k_spiral=(0,0,1/4)
- a FM calculation is also performed... it is the most stable.

* **example12**

Cr bcc: calculation of AF configuration by two different ways:
- standard collinear calculation with two atoms per unit cell
- spin spiral calculation with one atom per unit-cell and  k_spiral=(0,0,1/2) (xyz)=(-1/2,1/2,1/2) (direct)

* **example13**

Fe monatomic wire: spin spiral calculation of E(q) there is a shallow minimum for a given k_spiral this minimum disappears for larger lattice parameters.

* **example14**

Fe monatomic wire: calculation of spin spiral by two different ways:
- non-collinear calculation with 18 atoms per unit cell and theta_pen(i)=20*(i-1) i=1,18
- spin spiral calculation with one atom per unit-cell and  Q=(0,0,1/18)

* **example15**
DOS of Febcc with SOC
band structure of Febcc with SOC

* **example16**
SCF calculation of Cr trimer. Neel configuration is found (depending on the input magnetism)
verbose=.true. and the relaxation process is saved in scf->vizualization with ovito

* **example17**
This example shows how to use TBKOSTER.x to calculate a Fe cuboctahedron cluster
verbose=.true. and the relaxation process is saved in scf->vizualization with ovito"


# Atomic forces & molecular dynamics 

* **example01**
Basic test for forces

* **example02**
Testing the forces by computing finite differences in energy.

* **example03**
Testing the forces on a Pt trimer

* **example04**
Testing the forces by computing finite differences in energy for a 3D PBC system: moving the central atom bcc system. 

* **example05**
Testing the forces by computing finite differences in energy for a 2D PBC system: moving the surface layer of an Au(111) slab.

* **example06**
Test of the forces with respect to the finite difference for an Cr trimer (non coollinear)

* **example07**
Molecular dynamics of a Co trimer

# Spin dynamics examples

* **example01**

Single spin dynamics of a 5 atoms Fe chain.

* **example02**

Damped Spin Dynamics of a 5 atoms Fe chain with varying lattice parameter.

* **example03**

Damped Spin Dynamics of a Cr trimer-> Converge towards Néel structure.

* **example04**

Damped Spin Dynamics of a Fe trimer-> Converge towards FM structure.

* **example05**


