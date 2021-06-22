import f90nml
import matplotlib.pyplot as pyplot
import numpy
import os.path

# Read simulation data from txt files
if os.path.isfile('in_master.txt'):
    in_master_nml = f90nml.read('in_master.txt')
#else:
    # Read element

    # Read lattice

    # Read atom

    # Read mesh

    # Read hamtilonian_tb

    # Read energy

    # Read self_consistent_field

    # Read spin_dynamics

    # Read calculation

ns = in_master_nml['atom']['ns']
if ns == 1 or ns == 4:
    nsl = 1
else:
    nsl = 2

# Read units from out_units.txt
if os.path.isfile('out_units.txt'):
    out_units_nml = f90nml.read('out_units.txt')
else:
    print('plot_band(): out_units.txt not found')
units_energy = out_units_nml['units']['energy']
if units_energy == 'hau':
    units_energy = 'H.a.u.'
elif units_energy == 'rau':
    units_energy = 'R.a.u.'
elif units_energy == 'ev':
    units_energy = 'eV'

# Read Hamiltonian dimension from out_hamiltonian_tb.txt
if os.path.isfile('out_hamiltonian_tb.txt'):
    out_hamiltonian_tb_nml = f90nml.read('out_hamiltonian_tb.txt')
else:
    print('plot_band(): out_hamiltonian_tb.txt not found')
nh = out_hamiltonian_tb_nml['hamiltonian_tb']['nh']

# Read Fermi level from out_energy.txt
if os.path.isfile('out_energy.txt'):
    out_energy_nml = f90nml.read('out_energy.txt')
else:
    print('plot_band(): out_energy.txt not found')
en_f = out_energy_nml['energy']['en_f']

# Read high symmetry points from band/in_mesh.txt
if os.path.isfile('band/in_mesh.txt'):
    band_in_mesh_nml = f90nml.read('band/in_mesh.txt')
else:
    print('plot_band(): band/in_mesh.txt not found')
gxs = band_in_mesh_nml['mesh']['gxs']
nxs = band_in_mesh_nml['mesh']['nxs']
xs_label = band_in_mesh_nml['mesh']['xs_label']
xs_label = [xs_str.replace('G','$\Gamma$') for xs_str in xs_label]

# Read points from band/out_mesh.txt
if os.path.isfile('band/out_mesh.txt'):
    band_out_mesh_nml = f90nml.read('band/out_mesh.txt')
else:
    print('plot_band(): band/out_mesh.txt not found')
nx = band_out_mesh_nml['mesh']['nx']
x = numpy.asarray(band_out_mesh_nml['mesh']['x'])

# Calculate the cumulative distance between mesh points
xd = numpy.linalg.norm(x[:,1:nx]-x[:,0:nx-1],axis=0)
xd = numpy.insert(xd,0,0.0)
xd = numpy.cumsum(xd)

# Read energy bounds and energy eigenvalues from band/out_energy.txt
if os.path.isfile('band/out_energy.txt'):
    band_out_energy_nml = f90nml.read('band/out_energy.txt')
else:
    print('plot_band(): band/out_energy.txt not found')
en_min = band_out_energy_nml['energy']['en_min']
en_max = band_out_energy_nml['energy']['en_max']
en_k = numpy.asarray(band_out_energy_nml['energy']['en_k']) - en_f

print('test=',en_k)
# Plot
for isl in range(0,2):
    pyplot.figure()
    # Energy eigenvalues
    for ih in range(0,nh):
        pyplot.plot(xd,en_k[isl,:,ih],'k-')
    # Fermi level
    pyplot.hlines(0.0,xd[0],xd[-1],linestyles='dashed')
    # High symmetry points
    for ixs in range(1,nxs):
        pyplot.vlines(xd[ixs*gxs],en_min,en_max)
    # Axes properties
    pyplot.xlim(xd[0],xd[-1])
    pyplot.xticks(xd[range(0,nx,gxs)], xs_label)
    pyplot.ylabel('E (' + units_energy + ')')
    pyplot.ylim(en_min,en_max)
    # Save figure in .eps and .png format
    pyplot.savefig('band/band_' + str(isl+1) + '.eps', format='eps')
    pyplot.savefig('band/band_' + str(isl+1) + '.png', format='png')
