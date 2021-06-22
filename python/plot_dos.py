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

# Read Fermi level from out_energy.txt
if os.path.isfile('out_energy.txt'):
    out_energy_nml = f90nml.read('out_energy.txt')
else:
    print('plot_band(): out_energy.txt not found')
en_f = out_energy_nml['energy']['en_f']

# Read energy bounds from dos/in_dos.txt
if os.path.isfile('dos/in_dos.txt'):
    dos_out_dos_nml = f90nml.read('dos/out_dos.txt')
else:
    print('plot_band(): dos/out_dos.txt not found')
nen = dos_out_dos_nml['dos']['nen']
en_min = dos_out_dos_nml['dos']['en_min']
en_max = dos_out_dos_nml['dos']['en_max']
den = en_max-en_min/(nen+1)
na_dos = dos_out_dos_nml['dos']['na_dos']
dos = numpy.asarray(dos_out_dos_nml['dos']['dos'])
if na_dos > 0:
    ia = numpy.asarray(dos_out_dos_nml['dos']['ia'])
    dos_s = numpy.asarray(dos_out_dos_nml['dos']['dos_s'])
    dos_p  = numpy.asarray(dos_out_dos_nml['dos']['dos_p'])
    #dos_px = numpy.asarray(dos_out_dos_nml['dos']['dos_px'])
    #dos_py = numpy.asarray(dos_out_dos_nml['dos']['dos_py'])
    #dos_pz = numpy.asarray(dos_out_dos_nml['dos']['dos_pz'])
    dos_d     = numpy.asarray(dos_out_dos_nml['dos']['dos_d'])
    #dos_dxy   = numpy.asarray(os_out_dos_nml['dos']['dos_dxy'])
    #dos_dyz   = numpy.asarray(dos_out_dos_nml['dos']['dos_dyz'])
    #dos_dzx   = numpy.asarray(dos_out_dos_nml['dos']['dos_dzx'])
    #dos_dx2y2 = numpy.asarray(dos_out_dos_nml['dos']['dos_dx2y2'])
    #dos_dz2r2 = numpy.asarray(dos_out_dos_nml['dos']['dos_dz2r2'])
    if ns == 4:
        dos_s[1,:,:] = dos_s[1,:,:] + dos_s[2,:,:]
        dos_p[1,:,:] = dos_p[1,:,:] + dos_p[2,:,:]
        dos_d[1,:,:] = dos_d[1,:,:] + dos_d[2,:,:]
    dos_tot = dos_s + dos_p + dos_d

# Plot
en = numpy.linspace(en_min,en_max,num=nen)
dos_max = numpy.amax(dos)
for isl in range(0,2):
    pyplot.figure()
    # DOS
    pyplot.plot(en,dos[isl,:],'k-')
    # Fermi level
    pyplot.vlines(0.0,0.0,dos_max,linestyles='dashed')
    # Axes properties
    pyplot.xlabel('E (' + units_energy + ')')
    pyplot.xlim(en_min,en_max)
    pyplot.ylabel('DOS')
    pyplot.ylim(0.0,dos_max)
    # Save figure in .eps and .png format
    pyplot.savefig('dos/dos_isl' + str(isl+1) + '.eps', format='eps')
    pyplot.savefig('dos/dos_isl' + str(isl+1) + '.png', format='png')

    if na_dos > 0:
        for ia_dos in range(0,na_dos):
            pyplot.figure()
            # DOS
            pyplot.plot(en,dos_s[isl,:,ia_dos],'r-')
            pyplot.plot(en,dos_p[isl,:,ia_dos],'g-')
            pyplot.plot(en,dos_d[isl,:,ia_dos],'b-')
            pyplot.plot(en,dos_tot[isl,:,ia_dos],'k-')
            # Fermi level
            pyplot.vlines(0.0,0.0,dos_max,linestyles='dashed')
            # Axes properties
            pyplot.legend(['s','p','d','total'])
            pyplot.xlabel('E (' + units_energy + ')')
            pyplot.xlim(en_min,en_max)
            pyplot.ylabel('DOS')
            pyplot.ylim(0.0,dos_max)
            # Save figure in .eps and .png format
            pyplot.savefig('dos/dos_ia' + str(ia[ia_dos]) \
             + '_isl' + str(isl+1) + '.eps', format='eps')
            pyplot.savefig('dos/dos_ia' + str(ia[ia_dos]) \
             + '_isl' + str(isl+1) + '.png', format='png')
