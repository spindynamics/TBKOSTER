!
! Copyright (C) 2017
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
! Mathieu Cesar <mailto:mathieu.cesar@cea.fr>,
! Pascal Thibaudeau <mailto:pascal.thibaudeau@cea.fr>.
!
! This software is a computer program whose purpose is DyNaMol.
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and inRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
!
!  constant.f90
!  DyNaMol
module constant_mod
  use math_mod, only: pi
  use precision_mod, only: rp
  implicit none

  !> @defgroup Fundamental_constants Fundamental constants
  !> @{

  !> Vacuum permittivity /$ \epsilon_0 (\mathrm{F m}^{-1}]) /$ from
  !> https://en.wikipedia.org/wiki/Vacuum_permittivity
  real(rp),parameter :: epsilon_0 = 8.8541878128e-12_rp
  !> Electron spin g-factor /$ g_{\mathrm{e}} /$ from
  !> https://en.wikipedia.org/wiki/G-factor_(physics)#Electron_g-factors
  real(rp),parameter :: g_e = 2.002319304
  !> Planck constant \f$ \hbar (\mathrm{eV.fs}) \f$ from
  !> https://en.wikipedia.org/wiki/Planck_constant
  real(rp),parameter :: hbar = 6.582119569e-1_rp
  !> Boltzmann constant \f$ k_{\mathrm{B}} (\mathrm{eV.K}^{-1}) \f$ from
  !> https://en.wikipedia.org/wiki/Boltzmann_constant
  real(rp),parameter :: k_b = 8.617333262145e-5_rp
  !> Electron rest mass /$ m_{\mathrm{e}} (\mathrm{Da}) /$ from
  !> https://en.wikipedia.org/wiki/Electron_rest_mass
  real(rp),parameter :: m_e = 5.48579909065e-4_rp
  !> Elementary charge /$ m_{\mathrm{e}} (\mathrm{C}) /$ from
  !> https://en.wikipedia.org/wiki/Elementary_charge
  real(rp),parameter :: q_e = 1.60218e-19_rp
  !> @}

  !> @defgroup Derived_constants Derived constants
  !> @{

  !> Bohr radius \f$ a_0 (\AA) \f$ from
  !> https://en.wikipedia.org/wiki/Bohr_radius
  real(rp),parameter :: a_0 = 0.5291772107_rp
  !real(rp),parameter :: a_0 = 4*pi*epsilon_0*hbar**2/(m_e*q_e**2)
  !> Hartree energy \f$ E_{\mathrm{Ha}} (\mathrm{eV}) \f$ from
  !> https://en.wikipedia.org/wiki/Hartree
  real(rp),parameter :: e_ha = 27.211386245988_rp
  !real(rp),parameter :: e_ha = m_e*(q_e**2/(4*pi*epsilon_0*hbar))**2
  !> Rydberg energy \f$ E_{\mathrm{Ry}} (\mathrm{eV}) \f$ from
  !> https://en.wikipedia.org/wiki/Rydberg_constant
  !real(rp),parameter :: e_ry = 13.60569301_rp
  real(rp),parameter :: e_ry = e_ha/2.0_rp
  !> Electron gyromagnetic ratio \f$ \gamma_{\mathrm{e}}
  !> (\mathrm{fs}^{-1}.\mathrm{T}^{-1}) \f$ from
  !> https://en.wikipedia.org/wiki/Gyromagnetic_ratio
  real(rp),parameter :: gamma_e = 1.760859644e-4_rp
  !real(rp),parameter :: gamma_e = g_e*q_e/(2*m_e)
  !> Bohr magneton \f$ \mu_{\mathrm{B}} (\mathrm{eV.T}^{-1}) \f$ from
  !> https://en.wikipedia.org/wiki/Bohr_magneton
  real(rp),parameter :: mu_b = 5.78838180E-5_rp
  !real(rp),parameter :: mu_b = q_e*hbar/(2*m_e)
  !> Hartree atomic unit of time \f$ t_{\mathrm{Ha}} (fs) \f$
  real(rp),parameter :: t_ha = hbar/e_ha
  !> Rydberg atomic unit of time \f$ t_{\mathrm{Ry}} (fs) \f$
  real(rp),parameter :: t_ry = hbar/e_ry
  !> @}
end module constant_mod
