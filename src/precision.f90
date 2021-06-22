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
!  precision.f90
!  DyNaMol
module precision_mod
  implicit none
  !> Working precision of integer variables
  integer,parameter :: ip = 4
  !> Maximum integer length in decimal system
  integer,parameter :: il = 12 ! for ip=4
  !integer,parameter :: il = 21 ! for ip=8

  !> Working precision of real variables (4: single precision,
  ! 8: double precision)
  integer,parameter :: rp = 8
  !> Maximum real length in decimal system
  !integer,parameter :: rl = 17 ! for rp=4
  integer,parameter :: rl = 26 ! for rp=8

  ! Definition that guarantees to have at least 15 significant digits of
  ! precision and an exponent range of at least 307 (Lemmon and Schafer, 2005,
  ! p. 20–22)
  !integer, parameter :: wp = selected_real_kind(15, 307)
  ! In Fortran 95, the following constants are available:
  !integer, parameter :: sp = selected_real_kind(6, 37)
  !integer, parameter :: dp = selected_real_kind(15, 307)
  !integer, parameter :: qp = selected_real_kind(33, 4931)
end module precision_mod
