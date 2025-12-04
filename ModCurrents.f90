!
!Cheng, 02142020
!To calculate ionospheric currents from j=sigma dot E
!
!1. Define all variables in this module
!  1.1 Change Makefile to add ModCurrents.o to MODULES
!  1.2 Add usage of this module and 'call init_mod_currents' in init_grid.f90
!  1.3 Add usage of this module in ModSphereInterface.f90
!  1.4 Add usage of this module in calc_efield.f90
!  1.5 Add usage of this module in output_common.f90 for output
!2. Calculate currents in calc_currents.f90
!  2.1 Change Makefile to add calc_currents.o to OBJECTS
!  2.2 'call calc_currents' in calc_sources.f90 right after 'call calc_ion_v'
!
module ModCurrents

  use ModSizeGitm

  implicit none

  !This is the same as EField but with values on ghost cells
  real, allocatable :: IonosphericEField(:,:,:,:,:)

  !Horizontal Currents, in A/m2
  !J3D is the total ionospheric currents in geographical 
  !  coordinates
  !J3DUcB is the component due to UxB
  !J3DEF is the component due to external Electric Field
  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) :: &
       J3D,J3DUcB,J3DEF
  !Height-integrated J3D to obtain J2D
  real, dimension(-1:nLons+2, -1:nLats+2, 3) :: &
       J2D,J2DUcB,J2DEF
  !Expanding a 2D array to a 3D array with the same values
  !  on different altitudinal layers
  !The reason we have these parameters is that 
  !  UAM_Gradient_GC or UAM_Gradient works with 3D arrays 
  !  (lon,lat,alt)
  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) :: &
       J2D3D,J2D3DUcB,J2D3DEF

  !Field-aligned Currents, in nuA/m2=1e-6 A/m2
  !Refer to ionos_ch5 of OuluSpaceWiki or Lu1995
  !Jparallel=Jp
  !Jp1 = -1.*/sin(I)*alt_integral(divergence of SigmaDotE)
  !Jp2 = -1.*/sin(I)*divergence of (alt_integral(Sigma) Dot E)
  !So Jp1 is calculated from 3D currents, do divergence at 
  !  different altitude first, then do the integral along altitude
  !Jp2 is from 2D currents, from height-integrated conductance and 
  !  assuming E does not change with altitude, then get 
  !  height-integrated currents (2D currents). Divergence of 
  !  2D currents gives the parallel currents as well
  !These two methods can be identical under certain conditions, 
  !  google 'Leibniz'rule of Integral'
  real, dimension(-1:nLons+2, -1:nLats+2) :: &
       Jp1,Jp2, &
       Jp1UcB,Jp2UcB, &
       Jp1EF,Jp2EF

  !20231108, to add in reading and outputing Jp from SWMF outputs
  real, dimension(-1:nLons+2, -1:nLats+2) :: &
       JrSWMF,JpSWMF  

  !Divergence: EG = Gradient of East Component, NG = North Component
  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) ::           &
       J3DEG,J3DNG, &
       J3DUcBEG,J3DUcBNG, &
       J3DEFEG,J3DEFNG
  !
  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) ::           &
       J2D3DEG,J2D3DNG, &
       J2D3DUcBEG,J2D3DUcBNG, &
       J2D3DEFEG,J2D3DEFNG

  !Cheng, 20211107, smooth EFlux, 5-point
  !3-point test in get_potential.f90
  !ElectronEnergyFlux is a 2-dimensional variable, I do not know
  !  how to sync it with MPI because ModSphereInterface.f90 does 
  !  not include an example. So I made EFLux3D to be similar as 
  !  Rho, and follow the example of Rho.
  real, allocatable :: EFlux3D(:,:,:,:)

contains
  !=========================================================================
  subroutine init_mod_currents

    allocate(IonosphericEField(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3, nBlocks))
    !20211107, Cheng
    allocate(EFlux3D(-1:nLons+2, -1:nLats+2, -1:nAlts+2, nBlocks))

  end subroutine init_mod_currents
  !=========================================================================
  subroutine clean_mod_currents

    deallocate(IonosphericEField)
    !20211107, Cheng
    deallocate(EFlux3D)

  end subroutine clean_mod_currents
  !=========================================================================
end module ModCurrents
