!
!Cheng, 02142020
!
subroutine calc_currents(iBlock)

  use ModGITM
  !use ModInputs
  !use ModConstants
  !drop Sigma_0 in my current calculations, 08122020
  use ModElectrodynamics, only : Sigma_Pedersen,Sigma_Hall
  use ModCurrents
  use ModConstants, only : pi

  implicit none

  integer, intent(in) :: iBlock
  integer :: iLon, iLat, iAlt
  integer :: iDir

  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3) ::           &
       Epara,Eperp,EperpCrossBunit,UCrossB,Eprime
  real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2   ) ::           &
       EDotB
  !Cheng, 20211107
  !real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2   ) ::           &
  !     EFlux3Dtmp
  !real, dimension(-1:nLons+2, -1:nLats+2, -1:nAlts+2, 3   ) ::           &
  !     EFluxG

  !---------------------------------------------------------------------------

  call report("Current Calculation",1)
  call start_timing("Calculating Currents")

  UCrossB(:,:,:,iEast_)  =   (Velocity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iNorth_,iBlock) *&
                              B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iUp_   ,iBlock) - &
                              Velocity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iUp_   ,iBlock) *&
                              B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iNorth_,iBlock))
  UCrossB(:,:,:,iNorth_) = - (Velocity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iEast_ ,iBlock) *&
                              B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iUp_   ,iBlock) - &
                              Velocity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iUp_   ,iBlock) *&
                              B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iEast_ ,iBlock))
  UCrossB(:,:,:,iUp_)    =   (Velocity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iEast_ ,iBlock) *&
                              B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iNorth_,iBlock) - &
                              Velocity(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iNorth_,iBlock) *&
                              B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iEast_ ,iBlock))

  !I decided to do the same things for three times, the 
  !  different things will be different electric field 
  !  fed to Eprime. Should have a better way???
  !For J3D, Eprime=E+UxB
  !For J3DUcB, Eprime=UxB
  !For J3DEF, Eprime=E
  !
  !
  !------------------------------Total-------------------------
  Eprime(:,:,:,:)=0.
  EDotB=0.
  Epara(:,:,:,:)=0.
  Eperp(:,:,:,:)=0.
  !
  Eprime(:,:,:,:)=IonosphericEField(-1:nLons+2,-1:nLats+2,-1:nAlts+2,:,iBlock)+&
                  UCrossB(:,:,:,:)
  !
  EDotB=sum(Eprime(-1:nLons+2,-1:nLats+2,-1:nAlts+2,:) *&
            B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:3,iBlock), dim=4)
  do iDir=1,3
     Epara(:,:,:,iDir)=EDotB*B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iDir,iBlock)/&
          B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iMag_,iBlock)**2
     Eperp(:,:,:,iDir)=Eprime(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iDir)-&
          Epara(:,:,:,iDir)
  enddo
  !
  EperpCrossBunit(:,:,:,iEast_)  =   (Eperp(:,:,:,iNorth_) * B0(:,:,:,iUp_   ,iBlock) - &
                                      Eperp(:,:,:,iUp_   ) * B0(:,:,:,iNorth_,iBlock))/&
                                      B0(:,:,:,iMag_,iBlock)
  EperpCrossBunit(:,:,:,iNorth_) = - (Eperp(:,:,:,iEast_ ) * B0(:,:,:,iUp_   ,iBlock) - &
                                      Eperp(:,:,:,iUp_   ) * B0(:,:,:,iEast_ ,iBlock))/&
                                      B0(:,:,:,iMag_,iBlock)
  EperpCrossBunit(:,:,:,iUp_)    =   (Eperp(:,:,:,iEast_ ) * B0(:,:,:,iNorth_,iBlock) - &
                                      Eperp(:,:,:,iNorth_) * B0(:,:,:,iEast_ ,iBlock))/&
                                      B0(:,:,:,iMag_,iBlock)
  !
  J3D(:,:,:,:)=0.0
  J2D(:,:,:)=0.0
  do iDir = 1, 3
     do iLat = -1, nLats+2
     do iLon = -1, nLons+2
     do iAlt = -1, nAlts+2
        J3D(iLon,iLat,iAlt,iDir) = &
             Sigma_Pedersen(iLon,iLat,iAlt)*Eperp(iLon,iLat,iAlt,iDir)-&
             Sigma_Hall(iLon,iLat,iAlt)*EperpCrossBunit(iLon,iLat,iAlt,iDir)!+&
             !Sigma_0(iLon,iLat,iAlt)*Epara(iLon,iLat,iAlt,iDir)
        !J2D
        if (iAlt >= 2 .and. iAlt <= nAlts-1) then 
           J2D(iLon,iLat,iDir)=J2D(iLon,iLat,iDir) + &
                J3D(iLon,iLat,iAlt,iDir)* &
                ( Altitude_GB(iLon,iLat,iAlt+1,iBlock) - &
                Altitude_GB(iLon,iLat,iAlt-1,iBlock) )/2.0
        endif
        !
     enddo
     enddo
     enddo
  enddo
  !
  J3DEG(:,:,:,:)=0.0
  J3DNG(:,:,:,:)=0.0
  !What is the difference between UAM_Gradient and UAM_Gradient_GC???
  call UAM_Gradient_GC(J3D(:,:,:,iEast_),  J3DEG, iBlock)
  call UAM_Gradient_GC(J3D(:,:,:,iNorth_), J3DNG, iBlock)
  !
  Jp1(:,:)=0.0
  do iLat = -1, nLats+2
  do iLon = -1, nLons+2
  do iAlt = 2, nAlts-1
     Jp1(iLon,iLat)=Jp1(iLon,iLat) + &
          (-1.0)/sin( DipAngle(iLon,iLat,iAlt,iBlock) )* &
          ( J3DEG(iLon,iLat,iAlt,iEast_ ) + &
            J3DNG(iLon,iLat,iAlt,iNorth_) )* &
          ( Altitude_GB(iLon,iLat,iAlt+1,iBlock) - &
            Altitude_GB(iLon,iLat,iAlt-1,iBlock) )/2.0* &
          1.e6 !so that we get nuA/m2
  enddo
  enddo
  enddo
  !
  J2D3D(:,:,:,:)=0.0
  do iAlt=-1,nAlts+2
     J2D3D(:,:,iAlt,:)=J2D(:,:,:)
  enddo
  J2D3DEG(:,:,:,:)=0.0
  J2D3DNG(:,:,:,:)=0.0
  call UAM_GRadient_GC(J2D3D(:,:,:,iEast_ ),  J2D3DEG, iBlock)
  call UAM_GRadient_GC(J2D3D(:,:,:,iNorth_),  J2D3DNG, iBlock)
  Jp2(:,:)=0.0
  Jp2(:,:)=(-1.0)/sin( DipAngle(:,:,nAlts,iBlock) )* & 
       ( J2D3DEG(:,:,nAlts,iEast_ )+ &
       J2D3DNG(:,:,nAlts,iNorth_) )* &
       1.e6 !so that we get nuA/m2
  !
  !
  !
  !
  !
  !Cheng, 20231108, to pass JrSWMF to JpSWMF
  JpSWMF(:,:)=0.0
  JpSWMF(:,:)=(1.0)/sin( DipAngle(:,:,nAlts,iBlock) )* &
       JrSWMF(:,:)
  !
  !------------------------------Total-------------------------
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !------------------------------UcrossB-------------------------
  Eprime(:,:,:,:)=0.
  EDotB=0.
  Epara(:,:,:,:)=0.
  Eperp(:,:,:,:)=0.
  !
  Eprime(:,:,:,:)=UCrossB(:,:,:,:)
  !
  EDotB=sum(Eprime(-1:nLons+2,-1:nLats+2,-1:nAlts+2,:) *&
            B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:3,iBlock), dim=4)
  do iDir=1,3
     Epara(:,:,:,iDir)=EDotB*B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iDir,iBlock)/&
          B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iMag_,iBlock)**2
     Eperp(:,:,:,iDir)=Eprime(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iDir)-&
          Epara(:,:,:,iDir)
  enddo
  !
  EperpCrossBunit(:,:,:,iEast_)  =   (Eperp(:,:,:,iNorth_) * B0(:,:,:,iUp_   ,iBlock) - &
                                      Eperp(:,:,:,iUp_   ) * B0(:,:,:,iNorth_,iBlock))/&
                                      B0(:,:,:,iMag_,iBlock)
  EperpCrossBunit(:,:,:,iNorth_) = - (Eperp(:,:,:,iEast_ ) * B0(:,:,:,iUp_   ,iBlock) - &
                                      Eperp(:,:,:,iUp_   ) * B0(:,:,:,iEast_ ,iBlock))/&
                                      B0(:,:,:,iMag_,iBlock)
  EperpCrossBunit(:,:,:,iUp_)    =   (Eperp(:,:,:,iEast_ ) * B0(:,:,:,iNorth_,iBlock) - &
                                      Eperp(:,:,:,iNorth_) * B0(:,:,:,iEast_ ,iBlock))/&
                                      B0(:,:,:,iMag_,iBlock)
  !
  J3DUcB(:,:,:,:)=0.0
  J2DUcB(:,:,:)=0.0
  do iDir = 1, 3
     do iLat = -1, nLats+2
     do iLon = -1, nLons+2
     do iAlt = -1, nAlts+2
        J3DUcB(iLon,iLat,iAlt,iDir) = &
             Sigma_Pedersen(iLon,iLat,iAlt)*Eperp(iLon,iLat,iAlt,iDir)-&
             Sigma_Hall(iLon,iLat,iAlt)*EperpCrossBunit(iLon,iLat,iAlt,iDir)!+&
             !Sigma_0(iLon,iLat,iAlt)*Epara(iLon,iLat,iAlt,iDir)
        !J2DUcB
        if (iAlt >= 2 .and. iAlt <= nAlts-1) then 
           J2DUcB(iLon,iLat,iDir)=J2DUcB(iLon,iLat,iDir) + &
                J3DUcB(iLon,iLat,iAlt,iDir)* &
                ( Altitude_GB(iLon,iLat,iAlt+1,iBlock) - &
                Altitude_GB(iLon,iLat,iAlt-1,iBlock) )/2.0
        endif
        !
     enddo
     enddo
     enddo
  enddo
  !
  J3DUcBEG(:,:,:,:)=0.0
  J3DUcBNG(:,:,:,:)=0.0
  !What is the difference between UAM_Gradient and UAM_Gradient_GC???
  call UAM_Gradient_GC(J3DUcB(:,:,:,iEast_),  J3DUcBEG, iBlock)
  call UAM_Gradient_GC(J3DUcB(:,:,:,iNorth_), J3DUcBNG, iBlock)
  !
  Jp1UcB(:,:)=0.0
  do iLat = -1, nLats+2
  do iLon = -1, nLons+2
  do iAlt = 2, nAlts-1
     Jp1UcB(iLon,iLat)=Jp1UcB(iLon,iLat) + &
          (-1.0)/sin( DipAngle(iLon,iLat,iAlt,iBlock) )* &
          ( J3DUcBEG(iLon,iLat,iAlt,iEast_ ) + &
            J3DUcBNG(iLon,iLat,iAlt,iNorth_) )* &
          ( Altitude_GB(iLon,iLat,iAlt+1,iBlock) - &
            Altitude_GB(iLon,iLat,iAlt-1,iBlock) )/2.0* &
          1.e6 !so that we get nuA/m2
  enddo
  enddo
  enddo
  !
  J2D3DUcB(:,:,:,:)=0.0
  do iAlt=-1,nAlts+2
     J2D3DUcB(:,:,iAlt,:)=J2DUcB(:,:,:)
  enddo
  J2D3DUcBEG(:,:,:,:)=0.0
  J2D3DUcBNG(:,:,:,:)=0.0
  call UAM_GRadient_GC(J2D3DUcB(:,:,:,iEast_ ),  J2D3DUcBEG, iBlock)
  call UAM_GRadient_GC(J2D3DUcB(:,:,:,iNorth_),  J2D3DUcBNG, iBlock)
  Jp2UcB(:,:)=0.0
  Jp2UcB(:,:)=(-1.0)/sin( DipAngle(:,:,nAlts,iBlock) )* & 
       ( J2D3DUcBEG(:,:,nAlts,iEast_ )+ &
       J2D3DUcBNG(:,:,nAlts,iNorth_) )* &
       1.e6 !so that we get nuA/m2
  !------------------------------UcrossB-------------------------
  !
  !
  !
  !
  !
  !
  !
  !
  !
  !------------------------------EField-------------------------
  Eprime(:,:,:,:)=0.
  EDotB=0.
  Epara(:,:,:,:)=0.
  Eperp(:,:,:,:)=0.
  !
  Eprime(:,:,:,:)=IonosphericEField(-1:nLons+2,-1:nLats+2,-1:nAlts+2,:,iBlock)
  !
  EDotB=sum(Eprime(-1:nLons+2,-1:nLats+2,-1:nAlts+2,:) *&
            B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,1:3,iBlock), dim=4)
  do iDir=1,3
     Epara(:,:,:,iDir)=EDotB*B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iDir,iBlock)/&
          B0(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iMag_,iBlock)**2
     Eperp(:,:,:,iDir)=Eprime(-1:nLons+2,-1:nLats+2,-1:nAlts+2,iDir)-&
          Epara(:,:,:,iDir)
  enddo
  !
  EperpCrossBunit(:,:,:,iEast_)  =   (Eperp(:,:,:,iNorth_) * B0(:,:,:,iUp_   ,iBlock) - &
                                      Eperp(:,:,:,iUp_   ) * B0(:,:,:,iNorth_,iBlock))/&
                                      B0(:,:,:,iMag_,iBlock)
  EperpCrossBunit(:,:,:,iNorth_) = - (Eperp(:,:,:,iEast_ ) * B0(:,:,:,iUp_   ,iBlock) - &
                                      Eperp(:,:,:,iUp_   ) * B0(:,:,:,iEast_ ,iBlock))/&
                                      B0(:,:,:,iMag_,iBlock)
  EperpCrossBunit(:,:,:,iUp_)    =   (Eperp(:,:,:,iEast_ ) * B0(:,:,:,iNorth_,iBlock) - &
                                      Eperp(:,:,:,iNorth_) * B0(:,:,:,iEast_ ,iBlock))/&
                                      B0(:,:,:,iMag_,iBlock)
  !
  J3DEF(:,:,:,:)=0.0
  J2DEF(:,:,:)=0.0
  do iDir = 1, 3
     do iLat = -1, nLats+2
     do iLon = -1, nLons+2
     do iAlt = -1, nAlts+2
        J3DEF(iLon,iLat,iAlt,iDir) = &
             Sigma_Pedersen(iLon,iLat,iAlt)*Eperp(iLon,iLat,iAlt,iDir)-&
             Sigma_Hall(iLon,iLat,iAlt)*EperpCrossBunit(iLon,iLat,iAlt,iDir)!+&
             !Sigma_0(iLon,iLat,iAlt)*Epara(iLon,iLat,iAlt,iDir)
        !J2DEF
        if (iAlt >= 2 .and. iAlt <= nAlts-1) then 
           J2DEF(iLon,iLat,iDir)=J2DEF(iLon,iLat,iDir) + &
                J3DEF(iLon,iLat,iAlt,iDir)* &
                ( Altitude_GB(iLon,iLat,iAlt+1,iBlock) - &
                Altitude_GB(iLon,iLat,iAlt-1,iBlock) )/2.0
        endif
        !
     enddo
     enddo
     enddo
  enddo
  !
  J3DEFEG(:,:,:,:)=0.0
  J3DEFNG(:,:,:,:)=0.0
  !What is the difference between UAM_Gradient and UAM_Gradient_GC???
  call UAM_Gradient_GC(J3DEF(:,:,:,iEast_),  J3DEFEG, iBlock)
  call UAM_Gradient_GC(J3DEF(:,:,:,iNorth_), J3DEFNG, iBlock)
  !
  Jp1EF(:,:)=0.0
  do iLat = -1, nLats+2
  do iLon = -1, nLons+2
  do iAlt = 2, nAlts-1
     Jp1EF(iLon,iLat)=Jp1EF(iLon,iLat) + &
          (-1.0)/sin( DipAngle(iLon,iLat,iAlt,iBlock) )* &
          ( J3DEFEG(iLon,iLat,iAlt,iEast_ ) + &
            J3DEFNG(iLon,iLat,iAlt,iNorth_) )* &
          ( Altitude_GB(iLon,iLat,iAlt+1,iBlock) - &
            Altitude_GB(iLon,iLat,iAlt-1,iBlock) )/2.0* &
          1.e6 !so that we get nuA/m2
  enddo
  enddo
  enddo
  !
  J2D3DEF(:,:,:,:)=0.0
  do iAlt=-1,nAlts+2
     J2D3DEF(:,:,iAlt,:)=J2DEF(:,:,:)
  enddo
  J2D3DEFEG(:,:,:,:)=0.0
  J2D3DEFNG(:,:,:,:)=0.0
  call UAM_GRadient_GC(J2D3DEF(:,:,:,iEast_ ),  J2D3DEFEG, iBlock)
  call UAM_GRadient_GC(J2D3DEF(:,:,:,iNorth_),  J2D3DEFNG, iBlock)
  Jp2EF(:,:)=0.0
  Jp2EF(:,:)=(-1.0)/sin( DipAngle(:,:,nAlts,iBlock) )* & 
       ( J2D3DEFEG(:,:,nAlts,iEast_ )+ &
       J2D3DEFNG(:,:,nAlts,iNorth_) )* &
       1.e6 !so that we get nuA/m2
  !------------------------------UcrossB-------------------------

  !----------------------------debugging-------------------------
  !if ( abs(Longitude(9,iBlock)*180./pi-22.*15.) <= 15. .AND. &
  !     abs(Latitude(9,iBlock)*180./pi-65.) <= 15. )then
  !   write(*,*) "Longitude(9,iBlock)*180./pi,,Latitude(9,iBlock)*180./pi,J2DEF(:,9,1),AAA,J2D3DEFEG(:,9,nAlts,1)",&
  !        Longitude(9,iBlock)*180./pi,Latitude(9,iBlock)*180./pi,&
  !        J2DEF(:,9,1),"AAA",J2D3DEFEG(:,9,nAlts,1)
  !endif
  !----------------------------debugging-------------------------

  !----------------------------debugging-------------------------
  !do iAlt=-1,nAlts+2
  !   EFlux3Dtmp(:,:,iAlt)=ElectronEnergyFlux(:,:)
  !enddo
  !call UAM_GRadient_GC(EFlux3Dtmp(:,:,:),  EFluxG, iBlock)
  !if ( abs(Longitude(9,iBlock)*180./pi-312.5) <= 5. .AND. &
  !     abs(Latitude(0,iBlock)*180./pi-53.) <= 2. )then
  !   !write(*,*) "Longitude(9,iBlock)*180./pi,,Latitude(9,iBlock)*180./pi,J2DEF(:,9,1),AAA,J2D3DEFEG(:,9,nAlts,1)",&
  !   !     Longitude(9,iBlock)*180./pi,Latitude(9,iBlock)*180./pi,&
  !   !     J2DEF(:,9,1),"AAA",J2D3DEFEG(:,9,nAlts,1)
  !   do iLon=-1,nLons+2 
  !      !iAlt=8 for 110km
  !      write(*,*) iLon,Longitude(iLon,iBlock)*180./pi,Latitude(5,iBlock)*180./pi,ElectronEnergyFlux(iLon,5),&
  !           EFluxG(iLon,5,8,iEast_)
  !   enddo
  !endif
  !----------------------------debugging-------------------------

  call end_timing("Calculating Currents")

end subroutine calc_currents
