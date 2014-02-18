module Calc_Combo

use Calc

implicit none

contains

!This should give identical results to the CalcF12Correction routine,
!However, this should be much faster, as all is done in one loop, and so less i/o
    subroutine CalcF12Correction_Combo(Energy_F12)
        implicit none
        real(dp), intent(inout) :: Energy_F12
        integer :: unit_raaaa,unit_rabab,unit_rabba,ierr
        integer :: v,w,vw_pair,unit_r_read_ccoo,loop,ispin,unit_g_read,i,unit_r_read_oocc
        integer :: unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba,unit_g_oocc
        integer :: vw_pair_complete,vw_pair_orb,unit_PairE_write
        real(dp) :: DDOT,Temp,EnergyTemp_FTF(2)
        real(dp) :: EnergyTemp_B(2,13),EnergyTemp_V(2,3),EnergyTemp_X(2,3),EnergyTemp_Comm(2)
        real(dp) , allocatable :: FockGam(:),GamFockGam(:)
        real(dp) , allocatable :: TotOrbMat(:),OrbOrbMat(:),OrbCabsMat(:),OrbOrbMat_2(:)
        real(dp) , allocatable :: OrbCabsMat_2(:),OrbCabsMat_3(:),OrbCabsMat_4(:),OrbCabsMat_5(:)
        real(dp) , allocatable :: OrbCabsMat_6(:),OrbCabsMat_7(:),TotTotMat(:),TotTotMat_2(:)
        real(dp) , allocatable :: PairEnergies(:),TotTotMat_3(:)
        real(dp) , allocatable :: R_vw_oo(:),R_vw_cc(:),G_vw_oo(:),R_vw_ooT(:),G_vw_cc(:)
        real(dp) , allocatable :: AR1_vw_oo(:),AR1_vw_ooT(:),AG0_vw_oo(:),AR1_vw_cc(:)
        real(dp) , allocatable :: AR1_vw_ao(:),AR1_vw_oa(:),AR1_vw_ca(:),AG0_vw_ao(:) 
        real(dp) , allocatable :: RTLD_vw_oo(:),RTLD_vw_ao(:),RTLD_vw_oa(:),RTLD_vw_cc(:) 
        real(dp) , allocatable :: PhiTLD_vw_oo(:),PhiTLD_vw_oa(:)
        logical :: tFrozenPair
        character(*), parameter :: t_r="CalcF12Correction_Combo"
                    
        EnergyTemp_B(:,:) = 0.0_dp    !Contains the individual spin components for each term
        EnergyTemp_V(:,:) = 0.0_dp    !Contains the individual spin components for each term
        EnergyTemp_X(:,:) = 0.0_dp    !Contains the individual spin components for each term
        EnergyTemp_Comm(:) = 0.0_dp
        EnergyTemp_FTF(:) = 0.0_dp

        unit_r_read_ccoo=get_free_unit()
        open(unit_r_read_ccoo,file='R_Dir_cc.oo',status='old',form='unformatted',access='direct',action='read',recl=reclen_nmo_sq)
        unit_r_read_oocc=get_free_unit()
        open(unit_r_read_oocc,file='R_Dir_oo.cc',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)
        unit_g_read=get_free_unit()
        open(unit_g_read,file='G_Dir_ccoo',status='old',form='unformatted',access='direct',action='read',recl=reclen_nmo_sq)
        unit_g_oocc=get_free_unit()
        open(unit_g_oocc,file='G_Dir_oocc',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)
        
        unit_raaaa=get_free_unit()
        open(unit_raaaa,file='RTld_aaaa',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)
        unit_rabab=get_free_unit()
        open(unit_rabab,file='RTld_abab',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)
        unit_rabba=get_free_unit()
        open(unit_rabba,file='RTld_abba',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)

        unit_rtld_Phi_aaaa=get_free_unit()
        open(unit_rtld_Phi_aaaa,file='RTld_Phi_aaaa',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)
        unit_rtld_Phi_abab=get_free_unit()
        open(unit_rtld_Phi_abab,file='RTld_Phi_abab',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)
        unit_rtld_Phi_abba=get_free_unit()
        open(unit_rtld_Phi_abba,file='RTld_Phi_abba',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)

        allocate(PairEnergies(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Allocation err")
        PairEnergies(:) = 0.0_dp
        
        !Allocate temp arrays for partially contracted quantities
        allocate(FockGam(tnmo_sq),stat=ierr)
        allocate(GamFockGam(tnmo_sq),stat=ierr)
        allocate(TotOrbMat(tntmo_nmo),stat=ierr)
        allocate(OrbOrbMat(tnmo_sq),stat=ierr)
        allocate(OrbOrbMat_2(tnmo_sq),stat=ierr)
        allocate(OrbCabsMat(taux_nmo),stat=ierr)
        allocate(OrbCabsMat_2(taux_nmo),stat=ierr)
        allocate(OrbCabsMat_3(taux_nmo),stat=ierr)
        allocate(OrbCabsMat_4(taux_nmo),stat=ierr)
        allocate(OrbCabsMat_5(taux_nmo),stat=ierr)
        allocate(OrbCabsMat_6(taux_nmo),stat=ierr)
        allocate(OrbCabsMat_7(taux_nmo),stat=ierr)
        allocate(TotTotMat(tntmo_sq),stat=ierr)
        if(tExplicitCommTerm) then
            allocate(TotTotMat_2(tntmo_sq),stat=ierr)
            allocate(TotTotMat_3(tntmo_sq),stat=ierr)
        endif
        if(ierr.ne.0) call stop_all(t_r,"Allocation err")

        !Allocate memory for integrals
        allocate(R_vw_oo(tnmo_sq),stat=ierr)
        allocate(R_vw_cc(tntmo_sq),stat=ierr)
        allocate(G_vw_oo(tnmo_sq),stat=ierr)
        allocate(G_vw_cc(tntmo_sq),stat=ierr)
        allocate(R_vw_ooT(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Allocation err")

        !Allocate memory for antisymmetrised quantities
        allocate(AR1_vw_oo(tnmo_sq),stat=ierr)
        allocate(AR1_vw_ooT(tnmo_sq),stat=ierr)
        allocate(AG0_vw_oo(tnmo_sq),stat=ierr)
        allocate(AG0_vw_ao(taux_nmo),stat=ierr)
        allocate(AR1_vw_cc(tntmo_sq),stat=ierr)
        allocate(AR1_vw_ao(taux_nmo),stat=ierr)
        allocate(AR1_vw_oa(taux_nmo),stat=ierr)
        allocate(AR1_vw_ca(tntmo_nxmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Allocation err")

        !Allocate memory for RTLD arrays
        allocate(RTLD_vw_oo(tnmo_sq),stat=ierr)
        allocate(RTLD_vw_ao(taux_nmo),stat=ierr)
        allocate(RTLD_vw_oa(taux_nmo),stat=ierr)
        allocate(RTLD_vw_cc(tntmo_sq),stat=ierr)
        allocate(PhiTLD_vw_oa(taux_nmo),stat=ierr)
        allocate(PhiTLD_vw_oo(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Allocation err")

        !Construct quantities which can be precomputed...
        !Precomute f_qr Gam_rs
        call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,GenFockOrb(:,:),tnmo,OneRDM(:,:),tnmo,0.0_dp,FockGam,tnmo)
        !Now precomute contraction with Gam_pq - this will be used for B4
        call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,OneRDM(:,:),tnmo,FockGam,tnmo,0.0_dp,GamFockGam,tnmo)
        !This will be used for B6 - f_p'p Gam_pq
        call DGEMM('N','N',tntmo,tnmo,tnmo,1.0_dp,GenFockCompOrb(:,:),tntmo,OneRDM(:,:),tnmo,0.0_dp,TotOrbMat,tntmo)

        write(6,'(A)') "Entering main loop..."

        do v=1,tnmo

            !If doing MP2F12, only need to loop over occupied orbital pairs here...
            if(tHFRDMs.and.(.not.IsOrbOcc(v))) cycle

            if(int(real((v-1)*100,dp)/real(tnmo,dp)).gt.int(real((v-2)*100,dp)/real(tnmo,dp))) then
                write(6,"(I4,A)") int(real((v-1)*100,dp)/real(tnmo,dp)) , " % complete"
                call flush(6)
            endif
            do w=1,tnmo

                !If doing MP2F12, only need to loop over occupied orbital pairs here...
                if(tHFRDMs.and.(.not.IsOrbOcc(w))) cycle

                tFrozenPair = (IsOrbFrz(v).or.IsOrbFrz(w))

                vw_pair_complete = FindSqPairInd(v,w,tntmo)  
                vw_pair_orb = FindSqPairInd(v,w,tnmo)

                !Load into memory all rq
                if(.not.tFrozenPair) then
                    read(unit_r_read_ccoo,rec=vw_pair_complete) (R_vw_oo(loop),loop=1,tnmo_sq)   !q is fast index (opposite to matrix element)
                    read(unit_r_read_oocc,rec=vw_pair_orb) (R_vw_cc(loop),loop=1,tntmo_sq)       !q fast
                    call TransposeIntegralArr(R_vw_oo,R_vw_ooT,tnmo)    !Ensure t fast
                endif
                read(unit_g_oocc,rec=vw_pair_orb) (G_vw_cc(loop),loop=1,tntmo_sq)  !Read in (g_rs)^ta'
                read(unit_g_read,rec=vw_pair_complete) (G_vw_oo(loop),loop=1,tnmo_sq)  !Read in (g_rs)^tu with u fast.
                call TransposeIntegralArr_inplace(G_vw_oo,tnmo)    !Ensure t fast

                do ispin=1,2    !Loop over same spin and opposite spins, assuming abab = abba

                    !Antisymmetrise integrals as appropriate
                    !Combine these at the end
                    if(.not.tFrozenPair) then
                        call Antisym1Cusp(R_vw_oo,AR1_vw_oo,tnmo,ispin)
                        call Antisym1Cusp(R_vw_ooT,AR1_vw_ooT,tnmo,ispin)
                        call Antisym1Cusp(R_vw_cc,AR1_vw_cc,tntmo,ispin)
                        call Antisym1Cusp_ooao_from_oocc(R_vw_cc,AR1_vw_ao,ispin,.false.)
                        call Antisym1Cusp_ooao_from_oocc(R_vw_cc,AR1_vw_oa,ispin,.true.)    !This should be able to be removed by using 'T'
                        call Antisym1Cusp_ooca_from_oocc(R_vw_cc,AR1_vw_ca,ispin)
                    endif
                    call Antisym0Cusp(G_vw_oo,AG0_vw_oo,tnmo,ispin)
                    call Antisym0Cusp_ooao_from_oocc(G_vw_cc,AG0_vw_ao,ispin)

                    if(.not.tFrozenPair) then
                        call DGEMM('N','T',tnmo,tnmo,tnmo,1.0_dp,Genfockorb(1,1),tnmo,AR1_vw_oo,tnmo,0.0_dp,OrbOrbMat,tnmo) !For B1
                        call DGEMM('N','N',tnmo,tnmo,tnxmo,1.0_dp,GenFockOrbCABS(1,tnmo+1),tnmo,AR1_vw_ao,tnxmo,0.0_dp,OrbOrbMat_2,tnmo)  !For B2 
                        call DGEMM('N','N',tnxmo,tnmo,tnmo,1.0_dp,AR1_vw_ao,tnxmo,GamFockGam,tnmo,0.0_dp,OrbCabsMat_2,tnxmo)    !For B4
                        call DGEMM('N','N',tnxmo,tnmo,tnmo,1.0_dp,AR1_vw_oa,tnxmo,OneRDM(:,:),tnmo,0.0_dp,OrbCabsMat_3,tnxmo)  !For B5
                        call DGEMM('N','N',tnxmo,tnmo,tnxmo,1.0_dp,GenFockCABS(:,:),tnxmo,OrbCabsMat_3,tnxmo,0.0_dp,OrbCabsMat_4,tnxmo) !For B5
                        call DGEMM('T','N',tnxmo,tnmo,tntmo,1.0_dp,AR1_vw_ca,tntmo,TotOrbMat,tntmo,0.0_dp,OrbCabsMat_5,tnxmo)  !For B6
                        call DGEMM('T','T',tntmo,tntmo,tntmo,1.0_dp,GenExch,tntmo,AR1_vw_cc,tntmo,0.0_dp,TotTotMat,tntmo)  !For B13
                        call DGEMM('N','N',tnxmo,tnmo,tnmo,1.0_dp,AR1_vw_oa,tnxmo,OneRDM,tnmo,0.0_dp,OrbCabsMat_6,tnxmo)  !For X2
                        if(tExplicitCommTerm) then
                            !For explicit commutator terms (should not be added)
                            call DGEMM('T','T',tntmo,tntmo,tntmo,1.0_dp,GenFockComp(:,:),tntmo,AR1_vw_cc,tntmo,0.0_dp,TotTotMat_2,tntmo)    
                            !To calculate the FTF integrals by double RI...
                            call DGEMM('T','T',tntmo,tntmo,tntmo,1.0_dp,KinMat(:,:),tntmo,AR1_vw_cc,tntmo,0.0_dp,TotTotMat_3,tntmo)    
                        endif
                    endif
                    call DGEMM('N','N',tnxmo,tnmo,tnmo,1.0_dp,AG0_vw_ao,tnxmo,OneRDM,tnmo,0.0_dp,OrbCabsMat_7,tnxmo) !For V2

                    !Now load in R_TLD for the vw pair...
                    !Combine these at the end
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,RTLD_vw_oo,ispin,1,tnmo,1,tnmo) 
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,RTLD_vw_oa,ispin,1,tnmo,tnmo+1,tntmo) 
                    if(.not.tFrozenPair) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,RTLD_vw_ao,ispin,tnmo+1,tntmo,1,tnmo)
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,RTLD_vw_cc,ispin,1,tntmo,1,tntmo)

                        call GetAntisymR_TLD(unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba,v,w,PhiTLD_vw_oa,ispin,1,tnmo,tnmo+1,tntmo)
                        call GetAntisymR_TLD(unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba,v,w,PhiTLD_vw_oo,ispin,1,tnmo,1,tnmo)
                    endif
                    
!All contributions need to be x2, since we could flip all spins.
!This would need to be changed for open shell systems.
                    !V2/3 is needed even over frozen pairs.
                    Temp = DDOT(tnmo_sq,RTLD_vw_oo,1,AG0_vw_oo,1)
                    EnergyTemp_V(ispin,3) = EnergyTemp_V(ispin,3) - 0.5_dp * Temp
                    PairEnergies(vw_pair_orb) = PairEnergies(vw_pair_orb) - 0.5_dp * Temp

                    Temp = DDOT(taux_nmo,RTld_vw_oa,1,OrbCabsMat_7,1)
                    EnergyTemp_V(ispin,2) = EnergyTemp_V(ispin,2) - Temp
                    PairEnergies(vw_pair_orb) = PairEnergies(vw_pair_orb) - Temp
                    
                    if(.not.tFrozenPair) then
                        Temp = DDOT(tnmo_sq,RTLD_vw_oo,1,OrbOrbMat_2,1)
                        EnergyTemp_B(ispin,2) = EnergyTemp_B(ispin,2) - 0.5_dp * Temp 
                        PairEnergies(vw_pair_orb) = PairEnergies(vw_pair_orb) - 0.5_dp * Temp

                        Temp = DDOT(tnmo_sq,RTLD_vw_oo,1,OrbOrbMat,1)
                        EnergyTemp_B(ispin,1) = EnergyTemp_B(ispin,1) - 0.5_dp * Temp
                        PairEnergies(vw_pair_orb) = PairEnergies(vw_pair_orb) - 0.5_dp * Temp
                        
                        Temp = DDOT(taux_nmo,RTLD_vw_ao,1,OrbCabsMat,1)
                        EnergyTemp_B(ispin,3) = EnergyTemp_B(ispin,3) - 0.5_dp * Temp
                        PairEnergies(vw_pair_orb) = PairEnergies(vw_pair_orb) - 0.5_dp * Temp
                        
                        Temp = DDOT(taux_nmo,RTLD_vw_ao,1,OrbCabsMat_2,1)
                        EnergyTemp_B(ispin,4) = EnergyTemp_B(ispin,4) + 0.5_dp * Temp
                        PairEnergies(vw_pair_orb) = PairEnergies(vw_pair_orb) + 0.5_dp * Temp
                        
                        Temp = DDOT(taux_nmo,RTLD_vw_oa,1,OrbCabsMat_4,1)
                        EnergyTemp_B(ispin,5) = EnergyTemp_B(ispin,5) - 0.5_dp * Temp
                        PairEnergies(vw_pair_orb) = PairEnergies(vw_pair_orb) - 0.5_dp * Temp
                        
                        Temp = DDOT(taux_nmo,RTLD_vw_oa,1,OrbCabsMat_5,1)
                        EnergyTemp_B(ispin,6) = EnergyTemp_B(ispin,6) - 0.5_dp * Temp
                        PairEnergies(vw_pair_orb) = PairEnergies(vw_pair_orb) - 0.5_dp * Temp
                        
                        Temp = DDOT(tntmo_sq,RTLD_vw_cc,1,TotTotMat,1)
                        EnergyTemp_B(ispin,13) = EnergyTemp_B(ispin,13) - 0.5_dp * Temp
                        PairEnergies(vw_pair_orb) = PairEnergies(vw_pair_orb) - 0.5_dp * Temp
                        
                        Temp = DDOT(taux_nmo,PhiTLD_vw_oa,1,OrbCabsMat_6,1)
                        EnergyTemp_X(ispin,2) = EnergyTemp_X(ispin,2) + 0.5_dp * Temp
                        PairEnergies(vw_pair_orb) = PairEnergies(vw_pair_orb) + 0.5_dp * Temp
                        
                        Temp = DDOT(tnmo_sq,PhiTLD_vw_oo,1,AR1_vw_ooT,1)
                        EnergyTemp_X(ispin,3) = EnergyTemp_X(ispin,3) + 0.25_dp * Temp
                        PairEnergies(vw_pair_orb) = PairEnergies(vw_pair_orb) + 0.25_dp * Temp
                        
                        if(tExplicitCommTerm) then
                            Temp = DDOT(tntmo_sq,RTLD_vw_cc,1,TotTotMat_2,1)
                            EnergyTemp_Comm(ispin) = EnergyTemp_Comm(ispin) + 0.5_dp * Temp
                            
                            Temp = DDOT(tntmo_sq,RTLD_vw_cc,1,TotTotMat_3,1)
                            EnergyTemp_FTF(ispin) = EnergyTemp_FTF(ispin) + Temp
                        endif
                    endif

                enddo
            enddo
        enddo

        write(6,'(A)') "Main loop complete."
        write(6,*) ""
        call flush(6)

!        do i=1,tntmo
!            write(6,*) "orb, KE: ",i,tmat(i,i),Kinmat(i,i)
!        enddo

        EnergyTemp_B(:,7) = EnergyTemp_B(:,6)    !B7 is equivalent to B6
        EnergyTemp_B(:,3) = EnergyTemp_B(:,2)    !B3 is equivalent to B2
        
        write(6,"(A)")        "Energy breakdown from V term 2: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_V(1,2)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_V(2,2)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from V term 2:          ",EnergyTemp_V(1,2) + EnergyTemp_V(2,2)*2.0_dp
        Energy_F12 = Energy_F12 + EnergyTemp_V(1,2) + EnergyTemp_V(2,2)*2.0_dp
        
        write(6,"(A)")        "Energy breakdown from V term 3: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_V(1,3)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_V(2,3)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from V term 3:          ",EnergyTemp_V(1,3) + EnergyTemp_V(2,3)*2.0_dp
        Energy_F12 = Energy_F12 + EnergyTemp_V(1,3) + EnergyTemp_V(2,3)*2.0_dp

        write(6,"(A)")        "Energy breakdown from B term 1: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_B(1,1)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_B(2,1)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from B term 1:          ",EnergyTemp_B(1,1) + EnergyTemp_B(2,1)*2.0_dp
        Energy_F12 = Energy_F12 + EnergyTemp_B(1,1) + EnergyTemp_B(2,1)*2.0_dp
        
        write(6,"(A)")        "Energy breakdown from B term 2: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_B(1,2)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_B(2,2)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from B term 2:          ",EnergyTemp_B(1,2) + EnergyTemp_B(2,2)*2.0_dp
        Energy_F12 = Energy_F12 + EnergyTemp_B(1,2) + EnergyTemp_B(2,2)*2.0_dp
        
        write(6,"(A)")        "Energy breakdown from B term 3: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_B(1,3)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_B(2,3)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from B term 3:          ",EnergyTemp_B(1,3) + EnergyTemp_B(2,3)*2.0_dp
        Energy_F12 = Energy_F12 + EnergyTemp_B(1,3) + EnergyTemp_B(2,3)*2.0_dp
        
        write(6,"(A)")        "Energy breakdown from B term 4: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_B(1,4)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_B(2,4)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from B term 4:          ",EnergyTemp_B(1,4) + EnergyTemp_B(2,4)*2.0_dp
        Energy_F12 = Energy_F12 + EnergyTemp_B(1,4) + EnergyTemp_B(2,4)*2.0_dp
        
        write(6,"(A)")        "Energy breakdown from B term 5: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_B(1,5)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_B(2,5)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from B term 5:          ",EnergyTemp_B(1,5) + EnergyTemp_B(2,5)*2.0_dp
        Energy_F12 = Energy_F12 + EnergyTemp_B(1,5) + EnergyTemp_B(2,5)*2.0_dp
        
        write(6,"(A)")        "Energy breakdown from B term 6: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_B(1,6)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_B(2,6)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from B term 6:          ",EnergyTemp_B(1,6) + EnergyTemp_B(2,6)*2.0_dp
        Energy_F12 = Energy_F12 + EnergyTemp_B(1,6) + EnergyTemp_B(2,6)*2.0_dp
        
        write(6,"(A)")        "Energy breakdown from B term 7: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_B(1,7)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_B(2,7)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from B term 7:          ",EnergyTemp_B(1,7) + EnergyTemp_B(2,7)*2.0_dp
        Energy_F12 = Energy_F12 + EnergyTemp_B(1,7) + EnergyTemp_B(2,7)*2.0_dp
        
        write(6,"(A)")        "Energy breakdown from B term 13: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_B(1,13)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_B(2,13)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from B term 13:          ",EnergyTemp_B(1,13) + EnergyTemp_B(2,13)*2.0_dp
!        write(6,*) "EnergyTemp_13: ",EnergyTemp_B(:,13)
!        write(6,*) "EnergyTemp_FTF: ",EnergyTemp_FTF(:)
        Energy_F12 = Energy_F12 + EnergyTemp_B(1,13) + EnergyTemp_B(2,13)*2.0_dp
        
        write(6,"(A)")        "Energy breakdown from X term 2: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_X(1,2)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_X(2,2)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from X term 2:          ",EnergyTemp_X(1,2) + EnergyTemp_X(2,2)*2.0_dp
        Energy_F12 = Energy_F12 + EnergyTemp_X(1,2) + EnergyTemp_X(2,2)*2.0_dp
        
        write(6,"(A)")        "Energy breakdown from X term 3: "
        write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_X(1,3)
        write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_X(2,3)*2.0_dp  !*2 for abba term
        write(6,"(A,F25.12)") "Total energy from X term 3:          ",EnergyTemp_X(1,3) + EnergyTemp_X(2,3)*2.0_dp
        Energy_F12 = Energy_F12 + EnergyTemp_X(1,3) + EnergyTemp_X(2,3)*2.0_dp

        if(tExplicitCommTerm) then
            write(6,"(A)")        "Expected approximate energy breakdown from sum of commutator terms 8 -> 13: "
            write(6,"(A)")        "Warning: V slow CABS convergence, and energies not included in final result."
            write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_Comm(1)
            write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_Comm(2)*2.0_dp  !*2 for abba term
            write(6,"(A,F25.12)") "Total energy from Commutator terms:  ",EnergyTemp_Comm(1) + EnergyTemp_Comm(2)*2.0_dp
            deallocate(TotTotMat_2)
            deallocate(TotTotMat_3)
        endif
        
        !Deallocate memory
        deallocate(FockGam,GamFockGam,TotOrbMat,OrbOrbMat,OrbOrbMat_2,OrbCabsMat,OrbCabsMat_2,OrbCabsMat_3,OrbCabsMat_4)
        deallocate(OrbCabsMat_5,OrbCabsMat_6,OrbCabsMat_7,TotTotMat)
        deallocate(R_vw_oo,R_vw_cc,G_vw_oo,R_vw_ooT,G_vw_cc)
        deallocate(AR1_vw_oo,AR1_vw_ooT,AG0_vw_oo,AR1_vw_cc,AR1_vw_ao,AR1_vw_oa,AR1_vw_ca,AG0_vw_ao)
        deallocate(RTLD_vw_oo,RTLD_vw_ao,RTLD_vw_oa,RTLD_vw_cc,PhiTLD_vw_oa,PhiTLD_vw_oo)

        call CalcFGTerm(Energy_F12,PairEnergies)
        call CalcTauTerm(Energy_F12,PairEnergies,.false.)
        call CalcTauTerm(Energy_F12,PairEnergies,.true.)
        call CalcRRTerms(Energy_F12,PairEnergies,EnergyTemp_FTF,.false.)
        call CalcRRTerms(Energy_F12,PairEnergies,EnergyTemp_FTF,.true.)
        
        if(tWritePairEnergies) then
            unit_PairE_write=get_free_unit()
            open(unit_PairE_write,file='GemPairEnergy',status='unknown',form='formatted',action='write')

            do v=1,tnmo
                do w=1,tnmo
                    write(unit_PairE_write,"(2I6,G25.10)") v,w,PairEnergies(FindSqPairInd(v,w,tnmo))
                enddo
            enddo

            close(unit_PairE_write)
        endif
        
        deallocate(PairEnergies)
        
!        close(unit_g_read,status='delete')
        close(unit_g_read)
        close(unit_raaaa,status='delete')
        close(unit_rabab,status='delete')
        close(unit_rabba,status='delete')
        close(unit_rtld_Phi_aaaa,status='delete')
        close(unit_rtld_Phi_abab,status='delete')
        close(unit_rtld_Phi_abba,status='delete')
        close(unit_r_read_ccoo)
        close(unit_r_read_oocc)
        close(unit_g_oocc)
!        close(unit_r_read_ccoo,status='delete')
!        close(unit_r_read_oocc,status='delete')
!        close(unit_g_oocc,status='delete')

    end subroutine CalcF12Correction_Combo


end module Calc_Combo
