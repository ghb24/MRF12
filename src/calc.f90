module Calc

use utils, only: get_free_unit
use basis_data
use basis, only: FindSqPairInd,FindCabsOrbPairInd
use errors, only: stop_all,warning
use int_management, only: loadFockOrb,TransposeIntegralArr_inplace,TransposeIntegralArr 
use input_data, only: tHFRDMs,tWritePairEnergies,tCalcCombo,tDaltonFormat,tComplexOrbs,tReadBin,tExplicitCommTerm

implicit none

contains

!Each energy contribution from B_vw^xy goes as: 1/8 ( B Gamma_xy^vw - B Gamma_xy^wv )
!We use rplus and rminus, which are already contracted with r integrals.
    subroutine CalcF12Correction(Energy_F12)
        implicit none
        integer :: unit_raaaa,unit_rabab,unit_rabba,ierr,OrbPair,unit_PairE_write
        real(dp) , allocatable :: r_qr(:),Antisym_qr(:),S_pq(:)
        real(dp) , allocatable :: g_qr(:),Antisym_g_qr(:),OrbCabsMat(:),OrbCabsMat_2(:)
        real(dp) , allocatable :: NmoNmoMat(:),NmoNmoMat_2(:),TotOrbMat(:),TotCabsMat(:),TotCabsMat_2(:)
        real(dp) , allocatable :: TotTotMat(:),TotTotMat_2(:),NmoNmoMat_3(:),TotTotMat_3(:)
        real(dp) , allocatable :: PairEnergies(:)
        character(*), parameter :: t_r="CalcBCorrection"
        integer :: v,w,vw_pair,unit_r_read_ccoo,loop,ispin,unit_g_read,i,unit_r_read_oocc
        integer :: unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba
        real(dp), intent(inout) :: Energy_F12
        real(dp) :: DDOT,TermE,Temp 
        real(dp) :: EnergyTemp(3),EnergyTemp_FTF(2)
                    
        EnergyTemp(:)=0.0_dp    !Contains the individual spin components for each term
!        write(6,*) "Fock bounds: ",lbound(fockorb,1),ubound(fockorb,1),lbound(fockorb,2),ubound(fockorb,2)

        unit_r_read_ccoo=get_free_unit()
        open(unit_r_read_ccoo,file='R_Dir_cc.oo',status='old',form='unformatted',access='direct',action='read',recl=reclen_nmo_sq)
        unit_r_read_oocc=get_free_unit()
        open(unit_r_read_oocc,file='R_Dir_oo.cc',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)
        unit_g_read=get_free_unit()
        open(unit_g_read,file='G_Dir_ccoo',status='old',form='unformatted',access='direct',action='read',recl=reclen_nmo_sq)
        
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

        !Load (fpk) into memory!
        !Use fock matrix temporarily for testing purposes.
        allocate(PairEnergies(tnmo_sq),stat=ierr)
        allocate(r_qr(tnmo_sq),stat=ierr)
        allocate(g_qr(tnmo_sq),stat=ierr)
        allocate(OrbCabsMat(taux_nmo),stat=ierr)
        allocate(OrbCabsMat_2(taux_nmo),stat=ierr)
        allocate(TotOrbMat(tntmo_nmo),stat=ierr)
        allocate(TotCabsMat(tntmo_nxmo),stat=ierr)
        allocate(TotCabsMat_2(tntmo_nxmo),stat=ierr)
        allocate(NmoNmoMat(tnmo_sq),stat=ierr)
        allocate(NmoNmoMat_2(tnmo_sq),stat=ierr)
        allocate(NmoNmoMat_3(tnmo_sq),stat=ierr)
        allocate(TotTotMat(tntmo_sq),stat=ierr)
        allocate(TotTotMat_2(tntmo_sq),stat=ierr)
        allocate(TotTotMat_3(tntmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Allocation err 2")
        allocate(Antisym_g_qr(tnmo_sq),stat=ierr)
        allocate(Antisym_qr(tnmo_sq),stat=ierr)
        allocate(S_pq(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Allocation err 3")

        PairEnergies(:) = 0.0_dp

        !TERM 1: 
        !-f_{vw}^{rq} (fpk)_r^p r_{pq}^{xy}
!        open(24,file="FR-Full",status='unknown')
        !Term 1: Sum in -1/4 (S_vw)^p*q (r_tld^vw)_p*q
        do v=1,tnmo
            if(IsOrbFrz(v)) cycle
            do w=1,tnmo
                if(IsOrbFrz(w)) cycle

!Since these are stored for r integrals over the complete space (but not used here, we need to index over the complete space)
                vw_pair = FindSqPairInd(v,w,tntmo)  
                OrbPair = FindSqPairInd(v,w,tnmo)
                !Load into memory all rq
                read(unit_r_read_ccoo,rec=vw_pair) (r_qr(loop),loop=1,tnmo_sq)   !q is fast index (opposite to matrix element)

                do ispin=1,3    !Loop over all spins

                    !Antisymmetrise integrals as appropriate
                    call Antisym1Cusp(r_qr,Antisym_qr,tnmo,ispin)
!                    call write_matrix('Antir',Antisym_qr,tnmo)
!                    call write_matrix('FockMat',fockorb,tnmo)

!Don't transpose f (since symmetric), transpose r_qr, so that r is the 'fast' index.
!Next two arguments - resulting dimension of the array
!Then - dimension of r that we are contracting over
!Fast index in S will be second index in f (p)
!*** Fast index in S_pq is p ***
                    !If ispin=1, we have (S_va,wa)^pa,qa = (S_vb,wb)^pb,qb
                    !If ispin=2, we have (S_va,wb)^pa,qb = (S_vb,wa)^pb,qa
                    !If ispin=3, we have (S_va,wb)^pb,qa = (S_vb,wa)^pa,qb
                    call DGEMM('N','T',tnmo,tnmo,tnmo,1.0_dp,Genfockorb(1,1),tnmo,Antisym_qr,tnmo,0.0_dp,S_pq,tnmo)

                    !Now load in R_TLD for the vw pair...
!Antisym_qr is now antisym for R (Gamma^{vw})_pq (p fast). Always want the same spin as for S
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,Antisym_qr,ispin,1,tnmo,1,tnmo) 
                    
!All contributions need to be x2, since we could flip all spins.
!This would need to be changed for open shell systems.
!First, contract the same spin types of S and R_Tld
                    Temp = -0.25_dp * DDOT(tnmo_sq,Antisym_qr,1,S_pq,1)
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp
                    PairEnergies(OrbPair) = PairEnergies(OrbPair) + Temp
                    
!                    BContrib = -2.0_dp * DDOT(tnmo_sq,Antisym_qr,1,S_pq,1)

                    if(ispin.eq.3) then
!Antisym_qr is (r_tld)^wv_pq for spin = 2, and _3 for spin = 3 (i.e. abba) 
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,Antisym_qr,2,1,tnmo,1,tnmo) 
                        Temp = 0.25_dp * DDOT(tnmo_sq,Antisym_qr,1,S_pq,1)
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,Antisym_qr,3,1,tnmo,1,tnmo) 
                        Temp = 0.25_dp * DDOT(tnmo_sq,Antisym_qr,1,S_pq,1)
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,Antisym_qr,1,1,tnmo,1,tnmo)
                        Temp = 0.25_dp * DDOT(tnmo_sq,Antisym_qr,1,S_pq,1)
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp
                    PairEnergies(OrbPair) = PairEnergies(OrbPair) + Temp
                        
                enddo
            enddo
        enddo

        write(6,"(A)") "Energy breakdown from B term 1: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from term 1: ",TermE
        Energy_F12 = Energy_F12 + TermE
        EnergyTemp(:) = 0.0_dp

!Third V term - loop over all orbitals (inc. frozen)
        do v=1,tnmo
            do w=1,tnmo

!Since these are stored for r integrals over the complete space (but not used here, we need to index over the complete space)
                vw_pair = FindSqPairInd(v,w,tntmo)  
                OrbPair = FindSqPairInd(v,w,tnmo)  
                !Load into memory all rq
                read(unit_g_read,rec=vw_pair) (g_qr(loop),loop=1,tnmo_sq)  !Read in (g_rs)^tu with u fast.
                call TransposeIntegralArr_inplace(g_qr,tnmo)    !Ensure t fast

                do ispin=1,3    !Loop over all spins

                    !Antisymmetrise integrals as appropriate
                    call Antisym0Cusp(g_qr,Antisym_g_qr,tnmo,ispin)
!                    call write_matrix('Antir',Antisym_qr,tnmo)
!                    call write_matrix('FockMat',fockorb,tnmo)

                    !Now load in R_TLD for the vw pair...
!Antisym_qr is now antisym for R (Gamma^{vw})_pq (p fast). 
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,Antisym_qr,ispin,1,tnmo,1,tnmo) 
                    
!All contributions need to be x2, since we could flip all spins.
!This would need to be changed for open shell systems.
                    !Contrib from term 3 of V:
                    Temp = -0.25_dp * DDOT(tnmo_sq,Antisym_g_qr,1,Antisym_qr,1)
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp
                    PairEnergies(OrbPair) = PairEnergies(OrbPair) + Temp

                    if(ispin.eq.3) then
!Antisym_qr is (r_tld)^wv_pq for spin = 2, and _3 for spin = 3 (i.e. abba) 
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,Antisym_qr,2,1,tnmo,1,tnmo) 
                        Temp = 0.25_dp * DDOT(tnmo_sq,Antisym_qr,1,Antisym_g_qr,1)
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,Antisym_qr,3,1,tnmo,1,tnmo) 
                        Temp = 0.25_dp * DDOT(tnmo_sq,Antisym_qr,1,Antisym_g_qr,1)
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,Antisym_qr,1,1,tnmo,1,tnmo)
                        Temp = 0.25_dp * DDOT(tnmo_sq,Antisym_qr,1,Antisym_g_qr,1)
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(OrbPair) = PairEnergies(OrbPair) + Temp
                        
                enddo
            enddo
        enddo
        
        write(6,"(A)") "Energy breakdown from V term 3: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from term 3: ",TermE
        Energy_F12 = Energy_F12 + TermE
        EnergyTemp(:) = 0.0_dp
        
        !B Term 3
!       -1/4 ((r_vw)^pq f_p^a' (rtld^vw)_a'q 
        do v=1,tnmo
            if(IsOrbFrz(v)) cycle
            do w=1,tnmo
                if(IsOrbFrz(w)) cycle

                vw_Pair = FindSqPairInd(v,w,tntmo)
                OrbPair = FindSqPairInd(v,w,tnmo)  
                !Get (g_rs)^ta'
                read(unit_r_read_ccoo,rec=vw_pair) (r_qr(loop),loop=1,tnmo_sq)

                do ispin=1,3

                    !Antisymmetrise over p q (one geminal indices)
                    call Antisym1Cusp(r_qr,Antisym_qr,tnmo,ispin)

                    !Contract with f_p^a'
                    call DGEMM('T','T',tnxmo,tnmo,tnmo,1.0_dp,GenFockOrbCABS(1,tnmo+1),tnmo,Antisym_qr,tnmo,0.0_dp,OrbCabsMat,tnxmo)
                    !In OrbCabsMat, a' is fast

                    !Get (Rtld^vw)_a'q and (Rtld^wv)_a'q
                    !a' is returned as fast here...
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,OrbCabsMat_2,ispin,tnmo+1,tntmo,1,tnmo)

                    !Contract over a'q and sum in to energy
                    Temp = -0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp
                    PairEnergies(OrbPair) = PairEnergies(OrbPair) + Temp

                    if(ispin.eq.3) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat_2,2,tnmo+1,tntmo,1,tnmo)
                        Temp = 0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat_2,3,tnmo+1,tntmo,1,tnmo)
                        Temp = 0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat_2,1,tnmo+1,tntmo,1,tnmo)
                        Temp = 0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(OrbPair) = PairEnergies(OrbPair) + Temp

                enddo

            enddo
        enddo
        
        write(6,"(A)") "Energy breakdown from B term 3: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from B term 3: ",TermE
        Energy_F12 = Energy_F12 + TermE
        EnergyTemp(:) = 0.0_dp
        
!        deallocate(r_qr,Antisym_qr)
!        allocate(r_qr(tntmo_sq),stat=ierr)
!        allocate(Antisym_qr(taux_nmo),stat=ierr)
!        if(ierr.ne.0) call stop_all(t_r,"Alloc err")
        
        !TERM 2: 
        do v=1,tnmo
            if(IsOrbFrz(v)) cycle
            do w=1,tnmo
                if(IsOrbFrz(w)) cycle

                vw_pair = FindSqPairInd(v,w,tnmo)   !Now only indexed over orbital pair space  
                !Load into memory all rq, over all space, *q fast*
                read(unit_r_read_oocc,rec=vw_pair) (TotTotMat(loop),loop=1,tntmo_sq)

                do ispin=1,3    !Loop over all spins

                    !Antisymmetrise integrals as appropriate, and convert to proper ranges
                    call Antisym1Cusp_ooao_from_oocc(TotTotMat,OrbCabsMat,ispin,.false.)
                    !We now have (r_vw)^(a'q) with a' fast

!Swap around the order of the multiplication, so that there is no need to transpose S afterwards, and p is fast
                    call DGEMM('N','N',tnmo,tnmo,tnxmo,1.0_dp,GenFockOrbCABS(1,tnmo+1),tnmo,OrbCabsMat,tnxmo,0.0_dp,S_pq,tnmo)

                    !Now load in R_TLD for the vw pair...
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,NmoNmoMat,ispin,1,tnmo,1,tnmo) 
                    
!All contributions need to be x2, since we could flip all spins.
!This would need to be changed for open shell systems.
!First, contract the same spin types of S and R_Tld
                    Temp = -0.25_dp * DDOT(tnmo_sq,NmoNmoMat,1,S_pq,1)
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                    if(ispin.eq.3) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,NmoNmoMat,2,1,tnmo,1,tnmo) 
                        Temp = 0.25_dp * DDOT(tnmo_sq,NmoNmoMat,1,S_pq,1)
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,NmoNmoMat,3,1,tnmo,1,tnmo) 
                        Temp = 0.25_dp * DDOT(tnmo_sq,NmoNmoMat,1,S_pq,1)
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,NmoNmoMat,1,1,tnmo,1,tnmo)
                        Temp = 0.25_dp * DDOT(tnmo_sq,NmoNmoMat,1,S_pq,1)
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp
                        
                enddo
            enddo
        enddo

        write(6,"(A)") "Energy breakdown from B term 2: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from term 2: ",TermE
        Energy_F12 = Energy_F12 + TermE
        EnergyTemp(:) = 0.0_dp
        

        call CalcFGTerm(Energy_F12,PairEnergies)
        call CalcTauTerm(Energy_F12,PairEnergies,.false.)
        
        close(unit_raaaa)
        close(unit_rabab)
        close(unit_rabba)

        call CalcV2Term(Energy_F12,PairEnergies)
        
        unit_raaaa=get_free_unit()
        open(unit_raaaa,file='RTld_aaaa',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)
        unit_rabab=get_free_unit()
        open(unit_rabab,file='RTld_abab',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)
        unit_rabba=get_free_unit()
        open(unit_rabba,file='RTld_abba',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)

        !Term 4 B
        !1/4((r_vw)^pa' Gam_p^q f_q^r Gam_r^s r_sa'^xy Gam^vw_xy - permute)
        !Precomute f_qr Gam_rs
        call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,GenFockOrb(:,:),tnmo,OneRDM(:,:),tnmo,0.0_dp,NmoNmoMat,tnmo)
        !Now precomute contraction with Gam_pq
        call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,OneRDM(:,:),tnmo,NmoNmoMat(:),tnmo,0.0_dp,NmoNmoMat_2,tnmo)

        do v=1,tnmo
            if(IsOrbFrz(v)) cycle
            do w=1,tnmo
                if(IsOrbFrz(w)) cycle

                vw_pair = FindSqPairInd(v,w,tnmo)   !Now only indexed over orbital pair space  
                !Load into memory all rq, over all space, *q fast*
                read(unit_r_read_oocc,rec=vw_pair) (TotTotMat(loop),loop=1,tntmo_sq)

                do ispin=1,3    !Loop over all spins

                    !Antisymmetrise integrals as appropriate, and convert to proper ranges
                    call Antisym1Cusp_ooao_from_oocc(TotTotMat,OrbCabsMat,ispin,.false.)
                    !We now have (r_vw)^(a'p) with a' fast

                    !Contract over p
                    call DGEMM('N','N',tnxmo,tnmo,tnmo,1.0_dp,OrbCabsMat,tnxmo,NmoNmoMat_2(:),tnmo,0.0_dp,OrbCabsMat_2,tnxmo)

                    !Get (Rtld^vw)_a's and (Rtld^wv)_a's
                    !a' is returned as fast here...
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,OrbCabsMat,ispin,tnmo+1,tntmo,1,tnmo)

                    !Contract over a's and sum in to energy
                    Temp = 0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                    if(ispin.eq.3) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat,2,tnmo+1,tntmo,1,tnmo)
                        Temp = -0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat,3,tnmo+1,tntmo,1,tnmo)
                        Temp = -0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat,1,tnmo+1,tntmo,1,tnmo)
                        Temp = -0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                enddo

            enddo
        enddo
        
        write(6,"(A)") "Energy breakdown from B term 4: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from B term 4: ",TermE
        Energy_F12 = Energy_F12 + TermE
        EnergyTemp(:) = 0.0_dp

        !Term B 5
        !-1/4( (r_vw)^pa' Gam_p^q f_a'^b' (rtld^vw)_qb' - permute)
        do v=1,tnmo
            if(IsOrbFrz(v)) cycle
            do w=1,tnmo
                if(IsOrbFrz(w)) cycle

                vw_pair = FindSqPairInd(v,w,tnmo)   !Now only indexed over orbital pair space  
                !Load into memory all rq, over all space, *q fast*
                read(unit_r_read_oocc,rec=vw_pair) (TotTotMat(loop),loop=1,tntmo_sq)

                do ispin=1,3    !Loop over all spins

                    !Antisymmetrise integrals as appropriate, and convert to proper ranges
                    call Antisym1Cusp_ooao_from_oocc(TotTotMat,OrbCabsMat,ispin,.true.)
                    !We now have (r_vw)^(a'p) with a' fast
                    
                    !Now contract with OneRDM over p
                    call DGEMM('N','N',tnxmo,tnmo,tnmo,1.0_dp,OrbCabsMat,tnxmo,OneRDM(:,:),tnmo,0.0_dp,OrbCabsMat_2,tnxmo)
                    !We now have (X_vw)^qa' with a' fast

                    !Now contract with f_a'^b' over a'
!                    call DGEMM('T','N',tnmo,tnxmo,tnxmo,1.0_dp,OrbCabsMat_2,tnxmo,GenFockCABS(:,:),tnxmo,0.0_dp,OrbCabsMat,tnmo)
                    call DGEMM('N','N',tnxmo,tnmo,tnxmo,1.0_dp,GenFockCABS(:,:),tnxmo,OrbCabsMat_2,tnxmo,0.0_dp,OrbCabsMat,tnxmo)
                    !We now have (X_vw)^qb' with b' fast

                    !Get (Rtld^vw)_b'q and (Rtld^wv)_b'q
                    !b' is returned as fast here...
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,OrbCabsMat_2,ispin,1,tnmo,tnmo+1,tntmo)

                    !Contract over b'q and sum in to energy
                    Temp = -0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                    if(ispin.eq.3) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat_2,2,1,tnmo,tnmo+1,tntmo)
                        Temp = 0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat_2,3,1,tnmo,tnmo+1,tntmo)
                        Temp = 0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat_2,1,1,tnmo,tnmo+1,tntmo)
                        Temp = 0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                enddo

            enddo
        enddo
        
        write(6,"(A)") "Energy breakdown from B term 5: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from B term 5: ",TermE
        Energy_F12 = Energy_F12 + TermE
        EnergyTemp(:) = 0.0_dp


        !Term 6 B
        !-1/4((r_vw)^p'a' f_p'p Gam_p^q r_qa'^xy Gam^vw_xy - permute)
        !Precomute f_p'p Gam_pq
        call DGEMM('N','N',tntmo,tnmo,tnmo,1.0_dp,GenFockCompOrb(:,:),tntmo,OneRDM(:,:),tnmo,0.0_dp,TotOrbMat,tntmo)
        !X_p'q - p' fast

!        do i=1,tntmo_nmo
!            write(27,"(G25.10,2I5)") TotOrbMat(i)
!        enddo

        do v=1,tnmo
            if(IsOrbFrz(v)) cycle
            do w=1,tnmo
                if(IsOrbFrz(w)) cycle

                vw_pair = FindSqPairInd(v,w,tnmo)   !Now only indexed over orbital pair space  
                !Load into memory all rq, over all space, *q fast*
                read(unit_r_read_oocc,rec=vw_pair) (TotTotMat(loop),loop=1,tntmo_sq)

                do ispin=1,3    !Loop over all spins

                    !Antisymmetrise integrals as appropriate, and convert to proper ranges
                    call Antisym1Cusp_ooca_from_oocc(TotTotMat,TotCabsMat,ispin)
                    !We now have (r_vw)^(p'a') with p' fast
    
                    !Contract over p'
                    call DGEMM('T','N',tnxmo,tnmo,tntmo,1.0_dp,TotCabsMat,tntmo,TotOrbMat,tntmo,0.0_dp,OrbCabsMat_2,tnxmo)
                    !Returns (U_vw)^a'q with a' fast

                    !Get (Rtld^vw)_a'q and (Rtld^wv)_a'q
                    !a' is returned as fast here...
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,OrbCabsMat,ispin,1,tnmo,tnmo+1,tntmo)

                    !Contract over a's and sum in to energy
                    Temp = -0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                    if(ispin.eq.3) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat,2,1,tnmo,tnmo+1,tntmo)
                        Temp = 0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat,3,1,tnmo,tnmo+1,tntmo)
                        Temp = 0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,OrbCabsMat,1,1,tnmo,tnmo+1,tntmo)
                        Temp = 0.25_dp * DDOT(taux_nmo,OrbCabsMat_2,1,OrbCabsMat,1)
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                enddo

            enddo
        enddo
        
        write(6,"(A)") "Energy breakdown from B term 6: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from B term 6: ",TermE
        Energy_F12 = Energy_F12 + TermE
        EnergyTemp(:) = 0.0_dp
        
        !Term 7 B
        !-1/4((r_vw)^pa' Gam_p^q f_q^p' r_p'a'^xy Gam^vw_xy - permute)
        !Precomute Gam_pq f_qp'
        call DGEMM('N','T',tnmo,tntmo,tnmo,1.0_dp,OneRDM(:,:),tnmo,GenFockCompOrb(:,:),tntmo,0.0_dp,TotOrbMat,tnmo)
        !X_pp' - p fast

        do v=1,tnmo
            if(IsOrbFrz(v)) cycle
            do w=1,tnmo
                if(IsOrbFrz(w)) cycle

                vw_pair = FindSqPairInd(v,w,tnmo)   !Now only indexed over orbital pair space  
                !Load into memory all rq, over all space, *q fast*
                read(unit_r_read_oocc,rec=vw_pair) (TotTotMat(loop),loop=1,tntmo_sq)

                do ispin=1,3    !Loop over all spins

                    !Antisymmetrise integrals as appropriate, and convert to proper ranges
                    call Antisym1Cusp_ooao_from_oocc(TotTotMat,OrbCabsMat,ispin,.false.)
                    !We now have (r_vw)^(a'p) with a' fast
    
                    !Contract over p
                    call DGEMM('N','N',tnxmo,tntmo,tnmo,1.0_dp,OrbCabsMat,tnxmo,TotOrbMat,tnmo,0.0_dp,TotCabsMat,tnxmo)
                    !Returns (U_vw)^a'p' with a' fast

                    !Get (Rtld^vw)_p'a' and (Rtld^wv)_p'a'
                    !a' is returned as fast here...
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,TotCabsMat_2,ispin,1,tntmo,tnmo+1,tntmo)

                    !Contract over a's and sum in to energy
                    Temp = -0.25_dp * DDOT(tntmo_nxmo,TotCabsMat_2,1,TotCabsMat,1)
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                    if(ispin.eq.3) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,TotCabsMat_2,2,1,tntmo,tnmo+1,tntmo)
                        Temp = 0.25_dp * DDOT(tntmo_nxmo,TotCabsMat_2,1,TotCabsMat,1)
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,TotCabsMat_2,3,1,tntmo,tnmo+1,tntmo)
                        Temp = 0.25_dp * DDOT(tntmo_nxmo,TotCabsMat_2,1,TotCabsMat,1)
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,TotCabsMat_2,1,1,tntmo,tnmo+1,tntmo)
                        Temp = 0.25_dp * DDOT(tntmo_nxmo,TotCabsMat_2,1,TotCabsMat,1)
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                enddo

            enddo
        enddo
        
        write(6,"(A)") "Energy breakdown from B term 7: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from B term 7: ",TermE
        Energy_F12 = Energy_F12 + TermE
        EnergyTemp(:) = 0.0_dp


        !Term 13 B: This can be broken up into other terms if speedup is needed to avoid contracting over tntmo**2
        !-1/4((r_vw)^p'q' K_p'^r' r_r'q'^xy Gam^vw_xy - permute)
        do v=1,tnmo
            if(IsOrbFrz(v)) cycle
            do w=1,tnmo
                if(IsOrbFrz(w)) cycle

                vw_pair = FindSqPairInd(v,w,tnmo)   !Now only indexed over orbital pair space  
                !Load into memory all rq, over all space, *q fast*
                read(unit_r_read_oocc,rec=vw_pair) (TotTotMat_3(loop),loop=1,tntmo_sq)

                do ispin=1,3    !Loop over all spins

                    !Antisymmetrise integrals as appropriate, and convert to proper ranges
                    call Antisym1Cusp(TotTotMat_3,TotTotMat,tntmo,ispin)
                    !We now have (r_vw)^(p'q') with q' fast
    
                    !Contract over p'
                    call DGEMM('T','T',tntmo,tntmo,tntmo,1.0_dp,GenExch,tntmo,TotTotMat,tntmo,0.0_dp,TotTotMat_2,tntmo)
                    !The transposes are required to ensure that we end up with r' fast, so it is the right way round for r_tld
                    !Returns (X_vw)^r'q' with r' fast

                    !Get (Rtld^vw)_r'q' and (Rtld^wv)_r'q'
                    !r' is returned as fast here...
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,TotTotMat,ispin,1,tntmo,1,tntmo)

                    !Contract over r'q' and sum in to energy
                    Temp = -0.25_dp * DDOT(tntmo_sq,TotTotMat_2,1,TotTotMat,1)
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                    if(ispin.eq.3) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,TotTotMat,2,1,tntmo,1,tntmo)
                        Temp = 0.25_dp * DDOT(tntmo_sq,TotTotMat_2,1,TotTotMat,1)
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,TotTotMat,3,1,tntmo,1,tntmo)
                        Temp = 0.25_dp * DDOT(tntmo_sq,TotTotMat_2,1,TotTotMat,1)
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,w,v,TotTotMat,1,1,tntmo,1,tntmo)
                        Temp = 0.25_dp * DDOT(tntmo_sq,TotTotMat_2,1,TotTotMat,1)
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                enddo

            enddo
        enddo
        
        write(6,"(A)") "Energy breakdown from B term 13: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from B term 13: ",TermE
        Energy_F12 = Energy_F12 + TermE
        EnergyTemp(:) = 0.0_dp

        call CalcRRTerms(Energy_F12,PairEnergies,EnergyTemp_FTF,.false.)

        !Term 2 of X:
        !  -1/4 [ (r_vw)^ta' G_t^u r_ua'^xy phi_xy^vw - permute ]
        do v=1,tnmo
            if(IsOrbFrz(v)) cycle
            do w=1,tnmo
                if(IsOrbFrz(w)) cycle

                vw_pair = FindSqPairInd(v,w,tnmo)   !Now only indexed over orbital pair space  
                !Load into memory all rq, over all space, *q fast*
                read(unit_r_read_oocc,rec=vw_pair) (TotTotMat(loop),loop=1,tntmo_sq)

                do ispin=1,3    !Loop over all spins

                    !Antisymmetrise integrals as appropriate, and convert to proper ranges
                    call Antisym1Cusp_ooao_from_oocc(TotTotMat,OrbCabsMat,ispin,.true.)
                    !We now have (r_vw)^(a't) with a' fast

                    !Contract over t
                    call DGEMM('N','N',tnxmo,tnmo,tnmo,1.0_dp,OrbCabsMat,tnxmo,OneRDM,tnmo,0.0_dp,OrbCabsMat_2,tnxmo)
                    !(X_vw)^a'u with a' fast.

                    !Read RTld_Phi from disk
                    !a' is returned as fast here...
                    call GetAntisymR_TLD(unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba,v,w,OrbCabsMat,ispin, &
                        1,tnmo,tnmo+1,tntmo)

                    !Contract over ua' and sum in to energy
                    Temp = 0.25_dp * DDOT(taux_nmo,OrbCabsMat,1,OrbCabsMat_2,1)
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp

                    if(ispin.eq.3) then
                        call GetAntisymR_TLD(unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba,w,v,OrbCabsMat,2, &
                            1,tnmo,tnmo+1,tntmo)
                        Temp = -0.25_dp * DDOT(taux_nmo,OrbCabsMat,1,OrbCabsMat_2,1)
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba,w,v,OrbCabsMat,3, &
                            1,tnmo,tnmo+1,tntmo)
                        Temp = -0.25_dp * DDOT(taux_nmo,OrbCabsMat,1,OrbCabsMat_2,1)
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba,w,v,OrbCabsMat,1, &
                            1,tnmo,tnmo+1,tntmo)
                        Temp = -0.25_dp * DDOT(taux_nmo,OrbCabsMat,1,OrbCabsMat_2,1)
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(vw_pair) = PairEnergies(vw_pair) + Temp
                enddo
            enddo
        enddo

        write(6,"(A)") "Energy breakdown from X term 2: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from X term 2: ",TermE
        Energy_F12 = Energy_F12 + TermE
        EnergyTemp(:) = 0.0_dp


        !TERM 3 in X:
        !-1/8 [ (r_vw)^tu r^xy_tu Phi^vw_xy - permute ]
        do v=1,tnmo
            if(IsOrbFrz(v)) cycle
            do w=1,tnmo
                if(IsOrbFrz(w)) cycle

!Since these are stored for r integrals over the complete space (but not used here, we need to index over the complete space)
                vw_pair = FindSqPairInd(v,w,tntmo)  
                OrbPair = FindSqPairInd(v,w,tnmo)
                !Load into memory all rq
                read(unit_r_read_ccoo,rec=vw_pair) (NmoNmoMat_3(loop),loop=1,tnmo_sq)   !q is fast index (opposite to matrix element)
                call TransposeIntegralArr_inplace(NmoNmoMat_3,tnmo)    !Ensure t fast

                do ispin=1,3    !Loop over all spins

                    !Antisymmetrise integrals as appropriate
                    call Antisym1Cusp(NmoNmoMat_3,NmoNmoMat_2,tnmo,ispin)
                    
                    !Read RTld_Phi from disk
                    call GetAntisymR_TLD(unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba,v,w,NmoNmoMat,ispin,1,tnmo,1,tnmo)

                    !Contract over ua' and sum in to energy
                    Temp = DDOT(tnmo_sq,NmoNmoMat_2,1,NmoNmoMat,1) / 8.0_dp
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(OrbPair) = PairEnergies(OrbPair) + Temp

                    if(ispin.eq.3) then
                        call GetAntisymR_TLD(unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba,w,v,NmoNmoMat,2,1,tnmo,1,tnmo)
                        Temp = -DDOT(tnmo_sq,NmoNmoMat_2,1,NmoNmoMat,1) / 8.0_dp
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba,w,v,NmoNmoMat,3,1,tnmo,1,tnmo)
                        Temp = -DDOT(tnmo_sq,NmoNmoMat_2,1,NmoNmoMat,1) / 8.0_dp
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_rtld_Phi_aaaa,unit_rtld_Phi_abab,unit_rtld_Phi_abba,w,v,NmoNmoMat,1,1,tnmo,1,tnmo)
                        Temp = -DDOT(tnmo_sq,NmoNmoMat_2,1,NmoNmoMat,1) / 8.0_dp
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(OrbPair) = PairEnergies(OrbPair) + Temp
                enddo
            enddo
        enddo

        write(6,"(A)") "Energy breakdown from X term 3: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from X term 3: ",TermE
        Energy_F12 = Energy_F12 + TermE
        EnergyTemp(:) = 0.0_dp
        
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
        deallocate(r_qr)
        deallocate(Antisym_qr)
        deallocate(S_pq)
        deallocate(g_qr)
        deallocate(Antisym_g_qr)
        deallocate(OrbCabsMat)
        deallocate(OrbCabsMat_2)
        deallocate(NmoNmoMat)
        deallocate(NmoNmoMat_2)
        deallocate(NmoNmoMat_3)
        deallocate(TotOrbMat)
        deallocate(TotCabsMat)
        deallocate(TotCabsMat_2)
        deallocate(TotTotMat)
        deallocate(TotTotMat_2)
        deallocate(TotTotMat_3)

        close(unit_g_read,status='delete')
        close(unit_raaaa,status='delete')
        close(unit_rabab,status='delete')
        close(unit_rabba,status='delete')
        close(unit_rtld_Phi_aaaa,status='delete')
        close(unit_rtld_Phi_abab,status='delete')
        close(unit_rtld_Phi_abba,status='delete')
        close(unit_r_read_ccoo,status='delete')
        close(unit_r_read_oocc,status='delete')

    end subroutine CalcF12Correction

    subroutine CalcV2Term(Energy_F12,PairEnergies)
        implicit none
        real(dp) , intent(inout) :: Energy_F12,PairEnergies(tnmo_sq)
        real(dp) :: EnergyTemp(3),TermE,DDOT,Temp
        integer :: r,s,unit_raaaa,unit_rabab,unit_rabba,unit_g_oocc
        integer :: rsPair,loop,ispin,t,a_prime,i,ierr
        real(dp), allocatable :: r_temp(:),r_ta(:),S_ua(:),RTld_ua(:)
        character(*), parameter :: t_r='CalcV2Term'
        
        EnergyTemp(:) = 0.0_dp

        unit_raaaa=get_free_unit()
        open(unit_raaaa,file='RTld_aaaa',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)
        unit_rabab=get_free_unit()
        open(unit_rabab,file='RTld_abab',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)
        unit_rabba=get_free_unit()
        open(unit_rabba,file='RTld_abba',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)

        unit_g_oocc=get_free_unit()
        open(unit_g_oocc,file='G_Dir_oocc',status='old',form='unformatted',access='direct',action='read',recl=reclen_ntmo_sq)

        allocate(r_temp(tntmo_sq),stat=ierr)
        allocate(r_ta(taux_nmo),stat=ierr)
        allocate(S_ua(taux_nmo),stat=ierr)
        allocate(RTld_ua(taux_nmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc err")

        !Now for V term 2...
        do r=1,tnmo
            do s=1,tnmo

                rsPair = FindSqPairInd(r,s,tnmo)
                
                !Get (g_rs)^ta'
                read(unit_g_oocc,rec=rspair) (r_temp(loop),loop=1,tntmo_sq)

                do ispin=1,3

                    !Antisymmetrise over t a' (no geminal indices)
                    call Antisym0Cusp_ooao_from_oocc(r_temp,r_ta,ispin)
!                    r_ta(:) = 0.0_dp
!                    do t=1,tnmo
!                        do a_prime=tnmo+1,tntmo
!                            !a_prime is fast in r_ta
!                            if(ispin.eq.1) then
!                                r_ta(FindCABSOrbPairInd(a_prime,t)) = r_temp(FindSqPairInd(t,a_prime,tntmo)) &
!                                    - r_temp(FindSqPairInd(a_prime,t,tntmo))
!                            elseif(ispin.eq.2) then
!                                r_ta(FindCABSOrbPairInd(a_prime,t)) = r_temp(FindSqPairInd(t,a_prime,tntmo))
!                            elseif(ispin.eq.3) then
!                                r_ta(FindCABSOrbPairInd(a_prime,t)) = - r_temp(FindSqPairInd(a_prime,t,tntmo))
!                            endif
!                        enddo
!                    enddo

                    !Contract g with 1RDM_t^u
    !                call DGEMM('N','N',tnmo,tnmo,tnxmo,1.0_dp,GenFockOrbCABS(1,tnmo+1),tnmo,Antisym_qr(1),tnxmo,0.0_dp,S_pq,tnmo)
                    call DGEMM('N','N',tnxmo,tnmo,tnmo,1.0_dp,r_ta,tnxmo,OneRDM,tnmo,0.0_dp,S_ua,tnxmo)
    !               S_ua' - a' fast.

                    !Get (Rtld^rs)_ua' and (Rtld^sr)_ua'
                    !a' is returned as fast here...
                    call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,r,s,RTld_ua,ispin,1,tnmo,tnmo+1,tntmo)

                    !Contract over ua' and sum in to energy
                    Temp = -0.5_dp * DDOT(taux_nmo,RTld_ua,1,S_ua,1)
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(rsPair) = PairEnergies(rsPair) + Temp

                    if(ispin.eq.3) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,s,r,RTld_ua,2,1,tnmo,tnmo+1,tntmo)
                        Temp = 0.5_dp * DDOT(taux_nmo,RTld_ua,1,S_ua,1)
                    elseif(ispin.eq.2) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,s,r,RTld_ua,3,1,tnmo,tnmo+1,tntmo)
                        Temp = 0.5_dp * DDOT(taux_nmo,RTld_ua,1,S_ua,1)
                    elseif(ispin.eq.1) then
                        call GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,s,r,RTld_ua,1,1,tnmo,tnmo+1,tntmo)
                        Temp = 0.5_dp * DDOT(taux_nmo,RTld_ua,1,S_ua,1)
                    endif
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    PairEnergies(rsPair) = PairEnergies(rsPair) + Temp

                enddo

            enddo
        enddo

        close(unit_raaaa)
        close(unit_rabab)
        close(unit_rabba)
        close(unit_g_oocc,status='delete')
        
        write(6,"(A)") "Energy breakdown from V term 2: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from V term 2: ",TermE
        Energy_F12 = Energy_F12 + TermE

        deallocate(RTld_ua,S_ua,r_ta,r_temp)

    end subroutine CalcV2Term
    
    !1/2 (fpk_v^p' r2_p'w^xy Gamma_xy^vw - permute) : This is the FpK term
    !Also, the term in X: 1/4[ (RR_vw)^xy Gamma_xy^vw - permute ]
    subroutine CalcRRTerms(Energy_F12,PairEnergies,EnergyTemp_FTF,tUseFTFInts)
        use input_data, only: gamma_length
        implicit none
        real(dp) , intent(inout) :: Energy_F12,PairEnergies(tnmo_sq)
        real(dp) :: EnergyTemp(3),TermE,DDOT,EnergyTemp_X(3),Temp,Z,EnergyTemp_FTF(2)
        character(*), parameter :: t_r='CalcRRTerms'
        integer :: v,w,unit_2RDM_read_aaaa,unit_2RDM_read_abab,unit_2RDM_read_abba
        integer :: unit_rr,ierr,x,y,ios,ibuf,i_Ind,j_Ind,k_Ind,l_Ind,ispin,j,ijPair,jiPair
        integer :: p_prime,counter,xy_pair,wInd,vInd,loop,i,vwPair,wvPair
        real(dp) , allocatable :: rrints(:,:),Buf(:),CompOrbMat(:),NmoNmoMat(:),Gamma_vw(:),Gamma_wv(:)
        real(dp) , allocatable :: Phi_xy(:),OrbOrbMat(:)
        integer(i2) , allocatable :: Indices(:,:)
        integer(i8) :: maxlength,length
        integer :: unit_Phi_read_aaaa,unit_Phi_read_abab,unit_Phi_read_abba,vw_pair
        logical, intent(in) :: tUseFTFInts
        logical :: exists

        if(tUseFTFInts) then
            inquire(file='FTFDUMPBIN',exist=exists)
            if(.not.exists) return
        endif
        
        EnergyTemp(:) = 0.0_dp
        EnergyTemp_X(:) = 0.0_dp

        if(.not.tHFRDMs) then
            unit_2RDM_read_aaaa=get_free_unit()
            open(unit_2RDM_read_aaaa,file='TwoRDM_Dir_aaaa',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
            unit_2RDM_read_abab=get_free_unit()
            open(unit_2RDM_read_abab,file='TwoRDM_Dir_abab',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
            unit_2RDM_read_abba=get_free_unit()
            open(unit_2RDM_read_abba,file='TwoRDM_Dir_abba',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
        endif

        !Open files for reading in the Phi term
        unit_Phi_read_aaaa=get_free_unit()
        open(unit_Phi_read_aaaa,file='Phi_Dir_aaaa',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)
        unit_Phi_read_abab=get_free_unit()
        open(unit_Phi_read_abab,file='Phi_Dir_abab',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)
        unit_Phi_read_abba=get_free_unit()
        open(unit_Phi_read_abba,file='Phi_Dir_abba',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)

        !rrints stores all p',q',y in <p' q' | X y>
        allocate(rrints(tntmo_sq,tnmo),stat=ierr) !indices {p',q'},y
        if(ierr.ne.0) call stop_all(t_r,"Alloc err")

        allocate(CompOrbMat(tntmo_nmo),stat=ierr)
        allocate(NmoNmoMat(tnmo_sq),stat=ierr)
        allocate(Gamma_vw(tnmo_sq),stat=ierr)
        allocate(Gamma_wv(tnmo_sq),stat=ierr)
        allocate(OrbOrbMat(tnmo_sq),stat=ierr)
        allocate(Phi_xy(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc err")
        
        unit_rr = get_free_unit()
        if(tDaltonFormat) then
            open(unit_rr,file='MO_FF',status='old',form='unformatted',access='sequential',action='read')
        else
            if(tReadBin) then
                if(tUseFTFInts) then
                    open(unit_rr,file='FTFDUMPBIN',status='old',form='unformatted',action='read')
                else
                    open(unit_rr,file='FFDUMPBIN',status='old',form='unformatted',action='read')
                endif
            else
                open(unit_rr,file='FFDUMP',status='old',form='formatted',action='read')
            endif
        endif
        rewind(unit_rr)
        if(tDaltonFormat) then
            read(unit_rr) maxlength
            allocate(Indices(4,MaxLength),Buf(MaxLength),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err")
        else
            length = 1
        endif

        do x=1,tnmo !Loop over all x space
            !If doing MP2F12, only need to loop over occupied pairs here...
            if(tHFRDMs.and.(.not.IsOrbOcc(x))) cycle
            if(IsOrbFrz(x)) cycle

            rewind(unit_rr)
            rrints(:,:) = 0.0_dp
            if(tDaltonFormat) read(unit_rr) maxlength

            do while(.true.)
                if(tDaltonFormat) then
                    read(unit_rr,iostat=ios) length,Indices(1:4,1:length),Buf(1:length)
                else
                    if(tReadBin) then
                        read(unit_rr,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                    else
                        read(unit_rr,*,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                    endif
                endif
                if(ios.gt.0) call stop_all(t_r,"Error reading ints")
                if((length.le.0).or.(ios.lt.0)) exit

                do ibuf = 1,length
                    if(tDaltonFormat) then
                        i_Ind = int(Indices(1,ibuf),i4)
                        k_Ind = int(Indices(2,ibuf),i4)
                        j_Ind = int(Indices(3,ibuf),i4)
                        l_Ind = int(Indices(4,ibuf),i4)
                        Z = buf(ibuf)
!                        if(i_Ind.eq.1.and.i_Ind.eq.j_Ind.and.j_Ind.eq.k_Ind.and.k_Ind.eq.l_Ind) then
!                            write(6,*) "FF: ",i_Ind,j_Ind,k_Ind,l_Ind,Z
!                        endif
                    endif

                    !Check if at least one of the orbitals is 'x'
                    if((i_Ind.ne.x).and.(j_Ind.ne.x).and.(k_Ind.ne.x).and.(l_Ind.ne.x)) cycle
                    !Check if more than two indices outside orbital space
                    if(CountCabsInd(i_Ind,j_Ind,k_Ind,l_Ind).gt.2) cycle

                    if(i_Ind.eq.x) then
                        if(j_ind.le.tnmo) then
                            !<X y | p' q'>
                            rrints(FindSqPairInd(k_ind,l_ind,tntmo),j_ind) = Z 
                        endif
                        if(.not.tComplexOrbs) then
                            if(l_Ind.le.tnmo) then
                                !<X q' | p' y>
                                rrints(FindSqPairInd(k_Ind,j_Ind,tntmo),l_Ind) = Z 
                            endif
                        endif
                    endif
                    if(j_Ind.eq.x) then
                        if(i_ind.le.tnmo) then
                            !<y X | q' p'>
                            rrints(FindSqPairInd(l_ind,k_ind,tntmo),i_ind) = Z
                        endif
                        if(.not.tComplexOrbs) then
                            if(k_ind.le.tnmo) then
                                !<q' X | y p'>
                                rrints(FindSqPairInd(l_ind,i_ind,tntmo),k_ind) = Z
                            endif
                        endif
                    endif
                    if(k_Ind.eq.x) then
                        if(l_ind.le.tnmo) then
                            !<p' q' | X y>
                            rrints(FindSqPairInd(i_Ind,j_ind,tntmo),l_ind) = Z
                        endif
                        if(.not.tComplexOrbs) then
                            if(j_ind.le.tnmo) then
                                !<p' y | X q'>
                                rrints(FindSqPairInd(i_ind,l_ind,tntmo),j_ind) = Z
                            endif
                        endif
                    endif
                    if(l_Ind.eq.x) then
                        if(k_ind.le.tnmo) then
                            !<q' p' | y X>
                            rrints(FindSqPairInd(j_ind,i_ind,tntmo),k_ind) = Z
                        endif
                        if(.not.tComplexOrbs) then
                            if(i_ind.le.tnmo) then
                                !<y p' | q' X>
                                rrints(FindSqPairInd(j_ind,k_ind,tntmo),i_ind) = Z
                            endif
                        endif
                    endif

                enddo   !End loop over ibuf
            enddo
            !rrints now holds (rr^x)_p'q'^y  with order {p',q'},y  with q' fast in the first pair

            do y=1,tnmo
                !If doing MP2F12, only need to loop over occupied pairs here...
                if(tHFRDMs.and.(.not.IsOrbOcc(y))) cycle
                if(IsOrbFrz(y)) cycle

                do ispin=1,3

                    if(tCalcCombo.and.ispin.eq.3) cycle !Do not calculate abba term if tCalcCombo=.true.

                    !Construct the appropriatly antisymmetrised ((rr)^xy)_p'w
                    !Two Geminal pairs!
                    counter = 1
                    do w=1,tnmo
                        do p_prime=1,tntmo
                            !We want p_prime fast
                            if(ispin.eq.1) then
                                CompOrbMat(counter) = (rrints(FindSqPairInd(p_prime,w,tntmo),y) - &
                                    rrints(FindSqPairInd(w,p_prime,tntmo),y))/16.0_dp
                            elseif(ispin.eq.2) then
                                CompOrbMat(counter) = rrints(FindSqPairInd(p_prime,w,tntmo),y)*(5.0_dp/32.0_dp) + &
                                    rrints(FindSqPairInd(w,p_prime,tntmo),y)*(3.0_dp/32.0_dp)
                            else
                                CompOrbMat(counter) = rrints(FindSqPairInd(p_prime,w,tntmo),y)*(-3.0_dp/32.0_dp) - &
                                    rrints(FindSqPairInd(w,p_prime,tntmo),y)*(5.0_dp/32.0_dp)
!                                write(27,"(G20.10,4I5)") CompOrbMat(counter),p_prime,w,x,y
                            endif

                            !Also seperate out the antisymmetrised rr integrals for the X term
                            !Here, we only want antisymmetrised integrals over orbital:orbital space.
                            if(p_prime.le.tnmo) then
                                vw_pair = FindSqPairInd(p_prime,w,tnmo)
                                OrbOrbMat(vw_pair) = CompOrbMat(counter)
!                                write(6,*) OrbOrbMat(vw_pair),w,p_prime,counter
                            endif

                            counter = counter + 1
                        enddo
                    enddo

                    !Store (Phi_xy)^vw for the X term
                    xy_pair = FindSqPairInd(x,y,tnmo)
                    if(ispin.eq.1) then
                        read(unit_Phi_read_aaaa,rec=xy_pair) (Phi_xy(loop),loop=1,tnmo_sq)
                    elseif(ispin.eq.2) then
                        read(unit_Phi_read_abab,rec=xy_pair) (Phi_xy(loop),loop=1,tnmo_sq)
                    elseif(ispin.eq.3) then
                        read(unit_Phi_read_abba,rec=xy_pair) (Phi_xy(loop),loop=1,tnmo_sq)
                    endif
                    !Ensure that (Phi_xy)^vw only has vw over active pairs. Set other elements to zero.
                    !The xy pairs are already only read in over active pairs.
                    do i=1,tnfrz
                        do j=1,tnmo
                            ijPair = FindSqPairInd(FrzOrbs(i),j,tnmo)
                            jiPair = FindSqPairInd(j,FrzOrbs(i),tnmo)
                            Phi_xy(ijPair) = 0.0_dp
                            Phi_xy(jiPair) = 0.0_dp
                        enddo
                    enddo

                    !Contract with FpK over p'
                    call DGEMM('T','N',tnmo,tnmo,tntmo,1.0_dp,GenFpK_CompOrb,tntmo,CompOrbMat,tntmo,0.0_dp,NmoNmoMat,tnmo)
                    !NmoNmoMat now contains (X^xy)_vw with v fast

                    !Calculate (Gamma_xy)^vw with v fast
                    xy_pair = FindSqPairInd(x,y,tnmo)
                    if(tHFRDMs) then
                        Gamma_vw(:)=0.0_dp
                        Gamma_wv(:)=0.0_dp
                        do v=1,tnocc
                            vInd = OccOrbs(v)
                            if(IsOrbFrz(vInd)) cycle

                            do w=1,tnocc
                                wInd = OccOrbs(w)
                                if(IsOrbFrz(wInd)) cycle

                                vwPair = FindSqPairInd(wInd,vInd,tnmo)  !v fast
                                wvPair = FindSqPairInd(vInd,wInd,tnmo)  !w fast
                                if((x.eq.vInd).and.(y.eq.wInd).and.(x.eq.y)) then
                                    !All four indices the same
                                    if(ispin.eq.2) then
                                        Gamma_vw(vwPair) = 1.0_dp
                                        Gamma_wv(vwPair) = -1.0_dp
                                    elseif(ispin.eq.3) then
                                        Gamma_vw(vwPair) = -1.0_dp
                                        Gamma_wv(vwPair) = 1.0_dp
                                    endif
                                elseif((x.eq.vInd).and.(y.eq.wInd)) then
                                    if(ispin.eq.1) then
                                        Gamma_vw(vwPair) = 1.0_dp
                                    elseif(ispin.eq.2) then
                                        Gamma_vw(vwPair) = 1.0_dp
                                        Gamma_wv(vwPair) = -1.0_dp
                                    endif
                                elseif((y.eq.vInd).and.(x.eq.wInd)) then
                                    if(ispin.eq.1) then
                                        Gamma_vw(vwPair) = -1.0_dp
                                    elseif(ispin.eq.3) then
                                        Gamma_vw(vwPair) = -1.0_dp
                                        Gamma_wv(vwPair) = 1.0_dp
                                    endif
                                endif
                            enddo
                        enddo
                        if((ispin.eq.1).and.(.not.tCalcCombo)) then
                            !Calculate wv just by transposing vw
                            call TransposeIntegralArr(Gamma_vw(:),Gamma_wv(:),tnmo) !Find wv just by transposing vw
                        endif
                    else
                        !w is fast when reading in
                        if(ispin.eq.1) then
                            if(tCalcCombo) then
                                read(unit_2RDM_read_aaaa,rec=xy_pair) (Gamma_vw(loop),loop=1,tnmo_sq)
                                call TransposeIntegralArr_inplace(Gamma_vw,tnmo)
                            else
                                read(unit_2RDM_read_aaaa,rec=xy_pair) (Gamma_wv(loop),loop=1,tnmo_sq)
                                call TransposeIntegralArr(Gamma_wv,Gamma_vw,tnmo)
                            endif
                        elseif(ispin.eq.2) then
                            read(unit_2RDM_read_abab,rec=xy_pair) (Gamma_vw(loop),loop=1,tnmo_sq)
                            if(.not.tCalcCombo) read(unit_2RDM_read_abba,rec=xy_pair) (Gamma_wv(loop),loop=1,tnmo_sq)
                            call TransposeIntegralArr_inplace(Gamma_vw(:),tnmo) !We want v fast
                        elseif(ispin.eq.3) then
                            read(unit_2RDM_read_abba,rec=xy_pair) (Gamma_vw(loop),loop=1,tnmo_sq)
                            read(unit_2RDM_read_abab,rec=xy_pair) (Gamma_wv(loop),loop=1,tnmo_sq)
                            call TransposeIntegralArr_inplace(Gamma_vw(:),tnmo) !We want v fast
                        endif

                        !Ensure that (G_xy)^vw only has vw over active pairs. Set other elements to zero.
                        do i=1,tnfrz
                            do j=1,tnmo
                                ijPair = FindSqPairInd(FrzOrbs(i),j,tnmo)
                                jiPair = FindSqPairInd(j,FrzOrbs(i),tnmo)
                                Gamma_vw(ijPair) = 0.0_dp
                                Gamma_vw(jiPair) = 0.0_dp
                                Gamma_wv(ijPair) = 0.0_dp
                                Gamma_wv(jiPair) = 0.0_dp
                            enddo
                        enddo

                    endif
                    
                    !Contract over vw and sum in to energy

                    Temp = 0.5_dp * DDOT(tnmo_sq,Gamma_vw,1,NmoNmoMat,1)
                    if(tCalcCombo) Temp = Temp * 2.0_dp
                    EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                    if(tExplicitCommTerm) then
                        if(.not.tCalcCombo) call stop_all(t_r,"tExplicitCommTerm should not be set without tCalcCombo")
!                        EnergyTemp_FTF(ispin) = EnergyTemp_FTF(ispin) + Temp
                    endif
                    if(.not.tUseFTFInts) PairEnergies(xy_pair) = PairEnergies(xy_pair) + Temp

                    if(.not.tCalcCombo) then
                        if(ispin.eq.3) then
                            Temp = -0.5_dp * DDOT(tnmo_sq,Gamma_wv,1,NmoNmoMat,1)
                        elseif(ispin.eq.2) then
                            Temp = -0.5_dp * DDOT(tnmo_sq,Gamma_wv,1,NmoNmoMat,1)
                        elseif(ispin.eq.1) then
                            Temp = -0.5_dp * DDOT(tnmo_sq,Gamma_wv,1,NmoNmoMat,1)
                        endif
                        EnergyTemp(ispin) = EnergyTemp(ispin) + Temp 
                        if(.not.tUseFTFInts) PairEnergies(xy_pair) = PairEnergies(xy_pair) + Temp
                    endif

                    Temp = -0.5_dp * DDOT(tnmo_sq,OrbOrbMat,1,Phi_xy,1)
                    EnergyTemp_X(ispin) = EnergyTemp_X(ispin) + Temp 
                    if(.not.tUseFTFInts) PairEnergies(xy_pair) = PairEnergies(xy_pair) + Temp

                enddo   !End do over spins

            enddo   !End do running over y

        enddo   !Enddo over x

        if(tCalcCombo) then
            EnergyTemp(3) = EnergyTemp(2)
            EnergyTemp_X(3) = EnergyTemp_X(2)
        endif
        
        if(tUseFTFInts) then
            write(6,"(A)") "Energy breakdown from B term 9-12 with FTF Ints: "
            TermE = 0.0_dp
            do i=1,3
                TermE = TermE + EnergyTemp(i)
                write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
            enddo
            write(6,"(A,F25.12)") "Total energy from B term 9-12 with FTF Ints: ",TermE
        else
            write(6,"(A)") "Energy breakdown from B term 9-12 with FF Ints: "
            TermE = 0.0_dp
            do i=1,3
                TermE = TermE + EnergyTemp(i)
                write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
            enddo
            write(6,"(A,F25.12)") "Total energy from B term 9-12 with FF Ints: ",TermE
            write(6,"(A)") "(It is this contribution which is added to the final energy)"
            Energy_F12 = Energy_F12 + TermE
        endif
        
        if(tUseFTFInts) then
            write(6,"(A)") "Energy breakdown from X term 1 with FTF Ints: "
            TermE = 0.0_dp
            do i=1,3
                TermE = TermE + EnergyTemp_X(i)
                write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp_X(i)
            enddo
            write(6,"(A,F25.12)") "Total energy from X term 1 with FTF Ints: ",TermE
        else
            write(6,"(A)") "Energy breakdown from X term 1 with FF Ints: "
            TermE = 0.0_dp
            do i=1,3
                TermE = TermE + EnergyTemp_X(i)
                write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp_X(i)
            enddo
            write(6,"(A,F25.12)") "Total energy from X term 1 with FF Ints: ",TermE
            write(6,"(A)") "(It is this contribution which is added to the final energy)"
            Energy_F12 = Energy_F12 + TermE
        endif

        if(tExplicitCommTerm) then
            write(6,"(A)")        "Expected approximate energy breakdown from tau term: "
            write(6,"(A)")        "Warning: V slow CABS convergence, and energies not included in final result."
            write(6,"(A,F25.12)") "Same-spin energy contribution:       ",EnergyTemp_FTF(1)
            write(6,"(A,F25.12)") "Opposite-spin energy contribution:   ",EnergyTemp_FTF(2)*2.0_dp  !*2 for abba term
            write(6,"(A,F25.12)") "Total energy from approximate tau terms:  ",EnergyTemp_FTF(1) + EnergyTemp_FTF(2)*2.0_dp
        endif

        if(.not.tHFRDMs) then
            close(unit_2RDM_read_aaaa)
            close(unit_2RDM_read_abab)
            close(unit_2RDM_read_abba)
        endif
        close(unit_rr)
!        close(unit_Phi_read_aaaa,status='delete')
!        close(unit_Phi_read_abab,status='delete')
!        close(unit_Phi_read_abba,status='delete')
        close(unit_Phi_read_aaaa)
        close(unit_Phi_read_abab)
        close(unit_Phi_read_abba)
        
        deallocate(rrints,CompOrbMat,NmoNmoMat,Gamma_vw,Gamma_wv,Phi_xy,OrbOrbMat)
        if(tDaltonFormat) deallocate(Buf,Indices)

    end subroutine CalcRRTerms

    pure integer function CountCabsInd(i,j,k,l)
        implicit none
        integer, intent(in) :: i,j,k,l

        CountCabsInd = 0
        if(i.gt.tnmo) CountCabsInd = CountCabsInd + 1
        if(j.gt.tnmo) CountCabsInd = CountCabsInd + 1
        if(k.gt.tnmo) CountCabsInd = CountCabsInd + 1
        if(l.gt.tnmo) CountCabsInd = CountCabsInd + 1

    end function CountCabsInd
    
    ! 1/8 (tau_rs)^xy (2RDM^rs)_xy - 1/8 (tau_rs)^yx (2RDM^rs)_xy
    !If tUseFTFInts, use integrals from FTFDUMPBIN
    !Otherwise, use integrals from FFDUMPBIN
    !This only affects non dalton type calculations
    subroutine CalcTauTerm(Energy_F12,PairEnergies,tUseFTFInts)
        use input_data, only : gamma_length
        implicit none
        real(dp), intent(inout) :: Energy_F12,PairEnergies(tnmo_sq)
        real(dp) :: EnergyTemp(3),DDOT,TermE,Temp,Z
        integer :: unit_2RDM_read_aaaa,unit_2RDM_read_abab,unit_2RDM_read_abba,ierr,unit_tau,ios,j
        real(dp) , allocatable :: tauints(:,:),Gamma_xy_aaaa(:),Gamma_xy_abab(:),Gamma_xy_abba(:)
        real(dp) , allocatable :: Buf(:),tauints_permute(:),Antisym_tau(:)
        integer(i2), allocatable :: Indices(:,:)
        integer(i8) :: maxlength,length
        integer :: r,s,i_Ind,j_Ind,k_Ind,l_Ind,PairInd,ibuf,rs_pair,xInd,yInd,x,y,xyPair,loop,i
        integer :: ijPair,jiPair
        character(*), parameter :: t_r = 'CalcTauTerm'
        logical :: tFTFexists
        logical , intent(in) :: tUseFTFInts
        
        inquire(file='FTFDUMPBIN',exist=tFTFexists)
        if((.not.tFTFexists).and.tUseFTFInts) return

        EnergyTemp(:) = 0.0_dp

        if(.not.tHFRDMs) then
            unit_2RDM_read_aaaa=get_free_unit()
            open(unit_2RDM_read_aaaa,file='TwoRDM_Dir_aaaa',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
            unit_2RDM_read_abab=get_free_unit()
            open(unit_2RDM_read_abab,file='TwoRDM_Dir_abab',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
            unit_2RDM_read_abba=get_free_unit()
            open(unit_2RDM_read_abba,file='TwoRDM_Dir_abba',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
        endif

        !Allocate array to store all xy and S (second index) for a given R
        allocate(tauints(tnmo_sq,tnmo),stat=ierr)
        allocate(tauints_permute(tnmo_sq),stat=ierr)
        allocate(Antisym_tau(tnmo_sq),stat=ierr)
        allocate(Gamma_xy_aaaa(tnmo_sq),stat=ierr)
        allocate(Gamma_xy_abab(tnmo_sq),stat=ierr)
        allocate(Gamma_xy_abba(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc err")

        unit_tau = get_free_unit()
        if(tDaltonFormat) then
            open(unit_tau,file='MO_FTF',status='old',form='unformatted',access='sequential',action='read')
        else
            !The FTF integrals can be simply calculated from the FF integrals
            if(tReadBin) then
                if(tUseFTFInts) then
                    open(unit_tau,file='FTFDUMPBIN',status='old',form='unformatted',action='read')
                else
                    open(unit_tau,file='FFDUMPBIN',status='old',form='unformatted',action='read')
                endif
            else
                if(tUseFTFInts) then
                    open(unit_tau,file='FTFDUMP',status='old',form='formatted',action='read')
                else
                    open(unit_tau,file='FFDUMP',status='old',form='formatted',action='read')
                endif
            endif
        endif

        rewind(unit_tau)
        if(tDaltonFormat) then
            read(unit_tau) maxlength
            allocate(Indices(4,MaxLength),Buf(MaxLength),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err")
        else
            length = 1
        endif

        do r=1,tnmo !Loop over all r variables
            !If doing MP2F12, only need to loop over occupied pairs here...
            if(tHFRDMs.and.(.not.IsOrbOcc(r))) cycle
            if(IsOrbFrz(r)) cycle
            rewind(unit_tau)
            tauints(:,:) = 0.0_dp
            if(tDaltonFormat) read(unit_tau) maxlength

            do while(.true.)
                if(tDaltonFormat) then
                    read(unit_tau,iostat=ios) length,Indices(1:4,1:length),Buf(1:length)
                else
                    if(tReadBin) then
                        read(unit_tau,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                    else
                        read(unit_tau,*,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                    endif
                endif
                if(ios.gt.0) call stop_all(t_r,"Error reading ints")
                if((length.le.0).or.(ios.lt.0)) exit

                do ibuf = 1,length
                    if(tDaltonFormat) then
                        i_Ind = int(Indices(1,ibuf),i4)
                        k_Ind = int(Indices(2,ibuf),i4)
                        j_Ind = int(Indices(3,ibuf),i4)
                        l_Ind = int(Indices(4,ibuf),i4)
                        Z = Buf(ibuf)
!                        if(i_Ind.eq.j_Ind.and.k_Ind.eq.l_Ind.and.j_Ind.eq.k_Ind.and.i_Ind.eq.1) then
!                            write(6,*) "FTF: ",i_Ind,j_Ind,k_Ind,l_Ind,Z
!                        endif
                    else
                        !If not tDaltonformat, assume that we are using STGs, and so there is a simple transformation from FF to FTF
!                        Z = Z*4.0_dp*(gamma_length**(2.0_dp))
                        Z = Z*(gamma_length**(2.0_dp))
                    endif

                    !Check if Orbital indices outside orbital space
                    if(max(i_Ind,max(j_Ind,max(k_Ind,l_Ind))).gt.tnmo) cycle
                    !Check if at least one of the orbitals is 'r'
                    if((i_Ind.ne.r).and.(j_Ind.ne.r).and.(k_Ind.ne.r).and.(l_Ind.ne.r)) cycle

                    !Fill first index with kl (l fast)
                    if(i_Ind.eq.r) then
                        !<ij|kl>
                        PairInd = FindSqPairInd(k_Ind,l_Ind,tnmo)
                        tauints(PairInd,j_Ind) = Z 
                        if(.not.tComplexOrbs) then
                            !<il|kj>
                            PairInd = FindSqPairInd(k_Ind,j_Ind,tnmo)
                            tauints(PairInd,l_Ind) = Z
                        endif
                    endif
                    if(j_Ind.eq.r) then
                        !<ji|lk>
                        PairInd = FindSqPairInd(l_Ind,k_Ind,tnmo)
                        tauints(PairInd,i_Ind) = Z
                        if(.not.tComplexOrbs) then
                            !<jk|li>
                            PairInd = FindSqPairInd(l_Ind,i_Ind,tnmo)
                            tauints(PairInd,k_Ind) = Z
                        endif
                    endif
                    if(k_Ind.eq.r) then
                        !<kl|ij>
                        PairInd = FindSqPairInd(i_Ind,j_Ind,tnmo)
                        tauints(PairInd,l_Ind) = Z
                        if(.not.tComplexOrbs) then
                            !<kj|il>
                            PairInd = FindSqPairInd(i_Ind,l_Ind,tnmo)
                            tauints(PairInd,j_Ind) = Z
                        endif
                    endif
                    if(l_Ind.eq.r) then
                        !<lk|ji>
                        PairInd = FindSqPairInd(j_Ind,i_Ind,tnmo)
                        tauints(PairInd,k_Ind) = Z
                        if(.not.tComplexOrbs) then
                            !<li|jk>
                            PairInd = FindSqPairInd(j_Ind,k_Ind,tnmo)
                            tauints(PairInd,i_Ind) = Z
                        endif
                    endif

                enddo   !End loop over ibuf
            enddo

            !tauints now holds all xy and s for a given r

            do s=1,tnmo
                !If doing MP2F12, only need to loop over occupied pairs here...
                if(tHFRDMs.and.(.not.IsOrbOcc(s))) cycle
                if(IsOrbFrz(s)) cycle

                !Read in all (Gamma^rs)_xy
                rs_pair = FindSqPairInd(r,s,tnmo)

                if(tHFRDMs) then
                    Gamma_xy_aaaa(:)=0.0_dp
                    Gamma_xy_abab(:)=0.0_dp
                    if(.not.tCalcCombo) Gamma_xy_abba(:)=0.0_dp
                    do x=1,tnocc
                        xInd = OccOrbs(x)
                        if(IsOrbFrz(xInd)) cycle
                        do y=1,tnocc
                            yInd = OccOrbs(y)
                            if(IsOrbFrz(yInd)) cycle
                            if((r.eq.xInd).and.(s.eq.yInd).and.(xInd.eq.yInd)) then
                                !All four indices the same
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Gamma_xy_abab(xyPair) = 1.0_dp
                                Gamma_xy_abba(xyPair) = -1.0_dp
                            elseif((r.eq.xInd).and.(s.eq.yInd)) then
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Gamma_xy_aaaa(xyPair) = 1.0_dp
                                Gamma_xy_abab(xyPair) = 1.0_dp
                            elseif((r.eq.yInd).and.(s.eq.xInd)) then
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Gamma_xy_aaaa(xyPair) = -1.0_dp
                                Gamma_xy_abba(xyPair) = -1.0_dp
                            endif
                        enddo
                    enddo
                else
                    read(unit_2RDM_read_aaaa,rec=rs_pair) (Gamma_xy_aaaa(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abab,rec=rs_pair) (Gamma_xy_abab(loop),loop=1,tnmo_sq)
                    if(.not.tCalcCombo) read(unit_2RDM_read_abba,rec=rs_pair) (Gamma_xy_abba(loop),loop=1,tnmo_sq)

                    !Zero frozen index elements
                    do i=1,tnfrz
                        do j=1,tnmo 
                            ijPair = FindSqPairInd(FrzOrbs(i),j,tnmo)
                            jiPair = FindSqPairInd(j,FrzOrbs(i),tnmo)
                            Gamma_xy_aaaa(ijPair) = 0.0_dp
                            Gamma_xy_aaaa(jiPair) = 0.0_dp
                            Gamma_xy_abab(ijPair) = 0.0_dp
                            Gamma_xy_abab(jiPair) = 0.0_dp
                            if(.not.tCalcCombo) then
                                Gamma_xy_abba(ijPair) = 0.0_dp
                                Gamma_xy_abba(jiPair) = 0.0_dp
                            endif
                        enddo
                    enddo

                endif
                
                !Extract relevant tau integrals, and antisymmetrise
                call TransposeIntegralArr(tauints(:,s),tauints_permute,tnmo)

                !aaaa / bbbb
                Antisym_tau(:) = (tauints(:,s) - tauints_permute(:))/16.0_dp
                Temp = 0.25_dp * DDOT(tnmo_sq,Antisym_tau,1,Gamma_xy_aaaa,1)
                if(tCalcCombo) Temp = Temp * 2.0_dp !To account for the loss of the permutation
                EnergyTemp(1) = EnergyTemp(1) + Temp 
                if(tDaltonFormat.or.tUseFTFInts) PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp

                if(.not.tCalcCombo) then
                    !Calculate second term, -1/8 (tau_rs)^yx (Gamma^rs)_xy
                    call TransposeIntegralArr_inplace(Antisym_tau,tnmo)
                    Temp = -0.25_dp * DDOT(tnmo_sq,Antisym_tau,1,Gamma_xy_aaaa,1)
                    EnergyTemp(1) = EnergyTemp(1) + Temp 
                    if(tDaltonFormat.or.tUseFTFInts) PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp
                endif

                !mixed spin term - first calculate all {xy} terms
                Antisym_tau(:) = ((5.0_dp/32.0_dp)*tauints(:,s)) + ((3.0_dp/32.0_dp)*tauints_permute(:))
                Temp = 0.25_dp * DDOT(tnmo_sq,Antisym_tau,1,Gamma_xy_abab,1)
                if(tCalcCombo) Temp = Temp * 2.0_dp
                EnergyTemp(2) = EnergyTemp(2) + Temp 
                if(tDaltonFormat.or.tUseFTFInts) PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp

                if(.not.tCalcCombo) then
                    !Don't calculate permutation or abba term if using tCalcCombo
                    Antisym_tau(:) = (tauints(:,s)*(-3.0_dp/32.0_dp)) - ((5.0_dp/32.0_dp)*tauints_permute(:))
                    Temp = 0.25_dp * DDOT(tnmo_sq,Antisym_tau,1,Gamma_xy_abba,1)
                    EnergyTemp(3) = EnergyTemp(3) + Temp 
                    if(tDaltonFormat.or.tUseFTFInts) PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp

                    !Then calculate all {yx} terms
                    Antisym_tau(:) = ((5.0_dp/32.0_dp)*tauints(:,s)) + ((3.0_dp/32.0_dp)*tauints_permute(:))
                    call TransposeIntegralArr_inplace(Antisym_tau,tnmo)
                    Temp = -0.25_dp * DDOT(tnmo_sq,Antisym_tau,1,Gamma_xy_abba,1)
                    EnergyTemp(2) = EnergyTemp(2) + Temp 
                    if(tDaltonFormat.or.tUseFTFInts) PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp

                    Antisym_tau(:) = (tauints(:,s)*(-3.0_dp/32.0_dp)) - ((5.0_dp/32.0_dp)*tauints_permute(:))
                    call TransposeIntegralArr_inplace(Antisym_tau,tnmo)
                    Temp = -0.25_dp * DDOT(tnmo_sq,Antisym_tau,1,Gamma_xy_abab,1)
                    EnergyTemp(3) = EnergyTemp(3) + Temp 
                    if(tDaltonFormat.or.tUseFTFInts) PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp
                endif

            enddo   !End do running over s

        enddo   !Enddo over r

        if(tCalcCombo) EnergyTemp(3) = EnergyTemp(2)    !Set the abba term to equal the abab term

        if(.not.tHFRDMs) then
            close(unit_2RDM_read_aaaa)
            close(unit_2RDM_read_abab)
            close(unit_2RDM_read_abba)
        endif
        close(unit_tau)
        deallocate(tauints,tauints_permute,Antisym_tau)
        deallocate(Gamma_xy_aaaa,Gamma_xy_abab,Gamma_xy_abba)
        if(tDaltonFormat) deallocate(Indices,Buf)
        
        if(tFTFexists.and.(.not.tUseFTFInts)) then
            !We are not calculating with the FTF integrals here, but are likely to later on
            write(6,"(A)") "Energy breakdown from tau term in B with FF integrals: "
            TermE = 0.0_dp
            do i=1,3
                TermE = TermE + EnergyTemp(i)
                write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
            enddo
            write(6,"(A,F25.12)") "Total energy from tau term in B with FF integrals: ",TermE
        elseif(tUseFTFInts) then
            write(6,"(A)") "Energy breakdown from tau term in B with FTF integrals: "
            TermE = 0.0_dp
            do i=1,3
                TermE = TermE + EnergyTemp(i)
                write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
            enddo
            write(6,"(A,F25.12)") "Total energy from tau term in B with FTF integrals: ",TermE
            write(6,"(A)") "(This is the one which is summed into final energy)"
            Energy_F12 = Energy_F12 + TermE
        else
            !Molecular calculation?
            write(6,"(A)") "Energy breakdown from tau term in B: "
            TermE = 0.0_dp
            do i=1,3
                TermE = TermE + EnergyTemp(i)
                write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
            enddo
            write(6,"(A,F25.12)") "Total energy from tau term in B: ",TermE
            Energy_F12 = Energy_F12 + TermE
        endif


    end subroutine CalcTauTerm


    ! -1/4 (fg_rs)^xy (2RDM^rs)_xy + 1/4 (fg_rs)^yx (2RDM^rs)_xy
    subroutine CalcFGTerm(Energy_F12,PairEnergies)
        implicit none
        real(dp), intent(inout) :: Energy_F12,PairEnergies(tnmo_sq)
        real(dp) :: EnergyTemp(3),DDOT,TermE,Temp,Z
        integer :: unit_2RDM_read_aaaa,unit_2RDM_read_abab,unit_2RDM_read_abba,ierr,unit_fg,ios
        real(dp) , allocatable :: fgints(:,:),Gamma_xy_aaaa(:),Gamma_xy_abab(:),Gamma_xy_abba(:)
        real(dp) , allocatable :: Buf(:),fgints_permute(:),Antisym_fg(:)
        integer(i2), allocatable :: Indices(:,:)
        integer(i8) :: maxlength,length
        integer :: r,s,i_Ind,j_Ind,k_Ind,l_Ind,PairInd,ibuf,rs_pair,xInd,yInd,x,y,xyPair,loop,i
        character(*), parameter :: t_r = 'CalcFGTerm'

        EnergyTemp(:) = 0.0_dp

        if(.not.tHFRDMs) then
            unit_2RDM_read_aaaa=get_free_unit()
            open(unit_2RDM_read_aaaa,file='TwoRDM_Dir_aaaa',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
            unit_2RDM_read_abab=get_free_unit()
            open(unit_2RDM_read_abab,file='TwoRDM_Dir_abab',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
            unit_2RDM_read_abba=get_free_unit()
            open(unit_2RDM_read_abba,file='TwoRDM_Dir_abba',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
        endif

        !Allocate array to store all xy and S (second index) for a given R
        allocate(fgints(tnmo_sq,tnmo),stat=ierr)
        allocate(fgints_permute(tnmo_sq),stat=ierr)
        allocate(Antisym_fg(tnmo_sq),stat=ierr)
        allocate(Gamma_xy_aaaa(tnmo_sq),stat=ierr)
        allocate(Gamma_xy_abab(tnmo_sq),stat=ierr)
        allocate(Gamma_xy_abba(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc err")

        unit_fg = get_free_unit()
        if(tDaltonFormat) then
            open(unit_fg,file='MO_FG',status='old',form='unformatted',access='sequential',action='read')
        else
            if(tReadBin) then
                open(unit_fg,file='FGDUMPBIN',status='old',form='unformatted',action='read')
            else
                open(unit_fg,file='FGDUMP',status='old',form='formatted',action='read')
            endif
        endif
        rewind(unit_fg)
        if(tDaltonFormat) then
            read(unit_fg) maxlength
            allocate(Indices(4,MaxLength),Buf(MaxLength),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err")
        else
            length = 1
        endif

        do r=1,tnmo !Loop over all r variables
            !If doing MP2F12, only need to loop over occupied orbital pairs here.
            if(tHFRDMs.and.(.not.IsOrbOcc(r))) cycle
            if(IsOrbFrz(r)) cycle   !r isn't geminal, but take transpose of both matrices (ok since symmetric)
                                    !Then we are restricting r and s to be active, rather than x and y
            rewind(unit_fg)
            fgints(:,:) = 0.0_dp
            if(tDaltonFormat) read(unit_fg) maxlength

            do while(.true.)
                if(tDaltonFormat) then
                    read(unit_fg,iostat=ios) length,Indices(1:4,1:length),Buf(1:length)
                else
                    if(tReadBin) then
                        read(unit_fg,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                    else
                        read(unit_fg,*,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                    endif
                endif
                if(ios.gt.0) call stop_all(t_r,"Error reading ints")
                if((length.le.0).or.(ios.lt.0)) exit

                do ibuf = 1,length
                    if(tDaltonFormat) then
                        i_Ind = int(Indices(1,ibuf),i4)
                        k_Ind = int(Indices(2,ibuf),i4)
                        j_Ind = int(Indices(3,ibuf),i4)
                        l_Ind = int(Indices(4,ibuf),i4)
                        Z = Buf(ibuf)
!                        if(i_Ind.eq.j_Ind.and.k_Ind.eq.l_Ind.and.i_Ind.eq.k_Ind.and.i_Ind.eq.1) then
!                            write(6,*) i_Ind,j_Ind,k_Ind,l_Ind,Z
!                        endif
                    endif

                    !Check if Orbital indices outside orbital space
                    if(max(i_Ind,max(j_Ind,max(k_Ind,l_Ind))).gt.tnmo) cycle
                    !Check if at least one of the orbitals is 'r'
                    if((i_Ind.ne.r).and.(j_Ind.ne.r).and.(k_Ind.ne.r).and.(l_Ind.ne.r)) cycle

                    !Fill first index with kl (l fast)
                    if(i_Ind.eq.r) then
                        !<ij|kl>
                        PairInd = FindSqPairInd(k_Ind,l_Ind,tnmo)
                        fgints(PairInd,j_Ind) = Z 
                        if(.not.tComplexOrbs) then
                            !<il|kj>
                            PairInd = FindSqPairInd(k_Ind,j_Ind,tnmo)
                            fgints(PairInd,l_Ind) = Z
                        endif
                    endif
                    if(j_Ind.eq.r) then
                        !<ji|lk>
                        PairInd = FindSqPairInd(l_Ind,k_Ind,tnmo)
                        fgints(PairInd,i_Ind) = Z
                        if(.not.tComplexOrbs) then
                            !<jk|li>
                            PairInd = FindSqPairInd(l_Ind,i_Ind,tnmo)
                            fgints(PairInd,k_Ind) = Z
                        endif
                    endif
                    if(k_Ind.eq.r) then
                        !<kl|ij>
                        PairInd = FindSqPairInd(i_Ind,j_Ind,tnmo)
                        fgints(PairInd,l_Ind) = Z
                        if(.not.tComplexOrbs) then
                            !<kj|il>
                            PairInd = FindSqPairInd(i_Ind,l_Ind,tnmo)
                            fgints(PairInd,j_Ind) = Z
                        endif
                    endif
                    if(l_Ind.eq.r) then
                        !<lk|ji>
                        PairInd = FindSqPairInd(j_Ind,i_Ind,tnmo)
                        fgints(PairInd,k_Ind) = Z
                        if(.not.tComplexOrbs) then
                            !<li|jk>
                            PairInd = FindSqPairInd(j_Ind,k_Ind,tnmo)
                            fgints(PairInd,i_Ind) = Z
                        endif
                    endif

                enddo   !End loop over ibuf
            enddo

            !fgints now holds all xy and s for a given r

            do s=1,tnmo
                !If doing MP2F12, only need to loop over occupied orbital pairs here...
                if(tHFRDMs.and.(.not.IsOrbOcc(s))) cycle
                if(IsOrbFrz(s)) cycle

                !Read in all (Gamma^rs)_xy
                rs_pair = FindSqPairInd(r,s,tnmo)

                if(tHFRDMs) then
                    Gamma_xy_aaaa(:)=0.0_dp
                    Gamma_xy_abab(:)=0.0_dp
                    if(.not.tCalcCombo) Gamma_xy_abba(:)=0.0_dp
                    do x=1,tnocc
                        xInd = OccOrbs(x)
                        do y=1,tnocc
                            yInd = OccOrbs(y)
                            if((r.eq.xInd).and.(s.eq.yInd).and.(xInd.eq.yInd)) then
                                !All four indices the same
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Gamma_xy_abab(xyPair) = 1.0_dp
                                Gamma_xy_abba(xyPair) = -1.0_dp
                            elseif((r.eq.xInd).and.(s.eq.yInd)) then
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Gamma_xy_aaaa(xyPair) = 1.0_dp
                                Gamma_xy_abab(xyPair) = 1.0_dp
                            elseif((r.eq.yInd).and.(s.eq.xInd)) then
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Gamma_xy_aaaa(xyPair) = -1.0_dp
                                Gamma_xy_abba(xyPair) = -1.0_dp
                            endif
                        enddo
                    enddo
                else
                    read(unit_2RDM_read_aaaa,rec=rs_pair) (Gamma_xy_aaaa(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abab,rec=rs_pair) (Gamma_xy_abab(loop),loop=1,tnmo_sq)
                    if(.not.tCalcCombo) read(unit_2RDM_read_abba,rec=rs_pair) (Gamma_xy_abba(loop),loop=1,tnmo_sq)
                endif
                
                !Extract relevant fg integrals, and antisymmetrise
                call TransposeIntegralArr(fgints(:,s),fgints_permute,tnmo)

                !aaaa / bbbb
                Antisym_fg(:) = (fgints(:,s) - fgints_permute(:))/4.0_dp
                Temp = 0.5_dp * DDOT(tnmo_sq,Antisym_fg,1,Gamma_xy_aaaa,1)
                if(tCalcCombo) Temp = Temp * 2.0_dp !Account for permutation
                EnergyTemp(1) = EnergyTemp(1) + Temp 
                PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp

                if(.not.tCalcCombo) then
                    !Ignore permutation if using tCalcCombo
                    call TransposeIntegralArr_inplace(Antisym_fg,tnmo)
                    Temp = -0.5_dp * DDOT(tnmo_sq,Antisym_fg,1,Gamma_xy_aaaa,1)
                    EnergyTemp(1) = EnergyTemp(1) + Temp 
                    PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp
                endif

                !mixed spin term - first calculate all {xy} terms
                Antisym_fg(:) = ((3.0_dp/8.0_dp)*fgints(:,s)) + (fgints_permute(:)/8.0_dp)
                Temp = 0.5_dp * DDOT(tnmo_sq,Antisym_fg,1,Gamma_xy_abab,1)
                if(tCalcCombo) Temp = Temp * 2.0_dp !Account for permutation
                EnergyTemp(2) = EnergyTemp(2) + Temp 
                PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp

                if(.not.tCalcCombo) then
                    Antisym_fg(:) = (fgints(:,s)/(-8.0_dp)) - ((3.0_dp/8.0_dp)*fgints_permute(:))
                    Temp = 0.5_dp * DDOT(tnmo_sq,Antisym_fg,1,Gamma_xy_abba,1)
                    EnergyTemp(3) = EnergyTemp(3) + Temp 
                    PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp

                    !Then calculate all {yx} terms
                    Antisym_fg(:) = ((3.0_dp/8.0_dp)*fgints(:,s)) + (fgints_permute(:)/8.0_dp)
                    call TransposeIntegralArr_inplace(Antisym_fg,tnmo)
                    Temp = -0.5_dp * DDOT(tnmo_sq,Antisym_fg,1,Gamma_xy_abba,1)
                    EnergyTemp(2) = EnergyTemp(2) + Temp 
                    PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp

                    Antisym_fg(:) = (fgints(:,s)/(-8.0_dp)) - ((3.0_dp/8.0_dp)*fgints_permute(:))
                    call TransposeIntegralArr_inplace(Antisym_fg,tnmo)
                    Temp = -0.5_dp * DDOT(tnmo_sq,Antisym_fg,1,Gamma_xy_abab,1)
                    EnergyTemp(3) = EnergyTemp(3) + Temp 
                    PairEnergies(rs_pair) = PairEnergies(rs_pair) + Temp
                endif

            enddo   !End do running over s

        enddo   !Enddo over r

        if(tCalcCombo) EnergyTemp(3) = EnergyTemp(2)    !abba gives same energy as abab

        if(.not.tHFRDMs) then
            close(unit_2RDM_read_aaaa)
            close(unit_2RDM_read_abab)
            close(unit_2RDM_read_abba)
        endif
        close(unit_fg)
        deallocate(fgints,fgints_permute,Antisym_fg)
        deallocate(Gamma_xy_aaaa,Gamma_xy_abab,Gamma_xy_abba)
        if(tDaltonFormat) deallocate(Indices,Buf)
        
        write(6,"(A)") "Energy breakdown from fg term in V: "
        TermE = 0.0_dp
        do i=1,3
            TermE = TermE + EnergyTemp(i)
            write(6,"(A,I3,A,F25.12)") "Spin: ",i," Energy contribution: ",EnergyTemp(i)
        enddo
        write(6,"(A,F25.12)") "Total energy from fg term in V: ",TermE
        Energy_F12 = Energy_F12 + TermE

    end subroutine CalcFGTerm

    subroutine write_matrix(Filename,Arr,Dimen)
        implicit none
        character(*), intent(in) :: Filename
        integer, intent(in) :: Dimen
        real(dp), intent(in) :: Arr(Dimen,Dimen)
        integer :: iunit,i,j

        iunit=get_free_unit()
        open(iunit,file=trim(Filename),status='unknown')
        do i=1,Dimen
            do j=1,Dimen
                write(iunit,"(F13.7)",advance='no') Arr(j,i)
            enddo
            write(iunit,*)
        enddo

        close(iunit)
    end subroutine write_matrix

    
!For a given {v,w} pair, construct the appropriate antisymmetrised R_Tld array
!End up with (r_tld_vw)^pq, with **p** as the fast index (actually, fast index is dependant on space required)
!R_TLD arrays have two-fold permutational symmetry - (r_tld)_ij^kl = (r_tld)_ji^lk
!ispin = 1:    aaaa
!ispin = 2:    abab
!ispin = 3:    baab
!BEWARE OF ACCIDENTAL TRANSPOSES WHEN USING THIS ROUTINE! Care with index ordering needed. 
!*****   RTld is stored ' p fast ' for a given pq.
    subroutine GetAntisymR_TLD(unit_raaaa,unit_rabab,unit_rabba,v,w,Arr,ispin,MinOrbp,MaxOrbp,MinOrbq,MaxOrbq)
        implicit none
        integer , intent(in) :: v,w,ispin,MaxOrbp,MaxOrbq,MinOrbp,MinOrbq
        real(dp) , intent(out) :: Arr((MaxOrbp-MinOrbp+1)*(MaxOrbq-MinOrbq+1))
        integer , intent(in) :: unit_raaaa, unit_rabab,unit_rabba
        real(dp), allocatable :: TempInts(:)
        integer :: vwPair,loop,RTld_unit,ierr,p,q,counter,p_prime,a_prime
        character(*), parameter :: t_r='GetAntisymR_TLD'

        Arr(:)=0.0_dp
        vwPair = FindSqPairInd(v,w,tnmo)    !R_TLD only has orbital pairs on disk

        allocate(TempInts(tntmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error allocating temp array")

        if(ispin.eq.1) then
            RTld_unit = unit_raaaa
        elseif(ispin.eq.2) then
            RTld_unit = unit_rabab
        elseif(ispin.eq.3) then
            RTld_unit = unit_rabba
        else
            call stop_all(t_r,"Invalid spin")
        endif

        read(RTld_unit,rec=vwPair) (TempInts(loop),loop=1,tntmo_sq)
        if((MaxOrbp.eq.MaxOrbq).and.(MinOrbp.eq.MinOrbq).and.(MinOrbp.eq.1)) then
            !Only interested in pq both over orbital space
            do q=MinOrbq,MaxOrbq
                do p=MinOrbp,MaxOrbp
                    !Remember, p fast.
                    Arr(FindSqPairInd(q,p,MaxOrbp)) = TempInts(FindSqPairInd(q,p,tntmo))
                enddo
            enddo
        elseif((MinOrbp.eq.1).and.(MaxOrbp.eq.tnmo).and.(MinOrbq.eq.(tnmo+1)).and.(MaxOrbq.eq.tntmo)) then
            !We are interested in p over orbital space, and q over CABS space.
            do q=MinOrbq,MaxOrbq
                do p=MinOrbp,MaxOrbp
                    !Actually, a_prime (cabs) is fast here....
                    Arr(FindCabsOrbPairInd(q,p)) = TempInts(FindSqPairInd(q,p,tntmo))
                enddo
            enddo
        elseif((MinOrbq.eq.1).and.(MaxOrbq.eq.tnmo).and.(MinOrbp.eq.(tnmo+1)).and.(MaxOrbp.eq.tntmo)) then
            !We are interested in q over orbital space, and p over CABS space.
            do q=MinOrbq,MaxOrbq
                do p=MinOrbp,MaxOrbp    !p is CABS here
                    !CABS space always fast here
                    !TempInts gets RTld_a'q since first index is fast in RTld!!
                    Arr(FindCabsOrbPairInd(p,q)) = TempInts(FindSqPairInd(q,p,tntmo))
                enddo
            enddo
        elseif((MinOrbp.eq.1).and.(MaxOrbp.eq.tntmo).and.(MinOrbq.eq.(tnmo+1)).and.(MaxOrbq.eq.(tntmo))) then
            !First index over complete space, second index over aux space.
            !However, we want the second index to be fast.
            counter=1
            do p_prime=1,tntmo
                do a_prime=(tnmo+1),tntmo
                    Arr(counter) = TempInts(FindSqPairInd(p_prime,a_prime,tntmo))
                    counter=counter+1
                enddo
            enddo
        elseif((min(MaxOrbp,MaxOrbq).eq.tntmo).and.(max(MinOrbp,MinOrbq).eq.1)) then
            !Both indices over complete space.
            !We want second index to be fast.
            !This is already done!
            Arr(:) = TempInts(:)
        else
            call stop_all(t_r,"Other RTld shapes not coded up yet")
        endif

        deallocate(TempInts)

    end subroutine GetAntisymR_TLD
                    
!Arr runs over integrals from all tntmo pairs.
!However, we only want antisymmetrised pairs over {a_prime,p} with a_prime fast!
!Used for V2
    subroutine Antisym0Cusp_ooao_from_oocc(Arr,AntisymArr,ispin)
        implicit none
        real(dp) , intent(in) :: Arr(tntmo_sq)
        real(dp) , intent(out) :: AntisymArr(taux_nmo)
        integer , intent(in) :: ispin
        integer :: loop,t,a_prime
    
        !Antisymmetrise over t a' (no geminal indices)
        AntisymArr(:) = 0.0_dp
        loop = 1
        if(ispin.eq.1) then
            do t=1,tnmo
                do a_prime=tnmo+1,tntmo
                    !a_prime is fast in r_ta
                    AntisymArr(loop) = Arr(FindSqPairInd(t,a_prime,tntmo)) &
                        - Arr(FindSqPairInd(a_prime,t,tntmo))
                    loop = loop + 1    !Should give the same index as FindCABSOrbPairInd(a_prime,t)
                enddo
            enddo
        elseif(ispin.eq.2) then
            do t=1,tnmo
                do a_prime=tnmo+1,tntmo
                    !a_prime is fast in r_ta
                    AntisymArr(loop) = Arr(FindSqPairInd(t,a_prime,tntmo))
                    loop = loop + 1    !Should give the same index as FindCABSOrbPairInd(a_prime,t)
                enddo
            enddo
        elseif(ispin.eq.3) then
            do t=1,tnmo
                do a_prime=tnmo+1,tntmo
                    AntisymArr(loop) = -Arr(FindSqPairInd(a_prime,t,tntmo))
                    !a_prime is fast in r_ta
                    loop = loop + 1    !Should give the same index as FindCABSOrbPairInd(a_prime,t)
                enddo
            enddo
        endif

    end subroutine Antisym0Cusp_ooao_from_oocc

!Arr runs over integrals from all tntmo pairs.
!However, we only want antisymmetrised pairs over {p',p} with p_prime fast!
    subroutine Antisym1Cusp_ooca_from_oocc(Arr,AntisymArr,ispin)
        implicit none
        real(dp) , intent(in) :: Arr(tntmo_sq)
        real(dp) , intent(out) :: AntisymArr(tntmo_nxmo)
        integer , intent(in) :: ispin
        integer :: p_prime,a_prime,ijPair,jiPair,loop

        AntisymArr(:)=0.0_dp
        loop=1
        if(ispin.eq.1) then
            !alpha alpha alpha alpha spin: <aa||aa> = 1/4 [ <ij|kl> - <ij|lk> ]
            do a_prime=tnmo+1,tntmo
                do p_prime=1,tntmo
                    ijPair = FindSqPairInd(p_prime,a_prime,tntmo)
                    jiPair = FindSqPairInd(a_prime,p_prime,tntmo)

                    AntisymArr(loop)=(Arr(ijPair)-Arr(jiPair))/4.0_dp
                    loop=loop+1
                enddo
            enddo
        elseif(ispin.eq.2) then
            !a b a b: <ab||ab> = 3/8 <ij|kl> + 1/8 <ij|lk>
            do a_prime=tnmo+1,tntmo
                do p_prime=1,tntmo
                    ijPair = FindSqPairInd(p_prime,a_prime,tntmo)
                    jiPair = FindSqPairInd(a_prime,p_prime,tntmo)

                    AntisymArr(loop)=((3.0_dp/8.0_dp)*Arr(ijPair)) + (Arr(jiPair)/8.0_dp)
                    loop=loop+1
                enddo
            enddo
        elseif(ispin.eq.3) then
            !a b b a: <ab|ba> = -1/8<ij|kl> - 3/8<ij|lk>
            do a_prime=tnmo+1,tntmo
                do p_prime=1,tntmo
                    ijPair = FindSqPairInd(p_prime,a_prime,tntmo)
                    jiPair = FindSqPairInd(a_prime,p_prime,tntmo)

                    AntisymArr(loop)=(Arr(ijPair)/(-8.0_dp)) - ((3.0_dp/8.0_dp)*Arr(jiPair))
                    loop=loop+1
                enddo
            enddo
        endif
    end subroutine Antisym1Cusp_ooca_from_oocc
    
!Arr runs over integrals from all tntmo pairs.
!However, we only want antisymmetrised pairs over {a_prime,p} with a_prime fast!
!tPermute_ao means that we are actually interested in the {p,a_prime} pair (still with a_prime fast)
    subroutine Antisym1Cusp_ooao_from_oocc(Arr,AntisymArr,ispin,tPermute_ao)
        implicit none
        real(dp) , intent(in) :: Arr(tntmo_sq)
        real(dp) , intent(out) :: AntisymArr(taux_nmo)
        integer , intent(in) :: ispin
        logical , intent(in) :: tPermute_ao
        integer :: a_prime,p,ijPair,jiPair,AntisymPair

        AntisymArr(:)=0.0_dp
        if(ispin.eq.1) then
            !alpha alpha alpha alpha spin: <aa||aa> = 1/4 [ <ij|kl> - <ij|lk> ]
            do p=1,tnmo
                do a_prime=tnmo+1,tntmo
                    ijPair = FindSqPairInd(a_prime,p,tntmo)
                    jiPair = FindSqPairInd(p,a_prime,tntmo)

                    AntisymPair = FindCabsOrbPairInd(a_prime,p)
                    if(tPermute_ao) then
                        AntisymArr(AntisymPair)=(Arr(jiPair)-Arr(ijPair))/4.0_dp
                    else
                        AntisymArr(AntisymPair)=(Arr(ijPair)-Arr(jiPair))/4.0_dp
                    endif
!                    write(6,*) i,j,AntisymArr(ijPair)
                enddo
            enddo
        elseif(ispin.eq.2) then
            !a b a b: <ab||ab> = 3/8 <ij|kl> + 1/8 <ij|lk>
            do p=1,tnmo
                do a_prime=tnmo+1,tntmo
                    ijPair = FindSqPairInd(a_prime,p,tntmo)
                    jiPair = FindSqPairInd(p,a_prime,tntmo)

                    AntisymPair = FindCabsOrbPairInd(a_prime,p)
                    if(tPermute_ao) then
                        AntisymArr(AntisymPair)=((3.0_dp/8.0_dp)*Arr(jiPair)) + (Arr(ijPair)/8.0_dp)
                    else
                        AntisymArr(AntisymPair)=((3.0_dp/8.0_dp)*Arr(ijPair)) + (Arr(jiPair)/8.0_dp)
                    endif
                enddo
            enddo
        elseif(ispin.eq.3) then
            !a b b a: <ab|ba> = -1/8<ij|kl> - 3/8<ij|lk>
            do p=1,tnmo
                do a_prime=tnmo+1,tntmo
                    ijPair = FindSqPairInd(a_prime,p,tntmo)
                    jiPair = FindSqPairInd(p,a_prime,tntmo)

                    AntisymPair = FindCabsOrbPairInd(a_prime,p)
                    if(tPermute_ao) then
                        AntisymArr(AntisymPair)=(Arr(jiPair)/(-8.0_dp)) - ((3.0_dp/8.0_dp)*Arr(ijPair))
                    else
                        AntisymArr(AntisymPair)=(Arr(ijPair)/(-8.0_dp)) - ((3.0_dp/8.0_dp)*Arr(jiPair))
                    endif
                enddo
            enddo
        endif
    end subroutine Antisym1Cusp_ooao_from_oocc

    subroutine Antisym1Cusp(Arr,AntisymArr,nOrbs,ispin)
        implicit none
        real(dp) , intent(in) :: Arr(nOrbs*nOrbs)
        real(dp) , intent(out) :: AntisymArr(nOrbs*nOrbs)
        integer , intent(in) :: nOrbs,ispin
        integer :: i,j,ijPair,jiPair

        AntisymArr(:)=0.0_dp
        if(ispin.eq.1) then
            !alpha alpha alpha alpha spin: <aa||aa> = 1/4 [ <ij|kl> - <ij|lk> ]
            do i=1,nOrbs
                do j=1,nOrbs
                    ijPair = FindSqPairInd(i,j,nOrbs)
                    jiPair = FindSqPairInd(j,i,nOrbs)
                    AntisymArr(ijPair)=(Arr(ijPair)-Arr(jiPair))/4.0_dp
!                    write(6,*) i,j,AntisymArr(ijPair)
                enddo
            enddo
        elseif(ispin.eq.2) then
            !a b a b: <ab||ab> = 3/8 <ij|kl> + 1/8 <ij|lk>
            do i=1,nOrbs
                do j=1,nOrbs
                    ijPair = FindSqPairInd(i,j,nOrbs)
                    jiPair = FindSqPairInd(j,i,nOrbs)
                    AntisymArr(ijPair)=((3.0_dp/8.0_dp)*Arr(ijPair)) + (Arr(jiPair)/8.0_dp)
                enddo
            enddo
        elseif(ispin.eq.3) then
            !a b b a: <ab|ba> = -1/8<ij|kl> - 3/8<ij|lk>
            do i=1,nOrbs
                do j=1,nOrbs
                    ijPair = FindSqPairInd(i,j,nOrbs)
                    jiPair = FindSqPairInd(j,i,nOrbs)
                    AntisymArr(ijPair)=(Arr(ijPair)/(-8.0_dp)) - ((3.0_dp/8.0_dp)*Arr(jiPair))
                enddo
            enddo
        endif
    end subroutine Antisym1Cusp
    
    subroutine Antisym0Cusp(Arr,AntisymArr,nOrbs,ispin)
        implicit none
        real(dp) , intent(in) :: Arr(nOrbs*nOrbs)
        real(dp) , intent(out) :: AntisymArr(nOrbs*nOrbs)
        integer , intent(in) :: nOrbs,ispin

        AntisymArr(:)=0.0_dp
        if(ispin.eq.1) then
            !alpha alpha alpha alpha spin: <aa||aa> = <ij|kl> - <ij|lk> 
            call TransposeIntegralArr(Arr,AntisymArr,nOrbs) !This will get <ij|lk>
            AntisymArr(:) = - AntisymArr(:) + Arr(:)
!            do i=1,nOrbs
!                do j=1,nOrbs
!                    ijPair = FindSqPairInd(i,j,nOrbs)
!                    jiPair = FindSqPairInd(j,i,nOrbs)
!                    AntisymArr(ijPair)=(Arr(ijPair)-Arr(jiPair))
!!                    write(6,*) i,j,AntisymArr(ijPair)
!                enddo
!            enddo
        elseif(ispin.eq.2) then
            !a b a b: <ab||ab> = same as spatial! 
            AntisymArr(:) = Arr(:)
        elseif(ispin.eq.3) then
            !a b b a: <ab|ba> = - spatial integrals
            call TransposeIntegralArr(Arr,AntisymArr,nOrbs) !This will get <ij|lk>
            AntisymArr(:) = -AntisymArr(:)
        endif
    end subroutine Antisym0Cusp

    subroutine CheckRDMProperties()
        use input_data, only : tIntegrateTwoRDM,tUEG
        implicit none
        integer :: ierr,lWork,info,k,x,y,xInd,yInd,xyPair,in_pair
        real(dp), allocatable :: Two_Gamma_aaaa(:), Two_Gamma_abab(:), Two_Gamma_abba(:)
        real(dp), allocatable :: Two_Gamma_ji_aaaa(:), Two_Gamma_ji_abab(:), Two_Gamma_ji_abba(:)
        real(dp), allocatable :: OneRDMCopy(:,:),W(:),Work(:),g_in(:)
        real(dp) :: OneRDMTrace,TwoRDMTrace,DDOT,SpinFac,SpinA,SpinB,LinearRelDiff,AntisymDiff
        character(*), parameter :: t_r="CheckRDMProperties"
        integer :: unit_g_stored,i,j,unit_2RDM_read_aaaa,unit_2RDM_read_abab,unit_2RDM_read_abba,ijInd,loop,jiInd
        logical :: tNewErr
        
        if(.not.tHFRDMs) then
            write(6,"(A)") "Checking RDMs for permutational symmetry, hermiticity, N-representability and trace conditions..."
        endif
        
        unit_g_stored=get_free_unit()
        open(unit_g_stored,file='G_Dir_ccoo',status='old',form='unformatted',access='direct',action='read',recl=reclen_nmo_sq)
        allocate(g_in(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"allocation err")

!        write(6,"(A)") "Calculating variational energy from density matrices..."

        !Now calculate variational energy from RDM
        VarE = 0.0_dp
        OneRDMTrace = 0.0_dp    !Check that trace of 1RDM = N

        !One electron part  (Only run over alpha/alpha, which we store, and *2)
        do i=1,tnmo
            do j=1,tnmo
                VarE = VarE + tmat(i,j)*OneRDM(j,i)*2
                if(abs(OneRDM(j,i)-OneRDM(i,j)).gt.1.D-8) then
                    call stop_all(t_r,"Hermiticity condition on 1RDM not satisfied")
                endif
                if(i.eq.j) then
                    OneRDMTrace = OneRDMTrace + OneRDM(i,j)*2
                else
                    if(tUEG.and.(abs(OneRDM(i,j).gt.1.D-8))) then
                        call stop_all(t_r,"OneRDM is not diagonal for the UEG system!")
                    endif
                endif
            enddo
        enddo
        if(.not.tHFRDMs) write(6,"(A)") "Hermiticity condition on 1RDM satisfied"
        write(6,"(A,G25.15)") "One electron energy from 1RDM: ",VarE

        if(abs(OneRDMTrace-real(NEl,dp)).gt.1.D-6) then
            write(6,*) OneRDMTrace,NEl,abs(OneRDMTrace-NEl)
            call stop_all(t_r,"Trace condition on 1RDM not satisfied")
        else
            if(.not.tHFRDMs) write(6,"(A)") "Trace condition on 1RDM satisfied"
        endif

        !Diagonalising OneRDM and I - OneRDM to check complete N-representability conditions on 1RDM
        !Create another copy of OneRDM, since it will be destroyed by diagonalisation procedure
        allocate(OneRDMCopy(tnmo,tnmo),stat=ierr)
        allocate(W(tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Allocation Error")
        OneRDMCopy = OneRDM

        lWork = -1 
        allocate(Work(1))
        !Worksize queiry
        call DSYEV('N','U',tnmo,OneRDMCopy,tnmo,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,"Workspace quiery failed")

        lWork = Work(1)
        deallocate(Work)
        allocate(Work(lWork),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Allocation Error")
        OneRDMCopy = OneRDM
        call DSYEV('N','U',tnmo,OneRDMCopy,tnmo,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,"1RDM diagonalisation failed")
            
!        write(6,*) W(:)
        do i=1,tnmo
            if(W(i).lt.-1.D-9) then
                write(6,*) i,W(i)
                call warning(t_r,"N-representability constraint on 1RDM not satisfied - negative eigenvalue found")
            endif
        enddo
!        write(6,"(A)") "No negative eigenvalues found for 1RDM"

        !Now check I-OneRDM
        OneRDMCopy = OneRDM
        do i=1,tnmo
            OneRDMCopy(i,i) = 1.0_dp - OneRDMCopy(i,i)
        enddo
        call DSYEV('N','U',tnmo,OneRDMCopy,tnmo,W,Work,lWork,info)
        if(info.ne.0) call stop_all(t_r,"1RDM diagonalisation failed")
        do i=1,tnmo
            if(W(i).lt.-1.D-9) then
                write(6,*) i,W(i)
                call warning(t_r,"N-representability constraint on I-1RDM not satisfied - negative eigenvalue found")
            endif
        enddo
!        write(6,"(A)") "No negative eigenvalues found for I - 1RDM" 
        if(.not.tHFRDMs) write(6,"(A)") "Complete N-representability conditions satisfied for 1RDM"
        deallocate(W,Work)

        !Allocate memory for each ij slice of 2RDM
        allocate(Two_Gamma_aaaa(tnmo_sq),stat=ierr) !This stores all ij in (2R)_kl^ij    
        allocate(Two_Gamma_abab(tnmo_sq),stat=ierr) !This stores all ij in (2R)_kl^ij    
        allocate(Two_Gamma_abba(tnmo_sq),stat=ierr) !This stores all ij in (2R)_kl^ij    
        allocate(Two_Gamma_ji_aaaa(tnmo_sq),stat=ierr) !This stores all ji in (2R)_kl^ji    
        allocate(Two_Gamma_ji_abab(tnmo_sq),stat=ierr) !This stores all ji in (2R)_kl^ji    
        allocate(Two_Gamma_ji_abba(tnmo_sq),stat=ierr) !This stores all ji in (2R)_kl^ji    
        if(ierr.ne.0) call stop_all(t_r,"allocation err")
        TwoRDMTrace = 0.0_dp

        !This now holds information for the linear relation:
        !sum_j  2R_ij^kj = (N-1) 1R_i^k
        OneRDMCopy = 0.0_dp

        if(.not.tHFRDMs) then
            unit_2RDM_read_aaaa=get_free_unit()
            open(unit_2RDM_read_aaaa,file='TwoRDM_Dir_aaaa',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
            unit_2RDM_read_abab=get_free_unit()
            open(unit_2RDM_read_abab,file='TwoRDM_Dir_abab',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
            unit_2RDM_read_abba=get_free_unit()
            open(unit_2RDM_read_abba,file='TwoRDM_Dir_abba',status='old',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
        endif

        !Two electron part
        !Trace condition on 2-electron part: sum_ij 2RDM_ij^ij = N(N-1)
        !Check antisymmetry condition of 2RDMs: 2R_ij^kl = -2R_ij^lk = -2R_ji^kl = 2R_ji^lk = 2R_kl^ij (last one hermiticity - not checked)
        !Can check linear equality which has to be satisfied: sum_k 2R_ij^jk = (N-1)1RDM(i,j) for all i,j
        !Check hermiticity of 1&2RDM - this is difficult - easier to do this within NECI (if not explicitly done)
        !More complicated N-representability conditions possible to check - involve diagonalisation of N^4 matrix
        SpinFac = 0.0_dp
        AntisymDiff = 0.0_dp
        do i=1,tnmo
            do j=1,tnmo
                !Read in / calculate all {k,l} for the 2RDM
                !This also wants to refer to spatial orbitals?
                if(tHFRDMs) then
                    Two_Gamma_aaaa(:)=0.0_dp
                    Two_Gamma_abab(:)=0.0_dp
                    Two_Gamma_abba(:)=0.0_dp
                    do x=1,tnocc
                        xInd = OccOrbs(x)
                        do y=1,tnocc
                            yInd = OccOrbs(y)
                            if((i.eq.xInd).and.(j.eq.yInd).and.(xInd.eq.yInd)) then
                                !All four indices the same = 0
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Two_Gamma_abab(xyPair) = 1.0_dp
                                Two_Gamma_abba(xyPair) = -1.0_dp
                            elseif((i.eq.xInd).and.(j.eq.yInd)) then
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Two_Gamma_aaaa(xyPair) = 1.0_dp
                                Two_Gamma_abab(xyPair) = 1.0_dp
                            elseif((i.eq.yInd).and.(j.eq.xInd)) then
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Two_Gamma_aaaa(xyPair) = -1.0_dp
                                Two_Gamma_abba(xyPair) = -1.0_dp
                            endif
                        enddo
                    enddo
                else
                    ijInd = FindSqPairInd(i,j,tnmo)
                    read(unit_2RDM_read_aaaa,rec=ijInd) (Two_Gamma_aaaa(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abab,rec=ijInd) (Two_Gamma_abab(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abba,rec=ijInd) (Two_Gamma_abba(loop),loop=1,tnmo_sq)
                endif

                !Check antisymmetry requirements...
                if(tHFRDMs) then
                    Two_Gamma_ji_aaaa(:)=0.0_dp
                    Two_Gamma_ji_abab(:)=0.0_dp
                    Two_Gamma_ji_abba(:)=0.0_dp
                    do x=1,tnocc
                        xInd = OccOrbs(x)
                        do y=1,tnocc
                            yInd = OccOrbs(y)
                            if((j.eq.xInd).and.(i.eq.yInd).and.(xInd.eq.yInd)) then
                                !All four indices the same = 0
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Two_Gamma_ji_abab(xyPair) = 1.0_dp
                                Two_Gamma_ji_abba(xyPair) = -1.0_dp
                            elseif((j.eq.xInd).and.(i.eq.yInd)) then
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Two_Gamma_ji_aaaa(xyPair) = 1.0_dp
                                Two_Gamma_ji_abab(xyPair) = 1.0_dp
                            elseif((j.eq.yInd).and.(i.eq.xInd)) then
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Two_Gamma_ji_aaaa(xyPair) = -1.0_dp
                                Two_Gamma_ji_abba(xyPair) = -1.0_dp
                            endif
                        enddo
                    enddo
                else
                    !Read in 2R_ji^kl
                    jiInd = FindSqPairInd(j,i,tnmo)
                    read(unit_2RDM_read_aaaa,rec=jiInd) (Two_Gamma_ji_aaaa(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abab,rec=jiInd) (Two_Gamma_ji_abab(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abba,rec=jiInd) (Two_Gamma_ji_abba(loop),loop=1,tnmo_sq)
                endif

                ijInd = FindSqPairInd(i,j,tnmo)
                jiInd = FindSqPairInd(j,i,tnmo)
                SpinFac = SpinFac + 2.0_dp*Two_Gamma_aaaa(ijInd) - 2.0_dp*Two_Gamma_abab(ijInd) &
                    - 4.0_dp*Two_Gamma_abab(jiInd)

                !Check that 2R_ij^kl = - 2R_ji^kl (Remember to take spin component with it!)
                if(i.ne.j) then
                    do x=1,tnmo_sq

                        if(abs(Two_Gamma_ji_aaaa(x)+Two_Gamma_aaaa(x)).gt.AntisymDiff) then
                            AntisymDiff = abs(Two_Gamma_ji_aaaa(x)+Two_Gamma_aaaa(x))
                        endif
                        if(abs(Two_Gamma_ji_abab(x)+Two_Gamma_abba(x)).gt.AntisymDiff) then
                            AntisymDiff = abs(Two_Gamma_ji_abab(x)+Two_Gamma_abba(x))
                        endif
                        if(abs(Two_Gamma_ji_abba(x)+Two_Gamma_abab(x)).gt.1.D-3) then
                            AntisymDiff = abs(Two_Gamma_ji_abba(x)+Two_Gamma_abab(x))
                        endif

                        if(abs(Two_Gamma_ji_aaaa(x)+Two_Gamma_aaaa(x)).gt.1.D-3) then
                            write(6,*) "** WARNING ** Antisymmetry of 2RDM out by: ",   &
                                abs(Two_Gamma_aaaa(x)+Two_Gamma_aaaa(x))
                            if(abs(Two_Gamma_ji_aaaa(x)+Two_Gamma_aaaa(x)).gt.5.D-2) then
                                write(6,*) i,j,x,Two_Gamma_ji_aaaa(x),Two_Gamma_aaaa(x)
                                call warning(t_r,"Antisymmetry requirement of 2RDM not satisfied - aaaa")
                            endif
                        endif
                        if(abs(Two_Gamma_ji_abab(x)+Two_Gamma_abba(x)).gt.1.D-3) then
                            write(6,*) "** WARNING ** Antisymmetry of 2RDM out by: ",   &
                                abs(Two_Gamma_abab(x)+Two_Gamma_abba(x))
                            if(abs(Two_Gamma_ji_abab(x)+Two_Gamma_abba(x)).gt.5.D-2) then
                                write(6,*) i,j,x,Two_Gamma_ji_abab(x),Two_Gamma_abba(x)
                                call warning(t_r,"Antisymmetry requirement of 2RDM not satisfied - abab")
                            endif
                        endif
                        if(abs(Two_Gamma_ji_abba(x)+Two_Gamma_abab(x)).gt.1.D-3) then
                            write(6,*) "** WARNING ** Antisymmetry of 2RDM out by: ",   &
                                abs(Two_Gamma_abba(x)+Two_Gamma_abab(x))
                            if(abs(Two_Gamma_ji_abba(x)+Two_Gamma_abab(x)).gt.5.D-2) then
                                write(6,*) i,j,x,Two_Gamma_ji_abba(x),Two_Gamma_abab(x)
                                call warning(t_r,"Antisymmetry requirement of 2RDM not satisfied - abba")
                            endif
                        endif
                    enddo

                    !Also check that 2R_ij^kl = -2R_ij^lk (Remember again to take spins with you when swapping!)
                    do x=1,tnmo
                        do y=1,tnmo
                            if(x.eq.y) cycle
                            xyPair = FindSqPairInd(x,y,tnmo)
                            in_Pair = FindSqPairInd(y,x,tnmo)

                            if(abs(Two_Gamma_aaaa(xyPair)+Two_Gamma_aaaa(in_pair)).gt.AntisymDiff) then
                                AntisymDiff = abs(Two_Gamma_aaaa(xyPair)+Two_Gamma_aaaa(in_pair))
                            endif
                            if(abs(Two_Gamma_abab(xyPair)+Two_Gamma_abba(in_pair)).gt.AntisymDiff) then
                                AntisymDiff = abs(Two_Gamma_abab(xyPair)+Two_Gamma_abba(in_pair))
                            endif
                            if(abs(Two_Gamma_abba(xyPair)+Two_Gamma_abab(in_pair)).gt.AntisymDiff) then
                                AntisymDiff = abs(Two_Gamma_abba(xyPair)+Two_Gamma_abab(in_pair))
                            endif

                            if(abs(Two_Gamma_aaaa(xyPair)+Two_Gamma_aaaa(in_pair)).gt.1.D-3) then
                                write(6,*) "** WARNING ** Antisymmetry of 2RDM out by: ",   &
                                    abs(Two_Gamma_aaaa(xyPair)+Two_Gamma_aaaa(in_pair))
                                if(abs(Two_Gamma_aaaa(xyPair)+Two_Gamma_aaaa(in_pair)).gt.1.D-2) then
                                    call warning(t_r,"Antisymmetry requirement of 2RDM not satisfied - aaaa 2")
                                endif
                            endif
                            if(abs(Two_Gamma_abab(xyPair)+Two_Gamma_abba(in_pair)).gt.1.D-3) then
                                write(6,*) "** WARNING ** Antisymmetry of 2RDM out by: ",   &
                                    abs(Two_Gamma_abab(xyPair)+Two_Gamma_abba(in_pair))
                                if(abs(Two_Gamma_abab(xyPair)+Two_Gamma_abba(in_pair)).gt.1.D-2) then
                                    write(6,*) i,j,x,y,Two_Gamma_abab(xyPair),Two_Gamma_abba(in_pair)
                                    call warning(t_r,"Antisymmetry requirement of 2RDM not satisfied - abab 2")
                                endif
                            endif
                            if(abs(Two_Gamma_abba(xyPair)+Two_Gamma_abab(in_pair)).gt.1.D-3) then
                                write(6,*) "** WARNING ** Antisymmetry of 2RDM out by: ",   &
                                    abs(Two_Gamma_abba(xyPair)+Two_Gamma_abab(in_pair))
                                if(abs(Two_Gamma_abba(xyPair)+Two_Gamma_abab(in_pair)).gt.1.D-2) then
                                    call warning(t_r,"Antisymmetry requirement of 2RDM not satisfied - abba 2")
                                endif
                            endif
                        enddo
                    enddo
                endif

                !Check linear relation
                do k=1,tnmo
                    !We explicitly convert to spin-orbitals
                    x = FindSqPairInd(k,j,tnmo)
                    !aaaa part (Also same for bbbb)
                    OneRDMCopy(i,k) = OneRDMCopy(i,k) + Two_Gamma_aaaa(x)
                    OneRDMCopy(k,i) = OneRDMCopy(k,i) + Two_Gamma_aaaa(x)

                    !abab part
                    OneRDMCopy(i,k) = OneRDMCopy(i,k) + Two_Gamma_abab(x)
                    OneRDMCopy(k,i) = OneRDMCopy(k,i) + Two_Gamma_abab(x)
                    !abba part can not contribute due to spin symmetry
                enddo

                !Read in (G_ij)^kl
                in_pair = FindSqPairInd(i,j,tntmo)  !Complete space orbital pairs are stored
                read(unit_g_stored,rec=in_pair) (g_in(loop),loop=1,tnmo_sq)

                !Sum over kl, and EVar = <ij|kl> * 2RDM(k,l,i,j)
                !Remember that <ij|kl> is not antisymmetrised, but does have spin symmetry - hence no abba term
                VarE = VarE + DDOT(tnmo_sq,g_in,1,Two_Gamma_aaaa,1) 
                VarE = VarE + DDOT(tnmo_sq,g_in,1,Two_Gamma_abab,1) 
                
                !Check trace condition - sum in 2RDM_ij^ij - be careful about spin! abba term cannot contribute.
                !Other terms x 2 for bbbb and baba terms.
                xyPair = FindSqPairInd(i,j,tnmo)
                TwoRDMTrace = TwoRDMTrace + Two_Gamma_aaaa(xyPair) * 2.0_dp
                TwoRDMTrace = TwoRDMTrace + Two_Gamma_abab(xyPair) * 2.0_dp 

            enddo
        enddo

        !4S^2 + 4S = \sum_ij G_iaja^iaja + G_ibjb^ibjb - 2*G_iajb^iajb - 4*G_iajb^jaib + 3N
        SpinFac = SpinFac + 3.0_dp*NEl
        SpinFac = SpinFac / 4.0_dp
        !Solve quadratic equation for S
        if((1.0_dp + 4.0_dp*SpinFac).le.0.0_dp) then
            call warning(t_r,"Complex spin calculated from density matrices!")
        else
            SpinA = (-1.0_dp + sqrt(1.0_dp + 4.0_dp*SpinFac))/2.0_dp
            SpinB = (-1.0_dp - sqrt(1.0_dp + 4.0_dp*SpinFac))/2.0_dp
            if(max(SpinA,SpinB).lt.-1.D-7) then
                write(6,"(A,2F12.8)") "Spin: ",max(SpinA,SpinB),min(SpinA,SpinB)
                call warning(t_r,"Negative spin states found!")
            elseif(min(SpinA,SpinB).gt.0.D0) then
                write(6,"(A,2F10.5)") "Two possible spin states from density matrices? : ",SpinA,SpinB
            else
                write(6,"(A,F10.7)") "Spin state of wavefunction calculated from density matrices as ", &
                    max(SpinA,SpinB)
!                write(6,"(A,2F30.20)") "Spin state of wavefunction calculated from density matrices as ", &
!                    max(SpinA,SpinB),min(SpinA,SpinB)
            endif
        endif

        VarE = VarE + pot_nuc

        LinearRelDiff = 0.0_dp
        do i=1,tnmo
            do j=1,tnmo
                !elements are always 2x too large for a reason I haven't quite worked out.
                OneRDMCopy(i,j) = OneRDMCopy(i,j)/2.0_dp 

!                if(abs(OneRDM(i,j)).gt.1.D-8) then
!                    write(6,"(2I5,3F15.9)") i,j,(OneRDMCopy(i,j)/(real(NEl-1,dp)))/OneRDM(i,j),OneRDMCopy(i,j)/(real(NEl-1,dp)),OneRDM(i,j)
!                endif

                if((abs((OneRDMCopy(i,j)/(real(NEl-1,dp)))-OneRDM(i,j))).gt.LinearRelDiff) then
                    LinearRelDiff=abs((OneRDMCopy(i,j)/(real(NEl-1,dp)))-OneRDM(i,j))
                    if((OneRDM(i,j).gt.1.D-7).and.(abs(((OneRDMCopy(i,j)/(real(NEl-1,dp)))/OneRDM(i,j))-1.0_dp).gt.1.D-3)) then
                        call warning(t_r,"Integrated 2RDM and 1RDM not consistent")
                        write(6,"(A,2I4,A,F13.9,A,F13.9)") "Orbital pair: ",i,j,"  Integrated element: ", &
                            OneRDMCopy(i,j)/(real(NEl-1,dp)),"  Actual element: ",OneRDM(i,j)
                    endif
                endif
            enddo
        enddo

        if(LinearRelDiff.gt.1.D-3) then
            write(6,"(A,F20.9)") "** WARNING ** Integrated 2RDM differs from read in 1RDM by a maximum of: ", LinearRelDiff 
            write(6,"(A)") "RDMs may not be sufficiently converged"
        endif
        if(.not.tHFRDMs) write(6,"(A,F20.15)") "Linear relation between RDMs satisfied to accuracy of ",LinearRelDiff
        if(tIntegrateTwoRDM.and.(LinearRelDiff.gt.1.D-8)) then
            call stop_all(t_r,          &
         "1-RDM obtained from integration of 2-RDM, but linear relation between RDMs not satisfied - integration incorrect")
        endif
        if(abs(TwoRDMTrace-real(NEl*(NEl-1),dp)).gt.1.0D-7) then
            write(6,*) TwoRDMTrace, real(NEl*(NEl-1),dp)
            call warning(t_r,"Trace condition on 2RDM not satisfied")
        else
            if(.not.tHFRDMs) then
                write(6,"(A)") "Trace condition on 2RDM satisfied"
                write(6,"(A)") "Hermiticity of 2RDM not checked but assumed"
            endif
        endif
        if(AntisymDiff.gt.1.D-3) then
            write(6,"(A,F20.9)") "** WARNING ** Antisymmetry requirement of 2RDM out by: ",AntisymDiff
            write(6,"(A)") "RDMs may not be sufficiently converged"
        else
            if(.not.tHFRDMs) write(6,"(A,F20.9)") "Antisymmetry requirement of 2RDM satisfied to accuracy of ",AntisymDiff
        endif

        if(.not.tHFRDMs) then
            write(6,*) 
            write(6,"(A,F25.10)") "Variational energy of system calculated from density matrices to be: ",VarE
            write(6,*) ""
        else
            if(abs(VarE-ref_ene).gt.1.0d-8) then
                call stop_all(t_r,"Single determinant density matrix used, but energy from this is not the same as HF energy!")
            endif
        endif

        deallocate(g_in)
        if(.not.tExplicitCommTerm) then
            deallocate(tmat)    !I don't think we need this again?
        endif
        deallocate(OneRDMCopy)
        deallocate(Two_Gamma_aaaa,Two_Gamma_abab,Two_Gamma_abba)
        deallocate(Two_Gamma_ji_aaaa,Two_Gamma_ji_abab,Two_Gamma_ji_abba)
        close(unit_g_stored)
        if(.not.tHFRDMs) then
            close(unit_2RDM_read_aaaa)
            close(unit_2RDM_read_abab)
            close(unit_2RDM_read_abba)
        endif

    end subroutine CheckRDMProperties

    !Calculate the HF energy from two independant measures, and compare to the reference energy read in from DALTON
    !(if a dalton calculation).
    !In addition, check various properties of the density matrices, and calculate the variational energy from them.
    subroutine CalcHFEnergy()
        use input_data, only : tDaltonFormat
        implicit none
        integer :: unit_g_stored,ierr
        real(dp), allocatable :: g_in(:)
        character(*), parameter :: t_r="CalcHFEnergy"
        real(dp) :: SumFock,SumOneE,TwoElSum,FockHF,tmatHF
        integer :: i,Orb,j,Orbj,in_pair,loop,pair,pair2
        
        unit_g_stored=get_free_unit()
        open(unit_g_stored,file='G_Dir_ccoo',status='old',form='unformatted',access='direct',    &
            action='read',recl=reclen_nmo_sq)
        allocate(g_in(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"allocation err")

        SumFock=0.0_dp
        SumOneE=0.0_dp
        TwoElSum = 0.0_dp

        do i=1,tnocc
            Orb = OccOrbs(i)
            SumFock = SumFock + fockOrb(Orb,Orb)
            SumOneE = SumOneE + tmat(Orb,Orb)

            do j=1,tnocc
                Orbj = OccOrbs(j)

                in_pair = FindSqPairInd(Orb,Orbj,tntmo)  !Complete space orbital pairs are stored
                read(unit_g_stored,rec=in_pair) (g_in(loop),loop=1,tnmo_sq)

                pair = FindSqPairInd(Orb,Orbj,tnmo)
                pair2 = FindSqPairInd(Orbj,Orb,tnmo)

                TwoElSum = TwoElSum + 2.0_dp * g_in(pair) - g_in(pair2)
!                write(6,*) i,j,g_in(pair),g_in(pair2)
            enddo
        enddo
!        write(6,*) "TwoElSum: ",TwoElSum
        FockHF = 2.0_dp*SumFock - TwoElSum + pot_nuc   
        tmatHF = 2.0_dp*SumOneE + TwoElSum + pot_nuc

        if(abs(FockHF-tmatHF).gt.1.D-8) then
            write(6,*) FockHF,tmatHF,ref_ene
            call stop_all(t_r,"Two calculations of Hartree-Fock energy do not agree...")
        endif

        if(tDaltonFormat) then
            if(abs(FockHF-ref_ene).gt.1.D-8) then
                write(6,*) FockHF,tmatHF,ref_ene
                call stop_all(t_r,"Calculations of Hartree-Fock energy do not agree with value read in from DALTON...")
            endif
        else
            !For the solid-state system, the HF energy is not read in initially.
            ref_ene = tmatHF
        endif

        write(6,*) ""
        write(6,"(A,F20.10)") "Hartree-Fock energy correctly calculated to be: ",FockHF
        write(6,*) ""
!        write(6,*) "Fock bounds: ",lbound(fockorb,1),ubound(fockorb,1),lbound(fockorb,2),ubound(fockorb,2)

        close(unit_g_stored)
        deallocate(g_in)

    end subroutine CalcHFEnergy

!Calculate the MP2 energy! In terms of spatial orbitals:
! \sum_ijab [ 2*<ij|ab><ab|ij> - <ij|ab><ab|ji> ] / E_i + E_j - E_a - E_b
! = \sum_ijab [ 2*|<ab|ij>|^2 - <ab|ij><ab|ji> ] / E_i + E_j - E_a - E_b
    subroutine CalcMP2Energy()
        implicit none
        integer :: unit_g_stored,ierr,i,j,a,b,in_pair,loop,Orbi,Orbj,abijInd,abjiInd
        real(dp), allocatable :: g_in(:) 
        character(len=*), parameter :: t_r='CalcMP2Energy'
        
        unit_g_stored=get_free_unit()
        open(unit_g_stored,file='G_Dir_ccoo',status='old',form='unformatted',access='direct',    &
            action='read',recl=reclen_nmo_sq)
        allocate(g_in(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"allocation err")

        !Loop over virtual indices
        do a=1,tnmo
            if(IsOrbOcc(a)) cycle
            do b=1,tnmo
                if(IsOrbOcc(b)) cycle
                
                in_pair = FindSqPairInd(a,b,tntmo)  !Complete space orbital pairs are stored
                read(unit_g_stored,rec=in_pair) (g_in(loop),loop=1,tnmo_sq)

                !Now loop over occ indices

                do i=1,tnocc
                    Orbi = OccOrbs(i)
                    if(IsOrbFrz(Orbi)) cycle
                    do j=1,tnocc
                        Orbj = OccOrbs(j)
                        if(IsOrbFrz(Orbj)) cycle

                        !<ab|ij> index
                        abijInd = FindSqPairInd(Orbi,Orbj,tnmo)
                        abjiInd = FindSqPairInd(Orbj,Orbi,tnmo)

                        EMP2 = EMP2 + ((2.0_dp*(g_in(abijInd)**2.0_dp) - g_in(abijInd)*g_in(abjiInd)) /  &
                                (fockOrb(Orbi,Orbi) + fockOrb(Orbj,Orbj) - fockOrb(a,a) - fockOrb(b,b)))

                    enddo
                enddo
            enddo
        enddo
        deallocate(g_in)
        close(unit_g_stored)

        write(6,"(A,F25.10)") "MP2 energy calculated to be: ",EMP2 + ref_ene

    end subroutine CalcMP2Energy

end module Calc
