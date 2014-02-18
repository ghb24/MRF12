module int_management

use utils, only: get_free_unit
use basis_data
use basis, only: FindSqPairInd,FindCabsOrbPairInd
use errors, only: stop_all,warning
use input_data

    ! This module is concerned with the reading in of integrals,
    ! their storage on disk, and subsequent recall.

    implicit none

    contains

!Calculate the Phi term for the overlap matrix.
!This should reduce down to Phi_pq^rs = 2(e_p + e_q) \delta_pq^rs in the single reference case?
!This will be written out to disk for each spin type, just as the two RDM is
!It is written out as (Phi_pq)^rs, with each pq in a different record.
    subroutine CalcPhi()
        implicit none
        real(dp) , allocatable :: FGamma(:,:),GammaFGamma(:,:),GammaF(:,:),Z_pqrs_aaaa(:)
        real(dp) , allocatable :: Z_pqrs_abab(:),Z_pqrs_abba(:),TwoRDM_aaaa(:),TwoRDM_abab(:)
        real(dp) , allocatable :: TwoRDM_abba(:),TwoRDM_qp_aaaa(:),TwoRDM_qp_abab(:),TwoRDM_qp_abba(:)
        real(dp) , allocatable :: Z_qprs_aaaa(:),Z_qprs_abab(:),Z_qprs_abba(:),Phi_pqrs_aaaa(:)
        real(dp) , allocatable :: Phi_pqrs_abab(:),Phi_pqrs_abba(:)
        integer :: ierr,unit_Phi_write_aaaa,unit_Phi_write_abab,unit_Phi_write_abba
        integer :: unit_2RDM_read_aaaa,unit_2RDM_read_abab,unit_2RDM_read_abba,ispin
        integer :: p,q,pq_pair,qp_pair,r,s,rs_pair,rInd,sInd,loop,unit_Lambda_Tld_aaaa
        integer :: unit_Lambda_Tld_abab,unit_Lambda_Tld_abba,u
        character(*), parameter :: t_r='CalcPhi'

        !First, precalculate some of the intermediates in the terms.
        
        !FGamma: X_t^s = f_t^u G_u^s
        allocate(FGamma(tnmo,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Err")
        call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,GenFockOrb(:,:),tnmo,OneRDM,tnmo,0.0_dp,FGamma,tnmo)
        !X_t^s (t fast)

        !GammaFGamma: U_q^s = G_q^t f_t^u G_u^s
        allocate(GammaFGamma(tnmo,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Err")
        call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,OneRDM,tnmo,FGamma,tnmo,0.0_dp,GammaFGamma,tnmo)
        !U_q^s (q fast)

        !GammaF_1: X_u^s = G_t^s f_u^t
        allocate(GammaF(tnmo,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Err")
        call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,OneRDM,tnmo,GenFockOrb,tnmo,0.0_dp,GammaF,tnmo)
        !X_u^s (s fast).  This is the same as G_p^t f_t^u.

        !Open files to store Phi term
        unit_Phi_write_aaaa=get_free_unit()
        open(unit_Phi_write_aaaa,file='Phi_Dir_aaaa',status='unknown',form='unformatted',    &
            access='direct',action='write',recl=reclen_nmo_sq)
        unit_Phi_write_abab=get_free_unit()
        open(unit_Phi_write_abab,file='Phi_Dir_abab',status='unknown',form='unformatted',    &
            access='direct',action='write',recl=reclen_nmo_sq)
        unit_Phi_write_abba=get_free_unit()
        open(unit_Phi_write_abba,file='Phi_Dir_abba',status='unknown',form='unformatted',    &
            access='direct',action='write',recl=reclen_nmo_sq)

        !Open files to read in Lambda_Tld
        !Delete these files at the end as they will not need to be used again.
        unit_Lambda_Tld_aaaa=get_free_unit()
        open(unit_Lambda_Tld_aaaa,file='L_Tld_Dir_aaaa',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)
        unit_Lambda_Tld_abab=get_free_unit()
        open(unit_Lambda_Tld_abab,file='L_Tld_Dir_abab',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)
        unit_Lambda_Tld_abba=get_free_unit()
        open(unit_Lambda_Tld_abba,file='L_Tld_Dir_abba',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)

        if(.not.tHFRDMs) then
            !Open units to read 2RDM
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

        allocate(Z_pqrs_aaaa(tnmo_sq),stat=ierr)
        allocate(Z_pqrs_abab(tnmo_sq),stat=ierr)
        allocate(Z_pqrs_abba(tnmo_sq),stat=ierr)
        allocate(Z_qprs_aaaa(tnmo_sq),stat=ierr)
        allocate(Z_qprs_abab(tnmo_sq),stat=ierr)
        allocate(Z_qprs_abba(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Err")
        allocate(TwoRDM_aaaa(tnmo_sq),stat=ierr)
        allocate(TwoRDM_abab(tnmo_sq),stat=ierr)
        allocate(TwoRDM_abba(tnmo_sq),stat=ierr)
        allocate(TwoRDM_qp_aaaa(tnmo_sq),stat=ierr)
        allocate(TwoRDM_qp_abab(tnmo_sq),stat=ierr)
        allocate(TwoRDM_qp_abba(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Err")
        allocate(Phi_pqrs_aaaa(tnmo_sq),stat=ierr)
        allocate(Phi_pqrs_abab(tnmo_sq),stat=ierr)
        allocate(Phi_pqrs_abba(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Err")

        write(6,'(F10.3,A)') 3.0_dp*(real(tnmo,dp)**4.0_dp)*8.0_dp*9.31322575D-10,   &
            " Gb required to write out Phi_Dir_xxxx files."
        call flush(6)

        !Loop over p < q
        do p=1,tnmo
            do q=p,tnmo

                !q fast in records.
                pq_pair = FindSqPairInd(p,q,tnmo)   !This is the unit that this phi term is written out to.
                qp_pair = FindSqPairInd(q,p,tnmo)
                    
                !Read in / calculate 2RDM for eeach spin type.
                if(tHFRDMs) then
                    TwoRDM_aaaa(:) = 0.0_dp
                    TwoRDM_abab(:) = 0.0_dp
                    TwoRDM_abba(:) = 0.0_dp
                    TwoRDM_qp_aaaa(:) = 0.0_dp
                    TwoRDM_qp_abab(:) = 0.0_dp
                    TwoRDM_qp_abba(:) = 0.0_dp
                    if(IsOrbOcc(p).and.IsOrbOcc(q)) then
                        do r=1,tnocc
                            rInd=OccOrbs(r)
                            do s=1,tnocc
                                sInd=OccOrbs(s)
                                rs_pair = FindSqPairInd(rInd,sInd,tnmo)
                                !Fill up the 2RDM for the (G_pq)^rs array
                                if((p.eq.rInd).and.(q.eq.sInd).and.(p.ne.q)) then
                                    TwoRDM_aaaa(rs_pair) = 1.0_dp
                                    TwoRDM_abab(rs_pair) = 1.0_dp
                                elseif((p.eq.sInd).and.(q.eq.rInd).and.(p.ne.q)) then
                                    TwoRDM_aaaa(rs_pair) = -1.0_dp
                                    TwoRDM_abba(rs_pair) = -1.0_dp
                                elseif((p.eq.sInd).and.(q.eq.rInd).and.(p.eq.q)) then
                                    TwoRDM_abab(rs_pair) = 1.0_dp
                                    TwoRDM_abba(rs_pair) = -1.0_dp
                                endif
                                if((q.eq.rInd).and.(p.eq.sInd).and.(p.ne.q)) then
                                    TwoRDM_qp_aaaa(rs_pair) = 1.0_dp
                                    TwoRDM_qp_abab(rs_pair) = 1.0_dp
                                elseif((q.eq.sInd).and.(p.eq.rInd).and.(p.ne.q)) then
                                    TwoRDM_qp_aaaa(rs_pair) = -1.0_dp
                                    TwoRDM_qp_abba(rs_pair) = -1.0_dp
                                elseif((q.eq.sInd).and.(p.eq.rInd).and.(p.eq.q)) then
                                    TwoRDM_qp_abab(rs_pair) = 1.0_dp
                                    TwoRDM_qp_abba(rs_pair) = -1.0_dp
                                endif
                            enddo
                        enddo
                    endif
                else
                    read(unit_2RDM_read_aaaa,rec=pq_pair) (TwoRDM_aaaa(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_aaaa,rec=qp_pair) (TwoRDM_qp_aaaa(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abab,rec=pq_pair) (TwoRDM_abab(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abab,rec=qp_pair) (TwoRDM_qp_abab(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abba,rec=pq_pair) (TwoRDM_abba(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abba,rec=qp_pair) (TwoRDM_qp_abba(loop),loop=1,tnmo_sq)
                endif

!                do r=1,tnmo
!                    do s=1,tnmo
!                        rs_pair = FindSqPairInd(r,s,tnmo)
!                        write(47,"(4I5,2G20.10)") p,q,r,s,TwoRDM_aaaa(rs_pair),TwoRDM_qp_aaaa(rs_pair)
!                    enddo
!                enddo

                !calculate (Z_pq)^rs and (Z_qp)^rs at the same time
                do r=1,tnmo
                    do s=1,tnmo
                        rs_pair = FindSqPairInd(r,s,tnmo)   !s fast within record
                        !The abba term doesn't contribute here due to spin integration.
                        Z_pqrs_aaaa(rs_pair) = OneRDM(p,r)*GammaFGamma(q,s)
                        Z_pqrs_abab(rs_pair) = OneRDM(p,r)*GammaFGamma(q,s)
                        Z_pqrs_abba(rs_pair) = 0.0_dp

                        Z_qprs_aaaa(rs_pair) = OneRDM(q,r)*GammaFGamma(p,s)
                        Z_qprs_abab(rs_pair) = OneRDM(q,r)*GammaFGamma(p,s)
                        Z_qprs_abba(rs_pair) = 0.0_dp
                    enddo
                enddo

                !First term of Z_pqrs is correct for all spins.
                
                !Second term in (Z_pq)^rs = 1/2 Gamma_t^s f_u^t Lambda_pq^ru
                !First, the 2RDM part of the cumulant, and add it to Z
                !Remember - in the 2RDM, the second index is fast!
                call DGEMM('N','N',tnmo,tnmo,tnmo,0.5_dp,GammaF,tnmo,TwoRDM_aaaa,tnmo,1.0_dp,Z_pqrs_aaaa,tnmo)
                call DGEMM('N','N',tnmo,tnmo,tnmo,0.5_dp,GammaF,tnmo,TwoRDM_abab,tnmo,1.0_dp,Z_pqrs_abab,tnmo)
                call DGEMM('N','N',tnmo,tnmo,tnmo,0.5_dp,GammaF,tnmo,TwoRDM_abba,tnmo,1.0_dp,Z_pqrs_abba,tnmo)

                call DGEMM('N','N',tnmo,tnmo,tnmo,0.5_dp,GammaF,tnmo,TwoRDM_qp_aaaa,tnmo,1.0_dp,Z_qprs_aaaa,tnmo)
                call DGEMM('N','N',tnmo,tnmo,tnmo,0.5_dp,GammaF,tnmo,TwoRDM_qp_abab,tnmo,1.0_dp,Z_qprs_abab,tnmo)
                call DGEMM('N','N',tnmo,tnmo,tnmo,0.5_dp,GammaF,tnmo,TwoRDM_qp_abba,tnmo,1.0_dp,Z_qprs_abba,tnmo)
                !This means that s is fast in the output.
                
                !Now for the 1RDM parts of the second term.
                do r=1,tnmo
                    do s=1,tnmo
                        rs_pair = FindSqPairInd(r,s,tnmo)
                        Z_pqrs_aaaa(rs_pair) = Z_pqrs_aaaa(rs_pair) + (GammaFGamma(p,s)*OneRDM(q,r)/2.0_dp) - &
                            (GammaFGamma(q,s)*OneRDM(p,r)/2.0_dp)
                        !Some bits don't contribute for the abab and abba spins.
                        Z_pqrs_abab(rs_pair) = Z_pqrs_abab(rs_pair) - (GammaFGamma(q,s)*OneRDM(p,r)/2.0_dp)
                        Z_pqrs_abba(rs_pair) = Z_pqrs_abba(rs_pair) + (GammaFGamma(p,s)*OneRDM(q,r)/2.0_dp)

                        Z_qprs_aaaa(rs_pair) = Z_qprs_aaaa(rs_pair) + (GammaFGamma(q,s)*OneRDM(p,r)/2.0_dp) - &
                            (GammaFGamma(p,s)*OneRDM(q,r)/2.0_dp)
                        !Some bits don't contribute for the abab and abba spins.
                        Z_qprs_abab(rs_pair) = Z_qprs_abab(rs_pair) - (GammaFGamma(p,s)*OneRDM(q,r)/2.0_dp)
                        Z_qprs_abba(rs_pair) = Z_qprs_abba(rs_pair) + (GammaFGamma(q,s)*OneRDM(p,r)/2.0_dp)
                    enddo
                enddo

                !Now for the third term. This is stored on disk.
                !Overwrite the TwoRDM arrays - these aren't needed again.
                read(unit_Lambda_Tld_aaaa,rec=pq_pair) (TwoRDM_aaaa(loop),loop=1,tnmo_sq)
                read(unit_Lambda_Tld_abab,rec=pq_pair) (TwoRDM_abab(loop),loop=1,tnmo_sq)
                read(unit_Lambda_Tld_abba,rec=pq_pair) (TwoRDM_abba(loop),loop=1,tnmo_sq)
                read(unit_Lambda_Tld_aaaa,rec=qp_pair) (TwoRDM_qp_aaaa(loop),loop=1,tnmo_sq)
                read(unit_Lambda_Tld_abab,rec=qp_pair) (TwoRDM_qp_abab(loop),loop=1,tnmo_sq)
                read(unit_Lambda_Tld_abba,rec=qp_pair) (TwoRDM_qp_abba(loop),loop=1,tnmo_sq)

                Z_pqrs_aaaa(:) = Z_pqrs_aaaa(:) + TwoRDM_aaaa(:)
                Z_pqrs_abab(:) = Z_pqrs_abab(:) + TwoRDM_abab(:)
                Z_pqrs_abba(:) = Z_pqrs_abba(:) + TwoRDM_abba(:)
                
                Z_qprs_aaaa(:) = Z_qprs_aaaa(:) + TwoRDM_qp_aaaa(:)
                Z_qprs_abab(:) = Z_qprs_abab(:) + TwoRDM_qp_abab(:)
                Z_qprs_abba(:) = Z_qprs_abba(:) + TwoRDM_qp_abba(:)


                !Fourth term. -G_p^r LambdaBar_q^s
                !This doesn't contribute to the abba term!
                do r=1,tnmo
                    do s=1,tnmo
                        rs_pair = FindSqPairInd(r,s,tnmo)
                        Z_pqrs_aaaa(rs_pair) = Z_pqrs_aaaa(rs_pair) - OneRDM(p,r)*Lambda_Bar(q,s)
                        Z_pqrs_abab(rs_pair) = Z_pqrs_abab(rs_pair) - OneRDM(p,r)*Lambda_Bar(q,s)

                        Z_qprs_aaaa(rs_pair) = Z_qprs_aaaa(rs_pair) - OneRDM(q,r)*Lambda_Bar(p,s)
                        Z_qprs_abab(rs_pair) = Z_qprs_abab(rs_pair) - OneRDM(q,r)*Lambda_Bar(p,s)
                    enddo
                enddo
                
                !(Z_pq)^rs is now calculated correctly. However, we still need to take correct permutations to get phi.

                Phi_pqrs_aaaa(:) = Z_pqrs_aaaa(:) - Z_qprs_aaaa(:)
                call TransposeIntegralArr_inplace(Z_pqrs_aaaa,tnmo)
                call TransposeIntegralArr_inplace(Z_qprs_aaaa,tnmo)
                Phi_pqrs_aaaa(:) = Phi_pqrs_aaaa(:) - Z_pqrs_aaaa(:) + Z_qprs_aaaa(:)

                write(unit_Phi_write_aaaa,rec=pq_pair) (Phi_pqrs_aaaa(loop),loop=1,tnmo_sq)
                write(unit_Phi_write_aaaa,rec=qp_pair) (-Phi_pqrs_aaaa(loop),loop=1,tnmo_sq)  !(Phi_pq)^rs = -(Phi_qp)^rs

                !Now for mixed spin parts...
                Phi_pqrs_abab(:) = Z_pqrs_abab(:) - Z_qprs_abba(:)
                Phi_pqrs_abba(:) = Z_pqrs_abba(:) - Z_qprs_abab(:)
                call TransposeIntegralArr_inplace(Z_pqrs_abab,tnmo)
                call TransposeIntegralArr_inplace(Z_pqrs_abba,tnmo)
                call TransposeIntegralArr_inplace(Z_qprs_abab,tnmo)
                call TransposeIntegralArr_inplace(Z_qprs_abba,tnmo)
                Phi_pqrs_abab(:) = Phi_pqrs_abab(:) - Z_pqrs_abba(:) + Z_qprs_abab(:)
                Phi_pqrs_abba(:) = Phi_pqrs_abba(:) - Z_pqrs_abab(:) + Z_qprs_abba(:)
                
                write(unit_Phi_write_abab,rec=pq_pair) (Phi_pqrs_abab(loop),loop=1,tnmo_sq)
                write(unit_Phi_write_abba,rec=qp_pair) (-Phi_pqrs_abab(loop),loop=1,tnmo_sq)  !(Phi_pq)^rs = -(Phi_qp)^rs

                write(unit_Phi_write_abba,rec=pq_pair) (Phi_pqrs_abba(loop),loop=1,tnmo_sq)
                write(unit_Phi_write_abab,rec=qp_pair) (-Phi_pqrs_abba(loop),loop=1,tnmo_sq)  !(Phi_pq)^rs = -(Phi_qp)^rs

            enddo
        enddo

        close(unit_Lambda_Tld_aaaa,status='delete')
        close(unit_Lambda_Tld_abab,status='delete')
        close(unit_Lambda_Tld_abba,status='delete')

        close(unit_Phi_write_aaaa)
        close(unit_Phi_write_abab)
        close(unit_Phi_write_abba)
        if(.not.tHFRDMs) then
            close(unit_2RDM_read_aaaa)
            close(unit_2RDM_read_abab)
            close(unit_2RDM_read_abba)
        endif

!!Test!
!        unit_Phi_write_aaaa=get_free_unit()
!        open(unit_Phi_write_aaaa,file='Phi_Dir_aaaa',status='unknown',form='unformatted',    &
!            access='direct',action='read',recl=reclen_nmo_sq)
!        unit_Phi_write_abab=get_free_unit()
!        open(unit_Phi_write_abab,file='Phi_Dir_abab',status='unknown',form='unformatted',    &
!            access='direct',action='read',recl=reclen_nmo_sq)
!        unit_Phi_write_abba=get_free_unit()
!        open(unit_Phi_write_abba,file='Phi_Dir_abba',status='unknown',form='unformatted',    &
!            access='direct',action='read',recl=reclen_nmo_sq)
!        do p=1,tnmo
!            do q=1,tnmo
!                pq_pair = FindSqPairInd(p,q,tnmo)
!                read(unit_Phi_write_aaaa,rec=pq_pair) (Phi_pqrs_aaaa(loop),loop=1,tnmo_sq)
!                read(unit_Phi_write_abab,rec=pq_pair) (Phi_pqrs_abab(loop),loop=1,tnmo_sq)
!                read(unit_Phi_write_abba,rec=pq_pair) (Phi_pqrs_abba(loop),loop=1,tnmo_sq)
!                do r=1,tnmo
!                    do s=1,tnmo
!                        rs_pair = FindSqPairInd(r,s,tnmo)
!                        write(71,"(G25.10,4I5)") Phi_pqrs_aaaa(rs_pair),p,q,r,s
!                        write(72,"(G25.10,4I5)") Phi_pqrs_abab(rs_pair),p,q,r,s
!                        write(73,"(G25.10,4I5)") Phi_pqrs_abba(rs_pair),p,q,r,s
!                    enddo
!                enddo
!            enddo
!        enddo
!        close(unit_Phi_write_aaaa)
!        close(unit_Phi_write_abab)
!        close(unit_Phi_write_abba)

        deallocate(FGamma,GammaFGamma,GammaF)
        deallocate(Z_pqrs_aaaa,Z_pqrs_abab,Z_pqrs_abba)
        deallocate(Z_qprs_aaaa,Z_qprs_abab,Z_qprs_abba)
        deallocate(TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
        deallocate(TwoRDM_qp_aaaa,TwoRDM_qp_abab,TwoRDM_qp_abba)
        deallocate(Phi_pqrs_aaaa,Phi_pqrs_abab,Phi_pqrs_abba)

    end subroutine CalcPhi

!Calculate the generalised K and FpK matrices
!K matrix is one matrix over complete:complete space.
!fpk matrix is one matrix over complete:orbital space.
    subroutine CalcGenKandFpK()
        implicit none
        integer :: ierr,loop,in_pair,unit_g_stored,p_prime,r_prime
        real(dp) , allocatable :: r_in(:)
        real(dp) :: DDOT
        character(*), parameter :: t_r='CalcGenKandFpK'

        if(tHFRDMs) then
            write(6,"(A)") "Constructing exchange and fpk matrix..."
        else
            write(6,"(A)") "Constructing generalised exchange and fpk matrix..."
        endif

        !We store the alpha:alpha components of the matrix only
        allocate(GenExch(tntmo,tntmo),stat=ierr)
        allocate(GenFpK_CompOrb(tntmo,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc err")
        GenExch(:,:) = 0.0_dp
        !Just add the exchange part to fock matrix, rather than constructing h + J.
        if(.not.allocated(GenFockCompOrb)) call stop_all(t_r,"Gen fock matrix not constructed")
        GenFpK_CompOrb(:,:) = GenFockCompOrb(:,:) 
            
        allocate(r_in(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"allocation err")
            
        unit_g_stored=get_free_unit()
        open(unit_g_stored,file='G_Dir_ccoo',status='old',form='unformatted',access='direct',action='read',recl=reclen_nmo_sq)

        do p_prime=1,tntmo
            do r_prime=1,tntmo

                in_pair = FindSqPairInd(p_prime,r_prime,tntmo)  !Complete space orbital pairs are stored
                read(unit_g_stored,rec=in_pair) (r_in(loop),loop=1,tnmo_sq)
                !Get all (G_p'r')^qp with p fast

                !Contract with Gam_q^p
                GenExch(p_prime,r_prime) = DDOT(tnmo_sq,r_in,1,OneRDM(:,:),1)
                if(r_prime.le.tnmo) then
                    GenFpK_CompOrb(p_prime,r_prime) = GenFpK_CompOrb(p_prime,r_prime) + GenExch(p_prime,r_prime)
                endif
            enddo
        enddo

        deallocate(r_in)
        close(unit_g_stored)

    end subroutine CalcGenKandFpK

!Calculate the generalised fock matrices in the different orbital subspaces.
!If tHFRDMs is true, then just copy the fock matrix that was read in.
    subroutine CalcGenFock()
        implicit none
        integer :: unit_g_stored,unit_g_read,ierr,k,l,temp,qp_pair,ios,ibuf,pq_pair
        integer :: i_Ind,j_Ind,k_Ind,l_Ind,p,q,loop,in_pair,lq_pair,ql_pair,i,j
        integer(i2) , allocatable :: Indices(:,:)
        real(dp) , allocatable :: Buf(:),r_in(:),TempGInts(:,:,:)
        integer(i8) :: maxlength,length
        logical :: tCanonicalHF,tCanonicalVirtHF
        real(dp) :: MaxOccVirt,MaxOccCABSFock,MaxVirtCABSFock
        character(len=*), parameter :: t_r='CalcGenFock'
        integer :: nOrb,nAux,nElec,Ms2,nProp(3),PropBitLen,iSym
        real(dp) :: Z,GAM
        integer(i8) :: OrbSym(1000)
        logical :: UHF,UEG
        namelist /FCI/ nOrb,nAux,nElec,Ms2,OrbSym,iSym,UHF,UEG,GAM,PropBitLen,nProp

        if(tHFRDMs) then
            write(6,"(A)") "Constructing fock matrix..."
        else
            write(6,"(A)") "Constructing generalised fock matrix..."
        endif

        allocate(GenFockOrb(tnmo,tnmo),stat=ierr)
        allocate(GenFockOrbCABS(1:tnmo,tnmo+1:tntmo),stat=ierr)
        allocate(GenFockCABS(tnmo+1:tntmo,tnmo+1:tntmo),stat=ierr)
        allocate(GenFockCompOrb(tntmo,tnmo),stat=ierr)
        allocate(GenFockComp(tntmo,tntmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Gen Fock allocation err")

        if (tHFRDMs) then
            !Just copy the read in fock matrix to the generalised one
            GenFockOrb(:,:) = FockOrb(:,:)
            GenFockOrbCABS(:,:) = FockOrbCABS(:,:)
            GenFockCABS(:,:) = FockCABS(:,:)
            GenFockComp(1:tnmo,1:tnmo) = FockOrb(:,:)
            GenFockComp(tnmo+1:tntmo,tnmo+1:tntmo) = FockCABS(:,:)
            do i=1,tnmo
                do j=tnmo+1,tntmo
                    GenFockComp(i,j) = GenFockOrbCABS(i,j)
                    GenFockComp(j,i) = GenFockOrbCABS(i,j)
                enddo
            enddo
        else
            if(tExplicitCommTerm) then
                write(6,"(A)") "Cannot calculate the explicit commutator term, since not coded up for multireference case"
                write(6,"(A)") "However, this should be simple - just code up GenFockComp for general case"
                call stop_all(t_r,"Please disable calculation of the explicit commutator term")
            endif
            !Need to construct the generalised fock matrices
            GenFockOrb(:,:) = 0.0_dp 
            GenFockOrbCABS(:,:) = 0.0_dp
            GenFockCABS(:,:) = 0.0_dp

            unit_g_stored=get_free_unit()
            open(unit_g_stored,file='G_Dir_ccoo',status='old',form='unformatted',access='direct',action='read',recl=reclen_nmo_sq)
            
            unit_g_read=get_free_unit()
            if(tDaltonFormat) then
                open(unit_g_read,file='MO_G',status='old',form='unformatted',access='sequential',action='read')
                rewind(unit_g_read)
                read(unit_g_read) maxlength
                allocate(Indices(4,MaxLength),stat=ierr)
                if(ierr.ne.0) call stop_all(t_r,"Alloc err 1")
                allocate(Buf(MaxLength),stat=ierr)
                if(ierr.ne.0) call stop_all(t_r,"Alloc err 1a")
            else
                if(tReadBin) then
                    open(unit_g_read,file='GDUMPBIN',status='old',form='unformatted',action='read')
                else
                    open(unit_g_read,file='GDUMP',status='old',form='formatted',action='read')
                endif
                rewind(unit_g_read)
                length = 1
                if(.not.tReadBin) read(unit_g_read,fci)
            endif

            allocate(r_in(tnmo_sq),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"allocation err")

    !       f_kl = h_kl + [ <kp|lq> - <kp|ql> ] Gamma_qp
            !First, deal with Complete:Orbital space
            !K > L, therefore K will always be in CABS space if there is a cabs orbital
            do k=1,tntmo
                do l=1,tnmo

                    if((k.le.tnmo).and.(k.lt.l)) cycle

                    if(k.le.tnmo) then
                        GenFockOrb(k,l) = tmat(k,l)
                    else
                        !K is in CABS, and so must be the second index in GenFockOrbCABS
                        GenFockOrbCABS(l,k) = tmat(l,k)
                    endif

                    !Load needed g ints
                    !l is always in orbital space
                    r_in(:)=0.0_dp
                    do p=1,tnmo

                        in_pair = FindSqPairInd(k,p,tntmo)  !Complete space orbital pairs are stored
                        read(unit_g_stored,rec=in_pair) (r_in(loop),loop=1,tnmo_sq)

                        do q=1,tnmo
                            lq_pair = FindSqPairInd(l,q,tnmo)
                            ql_pair = FindSqPairInd(q,l,tnmo)

                            if(k.le.tnmo) then
                                GenFockOrb(k,l) = GenFockOrb(k,l) +     &
                                    (2.0_dp*r_in(lq_pair)*OneRDM(q,p)) - (r_in(ql_pair)*OneRDM(q,p))
                            else
                                !K is in CABS, and so must be the second index in GenFockOrbCABS
                                GenFockOrbCABS(l,k) = GenFockOrbCABS(l,k) + &
                                   (2.0_dp*r_in(lq_pair)*OneRDM(q,p)) - (r_in(ql_pair)*OneRDM(q,p))
                            endif

                        enddo
                    enddo

                    if(k.le.tnmo) then
                        GenFockOrb(l,k) = GenFockOrb(k,l)
                    endif

                enddo
            enddo

            !Now deal with Cabs:Cabs gen fock mat.
            !First, deal with exchange term and one-electron term
            do k=tnmo+1,tntmo
                do l=k,tntmo

                    GenFockCABS(k,l) = tmat(k,l)
                        
                    !Deal with exchange term first, for which we can use -<kl|qp> Gamma_qp
                    r_in(:)=0.0_dp
                    in_pair = FindSqPairInd(k,l,tntmo)  !Complete space orbital pairs are stored
                    read(unit_g_stored,rec=in_pair) (r_in(loop),loop=1,tnmo_sq)

                    do p=1,tnmo
                        do q=1,tnmo

                            qp_pair = FindSqPairInd(q,p,tnmo)

                            GenFockCABS(k,l) = GenFockCABS(k,l) - (r_in(qp_pair)*OneRDM(q,p))

                        enddo
                    enddo

                enddo
            enddo

            !Now for the difficult coulomb term
            !Store all < K p | l q >
            allocate(TempGInts(tnmo,tnmo,tnmo+1:tntmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err")

            do k=tnmo+1,tntmo

                !Store order is p,q,l
!                <Kp|lq>
!                <Kq|lp>
!                <pK|ql>
!                <qK|pl>
!                <lp|Kq>
!                <lq|Kp>
!                <pl|qK>
!                <ql|pK>

!However, if tComplexOrbs is true, then we only want to calculate:
!                <Kp|lq>
!                <pK|ql>
!                <lq|Kp>
!                <ql|pK>
                TempGInts(:,:,:) = 0.0_dp

                rewind(unit_g_read)
                if(tDaltonFormat) then
                    read(unit_g_read) maxlength
                else
                    if(.not.tReadBin) read(unit_g_read,fci)
                endif
                do while(.true.)
                    if(tDaltonFormat) then
                        read(unit_g_read,iostat=ios) length,Indices(1:4,1:length),Buf(1:length)
                    else
                        if(tReadBin) then
                            read(unit_g_read,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                        else
                            read(unit_g_read,*,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                        endif
                    endif
                    if(ios.gt.0) call stop_all(t_r,"Error reading integrals")
                    if((length.le.0).or.(ios.lt.0)) exit    !EOF

                    do ibuf = 1,length
                        if(tDaltonFormat) then
                            i_Ind=int(Indices(1,ibuf),i4)
                            j_Ind=int(Indices(3,ibuf),i4)
                            k_Ind=int(Indices(2,ibuf),i4)
                            l_Ind=int(Indices(4,ibuf),i4)
                            Z = Buf(ibuf)
                        else
                            if((i_ind.eq.0).or.(j_Ind.eq.0).or.(k_ind.eq.0).or.(l_ind.eq.0)) cycle  !We only want 2e-integrals
                        endif

                        if(i_Ind.eq.k) then
                            if((k_Ind.gt.tnmo).and.(max(j_Ind,l_Ind).le.tnmo)) then
                                TempGInts(j_Ind,l_Ind,k_Ind) = Z    !<Kp|lq>
                                if(.not.tComplexOrbs) TempGInts(l_Ind,j_Ind,k_Ind) = Z  !<Kq|lp>
                            endif
                        endif
                        if(j_Ind.eq.k) then
                            if((l_Ind.gt.tnmo).and.(max(i_Ind,k_Ind).le.tnmo)) then
                                TempGInts(i_Ind,k_Ind,l_Ind) = Z !<pK|ql>
                                if(.not.tComplexOrbs) TempGInts(k_Ind,i_Ind,l_Ind) = Z  !<qK|pl>
                            endif
                        endif
                        if(k_Ind.eq.k) then
                            if((i_Ind.gt.tnmo).and.(max(j_Ind,l_Ind).le.tnmo)) then
                                TempGInts(l_Ind,j_Ind,i_Ind) = Z !<lq|Kp>
                                if(.not.tComplexOrbs) TempGInts(j_Ind,l_Ind,i_Ind) = Z !<lp|Kq>
                            endif
                        endif
                        if(l_Ind.eq.k) then
                            if((j_Ind.gt.tnmo).and.(max(i_Ind,k_Ind).le.tnmo)) then
                                TempGInts(k_Ind,i_Ind,j_Ind) = Z !<ql|pK>
                                if(.not.tComplexOrbs) TempGInts(i_Ind,k_Ind,j_Ind) = Z !<pl|pK>
                            endif
                        endif

                    enddo
                enddo   !End do looping over buffers

                do l=k,tntmo    !Loop over l

                    do p=1,tnmo
                        do q=1,tnmo

                            !Store order is p,q,l
!                            <Kp|lq>
!                            GenFockCABS(k,l) = GenFockCABS(k,l) + 2.0_dp*TempGInts(q,p,l)*OneRDM(q,p)
                            GenFockCABS(k,l) = GenFockCABS(k,l) + 2.0_dp*TempGInts(p,q,l)*OneRDM(q,p)
                        enddo
                    enddo

                    GenFockCABS(l,k) = GenFockCABS(k,l)

                enddo

            enddo

!            !Commented out below is the old way of filling gen fock mat, which should still work
!            do k=1,tntmo
!                do l=1,tntmo
!
!                    if(min(k,l).le.tnmo) then
!                        !Orbital orbital space or one CABS index only.
!                        !The one CABS index will be fine, as long as we ensure that k is the CABS index. Swap
!                        !The indices if it is not
!                        if(l.gt.tnmo) then
!                            temp = k
!                            k_Ind = l
!                            l_Ind = temp
!                        else
!                            k_Ind = k
!                            l_Ind = l
!                        endif
!
!                        if(max(k_Ind,l_Ind).le.tnmo) then
!                            GenFockOrb(k_Ind,l_Ind) = tmat(k_Ind,l_Ind)
!                        else
!                            !K is in CABS, and so must be the second index in GenFockOrbCABS
!                            GenFockOrbCABS(l_Ind,k_Ind) = tmat(k_Ind,l_Ind)
!                        endif
!
!                        r_in(:)=0.0_dp
!                        do p=1,tnmo
!
!                            in_pair = FindSqPairInd(k_Ind,p,tntmo)  !Complete space orbital pairs are stored
!                            read(unit_g_stored,rec=in_pair) (r_in(loop),loop=1,tnmo_sq)
!
!                            do q=1,tnmo
!                                lq_pair = FindSqPairInd(l_Ind,q,tnmo)
!                                ql_pair = FindSqPairInd(q,l_Ind,tnmo)
!
!                                if(max(k_Ind,l_Ind).le.tnmo) then
!                                    GenFockOrb(k_Ind,l_Ind) = GenFockOrb(k_Ind,l_Ind) +     &
!                                        (2.0_dp*r_in(lq_pair)*OneRDM(q,p)) - (r_in(ql_pair)*OneRDM(q,p))
!                                else
!                                    !K is in CABS, and so must be the second index in GenFockOrbCABS
!                                    GenFockOrbCABS(l_Ind,k_Ind) = GenFockOrbCABS(l_Ind,k_Ind) + &
!                                       (2.0_dp*r_in(lq_pair)*OneRDM(q,p)) - (r_in(ql_pair)*OneRDM(q,p))
!                                endif
!
!                            enddo
!                        enddo
!
!                    else
!                        !CABS CABS space: Both k and l in CABS. Difficult coulomb term, since not stored in the G_Dir file.
!
!                        GenFockCABS(k,l) = tmat(k,l)
!
!                        !Deal with exchange term first, for which we can use -<kl|qp> Gamma_qp
!                        r_in(:)=0.0_dp
!                        in_pair = FindSqPairInd(k,l,tntmo)  !Complete space orbital pairs are stored
!                        read(unit_g_stored,rec=in_pair) (r_in(loop),loop=1,tnmo_sq)
!
!                        do p=1,tnmo
!                            do q=1,tnmo
!
!                                qp_pair = FindSqPairInd(q,p,tnmo)
!
!                                GenFockCABS(k,l) = GenFockCABS(k,l) - (r_in(qp_pair)*OneRDM(q,p))
!
!                            enddo
!                        enddo
!
!                        !Now for tricky coulomb term. <kp|lq>
!                        !Since 4traf only writes out < c c | c o > integrals for G, then there are only 4 possible permutations of how
!                        !the integrals needed could be read in.
!                        !Store all p,q for the given k,l
!
!                        rewind(unit_g_read)
!                        read(unit_g_read) maxlength
!                        do while(.true.)
!                            !Read all integrals for this k,l pair
!                            read(unit_g_read,iostat=ios) length,Indices(1:4,1:length),Buf(1:length)
!                            if(ios.ne.0) then
!                                call stop_all(t_r,"Error reading integrals")
!                            endif
!
!                            if(length.le.0) exit
!                            do ibuf = 1,length
!                                i_Ind=int(Indices(1,ibuf),i4)
!                                j_Ind=int(Indices(3,ibuf),i4)
!                                k_Ind=int(Indices(2,ibuf),i4)
!                                l_Ind=int(Indices(4,ibuf),i4)
!
!                                if(((i_Ind.eq.k).and.(k_Ind.eq.l)).or.((k_Ind.eq.k).and.(i_Ind.eq.l))) then
!                                    !We have a <kx|lx> or <lx|kx> integral
!                                    if(max(j_Ind,l_Ind).le.tnmo) then
!                                        !We want to store this integral.
!                                        pq_pair = FindSqPairInd(j_Ind,l_Ind,tnmo)
!                                        qp_pair = FindSqPairInd(l_Ind,j_Ind,tnmo)
!                                        r_in(pq_pair) = Buf(ibuf)
!                                        r_in(qp_pair) = Buf(ibuf)
!                                    endif
!                                endif
!                            enddo
!                        enddo   !End do looping over buffers
!
!                        !Now we have all pq for the given <k p | l q> integrals.
!                        !Sum in the correct coulomb contribution
!                        do p=1,tnmo
!                            do q=1,tnmo
!
!                                pq_pair = FindSqPairInd(p,q,tnmo)
!                                GenFockCABS(k,l) = GenFockCABS(k,l) + 2.0_dp*(r_in(pq_pair)*OneRDM(q,p))
!
!                            enddo
!                        enddo
!
!                    endif
!
!                enddo
!            enddo

            close(unit_g_stored)
            close(unit_g_read)
            if(tDaltonFormat) then
                deallocate(Indices,Buf)
            endif
            deallocate(r_in,TempGInts)
         
        endif

        !Construct the Generalised Fock matrix in complete:orbital space from others.
        do i=1,tntmo
            do j=1,tnmo
                if(i.le.tnmo) then
                    GenFockCompOrb(i,j) = GenFockOrb(i,j)
                else
                    GenFockCompOrb(i,j) = GenFockOrbCABS(j,i)
                endif
            enddo
        enddo


        !Now do tests on the generalised fock matrix for various forms of brillouins condition
        tCanonicalHF = .true.
        tCanonicalVirtHF = .true.
        MaxOccVirt = 0.0_dp
        MaxOccCABSFock = 0.0_dp
        MaxVirtCABSFock = 0.0_dp

        do i=1,tnmo
            do j=1,tnmo
                !Run over orbital:orbital block of the GFM
                if(IsOrbOcc(i).and.IsOrbOcc(j)) then
                    if(i.ne.j) then
                        !Are there non-zero off-diagonal elements in occ:occ block
                        if(abs(GenFockOrb(i,j)).gt.1.D-8) then
                            tCanonicalHF=.false.
                        endif
                    endif
                elseif((IsOrbOcc(i).and.(.not.IsOrbOcc(j))).or.(IsOrbOcc(j).and.(.not.IsOrbOcc(i)))) then
                    !In occ:virt block
                    !Brillouins condition - should be true if at HF soln
                    if(abs(GenFockOrb(i,j)).gt.MaxOccVirt) MaxOccVirt=abs(GenFockOrb(i,j))
                elseif((.not.IsOrbOcc(i)).and.(.not.IsOrbOcc(j))) then
                    if(i.ne.j) then
                        !In virt:virt block
                        if(abs(GenFockOrb(i,j)).gt.1.D-8) then
                            if(tHFRDMs) write(6,"(A,2I5,F20.10)") "Virt-Virt Fock element: ",i,j,GenFockOrb(i,j)
                            tCanonicalVirtHF=.false.
                        endif
                    endif
                else
                    !Should have covered all blocks in orb:orb space!
                    call stop_all(t_r,"Should not have got here")
                endif
            enddo
        enddo
        !Now run over orbital:CABS block
        do i=1,tnmo
            do j=tnmo+1,tntmo
                if(IsOrbOcc(i)) then
                    !In Occ:CABS block - will all be zero if GBC holds
                    if(abs(GenFockOrbCABS(i,j)).gt.MaxOccCABSFock) MaxOccCABSFock = abs(GenFockOrbCABS(i,j))
                else
                    !In virt:CABS block - will all be zero if EBC holds
                    if(abs(GenFockOrbCABS(i,j)).gt.MaxVirtCABSFock) MaxVirtCABSFock = abs(GenFockOrbCABS(i,j))
                endif
            enddo
        enddo
                
        !Here, we do some tests on the generalised fock matrix we have just constructed.
        !If we are doing a UEG system, then both the generalised, and original fock matrices should
        !both be strictly diagonal (since we are in a NO and HF basis). Check for this too.
        if(.not.tCanonicalHF) then
            write(6,"(A)") "!Warning! Occupied-Occupied off-diagonal generalised fock matrix elements found!"
            if(tUEG) call stop_all(t_r,"Occupied OBS not diagonal in generalised fock matrix.")
            if(tCanonicalVirtHF) then
                write(6,"(A)") "However, virtual-virtual generalised fock block *does* seem to be diagonal"
            endif
        endif
        if(.not.tCanonicalVirtHF) then
            write(6,"(A)") "!Warning! Virtual-Virtual off-diagonal generalised fock matrix elements found!"
            if(tUEG) call stop_all(t_r,"Virtual OBS not diagonal in generalised fock matrix.")
            if(tCanonicalHF) then
                write(6,"(A)") "However, occupied-occupied generalised fock block *does* seem to be diagonal"
            endif
        endif
        if((.not.tCanonicalHF).or.(.not.tCanonicalVirtHF)) then
            write(6,"(A)") "Generalised Orbital basis *not* diagonal."
        else
            write(6,"(A)") "Generalised Fock matrix diagonal in OBS basis."
        endif
        if(MaxOccVirt.gt.1.D-8) then
            write(6,"(A)") "Brillouin condition does not hold for generalised fock matrix."
            write(6,"(A,G20.10)") "Largest occupied-virtual generalised fock matrix element is: ",MaxOccVirt
            if(tUEG) call stop_all(t_r,"Non-zero off-diagonal matrix elements found in the GFM between occ and virt in OBS")
        else
            write(6,"(A)") "Brillouin's condition holds in orbital space for generalised fock matrix."
        endif
        if(MaxOccCABSFock.gt.1.D-8) then
            write(6,"(A)") "Generalised Brillouin condition does not hold for generalised fock matrix. Occupied space not closed "  &
            &   //"under generalised fock operator."
            write(6,"(A,G20.10)") "Largest Occupied-CABS generalised fock matrix element is: ",MaxOccCABSFock
            if(tUEG) call stop_all(t_r,"GBC not fulfilled in the generalised fock matrix (occ/CABS coupling)")
        else
            write(6,"(A)") "Generalised Brillouin condition holds. Occupied space closed under fock operator"
        endif
        if(MaxVirtCABSFock.gt.1.D-8) then
            write(6,"(A)") "Extended Brillouin condition does not hold. Virtual space not closed under "  &
            &   //"generalised fock operators."
            write(6,"(A,G20.10)") "Largest Virtual-CABS generalised fock matrix element is: ",MaxVirtCABSFock
            if(tUEG) call stop_all(t_r,"EBC not fulfilled in the generalised fock matrix (virt/CABS coupling)")
        else
            write(6,"(A)") "Extended Brillouin condition holds. Virtual space closed under fock operator"
        endif

    end subroutine CalcGenFock

!Calculate the one-electron hamiltonian matrix elements (over the whole space)
    subroutine CalcOneHamInts()
        integer :: ierr,unit_g_stored,unit_g_read,ios,ibuf,i_Ind,j_Ind,k_Ind,l_Ind
        integer(i8) :: MaxLength,length
        integer(i2) , allocatable :: Indices(:,:)
        integer :: n,m,i,OrbInd,in_pair,im_pair,mi_pair,mn_pair,ii_pair,loop
        integer :: n_Ind,m_Ind,j
        real(dp), allocatable :: TempGInts(:,:,:),Buf(:),r_in(:),TempGInts2(:,:)
        real(dp) :: HElem
        character(len=*), parameter :: t_r='CalcOneHamInts'
        integer :: nOrb,nAux,nElec,Ms2,nProp(3),PropBitLen,iSym
        real(dp) :: Z,GAM
        integer(i8) :: OrbSym(1000)
        logical :: UHF,UEG
        namelist /FCI/ nOrb,nAux,nElec,Ms2,OrbSym,iSym,UHF,UEG,GAM,PropBitLen,nProp

        if(tDaltonFormat) then
            !In FCIDUMP-type integral files, the tmat integrals are written out and already read in, therefore
            !we will only use this routine to check the values for self-consistency.
            allocate(tmat(tntmo,tntmo),stat=ierr)
            tmat(:,:)=0.D0
        endif

        unit_g_stored=get_free_unit()
        open(unit_g_stored,file='G_Dir_ccoo',status='old',form='unformatted',access='direct',action='read',recl=reclen_nmo_sq)

!Big pain: For the CABS:CABS fock matrix elements, we need <in|im> coulomb integrals, where i is occ, and n & m CABS - not stored!
!Instead, temporarily store a 3D matrix for this. TempGInts(i,n,m) = <in|im> integrals where n & m are both in CABS, and i is occ
!Incidently, for the UEG, *all* these integrals should be zero, since the momentum transfer vector = 0
!We do not need to worry about complex orbitals for these integrals, because due to the repeated index, there is only a four-fold
!permutational symmetry anyway...
!However, these are no longer zero, but should only approach zero in thermodynamic limit.
        allocate(TempGInts(1:tnocc,tnmo+1:tntmo,tnmo+1:tntmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error allocating large three-index temp array")
        TempGInts=0.0_dp

        if(tUEG) then
            !Hack to store <in|ni> integrals for UEG one-electron integrals
            !(will need to do something similar for general complex orbital systems in the future)
            allocate(TempGInts2(1:tnocc,tnmo+1:tntmo))
            TempGInts2(:,:)=0.0_dp
        endif

        unit_g_read=get_free_unit()
        if(tDaltonFormat) then
            open(unit_g_read,file='MO_G',status='old',form='unformatted',access='sequential',action='read')
        else
            if(tReadBin) then
                open(unit_g_read,file='GDUMPBIN',status='old',form='unformatted',action='read')
            else
                open(unit_g_read,file='GDUMP',status='old',form='formatted',action='read')
            endif
        endif
        rewind(unit_g_read)
        if(tDaltonFormat) then
            read(unit_g_read) maxlength
            allocate(Indices(4,MaxLength),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err 1")
            allocate(Buf(MaxLength),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err 1a")
        else
            if(.not.tReadBin) read(unit_g_read,fci)
            length = 1
        endif

        do while(.true.)
            if(tDaltonFormat) then
                read(unit_g_read,iostat=ios) length,Indices(1:4,1:length),Buf(1:length)
            else
                if(tReadBin) then
                    read(unit_g_read,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                else
                    read(unit_g_read,*,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                endif
            endif
            if(ios.gt.0) then
                call stop_all(t_r,"Error reading integrals")
            endif
            if((length.le.0).or.(ios.lt.0)) exit    !Finished reading all integrals.

            do ibuf=1,length
                !Indices are read in in chemical notation. We store in PHYSICAL
                !Find pair index - alternatively, we could back construct the ij pair, and compare to this
                if(tDaltonFormat) then
                    i_Ind=int(Indices(1,ibuf))
                    j_Ind=int(Indices(3,ibuf))
                    k_Ind=int(Indices(2,ibuf))
                    l_Ind=int(Indices(4,ibuf))
                    Z = Buf(ibuf)
                else
                    if((i_Ind.eq.0).or.(j_Ind.eq.0).or.(k_Ind.eq.0).or.(l_Ind.eq.0)) cycle
                endif

                !Permutational symmetry of 4traf integrals!
                if(i_Ind.le.tnmo) then
                    if(IsOrbOcc(i_Ind).and.(i_Ind.eq.k_Ind).and.(j_Ind.gt.tnmo).and.(l_Ind.gt.tnmo)) then
                        !<in|im> / <im|in>
                        TempGInts(OccOrbNum(i_Ind),j_Ind,l_Ind) = Z 
                        TempGInts(OccOrbNum(i_Ind),l_Ind,j_Ind) = Z
                    elseif(tUEG.and.IsOrbOcc(i_Ind).and.(i_Ind.eq.l_Ind).and.(j_Ind.gt.tnmo).and.(j_Ind.eq.k_Ind)) then
                        !<in|ni> with i occupied, and n cabs
                        TempGInts2(OccOrbNum(i_Ind),j_Ind) = Z
                    endif
                elseif(j_Ind.le.tnmo) then
                    if(IsOrbOcc(j_Ind).and.(j_Ind.eq.l_Ind).and.(i_Ind.gt.tnmo).and.(k_Ind.gt.tnmo)) then
                        !<ni|mi> / <mi|ni>
                        TempGInts(OccOrbNum(j_Ind),i_Ind,k_Ind) = Z 
                        TempGInts(OccOrbNum(j_Ind),k_Ind,i_Ind) = Z
                    elseif(tUEG.and.IsOrbOcc(j_Ind).and.(j_Ind.eq.k_Ind).and.(i_Ind.gt.tnmo).and.(i_Ind.eq.l_Ind)) then
                        !<ni|in> with i occupied and n cabs
                        TempGInts2(OccOrbNum(j_Ind),i_Ind) = Z
                    endif
                endif

            enddo

        enddo

        close(unit_g_read)
        if(tDaltonFormat) then
            deallocate(Indices,Buf)
        endif

!       H_nm = f_nm - \sum_{i \in HF} [ <in|im> - <in|mi> ]
        allocate(r_in(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"allocation err")
        do n=1,tntmo
            do m=n,tntmo

                if((n.le.tnmo).and.(m.le.tnmo)) then
                    HElem=fockOrb(n,m)
                elseif(((n.le.tnmo).and.(m.gt.tnmo)).or.((n.gt.tnmo).and.(m.le.tnmo))) then
                    !orb:CABS block
                    HElem=FockOrbCABS(min(n,m),max(n,m))
                else
                    HElem=FockCABS(n,m)
                endif

                if((n.le.tnmo).and.(m.le.tnmo)) then
                    !orbital, orbital space

                    r_in(:)=0.0_dp
                    do i=1,tnocc
                        OrbInd = OccOrbs(i)

                        in_pair = FindSqPairInd(OrbInd,n,tntmo)  !Complete space orbital pairs are stored
                        read(unit_g_stored,rec=in_pair) (r_in(loop),loop=1,tnmo_sq)

                        im_pair = FindSqPairInd(OrbInd,m,tnmo)
                        mi_pair = FindSqPairInd(m,OrbInd,tnmo)

                        !Hold on a sec - spin symmetry?
                        HElem = HElem - (2.0_dp*r_in(im_pair) - r_in(mi_pair))
                    enddo

                elseif((min(n,m)).le.tnmo) then
                    !Occupied:CABS

                    !Ensure that m < n (switch if necessary)
                    if(m.gt.n) then
                        !switch
                        n_Ind=m
                        m_Ind=n
                    else
                        n_Ind=n
                        m_Ind=m
                    endif
                    r_in(:)=0.0_dp
                    do i=1,tnocc
                        !Coulomb: <in|im>   (m occ)
                        OrbInd = OccOrbs(i)
                        in_pair = FindSqPairInd(OrbInd,n_Ind,tntmo)
                        read(unit_g_stored,rec=in_pair) (r_in(loop),loop=1,tnmo_sq)
                        im_pair = FindSqPairInd(OrbInd,m_Ind,tnmo)
                        HElem = HElem - (2.0_dp*r_in(im_pair))

                        !Exchange: <in|mi> (m occ)
                        mi_pair = FindSqPairInd(m_Ind,OrbInd,tnmo)
                        HElem = HElem + r_in(mi_pair)
                    enddo


!                    !Use permutational symmetry of integrals
!                    !Ensure that n is in orbital space, and then calculate <im|in> - <mn|ii>
!                    !WARNING: This permutational symmetry is not available with complex orbitals.
!                    if(n.gt.m) then
!                        n_Ind=m
!                        m_Ind=n
!                    else
!                        n_Ind=n
!                        m_Ind=m
!                    endif
!                    r_in(:)=0.0_dp
!                    do i=1,tnocc
!                        !Coulomb part - <im|in> with m in CABS and n in orbital
!                        OrbInd = OccOrbs(i)
!                        im_pair = FindSqPairInd(OrbInd,m_Ind,tntmo)
!                        read(unit_g_stored,rec=im_pair) (r_in(loop),loop=1,tnmo_sq)
!                        in_pair = FindSqPairInd(OrbInd,n_Ind,tnmo)
!                        HElem = HElem - (2.0_dp*r_in(in_pair))
!
!                        !Exchange part - <mn|ii> with m in CABS and n in orbital
!                        mn_pair = FindSqPairInd(m_Ind,n_Ind,tntmo)
!                        read(unit_g_stored,rec=mn_pair) (r_in(loop),loop=1,tnmo_sq)
!                        ii_pair = FindSqPairInd(OrbInd,OrbInd,tnmo)
!                        HElem = HElem + r_in(ii_pair)
!                    enddo

                else
                    !Cabs Cabs space
                    if(.not.tComplexOrbs) then
                        do i=1,tnocc
                            OrbInd = OccOrbs(i)
                            !Already explicitly stored the <im|in> integral in temp
                            HElem = HElem - (2.0_dp*TempGInts(i,m,n))

                            !Exchange part read from disk - <mn|ii>
                            mn_pair = FindSqPairInd(m,n,tntmo)
                            read(unit_g_stored,rec=mn_pair) (r_in(loop),loop=1,tnmo_sq)
                            ii_pair = FindSqPairInd(OrbInd,OrbInd,tnmo)
                            HElem = HElem + r_in(ii_pair)
                        enddo
                    else
                        !This is a lot more complicated to do. 
                        !We do not store <in|im> OR <in|mi> for complex orbitals. BOO
                        if(.not.tUEG) then
                            call stop_all(t_r,"Ugh - calculating one-electron integrals with general complex orbitals not implemented")
                        endif
                        if(m.eq.n) then
                            do i=1,tnocc
                                HElem = HElem - (2.0_dp*TempGInts(i,m,n))

                                !Exchange part in TempGInts2
                                HElem = HElem + TempGInts2(i,n)
                            enddo
                        endif
                    endif

                endif


                !Store the one-electron elements for calculating the HF energy later
                if(tDaltonFormat) then
                    tmat(n,m)=HElem
                    tmat(m,n)=HElem
                else
                    !We have already read in the tmat integrals directly from the GDUMP file.
                    !We are only calculating them here for self-consistency.
                    if(abs(tmat(n,m)-HElem).gt.1.D-7) then
                        write(6,*) "One electron integral error..."
                        write(6,*) n,m,tmat(n,m),HElem
!                        call stop_all(t_r,"One electron integral self-consistency check fail...")
                    elseif(abs(tmat(m,n)-HElem).gt.1.D-7) then
                        write(6,*) "One electron integral error..."
                        write(6,*) m,n,tmat(m,n),HElem
!                        call stop_all(t_r,"One electron integral self-consistency check fail 2...")
                    endif
                endif

            enddo
        enddo
        if(tUEG) deallocate(TempGInts2)
        deallocate(r_in,TempGInts)
        close(unit_g_stored)

        !write out one-electron integrals. 
        !This is more for ease of debugging more than anything else.
        open(unit_g_stored,file='TMAT',status='unknown')
        do i=1,tnmo
            do j=i,tnmo
                if(abs(tmat(i,j)).gt.1.D-8) then
                    write(unit_g_stored,*) i,j,tmat(i,j)
                endif
            enddo
        enddo
        close(unit_g_stored)

    end subroutine CalcOneHamInts


!This routine will read in the OneRDM from NECI or other...
!The OneRDM is stored explicitly, while the TwoRDM is read and written straight out to a direct access file
!1-Density matrix stored in memory in spin orbitals, but only alpha, alpha block - the beta beta block is assumed
!to be equivalent
    subroutine ReadOneRDM()
        implicit none
        character(*), parameter :: t_r='ReadOneRDM'
        integer :: ierr,i,OneRDM_read,j,ios,i_Ind,j_Ind,TwoRDM_read,pp,p,q,r,s
        integer :: pInd,qInd,rInd,sInd,ispin,Antisym_spin,pqPair,qpPair,iSpinNum
        real(dp) :: OneRDMVal,TwoRDMVal,temp
        logical :: exists
        real(dp) , allocatable :: TwoRDM_aaaa(:,:),TwoRDM_abab(:,:),TwoRDM_abba(:,:)

        allocate(OneRDM(tnmo,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
        OneRDM(:,:) = 0.0_dp

        if(tHFRDMs) then

            write(6,"(A)") "Using Hartree--Fock density matrix in calculation..."
            write(6,*) ""

            !Store OneRDM
            do i=1,tnmo
                if(IsOrbOcc(i)) OneRDM(i,i) = 1.0_dp
            enddo
        elseif(tIntegrateTwoRDM) then
            
            write(6,"(A)") "Integrating 2-RDM to calculate 1-RDM..."
            write(6,*) ""
            
            !Fill diagonal frozen orbital entries with the HF matrix.
            !Shouldn't need to do this, as we are using the 2RDM *with* frozen orbitals included.
            !However, it would be quicker to do it like this in the future
!            do i=1,tnfrz
!                OneRDM(FrzOrbs(i),FrzOrbs(i)) = 1.0_dp
!            enddo
            
            !Two RDM is more difficult - don't want to store whole thing at once.
            !Store all {rs},q in one go.
            allocate(TwoRDM_aaaa(tnmo_sq,tnmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
            allocate(TwoRDM_abab(tnmo_sq,tnmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
            allocate(TwoRDM_abba(tnmo_sq,tnmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Error in allocation")

            TwoRDM_read=get_free_unit()
            if(tSpatRDMs) then
                !Spatial orbitals. Start with reading in the aaaa file.
                open(TwoRDM_read,file='TwoRDM_aaaa',status='old',action='read')
                iSpinNum=1
            else
                open(TwoRDM_read,file='TwoRDM',status='old',action='read')
            endif
            
            do pp=1,tnmo
                TwoRDM_aaaa(:,:)=0.0_dp
                TwoRDM_abab(:,:)=0.0_dp
                TwoRDM_abba(:,:)=0.0_dp
                if(.not.IsOrbFrz(pp)) then
                    rewind(TwoRDM_read)
                    do while(.true.)
!                        read(TwoRDM_read,"(4I6,G25.17)",iostat=ios) p,q,r,s,TwoRDMVal
                        read(TwoRDM_read,*,iostat=ios) p,q,r,s,TwoRDMVal
                        if(ios.gt.0) call stop_all(t_r,"Error reading TwoRDM")
                        if(ios.lt.0) then
                            !EOF. What now.
                            if(.not.tSpatRDMs) then
                                !Spin orbital - should all be in one file
                                exit
                            else
                                !Spatial orbitals
                                !Close current spin file and open next one
                                call OpenNext2RDMFile(iSpinNum,TwoRDM_read)
                                if(iSpinNum.eq.1) then
                                    !Finished reading in all the files, and have reopened aaaa file for next read.
                                    !Exit.
                                    exit
                                else
                                    !We have opened the next file. Start reading again.
                                    cycle
                                endif
                            endif
                        endif

                        pInd = p
                        qInd = q
                        rInd = r
                        sInd = s
                        if(.not.tSpatRDMs) then
                            !Reading in spin-orbital.
                            !Work out the spin type of the element.
                            ispin = WhichSpin(pInd,qInd,rInd,sInd)  !ispin now tells us the spin of the element
                        else
                            !In spatial orbitals.
                            !The spin type of the element is then simply determined by which file we are reading in from.
                            ispin = iSpinNum
                        endif
                        !Antisym_spin is the spin of the antisymmetrised element of the 2RDM
                        if(ispin.eq.1) then
                            Antisym_spin = 1
                        elseif(ispin.eq.2) then
                            Antisym_spin = 3
                        elseif(ispin.eq.3) then
                            Antisym_spin = 2
                        elseif(ispin.eq.0) then
                            if(abs(TwoRDMVal).gt.1.D-8) then
                                write(6,*) p,q,r,s,TwoRDMVal
                                call stop_all(t_r,"Forbidden spin type in 2RDM values")
                            endif
                            cycle
                        endif
                        if(.not.tSpatRDMs) then
                            !Spin orbitals - convert to spatial
                            call Conv2Spat(pInd,qInd,rInd,sInd)
                        endif
                        !Convert to non-frozen basis indicies
                        pInd = OrbMapping(pInd)
                        qInd = OrbMapping(qInd)
                        rInd = OrbMapping(rInd)
                        sInd = OrbMapping(sInd)
                        if(pInd.eq.pp) then
                            !Store as p,q,r,s
                            call Fill2RDMArr(ispin,pInd,qInd,rInd,sInd,TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                            !Also store as p,q,s,r - remembering sign change, and potential spin change!
                            call Fill2RDMArr(Antisym_spin,pInd,qInd,sInd,rInd,-TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                        endif
                        if(qInd.eq.pp) then
                            !Store as q,p,s,r
                            call Fill2RDMArr(ispin,qInd,pInd,sInd,rInd,TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                            !Also store as q,p,r,s - remembering sign change, and potential spin change!
                            call Fill2RDMArr(Antisym_spin,qInd,pInd,rInd,sInd,-TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                        endif
                        if(rInd.eq.pp) then
                            !Store as r,s,p,q (Hermiticity)
                            call Fill2RDMArr(ispin,rInd,sInd,pInd,qInd,TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                            !Also store as r,s,q,p - remembering sign change, and potential spin change!
                            call Fill2RDMArr(Antisym_spin,rInd,sInd,qInd,pInd,-TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                        endif
                        if(sInd.eq.pp) then
                            !Store as s,r,q,p
                            call Fill2RDMArr(ispin,sInd,rInd,qInd,pInd,TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                            !Also store as s,r,p,q - remembering sign change, and potential spin change!
                            call Fill2RDMArr(Antisym_spin,sInd,rInd,pInd,qInd,-TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                        endif

                    enddo
                endif

                !Now integrate down the two-RDM
                !We don't actually need the abba matrix, and could write stuff out here to save time, and read in later
                !Even better would be to do this in NECI though...!
                do r=1,tnmo
                    if(IsOrbFrz(r)) cycle

!                    temp = 0.0_dp
                    do q=1,tnmo
                        if(IsOrbFrz(q)) cycle
                    
                        OneRDM(pp,r) = OneRDM(pp,r) + TwoRDM_aaaa(FindSqPairInd(r,q,tnmo),q)
                        OneRDM(pp,r) = OneRDM(pp,r) + TwoRDM_abab(FindSqPairInd(r,q,tnmo),q)
!                        temp = temp + TwoRDM_abba(FindSqPairInd(r,q,tnmo),q)

                    enddo

                    OneRDM(pp,r) = OneRDM(pp,r) / real((NEl-tnfrz*2)-1,dp)
!                    temp = temp / real(NEl-1,dp)
!                    write(6,*) r,pp,temp
                    if(abs(OneRDM(pp,r)).gt.1.D-8) then
                        write(47,"(2I6,G25.17)",iostat=ios) r,pp,OneRDM(pp,r)
                    endif
                enddo

            enddo

            close(TwoRDM_read)
            deallocate(TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)

            !Check that the frozen orbitals have come out correctly.
            do r=1,tnfrz
                OneRDM(FrzOrbs(r),FrzOrbs(r)) = 1.0_dp
            enddo

        else
            write(6,*) "Reading in arbitrary density matrices for calculation..."
            write(6,*) ""

            inquire(file='OneRDM',exist=exists)
            if(.not.exists) then
                call stop_all(t_r,"Cannot find file 'OneRDM' to read in One-body density matrix")
            endif

            !Fill diagonal frozen orbital entries with the HF matrix.
            do i=1,tnfrz
                OneRDM(FrzOrbs(i),FrzOrbs(i)) = 1.0_dp
            enddo

            OneRDM_read=get_free_unit()
            open(OneRDM_read,file='OneRDM',status='old',action='read')
            rewind(OneRDM_read)

            do while(.true.)
!                read(OneRDM_read,"(2I6,G25.17)",iostat=ios) i,j,OneRDMVal
                read(OneRDM_read,*,iostat=ios) i,j,OneRDMVal
                if(ios.gt.0) call stop_all(t_r,"Error reading OneRDM")
                if(ios.lt.0) exit
                !Read in in spin-orbital notation.
                if(.not.tSpatRDMs) then
                    !Convert to spatial orbitals to store
                    if(mod(i,2).eq.0) then
                        i_Ind = i/2
                    else
                        i_Ind = (i+1)/2
                    endif
                    if(mod(j,2).eq.0) then
                        j_Ind = j/2
                        if((mod(i,2).eq.1).and.(abs(OneRDMVal).gt.1.D-8)) then
                            call stop_all(t_r,"Spin not integrated correctly in OneRDM")
                        endif
                    else
                        j_Ind = (j+1)/2
                        if((mod(i,2).eq.0).and.(abs(OneRDMVal).gt.1.D-8)) then
                            call stop_all(t_r,"Spin not integrated correctly in OneRDM")
                        endif
                    endif
                else
                    !Orbitals read in already in spatial format
                    i_Ind = i
                    j_Ind = j
                endif
                !Now convert to the non-frozen basis
                i_Ind = OrbMapping(i_Ind) 
                j_Ind = OrbMapping(j_Ind)
                if((abs(OneRDM(i_Ind,j_Ind)-OneRDMVal).gt.1.D-4).and.(abs(OneRDM(i_Ind,j_Ind)).gt.1.D-8)) then
                    call stop_all(t_r,"Error in spin symmetry in read in RDM")
                endif
                OneRDM(i_Ind,j_Ind)=OneRDMVal
                OneRDM(j_Ind,i_Ind)=OneRDMVal
            enddo
            close(OneRDM_read)
        endif

!        do i=1,tnmo
!            do j=1,tnmo
!                write(6,"(f10.5)",advance='no') OneRDM(i,j)
!            enddo
!            write(6,*)
!        enddo

    end subroutine ReadOneRDM

    !Does what it says on the tin. Only for spatial orbitals, and iSpinNum should be 1 or 2, indicating
    !whether we have just finished reading in aaaa or abab files.
    subroutine OpenNext2RDMFile(iSpinNum,TwoRDM_read)
        implicit none
        integer, intent(inout) :: iSpinNum  !Indicating how far through the different spin files we are
        integer, intent(in) :: TwoRDM_read  !Unit number of current RDM file
        logical :: exists,exists_2
        character(*), parameter :: t_r='OpenNext2RDMFile'

        if(iSpinNum.eq.1) then
            !Just finished reading in aaaa file.
            !Now read abab file
            inquire(file='TwoRDM_abab',exist=exists)
            if(exists) then
                open(TwoRDM_read,file='TwoRDM_abab',status='old',action='read')
                rewind(TwoRDM_read)
                iSpinNum=2
            else
                !Seems to be no abab file. Could be ok if there is an abba file.
                inquire(file='TwoRDM_abba',exist=exists_2)
                if(exists_2) then
                    open(TwoRDM_read,file='TwoRDM_abba',status='old',action='read')
                    rewind(TwoRDM_read)
                    iSpinNum=3
                else
                    call stop_all(t_r,"Cannot seem to find the mixed spin spatial RDM files...")
                endif
            endif
        elseif(iSpinNum.eq.2) then
            !Just finished reading in abab file.
            !Now try to read in abba file.
            inquire(file='TwoRDM_abba',exist=exists)
            if(exists) then
                open(TwoRDM_read,file='TwoRDM_abba',status='old',action='read')
                rewind(TwoRDM_read)
                iSpinNum=3
            else
                !We've finished reading.
                !Reopen the aaaa file for the next loop
                inquire(file='TwoRDM_aaaa',exist=exists_2)
                if(exists_2) then
                    open(TwoRDM_read,file='TwoRDM_aaaa',status='old',action='read')
                    rewind(TwoRDM_read)
                    iSpinNum=1  !Indicate that we've finished reading everything
                else
                    call stop_all(t_r,"Cannot seem to find the same spin spatial RDM files...")
                endif
            endif
        elseif(iSpinNum.eq.3) then
            !Finished reading all the files. Return iSpinNum = 1 to indicate this, and reopen aaaa file.
            inquire(file='TwoRDM_aaaa',exist=exists)
            if(exists) then
                open(TwoRDM_read,file='TwoRDM_aaaa',status='old',action='read')
                rewind(TwoRDM_read)
                iSpinNum=1  !Indicate that we've finished reading everything
            else
                call stop_all(t_r,"Cannot seem to find the same spin spatial RDM files...")
            endif
        else
            call stop_all(t_r,"Should not get here with different iSpinNum")
        endif

    end subroutine OpenNext2RDMFile

!This routine will read in the TwoRDM from NECI or other...
!The TwoRDM is read and written straight out to a direct access file
!TwoRDM is stored in direct file, TwoRDM_Dir as (TwoRDM_ij)^kl  where each ij is another record.
!2-Density matrices stored on disk as spin orbitals!!
!We also calculate Lambda_Bar here
    subroutine ReadTwoRDM()
        implicit none
        integer :: ierr,ios,p,q,r,s,TwoRDM_read,counter,qpPair,j,x,i,xj_pair,g_ints_read
        integer :: unit_2RDM_write_aaaa,unit_2RDM_write_abab,unit_2RDM_write_abba,y,k,yk_pair
        integer :: loop,ispin,Antisym_spin,pInd,pp,pqPair,qInd,rInd,sInd,rspair,u,uInd,iSpinNum
        integer :: unit_Lambda_Tld_write_aaaa,unit_Lambda_Tld_write_abab,unit_Lambda_Tld_write_abba
        integer :: rs_pair,unit_Lambda_Tld_read_aaaa,unit_Lambda_Tld_read_abab,unit_Lambda_Tld_read_abba
        character(len=*), parameter :: t_r='ReadTwoRDM'
        real(dp) , allocatable :: TwoRDM_aaaa(:,:),TwoRDM_abab(:,:),TwoRDM_abba(:,:)
        real(dp) , allocatable :: TempArr_aaaa(:,:),TempArr_abab(:,:),TempArr_abba(:,:),OrbOrbMat(:,:)
        real(dp) , allocatable :: TempArr_aaaa_2(:),TempArr_abba_2(:),TempArr_abab_2(:)
        real(dp) , allocatable :: GammaF(:,:),Lambda_Tld_aaaa(:,:),Lambda_Tld_abab(:,:),Lambda_Tld_abba(:,:)
        real(dp) :: TwoRDMVal,fac,DDOT
        logical :: exists

        !In this routine, we construct Lambda_Bar
        !Here, we calculate LB_p^s = -f_r^q L_pq^rs
        !where L is the 2-cumulant
        allocate(Lambda_Bar(tnmo,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
        Lambda_Bar(:,:) = 0.0_dp

        !In this routine, we may also construct part of the X matrix for the 
        !MR-orbital relaxation too. In here, we calculate RelaxX_y^x = f_j^i 2RDM_yi^xj if non-dyall,
        !or RelaxX_y^x = G_yk^ij 2RDM_ij^xk if dyall.
        !Later, this will be combined with the rest of the cumulant.
        if(.not.tHFRDMs) then
            if(tDyall.or.tZerothRelaxBoth) then
                allocate(RelaxXDyall(tnmo,tnmo),stat=ierr)
                if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
                RelaxXDyall(:,:) = 0.0_dp
                !Need the g integrals to contract the 2RDM with.
                g_ints_read = get_free_unit()
                open(g_ints_read,file='G_Dir_ccoo',status='old',form='unformatted',   &
                            access='direct',action='read',recl=reclen_nmo_sq)
            endif
            if((.not.tDyall).or.tZerothRelaxBoth) then
                allocate(RelaxXFock(tnmo,tnmo),stat=ierr)
                if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
                RelaxXFock(:,:) = 0.0_dp
            endif
        endif

        !GammaF_1: U_p^u = 1/2 G_p^t f_t^u
        !This is used later in the construction of Lambda_Tld
        allocate(GammaF(tnmo,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Err")
        call DGEMM('N','N',tnmo,tnmo,tnmo,0.5_dp,OneRDM,tnmo,GenFockOrb,tnmo,0.0_dp,GammaF,tnmo)
        !U_p^u (p fast).
                
        !Contract U with 1RDM over u
        !This is again used later (in the OneRDM part) in Lambda_Tld
        allocate(OrbOrbMat(tnmo,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc Err")
        call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,GammaF,tnmo,OneRDM,tnmo,0.0_dp,OrbOrbMat,tnmo)
        !p fast again
        
        allocate(TempArr_aaaa_2(tnmo_sq),stat=ierr)
        allocate(TempArr_abab_2(tnmo_sq),stat=ierr)
        allocate(TempArr_abba_2(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc err")

        if(.not.tHFRDMs) then
            !Two RDM is more difficult - don't want to store whole thing at once.
            TwoRDM_read=get_free_unit()
            if(tSpatRDMs) then
                !Spatial orbitals. Start with reading in the aaaa file.
                open(TwoRDM_read,file='TwoRDM_aaaa',status='old',action='read')
                iSpinNum = 1
            else
                open(TwoRDM_read,file='TwoRDM',status='old',action='read')
            endif

            !The files to write the RDMs into...
            unit_2RDM_write_aaaa=get_free_unit()
            open(unit_2RDM_write_aaaa,file='TwoRDM_Dir_aaaa',status='unknown',form='unformatted',    &
                access='direct',action='write',recl=reclen_nmo_sq)
            unit_2RDM_write_abab=get_free_unit()
            open(unit_2RDM_write_abab,file='TwoRDM_Dir_abab',status='unknown',form='unformatted',    &
                access='direct',action='write',recl=reclen_nmo_sq)
            unit_2RDM_write_abba=get_free_unit()
            open(unit_2RDM_write_abba,file='TwoRDM_Dir_abba',status='unknown',form='unformatted',    &
                access='direct',action='write',recl=reclen_nmo_sq)
            
            !The following units are to write out intermediate quantities for the constuction of Lamda_Tld    
            unit_Lambda_Tld_write_aaaa=get_free_unit()
            open(unit_Lambda_Tld_write_aaaa,file='L_Tld_Tmp_aaaa',status='unknown',form='unformatted',    &
                access='direct',action='write',recl=reclen_nmo_sq)
            unit_Lambda_Tld_write_abab=get_free_unit()
            open(unit_Lambda_Tld_write_abab,file='L_Tld_Tmp_abab',status='unknown',form='unformatted',    &
                access='direct',action='write',recl=reclen_nmo_sq)
            unit_Lambda_Tld_write_abba=get_free_unit()
            open(unit_Lambda_Tld_write_abba,file='L_Tld_Tmp_abba',status='unknown',form='unformatted',    &
                access='direct',action='write',recl=reclen_nmo_sq)

            !Store all {rs},q in one go.
            allocate(TwoRDM_aaaa(tnmo_sq,tnmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
            allocate(TwoRDM_abab(tnmo_sq,tnmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
            allocate(TwoRDM_abba(tnmo_sq,tnmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
            allocate(Lambda_Tld_aaaa(tnmo_sq,tnmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
            allocate(Lambda_Tld_abab(tnmo_sq,tnmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
            allocate(Lambda_Tld_abba(tnmo_sq,tnmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
            allocate(TempArr_aaaa(tnmo,tnmo),stat=ierr)
            allocate(TempArr_abab(tnmo,tnmo),stat=ierr)
            allocate(TempArr_abba(tnmo,tnmo),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Error in allocation")

            write(6,'(F10.3,A)') 3.0_dp*(real(tnmo,dp)**4.0_dp)*8.0_dp*9.31322575D-10,   &
                " Gb required to write out TwoRDM_Dir files."
            call flush(6)

            do pp=1,tnmo
                TwoRDM_aaaa(:,:)=0.0_dp
                TwoRDM_abab(:,:)=0.0_dp
                TwoRDM_abba(:,:)=0.0_dp
                if(.not.IsOrbFrz(pp)) then
                    rewind(TwoRDM_read)
                    do while(.true.)
                        read(TwoRDM_read,*,iostat=ios) p,q,r,s,TwoRDMVal
!                        read(TwoRDM_read,"(4I6,G25.17)",iostat=ios) p,q,r,s,TwoRDMVal
                        if(ios.gt.0) call stop_all(t_r,"Error reading TwoRDM")
                        if(ios.lt.0) then
                            !EOF. What now?
                            if(.not.tSpatRDMs) then
                                !Spinorbitals - should all be in one file.
                                exit
                            else
                                !Spatial orbitals
                                !Close current spin file and open next one
                                call OpenNext2RDMFile(iSpinNum,TwoRDM_read)
                                if(iSpinNum.eq.1) then
                                    !Finished reading in all the files. Exit.
                                    exit
                                else
                                    cycle
                                endif
                            endif
                        endif

                        pInd = p
                        qInd = q
                        rInd = r
                        sInd = s
                        if(.not.tSpatRDMs) then
                            !Reading in spin-orbital
                            !Work out the spin type of the element
                            ispin = WhichSpin(pInd,qInd,rInd,sInd)  !ispin now tells us the spin of the element
                        else
                            !In spatial orbitals.
                            !The spin type of the element is then simply determined by which file we are reading in from.
                            ispin = iSpinNum
                        endif
                        !Antisym_spin is the spin of the antisymmetrised element of the 2RDM
                        if(ispin.eq.1) then
                            Antisym_spin = 1
                        elseif(ispin.eq.2) then
                            Antisym_spin = 3
                        elseif(ispin.eq.3) then
                            Antisym_spin = 2
                        elseif(ispin.eq.0) then
                            if(abs(TwoRDMVal).gt.1.D-8) then
                                write(6,*) p,q,r,s,TwoRDMVal
                                call stop_all(t_r,"Forbidden spin type in 2RDM values")
                            endif
                            cycle
                        endif
                        if(.not.tSpatRDMs) then
                            !Spin orbitals - convert to spatial
                            call Conv2Spat(pInd,qInd,rInd,sInd)
                        endif
                        !Convert to non-frozen basis indicies
                        pInd = OrbMapping(pInd)
                        qInd = OrbMapping(qInd)
                        rInd = OrbMapping(rInd)
                        sInd = OrbMapping(sInd)
                        if(pInd.eq.pp) then
                            !Store as p,q,r,s
                            call Fill2RDMArr(ispin,pInd,qInd,rInd,sInd,TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                            !Also store as p,q,s,r - remembering sign change, and potential spin change!
                            call Fill2RDMArr(Antisym_spin,pInd,qInd,sInd,rInd,-TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                        endif
                        if(qInd.eq.pp) then
                            !Store as q,p,s,r
                            call Fill2RDMArr(ispin,qInd,pInd,sInd,rInd,TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                            !Also store as q,p,r,s - remembering sign change, and potential spin change!
                            call Fill2RDMArr(Antisym_spin,qInd,pInd,rInd,sInd,-TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                        endif
                        if(rInd.eq.pp) then
                            !Store as r,s,p,q (Hermiticity)
                            call Fill2RDMArr(ispin,rInd,sInd,pInd,qInd,TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                            !Also store as r,s,q,p - remembering sign change, and potential spin change!
                            call Fill2RDMArr(Antisym_spin,rInd,sInd,qInd,pInd,-TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                        endif
                        if(sInd.eq.pp) then
                            !Store as s,r,q,p
                            call Fill2RDMArr(ispin,sInd,rInd,qInd,pInd,TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                            !Also store as s,r,p,q - remembering sign change, and potential spin change!
                            call Fill2RDMArr(Antisym_spin,sInd,rInd,pInd,qInd,-TwoRDMVal,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
                        endif

                    enddo
                endif

                do q=1,tnmo
                    
                    pqPair = FindSqPairInd(pp,q,tnmo)
                    qpPair = FindSqPairInd(q,pp,tnmo)

                    !Add elements for the neglected frozen orbitals
                    if(IsOrbFrz(pp).and.IsOrbFrz(q)) then
                        !If both pp and q orbitals are frozen, then only add diagonal elements.
                        if(pp.ne.q) then
                            if(TwoRDM_aaaa(pqPair,q).ne.0.0_dp) call stop_all(t_r,"Error in freezing 2RDM")
                            TwoRDM_aaaa(pqPair,q) = 1.0_dp

                            if(TwoRDM_aaaa(qpPair,q).ne.0.0_dp) call stop_all(t_r,"Error in freezing 2RDM")
                            TwoRDM_aaaa(qpPair,q) = -1.0_dp
                            
                            if(TwoRDM_abab(pqPair,q).ne.0.0_dp) call stop_all(t_r,"Error in freezing 2RDM")
                            TwoRDM_abab(pqPair,q) = 1.0_dp
                            
                            if(TwoRDM_abba(qpPair,q).ne.0.0_dp) call stop_all(t_r,"Error in freezing 2RDM")
                            TwoRDM_abba(qpPair,q) = -1.0_dp
                        else
                            if(TwoRDM_abab(pqPair,q).ne.0.0_dp) call stop_all(t_r,"Error in freezing 2RDM")
                            TwoRDM_abab(pqPair,q) = 1.0_dp  !All four indices same
                            
                            if(TwoRDM_abba(pqPair,q).ne.0.0_dp) call stop_all(t_r,"Error in freezing 2RDM")
                            TwoRDM_abba(pqPair,q) = -1.0_dp
                        endif
                    elseif(IsOrbFrz(pp)) then
                        !pp is frozen and q is not.
                        !Add the correct 1RDM element

                        do j=1,tnmo
                            if(IsOrbFrz(j)) cycle

                            TwoRDM_aaaa(FindSqPairInd(pp,j,tnmo),q) = OneRDM(q,j)
                            TwoRDM_aaaa(FindSqPairInd(j,pp,tnmo),q) = - OneRDM(q,j)

                            TwoRDM_abab(FindSqPairInd(pp,j,tnmo),q) = OneRDM(q,j)
                            TwoRDM_abba(FindSqPairInd(j,pp,tnmo),q) = - OneRDM(q,j)
                        enddo

                    elseif(IsOrbFrz(q)) then
                        !q is frozen and p is not
                        !Add the correct 1RDM element

                        do j=1,tnmo
                            if(IsOrbFrz(j)) cycle

                            TwoRDM_aaaa(FindSqPairInd(j,q,tnmo),q) = OneRDM(pp,j)
                            TwoRDM_aaaa(FindSqPairInd(q,j,tnmo),q) = - OneRDM(pp,j)

                            TwoRDM_abab(FindSqPairInd(j,q,tnmo),q) = OneRDM(pp,j)
                            TwoRDM_abba(FindSqPairInd(q,j,tnmo),q) = - OneRDM(pp,j)

                        enddo

                    endif

                    write(unit_2RDM_write_aaaa,rec=pqPair) (TwoRDM_aaaa(loop,q),loop=1,tnmo_sq) 
                    write(unit_2RDM_write_abab,rec=pqPair) (TwoRDM_abab(loop,q),loop=1,tnmo_sq) 
                    write(unit_2RDM_write_abba,rec=pqPair) (TwoRDM_abba(loop,q),loop=1,tnmo_sq) 
                enddo

                !Calculate Lambda_Bar_q^s = -f_r^p Lambda_qp^rs
                !We need then aaaa and abba parts
                do s=1,tnmo

                    TempArr_aaaa(:,:) = 0.0_dp
                    TempArr_abba(:,:) = 0.0_dp
                    do r=1,tnmo
                        do q=1,tnmo
                            rs_pair = FindSqPairInd(r,s,tnmo)
                            TempArr_aaaa(q,r) = TwoRDM_aaaa(rs_pair,q)
                            TempArr_abba(q,r) = TwoRDM_abba(rs_pair,q)
                        enddo
                    enddo

                    Lambda_Bar(pp,s) = -1.0_dp * DDOT(tnmo_sq,GenFockOrb,1,TempArr_aaaa,1)
                    Lambda_Bar(pp,s) = Lambda_Bar(pp,s) - DDOT(tnmo_sq,GenFockOrb,1,TempArr_abba,1)
                enddo

                if(.not.tHFRDMs) then
                    !Here we calculate the initial contraction over the two-rdm for the MR orbital relaxation.

                    if(tDyall.or.tZerothRelaxBoth) then
                        !Calculate G_yk^ij 2RDM_ij^xk = RelaxX_y^x (y fast)
                        !By permutational symmetry, this is the same as G_yk^ij 2RDM_xk^ij
                        !The loop over pp is looping over x
                        
                        do y=1,tnmo
                            !loop over y

                            do k=1,tnmo
                                !Contract over k, with contractions over i and j inside.

                                TempArr_aaaa_2(:) = 0.0_dp
                                TempArr_abab_2(:) = 0.0_dp
                                TempArr_abba_2(:) = 0.0_dp

                                !Get the antisymmetrised vector (g_yk)^ij for all spin types.
                                yk_pair = FindSqPairInd(y,k,tntbf)
                                read(g_ints_read,rec=yk_pair) (TempArr_abab_2(loop),loop=1,tnmo_sq)
                                call TransposeIntegralArr(TempArr_abab_2,TempArr_abba_2,tnmo)
                                TempArr_abba_2(:) = -TempArr_abba_2(:)
                                TempArr_aaaa_2(:) = TempArr_abab_2(:) + TempArr_abba_2(:)

                                RelaxXDyall(y,pp) = RelaxXDyall(y,pp) + DDOT(tnmo_sq,TwoRDM_aaaa(:,k),1,TempArr_aaaa_2(:),1)
                                RelaxXDyall(y,pp) = RelaxXDyall(y,pp) + DDOT(tnmo_sq,TwoRDM_abab(:,k),1,TempArr_abab_2(:),1)
                                RelaxXDyall(y,pp) = RelaxXDyall(y,pp) + DDOT(tnmo_sq,TwoRDM_abba(:,k),1,TempArr_abba_2(:),1)
!
!                                !Actually we only want the symmetrised versions...
!!                                RelaxXDyall(y,pp) = RelaxXDyall(y,pp) + DDOT(tnmo_sq,TwoRDM_aaaa(:,k),1,TempArr_aaaa_2(:),1)
!!                                RelaxXDyall(y,pp) = RelaxXDyall(y,pp) + DDOT(tnmo_sq,TwoRDM_abab(:,k),1,TempArr_abab_2(:),1)
!!                                RelaxXDyall(y,pp) = RelaxXDyall(y,pp) + DDOT(tnmo_sq,TwoRDM_abba(:,k),1,TempArr_abba_2(:),1)
!                                RelaxXDyall(y,pp) = RelaxXDyall(y,pp) + DDOT(tnmo_sq,TwoRDM_aaaa(:,k),1,TempArr_abab_2(:),1)
!!                                RelaxXDyall(y,pp) = RelaxXDyall(y,pp) + DDOT(tnmo_sq,TwoRDM_abab(:,k),1,TempArr_abab_2(:),1)

                            enddo
                        enddo

                    endif
                    if((.not.tDyall).or.tZerothRelaxBoth) then
                        !calculate f_j^i 2RDM_yi^xj
                        !This is part of the X_y^x (y fast) matrix for the orbital relaxation correction for the 
                        !non-Dyall zeroth order hamiltonian
                        !At the moment, for a given pp, we have all 2RDM_{pp,Q}^{RS}

                        !pp is looping over y
                        do x=1,tnmo

                            TempArr_aaaa(:,:) = 0.0_dp
                            TempArr_abab(:,:) = 0.0_dp

                            do i=1,tnmo
                                do j=1,tnmo
                                    xj_pair = FindSqPairInd(x,j,tnmo)
                                    TempArr_aaaa(i,j) = TwoRDM_aaaa(xj_pair,i)  !i is fast. I don't think this will matter
                                    TempArr_abab(i,j) = TwoRDM_abab(xj_pair,i)  !since GenFockOrb is symmetric
                                enddo
                            enddo

                            RelaxXFock(pp,x) = DDOT(tnmo_sq,GenFockOrb,1,TempArr_aaaa,1)
                            RelaxXFock(pp,x) = RelaxXFock(pp,x) + DDOT(tnmo_sq,GenFockOrb,1,TempArr_abab,1)

                        enddo
                    endif

                endif

                !Calculate the 2RDM part of Lambda_Tld.
                !Perform contraction of U_q^u G_pu^sr
                !Nasty DGEMMs!
                call DGEMM('N','T',tnmo_sq,tnmo,tnmo,1.0_dp,TwoRDM_aaaa,tnmo_sq,GammaF,tnmo,0.0_dp,Lambda_Tld_aaaa,tnmo_sq)
                call DGEMM('N','T',tnmo_sq,tnmo,tnmo,1.0_dp,TwoRDM_abab,tnmo_sq,GammaF,tnmo,0.0_dp,Lambda_Tld_abab,tnmo_sq)
                call DGEMM('N','T',tnmo_sq,tnmo,tnmo,1.0_dp,TwoRDM_abba,tnmo_sq,GammaF,tnmo,0.0_dp,Lambda_Tld_abba,tnmo_sq)

                do q=1,tnmo
                    pqPair = FindSqPairInd(q,pp,tnmo)

                    call TransposeIntegralArr_inplace(Lambda_Tld_aaaa(:,q),tnmo)
                    call TransposeIntegralArr_inplace(Lambda_Tld_abab(:,q),tnmo)
                    call TransposeIntegralArr_inplace(Lambda_Tld_abba(:,q),tnmo)

                    write(unit_Lambda_Tld_write_aaaa,rec=pqPair) (Lambda_Tld_aaaa(loop,q),loop=1,tnmo_sq) 
                    write(unit_Lambda_Tld_write_abab,rec=pqPair) (Lambda_Tld_abab(loop,q),loop=1,tnmo_sq) 
                    write(unit_Lambda_Tld_write_abba,rec=pqPair) (Lambda_Tld_abba(loop,q),loop=1,tnmo_sq) 
                enddo

            enddo   !End do over pp

            close(TwoRDM_read)
            deallocate(TwoRDM_aaaa)
            deallocate(TwoRDM_abab)
            deallocate(TwoRDM_abba)
            deallocate(TempArr_aaaa)
            deallocate(TempArr_abab)
            deallocate(TempArr_abba)
            deallocate(Lambda_Tld_aaaa,Lambda_Tld_abab,Lambda_Tld_abba)

!**********! This is the original - much slower - way of reading in the 2RDM - Also no perm sym in read in **********!
!            allocate(TwoRDM(tnmo_sq),stat=ierr)
!            if(ierr.ne.0) call stop_all(t_r,"Error in allocation")
!
!            do ispin=1,3
!                !Run over aaaa, abab, abba
!                if(ispin.eq.1) then
!                    unit_write = unit_2RDM_write_aaaa
!                elseif(ispin.eq.2) then
!                    unit_write = unit_2RDM_write_abab
!                elseif(ispin.eq.3) then
!                    unit_write = unit_2RDM_write_abba
!                endif
!
!                do i=1,tnmo
!                    do j=1,tnmo
!                        rewind(TwoRDM_read)
!                        ijInd=FindSqPairInd(i,j,tnmo)
!
!                        TwoRDM(:) = 0.0_dp
!                        do while(.true.)
!                            read(TwoRDM_read,"(4I6,G25.17)",iostat=ios) p,q,r,s,TwoRDMVal
!                            if(ios.gt.0) call stop_all(t_r,"Error reading TwoRDM")
!                            if(ios.lt.0) exit
!
!                            if(IsCorrectSpin(p,q,r,s,ispin)) then
!                                call Conv2Spat(p,q,r,s)
!                                if(FindSqPairInd(p,q,tnmo).eq.ijInd) then
!                                    rsPair = FindSqPairInd(r,s,tnmo)
!
!                                    if((abs(TwoRDM(rsPair)-TwoRDMVal).gt.5.D-4).and.(abs(TwoRDM(rsPair)).gt.1.D-8)) then
!                                        write(6,*) "Overwriting 2RDM value with different value..."
!                                        write(6,*) TwoRDM(rsPair),TwoRDMVal,p,q,r,s,ispin
!                                        call warning(t_r,"Error in spin symmetry in read in Two-RDM")
!                                    endif
!                                    TwoRDM(rsPair) = TwoRDMVal
!                                endif
!                            endif
!                        enddo
!
!                        write(unit_write,rec=ijInd) (TwoRDM(loop),loop=1,tnmo_sq) 
!
!                    enddo
!                enddo
!            enddo
!            deallocate(TwoRDM)

!Test
!            allocate(TwoRDM(tnmo_sq))
!            do i=1,tnmo_sq
!                read(unit_2RDM_write_abab,rec=i) (TwoRDM(loop),loop=1,tnmo_sq)
!                do p=1,tnmo_sq
!                    write(27,*) i,p,TwoRDM(p)
!                enddo
!            enddo
!            deallocate(TwoRDM)

            close(unit_2RDM_write_aaaa)
            close(unit_2RDM_write_abab)
            close(unit_2RDM_write_abba)
            close(unit_Lambda_Tld_write_aaaa)
            close(unit_Lambda_Tld_write_abab)
            close(unit_Lambda_Tld_write_abba)
            if(tDyall) close(g_ints_read)

        else    !tHFRDMs

            !We need to construct the part of Lambda_Bar that we have already done when reading in the 2RDMs
            !Do this in a really inefficient way, since it should be quick anyway
            Lambda_Bar(:,:) = 0.0_dp
            do p=1,tnocc
                pInd = OccOrbs(p)
                do s=1,tnocc
                    sInd = OccOrbs(s)

                    do q=1,tnocc
                        qInd = OccOrbs(q)
                        do r=1,tnocc
                            rInd = OccOrbs(r)
                            fac = 0.0_dp
                            if((pInd.eq.rInd).and.(qInd.eq.sInd)) fac = fac + 1.0_dp
                            if((pInd.eq.sInd).and.(qInd.eq.rInd)) fac = fac - 2.0_dp
!                            if(fac.ne.0.D0) write(6,*) Lambda_Bar(pInd,sInd),GenFockOrb(rInd,qInd),fac

                            Lambda_Bar(pInd,sInd) = Lambda_Bar(pInd,sInd) - GenFockOrb(rInd,qInd)*fac
                        enddo
                    enddo
                enddo
            enddo

        endif   !Endif tHFRDMs

        !We now have the other bits of Lambda_Bar to deal with.
        !First, -G_p^s f_r^q G_q^r
        !Multiplied by 2 since {q,r} can be aa or bb
        fac = 2.0_dp * DDOT(tnmo_sq,GenFockOrb,1,OneRDM,1)
        do p=1,tnmo
            do s=1,tnmo
                Lambda_Bar(p,s) = Lambda_Bar(p,s) - OneRDM(p,s)*fac
            enddo
        enddo

        allocate(TempArr_aaaa(tnmo,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc err")
        !Final bit of Lambda_Bar: 
        ! + f_r^q G_p^r G_q^s
        call DGEMM('N','N',tnmo,tnmo,tnmo,1.0_dp,GenFockOrb,tnmo,OneRDM,tnmo,0.0_dp,TempArr_aaaa,tnmo)
        call DGEMM('T','N',tnmo,tnmo,tnmo,1.0_dp,TempArr_aaaa,tnmo,OneRDM,tnmo,1.0_dp,Lambda_Bar,tnmo)
        deallocate(TempArr_aaaa)

        !Now deal with Lambda_Tld.
        !Lambda_Tld_pq^rs = 1/2 G_p^t f_t^u Lambda_qu^sr
        !The twoRDM part of Lambda_Tld is written out to file - read this back in.
        !Alternatively, if this is tHFRDMs, then calc on fly.
        if(.not.tHFRDMs) then
            unit_Lambda_Tld_read_aaaa=get_free_unit()
            open(unit_Lambda_Tld_read_aaaa,file='L_Tld_Tmp_aaaa',status='unknown',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
            unit_Lambda_Tld_read_abab=get_free_unit()
            open(unit_Lambda_Tld_read_abab,file='L_Tld_Tmp_abab',status='unknown',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
            unit_Lambda_Tld_read_abba=get_free_unit()
            open(unit_Lambda_Tld_read_abba,file='L_Tld_Tmp_abba',status='unknown',form='unformatted',    &
                access='direct',action='read',recl=reclen_nmo_sq)
        endif

        !Also open new files to write out the final (LambdaTld_pq)^rs term
        unit_Lambda_Tld_write_aaaa=get_free_unit()
        open(unit_Lambda_Tld_write_aaaa,file='L_Tld_Dir_aaaa',status='unknown',form='unformatted',    &
            access='direct',action='write',recl=reclen_nmo_sq)
        unit_Lambda_Tld_write_abab=get_free_unit()
        open(unit_Lambda_Tld_write_abab,file='L_Tld_Dir_abab',status='unknown',form='unformatted',    &
            access='direct',action='write',recl=reclen_nmo_sq)
        unit_Lambda_Tld_write_abba=get_free_unit()
        open(unit_Lambda_Tld_write_abba,file='L_Tld_Dir_abba',status='unknown',form='unformatted',    &
            access='direct',action='write',recl=reclen_nmo_sq)

        TempArr_aaaa_2(:) = 0.0_dp
        TempArr_abab_2(:) = 0.0_dp
        TempArr_abba_2(:) = 0.0_dp

        do p=1,tnmo
            do q=1,tnmo
                pqpair = FindSqPairInd(p,q,tnmo)

                !Obtain (X_pq)^rs
                if(.not.tHFRDMs) then
                    read(unit_Lambda_Tld_read_aaaa,rec=pqPair) (TempArr_aaaa_2(loop),loop=1,tnmo_sq)
                    read(unit_Lambda_Tld_read_abab,rec=pqPair) (TempArr_abab_2(loop),loop=1,tnmo_sq)
                    read(unit_Lambda_Tld_read_abba,rec=pqPair) (TempArr_abba_2(loop),loop=1,tnmo_sq)
                else
                    TempArr_aaaa_2(:) = 0.0_dp
                    TempArr_abab_2(:) = 0.0_dp
                    TempArr_abba_2(:) = 0.0_dp

                    do r=1,tnocc
                        if(.not.IsOrbOcc(q)) exit       !This line isnt checked!
                        rInd = OccOrbs(r)
                        do s=1,tnocc
                            sInd = OccOrbs(s)
                            rspair = FindSqPairInd(rInd,sInd,tnmo)
                            do u=1,tnocc
                                uInd = OccOrbs(u) 

                                if((uInd.eq.rInd).and.(q.eq.sInd).and.(uInd.ne.q)) then
                                    !TwoRDM is 1
                                    TempArr_aaaa_2(rspair) = TempArr_aaaa_2(rspair) + GammaF(p,uInd)
                                    TempArr_abab_2(rspair) = TempArr_abab_2(rspair) + GammaF(p,uInd)
                                elseif((uInd.eq.sInd).and.(q.eq.rInd).and.(uInd.ne.q)) then
                                    !TwoRDM is -1
                                    TempArr_aaaa_2(rspair) = TempArr_aaaa_2(rspair) - GammaF(p,uInd)
                                    TempArr_abba_2(rspair) = TempArr_abba_2(rspair) - GammaF(p,uInd)
                                elseif((uInd.eq.rInd).and.(q.eq.sInd).and.(uInd.eq.q)) then
                                    !All orbs equal
                                    TempArr_abab_2(rspair) = TempArr_abab_2(rspair) + GammaF(p,uInd)
                                    TempArr_abba_2(rspair) = TempArr_abba_2(rspair) - GammaF(p,uInd)
                                endif

                            enddo
                        enddo
                    enddo
                endif

                !Now include the oneRDM parts of this.
                !No abba component from the 1RDM parts? Only aaaa and abab?
                do r=1,tnmo
                    do s=1,tnmo
                        rspair = FindSqPairInd(r,s,tnmo)
                        TempArr_aaaa_2(rspair) = TempArr_aaaa_2(rspair) + OrbOrbMat(p,s)*OneRDM(q,r) - OrbOrbMat(p,r)*OneRDM(q,s)
                        TempArr_abab_2(rspair) = TempArr_abab_2(rspair) - OrbOrbMat(p,r)*OneRDM(q,s)
                        TempArr_abba_2(rspair) = TempArr_abba_2(rspair) + OrbOrbMat(p,s)*OneRDM(q,r)
                    enddo
                enddo

                !Now write out.
                write(unit_Lambda_Tld_write_aaaa,rec=pqPair) (TempArr_aaaa_2(loop),loop=1,tnmo_sq) 
                write(unit_Lambda_Tld_write_abab,rec=pqPair) (TempArr_abab_2(loop),loop=1,tnmo_sq) 
                write(unit_Lambda_Tld_write_abba,rec=pqPair) (TempArr_abba_2(loop),loop=1,tnmo_sq) 

            enddo
        enddo

        deallocate(TempArr_aaaa_2,TempArr_abab_2,TempArr_abba_2)
        deallocate(OrbOrbMat,GammaF)
        if(.not.tHFRDMs) then
            !Delete temp storage.
            close(unit_Lambda_Tld_read_aaaa,status='delete')
            close(unit_Lambda_Tld_read_abab,status='delete')
            close(unit_Lambda_Tld_read_abba,status='delete')
        endif
        close(unit_lambda_Tld_write_aaaa)
        close(unit_lambda_Tld_write_abab)
        close(unit_lambda_Tld_write_abba)
            
!       !Test
!        unit_Lambda_Tld_read_aaaa=get_free_unit()
!        open(unit_Lambda_Tld_read_aaaa,file='L_Tld_Dir_aaaa',status='unknown',form='unformatted',    &
!            access='direct',action='read',recl=reclen_nmo_sq)
!        unit_Lambda_Tld_read_abab=get_free_unit()
!        open(unit_Lambda_Tld_read_abab,file='L_Tld_Dir_abab',status='unknown',form='unformatted',    &
!            access='direct',action='read',recl=reclen_nmo_sq)
!        unit_Lambda_Tld_read_abba=get_free_unit()
!        open(unit_Lambda_Tld_read_abba,file='L_Tld_Dir_abba',status='unknown',form='unformatted',    &
!            access='direct',action='read',recl=reclen_nmo_sq)
!
!        allocate(TempArr_aaaa_2(tnmo_sq))
!        allocate(TempArr_abab_2(tnmo_sq))
!        allocate(TempArr_abba_2(tnmo_sq))
!
!        do p=1,tnmo
!            do q=1,tnmo
!                pqpair = FindSqPairInd(p,q,tnmo)
!                read(unit_Lambda_Tld_read_aaaa,rec=pqpair) (TempArr_aaaa_2(loop),loop=1,tnmo_sq)
!                read(unit_Lambda_Tld_read_abab,rec=pqpair) (TempArr_abab_2(loop),loop=1,tnmo_sq)
!                read(unit_Lambda_Tld_read_abba,rec=pqpair) (TempArr_abba_2(loop),loop=1,tnmo_sq)
!                counter = 1
!                do r=1,tnmo
!                    do s=1,tnmo
!                        write(71,"(G25.10,4I5)") TempArr_aaaa_2(counter),p,q,r,s
!                        write(72,"(G25.10,4I5)") TempArr_abab_2(counter),p,q,r,s
!                        write(73,"(G25.10,4I5)") TempArr_abba_2(counter),p,q,r,s
!                        counter = counter + 1
!                    enddo
!                enddo
!            enddo
!        enddo
!
!
!
!        do p=1,tnmo
!            do s=1,tnmo
!                write(27,"(G25.10,2I5)") Lambda_Bar(p,s),p,s
!            enddo
!        enddo

    end subroutine ReadTwoRDM

    subroutine Fill2RDMArr(ispin,p,q,r,s,Val,TwoRDM_aaaa,TwoRDM_abab,TwoRDM_abba)
        implicit none
        integer, intent(in) :: ispin,p,q,r,s
        real(dp), intent(in) :: Val
        real(dp), intent(inout) :: TwoRDM_aaaa(tnmo_sq,tnmo),TwoRDM_abab(tnmo_sq,tnmo),TwoRDM_abba(tnmo_sq,tnmo)
        integer :: pair
        character(*), parameter :: t_r='Fill2RDMArr'

        if(ispin.eq.1) then
            pair = FindSqPairInd(r,s,tnmo)
            if((abs(TwoRDM_aaaa(pair,q)-Val).gt.5.D-4).and.(abs(TwoRDM_aaaa(pair,q)).gt.1.D-8)) then
                write(6,*) "Overwriting 2RDM value with different value..."
                write(6,*) TwoRDM_aaaa(pair,q),Val,p,q,r,s,ispin
                call warning(t_r,"Error in spin symmetry in read in TwoRDM")
            endif
            TwoRDM_aaaa(pair,q) = Val
        elseif(ispin.eq.2) then
            pair = FindSqPairInd(r,s,tnmo)
            if((abs(TwoRDM_abab(pair,q)-Val).gt.5.D-4).and.(abs(TwoRDM_abab(pair,q)).gt.1.D-8)) then
                write(6,*) "Overwriting 2RDM value with different value..."
                write(6,*) TwoRDM_abab(pair,q),Val,p,q,r,s,ispin
                call warning(t_r,"Error in spin symmetry in read in TwoRDM")
            endif
            TwoRDM_abab(pair,q) = Val
        elseif(ispin.eq.3) then
            pair = FindSqPairInd(r,s,tnmo)
            if((abs(TwoRDM_abba(pair,q)-Val).gt.5.D-4).and.(abs(TwoRDM_abba(pair,q)).gt.1.D-8)) then
                write(6,*) "Overwriting 2RDM value with different value..."
                write(6,*) TwoRDM_abba(pair,q),Val,p,q,r,s,ispin
                call warning(t_r,"Error in spin symmetry in read in TwoRDM")
            endif
            TwoRDM_abba(pair,q) = Val
        endif
    end subroutine Fill2RDMArr


    subroutine Conv2Spat(p,q,r,s)
        implicit none
        integer, intent(inout) :: p,q,r,s

        if(mod(p,2).eq.0) then
            p=p/2
        else
            p=(p+1)/2
        endif
        if(mod(q,2).eq.0) then
            q=q/2
        else
            q=(q+1)/2
        endif
        if(mod(r,2).eq.0) then
            r=r/2
        else
            r=(r+1)/2
        endif
        if(mod(s,2).eq.0) then
            s=s/2
        else
            s=(s+1)/2
        endif
    end subroutine Conv2Spat

    !Tells you which spin, the four index quantity X_pq^rs is
    !If ispin = 1, then p,q,r,s must all be alpha / beta orbitals
    !If ispin = 2, then p,r must be alpha and q,s beta (or vice versa)
    !If ispin = 3, then p,s must be alpha and q,r beta (or vice versa)
    pure integer function WhichSpin(p,q,r,s)
        integer, intent(in) :: p,q,r,s
        integer :: pspin,qspin,rspin,sspin

        pspin = mod(p,2)    !0 = alpha, 1 = beta
        qspin = mod(q,2)
        rspin = mod(r,2)
        sspin = mod(s,2)
        if((pspin.eq.qspin).and.(pspin.eq.rspin).and.(pspin.eq.sspin)) then
            !aaaa or bbbb
            WhichSpin = 1
        elseif((pspin.eq.rspin).and.(qspin.eq.sspin).and.(pspin.ne.qspin)) then
            !abab / baba
            WhichSpin = 2
        elseif((qspin.eq.rspin).and.(pspin.eq.sspin).and.(qspin.ne.sspin)) then
            !abba / baab
            WhichSpin = 3
        else
            WhichSpin = 0   !Forbidden spin - matrix element must be 0
        endif

    end function WhichSpin

    !Read in p,q,r,s as spin orbitals, and check spin
    !If ispin = 1, then p,q,r,s must all be alpha / beta orbitals
    !If ispin = 2, then p,r must be alpha and q,s beta (or vice versa)
    !If ispin = 3, then p,s must be alpha and q,r beta (or vice versa)
    pure logical function IsCorrectSpin(p,q,r,s,ispin)
        integer, intent(in) :: p,q,r,s,ispin
        integer :: pspin,qspin,rspin,sspin

        pspin = mod(p,2)    !0 = alpha, 1 = beta
        qspin = mod(q,2)
        rspin = mod(r,2)
        sspin = mod(s,2)

        if((pspin.eq.qspin).and.(pspin.eq.rspin).and.(pspin.eq.sspin)) then
            !aaaa or bbbb
            if(ispin.eq.1) then
                IsCorrectSpin = .true.
            else
                IsCorrectSpin = .false.
            endif
        elseif((pspin.eq.rspin).and.(qspin.eq.sspin).and.(pspin.ne.qspin)) then
            !abab / baba
            if(ispin.eq.2) then
                IsCorrectSpin = .true.
            else
                IsCorrectSpin = .false.
            endif
        elseif((qspin.eq.rspin).and.(pspin.eq.sspin).and.(qspin.ne.sspin)) then
            !abba / baab
            if(ispin.eq.3) then
                IsCorrectSpin = .true.
            else
                IsCorrectSpin = .false.
            endif
        else
            IsCorrectSpin = .false.
        endif
    end function IsCorrectSpin

    subroutine StoreRInts() 
        integer :: unit_r_read,unit_r_write_ccoo,ierr,unit_r_write_oocc,loop
        integer :: rdump_unit
        integer(i2) , allocatable :: Indices(:,:)
        real(dp) , allocatable :: Buf(:),r_pairs(:,:),r_pairs_oocc(:,:)
        character(*), parameter :: t_r='StoreRInts'
        integer :: i,j,ios,ibuf,ijInd,i_Ind,j_Ind,k_Ind,l_Ind,ijInd_oocc
        real(dp) :: Z
        integer(i8) :: Maxlength,length
        integer(i2) :: MaxFilei,MaxFilej,MaxFilek,MaxFilel 
        logical :: ccoo_exists,oocc_exists

        inquire(file='R_Dir_cc.oo',exist=ccoo_exists)
        inquire(file='R_Dir_oo.cc',exist=oocc_exists)
        if(ccoo_exists.and.oocc_exists) then
            write(6,"(A)") "Direct access F12 files found..."
            write(6,"(A)") "Skipping reading from F12DUMPBIN..."
            return
        endif
        
        if(tRDUMP) then
            rdump_unit=get_free_unit()
            open(rdump_unit,file='RDUMP',status='unknown')
        endif

        !First, read in all primitive g integrals.
        !We want the first two indices to be over entire space (+ CABS).
        unit_r_read=get_free_unit()
        if(tDaltonFormat) then
            open(unit_r_read,file='MO_F12',status='old',form='unformatted',access='sequential',action='read')
        else
            if(tReadBin) then
                open(unit_r_read,file='F12DUMPBIN',status='old',form='unformatted',action='read')
            else
                open(unit_r_read,file='F12DUMP',status='old',form='formatted',action='read')
            endif
        endif

        !For the r integrals, from 4traf, the indices should be (physical):
        !i = complete
        !j = complete
        !k = orbital 
        !l = orbital
 !       CALL FindMaxIndices(unit_r_read,MaxFilei,MaxFilej,MaxFilek,MaxFilel)
!        write(6,"(A,4I5)") "Index upper limits for the g file are (physical notation): ",MaxFilei,MaxFilej,MaxFilek,MaxFilel

        !Integrals written out in a square array - all kl for a given p'q' written to one record. No perm sym.
        !R_Dir_cc.oo writes out a seperate record for all pairs in the complete space. Each record holds all orbital:orbital pairs
        unit_r_write_ccoo=get_free_unit()
        open(unit_r_write_ccoo,file='R_Dir_cc.oo',status='unknown',form='unformatted',   &
            access='direct',action='write',recl=reclen_nmo_sq)

        !R_Dir_oo.cc writes out records for all orbital:orbital pairs, and each record has total orbital pairs integrals
        unit_r_write_oocc=get_free_unit()
        open(unit_r_write_oocc,file='R_Dir_oo.cc',status='unknown',form='unformatted',  &
            access='direct',action='write',recl=reclen_ntmo_sq)

        rewind(unit_r_read)
        if(tDaltonFormat) then
            read(unit_r_read) maxlength
!            write(6,*) "Maxlength = ",maxlength
            allocate(Indices(4,MaxLength),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err 1")
            allocate(Buf(MaxLength),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err 1a")
        else
            length = 1  !Consider it as 1 integral per 'buffer'
        endif

        allocate(r_pairs(tnmo_sq,tntmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc err 2")
        allocate(r_pairs_oocc(tntmo_sq,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc err 2")
        
        write(6,'(F10.3,A)') 2.0_dp*real(tnmo_sq,dp)*real(tntmo_sq,dp)*8.0_dp*9.31322575D-10, &
            " Gb required to write out R_Dir_ccoo and R_Dir_oocc files."
        call flush(6)
        
!        open(24,file="FR-Full",status='unknown')

        !Loop over all i to store
        do i=1,tntmo

            !fill r_pairs with all < i' J' | K L >
            !fill r_pairs_oocc with all < i J | K' L' >
            r_pairs(:,:) = 0.0_dp
            r_pairs_oocc(:,:) = 0.0_dp

            rewind(unit_r_read)
            if(tDaltonFormat) read(unit_r_read) maxlength
            do while(.true.)
                !Read all integrals for each orbital i
                if(tDaltonFormat) then
                    read(unit_r_read,iostat=ios) length,Indices(1:4,1:length),Buf(1:length)
                else
                    if(tReadBin) then
                        read(unit_r_read,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                    else
                        read(unit_r_read,*,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                    endif
                endif
                if(ios.gt.0) call stop_all(t_r,"Error reading integrals")
                if((length.le.0).or.(ios.lt.0)) exit    !Finished reading all integrals.

                do ibuf=1,length
                    !Indices are read in in chemical notation. We store in PHYSICAL
                    !Find pair index - alternatively, we could back construct the ij pair, and compare to this
                    if(tDaltonFormat) then
                        i_Ind=int(Indices(1,ibuf),i4)
                        j_Ind=int(Indices(3,ibuf),i4)
                        k_Ind=int(Indices(2,ibuf),i4)
                        l_Ind=int(Indices(4,ibuf),i4)
                        Z = Buf(ibuf)
                    endif
                    !No need to consider the possibility of 2 index integrals in this file.
                        
                    if(tRDUMP) then
                        if(i.eq.1) write(rdump_unit,*) Z,i_Ind,k_Ind,j_Ind,l_Ind
                    endif

                    if((i_Ind.ne.i).and.(j_Ind.ne.i).and.(k_Ind.ne.i).and.(l_Ind.ne.i)) cycle

                    !At least one of the indices is equal to i
                    !Fill arrays, taking into account permutational sym, and size of spaces
                    if(i_Ind.eq.i) then
                        !<ij|kl>
                        if(max(k_Ind,l_Ind).le.tnmo) then
                            r_pairs(FindSqPairInd(k_Ind,l_Ind,tnmo),j_Ind) = Z 
                        endif
                        if(max(i_Ind,j_Ind).le.tnmo) then
                            r_pairs_oocc(FindSqPairInd(k_Ind,l_Ind,tntmo),j_Ind) = Z
                        endif
                        if(.not.tComplexOrbs) then
                            !<il|kj>
                            if(max(k_Ind,j_Ind).le.tnmo) then
                                r_pairs(FindSqPairInd(k_Ind,j_Ind,tnmo),l_Ind) = Z 
                            endif
                            if(max(i_Ind,l_Ind).le.tnmo) then
                                r_pairs_oocc(FindSqPairInd(k_Ind,j_Ind,tntmo),l_Ind) = Z 
                            endif
                        endif
                    endif

                    if(j_Ind.eq.i) then
                        !<ji|lk>
                        if(max(l_Ind,k_Ind).le.tnmo) then
                            r_pairs(FindSqPairInd(l_Ind,k_Ind,tnmo),i_Ind) = Z 
                        endif
                        if(max(j_Ind,i_Ind).le.tnmo) then
                            r_pairs_oocc(FindSqPairInd(l_Ind,k_Ind,tntmo),i_Ind) = Z
                        endif
                        if(.not.tComplexOrbs) then
                            !<jk|li>
                            if(max(l_Ind,i_Ind).le.tnmo) then
                                r_pairs(FindSqPairInd(l_Ind,i_Ind,tnmo),k_Ind) = Z 
                            endif
                            if(max(j_Ind,k_Ind).le.tnmo) then
                                r_pairs_oocc(FindSqPairInd(l_Ind,i_Ind,tntmo),k_Ind) = Z 
                            endif
                        endif
                    endif

                    if(k_Ind.eq.i) then
                        !<kl|ij>
                        if(max(i_Ind,j_Ind).le.tnmo) then
                            r_pairs(FindSqPairInd(i_Ind,j_Ind,tnmo),l_Ind) = Z
                        endif
                        if(max(k_Ind,l_Ind).le.tnmo) then
                            r_pairs_oocc(FindSqPairInd(i_Ind,j_Ind,tntmo),l_Ind) = Z
                        endif
                        if(.not.tComplexOrbs) then
                            !<kj|il>
                            if(max(i_Ind,l_Ind).le.tnmo) then
                                r_pairs(FindSqPairInd(i_Ind,l_Ind,tnmo),j_Ind) = Z
                            endif
                            if(max(k_Ind,j_Ind).le.tnmo) then
                                r_pairs_oocc(FindSqPairInd(i_Ind,l_Ind,tntmo),j_Ind) = Z
                            endif
                        endif
                    endif

                    if(l_Ind.eq.i) then
                        !<lk|ji>
                        if(max(j_Ind,i_Ind).le.tnmo) then
                            r_pairs(FindSqPairInd(j_Ind,i_Ind,tnmo),k_Ind) = Z
                        endif
                        if(max(l_Ind,k_Ind).le.tnmo) then
                            r_pairs_oocc(FindSqPairInd(j_Ind,i_Ind,tntmo),k_Ind) = Z
                        endif
                        if(.not.tComplexOrbs) then
                            !<li|jk>
                            if(max(j_Ind,k_Ind).le.tnmo) then
                                r_pairs(FindSqPairInd(j_Ind,k_Ind,tnmo),i_Ind) = Z
                            endif
                            if(max(l_Ind,i_Ind).le.tnmo) then
                                r_pairs_oocc(FindSqPairInd(j_Ind,k_Ind,tntmo),i_Ind) = Z 
                            endif
                        endif
                    endif

                enddo

            enddo


            !All {k,l},j for the given i are now stored. Write these out.

            do j=1,tntmo

                ijInd = FindSqPairInd(i,j,tntmo)
                write(unit_r_write_ccoo,rec=ijInd) (r_pairs(loop,j),loop=1,tnmo_sq)

                if(max(i,j).le.tnmo) then
                    ijInd = FindSqPairInd(i,j,tnmo)
                    write(unit_r_write_oocc,rec=ijInd) (r_pairs_oocc(loop,j),loop=1,tntmo_sq)
                endif

!                do k=1,tnmo
!                    do l=1,tnmo
!                        loop=loop+1
!                        write(24,"(4I4,F15.8)") I_Pair_ind,J_Pair_Ind,k,l,r_pairs(loop)
!                    enddo
!                enddo

            enddo
        enddo

        !All ij pairs have now been written out in direct access file.
        close(unit_r_read)
        close(unit_r_write_oocc)
        close(unit_r_write_ccoo)
        if(tRDUMP) then
            close(rdump_unit)
        endif
!        close(24)
        if(tDaltonFormat) then
            deallocate(Indices,Buf)
        endif
        deallocate(r_pairs,r_pairs_oocc)

    end subroutine StoreRInts


    
    !Routine to fill r_pairs for a given ij pair which runs over all orbital space pairs
    !ijInd is the index of the ij pair we are trying to fill.
    !The entries of the r_pairs run over complete pairs.
    subroutine PutinRpairs_oocc(i,j,k,l,intgl,r_pairs,ijInd)
        implicit none
        integer, intent(in) :: i,j,k,l,ijInd
        real(dp), intent(in) :: intgl
        real(dp), intent(inout) :: r_pairs(tntmo_sq)

        if((max(i,j).le.tnmo).and.(FindSqPairInd(i,j,tnmo).eq.ijInd)) then
            !<ij|kl>
            r_pairs(FindSqPairInd(k,l,tntmo))=intgl
        endif
        if((max(j,i).le.tnmo).and.(FindSqPairInd(j,i,tnmo).eq.ijInd)) then
            !<ji|lk>
            r_pairs(FindSqPairInd(l,k,tntmo))=intgl
        endif
        if((max(k,l).le.tnmo).and.(FindSqPairInd(k,l,tnmo).eq.ijInd)) then
            !<kl|ij>
            r_pairs(FindSqPairInd(i,j,tntmo))=intgl
        endif
        if((max(l,k).le.tnmo).and.(FindSqPairInd(l,k,tnmo).eq.ijInd)) then
            !<lk|ji>
            r_pairs(FindSqPairInd(j,i,tntmo))=intgl
        endif
        if((max(k,j).le.tnmo).and.(FindSqPairInd(k,j,tnmo).eq.ijInd)) then
            !<kj|il>
            r_pairs(FindSqPairInd(i,l,tntmo))=intgl
        endif
        if((max(l,i).le.tnmo).and.(FindSqPairInd(l,i,tnmo).eq.ijInd)) then
            !<li|jk>
            r_pairs(FindSqPairInd(j,k,tntmo))=intgl
        endif
        if((max(i,l).le.tnmo).and.(FindSqPairInd(i,l,tnmo).eq.ijInd)) then
            !<il|kj>
            r_pairs(FindSqPairInd(k,j,tntmo))=intgl
        endif
        if((max(j,k).le.tnmo).and.(FindSqPairInd(j,k,tnmo).eq.ijInd)) then
            !<jk|li>
            r_pairs(FindSqPairInd(l,i,tntmo))=intgl
        endif 

    end subroutine PutinRpairs_oocc

    !Routine to fill r_pairs for a given ij pair which runs over all complete space pairs
    !ijInd is the index of the ij pair we are trying to fill.
    !The entries of the r_pairs run over all *orbtial* space pairs only.
    subroutine PutinRpairs_ccoo(i,j,k,l,intgl,r_pairs,ijInd)
        implicit none
        integer, intent(in) :: i,j,k,l,ijInd
        real(dp), intent(in) :: intgl
        real(dp), intent(inout) :: r_pairs(tnmo_sq)

        if(FindSqPairInd(i,j,tntmo).eq.ijInd) then
            !<ij|kl>
            if((k.le.tnmo).and.(l.le.tnmo)) then
                r_pairs(FindSqPairInd(k,l,tnmo))=intgl
            endif
        endif
        if(FindSqPairInd(j,i,tntmo).eq.ijInd) then
            !<ji|lk>
            if((l.le.tnmo).and.(k.le.tnmo)) then
                r_pairs(FindSqPairInd(l,k,tnmo))=intgl
            endif
        endif
        if(FindSqPairInd(k,l,tntmo).eq.ijInd) then
            !<kl|ij>
            if((i.le.tnmo).and.(j.le.tnmo)) then
                r_pairs(FindSqPairInd(i,j,tnmo))=intgl
            endif
        endif
        if(FindSqPairInd(l,k,tntmo).eq.ijInd) then
            !<lk|ji>
            if((j.le.tnmo).and.(i.le.tnmo)) then
                r_pairs(FindSqPairInd(j,i,tnmo))=intgl
            endif
        endif
        if(FindSqPairInd(k,j,tntmo).eq.ijInd) then
            !<kj|il>
            if((i.le.tnmo).and.(l.le.tnmo)) then
                r_pairs(FindSqPairInd(i,l,tnmo))=intgl
            endif
        endif
        if(FindSqPairInd(l,i,tntmo).eq.ijInd) then
            !<li|jk>
            if((j.le.tnmo).and.(k.le.tnmo)) then
                r_pairs(FindSqPairInd(j,k,tnmo))=intgl
            endif
        endif
        if(FindSqPairInd(i,l,tntmo).eq.ijInd) then
            !<il|kj>
            if((k.le.tnmo).and.(j.le.tnmo)) then
                r_pairs(FindSqPairInd(k,j,tnmo))=intgl
            endif
        endif
        if(FindSqPairInd(j,k,tntmo).eq.ijInd) then
            !<jk|li>
            if((l.le.tnmo).and.(i.le.tnmo)) then
                r_pairs(FindSqPairInd(l,i,tnmo))=intgl
            endif
        endif 

    end subroutine PutinRpairs_ccoo

!Combine the routines CreateR_TLD and CreateR_TLD_Phi
    subroutine CreateR_TLD_Both()
        use input_data
        use basis, only: FourIndSym
        implicit none
        integer :: unit_raaaa_write,unit_rabab_write,unit_rabba_write,unit_r_read_ccoo,ierr,p_prime,q_prime,i,j
        integer :: unit_2RDM_read_aaaa,unit_2RDM_read_abab,unit_2RDM_read_abba,ijPair,jiPair
        real(dp), allocatable :: r_xy(:),r_yx(:),raaaa(:),rabab(:),rabba(:),Gamma_xy_aaaa(:),Gamma_xy_abab(:),Gamma_xy_abba(:)
        real(dp), allocatable :: Antisym_r_xy_aa(:),Antisym_r_xy_ab(:),Antisym_r_xy_ba(:)
        integer :: pq_pair,loop,x,y,xyPair,v,w,vw_pair,pq_pair2,wv_pair
        integer :: unit_aaaa_Phi_write,unit_abab_Phi_write,unit_abba_Phi_write
        integer :: unit_Phi_read_aaaa,unit_Phi_read_abab,unit_Phi_read_abba 
        real(dp), allocatable :: Phi_xy_aaaa(:),Phi_xy_abab(:),Phi_xy_abba(:),Phi_aaaa(:),Phi_abab(:),Phi_abba(:)
        integer :: xInd,yInd,q_sym,vw_sym,CurrPercent,LastPercent
        real(dp) :: DDOT,Memory!,temp1,temp2,temp3,temp4
        logical :: tCalcPhiTld
        character(*), parameter :: t_r="CreateR_TLD"
        
!        write(6,*) "Entering R_TLD..."
!        call flush(6)
        
        !RTld Integrals written out in a square array
        !all kl (orbital pairs) for a given ij (all pairs) written to one record. No perm sym.
        unit_raaaa_write=get_free_unit()
        open(unit_raaaa_write,file='RTld_aaaa',status='unknown',form='unformatted',access='direct',action='write',recl=reclen_ntmo_sq)
        unit_rabab_write=get_free_unit()
        open(unit_rabab_write,file='RTld_abab',status='unknown',form='unformatted',access='direct',action='write',recl=reclen_ntmo_sq)
        unit_rabba_write=get_free_unit()
        open(unit_rabba_write,file='RTld_abba',status='unknown',form='unformatted',access='direct',action='write',recl=reclen_ntmo_sq)

        unit_aaaa_Phi_write = get_free_unit()
        open(unit_aaaa_Phi_write,file='RTld_Phi_aaaa',status='unknown',form='unformatted',access='direct', &
            action='write',recl=reclen_ntmo_sq)
        unit_abab_Phi_write = get_free_unit()
        open(unit_abab_Phi_write,file='RTld_Phi_abab',status='unknown',form='unformatted',access='direct', &
            action='write',recl=reclen_ntmo_sq)
        unit_abba_Phi_write=get_free_unit()
        open(unit_abba_Phi_write,file='RTld_Phi_abba',status='unknown',form='unformatted',access='direct', &
            action='write',recl=reclen_ntmo_sq)

        unit_r_read_ccoo=get_free_unit()
        open(unit_r_read_ccoo,file='R_Dir_cc.oo',status='old',form='unformatted',access='direct',action='read',recl=reclen_nmo_sq)

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
        
        unit_Phi_read_aaaa=get_free_unit()
        open(unit_Phi_read_aaaa,file='Phi_Dir_aaaa',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)
        unit_Phi_read_abab=get_free_unit()
        open(unit_Phi_read_abab,file='Phi_Dir_abab',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)
        unit_Phi_read_abba=get_free_unit()
        open(unit_Phi_read_abba,file='Phi_Dir_abba',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)

        allocate(r_xy(tnmo_sq),stat=ierr)
        allocate(r_yx(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error allocating")
        r_xy(:)=0.D0
        r_yx(:)=0.D0

        allocate(Gamma_xy_aaaa(tnmo_sq),stat=ierr)
        allocate(Gamma_xy_abab(tnmo_sq),stat=ierr)
        allocate(Gamma_xy_abba(tnmo_sq),stat=ierr)
        allocate(Phi_xy_aaaa(tnmo_sq),stat=ierr)
        allocate(Phi_xy_abab(tnmo_sq),stat=ierr)
        allocate(Phi_xy_abba(tnmo_sq),stat=ierr)
        allocate(raaaa(tntmo_sq),stat=ierr)
        allocate(rabab(tntmo_sq),stat=ierr)
        allocate(rabba(tntmo_sq),stat=ierr)
        allocate(Phi_aaaa(tntmo_sq),stat=ierr)
        allocate(Phi_abab(tntmo_sq),stat=ierr)
        allocate(Phi_abba(tntmo_sq),stat=ierr)
        allocate(Antisym_r_xy_aa(tnmo_sq),stat=ierr)
        allocate(Antisym_r_xy_ab(tnmo_sq),stat=ierr)
        allocate(Antisym_r_xy_ba(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error allocating")

        if(tCalcCombo) then
            Memory = 2.0_dp*(real(tnmo_sq,dp) + (real(tnmo-tnfrz,dp)**2.0_dp))*real(tntmo_sq,dp)*8.0_dp*9.31322575D-10
        else
            Memory = 3.0_dp*(real(tnmo_sq,dp) + (real(tnmo-tnfrz,dp)**2.0_dp))*real(tntmo_sq,dp)*8.0_dp*9.31322575D-10
        endif
        write(6,'(F10.3,A)') Memory," Gb required to write out RTld and RTld_Phi files."
        call flush(6)

        LastPercent = 0
        do v=1,tnmo
            !This calculates the percentage complete, taking into account the triangular nature of w
            CurrPercent = nint((real(200*(v-1),dp)/real(tnmo,dp))-(real(100*(v-1)*(v-1),dp)/real(tnmo_sq,dp)))
            if(CurrPercent.gt.LastPercent) then
                write(6,"(I4,A)") CurrPercent , " % complete"
                call flush(6)
                LastPercent = CurrPercent
            endif
            do w=v,tnmo

                if((.not.tHFRDMs).or.(IsOrbOcc(v).and.IsOrbOcc(w))) then
                    !Do normal contraction and write out as normal

                    vw_sym = ieor(OrbSyms(v)-1,OrbSyms(w)-1)

                    tCalcPhiTld = .not.(IsOrbFrz(v).or.IsOrbFrz(w)) !Logical indicating whether to calculate PhiTld

                    vw_pair = FindSqPairInd(v,w,tnmo)
                    wv_pair = FindSqPairInd(w,v,tnmo)
                    raaaa(:)=0.0_dp
                    rabab(:)=0.0_dp
                    rabba(:)=0.0_dp

                    Phi_aaaa(:)=0.0_dp
                    Phi_abab(:)=0.0_dp
                    Phi_abba(:)=0.0_dp

                    !Construct Gamma^xy
                    if(tHFRDMs) then
                        Gamma_xy_aaaa(:)=0.0_dp
                        Gamma_xy_abab(:)=0.0_dp
                        Gamma_xy_abba(:)=0.0_dp
                        do x=1,tnocc
                            xInd = OccOrbs(x) 
                            if(IsOrbFrz(xInd)) cycle
                            do y=1,tnocc
                                yInd = OccOrbs(y)
                                if(IsOrbFrz(yInd)) cycle
                                if((v.eq.xInd).and.(w.eq.yInd).and.(xInd.eq.yInd)) then
                                    !All four indices the same
                                    xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                    Gamma_xy_abab(xyPair) = 1.0_dp
                                    Gamma_xy_abba(xyPair) = -1.0_dp
                                elseif((v.eq.xInd).and.(w.eq.yInd)) then
                                    xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                    Gamma_xy_aaaa(xyPair) = 1.0_dp
                                    Gamma_xy_abab(xyPair) = 1.0_dp
                                elseif((v.eq.yInd).and.(w.eq.xInd)) then
                                    xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                    Gamma_xy_aaaa(xyPair) = -1.0_dp
                                    Gamma_xy_abba(xyPair) = -1.0_dp
                                endif
                            enddo
                        enddo
                    else
                        read(unit_2RDM_read_aaaa,rec=vw_pair) (Gamma_xy_aaaa(loop),loop=1,tnmo_sq)
                        read(unit_2RDM_read_abab,rec=vw_pair) (Gamma_xy_abab(loop),loop=1,tnmo_sq)
                        read(unit_2RDM_read_abba,rec=vw_pair) (Gamma_xy_abba(loop),loop=1,tnmo_sq)

                        !Zero frozen matrix elements
                        do i=1,tnfrz
                            do j=1,tnmo
                                ijPair = FindSqPairInd(FrzOrbs(i),j,tnmo)
                                jiPair = FindSqPairInd(j,FrzOrbs(i),tnmo)
                                Gamma_xy_aaaa(ijPair) = 0.0_dp
                                Gamma_xy_aaaa(jiPair) = 0.0_dp
                                Gamma_xy_abab(ijPair) = 0.0_dp
                                Gamma_xy_abab(jiPair) = 0.0_dp
                                Gamma_xy_abba(ijPair) = 0.0_dp
                                Gamma_xy_abba(jiPair) = 0.0_dp
                            enddo
                        enddo

                    endif
                    
                    !Construct Phi^xy
                    if(tCalcPhiTld) then
                        !Do not read in Phi if v/w is frozen, since do not need r_tld_Phi over frozen orbs
                        read(unit_Phi_read_aaaa,rec=vw_pair) (Phi_xy_aaaa(loop),loop=1,tnmo_sq)
                        read(unit_Phi_read_abab,rec=vw_pair) (Phi_xy_abab(loop),loop=1,tnmo_sq)
                        read(unit_Phi_read_abba,rec=vw_pair) (Phi_xy_abba(loop),loop=1,tnmo_sq)
                        !They also do not want to be contracted over frozen orbital, so set these to zero.
                        do i=1,tnfrz
                            do j=1,tnmo
                                ijPair = FindSqPairInd(FrzOrbs(i),j,tnmo)
                                jiPair = FindSqPairInd(j,FrzOrbs(i),tnmo)
                                Phi_xy_aaaa(ijPair) = 0.0_dp
                                Phi_xy_aaaa(jiPair) = 0.0_dp
                                Phi_xy_abab(ijPair) = 0.0_dp
                                Phi_xy_abab(jiPair) = 0.0_dp
                                Phi_xy_abba(ijPair) = 0.0_dp
                                Phi_xy_abba(jiPair) = 0.0_dp
                            enddo
                        enddo
                    endif

                    do p_prime=1,tntmo

                        if(p_prime.le.tnmo) then
                            q_sym = ieor(OrbSyms(p_prime)-1,vw_sym) + 1
                        else
                            q_sym = ieor(AuxOrbSyms(p_prime)-1,vw_sym) + 1
                        endif

                        !Loop over q_prime for the correct symmetry only, from StartSym_n(x)mo(q_sym) to EndSym_n(x)mo(q_sym)
                        !Also ensure q_prime >= p_prime
                        do q_prime = max(StartSym_nmo(q_sym),p_prime) , EndSym_nmo(q_sym)
    !                    do q_prime = p_prime,tntmo

                            !Read R_{pq}^{xy} - since pq here goes over entire space, use tntmo
                            pq_pair = FindSqPairInd(p_prime,q_prime,tntmo)
                            read(unit_r_read_ccoo,rec=pq_pair) (r_xy(loop),loop=1,tnmo_sq)

                            !Since we are storing p'q' over entire space 
                            !We want **p'** to be the fast index
                            pq_pair2 = FindSqPairInd(q_prime,p_prime,tntmo)

                            call TransposeIntegralArr(r_xy,r_yx,tnmo)

                            !Create Antisym_aa, ab and ba
                            Antisym_r_xy_aa(:)=(r_xy(:)-r_yx(:))/4.0_dp
                            Antisym_r_xy_ab(:)=((3.0_dp/8.0_dp)*r_xy(:)) + (r_yx(:)/8.0_dp)
                            Antisym_r_xy_ba(:)=(r_xy(:)/(-8.0_dp)) - ((3.0_dp/8.0_dp)*r_yx(:))

                            !Contract over {xy}
                            if((v.ne.w).and.(p_prime.ne.q_prime)) then  !otherwise it must = 0
                                raaaa(pq_pair2) = DDOT(tnmo_sq,Antisym_r_xy_aa,1,Gamma_xy_aaaa,1)    
                                if(tCalcPhiTld) then
                                    Phi_aaaa(pq_pair2) = DDOT(tnmo_sq,Antisym_r_xy_aa,1,Phi_xy_aaaa,1)
                                endif
                            endif

                            !For abab/abba terms, need to sum over contraction over {x,y} = alpha, beta *and* beta, alpha
                            !We use permutational symmetry to assume that they are the same
                            rabab(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Gamma_xy_abab,1)
                            if(tCalcPhiTld) then
                                Phi_abab(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Phi_xy_abab,1)
                            endif
    !                        rabab(pq_pair2) = rabab(pq_pair2) + DDOT(tnmo_sq,Antisym_r_xy_ba,1,Gamma_xy_abba,1)

                            rabba(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Gamma_xy_abba,1)
                            if(tCalcPhiTld) then
                                Phi_abba(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Phi_xy_abba,1)
                            endif
    !                        rabba(pq_pair2) = rabba(pq_pair2) + DDOT(tnmo_sq,Antisym_r_xy_ba,1,Gamma_xy_abab,1)

                            !More permutational antisymmetry
                            if((q_prime.ne.p_prime).or.(.not.tCalcCombo)) then
                                raaaa(pq_pair) = -raaaa(pq_pair2)
                                rabab(pq_pair) = -rabba(pq_pair2)
                                rabba(pq_pair) = -rabab(pq_pair2)
                                if(tCalcPhiTld) then
                                    Phi_aaaa(pq_pair) = -Phi_aaaa(pq_pair2)
                                    Phi_abba(pq_pair) = -Phi_abab(pq_pair2)
                                    Phi_abab(pq_pair) = -Phi_abba(pq_pair2)
                                endif
                            endif

                        enddo

                        !Now loop over q_prime for the correct symmetry, but now in the CABS space.
                        !Also ensure q_prime >= p_prime
                        do q_prime = max(StartSym_nxmo(q_sym),p_prime) , EndSym_nxmo(q_sym)

                            !Read R_{pq}^{xy} - since pq here goes over entire space, use tntmo
                            pq_pair = FindSqPairInd(p_prime,q_prime,tntmo)
                            read(unit_r_read_ccoo,rec=pq_pair) (r_xy(loop),loop=1,tnmo_sq)

                            !Since we are storing p'q' over entire space 
                            !We want **p'** to be the fast index
                            pq_pair2 = FindSqPairInd(q_prime,p_prime,tntmo)

                            call TransposeIntegralArr(r_xy,r_yx,tnmo)

                            !Create Antisym_aa, ab and ba
                            Antisym_r_xy_aa(:)=(r_xy(:)-r_yx(:))/4.0_dp
                            Antisym_r_xy_ab(:)=((3.0_dp/8.0_dp)*r_xy(:)) + (r_yx(:)/8.0_dp)
                            Antisym_r_xy_ba(:)=(r_xy(:)/(-8.0_dp)) - ((3.0_dp/8.0_dp)*r_yx(:))

                            if((v.ne.w).and.(p_prime.ne.q_prime)) then
                                !Contract over {xy}
                                raaaa(pq_pair2) = DDOT(tnmo_sq,Antisym_r_xy_aa,1,Gamma_xy_aaaa,1)    

                                if(tCalcPhiTld) then
                                    Phi_aaaa(pq_pair2) = DDOT(tnmo_sq,Antisym_r_xy_aa,1,Phi_xy_aaaa,1)
                                endif
                            endif

                            !For abab/abba terms, need to sum over contraction over {x,y} = alpha, beta *and* beta, alpha
                            !We use permutational symmetry to assume that they are the same
                            rabab(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Gamma_xy_abab,1)
                            if(tCalcPhiTld) then
                                Phi_abab(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Phi_xy_abab,1)
                            endif
    !                        rabab(pq_pair2) = rabab(pq_pair2) + DDOT(tnmo_sq,Antisym_r_xy_ba,1,Gamma_xy_abba,1)

                            rabba(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Gamma_xy_abba,1)
                            if(tCalcPhiTld) then
                                Phi_abba(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Phi_xy_abba,1)
                            endif
    !                        rabba(pq_pair2) = rabba(pq_pair2) + DDOT(tnmo_sq,Antisym_r_xy_ba,1,Gamma_xy_abab,1)

                            !More permutational antisymmetry
                            if((q_prime.ne.p_prime).or.(.not.tCalcCombo)) then
                                raaaa(pq_pair) = -raaaa(pq_pair2)
                                rabab(pq_pair) = -rabba(pq_pair2)
                                rabba(pq_pair) = -rabab(pq_pair2)
                                if(tCalcPhiTld) then
                                    Phi_aaaa(pq_pair) = -Phi_aaaa(pq_pair2)
                                    Phi_abba(pq_pair) = -Phi_abab(pq_pair2)
                                    Phi_abab(pq_pair) = -Phi_abba(pq_pair2)
                                endif
                            endif

                        enddo
                    enddo

                    !All kl for the given ij are now stored. Write these out
                    write(unit_raaaa_write,rec=vw_pair) (raaaa(i),i=1,tntmo_sq)
                    write(unit_rabab_write,rec=vw_pair) (rabab(i),i=1,tntmo_sq)
                    if(.not.tCalcCombo) write(unit_rabba_write,rec=vw_pair) (rabba(i),i=1,tntmo_sq)

                    write(unit_raaaa_write,rec=wv_pair) (-raaaa(i),i=1,tntmo_sq)
                    write(unit_rabab_write,rec=wv_pair) (-rabba(i),i=1,tntmo_sq)
                    if(.not.tCalcCombo) write(unit_rabba_write,rec=wv_pair) (-rabab(i),i=1,tntmo_sq)

                    if(tCalcPhiTld) then
                        write(unit_aaaa_Phi_write,rec=vw_pair) (Phi_aaaa(i),i=1,tntmo_sq)
                        write(unit_abab_Phi_write,rec=vw_pair) (Phi_abab(i),i=1,tntmo_sq)
                        if(.not.tCalcCombo) write(unit_abba_Phi_write,rec=vw_pair) (Phi_abba(i),i=1,tntmo_sq)

                        write(unit_aaaa_Phi_write,rec=wv_pair) (-Phi_aaaa(i),i=1,tntmo_sq)
                        write(unit_abab_Phi_write,rec=wv_pair) (-Phi_abba(i),i=1,tntmo_sq)
                        if(.not.tCalcCombo) write(unit_abba_Phi_write,rec=wv_pair) (-Phi_abab(i),i=1,tntmo_sq)
                    endif

                else
                    !Since we are using a single reference, and one or both of v and w are not occupied orbitals,
                    !we know that both the twoRDM and phi are zero. Therefore we can set rtld and rtld_phi to
                    !zero for all pq, without calculating anything.
                    vw_pair = FindSqPairInd(v,w,tnmo)
                    wv_pair = FindSqPairInd(w,v,tnmo)

                    raaaa(:) = 0.0_dp
                    rabab(:) = 0.0_dp
                    rabba(:) = 0.0_dp

                    write(unit_raaaa_write,rec=vw_pair) (raaaa(i),i=1,tntmo_sq)
                    write(unit_rabab_write,rec=vw_pair) (rabab(i),i=1,tntmo_sq)
                    if(.not.tCalcCombo) write(unit_rabba_write,rec=vw_pair) (rabba(i),i=1,tntmo_sq)

                    !These are zero, so we don't need to worry about the minus sign!
                    write(unit_raaaa_write,rec=wv_pair) (raaaa(i),i=1,tntmo_sq)
                    write(unit_rabab_write,rec=wv_pair) (rabba(i),i=1,tntmo_sq)
                    if(.not.tCalcCombo) write(unit_rabba_write,rec=wv_pair) (rabab(i),i=1,tntmo_sq)

                    if(IsOrbFrz(v).and.IsOrbFrz(w)) then
                        write(6,*) v,w
                        call stop_all(t_r,"Error here - frozen orbitals should be considered occupied")
                    endif

                    Phi_aaaa(:) = 0.0_dp
                    Phi_abab(:) = 0.0_dp
                    Phi_abba(:) = 0.0_dp

                    write(unit_aaaa_Phi_write,rec=vw_pair) (Phi_aaaa(i),i=1,tntmo_sq)
                    write(unit_abab_Phi_write,rec=vw_pair) (Phi_abab(i),i=1,tntmo_sq)
                    if(.not.tCalcCombo) write(unit_abba_Phi_write,rec=vw_pair) (Phi_abba(i),i=1,tntmo_sq)

                    !These are zero, so we don't need to worry about the minus sign!
                    write(unit_aaaa_Phi_write,rec=wv_pair) (Phi_aaaa(i),i=1,tntmo_sq)
                    write(unit_abab_Phi_write,rec=wv_pair) (Phi_abba(i),i=1,tntmo_sq)
                    if(.not.tCalcCombo) write(unit_abba_Phi_write,rec=wv_pair) (Phi_abab(i),i=1,tntmo_sq)

                endif   !Whether we actually need to do the contraction or not 

            enddo   !loop over v,w pair
        enddo

        close(unit_raaaa_write)
        close(unit_rabab_write)
        close(unit_rabba_write)
        if(.not.tHFRDMs) then
            close(unit_2RDM_read_aaaa)
            close(unit_2RDM_read_abab)
            close(unit_2RDM_read_abba)
        endif
        close(unit_r_read_ccoo)
        close(unit_aaaa_Phi_write)
        close(unit_abab_Phi_write)
        close(unit_abba_Phi_write)
        close(unit_Phi_read_aaaa)
        close(unit_Phi_read_abab)
        close(unit_Phi_read_abba)
        deallocate(raaaa,rabab,rabba,r_xy,r_yx,Gamma_xy_aaaa,Gamma_xy_abab,Gamma_xy_abba)
        deallocate(Antisym_r_xy_aa,Antisym_r_xy_ab,Antisym_r_xy_ba)
        deallocate(Phi_aaaa,Phi_abab,Phi_abba,Phi_xy_aaaa,Phi_xy_abab,Phi_xy_abba)


    end subroutine CreateR_TLD_Both

!Create the contracted intermediates resulting from: 
!R_{p'q'}^{xy} Gamma_{xy}^{rs} = Raaaa_{p'q'}^{rs}      spaces - {p'q'}:all, {rs}:orbital
!R_{p'q'}^{xy} Gamma_{xy}^{rs} = Rabab_{p'q'}^{rs}      spaces - {p'q'}:all, {rs}:orbital
!R_{p'q'}^{xy} Gamma_{xy}^{rs} = Rabba_{p'q'}^{rs}      spaces - {p'q'}:all, {rs}:orbital
!Fast index is the *first* (i.e. r)
!These are written to direct access disk.
    subroutine CreateR_TLD()
        use input_data
        implicit none
        integer :: unit_raaaa_write,unit_rabab_write,unit_rabba_write,unit_r_read_ccoo,ierr,p_prime,q_prime,i
        integer :: unit_2RDM_read_aaaa,unit_2RDM_read_abab,unit_2RDM_read_abba
        real(dp), allocatable :: r_xy(:),r_yx(:),raaaa(:),rabab(:),rabba(:),Gamma_xy_aaaa(:),Gamma_xy_abab(:),Gamma_xy_abba(:)
        real(dp), allocatable :: Antisym_r_xy_aa(:),Antisym_r_xy_ab(:),Antisym_r_xy_ba(:)
        integer :: pq_pair,loop,x,y,xyPair,v,w,vw_pair,pq_pair2,wv_pair
        integer :: xInd,yInd
        real(dp) :: DDOT!,temp1,temp2,temp3,temp4
        character(*), parameter :: t_r="CreateR_TLD"
        
!        write(6,*) "Entering R_TLD..."
!        call flush(6)
        
        !RTld Integrals written out in a square array
        !all kl (orbital pairs) for a given ij (all pairs) written to one record. No perm sym.
        unit_raaaa_write=get_free_unit()
        open(unit_raaaa_write,file='RTld_aaaa',status='unknown',form='unformatted',access='direct',action='write',recl=reclen_ntmo_sq)
        unit_rabab_write=get_free_unit()
        open(unit_rabab_write,file='RTld_abab',status='unknown',form='unformatted',access='direct',action='write',recl=reclen_ntmo_sq)
        unit_rabba_write=get_free_unit()
        open(unit_rabba_write,file='RTld_abba',status='unknown',form='unformatted',access='direct',action='write',recl=reclen_ntmo_sq)

        unit_r_read_ccoo=get_free_unit()
        open(unit_r_read_ccoo,file='R_Dir_cc.oo',status='old',form='unformatted',access='direct',action='read',recl=reclen_nmo_sq)

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

        allocate(r_xy(tnmo_sq),stat=ierr)
        allocate(r_yx(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error allocating")
        r_xy(:)=0.D0
        r_yx(:)=0.D0

        allocate(Gamma_xy_aaaa(tnmo_sq),stat=ierr)
        allocate(Gamma_xy_abab(tnmo_sq),stat=ierr)
        allocate(Gamma_xy_abba(tnmo_sq),stat=ierr)
        allocate(raaaa(tntmo_sq),stat=ierr)
        allocate(rabab(tntmo_sq),stat=ierr)
        allocate(rabba(tntmo_sq),stat=ierr)
        allocate(Antisym_r_xy_aa(tnmo_sq),stat=ierr)
        allocate(Antisym_r_xy_ab(tnmo_sq),stat=ierr)
        allocate(Antisym_r_xy_ba(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error allocating")

        do v=1,tnmo
            do w=v,tnmo

                vw_pair = FindSqPairInd(v,w,tnmo)
                wv_pair = FindSqPairInd(w,v,tnmo)
                raaaa(:)=0.0_dp
                rabab(:)=0.0_dp
                rabba(:)=0.0_dp

                !Construct Gamma^xy
                if(tHFRDMs) then
                    Gamma_xy_aaaa(:)=0.0_dp
                    Gamma_xy_abab(:)=0.0_dp
                    Gamma_xy_abba(:)=0.0_dp
                    do x=1,tnocc
                        xInd = OccOrbs(x) 
                        do y=1,tnocc
                            yInd = OccOrbs(y)
                            if((v.eq.xInd).and.(w.eq.yInd).and.(xInd.eq.yInd)) then
                                !All four indices the same
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Gamma_xy_abab(xyPair) = 1.0_dp
                                Gamma_xy_abba(xyPair) = -1.0_dp
                            elseif((v.eq.xInd).and.(w.eq.yInd)) then
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Gamma_xy_aaaa(xyPair) = 1.0_dp
                                Gamma_xy_abab(xyPair) = 1.0_dp
                            elseif((v.eq.yInd).and.(w.eq.xInd)) then
                                xyPair = FindSqPairInd(xInd,yInd,tnmo)
                                Gamma_xy_aaaa(xyPair) = -1.0_dp
                                Gamma_xy_abba(xyPair) = -1.0_dp
                            endif
                        enddo
                    enddo
                else
                    read(unit_2RDM_read_aaaa,rec=vw_pair) (Gamma_xy_aaaa(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abab,rec=vw_pair) (Gamma_xy_abab(loop),loop=1,tnmo_sq)
                    read(unit_2RDM_read_abba,rec=vw_pair) (Gamma_xy_abba(loop),loop=1,tnmo_sq)
                endif

                do p_prime=1,tntmo
                    do q_prime=p_prime,tntmo

                        !Read R_{pq}^{xy} - since pq here goes over entire space, use tntmo
                        pq_pair = FindSqPairInd(p_prime,q_prime,tntmo)
                        read(unit_r_read_ccoo,rec=pq_pair) (r_xy(loop),loop=1,tnmo_sq)

                        !Since we are storing p'q' over entire space 
                        !We want **p'** to be the fast index
                        pq_pair2 = FindSqPairInd(q_prime,p_prime,tntmo)

                        call TransposeIntegralArr(r_xy,r_yx,tnmo)

                        !Create Antisym_aa, ab and ba
                        Antisym_r_xy_aa(:)=(r_xy(:)-r_yx(:))/4.0_dp
                        Antisym_r_xy_ab(:)=((3.0_dp/8.0_dp)*r_xy(:)) + (r_yx(:)/8.0_dp)
                        Antisym_r_xy_ba(:)=(r_xy(:)/(-8.0_dp)) - ((3.0_dp/8.0_dp)*r_yx(:))

                        !Contract over {xy}
!                        do x=1,tnmo_sq
!                            write(6,*) x,r_xy(x),Gamma_xy(x)
!                        enddo
                        raaaa(pq_pair2) = DDOT(tnmo_sq,Antisym_r_xy_aa,1,Gamma_xy_aaaa,1)    

                        !For abab/abba terms, need to sum over contraction over {x,y} = alpha, beta *and* beta, alpha
                        !We use permutational symmetry to assume that they are the same
                        rabab(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Gamma_xy_abab,1)
!                        rabab(pq_pair2) = rabab(pq_pair2) + DDOT(tnmo_sq,Antisym_r_xy_ba,1,Gamma_xy_abba,1)

                        rabba(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Gamma_xy_abba,1)
!                        rabba(pq_pair2) = rabba(pq_pair2) + DDOT(tnmo_sq,Antisym_r_xy_ba,1,Gamma_xy_abab,1)

                        !More permutational antisymmetry
                        raaaa(pq_pair) = -raaaa(pq_pair2)
                        rabab(pq_pair) = -rabba(pq_pair2)
                        rabba(pq_pair) = -rabab(pq_pair2)

!                        temp1 = DDOT(tnmo_sq,Antisym_r_xy_ab,1,Gamma_xy_abab,1)
!                        temp2 = DDOT(tnmo_sq,Antisym_r_xy_ba,1,Gamma_xy_abba,1)
!                        temp3 = DDOT(tnmo_sq,Antisym_r_xy_ab,1,Gamma_xy_abba,1)
!                        temp4 = DDOT(tnmo_sq,Antisym_r_xy_ba,1,Gamma_xy_abab,1)
!
!                        if(abs(temp1-temp2).gt.1.D-9) call stop_all(t_r,"Not same")
!                        if(abs(temp3-temp4).gt.1.D-9) call stop_all(t_r,"Not same")

!                        write(6,*) DDOT(tnmo_sq,Antisym_r_xy_ab,1,Gamma_xy_abab,1), &
!                                   DDOT(tnmo_sq,Antisym_r_xy_ba,1,Gamma_xy_abba,1), &
!                                   DDOT(tnmo_sq,Antisym_r_xy_ab,1,Gamma_xy_abba,1), &
!                                   DDOT(tnmo_sq,Antisym_r_xy_ba,1,Gamma_xy_abab,1)

                    enddo
                enddo

                !All kl for the given ij are now stored. Write these out
                write(unit_raaaa_write,rec=vw_pair) (raaaa(i),i=1,tntmo_sq)
                write(unit_rabab_write,rec=vw_pair) (rabab(i),i=1,tntmo_sq)
                write(unit_rabba_write,rec=vw_pair) (rabba(i),i=1,tntmo_sq)

                write(unit_raaaa_write,rec=wv_pair) (-raaaa(i),i=1,tntmo_sq)
                write(unit_rabab_write,rec=wv_pair) (-rabba(i),i=1,tntmo_sq)
                write(unit_rabba_write,rec=wv_pair) (-rabab(i),i=1,tntmo_sq)

            enddo
        enddo

        close(unit_raaaa_write)
        close(unit_rabab_write)
        close(unit_rabba_write)
        if(.not.tHFRDMs) then
            close(unit_2RDM_read_aaaa)
            close(unit_2RDM_read_abab)
            close(unit_2RDM_read_abba)
        endif
        close(unit_r_read_ccoo)
        deallocate(raaaa,rabab,rabba,r_xy,r_yx,Gamma_xy_aaaa,Gamma_xy_abab,Gamma_xy_abba)
        deallocate(Antisym_r_xy_aa,Antisym_r_xy_ab,Antisym_r_xy_ba)

    end subroutine CreateR_TLD

    
!Create the contracted intermediates resulting from: 
!R_{p'q'}^{xy} Psi_{xy}^{rs} = Raaaa_{p'q'}^{rs}      spaces - {p'q'}:all, {rs}:orbital
!R_{p'q'}^{xy} Psi_{xy}^{rs} = Rabab_{p'q'}^{rs}      spaces - {p'q'}:all, {rs}:orbital
!R_{p'q'}^{xy} Psi_{xy}^{rs} = Rabba_{p'q'}^{rs}      spaces - {p'q'}:all, {rs}:orbital
!Fast index is the *first* (i.e. r)
!These are written to direct access disk.
    subroutine CreateR_TLD_Phi()
        use input_data
        implicit none
        integer :: unit_raaaa_write,unit_rabab_write,unit_rabba_write,unit_r_read_ccoo,ierr,p_prime,q_prime,i
        integer :: unit_Phi_read_aaaa,unit_Phi_read_abab,unit_Phi_read_abba
        real(dp), allocatable :: r_xy(:),r_yx(:),raaaa(:),rabab(:),rabba(:),Phi_xy_aaaa(:),Phi_xy_abab(:),Phi_xy_abba(:)
        real(dp), allocatable :: Antisym_r_xy_aa(:),Antisym_r_xy_ab(:),Antisym_r_xy_ba(:)
        integer :: pq_pair,loop,x,y,xyPair,v,w,vw_pair,pq_pair2,wv_pair
        integer :: xInd,yInd
        real(dp) :: DDOT!,temp1,temp2,temp3,temp4
        character(*), parameter :: t_r="CreateR_TLD_Psi"
!        real(dp) , allocatable :: raaaa_2(:),rabab_2(:),rabba_2(:)
!        integer :: yxpair,qppair,pqpair,p,q,wv_pair
        
!        write(6,*) "Entering R_TLD..."
!        call flush(6)
        
        !RTld Integrals written out in a square array
        !all kl (orbital pairs) for a given ij (all pairs) written to one record. No perm sym.
        unit_raaaa_write=get_free_unit()
        open(unit_raaaa_write,file='RTld_Phi_aaaa',status='unknown',form='unformatted',access='direct', &
            action='write',recl=reclen_ntmo_sq)
        unit_rabab_write=get_free_unit()
        open(unit_rabab_write,file='RTld_Phi_abab',status='unknown',form='unformatted',access='direct', &
            action='write',recl=reclen_ntmo_sq)
        unit_rabba_write=get_free_unit()
        open(unit_rabba_write,file='RTld_Phi_abba',status='unknown',form='unformatted',access='direct', &
            action='write',recl=reclen_ntmo_sq)

        unit_r_read_ccoo=get_free_unit()
        open(unit_r_read_ccoo,file='R_Dir_cc.oo',status='old',form='unformatted',access='direct', &
            action='read',recl=reclen_nmo_sq)

        unit_Phi_read_aaaa=get_free_unit()
        open(unit_Phi_read_aaaa,file='Phi_Dir_aaaa',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)
        unit_Phi_read_abab=get_free_unit()
        open(unit_Phi_read_abab,file='Phi_Dir_abab',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)
        unit_Phi_read_abba=get_free_unit()
        open(unit_Phi_read_abba,file='Phi_Dir_abba',status='unknown',form='unformatted',    &
            access='direct',action='read',recl=reclen_nmo_sq)

        allocate(r_xy(tnmo_sq),stat=ierr)
        allocate(r_yx(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error allocating")
        r_xy(:)=0.D0
        r_yx(:)=0.D0

        allocate(Phi_xy_aaaa(tnmo_sq),stat=ierr)
        allocate(Phi_xy_abab(tnmo_sq),stat=ierr)
        allocate(Phi_xy_abba(tnmo_sq),stat=ierr)
        allocate(raaaa(tntmo_sq),stat=ierr)
        allocate(rabab(tntmo_sq),stat=ierr)
        allocate(rabba(tntmo_sq),stat=ierr)
        allocate(Antisym_r_xy_aa(tnmo_sq),stat=ierr)
        allocate(Antisym_r_xy_ab(tnmo_sq),stat=ierr)
        allocate(Antisym_r_xy_ba(tnmo_sq),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Error allocating")

        do v=1,tnmo
            do w=v,tnmo

                vw_pair = FindSqPairInd(v,w,tnmo)
                wv_pair = FindSqPairInd(w,v,tnmo)
                raaaa(:)=0.0_dp
                rabab(:)=0.0_dp
                rabba(:)=0.0_dp

                !Construct Psi^xy
                read(unit_Phi_read_aaaa,rec=vw_pair) (Phi_xy_aaaa(loop),loop=1,tnmo_sq)
                read(unit_Phi_read_abab,rec=vw_pair) (Phi_xy_abab(loop),loop=1,tnmo_sq)
                read(unit_Phi_read_abba,rec=vw_pair) (Phi_xy_abba(loop),loop=1,tnmo_sq)

                do p_prime=1,tntmo
                    do q_prime=p_prime,tntmo

                        !Read R_{pq}^{xy} - since pq here goes over entire space, use tntmo
                        pq_pair = FindSqPairInd(p_prime,q_prime,tntmo)
                        read(unit_r_read_ccoo,rec=pq_pair) (r_xy(loop),loop=1,tnmo_sq)

                        !Since we are storing p'q' over entire space 
                        !We want **p'** to be the fast index
                        pq_pair2 = FindSqPairInd(q_prime,p_prime,tntmo)

                        call TransposeIntegralArr(r_xy,r_yx,tnmo)

                        !Create Antisym_aa, ab and ba
                        Antisym_r_xy_aa(:)=(r_xy(:)-r_yx(:))/4.0_dp
                        Antisym_r_xy_ab(:)=((3.0_dp/8.0_dp)*r_xy(:)) + (r_yx(:)/8.0_dp)
                        Antisym_r_xy_ba(:)=(r_xy(:)/(-8.0_dp)) - ((3.0_dp/8.0_dp)*r_yx(:))

                        !Contract over {xy}
!                        do x=1,tnmo_sq
!                            write(6,*) x,r_xy(x),Phi_xy(x)
!                        enddo
                        raaaa(pq_pair2) = DDOT(tnmo_sq,Antisym_r_xy_aa,1,Phi_xy_aaaa,1)    
                        raaaa(pq_pair) = -raaaa(pq_pair2)    

                        !For abab/abba terms, need to sum over contraction over {x,y} = alpha, beta *and* beta, alpha
                        !However, these terms should be the same by permutational symmetry.
                        !We assume that they are here...
                        rabab(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Phi_xy_abab,1)
!                        rabab(pq_pair2) = rabab(pq_pair2) + DDOT(tnmo_sq,Antisym_r_xy_ba,1,Phi_xy_abba,1)

                        rabba(pq_pair2) = 2.0_dp * DDOT(tnmo_sq,Antisym_r_xy_ab,1,Phi_xy_abba,1)

                        rabab(pq_pair) = -rabba(pq_pair2)
                        rabba(pq_pair) = -rabab(pq_pair2)
!                        rabba(pq_pair2) = rabba(pq_pair2) + DDOT(tnmo_sq,Antisym_r_xy_ba,1,Phi_xy_abab,1)
!                        
!                        temp1 = DDOT(tnmo_sq,Antisym_r_xy_ab,1,Phi_xy_abab,1)
!                        temp2 = DDOT(tnmo_sq,Antisym_r_xy_ba,1,Phi_xy_abba,1)
!                        temp3 = DDOT(tnmo_sq,Antisym_r_xy_ab,1,Phi_xy_abba,1)
!                        temp4 = DDOT(tnmo_sq,Antisym_r_xy_ba,1,Phi_xy_abab,1)
!
!                        if(abs(temp1-temp2).gt.1.D-9) call stop_all(t_r,"Not same")
!                        if(abs(temp3-temp4).gt.1.D-9) call stop_all(t_r,"Not same")

                    enddo
                enddo

                !All kl for the given ij are now stored. Write these out
                write(unit_raaaa_write,rec=vw_pair) (raaaa(i),i=1,tntmo_sq)
                write(unit_rabab_write,rec=vw_pair) (rabab(i),i=1,tntmo_sq)
                write(unit_rabba_write,rec=vw_pair) (rabba(i),i=1,tntmo_sq)
                
                !Use permutational symmetry: rtld_pq^vw = - rtld_pq^wv
                !Remember to also swap spins.
                write(unit_raaaa_write,rec=wv_pair) (-raaaa(i),i=1,tntmo_sq)
                write(unit_rabab_write,rec=wv_pair) (-rabba(i),i=1,tntmo_sq)
                write(unit_rabba_write,rec=wv_pair) (-rabab(i),i=1,tntmo_sq)

            enddo
        enddo

!        !Check that rtld^vw_pq = rtld^vw_qp
!        allocate(raaaa_2(tntmo_sq))
!        allocate(rabab_2(tntmo_sq))
!        allocate(rabba_2(tntmo_sq))
!        do x = 1,tnmo
!            do y = 1,tnmo
!                xypair = FindSqPairInd(x,y,tnmo)
!                yxpair = FindSqPairInd(y,x,tnmo)
!                read(unit_raaaa_write,rec=xypair) (raaaa(i),i=1,tntmo_sq)
!                read(unit_raaaa_write,rec=yxpair) (raaaa_2(i),i=1,tntmo_sq)
!                read(unit_rabab_write,rec=xypair) (rabab(i),i=1,tntmo_sq)
!                read(unit_rabab_write,rec=yxpair) (rabab_2(i),i=1,tntmo_sq)
!                read(unit_rabba_write,rec=xypair) (rabba(i),i=1,tntmo_sq)
!                read(unit_rabba_write,rec=yxpair) (rabba_2(i),i=1,tntmo_sq)
!                do p=1,tntmo
!                    do q=1,tntmo
!                        pqpair = FindSqPairInd(p,q,tntmo)
!                        qppair = FindSqPairInd(q,p,tntmo)
!
!                        if(abs(raaaa(pqpair)+raaaa(qppair)).gt.1.D-8) then
!                            call stop_all(t_r,"Not same")
!                        endif
!                        if(abs(rabab(pqpair)+rabba(qppair)).gt.1.D-8) then
!                            call stop_all(t_r,"Not same")
!                        endif
!                        if(abs(rabba(qppair)+rabab(pqpair)).gt.1.D-8) then
!                            call stop_all(t_r,"Not same")
!                        endif
!                    enddo
!                enddo
!            enddo
!        enddo
!        deallocate(raaaa_2,rabab_2,rabba_2)

        close(unit_raaaa_write)
        close(unit_rabab_write)
        close(unit_rabba_write)
        close(unit_Phi_read_aaaa)
        close(unit_Phi_read_abab)
        close(unit_Phi_read_abba)
        close(unit_r_read_ccoo)
        deallocate(raaaa,rabab,rabba,r_xy,r_yx,Phi_xy_aaaa,Phi_xy_abab,Phi_xy_abba)
        deallocate(Antisym_r_xy_aa,Antisym_r_xy_ab,Antisym_r_xy_ba)

    end subroutine CreateR_TLD_Phi
    
!Loads into memory the orbital:orbital block of the fock matrix
    subroutine LoadFockOrb()
        implicit none
        integer :: unit_fock,i,ii,jj,ierr
        integer(i8) :: maxlength,length
        integer(i2) , allocatable :: Indices(:,:)
        real(dp) , allocatable :: Buf(:)
        real(dp) :: xout(2)
        character(len=*), parameter :: t_r="LoadFockOrb"

        allocate(FockOrb(1:tnmo,1:tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Fock allocation err")
        !Orbital is first index, CABS second
        allocate(FockOrbCABS(1:tnmo,tnmo+1:tntmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Fock allocation err")
        allocate(FockCABS(tnmo+1:tntmo,tnmo+1:tntmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Fock allocation err")

        fockOrb(:,:)=0.0_dp
        fockOrbCABS(:,:)=0.0_dp
        fockCABS(:,:)=0.0_dp

        if(.not.tDaltonFormat) then
            !Assume that the *entire* basis is composed of fock eigenfunctions and so is diagonal
            !EBC and GBC satisfied exactly
            !fockOrbCABS matrix is strictly zero.
            do i=1,tntmo
                if(i.le.tnmo) then
                    fockOrb(i,i) = OrbEnergies(i)
                else
                    fockCABS(i,i) = OrbEnergies(i)
                endif
            enddo

        else
            unit_fock=get_free_unit()
            open(unit_fock,file='MO_F',status='old',form='unformatted',access='sequential',action='read')

            rewind(unit_fock)
            read(unit_fock) maxlength,xout(2)
            allocate(Indices(2,MaxLength),Buf(MaxLength))
            do while(.true.)
                read(unit_fock,iostat=ierr) length,Indices(1:2,1:length),Buf(1:length)
                IF(length.le.0) exit

                do i=1,length
                    ii=Indices(1,i)
                    jj=Indices(2,i)
                    if((ii.le.tnmo).and.(jj.le.tnmo)) then
                        !Orbital:Orbital block
                        fockOrb(ii,jj)=Buf(i)
                        fockOrb(jj,ii)=Buf(i)
                    elseif(((ii.le.tnmo).and.(jj.gt.tnmo)).or.((ii.gt.tnmo).and.(jj.le.tnmo))) then
                        !In the orbital,CABS block of the space
                        !Orbital is first index, CABS second
                        fockOrbCABS(min(ii,jj),max(ii,jj)) = Buf(i)
                        if(ii.eq.jj) call stop_all(t_r,"Error here in filling Fock matrix")
                    else
                        !In the CABS,CABS block of the space
                        fockCABS(ii,jj)=Buf(i)
                        fockCABS(jj,ii)=Buf(i)
                    endif
                enddo
            enddo


            deallocate(Buf,Indices)
            close(unit_fock)
        endif

    end subroutine LoadFockOrb

!From a slice of integrals, r^{kl}, return the transpose - r^{lk}
!Assumes 1D square array indexing
    subroutine TransposeIntegralArr_inplace(Arr,nOrb)
        implicit none
        real(dp), intent(inout) :: Arr(nOrb*nOrb)
        integer, intent(in) :: nOrb
        integer :: i,j,PairInd,PairIndSwap
        real(dp) :: temp

        do i=1,nOrb
            do j=i+1,nOrb
                !Run along triangle of array
                PairInd = FindSqPairInd(i,j,nOrb)
                PairIndSwap = FindSqPairInd(j,i,nOrb)

                temp=Arr(PairIndSwap)
                Arr(PairIndSwap)=Arr(PairInd)
                Arr(PairInd)=temp
            enddo
        enddo
    end subroutine TransposeIntegralArr_inplace

!From a slice of integrals, r^{kl}, return the transpose - r^{lk}
!Assumes 1D square array indexing
    subroutine TransposeIntegralArr(Arr,Arr2,nOrb)
        implicit none
        real(dp), intent(in) :: Arr(nOrb*nOrb)
        real(dp), intent(out) :: Arr2(nOrb*nOrb)
        integer, intent(in) :: nOrb
        integer :: i,j,PairInd,PairIndSwap
        real(dp) :: temp

        Arr2(:)=Arr(:)

        do i=1,nOrb
            do j=i+1,nOrb
                !Run along triangle of array
                PairInd = FindSqPairInd(i,j,nOrb)
                PairIndSwap = FindSqPairInd(j,i,nOrb)

                temp=Arr2(PairIndSwap)
                Arr2(PairIndSwap)=Arr2(PairInd)
                Arr2(PairInd)=temp
            enddo
        enddo
!
!        Arr2(:)=0.0_dp
!
!        do i=1,nOrb
!            do j=1,nOrb
!                PairInd = FindSqPairInd(i,j,nOrb)
!                PairIndSwap = FindSqPairInd(j,i,nOrb)
!
!                Arr2(PairIndSwap) = Arr(PairInd)
!            enddo
!        enddo

!        write(35,*) Arr(:)
!        write(36,*) Arr2(:)

    end subroutine TransposeIntegralArr
                
    !Write FCIDUMP
    subroutine CreateFCIDUMP()
        implicit none
        integer :: unit_dump,unit_g_read,i,ios,ibuf,n,m
        integer :: ierr
        integer(i8) :: maxlength,length
        integer(i2), allocatable :: Indices(:,:)
        real(dp), allocatable :: Buf(:)
!        integer(i2) :: Mini,Minj,Mink,Minl
        character(len=*), parameter :: t_r='CreateFCIDUMP'
        integer :: nOrb,nAux,nElec,Ms2,nProp(3),PropBitLen,iSym
        integer :: i_Ind,j_Ind,k_Ind,l_Ind
        real(dp) :: Z,GAM
        integer(i8) :: OrbSym(1000)
        logical :: UHF,UEG
        namelist /FCI/ nOrb,nAux,nElec,Ms2,OrbSym,iSym,UHF,UEG,GAM,PropBitLen,nProp

        unit_dump=get_free_unit()
        open(unit_dump,file='FCIDUMP',status='unknown',form='formatted')
        WRITE(unit_dump,'(2A6,I3,A7,I3,A5,I2,A)') '&FCI ','NORB=',tnmo,'NELEC=',NEl,',MS2=',0,','
        write(unit_dump,'(A9)',advance='no') 'ORBSYM='
        do i=1,tnmo
            write(unit_dump,'(I1,A1)',advance='no') OrbSyms(i),','
        enddo
        write(unit_dump,*)
        write(unit_dump,'(A7,I1,A12)') 'ISYM=',1,' UHF=.FALSE.'
        write(unit_dump,'(A5)') '&END'

        unit_g_read=get_free_unit()
        if(tDaltonFormat) then
            open(unit_g_read,file='MO_G',status='old',form='unformatted',access='sequential',action='read')
        else
            if(tReadBin) then
                open(unit_g_read,file='GDUMPBIN',status='old',form='unformatted',action='read')
            else
                open(unit_g_read,file='GDUMP',status='old',form='formatted',action='read')
            endif
        endif

!        call FindMinIndices(unit_g_read,Mini,Minj,Mink,Minl)
!        write(6,*) "Minimum indices of MO_G is (physical): ",Mini,Minj,Mink,Minl

        !For the g integrals, from 4traf, the indices should be (physical):
        !i = complete
        !j = complete
        !k = complete
        !l = orbital

        rewind(unit_g_read)
        if(tDaltonFormat) then
            read(unit_g_read) maxlength
!            write(6,*) "Maxlength = ",maxlength
            allocate(Indices(4,MaxLength),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err 1")
            allocate(Buf(MaxLength),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err 1a")
        else
            if(.not.tReadBin) read(unit_g_read,fci)
            length = 1
        endif

        do while(.true.)

            if(tDaltonFormat) then
                read(unit_g_read,iostat=ios) length,Indices(1:4,1:length),Buf(1:length)
            else
                if(tReadBin) then
                    read(unit_g_read,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                else
                    read(unit_g_read,*,iostat=ios) Z,i_Ind,k_Ind,j_Ind,l_Ind
                endif
            endif
            if(ios.gt.0) then
                call stop_all(t_r,"Error reading integrals")
            endif
            if((length.le.0).or.(ios.lt.0)) exit    !Finished reading all integrals.

            do ibuf=1,length
!                write(6,*) Indices(:,ibuf),Buf(ibuf)
                !Indices are read in in chemical notation (and written out to the FCIDUMP file). 
                if(tDaltonFormat) then
                    if((Indices(1,ibuf).le.tnmo).and.(Indices(2,ibuf).le.tnmo)      &
                        .and.(Indices(3,ibuf).le.tnmo).and.(Indices(4,ibuf).le.tnmo)) then
                        if(abs(Buf(ibuf)).gt.1.D-9) then
                            write(unit_dump,'(1X,G20.14,4I3)') Buf(ibuf),Indices(1,ibuf),Indices(2,ibuf),   &
                                Indices(3,ibuf),Indices(4,ibuf)
!                            write(6,'(1X,G20.14,4I3)') Buf(ibuf),Indices(1,ibuf),Indices(2,ibuf),Indices(3,ibuf),Indices(4,ibuf)
                        endif
                    endif
                else
                    if((i_Ind.ne.0).and.(j_Ind.ne.0).and.(k_Ind.ne.0).and.(l_Ind.ne.0)) then
                        if((i_Ind.le.tnmo).and.(j_Ind.le.tnmo).and.(k_Ind.le.tnmo).and.         &
                                (l_Ind.le.tnmo).and.(abs(Z).gt.1.0e-9_dp)) then
                            write(unit_dump,'(1X,G20.14,4I3)') Z,i_Ind,k_Ind,j_Ind,l_Ind
                        endif
                    endif
                endif
            enddo
        enddo
        if(tDaltonFormat) then
            deallocate(Indices)
            deallocate(Buf)
        endif
        close(unit_g_read)
        
        do n=1,tnmo
            do m=n,tnmo
                if(abs(tmat(n,m)).gt.1.D-9) then
                    write(unit_dump,'(1X,G20.14,4I3)') tmat(n,m),n,m,0,0
                endif
            enddo
        enddo

        !Now write out diagonal fock matrix elements
        do i=1,tnmo
            write(unit_dump,'(1X,G20.14,4I3)') fockOrb(i,i),i,0,0,0
        enddo
        write(unit_dump,'(1X,G20.14,4I3)') pot_nuc,0,0,0,0 

        close(unit_dump)

    end subroutine CreateFCIDUMP

!The minimum indices in a particular file, where the indices are returned in PHYSICAL NOTATION!
    subroutine FindMinIndices(iUnit,Mini,Minj,Mink,Minl)
        implicit none
        integer(i2) , intent(out) :: Mini,Minj,Mink,Minl
        integer , intent(in) :: iUnit
        integer(i8) :: MaxLength,Length
        integer(i2) , allocatable :: Indices(:,:)
        integer :: i,iBuf

        rewind(iUnit)
        read(iUnit) MaxLength
        allocate(Indices(4,MaxLength))
        Mini=30000
        Minj=30000
        Mink=30000
        Minl=30000
        iBuf=0
        do while(.true.)
            iBuf=iBuf+1
            write(6,*) "Attempting to read buffer: ",iBuf
            read(iUnit,err=11) length,Indices(1:4,1:length)
            if(length.gt.0) THEN
                do i=1,length
                    write(6,*) Indices(:,i)
                    IF(Indices(1,i).lt.Mini) Mini=Indices(1,i)
                    IF(Indices(2,i).lt.Mink) Mink=Indices(2,i)
                    IF(Indices(3,i).lt.Minj) Minj=Indices(3,i)
                    IF(Indices(4,i).lt.Minl) Minl=Indices(4,i)
                enddo
            else
                exit
            endif
        enddo
        deallocate(Indices)
        rewind(iUnit)
        return
11      stop 'Error when finding largest indices'
    end subroutine FindMinIndices

!Read in g integrals from 4traf file, and write them out to a direct access file,
!with each ij pair in a different index. At the moment, no triangular indexing or
!permutational symmetry at all.
!ij pairs will run over (complete, complete) basis pairs (each record)
!the kl pairs will run over (orbital, orbital) basis pairs (within record)
!The 4traf file MO_R runs over < complete, complete | orbital orbital >.
    subroutine StoreGInts() 
        integer :: unit_g_read,unit_g_write_ccoo,ierr,unit_g_write_oocc,loop
        integer(i2) , allocatable :: Indices(:,:)
        real(dp) , allocatable :: Buf(:),g_pairs(:,:),g_pairs_oocc(:,:)
        character(*), parameter :: t_r='StoreGInts'
        integer :: i,j,ios,ibuf,ijInd,i_Ind,j_Ind,k_Ind,l_Ind,ijInd_oocc,k,l,n
        integer(i8) :: Maxlength,length
        integer(i2) :: MaxFilei,MaxFilej,MaxFilek,MaxFilel 
        integer :: nOrb,nAux,nElec,Ms2,nProp(3),PropBitLen,iSym
        integer, allocatable :: Indices_long(:,:)
        real(dp) :: Z,GAM
        integer(i8) :: OrbSym(1000)
        logical :: UHF,UEG,ccoo_exists,oocc_exists
        namelist /FCI/ nOrb,nAux,nElec,Ms2,OrbSym,iSym,UHF,UEG,GAM,PropBitLen,nProp

        inquire(file='G_Dir_ccoo',exist=ccoo_exists)
        inquire(file='G_Dir_oocc',exist=oocc_exists)
        if(ccoo_exists.and.oocc_exists) then
            write(6,"(A)") "Direct access G files found..."
            write(6,"(A)") "Skipping reading from GDUMPBIN..."
            return
        endif

        !First, read in all primitive g integrals.
        !We want the first two indices to be over entire space (+ CABS).
        unit_g_read=get_free_unit()
        if(tDaltonFormat) then
            open(unit_g_read,file='MO_G',status='old',form='unformatted',access='sequential',action='read')
        else
            if(tReadBin) then
                open(unit_g_read,file='GDUMPBIN',status='old',form='unformatted',action='read')
            else
                open(unit_g_read,file='GDUMP',status='old',form='formatted',action='read')
            endif
        endif

        if(tDaltonFormat) then
            !For the g integrals, from 4traf, the indices should be (physical):
            !i = complete
            !j = complete
            !k = complete
            !l = orbital
            CALL FindMaxIndices(unit_g_read,MaxFilei,MaxFilej,MaxFilek,MaxFilel)
!            write(6,"(A,4I5)") "Index upper limits for the g file are (physical notation): ",MaxFilei,MaxFilej,MaxFilek,MaxFilel
        endif

        !Integrals written out in a square array - all kl for a given ij written to one record. No perm sym.
        unit_g_write_ccoo=get_free_unit()
        open(unit_g_write_ccoo,file='G_Dir_ccoo',status='unknown',form='unformatted',access='direct',action='write',recl=reclen_nmo_sq)
        unit_g_write_oocc=get_free_unit()
        open(unit_g_write_oocc,file='G_Dir_oocc',status='unknown',form='unformatted',access='direct',action='write',recl=reclen_ntmo_sq)

        rewind(unit_g_read)
        if(tDaltonFormat) then
            read(unit_g_read) maxlength
!            write(6,*) "Maxlength = ",maxlength
            allocate(Indices(4,MaxLength),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err 1")
            allocate(Buf(MaxLength),stat=ierr)
            if(ierr.ne.0) call stop_all(t_r,"Alloc err 1a")
        else
!            length = 1  !Consider it as 1 integral per 'buffer'
             MaxLength = 65536  !Actually, optimise by reading in 65536 integrals at a time
             allocate(Indices_long(4,MaxLength),stat=ierr)
             if(ierr.ne.0) call stop_all(t_r,"Alloc err 1")
             allocate(Buf(MaxLength),stat=ierr)
             if(ierr.ne.0) call stop_all(t_r,"Alloc err 1a")
        endif

        allocate(g_pairs(tnmo_sq,tntmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc err 2")
        allocate(g_pairs_oocc(tntmo_sq,tnmo),stat=ierr)
        if(ierr.ne.0) call stop_all(t_r,"Alloc err 2")

        write(6,'(F10.3,A)') 2.0_dp*real(tnmo_sq,dp)*real(tntmo_sq,dp)*8.0_dp*9.31322575D-10, &
            " Gb required to write out G_Dir_ccoo and G_Dir_oocc files."
        call flush(6)
        
!        open(24,file="FR-Full",status='unknown')

        !Loop over all i to store
        do i=1,tntmo

            !fill g_pairs with all < i' J' | K L >
            !fill g_pairs_oocc with all < i J | K' L' >
            g_pairs(:,:) = 0.0_dp
            g_pairs_oocc(:,:) = 0.0_dp

            rewind(unit_g_read)
            if(tDaltonFormat) then
                read(unit_g_read) maxlength
            elseif(.not.tReadBin) then
                read(unit_g_read,fci)
            endif
            do while(.true.)
                !Read all integrals for each ij pair
                !This is slow, and doesn't consider any special ordering of the orbitals
                if(tDaltonFormat) then
                    read(unit_g_read,iostat=ios) length,Indices(1:4,1:length),Buf(1:length)
                    if(ios.lt.0) call stop_all(t_r,"Error reading integrals - should not have reached eof")
                else
                    !Read from formatted FCIDUMP-type file. Again, we are converting from chemical notation.
                    if(tReadBin) then
                        do n=1,MaxLength
                            read(unit_g_read,iostat=ios) Buf(n),Indices_long(1:4,n)
                            if(ios.lt.0) exit
                        enddo
                    else
                        do n=1,MaxLength
                            read(unit_g_read,*,iostat=ios) Buf(n),Indices_long(1:4,n)
                            if(ios.lt.0) exit
                        enddo
                    endif
                    if(ios.lt.0) then
                        !The last buffer is incomplete.
                        length = n-1
                    else
                        length = MaxLength  !Full buffer
                    endif
                endif
                if(ios.gt.0) call stop_all(t_r,"Error reading integrals")
                if(length.le.0) exit    !Finished reading all integrals.

                do ibuf=1,length
                    !Indices are read in in chemical notation. We store in PHYSICAL
                    !Find pair index - alternatively, we could back construct the ij pair, and compare to this
                    if(tDaltonFormat) then
                        i_Ind=int(Indices(1,ibuf))
                        j_Ind=int(Indices(3,ibuf))
                        k_Ind=int(Indices(2,ibuf))
                        l_Ind=int(Indices(4,ibuf))
                    else
                        i_Ind=Indices_long(1,ibuf)
                        j_Ind=Indices_long(3,ibuf)
                        k_Ind=Indices_long(2,ibuf)
                        l_Ind=Indices_long(4,ibuf)
                        if((i_ind.eq.0).or.(j_Ind.eq.0).or.(k_ind.eq.0).or.(l_ind.eq.0)) cycle  !We only want 2e-integrals
                    endif

                    if((i_Ind.ne.i).and.(j_Ind.ne.i).and.(k_Ind.ne.i).and.(l_Ind.ne.i)) cycle

                    !At least one of the indices is equal to i
                    !Fill arrays, taking into account permutational sym, and size of spaces
                    if(i_Ind.eq.i) then
                        !<ij|kl>
                        if(max(k_Ind,l_Ind).le.tnmo) then
                            g_pairs(FindSqPairInd(k_Ind,l_Ind,tnmo),j_Ind) = Buf(ibuf)
                        endif
                        if(max(i_Ind,j_Ind).le.tnmo) then
                            g_pairs_oocc(FindSqPairInd(k_Ind,l_Ind,tntmo),j_Ind) = Buf(ibuf) 
                        endif
                        if(.not.tComplexOrbs) then
                            !<il|kj>
                            if(max(k_Ind,j_Ind).le.tnmo) then
                                g_pairs(FindSqPairInd(k_Ind,j_Ind,tnmo),l_Ind) = Buf(ibuf)  
                            endif
                            if(max(i_Ind,l_Ind).le.tnmo) then
                                g_pairs_oocc(FindSqPairInd(k_Ind,j_Ind,tntmo),l_Ind) = Buf(ibuf)  
                            endif
                        endif
                    endif

                    if(j_Ind.eq.i) then
                        !<ji|lk>
                        if(max(l_Ind,k_Ind).le.tnmo) then
                            g_pairs(FindSqPairInd(l_Ind,k_Ind,tnmo),i_Ind) = Buf(ibuf)  
                        endif
                        if(max(j_Ind,i_Ind).le.tnmo) then
                            g_pairs_oocc(FindSqPairInd(l_Ind,k_Ind,tntmo),i_Ind) = Buf(ibuf)  
                        endif
                        if(.not.tComplexOrbs) then
                            !<jk|li>
                            if(max(l_Ind,i_Ind).le.tnmo) then
                                g_pairs(FindSqPairInd(l_Ind,i_Ind,tnmo),k_Ind) = Buf(ibuf)  
                            endif
                            if(max(j_Ind,k_Ind).le.tnmo) then
                                g_pairs_oocc(FindSqPairInd(l_Ind,i_Ind,tntmo),k_Ind) = Buf(ibuf)  
                            endif
                        endif
                    endif

                    if(k_Ind.eq.i) then
                        !<kl|ij>
                        if(max(i_Ind,j_Ind).le.tnmo) then
                            g_pairs(FindSqPairInd(i_Ind,j_Ind,tnmo),l_Ind) = Buf(ibuf)  
                        endif
                        if(max(k_Ind,l_Ind).le.tnmo) then
                            g_pairs_oocc(FindSqPairInd(i_Ind,j_Ind,tntmo),l_Ind) = Buf(ibuf)  
                        endif
                        if(.not.tComplexOrbs) then
                            !<kj|il>
                            if(max(i_Ind,l_Ind).le.tnmo) then
                                g_pairs(FindSqPairInd(i_Ind,l_Ind,tnmo),j_Ind) = Buf(ibuf)  
                            endif
                            if(max(k_Ind,j_Ind).le.tnmo) then
                                g_pairs_oocc(FindSqPairInd(i_Ind,l_Ind,tntmo),j_Ind) = Buf(ibuf)  
                            endif
                        endif
                    endif

                    if(l_Ind.eq.i) then
                        !<lk|ji>
                        if(max(j_Ind,i_Ind).le.tnmo) then
                            g_pairs(FindSqPairInd(j_Ind,i_Ind,tnmo),k_Ind) = Buf(ibuf)  
                        endif
                        if(max(l_Ind,k_Ind).le.tnmo) then
                            g_pairs_oocc(FindSqPairInd(j_Ind,i_Ind,tntmo),k_Ind) = Buf(ibuf)  
                        endif
                        if(.not.tComplexOrbs) then
                            !<li|jk>
                            if(max(j_Ind,k_Ind).le.tnmo) then
                                g_pairs(FindSqPairInd(j_Ind,k_Ind,tnmo),i_Ind) = Buf(ibuf)  
                            endif
                            if(max(l_Ind,i_Ind).le.tnmo) then
                                g_pairs_oocc(FindSqPairInd(j_Ind,k_Ind,tntmo),i_Ind) = Buf(ibuf)  
                            endif
                        endif
                    endif

                enddo

            enddo


            !All {k,l},j for the given i are now stored. Write these out.

            do j=1,tntmo

                ijInd = FindSqPairInd(i,j,tntmo)
                write(unit_g_write_ccoo,rec=ijInd) (g_pairs(loop,j),loop=1,tnmo_sq)

                if(max(i,j).le.tnmo) then
                    ijInd = FindSqPairInd(i,j,tnmo)
                    write(unit_g_write_oocc,rec=ijInd) (g_pairs_oocc(loop,j),loop=1,tntmo_sq)
                endif

!                do k=1,tnmo
!                    do l=1,tnmo
!                        loop=loop+1
!                        write(24,"(4I4,F15.8)") I_Pair_ind,J_Pair_Ind,k,l,r_pairs(loop)
!                    enddo
!                enddo

            enddo
        enddo

        !All ij pairs have now been written out in direct access file.
        close(unit_g_read)
        close(unit_g_write_oocc)
        close(unit_g_write_ccoo)
!        close(24)
        if(tDaltonFormat) then
            deallocate(Indices)
        else
            deallocate(Indices_long)
        endif
        deallocate(g_pairs,g_pairs_oocc,Buf)

    end subroutine StoreGInts



!The maximum indices in a particular file, where the indices are returned in PHYSICAL NOTATION!
    SUBROUTINE FindMaxIndices(iUnit,Maxi,Maxj,Maxk,Maxl)
        IMPLICIT NONE
        INTEGER(2) , INTENT(OUT) :: Maxi,Maxj,Maxk,Maxl
        INTEGER , INTENT(IN) :: iUnit
        INTEGER(8) :: MaxLength,Length
        INTEGER(2) , ALLOCATABLE :: Indices(:,:)
        INTEGER :: i,iBuf

        rewind(iUnit)
        read(iUnit) MaxLength
        ALLOCATE(Indices(4,MaxLength))
        Maxi=0
        Maxj=0
        Maxk=0
        Maxl=0
        iBuf=0
        do while(.true.)
            iBuf=iBuf+1
!            WRITE(6,*) "Attempting to read buffer: ",iBuf
            read(iUnit,err=11) length,Indices(1:4,1:length)
            IF(length.gt.0) THEN
                do i=1,length
                    IF(Indices(1,i).gt.Maxi) Maxi=Indices(1,i)
                    IF(Indices(2,i).gt.Maxk) Maxk=Indices(2,i)
                    IF(Indices(3,i).gt.Maxj) Maxj=Indices(3,i)
                    IF(Indices(4,i).gt.Maxl) Maxl=Indices(4,i)
                enddo
            ELSE
                EXIT
            ENDIF
        enddo
        DEALLOCATE(Indices)
        REWIND(iUnit)
        RETURN
11      STOP 'Error when finding largest indices'
    END SUBROUTINE FindMaxIndices

!!Read in r integrals from 4traf file, and write them out to a direct access file,
!!with each ij pair in a different index. At the moment, no triangular indexing or
!!permutational symmetry at all.
!!ij pairs will run over (complete, complete) basis pairs (each record)
!!the kl pairs will run over (orbital, orbital) basis pairs (within record)
!!The 4traf file MO_R runs over < complete, complete | orbital orbital >.
!    subroutine StoreRInts() 
!        implicit none
!        integer :: unit_r_read,unit_r_write,ierr,rdump_unit,unit_r_write_oocc!,k,l,loop
!        integer(i2) , allocatable :: Indices(:,:)
!        real(dp) , allocatable :: Buf(:),r_pairs(:),r_pairs_oocc(:)
!        character(*), parameter :: t_r='StoreRInts'
!        integer :: i,ios,ibuf,i_Pair_Ind,j_Pair_Ind,ijInd,i_Ind,j_Ind,k_Ind,l_Ind
!        integer :: ijInd_oocc
!        integer(i8) :: Maxlength,length
!        integer(i2) :: MaxFilei,MaxFilej,MaxFilek,MaxFilel 
!
!        if(tRDUMP) then
!            rdump_unit=get_free_unit()
!            open(rdump_unit,file='RDUMP',status='unknown')
!        endif
!
!        !First, read in all primitive r integrals.
!        !We want the first two indices to be over entire space (+ CABS).
!        unit_r_read=get_free_unit()
!        open(unit_r_read,file='MO_F12',status='old',form='unformatted',access='sequential',action='read')
!
!        !For the r integrals, from 4traf, the indices should be (physical):
!        !i = complete
!        !j = complete
!        !k = orbital 
!        !l = orbital
!        CALL FindMaxIndices(unit_r_read,MaxFilei,MaxFilej,MaxFilek,MaxFilel)
!!        write(6,"(A,4I5)") "Index upper limits for the g file are (physical notation): ",MaxFilei,MaxFilej,MaxFilek,MaxFilel
!
!        !Integrals written out in a square array - all kl for a given p'q' written to one record. No perm sym.
!        !R_Dir_cc.oo writes out a seperate record for all pairs in the complete space. Each record holds all orbital:orbital pairs
!        unit_r_write=get_free_unit()
!        open(unit_r_write,file='R_Dir_cc.oo',status='unknown',form='unformatted',   &
!            access='direct',action='write',recl=reclen_nmo_sq)
!
!        !R_Dir_oo.cc writes out records for all orbital:orbital pairs, and each record has total orbital pairs integrals
!        unit_r_write_oocc=get_free_unit()
!        open(unit_r_write_oocc,file='R_Dir_oo.cc',status='unknown',form='unformatted',  &
!            access='direct',action='write',recl=reclen_ntmo_sq)
!
!        rewind(unit_r_read)
!        read(unit_r_read) maxlength
!!        write(6,*) "Maxlength = ",maxlength
!        allocate(Indices(4,MaxLength),stat=ierr)
!        if(ierr.ne.0) call stop_all(t_r,"Alloc err 1")
!        allocate(Buf(MaxLength),stat=ierr)
!        if(ierr.ne.0) call stop_all(t_r,"Alloc err 1a")
!
!        allocate(r_pairs(tnmo_sq),stat=ierr)
!        if(ierr.ne.0) call stop_all(t_r,"Alloc err 2")
!        allocate(r_pairs_oocc(tntmo_sq),stat=ierr)
!        if(ierr.ne.0) call stop_all(t_r,"Alloc err 2")
!        
!!        open(24,file="FR-Full",status='unknown')
!
!        !Loop over all pairs, since we want to store each ij pair in a direct access record.
!        do i_Pair_Ind=1,tntmo
!            do j_Pair_Ind=1,tntmo
!
!                ijInd=FindSqPairInd(i_Pair_Ind,j_Pair_Ind,tntmo)  !ijInd is the index of the given (ij) pair
!
!                if(max(i_Pair_Ind,j_Pair_Ind).le.tnmo) then
!                    !ijInd_oocc is the index of the given (ij) pair if ordering over orbtal basis
!                    ijInd_oocc=FindSqPairInd(i_Pair_Ind,j_Pair_Ind,tnmo)  
!                endif
!                
!                r_pairs(:)=0.0_dp
!                r_pairs_oocc(:)=0.0_dp
!                rewind(unit_r_read)
!                read(unit_r_read) maxlength
!                do while(.true.)
!                    !Read all integrals for each ij pair
!                    !This is slow, and doesn't consider any special ordering of the orbitals
!                    read(unit_r_read,iostat=ios) length,Indices(1:4,1:length),Buf(1:length)
!                    if(ios.ne.0) then
!                        call stop_all(t_r,"Error reading integrals")
!                    endif
!
!                    if(length.le.0) exit    !Finished reading all integrals.
!
!                    do ibuf=1,length
!                        !Indices are read in in chemical notation. We store in PHYSICAL
!                        !Find pair index - alternatively, we could back construct the ij pair, and compare to this
!                        i_Ind=int(Indices(1,ibuf),i4)
!                        j_Ind=int(Indices(3,ibuf),i4)
!                        k_Ind=int(Indices(2,ibuf),i4)
!                        l_Ind=int(Indices(4,ibuf),i4)
!
!                        if(tRDUMP) then
!                            if(ijInd.eq.1) write(rdump_unit,'(1X,G20.14,4I3)') Buf(ibuf),i_Ind,k_Ind,j_Ind,l_Ind
!                        endif
!!                        write(19,"(4I4,F15.8)") i_Ind,j_Ind,k_Ind,l_Ind,Buf(ibuf)
!
!                        !Permutational symmetry of 4traf integrals!
!                        call PutinRpairs_ccoo(i_Ind,j_Ind,k_Ind,l_Ind,Buf(ibuf),r_pairs,ijInd)
!
!                        if(max(i_Pair_Ind,j_Pair_Ind).le.tnmo) then
!                            !We only want to store orbital:orbital pairs for ij in the oocc file
!!                            write(6,*) "Storing pair: ",i_Pair_Ind,j_Pair_Ind,ijInd_oocc
!                            call PutinRPairs_oocc(i_Ind,j_Ind,k_Ind,l_Ind,Buf(ibuf),r_pairs_oocc,ijInd_oocc)
!                        endif
!
!                    enddo
!
!                enddo
!
!                !All kl for the given ij are now stored. Write these out.
!!                write(6,*) "ijInd: ",ijInd
!                write(unit_r_write,rec=ijInd) (r_pairs(i),i=1,tnmo_sq)
!                        
!                if(max(i_Pair_Ind,j_Pair_Ind).le.tnmo) then
!                    !Ensure that we are dealing with 
!!                    write(6,*) "*****",ijInd_oocc,i_Pair_Ind,j_Pair_Ind
!!                    if(ijInd_oocc.eq.3) write(6,*) r_pairs_oocc(:)
!                    write(unit_r_write_oocc,rec=ijInd_oocc) (r_pairs_oocc(i),i=1,tntmo_sq)
!                endif
!
!!                loop=0
!!                do k=1,tnmo
!!                    do l=1,tnmo
!!                        loop=loop+1
!!                        write(24,"(4I4,F15.8)") I_Pair_ind,J_Pair_Ind,k,l,r_pairs(loop)
!!                    enddo
!!                enddo
!
!            enddo
!        enddo
!
!        !All ij pairs have now been written out in direct access file.
!        close(unit_r_read)
!        close(unit_r_write)
!        close(unit_r_write_oocc)
!        if(tRDUMP) then
!            close(rdump_unit)
!        endif
!!        close(24)
!        deallocate(Indices,Buf,r_pairs,r_pairs_oocc)
!
!    end subroutine StoreRInts
!    
!!Read in g integrals from 4traf file, and write them out to a direct access file,
!!with each ij pair in a different index. At the moment, no triangular indexing or
!!permutational symmetry at all.
!!ij pairs will run over (complete, complete) basis pairs (each record)
!!the kl pairs will run over (orbital, orbital) basis pairs (within record)
!!The 4traf file MO_R runs over < complete, complete | orbital orbital >.
!    subroutine StoreGInts() 
!        integer :: unit_g_read,unit_g_write_ccoo,ierr,unit_g_write_oocc!,k,l,loop
!        integer(i2) , allocatable :: Indices(:,:)
!        real(dp) , allocatable :: Buf(:),g_pairs(:),g_pairs_oocc(:)
!        character(*), parameter :: t_r='StoreGInts'
!        integer :: i,ios,ibuf,i_Pair_Ind,j_Pair_Ind,ijInd,i_Ind,j_Ind,k_Ind,l_Ind,ijInd_oocc
!        integer(i8) :: Maxlength,length
!        integer(i2) :: MaxFilei,MaxFilej,MaxFilek,MaxFilel 
!
!        !First, read in all primitive g integrals.
!        !We want the first two indices to be over entire space (+ CABS).
!        unit_g_read=get_free_unit()
!        open(unit_g_read,file='MO_G',status='old',form='unformatted',access='sequential',action='read')
!
!        !For the g integrals, from 4traf, the indices should be (physical):
!        !i = complete
!        !j = complete
!        !k = complete
!        !l = orbital
!        CALL FindMaxIndices(unit_g_read,MaxFilei,MaxFilej,MaxFilek,MaxFilel)
!!        write(6,"(A,4I5)") "Index upper limits for the g file are (physical notation): ",MaxFilei,MaxFilej,MaxFilek,MaxFilel
!
!        !Integrals written out in a square array - all kl for a given ij written to one record. No perm sym.
!        unit_g_write_ccoo=get_free_unit()
!        open(unit_g_write_ccoo,file='G_Dir_ccoo',status='unknown',form='unformatted',access='direct',action='write',recl=reclen_nmo_sq)
!        unit_g_write_oocc=get_free_unit()
!        open(unit_g_write_oocc,file='G_Dir_oocc',status='unknown',form='unformatted',access='direct',action='write',recl=reclen_ntmo_sq)
!
!        rewind(unit_g_read)
!        read(unit_g_read) maxlength
!!        write(6,*) "Maxlength = ",maxlength
!        allocate(Indices(4,MaxLength),stat=ierr)
!        if(ierr.ne.0) call stop_all(t_r,"Alloc err 1")
!        allocate(Buf(MaxLength),stat=ierr)
!        if(ierr.ne.0) call stop_all(t_r,"Alloc err 1a")
!
!        allocate(g_pairs(tnmo_sq),stat=ierr)
!        if(ierr.ne.0) call stop_all(t_r,"Alloc err 2")
!        allocate(g_pairs_oocc(tntmo_sq),stat=ierr)
!        if(ierr.ne.0) call stop_all(t_r,"Alloc err 2")
!        
!!        open(24,file="FR-Full",status='unknown')
!
!        !Loop over all pairs, since we want to store each ij pair in a direct access record.
!        do i_Pair_Ind=1,tntmo
!            do j_Pair_Ind=1,tntmo
!
!                ijInd=FindSqPairInd(i_Pair_Ind,j_Pair_Ind,tntmo)  !ijInd is the index of the given (ij) pair
!
!                if(max(i_Pair_Ind,j_Pair_Ind).le.tnmo) then
!                    ijInd_oocc = FindSqPairInd(i_Pair_Ind,j_Pair_Ind,tnmo)  !IjInd_oocc is the pair for storage over occupied pairs
!                else
!                    ijInd_oocc = 0  !To get array bounds errors
!                endif
!                
!                g_pairs(:) = 0.0_dp
!                g_pairs_oocc(:) = 0.0_dp
!                rewind(unit_g_read)
!                read(unit_g_read) maxlength
!                do while(.true.)
!                    !Read all integrals for each ij pair
!                    !This is slow, and doesn't consider any special ordering of the orbitals
!                    read(unit_g_read,iostat=ios) length,Indices(1:4,1:length),Buf(1:length)
!                    if(ios.ne.0) then
!                        call stop_all(t_r,"Error reading integrals")
!                    endif
!
!                    if(length.le.0) exit    !Finished reading all integrals.
!
!                    do ibuf=1,length
!                        !Indices are read in in chemical notation. We store in PHYSICAL
!                        !Find pair index - alternatively, we could back construct the ij pair, and compare to this
!                        i_Ind=int(Indices(1,ibuf),i4)
!                        j_Ind=int(Indices(3,ibuf),i4)
!                        k_Ind=int(Indices(2,ibuf),i4)
!                        l_Ind=int(Indices(4,ibuf),i4)
!
!!                        write(19,"(4I4,F15.8)") i_Ind,j_Ind,k_Ind,l_Ind,Buf(ibuf)
!
!!Block for the (g_cc)^oo integrals ***********************************************************
!                        !Permutational symmetry of 4traf integrals!
!                        if(FindSqPairInd(i_Ind,j_Ind,tntmo).eq.ijInd) then
!                            !<ij|kl>
!                            if((k_Ind.le.tnmo).and.(l_Ind.le.tnmo)) then
!                                g_pairs(FindSqPairInd(k_Ind,l_Ind,tnmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(FindSqPairInd(j_Ind,i_Ind,tntmo).eq.ijInd) then
!                            !<ji|lk>
!                            if((l_Ind.le.tnmo).and.(k_Ind.le.tnmo)) then
!                                g_pairs(FindSqPairInd(l_Ind,k_Ind,tnmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(FindSqPairInd(k_Ind,l_Ind,tntmo).eq.ijInd) then
!                            !<kl|ij>
!                            if((i_Ind.le.tnmo).and.(j_Ind.le.tnmo)) then
!                                g_pairs(FindSqPairInd(i_Ind,j_Ind,tnmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(FindSqPairInd(l_Ind,k_Ind,tntmo).eq.ijInd) then
!                            !<lk|ji>
!                            if((j_Ind.le.tnmo).and.(i_Ind.le.tnmo)) then
!                                g_pairs(FindSqPairInd(j_Ind,i_Ind,tnmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(FindSqPairInd(k_Ind,j_Ind,tntmo).eq.ijInd) then
!                            !<kj|il>
!                            if((i_Ind.le.tnmo).and.(l_Ind.le.tnmo)) then
!                                g_pairs(FindSqPairInd(i_Ind,l_Ind,tnmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(FindSqPairInd(l_Ind,i_Ind,tntmo).eq.ijInd) then
!                            !<li|jk>
!                            if((j_Ind.le.tnmo).and.(k_Ind.le.tnmo)) then
!                                g_pairs(FindSqPairInd(j_Ind,k_Ind,tnmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(FindSqPairInd(i_Ind,l_Ind,tntmo).eq.ijInd) then
!                            !<il|kj>
!                            if((k_Ind.le.tnmo).and.(j_Ind.le.tnmo)) then
!                                g_pairs(FindSqPairInd(k_Ind,j_Ind,tnmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(FindSqPairInd(j_Ind,k_Ind,tntmo).eq.ijInd) then
!                            !<jk|li>
!                            if((l_Ind.le.tnmo).and.(i_Ind.le.tnmo)) then
!                                g_pairs(FindSqPairInd(l_Ind,i_Ind,tnmo))=Buf(ibuf)
!                            endif
!                        endif 
!
!!Block for the (g_oo)^cc integrals ***********************************************************
!                        if(max(i_Ind,j_Ind).le.tnmo) then
!                            if(FindSqPairInd(i_Ind,j_Ind,tnmo).eq.ijInd_oocc) then
!                                !<ij|kl>
!                                g_pairs_oocc(FindSqPairInd(k_Ind,l_Ind,tntmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(max(j_Ind,i_Ind).le.tnmo) then
!                            if(FindSqPairInd(j_Ind,i_Ind,tnmo).eq.ijInd_oocc) then
!                                !<ji|lk>
!                                g_pairs_oocc(FindSqPairInd(l_Ind,k_Ind,tntmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(max(k_Ind,l_Ind).le.tnmo) then
!                            if(FindSqPairInd(k_Ind,l_Ind,tnmo).eq.ijInd_oocc) then
!                                !<kl|ij>
!                                g_pairs_oocc(FindSqPairInd(i_Ind,j_Ind,tntmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(max(l_Ind,k_Ind).le.tnmo) then
!                            if(FindSqPairInd(l_Ind,k_Ind,tnmo).eq.ijInd_oocc) then
!                                !<lk|ji>
!                                g_pairs_oocc(FindSqPairInd(j_Ind,i_Ind,tntmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(max(k_Ind,j_Ind).le.tnmo) then
!                            if(FindSqPairInd(k_Ind,j_Ind,tnmo).eq.ijInd_oocc) then
!                                !<kj|il>
!                                g_pairs_oocc(FindSqPairInd(i_Ind,l_Ind,tntmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(max(l_Ind,i_Ind).le.tnmo) then
!                            if(FindSqPairInd(l_Ind,i_Ind,tnmo).eq.ijInd_oocc) then
!                                !<li|jk>
!                                g_pairs_oocc(FindSqPairInd(j_Ind,k_Ind,tntmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(max(i_Ind,l_Ind).le.tnmo) then
!                            if(FindSqPairInd(i_Ind,l_Ind,tnmo).eq.ijInd_oocc) then
!                                !<il|kj>
!                                g_pairs_oocc(FindSqPairInd(k_Ind,j_Ind,tntmo))=Buf(ibuf)
!                            endif
!                        endif
!                        if(max(j_Ind,k_Ind).le.tnmo) then
!                            if(FindSqPairInd(j_Ind,k_Ind,tnmo).eq.ijInd_oocc) then
!                                !<jk|li>
!                                g_pairs_oocc(FindSqPairInd(l_Ind,i_Ind,tntmo))=Buf(ibuf)
!                            endif
!                        endif
!
!                    enddo
!
!                enddo
!
!                !All kl for the given ij are now stored. Write these out.
!!                write(6,*) "ijInd: ",ijInd
!                write(unit_g_write_ccoo,rec=ijInd) (g_pairs(i),i=1,tnmo_sq)
!                
!                if(max(i_Pair_Ind,j_Pair_Ind).le.tnmo) then
!                    write(unit_g_write_oocc,rec=ijInd_oocc) (g_pairs_oocc(i),i=1,tntmo_sq)
!                endif
!
!!                loop=0
!!                do k=1,tnmo
!!                    do l=1,tnmo
!!                        loop=loop+1
!!                        write(24,"(4I4,F15.8)") I_Pair_ind,J_Pair_Ind,k,l,r_pairs(loop)
!!                    enddo
!!                enddo
!
!            enddo
!        enddo
!
!        !All ij pairs have now been written out in direct access file.
!        close(unit_g_read)
!        close(unit_g_write_oocc)
!        close(unit_g_write_ccoo)
!!        close(24)
!        deallocate(Indices,Buf,g_pairs,g_pairs_oocc)
!
!    end subroutine StoreGInts


end module int_management
