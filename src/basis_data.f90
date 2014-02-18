!Storage of global variables relating to the basis functions.
module basis_data

use const

implicit none
save

real(dp) :: pot_nuc !Nuclear potential energy
real(dp) :: ref_ene !Reference energy
real(dp) :: EMP2    !MP2 correlation energy - calculated, not read in
real(dp) :: VarE    !Variational total energy

integer :: nsym     !Number of abelian symmetry elements in D2h or subgroup.

!Symmetry dependant arrays (i.e. will go from 1:nsym to store orbital information per symmetry group)
!These arrays all start with a 'n'
!Frozen orbitals
integer, allocatable :: nfrz(:)     !number of orbitals frozen
!Occupied Orbitals
integer, allocatable :: nocc(:)     !number of orbitals occupied in HF 
!Orbital basis
integer, allocatable :: nmo(:)      !number of functions to describe orbitals in orbital basis
integer, allocatable :: nbf(:)      !number of orbitals in AO orbital basis (This
                                    !will generally be the same as nmo unless linear dependancies)
!Auxiliary CABS space
integer, allocatable :: nxmo(:)     !number of MO's in CABS space 
integer, allocatable :: nxbf(:)     !number of AO orbitals in CABS space (will be same as nxmo unless linear dependencies)
!Total space                                    
integer, allocatable :: ntbf(:)     !number of orbitals in complete space
integer, allocatable :: ntmo(:)     !number of MO orbitals in complete space

integer, allocatable :: StartSym_nmo(:) !Starting index of the orbitals with symmetry i
integer, allocatable :: EndSym_nmo(:)   !Ending index of the orbitals with symmetry i
integer, allocatable :: StartSym_nxmo(:) !Starting index of the aux orbitals with symmetry i
integer, allocatable :: EndSym_nxmo(:)   !Ending index of the aux orbitals with symmetry i

integer, allocatable :: OrbPhaseFlip(:) !For NOCI, if this is -1, then flip the sign of the contribution to this orbital in the RDMs.

!Total basis function numbers (i.e. sum over all symmetries)
!These integers all start with a 't' before the symmetry dependant array name.
!Number of electrons
integer :: NEl
!FrozenOrbitals
integer :: tnfrz
!Occupied Orbitals
integer :: tnocc    !number of orbitals occupied in HF
!Orbital basis
integer :: tnmo     !number of orbitals in orbital basis
integer :: tnbf     !number of AOs in orbital basis (This
                    !will generally be the same as nmo unless linear dependancies)
!Auxiliary CABS space
integer :: tnxmo     !number of orbitals in CABS space. 
integer :: tnxbf     !number of AO orbitals in CABS space (same as !tnxmo unless linear dependencies)
!Total space                                    
integer :: tntbf     !number of orbitals in complete space
integer :: tntmo     !number of MOs in complete space
!Pair space
integer :: tntmo_sq  !Total pairs of orbitals in complete space (no perm sym)
integer :: tnmo_sq   !Total pairs of orbitals in orbtial space (no perm sym)
integer :: taux_nmo  !Total pairs of made from CABS x orbital space (no perm sym)
integer :: tntmo_nmo    !Total pairs of orbitals from complete x orbital space.
integer :: tntmo_nxmo   !Total pairs of orbitals from complete x aux space.

!Record lengths
integer :: reclen_nmo_sq
integer :: reclen_ntmo_sq
                         
!Orbital fock energies in the orbital basis
real(dp), allocatable :: OrbEnergies(:)
!Symmetry characters of the orbitals in the orbital basis
integer, allocatable :: OrbSyms(:)
!Symmetry characters of the orbitals in the CABS basis
integer, allocatable :: AuxOrbSyms(:)

!List the indices of the HF occupied orbitals
integer, allocatable :: OccOrbs(:)
integer, allocatable :: OccOrbNum(:)
integer, allocatable :: FrzOrbs(:)
integer, allocatable :: OrbMapping(:) !If read in orbital i, then store it in OrbMapping(i)
!Logical list, which is true if the orbital is occupied in the HF det.
logical, allocatable :: IsOrbOcc(:)
logical, allocatable :: IsOrbFrz(:)

!Fock matrix in orbital:orbital space
real(dp) , allocatable :: FockOrb(:,:)
!Fock matrix in orbital:CABS space
!Orbital is first index, CABS second
real(dp) , allocatable :: FockOrbCABS(:,:)
!Fock matrix in CABS,CABS space
real(dp) , allocatable :: FockCABS(:,:) 

!Generalised Fock matrix in orbital:orbital space
real(dp) , allocatable :: GenFockOrb(:,:)
!Generalised Fock matrix in orbital:CABS space
!Orbital is first index, CABS second
real(dp) , allocatable :: GenFockOrbCABS(:,:)
!Generalised Fock matrix in CABS,CABS space
real(dp) , allocatable :: GenFockCABS(:,:) 
!Generalised Fock matrix in complete:orbital space.
real(dp) , allocatable :: GenFockCompOrb(:,:)
!Generalised Fock matrix in complete:complete space.
real(dp) , allocatable :: GenFockComp(:,:)

!Generalised Exchange matrix in complete:complete space.
real(dp) , allocatable :: GenExch(:,:)
!Generalised fock + exchange matrix in complete:orbital space.
real(dp) , allocatable :: GenFpK_CompOrb(:,:)

!one-electron matrix elements over complete space
real(dp) , allocatable :: tmat(:,:)
!Just kinetic energy elements over complete space
real(dp) , allocatable :: Kinmat(:,:)

!OneRDM over orbital space
real(dp) , allocatable :: OneRDM(:,:)

!Lambda_Bar over orbital space
!Eventually, this wants to be LB_q^s = f_t^u L_uq^ts where L is the 2-cumulant
real(dp) , allocatable :: Lambda_Bar(:,:)

!This array is the self-styled X matrix used for the MR-orbital relaxation.
!Its precise form depends on tDyall
real(dp) , allocatable :: RelaxXDyall(:,:)
real(dp) , allocatable :: RelaxXFock(:,:)

end module basis_data
