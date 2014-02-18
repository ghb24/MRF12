module dalt_funcs

!This module contains routines which are specific to reading in and manipulation dalton output files.

implicit none

contains

   subroutine find_label(a,lu,luerr)               
      use errors, only: stop_all
!                                                                      
!  adapted from hjaaj (DALTON) and mollab (4traf) 
!                                                                      
!  Purpose:                                                            
!     Search for MOLECULE labels on file LU                            
!                                                                      
      character*8, intent(in) :: a      !Label to search for.
      integer, intent(in) :: lu         !Unit to read from.
      integer, intent(inout) :: luerr   !Unit to write messages, but also returns as error code of routine.
      character*8 :: b(4), c
      data c/'********'/    
                                           
    1 read (lu,end=3,err=6) b                                          
      if (b(1).ne.c) goto 1                                           
      if (b(4).ne.a) goto 1                                           
      if (luerr.lt.0) luerr = 0                                        
      return                                                           
                                                                       
    3 if (luerr.lt.0) then                                             
         luerr = -1                                                    
         return                                                        
      else                                                             
         write(luerr,4)a,lu                                            
         write(luerr,*)'Molecule label not found on file'              
         stop                                                          
          call stop_all('find_label','molecule label not found on file')     
      end if                                                           
    4 format(/' MOLECULE label ',A8,' not found on unit',I4)                                 
                                                                       
    6 if (luerr.lt.0) then                                             
         luerr = -2                                                    
         return                                                        
      else                                                             
         write (luerr,7) lu,a                                          
         write(luerr,*)'error reading file'                            
      end if                                                           
    7 format(/' error reading unit',I4,/T22,'when searching for label ',A8)                      

      end subroutine find_label 

end module dalt_funcs
