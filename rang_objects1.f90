MODULE parameter
  IMPLICIT NONE
!  nnnp     -- maximum number of open shells
  INTEGER, PARAMETER :: nnnp = 600
!  mX      -- maximum number of open shells
  INTEGER, PARAMETER :: MX   = 4
! recsize  -- record size for linked lists
  INTEGER, PARAMETER :: recsize=10
END MODULE parameter

MODULE case_variables   ! needed for dimensions
  IMPLICIT NONE
! Z    - atomic number
! A    - mass number
  REAL(kind=8) :: Z, A
! ncore  -- number of core orbita;s
! npeel  -- number of peel orbitals
! nwf    -- total number of orbitals, nwf = ncore +nwf
  INTEGER  :: ncore, npeel, nwf
! ncsf   -- total number of CSFs
! nblk   -- number of blocks of given J and Parity
  Integer:: ncsf, nblk
END MODULE case_variables

MODULE  orb_def
  USE parameter
  USE case_variables,  ONLY: nwf
  IMPLICIT NONE 

  TYPE  orbital
     Character(LEN=4) :: el                    ! symbol for the orbital such as ‘1s’
     Integer(kind=1)  :: n, l, jj, kappa
  END TYPE orbital

  TYPE :: radial
     integer(kind=4)  :: npt
     REAL (kind=8) :: e, az, scf
     REAL(kind=8), DIMENSION(nnnp) :: P, Q
  END TYPE radial

    TYPE(orbital),  DIMENSION(:), allocatable :: orbs
    TYPE(radial),   DIMENSION(:), allocatable :: rads

CONTAINS
    SUBROUTINe  alloc_orbs
        allocate(orbs(nwf))
    END subroutine alloc_orbs
   
    SUBROUTINe  alloc_rads
        allocate(rads(nwf))
    END subroutine alloc_rads
END MODULE orb_def


MODULE csf_def
  USE parameter

  TYPE CSF
!  m         --number of subshells
  INTEGER(kind=1)     :: m
!  iel       --index of orbital
!  qa        --occuation number
!  J_sub     -- J for subshell (in 2J+1)representation
!  q         -- quasi-spin, Q
!  qM        -- Q_M (or other quantum number)
!  jcup      -- coupling left to right in (2J++1) representation
  INTEGER(kind=1),DIMENSION(mX) :: iel, qa, J_sub, q, qM, jcup
  END Type CSF

 
  Type csf_rec
     INTEGER :: N
     TYPE(csf), dimension(recsize) :: c_set
     TYPE(csf_rec), POINTER :: next
  END TYPE csf_rec

CONTAINS
  SUBROUTINE init(head, tail)
     ! Dummy arguments
     TYPE (csf_rec), POINTER :: head, tail
     NULLIFY(head, tail)
  END SUBROUTINE init

  SUBROUTINE add(new, head, tail)
     !Add a new rec to the list
     !Dummy arguments
     TYPE(csf_rec), POINTER :: new, head, tail
     
     ! Check if the list is empty
     IF (ASSOCIATED(head)) THEN
        ! list is NOT empty`
        tail%next =>new    ! (no need to allocate?)
        NULLIFY(new%next)
        tail => new
     ELSE 
        ! list is empty setup the first rec
        head => new   
        tail => new
        NULLIFY(tail%next)
     END IF
  END SUBROUTINE add

  SUBROUTINE write(head) ! output the contents of the list
     ! Dummy argument
     TYPE(CSF_rec), POINTER, INTENT(in):: head
     ! Local variable
     TYPE(CSF_rec), POINTER :: ptr
     INTEGER :: j 
     !Check if list is empty
     IF (.NOT. ASSOCIATED(head)) THEN
        Print *, 'List is empty'
     ELSE
        ! Set local pointer to head
        ptr => head
        DO
           Do j = 1,ptr%N
             print *, ptr%c_set(j)
           End Do
           ptr => ptr%next
           If (.NOT. Associated(ptr)) EXIT
        END DO
     END IF
  END SUBROUTINE
END MODULE CSF_def

MODULE ang_linked_list
   Use parameter
   IMPLICIT  NONE

   TYPE ang_rec
      INTEGER :: N  ! The number of items in the rec
      REAL(kind=8),     dimension(recsize) :: coeff
      INTEGER(kind=8),  dimension(recsize) :: intgrl
      TYPE(ang_rec), POINTER :: next
   END TYPE ang_rec

CONTAINS
   SUBROUTINE init(head, tail)
      ! Dummy arguments
      TYPE(ang_rec), POINTER :: head, tail
      NULLIFY(head, tail)
   END SUBROUTINE init

   SUBROUTINE add(new, head, tail)
      !Add a new rec to the list
      !Dummy arguments
      TYPE(ang_rec), POINTER :: new, head, tail
      
      ! Check if the list is empty
      IF (ASSOCIATED(head)) THEN
         ! list is NOT empty`
         tail%next =>new    ! (no need to allocate?)
         NULLIFY(new%next)
         tail => new
      ELSE 
         ! list is empty setup the first rec
         head => new   
         tail => new
         NULLIFY(tail%next)
      END IF
   END SUBROUTINE add

   SUBROUTINE write(head) ! output the contents of the list
      ! Dummy argument
      TYPE(ang_rec), POINTER, INTENT(in):: head
      ! Local variable
      TYPE(ang_rec), POINTER :: ptr
      INTEGER :: j 
      !Check if list is empty
      IF (.NOT. ASSOCIATED(head)) THEN
         Print *, 'List is empty'
      ELSE
         ! Set local pointer to head
         ptr => head
         DO
            Do j = 1,ptr%N
              print *, ptr%coeff(j), ptr%intgrl(j)
            End Do
            ptr => ptr%next
            If (.NOT. Associated(ptr)) EXIT
         END DO
      END IF
   END SUBROUTINE
END MODULE ang_linked_list



Program  test_CSF

  USE CSF_def
  IMPLICIT NONE

  ! generate the linked list
  TYPE(csf_rec), POINTER :: clist_head,clist_tail 
  TYPE(csf_rec), POINTER :: clist
  INTEGER :: i, j, nsize=1, ierr, iend=5
 
  CALL init(clist_head, clist_tail)
  i=1       ! counter for angular data
  Do
     allocate(clist, STAT=ierr) 
     If (ierr .ne. 0) then 

        print *, 'Something wrong with alloc'
        stop
     end if

     clist% N = 0
     Do j = 1, nsize
        clist % c_set(j)%iel(1:4)= (/ 1, 2, 3, 4 /)
        print *, clist%c_set(1) % iel
        clist % N = clist % N + 1
        i=i+1
        if ( i .gt. iend) exit
     END DO 
     print *, 'Calling add when i=',i
     call add(clist, clist_head, clist_tail)
     if ( i .gt. iend) exit
  END DO
  
  print *, 'Finished generating data'
  call write(clist_head)

  CALL test_ang_data

CONTAINS 
   SUBROUTINE test_ang_data
      USE ang_linked_list
      IMPLICIT NONE
      ! generate the linked list
      TYPE(ang_rec), POINTER :: ang_head, ang_tail
      TYPE(ang_rec), POINTER :: my_data
      INTEGER :: i, j, nsize=100, ierr, iend=10
     
      CALL init(ang_head, ang_tail)
      i=1       ! counter for angular data
      Do
         allocate(my_data, STAT=ierr) 
         If (ierr .ne. 0) then
            print *, 'Something wrong with alloc'
            stop
         end if
  
         my_data% N = 0
         Do j = 1, nsize
            my_data % coeff(j) = i*1.d0
            my_data % intgrl(j) = i
            my_data % N = my_data % N + 1
            i=i+1
            if ( i .gt. iend) exit
         END DO 
         print *, 'Calling add when i=',i
         call add(my_data, ang_head, ang_tail)
         if ( i .gt. iend) exit
      END DO
      
      print *, 'Finished generating data'
      call write(ang_head)

     END SUBROUTINE test_ang_data

END program test_CSF

