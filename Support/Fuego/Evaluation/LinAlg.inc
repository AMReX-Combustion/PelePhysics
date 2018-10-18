!
! F90 interface definitions for (selected) routines in LinAlg.f
!

interface
  subroutine dgefa(a, lda, n, ipvt, info)
    integer:: lda, n
    integer :: info, ipvt(*)
    double precision :: a(lda, *)
  end subroutine dgefa
  subroutine dgesl(a, lda, n, ipvt, b, job)
    integer :: lda, n, job, ipvt(*)
    double precision :: a(lda,n), b(*)
  end subroutine dgesl
end interface
