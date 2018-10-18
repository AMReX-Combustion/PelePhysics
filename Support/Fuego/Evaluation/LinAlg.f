      subroutine dscal (n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if ( n.le.0 .or. incx.le.0 ) return
      if (incx.eq.1) goto 20
c
c     code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        dx(i) = da*dx(i)
      enddo
      return
c
c     code for increment equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        dx(i) = da*dx(i)
      enddo
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
      enddo
      end

      subroutine vdaxpy (n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if (n.le.0) return
      if (da .eq. 0.0d0) return
      if (incx.eq.1.and.incy.eq.1) goto 20
c
c     code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
      enddo
      return
c
c     code for both increments equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,4)
      if ( m .eq. 0 ) goto 40
      do i = 1,m
        dy(i) = dy(i) + da*dx(i)
      enddo
      if ( n .lt. 4 ) return
   40 mp1 = m + 1
      do i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
      enddo
      end


      subroutine dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if (n.le.0) return
      if (incx.eq.1.and.incy.eq.1) goto 20
c
c     code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
      enddo
      return
c
c     code for both increments equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) goto 40
      do i = 1,m
        dy(i) = dx(i)
      enddo
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
      enddo
      end


      double precision function vddot (n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      vddot = 0.0d0
      dtemp = 0.0d0
      if (n.le.0) return
      if (incx.eq.1.and.incy.eq.1) goto 20
c
c     code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      enddo
      vddot = dtemp
      return
c
c     code for both increments equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,5)
      if ( m .eq. 0 ) goto 40
      do i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
      enddo
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
      enddo
   60 vddot = dtemp
      end

      double precision function ddot (n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if (n.le.0) return
      if (incx.eq.1.and.incy.eq.1) goto 20
c
c     code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if (incx.lt.0) ix = (-n+1)*incx + 1
      if (incy.lt.0) iy = (-n+1)*incy + 1
      do i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
      enddo
      ddot = dtemp
      return
c
c     code for both increments equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,5)
      if ( m .eq. 0 ) goto 40
      do i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
      enddo
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
      enddo
   60 ddot = dtemp
      end


      subroutine dgesl (a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(*),job
      double precision a(lda,n),b(*)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas vdaxpy,vddot
c
c     internal variables
c
      double precision vddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) goto 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) goto 30
         do k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) goto 10
               b(l) = b(k)
               b(k) = t
   10       continue
c            call vdaxpy(n-k,t,a(k+1,k),1,b(k+1),1)
            call ccse_daxpy(n-k,t,a(k+1,k),b(k+1))
         enddo
   30    continue
c
c        now solve  u*x = y
c
         do kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
c            call vdaxpy(k-1,t,a(1,k),1,b(1),1)
            call ccse_daxpy(k-1,t,a(1,k),b(1))
         enddo
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do k = 1, n
            t = vddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
         enddo
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) goto 90
         do kb = 1, nm1
            k = n - kb
            b(k) = b(k) + vddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
         enddo
   90    continue
  100 continue
      end

      subroutine dgbsl (abd,lda,n,ml,mu,ipvt,b,job)
      integer lda,n,ml,mu,ipvt(*),job
      double precision abd(lda,*),b(*)
c
c     dgbsl solves the double precision band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgbco or dgbfa.
c
c     on entry
c
c        abd     double precision(lda, n)
c                the output from dgbco or dgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from dgbco or dgbfa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgbco has set rcond .gt. 0.0
c        or dgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas vdaxpy,vddot
c     fortran min0
c
c     internal variables
c
      double precision vddot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) goto 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) goto 30
         if (nm1 .lt. 1) goto 30
            do k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) goto 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
c               call vdaxpy(lm,t,abd(m+1,k),1,b(k+1),1)
               call ccse_daxpy(lm,t,abd(m+1,k),b(k+1))
            enddo
   30    continue
c
c        now solve  u*x = y
c
         do kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
c            call vdaxpy(lm,t,abd(la,k),1,b(lb),1)
            call ccse_daxpy(lm,t,abd(la,k),b(lb))
         enddo
      goto 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = vddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
         enddo
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) goto 90
         if (nm1 .lt. 1) goto 90
            do kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + vddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) goto 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
            enddo
   90    continue
  100 continue
      end

      subroutine dgbfa (abd,lda,n,ml,mu,ipvt,info)
      integer lda,n,ml,mu,ipvt(*),info
      double precision abd(lda,*)
c
c     dgbfa factors a double precision band matrix by elimination.
c
c     dgbfa is usually called by dgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     double precision(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgbsl will divide by zero if
c                     called.  use  rcond  in dgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas vdaxpy,dscal,idamax
c     fortran max0,min0
c
c     internal variables
c
      double precision t
      integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) goto 30
      do jz = j0, j1
         i0 = m + 1 - jz
         do i = i0, ml
            abd(i,jz) = 0.0d0
         enddo
      enddo
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) goto 130
      do k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) goto 50
         if (ml .lt. 1) goto 50
            do i = 1, ml
               abd(i,jz) = 0.0d0
            enddo
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0d0) goto 100
c
c           interchange if necessary
c
            if (l .eq. m) goto 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -1.0d0/abd(m,k)
            call dscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) goto 90
            do j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) goto 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
c               call vdaxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
               call ccse_daxpy(lm,t,abd(m+1,k),abd(mm+1,j))
            enddo
   90       continue
         goto 110
  100    continue
            info = k
  110    continue
      enddo
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      end


      subroutine dgefa (a,lda,n,ipvt,info)
      integer lda,n,ipvt(*),info
      double precision a(lda,*)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas vdaxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) goto 70
      do k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) goto 40
c
c           interchange if necessary
c
            if (l .eq. k) goto 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do j = kp1, n
               t = a(l,j)
               if (l .eq. k) goto 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
c               call vdaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
               call ccse_daxpy(n-k,t,a(k+1,k),a(k+1,j))
            enddo
         goto 50
   40    continue
            info = k
   50    continue
      enddo
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      end


      integer function idamax (n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if ( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if (n.eq.1) return
      if (incx.eq.1) goto 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do i = 2,n
         if(dabs(dx(ix)).le.dmax) goto 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
      enddo
      return
c
c     code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) goto 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      end


      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  Purpose
*  =======
*
*     constant times a vector plus a vector.
*     uses unrolled loops for increments equal to one.
*     jack dongarra, linpack, 3/11/78.
*     modified 12/3/93, array(1) declarations changed to array(*)
*
*
*     .. Local Scalars ..
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) GO TO 20
*
*        code for unequal increments or equal increments
*          not equal to 1
*
      IX = 1
      IY = 1
      IF (INCX.LT.0) IX = (-N+1)*INCX + 1
      IF (INCY.LT.0) IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
   20 M = MOD(N,4)
      IF (M.EQ.0) GO TO 40
      DO 30 I = 1,M
          DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF (N.LT.4) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
          DY(I) = DY(I) + DA*DX(I)
          DY(I+1) = DY(I+1) + DA*DX(I+1)
          DY(I+2) = DY(I+2) + DA*DX(I+2)
          DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
      RETURN
      END

      subroutine ccse_daxpy (n,da,dx,dy)
      double precision dx(*),dy(*),da
      integer n, i
      do i=1,n
         dy(i) = dy(i) + da*dx(i)
      enddo
      end
