c
c     Numerical Analysis:
c     Mathematics of Scientific Computing
c     Third Edition
c     D.R. Kincaid & E.W. Cheney
c     Brooks/Cole Publ., 2002
c     Copyright (c) 1996
c
c     Section 5.3
c
c     Example of modified Gram-Schmidt algorithm
c
c
c     file: qrshif.f
c
      parameter (n=4,m=10)
      dimension a(n,n),q(n,n),r(n,n)
      data (a(1,j),j=1,n)/1.,2.,3.,4./
      data (a(2,j),j=1,n)/5.,6.,7.,8./
      data (a(3,j),j=1,n)/0.,9.,10.,11./
      data (a(4,j),j=1,n)/0.,0.,12.,13./
c
      print *
      print *,' Modified Gram-Schimdt example'
      print *,' Section 5.3, Kincaid-Cheney'
      print *
c
      do 4 k=1,m
         print *,' Matrix A'
         call prtmtx(n,a)
         z = a(n,n)
         do 2 i=1,n
            a(i,i)= a(i,i) - z
 2       continue
         call mgs(a,q,r,n,n)
         print *,' Matrix Q'
         call prtmtx(n,q)
         print *,' Matrix R'
         call prtmtx(n,r)
         call mult(n,r,q,a)
         do 3 i=1,n
            a(i,i)=a(i,i) + z
 3       continue
 4    continue
c
      stop
      end 
c
      subroutine prtmtx(n,a)
c
c     print array a
c
      dimension a(n,n)
c
      do 2 i=1,n
         print 3,(a(i,j),j=1,n)
 2    continue
      print *
c
      return
 3    format(2x,4(e13.6,2x))
      end
c 
      subroutine mgs(a,q,t,m,n)
c
c     modified Gram-Schimdt
c
      dimension a(m,n),q(m,n),t(n,n)
c
      do 3 j=1,n
         do 2 i=1,n
            t(i,j)=0.
 2       continue
         do 3 i=1,m
            q(i,j)=a(i,j) 
 3    continue
      do 7 k=1,n
         z=0.
         do 4 i=1,m
            z=z+q(i,k)**2 
 4       continue
         t(k,k)=sqrt(z)
         do 5 i=1,m
            q(i,k)=q(i,k)/t(k,k)
 5       continue
         do 7 j=k+1,n
            z=0.
            do 6 i=1,m
               z=z+q(i,j)*q(i,k)
 6          continue
            t(k,j)=z
            do 7 i=1,m
               q(i,j)=q(i,j)-t(k,j)*q(i,k)
 7     continue
c
      return
      end 
c 
      subroutine mult (n,a,b,c)
c
c     compute matrix product c=ab
c
      dimension a(n,n),b(n,n),c(n,n)
c
      do 3 i=1,n
         do 3 j=1,n
            x=0.0
            do 2 k=1,n
               x=x+a(i,k)*b(k,j) 
 2          continue
 3    c(i,j)=x
c
      return
      end 
