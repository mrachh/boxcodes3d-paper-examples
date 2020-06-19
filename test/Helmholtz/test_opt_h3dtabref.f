      implicit real *8 (a-h,o-z)
      complex *16 zk, im, zero, one
      complex *16, allocatable :: tab_ref(:,:)
      character type

      character *100 fname
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /            
      
      call prini(6,13)
      

      izk = 1

      if(izk.eq.1) zk = 1.6d0
      if(izk.eq.2) zk = 0.16d0
      if(izk.eq.3) zk = 5.0d0
      if(izk.eq.4) zk = ima*1.6d0

      tol = 1.0d-11


      ndeg = 3
      n = ndeg + 1

      type = 't'
      call legetens_npol_3d(ndeg,type,npol3)

      ntarg0 = 10*n**3

      print *, ndeg,ntarg0,npol3


      allocate(tab_ref(ntarg0,npol3))

      write(fname,'(a,i2.2,a,i1,a)') 'tabref_',n,'_izk_',izk,'.dat'

      open(unit=33,file=trim(fname))

c
c
c
c
      ifgen = 1
      if(ifgen.eq.1) then
         print *, ndeg,n,ntarg0,npol3
         call get_tab_ref(zk,ndeg,n,ntarg0,npol3,tab_ref)
         do i=1,2
           do j=1,ntarg0
              write(33,*) real(tab_ref(j,i)),imag(tab_ref(j,i))
           enddo
         enddo
         close(33)
      endif

      if(igen.ne.1) then
        stop

      endif

      
      stop
      end
c
c      

      subroutine test1(zk,ndeg,itest,ix,iy,iz,tol)
      implicit real *8 (a-h,o-z)
      integer, allocatable :: ip2ind(:,:,:)
      real *8 xq(200), u, v, w, xyzc(3), xtest(3), tol
      real *8, allocatable :: xyz(:,:)
      complex *16, allocatable :: tab(:,:)
      complex *16 zk, im, zero, one, cv, ccv
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /            
      character type

      type = 't'
      call legetens_npol_3d(ndeg,type,npol3)

      call prin2('zk is *',zk,2)
      call prinf('ndeg is *',ndeg,1)      

      n = ndeg+1

      allocate(ip2ind(n,n,n))
      call legetens_pow2ind_3d(ndeg,type,ip2ind)
      
c
c     get reference points
c

      itype = 0
      call legeexps(itype,n,xq,u,v,w)

      do i = 1,n
         xq(i) = xq(i)+1.0d0
      enddo

      
      ntarg0 = 10*n**3
      allocate(xyz(3,ntarg0))

      xyzc(1) = -1.0d0
      xyzc(2) = -1.0d0
      xyzc(3) = -1.0d0
      bs = 2.0d0

      istart1 = 4*n**3+1
      istart2 = 7*n**3+1
      
      call tensrefpts3d(xq,n,bs,xyzc,xyz,xyz(1,istart1),
     1        xyz(1,istart2))

c
c     call table generator
c      

      allocate(tab(ntarg0,npol3))
      call cpu_time(t1)
c$    t1 = omp_get_wtime()      
      call h3dtabp_ref(ndeg,zk,tol,tab,ntarg0)
      call cpu_time(t2)
c$    t2 = omp_get_wtime()            

      call prin2('time to generate table *',t2-t1,1)
      
c
c     brute force one entry
c     

      xtest(1) = xyz(1,itest)
      xtest(2) = xyz(2,itest)
      xtest(3) = xyz(3,itest)

      ipol = ip2ind(ix,iy,iz)

      call mksurhelm3dp(xtest(1),xtest(2),xtest(3),ix,iy,iz,
     1     zk,n,cv,ifail)

      write(*,*) itest, ipol

      ccv = tab(itest,ipol)

c
c     compare
c

      err1 = abs(cv-ccv)

      call prinf('itest is *',itest,1)
      call prinf('ix is *',ix,1)
      call prinf('iy is *',iy,1)
      call prinf('iz is *',iz,1)
      call prinf('ifail is *',ifail,1)            
      call prin2_long('abs error is *',err1,1)
      call prin2_long('rel error is *',err1/abs(cv),1)
      print *, cv, ccv

      return
      end
c
c
c      
c
      subroutine get_tab_ref(zk,ndeg,n,ntarg0,npol3,tab)
      implicit real *8 (a-h,o-z)
      integer, allocatable :: iind2p(:,:)
      integer ifail
      real *8 xq(200), u, v, w, xyzc(3), xtest(3), tol
      real *8, allocatable :: xyz(:,:)
      complex *16 tab(ntarg0,npol3)
      complex *16 zk, im, zero, one, cv, ccv
      data im / (0.0d0,1.0d0) /
      data zero / (0.0d0,0.0d0) /
      data one / (1.0d0,0.0d0) /            
      character type


      call prin2('zk is *',zk,2)
      call prinf('ndeg is *',ndeg,1)      

      allocate(iind2p(3,npol3))
      type = 't'
      call legetens_ind2pow_3d(ndeg,type,iind2p)
      
c
c     get reference points
c

      itype = 0
      call legeexps(itype,n,xq,u,v,w)

      do i = 1,n
         xq(i) = xq(i)+1.0d0
      enddo

      allocate(xyz(3,ntarg0))

      xyzc(1) = -1.0d0
      xyzc(2) = -1.0d0
      xyzc(3) = -1.0d0
      bs = 2.0d0

      istart1 = 4*n**3+1
      istart2 = 7*n**3+1

      call prinf('n=*',n,1)
      call prin2('bs=*',bs,1)
      call prinf('ntarg0=*',ntarg0,1)
      
      call tensrefpts3d(xq,n,bs,xyzc,xyz,xyz(1,istart1),
     1        xyz(1,istart2))
      
c
c     brute force one entry
c     

      do ipol=1,2
        ix = iind2p(1,ipol)
        iy = iind2p(2,ipol)
        iz = iind2p(3,ipol)

        print *, ipol,ix,iy,iz

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ifail,xtest,cv) 
C$OMP$SCHEDULE(DYNAMIC)
        do itest=1,ntarg0
          xtest(1) = xyz(1,itest)
          xtest(2) = xyz(2,itest)
          xtest(3) = xyz(3,itest)
          ifail = 0
          cv = 0
          call mksurhelm3dp(xtest(1),xtest(2),xtest(3),ix,iy,iz,
     1         zk,n,cv,ifail)
          tab(itest,ipol) = cv
        enddo
C$OMP END PARALLEL DO        
      enddo

      return
      end
