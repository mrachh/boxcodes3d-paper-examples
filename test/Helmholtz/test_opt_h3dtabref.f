      implicit real *8 (a-h,o-z)
      complex *16 zk, im, zero, one
      complex *16, allocatable :: tab_ref(:,:),tab(:,:)
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

      tol = 1.0d-6


      allocate(tab_ref(ntarg0,npol3),tab(ntarg0,npol3))

      write(fname,'(a,i2.2,a,i1,a)') 'tabref_',n,'_izk_',izk,'.dat'

      open(unit=33,file=trim(fname))

c
c
c
c
      ifgen = 0
      if(ifgen.eq.1) then
         call h3dtabp_ref_brute(ndeg,zk,tab_ref,ntarg0)
         do i=1,npol3
           do j=1,ntarg0
              write(33,*) real(tab_ref(j,i)),imag(tab_ref(j,i))
           enddo
         enddo
         close(33)
         stop
      endif

      if(igen.ne.1) then
        do i=1,npol3
          do j=1,ntarg0
            read(33,*) tmp1,tmp2
            tab_ref(j,i) = tmp1+im*tmp2
          enddo
        enddo
        call h3dtabp_ref(ndeg,zk,tol,tab,ntarg0)
        
        errmax = 0
        do i=1,npol3
          do j=1,ntarg0
             erra = abs(tab(j,i)-tab_ref(j,i))/abs(tab_ref(j,i))
             if(erra.gt.errmax) errmax = erra
          enddo
        enddo
        call prin2('max error=*',erra,1)
        stop
      endif



      
      stop
      end
