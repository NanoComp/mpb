C***********************************************************************
c     Compute operator w = OP*v for elastic band structure.
c     Arguments:
c       p:  number of bands
c       ldp:  number of bands leading dimension (number of stored bands)
c       v(ldp,3,Nk,Nj,Ni): input = velocity field u * sqrt(rho)
c       w(ldp,3,Nk,Nj,Ni): output
c       edata: elastic_data structure (pointer), from elastic.c
c       sqrtrhoinv(Nk,Nj,Ni): 1 / sqrt(density = rho)
c       rhoct2(Nk,Nj,Ni): transverse velocity: rho*ct^2
c       rhocl2(Nk,Nj,Ni): long. velocity: rho*cl^2
c       vkx,vky,vkz: k vector, Cartesian coordinates (including 2*pi)
c       Nk,Nj,Ni: spatial grid
c       b1x,...,b3z: G vectors, Cartesian coordinates (including 2*pi)
c       Datx,Daty,Datz,Dat1..18: complex(Nk,Nj,Ni) work arrays

      Subroutine av (p,ldp, v, w, edata, 
     +     sqrtrhoinv,rhoct2,rhocl2,vkx,vky,vkz,
     +     Ni,Nj,Nk,
     +     b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z,
     +     Datx,Daty,Datz,Dat1,Dat2,Dat3,Dat4,Dat5,Dat6,Dat7,
     +     Dat8,Dat9,Dat10,Dat11,Dat12,Dat13,Dat14,Dat15,
     +     Dat16,Dat17,Dat18)
      implicit none
      Integer n, p, ldp, Ni,Nj,Nk, edata
      Complex*16 v(Ni*Nj*Nk*3), w(Ni*Nj*Nk*3)
      real*8 vkx,vky,vkz
      real*8 b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z
      Complex*16 tempx,tempy,tempz,Ureal,Uimag,hmx,hmy,hmz,h1,h2
      Complex*16 temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
      Complex*16 temp9,temp10,temp11,temp12,temp13,temp14,temp15
      Complex*16 temp16,temp17,temp18
      Complex*16 hm1,hm2,hm3,hm4,hm5,hm6,hm7,hm8,hm9,hm10
      Complex*16 hm11,hm12,hm13,hm14,hm15,hm16,hm17,hm18

      Real*8 Datx(Ni*Nj*Nk*2),Daty(Ni*Nj*Nk*2),Datz(Ni*Nj*Nk*2)
      Real*8 Dat1(Ni*Nj*Nk*2)
      Real*8 Dat2(Ni*Nj*Nk*2)
      Real*8 Dat3(Ni*Nj*Nk*2)
      Real*8 Dat4(Ni*Nj*Nk*2)
      Real*8 Dat5(Ni*Nj*Nk*2)
      Real*8 Dat6(Ni*Nj*Nk*2)
      Real*8 Dat7(Ni*Nj*Nk*2)
      Real*8 Dat8(Ni*Nj*Nk*2)
      Real*8 Dat9(Ni*Nj*Nk*2)
      Real*8 Dat10(Ni*Nj*Nk*2)
      Real*8 Dat11(Ni*Nj*Nk*2)
      Real*8 Dat12(Ni*Nj*Nk*2)
      Real*8 Dat13(Ni*Nj*Nk*2)
      Real*8 Dat14(Ni*Nj*Nk*2)
      Real*8 Dat15(Ni*Nj*Nk*2)
      Real*8 Dat16(Ni*Nj*Nk*2)
      Real*8 Dat17(Ni*Nj*Nk*2)
      Real*8 Dat18(Ni*Nj*Nk*2)
      Real*8 Dat19(Ni*Nj*Nk*2)

      Real*8 sqrtrhoinv(Ni*Nj*Nk),Rhoct2(Ni*Nj*Nk),Rhocl2(Ni*Nj*Nk)
      real*8 gx,gy,gz
      integer i,j,k,ip, ii,jj,kk
      Integer nn(3), ndim, isign
c
c     Computes w <--- OP*v
c
      n = Ni*Nj*Nk*3

      if (p.ne.1 .or. ldp.ne.1) stop 'arrgghh, p not 1 not supported'

      Ureal=(1.0,0.0)
      Uimag=(0.0,1.0)

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     Reciprocal to Real - /sqrt(rho) - Back to Reciprocal
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      nn(1)=Ni
      nn(2)=Nj
      nn(3)=Nk
      ndim=3

      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            hmx=v(3*ip-2)
            hmy=v(3*ip-1)
            hmz=v(3*ip  )
            Datx (2*ip-1) = Dble( hmx )
            Datx (2*ip  ) = Imag( hmx )
            Daty (2*ip-1) = Dble( hmy )
            Daty (2*ip  ) = Imag( hmy )
            Datz (2*ip-1) = Dble( hmz )
            Datz (2*ip  ) = Imag( hmz )
          Enddo
        Enddo
      Enddo

      isign=1
      Call elasticFFT(isign, edata, datx,1,1,1)
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Datx (2*ip-1)=Datx(2*ip-1)*sqrtrhoinv(ip)           ! x
            Datx (2*ip  )=Datx(2*ip  )*sqrtrhoinv(ip)           ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, datx,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Datx(i)=Datx(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, daty,1,1,1)
      Do i=0,Ni-1                                               ! y
        Do j=0,Nj-1                                             ! y
          Do k=0,Nk-1                                           ! y
            ip= i*Nj*Nk + j*Nk + k+1                            ! y
            Daty (2*ip-1)=Daty(2*ip-1)*sqrtrhoinv(ip)           ! y
            Daty (2*ip  )=Daty(2*ip  )*sqrtrhoinv(ip)           ! y
          Enddo                                                 ! y
        Enddo                                                   ! y
      Enddo                                                     ! y
      isign=-1                              ! -1                ! y
      Call elasticFFT(isign, edata, daty,1,1,1)                 ! y
      Do i=1,2*Ni*Nj*Nk                                         ! y
        Daty(i)=Daty(i)/(Ni*Nj*Nk)                              ! y
      Enddo                                                     ! y

      isign=1                               ! +1
      Call elasticFFT(isign, edata, datz,1,1,1)
      Do i=0,Ni-1                                               ! z
        Do j=0,Nj-1                                             ! z
          Do k=0,Nk-1                                           ! z
            ip= i*Nj*Nk + j*Nk + k+1                            ! z
            Datz (2*ip-1)=Datz(2*ip-1)*sqrtrhoinv(ip)           ! z
            Datz (2*ip  )=Datz(2*ip  )*sqrtrhoinv(ip)           ! z
          Enddo                                                 ! z
        Enddo                                                   ! z
      Enddo                                                     ! z
      isign=-1                              ! -1                ! z
      Call elasticFFT(isign, edata, datz,1,1,1)                 ! z
      Do i=1,2*Ni*Nj*Nk                                         ! z
        Datz(i)=Datz(i)/(Ni*Nj*Nk)                              ! z
      Enddo                                                     ! z

      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            hmx= Ureal * Datx(2*ip-1) + Uimag * Datx(2*ip)
            hmy= Ureal * Daty(2*ip-1) + Uimag * Daty(2*ip)
            hmz= Ureal * Datz(2*ip-1) + Uimag * Datz(2*ip)
            w(3*ip-2)=hmx
            w(3*ip-1)=hmy
            w(3*ip  )=hmz
          Enddo
        Enddo
      Enddo

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     1st Reciprocal Calculation - To Real- *rhoct2 - Back to Reciprocal 
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            ii=i
            jj=j
            kk=k
            If(ii.gt.Ni/2) ii=ii-Ni
            If(jj.gt.Nj/2) jj=jj-Nj
            If(kk.gt.Nk/2) kk=kk-Nk
            Gx = ii*b1x + jj*b2x + kk*b3x     ! Reciprocal Vectors
            Gy = ii*b1y + jj*b2y + kk*b3y     ! Reciprocal Vectors
            Gz = ii*b1z + jj*b2z + kk*b3z     ! Reciprocal Vectors
            Gx = vkx + Gx                                    ! k+G
            Gy = vky + Gy                                    ! k+G
            Gz = vkz + Gz                                    ! k+G
            hmx=w(3*ip-2)
            hmy=w(3*ip-1)
            hmz=w(3*ip  )

            Dat1 (2*ip-1) = Dble( Gx*hmx )
            Dat1 (2*ip  ) = Imag( Gx*hmx )
            Dat2 (2*ip-1) = Dble( Gy*hmx )
            Dat2 (2*ip  ) = Imag( Gy*hmx )
            Dat3 (2*ip-1) = Dble( Gz*hmx )
            Dat3 (2*ip  ) = Imag( Gz*hmx )
            Dat4 (2*ip-1) = Dble( Gx*hmy )
            Dat4 (2*ip  ) = Imag( Gx*hmy )
            Dat5 (2*ip-1) = Dble( Gy*hmy )
            Dat5 (2*ip  ) = Imag( Gy*hmy )
            Dat6 (2*ip-1) = Dble( Gz*hmy )
            Dat6 (2*ip  ) = Imag( Gz*hmy )
            Dat7 (2*ip-1) = Dble( Gx*hmz )
            Dat7 (2*ip  ) = Imag( Gx*hmz )
            Dat8 (2*ip-1) = Dble( Gy*hmz )
            Dat8 (2*ip  ) = Imag( Gy*hmz )
            Dat9 (2*ip-1) = Dble( Gz*hmz )
            Dat9 (2*ip  ) = Imag( Gz*hmz )

            Dat10 (2*ip-1) = Dble( Gx*hmx )
            Dat10 (2*ip  ) = Imag( Gx*hmx )
c            Dat11 (2*ip-1) = Dble( Gy*hmx )
c            Dat11 (2*ip  ) = Imag( Gy*hmx )
c            Dat12 (2*ip-1) = Dble( Gz*hmx )
c            Dat12 (2*ip  ) = Imag( Gz*hmx )
c            Dat13 (2*ip-1) = Dble( Gx*hmy )
c            Dat13 (2*ip  ) = Imag( Gx*hmy )
            Dat14 (2*ip-1) = Dble( Gy*hmy )
            Dat14 (2*ip  ) = Imag( Gy*hmy )
c            Dat15 (2*ip-1) = Dble( Gz*hmy )
c            Dat15 (2*ip  ) = Imag( Gz*hmy )
c            Dat16 (2*ip-1) = Dble( Gx*hmz )
c            Dat16 (2*ip  ) = Imag( Gx*hmz )
c            Dat17 (2*ip-1) = Dble( Gy*hmz )
c            Dat17 (2*ip  ) = Imag( Gy*hmz )
            Dat18 (2*ip-1) = Dble( Gz*hmz )
            Dat18 (2*ip  ) = Imag( Gz*hmz )

          Enddo
        Enddo
      Enddo

      nn(1)=Ni
      nn(2)=Nj
      nn(3)=Nk
      ndim=3

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat1,1,1,1)                 ! 1
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat1 (2*ip-1) = Dat1 (2*ip-1) * rhoct2(ip)          ! x
            Dat1 (2*ip  ) = Dat1 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat1,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat1(i)=Dat1(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat2,1,1,1)                 ! 2
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat2 (2*ip-1) = Dat2 (2*ip-1) * rhoct2(ip)          ! x
            Dat2 (2*ip  ) = Dat2 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat2,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat2(i)=Dat2(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat3,1,1,1)                 ! 3
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat3 (2*ip-1) = Dat3 (2*ip-1) * rhoct2(ip)          ! x
            Dat3 (2*ip  ) = Dat3 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat3,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat3(i)=Dat3(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat4,1,1,1)                 ! 4
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat4 (2*ip-1) = Dat4 (2*ip-1) * rhoct2(ip)          ! x
            Dat4 (2*ip  ) = Dat4 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat4,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat4(i)=Dat4(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat5,1,1,1)                 ! 5
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat5 (2*ip-1) = Dat5 (2*ip-1) * rhoct2(ip)          ! x
            Dat5 (2*ip  ) = Dat5 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat5,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat5(i)=Dat5(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat6,1,1,1)                 ! 6
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat6 (2*ip-1) = Dat6 (2*ip-1) * rhoct2(ip)          ! x
            Dat6 (2*ip  ) = Dat6 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat6,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat6(i)=Dat6(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat7,1,1,1)                 ! 7
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat7 (2*ip-1) = Dat7 (2*ip-1) * rhoct2(ip)          ! x
            Dat7 (2*ip  ) = Dat7 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat7,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat7(i)=Dat7(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat8,1,1,1)                 ! 8
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat8 (2*ip-1) = Dat8 (2*ip-1) * rhoct2(ip)          ! x
            Dat8 (2*ip  ) = Dat8 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat8,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat8(i)=Dat8(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat9,1,1,1)                 ! 9
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat9 (2*ip-1) = Dat9 (2*ip-1) * rhoct2(ip)          ! x
            Dat9 (2*ip  ) = Dat9 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat9,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat9(i)=Dat9(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat10,1,1,1)                !10
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat10(2*ip-1)=Dat10(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat10(2*ip  )=Dat10(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat10,1,1,1)                ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat10(i)=Dat10(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

c      isign=1                               ! +1
c      Call elasticFFT(isign, edata, dat11,1,1,1)                !11
c      Do i=0,Ni-1                                               ! x
c        Do j=0,Nj-1                                             ! x
c          Do k=0,Nk-1                                           ! x
c            ip= i*Nj*Nk + j*Nk + k+1                            ! x
c           Dat11(2*ip-1)=Dat11(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
c           Dat11(2*ip  )=Dat11(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
c          Enddo                                                 ! x
c        Enddo                                                   ! x
c      Enddo                                                     ! x
c      isign=-1                              ! -1                ! x
c      Call elasticFFT(isign, edata, dat11,1,1,1)                ! x        
c      Do i=1,2*Ni*Nj*Nk                                         ! x
c        Dat11(i)=Dat11(i)/(Ni*Nj*Nk)                            ! x
c      Enddo                                                     ! x

c      isign=1                               ! +1
c      Call elasticFFT(isign, edata, dat12,1,1,1)                !12
c      Do i=0,Ni-1                                               ! x
c        Do j=0,Nj-1                                             ! x 
c          Do k=0,Nk-1                                           ! x
c            ip= i*Nj*Nk + j*Nk + k+1                            ! x
c           Dat12(2*ip-1)=Dat12(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
c           Dat12(2*ip  )=Dat12(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
c          Enddo                                                 ! x
c        Enddo                                                   ! x
c      Enddo                                                     ! x
c      isign=-1                              ! -1                ! x
c      Call elasticFFT(isign, edata, dat12,1,1,1)                ! x       
c      Do i=1,2*Ni*Nj*Nk                                         ! x
c        Dat12(i)=Dat12(i)/(Ni*Nj*Nk)                            ! x
c      Enddo                                                     ! x

c      isign=1                               ! +1
c      Call elasticFFT(isign, edata, dat13,1,1,1)                !13
c      Do i=0,Ni-1                                               ! x
c        Do j=0,Nj-1                                             ! x
c          Do k=0,Nk-1                                           ! x
c            ip= i*Nj*Nk + j*Nk + k+1                            ! x
c           Dat13(2*ip-1)=Dat13(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
c           Dat13(2*ip  )=Dat13(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
c          Enddo                                                 ! x
c        Enddo                                                   ! x
c      Enddo                                                     ! x
c      isign=-1                              ! -1                ! x
c      Call elasticFFT(isign, edata, dat13,1,1,1)                ! x     
c      Do i=1,2*Ni*Nj*Nk                                         ! x
c        Dat13(i)=Dat13(i)/(Ni*Nj*Nk)                            ! x
c      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat14,1,1,1)                !14
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat14(2*ip-1)=Dat14(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat14(2*ip  )=Dat14(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat14,1,1,1)                ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat14(i)=Dat14(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

c      isign=1                               ! +1
c      Call elasticFFT(isign, edata, dat15,1,1,1)                !15 
c      Do i=0,Ni-1                                               ! x
c        Do j=0,Nj-1                                             ! x
c          Do k=0,Nk-1                                           ! x
c            ip= i*Nj*Nk + j*Nk + k+1                            ! x
c           Dat15(2*ip-1)=Dat15(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
c           Dat15(2*ip  )=Dat15(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
c          Enddo                                                 ! x
c        Enddo                                                   ! x
c      Enddo                                                     ! x
c      isign=-1                              ! -1                ! x
c      Call elasticFFT(isign, edata, dat15,1,1,1)                ! x
c      Do i=1,2*Ni*Nj*Nk                                         ! x
c        Dat15(i)=Dat15(i)/(Ni*Nj*Nk)                            ! x
c      Enddo                                                     ! x

c      isign=1                               ! +1
c      Call elasticFFT(isign, edata, dat16,1,1,1)                !16
c      Do i=0,Ni-1                                               ! x
c        Do j=0,Nj-1                                             ! x
c          Do k=0,Nk-1                                           ! x
c            ip= i*Nj*Nk + j*Nk + k+1                            ! x
c           Dat16(2*ip-1)=Dat16(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
c           Dat16(2*ip  )=Dat16(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
c          Enddo                                                 ! x
c        Enddo                                                   ! x
c      Enddo                                                     ! x
c      isign=-1                              ! -1                ! x
c      Call elasticFFT(isign, edata, dat16,1,1,1)                ! x
c      Do i=1,2*Ni*Nj*Nk                                         ! x
c        Dat16(i)=Dat16(i)/(Ni*Nj*Nk)                            ! x
c      Enddo                                                     ! x

c      isign=1                               ! +1
c      Call elasticFFT(isign, edata, dat17,1,1,1)                !17
c      Do i=0,Ni-1                                               ! x
c        Do j=0,Nj-1                                             ! x
c          Do k=0,Nk-1                                           ! x
c            ip= i*Nj*Nk + j*Nk + k+1                            ! x
c           Dat17(2*ip-1)=Dat17(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
c           Dat17(2*ip  )=Dat17(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
c          Enddo                                                 ! x
c        Enddo                                                   ! x
c      Enddo                                                     ! x
c      isign=-1                              ! -1                ! x
c      Call elasticFFT(isign, edata, dat17,1,1,1)                ! x       
c      Do i=1,2*Ni*Nj*Nk                                         ! x
c        Dat17(i)=Dat17(i)/(Ni*Nj*Nk)                            ! x
c      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat18,1,1,1)                !18
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat18(2*ip-1)=Dat18(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat18(2*ip  )=Dat18(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat18,1,1,1)                ! x     
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat18(i)=Dat18(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     2nd Calculation in Reciprocal Space
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     
      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            ii=i
            jj=j
            kk=k
            If(ii.gt.Ni/2) ii=ii-Ni
            If(jj.gt.Nj/2) jj=jj-Nj
            If(kk.gt.Nk/2) kk=kk-Nk
            Gx = ii*b1x + jj*b2x + kk*b3x     ! Reciprocal Vectors
            Gy = ii*b1y + jj*b2y + kk*b3y     ! Reciprocal Vectors
            Gz = ii*b1z + jj*b2z + kk*b3z     ! Reciprocal Vectors
            Gx = vkx + Gx                                    ! k+G
            Gy = vky + Gy                                    ! k+G
            Gz = vkz + Gz                                    ! k+G
            hm1= Ureal * Dat1(2*ip-1) + Uimag * Dat1(2*ip)
            hm2= Ureal * Dat2(2*ip-1) + Uimag * Dat2(2*ip)
            hm3= Ureal * Dat3(2*ip-1) + Uimag * Dat3(2*ip)
            hm4= Ureal * Dat4(2*ip-1) + Uimag * Dat4(2*ip)
            hm5= Ureal * Dat5(2*ip-1) + Uimag * Dat5(2*ip)
            hm6= Ureal * Dat6(2*ip-1) + Uimag * Dat6(2*ip)
            hm7= Ureal * Dat7(2*ip-1) + Uimag * Dat7(2*ip)
            hm8= Ureal * Dat8(2*ip-1) + Uimag * Dat8(2*ip)
            hm9= Ureal * Dat9(2*ip-1) + Uimag * Dat9(2*ip)

  	    hm10= Ureal * Dat10(2*ip-1) + Uimag * Dat10(2*ip)
c            hm11= Ureal * Dat11(2*ip-1) + Uimag * Dat11(2*ip)
c            hm12= Ureal * Dat12(2*ip-1) + Uimag * Dat12(2*ip)
c            hm13= Ureal * Dat13(2*ip-1) + Uimag * Dat13(2*ip)
            hm14= Ureal * Dat14(2*ip-1) + Uimag * Dat14(2*ip)
c            hm15= Ureal * Dat15(2*ip-1) + Uimag * Dat15(2*ip)
c            hm16= Ureal * Dat16(2*ip-1) + Uimag * Dat16(2*ip)
c            hm17= Ureal * Dat17(2*ip-1) + Uimag * Dat17(2*ip)
            hm18= Ureal * Dat18(2*ip-1) + Uimag * Dat18(2*ip)

            tempx=(Gx*hm1+Gy*hm2+Gz*hm3)+
     +		  (Gx*hm1+Gy*hm4+Gz*hm7)+
     +	          (Gx*hm10+Gx*hm14+Gx*hm18)

  	    tempy=(Gx*hm4+Gy*hm5+Gz*hm6)+
     +            (Gx*hm2+Gy*hm5+Gz*hm8)+
     +	          (Gy*hm10+Gy*hm14+Gy*hm18)

	    tempz=(Gx*hm7+Gy*hm8+Gz*hm9)+
     +            (Gx*hm3+Gy*hm6+Gz*hm9)+
     +	          (Gz*hm10+Gz*hm14+Gz*hm18)

             w(3*ip-2)=tempx
             w(3*ip-1)=tempy
             w(3*ip  )=tempz

          Enddo
        Enddo
      Enddo

C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
C     To Real - /sqrt(rho) - Back to Reciprocal  
C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      nn(1)=Ni
      nn(2)=Nj
      nn(3)=Nk
      ndim=3

      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            hmx=w(3*ip-2)
            hmy=w(3*ip-1)
            hmz=w(3*ip  )
            Datx (2*ip-1) = Dble( hmx )
            Datx (2*ip  ) = Imag( hmx )
            Daty (2*ip-1) = Dble( hmy )
            Daty (2*ip  ) = Imag( hmy )
            Datz (2*ip-1) = Dble( hmz )
            Datz (2*ip  ) = Imag( hmz )
          Enddo
        Enddo
      Enddo

	isign=1                               ! +1
      Call elasticFFT(isign, edata, datx,1,1,1)
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Datx (2*ip-1)=Datx(2*ip-1)*sqrtrhoinv(ip)           ! x
            Datx (2*ip  )=Datx(2*ip  )*sqrtrhoinv(ip)           ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, datx,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Datx(i)=Datx(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, daty,1,1,1)
      Do i=0,Ni-1                                               ! y
        Do j=0,Nj-1                                             ! y
          Do k=0,Nk-1                                           ! y
            ip= i*Nj*Nk + j*Nk + k+1                            ! y
            Daty (2*ip-1)=Daty(2*ip-1)*sqrtrhoinv(ip)           ! y
            Daty (2*ip  )=Daty(2*ip  )*sqrtrhoinv(ip)           ! y
          Enddo                                                 ! y
        Enddo                                                   ! y
      Enddo                                                     ! y
      isign=-1                              ! -1                ! y
      Call elasticFFT(isign, edata, daty,1,1,1)                 ! y
      Do i=1,2*Ni*Nj*Nk                                         ! y
        Daty(i)=Daty(i)/(Ni*Nj*Nk)                              ! y
      Enddo                                                     ! y

      isign=1                               ! +1
      Call elasticFFT(isign, edata, datz,1,1,1)
      Do i=0,Ni-1                                               ! z
        Do j=0,Nj-1                                             ! z
          Do k=0,Nk-1                                           ! z
            ip= i*Nj*Nk + j*Nk + k+1                            ! z
            Datz (2*ip-1)=Datz(2*ip-1)*sqrtrhoinv(ip)           ! z
            Datz (2*ip  )=Datz(2*ip  )*sqrtrhoinv(ip)           ! z
          Enddo                                                 ! z
        Enddo                                                   ! z
      Enddo                                                     ! z
      isign=-1                              ! -1                ! z
      Call elasticFFT(isign, edata, datz,1,1,1)                 ! z
      Do i=1,2*Ni*Nj*Nk                                         ! z
        Datz(i)=Datz(i)/(Ni*Nj*Nk)                              ! z
      Enddo                                                     ! z

      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            hmx= Ureal * Datx(2*ip-1) + Uimag * Datx(2*ip)
            hmy= Ureal * Daty(2*ip-1) + Uimag * Daty(2*ip)
            hmz= Ureal * Datz(2*ip-1) + Uimag * Datz(2*ip)
            w(3*ip-2)=hmx
            w(3*ip-1)=hmy
            w(3*ip  )=hmz
          Enddo
        Enddo
      Enddo

      Return
      End

***********************************************************************
c     Compute preconditioner operator w = OP*v for elastic band structure, where
c     OP is an approximate inverse to the eigen-operator.
c     Arguments:
c       p:  number of bands
c       ldp:  number of bands leading dimension (number of stored bands)
c       v(ldp,3,Nk,Nj,Ni): input = velocity field u * sqrt(rho)
c       w(ldp,3,Nk,Nj,Ni): output
c       edata: elastic_data structure (pointer), from elastic.c
c       sqrtrhoinv(Nk,Nj,Ni): 1 / sqrt(density = rho)
c       rhoct2(Nk,Nj,Ni): transverse velocity: rho*ct^2
c       rhocl2(Nk,Nj,Ni): long. velocity: rho*cl^2
c       vkx,vky,vkz: k vector, Cartesian coordinates (including 2*pi)
c       Nk,Nj,Ni: spatial grid
c       b1x,...,b3z: G vectors, Cartesian coordinates (including 2*pi)
c       Datx,Daty,Datz,Dat1..18: complex(Nk,Nj,Ni) work arrays

      Subroutine precond(p,ldp, v, w, edata, 
     +     sqrtrhoinv,rhoct2,rhocl2,vkx,vky,vkz,
     +     Ni,Nj,Nk,
     +     b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z,
     +     Datx,Daty,Datz,Dat1,Dat2,Dat3,Dat4,Dat5,Dat6,Dat7,
     +     Dat8,Dat9,Dat10,Dat11,Dat12,Dat13,Dat14,Dat15,
     +     Dat16,Dat17,Dat18)
      implicit none
      Integer n, p, ldp, Ni,Nj,Nk, edata
ccc      Complex*16 v(ldp,Ni*Nj*Nk*3), w(ldp,Ni*Nj*Nk*3)
      Complex*16 v(Ni*Nj*Nk*3),w(Ni*Nj*Nk*3)
      real*8 vkx,vky,vkz
      real*8 b1x,b1y,b1z,b2x,b2y,b2z,b3x,b3y,b3z
      real*8 denom
      Complex*16 tempx,tempy,tempz,Ureal,Uimag,hmx,hmy,hmz,h1,h2
      Complex*16 temp1,temp2,temp3,temp4,temp5,temp6,temp7,temp8
      Complex*16 temp9,temp10,temp11,temp12,temp13,temp14,temp15
      Complex*16 temp16,temp17,temp18
      Complex*16 hm1,hm2,hm3,hm4,hm5,hm6,hm7,hm8,hm9,hm10
      Complex*16 hm11,hm12,hm13,hm14,hm15,hm16,hm17,hm18
      real*8 X1,X2,X3,X4,X5,X6,X7,X8,X9

      Real*8 Datx(Ni*Nj*Nk*2),Daty(Ni*Nj*Nk*2),Datz(Ni*Nj*Nk*2)
      Real*8 Dat1(Ni*Nj*Nk*2)
      Real*8 Dat2(Ni*Nj*Nk*2)
      Real*8 Dat3(Ni*Nj*Nk*2)
      Real*8 Dat4(Ni*Nj*Nk*2)
      Real*8 Dat5(Ni*Nj*Nk*2)
      Real*8 Dat6(Ni*Nj*Nk*2)
      Real*8 Dat7(Ni*Nj*Nk*2)
      Real*8 Dat8(Ni*Nj*Nk*2)
      Real*8 Dat9(Ni*Nj*Nk*2)
      Real*8 Dat10(Ni*Nj*Nk*2)
      Real*8 Dat11(Ni*Nj*Nk*2)
      Real*8 Dat12(Ni*Nj*Nk*2)
      Real*8 Dat13(Ni*Nj*Nk*2)
      Real*8 Dat14(Ni*Nj*Nk*2)
      Real*8 Dat15(Ni*Nj*Nk*2)
      Real*8 Dat16(Ni*Nj*Nk*2)
      Real*8 Dat17(Ni*Nj*Nk*2)
      Real*8 Dat18(Ni*Nj*Nk*2)
      Real*8 Dat19(Ni*Nj*Nk*2)

      Real*8 sqrtrhoinv(Ni*Nj*Nk),Rhoct2(Ni*Nj*Nk),Rhocl2(Ni*Nj*Nk)
      real*8 gx,gy,gz
      integer i,j,k,ip, ii,jj,kk
      Integer nn(3), ndim, isign
c
c     Computes w <--- OP*v
c
      n = Ni*Nj*Nk*3

      if (p.ne.1 .or. ldp.ne.1) stop 'arrgghh, p not 1 not supported'

c     stupid identity preconditioner
      Do i=0,Ni-1
         Do j=0,Nj-1
            Do k=0,Nk-1
               ip= i*Nj*Nk + j*Nk + k+1
cc               do ii=1,p
cc                  w(ii, 3*ip - 2) = v(ii, 3*ip - 2)
cc                  w(ii, 3*ip - 1) = v(ii, 3*ip - 1)
cc                  w(ii, 3*ip - 0) = v(ii, 3*ip - 0)
cc               enddo
            w(3*ip-2)=v(3*ip-2)
            w(3*ip-1)=v(3*ip-1)
            w(3*ip  )=v(3*ip  )
          Enddo
        Enddo
      Enddo

      Ureal=(1.0,0.0)
      Uimag=(0.0,1.0)

C-----------------------------------------------------------------
C     From Reciprocal To Real - *sqrt(rho) -Back to Reciprocal
C-----------------------------------------------------------------
 
      nn(1)=Ni
      nn(2)=Nj
      nn(3)=Nk
      ndim=3

      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            hmx=v(3*ip-2)
            hmy=v(3*ip-1)
            hmz=v(3*ip  )
            Datx (2*ip-1) = Dble( hmx )
            Datx (2*ip  ) = Imag( hmx )
            Daty (2*ip-1) = Dble( hmy )
            Daty (2*ip  ) = Imag( hmy )
            Datz (2*ip-1) = Dble( hmz )
            Datz (2*ip  ) = Imag( hmz )
          Enddo
        Enddo
      Enddo

      isign=1                               ! +1
      Call elasticFFT(isign, edata, datx,1,1,1)
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Datx (2*ip-1)=Datx(2*ip-1)/sqrtrhoinv(ip)           ! x ###
            Datx (2*ip  )=Datx(2*ip  )/sqrtrhoinv(ip)           ! x ###
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, datx,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Datx(i)=Datx(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, daty,1,1,1)
      Do i=0,Ni-1                                               ! y
        Do j=0,Nj-1                                             ! y
          Do k=0,Nk-1                                           ! y
            ip= i*Nj*Nk + j*Nk + k+1                            ! y
            Daty (2*ip-1)=Daty(2*ip-1)/sqrtrhoinv(ip)           ! y ###
            Daty (2*ip  )=Daty(2*ip  )/sqrtrhoinv(ip)           ! y ###
          Enddo                                                 ! y 
        Enddo                                                   ! y
      Enddo                                                     ! y
      isign=-1                              ! -1                ! y
      Call elasticFFT(isign, edata, daty,1,1,1)                 ! y
      Do i=1,2*Ni*Nj*Nk                                         ! y
        Daty(i)=Daty(i)/(Ni*Nj*Nk)                              ! y
      Enddo                                                     ! y

      isign=1                               ! +1
      Call elasticFFT(isign, edata, datz,1,1,1)
      Do i=0,Ni-1                                               ! z
        Do j=0,Nj-1                                             ! z
          Do k=0,Nk-1                                           ! z
            ip= i*Nj*Nk + j*Nk + k+1                            ! z
            Datz (2*ip-1)=Datz(2*ip-1)/sqrtrhoinv(ip)           ! z ###
            Datz (2*ip  )=Datz(2*ip  )/sqrtrhoinv(ip)           ! z ###
          Enddo                                                 ! z
        Enddo                                                   ! z
      Enddo                                                     ! z
      isign=-1                              ! -1                ! z
      Call elasticFFT(isign, edata, datz,1,1,1)                 ! z
      Do i=1,2*Ni*Nj*Nk                                         ! z
        Datz(i)=Datz(i)/(Ni*Nj*Nk)                              ! z
      Enddo                                                     ! z


      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            hmx= Ureal * Datx(2*ip-1) + Uimag * Datx(2*ip)
            hmy= Ureal * Daty(2*ip-1) + Uimag * Daty(2*ip)
            hmz= Ureal * Datz(2*ip-1) + Uimag * Datz(2*ip)
            w(3*ip-2)=hmx
            w(3*ip-1)=hmy
            w(3*ip  )=hmz
          Enddo
        Enddo
      Enddo

C-------------------------------------------------------------
C     All Reciprocal
C-------------------------------------------------------------

      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            ii=i
            jj=j
            kk=k
            If(ii.gt.Ni/2) ii=ii-Ni
            If(jj.gt.Nj/2) jj=jj-Nj
            If(kk.gt.Nk/2) kk=kk-Nk
            Gx = ii*b1x + jj*b2x + kk*b3x     ! Reciprocal Vectors
            Gy = ii*b1y + jj*b2y + kk*b3y     ! Reciprocal Vectors
            Gz = ii*b1z + jj*b2z + kk*b3z     ! Reciprocal Vectors
            Gx = vkx + Gx                                    ! k+G
            Gy = vky + Gy                                    ! k+G
            Gz = vkz + Gz                                    ! k+G
            hmx=w(3*ip-2)
            hmy=w(3*ip-1)
            hmz=w(3*ip  )

           
            X1=3.*Gx*Gx+Gy*Gy+Gz*Gz
            X2=Gy*Gx+Gx*Gy
            X3=Gz*Gx+Gx*Gz
            X4=Gx*Gy+Gy*Gx
            X5=Gx*Gx+3.*Gy*Gy+Gz*Gz
            X6=Gz*Gy+Gy*Gz
            X7=Gx*Gz+Gz*Gx
            X8=Gy*Gz+Gz*Gy
            X9=Gx*Gx+Gy*Gy+3.*Gz*Gz

            denom=-x3*x5*x7+x2*x6*x7+x3*x4*x8-x1*x6*x8-x2*x4*x9+x1*x5*x9 

c           add small quantity (5e-3) to denom so that we don't divide by zero
c           for k+G = 0 (only happens for Gamma point, i.e. k=0):
            denom = 1.0 / (denom + sign(5D-3, denom))


          tempx=(-X6*X8+X5*X9)*hmx+( X3*X8-X2*X9)*hmy+(-X3*X5+X2*X6)*hmz
          tempy=( X6*X7-X4*X9)*hmx+(-X3*X7+X1*X9)*hmy+( X3*X4-X1*X6)*hmz
          tempz=(-X5*X7+X4*X8)*hmx+( X2*X7-X1*X8)*hmy+(-X2*X4+X1*X5)*hmz 
 
            
            

            w(3*ip-2)=tempx*denom
            w(3*ip-1)=tempy*denom
            w(3*ip  )=tempz*denom


            Dat1 (2*ip-1) = Dble( Gx*hmx )
            Dat1 (2*ip  ) = Imag( Gx*hmx )
            Dat2 (2*ip-1) = Dble( Gy*hmx )
            Dat2 (2*ip  ) = Imag( Gy*hmx )
            Dat3 (2*ip-1) = Dble( Gz*hmx )
            Dat3 (2*ip  ) = Imag( Gz*hmx )
            Dat4 (2*ip-1) = Dble( Gx*hmy )
            Dat4 (2*ip  ) = Imag( Gx*hmy )
            Dat5 (2*ip-1) = Dble( Gy*hmy )
            Dat5 (2*ip  ) = Imag( Gy*hmy )
            Dat6 (2*ip-1) = Dble( Gz*hmy )
            Dat6 (2*ip  ) = Imag( Gz*hmy )
            Dat7 (2*ip-1) = Dble( Gx*hmz )
            Dat7 (2*ip  ) = Imag( Gx*hmz )
            Dat8 (2*ip-1) = Dble( Gy*hmz )
            Dat8 (2*ip  ) = Imag( Gy*hmz )
            Dat9 (2*ip-1) = Dble( Gz*hmz )
            Dat9 (2*ip  ) = Imag( Gz*hmz )

		
            Dat10 (2*ip-1) = Dble( Gx*hmx )
            Dat10 (2*ip  ) = Imag( Gx*hmx )
c            Dat11 (2*ip-1) = Dble( Gy*hmx )
c            Dat11 (2*ip  ) = Imag( Gy*hmx )
c            Dat12 (2*ip-1) = Dble( Gz*hmx )
c            Dat12 (2*ip  ) = Imag( Gz*hmx )
c            Dat13 (2*ip-1) = Dble( Gx*hmy )
c            Dat13 (2*ip  ) = Imag( Gx*hmy )
            Dat14 (2*ip-1) = Dble( Gy*hmy )
            Dat14 (2*ip  ) = Imag( Gy*hmy )
c            Dat15 (2*ip-1) = Dble( Gz*hmy )
c            Dat15 (2*ip  ) = Imag( Gz*hmy )
c            Dat16 (2*ip-1) = Dble( Gx*hmz )
c            Dat16 (2*ip  ) = Imag( Gx*hmz )
c            Dat17 (2*ip-1) = Dble( Gy*hmz )
c            Dat17 (2*ip  ) = Imag( Gy*hmz )
            Dat18 (2*ip-1) = Dble( Gz*hmz )
            Dat18 (2*ip  ) = Imag( Gz*hmz )

          Enddo
        Enddo
      Enddo

C/////SKIP////////////////////////////////////////////////////////

      Goto 111
      
      nn(1)=Ni
      nn(2)=Nj
      nn(3)=Nk
      ndim=3

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat1,1,1,1)                 ! 1
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat1 (2*ip-1) = Dat1 (2*ip-1) * rhoct2(ip)          ! x
            Dat1 (2*ip  ) = Dat1 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat1,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat1(i)=Dat1(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat2,1,1,1)                 ! 2
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat2 (2*ip-1) = Dat2 (2*ip-1) * rhoct2(ip)          ! x
            Dat2 (2*ip  ) = Dat2 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat2,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat2(i)=Dat2(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat3,1,1,1)                 ! 3
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat3 (2*ip-1) = Dat3 (2*ip-1) * rhoct2(ip)          ! x
            Dat3 (2*ip  ) = Dat3 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat3,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat3(i)=Dat3(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat4,1,1,1)                 ! 4
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat4 (2*ip-1) = Dat4 (2*ip-1) * rhoct2(ip)          ! x
            Dat4 (2*ip  ) = Dat4 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat4,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat4(i)=Dat4(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat5,1,1,1)                 ! 5
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat5 (2*ip-1) = Dat5 (2*ip-1) * rhoct2(ip)          ! x
            Dat5 (2*ip  ) = Dat5 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat5,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat5(i)=Dat5(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat6,1,1,1)                 ! 6
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat6 (2*ip-1) = Dat6 (2*ip-1) * rhoct2(ip)          ! x
            Dat6 (2*ip  ) = Dat6 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat6,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat6(i)=Dat6(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat7,1,1,1)                 ! 7
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat7 (2*ip-1) = Dat7 (2*ip-1) * rhoct2(ip)          ! x
            Dat7 (2*ip  ) = Dat7 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat7,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat7(i)=Dat7(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat8,1,1,1)                 ! 8 
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat8 (2*ip-1) = Dat8 (2*ip-1) * rhoct2(ip)          ! x
            Dat8 (2*ip  ) = Dat8 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat8,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat8(i)=Dat8(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat9,1,1,1)                 ! 9
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Dat9 (2*ip-1) = Dat9 (2*ip-1) * rhoct2(ip)          ! x
            Dat9 (2*ip  ) = Dat9 (2*ip  ) * rhoct2(ip)          ! x
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat9,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat9(i)=Dat9(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat10,1,1,1)                !10
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat10(2*ip-1)=Dat10(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat10(2*ip  )=Dat10(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat10,1,1,1)                ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat10(i)=Dat10(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat11,1,1,1)                !11
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat11(2*ip-1)=Dat11(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat11(2*ip  )=Dat11(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat11,1,1,1)                ! x 
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat11(i)=Dat11(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat12,1,1,1)                !12  
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat12(2*ip-1)=Dat12(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat12(2*ip  )=Dat12(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat12,1,1,1)                ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat12(i)=Dat12(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat13,1,1,1)                !13
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat13(2*ip-1)=Dat13(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat13(2*ip  )=Dat13(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat13,1,1,1)                ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat13(i)=Dat13(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat14,1,1,1)                !14
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat14(2*ip-1)=Dat14(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat14(2*ip  )=Dat14(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat14,1,1,1)                ! x 
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat14(i)=Dat14(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat15,1,1,1)                !15
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat15(2*ip-1)=Dat15(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat15(2*ip  )=Dat15(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat15,1,1,1)                ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat15(i)=Dat15(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat16,1,1,1)                !16
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat16(2*ip-1)=Dat16(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat16(2*ip  )=Dat16(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat16,1,1,1)                ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat16(i)=Dat16(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat17,1,1,1)                !17
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat17(2*ip-1)=Dat17(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat17(2*ip  )=Dat17(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat17,1,1,1)                ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat17(i)=Dat17(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, dat18,1,1,1)                !18
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
           Dat18(2*ip-1)=Dat18(2*ip-1)*(rhocl2(ip)-2.0*rhoct2(ip))
           Dat18(2*ip  )=Dat18(2*ip  )*(rhocl2(ip)-2.0*rhoct2(ip))
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, dat18,1,1,1)                ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Dat18(i)=Dat18(i)/(Ni*Nj*Nk)                            ! x
      Enddo                                                     ! x

      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            ii=i
            jj=j
            kk=k
            If(ii.gt.Ni/2) ii=ii-Ni
            If(jj.gt.Nj/2) jj=jj-Nj
            If(kk.gt.Nk/2) kk=kk-Nk
            Gx = ii*b1x + jj*b2x + kk*b3x     ! Reciprocal Vectors
            Gy = ii*b1y + jj*b2y + kk*b3y     ! Reciprocal Vectors
            Gz = ii*b1z + jj*b2z + kk*b3z     ! Reciprocal Vectors
            Gx = vkx + Gx                                    ! k+G
            Gy = vky + Gy                                    ! k+G
            Gz = vkz + Gz                                    ! k+G
            hm1= Ureal * Dat1(2*ip-1) + Uimag * Dat1(2*ip)
            hm2= Ureal * Dat2(2*ip-1) + Uimag * Dat2(2*ip)
            hm3= Ureal * Dat3(2*ip-1) + Uimag * Dat3(2*ip)
            hm4= Ureal * Dat4(2*ip-1) + Uimag * Dat4(2*ip)
            hm5= Ureal * Dat5(2*ip-1) + Uimag * Dat5(2*ip)
            hm6= Ureal * Dat6(2*ip-1) + Uimag * Dat6(2*ip)
            hm7= Ureal * Dat7(2*ip-1) + Uimag * Dat7(2*ip)
            hm8= Ureal * Dat8(2*ip-1) + Uimag * Dat8(2*ip)
            hm9= Ureal * Dat9(2*ip-1) + Uimag * Dat9(2*ip)
		
 	    hm10= Ureal * Dat10(2*ip-1)+ Uimag * Dat10(2*ip)
            hm11= Ureal * Dat11(2*ip-1) + Uimag * Dat11(2*ip)
            hm12= Ureal * Dat12(2*ip-1) + Uimag * Dat12(2*ip)
            hm13= Ureal * Dat13(2*ip-1) + Uimag * Dat13(2*ip)
            hm14= Ureal * Dat14(2*ip-1) + Uimag * Dat14(2*ip)
            hm15= Ureal * Dat15(2*ip-1) + Uimag * Dat15(2*ip)
            hm16= Ureal * Dat16(2*ip-1) + Uimag * Dat16(2*ip)
            hm17= Ureal * Dat17(2*ip-1) + Uimag * Dat17(2*ip)
            hm18= Ureal * Dat18(2*ip-1) + Uimag * Dat18(2*ip)

            tempx=(Gx*hm1+Gy*hm2+Gz*hm3)+
     +		  (Gx*hm1+Gy*hm4+Gz*hm7)+
     +	          (Gx*hm10+Gx*hm14+Gx*hm18)

	    tempy=(Gx*hm4+Gy*hm5+Gz*hm6)+
     +            (Gx*hm2+Gy*hm5+Gz*hm8)+
     +	          (Gy*hm10+Gy*hm14+Gy*hm18)

	    tempz=(Gx*hm7+Gy*hm8+Gz*hm9)+
     +            (Gx*hm3+Gy*hm6+Gz*hm9)+
     +	          (Gz*hm10+Gz*hm14+Gz*hm18)

             w(3*ip-2)=tempx
             w(3*ip-1)=tempy
             w(3*ip  )=tempz

          Enddo
        Enddo
      Enddo

C////END SKIP////////////////////////////////////////////////////

C
C     From Reciprocal to Real - *sqrt(rho) - Back to Reciprocal     
C

111      nn(1)=Ni
      nn(2)=Nj
      nn(3)=Nk
      ndim=3

      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            hmx=w(3*ip-2)
            hmy=w(3*ip-1)
            hmz=w(3*ip  )
            Datx (2*ip-1) = Dble( hmx )
            Datx (2*ip  ) = Imag( hmx )
            Daty (2*ip-1) = Dble( hmy )
            Daty (2*ip  ) = Imag( hmy )
            Datz (2*ip-1) = Dble( hmz )
            Datz (2*ip  ) = Imag( hmz )
          Enddo
        Enddo
      Enddo

      isign=1                               ! +1
      Call elasticFFT(isign, edata, datx,1,1,1)
      Do i=0,Ni-1                                               ! x
        Do j=0,Nj-1                                             ! x
          Do k=0,Nk-1                                           ! x
            ip= i*Nj*Nk + j*Nk + k+1                            ! x
            Datx (2*ip-1)=Datx(2*ip-1)/sqrtrhoinv(ip)! x ###
            Datx (2*ip  )=Datx(2*ip  )/sqrtrhoinv(ip)! x ###
          Enddo                                                 ! x
        Enddo                                                   ! x
      Enddo                                                     ! x
      isign=-1                              ! -1                ! x
      Call elasticFFT(isign, edata, datx,1,1,1)                 ! x
      Do i=1,2*Ni*Nj*Nk                                         ! x
        Datx(i)=Datx(i)/(Ni*Nj*Nk)                              ! x
      Enddo                                                     ! x

      isign=1                               ! +1
      Call elasticFFT(isign, edata, daty,1,1,1)
      Do i=0,Ni-1                                               ! y
        Do j=0,Nj-1                                             ! y
          Do k=0,Nk-1                                           ! y
            ip= i*Nj*Nk + j*Nk + k+1                            ! y
            Daty (2*ip-1)=Daty(2*ip-1)/sqrtrhoinv(ip)! y ###
            Daty (2*ip  )=Daty(2*ip  )/sqrtrhoinv(ip)! y ###
          Enddo                                                 ! y
        Enddo                                                   ! y
      Enddo                                                     ! y
      isign=-1                              ! -1                ! y
      Call elasticFFT(isign, edata, daty,1,1,1)                 ! y
      Do i=1,2*Ni*Nj*Nk                                         ! y
        Daty(i)=Daty(i)/(Ni*Nj*Nk)                              ! y
      Enddo                                                     ! y

      isign=1                               ! +1
      Call elasticFFT(isign, edata, datz,1,1,1)
      Do i=0,Ni-1                                               ! z
        Do j=0,Nj-1                                             ! z
          Do k=0,Nk-1                                           ! z
            ip= i*Nj*Nk + j*Nk + k+1                            ! z
            Datz (2*ip-1)=Datz(2*ip-1)/sqrtrhoinv(ip)! z ###
            Datz (2*ip  )=Datz(2*ip  )/sqrtrhoinv(ip)! z ###
          Enddo                                                 ! z
        Enddo                                                   ! z
      Enddo                                                     ! z
      isign=-1                              ! -1                ! z
      Call elasticFFT(isign, edata, datz,1,1,1)                 ! z
      Do i=1,2*Ni*Nj*Nk                                         ! z
        Datz(i)=Datz(i)/(Ni*Nj*Nk)                              ! z
      Enddo                                                     ! z

      Do i=0,Ni-1
        Do j=0,Nj-1
          Do k=0,Nk-1
            ip= i*Nj*Nk + j*Nk + k+1
            hmx= Ureal * Datx(2*ip-1) + Uimag * Datx(2*ip)
            hmy= Ureal * Daty(2*ip-1) + Uimag * Daty(2*ip)
            hmz= Ureal * Datz(2*ip-1) + Uimag * Datz(2*ip)
            w(3*ip-2)=hmx
            w(3*ip-1)=hmy
            w(3*ip  )=hmz
          Enddo
        Enddo
      Enddo

222   Return
      End



