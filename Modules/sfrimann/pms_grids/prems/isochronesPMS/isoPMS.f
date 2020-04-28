c     program bda/progf/iso.f
c
c     calcul des isochrones pour PMS
c
      common /wpl/ wxmv(5000),wbmv(5000),npl,ltitr(10)
      dimension dm(200),dlm(200),dl(200,200),dt(200,200)
     1,da(200,200),xm(5000)
      dimension xlm(5000),basl(200),bast(200),ym(2)
      dimension kl(2),dlog(200,200),va(200),yl(2),yt(2),xmas(200,200)
      dimension dmas(200,200),dda(200,200)
     1,ddl(200,200),ddt(200,200)
      dimension ddmas(200,200),gg(200,200),zlage(200)
      dimension x1des(5000),y1des(5000),x2des(5000),y2des(5000)
      dimension y3des(5000),xmasac(5000)
      dimension x3des(5000),x4des(5000),glog(5000)
      character ltitr*8
      character titf*30
      real masset
      data xmax,xmin,ymax,ymin/-100.,100.,-100.,100./
c*****
      data agl,agt/0.030,0.006/
      data xminl,xmint/0.010,0.002/
      data dseg/0.01/
      data iprim/0/
      data k1,k2/0,1/
c*****
      ipms=1
c*****    
c
c
      open (unit=10,file='modele.dat',
     .status='old')
     
      read(10,400) titf
      read(10,401)
      read(10,402) nma,jmp
c      print*,nma,jmp
      
  400 format(a30)
  401 format(1x)
  402 format(1x,i3,2x,i3)
  
      do 77 ki=1,nma
      read(10,401)
      read(10,403) masset
c      write(6,403) masset
      dm(ki)=masset
      read(10,401)
      do 78 km=1,jmp
      read(10,404) dda(km,ki),ddmas(km,ki),ddl(km,ki),ddt(km,ki)
c      write(6,404) dda(km,ki),ddmas(km,ki),ddl(km,ki),ddt(km,ki)
   78 continue
   77 continue
   
  403 format(1x,f8.3)
  404 format(2x,1pe14.7,0p,f9.4,2f6.3)
   
c
      do 80 m=1,nma
C      write(6,301) dm(m)
      dlm(m)=log10(dm(m))
      do 1 j=1,jmp
      dl(m,j)=ddl(j,m)
      dt(m,j)=ddt(j,m)
      xmas(m,j)=ddmas(j,m)
      da(m,j)=dda(j,m)
      dmas(m,j)=log10(xmas(m,j))
      r2=3.1472+33.583+ddl(j,m)-4.*ddt(j,m)
      gg(m,j)=-7.176+33.299+log10(ddmas(j,m))-r2
      if(j.ne.1) goto 1
      basl(m)=dl(m,j)
      bast(m)=dt(m,j)
    1 dlog(m,j)=log10(da(m,j))
c     print 302,(dl(m,j),dt(m,j),da(m,j),xmas(m,j),gg(m,j),j=1,jmp)
c
c     print 106
   80 continue
c     print 303
      print 106
c
c     Lecture interactive des parametres
c
    9 write(6,120)
  120 format(/,' Give Log age:  ',$)
      read (5,*) aglog
      if(aglog.gt.9.5) then
      agt=0.010
      xmint=0.003
      agl=0.050
      xminl=0.005
      endif
c choix de la masse superieure
      do 600 m=1,nma 
      zlage(m)=dlog(m,jmp)
  600 continue
       deumas=120.0
C      deumas=10**(fipoi(aglog,nma,zlage,dlm))
C      deumas=min(120.0,deumas)
C      deumas=deumas*(1+0.2)
c choix de la masse inferieure
      xaglog=aglog+1
      premas=10**(fipoi(xaglog,nma,zlage,dlm))
C      premas=max(0.8,premas)
      if(ipms.eq.1) premas=0.8
      write(6,602) premas,deumas
  602 format(1x,' inferior mass=',f7.3,' superior mass=',f7.3)
c*****
   
  733 xmpre=premas
      xmvpre=0.
      bmvpre=0.
      age=10.**aglog
      open(15,file='result',status='new')
      i=1
      delm=0.
      xm2w=premas
   65 xm2w=xm2w+delm
      xmm=xm2w  
      if(delm.lt.0) print*,' delm.lt.0'
      if(xmm.gt.deumas) go to 40
      if(i.gt.2000) go to 40
      xlm(i)=log10(xmm)
      xm(i)=xmm     
      m=indice(xmm,dm,nma)
      do 8 k=1,jmp
      va(k)=dlog(m,k)+(dlog(m+1,k)-dlog(m,k))*(xlm(i)-dlm(m))/
     1(dlm(m+1)-dlm(m))
      if(aglog-va(1))11,11,12
   11 yt(1)=dt(m,1)
      yt(2)=dt(m+1,1)
      yl(1)=dl(m,1)
      yl(2)=dl(m+1,1)
      ym(1)=dmas(m,1)
      ym(2)=dmas(m+1,1)
      kl(1)=1
      kl(2)=1
      goto 18
   12 if(k-jmp)14,13,13
   13 if(aglog-va(jmp))14,16,40
   16 diff=1.
      kl(1)=jmp-1
      kl(2)=jmp
      goto 15
   14 if(aglog-va(k))19,8,8
   19 vak=10.**va(k)
      vakm1=10.**va(k-1)
      diff=(age-vakm1)/(vak-vakm1)
      kl(1)=k-1
      kl(2)=k
      goto 15
    8 continue
   15 do 17 l=1,2
      ll=m+l-1
      ym(l)=xmas(ll,k-1)+((xmas(ll,k)-xmas(ll,k-1))*diff)
      ym(l)=log10(ym(l))
      yt(l)=dt(ll,k-1)+((dt(ll,k)-dt(ll,k-1))*diff)
   17 yl(l)=dl(ll,k-1)+((dl(ll,k)-dl(ll,k-1))*diff)
   18 ff=(xlm(i)-dlm(m))/(dlm(m+1)-dlm(m))
      xl=yl(1)+ (yl(2)-yl(1))*ff
      xt=yt(1)+ (yt(2)-yt(1))*ff
      xmlog=ym(1)+(ym(2)-ym(1))*ff
      xmm=10.**xmlog
      glog(i)=-10.6113+xmlog+4.*xt-xl
c***********
      if(ipms.eq.1) then
        call bvmv(xl,xt,xmv,bmv,b2v1,ub,xp)
        goto 51
      endif
      if(aglog.gt.8.8) then
        if(kl(2).le.13) then
        call bvmv(xl,xt,xmv,bmv,b2v1,ub,xp)
        goto 51
      endif
       if(kl(2).eq.14) then
         xt7=dt(m,13)+(dt(m+1,13)-dt(m,13))*ff
         xl7=dl(m,13)+(dl(m+1,13)-dl(m,13))*ff
         xt8=dt(m,14)+(dt(m+1,14)-dt(m,14))*ff
         xl8=dl(m,14)+(dl(m+1,14)-dl(m,14))*ff
         call bvmv(xl7,xt7,xmv7,bmv7,b2v17,ub7,xp7)
         call bvmiii(xl8,xt8,xmv8,bmv8,b2v18,ub8,xp8)
         xmv=xmv7+(xmv8-xmv7)*diff
         bmv=bmv7+(bmv8-bmv7)*diff
         b2v1=b2v17+(b2v18-b2v17)*diff
         ub=ub7+(ub8-ub7)*diff
         xp=xp7+(xp8-xp7)*diff
         go to 51
       endif
        if(kl(2).gt.14) then
         call bvmiii(xl,xt,xmv,bmv,b2v1,ub,xp)
         go to 51
        endif
      endif
c**********
      if(aglog.le.8.8) then
        if(kl(2).le.13) then
        call bvmv(xl,xt,xmv,bmv,b2v1,ub,xp)
        goto 51
      endif
       if(kl(2).eq.14) then
         xt7=dt(m,13)+(dt(m+1,13)-dt(m,13))*ff
         xl7=dl(m,13)+(dl(m+1,13)-dl(m,13))*ff
         xt8=dt(m,14)+(dt(m+1,14)-dt(m,14))*ff
         xl8=dl(m,14)+(dl(m+1,14)-dl(m,14))*ff
         call bvmv(xl7,xt7,xmv7,bmv7,b2v17,ub7,xp7)
         call bvmiii(xl8,xt8,xmv8,bmv8,b2v18,ub8,xp8)
         xmv=xmv7+(xmv8-xmv7)*diff
         bmv=bmv7+(bmv8-bmv7)*diff
         b2v1=b2v17+(b2v18-b2v17)*diff
         ub=ub7+(ub8-ub7)*diff
         xp=xp7+(xp8-xp7)*diff
         go to 51
       endif
        if(kl(2).lt.20) then
         call bvmiii(xl,xt,xmv,bmv,b2v1,ub,xp)
         go to 51
        endif
       if(kl(2).eq.20) then
         xt8=dt(m,19)+(dt(m+1,19)-dt(m,19))*ff
         xl8=dl(m,19)+(dl(m+1,19)-dl(m,19))*ff
         xt9=dt(m,20)+(dt(m+1,20)-dt(m,20))*ff
         xl9=dl(m,20)+(dl(m+1,20)-dl(m,20))*ff
         call bvmiii(xl8,xt8,xmv8,bmv8,b2v18,ub8,xp8)
         call bvmi(xl9,xt9,xmv9,bmv9,b2v19,ub9,xp9)
         xmv=xmv8+(xmv9-xmv8)*diff
         bmv=bmv8+(bmv9-bmv8)*diff
         b2v1=b2v18+(b2v19-b2v18)*diff
         ub=ub8+(ub9-ub8)*diff
         xp=xp8+(xp9-xp8)*diff
         go to 51
       endif
        if(kl(2).gt.20) then
        call bvmi(xl,xt,xmv,bmv,b2v1,ub,xp)
        go to 51
        endif 
      endif
c**********
   51 xmbol=-2.5*xl+4.75
      base=fipoi(xt,nma,bast,basl)
      dev=2.5*(xl-base)
      x2logr=3.1472+33.583+xl-4.*xt
      xlogr=0.5*x2logr
      rsol=(10.**xlogr)/6.9599e+10
      xlogg=-7.176+33.299+xmlog-x2logr
c
c     variables pour le dessin
c
      x1des(i)=xt
      y3des(i)=xl
      y1des(i)=xmbol
      y2des(i)=xmv
      x2des(i)=ub
      x3des(i)=bmv
      x4des(i)=b2v1
      xmasac(i)=10.**xmlog
C      write(6,110) i,xm(i),x1des(i),y1des(i),y2des(i),x2des(i),
C     1x3des(i),glog(i),y3des(i),xmasac(i)
c*****
      n1des=i
      nombre=i
      if(i.eq.1) then
       delm=(dm(m)-dm(m+1))/10
      else
       distl=abs(y3des(i)-y3des(i-1))
       distt=abs(x1des(i)-x1des(i-1))

       if(distl.gt.agl.or.distt.gt.agt) then
        k1=i
        xm2w=xm2w-delm
        delm=delm/1.2
C        if(k1.eq.k2) go to 40
        go to 65
       endif

       if(distl.lt.xminl.and.distt.lt.xmint) then
        k2=i 
        xm2w=xm2w-delm
        delm=1.01*delm
        if((xm2w+delm).gt.deumas) delm=2.*delm/10.
C        if(k1.eq.k2) go to 40
        go to 65
       endif

      endif
      i=i+1
      go to 65
c***** 
    7 continue
   40 print 106
c*****
      do 601 jj=1,nombre
      write(15,110) jj,xm(jj),x1des(jj),y1des(jj),y2des(jj),x2des(jj),
     1x3des(jj),glog(jj),y3des(jj),xmasac(jj)
  601 continue
c
c
  100 format(e11.4,1x,f6.2,1x,f6.2,1x,i5,i3)
  101 format(f6.2,2x,i3)
  102 format(i2,1x,f6.4,2f7.3,f8.3)
c 103 format(1x,i4,2f10.3,2f8.3,1x,2f8.3,2x,2f9.3,2x,i3,' -',i3,
c    13x,f8.3,3x,f8.3,1x,e11.4)
  104 format(1x,'1')
  105 format(1x,' age',2x,e14.7,2x,' nombre',2x,i5,2x,' masse no. 1',2x,
     1f7.3,2x,' masse no. 2 ',f7.3,//)
  106 format(1x,///)
  107 format(8x,' Masse',12x,' log L',3x,' log Te',3x,' M Bol',13x,
     1' log g',4x,' r',8x,' dev'/)
c    1 ' Teff',8x,' xp',6x,' U-B',15x,' Mbol',7x,' dev',3x,' densli'/)
  108 format(1x,i4,f10.3,4f9.3,2i5,3f9.3)
  109 format(1x,2i4)
  110 format(1x,i4,f10.4,8f9.4,f10.4)
  111 format(1x,' log age ',f8.3,//)
  114 format(1x,' bin ', 9x,f10.3,f8.3)
  209 format(1h1,' Donnees des traces evolutifs'/2i4)
  301 format(/,1x,f8.3,/)
  302 format(f8.3,f6.3,e15.7,2f8.3)
  303 format(1h1)
      stop
      end
c
c
      subroutine sortw ( a,n,j )
c
c
c     *******************************************************
c
c          sous - programme de mise en ordre de nombres
c
c     *******************************************************
c
c
c     ce sous - programme est une traduction fortran
c     de l algorithme 271  * quickersort * de l acm
c
c          paul bartholdi        -        mai 1969
c                  observatoire de geneve
c
c        ******************************************
c
      dimension a(j),n(j),lt(20),ut(20)
      integer p,q,ut
c
      jj = j
      i=1
      m=1
    1 if ( j-i-1 ) 12,12,2
    2 p = (i+j)/2
      t = a(p)
      a(p) = a(i)
      it = n(p)
      n(p) = n(i)
      q = j
      ip = i+1
      do 3 k=ip,q
      if ( a(k) - t ) 3,3,4
    4 if ( q-k ) 20,19,19
   19 if ( a(q) - t ) 6,5,5
    6 x = a(k)
      a(k) = a(q)
      a(q) = x
      is = n(k)
      n(k) = n(q)
      n(q) = is
      q = q-1
      go to 3
    5 q = q-1
      go to 4
   20 q = k-1
      go to 8
    3 continue
    8 a(i) = a(q)
      a(q) = t
      n(i) = n(q)
      n(q) = it
      if ( 2*q-i-j ) 10,10,9
    9 lt(m) = i
      ut(m) = q-1
      i = q+1
      go to 11
   10 lt(m) = q+1
      ut(m) = j
      j = q-1
   11 m = m+1
      go to 1
   12 if ( i-j ) 17,14,14
   17 if ( a(i) - a(j) ) 14,14,13
   13 x = a(i)
      a(i) = a(j)
      a(j) = x
      is = n(i)
      n(i) = n(j)
      n(j) = is
   14 m = m-1
      if ( m ) 16,16,15
   15 i = lt(m)
      j = ut(m)
      go to 1
   16 j = jj
      return
      end
c
      function fipoi(x,n,a,b)
c
c     utilise le sous-programme search.
c     sous-programme d interpolation.
      dimension a(n),b(n)
      k=indice(x,a,n)
c     interpolation lineaire.
      fipoi=b(k)+(b(k+1)-b(k))*(x-a(k))/(a(k+1)-a(k))
      return
      end
c
c
      function indice(x0,x,m)
c
c     recherche rapide de la position d une valeur x0 dans une table
c     monotone croissante ou decroissante de m nombres x(i)
c
c     si  k = indice(x0,x,m)  on aura   x0 compris entre x(k) et x(k+1)
c             ou x0 = x(k).
c     si x0 exterieur a la table indice=1 si x0 du cote de x(1)
c                                indice = m-1 si x0 du cote de x(m)
c
      dimension x(m)
c
      n=m
      k=1
c
    1 if(n-k-1) 2,5,2
    2 i=(k+n)/2
      if((x(i)-x0)*(x(n)-x(1))) 4,4,3
    3 n=i
      go to 1
    4 k=i
      go to 1
    5 indice = k
c
      return
c
      end 
c
c
c
      function ran1(idum)
c
      dimension r(97)
      parameter (m1=259200,ia1=7141,ic1=54773,rm1=3.8580247e-6)
      parameter (m2=134456,ia2=8121,ic2=28411,rm2=7.4373773e-6)
      parameter (m3=243000,ia3=4561,ic3=51349)
      data iff /0/
c
      if (idum.lt.0.or.iff.eq.0) then
        iff=1
        ix1=mod(ic1-idum,m1)
        ix1=mod(ia1*ix1+ic1,m1)
        ix2=mod(ix1,m2)
        ix1=mod(ia1*ix1+ic1,m1)
        ix3=mod(ix1,m3)
        do 11 j=1,97
          ix1=mod(ia1*ix1+ic1,m1)
          ix2=mod(ia2*ix2+ic2,m2)
          r(j)=(float(ix1)+float(ix2)*rm2)*rm1
11      continue
        idum=1
      endif
      ix1=mod(ia1*ix1+ic1,m1)
      ix2=mod(ia2*ix2+ic2,m2)
      ix3=mod(ia3*ix3+ic3,m3)
      j=1+(97*ix3)/m3
      if(j.gt.97.or.j.lt.1)pause
      ran1=r(j)
      r(j)=(float(ix1)+float(ix2)*rm2)*rm1
      return
      end
C
C
      SUBROUTINE BVMV(WL,WT,WMV,WBV,WB2V1,WU,WX)
C
C     RELATION Teff vs B-V 
C     D'APRES
C     BOHM-VITENSE:1981,ANN. REV. A&A,19,295
C     CORR. BOLO. D'APRES MALAGNINI ET AL:1986,162,140
C     RELATION UBV SCHMIDT-KALER:1982,LANDOLT & BORNSTEIN
C     verifiee 17-4-89
C  
      DIMENSION ZT(28),ZBC(28),ZBV(28),ZUB(28)
      DATA ZT/ 3.740,3.778,3.813,3.845,3.875,3.903,3.929,3.954,
     1         3.978,4.000,4.041,4.079,4.114,4.146,4.176,4.204,
     2         4.230,4.255,4.301,4.342,4.352,4.380,4.398,4.415,
     3         4.447,4.477,4.653,4.688/   
      DATA ZBC/-0.21,-0.08, 0.00, 0.04, 0.04, 0.02,-0.02,-0.09,
     1         -0.16,-0.25,-0.45,-0.66,-0.87,-1.07,-1.26,-1.43,
     2         -1.58,-1.72,-1.96,-2.17,-2.22,-2.37,-2.46,-2.56,
     3         -2.75,-2.91,-4.00,-4.30/
      DATA ZBV/ 0.70, 0.55, 0.420, 0.32, 0.22, 0.14, 0.10, 0.05,
     1          0.00,-0.03,-0.075,-0.11,-0.13,-0.15,-0.17,-0.18,
     2         -0.19,-0.20,-0.220,-0.24,-0.24,-0.26,-0.26,-0.27,
     3         -0.29,-0.30,-0.360,-0.37/
      DATA ZUB/ 0.23, 0.04,-0.016, 0.02, 0.10, 0.10, 0.09, 0.05,
     1          0.01,-0.06,-0.218,-0.34,-0.43,-0.50,-0.58,-0.62,
     2         -0.67,-0.71,-0.775,-0.84,-0.87,-0.92,-0.95,-0.98,
     3         -1.05,-1.08,-1.300,-1.34/             
      WBV=FIPOI(WT,28,ZT,ZBV)
      CB=FIPOI(WT,28,ZT,ZBC)
      WU=FIPOI(WT,28,ZT,ZUB)
      WMV=-2.5*WL+4.75-CB
      IF(WBV.LE.0.197) THEN
      WB2V1=0.712*WBV-0.144
      ELSE  
      WB2V1=-0.149+0.707*WBV+0.708*WBV*WBV-1.083*WBV*WBV*WBV+
     10.442*WBV*WBV*WBV*WBV
      ENDIF
      RETURN
      END                
C
C
C
      SUBROUTINE BVMIII(WL,WT,WMV,WBV,WB2V1,WU,WX) 
C
C     RELATION Teff vs B-V ET
C     CORR BOL D'APRES FLOWER A A 54,31 (1977)
C     UBV (CLASSE III) SCHMIDT KALER,1982,LANDOLT ET B.
C
      DIMENSION ZT(38),ZBC(38),ZBV(38),ZUB(38)
      DATA ZT/ 3.544,3.574,3.592,3.622,3.650,3.674,3.698,3.706,
     1         3.724,3.747,3.769,3.793,3.813,3.845,3.892,3.914,
     2         3.930,3.944,3.954,3.981,4.000,4.023,4.057,4.099,
     3         4.121,4.211,4.234,4.256,4.301,4.324,4.369,4.414,
     4         4.437,4.482,4.504,4.527,4.549,4.594/
      DATA ZBC/-1.66,-1.19,-0.92,-0.66,-0.49,-0.37,-0.26,-0.22,
     1         -0.16,-0.10,-0.06,-0.03,-0.01, 0.01, 0.00,-0.02,
     1         -0.03,-0.06,-0.09,-0.19,-0.24,-0.32,-0.46,-0.72,
     2         -0.84,-1.34,-1.46,-1.58,-1.84,-1.97,-2.22,-2.48,
     3         -2.60,-2.87,-2.99,-3.12,-3.25,-3.50/
      DATA ZBV/ 1.61, 1.54, 1.45, 1.30, 1.160, 1.04, 0.92, 0.88,
     1          0.80, 0.70, 0.60, 0.50, 0.430, 0.33, 0.20, 0.14,
     1          0.10, 0.07, 0.05, 0.00,-0.025,-0.05,-0.08,-0.11,
     2         -0.12,-0.16,-0.17,-0.18,-0.200,-0.21,-0.23,-0.25,
     3         -0.26,-0.28,-0.29,-0.30,-0.310,-0.31/
      DATA ZUB/ 1.88, 1.84, 1.72, 1.44, 1.160, 0.94, 0.66, 0.59,
     1          0.45, 0.28, 0.16, 0.10, 0.090, 0.08, 0.11, 0.11,
     2          0.10, 0.09, 0.06, 0.03,-0.053,-0.13,-0.24,-0.37,
     3         -0.40,-0.54,-0.58,-0.63,-0.740,-0.78,-0.87,-0.94,
     4         -0.97,-1.04,-1.08,-1.10,-1.120,-1.13/             
      WBV=FIPOI(WT,38,ZT,ZBV)
      CB=FIPOI(WT,38,ZT,ZBC)
      WU=FIPOI(WT,38,ZT,ZUB)
      WMV=-2.5*WL+4.75-CB
      WB2V1=-0.140+0.739*WBV+0.284*WBV*WBV-0.321*WBV*WBV*WBV
     1+0.115*WBV*WBV*WBV*WBV
      RETURN
      END
C
      SUBROUTINE BVMI(WL,WT,WMV,WBV,WB2V1,WU,WX) 
C
C     RELATION Teff vs B-V ET
C     CORR BOL D'APRES FLOWER A A 54,31 (1977)
C     UBV (CLASSE Ia) SCHMIDT KALER,1982,LANDOLT ET B.
C
      DIMENSION ZT(38),ZBC(38),ZBV(38),ZUB(38)
      DATA ZT/ 3.447,3.477,3.512,3.544,3.574,3.588,3.605,3.632,
     1         3.656,3.677,3.698,3.705,3.743,3.766,3.793,3.845,
     1         3.892,3.914,3.930,3.944,3.960,4.011,4.038,4.081,
     2         4.137,4.194,4.213,4.288,4.307,4.325,4.363,4.419,
     3         4.457,4.476,4.513,4.532,4.570,4.607/
      DATA ZBC/-3.36,-2.50,-1.72,-1.43,-1.00,-0.84,-0.67,-0.46,
     1         -0.35,-0.22,-0.14,-0.12,-0.01, 0.04, 0.08, 0.13,
     1          0.14, 0.09, 0.00,-0.10,-0.17,-0.38,-0.51,-0.64, 
     2         -0.82,-1.05,-1.16,-1.56,-1.67,-1.76,-2.02,-2.40,
     3         -2.68,-2.79,-3.06,-3.20,-3.46,-3.73/
      DATA ZBV/ 1.80, 1.76, 1.72, 1.70, 1.61, 1.54, 1.450, 1.30,
     1          1.16, 1.04, 0.92, 0.88, 0.70, 0.60, 0.500, 0.33,
     1          0.20, 0.14, 0.10, 0.07, 0.05, 0.00,-0.025,-0.05,  
     2         -0.08,-0.11,-0.12,-0.16,-0.17,-0.18,-0.200,-0.23, 
     3         -0.25,-0.26,-0.28,-0.29,-0.31,-0.32/
      DATA ZUB/ 2.11, 2.05, 1.99, 1.92, 1.81, 1.71, 1.572, 1.25,
     1          1.06, 0.84, 0.70, 0.66, 0.50, 0.45, 0.392, 0.28,
     2          0.17, 0.12,-0.04,-0.15,-0.20,-0.58,-0.635,-0.70,
     3         -0.76,-0.83,-0.85,-0.96,-0.97,-0.99,-1.013,-1.05,
     4         -1.09,-1.11,-1.13,-1.13,-1.16,-1.17/    
      WBV=FIPOI(WT,38,ZT,ZBV)
      CB=FIPOI(WT,38,ZT,ZBC)
      WU=FIPOI(WT,38,ZT,ZUB)
      WMV=-2.5*WL+4.75-CB
      WB2V1=-0.121+0.722*WBV+0.079*WBV*WBV
      RETURN
      END
