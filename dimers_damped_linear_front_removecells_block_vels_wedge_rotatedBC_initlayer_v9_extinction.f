      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!  Completely overdamped MD for dimers with growth
      !!
      !!
      !!  BC's = rotated according to ETA wedge
      !!  Cells (1) grow in growth layer, (2) are pushed in 
      !!    propagation layer, (3) are held fixed in boundary 
      !!    layer, & (4) are removed beyond boundary layer
      !!
      !!  Depths defined as distance to closest cell in front
      !!  
      !!  Options - Restart: T/F = use restart file/start from scratch
      !!              Movie: T/F = do/do not output movie
      !!             Bottom: T/F = keep/discard bottom
      !!
      !!  F = b*m*dr/dt (m=1 implicit)   
      !!  T = b*I*dth/dt (I=inertia, th=orientation angle)   
      !!
      !!  Carl Schreck
      !!  4/7/2016
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program main
      implicit none
      integer Ntot
      parameter(Ntot=2**18)
      double precision x(Ntot),vx(Ntot),ax(Ntot),bx(Ntot),fx(Ntot),xa(2)
      double precision y(Ntot),vy(Ntot),ay(Ntot),by(Ntot),fy(Ntot),ya(2)
      double precision th(Ntot),vth(Ntot),ath(Ntot),bth(Ntot),fth(Ntot)
      double precision xp(Ntot),yp(Ntot),D(Ntot),alpha(Ntot),rate(Ntot)
      double precision inert(Ntot),depth(Ntot),rate0(Ntot),tdiv,alpha0
      double precision b,Lx,exp,desync,kinetic,KE,V,ran2,cc,ss,corr,ddsq
      double precision layerwidth,layerdepth,dr(2),dk(2),maxdis,alphamax
      double precision width,propdepth,bounddepth,propdist,bounddist,dd
      double precision dt,att,rateWT,w0,s1,dwdleta,wedgedepth,initdist
      double precision angleBC,yBC,xarot(2,2),yarot(2,2),sina,cosa,angle
      double precision bounddepthprod,propdepthprod,minyfront
      integer N,seed,steps,i,j,k,countn(3),nl(3,12*Ntot,2),kstart
      integer restartexist,dataskip,prodskip,div,layerskip,restskip
      integer nrem,forcelist(Ntot),proplist(Ntot),nprop,nsum,Nc,Ncb
      integer Np,Npb,Nf,Nfb,Nu,Nub,Nc2,Ncb2,celltype(Ntot),burnsteps
      integer mutate,rot,l,cut,cutrelax,cutsteps,label(Ntot),nforce
      integer seedstart,nmuts,nlayer
      character file1*199,file2*199,file3*199,file4*199
      logical restart,movie,bottom,bottomprod,boundlogic
      common /f1com/ exp,alpha
      common /f2com/ nl,countn 
      common /f3com/ proplist
      common /f4com/ bottom
      common /f5com/ alphamax,alpha0

      ! read geometric parameters
      read(*,*) alpha0
      read(*,*) alphamax
      read(*,*) Lx
      read(*,*) att

      ! read rates
      read(*,*) rateWT
      read(*,*) b

      ! read steps
      read(*,*) steps
      read(*,*) burnsteps
      read(*,*) layerskip
      read(*,*) dataskip
      read(*,*) prodskip
      read(*,*) restskip
      read(*,*) dt
 
      ! read growth layer parameters
      read(*,*) layerwidth      
      read(*,*) layerdepth

      ! read layer parameters for force calc
      read(*,*) propdepth
      read(*,*) bounddepthprod

      ! read run parameters
      read(*,*) desync
      read(*,*) seed
      
      ! read output files
      read(*,*) file1
      read(*,*) file2
      read(*,*) file3
      read(*,*) file4

      ! read options
      read(*,*) movie
      read(*,*) restart
      read(*,*) bottomprod

      ! rescue parameterss
      read(*,*) s1
      read(*,*) w0

      ! time-steps for cutting out front
      read(*,*) cutrelax

      ! parameters
      exp=2d0     ! 2 = LS, 2.5 = Hertzian, >2.9 = RLJ
      width=0.1d0 ! width of neighborlist 

      ! calculate # steps until division
      tdiv=dlog10(2d0)/dlog10(1d0+dt*rateWT)*dt

      ! calc vertical extent of wedge
      dwdleta=2d0*sqrt(-s1*(2d0+s1))/(1d0+s1)
      wedgedepth=dwdleta*(Lx-w0)/4d0 ! don't need

      ! calc position & rotation angle of BC's
      angleBC=2d0*datan(dwdleta/2d0)
      yBC=wedgedepth+Lx/dwdleta

      ! initialize system from scratch or restart file
      call initialize(file1,file2,file3,file4,restart,movie,
     +     rateWT,rate0,desync,kstart,seed,b,Lx,att,N,rate,
     +     depth,inert,d,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +     w0,celltype,s1,mutate,dwdleta,angleBC,yBC,width,xp,yp)

      ! total # of cells (including those removed)
      nsum=N

      ! initialize boundary for cutting loop
      bottom=.true. ! initially keep bottom 
      bounddepth=bounddepthprod  
      propdepth=propdepthprod

      ! distances of propagation/boundary layer from front
      propdist=layerdepth+propdepth
      bounddist=layerdepth+propdepth+bounddepth
      
      ! loop over time-steps - initialize
      k=0
      cut=0
      do while(cut.le.5)
         k=k+1

         call calc_boundary(N,x,y,Lx,depth,angleBC,yBC,
     +        layerwidth,bounddist,boundlogic)
  
         if(cut.eq.0.and.boundlogic) then
            cut=1
            cutsteps=k+cutrelax
            call cut_wedge(N,x,y,th,vx,vy,vth,ax,ay,ath,
     +           bx,by,bth,d,alpha,depth,rate,rate0,celltype,label,
     +           inert,proplist,bounddist,dwdleta,w0,desync,rateWT,s1)
            call makelist(N,x,y,d,Lx,xp,yp,width,att,angleBC,yBC)
            call calcdepth_wedge(N,x,y,d,Lx,
     +           layerwidth,depth,angleBC,yBC)
         elseif(cut.ge.1.and.k.eq.cutsteps) then
            cut=cut+1
            cutsteps=k+cutrelax
            bottom=bottomprod ! set bottom for producion
            call cut_wedge(N,x,y,th,vx,vy,vth,ax,ay,ath,
     +           bx,by,bth,d,alpha,depth,rate,rate0,celltype,label,
     +           inert,proplist,bounddist,dwdleta,w0,desync,rateWT,s1)
            call makelist(N,x,y,d,Lx,xp,yp,width,att,angleBC,yBC)
            call calcdepth_wedge(N,x,y,d,Lx,
     +           layerwidth,depth,angleBC,yBC)
         endif

         ! grow/divide cells
         call grow(dt,N,nsum,depth,layerdepth,
     +        rate,rate0,rateWT,Lx,width,att,D,x,y,th,vx,vy,vth,
     +        ax,ay,ath,bx,by,bth,xp,yp,seed,desync,celltype,label,s1)

         ! calc propagation list
         call calc_proplist(N,nprop,depth,proplist,
     +        vx,vy,vth,ax,ay,ath,bx,by,bth,propdist)
         
         ! remove cells & make neighbor list
         call checklist(N,x,y,xp,yp,maxdis)
         if(maxdis.gt.width*d(1)) then
            call makelist(N,x,y,d,Lx,xp,yp,width,att,angleBC,yBC)
         endif         

         ! calculate inertia of each cell
         call calc_inert(N,inert,D)
         
         ! Gear precictor-corrector
         call predict(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth)
         call force(N,x,y,th,d,V,fx,fy,fth,Lx,att,angleBC,yBC)    
         call correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,
     +        bx,by,bth,fx,fy,fth,inert,b)
         KE=kinetic(N,vx,vy,vth,inert)
         
         ! calc distance to front     
         if(mod(k,layerskip).eq.0) then
            call calcdepth_wedge(N,x,y,d,Lx,
     +           layerwidth,depth,angleBC,yBC)
         endif
         
         ! output data to screen
         if(mod(k,dataskip).eq.0) then
            nmuts=0
            nlayer=0
            do i=1,N
               if(depth(i).lt.layerdepth) then 
                  nlayer=nlayer+1
                  if(celltype(i).eq.1) then
                     nmuts=nmuts+1
                  endif
               endif
            enddo
            write(*,'(I,ES16.8,4I)') cut,dble(k)*dt,nsum,N,nlayer,nmuts
         endif   
      enddo

      ! barcode cells
      do i=1,N
         label(i)=i
      enddo

      ! loop over time-steps - production
      k=kstart-1
      minyfront=0d0
      do while(minyfront.lt.yBC) ! until front surpases yBC
         k=k+1

         ! grow/divide cells
         call grow(dt,N,nsum,depth,layerdepth,rate,rate0,
     +        rateWT,Lx,width,att,D,x,y,th,vx,vy,vth,ax,ay,ath,
     +        bx,by,bth,xp,yp,seed,desync,celltype,label,s1)

         ! calc propagation list
         call calc_proplist(N,nprop,depth,proplist,
     +        vx,vy,vth,ax,ay,ath,bx,by,bth,propdist)
         
         ! remove cells & make neighbor list
         call checklist(N,x,y,xp,yp,maxdis)
         if(maxdis.gt.width*d(1)) then
            call remove(N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,d,
     +           alpha,depth,rate,rate0,celltype,label,inert,proplist,
     +           bounddist)
            call makelist(N,x,y,d,Lx,xp,yp,width,att,angleBC,yBC)
         endif
         
         ! calculate inertia of each cell
         call calc_inert(N,inert,D)
         
         ! Gear precictor-corrector
         call predict(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth)
         call force(N,x,y,th,d,V,fx,fy,fth,Lx,att,angleBC,yBC)    
         call correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,
     +        bx,by,bth,fx,fy,fth,inert,b)
         KE=kinetic(N,vx,vy,vth,inert)
         
         ! calc distance to front     
         if(mod(k,layerskip).eq.0) then
            call calcdepth_wedge(N,x,y,d,Lx,
     +           layerwidth,depth,angleBC,yBC)
         endif
         
         ! output data to screen
         if(mod(k,dataskip).eq.0) then
            nmuts=0
            nlayer=0
            do i=1,N
               if(depth(i).lt.layerdepth) then 
                  nlayer=nlayer+1
                  if(celltype(i).eq.1) then
                     nmuts=nmuts+1
                  endif
               endif
            enddo
            write(*,'(I,ES16.8,4I)') cut,dble(k)*dt,nsum,N,nlayer,nmuts
         endif         
 
         ! save config
         if(movie.and.mod(k,prodskip).eq.0) then
            write(1,*) 3*2*N
            do i=1,N
               cc=dcos(th(i))
               ss=dsin(th(i))
               dd=alpha(i)-1d0
               dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
               dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
               do j=1,2
                  xa(j)=x(i)+dr(j)*cc
                  ya(j)=y(i)+dr(j)*ss
               enddo
               dk(1)=D(i)
               dk(2)=dd*D(i)

               do rot=1,3
                  angle=angleBC*dble(rot-2)  
                  cosa=dcos(angle)
                  sina=dsin(angle)
                  xarot(rot,1)=xa(1)*cosa-(ya(1)-yBC)*sina
                  yarot(rot,1)=xa(1)*sina+(ya(1)-yBC)*cosa+yBC
                  xarot(rot,2)=xa(2)*cosa-(ya(2)-yBC)*sina
                  yarot(rot,2)=xa(2)*sina+(ya(2)-yBC)*cosa+yBC
                  write(1,'(4F,3I)') xarot(rot,1),yarot(rot,1),
     +                 dk(1),depth(i),rot,celltype(i),label(i)
                  write(1,'(4F,3I)') xarot(rot,2),yarot(rot,2),
     +                 dk(2),depth(i),rot,celltype(i),label(i)
               enddo               
            enddo
            flush(1)
         endif   

         ! save restart file
         if(restart.and.mod(k,restskip).eq.0) then
            open(unit=2,file=TRIM(file2))
            write(2,'(4I)') k, N, seed, mutate
            do i=1,N
               write(2,'(17E26.18,I)') x(i),y(i),th(i),vx(i),vy(i),
     +              vth(i),ax(i),ay(i),ath(i),bx(i),by(i),bth(i),d(i),
     +           alpha(i),depth(i),rate(i),rate0(i),celltype(i),label(i)
            enddo
            flush(2)
            close(2)
         endif   

         if(mod(k,dataskip).eq.0) then
            write(4,*) nprop, V/dble(nprop) ! replace w growth layer
            flush(4)
         endif

         ! calc min y position of cell in front
         minyfront=1d16
         do i=1,N
            if(depth(i).lt.1d0.and.y(i).lt.minyfront) then
               minyfront=y(i)
            endif
         enddo
      enddo
      
      end ! end main

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!   cut out wedge for initialization   !!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine cut_wedge(N,x,y,th,vx,vy,vth,ax,ay,ath,
     +     bx,by,bth,d,alpha,depth,rate,rate0,celltype,label,
     +     inert,proplist,bounddist,dwdleta,w0,desync,rateWT,s1)
      integer Ntot
      parameter(Ntot=2**18)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot)
      double precision by(Ntot),bth(Ntot),inert(Ntot),depth(Ntot)
      double precision d(Ntot),alpha(Ntot),rate(Ntot),w0,rate0(Ntot)
      double precision dwdleta,ywedge,desync,rateWT,s1
      integer N,seed,i,proplist(Ntot),celltype(Ntot),nrem,label(Ntot)
      character file1*199,file2*199,file3*199,file4*199
      logical restart,movie

      ! cut out wedge
      nrem=0
      do i=N,1,-1
         if(dabs(x(i)).lt.w0/2d0) then         
            ywedge=0d0
         else
            ywedge=(dabs(x(i))-w0/2d0)*dwdleta/2d0
         endif

         if(y(i).gt.ywedge) then
            nrem=nrem+1
            do j=i+1,N
               x(j-1)=x(j)
               y(j-1)=y(j)
               th(j-1)=th(j)
               vx(j-1)=vx(j)
               vy(j-1)=vy(j)
               vth(j-1)=vth(j)
               ax(j-1)=ax(j)
               ay(j-1)=ay(j)
               ath(j-1)=ath(j)
               bx(j-1)=bx(j)
               by(j-1)=by(j)
               bth(j-1)=bth(j)
               d(j-1)=d(j)
               alpha(j-1)=alpha(j)
               depth(j-1)=depth(j)
               rate(j-1)=rate(j)
               rate0(j-1)=rate0(j)
               celltype(j-1)=celltype(j)
               label(j-1)=label(j)
               inert(j-1)=inert(j)
               proplist(j-1)=proplist(j)                     
            enddo
         endif
      enddo
      N=N-nrem
      
      ! assign cells types consistent with wedge
      do i=1,N
         if (dabs(x(i)).lt.w0/2d0-dwdleta/2d0*y(i)-0.5d0) then
            celltype(i)=1
            rate0(i)=(1d0+s1)*rateWT
         else
            celltype(i)=0
            rate0(i)=rateWT
         endif
         rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i)
      enddo
      
      return
      end ! cut wedge


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!  initialize cell position & momenta  !!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine initialize(file1,file2,file3,file4,restart,movie,
     +     rateWT,rate0,desync,kstart,seed,b,Lx,att,N,rate,depth,
     +     inert,d,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,w0,celltype,
     +     s1,mutate,dwdleta,angleBC,yBC,width,xp,yp)
      integer Ntot
      parameter(Ntot=2**18)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision Lx,x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot)
      double precision by(Ntot),bth(Ntot),fx(Ntot),fy(Ntot),fth(Ntot)
      double precision inert(Ntot),depth(Ntot),d(Ntot),alpha(Ntot),V
      double precision rate(Ntot),exp,att,b,dd,ddsq,tmp,w0,rateWT,ran2
      double precision rate0(Ntot),alpha0,alphamax,desync,s1,dwdleta
      double precision angleBC,yBC,width,xp(Ntot),yp(Ntot)
      integer N,N1,N2,kstart,seed,restartexist,seedstart,proplist(Ntot)
      integer i,celltype(Ntot),mutate,countn(3),nl(3,12*Ntot,2)
      character file1*199,file2*199,file3*199,file4*199
      logical restart,movie
      common /f1com/ exp,alpha
      common /f2com/ nl,countn 
      common /f3com/ proplist
      common /f5com/ alphamax,alpha0
      
      ! check if restart file exists
      inquire(file=file2,exist=restartexist)
      if(restart.and.restartexist) then  
         ! open files
         if(movie) open(unit=1,file=TRIM(file1),ACCESS="APPEND")
         if(movie) open(unit=4,file=TRIM(file4),ACCESS="APPEND")
         if(restart) open(unit=2,file=TRIM(file2))

         ! read restart file
         read(2,*) kstart, N, seedstart, mutate
         do i=1,N
            read(2,*) x(i),y(i),th(i),vx(i),vy(i),vth(i),
     +           ax(i),ay(i),ath(i),bx(i),by(i),bth(i),d(i),
     +           alpha(i),depth(i),rate(i),rate0(i),celltype(i)
         enddo
         flush(2)
         close(2)

         ! calculate inertia of each cell
         call calc_inert(N,inert,D)
 
         ! burn seeds
         do while (seed.ne.seedstart)
            tmp=ran2(seed)
         enddo
      else ! no restart file exists
         kstart=0
         mutate=0

         ! open files
         open(unit=1,file=TRIM(file1))
         open(unit=3,file=TRIM(file3))
         open(unit=4,file=TRIM(file4))

         ! random initial config - inside wedge
         N1=floor(w0/2d0)
         do i=1,N1
            d(i)=1d0
            x(i)=2d0*dble(i)-dble(N1+1)
            y(i)=0d0
            th(i)=(ran2(seed)-0.5d0)*2d0*pi
            depth(i)=0d0
            proplist(i)=1
         enddo

         ! random initial config - outside wedge
         N2=floor((Lx-w0)/4d0*dsqrt(1d0+dwdleta**2/4d0))
         do i=N1+1,N1+N2 
            d(i)=1d0
            x(i)=-0.5d0+w0/2d0+2d0*dble(i-N1)/dsqrt(1d0+dwdleta**2/4d0)
            y(i)=2d0*dble(i-N1)*dwdleta/2d0/dsqrt(1d0+dwdleta**2/4d0)
            th(i)=(ran2(seed)-0.5d0)*2d0*pi
            depth(i)=0d0
            proplist(i)=1
         enddo
         do i=N1+N2+1,N1+2*N2+1 ! +1: add extra cell
            d(i)=1d0
           x(i)=0.5d0-w0/2d0-2d0*dble(i-N1-N2)/dsqrt(1d0+dwdleta**2/4d0)
            y(i)=2d0*dble(i-N1-N2)*dwdleta/2d0/dsqrt(1d0+dwdleta**2/4d0)
            th(i)=(ran2(seed)-0.5d0)*2d0*pi
            depth(i)=0d0
            proplist(i)=1
         enddo
         N=N1+2*N2+1 ! +1: add extra cell

         ! initial growth rates
         do i=1,N
            if(dabs(x(i)).gt.w0/2d0) then
               celltype(i)=0
               rate0(i)=rateWT               
            else
               celltype(i)=1
               rate0(i)=(1d0+s1)*rateWT
            endif
         enddo

         ! assign initial aspect ratios & rates
         do i=1,N     
            alpha(i)=alpha0*(1d0+dble(i-1)/2d0)
            rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i)
         enddo
         
         ! assign initial aspect ratios & rates
         do i=1,N     
            alpha(i)=alpha0*(1d0+ran2(seed))
         enddo         
         
         ! calculate inertia of each cell
         call calc_inert(N,inert,D)
         call makelist(N,x,y,d,Lx,xp,yp,width,att,angleBC,yBC)
         call force(N,x,y,th,d,V,fx,fy,fth,Lx,att,angleBC,yBC) 
         do i=1,N
            vx(i)=b*fx(i)
            vy(i)=b*fy(i)
            vth(i)=b*fth(i)/inert(i)
            ax(i)=0d0
            ay(i)=0d0
            ath(i)=0d0
            bx(i)=0d0
            by(i)=0d0
            bth(i)=0d0         
         enddo
      endif

      return
      end ! end initialize


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!  initialize cell position & momenta  !!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine grow(dt,N,nsum,depth,layerdepth,rate,
     +     rate0,rateWT,Lx,width,att,D,x,y,th,vx,vy,vth,
     +     ax,ay,ath,bx,by,bth,xp,yp,seed,desync,celltype,label,s1)
      integer Ntot
      parameter(Ntot=2**18)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision exp,Lx,x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot),bx(Ntot),dt
      double precision by(Ntot),bth(Ntot),depth(Ntot),d(Ntot),rate(Ntot)
      double precision alpha(Ntot),rate0(Ntot),att,ran2,alpha0,alphamax
      double precision layerdepth,corr,desync,s1,rateWT
      integer N,nsum,seed,i,celltype(Ntot),label(Ntot)
      common /f1com/ exp,alpha
      common /f5com/ alphamax,alpha0

      do i=1,N
         ! grow cell i
         if(depth(i).lt.layerdepth) then
            corr=(1d0+(alpha(i)-1d0)**2)/2d0/(alpha(i)-1d0)
            alpha(i)=alpha(i)+corr*dt*rate(i)            
         endif

         ! divide cell i
         if(alpha(i).gt.alphamax) then
            ! divide into 2 - N=current cels, nsum=total 
            N=N+1
            nsum=nsum+1

            ! mutate
            if(celltype(i).eq.0) then         
               celltype(i)=0           
               celltype(N)=0           
               rate0(i)=rateWT
               rate0(N)=rateWT
            else if(celltype(i).eq.1) then
               celltype(i)=1
               celltype(N)=1
               rate0(i)=(1+s1)*rateWT
               rate0(N)=(1+s1)*rateWT
            endif

            ! divide into 2 - 1st assigned index N+1
            D(N)=D(i)
            x(N)=x(i)+alpha0/2d0*dcos(th(i))
            y(N)=y(i)+alpha0/2d0*dsin(th(i))
            th(N)=th(i)
            rate(N)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(N)
            alpha(N)=alpha0
            vx(N)=vx(i)
            vy(N)=vy(i)
            vth(N)=vth(i)               
            ax(N)=ax(i)
            ay(N)=ay(i)
            ath(N)=ath(i)               
            bx(N)=bx(i)
            by(N)=by(i)
            bth(N)=bth(i)
            label(N)=label(i)

            ! divide into 2 - 2nd assigned index i
            x(i)=x(i)-alpha0/2d0*dcos(th(i))
            y(i)=y(i)-alpha0/2d0*dsin(th(i))
            rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0(i) 
            alpha(i)=alpha0

            th(N)=th(i)
            th(i)=th(i)+pi

            ! update neighbor list
            call makelistind(N,N,x,y,d,Lx,xp,yp,
     +           width,att,angleBC,yBC)

            ! update depth
            depth(N)=depth(i)+alpha0/2d0*dsin(th(i))
            depth(i)=depth(i)-alpha0/2d0*dsin(th(i))
         endif
      enddo

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!  remove cells that fall behind front  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine remove(N,x,y,th,vx,vy,vth,ax,ay,ath,
     +     bx,by,bth,d,alpha,depth,rate,rate0,celltype,label,
     +     inert,proplist,bounddist)
      integer Ntot
      parameter(Ntot=2**18)
      double precision x(Ntot),y(Ntot),th(Ntot),vx(Ntot),vy(Ntot)
      double precision vth(Ntot),ax(Ntot),ay(Ntot),ath(Ntot)
      double precision bx(Ntot),by(Ntot),bth(Ntot),depth(Ntot)
      double precision d(Ntot),rate(Ntot),alpha(Ntot),inert(Ntot)
      double precision bounddist,rate0(Ntot)
      integer N,nrem,i,j,proplist(Ntot),celltype(Ntot),label(Ntot)

      nrem=0
      do i=N,1,-1
         if(depth(i).gt.bounddist) then
            nrem=nrem+1
            do j=i+1,N
               x(j-1)=x(j)
               y(j-1)=y(j)
               th(j-1)=th(j)
               vx(j-1)=vx(j)
               vy(j-1)=vy(j)
               vth(j-1)=vth(j)
               ax(j-1)=ax(j)
               ay(j-1)=ay(j)
               ath(j-1)=ath(j)
               bx(j-1)=bx(j)
               by(j-1)=by(j)
               bth(j-1)=bth(j)
               d(j-1)=d(j)
               alpha(j-1)=alpha(j)
               depth(j-1)=depth(j)
               rate(j-1)=rate(j)
               rate0(j-1)=rate0(j)
               celltype(j-1)=celltype(j)
               label(j-1)=label(j)
               inert(j-1)=inert(j)
               proplist(j-1)=proplist(j)                     
            enddo
         endif
      enddo
      N=N-nrem
         
      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!  calc propagation list  !!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_proplist(N,nprop,depth,proplist,
     +        vx,vy,vth,ax,ay,ath,bx,by,bth,propdist)
      integer Ntot
      parameter(Ntot=2**18)
      double precision vx(Ntot),vy(Ntot),vth(Ntot),ax(Ntot),ay(Ntot)
      double precision ath(Ntot),bx(Ntot),by(Ntot),bth(Ntot),depth(Ntot)
      double precision propdist
      integer N,nprop,i,proplist(Ntot)

      nprop=0
      do i=1,N
         if(depth(i).lt.propdist) then
            proplist(i)=1
            nprop=nprop+1
         else
            proplist(i)=0
            vx(i)=0d0
            vy(i)=0d0
            vth(i)=0d0
            ax(i)=0d0
            ay(i)=0d0
            ath(i)=0d0
            bx(i)=0d0
            by(i)=0d0
            bth(i)=0d0
         endif
      enddo

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  check max displacement to update list  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_inert(N,inert,D)
      integer Ntot
      parameter(Ntot=2**18)
      double precision inert(Ntot),D(Ntot),alpha(Ntot),dd,ddsq,exp
      integer i,N
      common /f1com/ exp,alpha
      
      do i=1,N
         dd=alpha(i)-1d0
         ddsq=dd*dd
         inert(i)=((1d0+ddsq**2)/(1d0+ddsq)+2d0*ddsq*
     +        (1d0+dd)**2/(1d0+ddsq)**2)*d(i)**2/8d0
      enddo

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  check max displacement to update list  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine checklist(N,x,y,xp,yp,maxdis)
      integer Ntot
      parameter(Ntot=2**18)
      double precision maxdis,x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),df
      integer N

      df=2d0

      maxdis=0d0
      do i=1,N
	maxdis=max(dabs(x(i)-xp(i)),maxdis)
	maxdis=max(dabs(y(i)-yp(i)),maxdis)
      enddo
      maxdis=2d0*dsqrt(df*maxdis*maxdis)

      return
      end ! end checklist

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!   make neighbor list   !!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine makelist(N,x,y,d,Lx,xp,yp,width,att,angleBC,yBC)
      integer Ntot
      parameter(Ntot=2**18)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),d(Ntot)
      double precision rij,dij,rijsq,width,di_up(Ntot),alphamax
      double precision alpha0,att,Lx,Ly,xij,yij,dijlist,yBC,angle
      double precision angleBC,cosa,sina,xrot(Ntot,3),yrot(Ntot,3)
      integer countn(3),nl(3,12*Ntot,2),N,i,j,rot
      common /f2com/ nl,countn
      common /f5com/ alphamax,alpha0

      do rot=1,3
         angle=angleBC*dble(rot-2)            
         cosa=dcos(angle)
         sina=dsin(angle)
         do j=1,N
            xrot(j,rot)=x(j)*cosa-(y(j)-yBC)*sina
            yrot(j,rot)=x(j)*sina+(y(j)-yBC)*cosa+yBC
         enddo
      enddo

      do rot=1,3
         countn(rot)=0      
      enddo
      do i=1,N
         do j=1,i-1
            do rot=1,3
               xij=x(i)-xrot(j,rot)
               dij=alphamax*d(i) 
               dijlist=dij+(width+att)*d(1)
               if(dabs(xij).lt.dijlist) then
                  yij=y(i)-yrot(j,rot)
                  rijsq=xij*xij+yij*yij
                  if(rijsq.lt.dijlist**2) then
                     countn(rot)=countn(rot)+1
                     nl(rot,countn(rot),1)=i
                     nl(rot,countn(rot),2)=j
                     if(rot.eq.1) then
                        countn(3)=countn(3)+1
                        nl(3,countn(3),1)=j
                        nl(3,countn(3),2)=i
                     elseif(rot.eq.3) then
                        countn(1)=countn(1)+1
                        nl(1,countn(1),1)=j
                        nl(1,countn(1),2)=i
                     endif
                  endif
               endif
            enddo
         enddo
      enddo
            
      do i=1,N
         xp(i)=x(i)
         yp(i)=y(i)
      enddo      

      return
      end ! end makelist      

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!   make neighbor list only for cell i   !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine makelistind(i,N,x,y,d,Lx,xp,yp,
     +     width,att,angleBC,yBC)
      integer Ntot
      parameter(Ntot=2**18)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),d(Ntot)
      double precision rij,dij,rijsq,width,di_up(Ntot),att,Lx,Ly
      double precision xij,yij,dijlist,alphamax,alpha0,angleBC,yBC
      double precision xrot(Ntot,3),yrot(Ntot,3),angle,cosa,sina
      integer countn(3),nl(3,12*Ntot,2),N,i,j,rot
      common /f2com/ nl,countn
      common /f5com/ alphamax,alpha0

      do rot=1,3
         angle=angleBC*dble(rot-2)            
         cosa=dcos(angle)
         sina=dsin(angle)
         do j=1,N
            xrot(j,rot)=x(j)*cosa-(y(j)-yBC)*sina
            yrot(j,rot)=x(j)*sina+(y(j)-yBC)*cosa+yBC
         enddo
      enddo

      do j=1,i-1       
         do rot=1,3
            xij=x(i)-xrot(j,rot)
            dij=alphamax*d(i) 
            dijlist=dij+(width+att)*d(1)
            if(dabs(xij).lt.dijlist) then
               yij=y(i)-yrot(j,rot)
               rijsq=xij*xij+yij*yij
               if(rijsq.lt.dijlist**2) then
                  countn(rot)=countn(rot)+1
                  nl(rot,countn(rot),1)=i
                  nl(rot,countn(rot),2)=j
                  if(rot.eq.1) then
                     countn(3)=countn(3)+1
                     nl(3,countn(3),1)=j
                     nl(3,countn(3),2)=i
                  elseif(rot.eq.3) then
                     countn(1)=countn(1)+1
                     nl(1,countn(1),1)=j
                     nl(1,countn(1),2)=i
                  endif
               end if
            endif
         enddo
      enddo

      xp(i)=x(i)
      yp(i)=y(i)
      
      return
      end ! end makelistind
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!    predicts new positions and velocities    !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine predict(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth)     
      integer Ntot
      parameter(Ntot=2**18)
      integer N,i,proplist(Ntot)
      double precision x(Ntot),y(Ntot),vx(Ntot),vy(Ntot),ax(Ntot)
      double precision ay(Ntot),bx(Ntot),by(Ntot),th(Ntot),vth(Ntot)
      double precision ath(Ntot),bth(Ntot),dt,c1,c2,c3
      common /f3com/ proplist

      c1 = dt
      c2 = c1*dt/2d0
      c3 = c2*dt/3d0

      do i=1,N
         if(proplist(i).eq.1) then 
            x(i) = x(i) + c1*vx(i) + c2*ax(i) + c3*bx(i)
            y(i) = y(i) + c1*vy(i) + c2*ay(i) + c3*by(i)
            th(i) = th(i) + c1*vth(i) + c2*ath(i) + c3*bth(i)         
            vx(i) = vx(i) + c1*ax(i) + c2*bx(i)
            vy(i) = vy(i) + c1*ay(i) + c2*by(i)     
            vth(i) = vth(i) + c1*ath(i) + c2*bth(i)     
            ax(i) = ax(i) + c1*bx(i)
            ay(i) = ay(i) + c1*by(i)
            ath(i) = ath(i) + c1*bth(i)
         endif
      enddo

      end ! end prediction step


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!   corrects prediction   !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine correct(dt,N,x,y,th,vx,vy,vth,ax,ay,ath,bx,by,bth,
     +     fx,fy,fth,inert,b)
      integer Ntot
      parameter(Ntot=2**18)
      integer i,N,proplist(Ntot)
      double precision b,dt,x(Ntot),y(Ntot),vx(Ntot),vy(Ntot),ax(Ntot)
      double precision ay(Ntot),bx(Ntot),by(Ntot),th(Ntot),vth(Ntot)
      double precision ath(Ntot),bth(Ntot),fx(Ntot),fy(Ntot),fth(Ntot)
      double precision inert(Ntot),c1,c2,c3,cg0,cg2,cg3
      double precision gear0,gear2,gear3,corrx,corry,corrth
      common /f3com/ proplist

      gear0 = 3d0/8d0
      gear2 = 3d0/4d0
      gear3 = 1d0/6d0

      c1 = dt
      c2 = c1*dt/2d0
      c3 = c2*dt/2d0

      cg0 = gear0*c1
      cg2 = gear2*c1/c2
      cg3 = gear3*c1/c3

      do i=1,N
         if(proplist(i).eq.1) then 
            vxi = b*fx(i)
            vyi = b*fy(i)
            vthi = b*fth(i)/inert(i)
            corrx = vxi - vx(i)
            corry = vyi - vy(i)
            corrth = vthi - vth(i)
            x(i) = x(i) + cg0*corrx
            y(i) = y(i) + cg0*corry
            th(i) = th(i) + cg0*corrth        
            vx(i) = vxi
            vy(i) = vyi
            vth(i) = vthi
            ax(i) = ax(i) + cg2*corrx
            ay(i) = ay(i) + cg2*corry
            ath(i) = ath(i) + cg2*corrth
            bx(i) = bx(i) + cg3*corrx
            by(i) = by(i) + cg3*corry
            bth(i) = bth(i) + cg3*corrth
         endif
      enddo

      end ! end correction step
            

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!           force            !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine force(N,x,y,th,d,V,fx,fy,fth,Lx,att,angleBC,yBC)
      integer Ntot
      parameter(Ntot=2**18)
      double precision x(Ntot),y(Ntot),th(Ntot),alpha(Ntot),D(Ntot),Lx
      double precision radi_up(Ntot),fx(Ntot),fy(Ntot),fth(Ntot),V,Vij
      double precision f_x,f_y,fc,fr,LJ,LJ0,exp,dij,rij,xij,yij,dij_up
      double precision fact,att,fthi,fthj,rijsq,c(Ntot),s(Ntot),dd,dd2dd
      double precision xa(Ntot,2),ya(Ntot,2),dk(Ntot,2),dr(Ntot,2),di1j1
      double precision angle,angleBC,yBC,xrot(Ntot,3),yrot(Ntot,3)
      double precision xarot(Ntot,3,2),yarot(Ntot,3,2)
      integer i,countn(3),nl(3,12*Ntot,2),N,ki,kj,jj,up,down
      integer proplist(Ntot)
      logical forcelogic,forcej
      common /f1com/ exp,alpha
      common /f2com/ nl,countn
      common /f3com/ proplist     

      do i=1,N
         fx(i)=0d0
         fy(i)=0d0
         fth(i)=0d0
      enddo
      V=0d0

      ! convert to from molecules to atoms
      do i=1,N
         c(i)=dcos(th(i))
         s(i)=dsin(th(i))
         dd=alpha(i)-1d0
         dd2=(1d0+dd)/(1d0+dd**2)*D(i)/2d0
         dr(i,1)=dd2*dd**2
         dr(i,2)=-dd2
         do k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         enddo
         dk(i,1)=D(i)
         dk(i,2)=dd*D(i)
         radi_up(i)=(dk(i,2)-2d0*dr(i,2))/2d0
      enddo

      ! rotate by boundary condition angle
      do rot=1,3
         angle=angleBC*dble(rot-2)  
         cosa=dcos(angle)
         sina=dsin(angle)
         do j=1,N
            xrot(j,rot)=x(j)*cosa-(y(j)-yBC)*sina
            yrot(j,rot)=x(j)*sina+(y(j)-yBC)*cosa+yBC
            do k=1,2
               xarot(j,rot,k)=xa(j,k)*cosa-(ya(j,k)-yBC)*sina
               yarot(j,rot,k)=xa(j,k)*sina+(ya(j,k)-yBC)*cosa+yBC
            enddo
         enddo
      enddo

      ! inter-particle interactions      
      do rot=1,3       
         do k=1,countn(rot)
            i=nl(rot,k,1)
            j=nl(rot,k,2)      
            if(rot.eq.2) then
               if(proplist(i).eq.1.or.proplist(j).eq.1) then    
                  forcelogic=.true.
                  forcej=.true.
               else 
                  forcelogic=.false.
               endif
            else
               if(proplist(i).eq.1) then                
                  forcelogic=.true.
                  forcej=.false.
               else 
                  forcelogic=.false.
               endif
            endif
            if(forcelogic) then
               dij_up=radi_up(i)+radi_up(j)
               xij=x(i)-xrot(j,rot)
               if(dabs(xij).lt.dij_up) then 
                  yij=y(i)-yrot(j,rot)        
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.dij_up*dij_up) then
                     di1j1=(dk(i,1)+dk(j,1))/2d0
                     do ki=1,2
                        do kj=1,2
                           dij=(dk(i,ki)+dk(j,kj))/2d0
                           xij=xa(i,ki)-xarot(j,rot,kj)
                           yij=ya(i,ki)-yarot(j,rot,kj)
                           rijsq=xij**2+yij**2
                           if(rijsq.lt.(dij+att)**2) then
                              rij=dsqrt(rijsq)
                              if(exp.eq.2d0) then
                                 fc=(1d0-rij/dij)/dij     
                                 Vij=(1d0-rij/dij)**2/exp
     +                                -(att/dij)**2/exp
                                 fact=(dij/di1j1)**2
                              elseif(exp.lt.2.9) then
                                 fc=(1d0-rij/dij)/dij     
                                 Vij=(1d0-rij/dij)**2/exp
     +                                -(att/dij)**2/exp
                                 fact=(dij/di1j1)**exp
                              else
                                 LJ=(dij/rij)*(dij/rij)
                                 LJ=LJ*LJ*LJ
                                 LJ0=(dij/(dij+att))**6
                                 fc=1d0/rij*LJ*(LJ-1d0)
                                 Vij=(LJ-1d0)**2-(LJ0-1d0)**2
                                 fact=(dij/di1j1)**2
                              endif                 
                              fr=fc/rij*fact
                              f_x=fr*xij
                              f_y=fr*yij
                              if(proplist(i).eq.1) then
                                 fx(i)=fx(i)+f_x
                                 fy(i)=fy(i)+f_y
                                 fth(i)=fth(i)+dr(i,ki)*
     +                                (c(i)*f_y-s(i)*f_x)
                              endif
                              if(forcej) then
                                 if(proplist(j).eq.1) then
                                    fx(j)=fx(j)-f_x
                                    fy(j)=fy(j)-f_y
                                    fth(j)=fth(j)-dr(j,kj)*
     +                                   (c(j)*f_y-s(j)*f_x)
                                 endif
                              endif
                              V=V+Vij*fact                     
                           endif
                        enddo
                     enddo
                  endif
               endif
            endif
         enddo
      enddo
 
      if(exp .gt. 2.9) then
         do i=1,N
            fx(i)=fx(i)/6d0
            fy(i)=fy(i)/6d0
            fth(i)=fth(i)/6d0 
         enddo         
         V=V/72d0
      endif     

      return							
      end ! end force calc

      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!    calc kinetic energy    !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      function kinetic(N,vx,vy,vth,inert)
      integer Ntot
      parameter(Ntot=2**18)
      integer i,N,proplist(Ntot)
      double precision vx(Ntot),vy(Ntot),vth(Ntot),inert(Ntot),kinetic
      common /f3com/ proplist

      kinetic=0d0
      do i=1,N
         if(proplist(i).eq.1) then         
            kinetic=kinetic+vx(i)**2+vy(i)**2+inert(i)*vth(i)**2   
         endif
      enddo   
      kinetic=kinetic/2d0

      end ! end kinetic energy calc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!    random number generator    !!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END ! end ran2
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  calc num in boundary, binned by x pos  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calc_boundary(N,x,y,Lx,depth,
     +     angleBC,yBC,layerwidth,bounddist,boundlogic)
      integer Ntot
      parameter(Ntot=2**18,LMAX=5000,DEPTHMAX=1000)
      double precision x(Ntot),y(Ntot),Lx,layerwidth,xi,depth(Ntot)
      double precision xrot(Ntot,3),yrot(Ntot,3),cosa,sina,angle
      double precision angleBC,yBC,offset,bounddist
      integer nbound(LMAX),N,numbins,i,bin,rot,nlayer(Ntot)
      logical boundlogic

      do i=1,N
         do rot=1,3
            angle=angleBC*dble(rot-2)  
            cosa=dcos(angle)
            sina=dsin(angle)
            xrot(i,rot)=x(i)*cosa-(y(i)-yBC)*sina
            yrot(i,rot)=x(i)*sina+(y(i)-yBC)*cosa+yBC
         enddo
      enddo

      numbins=nint(2d0*Lx/layerwidth)
      
      ! calc vertical depths
      do bin=1,numbins
         nbound(i)=0
         nlayer(i)=0
      enddo         
      offset=floor(Lx/layerwidth)+1
      do i=1,N
         do rot=1,3
            xi=xrot(i,rot)
            bin=floor(xi/layerwidth)+offset
            if(bin.gt.1.and.bin.le.numbins) then
               if(depth(i).gt.bounddist) then
                  nbound(bin)=nbound(bin)+1
               else if(depth(i).lt.1d0) then
                  nlayer(bin)=nlayer(bin)+1
               endif
            endif
         enddo
      enddo

      boundlogic=.true.
      do bin=1,numbins
         if(nbound(bin).eq.0.and.nlayer(bin).gt.0) then
            boundlogic=.false.
         endif
      enddo

      return
      end ! end depth calc
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!  check max displacement to update list  !!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      subroutine calcdepth_wedge(N,x,y,d,Lx,
     +     layerwidth,depth,angleBC,yBC)
      integer Ntot
      parameter(Ntot=2**18,LMAX=5000,DEPTHMAX=1000)
      double precision x(Ntot),y(Ntot),d(Ntot),Lx,layerwidth,xi,xij,yij
      double precision ymin,ymax,depth(Ntot),depthD(Ntot),depthU(Ntot)
      double precision offset,dx,drsq,angleBC,yBC,angle,cosa,sina
      double precision xrot(Ntot,3),yrot(Ntot,3)
      integer binirot(Ntot,3),npartrot(LMAX,3),ibinrot(LMAX,DEPTHMAX,3)
      integer N,numbins,i,j,bin,dbin,nj,numfrontU,numfrontD,rot
      integer bini(Ntot),npart(LMAX),ibin(LMAX,DEPTHMAX)
      integer frontU(Ntot),frontD(Ntot)
      logical bottom
      common /f4com/ bottom

      do j=1,N
         do rot=1,3
            angle=angleBC*dble(rot-2)  
            cosa=dcos(angle)
            sina=dsin(angle)
            xrot(j,rot)=x(j)*cosa-(y(j)-yBC)*sina
            yrot(j,rot)=x(j)*sina+(y(j)-yBC)*cosa+yBC
         enddo
      enddo

      numbins=nint(2d0*Lx/layerwidth)
      
      ! calc vertical depths
      do bin=1,numbins
         do rot=1,3
            npartrot(bin,rot)=0
         enddo
      enddo         
      offset=floor(Lx/layerwidth)+1
      do i=1,N
         do rot=1,3
            xi=xrot(i,rot)!-dnint(xrot(i,rot)/Lx)*Lx
            bin=floor(xi/layerwidth)+offset
            if(bin.gt.1.and.bin.le.numbins) then
               binirot(i,rot)=bin                    ! x bin of cell i
               npartrot(bin,rot)=npartrot(bin,rot)+1 ! # cells in bin
               ibinrot(bin,npartrot(bin,rot),rot)=i ! ibin=index in bin 
            endif
         enddo
      enddo
      do i=1,N
         ymin=1d16
         ymax=-1d16
         do dbin=-1,1
            bin=mod(numbins+binirot(i,2)+dbin-1,numbins)+1
            do rot=1,3
               do nj=1,npartrot(bin,rot)                  
                  j=ibinrot(bin,nj,rot)               
                  xij=xrot(i,2)-xrot(j,rot)
                  if(dabs(xij)<layerwidth) then
                     if(yrot(j,rot).gt.ymax) then
                        ymax=yrot(j,rot)
                     elseif(yrot(j,rot).lt.ymin) then
                        ymin=yrot(j,rot)
                     endif                     
                  endif
               enddo
            enddo
         enddo
         depthU(i)=ymax-yrot(i,2)
         depthD(i)=yrot(i,2)-ymin 
      enddo

      ! assign cells near front to be at front
      numfrontU=0
      numfrontD=0
      do i=1,N
         if(depthU(i).lt.D(i)) then 
            numfrontU=numfrontU+1
            frontU(numfrontU)=i
            depthU(i)=0d0
         endif
         if(depthD(i).lt.D(i)) then
            numfrontD=numfrontD+1
            frontD(numfrontD)=i
            depthD(i)=0d0
         endif
      enddo
      
      ! calc distance to nearest cell at front
      do i=1,N
         do jj=1,numfrontU
            j=frontU(jj)
            do rot=1,3
               dx=xrot(i,2)-xrot(j,rot)
               if(dabs(dx).lt.depthU(i)) then
                  dy=yrot(i,2)-yrot(j,rot)
                  drsq=dx*dx+dy*dy
                  if(drsq.lt.depthU(i)**2) then
                     depthU(i)=dsqrt(drsq)
                  endif
               endif
            enddo
         enddo
         do jj=1,numfrontD
            j=frontD(jj)
            do rot=1,3
               dx=xrot(i,2)-xrot(j,rot)
               if(dabs(dx).lt.depthD(i)) then
                  dy=yrot(i,2)-yrot(j,rot)
                  drsq=dx*dx+dy*dy
                  if(drsq.lt.depthD(i)**2) then
                     depthD(i)=dsqrt(drsq)
                  endif
               endif
            enddo
         enddo
         if(bottom) then
            depth(i)=min(depthU(i),depthD(i))   
         else
            depth(i)=depthU(i)
         endif
      enddo

      return
      end ! end depth calc

