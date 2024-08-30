program armd
implicit none
!
! Defining variables
!
integer::natom,ntime,nprint,i,j,k,ncont,ilts,itime,idd,istart,NorS
  real,dimension(500)::x,y,z,vx,vy,vz,ax,ay,az  ! velocities, accelerations 
  real::dseed,temp,epsln,sigma,Penergy,Kenergy,Tenergy,expTenergy,xbox,ybox,zbox,mar,kb,rcut,dx,dy,dz,dd,xx,yy, &
        zz,epsln4,epsln24,sigma2,rcut2,Kenergyi,itemp,dt,dth,dt2h,time,scal,reqKenergy,avgPenergy, avgTenergy, &
        avgKenergy, avgtemp, avgPenergy2, avgKenergy2, avgTenergy2, avgtemp2, sdPenergy, sdKenergy, sdTenergy,&
        sdtemp,fqqp,ffqn,r1,ApCharge,AnCharge,fpiEpsln,Avgd
  character(len=20)::gout,eout,xout,ipdb,fpdb
!
  write(6,*)'==================================='
  write(6,*)'...................................'
  write(6,*)'     Calculation is in progess     '
  write(6,*)'...................................'
  write(6,*)'==================================='
! Open input file
!
  open(unit=10,file='armd.dat',status='old')
!
! Read data from the file
!
  read(10,*) istart             ! CDT generate (0) or from the last MD run (1)
  read(10,*) ilts               ! rigid (0) OR random (1)
  read(10,*) natom              ! No. of atoms in the system
  read(10,*) ntime              ! No. of MD steps
  read(10,*) nprint             ! printing interval
  read(10,*) epsln              ! epsilon
  read(10,*) sigma              ! sigma
  read(10,*) temp               ! Initial temperature
  read(10,*) expTenergy         ! Expected Total Energy (kJ/mol)
  read(10,*) mar                ! Mass of Argon (g/mol)
  read(10,*) kb                 ! Boltzmann Constant (kJ/mol/K)
  read(10,*) rcut               ! Cutoff distance
  read(10,*) dt                 ! MD time step (ps)
  read(10,*) xbox               !
  read(10,*) ybox               ! Simulation box dimensions
  read(10,*) zbox               !
  read(10,*) NorS               ! Normal or Salt
  read(10,'(a)')gout            ! Name of the general output
  read(10,'(a)')eout            ! Energy out file
  read(10,'(a)')xout            ! Coordinate output file(trajectory file)
  read(10,'(a)')ipdb            ! Initial CDT file in PDB format)
  read(10,*) fpdb               ! Fianl CDT file in PDB format
!
! Create output files
!
  open(unit=20,file=gout,status='unknown')
  open(unit=21,file=eout,status='unknown')
  open(unit=22,file=xout,status='unknown')
  open(unit=23,file=ipdb,status='unknown')
  open(unit=24,file=fpdb,status='unknown')
!
! write data
!
  write(20,500) istart,ilts,natom,ntime,nprint,epsln,sigma,temp,expTenergy,mar,kb,rcut,dt,xbox,ybox,zbox,NorS,gout,eout, &
                xout,ipdb,fpdb
 500 format(/,2x,' Start or restart option     =',i5,/,&
              2x,' Initial CDT option          =',i5,/,&
              2x,' Number of atoms             =',i5,/,&
              2x,' Number of MD steps          =',i7,/,&
              2x,' Printing interval           =',i5,/,&
              2x,' Ar Epsilon                  =',f8.4,' kJ/mol',/,&
              2x,' Ar Sigma                    =',f8.4,' nm',/,&
              2x,' Initial temperature         =',f8.4,' K',/,&
              2x,' Expected Total Energy       =',f9.2,' kJ/mol',/,&
              2x,' Mass of Ar                  =',f8.4,' g/mol',/,&
              2x,' Boltzmann Constant          =',f8.6,' kJ/mol/K',/,&
              2x,' Cutoff distance             =',f8.4,' nm',/,&
              2x,' MD time step                =',f8.4,' ps',/,&
              2x,' xbox                        =',f8.4,' nm',/,&
              2x,' ybox                        =',f8.4,' nm',/,&
              2x,' zbox                        =',f8.4,' nm',/,&
              2x,' Normal(1) or Salt(0)        =',i5,/,&
              2x,' General output file         =',a20,/,&
              2x,' Energy output file          =',a20,/,&
              2x,' Coordinate out file         =',a20,/,&
              2x,' Initial Coordinates         =',a20,/,&
              2x,' Final Coordinates           =',a20,/)
!
 if (istart == 0) then   ! for start option
!
! -----------Generating initial coordinates
! (generating latice using do loops)
 if (ilts == 0) then
 ncont =0
 do i = 1,6
   xx = (i-1)*0.39
   do j = 1,6
     yy = (j-1)*0.39
     do k = 1,6
       zz = (k-1)*0.39
       ncont = ncont + 1
       x(ncont) = xx
       y(ncont) = yy
       z(ncont) = zz
     end do
   end do
 end do
!
 else
!
! (generating latice randomly)
   ncont = 0
99 do i = 1,200                        ! 
   xx = (rand() - 0.5) * xbox
   yy = (rand() - 0.5) * ybox
   zz = (rand() - 0.5) * zbox
   if (ncont == 0) then
      ncont = ncont + 1
      x(ncont) = xx
      y(ncont) = yy
      z(ncont) = zz
   else
    do j = 1,ncont
     dx = x(j) - xx
     dy = y(j) - yy
     dz = z(j) - zz
     dd = sqrt(dx*dx + dy*dy + dz*dz)
     if (dd<0.38) go to 99              ! To prevent getting very close atoms
    end do
    ncont = ncont + 1
    x(ncont) = xx
    y(ncont) = yy
    z(ncont) = zz
    if (mod(ncont,25) == 0) write(6,'(i5)') ncont
    if (ncont == natom) go to 88
   end if
   end do
  88 continue
 end if
!
 else    ! for start option
!
   open(unit=11, file='finl.pdb',status='old')
   read(11,*)
   read(11,*)
   do i = 1,natom
     read(11,'(30x,3f8.3)') xx,yy,zz
     x(i) = xx*0.10 
     y(i) = yy*0.10  ! convert anstrons into nm 
     z(i) = zz*0.10 
   end do
!
 end if
!
!-------Initial CDT file
!
  do i = 1, natom
    xx = x(i)*10
    yy = y(i)*10  ! 10 is used to convert into anstrons
    zz = z(i)*10
    if(NorS == 1)then
      write(23,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0,0       ! When the system is neutral
    else
      if(i<11) then       
        write(23,501) i,'Ap  ','Ap  ',i,xx,yy,zz,0.0,0.0,1     ! For positive atoms of the salt  
      else if(i<21) then  
        write(23,501) i,'An  ','An  ',i,xx,yy,zz,0.0,0.0,-1    ! For negative atoms of the salt
      else  
        write(23,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0,0     ! For neutral atoms
      end if
    end if
  end do
 501 format('ATOM',2x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2f6.2,4x,2x,i2)
!
!--------- Initial potential energy and forces of the system
!
  epsln4   = 4.0*epsln          
  epsln24  = 24.0*epsln           ! To reduce the computer time
  sigma2   = sigma*sigma
  rcut2 = rcut*rcut
!
  call eng(natom,mar,epsln4,epsln24,sigma2,rcut2,x,y,z,xbox,ybox,zbox,Penergy,ax,ay,az)
!
!---------Initial potential energy and acceleration
!
  write(20,'(/,a,f12.5,a)')  'Initial Potential Energy =', Penergy,' kJ/mol'
  write(20,*) 'Initial Accelerations'
  do i = 1,natom
    write(20,'(i5,2x,3f12.5)') i,ax(i),ay(i),az(i)
  end do
!
!---------Initial Velocities
!
  Kenergyi = 0.0
  do i = 1,natom
    vx(i) = rand() - 0.5    
    vy(i) = rand() - 0.5    
    vz(i) = rand() - 0.5 
    Kenergyi = Kenergyi + mar*0.5*(vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))   
  end do
!
!=================================================================================
!                (Activate this part for temperature calcaulation)
!---------------Calculate the initial kinetic energy and temperature
!---------Rescaling to get desired temperature by temp reading from the input file
!
  itemp = 2.0*Kenergyi/(3*natom*kb)
  Kenergyi = 0.0
  do i = 1,natom
    vx(i) = vx(i) * sqrt(temp/itemp)
    vy(i) = vy(i) * sqrt(temp/itemp)
    vz(i) = vz(i) * sqrt(temp/itemp)
    Kenergyi = Kenergyi + mar*0.5*(vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))   
  end do
  itemp = 2.0*Kenergyi/(3*natom*kb)
!
!===============================================================================
!
!
  write(20,*)''
  write(20,*) ' Initial Kinetic Energy    =', Kenergyi,' kJ/mol'
  write(20,*) ' Initial Temperature       =', itemp,' K'
!
!
!================================================================================
!----------(Activate this part to for expected total energy calculations)
!---------Rescaling velocities to get 'expTenergy' (Expected Total Energy)
!
!  reqKenergy = expTenergy - Penergy       ! reqKenergy - Required Kinetic Energy
!  scal       = sqrt(reqKenergy/Kenergyi)
!  Kenergyi   = 0.0
!  do i = 1,natom
!    vx(i) = vx(i) * scal
!    vy(i) = vy(i) * scal
!    vz(i) = vz(i) * scal
!    Kenergyi = Kenergyi + mar*0.5*(vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))   
!  end do
!  Tenergy = Penergy + Kenergyi
!  itemp = 2.0*Kenergyi/(3*natom*kb)
!
!  write(20,*) ' Initial Potentail Energy =', Penergy,' kJ/mol'
!  write(20,*) ' Initial Kinetic Energy   =', Kenergyi,' kJ/mol'
!  write(20,*) ' Initial Total Energy     =', Tenergy,' kJ/mol'
!  write(20,*) ' Initial Temperature      =', itemp,' K'
!================================================================================
  dt2h        = dt*dt/2.0
  dth         = dt/2.0
  idd         = 0
  avgPenergy  = 0.0
  avgKenergy  = 0.0
  avgTenergy  = 0.0
  avgtemp     = 0.0
  avgPenergy2 = 0.0
  avgKenergy2 = 0.0
  avgTenergy2 = 0.0
  avgtemp2    = 0.0
!
!-------MD run starts here (Using Verlet Algorithm)
!
  do itime = 1,ntime
    do i = 1,natom
      !
      ! New positions
      !
      x(i) = x(i) + dt*vx(i) + dt2h*ax(i)
      y(i) = y(i) + dt*vy(i) + dt2h*ay(i)
      z(i) = z(i) + dt*vz(i) + dt2h*az(i)
      !
      ! Volocities
      !
      vx(i) = vx(i) + dth*ax(i)
      vy(i) = vy(i) + dth*ay(i)
      vz(i) = vz(i) + dth*az(i)
    end do
!
!--------- Periodic boundry conditions
!
    do i = 1,natom
      x(i) = x(i) - anint(x(i)/xbox)*xbox
      y(i) = y(i) - anint(y(i)/ybox)*ybox
      z(i) = z(i) - anint(z(i)/zbox)*zbox
    end do
!
!-----------Call energy
!
  call eng(natom,mar,epsln4,epsln24,sigma2,rcut2,x,y,z,xbox,ybox,zbox,Penergy,ax,ay,az)
!
! ---------Second half of the velocity update and the kinetic energy
!
    Kenergy = 0.0
    do i = 1,natom
      vx(i) = vx(i) + dth*ax(i)
      vy(i) = vy(i) + dth*ay(i)
      vz(i) = vz(i) + dth*az(i)
      Kenergy = Kenergy + mar*0.5*(vx(i)*vx(i) + vy(i)*vy(i) + vz(i)*vz(i))
    end do
!
! Calculating the temperature and averages
!
    itemp       = 2.0*Kenergy/(3.0*natom*kb)
    Tenergy     = Penergy + Kenergy
    time        = itime*dt
    avgPenergy  = avgPenergy + Penergy
    avgKenergy  = avgKenergy + Kenergy
    avgTenergy  = avgTenergy + Tenergy
    avgtemp     = avgtemp    + itemp
    avgPenergy2 = avgPenergy2 + Penergy*Penergy
    avgKenergy2 = avgKenergy2 + Kenergy*Kenergy
    avgTenergy2 = avgTenergy2 + Tenergy*Tenergy
    avgtemp2    = avgtemp2    + itemp*itemp
!        
    if(mod(itime,nprint) == 0) then
      write(6,'(a,i5)') 'Number of steps = ',itime
      idd = idd + 1
      write(21,'(5f13.4)') time,Penergy,Kenergy,Tenergy,itemp
      write(22,'(a,i8)')'MODEL',idd
      write(22,'(a,3f8.3,3f6.1)') 'REMARKS',xbox*10,ybox*10,zbox*10,90.0,90.0,90.0
      do i = 1,natom
        xx = x(i)*10.0
        yy = y(i)*10.0 ! write in anstron unit
        zz = z(i)*10.0
        if(NorS == 1)then
          write(22,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0,0
        else
          if(i<11) then       
            write(22,501) i,'Ap  ','Ap  ',i,xx,yy,zz,0.0,0.0,1
          else if(i<21) then   
            write(22,501) i,'An  ','An  ',i,xx,yy,zz,0.0,0.0,-1
          else  
            write(22,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0,0
          end if
        end if
      end do
      write(22,'(a)')'ENDMDL'
      write(20,503) itime,time,Penergy,Kenergy,Tenergy,itemp
    end if
  end do
    write(6,*)'=========================================='
    write(6,*)'..........................................'
    write(6,*)'        Calculation is done               '
    write(6,*)'    armd.out  ! General output file       '
    write(6,*)'    armd.eng  ! Energy output file        '
    write(6,*)'    armd.cdt  ! Coordinate output file    '
    write(6,*)'    init.pdb  ! Coordinate output file    '
    write(6,*)'    finl.pdb  ! Coordinate output file    '
    write(6,*)'..........................................'
    write(6,*)'=========================================='
!
! -----------------MD run is over now
!
503 format(/,2x,'No. of steps                  = ',i7,/,&
             2x,'Time                          = ',f10.4,' ps',/,&
             2x,'Potential Energy              = ',f12.4,' kJ/mol',/,&
             2x,'Kinetic Energy                = ',f12.4,' kJ/mol',/,&
             2x,'Total Energy                  = ',f12.4,' kJ/mol',/,&
             2x,'Temperature                   = ',f12.4,' K',/)
!
! ----------------Final configuration
!
    write(24,'(a)') 'Final configuration'
    write(24,'(a,3f8.3,3f6.1)') 'REMARKS',xbox*10,ybox*10,zbox*10,90.0,90.0,90.0 ! 10 to make anstron
    do i = 1,natom
      xx = x(i)*10.0
      yy = y(i)*10.0 ! write in anstron unit
      zz = z(i)*10.0
      if(NorS == 1)then
        write(24,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0,0
      else
        if(i<11) then       
          write(24,501) i,'Ap  ','Ap  ',i,xx,yy,zz,0.0,0.0,1
        else if(i<21) then  
          write(24,501) i,'An  ','An  ',i,xx,yy,zz,0.0,0.0,-1
        else  
          write(24,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0,0
        end if
      end if
    end do
    write(24,'(a)')'ENDMDL'
    close(unit=24)        ! if write again it start from the top
!
! ----------------Calculating averages
!

     avgPenergy2  = avgPenergy2/ntime
     avgKenergy2  = avgKenergy2/ntime
     avgTenergy2  = avgTenergy2/ntime
     avgtemp2     = avgtemp2/ntime
     sdPenergy    = sqrt((ntime*avgPenergy2 - avgPenergy*avgPenergy)/(ntime*ntime))
     sdKenergy    = sqrt((ntime*avgKenergy2 - avgKenergy*avgKenergy)/(ntime*ntime))
     sdtemp       = sqrt((ntime*avgtemp2 - avgtemp*avgtemp)/(ntime*ntime))
     avgPenergy   = avgPenergy/ntime
     avgKenergy   = avgKenergy/ntime
     avgTenergy   = avgTenergy/ntime
     avgtemp      = avgtemp/ntime
     sdTenergy    = abs(avgTenergy/ntime)*sqrt((sdPenergy/avgPenergy)**2 + (sdKenergy/avgKenergy)**2)
     write(20,505) avgPenergy,sdPenergy,avgKenergy,sdKenergy,avgTenergy,sdTenergy,avgtemp,sdtemp
505  format(//,2x, 'Average properties -----------------',//,&
               2x,'Potential energy  = ',f12.4,' +/- ',f8.4,' kJ/mol',/,&
               2x,'Kinetic energy    = ',f12.4,' +/- ',f8.4,' kJ/mol',/,&
               2x,'Total energy      = ',f12.4,' +/- ',f8.4,' kJ/mol',/,&
               2x,'Temperature       = ',f12.4,' +/- ',f8.4,' K',/)
!
end program armd
!
! ***************************************************************
!
! Subroutine
!
  subroutine eng(natom,mar,epsln4,epsln24,sigma2,rcut2,x,y,z,xbox,ybox,zbox,Penergy,ax,ay,az)
  implicit none
!
  integer::i,j,k,natom,ii
  real,dimension(500)::x,y,z,ax,ay,az
  real::epsln4,epsln24,sigma2,rcut2,xbox,ybox,zbox,Penergy,c2,c6,c12,r2,xi,yi,zi,fxi,fyi, &
        fzi,dx,dy,dz,ff,mar,NorS,fqqp,fqqn,r1,ApCharge,AnCharge,fpiEpsln,Avgd
!
! ------Initialize
!
! (Pairwise calculation is carried out)
!
  Penergy = 0.0
  do i = 1,natom
    ax(i) = 0.0
    ay(i) = 0.0
    az(i) = 0.0
  end do
!
  do i = 1,natom-1
    xi = x(i)
    yi = y(i)
    zi = z(i)
    fxi = 0.0
    fyi = 0.0
    fzi = 0.0
    ii = i + 1
    do j = ii,natom     ! To save computer time :: Paralell computing
      dx = x(j) - xi
      dy = y(j) - yi
      dz = z(j) - zi
!
!------Minimum image convention
!
      dx = dx - xbox*anint(dx/xbox)
      dy = dy - ybox*anint(dy/ybox)      ! anint :: Will take the nearest integer
      dz = dz - zbox*anint(dz/zbox)
!
      r2 = dx*dx + dy*dy + dz*dz   ! Calculate distance between two atoms(nm^2)            
      if (r2 <= rcut2) then
         c2  = sigma2/r2
         c6  = c2*c2*c2
         c12 = c6*c6
         Penergy  = Penergy + epsln4*(c12 - c6)
!
!------Calculating force due to charges
! 
!=================================================================
         ApCharge =  1.602e-19    ! Positive ion charge in Coulomb
         AnCharge = -1.602e-19    ! Negative ion charge in Coulomb
         fpiEpsln =  1.113e-19    ! 4*Pi*Epsilon0 (C^2/J/nm)-------> (1.113e-10 C^2/J/m)
         Avgd     =  6.022e23     ! Avagardro constant(1/mol)
         !
         ff  = epsln24*(2*c12 - c6)/r2       ! Derivative of Lennard-Jones function (No force from charges) (Calculate before to reduce the computer time)
         !
         if (NorS == 0) then                 ! When the system has a salt (Force from charges are there)
           !
           ! Calculate force due to charges befoe if statemet to reduce the computer time
           ! (Avagardro contant and 10^3 is taken for unit adjustment)
           !
           fqqp = ((ApCharge*ApCharge)/(fpiEpsln*r2))*(Avgd/10**3)     ! +/+ charges and -/- charges  (kJ/mol/nm)
           fqqn = ((ApCharge*AnCharge)/(fpiEpsln*r2))*(Avgd/10**3)     ! +/- charges and -/+ charges
           !
           !        
           if ( i<11 .and. j<11) then 
             ff  = ff+fqqp          
           else if (i<11 .and. j<21) then 
             ff  = ff+fqqn                  
           else if (i<21 .and. j<11) then 
             ff  = ff+fqqn          
           else
             ff  = ff+fqqp          
           end if        
!
         else 
           ff  = ff            
         end if
!         
!=================================================================
!
         fxi = fxi - ff*dx
         fyi = fyi - ff*dy                        ! Forces on atom i due to atom j
         fzi = fzi - ff*dz
         ax(j) = ax(j) + ff*dx
         ay(j) = ay(j) + ff*dy
         az(j) = az(j) + ff*dz
      end if
    end do
    ax(i) = ax(i) + fxi
    ay(i) = ay(i) + fyi
    az(i) = az(i) + fzi
  end do
!
! ----------Convert force to acceleration
!
  do i = 1,natom
    ax(i) = ax(i)/mar
    ay(i) = ay(i)/mar
    az(i) = az(i)/mar
  end do
!
  end subroutine
