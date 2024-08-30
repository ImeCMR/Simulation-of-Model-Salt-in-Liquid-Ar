 program radialdis
 implicit none
!
! Defining variables
!
 integer:: i,j,ii,ibin,nbin,nconfig,iconfig,natom
 real,dimension(500)::x,y,z,bin
 real::distance,rr,xbox,ybox,zbox,xi,yi,zi,dx,dy,dz,binSize,volume,rho,pi,dv,GofR
 character(len=20)::icdt,ordf
!
! Open input data file 
!
 open(unit=10,file='rdf.dat',status='old')
!
! Reading the file
!
 read(10,*) natom               ! Number of atoms
 read(10,*) nconfig             ! Number of configurations
 read(10,*) nbin                ! Bin size
 read(10,*) distance            ! Half of the box length or little bit less
 read(10,'(a20)') icdt          ! Input coordinate file (armd.cdt)
 read(10,'(a20)') ordf          ! Output file name
!
! Open input and output files
!
 open(unit=11,file=icdt,status='old')
 open(unit=20,file=ordf,status='unknown')
!
! Calculating the bin size
!
 binSize = distance/nbin
!
! Intialize all the bins
!
 do i = 1,nbin+5
   bin(i) = 0.0         ! bin is defined as a real array
 end do                 ! (finally we need to have average no of atoms and it can be real)
!
! Read the armd.cdt file
!
 do iconfig = 1,nconfig             !  going through all the configurations
!
! Reading iconfig_th configuration
!
   read(11,*)
   read(11,'(7x,3f8.3)') xbox,ybox,zbox
   do i = 1,natom
     read(11,'(30x,3f8.3)') x(i),y(i),z(i)
   end do
   read(11,*)
!------------------------- Reading is done (iconfig_th configuration)
!==============================================================================================
!------------RDF calculation starts here
!
   do i = 1,20          ! Runs for charged atoms(if all atoms are identical make this natom-1) ***
     xi = x(i)
     yi = y(i)          ! For salt atoms(Solute)
     zi = z(i)
     !ii = i+1          !(If atoms are identical activate this) 
     do j = 21,natom    ! (if all atoms are identical ii+natom)  ***
       dx = x(j) - xi
       dy = y(j) - yi   ! taking the distances from solute atoms to solvent atoms
       dz = z(j) - zi
!
!----------Minimum Image Convension
!
       dx = dx - anint(dx/xbox)*xbox
       dy = dy - anint(dy/ybox)*ybox
       dz = dz - anint(dz/zbox)*zbox
!
!----------Calculating the distances (rr)
!
       rr = sqrt(dx*dx + dy*dy + dz*dz)
       if (rr <= distance) then                ! When rr is less than half a box length
         ibin = int(rr/binSize) + 1            ! Finds what is the bin that should update
         bin(ibin) = bin(ibin) + 1             ! (If all atoms are identical bin(ibin) + 2)  ***
       end if  
     end do
   end do
 end do
!==============================================================================================
! Average number of atoms in the bin
!
 do i = 1,nbin
   bin(i) = bin(i)/nconfig/20
 end do
 volume = xbox*ybox*zbox
 rho = real(natom)/volume
 pi = acos(-1.000)                               ! Cos inverse
 write(6,*) ' PI = ',pi,' AGAIN PI = ',22.0/7.0
!
! Average number of Ar atoms in the bin
!
 do i = 1,nbin
   rr = (i-1)*binSize + 0.5*binSize  ! 0.5*binSize will take the mid point of the bin
   dv = 4.0*pi*rr*rr*binSize         ! Volume of the shell at distance rr
   GofR = bin(i)/dv/rho
   write(20,'(f8.4,f12.5)') rr,GofR
 end do
!
!
 end program radialdis
