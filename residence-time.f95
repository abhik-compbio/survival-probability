! Survival probability calculations for interfacial water molecules
! Water molecules near 6 angstorm from protein Backbone atom (C,N,CA) are considerd only
! Written by Abhik Ghosh Moulick, SNBNCBS
	
integer,parameter:: atomno = 369  ! change,check last atomno of complex from pdb file
integer,parameter:: solno =  13118 ! change, check lst solno from pdb file
integer,parameter:: maxframe =  1005 ! change, maximum frame to consider
real :: xc(maxframe,1:atomno,3),xw(maxframe,atomno+1:solno,3)
real :: num(atomno+1:solno)  
integer :: ichain,atnum,gap,p,no
character :: line*80
character :: atnam*4, junk*1
integer :: resnum,j,i,steps,nbin,resinum,resn,origin,delay
real :: r,Gs,dr,rGs,Pt1,P1t1,P1,P2,Sp(0:300)
real :: dist, dist_x, dist_y, dist_z, dist_gap,M
real :: dist_x1, dist_x2, dist_x3, Multiply,distance(maxframe,atomno+1:solno)
data delt/1.0/
open(10, file='Input-survivalprobability.pdb') ! Change file name
open(20, file='Sp-water-around-protein-free.dat')! Change file name

! Read pdb file
ichain=1
do while (ichain.ge.1)	
read (10,'(a)') line
if (line(1:4).eq.'ATOM') then
	read (line,1) resnum,atnam,junk,resn,x,y,z
	!write (22,1) resnum,atnam,junk,resn,x,y,z
1 	format (6x,i5,2x,a4,3x,a2,1x,i4,3x,3f8.3)
	if(atnam.ne.'OW') then
		xc(ichain,resnum,1) = x
		xc(ichain,resnum,2) = y
		xc(ichain,resnum,3) = z
	else if(atnam.eq.'OW') then
		xw(ichain,resnum,1) = x
		xw(ichain,resnum,2) = y
		xw(ichain,resnum,3) = z
	end if
end if				
	if (line(1:4).eq.'TER ') THEN
		ichain = ichain + 1
		maxchain = ichain
		!write(*,*) ichain
		if (maxchain .eq. 1002) exit
	end if
end do
 close(10)
      	
! Survival probability calculations
      	
steps = maxchain - 1 
origin = 800
delay = 200
do gap = 0,delay
	Sp(gap)=0.0
end do

do ichain = 1,steps
	do j = atomno+1,solno
	Pt1 = 0.0
		do i = 1,atomno
			dist_x = xc(ichain,i,1)-xw(ichain,j,1)
      			dist_y = xc(ichain,i,2)-xw(ichain,j,2)
      			dist_z = xc(ichain,i,3)-xw(ichain,j,3)
      			dist = sqrt(dist_x**2+dist_y**2+dist_z**2)
			if (dist .le. 6.0) then
      				Pt1 = Pt1 + 1.0
      			end if
      		end do
      		if (Pt1 .gt. 0.0) then
			distance(ichain,j) = 1.0
		else 
			distance(ichain,j) = 0.0
		end if
	end do
end do


do j = atomno+1,solno
	do i = 1,origin
		M = distance(i,j)
		do k = 1,delay
			M = M*distance(i+k,j)
			Sp(k) = Sp(k) + M
		end do
	end do
end do

do gap = 1,delay
	Sp(gap)=Sp(gap)/origin
end do


do gap = 1,delay
	write(20,6) float(gap),Sp(gap)/Sp(1)
6     	format(2(f10.4,' 			'))
end do
end
