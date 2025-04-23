program main 

!--------------
! the 2nd version
! Sun Mar  6 15:59:47 CST 2022 in UNIST
! please check the algorithm very carefully
!---------------

!--------------
! the Chenlu Wang version
! Mon Mar 7 21:22:47 CST 2022 in IPE
!---------------

implicit none
real :: rtip,rcut
real :: dx,dy,dz,dr,box(3) 
real,allocatable :: rtype(:),zlo(:),zhi(:) 
real,allocatable :: pos(:,:,:),grid(:,:)
real :: mesh,v0,xlo,ylo,xhi,yhi,zmin,zmax 
integer :: lx,ly,lz,ngrid,neigh_max
integer :: null1 
character :: atomname 
integer :: i,j,k,l 
integer :: nsteps,nskip,natoms 
integer :: ntype
integer :: nvol,invol
integer,allocatable :: neigh(:,:),at(:),occupy(:,:),inv(:,:),check(:)
integer :: id,mol
real :: x,y,z


open(10,file = "profile",status = "old")
read(10,*) nsteps,nskip,mesh
read(10,*) rtip 
read(10,*) ntype
allocate (rtype(ntype))
do i = 1,ntype
 read(10,*) null1,rtype(i)
 rtype(i) = rtype(i)*rtype(i) 
end do 
v0 = mesh**3 
rcut = rtip*rtip 


open(11,file = "dump.lammpstrj",status = "old")

do i = 1,nskip
 read(11,*)
 read(11,*)
 read(11,*)
 read(11,*) natoms
 read(11,*)
 read(11,*) xlo,xhi
 read(11,*) ylo,yhi
 read(11,*) 
 read(11,*) 
 do j = 1,natoms
  read(11,*)
 end do
end do 


allocate(pos(nsteps,natoms,3),at(natoms),zhi(nsteps),zlo(nsteps))
zmin = 100.0
zmax = -100.0
do i = 1,nsteps
 read(11,*)
 read(11,*)
 read(11,*)
 read(11,*) natoms
 read(11,*)
 read(11,*) xlo,xhi
 read(11,*) ylo,yhi
 read(11,*) 
 read(11,*)
 box(1) = xhi-xlo  
 box(2) = yhi-ylo

 do j = 1,natoms
  read(11,*) id,mol,at(j),atomname,pos(i,j,:)
   if (at(j) == 17) then 
     zhi(i) = pos(i,j,3)
   end if 
   
   if (at(j) == 20) then
     zlo(i) = pos(i,j,3)
   end if



   if (at(j) == 17 .and. pos(i,j,3) >= zmax) then 
     zmax =  pos(i,j,3)
   end if 

   if (at(j) == 20 .and. pos(i,j,3) <= zmin) then
     zmin =  pos(i,j,3)
   end if   
 end do
end do
box(3) = zmax-zmin


lx = int(box(1)/mesh)+1
ly = int(box(2)/mesh)+1
lz = int(box(3)/mesh)+1
allocate(grid(lx*ly*lz,3))
ngrid = lx*ly*lz


neigh_max = 200

l = 0 
do i = 1,lx
 do j = 1,ly
  do k = 1,lz
   l = l+1
   grid(l,1) = xlo+(i-1)*mesh
   grid(l,2) = ylo+(j-1)*mesh
   grid(l,3) = zmin+(k-1)*mesh 
  end do
 end do
end do 
write(*,*) "hello"

allocate (neigh(ngrid,neigh_max))
neigh = 0  
write(*,*) "good"
do i = 1,ngrid
 l = 0 
 do j = 1,ngrid
  if (i /= j) then
   dx = abs(grid(i,1)-grid(j,1))
   dy = abs(grid(i,2)-grid(j,2))
   dz = abs(grid(i,3)-grid(j,3))
   dx = abs(dx-ANINT(dx/box(1))*box(1)) 
   dy = abs(dy-ANINT(dy/box(2))*box(2))
   dr = dx*dx+dy*dy+dz*dz 
   if (dr <= rcut) then  
    l=l+1
    neigh(i,l) = j 
   end if 
  end if 
 end do 
end do 
write(*,*) "evening"

allocate(occupy(nsteps,ngrid),inv(nsteps,ngrid),check(ngrid))
occupy = 0 
inv = 0 

do i = 1,nsteps
 do j = 1,natoms
  if (pos(i,j,3) <= zhi(i)) then 
   do k = 1,ngrid 
     dx = abs(pos(i,j,1)-grid(k,1))
     dy = abs(pos(i,j,2)-grid(k,2))
     dz = abs(pos(i,j,3)-grid(k,3))
     dx = abs(dx-ANINT(dx/box(1))*box(1))
     dy = abs(dy-ANINT(dy/box(2))*box(2))
     dr = dx*dx+dy*dy+dz*dz
     if (dr <= rtype(at(j))) then 
      occupy(i,k) = 1 
      inv(i,k) = 1 
     end if 
   end do 
  end if 
 end do 
end do   


open(12,file="volume.data",status = "unknown")
open(13,file="cluster.lammpstrj",status = "unknown")

do i = 1,nsteps
check = 0 
nvol=0
invol=0
 do j = 1,ngrid 
  if (occupy(i,j) == 0 )then  
   do k = 1,neigh_max 
    if (neigh(j,k) > 0 ) then  
      if (occupy(i,neigh(j,k)) == 1) then 
        inv(i,j) = 1 
      end if 
    end if 
   end do
  end if 
 end do

 do j = 1,ngrid  
  if (grid(j,3) <= zhi(i)) then 
    nvol = nvol+1 
   if (inv(i,j) == 1) then 
    invol = invol+1 
    check(j) = 1
   end if 
  end if 
 end do 

 write(12,"(I5,3F12.5)") i,nvol*v0,invol*v0,(nvol-invol)*v0

 write(13,"(A)") "ITEM: TIMESTEP"
 write(13,"(A)") "0"
 write(13,"(A)") "ITEM: NUMBER OF ATOMS"
 write(13,"(I8)") nvol-invol
 write(13,"(A)") "ITEM: BOX BOUNDS pp pp pp"
 write(13,"(A)") "-2.6402515399999999e+01 2.7632267680000002e+01"
 write(13,"(A)") "-3.7576991999999997e+01 3.8997007330000002e+01"
 write(13,"(A)") "-1.5000000000000000e+02 1.5000000000000000e+02"
 write(13,"(A)") "ITEM: ATOMS id x y z"
 l = 0 

 do j = 1,ngrid
  if (check(j) == 0 .and. grid(j,3)<= zhi(i)) then
    l = l+1
    write(13,"(I5,3F12.5)") l,grid(j,:)
  end if 
 end do 

end do 

end 
