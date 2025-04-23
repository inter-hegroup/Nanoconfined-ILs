!------
!  program for caculating electrical potential distribution
!------
program elepotential

implicit none 

real :: lx,ly,lz,bin,sumc,summ,e2c,a2m,ep,sig,low,q,scale1
integer :: nsteps,nsite,i,j,k,l,nskip
integer :: ca,an,co2,nca,nan,nco2,ngo,ngra,nar,z,nff
real,allocatable ::mass1(:),mass2(:),mass3(:),cha1(:),cha2(:),cha3(:),chard(:,:),massd(:,:),zc(:),zm(:)
real,allocatable :: pos1(:,:,:,:),pos2(:,:,:,:),pos3(:,:,:,:),pos4(:,:,:)
character(len = 20):: ti,cc,filename 


open(11,file = "profile",status = "old")
read(11,*) nsteps,nskip,q
read(11,*)lx,ly,lz,bin,low
read(11,*)ca,an,co2,nca,nan,nco2,ngo,ngra,nar
read(11,*)
 
allocate (mass1(nca),cha1(nca),mass2(nan),cha2(nan),mass3(nco2),cha3(nco2))
 
do i = 1,nca
  read(11,*)ti,mass1(i),cha1(i)
end do 
 
read(11,*)
 
do i = 1,nan
  read(11,*)ti,mass2(i),cha2(i)
end do
 
read(11,*)

do i = 1,nco2 
  read(11,*)ti,mass3(i),cha3(i)
end do 
 
read(11,*)filename 
open(13,file = filename, status = "old")

!------
!  Parameter 
!------
nsite = int(lz/bin)
e2c = 1.6e-19      
a2m = 1e-10         
scale1 = 1.0/0.6022 

allocate(pos1(nsteps,ca,nca,3),pos2(nsteps,an,nan,3),pos3(nsteps,co2,nco2,3),pos4(nsteps,nff,3))
allocate(chard(nsteps,nsite),massd(nsteps,nsite),zc(nsite),zm(nsite))



nff = ngo+nar+ngra

do i = 1,nskip
   read(13,*)
   read(13,*)
   do j = 1,ca
     do k = 1,nca
       read(13,*)
     end do
   end do

   do j = 1,an
     do k = 1,nan
       read(13,*)
     end do
   end do

   do j = 1,co2
     do k = 1,nco2
       read(13,*)
     end do
   end do

   do j = 1,nff
      read(13,*)
   end do

end do

do i = 1,nsteps  
   read(13,*)
   read(13,*)
   do j = 1,ca  
     do k = 1,nca  
       read(13,*)cc,pos1(i,j,k,:)
     end do
   end do
   
   do j = 1,an
     do k = 1,nan
       read(13,*)cc,pos2(i,j,k,:)
     end do 
   end do 

   do j = 1,co2 
     do k = 1,nco2 
       read(13,*)cc,pos3(i,j,k,:)
     end do 
   end do 
   
   do j = 1,nff
       read(13,*)
   end do
end do


do i = 1,nsteps
  do j = 1,ca
    do k = 1,nca
       l = floor((abs(pos1(i,j,k,3) - low))/bin)+1
       chard(i,l) = chard(i,l) + cha1(k)
       massd(i,l) = massd(i,l) + mass1(k)
    end do
  end do
  
  do j =1,an
    do k = 1,nan
      l = floor((abs(pos2(i,j,k,3) - low))/bin)+1
      chard(i,l) = chard(i,l) + cha2(k)
      massd(i,l) = massd(i,l) + mass2(k)
    end do
  end do 

  do j = 1,co2
    do k = 1,nco2
      l = floor((abs(pos3(i,j,k,3) - low))/bin)+1
      chard(i,l) = chard(i,l) + cha3(k)
      massd(i,l) = massd(i,l) + mass3(k)
    end do 
  end do 
end do


do l = 1,nsite
  sumc = 0.0
  summ = 0.0
  do i = 1,nsteps
    sumc = sumc + chard(i,l)
    summ = summ + massd(i,l)
  end do
    zc(l) = (sumc/dble(nsteps)/bin/lx/ly)   
    zm(l) = (summ/dble(nsteps)/bin/lx/ly)*scale1   
end do

open(16,file = "dc.dat",status = "unknown")


write(16,*)  'distance density charge'

do i = 1,nsite
  write(16,"(3F12.5)") dble(i)*bin+low,zm(i),zc(i)*1000 
end do 

end 
