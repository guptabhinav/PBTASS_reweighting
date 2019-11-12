! PBMTD along phi1 and phi2 and TAMD along phi1,phi2,psi1,psi2
PROGRAM reweighting
implicit none
REAL*8, ALLOCATABLE :: hh(:), ss(:), grid_max(:), grid_min(:)
REAL*8, ALLOCATABLE :: grid_width(:), dum(:), ds2(:), dt(:), diff_s2(:)
REAL*8, ALLOCATABLE :: v(:,:), hill(:,:), ht(:,:), width(:,:)
REAL*8, ALLOCATABLE ::  grid(:,:), ct(:), num(:), den(:), cv(:,:)

REAL*8 :: bias_fact, norm, t_cv
REAL*8 :: beta_cv, addition, t_sys                 
REAL*8 :: gridmin1, gridmin2, griddif1, griddif2, rweight, s1, s2
REAL*8 :: dum1, gamma_

INTEGER, ALLOCATABLE :: nbin(:)

INTEGER :: i_mtd, mtd_steps, ns, ios, is, ig, i, j, i_md
INTEGER :: md_steps, mtd_max, dum2,  i_s1, i_s2, np
INTEGER :: w_cv, w_hill, k, t_min, t_max, index1, index2

REAL*8,PARAMETER :: pi=4.d0*atan(1.d0)
REAL*8, PARAMETER :: kb=8.314472e-3                                                                  ! kJ/(mol*K)
REAL*8, PARAMETER :: kj_to_kcal=1.d0/4.214d0

open(1,file='input')

read(1,*) ns
read(1,*) t_cv, t_sys, bias_fact                                   ! CVs_temp, system_temperature, bias_factor
read(1,*) w_cv, w_hill                                ! frequency of printing cvmdck file, frequency of hill added 
read(1,*) t_min, t_max

open(2,file='HILLS_phi1')
open(3,file='HILLS_phi3')
open(7,file='HILLS_phi4')
open(4,file='COLVAR')

open(20,file='ct.dat',status='unknown')


CALL step_count(2,mtd_steps)
print *, 'mtd steps from count =', mtd_steps
rewind(2)

CALL step_count(4,md_steps)
print *, 'md steps from count =', md_steps
rewind(4)


allocate(width(ns,mtd_steps))
allocate(ht(ns,mtd_steps))
allocate(hill(ns,mtd_steps))
allocate(hh(ns))
allocate(ss(ns))
allocate(grid_max(ns))
allocate(grid_min(ns))
allocate(grid_width(ns))
allocate(nbin(ns))
allocate(ct(mtd_steps))
allocate(ds2(ns))
allocate(diff_s2(ns))
allocate(dt(ns))
allocate(num(ns))
allocate(den(ns))


read(1,*) grid_max(1:ns)
read(1,*) grid_min(1:ns)
read(1,*) nbin(1:ns)

print *, 'bins used=', nbin(1:ns)

grid_width(1:ns) = (grid_max(1:ns) - grid_min(1:ns))/float(nbin(1:ns)-1)

print *, 'grid widths =', grid_width(1:ns)


do i=1,mtd_steps
 read(2,*)   dum1, hill(1,i), width(1,i), ht(1,i), dum2
          IF( hill(1,i) .gt. pi) hill(1,i) = hill(1,i) - 2.d0*pi
         IF( hill(1,i) .lt. -pi) hill(1,i) = hill(1,i) + 2.d0*pi
end do

do i=1,mtd_steps
 read(3,*)   dum1, hill(2,i), width(2,i), ht(2,i), dum2 
          IF( hill(2,i) .gt. pi) hill(2,i) = hill(2,i) - 2.d0*pi
         IF( hill(2,i) .lt. -pi) hill(2,i) = hill(2,i) + 2.d0*pi
end do

do i=1,mtd_steps
 read(7,*)   dum1, hill(3,i), width(3,i), ht(3,i), dum2
          IF( hill(3,i) .gt. pi) hill(3,i) = hill(3,i) - 2.d0*pi
         IF( hill(3,i) .lt. -pi) hill(3,i) = hill(3,i) + 2.d0*pi
end do

allocate(cv(ns,md_steps))

do k = 1, md_steps
  read(4,*)   dum1,cv(1,k),dum1,dum1,dum1,cv(2,k),dum1,cv(3,k),dum1
              if (cv(1,k) .gt. pi)  cv(1,k) = cv(1,k) - 2.d0*pi
              if (cv(1,k) .lt. -pi) cv(1,k) = cv(1,k) + 2.d0*pi
              if (cv(2,k) .gt. pi)  cv(2,k) = cv(2,k) - 2.d0*pi
              if (cv(2,k) .lt. -pi) cv(2,k) = cv(2,k) + 2.d0*pi
              if (cv(3,k) .gt. pi)  cv(3,k) = cv(3,k) - 2.d0*pi
              if (cv(3,k) .lt. -pi) cv(3,k) = cv(3,k) + 2.d0*pi
end do


beta_cv =  1.d0/(kb*t_cv)
gamma_ = (bias_fact - 1.0)/bias_fact

print *, 'gamma_ = ', gamma_
print *, 'beta_cv = 1/(kb*t_cv) =  ', beta_cv

allocate (grid(ns,nbin(1)))

do i=1,nbin(1)
grid(1,i) = grid_min(1) + float(i-1)*grid_width(1)
end do

do j=1,nbin(2)
grid(2,j) = grid_min(2) + float(j-1)*grid_width(2)
end do

do j=1,nbin(3)
grid(3,j) = grid_min(3) + float(j-1)*grid_width(3)
end do

allocate(v(ns,nbin(1)))

print *, 'calculating ct factor.'
      v=0.d0
      DO i_mtd=1,mtd_steps
         ds2(1:ns)=width(1:ns,i_mtd)*width(1:ns,i_mtd)
         ss(1:ns)=hill(1:ns,i_mtd)
         
         hh(1:ns)=ht(1:ns,i_mtd)                  
         
         num=0.D0
         den=0.D0
         DO is=1,ns

            DO ig=1,nbin(is)
              diff_s2(is)=grid(is,ig)-ss(is)

               if (diff_s2(is) .gt. pi ) diff_s2(is) =diff_s2(is) - 2.d0*pi  
               if (diff_s2(is) .lt.-pi ) diff_s2(is) =diff_s2(is) + 2.d0*pi 

              diff_s2(is)=diff_s2(is)*diff_s2(is)*0.5D0
              v(is,ig) = v(is,ig) + hh(is)*dexp(-diff_s2(is)/ds2(is))               
              num(is) = num(is) + dexp(beta_cv*v(is,ig) - beta_cv*gamma_*v(is,ig))                 
              den(is) = den(is) + dexp(beta_cv*v(is,ig))                                 

            END DO
         END DO

         addition=0.d0
         do is=1,ns
          addition = addition + num(is)/den(is)! + 1.d0
         end do       
       
         ct(i_mtd)= (1.d0/beta_cv)*dlog(addition)
            write(20,'(I10,F16.6)')  i_mtd,  ct(i_mtd)
         END DO
close (20)
print *, 'done calculating ct. => ct.dat'


print *, 'done.'
close (1)
close (2)
close (3)


deallocate(cv)
deallocate(width)
deallocate(ht)
deallocate(hill)
deallocate(hh)
deallocate(ss)
deallocate(nbin)
deallocate(grid_max)
deallocate(grid_min)
deallocate(grid_width)
deallocate(ct)
deallocate(ds2)
deallocate(diff_s2)
deallocate(grid)
deallocate(num)
deallocate(den)

end program


subroutine step_count(file_number,steps)
integer :: file_number, steps, ios,i
steps=0
do
 read(file_number,*,iostat=ios)
 if(ios.ne.0) exit
  steps=steps+1
end do
end subroutine
