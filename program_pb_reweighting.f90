PROGRAM reweighting
implicit none
real*8, allocatable :: grid_max(:), grid_min(:), grid_width(:)
real*8, allocatable :: ct(:), dum(:), ds2(:), dt(:), diff_s2(:)
real*8, allocatable :: prob(:,:,:,:), cv(:,:), vbias(:)

real*8 :: bias_fact, norm, t_cv, t_sys, alpha_cv, beta_cv, addition, dum2, dummy    
real*8 :: gridmin1, gridmin2, griddif1, griddif2, rweight, s1, s2, s3, s4, gamma_       

integer, allocatable :: nbin(:)

integer :: i_mtd, mtd_steps, ns, ios, is, ig, i, j, i_md, md_steps, mtd_max, dum3
integer :: w_cv, w_hill, k, t_min, t_max, index1, index2, index3, index4, i_s1, i_s2,i_s3, i_s4

REAL*8,PARAMETER :: pi=4.d0*atan(1.d0)
REAL*8, PARAMETER :: kb=8.314472e-3                                                                  ! kJ/(mol*K)
REAL*8, PARAMETER :: kj_to_kcal=1.d0/4.214d0

open(1,file='input')
open(2,file='HILLS_phi1')
open(4,file='COLVAR')
open(7,file='vbias.dat')
open(8,file='ct.dat')

open(50,file='PROB_4D.dat',form='unformatted',status='replace' )

read(1,*) ns
read(1,*) t_cv, t_sys, bias_fact                                   ! CVs_temp, system_temperature, bias_factor
read(1,*) w_cv, w_hill                                ! frequency of printing cvmdck file, frequency of hill added 
read(1,*) t_min, t_max

print *, 'cv freq=', w_cv
print *, 'hill freq=', w_hill

CALL step_count(2,mtd_steps)
print *, 'mtd steps from count =', mtd_steps
rewind(2)

CALL step_count(4,md_steps)
print *, 'md steps from count =', md_steps
rewind(4)

IF(t_max.gt.md_steps)STOP '!!ERROR: t_max > total MD steps'

allocate(grid_max(ns))
allocate(grid_min(ns))
allocate(grid_width(ns))
allocate(nbin(ns))
allocate(ct(mtd_steps))
allocate(dt(ns))


read(1,*) grid_max(1:ns)
read(1,*) grid_min(1:ns)
read(1,*) nbin(1:ns)

print *, 'bins used=', nbin(1:ns)

grid_width(1:ns) = (grid_max(1:ns) - grid_min(1:ns))/float(nbin(1:ns)-1)

print *, 'grid widths =', grid_width(1:ns)

allocate(cv(ns,md_steps))

do k = 1, md_steps
  read(4,*)   dum2, cv(3,k), cv(4,k) ,cv(1,k), cv(2,k),dum, dum, dum, dum, dum
!                             phi2, psi2, phi1, psi1
              if (cv(1,k) .gt. pi)  cv(1,k) = cv(1,k) - 2.d0*pi
              if (cv(1,k) .lt. -pi) cv(1,k) = cv(1,k) + 2.d0*pi
              if (cv(2,k) .gt. pi)  cv(2,k) = cv(2,k) - 2.d0*pi
              if (cv(2,k) .lt. -pi) cv(2,k) = cv(2,k) + 2.d0*pi
              if (cv(3,k) .gt. pi)  cv(3,k) = cv(3,k) - 2.d0*pi
              if (cv(3,k) .lt. -pi) cv(3,k) = cv(3,k) + 2.d0*pi
              if (cv(4,k) .gt. pi)  cv(4,k) = cv(4,k) - 2.d0*pi
              if (cv(4,k) .lt. -pi) cv(4,k) = cv(4,k) + 2.d0*pi
end do



gamma_ = (bias_fact - 1.0)/bias_fact
beta_cv =  1.d0/(kb*t_cv)

print *, 'gamma_ = ', gamma_
print *, 'beta_cv = 1/(kb*t_cv) =  ', beta_cv

print *, ' reading ct and vbias factor.'

allocate(vbias(md_steps))

do i=1,mtd_steps
read(8,'(I10,F16.6)') dum3, ct(i)
end do

do k=1,md_steps
read(7,'(I10,2F16.6)') dum3, vbias(k), dum2
end do



allocate(prob(nbin(1),nbin(2),nbin(3),nbin(4)))
norm=0.0
prob=0.0
print *, 'calculating unbiased probabilty density => free energy.dat'

DO i_md=1,md_steps
        IF((i_md.GE.t_min).AND.(i_md.LT.t_max))THEN                                      ! t_min = 1 not 0
          index1 = nint((cv(1,i_md)-grid_min(1))/grid_width(1))+1            ! US is applied
          index2 = nint((cv(2,i_md)-grid_min(2))/grid_width(2))+1            ! pbmtd applied
          index3 = nint((cv(3,i_md)-grid_min(3))/grid_width(3))+1            ! pbmtd applied
          index4 = nint((cv(4,i_md)-grid_min(4))/grid_width(4))+1            ! pbmtd applied
           if(index1.gt.0.and.index2.gt.0.and.index3.gt.0.and.index1.le.nbin(1).and.index2.le.&
              nbin(2).and.index3.le.nbin(3).and.index4.gt.0.and.index4.le.nbin(4))then
!            i_mtd =  (i_md-1)*w_cv/w_hill 

             i_mtd= ((i_md-1)*w_cv/w_hill) + 1

             dummy = vbias(i_md) + ct(i_mtd)      
                    
!             print *, 'i_md, i_mtd, vbiad, ct, dummy =', i_md, i_mtd, vbias(i_md), ct(i_md), dummy 
             rweight = dexp(beta_cv*dummy)
             prob(index1,index2,index3,index4) = prob(index1,index2,index3,index4) + rweight
             norm=norm+rweight
           end if
        END IF
      END DO



 print *,' norm =', norm
 do is=1,ns
 norm=norm*grid_width(is)
 end do                         
 print *, 'norm per unit area=', norm

    DO i_s1=1,nbin(1)
    s1= grid_min(1) + float(i_s1-1)*grid_width(1)
       DO i_s2=1,nbin(2)
       s2= grid_min(2) + float(i_s2-1)*grid_width(2)
          DO i_s3=1,nbin(3)
           s3= grid_min(3) + float(i_s3-1)*grid_width(3)
              DO i_s4=1,nbin(4)
              s4= grid_min(4) + float(i_s4-1)*grid_width(4)

               prob(i_s1,i_s2,i_s3,i_s4) = prob(i_s1,i_s2,i_s3,i_s4)*(1.d0/norm)                                    
              
               WRITE(50) prob(i_s1,i_s2,i_s3,i_s4)

              END DO
          END DO
       END DO
    END DO

print *, 'PROB_4D.dat =>  prob_density'
close (50)

print *, 'done.'
close (1)
close (2)
close (3)


deallocate(cv)
deallocate(nbin)
deallocate(grid_max)
deallocate(grid_min)
deallocate(grid_width)
deallocate(ct)
deallocate(vbias)
deallocate(prob)

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
