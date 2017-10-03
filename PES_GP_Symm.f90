module PES_details
  double precision :: gpRMax = 9.0   *   1.8897259885789
  double precision :: gpRMin = 2.5   *   1.8897259885789
  double precision :: gpEmax = 0.027031002
  !! Maximum energy in Hartree seen within the geometric constraint
  interface PES_GP 
     function PES_GP(xStar) 
       implicit none 
       double precision:: PES_GP
       double precision, dimension(:) ::  xStar
     end function PES_GP
  end interface PES_GP
end module PES_details


module GP_variables
  double precision, allocatable :: alpha (:), lScale(:), xTraining(:,:), xTrainingPerm(:,:)
  double precision expVar,NuggVar
  integer :: nDim=1
  integer :: nTraining=72 
end module GP_variables


 ! Test program
implicit none


integer k,i, choice


call load_GP_Data
call fixedAngleSlice

end
!

subroutine fixedAngleSlice
  use PES_details
  !use GP_variables
  implicit none
  double precision rab(1)
  integer i, itot
  double precision  r, tha, thb, ph, e, e_GP, asymp, PES,AngToBohr
  AngToBohr= 1.8897259885789
    
  itot=700

  open (unit=15, file="PES_Out.dat ", status='replace')
  
  do i=0, itot

     ! specify centre-to-centre separation in Bohr
     r = (  2.4 + 15.0*i/(1.0*itot) ) * AngToBohr
     rab(1)=r
     
     
     e=PES( rab)
     write(15,*) r/AngToBohr , e 
     
  enddo
  write(6,*)'Written to file: PES_Out.dat '
  close(15)

end subroutine fixedAngleSlice
  
  
  
double precision function asymp(rab)
! Work out the asymptotic energy for Ar-Ar
! Everything is in Atomic Units (Bohr, Hartree)
implicit none
double precision rab(1)

asymp= - 64.3 / rab(1)**6 - 1642 / rab(1)**8
!- 64.3 / r^6 - 1642 / r^8
end
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Gaussian Process Code!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!pom


  
subroutine load_GP_Data
  use GP_variables
  use PES_details
  implicit none
  
  double precision, allocatable::  xStar(:)
  integer i,j
  double precision :: dum
  character (len=90) :: filename

  allocate (alpha(nTraining), lScale(nDim), xTraining(nDim,nTraining),xTrainingPerm(nDim,nTraining), xStar(nDim))

  !====Load hyperparameters====
  write (filename,  '( "TrainingData/HyperParams_Symm", I2.2, ".dat" )' )  nTraining
  open (unit = 7, file = filename)
  !Only need to read some as others are tied.
  read (7,*) lScale(1),expVar,NuggVar,dum
  !Copy tied hyperparameter values
  !Not needed for 2Ar as no symmetries
  !print *,"HyperParams",lScale(1),expVar,NuggVar, dum
  close(7)
  
  !====Load alpha coefficients====
  write (filename, '( "TrainingData/alpha_Symm", I2.2, ".dat" )' )  nTraining
  open (unit = 7, file = filename)
  do i=1,nTraining
     read (7,*) alpha(i)
     !!print *,"alpha ",i, alpha(i)
  end do
  close(7)


  !====Load training data x values ====
  write (filename, '( "TrainingData/xTraining", I2.2, ".dat" )'  )  nTraining
  open (unit = 7, file = filename)
    
  do i=1,nTraining
     read (7,*) xTraining(1,i)
     !print *,"xTraining ",i,xTraining(1,i)
  end do
  close(7)
  
  ! Permute the training vectors (not needed for 2Ar as only one vector)
  !xTrainingPerm = xTraining
  !do i=1,nTraining
  !   xTrainingPerm(3,i)=xTraining(5,i)
   !  xTrainingPerm(5,i)=xTraining(3,i)
    ! xTrainingPerm(4,i)=xTraining(6,i)
    ! xTrainingPerm(6,i)=xTraining(4,i)
  !end do

end subroutine load_GP_Data
  
function PES_GP(xStar)
  use GP_variables
  implicit none
  double precision, dimension(:) :: xStar
  double precision:: PES_GP
  integer beta,i
  double precision kSqExp, kSqExpPerm, kKern

  kKern=0

  do i=1,nTraining
     kSqExp=1;
     do beta=1,nDim
        kSqExp  =  kSqExp * ( exp( - (xStar(beta)-xTraining(beta,i))**2 /2.0/lScale(beta)**2) )
     end do
     kKern = kKern + alpha(i) * kSqExp 
  end do
  
  PES_GP=kKern * expVar
end function PES_GP



function PES( rab)
  !! Takes in rab in Angstrom
  use PES_details
  implicit none
  double precision rab(1), xStar(1),  asymp
  double precision AngToBohr, PES

  AngToBohr= 1.8897259885789

  
  if( rab(1) > gpRMax) then !!Use asymptotic function
     PES = asymp(rab)
     
  else if (rab(1) < gpRMin  ) then !! Use repulsive approximation function
     PES=gpEmax* ( (gpRMin)/rab(1))**12 
     
  else !! Use the Guassian Process function
     xStar(:) = AngToBohr/rab(:)
     PES = PES_GP( xStar)
  end if

end function PES
