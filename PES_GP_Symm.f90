! Test program
implicit none


integer k,i, choice


call load_GP_2Ar_Data
call fixedAngleSlice

end



subroutine fixedAngleSlice
  use PES_2Ar_details
  !use GP_2Ar_variables
  implicit none
  double precision rab(1)
  integer i, itot
  double precision  r, tha, thb, ph, e, e_GP, asymp, PES_2Ar,AngToBohr
  AngToBohr= 1.8897259885789
    
  itot=700

  open (unit=15, file="PES_2Ar_Out.dat ", status='replace')
  
  do i=0, itot

     ! specify centre-to-centre separation in Bohr
     r = (  2.4 + 15.0*i/(1.0*itot) ) * AngToBohr
     rab(1)=r
     
     
     e=PES_2Ar( rab)
     write(15,*) r/AngToBohr , e 
     
  enddo
  write(6,*)'Written to file: PES_2Ar_Out.dat '
  close(15)

end subroutine fixedAngleSlice
  



!
