! Monte Carlo algorithm to simulate the Renyi ensemble probability distribution. Based on the Monte Carlo code in the ising.f90 file in https://github.com/carrasqu/2d_ising_model_square_triangular

!Modules are used for
!  Packaging subprograms, data and interface blocks.
!  Defining global data that can be used by more than one routine.
!  Declaring variables that can be made available within any routines you choose.
!  Importing a module entirely, for use, into another program or subroutine.

!Nx,Ny = self_explanatory
!nh is later defined as Nx * Ny (total nbumber of spins)
!which_lattice = 1 --> square lattice; which_lattice = 2 --> triangular lattice
!nn = number of nearest neighbors. For the square lattice, it will be 4. For the triangular lattice it will be 6.
!ivic = array that stores indices of nearest neighbors for each site
!ivict = array used for triangular lattice case, I did not investigate it for my model.
!ixpiy: It's needed for the calculation of the staggered magnetization in the event the coupling is antiferromagnetic (H = -J * sum... where J is negative --> antiferromagnetic). The idea is ixpiy = i+j (row index + column index) of a given site, so that (-1)**ixpiy alternates between -1 and 1 as you move from one site to the next (in any direction). ixpiy is then used for the staggered magnetization calculation
!spin(:) = array of size nh that will store the spin value associated with every site (nh=Nx*Ny) sites in total.
! ...
!drand1() produces a random floating point number between 0 and 1. How it's done doesn't really matter. See random.c file.
!temperature and beta should be self_explanatory. Jh I think is the nearest neightbor coupling constant (nearest neighbor interaction energy).

! First model: configuration
module configuration
save 
integer(4)Nx,Ny,nh,which_lattice,nn
integer(4), allocatable :: ivic(:,:),ivict(:,:),ixpiy(:) ! nearest neighbors
integer(4), allocatable :: spin(:) ! spin (-1 or 1) or (boson 0 or 1 respectively) configuration
real(8) temperature,beta,Jh,h,alpha !alpha = Renyi index... adding Tmin, Tmax and Tstep to be able to loop through temperatures.
real(4) Tmin,Tmax,Tstep
end module configuration
!The allocatable attribute provides a safe way for memory handling. In comparison to variables with pointer attribute the memory is managed automatically and will be deallocated
!automatically once the variable goes out-of-scope. Using allocatable variables removes the possibility to create memory leaks in an application.

!Next module: measurement
!lm is set to 2, and while I don't know what it stands for, it ends up defining the size of the mea array which contains the measurement results (mag and energy) I think.
!All the "Vm" variables are quite confusing right now. Unclear what they are. They contain measurement results but of what, I don't know.
!E's gotta be the energy. It is, but is it energy per spin, energy of a configuration at each Monte Carlo step, we shall see.
module measurement
save
integer(4)lm !number of quantities we are measuring during the Monte Carlo process
real(8), allocatable :: Vm(:),Vm2(:),Vmtemp(:),Vm2temp(:),Vm2temp_error(:),mea(:) !Vm2temp contains average of squares of mag, E, etc, Vm2temp_error contains corresponding errors
real(8)E !Energy of a configuration
real(8)E_sq !E^2.
real(8), allocatable :: Vm_sq(:), Vmtemp_sq(:), mea_sq(:) !arrays to store measurement values... this will become clear below
end module measurement

! Module with three subroutines, containing the key Monte Carlo subroutines
MODULE runbin_accmeasure

  implicit none

CONTAINS

!Perform Monte Carlo simulation for msteps Monte Carlo steps (msteps * nh total steps) for the Gibbs state.
subroutine runbin_gibbs(msteps)
use configuration !this module defined Nx,Ny,nh,which_lattice,nn,ivic(:,:),ivict(:,:),ixpiy(:),spin(:),temperature,beta,Jh                                                                          
use measurement !this module defined lm,Vm(:),Vm2(:),Vmtemp(:),Vm2temp(:),Vm2temp_error(:),mea(:),E   
implicit none
integer(4)i,j,k,kk,msteps
real(8) DE,Eest,r,drand1,mag,mag_abs,magt
real(8) mag_sq, mag_abs_sq,Esqest
!DE = delta_E, change in energy from present configuration to proposed configuration
mea=0.0 
Eest=0.0 !estimate of E, where E
mag=0.0 !initialization of magnetization per spin
mag_abs=0.0 !initialization of absolute value of magnetization per spin.

!equivalent initialization for squared quantities
mea_sq = 0.0
Esqest = 0.0
mag_sq = 0.0 !initialization of square of mag
mag_abs_sq = 0.0


!Total number of single-spin flip updates = msteps * nh. This will form one "bin".
!The configuration, spin, has already been initialized using the initial() subroutine. We can use "spin" because of the "use configuration" line above. We can also use "mea" because of the "use measurement" line of code as well.
do i=1,msteps !1 mstep = "1 Monte Carlo Step" which means "1 Monte Carlo step per lattice site" for a total of nh steps
    do j=1,nh !nh = total number of sites
        k=int(drand1()*nh+1) !this is the code that generates a random number between 1 and nh (the number k). Without the +1 you generate a random number between 0 and (nh-1).
        DE=0 ! initializing DE, i.e. delta_E
        do kk=1,nn
            DE=DE+spin(ivic(k,kk)) !sum the nearest neighbor spins to spin(k) - see Barkema & Newman equation (3.10) to understand why this is done
        end do
        DE=2*DE*Jh*spin(k) !equation (3.10) in Barkema & Newman- calculation of delta_E, or DE
        
        !#metropolis
        if (DE<0)then !flip the spin if the change in energy is negative. Actually it should be a flipped spin if it's less than or equal to zero. This is accomplished below (see my comment).
           spin(k)=-spin(k)
           E=E+DE
        else
           r=drand1() !random number between 0 and 1
           if (r<exp(-beta*DE))then !if DE = 0, then r which is [0,1) will always be less than exp(0)=1, so the spin will be flipped.
             spin(k)=-spin(k) !flip spin
             E=E+DE !update the energy
           end if
        end if
     end do
     !so here we have completed one Monte Carlo step (nh updates)

    !only calculating the magnetization after nh updates (1 mstep) and then accumulating all mstep magnetization values in "mag", then dividing by mag at the end
     if(which_lattice==1.and.Jh>0.0)then !if lattice is square and if the coupling is ferromagnetic, then...
        mag=mag+sum(dble(spin))/(Nx*Ny)
        mag_sq = mag_sq + (sum(dble(spin)))**2 !do not "per spin" this guy right now... in case we want to calculate magnetic susceptibility.
        mag_abs=mag_abs+abs(sum(dble(spin)))/(Nx*Ny) !this is absolute value of total magnetization per site. I think dble(spin) makes each element of spin a double precision float, sum(...) sums the elements, and then after that abs(...) takes the absolute value.                                         
        mag_abs_sq = mag_abs_sq + (abs(sum(dble(spin))))**2 !do not "per spin" this guy right now... in case we want to calculate magnetic susceptibility.
     else if(which_lattice==1.and.Jh<0.0)then !if the coupling is antiferromagnetic, then we care about the staggered magnetization I think.
        magt=0.0d0 !initialize object used to calculate staggered magnetization
        do k=1,nh
           magt=magt+dble(spin(k))*(-1.0d0)**ixpiy(k) !calculation of staggered magnetization- multiply every other spin by -1, which amounts to multiplying each spin by (-1)**(i+j).
           !write(6,*)'checking',spin(k),ixpiy(k),(-1.0d0)**ixpiy(k),dble(spin(k))*(-1.0d0)**ixpiy(k)
        end do
        mag=mag+magt/(Nx*Ny)
        mag_sq = mag_sq + magt**2 !also we are not "per spin-ing" this guy here
        mag_abs=mag_abs+abs(magt)/(Nx*Ny) !ok so take absolute value of "staggered magnetization", divide by number of spins
        mag_abs_sq = mag_abs_sq + (abs(magt))**2 !also we are not "per spin-ing" this guy here
     else if(which_lattice==2)then
        mag=mag+sum(dble(spin))/(Nx*Ny)
        mag_sq = mag_sq + (sum(dble(spin)))**2 !we are not "per spin-ing" this guy here
        mag_abs=mag_abs+abs(sum(dble(spin)))/(Nx*Ny)
        mag_abs_sq = mag_abs_sq + (abs(sum(dble(spin))))**2 !we are not "per spin-ing" this guy here
     end if
    !We are only including the magnetization (absolute value) at the end of a full "Monte Carlo step" (nh steps) in the calculation of the average magnetization for that bin. Same with the energy. You sum all the abs(mag) and E values at the end of each mstep, and then divide by msteps, and store the values in mea(1), mea(2) and mea(3). Those are your values for a given bin.
    Eest=Eest+E/(Nx*Ny)
 end do
!Store the bin averages in the mea array
 mea(1)=mag/msteps     !regular mag per spin
 mea(2)=mag_abs/msteps !abs mag per spin
 mea(3)=Eest/msteps    !energy per spin

 !bin average for squared values. Needed for error calculations.
 mea_sq(1) = mag_sq/msteps
 mea_sq(2) = mag_abs_sq/msteps
 mea_sq(3) = Esqest/msteps

 
 
!write(6,*)"mea",mea

end subroutine runbin_gibbs



!SUBROUTINE 2 INSIDE "CONTAINS"- runbin for the Renyi ensemble
subroutine runbin_renyi(msteps,E_bar)
use configuration !this module defined Nx,Ny,nh,which_lattice,nn,ivic(:,:),ivict(:,:),ixpiy(:),spin(:),temperature,beta,Jh                                                                          
use measurement !this module defined lm,Vm(:),Vm2(:),Vmtemp(:),Vm2temp(:),Vm2temp_error(:),mea(:),E    - I think E is the energy of the current configuration. Fair enough.
implicit none
integer(4)i,j,k,kk,msteps
real(8) DE,Eest,r,drand1,mag,mag_abs,magt,E_mu,E_nu,numer,denom,E_bar !E_mu is the energy of the current configuration, E_nu is the energy of the proposed configuration,numer and denom are the numerator and denominator of the acceptance ratio
real(8) mag_sq, mag_abs_sq,Esqest
!DE = delta_E, change in energy from present configuration to proposed configuration
mea=0.0 !array that will store a bin average of mag, abs(mag) and E.
Eest=0.0
mag=0.0
mag_abs=0.0

!equivalent initialization for squared quantities
mea_sq = 0.0 !testing this for Cv and Chi calculations. squared values...
Esqest = 0.0
mag_sq = 0.0 !initialization of square of mag. Defining this to test Cv and Chi calculations, to try something different.
mag_abs_sq = 0.0

!Total number of single-spin flip updates = msteps * nh. This will form one "bin"..
!The configuration, spin, has already been initialized using the initial() subroutine. We can use "spin" because of the "use configuration" line above. We can also use "mea" because of the "use measurement" line of code as well.
do i=1,msteps !OK 1 mstep = "1 Monte Carlo Step" which means "1 Monte Carlo step per lattice site" for a total of nh steps
    do j=1,nh !nh = total number of sites
        k=int(drand1()*nh+1) !this is the code that generates a random number between 1 and nh (the number k). Makes perfect sense. Without the +1 you generate a random number between 0 and (nh-1).
        DE=0 ! initializing DE for every single update
        do kk=1,nn
            DE=DE+spin(ivic(k,kk)) !sum the nearest neighbor spins to spin(k) - see Barkema & Newman equation (3.10) to understand why this is done
        end do
        DE=2*DE*Jh*spin(k) !equation (3.10) in Barkema & Newman- makes perfect sense- calculation of delta_E, or DE
        E_mu = E !energy of current configuration
        !#metropolis
        if (DE<0)then !flip the spin if the change in energy is negative. Actually it should be a flipped spin if it's less than or equal to zero. This is accomplished below (see my comment).
           spin(k)=-spin(k)
           E=E+DE
        else !now we need the acceptance ratio
           E_nu = E+DE !energy of proposed configuration
           !we need to check if E_nu is allowed as per the Renyi constraint. Let's implement the code I wrote in Julia here
           if (E_nu <= (alpha/(beta*(alpha-1))+E_bar))then
              !Now we can proceed if we are here
              numer = (abs(1 - beta * ((alpha-1)/alpha) * (E_nu-E_bar)))**(1/(alpha-1))
              denom = (abs(1 - beta * ((alpha-1)/alpha) * (E_mu-E_bar)))**(1/(alpha-1))
              r=drand1() !random number between 0 and 1
              if (r<(numer/denom))then !if DE = 0, then r which is [0,1) will always be less than exp(0)=1, so the spin will be flipped.
                 spin(k)=-spin(k) !flip spin
                 E=E+DE !update the energy
              end if
           end if
           !Recall, we're only updating E if a spin is flipped.
        end if
     end do
     !so here we have completed one Monte Carlo step (nh updates)

     !only calculating the magnetization after nh updates (1 mstep) and then accumulating all mstep magnetization values in "mag", then dividing by mag at the end
     if(which_lattice==1.and.Jh>0.0)then !if lattice is square and if the coupling is ferromagnetic, then...
        mag=mag+sum(dble(spin))/(Nx*Ny)
        mag_sq = mag_sq + (sum(dble(spin)))**2 !do not "per spin" this guy right now
        mag_abs=mag_abs+abs(sum(dble(spin)))/(Nx*Ny) !this is absolute value of total magnetization per site. I think dble(spin) makes each element of spin a double precision float, sum(...) sums the elements, and then after that abs(...) takes the absolute value                                         
        mag_abs_sq = mag_abs_sq + (abs(sum(dble(spin))))**2 !do not "per spin" this guy right now
     else if(which_lattice==1.and.Jh<0.0)then !if the coupling is antiferromagnetic, then we care about the staggered magnetization.
        magt=0.0d0 !initialize object used to calculate staggered magnetization
        do k=1,nh
           magt=magt+dble(spin(k))*(-1.0d0)**ixpiy(k) !calculation of staggered magnetization- multiply every other spin by -1, which amounts to multiplying each spin by (-1)**(i+j).
           !write(6,*)'checking',spin(k),ixpiy(k),(-1.0d0)**ixpiy(k),dble(spin(k))*(-1.0d0)**ixpiy(k)
        end do
        mag=mag+magt/(Nx*Ny)
        mag_sq = mag_sq + magt**2 !also we are not "per spin-ing" this guy here
        mag_abs=mag_abs+abs(magt)/(Nx*Ny) !ok so take absolute value of "staggered magnetization", divide by number of spins
        mag_abs_sq = mag_abs_sq + (abs(magt))**2 !also we are not "per spin-ing" this guy here
     else if(which_lattice==2)then
        mag=mag+sum(dble(spin))/(Nx*Ny)
        mag_sq = mag_sq + (sum(dble(spin)))**2 !we are not "per spin-ing" this guy here
        mag_abs=mag_abs+abs(sum(dble(spin)))/(Nx*Ny)
        mag_abs_sq = mag_abs_sq + (abs(sum(dble(spin))))**2 !we are not "per spin-ing" this guy here
     end if
     !We are only including the magnetization (absolute value) at the end of a full "Monte Carlo step" (nh steps) in the calculation of the average magnetization for that bin. Same with the energy. You sum all the mag, abs(mag) and E values at the end of each mstep, and then divide by msteps, and store the values in mea(1), mea(2) and mea(3).
     Eest=Eest+E/(Nx*Ny) !aggregate the energies per spin at the end of every mstep into Eest
     Esqest=Esqest+(E**2) !do not per spin this guy right now
  end do
  
!Store the bin averages in the mea array
mea(1)=mag/msteps
mea(2)=mag_abs/msteps
mea(3)=Eest/msteps
!Every time this subroutine is called, mea(1), mea(2) and mea(3) get re-assigned. They essentially represent the average values from the latest Monte Carlo process. Recall that to find the fixed point, we have to repeat the Monte Carlo process multiple times, and update the equilibrium probability distribution to be simulated with the latest E_bar value. 

!bin average for squared values.
mea_sq(1) = mag_sq/msteps
mea_sq(2) = mag_abs_sq/msteps
mea_sq(3) = Esqest/msteps

!write(6,*)"mea",mea

end subroutine runbin_renyi



!SUBROUTINE 3 INSIDE "CONTAINS"- measurements
subroutine accmeasure(i)
use configuration !this module defined Nx,Ny,nh,which_lattice,nn,ivic(:,:),ivict(:,:),ixpiy(:),spin(:),temperature,beta,Jh                                                                          
use measurement !this module defined lm,Vm(:),Vm2(:),Vmtemp(:),Vm2temp(:),Vm2temp_error(:),mea(:),E
implicit none
integer(4)i,j,k !is this the same "i" as the "i" that was input into this function?
real(8)magt,magt_abs

!Recall that Vm and Vm2 are 3-component arrays, 3 elements, and were initialized to zero in the initial() function. Ok, accmeasure(i) is called after each runbin(msteps) call. mea contains the bin averages from the latest runbin(msteps) call, i.e. the latest bin averages computed. Vm2 contains the square of these values, which are used for error calculation. Also, recall, the bin averages in mea are "per spin" average values.
Vm=Vm+mea
Vm2=Vm2+mea**2 !element wise
Vm_sq=Vm_sq+mea_sq !new addition



!write(15,*)mea(1),mea(2),mea(3)

if(which_lattice==1.and.Jh>0.0)then !ferromagnetic coupling
   magt=sum(dble(spin))/(Nx*Ny)
   magt_abs=abs(sum(dble(spin)))/(Nx*Ny) !"spin" refers to the spin configuration at the end of the latest bin (msteps * nh total updates corresponding to one bin)
   !magt must be the abs(mag) per spin value for the latest bin
else if(which_lattice==1.and.Jh<0.0)then !forget about the antiferromagnetic coupling for now
   magt=0.0d0
   do k=1,nh
      magt=magt+dble(spin(k))*(-1.0d0)**ixpiy(k)
      !write(6,*)'checking',spin(k),ixpiy(k),(-1.0d0)**ixpiy(k),dble(spin(k))*(-1.0d0)**ixpiy(k)
   end do
   magt=magt/(Nx*Ny)
   magt_abs=abs(magt)/(Nx*Ny)
else if(which_lattice==2)then
   magt=sum(dble(spin))/(Nx*Ny)
   magt_abs=abs(sum(dble(spin)))/(Nx*Ny)
end if

!For now we will not write the equivalent block of code for the squared averages...

!Recall 15 and 20 were created and opened in the ising program
!open(15,file='measurement.dat',status='unknown',form='formatted')
!open(20,file='configuration.dat',status='unknown',form='formatted')
write(15,*)E/(Nx*Ny),magt !E is the energy of the spin configuration after the latest round of msteps * nh updates. measurement.dat is just about keeping track of the values.



!#Current estimates
Vmtemp=Vm/i !temp stands for temporary. Vm sums the bin averages (i is the number of bins completed) for magnetization (Vm(1)), abs(magnetization) (Vm(2)) and energy (Vm(3)), all per spin
Vm2temp=Vm2/i !This is only useful for the calculation of the errors, NOT the squared averages
Vm2temp_error=sqrt(abs(Vm2temp-Vmtemp**2)/dble(i))
Vmtemp_sq=Vm_sq/i ! new addition for squared averages

open(10,file='results.txt',status='replace') !replace means start from scratch, every time accmeasure(i) is called, the results.txt file is rewritten with the latest values
write(10,*)'Current temperature :', temperature
write(10,*)' Number of bins completed : ',i
write(10,*)' ========================================='
write(10,10)'  M         : ',Vmtemp(1),Vm2temp_error(1) !Vmtemp(1), Vmtemp(2) and Vmtemp(3) are the averages of the bin values.
write(10,10)'  M_ABS     : ',Vmtemp(2),Vm2temp_error(2)     
write(10,10)'  E         : ',Vmtemp(3),Vm2temp_error(3)
write(10,10)'  M_SQ      : ',Vmtemp_sq(1)  
write(10,10)'  M_ABS_SQ  : ',Vmtemp_sq(2)
write(10,10)'  E_SQ      : ',Vmtemp_sq(3)
write(10,*)' ========================================='
10 format(1x,a,2f14.8)
close(10)


write(20,*)spin
!write(6,*)sum(dble(spin)/nh)

end subroutine accmeasure
  
END MODULE runbin_accmeasure





!###################################################################################################################################################
!############################################################### FINISHED WITH MODULES #############################################################
!###################################################################################################################################################







!Ok the main program. It calls all the subroutines that we will then need to individually understand. Some of the subroutines were included in the runbin_accmeasure module, because the find_fixed_point subroutine for example uses those subroutines themselves.
program ising
use configuration !so now we can use all the global variables defined within the configuration module
use measurement
use runbin_accmeasure !module that contains the runbin and accmeasure subroutines
implicit none !enforces explicit declaration of all variables- Fortran 90 implicitly treats all variables starting with i, j, etc. as integers for example. "implicit none" avoids this issue.
integer(4) i,j,k,N,iseedd,nther,nrun,msteps,nther_fps,nrun_fps,msteps_fps,upOrDown,num_of_oscillations ! added a bunch of variables needed for the fixed point search.
real(8) :: drand1, rate !I added rate, for the measuring of execution time using system_clock. 
real(8) :: E_bar_init, E_bar_fp !E_bar_init = Initial E_bar value that will parametrize our Renyi ensemble probability distribution, E_bar_fp is the fixed point
integer(8) :: beginning, end !to measure execution time
real(8) m_avg,m_sq_avg,m_error,m_abs_avg,m_abs_sq_avg,m_abs_error,E_avg,E_sq_avg,E_error
real(8) m_sq_avg_CvChi, E_sq_avg_CvChi, E_error_CvChi, m_error_CvChi !new additions for specific heat and magnetic susceptibility calculation tests... testing error calculation as well
character(1000) f0, f1, f2, f3, f4, f5, f6, f7, f8, f9 !strings to be used for file naming purposes...
logical :: exist
real(4) Tmin_forFileWriting, E_bar_init_forFileWriting

!i, j, k should specify site indices, or at least elements of a given array (e.g. num of rows, num of columns in the 2D lattice)
!iseedd is a seed for the random number generator. It is an argument of rand_init which is a function of subroutine that initializes the random number generator drand1. rand_init must be called before drand1 is called.
!nther. nther is the number of times runbin(msteps) is called. It is the number of bins that are being run in the initial "thermalization" phase, i.e. the number of bins that we are running to ensure thermalization occurs. It is the number of bins that are discarded for the averaging
!nrun = number of bins being included in the thermal averaging. We discard nther bins and consider nrun bins once thermalization has (hopefully) been achieved
!msteps = # of Monte Carlo steps per lattice site (total steps = msteps * nh) in a given bin.
!upOrDown: if = -1, all spins down for the initial state. if = +1, all spins up for the initial state.
!drand1 comes from the random.c file. Produces a random floating number between 0 (inclusive) and 1 (non-inclusive).

!Inputs that must be manually input, but we do not manually input anything for our Renyi ensemble simulations. So comment this out.
!write(6,*)'Nx,Ny,which,alpha,temperature,Jh,upOrDown' !Nx,Ny, self-explanatory, which = 1(2) means square(triangular) lattice.
!read(5,*)Nx,Ny,which_lattice,alpha,temperature,Jh,upOrDown
!write(6,*)'num_of_oscillations,nther_fps,nrun_fps,msteps_fps,E_bar_init' !additional shit needed for the fixed point search. "_fps" specifies variables used in the fixed point search
!read(5,*)num_of_oscillations,nther_fps,nrun_fps,msteps_fps,E_bar_init
!write(6,*)'iseedd,nther,nrun,msteps'
!read(5,*)iseedd,nther,nrun,msteps

Nx = 40
Ny = 40
nh = Nx * Ny ! self-explanatory
which_lattice = 1
alpha = 2.0
!temperature = 2.8
Jh = 1
h = 0.0 !including transverse field value of h = 0.0 in file name for consistency purposes
upOrDown = -1
num_of_oscillations = 30 !We will start with this to see if we can reproduce the results we produced using Julia. Moving to num_of_oscillations = 30. See paper, where this is called N_{osc}
nther_fps = 50
nrun_fps = 50
msteps_fps = 50
E_bar_init_forFileWriting = -1.60 * (Nx*Ny)
!E_bar_init = -1.95 * (Nx*Ny) !Trying something new: annealing. We will start with an energy close to the ground state energy for low temperatures, then the annealing begins.
iseedd = 111
nther = 200
nrun = 200
msteps = 600


Tmin_forFileWriting = 0.001
Tmax = 5.000
Tstep = 0.001
 

!Tmin and E_bar_init                                                                                                                                                                                     
write(f0,"('last_tempAndFp_completed.txt')")
inquire(file=f0, exist=exist)
if (exist) then
   open(34, file=f0, status="old",  action="read")
   read(34,*) Tmin, E_bar_init
   Tmin = Tmin + Tstep !next temperature                                                                                                                                                                 
   E_bar_init = E_bar_init + 9.0
   close(34)
else ! initialize from scratch if f0 doesn't exist
   !open(34, file=f0, status="new", action="write")                                                                                                                                                       
   Tmin = 0.001
   E_bar_init = -1.95 * (Nx*Ny)
end if



!From here I want to measure execution time. Trial system_clock
call system_clock(beginning, rate)


!recall, ivic stores the indices associated with all 4 nearest neighbors of each of the nh sites. The first site is the bottom left site, and the site indices are read out in a snaking pattern from left to right, one row at a time.
if(which_lattice==1)then
   nn=4
   allocate(ivic(nh,nn)) !square
   allocate(ixpiy(nh))
   call square() !essentially fills the ivic and ixpiy arrays
elseif(which_lattice==2)then !we don't care about the triangular lattice for now
   nn=6
   allocate(ivic(nh,nn),ivict(nh,4)) !triangular
   call triangular()  
end if !this if statement doesn't need to be in the temperature loop


temperature = Tmin-Tstep !initialization of temperature
!num_of_temperatures = int((Tmax-Tmin)/Tstep)+1


!allocation of variables- we do it here, before the loop over temperatures, because we are having issues
allocate(spin(nh)) !the spin array should hold the spin values associated with each site
lm=3
allocate(mea(lm),mea_sq(lm),Vm(lm),Vm2(lm),Vmtemp(lm),Vm2temp(lm),Vm2temp_error(lm),Vm_sq(lm),Vmtemp_sq(lm))

!Filename creation (i.e. creation of the filename strings) for the various quantities we will calculate as a function of temperature.
! I will add leading 0's ahead of the inclusion of the h, E_bar_init, Tmin, and Tstep values, because so far we have used values for those quantities of 0.___ and Fortran doesn't show the leading 0 to the left of the decimal point. I need it to because I need it to. Now we may change the E_bar_init value to -2.0 or to the result of some annealing process and if we do, we'll have to amend the below code as well.
write(f1,"('m_avg_vs_T_alpha',F0.2,'_',I0,'x',I0,'_h0',F0.1,'_J',F0.1,'_Tmin0',F0.3,'_Tmax',F0.3,'_Tstep0',F0.3,'_numOsc', &
I0,'_EbarInit',F0.2,'_nTherFps',I0,'_nRunFps',I0,'_mStepsFps',I0,'_nTher',I0,'_nRun',I0,'_mSteps',I0,'_iSeedd',I0, &
'_Periodic_Annealing_MonteCarlo_Ising2D_Renyi.txt')") alpha, Nx, Ny, h, Jh, Tmin_forFileWriting, Tmax, Tstep, &
     num_of_oscillations,E_bar_init_forFileWriting/nh,nther_fps,nrun_fps,msteps_fps, nther,nrun,msteps,iseedd

write(f2,"('m_sq_avg_vs_T_alpha',F0.2,'_',I0,'x',I0,'_h0',F0.1,'_J',F0.1,'_Tmin0',F0.3,'_Tmax',F0.3,'_Tstep0',F0.3,'_numOsc',I0, &
'_EbarInit',F0.2,'_nTherFps',I0,'_nRunFps',I0,'_mStepsFps',I0,'_nTher',I0,'_nRun',I0,'_mSteps',I0,'_iSeedd',I0, &
'_Periodic_Annealing_MonteCarlo_Ising2D_Renyi.txt')") alpha, Nx, Ny, h, Jh, Tmin_forFileWriting, Tmax, Tstep, &
     num_of_oscillations,E_bar_init_forFileWriting/nh,nther_fps,nrun_fps,msteps_fps, nther,nrun,msteps,iseedd

write(f3,"('m_error_vs_T_alpha',F0.2,'_',I0,'x',I0,'_h0',F0.1,'_J',F0.1,'_Tmin0',F0.3,'_Tmax',F0.3,'_Tstep0',F0.3,'_numOsc',I0, &
'_EbarInit',F0.2,'_nTherFps',I0,'_nRunFps',I0,'_mStepsFps',I0,'_nTher',I0,'_nRun',I0,'_mSteps',I0,'_iSeedd',I0, &
'_Periodic_Annealing_MonteCarlo_Ising2D_Renyi.txt')") alpha, Nx, Ny, h, Jh, Tmin_forFileWriting, Tmax, Tstep, &
     num_of_oscillations,E_bar_init_forFileWriting/nh,nther_fps,nrun_fps,msteps_fps, nther,nrun,msteps,iseedd

write(f4,"('m_abs_avg_vs_T_alpha',F0.2,'_',I0,'x',I0,'_h0',F0.1,'_J',F0.1,'_Tmin0',F0.3,'_Tmax',F0.3,'_Tstep0',F0.3,'_numOsc',I0, &
'_EbarInit',F0.2,'_nTherFps',I0,'_nRunFps',I0,'_mStepsFps',I0,'_nTher',I0,'_nRun',I0,'_mSteps',I0,'_iSeedd',I0, &
'_Periodic_Annealing_MonteCarlo_Ising2D_Renyi.txt')") alpha, Nx, Ny, h, Jh, Tmin_forFileWriting, Tmax, Tstep, &
     num_of_oscillations,E_bar_init_forFileWriting/nh,nther_fps,nrun_fps,msteps_fps, nther,nrun,msteps,iseedd

write(f5,"('m_abs_sq_avg_vs_T_alpha',F0.2,'_',I0,'x',I0,'_h0',F0.1,'_J',F0.1,'_Tmin0',F0.3,'_Tmax',F0.3,'_Tstep0',F0.3,'_numOsc', &
I0,'_EbarInit',F0.2,'_nTherFps',I0,'_nRunFps',I0,'_mStepsFps',I0,'_nTher',I0,'_nRun',I0,'_mSteps',I0,'_iSeedd',I0, &
'_Periodic_Annealing_MonteCarlo_Ising2D_Renyi.txt')") alpha, Nx, Ny, h, Jh, Tmin_forFileWriting, Tmax, Tstep, &
     num_of_oscillations,E_bar_init_forFileWriting/nh,nther_fps,nrun_fps,msteps_fps, nther,nrun,msteps,iseedd

write(f6,"('m_abs_error_vs_T_alpha',F0.2,'_',I0,'x',I0,'_h0',F0.1,'_J',F0.1,'_Tmin0',F0.3,'_Tmax',F0.3,'_Tstep0',F0.3,'_numOsc', &
I0,'_EbarInit',F0.2,'_nTherFps',I0,'_nRunFps',I0,'_mStepsFps',I0,'_nTher',I0,'_nRun',I0,'_mSteps',I0,'_iSeedd',I0, &
'_Periodic_Annealing_MonteCarlo_Ising2D_Renyi.txt')") alpha, Nx, Ny, h, Jh, Tmin_forFileWriting, Tmax, Tstep, &
     num_of_oscillations,E_bar_init_forFileWriting/nh,nther_fps,nrun_fps,msteps_fps, nther,nrun,msteps,iseedd

write(f7,"('E_avg_vs_T_alpha',F0.2,'_',I0,'x',I0,'_h0',F0.1,'_J',F0.1,'_Tmin0',F0.3,'_Tmax',F0.3,'_Tstep0',F0.3,'_numOsc',I0, &
'_EbarInit',F0.2,'_nTherFps',I0,'_nRunFps',I0,'_mStepsFps',I0,'_nTher',I0,'_nRun',I0,'_mSteps',I0,'_iSeedd',I0, &
'_Periodic_Annealing_MonteCarlo_Ising2D_Renyi.txt')") alpha, Nx, Ny, h, Jh, Tmin_forFileWriting, Tmax, Tstep, &
     num_of_oscillations,E_bar_init_forFileWriting/nh,nther_fps,nrun_fps,msteps_fps, nther,nrun,msteps,iseedd

write(f8,"('E_sq_avg_vs_T_alpha',F0.2,'_',I0,'x',I0,'_h0',F0.1,'_J',F0.1,'_Tmin0',F0.3,'_Tmax',F0.3,'_Tstep0',F0.3,'_numOsc',I0, &
'_EbarInit',F0.2,'_nTherFps',I0,'_nRunFps',I0,'_mStepsFps',I0,'_nTher',I0,'_nRun',I0,'_mSteps',I0,'_iSeedd',I0, &
'_Periodic_Annealing_MonteCarlo_Ising2D_Renyi.txt')") alpha, Nx, Ny, h, Jh, Tmin_forFileWriting, Tmax, Tstep, &
     num_of_oscillations,E_bar_init_forFileWriting/nh,nther_fps,nrun_fps,msteps_fps, nther,nrun,msteps,iseedd

write(f9,"('E_error_vs_T_alpha',F0.2,'_',I0,'x',I0,'_h0',F0.1,'_J',F0.1,'_Tmin0',F0.3,'_Tmax',F0.3,'_Tstep0',F0.3,'_numOsc',I0, &
'_EbarInit',F0.2,'_nTherFps',I0,'_nRunFps',I0,'_mStepsFps',I0,'_nTher',I0,'_nRun',I0,'_mSteps',I0,'_iSeedd',I0, &
'_Periodic_Annealing_MonteCarlo_Ising2D_Renyi.txt')") alpha, Nx, Ny, h, Jh, Tmin_forFileWriting, Tmax, Tstep, &
     num_of_oscillations,E_bar_init_forFileWriting/nh,nther_fps,nrun_fps,msteps_fps, nther,nrun,msteps,iseedd

!Open / Create files
! Testing something new.
!open(25,file=f1,status='unknown')
!open(26,file=f2,status='unknown')
!open(27,file=f3,status='unknown')
!open(28,file=f4,status='unknown')
!open(29,file=f5,status='unknown')
!open(30,file=f6,status='unknown')
!open(31,file=f7,status='unknown')
!open(32,file=f8,status='unknown')
!open(33,file=f9,status='unknown')

print*, "num of temperatures = ", nint((Tmax-Tmin)/Tstep)+1 !Used this line for debugging purposes
do i = 1, nint((Tmax-Tmin)/Tstep)+1
   !OPEN FILES FIRST. We are doing this in the loop from now on. I don't know if this will make things slower or not. But we have to do it here because of the threat of job termination, and for this ridiculous checkpointing that I'm being forced to do. Otherwise, the files that are saved won't contain all the latest data that was produced, I tested this. We need to open and close those files every iteration.
   inquire(file=f1, exist=exist)
   if (exist) then
      open(25, file=f1, status="old", position="append", action="write")
   else
      open(25, file=f1, status="new", action="write")
   end if
   
   inquire(file=f2, exist=exist)
   if (exist) then
      open(26, file=f2, status="old", position="append", action="write")
   else
      open(26, file=f2, status="new", action="write")
   end if
   
   inquire(file=f3, exist=exist)
   if (exist) then
      open(27, file=f3, status="old", position="append", action="write")
   else
      open(27, file=f3, status="new", action="write")
   end if
   
   inquire(file=f4, exist=exist)
   if (exist) then
      open(28, file=f4, status="old", position="append", action="write")
   else
      open(28, file=f4, status="new", action="write")
   end if
   
   inquire(file=f5, exist=exist)
   if (exist) then
      open(29, file=f5, status="old", position="append", action="write")
   else
      open(29, file=f5, status="new", action="write")
   end if
   
   inquire(file=f6, exist=exist)
   if (exist) then
      open(30, file=f6, status="old", position="append", action="write")
   else
      open(30, file=f6, status="new", action="write")
   end if
   
   !E_avg file that we opened earlier (if it already existed) in order to read from it.
   inquire(file=f7, exist=exist)
   if (exist) then
      open(31, file=f7, status="old", position="append", action="write")
   else
      open(31, file=f7, status="new", action="write")
   end if
   
   inquire(file=f8, exist=exist)
   if (exist) then
      open(32, file=f8, status="old", position="append", action="write")
   else
      open(32, file=f8, status="new", action="write")
   end if
   
   inquire(file=f9, exist=exist)
   if (exist) then
      open(33, file=f9, status="old", position="append", action="write")
   else
      open(33, file=f9, status="new", action="write")
   end if
   
   



   temperature = temperature + Tstep 
   beta=1.0d0/temperature
   
   
   call rand_init(iseedd) !Seeding random number generator, see random.c file.
   call initial_groundstate_FM(upOrDown)
   call find_fixed_point(nther_fps,nrun_fps,msteps_fps,num_of_oscillations,E_bar_init)

   !call find_fixed_point(nther_fps,nrun_fps,msteps_fps,num_of_oscillations,E_bar_init)
   !OK, fixed point has been found here. I think Vmtemp(3) contains the fixed point right now. The latest Vmtemp(3) value contains the desired E_bar value

   E_bar_fp = Vmtemp(3) * (Nx*Ny) !this is the fixed point up to an error
   !m_avg = Vmtemp(1) !Wrote this line to debug
   
   !Once the fixed point is found, run the MC process one last time with a larger number of bins and a larger number of steps.



   !Now. We only want to run the Monte Carlo process again if E_bar_fp != -2.0*N

   if (E_bar_fp == -2.0 * (Nx*Ny)) then
      !print*, "SYSTEM STAYS IN GROUND STATE"
      m_avg = upOrDown * (Nx*Ny)
      m_sq_avg = (Nx*Ny)**2
      m_error = 0.0 !we know this will be true without having to do the averaging explicitly 
      m_abs_avg = Nx*Ny
      m_abs_sq_avg = (Nx*Ny)**2
      m_abs_error = 0.0 !we know this will be true without having to do the averaging explicitly 
      E_avg = E_bar_fp
      E_sq_avg = E_avg**2 !it's fine to do this here since we are stuck in the ground state
      E_error = 0.0 !we know this will be true
      ! The above values are not "per spin" values
      
      !call accmeasure(nrun_fps) !just to update output files
      
      !print*, "##### AFTER MONTE CARLO SIMULATION #####"
      !print*, "E_bar per spin, final = ", E_avg / (Nx*Ny)
      !print*, "mag per spin = ", m_avg / (Nx*Ny)
      !print*, "abs(mag) per spin = ", m_abs_avg / (Nx*Ny)
      !print*, "energy error per spin = ", E_error / (Nx*Ny)
      !print*, "abs(mag) error per spin = ", m_abs_error / (Nx*Ny)
      !print*, "Cv = 0.0"
      !print*, "Chi = 0.0"

      !Tracking progress of loop
      open(16,file='trackProgress.txt',status='replace')
      write(16,*) 'T = ', temperature, '   E_bar_fp / N  = ', E_avg / (Nx*Ny)
      close(16)

      !write(16,*) E_avg / (Nx*Ny)
      !16    format(F4.2)
      !10 format(1x,a,2f14.8)
      !16    format(1x,a,2f14.8)
   else
      ! Now do the MC process: thermalization
      ! nther bins needed to achieve thermalization. For each bin, you do msteps * nh local updates (or proposed updates, since not all updates are accepted)
   
      do j = 1,nther
         call runbin_renyi(msteps,E_bar_fp) !one bin average is computed every time runbin(msteps) is called. The bin has msteps * nh local updates (where mstep = Monte Carlo step). The nther bins will be discarded and and will not be used to calculate the final averages, the assumption is that nther bins are needed to achieve thermalization.
      end do

      !write(6,*)"Thermalization ready..."
      ! Once we are here, we assume thermalization has been realized, and we are ready to calculate the final averages.

      !nrun bins
      do j = 1,nrun !nrun bins. Each time runbin(msteps) is called, we get one bin average (for mag and another for E) 
         call runbin_renyi(msteps,E_bar_fp)
         call accmeasure(j) !This will do the averaging of the nrun bin averages.
      end do
      !Everything we calculate and store in Vmtemp, Vm2temp and Vm2temp_error is a per spin value, or in the case of squared averages, a "per spin squared" value.
      m_avg = Vmtemp(1) * nh
      m_sq_avg = Vmtemp_sq(1) !we did not "per spin" the squared averages stored in Vmtemp_sq, so no need to multiply by nh
      m_error = Vm2temp_error(1) * nh !times nh  rather than times nh**2, because m_error is the standard deviation, not the variance. And Vm2temp_error represents the error per spin (sqrt used)
      m_abs_avg = Vmtemp(2) * nh
      m_abs_sq_avg = Vmtemp_sq(2) !we did not "per spin" the squared averages stored in Vmtemp_sq, so no need to multiply by nh
      m_abs_error = Vm2temp_error(2) * nh !see comment I wrote above on m_error line
      E_avg = Vmtemp(3) * nh
      E_sq_avg = Vmtemp_sq(3) !we did not "per spin" the squared averages stored in Vmtemp_sq, so no need to multiply by nh
      E_error = Vm2temp_error(3) * nh !see comment I wrote above on m_error line
      
      !m_avg = Vmtemp(1)
      !m_sq_avg = Vmtemp_sq(1)
      !m_error = Vm2temp_error(1)
      !m_abs_avg = Vmtemp(2)
      !m_abs_sq_avg = Vmtemp_sq(2)
      !m_abs_error = Vm2temp_error(2)
      !E_avg = Vmtemp(3)
      !E_sq_avg = Vmtemp_sq(3)
      !E_error = Vm2temp_error(3)

      
      ! Remember, Vmtemp already contains "per spin" values, so no need to divide by (Nx*Ny) here
      !print*, "##### AFTER MONTE CARLO SIMULATION #####"
      !print*, "E_bar per spin, final = ", E_avg / (Nx*Ny)
      !print*, "abs(mag) per spin = ", m_abs_avg / (Nx*Ny)
      !print*, "energy error per spin = ", E_error / (Nx*Ny)
      !print*, "abs(mag) error per spin = ", m_abs_error / (Nx*Ny)
      !print*, "Cv = ", (E_sq_avg - E_avg**2) / temperature**2
      !print*, "Chi = ", (m_sq_avg - m_avg**2) / temperature
      !Tracking progress of loop
      open(16,file='trackProgress.txt',status='replace')
      write(16,*) 'T = ', temperature, '   E_bar_fp / N  = ', E_avg / (Nx*Ny)
      close(16)
      !write(16,'(f2.3)') 'T = ', temperature, 'E_bar_fp / N  = ', E_avg / (Nx*Ny)
   end if

   
   !OUTPUT FILES WRITTEN HERE AND UPDATED AT EVERY TEMPERATURE
   !CHECK THAT I AM CALCULATING <M^2> and squares of averages correctly
   write(25,*) m_avg
   write(26,*) m_sq_avg
   write(27,*) m_error
   write(28,*) m_abs_avg
   write(29,*) m_abs_sq_avg
   write(30,*) m_abs_error
   write(31,*) E_avg
   write(32,*) E_sq_avg
   write(33,*) E_error

   !close(15)
   !close(20)

   !print*, ""

   !ANNEALING. We recognize that the "next energy level" to be accessed is one that has a single spin flipped compared to the previous state
   E_bar_init = E_bar_fp + 9.0 !This should be good enough to get the correct fixed point for all system sizes, because the "next" energy state to be accessed by the system is at most greater in energy than the last energy state by 8.0 J, by virtue of the energy cost of flipping one spin / forming a minority droplet of size 1.
   

   !if we are here we have completed the current temperature. Update file 24.
   !temperature and E_avg for last temperature that was run (file). If it already existed, then we opened it earlier to read the last temperature and fixed point from it.
   open(24, file='last_tempAndFp_completed.txt', status="replace",action='write') !open this file if it exists, create it if it doesn't
   write(24,*) temperature, E_avg
   close(24)

   !Closing all files here in order to ensure the new data is saved in the event of job termination.
   close(25)
   close(26)
   close(27)
   close(28)
   close(29)
   close(30)
   close(31)
   close(32)
   close(33)

end do

!Close files - maybe we should be opening and closing these files within the loop. Or maybe the files are automatically closed when the program terminates. We can see.
!close(16) !close trackProgress.txt file
!close(24) !close last_temperature_completed.txt file
!close(34) ! closing current temperature file

!Measure execution time, measure it in seconds begin with:
call system_clock(end)
print *, "elapsed time: ", real(end - beginning) / real(rate)

end program ising


!This subroutine stores the indices of nearest neighbors associated with every single site. The sites go from 1 to Nx*Ny = nh. They are counted from left to right starting
!from the bottom row, and then you proceed in a snake pattern as you go from one row to the next. Bottom left corner of the grid is the first site, index = 1 (ii = 1).
!Anyways, ivic is an array of shape nh x nn, or Nx*Ny x nn. Each site has 4 nearest neighbors and the goal is to store the indices of all 4 nearest neighbors associated with a given site.
!So nn=1 refers to the right nn, 2 refers to the up nn, 3 refers to the left nn, 4 refers to the bottom nn. I understand the math below except the line: if (j==Ny) then ivic(ii,2)=ii-Ny*(Nx-1);
!If Ny refers to the vertical coordinate (number of rows) and Nx refers to the horizontal coordinate (number of columns), then for me this should be ii-(Ny-1)*Nx, not ii-Ny*(Nx-1).
!It doesn't matter if the lattice is a square lattice, because Nx=Ny in that case, but if it's not a square lattice I think this line is wrong. I may be thinking about this wrongly but I don't think it matters, we're only ever going to be dealing with square lattices in our analyses I think.
! ....
!ivic stores the indices of the 4 nearest neighbors associated with each site (each site being specified by an index)
subroutine square()
use configuration
implicit none
integer(4)i,j,k,ii
ii=1

do j=1,Ny
   do i=1,Nx
   
     !#right  (1)
       if (i==Nx)then
          ivic(ii,1)=ii-(Nx-1)
       else
          ivic(ii,1)=ii+1
       end if
      ! # up (2)
       if (j==Ny)then
          ivic(ii,2)=ii-Ny*(Nx-1)
       else
          ivic(ii,2)=ii+Nx
       end if
       !#left (3)
       if (i==1)then
          ivic(ii,3)=ii+(Nx-1)
       else
        ivic(ii,3)=ii-1
       end if
       !#down (4)
       if (j==1)then
          ivic(ii,4)=ii+(Ny-1)*(Nx)
       else
          ivic(ii,4)=ii-Nx
       end if


       ixpiy(ii)=i+j
  
       ii=ii+1      
   

   end do  
end do

!do i=1,nh
!write(6,*)i,ivic(i,:)
!end do

end subroutine square

subroutine triangular()
use configuration
implicit none
integer(4)i,j,k,ii
ii=1

do j=1,Ny
   do i=1,Nx
   
     !#right  (1)
       if (i==Nx)then
          ivict(ii,1)=ii-(Nx-1)
       else
          ivict(ii,1)=ii+1
       end if
      ! # up (2)
       if (j==Ny)then
          ivict(ii,2)=ii-Ny*(Nx-1)
       else
          ivict(ii,2)=ii+Nx
       end if
       !#left (3)
       if (i==1)then
          ivict(ii,3)=ii+(Nx-1)
       else
          ivict(ii,3)=ii-1
       end if
       !#down (4)
       if (j==1)then
          ivict(ii,4)=ii+(Ny-1)*(Nx)
       else
          ivict(ii,4)=ii-Nx
       end if
       ii=ii+1       

   end do  
end do

ivic=0

ii=1
do j=1,Ny
    do i=1,Nx

        ivic(ii,1)=ivict(ii,1)
        ivic(ii,3)=ivict(ii,2)
        ivic(ii,4)=ivict(ii,3)
        ivic(ii,6)=ivict(ii,4)
        ivic(ii,2)=ivict(ivict(ii,1),2)
        ivic(ii,5)=ivict(ivict(ii,3),4)
        ii=ii+1
    end do
end do

!do i=1,nh
!write(6,*)i,ivic(i,:)
!end do

end subroutine triangular


!------------------------
!initial configuration- infinite temperature state (random spins)
!-----------------------
subroutine initial_random()

use configuration !this module defined Nx,Ny,nh,which_lattice,nn,ivic(:,:),ivict(:,:),ixpiy(:),spin(:),temperature,beta,Jh
use measurement !this module defined lm,Vm(:),Vm2(:),Vmtemp(:),Vm2temp(:),Vm2temp_error(:),mea(:),E !Vm2temp contains the average of the squares of mag, E, etc. Vm2temp_error contains corresponding errors
implicit none
integer :: i,j 
real(8) :: drand1,bu !unsure about bu for now.

allocate(spin(nh)) !the spin array should hold the spin values associated with each site

do i=1,nh
   ! random initialization at each site
   !Andrew's comment: I get it now. drand1() produces a random float between 0 (inclusive) and 1 (non-inclusive). 2*drand1() will produce a number between 0 (inclusive) and 2 (non-inclusive), 50-50 chance that the number will be between 0 and 1 or between 1 and 2. then int() will cut off the decimal places and you get a 50-50 chance of producing the integer 0 or the integer 1. 2*int(...)-1 produces a -1 with 50% probability and a +1 with 50% probability, clearly.
  spin(i)=2*int(2.*drand1())-1
end do
!so now we have our initial spin configuration. I will have to change it for my case, because Juan chooses an infinite temperature configuration as the initial configuration. I cannot do that for the Renyi ensemble of course, I have to start with a ground state configuration, either all spins up or all spins down.

beta=1/temperature !not sure why this is done here. Don't we do this elsewhere? it is indeed also done in the ising program.

!#initial energy
!I guess this calculates the energy of the initial configuration
E=0.0
do i=1,nh
    bu=0.0
    do j=1,nn/2 ! just need to consider the energy of 2 bonds associated with each site in order to not double count any bonds. Juan goes with 1 (right) and 2 (up) as the bonds of choice.
       bu=bu+spin(ivic(i,j))
       !I get it. bu sums the spin values of the spins to the right and up relative to the desired spin. if bu = 0 at the end, whatever spin(i) is, the contribution to the energy will be 0. If bu = -2 or +2, the contribution to the energy will be nonzero
    end do
    E=E-bu*Jh*spin(i) ! add the energy contribution due to bonds between spin(i) and the right spin and spin(i) and the up spin.
 end do
!first value of E, for the initial spin configuration, successfully calculated. 

!#initialize arrays to accumulate measurements - set all the elements equal to 0. One element will be associated with magnetization, the other with energy
!...
!UPDATE: Let's set lm to 3. I also want regular magnetization in the mix here, not just absolute value. We'll go
 !   element(1) = regular magnetization per spin
 !   element(2) = abs magnetization per spin
 !   element(3) = energy per spin
!lm=3
!allocate(mea(lm),Vm(lm),Vm2(lm),Vmtemp(lm),Vm2temp(lm),Vm2temp_error(lm))
Vm=0.0d0
Vm2=0.0d0
mea=0.0d0
Vm_sq=0.0d0
mea_sq=0.0d0

!write(6,*)'initial energy', E/nh
end subroutine initial_random


!------------------------
!initial configuration- ground state (all spins up or all spins down) for the ferromagnet, I won't write any code for now for the antiferromagnet
!-----------------------
subroutine initial_groundstate_FM(upOrDown)

use configuration !this module defined Nx,Ny,nh,which_lattice,nn,ivic(:,:),ivict(:,:),ixpiy(:),spin(:),temperature,beta,Jh
use measurement !this module defined lm,Vm(:),Vm2(:),Vmtemp(:),Vm2temp(:),Vm2temp_error(:),mea(:),E
implicit none
integer :: i,j,upOrDown !if upOrDown = -1, select all-down ground state as initial state (one of the two symmetry broken ground states). If == +1, select all-up ground state as initial state.
real(8) :: drand1,bu !bu unused for the square lattice.

if (upOrDown==-1) then
   spin = -1 ! all down state
else if (upOrDown==1) then
   spin = 1 ! all up state
end if

!#initial energy- no need to calculate it directly from the configuration in this case since we are working with the ground state. We are assuming Jh is positive (ferromagnet).
E = -Jh * 2 * nh
E_sq = (-Jh * 2 * nh)**2 ! initialize E_sq
!first value of E, for the initial spin configuration, successfully calculated. Henceforth, whenever the measurement module is used, the latest value of E will be the value of E.

!set all the elements equal to 0. One element will be associated with magnetization, the next with the absolute value of mag, the third with energy.
Vm=0.0d0
Vm2=0.0d0
mea=0.0d0
Vm_sq=0.0d0
mea_sq=0.0d0 !new addition for Cv and Chi testing
!write(6,*)'initial energy', E/nh
end subroutine initial_groundstate_FM



subroutine find_fixed_point(nther_fps,nrun_fps,msteps_fps,num_of_oscillations,E_bar_init)
  use configuration !this module defined Nx,Ny,nh,which_lattice,nn,ivic(:,:),ivict(:,:),ixpiy(:),spin(:),temperature,beta,Jh
  use measurement !this module defined lm,Vm(:),Vm2(:),Vmtemp(:),Vm2temp(:),Vm2temp_error(:),mea(:),E
  use runbin_accmeasure !Now we should be able to use the subroutines in the runbin_accmeasure module inside the find_fixed_point subroutine
  implicit none
  integer(4) nther_fps, nrun_fps, msteps_fps, num_of_oscillations, i, j, k
  integer(4) didEbarOscillate, sign_last, sign_next
  real(8) E_bar_init, E_bar_last, E_bar_next
  
  didEbarOscillate = 0 !initialize
  
  !Now do one round of dynamics.
  !Achieve thermalization
  do i = 1,nther_fps
     call runbin_renyi(msteps_fps,E_bar_init)
  end do
  !Bins for averaging
  do i = 1,nrun_fps
     call runbin_renyi(msteps_fps,E_bar_init)
     Vm=Vm+mea
     Vm2=Vm2+mea**2
  end do
  !Calculate the averages and associated errors in case I need them
  Vmtemp=Vm/nrun_fps
  !print*, "Vmtemp = ", Vmtemp
  Vm2temp=Vm2/nrun_fps 
  Vm2temp_error=sqrt(abs(Vm2temp-Vmtemp**2)/dble(nrun_fps))

  E_bar_last = E_bar_init
  E_bar_next = mea(3) * (Nx * Ny) !we need total energy, not energy per spin. Average total energy
  
  sign_next = sign(1.0d0,(E_bar_next-E_bar_last))
  sign_last = sign_next !assume first "sign_last" value is the same as sign_next value. This is just a little trick to proceed.
  
  !Initial E_bar_init value used for MC simulation
  do while (didEbarOscillate < num_of_oscillations)
     !print*, "T = ", temperature, "E_bar_next = ", E_bar_next ! used for debugging
     !These need to be reset to 0 at the start of every MC simulation.
     Vm=0.0d0
     Vm2=0.0d0
     Vm_sq=0.0d0
     !mea=0.0d0 - mea is set to 0 every time runbin_renyi or runbin_gibbs is called. No need to set it to 0 explicitly.
     
     E_bar_last = E_bar_next
     sign_last = sign_next
     !Achieve thermalization first
     do i = 1,nther_fps
        call runbin_renyi(msteps_fps,E_bar_last)
     end do
     !Once equilibrium is realized, run for nrun_fps bins and then proceed
     do i = 1,nrun_fps
        call runbin_renyi(msteps_fps,E_bar_last) !msteps_fps per bin, we will just alter the number of bins
        Vm=Vm+mea
        Vm2=Vm2+mea**2 !element wise
        !call accmeasure(i)
     end do

     !Calculate the averages and associated errors in case I need them
     Vmtemp=Vm/nrun_fps
     Vm2temp=Vm2/nrun_fps 
     Vm2temp_error=sqrt(abs(Vm2temp-Vmtemp**2)/dble(nrun_fps))

     !E_bar_next = mea(3) * (Nx * Ny) ! mea just stored the averages for the last bin... so it was wrong to write use mea(3). Which is why we use Vmtemp(3) for E_bar_next below.

     E_bar_next = Vmtemp(3) * (Nx * Ny)
     sign_next = sign(1.0d0,(E_bar_next-E_bar_last))

     if (E_bar_next-E_bar_last == 0.)then
        !GROUND STATE FOUND OR CORRECT FIXED POINT FOUND, UP TO MONTE CARLO ERROR. This will occur at low temperatures when the fixed point is exactly the ground state.
        exit !break out of "do while" loop
     else if (sign_next * sign_last == -1.)then
        didEbarOscillate = didEbarOscillate + 1
        !print*, "E_bar_next / N = ", E_bar_next / (Nx*Ny)
        !print*, "didEbarOscillate = ", didEbarOscillate
     end if

     

     
     !ANOTHER WAY TO CALCULATE SIGNS
     !Calculate sign_next (assume that the first sign_last has already been determined outside the loop)
     !if (Ebar_next - E_bar_last < 0) then
     !   sign_next = -1
     !else if (E_bar_next - E_bar_last > 0) then
     !   sign_next = 1
     !else if (E_bar_next - E_bar_last == 0) then
        !more code will be needed here I'm sure
     !   sign_next = 0
     !end if
     
     !if (sign_last * sign_next == 0) then !i.e. if the fixed point is at the ground state value of -2.0 * Nx * Ny
     !   exit !no need to proceed in this case
     !   print*, "E_bar_next = ", E_bar_next / (Nx*Ny)
     !else if (sign_last * sign_next == -1)
     !   didEbarOscillate = didEbarOscillate + 1 !increment didEbarOscillate
     !   print*, "E_bar_next = ", E_bar_next / (Nx*Ny)
     !   print*, "didEbarOscillate = ", didEbarOscillate
     !end if


  end do

  
  
  !Once fixed point is found, reset Vm and Vm2 to 0.0 for the final simulation to be done in the main program.
  Vm=0.0d0
  Vm2=0.0d0
  Vm_sq=0.0d0
  !mea already reset to 0 every time runbin_renyi or runbin_gibbs is run. Vm set to 0 only when the initial state is prepared. That's all well and good if you run one MC simulation only. In our case, we run a number of MC simulations just to find the fixed point, then we run a final MC simulation. So we need to reset Vm and Vm2 to 0 at the start of every new MC simulation. 
  
end subroutine find_fixed_point
