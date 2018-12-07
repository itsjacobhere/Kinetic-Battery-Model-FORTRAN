!Source code written by Jacob Turner finished on 10/28/2018
!This program uses runge kutta fehlberg algorithm to evaluate the system of differential eqns
!that describe the movement of charge through a Battery. The system is referenced via a function
!and is derived from the Kinetic Battery Model (KiBaM). The model's charge in the wells is written
!to a CSV file for reading from an R program so that it can be plotted. Depending on user preference
!a sensitivity analysis can also be performed on the charge distribution, potential difference
!(voltage), internal resistence, and the battery run-time required from the battery. The results
!for all of these analyses are written to CSV files to be read and plotted in R studio.
!___________________________________________________________________________________________________
!This module stores and initializes the global variables needed throughout the program.
module var
    !this module stores the variables for use throughout the program.
    implicit none
    integer, parameter :: wp = selected_real_kind(15) !precision
    !subroutine variables:
    real(wp) :: tstart, tend ,h, hmin, hmax
    real(wp), parameter :: epsi1 = 0.000001_wp, epsi2 = 0.00001_wp
    integer :: l, j, exitflag, neqn, w
    integer, parameter :: nsteps = 11000 !size of array to store model results.
    !function/program variables:
    real(wp) :: h2 ! height of the bound charge well.
    real(wp) :: h1 !height of the available charge well.
    real(wp) :: c, I, R, k, v, i0, tgoalhrs
    ! c charge distribution between the wells.
    ! k is fixed conductance equal to 1/resistence.
    ! v is the voltage.
    ! i is current as a function of time constant discharge current i(t) = I at time t
    ! tgoalhrs is the time required or the battery run-time needed.
end module var
!___________________________________________________________________________________________________
!This module contains the RKF algorithm used to solve ODE's
Module Rung_Kutta
    use var
    implicit none
    contains
    subroutine rungk(tstart, tend, y, h, hmin, hmax, D, exitflag, neqn)
        implicit none
        integer, intent(in) :: neqn
        integer, intent(out) :: exitflag
        real(wp), intent(in) :: hmin, hmax, tstart, tend
        real(wp), intent(inout) :: h
        !local vars used only within subroutine
        real(wp) :: t, e, hsave
        real(wp), dimension(neqn) :: k1,k2,k3,k4,k5,k6,y4,y5,ysave
        !see page 747 Chapra & canale, for Runge-Kutta Fehlberg walkthrough and values.
        real(wp), parameter :: &
        c1 = 37._wp/378._wp,       c2 = 250._wp/621._wp,    c3 = 125._wp/594._wp,       &
        c4 = 512._wp/1771._wp,     c5 = 2825._wp/27648._wp, c6 = 18575._wp/48384._wp,   &
        c7 = 13525._wp/55296._wp,  c8 = 277._wp/14336._wp,  c9 = 1._wp/4._wp,           &
        a1 = 1._wp/5._wp,          a2 = 3._wp/10._wp,       a22= 3._wp/40._wp,          &
        a3 = 9._wp/40._wp,         a4 = 3._wp/5._wp,        a5 = -9._wp/10._wp,         &
        a6 = 6._wp/5._wp,          a7 = -11._wp/54._wp,     a8 = 5._wp/2._wp,           &
        a9 = -70._wp/27._wp,       a10= 35._wp/27._wp,      a11= 7._wp/8._wp,           &
        a12= 1631._wp/55296._wp,   a13= 175._wp/512._wp,    a14= 575._wp/13824._wp,     &
        a15= 44275._wp/110592._wp, a16= 253._wp/4096._wp
        real(wp), dimension(:), intent(inout) :: y
        !interface to function that contains the system of diffs
        interface
            Function D(t,y,neqns)
                use var
                implicit none
                integer, intent(in) :: neqns
                real(wp), intent(in) :: t
                real(wp), dimension(:), intent(in) :: y
                real(wp), dimension(neqns) :: D
            end function D
        end interface
        !checking to make sure answer makes sense
        !if (tstart>tend) .or. (epsi1>epsi2) .or. (hmin>hmax) exit
        if (tstart>tend) then
            exitflag=1
            return
        else if (epsi1>epsi2) then
            exitflag=2
            return
        else if (hmin>hmax) then
            exitflag=3
            return
        end if 
        if (h>(tend-tstart)) then
            h = tend-tstart
        end if
        hsave = h
        !t = tstart !time begins at zero
        do 
            l=l+1
            ysave=y(:neqn)!!
            !RFK fehlberg eqns:
            k1 = D(t,y,neqn)
            k2 = D((t+(a1*h)), (y+(a1*k1*h)), neqn)
            k3 = D((t+(a2*h)), (y+(a22*k1*h)+(a3*k2*h)), neqn)
            k4 = D((t+(a4*h)), (y+(a2*k1*h)+(a5*k2*h)+(a6*k3*h)), neqn)
            k5 = D((t+h), (y+(a7*k1*h)+(a8*k2*h)+(a9*k3*h)+(a10*k4*h)), neqn)
            k6 = D(t+(a11*h), y+(a12*k1*h)+(a13*k2*h)+(a14*k3*h)+(a15*k4*h)+(a16*k5*h), neqn)
            !4th and 5th order estimate equations
            y4 = y+(c1*k1*h)+((c2*k3*h)+(c3+k4*h)+(c4*k6*h))
            y5 = y+((c5*k1*h)+(c6*k3*h)+(c7*k4*h)+(c8*k5*h)+(c9*k6*h))
            !compute error estimate differenece between 4th and 5th order estimate.
            e = maxval((abs(y5-y4))/y5) !!
            !controlling error in RKF algorithm
                if (e>epsi2) then
                    if (abs(hmin-h)<0.000001) then
                        exitflag=4
                        return
                    !end if
                    h = h/2._wp
                    if (h<hmin) then 
                        h=hmin
                    end if
                    y = ysave
                else
                    t = t + h
                    y = y5
                    if(e<epsi1) then
                        h = h*2 !double step size
                        if (h>hmax) h = hmax
                    end if
                    !within this range calculations are no longer necessary
                    if(abs(t-tend)<1.e-6) then
                        exitflag = 10
                        return
                    end if
                    !t cannot exceed tend, duh
                    if (t+h>tend) then
                        hsave = h
                        h = tend-t
                        exit
                    end if
                end if
            end if
        end do
    end subroutine rungk
end module Rung_Kutta
!_________________________________________________________________________________________________________________________________________
program charge
    use var
    use Rung_Kutta
    implicit none
    interface
        Function D(t,y,neqns)
            use var
            implicit none
            integer, intent(in) :: neqns
            real(wp), intent(in) :: t
            real(wp), dimension(:), intent(in) :: y
            real(wp), dimension(neqns) :: D
        end function D
    end interface
    real(wp), dimension(:), allocatable :: y
    real(wp), dimension(:,:), allocatable :: f 
    real(wp) :: capacity, batterysize,timeinterval
    real(wp) :: initial_available, initial_bound, tgoal
    integer :: ierror
    character(len=40) :: oldcap="output_original_capacity.csv",afmt &
    ,newcap= "output_NEW_capacity.csv", search = "search.csv", answer
    character (len=15) :: c1 = "c1.csv", r1 = "r1.csv", v1 = "v1.csv",t1 = "t1.csv"
    character(len=3) :: sensitivity
!__________________________________________________________________________________________________
    !create output files to write results to.
    open(unit=20,file=oldcap,status="replace",iostat=ierror)
    open(unit=30,file=newcap,status="replace",iostat=ierror)
    open(unit=40,file=search,status="replace",iostat=ierror)
    !initialization of array for storing model results
    allocate(f(nsteps+1,4)) !allocate final array used for writing to file.
    neqn = 2 !number of equations in system
    !which variable do you want to perform the sensitivity analysis on?
    do
        if (w /= 0) then
            write(*,*) "sensitivity on, ",sensitivity, "completed."
            write(*,*) "want to repeat?"
            read(*,*) answer
            if (answer == "no") then
                write(*,*) "stopping"
                STOP
            end if
        end if
        !user input
        write(*,*) "Do you want to perform a sensitivity analysis on c, V, R or t?"
        write(*,*) "('no' will run the program with all parameters unchanged.)"
        write(*,*)
        read(*,*) sensitivity
        !ask the user if they want to perform a sensitivity analysis
        if (sensitivity == "c") then
            open(unit=50,file = c1 ,status="replace",iostat=ierror)
            write(50,*) "percent,capacity,c"
        else if (sensitivity == "V" .or. sensitivity == "v") then
            open(unit=60,file = v1 ,status="replace",iostat=ierror)
            write(60,*) "percent,capacity,voltage"
        else if (sensitivity == "R" .or. sensitivity == "r") then
            open(unit=70,file = r1 ,status="replace",iostat=ierror)
            write(70,*) "percent,capacity,resistance"
        else if (sensitivity == "t") then
            open(unit=80,file = t1 ,status="replace",iostat=ierror)
            write(80,*) "percent,capacity,tgoal"
        end if
        
        if (sensitivity == "no") then
            w = 0
            go to 100
        end if
        
        do w = 50, 150, 1 !percent range to of original parameter for sensitivity analysis.
            !forumulas:
            !capacity = total battery capacity
            !typical phone battery is 1500 - 2800 mAh will test various values, 
            !this is initial guess for capacity
            if (sensitivity /= "no") then
                write (*,*) w, "%"
            end if
            100 continue
            batterysize = 1000 !mAh
            !parameters: 
            !c, R, v, I
            !c = ratio of charge in available well to bound well
            !R = internal resistence of the battery
            !v = voltage provided from the battery, used to determine constant discharge current
            !********************************************************************************
            !the sensitivity analysis will be on vars: c, R, v and also the time required
            c = 0.65! charge distrubution ratio between the wells
            !c changes the steepness of the initial exponential curve in the R Plot
            R =100!mOHM !increases the battery life and lowers the height of the curve
            !increased voltage increases the length of the initial drop in available charge
            v = 4 ! 4.2 reasonable range for phones is typically between 3.7 to 4.2, split the difference
            !how long must the battery last with  these parameters?
            tgoalhrs = 14 !hrs
            tgoal = 60*tgoalhrs ! 840 minutes which is 14 hours
            !calc for sens analysis
            !if statements before to select variable
            if (sensitivity == "c") then
                c = c*((real(w))/100)
            else if (sensitivity == "V" .or. sensitivity == "v") then
                v = v*((real(w))/100)
            else if (sensitivity == "R" .or. sensitivity == "r") then
                R = R*((real(w))/100)
            else if (sensitivity == "t") then
                tgoal = tgoal*((real(w))/100)
            end if
            k = 1._wp/R !fixed conductance; defined as inverse resistance. (internal resistence of battery)
            i0 = (v/R) ! initial current, use this for a decaying current with time
            !notes:
            !i = (v/R) !*60!1000 !use this for constant current discharge
            !1 Ah = 3600 coulombs
            !C = Q/v ! capacitance = charge/voltage (farads)
            !I = Q/t ! current = charge/time (amps)
            !converts batterysize in mAh to a charge in coulombs
            capacity = (batterysize/1000)*3600!/1200 !convert to coulomb from mAh
            !more notes:
            !y1(0) = c*C 
            !y2(0) = (1-c)*C
            !initial_available = c*capacity
            !Batterysize = ((initial_available/c)/3600)*1000
            !Bound->available->load i(t)
            y = [c*capacity ,(1-c)*capacity] !initial conditions coulombs of charge in well eavh well
            !(initial_available/c)*(1-c)
            initial_available = y(1)
            initial_bound = y(2)
            tstart = 0._wp !start at time zero
            tend = nsteps ! minutes
            timeinterval = 1._wp !(tend-tstart)/nsteps !time interval time advances by 1 unit (minute)
            !setting step size for subroutine
            hmax = timeinterval
            hmin = hmax/100._wp !tminimum step size for subroutine
            h = hmax !start at largest step size for RKF
            !write out parameters for user only if sensitivity analysis is not performed.
            if (sensitivity == "no") then
                write(*,*) "-------------------------------------------------------------------------------------"
                write(*,*) "Available and Bound charge well over time"
                write(*,*) "initial charge in available well= ", initial_available !coulombs
                write(*,*) "initial charge bound well = ", initial_bound !coulombs
                write(*,*) "capacity = ", capacity, " coulombs"
                write(*,*) "k=",k !1/ohms
                write(*,*) "R=",R !ohms
                write(*,*) "c=",c !ratio
                write(*,*) "V=",v !volts
                write(*,*) "I0=",i0 !Amps = coulomb / second
                write(*,*) "The battery must run for, ", tgoalhrs, "hours"
            end if
            
            do j = 1,nsteps
                if (y(1) < 0.15) then
                    exitflag=6
                    exit
                    return
                end if
                tend = tstart + timeinterval
                call rungk(tstart,tend,y,h, hmin, hmax, D, exitflag, neqn)
                f(j,:) = [tend,y(1),y(2),y(1)+y(2)] !store results
                tstart = tend
            end do
            
            write(20,*) "Time (minutes),Available (Coulombs),Bound (Coulombs),Total Charge (Coulombs) "
            afmt = '(F7.2,",",F15.4,",",F15.4,",",F15.4)'
            
            if ((w-100)==0 .or. (w == 0)) then !(w-100)=0 when nothing is changed in sensitivity
                do j =1,nsteps
                    write(20,afmt) f(j,:)
                end do
            end if
            
            if (sensitivity == "no") then
                write(*,*) "-------------------------------------------------------------------------------------"
                write(*,*)
                write(*,*) "Report of initial guess run:"
                write(*,*)
                write(*,*) "Battery Length: ", tend," minutes or,", tend/60,"Hours." 
                write(*,*)
                write(*,*) "Minimum battery size:", batterysize," mAh."
                write(*,*)
                write(*,*) "Initial available charge:", initial_available," Coulombs."
                write(*,*)
                write(*,*) "Initial bound charge:" , initial_bound," Coulombs."
                write(*,*) 
                write(*,*) "Total initial charge:, ",initial_available+initial_bound," Coulombs."
                write(*,*)
                WRITE(*,*) "program finished writing to external file."    
                write(*,*) "-------------------------------------------------------------------------------------"
                write(*,*)
                write(*,*) "Begin search algorithm for charge needed in available well"
                write(*,*)
            end if
            
            searchalgo :do !searches for the battery capacity until conditions are met.
                if ( Abs(tend - tgoal) <= 1 ) exit !5
                !searching for when the battery runs out at the desired time.
                if (tend > tgoal) then
                    initial_available = initial_available - 11
                    write(40,*) initial_available, tend, "above time needed"
                else if ( tend < tgoal) then
                    initial_available = initial_available + 8.1
                    write(40,*) initial_available, tend, "below time needed"
                end if
                
                tstart=0
                y = [initial_available,initial_bound]
                
                do j = 1,nsteps
                    if (y(1)<0.15) then
                        exitflag=6
                        exit
                        return
                    end if
                    tend = tstart + timeinterval
                    call rungk(tstart,tend,y, h, hmin, hmax, D, exitflag, neqn)
                    f(j,:) = [tend,y(1),y(2),y(1)+y(2)]
                    tstart = tend
                end do
            end do searchalgo
            
            if (sensitivity == "no") then
                write(*,*)
                Write(*,*) " Success! results written to CSV file."
                write(*,*)
            end if
            
            if ((w-100)==0 .or. (w == 0)) then
                write(30,*) "Time(minutes),Available(Coulombs),Bound(Coulombs),Total Charge(Coulombs) "
                do j =1,nsteps
                    write(30,afmt) f(j,:) !writing out model results of charge over time.
                end do
            end if
            
            batterysize = ((initial_available/c)/3600)*1000
            initial_bound = (initial_available/c)*(1-c)
            
            if (sensitivity == "no") then
                write(*,*) "Final Report:"
                write(*,*)
                write(*,*) "Battery Length: ", tend," minutes or,", tend/60,"Hours." 
                write(*,*)
                write(*,*) "Minimum battery size:", batterysize," mAh."
                write(*,*)
                write(*,*) "Initial available charge:", initial_available," Coulombs."
                write(*,*)
                write(*,*) "Initial bound charge:" , initial_bound," Coulombs."
                write(*,*) 
                write(*,*) "Total initial charge:, ",initial_available+initial_bound," Coulombs."
                write(*,*) "-------------------------------------------------------------------------------------"
            end if
            
            if (sensitivity =="c") then
                write(50,'(I0,",",f10.4,",",f10.4)') w-100, batterysize ,c
            else if (sensitivity == "V".or. sensitivity == "v") then
                write(60,'(I0,",",f10.4,",",f10.4)') w-100, batterysize ,V
            else if (sensitivity == "R".or. sensitivity == "r") then
                write(70,'(I0,",",f10.4,",",f10.4)') w-100, batterysize ,R
            else if (sensitivity == "t") then
                write(80,'(I0,",",f10.4,",",f10.4)') w-100, batterysize ,tgoal
            end if
            if (sensitivity == "no") go to 200
        end do
    end do
    
    if (sensitivity == "no") then
        200 continue
        write(*,*) "sensitivity analysis not performed." 
        write(*,*) 
    end if
    
    write(*,*) "new capacity values written in: ", newcap
    write(*,*) 
    !report the reason for exiting
    if (exitflag==1) write(*,*) "tstart cannot be greater than tend"
    if (exitflag==2) write(*,*) "epsi1 must be smaller than epsi2"
    if (exitflag==3) write(*,*) "hmin should be less than hmax"
    if (exitflag==4) write(*,*) "step size (h) is too small"
    if (exitflag==5) write(*,*) "too many iterations in subroutine"
    if (exitflag==6) write(*,*) "Program is finished running."
    stop
end program charge
!__________________________________________________________________________________________________
Function D(t,y,neqns) 
    use var
    implicit none
    integer, intent(in) :: neqns
    real(wp), intent(in) :: t
    real(wp), dimension(:), intent(in) :: y
    real(wp), dimension(neqns) :: D
    !function relationships:
    h1 = y(1)/c
    h2 = y(2)/(1-c)
    !decaying current given initial current
    i = i0*t*EXP(-t/(tgoalhrs*60)) ! = tgoalhrs to minutes
    !System of Differential Equations:
    D(1) = -i + k*(h2 - h1)
    D(2) = -k*(h2 - h1)
end function D