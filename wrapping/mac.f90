!> @brief This module holds parameter data for the solver.
!!
!! We store all data not specific to the conserved variables in this module.
!! The variables may be set through subroutines in the solver code.
module dim

    implicit none
    
    integer, parameter  :: nx = 90  !< Number of cells in x-direction.
    integer, parameter  :: ny = 90  !< Number of cells in y-direction.
    
    double precision, parameter     :: dx = 1.0d0/(nx-1.0d0) !< Calculated spatial step size, \f$ \Delta x \f$.
    double precision, parameter     :: dy = 1.0d0/(ny-1.0d0) !< Calculated spatial step size, \f$ \Delta y \f$.
    
    double precision, parameter     :: wconst = 1.0d0/(4.0d0*dx) !< Carried out constant from the calculation of vorticity.
    
    integer             :: t        !< Current timestep value
    double precision    :: dt       !< Current \f$ \Delta t \f$ for timestepping
    
    double precision    :: re       !< Reynolds number of flow, \f$ Re \f$
    double precision    :: r        !< Additional convergence criteria for CFL condition
end module dim

!> @brief The MAC Method solver and it's persistent data.
!! 
module solver
    use dim
    implicit none
    double precision, dimension(nx,ny)  :: u !< Instantaneous velocity in x-direction, \f$ u_{i,j}^{n+1} \f$.
    double precision, dimension(nx,ny)  :: v !< Instantaneous velocity in y-direction, \f$ v_{i,j}^{n+1} \f$.
    double precision, dimension(nx,ny)  :: p !< Scalar pressure at cell center \f$ P_{i,j} \f$
    double precision, dimension(nx,ny)  :: Fn !< See Harlow & Welch 1965
    double precision, dimension(nx,ny)  :: Gn !< See Harlow & Welch 1965
    double precision, dimension(nx,ny)  :: Q !< Pressure source term for Poisson
    
    contains
    
    !> Initialize all variables with zero values
    subroutine initZero()
        u(:,:)  = 0.0d0
        v(:,:)  = 0.0d0
        p(:,:)  = 10.0d0
        Fn(:,:) = 0.0d0
        Gn(:,:) = 0.0d0
        Q(:,:)  = 0.0d0
    end subroutine initZero
    
    !> Accessor for U-Velocity.
    subroutine returnU(array)
        double precision, dimension(nx,ny), intent(out)   :: array !< Function returns a Numpy array of U
        array(:,:) = u(:,:)
    end subroutine returnU
    
    !> Accessor for V-Velocity.
    subroutine returnV(array)
        double precision, dimension(nx,ny), intent(out)   :: array !< Function returns a Numpy array of V
        array(:,:) = v(:,:)
    end subroutine returnV
    
    !> Accessor for Pressure.
    subroutine returnP(array)
        double precision, dimension(nx,ny), intent(out)   :: array !< Function returns a Numpy array of P
        array(:,:) = p(:,:)
    end subroutine returnP
    
    !> @brief Applies ghost cell boundary conditions at the walls. 
    !!
    !! For a rectanglular shape, we see that:
    !! \f[
    !! u|_{L_x=0} = 0\f]
    !! \f[u|_{L_x=1} = 0\f]
    !! \f[v|_{L_y=0} = 0\f]
    !! \f[v|_{L_y=1} = 0
    !! \f]
    !! This mathematically states that a fluid element cannot pass through a wall.
    subroutine ghostCondition()
        use omp_lib
        use dim
        integer     :: i,j

        !X Boundary condition
        !$omp parallel do shared(u) private(i) schedule(dynamic)
        do j=1,ny
            u(1,j)      = 0.0d0
            u(nx,j)     = 0.0d0
        end do
        !$omp end parallel do

        !Y Boundary condition
        !$omp parallel do shared(v)
        do i=1,nx
            v(i,1)      = 0.0d0
            v(i,ny)     = 0.0d0
        end do
        !$omp end parallel do
    end subroutine ghostCondition
    
    !> @brief Applies a moving lid to the top of the cavity, and a no slip condition at every other wall.
    !!
    !! The no slip condition discretized for MAC Method is posed as such:
    !! \f[
    !! u|_{L_y=0} = 0\f]
    !! \f[u|_{L_y=1} = 1\f]
    !! \f[v|_{L_x=0} = 0\f]
    !! \f[v|_{L_x=1} = 0
    !! \f]
    subroutine lidCondition()
        use omp_lib
        use dim
        integer     :: i, j

        !U velocity condition
        !$omp parallel do private(i) shared(u) schedule(dynamic)
        do i=2,nx-1
            u(i,1)      =       - u(i,2)
            u(i,ny)     = 2.0d0 - u(i,ny-1)
        end do
        !$omp end parallel do

        !V Velocity condition
        !$omp parallel do private(j) shared(v) schedule(dynamic)
        do j=2,ny-1
            v(1,j)      =       - v(2,j)
            v(nx,j)     =       - v(nx-1,j)
        end do
        !$omp end parallel do
    end subroutine lidCondition
    
    !> @brief Calculates are maximum allowable timestep, \f$ \Delta t \f$, based on our stability conditions.
    !!
    !! The first stability condition to be evaluated is:
    !! \f[ \Delta t = Re \Delta x \Delta y * r \f] Where \f$ r \f$ must be less than 0.25.
    !! The second stability condition, called the CFL condition, is shown to be:
    !! \f[ \Delta t = \frac{4}{Re (|u|+|v|)^2} \f]
    !! We take the minimum of these two conditions to ensure stability.
    subroutine calcTStep()
        use omp_lib
        use dim
        implicit none
        double precision                :: umax1(2),umax2,vmax1(2),vmax2,dtC,dtR

        !We use different step sizes depending on how long we've been running
        if(t .gt. 10) then
            !Calculate first stability condition
            dtR = r*dx*dy*re

            !get the maximum value of u and v using worksharing
            !$omp workshare
                umax1 = maxval(u)
                umax2 = maxval(umax1)
                vmax1 = maxval(v)
                vmax2 = maxval(vmax1)
            !$omp end workshare

            !Calculate second stability condition
            dtC = (1.0d0 / r) / (re * (abs(umax2) + abs(vmax2))**2 )

            !use the smaller of the two
            dt = min(dtR,dtC)
        else
            !If first few steps then take a very small step
            dt = (r/25.0d0)*dx*dy*re
        end if
    end subroutine calcTStep
    
    !> @brief Calculates intermediate step of \f$ F_n \f$ and \f$ G_n \f$. Reference Harlow & Welch (1965) for explanation.
    !!
    !! These variables are calculated based on the previous time-step's velocity. \f$ F_n \f$ is:
    !! \f[
    !!  F_{i+\frac{1}{2},j}^n = u_{i+\frac{1}{2},j} + \Delta t \left [ 
    !! \frac{u^n_{i+\frac{3}{2},j} - 2u^n_{i+\frac{1}{2},j} + u^n_{i-\frac{1}{2},j}}{Re(\Delta x)^2}
    !! +
    !! \frac{u^n_{i+\frac{1}{2},j-1} - 2u^n_{i+\frac{1}{2},j} + u^n_{i+\frac{1}{2},j-1}}{Re(\Delta y)^2} 
    !! -
    !! \frac{(uv)^n_{i+\frac{1}{2},j+\frac{1}{2}} - (uv)^n_{i-\frac{1}{2},j-\frac{1}{2}}}{\Delta y}
    !! \right ]
    !! \f]
    !! \f$ G_n \f$ is defined as:
    !! \f[
    !!  G_{i,j+\frac{1}{2}}^n = v_{i,j+\frac{1}{2}} + \Delta t \left [ 
    !!  \frac{v^n_{i+1,j+\frac{1}{2}} - 2v^n_{i,j+\frac{1}{2}} + v^n_{i-1,j+\frac{1}{2}}}{Re(\Delta x)^2}
    !!  +
    !!  \frac{v^n_{i,j+\frac{3}{2}} - 2v^n_{i,j+\frac{1}{2}} + u^n_{i,j-\frac{1}{2}}}{Re(\Delta y)^2} 
    !!  -
    !!  \frac{(uv)^n_{i+\frac{1}{2},j+\frac{1}{2}} - (uv)^n_{i-\frac{1}{2},j-\frac{1}{2}}}{\Delta y}
    !!  -
    !! \frac{(v^n_{i,j+1})^2 - (v^n_{i,j})^2}{\Delta y}
    !!  \right ]
    !! \f]
    subroutine calcFnGn()
        use dim
        use omp_lib
        implicit none
        integer         :: i, j, err

        !$omp parallel do private(i,j) schedule(dynamic)
        do i=2,nx-1
            do j=2,ny-1
                Fn(i,j) = u(i,j)+dt*(((u(i+1,j)-2*u(i,j)+u(i-1,j))/(re*dx*dx))+((u(i,j-1)               &
                          -2.0d0*u(i,j)+u(i,j+1))/(re*dy*dy))-(((u(i,j)+u(i+1,j))*(u(i,j)+u(i+1,j))         &
                          -(u(i-1,j)+u(i,j))*(u(i-1,j)+u(i,j)))/(4.0d0*dx))-((((u(i,j)+u(i,j+1))*(v(i+1,j)  &
                          +v(i,j)))-((u(i,j-1)+u(i,j))*(v(i+1,j-1)+v(i,j-1))))/(4.0d0*dy)))
                Gn(i,j) = v(i,j)+dt*(((v(i+1,j)-2.0d0*v(i,j)+v(i-1,j))/(re*dx*dx))+((v(i,j-1)               &
                          -2.0d0*v(i,j)+v(i,j+1))/(re*dy*dy))-(((v(i,j+1)+v(i,j))*(v(i,j+1)+v(i,j))         &
                          -(v(i,j-1)+v(i,j))*(v(i,j-1)+v(i,j)))/(4.0d0*dy))-((((u(i,j)+u(i,j+1))*(v(i+1,j)  &
                          +v(i,j)))-((u(i-1,j)+u(i-1,j+1))*(v(i,j)+v(i-1,j))))/(4.0d0*dx)))
            end do
        end do
        !$omp end parallel do
    end subroutine calcFnGn
    
    !> @brief \f$ Q_n \f$ is calculated as a source for the RHS of the pressure Poisson equation.
    !!
    !! It is calculated based on the values of \f$ F_n \f$ and \f$ Q_n \f$ calculated in the routine calcfngn().
    !! \f[
    !! Q^n_{i,j} = \left [ \frac{P_{i-1,j} -2P_{i+1,j}}{(\Delta x)^2}
    !!    +
    !!  \frac{P_{i,j-1} - 2P_{i,j} + P_{i,j+1}}{(\Delta y)^2}\right ]^{n+1}
    !! \f]
    !! which may be discretized to show:
    !! \f[
    !! Q^n_{i,j} = \frac{1}{\Delta t} \left [ 
    !!  \frac{F^n_{i+\frac{1}{2},j} - F^n_{i-\frac{1}{2},j}}{\Delta x}
    !!  +
    !!  \frac{G^n_{i,j+\frac{1}{2}} - G^n_{i,j-\frac{1}{2}}}{\Delta y}
    !!  \right ].
    !! \f]
    subroutine calcQn()
        use dim
        use omp_lib
        integer                                         :: i, j, err

        !$omp parallel do private(i,j) schedule(dynamic)
        do i=2,nx-1
            do j=2,ny-1
                Q(i,j) =((Fn(i,j)-Fn(i-1,j)+Gn(i,j)-Gn(i,j-1))/(dt*dx))
            end do
        end do
        !$omp end parallel do
    end subroutine calcQn
    
    !> @brief This subroutine calculates velocity of the field and stores the data in the module.
    !!
    !! Values of the velocity may be accessed in Python from the routines returnu() and returnv(). 
    !! The velocity is calculated by the equations:
    !! \f[
    !! u_{i+\frac{1}{2},j}^{n+1} = F_{i+\frac{1}{2},j}^n - \frac{\Delta t}{\Delta x}(P_{i+1,j}^{n+1} - P_{i,j}^{n+1})
    !! \f]
    !! \f[
    !!  v_{i+,j+\frac{1}{2}
    !!    }^{n+1} = G_{i,j+\frac{1}{2}
    !!    }^n - \frac{\Delta t}{\Delta x}(P_{i,j+1}^{n+1} - P_{i,j}^{n+1})
    !! \f]
    subroutine calcVel()
        use dim
        use omp_lib
        integer                                         :: i, j, err

        !$omp parallel do private(i,j) schedule(dynamic)
        do i=2,nx-2
            do j=2,ny-1
                u(i,j) = Fn(i,j) - ( p(i+1,j) - p(i,j) ) * (dt/dx)
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i,j) schedule(dynamic)
        do i=2,nx-1
            do j=2,ny-2
                v(i,j) = Gn(i,j) - ( p(i,j+1) - p(i,j) ) * (dt/dy)
            end do
        end do
        !$omp end parallel do

    end subroutine calcVel
    
    !> @brief Calculates vorticity based on velocity field
    !!
    !! Vorticity is mathematically defined as half the curl of the velocity field.
    !! \f[
    !! \omega = \frac{1}{2} \nabla \times \bar U
    !! \f]
    subroutine vort(w)
        use dim
        implicit none
        double precision, dimension(nx,ny), intent(out) :: w
        integer                                         :: i,j
        
        w(:,:) = 0.0d0
        do j=2,ny-1
            do i=2,nx-1
                w(i,j) = wconst * (v(i+1,j) - v(i-1,j) - u(i,j+1) + u(i,j-1))
            end do
        end do
    end subroutine vort   
    
    !> @brief Poisson solver for the pressure field
    !! 
    !! The equation for the pressure field is defined to be:
    !! \f[
    !!  \nabla^2 P = Q_n
    !! \f]
    !! For two dimensions the equation becomes
    !! \f[
    !!  \frac{\partial P}{\partial x} + \frac{\partial P}{\partial y} - Q_n = 0
    !! \f]
    !! In discretized form, using an explicit Gauss-Seidel method, we see
    !! \f[
    !! P^{n+1}_{i,j} = \frac{1}{4} \left [ P^{n+1}_{i-1,j} + P^{n+1}_{i,j-1} 
    !! + P^n_{i+1,j} + P^n_{i,j+1} - (\Delta x)^2 Q^n_{i,j}\right ]
    !! \f]
    subroutine poisson()
        use dim
        implicit none
        double precision, parameter                     :: rf = 1.60d0
        double precision, parameter                     :: conv = 0.00010d0
        integer, parameter                              :: itmax = 200
        double precision                                :: change, pold, po, ch
        integer                                         :: iter, im, jm, k
        integer                                         :: i,j

        iter = 0
        change = 1.0d0
        do while((change .ge. conv))
            iter = iter + 1
            change = 0.0d0

            do i=2,nx-1
                do j=2,ny-1
                    pold=p(i,j)
                    p(i,j) = 0.250d0*((p(i-1,j)+p(i,j-1)+p(i+1,j)+p(i,j+1))-(q(i,j)*dx*dx))
                    p(i,j) = pold + rf*(p(i,j)-pold);

                    !Calculates change
                    if(pold .ne. 0) ch = abs((p(i,j)-pold)/pold)

                    !Stores greatest change value
                    if(ch .gt. change) then
                        change = ch
                        im = i
                        jm = j
                        po = pold
                    end if
                end do
            end do

            !update boundaries on poisson solver
            do k=1,nx
                p(k,1)      = p(k,2)       - ((2.0 * v(k,2))       / (re*dx))
                p(k,ny)     = p(k,ny-1) + ((2.0 * v(k,ny-2)) / (re*dx))
                p(1,k)      = p(2,k)       - ((2.0 * u(2,k))       / (re*dx))
                p(nx,k)     = p(nx-1,k) + ((2.0 * u(nx-2,k)) / (re*dx))
            end do

            if(iter .gt. itmax) then
                goto 100
            end if
        end do
100     continue        
    end subroutine poisson
    
    !> @brief Parallelized Poisson solver
    !!
    !! @todo Parallel Poisson solver does not currently work due to F2Py bug with OpenMP barriers.
    !! Use the function poisson()
    subroutine parPoisson()
        use omp_lib
        use dim
        implicit none

        double precision, parameter                     :: rf = 1.60d0
        double precision, parameter                     :: conv = 0.00010d0
        integer, parameter                              :: itmax = 500
        double precision                                :: pold, po, ch
        integer                                         :: iter, im, jm, k
        integer                                         :: i,j,tstep
        integer                                         :: rowPer, colPer, numThreads, istart, iend, myid, kstart, kend
        double precision                                :: change, changeThread

        iter = 0
        change = 1.0d0

        !For OpenMP
        numThreads = 4
        call omp_set_num_threads(numThreads)

        !calculate sizes of blocks
        rowPer = int(ny/numThreads)
        colPer = int(nx/numThreads)

        !Fork the threads
        !$omp parallel private(myid, istart, iend, i, j, k,t, iter, pold, changeThread) shared(change)
        myid = omp_get_thread_num()

        istart = myid*rowPer + 1
        if(myid == 0) istart = istart + 1
        iend   = (myid*rowPer) + rowPer
        if(myid == numThreads-1) iend = iend - 1

        !Main iterative loop
        do while(iter .lt. itmax)
            iter = iter + 1

            !$omp single
            change = 0.0d0
            !$omp end single
            !$omp barrier
            !$omp flush

            changeThread = 0.0d0
            do j=istart,iend
                do i=2,nx-1
                    pold=p(i,j)
                    p(i,j) = 0.250d0*((p(i-1,j)+p(i,j-1)+p(i+1,j)+p(i,j+1))-(q(i,j)*dx*dx))
                    p(i,j) = pold + rf*(p(i,j)-pold);

                    !Calculates change
                    if(pold .ne. 0) changeThread = max(changeThread,abs((p(i,j)-pold)/pold))
                end do
            end do

            !calculates change
            !$omp critical
            change = max(change,changeThread)
            !$omp end critical
            !$omp barrier

            if(change .lt. conv) exit

            !$omp barrier

            !update boundaries on poisson solver
            !$omp master
            do k=1,nx
                p(k,1)     = p(k,2)       - ((2.0d0 * v(k,2))       / (re*dx))
                p(k,ny) = p(k,ny-1) + ((2.0d0 * v(k,ny-2)) / (re*dx))
                p(1,k)     = p(2,k)       - ((2.0d0 * u(2,k))       / (re*dx))
                p(nx,k) = p(nx-1,k) + ((2.0d0 * u(nx-2,k)) / (re*dx))
            end do
            !$omp end master

        end do
100     continue
        !$omp barrier
        !$omp end parallel
    end subroutine parPoisson
    
end module solver