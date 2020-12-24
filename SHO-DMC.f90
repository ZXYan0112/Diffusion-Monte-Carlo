program dmc_ho

   implicit none

   integer, parameter :: num = 300, ndim = 1
   integer, parameter :: nwalker_aim = 300
   integer, parameter :: mc_steps = 4000
   real(8), parameter :: x_max = 1.0d0
   real(8), parameter :: alpha = 0.1d0
   real(8), parameter :: gamma1 = 0.5d0
   real(8), parameter :: tau = 0.05d0
   real(8), parameter :: E0 = 0.5d0
 
   integer :: nw(mc_steps) = 0
   real(8) :: walkerV(mc_steps) = 0.0d0
   integer:: nwalker, nwalker_old, n_die, n_birth
   integer :: istep, iwalker, icopy, is, idim1
   real(8) :: time_begin, time_end
   real(8) :: r, force, el1, el2, tmp, E_trial, tmpg(ndim), vtot, ktot
   real(8), allocatable :: x(:,:), x_store(:,:)
   real(8), allocatable :: q(:), v_vec(:)
   integer, allocatable :: s(:)   

   call CPU_TIME(time_begin)
   call random_seed ()

   E_trial = E0

   allocate(x(num,ndim))
   allocate(x_store(num,ndim))

   !--- Put the walkers at random positions in configurational space ---!
   call random_number(x)
   x(:,:)=(x(:,:)*2.0-1.0)*x_max
   x_store = x   
   print *, "Path initialization has been done."

   nwalker = num

   do istep = 1, mc_steps
    allocate(s(nwalker))
    allocate(q(nwalker))
    nwalker_old = nwalker

    do iwalker = 1, nwalker     
      !--- Shift walker from its position R to a new position R' according to the Gaussian transition probability ---!
      do idim1 = 1, ndim
       call gauss_rand(tmpg(idim1))
      enddo
      x(iwalker,:) = x(iwalker,:) + tmpg(:)*sqrt(2.0d0*gamma1*tau)  
  

      !--- Evaluate the q and s = q + r ---!
      vtot = 0.0d0
      do idim1 = 1, ndim
          vtot = vtot + v(x(iwalker,idim1))
      enddo
      call random_number(tmp)
      q(iwalker) = exp(-tau*(vtot-E_trial))

      if (q(iwalker) < tmp) then
        s(iwalker) = 0  
      else
        if (q(iwalker)-1.0d0 > tmp) then
          s(iwalker) = 2
        else
          s(iwalker) = 1
        endif
      endif
    enddo
    print *, "Shifting walkers and Evaluations of q(and s) have been done at istep =", istep

    x_store = x

    nwalker = sum(s)
    nw(istep) = nwalker

    print *, "Current number of walkers is", nwalker
 
    x = 0.0d0   
    deallocate(x)
    allocate(x(nwalker,ndim))    
    x = 0.0d0

    !--- Eliminate the walker or create new ones at R' depending on s ---!
    icopy = 1
    n_die = 0
    n_birth = 0
    do iwalker = 1, nwalker_old
      if ((s(iwalker)) > 0) then     
        do is = 1, s(iwalker) 
          x(icopy,:) = x_store(iwalker,:)
          icopy = icopy + 1
          n_birth = n_birth + 1
        enddo
        n_birth = n_birth - 1
      else
        n_die = n_die +1
      endif
    enddo

    print *, "Birth or Death process of walkers have been done at istep =", istep
    print *, "Number of Birth =", n_birth
    print *, "Number of Death =", n_die
   
    if ((icopy-1) .ne. nwalker) then
      print *, "Current copies of walkers are", icopy-1
      STOP 'New number of walkers might have some problems. Stop to check.'
    endif

    x_store = 0.0d0
    deallocate(x_store)
    allocate(x_store(nwalker,ndim))
    x_store = x

    allocate(v_vec(nwalker))
    v_vec = 0.0d0
    do iwalker = 1, nwalker
     do idim1 = 1, ndim
      v_vec(iwalker) = v_vec(iwalker) + v(x(iwalker,idim1))
     enddo
    enddo

    !Update E_trial
    !E_trial = E0 + alpha*log((nwalker_aim*1.0d0/nwalker))

    E_trial = sum(v_vec)/nwalker - alpha*(nwalker-nwalker_aim)*1.0d0/(nwalker_aim)
    print *, "Current Trial Energy is", E_trial
    walkerV(istep) = E_trial

    write(77,*), istep, nwalker, E_trial
    deallocate(s)
    deallocate(q)
    deallocate(v_vec)

   enddo

   deallocate(x)
   deallocate(x_store)
  

   call CPU_TIME(time_end)

   !--------The results--------
   print *, "The mc_steps is", mc_steps
   print *, "Average number of walkers are", sum(nw)/mc_steps
   print *, "Average Trial Energy is", sum(walkerV)/mc_steps
   print *, "The used time is (in seconds)", time_end-time_begin
   !--------The results--------  
   

 contains
   subroutine gauss_rand(randx)

     implicit none
     integer, save :: flag = 1
     real(8), save :: ss
     real(8) :: x, y, sq, fac, v1, v2
     real(8), intent(out) :: randx
     if (flag == 1) then
      do
       call random_number(x)
       call random_number(y)
        v1 = 2.0*x-1.0
        v2 = 2.0*y-1.0
        sq = v1*v1+v2*v2
       if (sq <= 1.0d0) exit
      enddo
      fac = sqrt(-2.0d0*log(sq)/sq)
      ss = v1*fac
      randx = v2*fac
      flag = 0
     else
      flag = 1
      randx = ss
     endif
   end subroutine

   function v(x)
     real(8) :: x, v
     v = 0.5d0*x*x
   end function

   function dvdx(x)
     real(8) :: dvdx, x
     dvdx = x
   end function

   function f(x,alpha1)
     real(8) :: x, f, alpha1
     f = -4.0d0*x*alpha1
   end function

   function el(x,alpha1)
     real(8) :: el, x, alpha1
     el = alpha1 + x*x*(0.5d0-2.0d0*alpha1*alpha1)
   end function
end program



