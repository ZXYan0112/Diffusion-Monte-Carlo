program dmc_ho_gf

   implicit none

   integer, parameter :: num = 3000, ndim = 1
   integer, parameter :: nwalker_aim = 3000
   integer, parameter :: mc_steps = 20000
   real(8), parameter :: x_max = 1.0d0
   real(8), parameter :: alpha = 0.1d0
   real(8), parameter :: gamma1 = 0.5d0
   real(8), parameter :: tau = 0.05d0
   real(8), parameter :: E0 = 0.5d0
  
   !---Parameters of Guide Function---!
   real(8), parameter :: alpha_t = 0.45d0
   real(8) :: x_old(ndim), delta_new(ndim), delta_old(ndim)
   real(8) :: Tran_new, Tran_old, rho_new, rho_old, acc_rate, vtot_old

   !---Use for plotting wave function---!
   real(8), parameter :: p_range = 5.0d0
   integer, parameter :: p_bin = 500
   integer :: p_dist(p_bin) = 0
   real(8) :: dist = 0.0d0
   integer :: ibin
 
   integer :: nw(mc_steps) = 0
   real(8) :: walkerV(mc_steps) = 0.0d0
   integer :: nwalker, nwalker_old, n_die, n_birth
   integer :: istep, iwalker, icopy, is, idim1
   real(8) :: time_begin, time_end
   real(8) :: r, el1, el2, tmp, E_trial, tmpg(ndim), vtot, force(ndim), force_new(ndim)
   real(8), allocatable :: x(:,:), x_store(:,:)
   real(8), allocatable :: q(:), v_vec(:)
   integer, allocatable :: s(:)   

   call CPU_TIME(time_begin)
   call random_seed ()

   open (77, File = 'SHO_DMC_GF_results.log')
   write(77,*), "System Parameters:"
   write(77,*), "Aiming Numbers of walker = ", nwalker_aim
   write(77,*), "Initial Numbers of walker = ", num
   write(77,*), "Diffusion parameter = ", gamma1
   write(77,*), "Time step = ", tau
   write(77,*), "Gaussed Trial Energy =", E0
   write(77,*), "Trial Energy adjustable parameter = ", alpha
   write(77,*), "Guide Function's adjustable parameter = ", alpha_t
   
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
      x_old = x(iwalker,:)
      do idim1 = 1, ndim
       call gauss_rand(tmpg(idim1))
       force(idim1) = f(x(iwalker,idim1),alpha_t)
      enddo
      x(iwalker,:) = x(iwalker,:) + tmpg(:)*sqrt(2.0d0*gamma1*tau) + force(:)*tau*0.5d0 

      do idim1 = 1, ndim
       force_new(idim1) = f(x(iwalker,idim1),alpha_t)
      enddo

      vtot = 0.0d0
      do idim1 = 1, ndim
          vtot = vtot + v(x(iwalker,idim1))
      enddo

      vtot_old = 0.0d0
      do idim1 = 1, ndim
          vtot_old = vtot_old + v(x_old(idim1))
      enddo

      !--- Accept the move with a probability of min{1, (Tran_new*rho_new)/(Tran_old*rho_old)} ---!
      !delta_new(:) = x_old(:)-x(iwalker,:)-force_new(:)*tau*0.5d0
      !Tran_new = exp(-dot_product(delta_new, delta_new)/(2.0d0*tau))
      !delta_old(:) = x(iwalker,:)-x_old(:)-force(:)*tau*0.5d0
      !Tran_old = exp(-dot_product(delta_old, delta_old)/(2.0d0*tau))
      !rho_new = exp(-alpha_t*dot_product(x(iwalker,:),x(iwalker,:)))**2.0d0
      !rho_old = exp(-alpha_t*dot_product(x_old,x_old))**2.0d0
      !acc_rate = rho_new*Tran_new/(rho_old*Tran_old)


      acc_rate = 0.0d0
      do idim1 = 1, ndim
        acc_rate = acc_rate + 0.5d0*(force(idim1)+force_new(idim1))*&
                            &(0.25d0*tau*(force(idim1)-force_new(idim1))-x(iwalker,idim1)+x_old(idim1))
      enddo

      acc_rate = acc_rate - 2.0d0*alpha_t*(dot_product(x(iwalker,:),x(iwalker,:)) - dot_product(x_old,x_old))
      acc_rate = exp(acc_rate)

      !print *, "Accepted ratio =", acc_rate

      call random_number(tmp)
      if (acc_rate < tmp) then
        x(iwalker,:) = x_old
        s(iwalker) = 1
      else

       !--- Evaluate the q and s = q + r ---!
       call random_number(tmp)
       q(iwalker) = exp(-tau*(0.5d0*(el(dot_product(x(iwalker,:),x(iwalker,:)),alpha_t,ndim)+el(dot_product(x_old,x_old),alpha_t,ndim))-E_trial))
       !print *, "E_local(R') =", el(dot_product(x(iwalker,:),x(iwalker,:)),alpha_t,ndim)
       !print *, "E_local(R) =", el(dot_product(x_old,x_old),alpha_t,ndim)


       if (q(iwalker) < tmp) then
        s(iwalker) = 0  
       else
        if (q(iwalker)-1.0d0 > tmp) then
          s(iwalker) = 2
        else
          s(iwalker) = 1
        endif
       endif
   
       !---Produce more walkers---!
       s(iwalker) = int(q(iwalker)+tmp)
      endif
      !print *, "(iwalker,s,q)=", iwalker, s(iwalker), q(iwalker)
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
    E_trial = E_trial + alpha*log((nwalker_aim*1.0d0/nwalker))
    !E_trial = sum(v_vec)/nwalker - alpha*(nwalker-nwalker_aim)*1.0d0/(nwalker_aim)
   
    print *, "Current Trial Energy is", E_trial
    walkerV(istep) = E_trial

    write(77,*), istep, nwalker, E_trial
    deallocate(s)
    deallocate(q)
    deallocate(v_vec)

   enddo

   p_dist(:) = 0
   do iwalker = 1, nwalker
     dist = sqrt(dot_product(x(iwalker,:),x(iwalker,:)))

     ibin = int(dist/(p_range/p_bin*1.0)) + 1
     if (ibin > p_bin) ibin = p_bin

     p_dist(ibin) = p_dist(ibin) + 1
   enddo

   deallocate(x)
   deallocate(x_store)
 
   close(77) 

   open (88, File = 'SHO_DMC_wavefunction.log')
   do ibin = 1, p_bin
     write(88,*), ibin*(p_range/p_bin), p_dist(ibin)
   enddo
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

   function el(x2,alpha1,ndim)
     real(8) :: el, x2, alpha1
     integer :: ndim
     el = ndim*alpha1 + x2*(0.5d0-2.0d0*alpha1*alpha1)
   end function
end program



