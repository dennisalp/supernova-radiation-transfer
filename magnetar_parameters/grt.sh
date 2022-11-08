#!/bin/bash

function print_script {
    echo "c     h5fc -fopenmp -march=native -O3 grt.f; time ./a.out
c     
c     dalp the version i have tried to modify
c     dalp profiling on mac: 'instruments -t \"Time Profiler\" a.out'
c     dalp then just launch from finder gui (Instruments.app)

      module constants
      implicit none

!     I/O
! Remember to change lengths of CHARACTER, grid size, t_mod, t_min, t_max, sim_gam
      CHARACTER(LEN=${buffer_size}), PARAMETER :: filename =
     &     \"${model_path}\"
      real*8, parameter :: t_mod=${t_mod} ! seconds
      integer, parameter :: n_rad = ${n_rad}, n_tht = ${n_tht}, n_phi = ${n_phi}
      
      real*8, parameter :: tmin=${tmin}, tmax=${tmax} ! days
      real*8, parameter :: sim_gam=${sim_gam}
      integer, parameter :: DBG=2, n_pho=${n_pho}

!     Model
      integer, parameter :: n_rad1=n_rad+1,n_tht1=n_tht+1,n_phi1=n_phi+1
      integer, parameter :: n_tot = n_rad*n_tht*n_phi, n_abu = 19
      real*8, DIMENSION(:,:,:,:), ALLOCATABLE :: abu
      real*8 rad(n_rad1), tht(n_tht1), phi(n_phi1),
     &     den(n_rad, n_tht, n_phi), dep(n_rad, n_tht, n_phi),
     &     vol(n_rad, n_tht, n_phi), num_cum(n_tot), dif_tht, dif_phi,
     &     ne_tot(n_rad, n_tht, n_phi), vex(n_rad, n_tht, n_phi)

!     Computational
      integer, parameter :: n_ene=301
      
      real*8, parameter :: e_min=${e_min}, e_max=${e_max}
      real*8, parameter :: RAM_LIM=10. ! GB
      real*8, parameter :: TOL=1.d-2 ! (fraction of radius)
      real*8 ene(n_ene), sig(n_abu, n_ene)

!     Physics
      real*8, parameter :: pi=3.14159265d0, cc=2.99792458d10
      real*8, parameter :: m_e=9.10938356d-28, qq=1.60217662d-12
      real*8, parameter :: uu=1.66053904d-24, sig_th=6.652461528d-25
      integer, parameter :: doppler=1
      real*8, parameter :: t_e=1.d0
      real*8 aw(n_abu)
      data aw/1.00782503207d0, 4.00260325415d0, 12.d0, 15.99491461956d0,
     &     19.99244017540d0, 23.9850417d0, 27.9769265325d0, 31.972071d0,
     &     35.967545106d0, 39.96259098d0, 43.955481754d0,43.959402752d0,
     &     43.959690069d0, 47.954032d0, 51.948114d0, 55.934937475d0,
     &     55.939839278d0, 55.942132022d0, 55.934937475d0/
      integer nz(n_abu)
      data nz/1,2,6,8,10,12,14,16,18,20,20,21,22,24,26,26,27,28,26/

      end module constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      use constants
      use omp_lib
      implicit none
      real*8 coord(3), om_i(3), v_i(3), om_f(3), v_f(3), omabs, vabs
      real*8 dl, dl_tot, radius, rn1, rn2, rn1s, rn2s, rn, rn3, chi
      real*8 tau_phabs, tau_lim, x_ph, cum_tau_phabs, cum_tau_scatt
      real*8 opac_phabs, opac_scatt, kn_cross, pp_cross, x_i, phi_f
      real*8 opac_pairp, cum_tau_pairp, opac_tot, cum_tau_tot
      real*8 e_ph_i_ev, e_ph_i, e_ph_ev, e_ph_f, delta_e
      real*8 init_cell, init_step, init_rad, init_tht, init_phi
      integer below, above, mid, itmp, iabs, iesc, ipai
      real*8 rtmp, beta_i, gamma_i, beta_f, gamma_f, mu_i, p_e, mu_f
      real*8 mu_ph, r_xf_xi, rho, overstep, tt, ttmp, ti, tmp
      integer i,ii,jj,kk,nn,ns,np,ene_i,rad_i,tht_i,phi_i
      integer ns_tot
      integer bin_sea_idx, get_ene_idx, get_rad_idx, get_tht_idx,
     &     get_phi_idx
      integer ktot, n_abs, n_esc
      integer,parameter :: seed = 84456
      integer start, count_rate, count_max
      integer thread
      character*18 :: f_dat, f_abs, f_xyz

      call init()
      call get_mod()
      call get_sig()
      call srand(seed)
      call system_clock(start, count_rate, count_max)
      write(6,2349) 'Starting'
 2349 format(A)
      
!\$omp parallel num_threads(3) default(firstprivate) shared(abu, 
!\$omp& rad, tht, phi, den, vol, ne_tot)
      thread = omp_get_thread_num()
      write(f_dat,'(\"../out/tmp_\",I3.3,\".esc\")')thread
      open (unit = 1024+thread, file = f_dat)
!abs      write(f_abs,'(\"../out/tmp_\",I3.3,\".abs\")')thread
!abs      open (unit = 4855+thread, file = f_abs)
!xyz      write(f_xyz,'(\"../out/tmp_\",I3.3,\".xyz\")')thread
!xyz      open (unit = 4895+thread, file = f_xyz)

      n_esc=0
      n_abs=0
      ns_tot=0

c     dalp do the scattering
!\$omp do
      do nn=1,n_pho
         if (thread.eq.0.and.mod(nn,100000).eq.10000) then
            call timing(start,(nn-1)*omp_get_num_threads()/
     &           real(n_pho,8))
         endif
         
!     seed energy and position
         call ini_pow(e_ph_i_ev)
         rad_i=1
         tmp=acos(2*rand()-1)/pi
         tht_i=int(tmp*n_tht)+1
         phi_i=int(rand()*n_phi)+1
         
!     get coordinates
         init_rad = (rad(rad_i)+rad(rad_i+1))/2.
         init_tht = (tht(tht_i)+tht(tht_i+1))/2.
         init_phi = (phi(phi_i)+phi(phi_i+1))/2.
         coord(1)=init_rad*cos(init_phi)*sin(init_tht)
         coord(2)=init_rad*sin(init_phi)*sin(init_tht)
         coord(3)=init_rad*cos(init_tht)

c     photon initialization
         e_ph_ev=e_ph_i_ev
         call ini_pho(e_ph_ev, om_i, rad_i, tht_i, phi_i, coord,
     &        tt, nn)
         ti = tt
         coord = coord*tt/t_mod
         e_ph_i=e_ph_ev*qq
         
         iesc=0
         iabs=0
         ipai=0
         dl_tot=0.
c     scattering loop
         do ns=1,9999
c     move the photon
 5934       x_ph=2*e_ph_i/(m_e*cc**2)
            ene_i=get_ene_idx(e_ph_i/qq)
            e_ph_ev=e_ph_i/qq
c     integrate optical depth along photon path
            tau_lim=-log(rand())
c     tau_lim=1000. ! dbg
            cum_tau_scatt=0.
            cum_tau_phabs=0.
            cum_tau_pairp=0.
            
            do np=1,9999
               opac_scatt=0.
               opac_phabs=0.
               opac_pairp=0.
               ttmp = (t_mod/tt)**3
               do i=1,n_abu
                  rtmp=log10(sig(i,ene_i))+
     &                 (log10(e_ph_ev)-log10(ene(ene_i))) *
     &                 (log10(sig(i,ene_i+1))-log10(sig(i,ene_i))) /
     &                 (log10(ene(ene_i+1))-log10(ene(ene_i)))
                  rtmp=10.**rtmp
                  opac_phabs=opac_phabs+abu(rad_i,tht_i,phi_i,i)*rtmp
                  if (DBG.gt.2) print *,'opac',sig(i, ene_i)/rtmp
                  opac_pairp=opac_pairp+
     &                 abu(rad_i,tht_i,phi_i,i)*nz(i)**2 !swartz95 eq(1)
               enddo
               
               opac_pairp=pp_cross(e_ph_ev)*opac_pairp*ttmp
               opac_scatt=ne_tot(rad_i, tht_i, phi_i)*ttmp*
     &              kn_cross(x_ph)
               opac_phabs=opac_phabs*ttmp
               dl=TOL*sqrt(sum(coord**2))
               cum_tau_scatt=cum_tau_scatt+opac_scatt*dl
               cum_tau_phabs=cum_tau_phabs+opac_phabs*dl
               cum_tau_pairp=cum_tau_pairp+opac_pairp*dl
c     print *, dl, opac_phabs, opac_scatt
c     print *, ne_tot(rad_i, tht_i, phi_i), ttmp
c     print *, kn_cross(x_ph)
               
               cum_tau_tot = cum_tau_scatt+cum_tau_phabs+cum_tau_pairp
               opac_tot = opac_phabs+opac_scatt+opac_pairp
               if (cum_tau_tot.gt.tau_lim) then
                  overstep = 1-(cum_tau_tot-tau_lim)/(opac_tot*dl)
                  call make_step(coord, om_i, overstep*dl, tt)
                  radius=sqrt(sum(coord**2))
                  tmp = rand()
!xyz                  write(4895+thread, '(4Es20.10)') coord, e_ph_ev
                  if (radius.gt.rad(n_rad1)*tt/t_mod) then
!xyz                     write(4895+thread, '(A)') 'ESC'
                     iesc=1
                     goto 99
                  else if (tmp.lt.cum_tau_phabs/cum_tau_tot) then
!xyz                     write(4895+thread, '(A)') 'ABS'
cccp                      dep(rad_i, tht_i, phi_i) = e_ph_i +
cccp      &                    dep(rad_i, tht_i, phi_i)
                     iabs=1
                     goto 99
                  else if (tmp.lt.(cum_tau_phabs+cum_tau_pairp)/
     &                    cum_tau_tot) then
                     ipai=1
                     e_ph_ev=m_e*cc**2/qq
                     call ini_ppp(e_ph_ev, om_i,
     &                    rad_i, tht_i, phi_i, coord)
                     e_ph_i=e_ph_ev*qq
                     goto 5934
                  endif
                  goto 7171
               endif

               call make_step(coord, om_i, dl, tt)
               radius=sqrt(sum(coord**2))
               if (radius.gt.rad(n_rad1)*tt/t_mod) then
!xyz                  write(4895+thread, '(4Es20.10)') coord, e_ph_ev
!xyz                  write(4895+thread, '(A)') 'ESC'
                  iesc=1
                  goto 99
               endif
               
               dl_tot=dl_tot+dl
               rad_i=get_rad_idx(radius, tt)
               tht_i=get_tht_idx(acos(coord(3)/radius))
               phi_i=get_phi_idx(atan2(coord(2), coord(1)))
               if (DBG.gt.2) print *,'path',np,dl,
     &              t_mod/tt*radius/rad(n_rad1),
     &              cum_tau_scatt, opac_scatt, tau_lim
            enddo
            print *, 'ERROR PATH INTEGRATION FAILED'
            stop
 7171       continue

            beta_i = rad(rad_i)/(t_mod*cc)
            gamma_i = 1/sqrt(1-beta_i**2)
            p_e=m_e*cc*beta_i*gamma_i
            v_i = coord/sqrt(sum(coord**2))
            
c     angle between photon direction and electron velocity
            mu_i=0.
            do i=1,3
               mu_i=v_i(i)*om_i(i)+mu_i
            enddo
            
c     dimensionless photon energy in electron rest frame
            x_i=2.*gamma_i*e_ph_i*(1.-mu_i*beta_i)/(m_e*cc**2)
            
c     calculate direction and energy of photon after scattering
 55         rn1=rand()
            mu_f=(beta_i+2.*rn1-1.)/(1.+beta_i*(2.*rn1-1.))
            rn2=rand()
            phi_f=2.*pi*rn2
c     new photon direction
            rho=sqrt(v_i(1)**2+v_i(2)**2)
c     error in Santana
            om_f(1)=mu_f*v_i(1)+sqrt(1.-mu_f**2)*(v_i(2)*cos(phi_f)+
     &           v_i(1)*v_i(3)*sin(phi_f))/rho
            om_f(2)=mu_f*v_i(2)+sqrt(1.-mu_f**2)*(-v_i(1)*cos(phi_f)+
     &           v_i(2)*v_i(3)*sin(phi_f))/rho
c     note that santana has the wrong sign here!
            om_f(3)=mu_f*v_i(3)-sqrt(1.-mu_f**2)*rho*sin(phi_f)
            
            omabs=sqrt(om_f(1)**2+om_f(2)**2+om_f(3)**2)
            vabs=sqrt(v_i(1)**2+v_i(2)**2+v_i(3)**2)
            
c     angle between initial and final photon direction
            mu_ph=0.
            do i=1,3
               mu_ph=mu_ph+om_i(i)*om_f(i)
            enddo
            
            r_xf_xi = 1./(1.+e_ph_i*(1.-mu_ph)/
     &           (gamma_i*m_e*cc**2*(1.-mu_f*beta_i)))
            if(r_xf_xi<0.) then
               write(6,9238)e_ph_i,mu_ph,gamma_i,mu_f,beta_i,r_xf_xi
 9238          format('ERROR r_xf_xi < 0 ',1pe12.3,10e12.3)
               stop
            endif

 34         rn3=rand()
            chi=1./r_xf_xi + r_xf_xi + (4./x_i)*(1.-1./r_xf_xi) +
     &           (4./x_i**2)*(1.-1./r_xf_xi)**2
            if(2.*rn3 < r_xf_xi**2*chi) then
               goto 32
            else
               goto 55
            endif
 32         e_ph_f=r_xf_xi*x_i*m_e*cc**2/
     &           (2.*gamma_i*(1.-mu_f*beta_i))
            
c     calculate new electron energy and direction
            gamma_f = (e_ph_i-e_ph_f+m_e*cc**2*gamma_i)/(m_e*cc**2)
            beta_f=1.-1./gamma_f**2
            
c     deposition to electrons
            delta_e=(gamma_f-gamma_i)*m_e*cc**2
cccp             dep(rad_i, tht_i, phi_i)=dep(rad_i, tht_i, phi_i)+delta_e
            
c     new direction
            v_f(1)=(e_ph_i*om_i(1)/cc+m_e*cc*beta_i*gamma_i*v_i(1) -
     &           e_ph_f*om_f(1)/cc)/(m_e*cc*beta_f)
            v_f(2)=(e_ph_i*om_i(2)/cc+m_e*cc*beta_i*gamma_i*v_i(2) -
     &           e_ph_f*om_f(2)/cc)/(m_e*cc*beta_f)
            v_f(3)=(e_ph_i*om_i(3)/cc+m_e*cc*beta_i*gamma_i*v_i(3) -
     &           e_ph_f*om_f(3)/cc)/(m_e*cc*beta_f)
            
            do i=1,3
               om_i(i)=om_f(i)
            enddo
            e_ph_i=e_ph_f
         enddo
         print *, 'ERROR SCATTERING FAILED'
         stop
 99      continue

         ns=ns-1
c     dalp tally the photons
         if(iesc.eq.1) then
            n_esc=n_esc+1
            if(ipai.eq.1) then
               ns = 9999
               call print_event(e_ph_i_ev, e_ph_ev,
     &              om_i, ns, coord, tt, ti, thread)
            endif
            call print_event(e_ph_i_ev, e_ph_ev,
     &           om_i, ns, coord, tt, ti, thread)
            
         elseif(iabs.eq.1) then
            n_abs=n_abs+1
!abs            write(4855+thread, '(4Es20.13)') e_ph_i_ev,
!abs     &           e_ph_ev,ti,tt-sum(coord*om_i)/cc
         endif
         ns_tot=ns_tot+ns
         e_ph_ev=e_ph_f/qq
         
c     dalp this closes the scattering loop
      enddo
!\$omp end do
      write(6,*)'n_pho, n_esc, n_abs ',n_pho, n_esc, n_abs
!\$omp end parallel
      end

      real*8 function kn_cross(x)
      implicit none
      real*8 x,sig
      if(x<=0.5) then
         sig=1./3.+0.141*x-0.12*x**2+(1.+0.5*x)*(1.+x)**(-2)
      elseif(x>0.5.and.x<=3.5) then
         sig=(log(1.+x)+0.06)/x
      elseif(x>=3.5) then
         sig=(log(1.+x)+0.5-1./(2.+0.076*x))/x
      endif
      kn_cross=sig*4.989346146d-25
      return
      end

      real*8 function pp_cross(ee)
      use constants
      implicit none
      real*8 ee, e2, hlp
      e2 = ee/1.d6
      hlp = 2*m_e*cc**2/qq/1.d6
      if(e2.le.hlp) then
         pp_cross=0.d0
      elseif(e2.le.1.5) then
         pp_cross=0.10063*(e2-hlp)
      else
         pp_cross=0.0481+0.301*(e2-1.5)
      endif
      pp_cross=pp_cross*1.d-27
      return
      end
     
      subroutine make_step(coord, om_i, dl, tt)
      use constants
      implicit none
      real*8 coord(3), om_i(3), dl, tt
      coord = coord+om_i*dl
      tt = tt + dl/cc
      end

      subroutine print_event(e_ph_i_ev, e_ph_ev,
     &     om_i, ns, coord, tt, ti, thread)
      use constants
      implicit none
      real*8 e_ph_i_ev, e_ph_ev, om_i(3), coord(3), tt, t2, ti
      real*8 orth(3) ! orthogonal to line-of-sight
      integer ns, thread

c     project tt into observer's time
      t2=tt-sum(coord*om_i)/cc

c     print to file
      write(1024+thread, '(Es19.13, Es20.13, 3F12.7, I5, 2Es20.13)')
     &     e_ph_i_ev,e_ph_ev,om_i,ns,t2,ti
      end

      subroutine ini_ppp(e_ph_ev, om_i, rad_i, tht_i, phi_i, coord)
      use constants
      implicit none
      real*8 e_ph_ev, e_tmp, om_i(3), om_tmp(3), coord(3), rn1, rn2
      real*8 gamma, beta(3), babs
      integer rad_i, tht_i, phi_i

c     input energy and random direction in electron frame
      e_tmp = e_ph_ev
      rn1=rand()
      rn2=rand()
      om_tmp(3)=2*rn1-1.
      om_tmp(2)=sqrt(1.-om_tmp(3)**2)*sin(2*pi*rn2)
      om_tmp(1)=sqrt(1.-om_tmp(3)**2)*cos(2*pi*rn2)
      om_tmp = e_ph_ev*om_tmp

c     electron velocity in observer frame
      beta = rad(rad_i)/(t_mod*cc)*coord/sqrt(sum(coord**2))
      babs = sqrt(sum(beta**2))
      gamma = 1/sqrt(1-babs**2)

c     reverse boosts from electron to observer frame
      e_ph_ev = gamma*e_tmp-
     &     gamma*beta(1)*om_tmp(1)-
     &     gamma*beta(2)*om_tmp(2)-
     &     gamma*beta(3)*om_tmp(3)
      om_i(1) = -gamma*beta(1)*e_tmp+
     &     (1+(gamma-1)*beta(1)**2/babs**2)*om_tmp(1)+
     &     ((gamma-1)*beta(1)*beta(2)/babs**2)*om_tmp(2)+
     &     ((gamma-1)*beta(1)*beta(3)/babs**2)*om_tmp(3)
      om_i(2) = -gamma*beta(2)*e_tmp+
     &     ((gamma-1)*beta(1)*beta(2)/babs**2)*om_tmp(1)+
     &     (1+(gamma-1)*beta(2)**2/babs**2)*om_tmp(2)+
     &     ((gamma-1)*beta(3)*beta(2)/babs**2)*om_tmp(3)
      om_i(3) = -gamma*beta(3)*e_tmp+
     &     ((gamma-1)*beta(1)*beta(3)/babs**2)*om_tmp(1)+
     &     ((gamma-1)*beta(2)*beta(3)/babs**2)*om_tmp(2)+
     &     (1+(gamma-1)*beta(3)**2/babs**2)*om_tmp(3)
      om_i = om_i/sqrt(sum(om_i**2))
      end
      
      subroutine ini_pho(e_ph_ev, om_i, rad_i, tht_i, phi_i, coord, tt,
     &     nn)
      use constants
      implicit none
      real*8 e_ph_ev, e_tmp, om_i(3), om_tmp(3), coord(3), rn1, rn2
      real*8 gamma, beta(3), babs, tt
      integer rad_i, tht_i, phi_i, nn

c     time of emission
      tt=(tmax-tmin)*rand()+tmin
      tt=tt*24.d0*60.d0*60.d0

c     input energy and random direction in electron frame
      e_tmp = e_ph_ev
      rn1=rand()
      rn2=rand()
      om_tmp(3)=2*rn1-1.
      om_tmp(2)=sqrt(1.-om_tmp(3)**2)*sin(2*pi*rn2)
      om_tmp(1)=sqrt(1.-om_tmp(3)**2)*cos(2*pi*rn2)
c      om_tmp(3)=0.; om_tmp(2)=0.6; om_tmp(1)=0.8 ! dbg
      if (doppler.eq.0) then
         om_i = om_tmp
         goto 7391
      endif
      om_tmp = e_ph_ev*om_tmp

c     electron velocity in observer frame
      beta = rad(rad_i)/(t_mod*cc)*coord/sqrt(sum(coord**2))
      babs = sqrt(sum(beta**2))
      gamma = 1/sqrt(1-babs**2)

c     reverse boosts from electron to observer frame
      e_ph_ev = gamma*e_tmp-
     &     gamma*beta(1)*om_tmp(1)-
     &     gamma*beta(2)*om_tmp(2)-
     &     gamma*beta(3)*om_tmp(3)
      om_i(1) = -gamma*beta(1)*e_tmp+
     &     (1+(gamma-1)*beta(1)**2/babs**2)*om_tmp(1)+
     &     ((gamma-1)*beta(1)*beta(2)/babs**2)*om_tmp(2)+
     &     ((gamma-1)*beta(1)*beta(3)/babs**2)*om_tmp(3)
      om_i(2) = -gamma*beta(2)*e_tmp+
     &     ((gamma-1)*beta(1)*beta(2)/babs**2)*om_tmp(1)+
     &     (1+(gamma-1)*beta(2)**2/babs**2)*om_tmp(2)+
     &     ((gamma-1)*beta(3)*beta(2)/babs**2)*om_tmp(3)
      om_i(3) = -gamma*beta(3)*e_tmp+
     &     ((gamma-1)*beta(1)*beta(3)/babs**2)*om_tmp(1)+
     &     ((gamma-1)*beta(2)*beta(3)/babs**2)*om_tmp(2)+
     &     (1+(gamma-1)*beta(3)**2/babs**2)*om_tmp(3)
      om_i = om_i/sqrt(sum(om_i**2))
 7391 end

      integer function bin_sea_idx(val, arr, len)
      use constants
      implicit none
      integer len, below, above, mid, ii
      real*8 val, arr(len)

      below=1
      above=len
      do ii=1,99999
         mid=(below+above)/2
         if (mid.eq.below) goto 9137
         if (arr(mid).gt.val) then
            above = mid
         else
            below = mid
         endif
      enddo
      print *, 'ERROR BINARY SEARCH FOR INITIAL POSITION FAILED'
      stop

 9137 bin_sea_idx=mid
      if (DBG.gt.2) print *, 'bin_sea_idx', bin_sea_idx
      return
      end

c     this assumes ene_in in units of ev
      integer function get_ene_idx(ene_in)
      use constants
      implicit none
      real*8 de, ene_in
      de=log10(e_max/e_min)/real(n_ene-1)

      get_ene_idx = int(log10(ene_in/e_min)/de+1)
      get_ene_idx = max(1, get_ene_idx)
      get_ene_idx = min(n_ene, get_ene_idx)
      if (DBG.gt.2) print *, 'get_ene_idx', get_ene_idx,
     &     ene(get_ene_idx), ene_in, ene(get_ene_idx)/ene_in
      return
      end

      integer function get_rad_idx(rad_in, tt)
      use constants
      implicit none
      integer bin_sea_idx
      real*8 rad_in, tt
      get_rad_idx = bin_sea_idx(rad_in, tt/t_mod*rad, n_rad1)
      if (DBG.gt.2) print *, 'get_rad_idx', get_rad_idx,
     &     rad(get_rad_idx), rad_in, rad(get_rad_idx)/rad_in
      return
      end
      integer function get_tht_idx(tht_in)
      use constants
      implicit none
      real*8 tht_in
      get_tht_idx = int(tht_in/dif_tht)+1
      if (DBG.gt.2) print *, 'get_tht_idx', get_tht_idx,
     &     tht(get_tht_idx), tht_in, tht(get_tht_idx)-tht_in
      return
      end
      integer function get_phi_idx(phi_in)
      use constants
      implicit none
      real*8 phi_in
      get_phi_idx = int((phi_in+pi)/dif_phi)+1
      if (DBG.gt.2) print *, 'get_phi_idx', get_phi_idx,
     &     phi(get_phi_idx), phi_in, phi(get_phi_idx)-phi_in
      return
      end

      subroutine init()
      use constants
      implicit none
      integer ii,jj,kk,ll,ne,idx
      real*8 de, ram

      ram = real(n_rad*n_tht*n_phi,8)*(n_abu+5)*8/1e9
      if (ram.gt.RAM_LIM) then
         write(6,*) 'ERROR, RAM LIMIT EXCEEDED'
         stop
      endif
      write(6, 2981) 'Allocating approximately ', ram, ' GB of memory'
 2981 format(A, F3.1, A)

      de=log10(e_max/e_min)/real(n_ene-1)
      do ne=1, n_ene
         ene(ne) = 10.**(de*(ne-1))*e_min
      enddo

      ALLOCATE(abu(n_rad, n_tht, n_phi, n_abu))
      do ii = 1, n_rad
         rad(ii)  = 0.
         do jj = 1, n_tht
            if (ii.eq.1) tht(jj)  = 0.
            do kk = 1, n_phi
               if (ii.eq.1) phi(kk)  = 0.
               den(ii,jj,kk)  = 0.
               dep(ii,jj,kk)  = 0.
               vol(ii,jj,kk)  = 0.
               ne_tot(ii,jj,kk) = 0.
               idx = (ii-1)*(n_tht*n_phi)+(jj-1)*n_phi+kk
               do ll = 1, n_abu
                  abu(ii,jj,kk,ll) = 0.
               end do
            end do
         end do
      end do
      end

      subroutine get_mod()
      use constants
      USE HDF5
      IMPLICIT NONE
!     File
      character(len=4),dimension(19),parameter :: ele = ['p   ',
     &     'he4 ', 'c12 ', 'o16 ', 'ne20', 'mg24', 'si28', 's32 ',
     &     'ar36', 'ca40', 'ca44', 'sc44', 'ti44', 'cr48', 'fe52',
     &     'fe56', 'co56', 'ni56', 'x56 ']

!     HDF5
      INTEGER(HID_T) :: file_id, temp_id
      INTEGER error, ii, jj, kk, ll, idx
      INTEGER(HSIZE_T), DIMENSION(3) :: dat_dim
      data dat_dim/n_rad,n_tht,n_phi/
      INTEGER(HSIZE_T), DIMENSION(1) :: rad_dim, tht_dim, phi_dim
      data rad_dim/n_rad1/,tht_dim/n_tht1/,phi_dim/n_phi1/
      
!     Initialize FORTRAN interface and open an existing file.
      CALL h5open_f(error)
      CALL h5fopen_f (filename, H5F_ACC_RDONLY_F, file_id, error)

!     Open, read, and close an existing dataset
      write(6,4220) 'Loading the model'
 4220 format(A)
      CALL h5dopen_f(file_id, 'den', temp_id, error)
      CALL h5dread_f(temp_id, H5T_NATIVE_DOUBLE, den, dat_dim, error)
      if (error.ne.0) write(*,*) 'ERROR LOADING DATA', error, temp_id
      CALL h5dclose_f(temp_id, error)

      CALL h5dopen_f(file_id, 'radius', temp_id, error)
      CALL h5dread_f(temp_id, H5T_NATIVE_DOUBLE, rad, rad_dim, error)
      if (error.ne.0) write(*,*) 'ERROR LOADING DATA', error, temp_id
      CALL h5dclose_f(temp_id, error)

      CALL h5dopen_f(file_id, 'theta', temp_id, error)
      CALL h5dread_f(temp_id, H5T_NATIVE_DOUBLE, tht, tht_dim, error)
      if (error.ne.0) write(*,*) 'ERROR LOADING DATA', error, temp_id
      CALL h5dclose_f(temp_id, error)
      
      CALL h5dopen_f(file_id, 'phi', temp_id, error)
      CALL h5dread_f(temp_id, H5T_NATIVE_DOUBLE, phi, phi_dim, error)
      if (error.ne.0) write(*,*) 'ERROR LOADING DATA', error, temp_id
      CALL h5dclose_f(temp_id, error)

      do ii = 1, n_abu
         CALL h5dopen_f(file_id, trim(ele(ii)), temp_id, error)
         CALL h5dread_f(temp_id, H5T_NATIVE_DOUBLE, abu(:,:,:,ii),
     &        dat_dim, error)
         if (error.ne.0) write(*,*) 'ERROR LOADING DATA', error, temp_id
         CALL h5dclose_f(temp_id, error)
      end do
      
!     Close the file and FORTRAN interface.
      CALL h5fclose_f(file_id, error)
      CALL h5close_f(error)
      
!     Convert to number densities, expand homologously, and compute
!     electron number density
      do ii = 1, n_abu
         abu(:,:,:,ii) = den*abu(:,:,:,ii)/(uu*aw(ii))
         ne_tot = ne_tot + abu(:,:,:,ii)*nz(ii)
      end do
      
      dif_tht = tht(2)-tht(1)
      dif_phi = phi(2)-phi(1)

!     dbg print *, (rad(ii), ii = 1, n_rad1)
!     dbg print *, (tht(ii), ii = 1, n_tht1)
!     dbg print *, (phi(ii), ii = 1, n_phi1)
!     dbg      do jj = 10, 10
!     dbg         do kk = 100, 100
!     dbg            print *, (abu(ii,jj,kk,9), ii = 1, n_rad)
!     dbg         end do
!     dbg      end do
      end

      subroutine ini_pow(e_ph_i_ev)
      use constants
      implicit none
      real*8 e_ph_i_ev, tmp, sg1
      
!     Initialize energy as a power law
!     https://stackoverflow.com/questions/918736/random-number-generator-that-produces-a-power-law-distribution
!     https://mathworld.wolfram.com/RandomNumber.html
!     x = [(x1^(n+1) - x0^(n+1))*y + x0^(n+1)]^(1/(n+1))

      sg1 = -sim_gam+1
      if (abs(sg1).lt.1.d-6) then
         e_ph_i_ev = (e_max/e_min)**rand()*e_min
      else
         tmp = (e_max**sg1 - e_min**sg1)*rand()
         tmp = tmp + e_min**sg1
         e_ph_i_ev = tmp**(1/sg1)
      endif
      end
 

      
      subroutine get_sig()
      use constants
      implicit none
      integer ne,na,ns,nel,nsh
      real*8 sig_is
c      open(unit=2299, file='sig_verner96.txt')
      write(6,3322) 'Computing cross sections'
 3322 format(A)
      do ne=1,n_ene
         do na=1,n_abu
            sig(na,ne) = 0.
c     total cross sections by summing over shells
            nel = nz(na)
c     dalp number of shells for given atomic number nel
            if(nel.ge.20) then
               nsh = 7
            elseif(nel.ge.13) then
               nsh = 5
            elseif(nel.ge.11) then
               nsh = 4
            elseif(nel.ge.6) then
               nsh = 3
            elseif(nel.ge.1) then
               nsh = 1
            endif
            do ns = 1,nsh
               call phfit2(nel,nel,ns,ene(ne),sig_is)
               sig(na,ne) = sig(na,ne) + sig_is
            enddo
         enddo
      enddo
c      do na=1,n_abu
c         write(2299,*) sig(na,:)
c      enddo
      end
      
      subroutine timing(start, frac)
      implicit none
      real now, elapsed
      real*8 frac
      integer start, count, count_rate, count_max

      call system_clock(count, count_rate, count_max)
      elapsed = (count-start)/count_rate
      write(6,4105) 'Progress:', 100*frac, '%, time elapsed:', elapsed,
     &     ' s, time left:', elapsed*(1/frac-1), ' s'
 4105 format(A, F5.1, A, F8.1, A, F8.1, A)
      end




      subroutine phfit2(nz,ne,is,e,s)
*** Version 2. March 25, 1996.
*** Written by D. A. Verner, verner@pa.uky.edu
*** Inner-shell ionization energies of some low-ionized species are slightly
*** improved to fit smoothly the experimental inner-shell ionization energies
*** of neutral atoms.
******************************************************************************
*** This subroutine calculates partial photoionization cross sections
*** for all ionization stages of all atoms from H to Zn (Z=30) by use of
*** the following fit parameters:
*** Outer shells of the Opacity Project (OP) elements:
***    Verner, Ferland, Korista, Yakovlev, 1996, ApJ, in press.
*** Inner shells of all elements, and outer shells of the non-OP elements:
***    Verner and Yakovlev, 1995, A&AS, 109, 125
*** Input parameters:  nz - atomic number from 1 to 30 (integer)
***                    ne - number of electrons from 1 to iz (integer)
***                    is - shell number (integer)
***                    e - photon energy, eV
*** Output parameter:  s - photoionization cross section, Mb
*** Shell numbers:
*** 1 - 1s, 2 - 2s, 3 - 2p, 4 - 3s, 5 - 3p, 6 - 3d, 7 - 4s.
*** If a species in the ground state has no electrons on the given shell,
*** the subroutine returns s=0.
******************************************************************************
      implicit real*8(a-h,o-z)
      integer nz,ne,is
      real*8 e,s
      common/l/l(7)
      common/ninn/ninn(30)
      common/ntot/ntot(30)
      common/ph1/ph1(6,30,30,7)
      common/ph2/ph2(7,30,30)

c      write(6,*)' phfit ',nz,ne,is,e,s

      s=0.0
      if(nz.lt.1.or.nz.gt.30)return
      if(ne.lt.1.or.ne.gt.nz)return
      nout=ntot(ne)
      if(nz.eq.ne.and.nz.gt.18)nout=7
      if(nz.eq.(ne+1).and.(nz.eq.20.or.nz.eq.21.or.nz.eq.22.or.nz.
     &eq.25.or.nz.eq.26))nout=7
      if(is.gt.nout)return
      if(e.lt.ph1(1,nz,ne,is))return
      nint=ninn(ne)
      if(nz.eq.15.or.nz.eq.17.or.nz.eq.19.or.
     &(nz.gt.20.and.nz.ne.26))then
         einn=0.0
      else
         if(ne.lt.3)then
            einn=1.0e+30
         else
            einn=ph1(1,nz,ne,nint)
         endif
      endif
      if(is.lt.nout.and.is.gt.nint.and.e.lt.einn)return
      if(is.le.nint.or.e.ge.einn)then
         p1=-ph1(5,nz,ne,is)
         y=e/ph1(2,nz,ne,is)
         q=-0.5*p1-l(is)-5.5
         a=ph1(3,nz,ne,is)*((y-1.0)**2+ph1(6,nz,ne,is)**2)
         b=sqrt(y/ph1(4,nz,ne,is))+1.0
         s=a*y**q*b**p1
      else
         p1=-ph2(4,nz,ne)
         q=-0.5*p1-5.5
         x=e/ph2(1,nz,ne)-ph2(6,nz,ne)
         z=sqrt(x*x+ph2(7,nz,ne)**2)
         a=ph2(2,nz,ne)*((x-1.0)**2+ph2(5,nz,ne)**2)
         b=1.0+sqrt(z/ph2(3,nz,ne))
         s=a*z**q*b**p1
      endif
      s = s*1.e-18

c      write(6,*)' phfit ',nz,ne,is,e,s

      return
      end

      BLOCK DATA BDATA
      implicit real*8(a-h,o-z)
      COMMON/L/L(7)
      COMMON/NINN/NINN(30)
      COMMON/NTOT/NTOT(30)
      COMMON/PH1/PH1(6,30,30,7)
      COMMON/PH2/PH2(7,30,30)
      DATA (L(I),I=1,7) /0,0,1,0,1,2,0/
      DATA (NINN(I),I=1,30) /0,0,1,1,1,1,1,1,1,1,3,3,
     & 3,3,3,3,3,3,5,5,5,5,5,5,5,5,5,5,5,5/
      DATA (NTOT(I),I=1,30) /1,1,2,2,3,3,3,3,3,3,4,4,
     & 5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7/
      DATA (PH1(I, 1, 1, 1),I=1,6) /1.360E+01, 4.298E-01,
TRUNCATED
     1 2.671E+07, 7.923E+00, 2.069E+01, 1.382E+02, 2.481E-01/
      END

      DOUBLE PRECISION FUNCTION EXP1(X)
      implicit none
      REAL*8 x,e1
      IF(X.LT.1.) THEN
            E1=-LOG(X)-.577216+X-.24991*X*X+.0552*X*X*X-.00976
     &            *X*X*X*X
      ELSE
            E1=EXP(-X)*(X*X+2.334733*X+.250621)/((X*X+3.330657*X
     &            +1.681534)*X)
      ENDIF
      EXP1=E1
      RETURN
      END" > grt.f
}

function print_meta {
    echo "id	${ID}
model	${model}
tmin	${tmin/d/e}
tmax	${tmax/d/e}
e_min	${e_min/d/e}
e_max	${e_max/d/e}
sim_gam	${sim_gam/d/e}
n_pho	${n_pho}" > ../out/${output_root}_${ID}.txt
}

function run_script {
    h5fc -fopenmp -march=native -O3 grt.f -o grt.out
    time ./grt.out
    cat ../out/tmp*.esc > ../out/${output_root}_${ID}.esc
    rm ../out/tmp*.esc
}

function run_mod {
    model_path="/Users/silver/dat/sne/${model}.h5"
    buffer_size=${#model_path}
    output_root=$(echo "${model}" | tr '[:upper:]' '[:lower:]')
    
    # ID="r00"
    # tmin=1.d0
    # tmax=1030.d0
    # e_min=1.d4
    # e_max=1.d7
    # sim_gam=1.d0
    # print_script
    # print_meta
    # run_script

    ID="r00"
    tmin=1.d0
    tmax=10300.d0
    e_min=1.d2
    e_max=1.d7
    sim_gam=1.d0
    print_script
    print_meta
    run_script

    # ID="r99"
    # tmin=9700.d0
    # tmax=10300.d0
    # e_min=1.d6
    # e_max=1.d7
    # sim_gam=1.d0
    # print_script
    # print_meta
    # run_script
}

################################################################
# Specific for the model

model="SLSN-I_000"
t_mod=762068.30849
n_rad=1169
n_tht=90
n_phi=180
n_pho=100000000
run_mod

# post processing time python kicks/kicks.py /Users/silver/dat/sne/ M178a
# pp ./run_all.sh m178a+
#python kicks/kicks.py /Users/silver/dat/sne/ M157b+
# pp ./run_all_stripped.sh ic
