#!/bin/bash

function print_script {
    echo "c     h5fc -fopenmp -march=native -O3 src/grt.f; time ./a.out
c     
c     dalp the version i have tried to modify
c     dalp profiling on mac: 'instruments -t \"Time Profiler\" a.out'
c     dalp then just launch from finder gui (Instruments.app)

      module constants
      implicit none

!     I/O
! Remember to change lengths of CHARACTER, grid size, t_mod, t_min, t_max, src_abu
      CHARACTER(LEN=${buffer_size}), PARAMETER :: filename =
     &     \"${model_path}\"
      real*8, parameter :: t_mod=${t_mod} ! seconds
      integer, parameter :: n_rad = ${n_rad}, n_tht = ${n_tht}, n_phi = ${n_phi}
      
      real*8, parameter :: tmin=${tmin}, tmax=${tmax} ! days
      integer, parameter :: DBG=2, src_abu=${src_abu}, n_pho=${n_pho}

!     Model
      integer, parameter :: n_rad1=n_rad+1,n_tht1=n_tht+1,n_phi1=n_phi+1
      integer, parameter :: n_tot = n_rad*n_tht*n_phi, n_abu = 19
      real*8, DIMENSION(:,:,:,:), ALLOCATABLE :: abu
      real*8 rad(n_rad1), tht(n_tht1), phi(n_phi1),
     &     den(n_rad, n_tht, n_phi), dep(n_rad, n_tht, n_phi),
     &     vol(n_rad, n_tht, n_phi), num_cum(n_tot), dif_tht, dif_phi,
     &     ne_tot(n_rad, n_tht, n_phi), vex(n_rad, n_tht, n_phi)

!     Computational
      integer, parameter :: n_ene=501
      
      real*8, parameter :: e_min=1., e_max=4.d6
      real*8, parameter :: RAM_LIM=10. ! GB
      real*8, parameter :: TOL=1.d-2 ! (fraction of radius)
      real*8 ene(n_ene), sig(n_abu, n_ene)

!     Physics
      real*8, parameter :: pi=3.14159265d0, cc=2.99792458d10
      real*8, parameter :: m_e=9.10938356d-28, qq=1.60217662d-12
      real*8, parameter :: uu=1.66053904d-24, sig_th=6.652461528d-25
      integer, parameter :: maxwell=0, doppler=1
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
      character*15 :: f_dat, f_abs, f_xyz

      call init()
      call get_mod()
      call get_sig_hoeflich92()
      call srand(seed)
      call system_clock(start, count_rate, count_max)
      write(6,2349) 'Starting'
 2349 format(A)
      
!\$omp parallel num_threads(3) default(firstprivate) shared(abu, 
!\$omp& rad, tht, phi, den, vol, num_cum, ne_tot, vex)
      thread = omp_get_thread_num()
      write(f_dat,'(\"out/tmp_\",I3.3,\".esc\")')thread
      open (unit = 1024+thread, file = f_dat)
!abs      write(f_abs,'(\"out/tmp_\",I3.3,\".abs\")')thread
!abs      open (unit = 4855+thread, file = f_abs)
!xyz      write(f_xyz,'(\"out/tmp_\",I3.3,\".xyz\")')thread
!xyz      open (unit = 4895+thread, file = f_xyz)

      n_esc=0
      n_abs=0
      ns_tot=0
      
c     dalp do the scattering
!\$omp do
      do nn=1,n_pho
         if(src_abu.eq.0) then
            call ini_44(e_ph_i_ev, nn)
         elseif(src_abu.eq.1) then
            call ini_56ni(e_ph_i_ev, nn)
         elseif(src_abu.eq.2) then
            call ini_56co(e_ph_i_ev, nn)
         elseif(src_abu.eq.3) then
            call ini_57co(e_ph_i_ev, nn)
         endif
         if (thread.eq.0.and.mod(nn,100000).eq.10000) then
            call timing(start,(nn-1)*omp_get_num_threads()/
     &           real(n_pho,8))
         endif
         
         init_cell=rand()*num_cum(n_tot)
c     init_cell=0.30235047009401878*num_cum(n_tot) ! dbg
         itmp=bin_sea_idx(init_cell, num_cum, n_tot)
         
!     get indices
         rad_i=(itmp-1)/(n_tht*n_phi)+1
         itmp=itmp-(rad_i-1)*(n_tht*n_phi)
         tht_i=(itmp-1)/n_phi+1
         itmp=itmp-(tht_i-1)*n_phi
         phi_i=itmp
c     rad_i = 1 ! dbg
c     tht_i = 1 ! dbg
c     phi_i = 1 ! dbg
         
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
c     om_i(3)=0. ! dbg
c     om_i(2)=0. ! dbg
c     om_i(1)=1. ! dbg
         
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
                  if (DBG.gt.3) print *,'opac',sig(i, ene_i)/rtmp
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

            if (maxwell.eq.1) then
c     draw the momentum and energy of the electron for a maxwellian distribution
               call el_en(t_e,beta_i,gamma_i)
               p_e=m_e*cc*beta_i*gamma_i
c     electron velocity direction
 76            rn1=rand()
               rn2=rand()
               rn1s=rn1
               rn2s=rn2
               if(rn1==0..or.rn2==0.) then
                  goto 76
               endif
               v_i(3)=2*rn1-1.
               v_i(2)=sqrt(1.-v_i(3)**2)*sin(2*pi*rn2)
               v_i(1)=sqrt(1.-v_i(3)**2)*cos(2*pi*rn2)
            else
               beta_i = vex(rad_i, tht_i, phi_i)/cc
               gamma_i = 1/sqrt(1-beta_i**2)
               p_e=m_e*cc*beta_i*gamma_i
               v_i = coord/sqrt(sum(coord**2))
            endif
            
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
      beta = vex(rad_i, tht_i, phi_i)/cc*coord/sqrt(sum(coord**2))
      babs = sqrt(sum(beta**2))
      gamma = 1/sqrt(1-babs**2)

c     reverse boosts from electron to observer frame
      e_ph_ev = gamma*e_tmp+
     &     gamma*beta(1)*om_tmp(1)+
     &     gamma*beta(2)*om_tmp(2)+
     &     gamma*beta(3)*om_tmp(3)
      om_i(1) = gamma*beta(1)*e_tmp+
     &     (1+(gamma-1)*beta(1)**2/babs**2)*om_tmp(1)+
     &     ((gamma-1)*beta(1)*beta(2)/babs**2)*om_tmp(2)+
     &     ((gamma-1)*beta(1)*beta(3)/babs**2)*om_tmp(3)
      om_i(2) = gamma*beta(2)*e_tmp+
     &     ((gamma-1)*beta(1)*beta(2)/babs**2)*om_tmp(1)+
     &     (1+(gamma-1)*beta(2)**2/babs**2)*om_tmp(2)+
     &     ((gamma-1)*beta(3)*beta(2)/babs**2)*om_tmp(3)
      om_i(3) = gamma*beta(3)*e_tmp+
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
      beta = vex(rad_i, tht_i, phi_i)/cc*coord/sqrt(sum(coord**2))
      babs = sqrt(sum(beta**2))
      gamma = 1/sqrt(1-babs**2)

c     reverse boosts from electron to observer frame
      e_ph_ev = gamma*e_tmp+
     &     gamma*beta(1)*om_tmp(1)+
     &     gamma*beta(2)*om_tmp(2)+
     &     gamma*beta(3)*om_tmp(3)
      om_i(1) = gamma*beta(1)*e_tmp+
     &     (1+(gamma-1)*beta(1)**2/babs**2)*om_tmp(1)+
     &     ((gamma-1)*beta(1)*beta(2)/babs**2)*om_tmp(2)+
     &     ((gamma-1)*beta(1)*beta(3)/babs**2)*om_tmp(3)
      om_i(2) = gamma*beta(2)*e_tmp+
     &     ((gamma-1)*beta(1)*beta(2)/babs**2)*om_tmp(1)+
     &     (1+(gamma-1)*beta(2)**2/babs**2)*om_tmp(2)+
     &     ((gamma-1)*beta(3)*beta(2)/babs**2)*om_tmp(3)
      om_i(3) = gamma*beta(3)*e_tmp+
     &     ((gamma-1)*beta(1)*beta(3)/babs**2)*om_tmp(1)+
     &     ((gamma-1)*beta(2)*beta(3)/babs**2)*om_tmp(2)+
     &     (1+(gamma-1)*beta(3)**2/babs**2)*om_tmp(3)
      om_i = om_i/sqrt(sum(om_i**2))
 7391 end

      subroutine el_en(t_e,beta,gamma)
      implicit none
      real*8 t_e,beta,gamma,t_e_ev,x_e,r1,r2,r3,r4,z,z1,eta_p,eta_b
      real*8 eta
      integer ii
      t_e_ev=t_e/1.1609e4
      x_e=t_e_ev/5.11e5
      if(t_e_ev < 150 ) then
         ii=0
 11      r1=rand()
         r2=rand()
         z1=-3.*log(r1)
         ii=ii+1
         z=0.151*(1.+x_e*z1)**2*z1*(2.+x_e*z1)*r1
         if(r2**2 < z ) then
            eta=sqrt(x_e*z1*(2.+x_e*z1))
            goto 22
         else
            if(ii>500) then
               write(0,*)'zz ',ii,x_e,r1,r2,z1,z
               stop
            endif
            goto 11
         endif
 22      continue
      else
 44      r1=rand()
         r2=rand()
         r3=rand()
         r4=rand()
         eta_p=-x_e*log(r1*r2*r3)
         eta_b=-x_e*log(r1*r2*r3*r4)
         if(eta_b**2-eta_p**2 > 1.) then
            eta=eta_p
            goto 33
         else
            goto 44
         endif
 33      continue
      endif
      gamma=sqrt(eta**2+1.)
      beta=eta/gamma
      return
      end

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

      get_ene_idx = int(log10(ene_in)/de+1)
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
         ene(ne) = 10.**(de*(ne-1))
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
               num_cum(idx) = 0.
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

      CALL h5dopen_f(file_id, 'vex', temp_id, error)
      CALL h5dread_f(temp_id, H5T_NATIVE_DOUBLE, vex, dat_dim, error)
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

!     Clear NaNs
      abu(n_rad,:,:,:) = abu(n_rad-1,:,:,:)
      den(n_rad,:,:) = den(n_rad-1,:,:)
      vex(n_rad,:,:) = vex(n_rad-1,:,:)
      
!     Convert to number densities, expand homologously, and compute
!     electron number density
      do ii = 1, n_abu
         abu(:,:,:,ii) = den*abu(:,:,:,ii)/(uu*aw(ii))
         abu(:,:,:,ii) = abu(:,:,:,ii)
         ne_tot = ne_tot + abu(:,:,:,ii)*nz(ii)
! REDUCE HYDROGEN SCATTERING !         if (ii.eq.1) then
! REDUCE HYDROGEN SCATTERING !            ne_tot = ne_tot + abu(:,:,:,ii)*0.5*nz(ii)
! REDUCE HYDROGEN SCATTERING !         else
! REDUCE HYDROGEN SCATTERING !            ne_tot = ne_tot + abu(:,:,:,ii)*nz(ii)
! REDUCE HYDROGEN SCATTERING !         endif         
      end do
      
      dif_tht = tht(2)-tht(1)
      dif_phi = phi(2)-phi(1)
      do ii = 1,n_rad
         do jj = 1,n_tht
            do kk = 1,n_phi
               vol(ii,jj,kk) = ((rad(ii)+rad(ii+1))/2.)**2*
     &              sin((tht(jj)+tht(jj+1))/2.)*(rad(ii+1)-rad(ii))*
     &              dif_tht*dif_phi

               idx = (ii-1)*(n_tht*n_phi)+(jj-1)*n_phi+kk
c     44Ti (44Ti)
               if (src_abu.eq.0) then
                  if (idx.eq.1) then
                     num_cum(idx) = vol(ii,jj,kk)*abu(ii,jj,kk,13)
                  else
                     num_cum(idx) = num_cum(idx-1) +
     &                    vol(ii,jj,kk)*abu(ii,jj,kk,13)
                  endif
c     56Ni & 56Co (56Fe+56Co+56Ni+0.5*X)
               else if (src_abu.eq.1.or.src_abu.eq.2) then
                  if (idx.eq.1) then
                     num_cum(idx) = vol(ii,jj,kk)*(abu(ii,jj,kk,16)+
     &           abu(ii,jj,kk,17)+abu(ii,jj,kk,18)+0.5*abu(ii,jj,kk,19))
                  else
                     num_cum(idx) = num_cum(idx-1) +
     &                    vol(ii,jj,kk)*(abu(ii,jj,kk,16)+
     &           abu(ii,jj,kk,17)+abu(ii,jj,kk,18)+0.5*abu(ii,jj,kk,19))
                  endif
c     57Co (X)
               else if (src_abu.eq.3) then
                  if (idx.eq.1) then
                     num_cum(idx) = vol(ii,jj,kk)*abu(ii,jj,kk,19)
                  else
                     num_cum(idx) = num_cum(idx-1) +
     &                    vol(ii,jj,kk)*abu(ii,jj,kk,19)
                  endif
               endif
            end do
         end do
      end do
! dbg      print *,n_src, SUM(vol(:n_rad-1,:,:)*abu(:n_rad-1,:,:,9))
! dbg     ^ *uu*36/MSun

!     dbg print *, (rad(ii), ii = 1, n_rad1)
!     dbg print *, (tht(ii), ii = 1, n_tht1)
!     dbg print *, (phi(ii), ii = 1, n_phi1)
!     dbg      do jj = 10, 10
!     dbg         do kk = 100, 100
!     dbg            print *, (abu(ii,jj,kk,9), ii = 1, n_rad)
!     dbg         end do
!     dbg      end do
      end

      subroutine ini_44(e_ph_i_ev, ne)
      use constants
      implicit none
      real*8 e_ph_i_ev, tmp, cum_eff, ee(8), eff(8)
      integer ne, ii
c     44Ti AND 44Sc, from table of isotopes, endt90, ahmad06, chen11
c     http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=220044
c     http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=210044
      data ee/67.875d3, 78.337d3, 146.212d3, 1157.031d3, 1499.43d3,
     &     2144.2d3, 2656.41d3, 3301.2d3/
      data eff/94.4d-2, 96d-2, 0.089d-2, 99.9d-2, 0.912d-2, 0.0069d-2,
     &     0.115d-2, 0.0031d-2/

c      print *, 'Photons per decay:', sum(eff)
      tmp = real(ne,8)/real(n_pho,8)*sum(eff)
      ii = 0
      cum_eff = 0.d0
      do while (cum_eff+1e-13.le.tmp)
         ii = ii+1
         cum_eff = cum_eff + eff(ii)
      enddo
      e_ph_i_ev = ee(ii)
      end
 
      subroutine ini_56ni(e_ph_i_ev, ne)
      use constants
      implicit none
      real*8 e_ph_i_ev, tmp, cum_eff, ee(6), eff(6)
      integer ne, ii
c     from table of isotopes, junde92, junde99, junde11
c     http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=270056
      data ee/158.38d3, 269.50d3, 480.44d3, 749.95d3,811.85d3,1561.80d3/
      data eff/98.8d-2, 36.5d-2, 36.5d-2, 49.5d-2, 86.0d-2, 14.0d-2/

c      print *, 'Photons per decay:', sum(eff)
      tmp = real(ne,8)/real(n_pho,8)*sum(eff)
      ii = 0
      cum_eff = 0.d0
      do while (cum_eff+1e-13.le.tmp)
         ii = ii+1
         cum_eff = cum_eff + eff(ii)
      enddo
      e_ph_i_ev = ee(ii)
      end
     
      subroutine ini_56co(e_ph_i_ev, ne)
      use constants
      implicit none
      real*8 e_ph_i_ev, tmp, cum_eff, ee(45), eff(45)
      integer ne, ii
c     from table of isotopes, junde92, junde99, junde11
c     http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=270056
      data ee/263.410d3, 411.380d3, 486.540d3, 655.000d3, 674.700d3,
     & 733.511d3, 787.742d3, 846.771d3, 852.780d3, 896.531d3, 977.373d3,
     & 997.330d3, 1037.840d3, 1089.030d3, 1140.404d3, 1159.920d3,
     & 1175.102d3, 1198.780d3, 1238.282d3, 1272.200d3, 1335.589d3,
     & 1360.215d3, 1442.750d3, 1462.340d3, 1640.404d3, 1771.351d3,
     & 1810.772d3, 1963.714d3, 2015.181d3, 2034.755d3, 2113.123d3, 
     & 2212.933d3, 2276.360d3, 2373.700d3, 2598.459d3, 2657.450d3,
     & 3009.596d3, 3201.962d3, 3253.416d3, 3272.990d3, 3451.152d3,
     & 3547.930d3, 3600.700d3, 3611.800d3, 510.9989461d3/
      data eff/0.022d-2, 0.025d-2, 0.055d-2, 0.038d-2, 0.038d-2,
     & 0.200d-2, 0.310d-2, 100.000d-2, 0.050d-2, 0.070d-2, 1.430d-2,
     & 0.112d-2, 13.990d-2, 0.050d-2, 0.150d-2, 0.100d-2, 2.270d-2,
     & 0.050d-2, 67.600d-2, 0.020d-2, 0.120d-2, 4.330d-2, 0.200d-2,
     & 0.077d-2, 0.060d-2, 15.600d-2, 0.640d-2, 0.720d-2, 3.080d-2,
     & 7.880d-2, 0.380d-2, 0.350d-2, 0.110d-2, 0.080d-2, 17.280d-2,
     & 0.025d-2, 1.049d-2, 3.240d-2, 7.930d-2, 1.889d-2, 0.953d-2,
     & 0.198d-2, 0.018d-2, 0.010d-2, 38.00d-2/

c      print *, 'Photons per decay:', sum(eff)
      tmp = real(ne,8)/real(n_pho,8)*sum(eff)
      ii = 0
      cum_eff = 0.d0
      do while (cum_eff+1e-13.le.tmp)
         ii = ii+1
         cum_eff = cum_eff + eff(ii)
      enddo
      e_ph_i_ev = ee(ii)
      end

      subroutine ini_57co(e_ph_i_ev, ne)
      use constants
      implicit none
      real*8 e_ph_i_ev, tmp, cum_eff, ee(10), eff(10)
      integer ne, ii
c     44Ti AND 44Sc, from table of isotopes, endt90, ahmad06, chen11
c     http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=220044
c     http://nucleardata.nuclear.lu.se/toi/nuclide.asp?iZA=210044
      data ee/14.41300d3, 122.0614d3, 136.4743d3, 230.29d3, 339.54d3,
     &     352.36d3, 366.75d3, 569.92d3, 692.03d3, 706.40d3/
      data eff/9.16d-2, 85.60d-2, 10.68d-2, 0.0004d-2, 0.0139d-2,
     &     0.0132d-2, 0.0013d-2, 0.017d-2, 0.157d-2, 0.0253d-2/

c      print *, 'Photons per decay:', sum(eff)
      tmp = real(ne,8)/real(n_pho,8)*sum(eff)
      ii = 0
      cum_eff = 0.d0
      do while (cum_eff+1e-13.le.tmp)
         ii = ii+1
         cum_eff = cum_eff + eff(ii)
      enddo
      e_ph_i_ev = ee(ii)
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

! H and He are switched in Table 3 of hoeflich92!
      subroutine get_sig_hoeflich92()
      use constants
      implicit none
      integer ne,na,ns,nel
      real*8 sig_is
      open(unit=2299, file='sig_hoeflich92.txt')
      write(6,3322) 'Computing cross sections'
 3322 format(A)
      do ne=1,n_ene
         do na=1,n_abu
            sig(na,ne) = 0.
            nel = nz(na)
            call hoeflich92(nel,ene(ne),sig_is)
            sig(na,ne) = sig_is
         enddo
      enddo
      do na=1,n_abu
         write(2299,*) sig(na,:)
      enddo
      end
      subroutine hoeflich92(nel,e,s)
      use constants
      implicit none
      integer nel
      real*8 e,s
      if(nel.eq.1) then
         if(e.gt.1.d4) then
            s = 7.25d-5*(e/1.d5)**(-3.42)
         elseif(e.gt.1d3) then
            s = 1.92d-1*(e/1.d4)**(-3.33)
         elseif(e.gt.1d2) then
            s = 4.13d2*(e/1.d3)**(-2.94)
         else
            s = 3.59d5*(e/1.d2)**(-3.00)
         endif
      elseif(nel.eq.2) then
         if(e.gt.1.d4) then
            s = 1.55d-6*(e/1.d5)**(-3.47)
         elseif(e.gt.1d3) then
            s = 4.55d-3*(e/1.d4)**(-3.40)
         elseif(e.gt.1d2) then
            s = 1.15d1*(e/1.d3)**(-3.22)
         else
            s = 1.93d4*(e/1.d2)**(-3.00)
         endif
      elseif(nel.eq.6) then
         if(e.gt.1.d4) then
            s = 2.13d-2*(e/1.d5)**(-3.28)
         elseif(e.gt.1d3) then
            s = 4.03d1*(e/1.d4)**(-3.04)
         elseif(e.gt.2.83d2) then
            s = 4.45d4*(e/1.d3)**(-2.46)
         elseif(e.gt.1.d2) then
            s = 5.05d4*(e/2.83d2)**(-2.26)
         else
            s = 5.29d5*(e/1.d2)**(-3.00)
         endif
      elseif(nel.eq.8) then
         if(e.gt.1.d4) then
            s = 8.53d-2*(e/1.d5)**(-3.24)
         elseif(e.gt.1d3) then
            s = 1.47d2*(e/1.d4)**(-2.91)
         elseif(e.gt.5.31d2) then
            s = 1.19d5*(e/1.d3)**(-2.30)
         elseif(e.gt.1.d2) then
            s = 3.26d4*(e/5.31d2)**(-2.35)
         else
            s = 1.65d6*(e/1.d2)**(-3.00)
         endif
      elseif(nel.eq.10) then
         if(e.gt.1.d4) then
            s = 2.48d-1*(e/1.d5)**(-3.19)
         elseif(e.gt.1d3) then
            s = 3.80d2*(e/1.d4)**(-2.83)
         elseif(e.gt.8.66d2) then
            s = 2.57d5*(e/1.d3)**(-1.53)
         elseif(e.gt.1.d2) then
            s = 2.37d4*(e/8.66d2)**(-2.39)
         else
            s = 4.11d6*(e/1.d2)**(-3.00)
         endif
      elseif(nel.eq.12) then
         if(e.gt.1.d5) then
            s = 3.93d-4*(e/1.d6)**(-3.17)
         elseif(e.gt.1d4) then
            s = 5.76d-1*(e/1.d5)**(-3.15)
         elseif(e.gt.1.3d3) then
            s = 8.06d2*(e/1.d4)**(-2.81)
         elseif(e.gt.1.d3) then
            s = 1.79d4*(e/1.3d3)**(-2.88)
         elseif(e.gt.1.d2) then
            s = 3.85d4*(e/1.d3)**(-2.19)
         else
            s = 5.97d6*(e/1.d2)**(-3.00)
         endif
      elseif(nel.eq.14) then
         if(e.gt.1.d5) then
            s = 8.64d-4*(e/1.d6)**(-3.13)
         elseif(e.gt.1d4) then
            s = 1.16d0*(e/1.d5)**(-3.11)
         elseif(e.gt.1.83d3) then
            s = 1.49d3*(e/1.d4)**(-2.81)
         elseif(e.gt.1.d3) then
            s = 1.43d4*(e/1.83d3)**(-2.71)
         elseif(e.gt.1.d2) then
            s = 7.40d4*(e/1.d3)**(-1.56)
         else
            s = 2.66d6*(e/1.d2)**(-3.00)
         endif
      elseif(nel.eq.16) then
         if(e.gt.1.d5) then
            s = 1.63d-3*(e/1.d6)**(-3.13)
         elseif(e.gt.1d4) then
            s = 2.18d0*(e/1.d5)**(-3.07)
         elseif(e.gt.2.47d3) then
            s = 2.57d3*(e/1.d4)**(-2.79)
         elseif(e.gt.1.d3) then
            s = 1.16d4*(e/2.47d3)**(-2.68)
         elseif(e.gt.1.d2) then
            s = 1.31d5*(e/1.d3)**(-0.756)
         else
            s = 7.47d5*(e/1.d2)**(-3.00)
         endif
      elseif(nel.le.23) then
         if(e.gt.1.d5) then
            s = 5.13d-3*(e/1.d6)**(-3.06)
         elseif(e.gt.1d4) then
            s = 5.93d0*(e/1.d5)**(-3.01)
         elseif(e.gt.4.04d3) then
            s = 6.10d3*(e/1.d4)**(-2.75)
         elseif(e.gt.1.d3) then
            s = 8.09d3*(e/4.04d3)**(-2.51)
         elseif(e.gt.3.77d2) then
            s = 2.70d5*(e/1.d3)**(-2.25)
         elseif(e.gt.1.d2) then
            s = 2.53d6*(e/3.77d2)**0.364
         else
            s = 1.56d6*(e/1.d2)**(-3.00)
         endif
      else
         if(e.gt.1.d5) then
            s = 1.87d-2*(e/1.d6)**(-3.00)
         elseif(e.gt.1d4) then
            s = 1.87d1*(e/1.d5)**(-2.94)
         elseif(e.gt.7.11d3) then
            s = 1.62d4*(e/1.d4)**(-2.68)
         elseif(e.gt.1.d3) then
            s = 4.98d3*(e/7.11d3)**(-2.71)
         elseif(e.gt.8.47d2) then
            s = 1.01d6*(e/1.d3)**(-1.14)
         elseif(e.gt.7.12d2) then
            s = 1.05d6*(e/8.47d2)**1.47
         elseif(e.gt.1.d2) then
            s = 1.92d5*(e/7.12d2)**(-1.64)
         else
            s = 4.83d6*(e/1.d2)**(-3.00)
         endif
      endif
      s = s*1.d-24
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
     1 5.475E+04, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 2, 1, 1),I=1,6) /5.442E+01, 1.720E+00,
     1 1.369E+04, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 2, 2, 1),I=1,6) /2.459E+01, 5.996E+00,
     1 4.470E+03, 2.199E+00, 6.098E+00, 0.000E+00/
      DATA (PH1(I, 3, 1, 1),I=1,6) /1.225E+02, 3.871E+00,
     1 6.083E+03, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 3, 2, 1),I=1,6) /7.564E+01, 2.006E+01,
     1 3.201E+02, 7.391E+00, 2.916E+00, 0.000E+00/
      DATA (PH1(I, 3, 3, 1),I=1,6) /6.439E+01, 2.740E+01,
     1 1.564E+02, 3.382E+01, 1.490E+00, 0.000E+00/
      DATA (PH1(I, 3, 3, 2),I=1,6) /5.392E+00, 3.466E+00,
     1 4.774E+01, 2.035E+01, 4.423E+00, 0.000E+00/
      DATA (PH1(I, 4, 1, 1),I=1,6) /2.177E+02, 6.879E+00,
     1 3.422E+03, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 4, 2, 1),I=1,6) /1.539E+02, 1.760E+01,
     1 5.458E+02, 1.719E+01, 3.157E+00, 0.000E+00/
      DATA (PH1(I, 4, 3, 1),I=1,6) /1.299E+02, 4.759E+01,
     1 9.796E+01, 1.166E+03, 1.022E+00, 0.000E+00/
      DATA (PH1(I, 4, 3, 2),I=1,6) /1.821E+01, 6.685E+00,
     1 4.904E+01, 3.419E+01, 3.738E+00, 0.000E+00/
      DATA (PH1(I, 4, 4, 1),I=1,6) /1.193E+02, 4.093E+01,
     1 1.306E+02, 1.212E+02, 1.348E+00, 0.000E+00/
      DATA (PH1(I, 4, 4, 2),I=1,6) /9.323E+00, 3.427E+00,
     1 1.423E+03, 2.708E+00, 1.064E+01, 0.000E+00/
      DATA (PH1(I, 5, 1, 1),I=1,6) /3.402E+02, 1.075E+01,
     1 2.190E+03, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 5, 2, 1),I=1,6) /2.594E+02, 3.336E+01,
     1 2.846E+02, 2.163E+01, 2.624E+00, 0.000E+00/
      DATA (PH1(I, 5, 3, 1),I=1,6) /2.274E+02, 6.592E+01,
     1 8.605E+01, 1.906E+02, 1.210E+00, 0.000E+00/
      DATA (PH1(I, 5, 3, 2),I=1,6) /3.793E+01, 1.054E+01,
     1 4.273E+01, 4.433E+01, 3.501E+00, 0.000E+00/
      DATA (PH1(I, 5, 4, 1),I=1,6) /2.098E+02, 5.984E+01,
     1 1.037E+02, 7.915E+01, 1.436E+00, 0.000E+00/
      DATA (PH1(I, 5, 4, 2),I=1,6) /2.516E+01, 2.403E+00,
     1 1.530E+02, 9.203E+00, 9.374E+00, 0.000E+00/
      DATA (PH1(I, 5, 5, 1),I=1,6) /1.940E+02, 6.155E+01,
     1 9.698E+01, 7.354E+01, 1.438E+00, 0.000E+00/
      DATA (PH1(I, 5, 5, 2),I=1,6) /1.405E+01, 6.495E+00,
     1 1.525E+04, 1.352E+00, 1.218E+01, 0.000E+00/
      DATA (PH1(I, 5, 5, 3),I=1,6) /8.298E+00, 6.658E+00,
     1 3.643E+02, 9.500E+00, 5.469E+00, 5.608E-01/
      DATA (PH1(I, 6, 1, 1),I=1,6) /4.900E+02, 1.548E+01,
     1 1.521E+03, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 6, 2, 1),I=1,6) /3.921E+02, 4.624E+01,
     1 2.344E+02, 2.183E+01, 2.581E+00, 0.000E+00/
      DATA (PH1(I, 6, 3, 1),I=1,6) /3.522E+02, 8.412E+01,
     1 8.111E+01, 7.459E+01, 1.428E+00, 0.000E+00/
      DATA (PH1(I, 6, 3, 2),I=1,6) /6.449E+01, 7.843E+00,
     1 7.109E+01, 2.828E+01, 4.754E+00, 0.000E+00/
      DATA (PH1(I, 6, 4, 1),I=1,6) /3.289E+02, 8.370E+01,
     1 8.067E+01, 7.471E+01, 1.442E+00, 0.000E+00/
      DATA (PH1(I, 6, 4, 2),I=1,6) /4.789E+01, 3.942E+00,
     1 1.219E+02, 1.499E+01, 7.489E+00, 0.000E+00/
      DATA (PH1(I, 6, 5, 1),I=1,6) /3.076E+02, 9.113E+01,
     1 6.649E+01, 9.609E+01, 1.338E+00, 0.000E+00/
      DATA (PH1(I, 6, 5, 2),I=1,6) /3.047E+01, 2.991E+00,
     1 1.184E+03, 3.085E+00, 1.480E+01, 0.000E+00/
      DATA (PH1(I, 6, 5, 3),I=1,6) /2.438E+01, 1.094E+01,
     1 1.792E+02, 3.308E+01, 4.150E+00, 5.276E-01/
      DATA (PH1(I, 6, 6, 1),I=1,6) /2.910E+02, 8.655E+01,
     1 7.421E+01, 5.498E+01, 1.503E+00, 0.000E+00/
      DATA (PH1(I, 6, 6, 2),I=1,6) /1.939E+01, 1.026E+01,
     1 4.564E+03, 1.568E+00, 1.085E+01, 0.000E+00/
      DATA (PH1(I, 6, 6, 3),I=1,6) /1.126E+01, 9.435E+00,
     1 1.152E+03, 5.687E+00, 6.336E+00, 4.474E-01/
      DATA (PH1(I, 7, 1, 1),I=1,6) /6.671E+02, 2.108E+01,
     1 1.117E+03, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 7, 2, 1),I=1,6) /5.521E+02, 6.943E+01,
     1 1.519E+02, 2.627E+01, 2.315E+00, 0.000E+00/
      DATA (PH1(I, 7, 3, 1),I=1,6) /5.043E+02, 1.060E+02,
     1 7.304E+01, 5.547E+01, 1.538E+00, 0.000E+00/
      DATA (PH1(I, 7, 3, 2),I=1,6) /9.789E+01, 1.862E+01,
     1 3.447E+01, 4.231E+01, 3.606E+00, 0.000E+00/
      DATA (PH1(I, 7, 4, 1),I=1,6) /4.753E+02, 1.070E+02,
     1 7.046E+01, 5.342E+01, 1.552E+00, 0.000E+00/
      DATA (PH1(I, 7, 4, 2),I=1,6) /7.747E+01, 6.225E+00,
     1 1.110E+02, 1.733E+01, 6.719E+00, 0.000E+00/
      DATA (PH1(I, 7, 5, 1),I=1,6) /4.473E+02, 1.220E+02,
     1 5.235E+01, 9.428E+01, 1.335E+00, 0.000E+00/
      DATA (PH1(I, 7, 5, 2),I=1,6) /5.545E+01, 5.853E+00,
     1 1.908E+02, 6.264E+00, 9.711E+00, 0.000E+00/
      DATA (PH1(I, 7, 5, 3),I=1,6) /4.745E+01, 1.925E+01,
     1 9.400E+01, 1.152E+02, 3.194E+00, 5.496E-01/
      DATA (PH1(I, 7, 6, 1),I=1,6) /4.236E+02, 1.242E+02,
     1 5.002E+01, 9.100E+01, 1.335E+00, 0.000E+00/
      DATA (PH1(I, 7, 6, 2),I=1,6) /3.796E+01, 1.094E+01,
     1 7.483E+02, 2.793E+00, 9.956E+00, 0.000E+00/
      DATA (PH1(I, 7, 6, 3),I=1,6) /2.960E+01, 1.827E+01,
     1 1.724E+02, 8.893E+01, 3.348E+00, 4.209E-01/
      DATA (PH1(I, 7, 7, 1),I=1,6) /4.048E+02, 1.270E+02,
     1 4.748E+01, 1.380E+02, 1.252E+00, 0.000E+00/
      DATA (PH1(I, 7, 7, 2),I=1,6) /2.541E+01, 1.482E+01,
     1 7.722E+02, 2.306E+00, 9.139E+00, 0.000E+00/
      DATA (PH1(I, 7, 7, 3),I=1,6) /1.453E+01, 1.164E+01,
     1 1.029E+04, 2.361E+00, 8.821E+00, 4.239E-01/
      DATA (PH1(I, 8, 1, 1),I=1,6) /8.714E+02, 2.754E+01,
     1 8.554E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 8, 2, 1),I=1,6) /7.393E+02, 8.709E+01,
     1 1.329E+02, 2.535E+01, 2.336E+00, 0.000E+00/
      DATA (PH1(I, 8, 3, 1),I=1,6) /6.837E+02, 1.354E+02,
     1 6.029E+01, 5.682E+01, 1.533E+00, 0.000E+00/
      DATA (PH1(I, 8, 3, 2),I=1,6) /1.381E+02, 9.141E+00,
     1 6.896E+01, 3.896E+01, 4.943E+00, 0.000E+00/
      DATA (PH1(I, 8, 4, 1),I=1,6) /6.491E+02, 1.377E+02,
     1 5.735E+01, 5.486E+01, 1.540E+00, 0.000E+00/
      DATA (PH1(I, 8, 4, 2),I=1,6) /1.139E+02, 8.467E+00,
     1 9.641E+01, 2.287E+01, 6.061E+00, 0.000E+00/
      DATA (PH1(I, 8, 5, 1),I=1,6) /6.144E+02, 1.593E+02,
     1 4.123E+01, 1.141E+02, 1.287E+00, 0.000E+00/
      DATA (PH1(I, 8, 5, 2),I=1,6) /8.737E+01, 7.942E+00,
     1 1.063E+02, 9.708E+00, 8.183E+00, 0.000E+00/
      DATA (PH1(I, 8, 5, 3),I=1,6) /7.741E+01, 1.937E+01,
     1 1.619E+02, 7.312E+01, 3.648E+00, 3.760E-02/
      DATA (PH1(I, 8, 6, 1),I=1,6) /5.840E+02, 1.620E+02,
     1 3.939E+01, 1.104E+02, 1.289E+00, 0.000E+00/
      DATA (PH1(I, 8, 6, 2),I=1,6) /6.551E+01, 1.594E+01,
     1 1.217E+02, 6.156E+00, 7.271E+00, 0.000E+00/
      DATA (PH1(I, 8, 6, 3),I=1,6) /5.494E+01, 2.067E+01,
     1 2.318E+02, 7.136E+01, 3.618E+00, 5.538E-02/
      DATA (PH1(I, 8, 7, 1),I=1,6) /5.581E+02, 1.690E+02,
     1 3.584E+01, 1.894E+02, 1.185E+00, 0.000E+00/
      DATA (PH1(I, 8, 7, 2),I=1,6) /4.599E+01, 1.759E+01,
     1 1.962E+02, 4.020E+00, 7.999E+00, 0.000E+00/
      DATA (PH1(I, 8, 7, 3),I=1,6) /3.512E+01, 1.745E+01,
     1 5.186E+02, 1.728E+01, 4.995E+00, 2.182E-02/
      DATA (PH1(I, 8, 8, 1),I=1,6) /5.380E+02, 1.774E+02,
     1 3.237E+01, 3.812E+02, 1.083E+00, 0.000E+00/
      DATA (PH1(I, 8, 8, 2),I=1,6) /2.848E+01, 1.994E+01,
     1 2.415E+02, 3.241E+00, 8.037E+00, 0.000E+00/
      DATA (PH1(I, 8, 8, 3),I=1,6) /1.362E+01, 1.391E+01,
     1 1.220E+05, 1.364E+00, 1.140E+01, 4.103E-01/
      DATA (PH1(I, 9, 1, 1),I=1,6) /1.103E+03, 3.485E+01,
     1 6.759E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I, 9, 2, 1),I=1,6) /9.539E+02, 1.131E+02,
     1 1.039E+02, 2.657E+01, 2.255E+00, 0.000E+00/
      DATA (PH1(I, 9, 3, 1),I=1,6) /8.905E+02, 1.711E+02,
     1 4.890E+01, 6.137E+01, 1.501E+00, 0.000E+00/
      DATA (PH1(I, 9, 3, 2),I=1,6) /1.852E+02, 6.423E+00,
     1 8.497E+01, 6.075E+01, 5.214E+00, 0.000E+00/
      DATA (PH1(I, 9, 4, 1),I=1,6) /8.502E+02, 1.737E+02,
     1 4.668E+01, 5.876E+01, 1.511E+00, 0.000E+00/
      DATA (PH1(I, 9, 4, 2),I=1,6) /1.572E+02, 1.146E+01,
     1 8.453E+01, 2.457E+01, 5.771E+00, 0.000E+00/
      DATA (PH1(I, 9, 5, 1),I=1,6) /8.091E+02, 1.925E+02,
     1 3.680E+01, 7.933E+01, 1.377E+00, 0.000E+00/
      DATA (PH1(I, 9, 5, 2),I=1,6) /1.262E+02, 1.889E+01,
     1 7.149E+01, 1.194E+01, 6.030E+00, 0.000E+00/
      DATA (PH1(I, 9, 5, 3),I=1,6) /1.142E+02, 2.528E+01,
     1 1.419E+02, 7.089E+01, 3.628E+00, 3.418E-01/
      DATA (PH1(I, 9, 6, 1),I=1,6) /7.709E+02, 1.603E+02,
     1 5.302E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I, 9, 6, 2),I=1,6) /9.957E+01, 2.629E+01,
     1 3.563E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I, 9, 6, 3),I=1,6) /8.714E+01, 2.292E+01,
     1 3.005E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I, 9, 7, 1),I=1,6) /7.392E+02, 2.055E+02,
     1 3.144E+01, 1.230E+02, 1.263E+00, 0.000E+00/
      DATA (PH1(I, 9, 7, 2),I=1,6) /7.610E+01, 2.307E+01,
     1 7.140E+01, 7.282E+00, 6.543E+00, 0.000E+00/
      DATA (PH1(I, 9, 7, 3),I=1,6) /6.271E+01, 2.744E+01,
     1 2.458E+02, 7.329E+01, 3.596E+00, 1.390E-01/
      DATA (PH1(I, 9, 8, 1),I=1,6) /7.122E+02, 1.660E+02,
     1 4.798E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I, 9, 8, 2),I=1,6) /5.459E+01, 2.882E+01,
     1 2.439E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I, 9, 8, 3),I=1,6) /3.497E+01, 2.658E+01,
     1 2.649E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I, 9, 9, 1),I=1,6) /6.940E+02, 2.390E+02,
     1 2.295E+01, 1.257E+03, 9.638E-01, 0.000E+00/
      DATA (PH1(I, 9, 9, 2),I=1,6) /3.786E+01, 2.568E+01,
     1 1.097E+02, 4.297E+00, 7.303E+00, 0.000E+00/
      DATA (PH1(I, 9, 9, 3),I=1,6) /1.742E+01, 1.658E+01,
     1 2.775E+05, 1.242E+00, 1.249E+01, 3.857E-01/
      DATA (PH1(I,10, 1, 1),I=1,6) /1.362E+03, 4.304E+01,
     1 5.475E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,10, 2, 1),I=1,6) /1.196E+03, 1.586E+02,
     1 6.695E+01, 3.352E+01, 2.002E+00, 0.000E+00/
      DATA (PH1(I,10, 3, 1),I=1,6) /1.125E+03, 2.193E+02,
     1 3.719E+01, 8.181E+01, 1.396E+00, 0.000E+00/
      DATA (PH1(I,10, 3, 2),I=1,6) /2.391E+02, 2.859E+01,
     1 2.897E+01, 3.836E+01, 3.992E+00, 0.000E+00/
      DATA (PH1(I,10, 4, 1),I=1,6) /1.078E+03, 1.886E+02,
     1 5.011E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,10, 4, 2),I=1,6) /2.073E+02, 2.822E+01,
     1 4.978E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,10, 5, 1),I=1,6) /1.031E+03, 2.308E+02,
     1 3.232E+01, 7.167E+01, 1.411E+00, 0.000E+00/
      DATA (PH1(I,10, 5, 2),I=1,6) /1.719E+02, 2.496E+01,
     1 5.575E+01, 1.427E+01, 5.562E+00, 0.000E+00/
      DATA (PH1(I,10, 5, 3),I=1,6) /1.579E+02, 2.740E+01,
     1 1.717E+02, 6.423E+01, 3.819E+00, 1.680E-01/
      DATA (PH1(I,10, 6, 1),I=1,6) /9.873E+02, 2.343E+02,
     1 3.097E+01, 6.907E+01, 1.416E+00, 0.000E+00/
      DATA (PH1(I,10, 6, 2),I=1,6) /1.415E+02, 2.615E+01,
     1 5.421E+01, 1.117E+01, 5.907E+00, 0.000E+00/
      DATA (PH1(I,10, 6, 3),I=1,6) /1.262E+02, 3.316E+01,
     1 1.962E+02, 7.397E+01, 3.563E+00, 2.895E-01/
      DATA (PH1(I,10, 7, 1),I=1,6) /9.480E+02, 2.360E+02,
     1 3.027E+01, 6.734E+01, 1.422E+00, 0.000E+00/
      DATA (PH1(I,10, 7, 2),I=1,6) /1.132E+02, 3.092E+01,
     1 4.190E+01, 1.139E+01, 5.554E+00, 0.000E+00/
      DATA (PH1(I,10, 7, 3),I=1,6) /9.712E+01, 3.494E+01,
     1 2.214E+02, 7.203E+01, 3.563E+00, 2.573E-01/
      DATA (PH1(I,10, 8, 1),I=1,6) /9.131E+02, 2.387E+02,
     1 2.943E+01, 6.525E+01, 1.424E+00, 0.000E+00/
      DATA (PH1(I,10, 8, 2),I=1,6) /8.721E+01, 3.192E+01,
     1 3.797E+01, 1.054E+01, 5.652E+00, 0.000E+00/
      DATA (PH1(I,10, 8, 3),I=1,6) /6.346E+01, 3.485E+01,
     1 2.451E+02, 6.937E+01, 3.645E+00, 2.042E-01/
      DATA (PH1(I,10, 9, 1),I=1,6) /8.831E+02, 2.452E+02,
     1 2.783E+01, 6.075E+01, 1.420E+00, 0.000E+00/
      DATA (PH1(I,10, 9, 2),I=1,6) /6.374E+01, 3.187E+01,
     1 4.025E+01, 8.495E+00, 6.038E+00, 0.000E+00/
      DATA (PH1(I,10, 9, 3),I=1,6) /4.096E+01, 3.428E+01,
     1 2.766E+02, 4.179E+01, 4.029E+00, 3.052E-01/
      DATA (PH1(I,10,10, 1),I=1,6) /8.701E+02, 3.144E+02,
     1 1.664E+01, 2.042E+05, 8.450E-01, 0.000E+00/
      DATA (PH1(I,10,10, 2),I=1,6) /4.847E+01, 3.204E+01,
     1 5.615E+01, 5.808E+00, 6.678E+00, 0.000E+00/
      DATA (PH1(I,10,10, 3),I=1,6) /2.156E+01, 2.000E+01,
     1 1.691E+04, 2.442E+00, 1.043E+01, 3.345E-01/
      DATA (PH1(I,11, 1, 1),I=1,6) /1.649E+03, 5.211E+01,
     1 4.525E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,11, 2, 1),I=1,6) /1.465E+03, 2.268E+02,
     1 3.995E+01, 5.315E+01, 1.678E+00, 0.000E+00/
      DATA (PH1(I,11, 3, 1),I=1,6) /1.386E+03, 2.406E+02,
     1 3.850E+01, 5.198E+01, 1.575E+00, 0.000E+00/
      DATA (PH1(I,11, 3, 2),I=1,6) /2.999E+02, 1.758E+01,
     1 4.531E+01, 4.549E+01, 4.689E+00, 0.000E+00/
      DATA (PH1(I,11, 4, 1),I=1,6) /1.335E+03, 2.464E+02,
     1 3.613E+01, 4.968E+01, 1.579E+00, 0.000E+00/
      DATA (PH1(I,11, 4, 2),I=1,6) /2.642E+02, 2.097E+01,
     1 6.127E+01, 2.644E+01, 5.246E+00, 0.000E+00/
      DATA (PH1(I,11, 5, 1),I=1,6) /1.281E+03, 2.777E+02,
     1 2.745E+01, 7.512E+01, 1.397E+00, 0.000E+00/
      DATA (PH1(I,11, 5, 2),I=1,6) /2.244E+02, 2.947E+01,
     1 4.802E+01, 1.621E+01, 5.395E+00, 0.000E+00/
      DATA (PH1(I,11, 5, 3),I=1,6) /2.085E+02, 3.552E+01,
     1 1.374E+02, 7.062E+01, 3.675E+00, 1.613E-01/
      DATA (PH1(I,11, 6, 1),I=1,6) /1.230E+03, 2.806E+02,
     1 2.654E+01, 7.237E+01, 1.405E+00, 0.000E+00/
      DATA (PH1(I,11, 6, 2),I=1,6) /1.899E+02, 3.303E+01,
     1 4.218E+01, 1.377E+01, 5.458E+00, 0.000E+00/
      DATA (PH1(I,11, 6, 3),I=1,6) /1.722E+02, 4.212E+01,
     1 1.650E+02, 9.649E+01, 3.367E+00, 2.182E-01/
      DATA (PH1(I,11, 7, 1),I=1,6) /1.185E+03, 2.360E+02,
     1 3.753E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,11, 7, 2),I=1,6) /1.567E+02, 3.868E+01,
     1 2.607E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,11, 7, 3),I=1,6) /1.384E+02, 3.523E+01,
     1 3.028E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,11, 8, 1),I=1,6) /1.143E+03, 2.841E+02,
     1 2.551E+01, 6.977E+01, 1.414E+00, 0.000E+00/
      DATA (PH1(I,11, 8, 2),I=1,6) /1.269E+02, 3.807E+01,
     1 3.144E+01, 1.250E+01, 5.405E+00, 0.000E+00/
      DATA (PH1(I,11, 8, 3),I=1,6) /9.892E+01, 4.314E+01,
     1 2.232E+02, 8.352E+01, 3.511E+00, 1.644E-01/
      DATA (PH1(I,11, 9, 1),I=1,6) /1.118E+03, 2.426E+02,
     1 3.466E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,11, 9, 2),I=1,6) /9.945E+01, 4.130E+01,
     1 1.933E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,11, 9, 3),I=1,6) /7.162E+01, 3.957E+01,
     1 2.748E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,11,10, 1),I=1,6) /1.097E+03, 3.266E+02,
     1 1.889E+01, 2.527E+02, 1.120E+00, 0.000E+00/
      DATA (PH1(I,11,10, 2),I=1,6) /7.347E+01, 3.939E+01,
     1 2.629E+01, 1.103E+01, 5.618E+00, 0.000E+00/
      DATA (PH1(I,11,10, 3),I=1,6) /4.729E+01, 3.911E+01,
     1 3.551E+02, 2.257E+01, 4.643E+00, 2.900E-01/
      DATA (PH1(I,11,11, 1),I=1,6) /1.079E+03, 4.216E+02,
     1 1.119E+01, 5.642E+07, 7.736E-01, 0.000E+00/
      DATA (PH1(I,11,11, 2),I=1,6) /7.084E+01, 4.537E+01,
     1 1.142E+01, 2.395E+02, 3.380E+00, 0.000E+00/
      DATA (PH1(I,11,11, 3),I=1,6) /3.814E+01, 3.655E+01,
     1 2.486E+02, 3.222E+02, 3.570E+00, 1.465E-01/
      DATA (PH1(I,11,11, 4),I=1,6) /5.139E+00, 5.968E+00,
     1 1.460E+00, 2.557E+07, 3.789E+00, 0.000E+00/
      DATA (PH1(I,12, 1, 1),I=1,6) /1.963E+03, 6.203E+01,
     1 3.802E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,12, 2, 1),I=1,6) /1.762E+03, 2.042E+02,
     1 6.140E+01, 2.778E+01, 2.161E+00, 0.000E+00/
      DATA (PH1(I,12, 3, 1),I=1,6) /1.675E+03, 2.858E+02,
     1 3.290E+01, 5.384E+01, 1.560E+00, 0.000E+00/
      DATA (PH1(I,12, 3, 2),I=1,6) /3.675E+02, 4.041E+01,
     1 2.228E+01, 3.904E+01, 3.986E+00, 0.000E+00/
      DATA (PH1(I,12, 4, 1),I=1,6) /1.618E+03, 2.671E+02,
     1 3.699E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,12, 4, 2),I=1,6) /3.282E+02, 3.981E+01,
     1 3.936E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,12, 5, 1),I=1,6) /1.558E+03, 3.082E+02,
     1 2.721E+01, 5.108E+01, 1.539E+00, 0.000E+00/
      DATA (PH1(I,12, 5, 2),I=1,6) /2.839E+02, 3.967E+01,
     1 3.755E+01, 1.873E+01, 4.950E+00, 0.000E+00/
      DATA (PH1(I,12, 5, 3),I=1,6) /2.660E+02, 2.461E+01,
     1 3.546E+02, 6.160E+01, 4.366E+00, 6.452E-03/
      DATA (PH1(I,12, 6, 1),I=1,6) /1.501E+03, 2.742E+02,
     1 3.407E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,12, 6, 2),I=1,6) /2.444E+02, 4.341E+01,
     1 2.832E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,12, 6, 3),I=1,6) /2.249E+02, 4.016E+01,
     1 2.442E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,12, 7, 1),I=1,6) /1.449E+03, 2.777E+02,
     1 3.278E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,12, 7, 2),I=1,6) /2.076E+02, 4.500E+01,
     1 2.434E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,12, 7, 3),I=1,6) /1.865E+02, 4.155E+01,
     1 2.901E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,12, 8, 1),I=1,6) /1.400E+03, 3.259E+02,
     1 2.350E+01, 6.173E+01, 1.464E+00, 0.000E+00/
      DATA (PH1(I,12, 8, 2),I=1,6) /1.735E+02, 4.710E+01,
     1 2.448E+01, 1.636E+01, 4.951E+00, 0.000E+00/
      DATA (PH1(I,12, 8, 3),I=1,6) /1.413E+02, 5.114E+01,
     1 2.134E+02, 9.563E+01, 3.438E+00, 5.419E-03/
      DATA (PH1(I,12, 9, 1),I=1,6) /1.356E+03, 3.345E+02,
     1 2.212E+01, 7.147E+01, 1.413E+00, 0.000E+00/
      DATA (PH1(I,12, 9, 2),I=1,6) /1.411E+02, 4.810E+01,
     1 2.118E+01, 1.675E+01, 4.941E+00, 0.000E+00/
      DATA (PH1(I,12, 9, 3),I=1,6) /1.093E+02, 5.137E+01,
     1 2.223E+02, 9.301E+01, 3.513E+00, 5.276E-03/
      DATA (PH1(I,12,10, 1),I=1,6) /1.336E+03, 3.607E+02,
     1 1.877E+01, 1.401E+02, 1.238E+00, 0.000E+00/
      DATA (PH1(I,12,10, 2),I=1,6) /1.111E+02, 4.804E+01,
     1 1.914E+01, 1.613E+01, 5.050E+00, 0.000E+00/
      DATA (PH1(I,12,10, 3),I=1,6) /8.014E+01, 5.319E+01,
     1 2.104E+02, 7.813E+01, 3.630E+00, 2.851E-01/
      DATA (PH1(I,12,11, 1),I=1,6) /1.320E+03, 3.709E+02,
     1 1.767E+01, 1.346E+02, 1.225E+00, 0.000E+00/
      DATA (PH1(I,12,11, 2),I=1,6) /9.881E+01, 5.142E+01,
     1 1.290E+01, 5.775E+01, 3.903E+00, 0.000E+00/
      DATA (PH1(I,12,11, 3),I=1,6) /6.569E+01, 4.940E+01,
     1 2.049E+02, 4.112E+03, 2.995E+00, 2.223E-05/
      DATA (PH1(I,12,11, 4),I=1,6) /1.504E+01, 8.139E+00,
     1 3.278E+00, 4.341E+07, 3.610E+00, 0.000E+00/
      DATA (PH1(I,12,12, 1),I=1,6) /1.311E+03, 2.711E+02,
     1 3.561E+01, 2.374E+01, 1.952E+00, 0.000E+00/
      DATA (PH1(I,12,12, 2),I=1,6) /9.400E+01, 4.587E+01,
     1 1.671E+01, 2.389E+01, 4.742E+00, 0.000E+00/
      DATA (PH1(I,12,12, 3),I=1,6) /5.490E+01, 4.937E+01,
     1 2.023E+02, 1.079E+04, 2.960E+00, 1.463E-02/
      DATA (PH1(I,12,12, 4),I=1,6) /7.646E+00, 9.393E+00,
     1 3.034E+00, 2.625E+07, 3.923E+00, 0.000E+00/
      DATA (PH1(I,13, 1, 1),I=1,6) /2.304E+03, 7.281E+01,
     1 3.239E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,13, 2, 1),I=1,6) /2.086E+03, 2.738E+02,
     1 4.036E+01, 3.567E+01, 1.915E+00, 0.000E+00/
      DATA (PH1(I,13, 3, 1),I=1,6) /1.992E+03, 3.137E+02,
     1 3.272E+01, 4.295E+01, 1.676E+00, 0.000E+00/
      DATA (PH1(I,13, 3, 2),I=1,6) /4.421E+02, 4.434E+01,
     1 2.098E+01, 3.913E+01, 4.062E+00, 0.000E+00/
      DATA (PH1(I,13, 4, 1),I=1,6) /1.929E+03, 3.117E+02,
     1 3.223E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,13, 4, 2),I=1,6) /3.994E+02, 4.639E+01,
     1 3.516E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,13, 5, 1),I=1,6) /1.862E+03, 3.305E+02,
     1 2.857E+01, 3.660E+01, 1.712E+00, 0.000E+00/
      DATA (PH1(I,13, 5, 2),I=1,6) /3.502E+02, 4.550E+01,
     1 3.393E+01, 1.945E+01, 4.926E+00, 0.000E+00/
      DATA (PH1(I,13, 5, 3),I=1,6) /3.301E+02, 4.084E+01,
     1 1.706E+02, 6.480E+01, 3.927E+00, 1.386E-01/
      DATA (PH1(I,13, 6, 1),I=1,6) /1.799E+03, 3.193E+02,
     1 2.989E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,13, 6, 2),I=1,6) /3.065E+02, 5.016E+01,
     1 2.604E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,13, 6, 3),I=1,6) /2.846E+02, 4.719E+01,
     1 2.229E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,13, 7, 1),I=1,6) /1.739E+03, 3.396E+02,
     1 2.653E+01, 3.478E+01, 1.721E+00, 0.000E+00/
      DATA (PH1(I,13, 7, 2),I=1,6) /2.662E+02, 5.726E+01,
     1 2.274E+01, 1.983E+01, 4.594E+00, 0.000E+00/
      DATA (PH1(I,13, 7, 3),I=1,6) /2.414E+02, 5.746E+01,
     1 1.941E+02, 8.747E+01, 3.460E+00, 1.206E-01/
      DATA (PH1(I,13, 8, 1),I=1,6) /1.688E+03, 3.268E+02,
     1 2.786E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,13, 8, 2),I=1,6) /2.268E+02, 5.340E+01,
     1 1.982E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,13, 8, 3),I=1,6) /1.905E+02, 5.032E+01,
     1 2.888E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,13, 9, 1),I=1,6) /1.634E+03, 3.432E+02,
     1 2.572E+01, 3.376E+01, 1.730E+00, 0.000E+00/
      DATA (PH1(I,13, 9, 2),I=1,6) /1.903E+02, 5.651E+01,
     1 1.865E+01, 1.904E+01, 4.763E+00, 0.000E+00/
      DATA (PH1(I,13, 9, 3),I=1,6) /1.538E+02, 6.049E+01,
     1 2.132E+02, 8.633E+01, 3.522E+00, 1.041E-01/
      DATA (PH1(I,13,10, 1),I=1,6) /1.604E+03, 3.482E+02,
     1 2.494E+01, 3.260E+01, 1.732E+00, 0.000E+00/
      DATA (PH1(I,13,10, 2),I=1,6) /1.558E+02, 5.735E+01,
     1 1.619E+01, 1.989E+01, 4.749E+00, 0.000E+00/
      DATA (PH1(I,13,10, 3),I=1,6) /1.200E+02, 6.033E+01,
     1 2.186E+02, 8.353E+01, 3.609E+00, 9.152E-02/
      DATA (PH1(I,13,11, 1),I=1,6) /1.583E+03, 3.200E+02,
     1 2.994E+01, 2.757E+01, 1.876E+00, 0.000E+00/
      DATA (PH1(I,13,11, 2),I=1,6) /1.407E+02, 5.629E+01,
     1 1.398E+01, 3.550E+01, 4.267E+00, 0.000E+00/
      DATA (PH1(I,13,11, 3),I=1,6) /1.026E+02, 6.016E+01,
     1 1.990E+02, 5.082E+02, 3.110E+00, 2.016E-02/
      DATA (PH1(I,13,11, 4),I=1,6) /2.845E+01, 1.027E+01,
     1 4.915E+00, 1.990E+06, 3.477E+00, 0.000E+00/
      DATA (PH1(I,13,12, 1),I=1,6) /1.571E+03, 3.049E+02,
     1 3.295E+01, 2.905E+01, 1.898E+00, 0.000E+00/
      DATA (PH1(I,13,12, 2),I=1,6) /1.281E+02, 5.171E+01,
     1 1.642E+01, 2.275E+01, 4.810E+00, 0.000E+00/
      DATA (PH1(I,13,12, 3),I=1,6) /8.997E+01, 6.154E+01,
     1 1.842E+02, 2.404E+03, 2.920E+00, 7.839E-03/
      DATA (PH1(I,13,12, 4),I=1,6) /1.883E+01, 1.012E+01,
     1 6.324E+00, 2.195E+02, 4.481E+00, 0.000E+00/
      DATA (PH1(I,13,13, 1),I=1,6) /1.567E+03, 3.670E+02,
     1 2.206E+01, 4.405E+01, 1.588E+00, 0.000E+00/
      DATA (PH1(I,13,13, 2),I=1,6) /1.256E+02, 5.594E+01,
     1 1.425E+01, 3.094E+01, 4.399E+00, 0.000E+00/
      DATA (PH1(I,13,13, 3),I=1,6) /8.040E+01, 6.445E+01,
     1 1.735E+02, 1.131E+04, 2.762E+00, 2.337E-02/
      DATA (PH1(I,13,13, 4),I=1,6) /1.133E+01, 1.204E+01,
     1 5.384E+00, 4.341E+02, 4.088E+00, 0.000E+00/
      DATA (PH1(I,13,13, 5),I=1,6) /5.986E+00, 1.860E+01,
     1 1.828E+02, 2.797E+00, 1.084E+01, 3.076E-01/
      DATA (PH1(I,14, 1, 1),I=1,6) /2.673E+03, 8.447E+01,
     1 2.793E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,14, 2, 1),I=1,6) /2.438E+03, 2.752E+02,
     1 4.754E+01, 2.848E+01, 2.135E+00, 0.000E+00/
      DATA (PH1(I,14, 3, 1),I=1,6) /2.336E+03, 3.150E+02,
     1 3.880E+01, 3.034E+01, 1.925E+00, 0.000E+00/
      DATA (PH1(I,14, 3, 2),I=1,6) /5.235E+02, 4.847E+01,
     1 1.966E+01, 3.851E+01, 4.148E+00, 0.000E+00/
      DATA (PH1(I,14, 4, 1),I=1,6) /2.268E+03, 3.599E+02,
     1 2.832E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,14, 4, 2),I=1,6) /4.761E+02, 5.350E+01,
     1 3.155E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,14, 5, 1),I=1,6) /2.194E+03, 3.356E+02,
     1 3.343E+01, 2.487E+01, 1.982E+00, 0.000E+00/
      DATA (PH1(I,14, 5, 2),I=1,6) /4.234E+02, 6.098E+01,
     1 2.624E+01, 2.213E+01, 4.520E+00, 0.000E+00/
      DATA (PH1(I,14, 5, 3),I=1,6) /4.014E+02, 3.950E+01,
     1 2.223E+02, 6.521E+01, 4.118E+00, 4.760E-05/
      DATA (PH1(I,14, 6, 1),I=1,6) /2.125E+03, 3.679E+02,
     1 2.641E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,14, 6, 2),I=1,6) /3.756E+02, 5.745E+01,
     1 2.394E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,14, 6, 3),I=1,6) /3.511E+02, 5.486E+01,
     1 2.031E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,14, 7, 1),I=1,6) /2.058E+03, 3.442E+02,
     1 3.135E+01, 2.260E+01, 2.020E+00, 0.000E+00/
      DATA (PH1(I,14, 7, 2),I=1,6) /3.310E+02, 6.593E+01,
     1 2.049E+01, 2.121E+01, 4.515E+00, 0.000E+00/
      DATA (PH1(I,14, 7, 3),I=1,6) /3.032E+02, 5.911E+01,
     1 2.326E+02, 6.806E+01, 3.718E+00, 2.813E-05/
      DATA (PH1(I,14, 8, 1),I=1,6) /2.001E+03, 3.759E+02,
     1 2.474E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,14, 8, 2),I=1,6) /2.872E+02, 6.086E+01,
     1 1.857E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,14, 8, 3),I=1,6) /2.465E+02, 5.786E+01,
     1 2.777E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,14, 9, 1),I=1,6) /1.946E+03, 3.799E+02,
     1 2.398E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,14, 9, 2),I=1,6) /2.468E+02, 6.237E+01,
     1 1.648E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,14, 9, 3),I=1,6) /2.051E+02, 6.007E+01,
     1 2.789E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,14,10, 1),I=1,6) /1.887E+03, 3.460E+02,
     1 3.100E+01, 1.979E+01, 2.094E+00, 0.000E+00/
      DATA (PH1(I,14,10, 2),I=1,6) /2.076E+02, 6.653E+01,
     1 1.463E+01, 2.200E+01, 4.611E+00, 0.000E+00/
      DATA (PH1(I,14,10, 3),I=1,6) /1.668E+02, 7.021E+01,
     1 2.109E+02, 8.199E+01, 3.588E+00, 6.329E-05/
      DATA (PH1(I,14,11, 1),I=1,6) /1.868E+03, 3.131E+02,
     1 3.784E+01, 2.017E+01, 2.186E+00, 0.000E+00/
      DATA (PH1(I,14,11, 2),I=1,6) /1.899E+02, 5.977E+01,
     1 1.486E+01, 2.786E+01, 4.565E+00, 0.000E+00/
      DATA (PH1(I,14,11, 3),I=1,6) /1.466E+02, 6.482E+01,
     1 2.326E+02, 1.188E+02, 3.547E+00, 8.417E-04/
      DATA (PH1(I,14,11, 4),I=1,6) /4.514E+01, 1.288E+01,
     1 6.083E+00, 1.356E+06, 3.353E+00, 0.000E+00/
      DATA (PH1(I,14,12, 1),I=1,6) /1.852E+03, 3.201E+02,
     1 3.591E+01, 2.174E+01, 2.123E+00, 0.000E+00/
      DATA (PH1(I,14,12, 2),I=1,6) /1.744E+02, 5.843E+01,
     1 1.569E+01, 2.296E+01, 4.803E+00, 0.000E+00/
      DATA (PH1(I,14,12, 3),I=1,6) /1.311E+02, 6.652E+01,
     1 2.197E+02, 1.169E+02, 3.529E+00, 8.658E-04/
      DATA (PH1(I,14,12, 4),I=1,6) /3.349E+01, 1.182E+01,
     1 8.666E+00, 3.532E+02, 4.208E+00, 0.000E+00/
      DATA (PH1(I,14,13, 1),I=1,6) /1.848E+03, 4.580E+02,
     1 1.628E+01, 7.706E+01, 1.385E+00, 0.000E+00/
      DATA (PH1(I,14,13, 2),I=1,6) /1.619E+02, 6.738E+01,
     1 1.263E+01, 3.623E+01, 4.172E+00, 0.000E+00/
      DATA (PH1(I,14,13, 3),I=1,6) /1.186E+02, 7.154E+01,
     1 1.832E+02, 3.537E+02, 3.133E+00, 2.870E-04/
      DATA (PH1(I,14,13, 4),I=1,6) /2.240E+01, 1.364E+01,
     1 8.454E+00, 7.489E+01, 4.676E+00, 0.000E+00/
      DATA (PH1(I,14,13, 5),I=1,6) /1.635E+01, 2.123E+01,
     1 6.975E+01, 4.907E+00, 9.525E+00, 3.169E-01/
      DATA (PH1(I,14,14, 1),I=1,6) /1.846E+03, 5.322E+02,
     1 1.184E+01, 2.580E+02, 1.102E+00, 0.000E+00/
      DATA (PH1(I,14,14, 2),I=1,6) /1.560E+02, 7.017E+01,
     1 1.166E+01, 4.742E+01, 3.933E+00, 0.000E+00/
      DATA (PH1(I,14,14, 3),I=1,6) /1.060E+02, 7.808E+01,
     1 1.532E+02, 5.765E+06, 2.639E+00, 2.774E-04/
      DATA (PH1(I,14,14, 4),I=1,6) /1.517E+01, 1.413E+01,
     1 1.166E+01, 2.288E+01, 5.334E+00, 0.000E+00/
      DATA (PH1(I,14,14, 5),I=1,6) /8.152E+00, 2.212E+01,
     1 1.845E+02, 3.849E+00, 9.721E+00, 2.921E-01/
      DATA (PH1(I,15, 1, 1),I=1,6) /3.070E+03, 9.701E+01,
     1 2.433E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,15, 2, 1),I=1,6) /2.817E+03, 3.381E+02,
     1 3.646E+01, 3.254E+01, 2.002E+00, 0.000E+00/
      DATA (PH1(I,15, 3, 1),I=1,6) /2.707E+03, 3.711E+02,
     1 3.220E+01, 3.257E+01, 1.869E+00, 0.000E+00/
      DATA (PH1(I,15, 3, 2),I=1,6) /6.119E+02, 6.401E+01,
     1 1.519E+01, 3.741E+01, 3.995E+00, 0.000E+00/
      DATA (PH1(I,15, 4, 1),I=1,6) /2.633E+03, 3.790E+02,
     1 3.053E+01, 3.068E+01, 1.885E+00, 0.000E+00/
      DATA (PH1(I,15, 4, 2),I=1,6) /5.604E+02, 5.712E+01,
     1 2.980E+01, 2.988E+01, 4.450E+00, 0.000E+00/
      DATA (PH1(I,15, 5, 1),I=1,6) /2.553E+03, 3.860E+02,
     1 2.914E+01, 2.656E+01, 1.945E+00, 0.000E+00/
      DATA (PH1(I,15, 5, 2),I=1,6) /5.035E+02, 7.791E+01,
     1 2.079E+01, 2.573E+01, 4.191E+00, 0.000E+00/
      DATA (PH1(I,15, 5, 3),I=1,6) /4.796E+02, 7.305E+01,
     1 7.867E+01, 7.330E+01, 3.527E+00, 1.286E-03/
      DATA (PH1(I,15, 6, 1),I=1,6) /2.477E+03, 3.824E+02,
     1 2.955E+01, 2.402E+01, 2.010E+00, 0.000E+00/
      DATA (PH1(I,15, 6, 2),I=1,6) /4.522E+02, 7.664E+01,
     1 1.966E+01, 2.403E+01, 4.317E+00, 0.000E+00/
      DATA (PH1(I,15, 6, 3),I=1,6) /4.245E+02, 6.502E+01,
     1 1.792E+02, 6.744E+01, 3.727E+00, 1.364E-03/
      DATA (PH1(I,15, 7, 1),I=1,6) /2.407E+03, 4.243E+02,
     1 2.277E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,15, 7, 2),I=1,6) /4.018E+02, 6.712E+01,
     1 1.952E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,15, 7, 3),I=1,6) /3.717E+02, 6.437E+01,
     1 2.384E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,15, 8, 1),I=1,6) /2.337E+03, 3.811E+02,
     1 2.964E+01, 2.122E+01, 2.087E+00, 0.000E+00/
      DATA (PH1(I,15, 8, 2),I=1,6) /3.552E+02, 7.217E+01,
     1 1.771E+01, 2.129E+01, 4.617E+00, 0.000E+00/
      DATA (PH1(I,15, 8, 3),I=1,6) /3.094E+02, 5.587E+01,
     1 3.696E+02, 4.738E+01, 4.244E+00, 1.094E-03/
      DATA (PH1(I,15, 9, 1),I=1,6) /2.280E+03, 4.327E+02,
     1 2.146E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,15, 9, 2),I=1,6) /3.098E+02, 7.044E+01,
     1 1.553E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,15, 9, 3),I=1,6) /2.632E+02, 6.819E+01,
     1 2.712E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,15,10, 1),I=1,6) /2.214E+03, 3.828E+02,
     1 2.932E+01, 1.991E+01, 2.123E+00, 0.000E+00/
      DATA (PH1(I,15,10, 2),I=1,6) /2.663E+02, 6.951E+01,
     1 1.560E+01, 1.845E+01, 4.938E+00, 0.000E+00/
      DATA (PH1(I,15,10, 3),I=1,6) /2.204E+02, 6.893E+01,
     1 2.980E+02, 3.497E+01, 4.307E+00, 1.008E-03/
      DATA (PH1(I,15,11, 1),I=1,6) /2.192E+03, 3.423E+02,
     1 3.671E+01, 2.010E+01, 2.232E+00, 0.000E+00/
      DATA (PH1(I,15,11, 2),I=1,6) /2.461E+02, 6.270E+01,
     1 1.535E+01, 2.486E+01, 4.787E+00, 0.000E+00/
      DATA (PH1(I,15,11, 3),I=1,6) /1.977E+02, 6.186E+01,
     1 3.368E+02, 4.700E+01, 4.256E+00, 1.118E-03/
      DATA (PH1(I,15,11, 4),I=1,6) /6.503E+01, 1.526E+00,
     1 3.688E+01, 1.111E+03, 4.890E+00, 0.000E+00/
      DATA (PH1(I,15,12, 1),I=1,6) /2.173E+03, 3.969E+02,
     1 2.640E+01, 2.945E+01, 1.901E+00, 0.000E+00/
      DATA (PH1(I,15,12, 2),I=1,6) /2.280E+02, 6.564E+01,
     1 1.491E+01, 2.333E+01, 4.787E+00, 0.000E+00/
      DATA (PH1(I,15,12, 3),I=1,6) /1.795E+02, 6.946E+01,
     1 2.661E+02, 5.760E+01, 3.967E+00, 4.482E-03/
      DATA (PH1(I,15,12, 4),I=1,6) /5.144E+01, 1.151E+01,
     1 1.153E+01, 2.303E+02, 4.447E+00, 0.000E+00/
      DATA (PH1(I,15,13, 1),I=1,6) /2.157E+03, 5.025E+02,
     1 1.577E+01, 6.163E+01, 1.467E+00, 0.000E+00/
      DATA (PH1(I,15,13, 2),I=1,6) /2.124E+02, 7.708E+01,
     1 1.206E+01, 3.388E+01, 4.192E+00, 0.000E+00/
      DATA (PH1(I,15,13, 3),I=1,6) /1.639E+02, 7.943E+01,
     1 1.949E+02, 1.392E+02, 3.391E+00, 5.039E-03/
      DATA (PH1(I,15,13, 4),I=1,6) /3.828E+01, 1.832E+01,
     1 8.260E+00, 7.058E+02, 3.682E+00, 0.000E+00/
      DATA (PH1(I,15,13, 5),I=1,6) /3.020E+01, 2.401E+01,
     1 2.851E+01, 9.795E+00, 8.210E+00, 3.290E-01/
      DATA (PH1(I,15,14, 1),I=1,6) /2.155E+03, 4.445E+02,
     1 1.990E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,15,14, 2),I=1,6) /1.987E+02, 7.299E+01,
     1 1.328E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,15,14, 3),I=1,6) /1.495E+02, 7.517E+01,
     1 2.347E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,15,14, 4),I=1,6) /2.709E+01, 1.666E+01,
     1 7.930E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,15,14, 5),I=1,6) /1.973E+01, 2.511E+01,
     1 2.564E+01, 1.800E+01, 7.200E+00, 2.600E-01/
      DATA (PH1(I,15,15, 1),I=1,6) /2.154E+03, 6.472E+02,
     1 9.167E+00, 2.562E+02, 1.063E+00, 0.000E+00/
      DATA (PH1(I,15,15, 2),I=1,6) /1.940E+02, 8.632E+01,
     1 9.931E+00, 6.594E+01, 3.617E+00, 0.000E+00/
      DATA (PH1(I,15,15, 3),I=1,6) /1.400E+02, 8.812E+01,
     1 1.512E+02, 2.230E+03, 2.795E+00, 3.422E-04/
      DATA (PH1(I,15,15, 4),I=1,6) /2.017E+01, 1.658E+01,
     1 1.125E+01, 2.613E+01, 5.205E+00, 0.000E+00/
      DATA (PH1(I,15,15, 5),I=1,6) /1.049E+01, 2.580E+01,
     1 9.925E+01, 6.712E+00, 8.516E+00, 2.765E-01/
      DATA (PH1(I,16, 1, 1),I=1,6) /3.494E+03, 1.104E+02,
     1 2.139E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,16, 2, 1),I=1,6) /3.224E+03, 4.390E+02,
     1 2.453E+01, 4.405E+01, 1.765E+00, 0.000E+00/
      DATA (PH1(I,16, 3, 1),I=1,6) /3.107E+03, 4.667E+02,
     1 2.293E+01, 4.459E+01, 1.668E+00, 0.000E+00/
      DATA (PH1(I,16, 3, 2),I=1,6) /7.072E+02, 4.414E+01,
     1 2.178E+01, 4.365E+01, 4.480E+00, 0.000E+00/
      DATA (PH1(I,16, 4, 1),I=1,6) /3.029E+03, 4.669E+02,
     1 2.233E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,16, 4, 2),I=1,6) /6.517E+02, 6.930E+01,
     1 2.571E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,16, 5, 1),I=1,6) /2.941E+03, 4.828E+02,
     1 2.087E+01, 3.742E+01, 1.720E+00, 0.000E+00/
      DATA (PH1(I,16, 5, 2),I=1,6) /5.906E+02, 6.791E+01,
     1 2.451E+01, 2.306E+01, 4.707E+00, 0.000E+00/
      DATA (PH1(I,16, 5, 3),I=1,6) /5.647E+02, 5.542E+01,
     1 1.641E+02, 6.822E+01, 4.000E+00, 4.183E-02/
      DATA (PH1(I,16, 6, 1),I=1,6) /2.859E+03, 4.758E+02,
     1 2.101E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,16, 6, 2),I=1,6) /5.346E+02, 7.360E+01,
     1 2.030E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,16, 6, 3),I=1,6) /5.048E+02, 7.213E+01,
     1 1.689E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,16, 7, 1),I=1,6) /2.782E+03, 4.802E+02,
     1 2.041E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,16, 7, 2),I=1,6) /4.804E+02, 7.555E+01,
     1 1.812E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,16, 7, 3),I=1,6) /4.471E+02, 7.326E+01,
     1 2.213E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,16, 8, 1),I=1,6) /2.705E+03, 4.937E+02,
     1 1.946E+01, 3.568E+01, 1.737E+00, 0.000E+00/
      DATA (PH1(I,16, 8, 2),I=1,6) /4.296E+02, 8.031E+01,
     1 1.659E+01, 2.170E+01, 4.615E+00, 0.000E+00/
      DATA (PH1(I,16, 8, 3),I=1,6) /3.791E+02, 6.911E+01,
     1 2.970E+02, 5.519E+01, 4.011E+00, 2.178E-02/
      DATA (PH1(I,16, 9, 1),I=1,6) /2.641E+03, 4.891E+02,
     1 1.931E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,16, 9, 2),I=1,6) /3.797E+02, 7.905E+01,
     1 1.462E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,16, 9, 3),I=1,6) /3.282E+02, 7.695E+01,
     1 2.615E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,16,10, 1),I=1,6) /2.569E+03, 4.952E+02,
     1 1.916E+01, 3.555E+01, 1.742E+00, 0.000E+00/
      DATA (PH1(I,16,10, 2),I=1,6) /3.321E+02, 8.439E+01,
     1 1.291E+01, 2.335E+01, 4.545E+00, 0.000E+00/
      DATA (PH1(I,16,10, 3),I=1,6) /2.809E+02, 7.813E+01,
     1 2.751E+02, 4.495E+01, 4.110E+00, 2.108E-02/
      DATA (PH1(I,16,11, 1),I=1,6) /2.544E+03, 4.681E+02,
     1 2.148E+01, 4.092E+01, 1.738E+00, 0.000E+00/
      DATA (PH1(I,16,11, 2),I=1,6) /3.094E+02, 6.780E+01,
     1 1.495E+01, 2.519E+01, 4.826E+00, 0.000E+00/
      DATA (PH1(I,16,11, 3),I=1,6) /2.557E+02, 5.747E+01,
     1 4.712E+02, 3.610E+01, 4.742E+00, 2.480E-02/
      DATA (PH1(I,16,11, 4),I=1,6) /8.805E+01, 1.413E+01,
     1 9.139E+00, 1.656E+03, 3.626E+00, 0.000E+00/
      DATA (PH1(I,16,12, 1),I=1,6) /2.522E+03, 4.813E+02,
     1 2.029E+01, 3.854E+01, 1.736E+00, 0.000E+00/
      DATA (PH1(I,16,12, 2),I=1,6) /2.888E+02, 7.355E+01,
     1 1.407E+01, 2.397E+01, 4.754E+00, 0.000E+00/
      DATA (PH1(I,16,12, 3),I=1,6) /2.350E+02, 7.411E+01,
     1 2.919E+02, 4.864E+01, 4.142E+00, 2.785E-02/
      DATA (PH1(I,16,12, 4),I=1,6) /7.268E+01, 1.271E+01,
     1 1.363E+01, 2.570E+02, 4.361E+00, 0.000E+00/
      DATA (PH1(I,16,13, 1),I=1,6) /2.502E+03, 5.523E+02,
     1 1.504E+01, 5.527E+01, 1.517E+00, 0.000E+00/
      DATA (PH1(I,16,13, 2),I=1,6) /2.703E+02, 8.665E+01,
     1 1.155E+01, 3.202E+01, 4.227E+00, 0.000E+00/
      DATA (PH1(I,16,13, 3),I=1,6) /2.164E+02, 8.744E+01,
     1 2.037E+02, 9.310E+01, 3.565E+00, 9.497E-03/
      DATA (PH1(I,16,13, 4),I=1,6) /5.750E+01, 2.143E+01,
     1 9.080E+00, 7.224E+02, 3.618E+00, 0.000E+00/
      DATA (PH1(I,16,13, 5),I=1,6) /4.731E+01, 3.325E+01,
     1 1.099E+01, 3.639E+01, 5.977E+00, 3.761E-01/
      DATA (PH1(I,16,14, 1),I=1,6) /2.486E+03, 6.434E+02,
     1 1.082E+01, 1.410E+02, 1.221E+00, 0.000E+00/
      DATA (PH1(I,16,14, 2),I=1,6) /2.536E+02, 9.233E+01,
     1 1.048E+01, 4.007E+01, 3.968E+00, 0.000E+00/
      DATA (PH1(I,16,14, 3),I=1,6) /1.995E+02, 9.246E+01,
     1 1.780E+02, 1.498E+02, 3.319E+00, 2.142E-02/
      DATA (PH1(I,16,14, 4),I=1,6) /4.415E+01, 2.175E+01,
     1 8.208E+00, 7.941E+02, 3.591E+00, 0.000E+00/
      DATA (PH1(I,16,14, 5),I=1,6) /3.483E+01, 2.934E+01,
     1 3.702E+01, 1.445E+01, 7.321E+00, 2.989E-01/
      DATA (PH1(I,16,15, 1),I=1,6) /2.478E+03, 7.306E+02,
     1 8.303E+00, 7.483E+02, 9.844E-01, 0.000E+00/
      DATA (PH1(I,16,15, 2),I=1,6) /2.387E+02, 9.694E+01,
     1 9.659E+00, 5.188E+01, 3.738E+00, 0.000E+00/
      DATA (PH1(I,16,15, 3),I=1,6) /1.846E+02, 9.058E+01,
     1 1.896E+02, 7.538E+01, 3.635E+00, 2.934E-01/
      DATA (PH1(I,16,15, 4),I=1,6) /3.190E+01, 2.047E+01,
     1 7.824E+00, 2.396E+02, 3.950E+00, 0.000E+00/
      DATA (PH1(I,16,15, 5),I=1,6) /2.333E+01, 2.890E+01,
     1 1.054E+02, 7.474E+00, 8.421E+00, 2.840E-01/
      DATA (PH1(I,16,16, 1),I=1,6) /2.477E+03, 8.114E+02,
     1 6.649E+00, 3.734E+03, 8.646E-01, 0.000E+00/
      DATA (PH1(I,16,16, 2),I=1,6) /2.350E+02, 1.047E+02,
     1 8.520E+00, 9.469E+01, 3.346E+00, 0.000E+00/
      DATA (PH1(I,16,16, 3),I=1,6) /1.700E+02, 9.152E+01,
     1 1.883E+02, 7.193E+01, 3.633E+00, 2.485E-01/
      DATA (PH1(I,16,16, 4),I=1,6) /2.130E+01, 1.916E+01,
     1 1.003E+01, 3.296E+01, 5.038E+00, 0.000E+00/
      DATA (PH1(I,16,16, 5),I=1,6) /1.036E+01, 2.975E+01,
     1 5.644E+01, 1.321E+01, 7.513E+00, 2.621E-01/
      DATA (PH1(I,17, 1, 1),I=1,6) /3.946E+03, 1.247E+02,
     1 1.894E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,17, 2, 1),I=1,6) /3.659E+03, 4.174E+02,
     1 3.143E+01, 3.162E+01, 2.037E+00, 0.000E+00/
      DATA (PH1(I,17, 3, 1),I=1,6) /3.534E+03, 4.432E+02,
     1 2.956E+01, 3.058E+01, 1.953E+00, 0.000E+00/
      DATA (PH1(I,17, 3, 2),I=1,6) /8.094E+02, 7.444E+01,
     1 1.356E+01, 3.797E+01, 4.103E+00, 0.000E+00/
      DATA (PH1(I,17, 4, 1),I=1,6) /3.448E+03, 4.473E+02,
     1 2.878E+01, 2.772E+01, 1.999E+00, 0.000E+00/
      DATA (PH1(I,17, 4, 2),I=1,6) /7.498E+02, 7.310E+01,
     1 2.449E+01, 3.044E+01, 4.421E+00, 0.000E+00/
      DATA (PH1(I,17, 5, 1),I=1,6) /3.356E+03, 4.633E+02,
     1 2.651E+01, 2.523E+01, 2.025E+00, 0.000E+00/
      DATA (PH1(I,17, 5, 2),I=1,6) /6.846E+02, 1.049E+02,
     1 1.622E+01, 2.683E+01, 4.067E+00, 0.000E+00/
      DATA (PH1(I,17, 5, 3),I=1,6) /6.567E+02, 3.905E+01,
     1 3.845E+02, 8.984E+01, 4.272E+00, 3.017E-01/
      DATA (PH1(I,17, 6, 1),I=1,6) /3.266E+03, 4.610E+02,
     1 2.666E+01, 2.298E+01, 2.085E+00, 0.000E+00/
      DATA (PH1(I,17, 6, 2),I=1,6) /6.249E+02, 9.701E+01,
     1 1.646E+01, 2.531E+01, 4.266E+00, 0.000E+00/
      DATA (PH1(I,17, 6, 3),I=1,6) /5.920E+02, 7.621E+01,
     1 1.831E+02, 6.513E+01, 3.852E+00, 3.987E-01/
      DATA (PH1(I,17, 7, 1),I=1,6) /3.184E+03, 5.397E+02,
     1 1.840E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,17, 7, 2),I=1,6) /5.660E+02, 8.450E+01,
     1 1.684E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,17, 7, 3),I=1,6) /5.293E+02, 8.280E+01,
     1 2.053E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,17, 8, 1),I=1,6) /3.100E+03, 4.593E+02,
     1 2.668E+01, 2.099E+01, 2.148E+00, 0.000E+00/
      DATA (PH1(I,17, 8, 2),I=1,6) /5.109E+02, 9.134E+01,
     1 1.502E+01, 2.321E+01, 4.515E+00, 0.000E+00/
      DATA (PH1(I,17, 8, 3),I=1,6) /4.556E+02, 8.470E+01,
     1 2.378E+02, 5.655E+01, 3.889E+00, 3.329E-01/
      DATA (PH1(I,17, 9, 1),I=1,6) /3.030E+03, 5.491E+02,
     1 1.746E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,17, 9, 2),I=1,6) /4.565E+02, 8.818E+01,
     1 1.375E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,17, 9, 3),I=1,6) /4.001E+02, 8.635E+01,
     1 2.505E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,17,10, 1),I=1,6) /2.951E+03, 4.603E+02,
     1 2.645E+01, 1.991E+01, 2.182E+00, 0.000E+00/
      DATA (PH1(I,17,10, 2),I=1,6) /4.048E+02, 8.397E+01,
     1 1.392E+01, 2.002E+01, 4.904E+00, 0.000E+00/
      DATA (PH1(I,17,10, 3),I=1,6) /3.483E+02, 8.434E+01,
     1 2.880E+02, 4.019E+01, 4.242E+00, 2.672E-01/
      DATA (PH1(I,17,11, 1),I=1,6) /2.923E+03, 4.414E+02,
     1 2.836E+01, 2.402E+01, 2.121E+00, 0.000E+00/
      DATA (PH1(I,17,11, 2),I=1,6) /3.797E+02, 6.106E+01,
     1 1.676E+01, 2.195E+01, 5.347E+00, 0.000E+00/
      DATA (PH1(I,17,11, 3),I=1,6) /3.207E+02, 6.471E+01,
     1 4.512E+02, 3.950E+01, 4.637E+00, 2.356E-01/
      DATA (PH1(I,17,11, 4),I=1,6) /1.142E+02, 1.534E+00,
     1 6.782E+01, 1.497E+03, 4.756E+00, 0.000E+00/
      DATA (PH1(I,17,12, 1),I=1,6) /2.898E+03, 5.111E+02,
     1 2.068E+01, 3.167E+01, 1.860E+00, 0.000E+00/
      DATA (PH1(I,17,12, 2),I=1,6) /3.566E+02, 7.233E+01,
     1 1.520E+01, 2.042E+01, 5.148E+00, 0.000E+00/
      DATA (PH1(I,17,12, 3),I=1,6) /2.974E+02, 7.682E+01,
     1 3.345E+02, 3.803E+01, 4.431E+00, 2.396E-01/
      DATA (PH1(I,17,12, 4),I=1,6) /9.703E+01, 8.603E+00,
     1 2.080E+01, 2.738E+02, 4.689E+00, 0.000E+00/
      DATA (PH1(I,17,13, 1),I=1,6) /2.875E+03, 6.130E+02,
     1 1.391E+01, 5.524E+01, 1.527E+00, 0.000E+00/
      DATA (PH1(I,17,13, 2),I=1,6) /3.354E+02, 8.204E+01,
     1 1.391E+01, 2.025E+01, 4.934E+00, 0.000E+00/
      DATA (PH1(I,17,13, 3),I=1,6) /2.760E+02, 8.644E+01,
     1 2.701E+02, 3.890E+01, 4.239E+00, 2.462E-01/
      DATA (PH1(I,17,13, 4),I=1,6) /7.997E+01, 1.825E+01,
     1 1.148E+01, 1.976E+02, 4.245E+00, 0.000E+00/
      DATA (PH1(I,17,13, 5),I=1,6) /6.782E+01, 3.654E+01,
     1 1.195E+01, 3.759E+01, 6.009E+00, 4.036E-01/
      DATA (PH1(I,17,14, 1),I=1,6) /2.851E+03, 5.520E+02,
     1 1.713E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,17,14, 2),I=1,6) /3.151E+02, 9.076E+01,
     1 1.198E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,17,14, 3),I=1,6) /2.552E+02, 9.562E+01,
     1 2.144E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,17,14, 4),I=1,6) /6.470E+01, 2.019E+01,
     1 9.833E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,17,14, 5),I=1,6) /5.347E+01, 3.505E+01,
     1 2.148E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,17,15, 1),I=1,6) /2.838E+03, 7.650E+02,
     1 8.653E+00, 2.714E+02, 1.105E+00, 0.000E+00/
      DATA (PH1(I,17,15, 2),I=1,6) /2.979E+02, 9.151E+01,
     1 1.265E+01, 2.088E+01, 4.703E+00, 0.000E+00/
      DATA (PH1(I,17,15, 3),I=1,6) /2.382E+02, 9.580E+01,
     1 2.222E+02, 4.214E+01, 4.019E+00, 2.237E-01/
      DATA (PH1(I,17,15, 4),I=1,6) /5.019E+01, 2.307E+01,
     1 8.302E+00, 2.450E+02, 3.950E+00, 0.000E+00/
      DATA (PH1(I,17,15, 5),I=1,6) /3.961E+01, 3.460E+01,
     1 3.044E+01, 3.383E+01, 6.261E+00, 2.707E-01/
      DATA (PH1(I,17,16, 1),I=1,6) /2.832E+03, 5.619E+02,
     1 1.634E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,17,16, 2),I=1,6) /2.837E+02, 9.168E+01,
     1 1.167E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,17,16, 3),I=1,6) /2.236E+02, 9.695E+01,
     1 2.052E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,17,16, 4),I=1,6) /3.686E+01, 2.191E+01,
     1 7.493E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,17,16, 5),I=1,6) /2.381E+01, 3.339E+01,
     1 4.993E+01, 1.800E+01, 7.200E+00, 2.600E-01/
      DATA (PH1(I,17,17, 1),I=1,6) /2.830E+03, 9.700E+02,
     1 5.255E+00, 1.856E+06, 7.888E-01, 0.000E+00/
      DATA (PH1(I,17,17, 2),I=1,6) /2.780E+02, 1.092E+02,
     1 1.059E+01, 2.491E+01, 4.205E+00, 0.000E+00/
      DATA (PH1(I,17,17, 3),I=1,6) /2.090E+02, 8.004E+01,
     1 3.053E+02, 3.498E+01, 4.457E+00, 2.017E-01/
      DATA (PH1(I,17,17, 4),I=1,6) /2.531E+01, 2.231E+01,
     1 6.628E+00, 1.843E+02, 4.196E+00, 0.000E+00/
      DATA (PH1(I,17,17, 5),I=1,6) /1.297E+01, 3.398E+01,
     1 4.539E+01, 2.232E+01, 6.896E+00, 2.479E-01/
      DATA (PH1(I,18, 1, 1),I=1,6) /4.426E+03, 1.399E+02,
     1 1.690E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,18, 2, 1),I=1,6) /4.121E+03, 4.468E+02,
     1 3.108E+01, 3.039E+01, 2.092E+00, 0.000E+00/
      DATA (PH1(I,18, 3, 1),I=1,6) /3.988E+03, 4.680E+02,
     1 3.003E+01, 2.854E+01, 2.037E+00, 0.000E+00/
      DATA (PH1(I,18, 3, 2),I=1,6) /9.180E+02, 5.440E+01,
     1 1.838E+01, 4.656E+01, 4.439E+00, 0.000E+00/
      DATA (PH1(I,18, 4, 1),I=1,6) /3.898E+03, 4.756E+02,
     1 2.883E+01, 2.615E+01, 2.074E+00, 0.000E+00/
      DATA (PH1(I,18, 4, 2),I=1,6) /8.548E+02, 8.162E+01,
     1 2.235E+01, 3.057E+01, 4.418E+00, 0.000E+00/
      DATA (PH1(I,18, 5, 1),I=1,6) /3.798E+03, 4.749E+02,
     1 2.874E+01, 2.235E+01, 2.171E+00, 0.000E+00/
      DATA (PH1(I,18, 5, 2),I=1,6) /7.856E+02, 1.086E+02,
     1 1.606E+01, 2.688E+01, 4.173E+00, 0.000E+00/
      DATA (PH1(I,18, 5, 3),I=1,6) /7.558E+02, 8.270E+01,
     1 1.006E+02, 7.201E+01, 3.789E+00, 5.509E-03/
      DATA (PH1(I,18, 6, 1),I=1,6) /3.702E+03, 4.731E+02,
     1 2.888E+01, 2.042E+01, 2.234E+00, 0.000E+00/
      DATA (PH1(I,18, 6, 2),I=1,6) /7.217E+02, 1.081E+02,
     1 1.515E+01, 2.572E+01, 4.251E+00, 0.000E+00/
      DATA (PH1(I,18, 6, 3),I=1,6) /6.861E+02, 5.722E+01,
     1 3.784E+02, 7.663E+01, 4.151E+00, 1.194E-02/
      DATA (PH1(I,18, 7, 1),I=1,6) /3.613E+03, 6.028E+02,
     1 1.666E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,18, 7, 2),I=1,6) /6.584E+02, 9.398E+01,
     1 1.567E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,18, 7, 3),I=1,6) /6.183E+02, 9.297E+01,
     1 1.904E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,18, 8, 1),I=1,6) /3.523E+03, 4.789E+02,
     1 2.796E+01, 1.917E+01, 2.271E+00, 0.000E+00/
      DATA (PH1(I,18, 8, 2),I=1,6) /5.992E+02, 9.884E+01,
     1 1.431E+01, 2.351E+01, 4.545E+00, 0.000E+00/
      DATA (PH1(I,18, 8, 3),I=1,6) /5.390E+02, 7.562E+01,
     1 3.466E+02, 5.329E+01, 4.195E+00, 1.392E-02/
      DATA (PH1(I,18, 9, 1),I=1,6) /3.446E+03, 6.126E+02,
     1 1.585E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,18, 9, 2),I=1,6) /5.403E+02, 9.783E+01,
     1 1.294E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,18, 9, 3),I=1,6) /4.787E+02, 9.639E+01,
     1 2.389E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,18,10, 1),I=1,6) /3.361E+03, 4.679E+02,
     1 2.931E+01, 1.744E+01, 2.362E+00, 0.000E+00/
      DATA (PH1(I,18,10, 2),I=1,6) /4.845E+02, 9.273E+01,
     1 1.306E+01, 2.086E+01, 4.862E+00, 0.000E+00/
      DATA (PH1(I,18,10, 3),I=1,6) /4.225E+02, 9.230E+01,
     1 2.855E+02, 4.508E+01, 4.165E+00, 8.883E-03/
      DATA (PH1(I,18,11, 1),I=1,6) /3.331E+03, 4.486E+02,
     1 3.143E+01, 2.008E+01, 2.315E+00, 0.000E+00/
      DATA (PH1(I,18,11, 2),I=1,6) /4.570E+02, 6.529E+01,
     1 1.602E+01, 2.336E+01, 5.323E+00, 0.000E+00/
      DATA (PH1(I,18,11, 3),I=1,6) /3.925E+02, 7.368E+01,
     1 4.198E+02, 4.419E+01, 4.492E+00, 7.712E-03/
      DATA (PH1(I,18,11, 4),I=1,6) /1.435E+02, 3.884E+00,
     1 3.295E+01, 7.082E+02, 4.645E+00, 0.000E+00/
      DATA (PH1(I,18,12, 1),I=1,6) /3.303E+03, 6.199E+02,
     1 1.559E+01, 4.503E+01, 1.664E+00, 0.000E+00/
      DATA (PH1(I,18,12, 2),I=1,6) /4.314E+02, 8.000E+01,
     1 1.424E+01, 2.136E+01, 5.094E+00, 0.000E+00/
      DATA (PH1(I,18,12, 3),I=1,6) /3.667E+02, 8.563E+01,
     1 3.204E+02, 4.230E+01, 4.329E+00, 7.258E-03/
      DATA (PH1(I,18,12, 4),I=1,6) /1.243E+02, 8.091E+00,
     1 2.596E+01, 3.212E+02, 4.685E+00, 0.000E+00/
      DATA (PH1(I,18,13, 1),I=1,6) /3.277E+03, 6.427E+02,
     1 1.444E+01, 4.208E+01, 1.659E+00, 0.000E+00/
      DATA (PH1(I,18,13, 2),I=1,6) /4.076E+02, 8.974E+01,
     1 1.313E+01, 2.094E+01, 4.919E+00, 0.000E+00/
      DATA (PH1(I,18,13, 3),I=1,6) /3.426E+02, 9.233E+01,
     1 2.786E+02, 4.220E+01, 4.227E+00, 7.408E-03/
      DATA (PH1(I,18,13, 4),I=1,6) /1.056E+02, 1.839E+01,
     1 1.304E+01, 2.002E+02, 4.303E+00, 0.000E+00/
      DATA (PH1(I,18,13, 5),I=1,6) /9.101E+01, 3.564E+01,
     1 1.418E+01, 3.945E+01, 6.194E+00, 7.398E-02/
      DATA (PH1(I,18,14, 1),I=1,6) /3.253E+03, 6.819E+02,
     1 1.269E+01, 4.958E+01, 1.562E+00, 0.000E+00/
      DATA (PH1(I,18,14, 2),I=1,6) /3.852E+02, 9.393E+01,
     1 1.265E+01, 2.082E+01, 4.848E+00, 0.000E+00/
      DATA (PH1(I,18,14, 3),I=1,6) /3.200E+02, 9.835E+01,
     1 2.474E+02, 4.284E+01, 4.125E+00, 7.283E-03/
      DATA (PH1(I,18,14, 4),I=1,6) /8.828E+01, 2.249E+01,
     1 1.048E+01, 2.022E+02, 4.127E+00, 0.000E+00/
      DATA (PH1(I,18,14, 5),I=1,6) /7.502E+01, 3.732E+01,
     1 2.476E+01, 3.876E+01, 6.142E+00, 1.882E-01/
      DATA (PH1(I,18,15, 1),I=1,6) /3.228E+03, 6.171E+02,
     1 1.548E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,18,15, 2),I=1,6) /3.642E+02, 1.010E+02,
     1 1.120E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,18,15, 3),I=1,6) /2.987E+02, 1.073E+02,
     1 2.007E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,18,15, 4),I=1,6) /7.174E+01, 2.269E+01,
     1 9.306E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,18,15, 5),I=1,6) /5.981E+01, 3.958E+01,
     1 3.058E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,18,16, 1),I=1,6) /3.216E+03, 8.408E+02,
     1 8.071E+00, 1.847E+02, 1.160E+00, 0.000E+00/
      DATA (PH1(I,18,16, 2),I=1,6) /3.455E+02, 1.031E+02,
     1 1.164E+01, 2.134E+01, 4.657E+00, 0.000E+00/
      DATA (PH1(I,18,16, 3),I=1,6) /2.801E+02, 1.025E+02,
     1 2.281E+02, 4.380E+01, 4.046E+00, 7.167E-03/
      DATA (PH1(I,18,16, 4),I=1,6) /5.637E+01, 2.612E+01,
     1 7.899E+00, 2.318E+02, 3.961E+00, 0.000E+00/
      DATA (PH1(I,18,16, 5),I=1,6) /4.074E+01, 3.862E+01,
     1 3.973E+01, 3.208E+01, 6.347E+00, 2.643E-01/
      DATA (PH1(I,18,17, 1),I=1,6) /3.208E+03, 6.260E+02,
     1 1.489E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,18,17, 2),I=1,6) /3.317E+02, 1.018E+02,
     1 1.096E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,18,17, 3),I=1,6) /2.662E+02, 1.085E+02,
     1 1.937E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,18,17, 4),I=1,6) /4.198E+01, 2.473E+01,
     1 7.104E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,18,17, 5),I=1,6) /2.763E+01, 3.780E+01,
     1 6.006E+01, 1.800E+01, 7.200E+00, 2.600E-01/
      DATA (PH1(I,18,18, 1),I=1,6) /3.203E+03, 1.135E+03,
     1 4.280E+00, 3.285E+07, 7.631E-01, 0.000E+00/
      DATA (PH1(I,18,18, 2),I=1,6) /3.260E+02, 1.302E+02,
     1 9.185E+00, 2.693E+01, 4.021E+00, 0.000E+00/
      DATA (PH1(I,18,18, 3),I=1,6) /2.492E+02, 1.647E+02,
     1 8.372E+01, 5.452E+01, 3.328E+00, 6.270E-01/
      DATA (PH1(I,18,18, 4),I=1,6) /2.892E+01, 2.525E+01,
     1 6.394E+00, 1.700E+02, 4.223E+00, 0.000E+00/
      DATA (PH1(I,18,18, 5),I=1,6) /1.576E+01, 3.854E+01,
     1 4.872E+01, 2.640E+01, 6.662E+00, 2.355E-01/
      DATA (PH1(I,19, 1, 1),I=1,6) /4.934E+03, 1.559E+02,
     1 1.517E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,19, 2, 1),I=1,6) /4.611E+03, 5.177E+02,
     1 2.589E+01, 3.282E+01, 2.018E+00, 0.000E+00/
      DATA (PH1(I,19, 3, 1),I=1,6) /4.471E+03, 5.376E+02,
     1 2.534E+01, 3.129E+01, 1.966E+00, 0.000E+00/
      DATA (PH1(I,19, 3, 2),I=1,6) /1.035E+03, 1.155E+02,
     1 8.875E+00, 3.732E+01, 3.844E+00, 0.000E+00/
      DATA (PH1(I,19, 4, 1),I=1,6) /4.375E+03, 5.442E+02,
     1 2.451E+01, 2.858E+01, 2.005E+00, 0.000E+00/
      DATA (PH1(I,19, 4, 2),I=1,6) /9.680E+02, 8.979E+01,
     1 2.067E+01, 3.141E+01, 4.401E+00, 0.000E+00/
      DATA (PH1(I,19, 5, 1),I=1,6) /4.269E+03, 5.279E+02,
     1 2.596E+01, 2.351E+01, 2.144E+00, 0.000E+00/
      DATA (PH1(I,19, 5, 2),I=1,6) /8.935E+02, 1.199E+02,
     1 1.483E+01, 2.724E+01, 4.169E+00, 0.000E+00/
      DATA (PH1(I,19, 5, 3),I=1,6) /8.611E+02, 9.596E+01,
     1 7.963E+01, 6.368E+01, 3.884E+00, 1.487E+00/
      DATA (PH1(I,19, 6, 1),I=1,6) /4.166E+03, 5.297E+02,
     1 2.562E+01, 2.213E+01, 2.182E+00, 0.000E+00/
      DATA (PH1(I,19, 6, 2),I=1,6) /8.255E+02, 1.179E+02,
     1 1.420E+01, 2.609E+01, 4.262E+00, 0.000E+00/
      DATA (PH1(I,19, 6, 3),I=1,6) /7.867E+02, 1.054E+02,
     1 1.244E+02, 6.541E+01, 3.764E+00, 9.280E-01/
      DATA (PH1(I,19, 7, 1),I=1,6) /4.070E+03, 6.694E+02,
     1 1.515E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,19, 7, 2),I=1,6) /7.578E+02, 1.040E+02,
     1 1.460E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,19, 7, 3),I=1,6) /7.147E+02, 1.038E+02,
     1 1.768E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,19, 8, 1),I=1,6) /3.974E+03, 5.204E+02,
     1 2.646E+01, 1.958E+01, 2.282E+00, 0.000E+00/
      DATA (PH1(I,19, 8, 2),I=1,6) /6.945E+02, 1.112E+02,
     1 1.306E+01, 2.458E+01, 4.472E+00, 0.000E+00/
      DATA (PH1(I,19, 8, 3),I=1,6) /6.295E+02, 1.376E+02,
     1 1.164E+02, 5.816E+01, 3.609E+00, 8.872E-01/
      DATA (PH1(I,19, 9, 1),I=1,6) /3.890E+03, 6.796E+02,
     1 1.446E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,19, 9, 2),I=1,6) /6.310E+02, 1.080E+02,
     1 1.218E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,19, 9, 3),I=1,6) /5.647E+02, 1.071E+02,
     1 2.272E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,19,10, 1),I=1,6) /3.799E+03, 5.158E+02,
     1 2.685E+01, 1.837E+01, 2.338E+00, 0.000E+00/
      DATA (PH1(I,19,10, 2),I=1,6) /5.712E+02, 1.100E+02,
     1 1.119E+01, 2.456E+01, 4.565E+00, 0.000E+00/
      DATA (PH1(I,19,10, 3),I=1,6) /5.038E+02, 1.558E+02,
     1 1.104E+02, 5.635E+01, 3.554E+00, 7.774E-01/
      DATA (PH1(I,19,11, 1),I=1,6) /3.766E+03, 5.034E+02,
     1 2.785E+01, 2.112E+01, 2.275E+00, 0.000E+00/
      DATA (PH1(I,19,11, 2),I=1,6) /5.413E+02, 8.084E+01,
     1 1.395E+01, 2.604E+01, 4.988E+00, 0.000E+00/
      DATA (PH1(I,19,11, 3),I=1,6) /4.713E+02, 1.737E+02,
     1 8.545E+01, 4.865E+01, 3.547E+00, 8.892E-01/
      DATA (PH1(I,19,11, 4),I=1,6) /1.758E+02, 1.352E+01,
     1 1.538E+01, 7.752E+02, 3.884E+00, 0.000E+00/
      DATA (PH1(I,19,12, 1),I=1,6) /3.735E+03, 6.944E+02,
     1 1.389E+01, 4.826E+01, 1.635E+00, 0.000E+00/
      DATA (PH1(I,19,12, 2),I=1,6) /5.133E+02, 1.060E+02,
     1 1.113E+01, 2.784E+01, 4.499E+00, 0.000E+00/
      DATA (PH1(I,19,12, 3),I=1,6) /4.430E+02, 1.414E+02,
     1 1.332E+02, 5.593E+01, 3.674E+00, 7.309E-01/
      DATA (PH1(I,19,12, 4),I=1,6) /1.547E+02, 1.930E+01,
     1 1.619E+01, 3.223E+02, 4.079E+00, 0.000E+00/
      DATA (PH1(I,19,13, 1),I=1,6) /3.706E+03, 7.193E+02,
     1 1.288E+01, 4.554E+01, 1.628E+00, 0.000E+00/
      DATA (PH1(I,19,13, 2),I=1,6) /4.868E+02, 1.155E+02,
     1 1.017E+01, 2.932E+01, 4.325E+00, 0.000E+00/
      DATA (PH1(I,19,13, 3),I=1,6) /4.162E+02, 1.360E+02,
     1 1.449E+02, 5.618E+01, 3.706E+00, 6.504E-01/
      DATA (PH1(I,19,13, 4),I=1,6) /1.344E+02, 3.470E+01,
     1 9.467E+00, 1.032E+03, 3.365E+00, 0.000E+00/
      DATA (PH1(I,19,13, 5),I=1,6) /1.176E+02, 4.576E+01,
     1 1.153E+01, 6.308E+01, 5.520E+00, 2.470E-01/
      DATA (PH1(I,19,14, 1),I=1,6) /3.679E+03, 8.147E+02,
     1 9.835E+00, 9.408E+01, 1.355E+00, 0.000E+00/
      DATA (PH1(I,19,14, 2),I=1,6) /4.618E+02, 1.184E+02,
     1 9.849E+00, 2.972E+01, 4.280E+00, 0.000E+00/
      DATA (PH1(I,19,14, 3),I=1,6) /3.909E+02, 2.054E+02,
     1 5.562E+01, 1.347E+04, 2.349E+00, 5.284E-01/
      DATA (PH1(I,19,14, 4),I=1,6) /1.152E+02, 3.602E+01,
     1 8.591E+00, 2.262E+03, 3.242E+00, 0.000E+00/
      DATA (PH1(I,19,14, 5),I=1,6) /9.944E+01, 5.453E+01,
     1 1.638E+01, 9.114E+01, 5.035E+00, 4.138E-01/
      DATA (PH1(I,19,15, 1),I=1,6) /3.651E+03, 6.817E+02,
     1 1.428E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,19,15, 2),I=1,6) /4.374E+02, 1.113E+02,
     1 1.061E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,19,15, 3),I=1,6) /3.665E+02, 1.188E+02,
     1 1.922E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,19,15, 4),I=1,6) /9.657E+01, 2.466E+01,
     1 9.855E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,19,15, 5),I=1,6) /8.266E+01, 4.382E+01,
     1 3.214E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,19,16, 1),I=1,6) /3.633E+03, 8.442E+02,
     1 9.086E+00, 8.929E+01, 1.343E+00, 0.000E+00/
      DATA (PH1(I,19,16, 2),I=1,6) /4.165E+02, 1.249E+02,
     1 9.136E+00, 3.259E+01, 4.130E+00, 0.000E+00/
      DATA (PH1(I,19,16, 3),I=1,6) /3.452E+02, 1.876E+02,
     1 6.770E+01, 1.658E+07, 2.363E+00, 4.546E-01/
      DATA (PH1(I,19,16, 4),I=1,6) /7.925E+01, 3.296E+01,
     1 7.619E+00, 1.323E+03, 3.423E+00, 0.000E+00/
      DATA (PH1(I,19,16, 5),I=1,6) /6.091E+01, 4.791E+01,
     1 3.344E+01, 5.426E+01, 5.626E+00, 3.140E-01/
      DATA (PH1(I,19,17, 1),I=1,6) /3.623E+03, 6.897E+02,
     1 1.383E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,19,17, 2),I=1,6) /3.990E+02, 1.121E+02,
     1 1.039E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,19,17, 3),I=1,6) /3.279E+02, 1.199E+02,
     1 1.862E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,19,17, 4),I=1,6) /6.269E+01, 2.670E+01,
     1 7.455E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,19,17, 5),I=1,6) /4.581E+01, 4.055E+01,
     1 6.784E+01, 1.800E+01, 7.200E+00, 2.600E-01/
      DATA (PH1(I,19,18, 1),I=1,6) /3.617E+03, 1.103E+03,
     1 5.173E+00, 2.564E+03, 8.917E-01, 0.000E+00/
      DATA (PH1(I,19,18, 2),I=1,6) /3.867E+02, 1.488E+02,
     1 7.085E+00, 6.561E+01, 3.449E+00, 0.000E+00/
      DATA (PH1(I,19,18, 3),I=1,6) /3.067E+02, 2.286E+02,
     1 4.199E+01, 8.731E+06, 2.239E+00, 5.560E-01/
      DATA (PH1(I,19,18, 4),I=1,6) /4.728E+01, 2.985E+01,
     1 6.625E+00, 3.001E+02, 3.891E+00, 0.000E+00/
      DATA (PH1(I,19,18, 5),I=1,6) /3.163E+01, 4.195E+01,
     1 7.392E+01, 1.697E+01, 7.208E+00, 2.531E-01/
      DATA (PH1(I,19,19, 1),I=1,6) /3.614E+03, 1.171E+03,
     1 4.540E+00, 6.165E+03, 8.392E-01, 0.000E+00/
      DATA (PH1(I,19,19, 2),I=1,6) /3.843E+02, 1.602E+02,
     1 6.389E+00, 1.044E+02, 3.159E+00, 0.000E+00/
      DATA (PH1(I,19,19, 3),I=1,6) /3.014E+02, 2.666E+02,
     1 3.107E+01, 7.187E+06, 2.067E+00, 5.274E-01/
      DATA (PH1(I,19,19, 4),I=1,6) /4.080E+01, 2.910E+01,
     1 6.377E+00, 2.229E+03, 3.587E+00, 0.000E+00/
      DATA (PH1(I,19,19, 5),I=1,6) /2.466E+01, 4.138E+01,
     1 2.614E+01, 2.143E+02, 5.631E+00, 2.437E-01/
      DATA (PH1(I,19,19, 6),I=1,6) /1.000E+00, 1.000E+00,
     1 0.000E+00, 1.000E+00, 1.000E+00, 1.000E+00/
      DATA (PH1(I,19,19, 7),I=1,6) /4.341E+00, 3.824E+00,
     1 7.363E-01, 2.410E+07, 4.427E+00, 2.049E-04/
      DATA (PH1(I,20, 1, 1),I=1,6) /5.470E+03, 1.729E+02,
     1 1.369E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,20, 2, 1),I=1,6) /5.129E+03, 6.297E+02,
     1 1.936E+01, 3.921E+01, 1.862E+00, 0.000E+00/
      DATA (PH1(I,20, 3, 1),I=1,6) /4.982E+03, 6.640E+02,
     1 1.820E+01, 3.979E+01, 1.777E+00, 0.000E+00/
      DATA (PH1(I,20, 3, 2),I=1,6) /1.157E+03, 9.395E+01,
     1 1.115E+01, 3.959E+01, 4.176E+00, 0.000E+00/
      DATA (PH1(I,20, 4, 1),I=1,6) /4.880E+03, 7.235E+02,
     1 1.486E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,20, 4, 2),I=1,6) /1.087E+03, 1.072E+02,
     1 1.788E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,20, 5, 1),I=1,6) /4.767E+03, 6.862E+02,
     1 1.665E+01, 3.443E+01, 1.823E+00, 0.000E+00/
      DATA (PH1(I,20, 5, 2),I=1,6) /1.008E+03, 1.169E+02,
     1 1.552E+01, 2.618E+01, 4.389E+00, 0.000E+00/
      DATA (PH1(I,20, 5, 3),I=1,6) /9.745E+02, 4.041E+01,
     1 5.215E+02, 1.044E+02, 4.452E+00, 1.164E-05/
      DATA (PH1(I,20, 6, 1),I=1,6) /4.659E+03, 7.342E+02,
     1 1.417E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,20, 6, 2),I=1,6) /9.357E+02, 1.122E+02,
     1 1.490E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,20, 6, 3),I=1,6) /8.946E+02, 1.144E+02,
     1 1.199E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,20, 7, 1),I=1,6) /4.555E+03, 7.396E+02,
     1 1.384E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,20, 7, 2),I=1,6) /8.642E+02, 1.145E+02,
     1 1.363E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,20, 7, 3),I=1,6) /8.177E+02, 1.152E+02,
     1 1.643E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,20, 8, 1),I=1,6) /4.453E+03, 6.989E+02,
     1 1.570E+01, 3.218E+01, 1.851E+00, 0.000E+00/
      DATA (PH1(I,20, 8, 2),I=1,6) /7.968E+02, 1.154E+02,
     1 1.287E+01, 2.461E+01, 4.561E+00, 0.000E+00/
      DATA (PH1(I,20, 8, 3),I=1,6) /7.267E+02, 8.942E+01,
     1 3.336E+02, 5.766E+01, 4.171E+00, 1.104E-05/
      DATA (PH1(I,20, 9, 1),I=1,6) /4.362E+03, 7.503E+02,
     1 1.324E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,20, 9, 2),I=1,6) /7.287E+02, 1.187E+02,
     1 1.147E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,20, 9, 3),I=1,6) /6.572E+02, 1.184E+02,
     1 2.157E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,20,10, 1),I=1,6) /4.265E+03, 7.010E+02,
     1 1.547E+01, 3.197E+01, 1.858E+00, 0.000E+00/
      DATA (PH1(I,20,10, 2),I=1,6) /6.649E+02, 9.384E+01,
     1 1.413E+01, 1.780E+01, 5.361E+00, 0.000E+00/
      DATA (PH1(I,20,10, 3),I=1,6) /5.919E+02, 7.052E+01,
     1 5.909E+02, 4.129E+01, 4.878E+00, 7.117E-06/
      DATA (PH1(I,20,11, 1),I=1,6) /4.229E+03, 6.805E+02,
     1 1.642E+01, 3.522E+01, 1.842E+00, 0.000E+00/
      DATA (PH1(I,20,11, 2),I=1,6) /6.326E+02, 8.448E+01,
     1 1.359E+01, 2.702E+01, 5.021E+00, 0.000E+00/
      DATA (PH1(I,20,11, 3),I=1,6) /5.569E+02, 6.511E+01,
     1 6.616E+02, 4.371E+01, 4.937E+00, 7.881E-06/
      DATA (PH1(I,20,11, 4),I=1,6) /2.113E+02, 1.605E+01,
     1 1.437E+01, 6.989E+02, 3.857E+00, 0.000E+00/
      DATA (PH1(I,20,12, 1),I=1,6) /4.198E+03, 7.401E+02,
     1 1.369E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,20,12, 2),I=1,6) /6.018E+02, 1.203E+02,
     1 1.052E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,20,12, 3),I=1,6) /5.270E+02, 1.293E+02,
     1 1.934E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,20,12, 4),I=1,6) /1.883E+02, 2.552E+01,
     1 1.335E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,20,13, 1),I=1,6) /4.163E+03, 7.816E+02,
     1 1.217E+01, 4.492E+01, 1.646E+00, 0.000E+00/
      DATA (PH1(I,20,13, 2),I=1,6) /5.732E+02, 1.243E+02,
     1 9.859E+00, 2.823E+01, 4.388E+00, 0.000E+00/
      DATA (PH1(I,20,13, 3),I=1,6) /4.967E+02, 1.070E+02,
     1 2.788E+02, 4.768E+01, 4.200E+00, 4.591E-06/
      DATA (PH1(I,20,13, 4),I=1,6) /1.664E+02, 3.989E+01,
     1 9.257E+00, 9.884E+02, 3.326E+00, 0.000E+00/
      DATA (PH1(I,20,13, 5),I=1,6) /1.472E+02, 5.202E+01,
     1 1.127E+01, 6.821E+01, 5.408E+00, 2.204E-01/
      DATA (PH1(I,20,14, 1),I=1,6) /4.133E+03, 8.261E+02,
     1 1.079E+01, 5.751E+01, 1.531E+00, 0.000E+00/
      DATA (PH1(I,20,14, 2),I=1,6) /5.455E+02, 1.268E+02,
     1 9.589E+00, 2.852E+01, 4.353E+00, 0.000E+00/
      DATA (PH1(I,20,14, 3),I=1,6) /4.687E+02, 1.208E+02,
     1 2.181E+02, 5.816E+01, 3.907E+00, 4.346E-06/
      DATA (PH1(I,20,14, 4),I=1,6) /1.452E+02, 4.100E+01,
     1 8.521E+00, 1.977E+03, 3.218E+00, 0.000E+00/
      DATA (PH1(I,20,14, 5),I=1,6) /1.272E+02, 6.097E+01,
     1 1.648E+01, 9.236E+01, 5.007E+00, 4.261E-01/
      DATA (PH1(I,20,15, 1),I=1,6) /4.105E+03, 8.505E+02,
     1 1.012E+01, 6.416E+01, 1.482E+00, 0.000E+00/
      DATA (PH1(I,20,15, 2),I=1,6) /5.193E+02, 1.297E+02,
     1 9.293E+00, 2.927E+01, 4.300E+00, 0.000E+00/
      DATA (PH1(I,20,15, 3),I=1,6) /4.423E+02, 1.303E+02,
     1 1.855E+02, 6.993E+01, 3.707E+00, 1.400E-05/
      DATA (PH1(I,20,15, 4),I=1,6) /1.249E+02, 3.994E+01,
     1 8.068E+00, 2.341E+03, 3.235E+00, 0.000E+00/
      DATA (PH1(I,20,15, 5),I=1,6) /1.088E+02, 5.749E+01,
     1 2.483E+01, 7.819E+01, 5.203E+00, 3.800E-01/
      DATA (PH1(I,20,16, 1),I=1,6) /4.078E+03, 8.635E+02,
     1 9.791E+00, 6.215E+01, 1.479E+00, 0.000E+00/
      DATA (PH1(I,20,16, 2),I=1,6) /4.948E+02, 1.315E+02,
     1 9.126E+00, 2.928E+01, 4.281E+00, 0.000E+00/
      DATA (PH1(I,20,16, 3),I=1,6) /4.175E+02, 1.384E+02,
     1 1.622E+02, 8.811E+01, 3.521E+00, 4.384E-04/
      DATA (PH1(I,20,16, 4),I=1,6) /1.054E+02, 3.722E+01,
     1 7.757E+00, 1.194E+03, 3.398E+00, 0.000E+00/
      DATA (PH1(I,20,16, 5),I=1,6) /8.451E+01, 5.673E+01,
     1 3.044E+01, 7.647E+01, 5.265E+00, 3.557E-01/
      DATA (PH1(I,20,17, 1),I=1,6) /4.063E+03, 9.606E+02,
     1 7.782E+00, 1.191E+02, 1.271E+00, 0.000E+00/
      DATA (PH1(I,20,17, 2),I=1,6) /4.719E+02, 1.373E+02,
     1 8.622E+00, 3.169E+01, 4.152E+00, 0.000E+00/
      DATA (PH1(I,20,17, 3),I=1,6) /3.944E+02, 1.413E+02,
     1 1.542E+02, 9.906E+01, 3.446E+00, 1.107E-03/
      DATA (PH1(I,20,17, 4),I=1,6) /8.680E+01, 3.645E+01,
     1 7.134E+00, 1.200E+03, 3.447E+00, 0.000E+00/
      DATA (PH1(I,20,17, 5),I=1,6) /6.727E+01, 5.169E+01,
     1 4.224E+01, 4.415E+01, 5.864E+00, 3.052E-01/
      DATA (PH1(I,20,18, 1),I=1,6) /4.053E+03, 1.003E+03,
     1 7.079E+00, 1.108E+02, 1.256E+00, 0.000E+00/
      DATA (PH1(I,20,18, 2),I=1,6) /4.542E+02, 1.483E+02,
     1 7.770E+00, 3.866E+01, 3.891E+00, 0.000E+00/
      DATA (PH1(I,20,18, 3),I=1,6) /3.731E+02, 1.260E+02,
     1 1.945E+02, 6.819E+01, 3.770E+00, 4.791E-04/
      DATA (PH1(I,20,18, 4),I=1,6) /6.920E+01, 3.504E+01,
     1 6.596E+00, 7.909E+02, 3.588E+00, 0.000E+00/
      DATA (PH1(I,20,18, 5),I=1,6) /5.091E+01, 4.617E+01,
     1 7.930E+01, 1.711E+01, 7.186E+00, 2.658E-01/
      DATA (PH1(I,20,19, 1),I=1,6) /4.047E+03, 7.997E+02,
     1 1.168E+01, 3.233E+01, 1.744E+00, 0.000E+00/
      DATA (PH1(I,20,19, 2),I=1,6) /4.445E+02, 1.335E+02,
     1 8.939E+00, 2.914E+01, 4.269E+00, 0.000E+00/
      DATA (PH1(I,20,19, 3),I=1,6) /3.638E+02, 1.448E+02,
     1 1.446E+02, 1.349E+02, 3.300E+00, 3.358E-04/
      DATA (PH1(I,20,19, 4),I=1,6) /6.037E+01, 3.176E+01,
     1 6.924E+00, 5.246E+02, 3.771E+00, 0.000E+00/
      DATA (PH1(I,20,19, 5),I=1,6) /4.090E+01, 4.498E+01,
     1 7.314E+01, 1.898E+01, 7.152E+00, 2.735E-01/
      DATA (PH1(I,20,19, 6),I=1,6) /1.000E+00, 1.000E+00,
     1 0.000E+00, 1.000E+00, 1.000E+00, 1.000E+00/
      DATA (PH1(I,20,19, 7),I=1,6) /1.187E+01, 4.155E+00,
     1 2.235E+00, 1.595E+04, 4.313E+00, 3.539E-01/
      DATA (PH1(I,20,20, 1),I=1,6) /4.043E+03, 6.947E+02,
     1 1.586E+01, 2.563E+01, 1.966E+00, 0.000E+00/
      DATA (PH1(I,20,20, 2),I=1,6) /4.425E+02, 1.201E+02,
     1 1.010E+01, 2.468E+01, 4.592E+00, 0.000E+00/
      DATA (PH1(I,20,20, 3),I=1,6) /3.523E+02, 1.529E+02,
     1 1.282E+02, 2.217E+02, 3.087E+00, 3.343E-03/
      DATA (PH1(I,20,20, 4),I=1,6) /4.830E+01, 3.012E+01,
     1 7.227E+00, 1.736E+02, 4.165E+00, 0.000E+00/
      DATA (PH1(I,20,20, 5),I=1,6) /3.443E+01, 4.487E+01,
     1 9.017E+01, 1.465E+01, 7.498E+00, 2.754E-01/
      DATA (PH1(I,20,20, 6),I=1,6) /1.000E+00, 1.000E+00,
     1 0.000E+00, 1.000E+00, 1.000E+00, 1.000E+00/
      DATA (PH1(I,20,20, 7),I=1,6) /6.113E+00, 7.366E+00,
     1 2.373E+00, 2.082E+02, 4.841E+00, 5.841E-04/
      DATA (PH1(I,21, 1, 1),I=1,6) /6.034E+03, 1.907E+02,
     1 1.241E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,21, 2, 1),I=1,6) /5.675E+03, 6.267E+02,
     1 2.182E+01, 3.452E+01, 1.998E+00, 0.000E+00/
      DATA (PH1(I,21, 3, 1),I=1,6) /5.520E+03, 6.439E+02,
     1 2.169E+01, 3.290E+01, 1.960E+00, 0.000E+00/
      DATA (PH1(I,21, 3, 2),I=1,6) /1.288E+03, 1.197E+02,
     1 8.820E+00, 3.718E+01, 4.054E+00, 0.000E+00/
      DATA (PH1(I,21, 4, 1),I=1,6) /5.413E+03, 6.607E+02,
     1 2.038E+01, 3.109E+01, 1.971E+00, 0.000E+00/
      DATA (PH1(I,21, 4, 2),I=1,6) /1.213E+03, 1.109E+02,
     1 1.725E+01, 3.167E+01, 4.371E+00, 0.000E+00/
      DATA (PH1(I,21, 5, 1),I=1,6) /5.294E+03, 6.684E+02,
     1 1.971E+01, 2.792E+01, 2.022E+00, 0.000E+00/
      DATA (PH1(I,21, 5, 2),I=1,6) /1.130E+03, 1.386E+02,
     1 1.325E+01, 2.785E+01, 4.218E+00, 0.000E+00/
      DATA (PH1(I,21, 5, 3),I=1,6) /1.094E+03, 1.435E+02,
     1 4.807E+01, 8.557E+01, 3.452E+00, 4.096E-02/
      DATA (PH1(I,21, 6, 1),I=1,6) /5.178E+03, 6.752E+02,
     1 1.918E+01, 2.616E+01, 2.052E+00, 0.000E+00/
      DATA (PH1(I,21, 6, 2),I=1,6) /1.054E+03, 1.367E+02,
     1 1.270E+01, 2.670E+01, 4.303E+00, 0.000E+00/
      DATA (PH1(I,21, 6, 3),I=1,6) /1.009E+03, 1.269E+02,
     1 1.155E+02, 8.134E+01, 3.618E+00, 1.480E-02/
      DATA (PH1(I,21, 7, 1),I=1,6) /5.067E+03, 8.133E+02,
     1 1.269E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,21, 7, 2),I=1,6) /9.775E+02, 1.256E+02,
     1 1.274E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,21, 7, 3),I=1,6) /9.275E+02, 1.273E+02,
     1 1.529E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,21, 8, 1),I=1,6) /4.960E+03, 6.739E+02,
     1 1.909E+01, 2.432E+01, 2.101E+00, 0.000E+00/
      DATA (PH1(I,21, 8, 2),I=1,6) /9.062E+02, 1.274E+02,
     1 1.195E+01, 2.504E+01, 4.532E+00, 0.000E+00/
      DATA (PH1(I,21, 8, 3),I=1,6) /8.308E+02, 7.899E+01,
     1 4.780E+02, 6.650E+01, 4.304E+00, 9.238E-03/
      DATA (PH1(I,21, 9, 1),I=1,6) /4.861E+03, 8.244E+02,
     1 1.216E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,21, 9, 2),I=1,6) /8.333E+02, 1.300E+02,
     1 1.081E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,21, 9, 3),I=1,6) /7.567E+02, 1.304E+02,
     1 2.045E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,21,10, 1),I=1,6) /4.759E+03, 6.781E+02,
     1 1.869E+01, 2.393E+01, 2.110E+00, 0.000E+00/
      DATA (PH1(I,21,10, 2),I=1,6) /7.657E+02, 1.156E+02,
     1 1.130E+01, 2.296E+01, 4.853E+00, 0.000E+00/
      DATA (PH1(I,21,10, 3),I=1,6) /6.874E+02, 8.313E+01,
     1 5.075E+02, 4.739E+01, 4.637E+00, 3.716E-03/
      DATA (PH1(I,21,11, 1),I=1,6) /4.720E+03, 6.682E+02,
     1 1.915E+01, 2.701E+01, 2.060E+00, 0.000E+00/
      DATA (PH1(I,21,11, 2),I=1,6) /7.309E+02, 7.807E+01,
     1 1.420E+01, 2.694E+01, 5.297E+00, 0.000E+00/
      DATA (PH1(I,21,11, 3),I=1,6) /6.505E+02, 8.091E+01,
     1 5.240E+02, 4.759E+01, 4.674E+00, 3.660E-03/
      DATA (PH1(I,21,11, 4),I=1,6) /2.498E+02, 3.414E+00,
     1 5.976E+01, 1.109E+03, 4.577E+00, 0.000E+00/
      DATA (PH1(I,21,12, 1),I=1,6) /4.684E+03, 7.696E+02,
     1 1.415E+01, 3.494E+01, 1.825E+00, 0.000E+00/
      DATA (PH1(I,21,12, 2),I=1,6) /6.982E+02, 1.019E+02,
     1 1.215E+01, 2.361E+01, 5.033E+00, 0.000E+00/
      DATA (PH1(I,21,12, 3),I=1,6) /6.173E+02, 1.007E+02,
     1 3.568E+02, 4.400E+01, 4.463E+00, 3.640E-03/
      DATA (PH1(I,21,12, 4),I=1,6) /2.251E+02, 1.104E+01,
     1 2.980E+01, 3.451E+02, 4.545E+00, 0.000E+00/
      DATA (PH1(I,21,13, 1),I=1,6) /4.649E+03, 8.283E+02,
     1 1.206E+01, 4.105E+01, 1.707E+00, 0.000E+00/
      DATA (PH1(I,21,13, 2),I=1,6) /6.666E+02, 1.131E+02,
     1 1.124E+01, 2.303E+01, 4.893E+00, 0.000E+00/
      DATA (PH1(I,21,13, 3),I=1,6) /5.852E+02, 1.050E+02,
     1 3.290E+02, 4.379E+01, 4.411E+00, 3.675E-03/
      DATA (PH1(I,21,13, 4),I=1,6) /2.015E+02, 1.499E+01,
     1 2.155E+01, 2.761E+02, 4.455E+00, 0.000E+00/
      DATA (PH1(I,21,13, 5),I=1,6) /1.800E+02, 4.763E+01,
     1 1.441E+01, 4.454E+01, 6.085E+00, 5.576E-03/
      DATA (PH1(I,21,14, 1),I=1,6) /4.615E+03, 9.428E+02,
     1 9.117E+00, 7.840E+01, 1.428E+00, 0.000E+00/
      DATA (PH1(I,21,14, 2),I=1,6) /6.362E+02, 1.182E+02,
     1 1.081E+01, 2.278E+01, 4.836E+00, 0.000E+00/
      DATA (PH1(I,21,14, 3),I=1,6) /5.544E+02, 1.106E+02,
     1 2.980E+02, 4.355E+01, 4.348E+00, 3.645E-03/
      DATA (PH1(I,21,14, 4),I=1,6) /1.784E+02, 2.031E+01,
     1 1.547E+01, 2.264E+02, 4.346E+00, 0.000E+00/
      DATA (PH1(I,21,14, 5),I=1,6) /1.581E+02, 4.941E+01,
     1 2.617E+01, 4.344E+01, 6.056E+00, 6.050E-03/
      DATA (PH1(I,21,15, 1),I=1,6) /4.582E+03, 8.245E+02,
     1 1.210E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,21,15, 2),I=1,6) /6.052E+02, 1.335E+02,
     1 9.523E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,21,15, 3),I=1,6) /5.236E+02, 1.434E+02,
     1 1.765E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,21,15, 4),I=1,6) /1.557E+02, 2.891E+01,
     1 1.052E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,21,15, 5),I=1,6) /1.380E+02, 5.322E+01,
     1 3.348E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,21,16, 1),I=1,6) /4.554E+03, 9.555E+02,
     1 8.840E+00, 7.637E+01, 1.426E+00, 0.000E+00/
      DATA (PH1(I,21,16, 2),I=1,6) /5.803E+02, 1.249E+02,
     1 1.026E+01, 2.249E+01, 4.764E+00, 0.000E+00/
      DATA (PH1(I,21,16, 3),I=1,6) /4.977E+02, 1.266E+02,
     1 2.303E+02, 4.517E+01, 4.139E+00, 3.408E-03/
      DATA (PH1(I,21,16, 4),I=1,6) /1.348E+02, 2.899E+01,
     1 9.855E+00, 1.921E+02, 4.175E+00, 0.000E+00/
      DATA (PH1(I,21,16, 5),I=1,6) /1.107E+02, 5.034E+01,
     1 4.505E+01, 3.882E+01, 6.159E+00, 6.694E-03/
      DATA (PH1(I,21,17, 1),I=1,6) /4.531E+03, 8.307E+02,
     1 1.185E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,21,17, 2),I=1,6) /5.552E+02, 1.344E+02,
     1 9.326E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,21,17, 3),I=1,6) /4.730E+02, 1.444E+02,
     1 1.721E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,21,17, 4),I=1,6) /1.141E+02, 3.036E+01,
     1 8.550E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,21,17, 5),I=1,6) /9.187E+01, 5.439E+01,
     1 4.585E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,21,18, 1),I=1,6) /4.517E+03, 9.957E+02,
     1 8.077E+00, 7.124E+01, 1.415E+00, 0.000E+00/
      DATA (PH1(I,21,18, 2),I=1,6) /5.306E+02, 1.371E+02,
     1 9.401E+00, 2.320E+01, 4.579E+00, 0.000E+00/
      DATA (PH1(I,21,18, 3),I=1,6) /4.476E+02, 1.422E+02,
     1 1.820E+02, 5.428E+01, 3.854E+00, 3.269E-03/
      DATA (PH1(I,21,18, 4),I=1,6) /9.451E+01, 3.451E+01,
     1 7.316E+00, 1.881E+02, 4.071E+00, 0.000E+00/
      DATA (PH1(I,21,18, 5),I=1,6) /7.349E+01, 5.468E+01,
     1 5.377E+01, 3.321E+01, 6.239E+00, 2.941E-01/
      DATA (PH1(I,21,19, 1),I=1,6) /4.508E+03, 1.030E+03,
     1 7.512E+00, 6.545E+01, 1.410E+00, 0.000E+00/
      DATA (PH1(I,21,19, 2),I=1,6) /5.148E+02, 1.443E+02,
     1 8.939E+00, 2.406E+01, 4.456E+00, 0.000E+00/
      DATA (PH1(I,21,19, 3),I=1,6) /4.285E+02, 1.348E+02,
     1 2.014E+02, 5.215E+01, 3.952E+00, 3.368E-03/
      DATA (PH1(I,21,19, 4),I=1,6) /7.762E+01, 3.513E+01,
     1 6.785E+00, 1.898E+02, 4.068E+00, 0.000E+00/
      DATA (PH1(I,21,19, 5),I=1,6) /5.591E+01, 5.169E+01,
     1 5.465E+01, 2.788E+01, 6.613E+00, 2.567E-01/
      DATA (PH1(I,21,19, 6),I=1,6) /2.476E+01, 2.183E+01,
     1 8.459E+01, 6.603E+01, 6.452E+00, 3.890E-01/
      DATA (PH1(I,21,20, 1),I=1,6) /4.500E+03, 7.669E+02,
     1 1.437E+01, 2.629E+01, 1.953E+00, 0.000E+00/
      DATA (PH1(I,21,20, 2),I=1,6) /5.054E+02, 1.274E+02,
     1 1.003E+01, 2.215E+01, 4.755E+00, 0.000E+00/
      DATA (PH1(I,21,20, 3),I=1,6) /4.190E+02, 1.472E+02,
     1 1.711E+02, 5.517E+01, 3.796E+00, 3.289E-03/
      DATA (PH1(I,21,20, 4),I=1,6) /6.848E+01, 3.246E+01,
     1 6.970E+00, 1.757E+02, 4.198E+00, 0.000E+00/
      DATA (PH1(I,21,20, 5),I=1,6) /4.678E+01, 4.990E+01,
     1 5.795E+01, 2.476E+01, 6.870E+00, 2.658E-01/
      DATA (PH1(I,21,20, 6),I=1,6) /1.444E+01, 1.611E+01,
     1 1.269E+02, 5.229E+01, 7.300E+00, 3.103E-01/
      DATA (PH1(I,21,20, 7),I=1,6) /1.280E+01, 3.837E+00,
     1 1.909E+00, 5.825E+02, 5.018E+00, 1.415E-02/
      DATA (PH1(I,21,21, 1),I=1,6) /4.494E+03, 7.367E+02,
     1 1.554E+01, 2.940E+01, 1.937E+00, 0.000E+00/
      DATA (PH1(I,21,21, 2),I=1,6) /5.032E+02, 1.180E+02,
     1 1.065E+01, 2.172E+01, 4.911E+00, 0.000E+00/
      DATA (PH1(I,21,21, 3),I=1,6) /4.054E+02, 1.593E+02,
     1 1.473E+02, 6.007E+01, 3.635E+00, 3.208E-03/
      DATA (PH1(I,21,21, 4),I=1,6) /5.640E+01, 3.284E+01,
     1 6.849E+00, 1.709E+02, 4.207E+00, 0.000E+00/
      DATA (PH1(I,21,21, 5),I=1,6) /3.360E+01, 4.936E+01,
     1 6.176E+01, 2.186E+01, 7.081E+00, 2.671E-01/
      DATA (PH1(I,21,21, 6),I=1,6) /8.010E+00, 9.733E+00,
     1 2.488E+02, 2.066E+01, 1.022E+01, 3.104E-01/
      DATA (PH1(I,21,21, 7),I=1,6) /7.342E+00, 8.324E+00,
     1 2.252E+00, 4.118E+02, 4.588E+00, 1.441E-04/
      DATA (PH1(I,22, 1, 1),I=1,6) /6.626E+03, 2.094E+02,
     1 1.131E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,22, 2, 1),I=1,6) /6.249E+03, 8.656E+02,
     1 1.234E+01, 5.893E+01, 1.622E+00, 0.000E+00/
      DATA (PH1(I,22, 3, 1),I=1,6) /6.087E+03, 8.854E+02,
     1 1.240E+01, 8.225E+01, 1.477E+00, 0.000E+00/
      DATA (PH1(I,22, 3, 2),I=1,6) /1.425E+03, 1.025E+02,
     1 1.038E+01, 4.131E+01, 4.265E+00, 0.000E+00/
      DATA (PH1(I,22, 4, 1),I=1,6) /5.974E+03, 7.415E+02,
     1 1.773E+01, 3.385E+01, 1.915E+00, 0.000E+00/
      DATA (PH1(I,22, 4, 2),I=1,6) /1.346E+03, 1.388E+02,
     1 1.408E+01, 3.133E+01, 4.198E+00, 0.000E+00/
      DATA (PH1(I,22, 5, 1),I=1,6) /5.848E+03, 7.545E+02,
     1 1.694E+01, 3.060E+01, 1.955E+00, 0.000E+00/
      DATA (PH1(I,22, 5, 2),I=1,6) /1.260E+03, 1.433E+02,
     1 1.296E+01, 2.842E+01, 4.278E+00, 0.000E+00/
      DATA (PH1(I,22, 5, 3),I=1,6) /1.221E+03, 1.686E+02,
     1 3.856E+01, 9.198E+01, 3.355E+00, 3.372E-05/
      DATA (PH1(I,22, 6, 1),I=1,6) /5.726E+03, 7.617E+02,
     1 1.648E+01, 2.934E+01, 1.972E+00, 0.000E+00/
      DATA (PH1(I,22, 6, 2),I=1,6) /1.179E+03, 1.419E+02,
     1 1.240E+01, 2.693E+01, 4.373E+00, 0.000E+00/
      DATA (PH1(I,22, 6, 3),I=1,6) /1.131E+03, 1.291E+02,
     1 1.268E+02, 8.686E+01, 3.655E+00, 1.321E-05/
      DATA (PH1(I,22, 7, 1),I=1,6) /5.607E+03, 8.905E+02,
     1 1.167E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,22, 7, 2),I=1,6) /1.098E+03, 1.372E+02,
     1 1.192E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,22, 7, 3),I=1,6) /1.044E+03, 1.401E+02,
     1 1.425E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,22, 8, 1),I=1,6) /5.495E+03, 7.690E+02,
     1 1.599E+01, 2.789E+01, 1.997E+00, 0.000E+00/
      DATA (PH1(I,22, 8, 2),I=1,6) /1.022E+03, 1.341E+02,
     1 1.154E+01, 2.533E+01, 4.576E+00, 0.000E+00/
      DATA (PH1(I,22, 8, 3),I=1,6) /9.419E+02, 1.084E+02,
     1 2.972E+02, 6.609E+01, 4.056E+00, 1.465E-05/
      DATA (PH1(I,22, 9, 1),I=1,6) /5.387E+03, 9.021E+02,
     1 1.121E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,22, 9, 2),I=1,6) /9.448E+02, 1.417E+02,
     1 1.019E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,22, 9, 3),I=1,6) /8.631E+02, 1.430E+02,
     1 1.937E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,22,10, 1),I=1,6) /5.281E+03, 7.712E+02,
     1 1.576E+01, 2.765E+01, 2.005E+00, 0.000E+00/
      DATA (PH1(I,22,10, 2),I=1,6) /8.734E+02, 1.250E+02,
     1 1.070E+01, 2.364E+01, 4.830E+00, 0.000E+00/
      DATA (PH1(I,22,10, 3),I=1,6) /7.878E+02, 9.852E+01,
     1 4.256E+02, 4.955E+01, 4.485E+00, 1.750E-05/
      DATA (PH1(I,22,11, 1),I=1,6) /5.239E+03, 7.394E+02,
     1 1.717E+01, 2.884E+01, 2.020E+00, 0.000E+00/
      DATA (PH1(I,22,11, 2),I=1,6) /8.363E+02, 9.240E+01,
     1 1.291E+01, 2.687E+01, 5.160E+00, 0.000E+00/
      DATA (PH1(I,22,11, 3),I=1,6) /7.511E+02, 9.404E+01,
     1 4.550E+02, 5.003E+01, 4.543E+00, 1.753E-05/
      DATA (PH1(I,22,11, 4),I=1,6) /2.915E+02, 3.359E+00,
     1 6.875E+01, 1.244E+03, 4.559E+00, 0.000E+00/
      DATA (PH1(I,22,12, 1),I=1,6) /5.200E+03, 9.318E+02,
     1 1.045E+01, 5.288E+01, 1.604E+00, 0.000E+00/
      DATA (PH1(I,22,12, 2),I=1,6) /8.013E+02, 1.259E+02,
     1 1.049E+01, 2.445E+01, 4.779E+00, 0.000E+00/
      DATA (PH1(I,22,12, 3),I=1,6) /7.156E+02, 1.122E+02,
     1 3.306E+02, 4.701E+01, 4.376E+00, 1.686E-05/
      DATA (PH1(I,22,12, 4),I=1,6) /2.650E+02, 9.778E+00,
     1 3.726E+01, 4.225E+02, 4.552E+00, 0.000E+00/
      DATA (PH1(I,22,13, 1),I=1,6) /5.162E+03, 9.469E+02,
     1 1.008E+01, 5.127E+01, 1.603E+00, 0.000E+00/
      DATA (PH1(I,22,13, 2),I=1,6) /7.670E+02, 1.254E+02,
     1 1.046E+01, 2.368E+01, 4.830E+00, 0.000E+00/
      DATA (PH1(I,22,13, 3),I=1,6) /6.807E+02, 1.133E+02,
     1 3.222E+02, 4.654E+01, 4.373E+00, 1.686E-05/
      DATA (PH1(I,22,13, 4),I=1,6) /2.397E+02, 1.712E+01,
     1 2.117E+01, 2.738E+02, 4.407E+00, 0.000E+00/
      DATA (PH1(I,22,13, 5),I=1,6) /2.159E+02, 5.265E+01,
     1 1.405E+01, 4.648E+01, 6.025E+00, 8.018E-03/
      DATA (PH1(I,22,14, 1),I=1,6) /5.126E+03, 9.454E+02,
     1 1.011E+01, 5.151E+01, 1.603E+00, 0.000E+00/
      DATA (PH1(I,22,14, 2),I=1,6) /7.342E+02, 1.245E+02,
     1 1.045E+01, 2.319E+01, 4.873E+00, 0.000E+00/
      DATA (PH1(I,22,14, 3),I=1,6) /6.473E+02, 1.160E+02,
     1 3.065E+02, 4.592E+01, 4.355E+00, 1.680E-05/
      DATA (PH1(I,22,14, 4),I=1,6) /2.148E+02, 1.736E+01,
     1 1.941E+01, 2.691E+02, 4.420E+00, 0.000E+00/
      DATA (PH1(I,22,14, 5),I=1,6) /1.921E+02, 5.337E+01,
     1 2.611E+01, 4.514E+01, 6.042E+00, 7.441E-03/
      DATA (PH1(I,22,15, 1),I=1,6) /5.089E+03, 9.027E+02,
     1 1.114E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,22,15, 2),I=1,6) /6.999E+02, 1.454E+02,
     1 9.019E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,22,15, 3),I=1,6) /6.129E+02, 1.565E+02,
     1 1.694E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,22,15, 4),I=1,6) /1.900E+02, 3.119E+01,
     1 1.069E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,22,15, 5),I=1,6) /1.704E+02, 5.839E+01,
     1 3.348E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,22,16, 1),I=1,6) /5.059E+03, 9.613E+02,
     1 9.746E+00, 5.327E+01, 1.581E+00, 0.000E+00/
      DATA (PH1(I,22,16, 2),I=1,6) /6.729E+02, 1.313E+02,
     1 9.948E+00, 2.269E+01, 4.817E+00, 0.000E+00/
      DATA (PH1(I,22,16, 3),I=1,6) /5.853E+02, 1.296E+02,
     1 2.477E+02, 4.576E+01, 4.216E+00, 1.795E-05/
      DATA (PH1(I,22,16, 4),I=1,6) /1.673E+02, 2.786E+01,
     1 1.099E+01, 1.974E+02, 4.259E+00, 0.000E+00/
      DATA (PH1(I,22,16, 5),I=1,6) /1.408E+02, 5.484E+01,
     1 4.527E+01, 3.999E+01, 6.137E+00, 9.213E-03/
      DATA (PH1(I,22,17, 1),I=1,6) /5.027E+03, 9.080E+02,
     1 1.096E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,22,17, 2),I=1,6) /6.441E+02, 1.464E+02,
     1 8.835E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,22,17, 3),I=1,6) /5.563E+02, 1.574E+02,
     1 1.656E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,22,17, 4),I=1,6) /1.445E+02, 3.258E+01,
     1 8.817E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,22,17, 5),I=1,6) /1.195E+02, 5.940E+01,
     1 4.665E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,22,18, 1),I=1,6) /5.011E+03, 1.031E+03,
     1 8.364E+00, 6.468E+01, 1.478E+00, 0.000E+00/
      DATA (PH1(I,22,18, 2),I=1,6) /6.178E+02, 1.446E+02,
     1 9.107E+00, 2.319E+01, 4.635E+00, 0.000E+00/
      DATA (PH1(I,22,18, 3),I=1,6) /5.296E+02, 1.466E+02,
     1 1.952E+02, 4.874E+01, 4.006E+00, 1.732E-05/
      DATA (PH1(I,22,18, 4),I=1,6) /1.231E+02, 3.545E+01,
     1 7.729E+00, 1.792E+02, 4.142E+00, 0.000E+00/
      DATA (PH1(I,22,18, 5),I=1,6) /9.930E+01, 5.902E+01,
     1 5.520E+01, 3.434E+01, 6.226E+00, 2.921E-01/
      DATA (PH1(I,22,19, 1),I=1,6) /4.998E+03, 1.048E+03,
     1 8.085E+00, 6.252E+01, 1.476E+00, 0.000E+00/
      DATA (PH1(I,22,19, 2),I=1,6) /5.949E+02, 1.452E+02,
     1 9.082E+00, 2.289E+01, 4.645E+00, 0.000E+00/
      DATA (PH1(I,22,19, 3),I=1,6) /5.066E+02, 1.497E+02,
     1 1.878E+02, 4.955E+01, 3.967E+00, 1.717E-05/
      DATA (PH1(I,22,19, 4),I=1,6) /1.035E+02, 3.700E+01,
     1 6.988E+00, 1.778E+02, 4.126E+00, 0.000E+00/
      DATA (PH1(I,22,19, 5),I=1,6) /7.917E+01, 5.872E+01,
     1 5.237E+01, 3.131E+01, 6.379E+00, 2.839E-01/
      DATA (PH1(I,22,19, 6),I=1,6) /4.327E+01, 2.624E+01,
     1 9.279E+01, 6.851E+01, 6.300E+00, 3.297E-01/
      DATA (PH1(I,22,20, 1),I=1,6) /4.988E+03, 1.086E+03,
     1 7.494E+00, 5.795E+01, 1.469E+00, 0.000E+00/
      DATA (PH1(I,22,20, 2),I=1,6) /5.792E+02, 1.542E+02,
     1 8.585E+00, 2.370E+01, 4.508E+00, 0.000E+00/
      DATA (PH1(I,22,20, 3),I=1,6) /4.870E+02, 1.508E+02,
     1 1.834E+02, 5.660E+01, 3.868E+00, 1.691E-05/
      DATA (PH1(I,22,20, 4),I=1,6) /8.603E+01, 3.786E+01,
     1 6.489E+00, 1.768E+02, 4.124E+00, 0.000E+00/
      DATA (PH1(I,22,20, 5),I=1,6) /6.201E+01, 5.643E+01,
     1 5.284E+01, 2.686E+01, 6.696E+00, 2.571E-01/
      DATA (PH1(I,22,20, 6),I=1,6) /2.749E+01, 2.399E+01,
     1 1.731E+02, 5.474E+01, 6.742E+00, 3.650E-01/
      DATA (PH1(I,22,21, 1),I=1,6) /4.980E+03, 6.950E+02,
     1 1.983E+01, 1.966E+01, 2.290E+00, 0.000E+00/
      DATA (PH1(I,22,21, 2),I=1,6) /5.704E+02, 1.312E+02,
     1 9.886E+00, 2.185E+01, 4.872E+00, 0.000E+00/
      DATA (PH1(I,22,21, 3),I=1,6) /4.772E+02, 1.630E+02,
     1 1.584E+02, 5.914E+01, 3.739E+00, 1.643E-05/
      DATA (PH1(I,22,21, 4),I=1,6) /7.654E+01, 3.518E+01,
     1 6.620E+00, 1.617E+02, 4.260E+00, 0.000E+00/
      DATA (PH1(I,22,21, 5),I=1,6) /5.253E+01, 5.494E+01,
     1 5.561E+01, 2.416E+01, 6.917E+00, 2.642E-01/
      DATA (PH1(I,22,21, 6),I=1,6) /1.613E+01, 1.771E+01,
     1 2.597E+02, 3.909E+01, 7.825E+00, 2.982E-01/
      DATA (PH1(I,22,21, 7),I=1,6) /1.358E+01, 5.011E+00,
     1 1.639E+00, 5.095E+02, 4.958E+00, 1.633E-03/
      DATA (PH1(I,22,22, 1),I=1,6) /4.972E+03, 6.826E+02,
     1 2.029E+01, 2.415E+01, 2.187E+00, 0.000E+00/
      DATA (PH1(I,22,22, 2),I=1,6) /5.690E+02, 1.177E+02,
     1 1.068E+01, 2.144E+01, 5.085E+00, 0.000E+00/
      DATA (PH1(I,22,22, 3),I=1,6) /4.640E+02, 1.814E+02,
     1 1.291E+02, 6.629E+01, 3.531E+00, 5.519E-05/
      DATA (PH1(I,22,22, 4),I=1,6) /6.500E+01, 3.562E+01,
     1 6.517E+00, 1.572E+02, 4.268E+00, 0.000E+00/
      DATA (PH1(I,22,22, 5),I=1,6) /4.000E+01, 5.424E+01,
     1 5.924E+01, 2.151E+01, 7.123E+00, 2.645E-01/
      DATA (PH1(I,22,22, 6),I=1,6) /9.940E+00, 1.102E+01,
     1 5.478E+02, 1.610E+01, 1.096E+01, 2.972E-01/
      DATA (PH1(I,22,22, 7),I=1,6) /6.820E+00, 9.184E+00,
     1 2.167E+00, 4.297E+02, 4.552E+00, 3.612E-03/
      DATA (PH1(I,23, 1, 1),I=1,6) /7.246E+03, 2.290E+02,
     1 1.035E+02, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,23, 2, 1),I=1,6) /6.852E+03, 7.669E+02,
     1 1.757E+01, 3.884E+01, 1.928E+00, 0.000E+00/
      DATA (PH1(I,23, 3, 1),I=1,6) /6.682E+03, 7.973E+02,
     1 1.694E+01, 3.805E+01, 1.874E+00, 0.000E+00/
      DATA (PH1(I,23, 3, 2),I=1,6) /1.570E+03, 4.820E+01,
     1 1.999E+01, 7.258E+01, 4.724E+00, 0.000E+00/
      DATA (PH1(I,23, 4, 1),I=1,6) /6.563E+03, 8.146E+02,
     1 1.606E+01, 3.587E+01, 1.887E+00, 0.000E+00/
      DATA (PH1(I,23, 4, 2),I=1,6) /1.487E+03, 1.341E+02,
     1 1.463E+01, 3.226E+01, 4.336E+00, 0.000E+00/
      DATA (PH1(I,23, 5, 1),I=1,6) /6.431E+03, 8.301E+02,
     1 1.528E+01, 3.284E+01, 1.918E+00, 0.000E+00/
      DATA (PH1(I,23, 5, 2),I=1,6) /1.396E+03, 1.597E+02,
     1 1.178E+01, 2.849E+01, 4.246E+00, 0.000E+00/
      DATA (PH1(I,23, 5, 3),I=1,6) /1.355E+03, 1.329E+02,
     1 7.374E+01, 1.017E+02, 3.596E+00, 2.681E-02/
      DATA (PH1(I,23, 6, 1),I=1,6) /6.302E+03, 8.379E+02,
     1 1.487E+01, 3.151E+01, 1.935E+00, 0.000E+00/
      DATA (PH1(I,23, 6, 2),I=1,6) /1.311E+03, 1.437E+02,
     1 1.233E+01, 2.743E+01, 4.461E+00, 0.000E+00/
      DATA (PH1(I,23, 6, 3),I=1,6) /1.260E+03, 1.761E+02,
     1 7.417E+01, 9.459E+01, 3.402E+00, 1.814E-02/
      DATA (PH1(I,23, 7, 1),I=1,6) /6.174E+03, 9.714E+02,
     1 1.078E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,23, 7, 2),I=1,6) /1.225E+03, 1.493E+02,
     1 1.118E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,23, 7, 3),I=1,6) /1.168E+03, 1.535E+02,
     1 1.330E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,23, 8, 1),I=1,6) /6.058E+03, 8.454E+02,
     1 1.444E+01, 3.020E+01, 1.955E+00, 0.000E+00/
      DATA (PH1(I,23, 8, 2),I=1,6) /1.146E+03, 1.434E+02,
     1 1.098E+01, 2.603E+01, 4.573E+00, 0.000E+00/
      DATA (PH1(I,23, 8, 3),I=1,6) /1.060E+03, 1.412E+02,
     1 1.992E+02, 6.905E+01, 3.834E+00, 1.640E-02/
      DATA (PH1(I,23, 9, 1),I=1,6) /5.941E+03, 9.834E+02,
     1 1.036E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,23, 9, 2),I=1,6) /1.063E+03, 1.540E+02,
     1 9.623E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,23, 9, 3),I=1,6) /9.758E+02, 1.562E+02,
     1 1.835E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,23,10, 1),I=1,6) /5.831E+03, 8.491E+02,
     1 1.420E+01, 2.950E+01, 1.968E+00, 0.000E+00/
      DATA (PH1(I,23,10, 2),I=1,6) /9.883E+02, 1.272E+02,
     1 1.058E+01, 2.434E+01, 4.902E+00, 0.000E+00/
      DATA (PH1(I,23,10, 3),I=1,6) /8.960E+02, 1.198E+02,
     1 3.365E+02, 5.275E+01, 4.289E+00, 1.591E-02/
      DATA (PH1(I,23,11, 1),I=1,6) /5.787E+03, 8.254E+02,
     1 1.504E+01, 3.096E+01, 1.969E+00, 0.000E+00/
      DATA (PH1(I,23,11, 2),I=1,6) /9.488E+02, 9.010E+01,
     1 1.292E+01, 2.883E+01, 5.244E+00, 0.000E+00/
      DATA (PH1(I,23,11, 3),I=1,6) /8.589E+02, 1.178E+02,
     1 3.418E+02, 5.198E+01, 4.329E+00, 1.528E-02/
      DATA (PH1(I,23,11, 4),I=1,6) /3.363E+02, 4.802E+00,
     1 5.149E+01, 9.708E+02, 4.521E+00, 0.000E+00/
      DATA (PH1(I,23,12, 1),I=1,6) /5.745E+03, 8.791E+02,
     1 1.315E+01, 3.287E+01, 1.891E+00, 0.000E+00/
      DATA (PH1(I,23,12, 2),I=1,6) /9.114E+02, 1.194E+02,
     1 1.089E+01, 2.492E+01, 4.980E+00, 0.000E+00/
      DATA (PH1(I,23,12, 3),I=1,6) /8.209E+02, 1.210E+02,
     1 3.219E+02, 4.969E+01, 4.340E+00, 1.495E-02/
      DATA (PH1(I,23,12, 4),I=1,6) /3.081E+02, 1.475E+01,
     1 2.755E+01, 3.277E+02, 4.457E+00, 0.000E+00/
      DATA (PH1(I,23,13, 1),I=1,6) /5.704E+03, 9.381E+02,
     1 1.144E+01, 3.678E+01, 1.793E+00, 0.000E+00/
      DATA (PH1(I,23,13, 2),I=1,6) /8.746E+02, 1.297E+02,
     1 1.022E+01, 2.418E+01, 4.886E+00, 0.000E+00/
      DATA (PH1(I,23,13, 3),I=1,6) /7.834E+02, 1.219E+02,
     1 3.155E+02, 4.959E+01, 4.332E+00, 1.503E-02/
      DATA (PH1(I,23,13, 4),I=1,6) /2.810E+02, 1.396E+01,
     1 2.795E+01, 3.475E+02, 4.463E+00, 0.000E+00/
      DATA (PH1(I,23,13, 5),I=1,6) /2.557E+02, 6.549E+01,
     1 1.189E+01, 4.523E+01, 5.865E+00, 5.604E-01/
      DATA (PH1(I,23,14, 1),I=1,6) /5.664E+03, 9.448E+02,
     1 1.126E+01, 3.775E+01, 1.777E+00, 0.000E+00/
      DATA (PH1(I,23,14, 2),I=1,6) /8.391E+02, 1.492E+02,
     1 8.956E+00, 2.627E+01, 4.576E+00, 0.000E+00/
      DATA (PH1(I,23,14, 3),I=1,6) /7.474E+02, 1.218E+02,
     1 3.136E+02, 4.700E+01, 4.379E+00, 1.470E-02/
      DATA (PH1(I,23,14, 4),I=1,6) /2.543E+02, 5.439E+01,
     1 8.345E+00, 9.136E+02, 3.280E+00, 0.000E+00/
      DATA (PH1(I,23,14, 5),I=1,6) /2.305E+02, 7.687E+01,
     1 1.769E+01, 7.556E+01, 5.183E+00, 4.070E-01/
      DATA (PH1(I,23,15, 1),I=1,6) /5.624E+03, 9.854E+02,
     1 1.025E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,23,15, 2),I=1,6) /8.017E+02, 1.579E+02,
     1 8.542E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,23,15, 3),I=1,6) /7.093E+02, 1.701E+02,
     1 1.628E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,23,15, 4),I=1,6) /2.274E+02, 3.357E+01,
     1 1.077E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,23,15, 5),I=1,6) /2.058E+02, 6.386E+01,
     1 3.315E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,23,16, 1),I=1,6) /5.591E+03, 9.640E+02,
     1 1.077E+01, 4.005E+01, 1.740E+00, 0.000E+00/
      DATA (PH1(I,23,16, 2),I=1,6) /7.727E+02, 1.513E+02,
     1 8.749E+00, 2.570E+01, 4.587E+00, 0.000E+00/
      DATA (PH1(I,23,16, 3),I=1,6) /6.800E+02, 1.400E+02,
     1 2.393E+02, 5.044E+01, 4.150E+00, 1.006E-02/
      DATA (PH1(I,23,16, 4),I=1,6) /2.031E+02, 5.013E+01,
     1 7.844E+00, 8.057E+02, 3.393E+00, 0.000E+00/
      DATA (PH1(I,23,16, 5),I=1,6) /1.735E+02, 6.959E+01,
     1 3.541E+01, 6.440E+01, 5.460E+00, 2.830E-01/
      DATA (PH1(I,23,17, 1),I=1,6) /5.557E+03, 9.898E+02,
     1 1.013E+01, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,23,17, 2),I=1,6) /7.401E+02, 1.589E+02,
     1 8.371E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,23,17, 3),I=1,6) /6.468E+02, 1.710E+02,
     1 1.594E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,23,17, 4),I=1,6) /1.781E+02, 3.490E+01,
     1 9.007E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,23,17, 5),I=1,6) /1.506E+02, 6.472E+01,
     1 4.696E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,23,18, 1),I=1,6) /5.531E+03, 9.782E+02,
     1 1.044E+01, 3.841E+01, 1.745E+00, 0.000E+00/
      DATA (PH1(I,23,18, 2),I=1,6) /7.123E+02, 1.580E+02,
     1 8.337E+00, 2.583E+01, 4.520E+00, 0.000E+00/
      DATA (PH1(I,23,18, 3),I=1,6) /6.188E+02, 1.692E+02,
     1 1.626E+02, 6.988E+01, 3.712E+00, 1.180E-02/
      DATA (PH1(I,23,18, 4),I=1,6) /1.549E+02, 4.722E+01,
     1 6.956E+00, 5.988E+02, 3.553E+00, 0.000E+00/
      DATA (PH1(I,23,18, 5),I=1,6) /1.281E+02, 6.736E+01,
     1 4.931E+01, 4.579E+01, 5.859E+00, 3.123E-01/
      DATA (PH1(I,23,19, 1),I=1,6) /5.516E+03, 9.875E+02,
     1 1.024E+01, 3.726E+01, 1.748E+00, 0.000E+00/
      DATA (PH1(I,23,19, 2),I=1,6) /6.859E+02, 1.633E+02,
     1 8.009E+00, 2.702E+01, 4.427E+00, 0.000E+00/
      DATA (PH1(I,23,19, 3),I=1,6) /5.923E+02, 1.721E+02,
     1 1.571E+02, 6.951E+01, 3.697E+00, 1.225E-02/
      DATA (PH1(I,23,19, 4),I=1,6) /1.329E+02, 4.591E+01,
     1 6.536E+00, 5.301E+02, 3.633E+00, 0.000E+00/
      DATA (PH1(I,23,19, 5),I=1,6) /1.059E+02, 6.053E+01,
     1 5.762E+01, 2.883E+01, 6.578E+00, 2.283E-01/
      DATA (PH1(I,23,19, 6),I=1,6) /6.528E+01, 3.190E+01,
     1 8.329E+01, 1.369E+02, 5.672E+00, 1.237E-01/
      DATA (PH1(I,23,20, 1),I=1,6) /5.504E+03, 9.965E+02,
     1 9.948E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,23,20, 2),I=1,6) /6.654E+02, 1.599E+02,
     1 8.204E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,23,20, 3),I=1,6) /5.715E+02, 1.723E+02,
     1 1.553E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,23,20, 4),I=1,6) /1.118E+02, 3.793E+01,
     1 6.694E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,23,20, 5),I=1,6) /8.555E+01, 6.751E+01,
     1 4.255E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,23,20, 6),I=1,6) /4.671E+01, 2.794E+01,
     1 2.071E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,23,21, 1),I=1,6) /5.493E+03, 1.056E+03,
     1 8.926E+00, 3.290E+01, 1.738E+00, 0.000E+00/
      DATA (PH1(I,23,21, 2),I=1,6) /6.421E+02, 1.795E+02,
     1 7.208E+00, 3.031E+01, 4.183E+00, 0.000E+00/
      DATA (PH1(I,23,21, 3),I=1,6) /5.487E+02, 1.685E+02,
     1 1.634E+02, 6.966E+01, 3.723E+00, 1.158E-02/
      DATA (PH1(I,23,21, 4),I=1,6) /9.452E+01, 4.339E+01,
     1 5.872E+00, 3.580E+02, 3.847E+00, 0.000E+00/
      DATA (PH1(I,23,21, 5),I=1,6) /6.813E+01, 6.101E+01,
     1 5.383E+01, 2.353E+01, 6.904E+00, 2.586E-01/
      DATA (PH1(I,23,21, 6),I=1,6) /2.931E+01, 1.913E+01,
     1 3.556E+03, 5.763E+00, 1.187E+01, 5.113E-01/
      DATA (PH1(I,23,22, 1),I=1,6) /5.484E+03, 1.467E+03,
     1 4.369E+00, 3.858E+02, 1.042E+00, 0.000E+00/
      DATA (PH1(I,23,22, 2),I=1,6) /6.389E+02, 1.378E+02,
     1 9.535E+00, 2.195E+01, 4.927E+00, 0.000E+00/
      DATA (PH1(I,23,22, 3),I=1,6) /5.292E+02, 1.302E+02,
     1 2.706E+02, 4.385E+01, 4.365E+00, 6.490E-05/
      DATA (PH1(I,23,22, 4),I=1,6) /7.942E+01, 4.123E+01,
     1 5.887E+00, 1.634E+02, 4.188E+00, 0.000E+00/
      DATA (PH1(I,23,22, 5),I=1,6) /5.323E+01, 6.049E+01,
     1 5.264E+01, 2.216E+01, 7.058E+00, 2.583E-01/
      DATA (PH1(I,23,22, 6),I=1,6) /1.466E+01, 1.451E+01,
     1 7.154E+02, 1.749E+01, 1.038E+01, 4.665E-01/
      DATA (PH1(I,23,23, 1),I=1,6) /5.475E+03, 6.550E+02,
     1 2.420E+01, 2.343E+01, 2.326E+00, 0.000E+00/
      DATA (PH1(I,23,23, 2),I=1,6) /6.380E+02, 9.688E+01,
     1 1.212E+01, 1.905E+01, 5.757E+00, 0.000E+00/
      DATA (PH1(I,23,23, 3),I=1,6) /5.270E+02, 2.290E+02,
     1 8.513E+01, 5.037E+02, 2.767E+00, 9.964E-04/
      DATA (PH1(I,23,23, 4),I=1,6) /7.700E+01, 3.668E+01,
     1 6.690E+00, 7.759E+01, 4.706E+00, 0.000E+00/
      DATA (PH1(I,23,23, 5),I=1,6) /4.700E+01, 5.937E+01,
     1 7.599E+01, 1.488E+01, 7.586E+00, 2.690E-01/
      DATA (PH1(I,23,23, 6),I=1,6) /1.200E+01, 1.146E+01,
     1 3.134E+03, 7.037E+00, 1.377E+01, 3.417E-01/
      DATA (PH1(I,23,23, 7),I=1,6) /6.740E+00, 1.002E+01,
     1 2.059E+00, 3.914E+02, 4.565E+00, 5.057E-04/
      DATA (PH1(I,24, 1, 1),I=1,6) /7.895E+03, 2.495E+02,
     1 9.505E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,24, 2, 1),I=1,6) /7.482E+03, 9.695E+02,
     1 1.185E+01, 5.351E+01, 1.696E+00, 0.000E+00/
      DATA (PH1(I,24, 3, 1),I=1,6) /7.306E+03, 9.824E+02,
     1 1.199E+01, 5.185E+01, 1.668E+00, 0.000E+00/
      DATA (PH1(I,24, 3, 2),I=1,6) /1.721E+03, 1.480E+02,
     1 7.312E+00, 3.880E+01, 4.083E+00, 0.000E+00/
      DATA (PH1(I,24, 4, 1),I=1,6) /7.181E+03, 9.952E+02,
     1 1.156E+01, 4.841E+01, 1.687E+00, 0.000E+00/
      DATA (PH1(I,24, 4, 2),I=1,6) /1.634E+03, 1.488E+02,
     1 1.333E+01, 3.225E+01, 4.310E+00, 0.000E+00/
      DATA (PH1(I,24, 5, 1),I=1,6) /7.042E+03, 1.003E+03,
     1 1.125E+01, 4.457E+01, 1.717E+00, 0.000E+00/
      DATA (PH1(I,24, 5, 2),I=1,6) /1.539E+03, 1.748E+02,
     1 1.089E+01, 2.892E+01, 4.221E+00, 0.000E+00/
      DATA (PH1(I,24, 5, 3),I=1,6) /1.497E+03, 2.216E+02,
     1 2.692E+01, 1.091E+02, 3.193E+00, 3.540E-03/
      DATA (PH1(I,24, 6, 1),I=1,6) /6.907E+03, 1.004E+03,
     1 1.114E+01, 4.238E+01, 1.740E+00, 0.000E+00/
      DATA (PH1(I,24, 6, 2),I=1,6) /1.450E+03, 1.749E+02,
     1 1.035E+01, 2.790E+01, 4.278E+00, 0.000E+00/
      DATA (PH1(I,24, 6, 3),I=1,6) /1.396E+03, 2.006E+02,
     1 6.290E+01, 1.016E+02, 3.330E+00, 1.634E-03/
      DATA (PH1(I,24, 7, 1),I=1,6) /6.769E+03, 1.056E+03,
     1 9.975E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24, 7, 2),I=1,6) /1.359E+03, 1.619E+02,
     1 1.050E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24, 7, 3),I=1,6) /1.299E+03, 1.675E+02,
     1 1.244E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,24, 8, 1),I=1,6) /6.650E+03, 1.010E+03,
     1 1.086E+01, 4.093E+01, 1.755E+00, 0.000E+00/
      DATA (PH1(I,24, 8, 2),I=1,6) /1.276E+03, 1.514E+02,
     1 1.053E+01, 2.659E+01, 4.591E+00, 0.000E+00/
      DATA (PH1(I,24, 8, 3),I=1,6) /1.185E+03, 1.297E+02,
     1 2.646E+02, 7.449E+01, 3.960E+00, 6.187E-03/
      DATA (PH1(I,24, 9, 1),I=1,6) /6.523E+03, 1.068E+03,
     1 9.610E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24, 9, 2),I=1,6) /1.189E+03, 1.668E+02,
     1 9.095E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24, 9, 3),I=1,6) /1.097E+03, 1.701E+02,
     1 1.738E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,24,10, 1),I=1,6) /6.409E+03, 1.017E+03,
     1 1.063E+01, 4.093E+01, 1.756E+00, 0.000E+00/
      DATA (PH1(I,24,10, 2),I=1,6) /1.110E+03, 1.388E+02,
     1 9.935E+00, 2.496E+01, 4.858E+00, 0.000E+00/
      DATA (PH1(I,24,10, 3),I=1,6) /1.011E+03, 1.397E+02,
     1 2.826E+02, 5.599E+01, 4.153E+00, 5.775E-03/
      DATA (PH1(I,24,11, 1),I=1,6) /6.362E+03, 9.927E+02,
     1 1.118E+01, 4.335E+01, 1.752E+00, 0.000E+00/
      DATA (PH1(I,24,11, 2),I=1,6) /1.068E+03, 9.748E+01,
     1 1.227E+01, 2.968E+01, 5.204E+00, 0.000E+00/
      DATA (PH1(I,24,11, 3),I=1,6) /9.738E+02, 1.312E+02,
     1 3.124E+02, 5.494E+01, 4.251E+00, 6.239E-03/
      DATA (PH1(I,24,11, 4),I=1,6) /3.842E+02, 5.856E+00,
     1 4.533E+01, 8.787E+02, 4.495E+00, 0.000E+00/
      DATA (PH1(I,24,12, 1),I=1,6) /6.318E+03, 1.103E+03,
     1 8.914E+00, 5.657E+01, 1.588E+00, 0.000E+00/
      DATA (PH1(I,24,12, 2),I=1,6) /1.029E+03, 1.295E+02,
     1 1.028E+01, 2.555E+01, 4.945E+00, 0.000E+00/
      DATA (PH1(I,24,12, 3),I=1,6) /9.334E+02, 1.427E+02,
     1 2.652E+02, 5.321E+01, 4.178E+00, 5.905E-03/
      DATA (PH1(I,24,12, 4),I=1,6) /3.548E+02, 1.135E+01,
     1 3.894E+01, 4.475E+02, 4.496E+00, 0.000E+00/
      DATA (PH1(I,24,13, 1),I=1,6) /6.274E+03, 1.115E+03,
     1 8.708E+00, 5.540E+01, 1.588E+00, 0.000E+00/
      DATA (PH1(I,24,13, 2),I=1,6) /9.893E+02, 1.454E+02,
     1 9.421E+00, 2.472E+01, 4.803E+00, 0.000E+00/
      DATA (PH1(I,24,13, 3),I=1,6) /8.933E+02, 1.424E+02,
     1 2.644E+02, 5.317E+01, 4.181E+00, 5.962E-03/
      DATA (PH1(I,24,13, 4),I=1,6) /3.255E+02, 1.164E+01,
     1 3.651E+01, 4.413E+02, 4.490E+00, 0.000E+00/
      DATA (PH1(I,24,13, 5),I=1,6) /2.981E+02, 6.675E+01,
     1 1.246E+01, 4.963E+01, 5.860E+00, 3.057E-01/
      DATA (PH1(I,24,14, 1),I=1,6) /6.231E+03, 1.111E+03,
     1 8.763E+00, 5.576E+01, 1.589E+00, 0.000E+00/
      DATA (PH1(I,24,14, 2),I=1,6) /9.512E+02, 1.447E+02,
     1 9.388E+00, 2.428E+01, 4.837E+00, 0.000E+00/
      DATA (PH1(I,24,14, 3),I=1,6) /8.546E+02, 1.426E+02,
     1 2.618E+02, 5.270E+01, 4.188E+00, 5.989E-03/
      DATA (PH1(I,24,14, 4),I=1,6) /2.969E+02, 2.539E+01,
     1 1.632E+01, 2.419E+02, 4.286E+00, 0.000E+00/
      DATA (PH1(I,24,14, 5),I=1,6) /2.708E+02, 6.544E+01,
     1 2.411E+01, 4.849E+01, 5.923E+00, 2.753E-01/
      DATA (PH1(I,24,15, 1),I=1,6) /6.187E+03, 1.073E+03,
     1 9.448E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24,15, 2),I=1,6) /9.108E+02, 1.710E+02,
     1 8.093E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24,15, 3),I=1,6) /8.130E+02, 1.843E+02,
     1 1.566E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,24,15, 4),I=1,6) /2.680E+02, 3.606E+01,
     1 1.078E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,24,15, 5),I=1,6) /2.444E+02, 6.965E+01,
     1 3.257E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,24,16, 1),I=1,6) /6.152E+03, 1.110E+03,
     1 8.774E+00, 5.647E+01, 1.586E+00, 0.000E+00/
      DATA (PH1(I,24,16, 2),I=1,6) /8.797E+02, 1.476E+02,
     1 9.146E+00, 2.358E+01, 4.847E+00, 0.000E+00/
      DATA (PH1(I,24,16, 3),I=1,6) /7.819E+02, 1.427E+02,
     1 2.573E+02, 5.042E+01, 4.228E+00, 5.960E-03/
      DATA (PH1(I,24,16, 4),I=1,6) /2.419E+02, 2.465E+01,
     1 1.419E+01, 2.341E+02, 4.369E+00, 0.000E+00/
      DATA (PH1(I,24,16, 5),I=1,6) /2.093E+02, 6.644E+01,
     1 4.311E+01, 4.198E+01, 6.057E+00, 3.135E-01/
      DATA (PH1(I,24,17, 1),I=1,6) /6.114E+03, 1.076E+03,
     1 9.366E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24,17, 2),I=1,6) /8.432E+02, 1.720E+02,
     1 7.934E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24,17, 3),I=1,6) /7.444E+02, 1.850E+02,
     1 1.537E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,24,17, 4),I=1,6) /2.148E+02, 3.733E+01,
     1 9.131E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,24,17, 5),I=1,6) /1.847E+02, 7.035E+01,
     1 4.687E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,24,18, 1),I=1,6) /6.080E+03, 1.122E+03,
     1 8.562E+00, 5.603E+01, 1.582E+00, 0.000E+00/
      DATA (PH1(I,24,18, 2),I=1,6) /8.140E+02, 1.553E+02,
     1 8.726E+00, 2.320E+01, 4.789E+00, 0.000E+00/
      DATA (PH1(I,24,18, 3),I=1,6) /7.153E+02, 1.571E+02,
     1 2.130E+02, 5.017E+01, 4.114E+00, 5.706E-03/
      DATA (PH1(I,24,18, 4),I=1,6) /1.900E+02, 3.531E+01,
     1 8.881E+00, 1.815E+02, 4.264E+00, 0.000E+00/
      DATA (PH1(I,24,18, 5),I=1,6) /1.602E+02, 6.754E+01,
     1 5.675E+01, 3.729E+01, 6.191E+00, 2.358E-01/
      DATA (PH1(I,24,19, 1),I=1,6) /6.062E+03, 1.128E+03,
     1 8.464E+00, 5.569E+01, 1.580E+00, 0.000E+00/
      DATA (PH1(I,24,19, 2),I=1,6) /7.842E+02, 1.551E+02,
     1 8.723E+00, 2.287E+01, 4.810E+00, 0.000E+00/
      DATA (PH1(I,24,19, 3),I=1,6) /6.855E+02, 1.640E+02,
     1 1.963E+02, 5.019E+01, 4.060E+00, 5.648E-03/
      DATA (PH1(I,24,19, 4),I=1,6) /1.656E+02, 3.830E+01,
     1 7.706E+00, 1.719E+02, 4.245E+00, 0.000E+00/
      DATA (PH1(I,24,19, 5),I=1,6) /1.359E+02, 6.693E+01,
     1 5.368E+01, 3.420E+01, 6.344E+00, 2.242E-01/
      DATA (PH1(I,24,19, 6),I=1,6) /9.064E+01, 3.098E+01,
     1 1.201E+02, 6.811E+01, 6.392E+00, 3.113E-03/
      DATA (PH1(I,24,20, 1),I=1,6) /6.046E+03, 1.081E+03,
     1 9.244E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24,20, 2),I=1,6) /7.598E+02, 1.731E+02,
     1 7.776E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24,20, 3),I=1,6) /6.602E+02, 1.862E+02,
     1 1.500E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,24,20, 4),I=1,6) /1.427E+02, 4.027E+01,
     1 6.892E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,24,20, 5),I=1,6) /1.134E+02, 7.290E+01,
     1 4.327E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,24,20, 6),I=1,6) /6.946E+01, 3.105E+01,
     1 2.269E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,24,21, 1),I=1,6) /6.033E+03, 1.224E+03,
     1 7.106E+00, 7.179E+01, 1.455E+00, 0.000E+00/
      DATA (PH1(I,24,21, 2),I=1,6) /7.331E+02, 1.640E+02,
     1 8.333E+00, 2.286E+01, 4.721E+00, 0.000E+00/
      DATA (PH1(I,24,21, 3),I=1,6) /6.344E+02, 1.721E+02,
     1 1.795E+02, 5.067E+01, 3.993E+00, 5.538E-03/
      DATA (PH1(I,24,21, 4),I=1,6) /1.219E+02, 4.256E+01,
     1 6.252E+00, 1.605E+02, 4.224E+00, 0.000E+00/
      DATA (PH1(I,24,21, 5),I=1,6) /9.275E+01, 6.776E+01,
     1 4.856E+01, 2.817E+01, 6.628E+00, 2.673E-01/
      DATA (PH1(I,24,21, 6),I=1,6) /4.916E+01, 3.483E+01,
     1 2.224E+02, 5.448E+01, 6.484E+00, 3.436E-01/
      DATA (PH1(I,24,22, 1),I=1,6) /6.021E+03, 1.085E+03,
     1 9.165E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,24,22, 2),I=1,6) /7.161E+02, 1.735E+02,
     1 7.720E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,24,22, 3),I=1,6) /6.164E+02, 1.871E+02,
     1 1.480E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,24,22, 4),I=1,6) /1.032E+02, 4.293E+01,
     1 5.797E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,24,22, 5),I=1,6) /7.436E+01, 6.450E+01,
     1 5.388E+01, 2.100E+01, 7.200E+00, 2.600E-01/
      DATA (PH1(I,24,22, 6),I=1,6) /3.096E+01, 9.947E+00,
     1 1.924E+03, 9.000E+00, 1.450E+01, 3.500E-01/
      DATA (PH1(I,24,23, 1),I=1,6) /6.009E+03, 1.497E+03,
     1 4.593E+00, 1.465E+02, 1.189E+00, 0.000E+00/
      DATA (PH1(I,24,23, 2),I=1,6) /7.074E+02, 2.013E+02,
     1 6.869E+00, 2.572E+01, 4.266E+00, 0.000E+00/
      DATA (PH1(I,24,23, 3),I=1,6) /5.970E+02, 1.477E+02,
     1 2.382E+02, 4.633E+01, 4.261E+00, 6.156E-03/
      DATA (PH1(I,24,23, 4),I=1,6) /8.745E+01, 4.462E+01,
     1 5.566E+00, 1.542E+02, 4.230E+00, 0.000E+00/
      DATA (PH1(I,24,23, 5),I=1,6) /5.879E+01, 6.599E+01,
     1 5.021E+01, 2.175E+01, 7.102E+00, 2.576E-01/
      DATA (PH1(I,24,23, 6),I=1,6) /1.650E+01, 1.521E+01,
     1 1.031E+03, 1.415E+01, 1.123E+01, 4.827E-01/
      DATA (PH1(I,24,24, 1),I=1,6) /5.996E+03, 1.035E+03,
     1 1.025E+01, 3.343E+01, 1.822E+00, 0.000E+00/
      DATA (PH1(I,24,24, 2),I=1,6) /7.030E+02, 1.588E+02,
     1 8.555E+00, 2.258E+01, 4.789E+00, 0.000E+00/
      DATA (PH1(I,24,24, 3),I=1,6) /5.850E+02, 1.865E+02,
     1 1.540E+02, 5.560E+01, 3.823E+00, 5.785E-03/
      DATA (PH1(I,24,24, 4),I=1,6) /7.900E+01, 4.234E+01,
     1 5.602E+00, 1.356E+02, 4.374E+00, 0.000E+00/
      DATA (PH1(I,24,24, 5),I=1,6) /4.900E+01, 6.567E+01,
     1 5.313E+01, 1.981E+01, 7.258E+00, 2.601E-01/
      DATA (PH1(I,24,24, 6),I=1,6) /8.660E+00, 7.244E+00,
     1 1.485E+03, 9.671E+00, 1.575E+01, 7.760E-01/
      DATA (PH1(I,24,24, 7),I=1,6) /6.767E+00, 9.636E+00,
     1 6.532E-01, 5.232E+02, 4.641E+00, 9.332E-05/
      DATA (PH1(I,25, 1, 1),I=1,6) /8.572E+03, 2.709E+02,
     1 8.759E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,25, 2, 1),I=1,6) /8.141E+03, 9.661E+02,
     1 1.309E+01, 4.829E+01, 1.794E+00, 0.000E+00/
      DATA (PH1(I,25, 3, 1),I=1,6) /7.957E+03, 9.997E+02,
     1 1.266E+01, 4.727E+01, 1.749E+00, 0.000E+00/
      DATA (PH1(I,25, 3, 2),I=1,6) /1.880E+03, 2.004E+02,
     1 5.352E+00, 3.934E+01, 3.800E+00, 0.000E+00/
      DATA (PH1(I,25, 4, 1),I=1,6) /7.827E+03, 9.916E+02,
     1 1.277E+01, 4.256E+01, 1.799E+00, 0.000E+00/
      DATA (PH1(I,25, 4, 2),I=1,6) /1.788E+03, 1.631E+02,
     1 1.229E+01, 3.227E+01, 4.294E+00, 0.000E+00/
      DATA (PH1(I,25, 5, 1),I=1,6) /7.682E+03, 9.946E+02,
     1 1.255E+01, 3.824E+01, 1.847E+00, 0.000E+00/
      DATA (PH1(I,25, 5, 2),I=1,6) /1.689E+03, 1.952E+02,
     1 9.881E+00, 2.837E+01, 4.199E+00, 0.000E+00/
      DATA (PH1(I,25, 5, 3),I=1,6) /1.644E+03, 2.576E+02,
     1 2.159E+01, 1.234E+02, 3.086E+00, 1.723E-02/
      DATA (PH1(I,25, 6, 1),I=1,6) /7.539E+03, 9.973E+02,
     1 1.238E+01, 3.663E+01, 1.868E+00, 0.000E+00/
      DATA (PH1(I,25, 6, 2),I=1,6) /1.596E+03, 1.879E+02,
     1 9.744E+00, 2.825E+01, 4.277E+00, 0.000E+00/
      DATA (PH1(I,25, 6, 3),I=1,6) /1.539E+03, 2.403E+02,
     1 4.744E+01, 1.147E+02, 3.190E+00, 1.237E-02/
      DATA (PH1(I,25, 7, 1),I=1,6) /7.391E+03, 1.144E+03,
     1 9.260E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25, 7, 2),I=1,6) /1.500E+03, 1.751E+02,
     1 9.871E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25, 7, 3),I=1,6) /1.437E+03, 1.822E+02,
     1 1.165E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,25, 8, 1),I=1,6) /7.270E+03, 1.005E+03,
     1 1.205E+01, 3.464E+01, 1.893E+00, 0.000E+00/
      DATA (PH1(I,25, 8, 2),I=1,6) /1.414E+03, 1.452E+02,
     1 1.087E+01, 2.737E+01, 4.741E+00, 0.000E+00/
      DATA (PH1(I,25, 8, 3),I=1,6) /1.317E+03, 1.889E+02,
     1 1.381E+02, 7.915E+01, 3.619E+00, 1.713E-02/
      DATA (PH1(I,25, 9, 1),I=1,6) /7.132E+03, 1.157E+03,
     1 8.934E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25, 9, 2),I=1,6) /1.321E+03, 1.802E+02,
     1 8.605E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25, 9, 3),I=1,6) /1.224E+03, 1.847E+02,
     1 1.647E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,25,10, 1),I=1,6) /7.016E+03, 1.007E+03,
     1 1.191E+01, 3.374E+01, 1.907E+00, 0.000E+00/
      DATA (PH1(I,25,10, 2),I=1,6) /1.239E+03, 1.468E+02,
     1 9.533E+00, 2.564E+01, 4.859E+00, 0.000E+00/
      DATA (PH1(I,25,10, 3),I=1,6) /1.133E+03, 1.588E+02,
     1 2.466E+02, 5.979E+01, 4.044E+00, 1.415E-02/
      DATA (PH1(I,25,11, 1),I=1,6) /6.966E+03, 9.952E+02,
     1 1.220E+01, 3.635E+01, 1.885E+00, 0.000E+00/
      DATA (PH1(I,25,11, 2),I=1,6) /1.195E+03, 1.084E+02,
     1 1.148E+01, 3.003E+01, 5.139E+00, 0.000E+00/
      DATA (PH1(I,25,11, 3),I=1,6) /1.096E+03, 1.576E+02,
     1 2.464E+02, 5.814E+01, 4.081E+00, 1.417E-02/
      DATA (PH1(I,25,11, 4),I=1,6) /4.352E+02, 5.053E+00,
     1 5.856E+01, 1.103E+03, 4.490E+00, 0.000E+00/
      DATA (PH1(I,25,12, 1),I=1,6) /6.919E+03, 1.016E+03,
     1 1.168E+01, 3.428E+01, 1.894E+00, 0.000E+00/
      DATA (PH1(I,25,12, 2),I=1,6) /1.153E+03, 1.412E+02,
     1 9.670E+00, 2.617E+01, 4.900E+00, 0.000E+00/
      DATA (PH1(I,25,12, 3),I=1,6) /1.053E+03, 1.572E+02,
     1 2.455E+02, 5.660E+01, 4.109E+00, 1.411E-02/
      DATA (PH1(I,25,12, 4),I=1,6) /4.030E+02, 1.345E+01,
     1 3.535E+01, 4.180E+02, 4.462E+00, 0.000E+00/
      DATA (PH1(I,25,13, 1),I=1,6) /6.872E+03, 1.025E+03,
     1 1.145E+01, 3.313E+01, 1.902E+00, 0.000E+00/
      DATA (PH1(I,25,13, 2),I=1,6) /1.111E+03, 1.437E+02,
     1 9.475E+00, 2.535E+01, 4.916E+00, 0.000E+00/
      DATA (PH1(I,25,13, 3),I=1,6) /1.010E+03, 1.534E+02,
     1 2.551E+02, 5.638E+01, 4.140E+00, 1.436E-02/
      DATA (PH1(I,25,13, 4),I=1,6) /3.731E+02, 2.286E+01,
     1 2.034E+01, 2.772E+02, 4.326E+00, 0.000E+00/
      DATA (PH1(I,25,13, 5),I=1,6) /3.436E+02, 7.288E+01,
     1 1.205E+01, 5.254E+01, 5.787E+00, 1.633E-04/
      DATA (PH1(I,25,14, 1),I=1,6) /6.826E+03, 1.018E+03,
     1 1.161E+01, 3.333E+01, 1.905E+00, 0.000E+00/
      DATA (PH1(I,25,14, 2),I=1,6) /1.070E+03, 1.522E+02,
     1 9.043E+00, 2.474E+01, 4.856E+00, 0.000E+00/
      DATA (PH1(I,25,14, 3),I=1,6) /9.690E+02, 1.541E+02,
     1 2.509E+02, 5.596E+01, 4.141E+00, 1.435E-02/
      DATA (PH1(I,25,14, 4),I=1,6) /3.427E+02, 2.273E+01,
     1 1.932E+01, 2.771E+02, 4.340E+00, 0.000E+00/
      DATA (PH1(I,25,14, 5),I=1,6) /3.144E+02, 7.442E+01,
     1 2.225E+01, 4.974E+01, 5.829E+00, 4.047E-01/
      DATA (PH1(I,25,15, 1),I=1,6) /6.778E+03, 1.164E+03,
     1 8.718E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25,15, 2),I=1,6) /1.027E+03, 1.846E+02,
     1 7.670E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25,15, 3),I=1,6) /9.238E+02, 1.990E+02,
     1 1.508E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,25,15, 4),I=1,6) /3.118E+02, 3.864E+01,
     1 1.074E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,25,15, 5),I=1,6) /2.860E+02, 7.574E+01,
     1 3.181E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,25,16, 1),I=1,6) /6.741E+03, 1.013E+03,
     1 1.172E+01, 3.345E+01, 1.909E+00, 0.000E+00/
      DATA (PH1(I,25,16, 2),I=1,6) /9.939E+02, 1.554E+02,
     1 8.799E+00, 2.403E+01, 4.866E+00, 0.000E+00/
      DATA (PH1(I,25,16, 3),I=1,6) /8.910E+02, 1.555E+02,
     1 2.430E+02, 5.379E+01, 4.166E+00, 1.402E-02/
      DATA (PH1(I,25,16, 4),I=1,6) /2.840E+02, 2.364E+01,
     1 1.585E+01, 2.573E+02, 4.393E+00, 0.000E+00/
      DATA (PH1(I,25,16, 5),I=1,6) /2.483E+02, 7.277E+01,
     1 4.169E+01, 4.324E+01, 6.012E+00, 3.835E-01/
      DATA (PH1(I,25,17, 1),I=1,6) /6.698E+03, 1.167E+03,
     1 8.667E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25,17, 2),I=1,6) /9.536E+02, 1.856E+02,
     1 7.522E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25,17, 3),I=1,6) /8.492E+02, 1.997E+02,
     1 1.482E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,25,17, 4),I=1,6) /2.547E+02, 3.986E+01,
     1 9.199E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,25,17, 5),I=1,6) /2.218E+02, 7.629E+01,
     1 4.645E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,25,18, 1),I=1,6) /6.664E+03, 1.021E+03,
     1 1.152E+01, 3.256E+01, 1.914E+00, 0.000E+00/
      DATA (PH1(I,25,18, 2),I=1,6) /9.230E+02, 1.679E+02,
     1 8.168E+00, 2.438E+01, 4.732E+00, 0.000E+00/
      DATA (PH1(I,25,18, 3),I=1,6) /8.191E+02, 1.748E+02,
     1 1.921E+02, 5.809E+01, 3.973E+00, 1.717E-02/
      DATA (PH1(I,25,18, 4),I=1,6) /2.282E+02, 5.334E+01,
     1 7.236E+00, 4.454E+02, 3.624E+00, 0.000E+00/
      DATA (PH1(I,25,18, 5),I=1,6) /1.945E+02, 7.818E+01,
     1 4.992E+01, 4.763E+01, 5.832E+00, 3.093E-01/
      DATA (PH1(I,25,19, 1),I=1,6) /6.635E+03, 1.021E+03,
     1 1.154E+01, 3.219E+01, 1.919E+00, 0.000E+00/
      DATA (PH1(I,25,19, 2),I=1,6) /8.900E+02, 1.668E+02,
     1 8.196E+00, 2.394E+01, 4.767E+00, 0.000E+00/
      DATA (PH1(I,25,19, 3),I=1,6) /7.860E+02, 1.739E+02,
     1 1.939E+02, 5.345E+01, 4.041E+00, 1.659E-02/
      DATA (PH1(I,25,19, 4),I=1,6) /2.016E+02, 5.184E+01,
     1 6.791E+00, 3.811E+02, 3.719E+00, 0.000E+00/
      DATA (PH1(I,25,19, 5),I=1,6) /1.691E+02, 7.618E+01,
     1 4.935E+01, 3.879E+01, 6.116E+00, 3.087E-01/
      DATA (PH1(I,25,19, 6),I=1,6) /1.193E+02, 3.865E+01,
     1 9.713E+01, 1.133E+02, 5.798E+00, 1.437E-01/
      DATA (PH1(I,25,20, 1),I=1,6) /6.617E+03, 1.171E+03,
     1 8.591E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25,20, 2),I=1,6) /8.614E+02, 1.867E+02,
     1 7.373E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25,20, 3),I=1,6) /7.560E+02, 2.007E+02,
     1 1.450E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,25,20, 4),I=1,6) /1.768E+02, 4.271E+01,
     1 7.047E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,25,20, 5),I=1,6) /1.444E+02, 7.861E+01,
     1 4.367E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,25,20, 6),I=1,6) /9.575E+01, 3.452E+01,
     1 2.355E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,25,21, 1),I=1,6) /6.602E+03, 1.045E+03,
     1 1.101E+01, 3.201E+01, 1.902E+00, 0.000E+00/
      DATA (PH1(I,25,21, 2),I=1,6) /8.319E+02, 1.782E+02,
     1 7.660E+00, 2.515E+01, 4.612E+00, 0.000E+00/
      DATA (PH1(I,25,21, 3),I=1,6) /7.280E+02, 1.909E+02,
     1 1.610E+02, 6.064E+01, 3.847E+00, 2.553E-02/
      DATA (PH1(I,25,21, 4),I=1,6) /1.530E+02, 5.019E+01,
     1 5.951E+00, 3.117E+02, 3.866E+00, 0.000E+00/
      DATA (PH1(I,25,21, 5),I=1,6) /1.210E+02, 6.892E+01,
     1 5.293E+01, 2.594E+01, 6.855E+00, 1.836E-01/
      DATA (PH1(I,25,21, 6),I=1,6) /7.240E+01, 3.798E+01,
     1 2.289E+02, 7.678E+01, 6.206E+00, 6.108E-02/
      DATA (PH1(I,25,22, 1),I=1,6) /6.589E+03, 1.174E+03,
     1 8.541E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,25,22, 2),I=1,6) /8.118E+02, 1.872E+02,
     1 7.318E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,25,22, 3),I=1,6) /7.063E+02, 2.015E+02,
     1 1.433E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,25,22, 4),I=1,6) /1.300E+02, 4.531E+01,
     1 5.863E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,25,22, 5),I=1,6) /9.879E+01, 8.115E+01,
     1 3.665E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,25,22, 6),I=1,6) /5.120E+01, 3.615E+01,
     1 3.188E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,25,23, 1),I=1,6) /6.576E+03, 1.154E+03,
     1 8.979E+00, 2.899E+01, 1.859E+00, 0.000E+00/
      DATA (PH1(I,25,23, 2),I=1,6) /7.947E+02, 2.016E+02,
     1 6.777E+00, 2.756E+01, 4.332E+00, 0.000E+00/
      DATA (PH1(I,25,23, 3),I=1,6) /6.822E+02, 1.908E+02,
     1 1.609E+02, 6.198E+01, 3.835E+00, 2.549E-02/
      DATA (PH1(I,25,23, 4),I=1,6) /1.120E+02, 4.956E+01,
     1 5.285E+00, 2.544E+02, 4.005E+00, 0.000E+00/
      DATA (PH1(I,25,23, 5),I=1,6) /8.062E+01, 7.190E+01,
     1 4.707E+01, 2.356E+01, 6.977E+00, 2.565E-01/
      DATA (PH1(I,25,23, 6),I=1,6) /3.367E+01, 1.604E+01,
     1 1.991E+05, 2.538E+00, 1.808E+01, 6.536E-01/
      DATA (PH1(I,25,24, 1),I=1,6) /6.564E+03, 8.684E+02,
     1 1.627E+01, 2.405E+01, 2.207E+00, 0.000E+00/
      DATA (PH1(I,25,24, 2),I=1,6) /7.849E+02, 1.346E+02,
     1 9.823E+00, 1.970E+01, 5.386E+00, 0.000E+00/
      DATA (PH1(I,25,24, 3),I=1,6) /6.714E+02, 2.207E+02,
     1 1.185E+02, 9.212E+01, 3.460E+00, 1.974E-01/
      DATA (PH1(I,25,24, 4),I=1,6) /1.015E+02, 4.359E+01,
     1 5.735E+00, 1.016E+02, 4.563E+00, 0.000E+00/
      DATA (PH1(I,25,24, 5),I=1,6) /7.009E+01, 7.141E+01,
     1 4.984E+01, 2.130E+01, 7.133E+00, 2.616E-01/
      DATA (PH1(I,25,24, 6),I=1,6) /2.058E+01, 1.624E+01,
     1 2.592E+05, 2.411E+00, 1.797E+01, 4.774E-01/
      DATA (PH1(I,25,24, 7),I=1,6) /1.564E+01, 9.306E+00,
     1 1.189E+00, 1.439E+09, 4.095E+00, 8.740E-04/
      DATA (PH1(I,25,25, 1),I=1,6) /6.550E+03, 8.758E+02,
     1 1.592E+01, 3.965E+01, 1.947E+00, 0.000E+00/
      DATA (PH1(I,25,25, 2),I=1,6) /7.816E+02, 8.316E+01,
     1 1.156E+01, 2.187E+01, 6.149E+00, 0.000E+00/
      DATA (PH1(I,25,25, 3),I=1,6) /6.554E+02, 2.737E+02,
     1 7.498E+01, 3.952E+02, 2.793E+00, 4.661E-02/
      DATA (PH1(I,25,25, 4),I=1,6) /9.460E+01, 4.114E+01,
     1 6.172E+00, 5.928E+01, 4.989E+00, 0.000E+00/
      DATA (PH1(I,25,25, 5),I=1,6) /5.940E+01, 7.040E+01,
     1 6.697E+01, 1.485E+01, 7.643E+00, 2.659E-01/
      DATA (PH1(I,25,25, 6),I=1,6) /1.430E+01, 1.311E+01,
     1 1.668E+04, 4.497E+00, 1.646E+01, 3.881E-01/
      DATA (PH1(I,25,25, 7),I=1,6) /7.434E+00, 1.183E+01,
     1 1.549E+00, 2.920E+06, 4.113E+00, 3.256E-02/
      DATA (PH1(I,26, 1, 1),I=1,6) /9.278E+03, 2.932E+02,
     1 8.099E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,26, 2, 1),I=1,6) /8.829E+03, 1.057E+03,
     1 1.195E+01, 5.769E+01, 1.718E+00, 0.000E+00/
      DATA (PH1(I,26, 3, 1),I=1,6) /8.638E+03, 1.087E+03,
     1 1.157E+01, 5.086E+01, 1.722E+00, 0.000E+00/
      DATA (PH1(I,26, 3, 2),I=1,6) /2.046E+03, 1.873E+02,
     1 5.833E+00, 3.849E+01, 3.998E+00, 0.000E+00/
      DATA (PH1(I,26, 4, 1),I=1,6) /8.484E+03, 1.215E+03,
     1 9.098E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,26, 4, 2),I=1,6) /1.950E+03, 1.799E+02,
     1 1.138E+01, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26, 5, 1),I=1,6) /8.350E+03, 1.066E+03,
     1 1.181E+01, 4.116E+01, 1.827E+00, 0.000E+00/
      DATA (PH1(I,26, 5, 2),I=1,6) /1.847E+03, 1.708E+02,
     1 1.125E+01, 2.839E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26, 5, 3),I=1,6) /1.799E+03, 8.070E+01,
     1 2.675E+02, 1.191E+02, 4.194E+00, 1.901E-05/
      DATA (PH1(I,26, 6, 1),I=1,6) /8.184E+03, 1.228E+03,
     1 8.773E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,26, 6, 2),I=1,6) /1.745E+03, 1.860E+02,
     1 9.942E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26, 6, 3),I=1,6) /1.689E+03, 1.970E+02,
     1 7.720E+01, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,26, 7, 1),I=1,6) /8.041E+03, 1.235E+03,
     1 8.619E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,26, 7, 2),I=1,6) /1.648E+03, 1.888E+02,
     1 9.298E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26, 7, 3),I=1,6) /1.582E+03, 1.975E+02,
     1 1.093E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,26, 8, 1),I=1,6) /7.918E+03, 1.068E+03,
     1 1.155E+01, 3.578E+01, 1.895E+00, 0.000E+00/
      DATA (PH1(I,26, 8, 2),I=1,6) /1.559E+03, 1.623E+02,
     1 1.000E+01, 2.675E+01, 4.711E+00, 0.000E+00/
      DATA (PH1(I,26, 8, 3),I=1,6) /1.456E+03, 1.338E+02,
     1 3.000E+02, 7.010E+01, 4.143E+00, 1.631E-05/
      DATA (PH1(I,26, 9, 1),I=1,6) /7.769E+03, 1.249E+03,
     1 8.327E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,26, 9, 2),I=1,6) /1.460E+03, 1.940E+02,
     1 8.150E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26, 9, 3),I=1,6) /1.358E+03, 1.999E+02,
     1 1.561E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,26,10, 1),I=1,6) /7.651E+03, 1.067E+03,
     1 1.150E+01, 3.412E+01, 1.922E+00, 0.000E+00/
      DATA (PH1(I,26,10, 2),I=1,6) /1.375E+03, 1.534E+02,
     1 9.191E+00, 2.691E+01, 4.849E+00, 0.000E+00/
      DATA (PH1(I,26,10, 3),I=1,6) /1.262E+03, 1.412E+02,
     1 3.368E+02, 5.569E+01, 4.328E+00, 1.177E-05/
      DATA (PH1(I,26,11, 1),I=1,6) /7.599E+03, 1.061E+03,
     1 1.161E+01, 3.713E+01, 1.889E+00, 0.000E+00/
      DATA (PH1(I,26,11, 2),I=1,6) /1.329E+03, 1.363E+02,
     1 9.940E+00, 3.026E+01, 4.885E+00, 0.000E+00/
      DATA (PH1(I,26,11, 3),I=1,6) /1.216E+03, 1.396E+02,
     1 3.383E+02, 5.459E+01, 4.366E+00, 1.179E-05/
      DATA (PH1(I,26,11, 4),I=1,6) /4.893E+02, 2.873E+01,
     1 1.207E+01, 5.150E+02, 3.846E+00, 0.000E+00/
      DATA (PH1(I,26,12, 1),I=1,6) /7.553E+03, 1.258E+03,
     1 8.097E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,26,12, 2),I=1,6) /1.287E+03, 1.967E+02,
     1 7.555E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,26,12, 3),I=1,6) /1.181E+03, 2.136E+02,
     1 1.496E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,26,12, 4),I=1,6) /4.570E+02, 4.064E+01,
     1 1.228E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,26,13, 1),I=1,6) /7.499E+03, 1.112E+03,
     1 1.052E+01, 3.562E+01, 1.871E+00, 0.000E+00/
      DATA (PH1(I,26,13, 2),I=1,6) /1.240E+03, 1.724E+02,
     1 8.240E+00, 2.693E+01, 4.675E+00, 0.000E+00/
      DATA (PH1(I,26,13, 3),I=1,6) /1.125E+03, 1.387E+02,
     1 3.367E+02, 5.279E+01, 4.407E+00, 1.111E-05/
      DATA (PH1(I,26,13, 4),I=1,6) /4.238E+02, 7.077E+01,
     1 7.969E+00, 5.650E+02, 3.329E+00, 0.000E+00/
      DATA (PH1(I,26,13, 5),I=1,6) /3.922E+02, 7.993E+01,
     1 1.135E+01, 5.807E+01, 5.676E+00, 6.512E-02/
      DATA (PH1(I,26,14, 1),I=1,6) /7.450E+03, 1.103E+03,
     1 1.069E+01, 3.585E+01, 1.875E+00, 0.000E+00/
      DATA (PH1(I,26,14, 2),I=1,6) /1.197E+03, 1.722E+02,
     1 8.199E+00, 2.626E+01, 4.710E+00, 0.000E+00/
      DATA (PH1(I,26,14, 3),I=1,6) /1.081E+03, 1.392E+02,
     1 3.316E+02, 5.237E+01, 4.410E+00, 1.107E-05/
      DATA (PH1(I,26,14, 4),I=1,6) /3.916E+02, 6.712E+01,
     1 8.093E+00, 5.886E+02, 3.360E+00, 0.000E+00/
      DATA (PH1(I,26,14, 5),I=1,6) /3.610E+02, 8.203E+01,
     1 2.092E+01, 6.128E+01, 5.596E+00, 3.441E-02/
      DATA (PH1(I,26,15, 1),I=1,6) /7.403E+03, 1.096E+03,
     1 1.083E+01, 3.612E+01, 1.878E+00, 0.000E+00/
      DATA (PH1(I,26,15, 2),I=1,6) /1.155E+03, 1.720E+02,
     1 8.161E+00, 2.570E+01, 4.741E+00, 0.000E+00/
      DATA (PH1(I,26,15, 3),I=1,6) /1.039E+03, 1.383E+02,
     1 3.325E+02, 5.090E+01, 4.446E+00, 1.114E-05/
      DATA (PH1(I,26,15, 4),I=1,6) /3.600E+02, 6.329E+01,
     1 8.143E+00, 5.637E+02, 3.421E+00, 0.000E+00/
      DATA (PH1(I,26,15, 5),I=1,6) /3.308E+02, 8.150E+01,
     1 3.022E+01, 6.039E+01, 5.612E+00, 2.772E-02/
      DATA (PH1(I,26,16, 1),I=1,6) /7.359E+03, 1.095E+03,
     1 1.084E+01, 3.598E+01, 1.880E+00, 0.000E+00/
      DATA (PH1(I,26,16, 2),I=1,6) /1.115E+03, 1.714E+02,
     1 8.138E+00, 2.516E+01, 4.774E+00, 0.000E+00/
      DATA (PH1(I,26,16, 3),I=1,6) /9.983E+02, 1.632E+02,
     1 2.449E+02, 5.452E+01, 4.187E+00, 1.594E-05/
      DATA (PH1(I,26,16, 4),I=1,6) /3.292E+02, 6.060E+01,
     1 7.963E+00, 5.214E+02, 3.487E+00, 0.000E+00/
      DATA (PH1(I,26,16, 5),I=1,6) /2.902E+02, 7.972E+01,
     1 3.922E+01, 5.574E+01, 5.729E+00, 2.828E-02/
      DATA (PH1(I,26,17, 1),I=1,6) /7.316E+03, 1.098E+03,
     1 1.078E+01, 3.596E+01, 1.878E+00, 0.000E+00/
      DATA (PH1(I,26,17, 2),I=1,6) /1.076E+03, 1.710E+02,
     1 8.135E+00, 2.440E+01, 4.817E+00, 0.000E+00/
      DATA (PH1(I,26,17, 3),I=1,6) /9.590E+02, 1.647E+02,
     1 2.392E+02, 5.276E+01, 4.204E+00, 1.574E-05/
      DATA (PH1(I,26,17, 4),I=1,6) /2.990E+02, 5.854E+01,
     1 7.650E+00, 4.679E+02, 3.561E+00, 0.000E+00/
      DATA (PH1(I,26,17, 5),I=1,6) /2.621E+02, 7.914E+01,
     1 4.686E+01, 5.125E+01, 5.835E+00, 2.574E-02/
      DATA (PH1(I,26,18, 1),I=1,6) /7.275E+03, 1.102E+03,
     1 1.068E+01, 3.576E+01, 1.878E+00, 0.000E+00/
      DATA (PH1(I,26,18, 2),I=1,6) /1.039E+03, 1.706E+02,
     1 8.130E+00, 2.387E+01, 4.850E+00, 0.000E+00/
      DATA (PH1(I,26,18, 3),I=1,6) /9.211E+02, 1.715E+02,
     1 2.205E+02, 5.298E+01, 4.154E+00, 1.508E-05/
      DATA (PH1(I,26,18, 4),I=1,6) /2.696E+02, 5.701E+01,
     1 7.237E+00, 4.092E+02, 3.643E+00, 0.000E+00/
      DATA (PH1(I,26,18, 5),I=1,6) /2.336E+02, 7.704E+01,
     1 5.530E+01, 4.380E+01, 6.059E+00, 2.673E-02/
      DATA (PH1(I,26,19, 1),I=1,6) /7.237E+03, 1.103E+03,
     1 1.067E+01, 3.584E+01, 1.876E+00, 0.000E+00/
      DATA (PH1(I,26,19, 2),I=1,6) /1.003E+03, 1.694E+02,
     1 8.158E+00, 2.345E+01, 4.887E+00, 0.000E+00/
      DATA (PH1(I,26,19, 3),I=1,6) /8.849E+02, 1.788E+02,
     1 2.032E+02, 5.293E+01, 4.107E+00, 1.398E-05/
      DATA (PH1(I,26,19, 4),I=1,6) /2.409E+02, 5.450E+01,
     1 6.938E+00, 3.605E+02, 3.741E+00, 0.000E+00/
      DATA (PH1(I,26,19, 5),I=1,6) /2.055E+02, 7.556E+01,
     1 5.336E+01, 3.855E+01, 6.265E+00, 2.667E-02/
      DATA (PH1(I,26,19, 6),I=1,6) /1.511E+02, 3.652E+01,
     1 1.323E+02, 9.277E+01, 6.150E+00, 2.051E-02/
      DATA (PH1(I,26,20, 1),I=1,6) /7.217E+03, 1.110E+03,
     1 1.055E+01, 3.563E+01, 1.873E+00, 0.000E+00/
      DATA (PH1(I,26,20, 2),I=1,6) /9.693E+02, 1.694E+02,
     1 8.138E+00, 2.316E+01, 4.905E+00, 0.000E+00/
      DATA (PH1(I,26,20, 3),I=1,6) /8.512E+02, 1.840E+02,
     1 1.920E+02, 5.238E+01, 4.082E+00, 1.377E-05/
      DATA (PH1(I,26,20, 4),I=1,6) /2.135E+02, 5.365E+01,
     1 6.463E+00, 3.169E+02, 3.824E+00, 0.000E+00/
      DATA (PH1(I,26,20, 5),I=1,6) /1.783E+02, 7.433E+01,
     1 5.210E+01, 3.306E+01, 6.509E+00, 2.404E-02/
      DATA (PH1(I,26,20, 6),I=1,6) /1.250E+02, 4.332E+01,
     1 1.701E+02, 9.791E+01, 5.905E+00, 2.604E-02/
      DATA (PH1(I,26,21, 1),I=1,6) /7.199E+03, 1.113E+03,
     1 1.049E+01, 3.537E+01, 1.874E+00, 0.000E+00/
      DATA (PH1(I,26,21, 2),I=1,6) /9.383E+02, 1.803E+02,
     1 7.682E+00, 2.415E+01, 4.757E+00, 0.000E+00/
      DATA (PH1(I,26,21, 3),I=1,6) /8.202E+02, 1.952E+02,
     1 1.709E+02, 5.690E+01, 3.953E+00, 1.330E-05/
      DATA (PH1(I,26,21, 4),I=1,6) /1.876E+02, 5.277E+01,
     1 6.047E+00, 2.797E+02, 3.910E+00, 0.000E+00/
      DATA (PH1(I,26,21, 5),I=1,6) /1.527E+02, 7.260E+01,
     1 5.246E+01, 2.751E+01, 6.823E+00, 2.105E-02/
      DATA (PH1(I,26,21, 6),I=1,6) /9.906E+01, 4.487E+01,
     1 2.082E+02, 9.579E+01, 5.914E+00, 2.368E-02/
      DATA (PH1(I,26,22, 1),I=1,6) /7.184E+03, 1.124E+03,
     1 1.027E+01, 3.534E+01, 1.866E+00, 0.000E+00/
      DATA (PH1(I,26,22, 2),I=1,6) /9.101E+02, 1.818E+02,
     1 7.623E+00, 2.401E+01, 4.752E+00, 0.000E+00/
      DATA (PH1(I,26,22, 3),I=1,6) /7.920E+02, 1.996E+02,
     1 1.636E+02, 5.719E+01, 3.924E+00, 1.320E-05/
      DATA (PH1(I,26,22, 4),I=1,6) /1.633E+02, 5.218E+01,
     1 5.672E+00, 2.482E+02, 3.991E+00, 0.000E+00/
      DATA (PH1(I,26,22, 5),I=1,6) /1.288E+02, 7.273E+01,
     1 5.143E+01, 2.428E+01, 7.028E+00, 1.365E-01/
      DATA (PH1(I,26,22, 6),I=1,6) /7.501E+01, 4.369E+01,
     1 2.661E+02, 6.273E+01, 6.344E+00, 2.402E-01/
      DATA (PH1(I,26,23, 1),I=1,6) /7.169E+03, 1.155E+03,
     1 9.711E+00, 3.329E+01, 1.869E+00, 0.000E+00/
      DATA (PH1(I,26,23, 2),I=1,6) /8.871E+02, 2.059E+02,
     1 6.695E+00, 2.772E+01, 4.410E+00, 0.000E+00/
      DATA (PH1(I,26,23, 3),I=1,6) /7.669E+02, 2.032E+02,
     1 1.578E+02, 5.933E+01, 3.879E+00, 1.308E-05/
      DATA (PH1(I,26,23, 4),I=1,6) /1.411E+02, 5.256E+01,
     1 5.295E+00, 2.332E+02, 4.036E+00, 0.000E+00/
      DATA (PH1(I,26,23, 5),I=1,6) /1.067E+02, 7.625E+01,
     1 4.810E+01, 2.286E+01, 7.043E+00, 2.449E-01/
      DATA (PH1(I,26,23, 6),I=1,6) /5.480E+01, 4.154E+01,
     1 3.683E+02, 3.529E+01, 7.056E+00, 3.316E-01/
      DATA (PH1(I,26,24, 1),I=1,6) /7.155E+03, 1.245E+03,
     1 8.309E+00, 3.170E+01, 1.826E+00, 0.000E+00/
      DATA (PH1(I,26,24, 2),I=1,6) /8.710E+02, 2.158E+02,
     1 6.454E+00, 2.722E+01, 4.355E+00, 0.000E+00/
      DATA (PH1(I,26,24, 3),I=1,6) /7.451E+02, 2.023E+02,
     1 1.591E+02, 5.990E+01, 3.878E+00, 1.309E-05/
      DATA (PH1(I,26,24, 4),I=1,6) /1.211E+02, 5.277E+01,
     1 5.015E+00, 2.131E+02, 4.093E+00, 0.000E+00/
      DATA (PH1(I,26,24, 5),I=1,6) /8.705E+01, 7.777E+01,
     1 4.422E+01, 2.336E+01, 7.017E+00, 2.557E-01/
      DATA (PH1(I,26,24, 6),I=1,6) /3.065E+01, 2.670E+01,
     1 6.301E+03, 5.385E+00, 1.232E+01, 4.407E-01/
      DATA (PH1(I,26,25, 1),I=1,6) /7.140E+03, 8.931E+02,
     1 1.666E+01, 2.381E+01, 2.262E+00, 0.000E+00/
      DATA (PH1(I,26,25, 2),I=1,6) /8.608E+02, 1.431E+02,
     1 9.289E+00, 2.026E+01, 5.376E+00, 0.000E+00/
      DATA (PH1(I,26,25, 3),I=1,6) /7.341E+02, 2.281E+02,
     1 1.241E+02, 8.058E+01, 3.572E+00, 1.223E-01/
      DATA (PH1(I,26,25, 4),I=1,6) /1.102E+02, 4.663E+01,
     1 5.434E+00, 9.271E+01, 4.640E+00, 0.000E+00/
      DATA (PH1(I,26,25, 5),I=1,6) /7.617E+01, 7.750E+01,
     1 4.624E+01, 2.155E+01, 7.138E+00, 2.599E-01/
      DATA (PH1(I,26,25, 6),I=1,6) /2.193E+01, 1.933E+01,
     1 3.679E+04, 3.564E+00, 1.566E+01, 4.144E-01/
      DATA (PH1(I,26,25, 7),I=1,6) /1.619E+01, 1.014E+01,
     1 1.084E+00, 2.562E+04, 4.167E+00, 1.598E-02/
      DATA (PH1(I,26,26, 1),I=1,6) /7.124E+03, 8.044E+02,
     1 2.055E+01, 3.633E+01, 2.118E+00, 0.000E+00/
      DATA (PH1(I,26,26, 2),I=1,6) /8.570E+02, 5.727E+01,
     1 1.076E+01, 2.785E+01, 6.635E+00, 0.000E+00/
      DATA (PH1(I,26,26, 3),I=1,6) /7.240E+02, 2.948E+02,
     1 7.191E+01, 3.219E+02, 2.837E+00, 6.314E-02/
      DATA (PH1(I,26,26, 4),I=1,6) /1.040E+02, 4.334E+01,
     1 5.921E+00, 5.293E+01, 5.129E+00, 0.000E+00/
      DATA (PH1(I,26,26, 5),I=1,6) /6.600E+01, 7.630E+01,
     1 6.298E+01, 1.479E+01, 7.672E+00, 2.646E-01/
      DATA (PH1(I,26,26, 6),I=1,6) /1.470E+01, 1.407E+01,
     1 1.850E+04, 4.458E+00, 1.691E+01, 4.039E-01/
      DATA (PH1(I,26,26, 7),I=1,6) /7.902E+00, 1.277E+01,
     1 1.468E+00, 1.116E+05, 4.112E+00, 3.238E-02/
      DATA (PH1(I,27, 1, 1),I=1,6) /1.001E+04, 3.163E+02,
     1 7.510E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,27, 2, 1),I=1,6) /9.545E+03, 1.251E+03,
     1 9.046E+00, 6.786E+01, 1.613E+00, 0.000E+00/
      DATA (PH1(I,27, 3, 1),I=1,6) /9.347E+03, 1.249E+03,
     1 9.384E+00, 6.378E+01, 1.609E+00, 0.000E+00/
      DATA (PH1(I,27, 3, 2),I=1,6) /2.219E+03, 2.352E+02,
     1 4.599E+00, 3.873E+01, 3.808E+00, 0.000E+00/
      DATA (PH1(I,27, 4, 1),I=1,6) /9.205E+03, 1.267E+03,
     1 9.026E+00, 6.059E+01, 1.619E+00, 0.000E+00/
      DATA (PH1(I,27, 4, 2),I=1,6) /2.119E+03, 1.593E+02,
     1 1.253E+01, 3.366E+01, 4.493E+00, 0.000E+00/
      DATA (PH1(I,27, 5, 1),I=1,6) /9.046E+03, 1.275E+03,
     1 8.806E+00, 5.513E+01, 1.651E+00, 0.000E+00/
      DATA (PH1(I,27, 5, 2),I=1,6) /2.012E+03, 1.799E+02,
     1 1.074E+01, 2.878E+01, 4.518E+00, 0.000E+00/
      DATA (PH1(I,27, 5, 3),I=1,6) /1.961E+03, 1.923E+02,
     1 4.987E+01, 9.392E+01, 3.631E+00, 2.239E-05/
      DATA (PH1(I,27, 6, 1),I=1,6) /8.890E+03, 1.277E+03,
     1 8.703E+00, 5.267E+01, 1.668E+00, 0.000E+00/
      DATA (PH1(I,27, 6, 2),I=1,6) /1.910E+03, 1.760E+02,
     1 1.043E+01, 2.775E+01, 4.612E+00, 0.000E+00/
      DATA (PH1(I,27, 6, 3),I=1,6) /1.846E+03, 2.159E+02,
     1 7.348E+01, 9.457E+01, 3.535E+00, 1.182E-05/
      DATA (PH1(I,27, 7, 1),I=1,6) /8.718E+03, 1.330E+03,
     1 8.042E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27, 7, 2),I=1,6) /1.803E+03, 2.030E+02,
     1 8.771E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27, 7, 3),I=1,6) /1.735E+03, 2.134E+02,
     1 1.027E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,27, 8, 1),I=1,6) /8.595E+03, 1.260E+03,
     1 8.868E+00, 4.790E+01, 1.717E+00, 0.000E+00/
      DATA (PH1(I,27, 8, 2),I=1,6) /1.711E+03, 1.706E+02,
     1 9.597E+00, 2.731E+01, 4.722E+00, 0.000E+00/
      DATA (PH1(I,27, 8, 3),I=1,6) /1.603E+03, 2.448E+02,
     1 9.893E+01, 9.263E+01, 3.432E+00, 9.870E-06/
      DATA (PH1(I,27, 9, 1),I=1,6) /8.433E+03, 1.344E+03,
     1 7.779E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27, 9, 2),I=1,6) /1.606E+03, 2.084E+02,
     1 7.728E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27, 9, 3),I=1,6) /1.505E+03, 2.157E+02,
     1 1.481E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,27,10, 1),I=1,6) /8.315E+03, 1.252E+03,
     1 8.903E+00, 4.635E+01, 1.738E+00, 0.000E+00/
      DATA (PH1(I,27,10, 2),I=1,6) /1.519E+03, 1.643E+02,
     1 8.754E+00, 2.695E+01, 4.851E+00, 0.000E+00/
      DATA (PH1(I,27,10, 3),I=1,6) /1.397E+03, 2.062E+02,
     1 1.816E+02, 6.860E+01, 3.824E+00, 9.070E-06/
      DATA (PH1(I,27,11, 1),I=1,6) /8.260E+03, 1.249E+03,
     1 8.967E+00, 4.909E+01, 1.719E+00, 0.000E+00/
      DATA (PH1(I,27,11, 2),I=1,6) /1.470E+03, 1.433E+02,
     1 9.597E+00, 2.996E+01, 4.928E+00, 0.000E+00/
      DATA (PH1(I,27,11, 3),I=1,6) /1.362E+03, 2.054E+02,
     1 1.804E+02, 6.660E+01, 3.855E+00, 8.999E-06/
      DATA (PH1(I,27,11, 4),I=1,6) /5.466E+02, 1.262E+01,
     1 2.734E+01, 6.463E+02, 4.280E+00, 0.000E+00/
      DATA (PH1(I,27,12, 1),I=1,6) /8.207E+03, 1.277E+03,
     1 8.536E+00, 4.760E+01, 1.714E+00, 0.000E+00/
      DATA (PH1(I,27,12, 2),I=1,6) /1.423E+03, 1.673E+02,
     1 8.543E+00, 2.705E+01, 4.820E+00, 0.000E+00/
      DATA (PH1(I,27,12, 3),I=1,6) /1.314E+03, 2.058E+02,
     1 1.783E+02, 6.530E+01, 3.868E+00, 8.972E-06/
      DATA (PH1(I,27,12, 4),I=1,6) /5.120E+02, 2.690E+01,
     1 2.122E+01, 3.166E+02, 4.210E+00, 0.000E+00/
      DATA (PH1(I,27,13, 1),I=1,6) /8.154E+03, 1.277E+03,
     1 8.527E+00, 4.599E+01, 1.727E+00, 0.000E+00/
      DATA (PH1(I,27,13, 2),I=1,6) /1.376E+03, 1.726E+02,
     1 8.292E+00, 2.627E+01, 4.812E+00, 0.000E+00/
      DATA (PH1(I,27,13, 3),I=1,6) /1.266E+03, 1.933E+02,
     1 2.002E+02, 6.453E+01, 3.945E+00, 9.099E-06/
      DATA (PH1(I,27,13, 4),I=1,6) /4.777E+02, 2.652E+01,
     1 1.977E+01, 2.830E+02, 4.297E+00, 0.000E+00/
      DATA (PH1(I,27,13, 5),I=1,6) /4.440E+02, 1.358E+02,
     1 6.171E+00, 4.715E+01, 5.247E+00, 8.240E-01/
      DATA (PH1(I,27,14, 1),I=1,6) /8.102E+03, 1.272E+03,
     1 8.593E+00, 4.604E+01, 1.730E+00, 0.000E+00/
      DATA (PH1(I,27,14, 2),I=1,6) /1.331E+03, 1.718E+02,
     1 8.247E+00, 2.592E+01, 4.840E+00, 0.000E+00/
      DATA (PH1(I,27,14, 3),I=1,6) /1.219E+03, 1.931E+02,
     1 1.993E+02, 6.407E+01, 3.953E+00, 9.111E-06/
      DATA (PH1(I,27,14, 4),I=1,6) /4.436E+02, 2.529E+01,
     1 1.979E+01, 2.943E+02, 4.315E+00, 0.000E+00/
      DATA (PH1(I,27,14, 5),I=1,6) /4.110E+02, 9.918E+01,
     1 1.786E+01, 5.216E+01, 5.592E+00, 5.707E-01/
      DATA (PH1(I,27,15, 1),I=1,6) /8.045E+03, 1.361E+03,
     1 7.456E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27,15, 2),I=1,6) /1.281E+03, 2.135E+02,
     1 6.900E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27,15, 3),I=1,6) /1.167E+03, 2.301E+02,
     1 1.403E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,27,15, 4),I=1,6) /4.088E+02, 4.412E+01,
     1 1.052E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,27,15, 5),I=1,6) /3.790E+02, 8.885E+01,
     1 2.993E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,27,16, 1),I=1,6) /8.005E+03, 1.264E+03,
     1 8.703E+00, 4.638E+01, 1.732E+00, 0.000E+00/
      DATA (PH1(I,27,16, 2),I=1,6) /1.244E+03, 1.755E+02,
     1 8.017E+00, 2.508E+01, 4.854E+00, 0.000E+00/
      DATA (PH1(I,27,16, 3),I=1,6) /1.131E+03, 1.925E+02,
     1 1.972E+02, 6.228E+01, 3.980E+00, 9.096E-06/
      DATA (PH1(I,27,16, 4),I=1,6) /3.775E+02, 2.185E+01,
     1 1.954E+01, 3.125E+02, 4.423E+00, 0.000E+00/
      DATA (PH1(I,27,16, 5),I=1,6) /3.360E+02, 9.235E+01,
     1 3.595E+01, 4.513E+01, 5.839E+00, 5.704E-01/
      DATA (PH1(I,27,17, 1),I=1,6) /7.951E+03, 1.362E+03,
     1 7.447E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27,17, 2),I=1,6) /1.196E+03, 2.146E+02,
     1 6.772E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27,17, 3),I=1,6) /1.080E+03, 2.305E+02,
     1 1.383E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,27,17, 4),I=1,6) /3.439E+02, 4.522E+01,
     1 9.201E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,27,17, 5),I=1,6) /3.053E+02, 8.909E+01,
     1 4.487E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,27,18, 1),I=1,6) /7.915E+03, 1.266E+03,
     1 8.671E+00, 4.646E+01, 1.731E+00, 0.000E+00/
      DATA (PH1(I,27,18, 2),I=1,6) /1.163E+03, 1.808E+02,
     1 7.771E+00, 2.435E+01, 4.849E+00, 0.000E+00/
      DATA (PH1(I,27,18, 3),I=1,6) /1.048E+03, 1.922E+02,
     1 1.949E+02, 5.995E+01, 4.013E+00, 9.089E-06/
      DATA (PH1(I,27,18, 4),I=1,6) /3.142E+02, 3.099E+01,
     1 1.193E+01, 2.258E+02, 4.398E+00, 0.000E+00/
      DATA (PH1(I,27,18, 5),I=1,6) /2.754E+02, 9.052E+01,
     1 5.045E+01, 4.052E+01, 6.001E+00, 4.620E-01/
      DATA (PH1(I,27,19, 1),I=1,6) /7.877E+03, 1.287E+03,
     1 8.375E+00, 4.957E+01, 1.697E+00, 0.000E+00/
      DATA (PH1(I,27,19, 2),I=1,6) /1.123E+03, 1.802E+02,
     1 7.760E+00, 2.395E+01, 4.878E+00, 0.000E+00/
      DATA (PH1(I,27,19, 3),I=1,6) /1.009E+03, 1.943E+02,
     1 1.901E+02, 5.795E+01, 4.029E+00, 8.959E-06/
      DATA (PH1(I,27,19, 4),I=1,6) /2.833E+02, 3.765E+01,
     1 9.284E+00, 1.936E+02, 4.355E+00, 0.000E+00/
      DATA (PH1(I,27,19, 5),I=1,6) /2.450E+02, 9.011E+01,
     1 4.789E+01, 3.807E+01, 6.102E+00, 4.132E-01/
      DATA (PH1(I,27,19, 6),I=1,6) /1.861E+02, 3.461E+01,
     1 1.735E+02, 7.681E+01, 6.518E+00, 1.706E-04/
      DATA (PH1(I,27,20, 1),I=1,6) /7.840E+03, 1.363E+03,
     1 7.433E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27,20, 2),I=1,6) /1.086E+03, 2.158E+02,
     1 6.639E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27,20, 3),I=1,6) /9.692E+02, 2.313E+02,
     1 1.357E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,27,20, 4),I=1,6) /2.544E+02, 4.790E+01,
     1 7.244E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,27,20, 5),I=1,6) /2.158E+02, 9.095E+01,
     1 4.368E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,27,20, 6),I=1,6) /1.578E+02, 4.257E+01,
     1 2.299E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,27,21, 1),I=1,6) /7.820E+03, 1.295E+03,
     1 8.260E+00, 5.104E+01, 1.682E+00, 0.000E+00/
      DATA (PH1(I,27,21, 2),I=1,6) /1.052E+03, 1.797E+02,
     1 7.740E+00, 2.344E+01, 4.914E+00, 0.000E+00/
      DATA (PH1(I,27,21, 3),I=1,6) /9.379E+02, 2.015E+02,
     1 1.767E+02, 5.582E+01, 4.018E+00, 8.773E-06/
      DATA (PH1(I,27,21, 4),I=1,6) /2.255E+02, 4.430E+01,
     1 6.990E+00, 1.659E+02, 4.350E+00, 0.000E+00/
      DATA (PH1(I,27,21, 5),I=1,6) /1.876E+02, 8.695E+01,
     1 4.460E+01, 3.324E+01, 6.378E+00, 3.125E-01/
      DATA (PH1(I,27,21, 6),I=1,6) /1.290E+02, 4.433E+01,
     1 2.748E+02, 5.829E+01, 6.482E+00, 1.086E-04/
      DATA (PH1(I,27,22, 1),I=1,6) /7.803E+03, 1.364E+03,
     1 7.423E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27,22, 2),I=1,6) /1.025E+03, 2.163E+02,
     1 6.587E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27,22, 3),I=1,6) /9.074E+02, 2.319E+02,
     1 1.344E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,27,22, 4),I=1,6) /1.999E+02, 5.038E+01,
     1 6.114E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,27,22, 5),I=1,6) /1.619E+02, 9.318E+01,
     1 3.733E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,27,22, 6),I=1,6) /1.020E+02, 4.352E+01,
     1 3.483E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,27,23, 1),I=1,6) /7.788E+03, 1.313E+03,
     1 8.026E+00, 5.052E+01, 1.675E+00, 0.000E+00/
      DATA (PH1(I,27,23, 2),I=1,6) /9.918E+02, 1.840E+02,
     1 7.606E+00, 2.298E+01, 4.903E+00, 0.000E+00/
      DATA (PH1(I,27,23, 3),I=1,6) /8.775E+02, 2.077E+02,
     1 1.666E+02, 5.503E+01, 3.993E+00, 8.505E-06/
      DATA (PH1(I,27,23, 4),I=1,6) /1.740E+02, 5.017E+01,
     1 5.619E+00, 1.481E+02, 4.342E+00, 0.000E+00/
      DATA (PH1(I,27,23, 5),I=1,6) /1.366E+02, 8.377E+01,
     1 4.275E+01, 2.773E+01, 6.744E+00, 2.457E-01/
      DATA (PH1(I,27,23, 6),I=1,6) /7.950E+01, 4.369E+01,
     1 3.792E+02, 4.396E+01, 6.895E+00, 1.079E-04/
      DATA (PH1(I,27,24, 1),I=1,6) /7.773E+03, 1.365E+03,
     1 7.414E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,27,24, 2),I=1,6) /9.708E+02, 2.166E+02,
     1 6.563E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,27,24, 3),I=1,6) /8.556E+02, 2.326E+02,
     1 1.332E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,27,24, 4),I=1,6) /1.495E+02, 5.342E+01,
     1 5.159E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,27,24, 5),I=1,6) /1.128E+02, 9.621E+01,
     1 3.175E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,27,24, 6),I=1,6) /5.127E+01, 4.530E+01,
     1 3.755E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,27,25, 1),I=1,6) /7.758E+03, 1.487E+03,
     1 6.149E+00, 5.282E+01, 1.571E+00, 0.000E+00/
      DATA (PH1(I,27,25, 2),I=1,6) /9.547E+02, 2.172E+02,
     1 6.667E+00, 2.353E+01, 4.602E+00, 0.000E+00/
      DATA (PH1(I,27,25, 3),I=1,6) /8.295E+02, 2.089E+02,
     1 1.651E+02, 5.457E+01, 3.992E+00, 8.474E-06/
      DATA (PH1(I,27,25, 4),I=1,6) /1.305E+02, 5.435E+01,
     1 4.874E+00, 1.379E+02, 4.336E+00, 0.000E+00/
      DATA (PH1(I,27,25, 5),I=1,6) /9.362E+01, 8.419E+01,
     1 4.192E+01, 2.297E+01, 7.043E+00, 2.547E-01/
      DATA (PH1(I,27,25, 6),I=1,6) /3.350E+01, 2.181E+01,
     1 1.205E+03, 1.434E+01, 1.114E+01, 1.266E-02/
      DATA (PH1(I,27,26, 1),I=1,6) /7.742E+03, 2.079E+03,
     1 3.020E+00, 1.067E+03, 9.653E-01, 0.000E+00/
      DATA (PH1(I,27,26, 2),I=1,6) /9.443E+02, 1.709E+02,
     1 8.001E+00, 2.274E+01, 5.044E+00, 0.000E+00/
      DATA (PH1(I,27,26, 3),I=1,6) /8.130E+02, 1.989E+02,
     1 1.810E+02, 5.274E+01, 4.078E+00, 1.431E-04/
      DATA (PH1(I,27,26, 4),I=1,6) /1.213E+02, 5.602E+01,
     1 4.666E+00, 1.358E+02, 4.323E+00, 0.000E+00/
      DATA (PH1(I,27,26, 5),I=1,6) /7.621E+01, 8.408E+01,
     1 4.300E+01, 2.074E+01, 7.214E+00, 2.553E-01/
      DATA (PH1(I,27,26, 6),I=1,6) /1.708E+01, 2.419E+00,
     1 1.275E+01, 2.823E+01, 1.843E+01, 4.806E-05/
      DATA (PH1(I,27,27, 1),I=1,6) /7.725E+03, 6.269E+02,
     1 3.582E+01, 3.161E+01, 2.476E+00, 0.000E+00/
      DATA (PH1(I,27,27, 2),I=1,6) /9.400E+02, 8.888E+01,
     1 1.030E+01, 2.797E+01, 5.913E+00, 0.000E+00/
      DATA (PH1(I,27,27, 3),I=1,6) /8.000E+02, 2.832E+02,
     1 9.075E+01, 7.686E+01, 3.416E+00, 4.833E-05/
      DATA (PH1(I,27,27, 4),I=1,6) /1.150E+02, 5.097E+01,
     1 4.896E+00, 1.198E+02, 4.513E+00, 0.000E+00/
      DATA (PH1(I,27,27, 5),I=1,6) /7.300E+01, 8.256E+01,
     1 4.587E+01, 1.973E+01, 7.331E+00, 2.573E-01/
      DATA (PH1(I,27,27, 6),I=1,6) /1.580E+01, 1.581E+01,
     1 4.931E+03, 6.607E+00, 1.532E+01, 3.676E-01/
      DATA (PH1(I,27,27, 7),I=1,6) /7.864E+00, 1.370E+01,
     1 1.555E+00, 7.559E+02, 4.337E+00, 3.355E-02/
      DATA (PH1(I,28, 1, 1),I=1,6) /1.078E+04, 3.406E+02,
     1 6.983E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,28, 2, 1),I=1,6) /1.029E+04, 1.349E+03,
     1 8.382E+00, 7.431E+01, 1.588E+00, 0.000E+00/
      DATA (PH1(I,28, 3, 1),I=1,6) /1.008E+04, 1.373E+03,
     1 8.336E+00, 7.347E+01, 1.558E+00, 0.000E+00/
      DATA (PH1(I,28, 3, 2),I=1,6) /2.399E+03, 1.314E+02,
     1 8.142E+00, 4.522E+01, 4.477E+00, 0.000E+00/
      DATA (PH1(I,28, 4, 1),I=1,6) /9.937E+03, 1.382E+03,
     1 8.149E+00, 6.782E+01, 1.580E+00, 0.000E+00/
      DATA (PH1(I,28, 4, 2),I=1,6) /2.295E+03, 1.813E+02,
     1 1.117E+01, 3.332E+01, 4.426E+00, 0.000E+00/
      DATA (PH1(I,28, 5, 1),I=1,6) /9.771E+03, 1.404E+03,
     1 7.788E+00, 6.396E+01, 1.592E+00, 0.000E+00/
      DATA (PH1(I,28, 5, 2),I=1,6) /2.184E+03, 1.649E+02,
     1 1.145E+01, 2.997E+01, 4.706E+00, 0.000E+00/
      DATA (PH1(I,28, 5, 3),I=1,6) /2.131E+03, 2.473E+02,
     1 3.198E+01, 1.014E+02, 3.435E+00, 2.062E-06/
      DATA (PH1(I,28, 6, 1),I=1,6) /9.609E+03, 1.415E+03,
     1 7.600E+00, 6.229E+01, 1.599E+00, 0.000E+00/
      DATA (PH1(I,28, 6, 2),I=1,6) /2.077E+03, 1.726E+02,
     1 1.052E+01, 2.738E+01, 4.768E+00, 0.000E+00/
      DATA (PH1(I,28, 6, 3),I=1,6) /2.011E+03, 2.139E+02,
     1 8.219E+01, 9.511E+01, 3.615E+00, 6.838E-07/
      DATA (PH1(I,28, 7, 1),I=1,6) /9.423E+03, 1.429E+03,
     1 7.520E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28, 7, 2),I=1,6) /1.965E+03, 2.177E+02,
     1 8.285E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28, 7, 3),I=1,6) /1.894E+03, 2.300E+02,
     1 9.666E+01, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,28, 8, 1),I=1,6) /9.300E+03, 1.401E+03,
     1 7.673E+00, 5.681E+01, 1.640E+00, 0.000E+00/
      DATA (PH1(I,28, 8, 2),I=1,6) /1.870E+03, 1.608E+02,
     1 9.953E+00, 2.835E+01, 4.878E+00, 0.000E+00/
      DATA (PH1(I,28, 8, 3),I=1,6) /1.756E+03, 2.140E+02,
     1 1.434E+02, 7.628E+01, 3.752E+00, 5.528E-07/
      DATA (PH1(I,28, 9, 1),I=1,6) /9.125E+03, 1.443E+03,
     1 7.283E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28, 9, 2),I=1,6) /1.760E+03, 2.233E+02,
     1 7.335E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28, 9, 3),I=1,6) /1.648E+03, 2.321E+02,
     1 1.405E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,28,10, 1),I=1,6) /9.007E+03, 1.400E+03,
     1 7.626E+00, 5.528E+01, 1.653E+00, 0.000E+00/
      DATA (PH1(I,28,10, 2),I=1,6) /1.669E+03, 1.441E+02,
     1 9.448E+00, 2.779E+01, 5.132E+00, 0.000E+00/
      DATA (PH1(I,28,10, 3),I=1,6) /1.541E+03, 2.358E+02,
     1 1.530E+02, 7.451E+01, 3.705E+00, 4.687E-07/
      DATA (PH1(I,28,11, 1),I=1,6) /8.950E+03, 1.383E+03,
     1 7.837E+00, 5.839E+01, 1.644E+00, 0.000E+00/
      DATA (PH1(I,28,11, 2),I=1,6) /1.618E+03, 1.444E+02,
     1 9.470E+00, 3.104E+01, 4.981E+00, 0.000E+00/
      DATA (PH1(I,28,11, 3),I=1,6) /1.506E+03, 2.334E+02,
     1 1.540E+02, 7.184E+01, 3.746E+00, 4.400E-07/
      DATA (PH1(I,28,11, 4),I=1,6) /6.071E+02, 1.032E+01,
     1 3.674E+01, 8.252E+02, 4.304E+00, 0.000E+00/
      DATA (PH1(I,28,12, 1),I=1,6) /8.894E+03, 1.442E+03,
     1 7.159E+00, 6.122E+01, 1.602E+00, 0.000E+00/
      DATA (PH1(I,28,12, 2),I=1,6) /1.569E+03, 1.824E+02,
     1 8.013E+00, 2.738E+01, 4.780E+00, 0.000E+00/
      DATA (PH1(I,28,12, 3),I=1,6) /1.456E+03, 2.327E+02,
     1 1.539E+02, 7.045E+01, 3.763E+00, 4.221E-07/
      DATA (PH1(I,28,12, 4),I=1,6) /5.713E+02, 2.577E+01,
     1 2.339E+01, 3.455E+02, 4.232E+00, 0.000E+00/
      DATA (PH1(I,28,13, 1),I=1,6) /8.838E+03, 1.452E+03,
     1 7.048E+00, 6.036E+01, 1.602E+00, 0.000E+00/
      DATA (PH1(I,28,13, 2),I=1,6) /1.519E+03, 1.853E+02,
     1 7.851E+00, 2.679E+01, 4.785E+00, 0.000E+00/
      DATA (PH1(I,28,13, 3),I=1,6) /1.405E+03, 2.180E+02,
     1 1.742E+02, 6.951E+01, 3.842E+00, 4.401E-07/
      DATA (PH1(I,28,13, 4),I=1,6) /5.347E+02, 2.505E+01,
     1 2.361E+01, 3.622E+02, 4.225E+00, 0.000E+00/
      DATA (PH1(I,28,13, 5),I=1,6) /4.984E+02, 1.711E+02,
     1 4.708E+00, 4.790E+01, 5.019E+00, 8.478E-01/
      DATA (PH1(I,28,14, 1),I=1,6) /8.783E+03, 1.448E+03,
     1 7.081E+00, 6.137E+01, 1.599E+00, 0.000E+00/
      DATA (PH1(I,28,14, 2),I=1,6) /1.471E+03, 1.841E+02,
     1 7.825E+00, 2.644E+01, 4.816E+00, 0.000E+00/
      DATA (PH1(I,28,14, 3),I=1,6) /1.356E+03, 2.108E+02,
     1 1.846E+02, 6.874E+01, 3.886E+00, 4.288E-07/
      DATA (PH1(I,28,14, 4),I=1,6) /4.988E+02, 2.375E+01,
     1 2.383E+01, 3.786E+02, 4.245E+00, 0.000E+00/
      DATA (PH1(I,28,14, 5),I=1,6) /4.637E+02, 1.543E+02,
     1 1.053E+01, 4.743E+01, 5.176E+00, 8.373E-01/
      DATA (PH1(I,28,15, 1),I=1,6) /8.720E+03, 1.467E+03,
     1 6.912E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28,15, 2),I=1,6) /1.419E+03, 2.287E+02,
     1 6.549E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28,15, 3),I=1,6) /1.299E+03, 2.464E+02,
     1 1.356E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,28,15, 4),I=1,6) /4.620E+02, 4.702E+01,
     1 1.037E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,28,15, 5),I=1,6) /4.302E+02, 9.587E+01,
     1 2.890E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,28,16, 1),I=1,6) /8.680E+03, 1.439E+03,
     1 7.182E+00, 6.248E+01, 1.598E+00, 0.000E+00/
      DATA (PH1(I,28,16, 2),I=1,6) /1.379E+03, 1.886E+02,
     1 7.589E+00, 2.560E+01, 4.825E+00, 0.000E+00/
      DATA (PH1(I,28,16, 3),I=1,6) /1.262E+03, 2.088E+02,
     1 1.852E+02, 6.640E+01, 3.924E+00, 4.306E-07/
      DATA (PH1(I,28,16, 4),I=1,6) /4.291E+02, 2.003E+01,
     1 2.274E+01, 3.611E+02, 4.438E+00, 0.000E+00/
      DATA (PH1(I,28,16, 5),I=1,6) /3.840E+02, 1.158E+02,
     1 2.906E+01, 4.380E+01, 5.646E+00, 7.430E-01/
      DATA (PH1(I,28,17, 1),I=1,6) /8.620E+03, 1.467E+03,
     1 6.916E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28,17, 2),I=1,6) /1.328E+03, 2.299E+02,
     1 6.430E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28,17, 3),I=1,6) /1.207E+03, 2.467E+02,
     1 1.337E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,28,17, 4),I=1,6) /3.933E+02, 4.806E+01,
     1 9.150E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,28,17, 5),I=1,6) /3.521E+02, 9.596E+01,
     1 4.382E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,28,18, 1),I=1,6) /8.584E+03, 1.440E+03,
     1 7.163E+00, 6.246E+01, 1.599E+00, 0.000E+00/
      DATA (PH1(I,28,18, 2),I=1,6) /1.293E+03, 1.925E+02,
     1 7.399E+00, 2.482E+01, 4.836E+00, 0.000E+00/
      DATA (PH1(I,28,18, 3),I=1,6) /1.174E+03, 2.087E+02,
     1 1.825E+02, 6.397E+01, 3.955E+00, 4.274E-07/
      DATA (PH1(I,28,18, 4),I=1,6) /3.620E+02, 2.901E+01,
     1 1.340E+01, 2.515E+02, 4.429E+00, 0.000E+00/
      DATA (PH1(I,28,18, 5),I=1,6) /3.210E+02, 1.005E+02,
     1 4.742E+01, 4.168E+01, 5.919E+00, 5.228E-01/
      DATA (PH1(I,28,19, 1),I=1,6) /8.542E+03, 1.439E+03,
     1 7.170E+00, 6.247E+01, 1.599E+00, 0.000E+00/
      DATA (PH1(I,28,19, 2),I=1,6) /1.251E+03, 1.931E+02,
     1 7.348E+00, 2.450E+01, 4.851E+00, 0.000E+00/
      DATA (PH1(I,28,19, 3),I=1,6) /1.132E+03, 2.109E+02,
     1 1.780E+02, 6.201E+01, 3.969E+00, 4.273E-07/
      DATA (PH1(I,28,19, 4),I=1,6) /3.290E+02, 3.448E+01,
     1 1.054E+01, 2.150E+02, 4.412E+00, 0.000E+00/
      DATA (PH1(I,28,19, 5),I=1,6) /2.877E+02, 1.033E+02,
     1 4.372E+01, 3.878E+01, 5.977E+00, 5.276E-01/
      DATA (PH1(I,28,19, 6),I=1,6) /2.246E+02, 3.839E+01,
     1 1.714E+02, 7.972E+01, 6.440E+00, 1.220E-04/
      DATA (PH1(I,28,20, 1),I=1,6) /8.503E+03, 1.439E+03,
     1 7.178E+00, 6.253E+01, 1.598E+00, 0.000E+00/
      DATA (PH1(I,28,20, 2),I=1,6) /1.211E+03, 1.935E+02,
     1 7.316E+00, 2.416E+01, 4.868E+00, 0.000E+00/
      DATA (PH1(I,28,20, 3),I=1,6) /1.092E+03, 2.189E+02,
     1 1.653E+02, 6.160E+01, 3.932E+00, 4.163E-07/
      DATA (PH1(I,28,20, 4),I=1,6) /2.972E+02, 3.985E+01,
     1 8.589E+00, 1.892E+02, 4.394E+00, 0.000E+00/
      DATA (PH1(I,28,20, 5),I=1,6) /2.561E+02, 1.023E+02,
     1 4.197E+01, 3.659E+01, 6.079E+00, 4.805E-01/
      DATA (PH1(I,28,20, 6),I=1,6) /1.930E+02, 4.212E+01,
     1 2.596E+02, 6.923E+01, 6.485E+00, 1.416E-04/
      DATA (PH1(I,28,21, 1),I=1,6) /8.472E+03, 1.466E+03,
     1 6.924E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28,21, 2),I=1,6) /1.174E+03, 2.314E+02,
     1 6.277E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28,21, 3),I=1,6) /1.051E+03, 2.476E+02,
     1 1.308E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,28,21, 4),I=1,6) /2.682E+02, 5.179E+01,
     1 6.728E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,28,21, 5),I=1,6) /2.266E+02, 9.852E+01,
     1 4.031E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,28,21, 6),I=1,6) /1.620E+02, 4.735E+01,
     1 2.988E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,28,22, 1),I=1,6) /8.452E+03, 1.430E+03,
     1 7.274E+00, 6.333E+01, 1.599E+00, 0.000E+00/
      DATA (PH1(I,28,22, 2),I=1,6) /1.139E+03, 1.947E+02,
     1 7.253E+00, 2.382E+01, 4.880E+00, 0.000E+00/
      DATA (PH1(I,28,22, 3),I=1,6) /1.019E+03, 2.193E+02,
     1 1.639E+02, 5.860E+01, 3.969E+00, 4.138E-07/
      DATA (PH1(I,28,22, 4),I=1,6) /2.378E+02, 4.756E+01,
     1 6.461E+00, 1.609E+02, 4.380E+00, 0.000E+00/
      DATA (PH1(I,28,22, 5),I=1,6) /1.971E+02, 9.385E+01,
     1 4.144E+01, 3.216E+01, 6.437E+00, 3.229E-01/
      DATA (PH1(I,28,22, 6),I=1,6) /1.330E+02, 5.059E+01,
     1 3.124E+02, 5.417E+01, 6.505E+00, 1.508E-04/
      DATA (PH1(I,28,23, 1),I=1,6) /8.434E+03, 1.466E+03,
     1 6.929E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28,23, 2),I=1,6) /1.112E+03, 2.318E+02,
     1 6.240E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28,23, 3),I=1,6) /9.886E+02, 2.482E+02,
     1 1.297E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,28,23, 4),I=1,6) /2.119E+02, 5.449E+01,
     1 5.709E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,28,23, 5),I=1,6) /1.709E+02, 1.010E+02,
     1 3.462E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,28,23, 6),I=1,6) /1.080E+02, 4.837E+01,
     1 3.797E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,28,24, 1),I=1,6) /8.418E+03, 1.463E+03,
     1 6.925E+00, 6.343E+01, 1.583E+00, 0.000E+00/
      DATA (PH1(I,28,24, 2),I=1,6) /1.077E+03, 1.956E+02,
     1 7.233E+00, 2.315E+01, 4.910E+00, 0.000E+00/
      DATA (PH1(I,28,24, 3),I=1,6) /9.576E+02, 2.215E+02,
     1 1.608E+02, 5.699E+01, 3.980E+00, 4.065E-07/
      DATA (PH1(I,28,24, 4),I=1,6) /1.849E+02, 5.389E+01,
     1 5.242E+00, 1.434E+02, 4.370E+00, 0.000E+00/
      DATA (PH1(I,28,24, 5),I=1,6) /1.447E+02, 9.036E+01,
     1 4.016E+01, 2.699E+01, 6.796E+00, 2.509E-01/
      DATA (PH1(I,28,24, 6),I=1,6) /7.610E+01, 4.694E+01,
     1 4.360E+02, 3.874E+01, 7.119E+00, 3.263E-03/
      DATA (PH1(I,28,25, 1),I=1,6) /8.402E+03, 1.466E+03,
     1 6.933E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,28,25, 2),I=1,6) /1.059E+03, 2.320E+02,
     1 6.226E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,28,25, 3),I=1,6) /9.358E+02, 2.488E+02,
     1 1.288E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,28,25, 4),I=1,6) /1.597E+02, 5.774E+01,
     1 4.849E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,28,25, 5),I=1,6) /1.201E+02, 1.043E+02,
     1 2.963E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,28,25, 6),I=1,6) /5.490E+01, 5.023E+01,
     1 3.918E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,28,26, 1),I=1,6) /8.386E+03, 1.647E+03,
     1 5.367E+00, 6.758E+01, 1.485E+00, 0.000E+00/
      DATA (PH1(I,28,26, 2),I=1,6) /1.040E+03, 2.324E+02,
     1 6.314E+00, 2.361E+01, 4.604E+00, 0.000E+00/
      DATA (PH1(I,28,26, 3),I=1,6) /9.084E+02, 2.230E+02,
     1 1.588E+02, 5.650E+01, 3.977E+00, 4.038E-07/
      DATA (PH1(I,28,26, 4),I=1,6) /1.401E+02, 5.823E+01,
     1 4.586E+00, 1.330E+02, 4.367E+00, 0.000E+00/
      DATA (PH1(I,28,26, 5),I=1,6) /1.003E+02, 9.061E+01,
     1 3.980E+01, 2.247E+01, 7.093E+00, 2.543E-01/
      DATA (PH1(I,28,26, 6),I=1,6) /3.532E+01, 2.971E+01,
     1 1.965E+03, 9.751E+00, 1.140E+01, 4.019E-01/
      DATA (PH1(I,28,27, 1),I=1,6) /8.368E+03, 2.249E+03,
     1 2.776E+00, 1.325E+03, 9.540E-01, 0.000E+00/
      DATA (PH1(I,28,27, 2),I=1,6) /1.029E+03, 1.808E+02,
     1 7.632E+00, 2.289E+01, 5.060E+00, 0.000E+00/
      DATA (PH1(I,28,27, 3),I=1,6) /8.903E+02, 2.205E+02,
     1 1.625E+02, 5.653E+01, 3.990E+00, 1.246E-04/
      DATA (PH1(I,28,27, 4),I=1,6) /1.315E+02, 6.024E+01,
     1 4.394E+00, 1.314E+02, 4.347E+00, 0.000E+00/
      DATA (PH1(I,28,27, 5),I=1,6) /8.232E+01, 9.066E+01,
     1 4.075E+01, 2.048E+01, 7.246E+00, 2.545E-01/
      DATA (PH1(I,28,27, 6),I=1,6) /1.817E+01, 1.501E+00,
     1 9.258E-01, 4.442E+01, 1.928E+01, 3.708E-02/
      DATA (PH1(I,28,28, 1),I=1,6) /8.348E+03, 7.366E+02,
     1 2.836E+01, 3.622E+01, 2.316E+00, 0.000E+00/
      DATA (PH1(I,28,28, 2),I=1,6) /1.024E+03, 1.132E+02,
     1 9.424E+00, 2.712E+01, 5.643E+00, 0.000E+00/
      DATA (PH1(I,28,28, 3),I=1,6) /8.760E+02, 3.043E+02,
     1 8.611E+01, 7.868E+01, 3.408E+00, 1.680E-05/
      DATA (PH1(I,28,28, 4),I=1,6) /1.250E+02, 5.448E+01,
     1 4.611E+00, 1.157E+02, 4.548E+00, 0.000E+00/
      DATA (PH1(I,28,28, 5),I=1,6) /8.200E+01, 8.896E+01,
     1 4.351E+01, 1.942E+01, 7.372E+00, 2.566E-01/
      DATA (PH1(I,28,28, 6),I=1,6) /1.700E+01, 6.063E+00,
     1 1.186E+03, 6.823E+00, 2.223E+01, 6.227E-03/
      DATA (PH1(I,28,28, 7),I=1,6) /7.637E+00, 1.468E+01,
     1 1.437E+00, 7.411E+02, 4.342E+00, 3.908E-02/
      DATA (PH1(I,29, 1, 1),I=1,6) /1.157E+04, 3.656E+02,
     1 6.510E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,29, 2, 1),I=1,6) /1.106E+04, 1.427E+03,
     1 8.057E+00, 7.912E+01, 1.583E+00, 0.000E+00/
      DATA (PH1(I,29, 3, 1),I=1,6) /1.085E+04, 1.435E+03,
     1 8.196E+00, 7.379E+01, 1.577E+00, 0.000E+00/
      DATA (PH1(I,29, 3, 2),I=1,6) /2.585E+03, 1.732E+02,
     1 6.336E+00, 4.099E+01, 4.316E+00, 0.000E+00/
      DATA (PH1(I,29, 4, 1),I=1,6) /1.070E+04, 1.476E+03,
     1 7.658E+00, 7.281E+01, 1.567E+00, 0.000E+00/
      DATA (PH1(I,29, 4, 2),I=1,6) /2.459E+03, 1.822E+02,
     1 1.108E+01, 3.417E+01, 4.489E+00, 0.000E+00/
      DATA (PH1(I,29, 5, 1),I=1,6) /1.053E+04, 1.516E+03,
     1 7.156E+00, 7.157E+01, 1.560E+00, 0.000E+00/
      DATA (PH1(I,29, 5, 2),I=1,6) /2.363E+03, 2.174E+02,
     1 9.074E+00, 2.892E+01, 4.445E+00, 0.000E+00/
      DATA (PH1(I,29, 5, 3),I=1,6) /2.298E+03, 1.099E+02,
     1 1.898E+02, 1.272E+02, 4.081E+00, 2.821E-08/
      DATA (PH1(I,29, 6, 1),I=1,6) /1.036E+04, 1.537E+03,
     1 6.899E+00, 7.045E+01, 1.560E+00, 0.000E+00/
      DATA (PH1(I,29, 6, 2),I=1,6) /2.253E+03, 2.118E+02,
     1 8.876E+00, 2.857E+01, 4.513E+00, 0.000E+00/
      DATA (PH1(I,29, 6, 3),I=1,6) /2.173E+03, 1.503E+02,
     1 1.866E+02, 1.055E+02, 3.944E+00, 1.972E-08/
      DATA (PH1(I,29, 7, 1),I=1,6) /1.016E+04, 1.531E+03,
     1 7.047E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29, 7, 2),I=1,6) /2.133E+03, 2.330E+02,
     1 7.838E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29, 7, 3),I=1,6) /2.045E+03, 2.473E+02,
     1 9.111E+01, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,29, 8, 1),I=1,6) /1.004E+04, 1.556E+03,
     1 6.650E+00, 6.987E+01, 1.559E+00, 0.000E+00/
      DATA (PH1(I,29, 8, 2),I=1,6) /2.037E+03, 2.003E+02,
     1 8.422E+00, 2.803E+01, 4.659E+00, 0.000E+00/
      DATA (PH1(I,29, 8, 3),I=1,6) /1.905E+03, 1.999E+02,
     1 1.798E+02, 7.715E+01, 3.888E+00, 1.496E-08/
      DATA (PH1(I,29, 9, 1),I=1,6) /9.844E+03, 1.546E+03,
     1 6.833E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29, 9, 2),I=1,6) /1.920E+03, 2.388E+02,
     1 6.970E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29, 9, 3),I=1,6) /1.793E+03, 2.493E+02,
     1 1.335E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,29,10, 1),I=1,6) /9.729E+03, 1.540E+03,
     1 6.748E+00, 6.512E+01, 1.590E+00, 0.000E+00/
      DATA (PH1(I,29,10, 2),I=1,6) /1.827E+03, 1.901E+02,
     1 7.838E+00, 2.795E+01, 4.791E+00, 0.000E+00/
      DATA (PH1(I,29,10, 3),I=1,6) /1.690E+03, 2.206E+02,
     1 1.906E+02, 6.390E+01, 3.955E+00, 1.558E-08/
      DATA (PH1(I,29,11, 1),I=1,6) /9.668E+03, 1.523E+03,
     1 6.915E+00, 6.870E+01, 1.581E+00, 0.000E+00/
      DATA (PH1(I,29,11, 2),I=1,6) /1.773E+03, 1.626E+02,
     1 8.723E+00, 3.155E+01, 4.881E+00, 0.000E+00/
      DATA (PH1(I,29,11, 3),I=1,6) /1.657E+03, 2.651E+02,
     1 1.307E+02, 7.793E+01, 3.637E+00, 1.894E-08/
      DATA (PH1(I,29,11, 4),I=1,6) /6.706E+02, 1.529E+01,
     1 2.453E+01, 6.172E+02, 4.260E+00, 0.000E+00/
      DATA (PH1(I,29,12, 1),I=1,6) /9.610E+03, 1.558E+03,
     1 6.567E+00, 6.667E+01, 1.576E+00, 0.000E+00/
      DATA (PH1(I,29,12, 2),I=1,6) /1.722E+03, 1.953E+02,
     1 7.587E+00, 2.834E+01, 4.738E+00, 0.000E+00/
      DATA (PH1(I,29,12, 3),I=1,6) /1.604E+03, 2.619E+02,
     1 1.331E+02, 7.654E+01, 3.661E+00, 1.879E-08/
      DATA (PH1(I,29,12, 4),I=1,6) /6.330E+02, 2.604E+01,
     1 2.426E+01, 3.618E+02, 4.240E+00, 0.000E+00/
      DATA (PH1(I,29,13, 1),I=1,6) /9.550E+03, 1.563E+03,
     1 6.515E+00, 6.616E+01, 1.577E+00, 0.000E+00/
      DATA (PH1(I,29,13, 2),I=1,6) /1.670E+03, 1.993E+02,
     1 7.415E+00, 2.757E+01, 4.743E+00, 0.000E+00/
      DATA (PH1(I,29,13, 3),I=1,6) /1.551E+03, 2.401E+02,
     1 1.577E+02, 7.468E+01, 3.767E+00, 2.341E-08/
      DATA (PH1(I,29,13, 4),I=1,6) /5.949E+02, 3.706E+01,
     1 1.653E+01, 2.918E+02, 4.111E+00, 0.000E+00/
      DATA (PH1(I,29,13, 5),I=1,6) /5.570E+02, 2.011E+02,
     1 3.929E+00, 5.061E+01, 4.854E+00, 8.398E-01/
      DATA (PH1(I,29,14, 1),I=1,6) /9.493E+03, 1.558E+03,
     1 6.561E+00, 6.691E+01, 1.575E+00, 0.000E+00/
      DATA (PH1(I,29,14, 2),I=1,6) /1.619E+03, 1.974E+02,
     1 7.417E+00, 2.691E+01, 4.791E+00, 0.000E+00/
      DATA (PH1(I,29,14, 3),I=1,6) /1.499E+03, 2.384E+02,
     1 1.588E+02, 7.405E+01, 3.780E+00, 2.405E-08/
      DATA (PH1(I,29,14, 4),I=1,6) /5.572E+02, 3.531E+01,
     1 1.668E+01, 3.042E+02, 4.130E+00, 0.000E+00/
      DATA (PH1(I,29,14, 5),I=1,6) /5.200E+02, 1.207E+02,
     1 1.569E+01, 5.849E+01, 5.405E+00, 4.555E-01/
      DATA (PH1(I,29,15, 1),I=1,6) /9.423E+03, 1.577E+03,
     1 6.419E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,15, 2),I=1,6) /1.563E+03, 2.446E+02,
     1 6.221E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,15, 3),I=1,6) /1.439E+03, 2.632E+02,
     1 1.311E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,15, 4),I=1,6) /5.184E+02, 5.001E+01,
     1 1.019E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,15, 5),I=1,6) /4.840E+02, 1.032E+02,
     1 2.783E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,16, 1),I=1,6) /9.384E+03, 1.553E+03,
     1 6.609E+00, 6.923E+01, 1.568E+00, 0.000E+00/
      DATA (PH1(I,29,16, 2),I=1,6) /1.522E+03, 1.990E+02,
     1 7.271E+00, 2.604E+01, 4.826E+00, 0.000E+00/
      DATA (PH1(I,29,16, 3),I=1,6) /1.400E+03, 2.302E+02,
     1 1.675E+02, 7.119E+01, 3.847E+00, 2.441E-08/
      DATA (PH1(I,29,16, 4),I=1,6) /4.838E+02, 3.424E+01,
     1 1.441E+01, 2.579E+02, 4.283E+00, 0.000E+00/
      DATA (PH1(I,29,16, 5),I=1,6) /4.350E+02, 1.064E+02,
     1 3.399E+01, 5.025E+01, 5.723E+00, 4.731E-01/
      DATA (PH1(I,29,17, 1),I=1,6) /9.317E+03, 1.576E+03,
     1 6.433E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,17, 2),I=1,6) /1.467E+03, 2.457E+02,
     1 6.110E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,17, 3),I=1,6) /1.340E+03, 2.635E+02,
     1 1.295E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,17, 4),I=1,6) /4.458E+02, 5.099E+01,
     1 9.072E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,17, 5),I=1,6) /4.010E+02, 1.031E+02,
     1 4.266E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,18, 1),I=1,6) /9.282E+03, 1.550E+03,
     1 6.633E+00, 6.922E+01, 1.570E+00, 0.000E+00/
      DATA (PH1(I,29,18, 2),I=1,6) /1.431E+03, 1.960E+02,
     1 7.274E+00, 2.502E+01, 4.908E+00, 0.000E+00/
      DATA (PH1(I,29,18, 3),I=1,6) /1.307E+03, 2.291E+02,
     1 1.666E+02, 6.860E+01, 3.881E+00, 2.421E-08/
      DATA (PH1(I,29,18, 4),I=1,6) /4.129E+02, 3.356E+01,
     1 1.254E+01, 2.441E+02, 4.373E+00, 0.000E+00/
      DATA (PH1(I,29,18, 5),I=1,6) /3.688E+02, 1.023E+02,
     1 4.888E+01, 4.483E+01, 5.919E+00, 3.996E-01/
      DATA (PH1(I,29,19, 1),I=1,6) /9.237E+03, 1.550E+03,
     1 6.630E+00, 6.917E+01, 1.570E+00, 0.000E+00/
      DATA (PH1(I,29,19, 2),I=1,6) /1.386E+03, 1.983E+02,
     1 7.172E+00, 2.486E+01, 4.901E+00, 0.000E+00/
      DATA (PH1(I,29,19, 3),I=1,6) /1.262E+03, 2.322E+02,
     1 1.615E+02, 6.669E+01, 3.889E+00, 2.403E-08/
      DATA (PH1(I,29,19, 4),I=1,6) /3.778E+02, 3.308E+01,
     1 1.152E+01, 2.347E+02, 4.436E+00, 0.000E+00/
      DATA (PH1(I,29,19, 5),I=1,6) /3.335E+02, 1.017E+02,
     1 4.646E+01, 4.207E+01, 6.020E+00, 3.632E-01/
      DATA (PH1(I,29,19, 6),I=1,6) /2.661E+02, 4.605E+01,
     1 1.458E+02, 7.967E+01, 6.285E+00, 4.503E-03/
      DATA (PH1(I,29,20, 1),I=1,6) /9.181E+03, 1.574E+03,
     1 6.454E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,20, 2),I=1,6) /1.339E+03, 2.470E+02,
     1 5.993E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,20, 3),I=1,6) /1.211E+03, 2.640E+02,
     1 1.274E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,20, 4),I=1,6) /3.447E+02, 5.350E+01,
     1 7.318E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,20, 5),I=1,6) /2.999E+02, 1.045E+02,
     1 4.286E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,20, 6),I=1,6) /2.320E+02, 5.210E+01,
     1 2.093E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,29,21, 1),I=1,6) /9.155E+03, 1.546E+03,
     1 6.671E+00, 6.964E+01, 1.570E+00, 0.000E+00/
      DATA (PH1(I,29,21, 2),I=1,6) /1.302E+03, 1.958E+02,
     1 7.186E+00, 2.423E+01, 4.960E+00, 0.000E+00/
      DATA (PH1(I,29,21, 3),I=1,6) /1.178E+03, 2.383E+02,
     1 1.525E+02, 6.417E+01, 3.890E+00, 2.517E-08/
      DATA (PH1(I,29,21, 4),I=1,6) /3.113E+02, 4.346E+01,
     1 7.773E+00, 1.823E+02, 4.413E+00, 0.000E+00/
      DATA (PH1(I,29,21, 5),I=1,6) /2.672E+02, 9.821E+01,
     1 4.314E+01, 3.667E+01, 6.291E+00, 2.759E-01/
      DATA (PH1(I,29,21, 6),I=1,6) /1.990E+02, 5.000E+01,
     1 3.041E+02, 6.355E+01, 6.458E+00, 2.040E-03/
      DATA (PH1(I,29,22, 1),I=1,6) /9.127E+03, 1.573E+03,
     1 6.468E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,22, 2),I=1,6) /1.267E+03, 2.476E+02,
     1 5.944E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,22, 3),I=1,6) /1.137E+03, 2.644E+02,
     1 1.263E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,22, 4),I=1,6) /2.824E+02, 5.587E+01,
     1 6.260E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,22, 5),I=1,6) /2.375E+02, 1.064E+02,
     1 3.728E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,22, 6),I=1,6) /1.670E+02, 5.236E+01,
     1 3.470E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,29,23, 1),I=1,6) /9.107E+03, 1.541E+03,
     1 6.731E+00, 7.196E+01, 1.562E+00, 0.000E+00/
      DATA (PH1(I,29,23, 2),I=1,6) /1.229E+03, 1.968E+02,
     1 7.138E+00, 2.379E+01, 4.979E+00, 0.000E+00/
      DATA (PH1(I,29,23, 3),I=1,6) /1.104E+03, 2.376E+02,
     1 1.529E+02, 6.163E+01, 3.925E+00, 2.511E-08/
      DATA (PH1(I,29,23, 4),I=1,6) /2.504E+02, 5.085E+01,
     1 5.992E+00, 1.563E+02, 4.410E+00, 0.000E+00/
      DATA (PH1(I,29,23, 5),I=1,6) /2.067E+02, 1.010E+02,
     1 3.860E+01, 3.121E+01, 6.491E+00, 3.331E-01/
      DATA (PH1(I,29,23, 6),I=1,6) /1.390E+02, 5.652E+01,
     1 3.452E+02, 5.007E+01, 6.558E+00, 2.171E-03/
      DATA (PH1(I,29,24, 1),I=1,6) /9.088E+03, 1.572E+03,
     1 6.482E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,24, 2),I=1,6) /1.203E+03, 2.479E+02,
     1 5.918E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,24, 3),I=1,6) /1.073E+03, 2.649E+02,
     1 1.255E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,24, 4),I=1,6) /2.243E+02, 5.878E+01,
     1 5.342E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,24, 5),I=1,6) /1.801E+02, 1.092E+02,
     1 3.217E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,24, 6),I=1,6) /1.030E+02, 5.346E+01,
     1 4.010E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,29,25, 1),I=1,6) /9.070E+03, 1.572E+03,
     1 6.440E+00, 6.913E+01, 1.560E+00, 0.000E+00/
      DATA (PH1(I,29,25, 2),I=1,6) /1.166E+03, 2.072E+02,
     1 6.884E+00, 2.337E+01, 4.917E+00, 0.000E+00/
      DATA (PH1(I,29,25, 3),I=1,6) /1.041E+03, 2.368E+02,
     1 1.537E+02, 5.951E+01, 3.955E+00, 2.525E-08/
      DATA (PH1(I,29,25, 4),I=1,6) /1.960E+02, 5.776E+01,
     1 4.896E+00, 1.391E+02, 4.397E+00, 0.000E+00/
      DATA (PH1(I,29,25, 5),I=1,6) /1.528E+02, 9.650E+01,
     1 3.804E+01, 2.618E+01, 6.870E+00, 2.460E-01/
      DATA (PH1(I,29,25, 6),I=1,6) /7.990E+01, 5.875E+01,
     1 3.868E+02, 3.774E+01, 6.878E+00, 3.091E-01/
      DATA (PH1(I,29,26, 1),I=1,6) /9.052E+03, 1.571E+03,
     1 6.497E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,29,26, 2),I=1,6) /1.146E+03, 2.480E+02,
     1 5.914E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,29,26, 3),I=1,6) /1.020E+03, 2.654E+02,
     1 1.248E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,29,26, 4),I=1,6) /1.703E+02, 6.225E+01,
     1 4.565E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,29,26, 5),I=1,6) /1.277E+02, 1.127E+02,
     1 2.769E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,29,26, 6),I=1,6) /5.738E+01, 5.539E+01,
     1 4.028E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,29,27, 1),I=1,6) /9.032E+03, 2.025E+03,
     1 3.761E+00, 2.259E+02, 1.184E+00, 0.000E+00/
      DATA (PH1(I,29,27, 2),I=1,6) /1.128E+03, 2.508E+02,
     1 5.929E+00, 2.384E+01, 4.581E+00, 0.000E+00/
      DATA (PH1(I,29,27, 3),I=1,6) /9.911E+02, 2.369E+02,
     1 1.537E+02, 5.847E+01, 3.967E+00, 2.525E-08/
      DATA (PH1(I,29,27, 4),I=1,6) /1.499E+02, 6.247E+01,
     1 4.311E+00, 1.291E+02, 4.390E+00, 0.000E+00/
      DATA (PH1(I,29,27, 5),I=1,6) /1.072E+02, 9.745E+01,
     1 3.768E+01, 2.211E+01, 7.130E+00, 2.537E-01/
      DATA (PH1(I,29,27, 6),I=1,6) /3.684E+01, 3.541E+01,
     1 1.181E+03, 1.390E+01, 1.017E+01, 3.416E-01/
      DATA (PH1(I,29,28, 1),I=1,6) /9.012E+03, 2.378E+03,
     1 2.659E+00, 7.489E+02, 9.937E-01, 0.000E+00/
      DATA (PH1(I,29,28, 2),I=1,6) /1.114E+03, 3.054E+02,
     1 4.889E+00, 2.675E+01, 4.155E+00, 0.000E+00/
      DATA (PH1(I,29,28, 3),I=1,6) /9.715E+02, 2.328E+02,
     1 1.588E+02, 5.776E+01, 3.997E+00, 2.543E-08/
      DATA (PH1(I,29,28, 4),I=1,6) /1.312E+02, 6.472E+01,
     1 4.139E+00, 1.278E+02, 4.366E+00, 0.000E+00/
      DATA (PH1(I,29,28, 5),I=1,6) /8.861E+01, 9.750E+01,
     1 3.865E+01, 2.022E+01, 7.277E+00, 2.537E-01/
      DATA (PH1(I,29,28, 6),I=1,6) /2.029E+01, 1.590E+01,
     1 3.992E+03, 7.098E+00, 1.616E+01, 7.439E-01/
      DATA (PH1(I,29,29, 1),I=1,6) /8.988E+03, 1.788E+03,
     1 4.870E+00, 7.645E+01, 1.451E+00, 0.000E+00/
      DATA (PH1(I,29,29, 2),I=1,6) /1.106E+03, 2.502E+02,
     1 5.938E+00, 2.402E+01, 4.576E+00, 0.000E+00/
      DATA (PH1(I,29,29, 3),I=1,6) /9.470E+02, 2.826E+02,
     1 1.093E+02, 6.688E+01, 3.668E+00, 9.174E-07/
      DATA (PH1(I,29,29, 4),I=1,6) /1.288E+02, 6.198E+01,
     1 4.164E+00, 1.158E+02, 4.488E+00, 0.000E+00/
      DATA (PH1(I,29,29, 5),I=1,6) /8.300E+01, 9.617E+01,
     1 4.275E+01, 1.747E+01, 7.555E+00, 2.599E-01/
      DATA (PH1(I,29,29, 6),I=1,6) /1.064E+01, 7.279E+00,
     1 1.027E+03, 7.988E+00, 2.033E+01, 1.582E+00/
      DATA (PH1(I,29,29, 7),I=1,6) /7.726E+00, 1.436E+01,
     1 4.681E-01, 2.383E+03, 4.224E+00, 3.736E-02/
      DATA (PH1(I,30, 1, 1),I=1,6) /1.239E+04, 3.915E+02,
     1 6.083E+01, 3.288E+01, 2.963E+00, 0.000E+00/
      DATA (PH1(I,30, 2, 1),I=1,6) /1.187E+04, 1.564E+03,
     1 7.165E+00, 9.095E+01, 1.535E+00, 0.000E+00/
      DATA (PH1(I,30, 3, 1),I=1,6) /1.165E+04, 1.600E+03,
     1 7.032E+00, 9.272E+01, 1.497E+00, 0.000E+00/
      DATA (PH1(I,30, 3, 2),I=1,6) /2.780E+03, 3.219E+02,
     1 3.342E+00, 3.904E+01, 3.684E+00, 0.000E+00/
      DATA (PH1(I,30, 4, 1),I=1,6) /1.149E+04, 1.565E+03,
     1 7.294E+00, 7.835E+01, 1.557E+00, 0.000E+00/
      DATA (PH1(I,30, 4, 2),I=1,6) /2.647E+03, 1.989E+02,
     1 1.025E+01, 3.502E+01, 4.434E+00, 0.000E+00/
      DATA (PH1(I,30, 5, 1),I=1,6) /1.131E+04, 1.582E+03,
     1 7.052E+00, 7.246E+01, 1.577E+00, 0.000E+00/
      DATA (PH1(I,30, 5, 2),I=1,6) /2.550E+03, 2.168E+02,
     1 9.068E+00, 2.977E+01, 4.514E+00, 0.000E+00/
      DATA (PH1(I,30, 5, 3),I=1,6) /2.479E+03, 1.230E+02,
     1 1.663E+02, 1.330E+02, 4.013E+00, 2.250E-08/
      DATA (PH1(I,30, 6, 1),I=1,6) /1.113E+04, 1.607E+03,
     1 6.765E+00, 7.173E+01, 1.574E+00, 0.000E+00/
      DATA (PH1(I,30, 6, 2),I=1,6) /2.435E+03, 2.335E+02,
     1 8.158E+00, 2.769E+01, 4.506E+00, 0.000E+00/
      DATA (PH1(I,30, 6, 3),I=1,6) /2.363E+03, 1.386E+02,
     1 2.396E+02, 1.169E+02, 4.016E+00, 1.579E-08/
      DATA (PH1(I,30, 7, 1),I=1,6) /1.092E+04, 1.637E+03,
     1 6.617E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30, 7, 2),I=1,6) /2.309E+03, 2.488E+02,
     1 7.424E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30, 7, 3),I=1,6) /2.216E+03, 2.651E+02,
     1 8.601E+01, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,30, 8, 1),I=1,6) /1.080E+04, 1.629E+03,
     1 6.509E+00, 7.085E+01, 1.574E+00, 0.000E+00/
      DATA (PH1(I,30, 8, 2),I=1,6) /2.210E+03, 2.179E+02,
     1 7.856E+00, 2.852E+01, 4.612E+00, 0.000E+00/
      DATA (PH1(I,30, 8, 3),I=1,6) /2.070E+03, 2.034E+02,
     1 1.893E+02, 8.099E+01, 3.908E+00, 6.713E-09/
      DATA (PH1(I,30, 9, 1),I=1,6) /1.059E+04, 1.652E+03,
     1 6.423E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30, 9, 2),I=1,6) /2.087E+03, 2.548E+02,
     1 6.629E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30, 9, 3),I=1,6) /1.953E+03, 2.670E+02,
     1 1.269E+02, 5.000E+01, 4.000E+00, 3.000E-01/
      DATA (PH1(I,30,10, 1),I=1,6) /1.048E+04, 1.636E+03,
     1 6.400E+00, 6.864E+01, 1.584E+00, 0.000E+00/
      DATA (PH1(I,30,10, 2),I=1,6) /1.992E+03, 2.104E+02,
     1 7.257E+00, 2.813E+01, 4.728E+00, 0.000E+00/
      DATA (PH1(I,30,10, 3),I=1,6) /1.846E+03, 2.254E+02,
     1 1.992E+02, 6.645E+01, 3.978E+00, 6.413E-09/
      DATA (PH1(I,30,11, 1),I=1,6) /1.042E+04, 1.605E+03,
     1 6.675E+00, 7.053E+01, 1.588E+00, 0.000E+00/
      DATA (PH1(I,30,11, 2),I=1,6) /1.936E+03, 1.953E+02,
     1 7.654E+00, 3.173E+01, 4.697E+00, 0.000E+00/
      DATA (PH1(I,30,11, 3),I=1,6) /1.815E+03, 2.273E+02,
     1 1.928E+02, 6.347E+01, 4.011E+00, 6.017E-09/
      DATA (PH1(I,30,11, 4),I=1,6) /7.374E+02, 1.960E+01,
     1 1.939E+01, 5.293E+02, 4.225E+00, 0.000E+00/
      DATA (PH1(I,30,12, 1),I=1,6) /1.035E+04, 1.631E+03,
     1 6.432E+00, 6.791E+01, 1.590E+00, 0.000E+00/
      DATA (PH1(I,30,12, 2),I=1,6) /1.882E+03, 2.070E+02,
     1 7.256E+00, 2.865E+01, 4.736E+00, 0.000E+00/
      DATA (PH1(I,30,12, 3),I=1,6) /1.760E+03, 2.265E+02,
     1 1.930E+02, 6.308E+01, 4.019E+00, 6.018E-09/
      DATA (PH1(I,30,12, 4),I=1,6) /6.980E+02, 3.115E+01,
     1 2.108E+01, 3.355E+02, 4.197E+00, 0.000E+00/
      DATA (PH1(I,30,13, 1),I=1,6) /1.029E+04, 1.636E+03,
     1 6.381E+00, 6.693E+01, 1.592E+00, 0.000E+00/
      DATA (PH1(I,30,13, 2),I=1,6) /1.828E+03, 2.085E+02,
     1 7.155E+00, 2.796E+01, 4.758E+00, 0.000E+00/
      DATA (PH1(I,30,13, 3),I=1,6) /1.704E+03, 2.296E+02,
     1 1.871E+02, 6.376E+01, 3.995E+00, 5.868E-09/
      DATA (PH1(I,30,13, 4),I=1,6) /6.583E+02, 3.116E+01,
     1 2.073E+01, 3.448E+02, 4.181E+00, 0.000E+00/
      DATA (PH1(I,30,13, 5),I=1,6) /6.190E+02, 2.539E+02,
     1 2.886E+00, 5.384E+01, 4.592E+00, 8.255E-01/
      DATA (PH1(I,30,14, 1),I=1,6) /1.023E+04, 1.633E+03,
     1 6.399E+00, 6.753E+01, 1.591E+00, 0.000E+00/
      DATA (PH1(I,30,14, 2),I=1,6) /1.775E+03, 2.101E+02,
     1 7.058E+00, 2.748E+01, 4.771E+00, 0.000E+00/
      DATA (PH1(I,30,14, 3),I=1,6) /1.650E+03, 2.720E+02,
     1 1.333E+02, 8.071E+01, 3.661E+00, 4.757E-09/
      DATA (PH1(I,30,14, 4),I=1,6) /6.187E+02, 2.985E+01,
     1 2.083E+01, 3.577E+02, 4.196E+00, 0.000E+00/
      DATA (PH1(I,30,14, 5),I=1,6) /5.790E+02, 1.403E+02,
     1 1.361E+01, 6.070E+01, 5.261E+00, 5.034E-01/
      DATA (PH1(I,30,15, 1),I=1,6) /1.015E+04, 1.691E+03,
     1 5.972E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30,15, 2),I=1,6) /1.716E+03, 2.609E+02,
     1 5.912E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30,15, 3),I=1,6) /1.586E+03, 2.806E+02,
     1 1.269E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,30,15, 4),I=1,6) /5.780E+02, 5.311E+01,
     1 1.000E+01, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,30,15, 5),I=1,6) /5.420E+02, 1.108E+02,
     1 2.676E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,30,16, 1),I=1,6) /1.012E+04, 1.617E+03,
     1 6.536E+00, 6.839E+01, 1.594E+00, 0.000E+00/
      DATA (PH1(I,30,16, 2),I=1,6) /1.673E+03, 2.142E+02,
     1 6.875E+00, 2.631E+01, 4.800E+00, 0.000E+00/
      DATA (PH1(I,30,16, 3),I=1,6) /1.546E+03, 2.596E+02,
     1 1.443E+02, 7.722E+01, 3.739E+00, 5.496E-09/
      DATA (PH1(I,30,16, 4),I=1,6) /5.416E+02, 3.912E+01,
     1 1.395E+01, 2.883E+02, 4.155E+00, 0.000E+00/
      DATA (PH1(I,30,16, 5),I=1,6) /4.900E+02, 1.153E+02,
     1 3.254E+01, 5.308E+01, 5.645E+00, 4.056E-01/
      DATA (PH1(I,30,17, 1),I=1,6) /1.004E+04, 1.689E+03,
     1 5.993E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30,17, 2),I=1,6) /1.613E+03, 2.621E+02,
     1 5.810E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30,17, 3),I=1,6) /1.481E+03, 2.807E+02,
     1 1.254E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,30,17, 4),I=1,6) /5.015E+02, 5.403E+01,
     1 8.971E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,30,17, 5),I=1,6) /4.540E+02, 1.106E+02,
     1 4.142E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,30,18, 1),I=1,6) /1.001E+04, 1.617E+03,
     1 6.532E+00, 6.768E+01, 1.597E+00, 0.000E+00/
      DATA (PH1(I,30,18, 2),I=1,6) /1.576E+03, 2.135E+02,
     1 6.802E+00, 2.575E+01, 4.839E+00, 0.000E+00/
      DATA (PH1(I,30,18, 3),I=1,6) /1.448E+03, 2.538E+02,
     1 1.487E+02, 7.400E+01, 3.794E+00, 5.551E-09/
      DATA (PH1(I,30,18, 4),I=1,6) /4.671E+02, 3.780E+01,
     1 1.186E+01, 2.389E+02, 4.334E+00, 0.000E+00/
      DATA (PH1(I,30,18, 5),I=1,6) /4.197E+02, 1.073E+02,
     1 4.863E+01, 4.776E+01, 5.875E+00, 2.432E-01/
      DATA (PH1(I,30,19, 1),I=1,6) /9.960E+03, 1.610E+03,
     1 6.593E+00, 6.774E+01, 1.600E+00, 0.000E+00/
      DATA (PH1(I,30,19, 2),I=1,6) /1.528E+03, 2.118E+02,
     1 6.804E+00, 2.536E+01, 4.875E+00, 0.000E+00/
      DATA (PH1(I,30,19, 3),I=1,6) /1.399E+03, 2.529E+02,
     1 1.489E+02, 7.156E+01, 3.823E+00, 5.330E-09/
      DATA (PH1(I,30,19, 4),I=1,6) /4.299E+02, 3.741E+01,
     1 1.098E+01, 2.313E+02, 4.386E+00, 0.000E+00/
      DATA (PH1(I,30,19, 5),I=1,6) /3.828E+02, 1.057E+02,
     1 4.664E+01, 4.477E+01, 5.995E+00, 2.234E-01/
      DATA (PH1(I,30,19, 6),I=1,6) /3.108E+02, 5.749E+01,
     1 1.122E+02, 7.804E+01, 6.097E+00, 3.731E-01/
      DATA (PH1(I,30,20, 1),I=1,6) /9.915E+03, 1.611E+03,
     1 6.591E+00, 6.799E+01, 1.598E+00, 0.000E+00/
      DATA (PH1(I,30,20, 2),I=1,6) /1.482E+03, 2.023E+02,
     1 6.987E+00, 2.500E+01, 4.972E+00, 0.000E+00/
      DATA (PH1(I,30,20, 3),I=1,6) /1.353E+03, 2.543E+02,
     1 1.465E+02, 6.961E+01, 3.838E+00, 5.027E-09/
      DATA (PH1(I,30,20, 4),I=1,6) /3.939E+02, 3.705E+01,
     1 1.008E+01, 2.216E+02, 4.447E+00, 0.000E+00/
      DATA (PH1(I,30,20, 5),I=1,6) /3.467E+02, 1.046E+02,
     1 4.454E+01, 4.180E+01, 6.115E+00, 1.997E-01/
      DATA (PH1(I,30,20, 6),I=1,6) /2.740E+02, 5.735E+01,
     1 2.033E+02, 7.143E+01, 6.214E+00, 3.888E-01/
      DATA (PH1(I,30,21, 1),I=1,6) /9.854E+03, 1.685E+03,
     1 6.036E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30,21, 2),I=1,6) /1.436E+03, 2.638E+02,
     1 5.673E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30,21, 3),I=1,6) /1.302E+03, 2.813E+02,
     1 1.231E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,30,21, 4),I=1,6) /3.609E+02, 5.754E+01,
     1 6.797E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,30,21, 5),I=1,6) /3.128E+02, 1.126E+02,
     1 3.961E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,30,21, 6),I=1,6) /2.380E+02, 5.727E+01,
     1 2.777E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,30,22, 1),I=1,6) /9.831E+03, 1.605E+03,
     1 6.653E+00, 6.874E+01, 1.597E+00, 0.000E+00/
      DATA (PH1(I,30,22, 2),I=1,6) /1.397E+03, 2.114E+02,
     1 6.748E+00, 2.443E+01, 4.934E+00, 0.000E+00/
      DATA (PH1(I,30,22, 3),I=1,6) /1.267E+03, 2.512E+02,
     1 1.490E+02, 6.562E+01, 3.896E+00, 5.943E-09/
      DATA (PH1(I,30,22, 4),I=1,6) /3.257E+02, 4.721E+01,
     1 7.066E+00, 1.761E+02, 4.431E+00, 0.000E+00/
      DATA (PH1(I,30,22, 5),I=1,6) /2.785E+02, 1.016E+02,
     1 4.118E+01, 3.606E+01, 6.396E+00, 1.450E-01/
      DATA (PH1(I,30,22, 6),I=1,6) /2.030E+02, 5.906E+01,
     1 3.187E+02, 5.850E+01, 6.430E+00, 2.472E-01/
      DATA (PH1(I,30,23, 1),I=1,6) /9.809E+03, 1.683E+03,
     1 6.057E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30,23, 2),I=1,6) /1.362E+03, 2.642E+02,
     1 5.637E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30,23, 3),I=1,6) /1.227E+03, 2.816E+02,
     1 1.222E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,30,23, 4),I=1,6) /2.969E+02, 6.012E+01,
     1 5.837E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,30,23, 5),I=1,6) /2.487E+02, 1.147E+02,
     1 3.456E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,30,23, 6),I=1,6) /1.750E+02, 5.760E+01,
     1 3.804E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,30,24, 1),I=1,6) /9.788E+03, 1.697E+03,
     1 5.935E+00, 9.066E+01, 1.488E+00, 0.000E+00/
      DATA (PH1(I,30,24, 2),I=1,6) /1.323E+03, 2.094E+02,
     1 6.772E+00, 2.401E+01, 4.975E+00, 0.000E+00/
      DATA (PH1(I,30,24, 3),I=1,6) /1.193E+03, 2.522E+02,
     1 1.476E+02, 6.376E+01, 3.914E+00, 6.079E-09/
      DATA (PH1(I,30,24, 4),I=1,6) /2.633E+02, 5.492E+01,
     1 5.518E+00, 1.512E+02, 4.431E+00, 0.000E+00/
      DATA (PH1(I,30,24, 5),I=1,6) /2.165E+02, 1.096E+02,
     1 3.576E+01, 3.034E+01, 6.524E+00, 3.502E-01/
      DATA (PH1(I,30,24, 6),I=1,6) /1.360E+02, 6.530E+01,
     1 3.456E+02, 4.608E+01, 6.555E+00, 2.476E-01/
      DATA (PH1(I,30,25, 1),I=1,6) /9.770E+03, 1.682E+03,
     1 6.079E+00, 5.000E+01, 1.650E+00, 0.000E+00/
      DATA (PH1(I,30,25, 2),I=1,6) /1.298E+03, 2.645E+02,
     1 5.621E+00, 2.600E+01, 4.500E+00, 0.000E+00/
      DATA (PH1(I,30,25, 3),I=1,6) /1.162E+03, 2.821E+02,
     1 1.215E+02, 7.000E+01, 3.700E+00, 1.000E-01/
      DATA (PH1(I,30,25, 4),I=1,6) /2.369E+02, 6.326E+01,
     1 5.006E+00, 1.800E+02, 4.200E+00, 0.000E+00/
      DATA (PH1(I,30,25, 5),I=1,6) /1.896E+02, 1.177E+02,
     1 2.996E+01, 4.500E+01, 5.900E+00, 2.600E-01/
      DATA (PH1(I,30,25, 6),I=1,6) /1.080E+02, 5.878E+01,
     1 4.151E+02, 4.000E+01, 7.000E+00, 3.000E-01/
      DATA (PH1(I,30,26, 1),I=1,6) /9.751E+03, 1.732E+03,
     1 5.666E+00, 8.751E+01, 1.484E+00, 0.000E+00/
      DATA (PH1(I,30,26, 2),I=1,6) /1.259E+03, 2.231E+02,
     1 6.480E+00, 2.364E+01, 4.891E+00, 0.000E+00/
      DATA (PH1(I,30,26, 3),I=1,6) /1.129E+03, 2.523E+02,
     1 1.473E+02, 6.205E+01, 3.933E+00, 6.078E-09/
      DATA (PH1(I,30,26, 4),I=1,6) /2.075E+02, 6.183E+01,
     1 4.578E+00, 1.352E+02, 4.422E+00, 0.000E+00/
      DATA (PH1(I,30,26, 5),I=1,6) /1.612E+02, 1.030E+02,
     1 3.598E+01, 2.554E+01, 6.932E+00, 2.426E-01/
      DATA (PH1(I,30,26, 6),I=1,6) /8.260E+01, 6.378E+01,
     1 4.172E+02, 3.462E+01, 7.005E+00, 3.108E-01/
      DATA (PH1(I,30,27, 1),I=1,6) /9.733E+03, 1.890E+03,
     1 4.708E+00, 1.202E+02, 1.363E+00, 0.000E+00/
      DATA (PH1(I,30,27, 2),I=1,6) /1.239E+03, 2.366E+02,
     1 6.210E+00, 2.347E+01, 4.803E+00, 0.000E+00/
      DATA (PH1(I,30,27, 3),I=1,6) /1.101E+03, 2.528E+02,
     1 1.467E+02, 6.144E+01, 3.938E+00, 6.071E-09/
      DATA (PH1(I,30,27, 4),I=1,6) /1.826E+02, 6.427E+01,
     1 4.285E+00, 1.295E+02, 4.424E+00, 0.000E+00/
      DATA (PH1(I,30,27, 5),I=1,6) /1.365E+02, 1.033E+02,
     1 3.582E+01, 2.337E+01, 7.076E+00, 2.470E-01/
      DATA (PH1(I,30,27, 6),I=1,6) /5.940E+01, 5.376E+01,
     1 6.007E+02, 2.429E+01, 7.985E+00, 3.026E-01/
      DATA (PH1(I,30,28, 1),I=1,6) /9.713E+03, 1.981E+03,
     1 4.228E+00, 1.057E+02, 1.361E+00, 0.000E+00/
      DATA (PH1(I,30,28, 2),I=1,6) /1.222E+03, 2.686E+02,
     1 5.600E+00, 2.394E+01, 4.574E+00, 0.000E+00/
      DATA (PH1(I,30,28, 3),I=1,6) /1.077E+03, 2.519E+02,
     1 1.478E+02, 6.090E+01, 3.949E+00, 6.087E-09/
      DATA (PH1(I,30,28, 4),I=1,6) /1.601E+02, 6.684E+01,
     1 4.057E+00, 1.255E+02, 4.414E+00, 0.000E+00/
      DATA (PH1(I,30,28, 5),I=1,6) /1.143E+02, 1.045E+02,
     1 3.568E+01, 2.181E+01, 7.163E+00, 2.533E-01/
      DATA (PH1(I,30,28, 6),I=1,6) /3.972E+01, 3.640E+01,
     1 1.498E+03, 1.188E+01, 1.082E+01, 3.519E-01/
      DATA (PH1(I,30,29, 1),I=1,6) /9.691E+03, 1.125E+03,
     1 1.382E+01, 2.765E+01, 2.242E+00, 0.000E+00/
      DATA (PH1(I,30,29, 2),I=1,6) /1.209E+03, 2.121E+02,
     1 6.713E+00, 2.342E+01, 4.988E+00, 0.000E+00/
      DATA (PH1(I,30,29, 3),I=1,6) /1.065E+03, 2.664E+02,
     1 1.329E+02, 6.260E+01, 3.864E+00, 5.913E-09/
      DATA (PH1(I,30,29, 4),I=1,6) /1.479E+02, 6.443E+01,
     1 4.075E+00, 1.175E+02, 4.506E+00, 0.000E+00/
      DATA (PH1(I,30,29, 5),I=1,6) /1.021E+02, 1.047E+02,
     1 3.618E+01, 2.122E+01, 7.201E+00, 2.534E-01/
      DATA (PH1(I,30,29, 6),I=1,6) /2.694E+01, 2.593E+01,
     1 3.408E+03, 7.832E+00, 1.363E+01, 3.562E-01/
      DATA (PH1(I,30,29, 7),I=1,6) /1.796E+01, 1.411E+01,
     1 8.168E-01, 8.775E+02, 4.393E+00, 4.021E-04/
      DATA (PH1(I,30,30, 1),I=1,6) /9.667E+03, 8.320E+02,
     1 2.586E+01, 4.497E+01, 2.215E+00, 0.000E+00/
      DATA (PH1(I,30,30, 2),I=1,6) /1.203E+03, 9.755E+01,
     1 9.077E+00, 3.219E+01, 5.888E+00, 0.000E+00/
      DATA (PH1(I,30,30, 3),I=1,6) /1.037E+03, 3.486E+02,
     1 7.784E+01, 8.298E+01, 3.391E+00, 1.125E-08/
      DATA (PH1(I,30,30, 4),I=1,6) /1.450E+02, 6.195E+01,
     1 4.094E+00, 1.086E+02, 4.614E+00, 0.000E+00/
      DATA (PH1(I,30,30, 5),I=1,6) /9.700E+01, 1.025E+02,
     1 3.930E+01, 1.876E+01, 7.456E+00, 2.559E-01/
      DATA (PH1(I,30,30, 6),I=1,6) /1.730E+01, 1.818E+01,
     1 1.017E+04, 5.288E+00, 1.736E+01, 4.667E-01/
      DATA (PH1(I,30,30, 7),I=1,6) /9.394E+00, 1.673E+01,
     1 1.236E+00, 1.029E+03, 4.259E+00, 3.962E-02/
      DATA (PH2(I, 1, 1),I=1,7) /4.298E-01, 5.475E+04,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 2, 1),I=1,7) /1.720E+00, 1.369E+04,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 2, 2),I=1,7) /1.361E+01, 9.492E+02,
     1 1.469E+00, 3.188E+00, 2.039E+00, 4.434E-01, 2.136E+00/
      DATA (PH2(I, 3, 1),I=1,7) /3.871E+00, 6.083E+03,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 3, 2),I=1,7) /2.006E+01, 3.201E+02,
     1 7.391E+00, 2.916E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 3, 3),I=1,7) /3.107E+00, 6.245E+01,
     1 1.501E+01, 4.895E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 4, 1),I=1,7) /6.879E+00, 3.422E+03,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 4, 2),I=1,7) /1.760E+01, 5.458E+02,
     1 1.719E+01, 3.157E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 4, 3),I=1,7) /1.181E+00, 2.678E+02,
     1 5.645E+00, 1.170E+01, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 4, 4),I=1,7) /9.539E+00, 2.932E+05,
     1 4.301E-01, 1.052E+01, 3.655E-01, 8.278E-04, 1.269E-02/
      DATA (PH2(I, 5, 1),I=1,7) /1.075E+01, 2.190E+03,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 5, 2),I=1,7) /3.336E+01, 2.846E+02,
     1 2.163E+01, 2.624E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 5, 3),I=1,7) /1.041E+00, 5.393E+01,
     1 1.767E+01, 9.540E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 5, 4),I=1,7) /2.869E+00, 1.859E+04,
     1 1.783E+00, 1.618E+01, 3.503E+00, 4.960E-03, 3.400E-02/
      DATA (PH2(I, 5, 5),I=1,7) /5.213E-01, 5.466E+00,
     1 8.618E+00, 1.728E+01, 1.887E+01, 1.319E+01, 4.556E+00/
      DATA (PH2(I, 6, 1),I=1,7) /1.548E+01, 1.521E+03,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 6, 2),I=1,7) /4.624E+01, 2.344E+02,
     1 2.183E+01, 2.581E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 6, 3),I=1,7) /3.506E+00, 1.068E+02,
     1 1.436E+01, 7.457E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 6, 4),I=1,7) /4.614E+00, 1.539E+04,
     1 1.737E+00, 1.593E+01, 5.922E+00, 4.378E-03, 2.528E-02/
      DATA (PH2(I, 6, 5),I=1,7) /4.058E-01, 8.709E+00,
     1 1.261E+02, 8.578E+00, 2.093E+00, 4.929E+01, 3.234E+00/
      DATA (PH2(I, 6, 6),I=1,7) /2.144E+00, 5.027E+02,
     1 6.216E+01, 5.101E+00, 9.157E-02, 1.133E+00, 1.607E+00/
      DATA (PH2(I, 7, 1),I=1,7) /2.108E+01, 1.117E+03,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 7, 2),I=1,7) /6.943E+01, 1.519E+02,
     1 2.627E+01, 2.315E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 7, 3),I=1,7) /4.471E+00, 8.376E+01,
     1 3.297E+01, 6.003E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 7, 4),I=1,7) /5.494E+00, 1.690E+04,
     1 1.714E+00, 1.706E+01, 7.904E+00, 6.415E-03, 1.937E-02/
      DATA (PH2(I, 7, 5),I=1,7) /2.420E-01, 9.375E-01,
     1 2.788E+02, 9.156E+00, 1.850E+00, 1.877E+02, 3.999E+00/
      DATA (PH2(I, 7, 6),I=1,7) /6.128E-02, 1.944E+00,
     1 8.163E+02, 8.773E+00, 1.043E+01, 4.280E+02, 2.030E+01/
      DATA (PH2(I, 7, 7),I=1,7) /4.034E+00, 8.235E+02,
     1 8.033E+01, 3.928E+00, 9.097E-02, 8.598E-01, 2.325E+00/
      DATA (PH2(I, 8, 1),I=1,7) /2.754E+01, 8.554E+02,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 8, 2),I=1,7) /8.709E+01, 1.329E+02,
     1 2.535E+01, 2.336E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 8, 3),I=1,7) /7.824E+00, 6.864E+01,
     1 3.210E+01, 5.495E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 8, 4),I=1,7) /2.854E+00, 1.642E+04,
     1 1.792E+00, 2.647E+01, 2.836E+01, 3.036E-02, 5.554E-02/
      DATA (PH2(I, 8, 5),I=1,7) /2.044E-01, 8.659E-01,
     1 4.931E+02, 8.785E+00, 3.143E+00, 3.328E+02, 4.285E+01/
      DATA (PH2(I, 8, 6),I=1,7) /1.723E-01, 6.753E+02,
     1 3.852E+02, 6.822E+00, 1.191E-01, 3.839E-03, 4.569E-01/
      DATA (PH2(I, 8, 7),I=1,7) /1.386E+00, 5.967E+01,
     1 3.175E+01, 8.943E+00, 1.934E-02, 2.131E+01, 1.503E-02/
      DATA (PH2(I, 8, 8),I=1,7) /1.240E+00, 1.745E+03,
     1 3.784E+00, 1.764E+01, 7.589E-02, 8.698E+00, 1.271E-01/
      DATA (PH2(I, 9, 1),I=1,7) /3.485E+01, 6.759E+02,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 9, 2),I=1,7) /1.131E+02, 1.039E+02,
     1 2.657E+01, 2.255E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 9, 3),I=1,7) /2.563E+00, 6.930E+01,
     1 7.547E+01, 6.448E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I, 9, 4),I=1,7) /4.008E+00, 1.157E+04,
     1 1.848E+00, 2.446E+01, 2.411E+01, 2.071E-02, 3.998E-02/
      DATA (PH2(I, 9, 5),I=1,7) /7.286E-01, 4.690E-01,
     1 1.400E+02, 9.718E+00, 2.570E-01, 1.506E+02, 2.574E-01/
      DATA (PH2(I, 9, 6),I=1,7) /7.744E-01, 3.165E+00,
     1 1.099E+02, 9.203E+00, 6.812E+00, 9.531E+01, 9.781E+00/
      DATA (PH2(I, 9, 7),I=1,7) /2.542E+00, 1.541E+02,
     1 5.742E+01, 6.614E+00, 1.115E+00, 1.641E+01, 5.124E+00/
      DATA (PH2(I, 9, 8),I=1,7) /1.763E+00, 8.013E+01,
     1 1.667E+01, 1.050E+01, 5.103E-01, 1.715E+01, 7.724E-01/
      DATA (PH2(I, 9, 9),I=1,7) /1.297E+01, 3.803E+03,
     1 2.587E+00, 7.275E+00, 2.170E-03, 1.701E-04, 1.345E-02/
      DATA (PH2(I,10, 1),I=1,7) /4.304E+01, 5.475E+02,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,10, 2),I=1,7) /1.586E+02, 6.695E+01,
     1 3.352E+01, 2.002E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,10, 3),I=1,7) /1.003E+01, 5.631E+01,
     1 3.628E+01, 5.585E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,10, 4),I=1,7) /4.888E+00, 1.198E+04,
     1 1.788E+00, 2.550E+01, 2.811E+01, 2.536E-02, 4.417E-02/
      DATA (PH2(I,10, 5),I=1,7) /1.499E+00, 9.854E-01,
     1 1.350E+02, 8.836E+00, 1.656E+00, 1.042E+02, 1.435E+00/
      DATA (PH2(I,10, 6),I=1,7) /1.248E+00, 2.430E+00,
     1 1.066E+02, 8.999E+00, 6.855E-01, 9.169E+01, 3.702E-01/
      DATA (PH2(I,10, 7),I=1,7) /5.566E+00, 1.685E+03,
     1 6.409E+02, 3.056E+00, 8.290E-03, 5.149E+00, 6.687E+00/
      DATA (PH2(I,10, 8),I=1,7) /7.753E-01, 5.708E+00,
     1 6.725E+01, 1.005E+01, 4.633E-01, 7.654E+01, 2.023E+00/
      DATA (PH2(I,10, 9),I=1,7) /1.247E+01, 1.583E+03,
     1 3.935E+00, 7.810E+00, 6.558E-02, 1.520E+00, 1.084E-01/
      DATA (PH2(I,10,10),I=1,7) /4.870E+00, 4.287E+03,
     1 5.798E+00, 8.355E+00, 2.434E-01, 4.236E-02, 5.873E+00/
      DATA (PH2(I,11, 1),I=1,7) /5.211E+01, 4.525E+02,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,11, 2),I=1,7) /2.268E+02, 3.995E+01,
     1 5.315E+01, 1.678E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,11, 3),I=1,7) /1.391E+01, 4.729E+01,
     1 3.889E+01, 5.265E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,11, 4),I=1,7) /1.535E+02, 7.215E+03,
     1 3.886E-01, 8.476E+00, 9.121E-01, 1.667E-01, 1.766E-02/
      DATA (PH2(I,11, 5),I=1,7) /2.096E+00, 1.609E+00,
     1 2.473E+02, 7.681E+00, 1.895E+00, 9.940E+01, 3.278E+00/
      DATA (PH2(I,11, 6),I=1,7) /4.846E+01, 7.101E+01,
     1 3.945E+01, 2.832E+00, 1.285E-02, 9.603E-04, 6.378E-03/
      DATA (PH2(I,11, 7),I=1,7) /5.408E+00, 2.346E+01,
     1 2.913E+01, 8.260E+00, 9.275E-01, 2.204E+01, 7.577E-01/
      DATA (PH2(I,11, 8),I=1,7) /6.690E-01, 2.330E+00,
     1 1.205E+02, 9.714E+00, 7.365E-01, 1.383E+02, 4.260E+00/
      DATA (PH2(I,11, 9),I=1,7) /1.069E+01, 1.885E+03,
     1 3.613E+00, 9.803E+00, 8.579E-02, 3.725E+00, 2.279E-01/
      DATA (PH2(I,11,10),I=1,7) /8.203E+00, 1.040E+03,
     1 8.259E+00, 7.362E+00, 2.328E+00, 3.375E+00, 4.010E+00/
      DATA (PH2(I,11,11),I=1,7) /6.139E+00, 1.601E+00,
     1 6.148E+03, 3.839E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,12, 1),I=1,7) /6.203E+01, 3.802E+02,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,12, 2),I=1,7) /2.042E+02, 6.140E+01,
     1 2.778E+01, 2.161E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,12, 3),I=1,7) /1.452E+01, 4.427E+01,
     1 3.826E+01, 5.460E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,12, 4),I=1,7) /3.482E+01, 9.008E+02,
     1 1.823E+00, 1.444E+01, 2.751E+00, 5.444E+00, 7.918E-02/
      DATA (PH2(I,12, 5),I=1,7) /4.884E-01, 6.344E-02,
     1 5.085E+02, 9.385E+00, 6.666E-01, 5.348E+02, 3.997E-03/
      DATA (PH2(I,12, 6),I=1,7) /3.570E+00, 3.104E+00,
     1 6.060E+01, 8.857E+00, 1.422E+00, 5.452E+01, 2.078E+00/
      DATA (PH2(I,12, 7),I=1,7) /1.711E+00, 2.185E+00,
     1 9.350E+01, 9.202E+00, 6.325E-01, 1.007E+02, 1.729E+00/
      DATA (PH2(I,12, 8),I=1,7) /9.762E-01, 1.728E+00,
     1 9.184E+01, 1.006E+01, 8.090E-01, 1.276E+02, 3.979E+00/
      DATA (PH2(I,12, 9),I=1,7) /2.912E+01, 1.394E+03,
     1 2.895E+00, 6.487E+00, 4.326E-02, 9.402E-01, 1.135E-01/
      DATA (PH2(I,12,10),I=1,7) /1.086E+01, 5.377E+02,
     1 9.779E+00, 7.117E+00, 2.604E+00, 4.860E+00, 3.722E+00/
      DATA (PH2(I,12,11),I=1,7) /8.139E+00, 3.278E+00,
     1 4.341E+07, 3.610E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,12,12),I=1,7) /1.197E+01, 1.372E+08,
     1 2.228E-01, 1.574E+01, 2.805E-01, 0.000E+00, 0.000E+00/
      DATA (PH2(I,13, 1),I=1,7) /7.281E+01, 3.239E+02,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,13, 2),I=1,7) /2.738E+02, 4.036E+01,
     1 3.567E+01, 1.915E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,13, 3),I=1,7) /2.355E+01, 3.388E+01,
     1 3.432E+01, 5.085E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,13, 4),I=1,7) /8.044E+00, 1.774E+04,
     1 1.653E+00, 2.655E+01, 2.953E+01, 2.538E-02, 1.203E-02/
      DATA (PH2(I,13, 5),I=1,7) /1.842E+00, 4.982E-01,
     1 2.568E+02, 8.406E+00, 6.945E-01, 1.719E+02, 6.595E+00/
      DATA (PH2(I,13, 6),I=1,7) /4.866E-01, 2.350E-01,
     1 7.216E+02, 8.659E+00, 2.773E-01, 5.704E+02, 1.580E-01/
      DATA (PH2(I,13, 7),I=1,7) /2.636E+00, 1.889E+02,
     1 1.338E+02, 6.204E+00, 1.836E+00, 3.552E+01, 8.223E-03/
      DATA (PH2(I,13, 8),I=1,7) /3.483E-01, 1.962E-02,
     1 1.856E+01, 2.084E+01, 8.839E+00, 5.675E-02, 2.768E-01/
      DATA (PH2(I,13, 9),I=1,7) /2.414E+01, 2.925E+02,
     1 6.973E+00, 6.724E+00, 1.000E-01, 3.495E+00, 2.701E-01/
      DATA (PH2(I,13,10),I=1,7) /3.130E+00, 1.513E+01,
     1 1.674E+01, 1.180E+01, 5.342E+00, 3.994E+01, 4.803E+00/
      DATA (PH2(I,13,11),I=1,7) /1.027E+01, 4.915E+00,
     1 1.990E+06, 3.477E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,13,12),I=1,7) /2.048E-01, 6.948E-02,
     1 5.675E+02, 9.049E+00, 4.615E-01, 9.149E+01, 6.565E-01/
      DATA (PH2(I,13,13),I=1,7) /1.381E+01, 7.195E+00,
     1 1.621E+03, 3.642E+00, 3.166E-01, 2.041E-01, 4.753E-01/
      DATA (PH2(I,14, 1),I=1,7) /8.447E+01, 2.793E+02,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,14, 2),I=1,7) /2.752E+02, 4.754E+01,
     1 2.848E+01, 2.135E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,14, 3),I=1,7) /3.560E+01, 2.539E+01,
     1 3.307E+01, 4.728E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,14, 4),I=1,7) /1.205E+01, 1.992E+04,
     1 1.582E+00, 2.425E+01, 2.392E+01, 1.990E-02, 1.007E-02/
      DATA (PH2(I,14, 5),I=1,7) /8.787E-01, 1.950E-01,
     1 7.461E+02, 8.302E+00, 4.489E-01, 4.528E+02, 1.015E+00/
      DATA (PH2(I,14, 6),I=1,7) /3.343E-01, 1.465E-01,
     1 1.404E+03, 8.503E+00, 1.646E+00, 1.036E+03, 2.936E-01/
      DATA (PH2(I,14, 7),I=1,7) /7.655E-01, 3.477E-01,
     1 3.733E+02, 8.986E+00, 1.476E-03, 3.850E+02, 8.999E-02/
      DATA (PH2(I,14, 8),I=1,7) /3.277E-01, 6.680E-02,
     1 4.132E+01, 1.606E+01, 3.280E+00, 1.149E-02, 6.396E-01/
      DATA (PH2(I,14, 9),I=1,7) /6.305E+01, 7.293E+01,
     1 1.558E+02, 2.400E+00, 2.989E-03, 1.115E+00, 8.051E-02/
      DATA (PH2(I,14,10),I=1,7) /7.761E-01, 8.863E-01,
     1 1.541E+02, 9.980E+00, 1.303E+00, 2.009E+02, 4.537E+00/
      DATA (PH2(I,14,11),I=1,7) /1.288E+01, 6.083E+00,
     1 1.356E+06, 3.353E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,14,12),I=1,7) /1.659E-01, 5.790E-04,
     1 1.474E+02, 1.336E+01, 8.626E-01, 9.613E+01, 6.442E-01/
      DATA (PH2(I,14,13),I=1,7) /2.556E+00, 4.140E+00,
     1 1.337E+01, 1.191E+01, 1.570E+00, 6.634E+00, 1.272E-01/
      DATA (PH2(I,14,14),I=1,7) /2.317E+01, 2.506E+01,
     1 2.057E+01, 3.546E+00, 2.837E-01, 1.672E-05, 4.207E-01/
      DATA (PH2(I,16, 1),I=1,7) /1.104E+02, 2.139E+02,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,16, 2),I=1,7) /4.390E+02, 2.453E+01,
     1 4.405E+01, 1.765E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,16, 3),I=1,7) /3.310E+01, 2.555E+01,
     1 3.821E+01, 5.037E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,16, 4),I=1,7) /1.474E+01, 2.294E+04,
     1 1.529E+00, 2.568E+01, 2.738E+01, 2.203E-02, 1.073E-02/
      DATA (PH2(I,16, 5),I=1,7) /2.443E+00, 3.490E-01,
     1 5.411E+02, 7.769E+00, 7.033E-01, 2.279E+02, 1.172E+00/
      DATA (PH2(I,16, 6),I=1,7) /6.485E+00, 1.275E+01,
     1 6.583E+01, 7.692E+00, 1.678E+00, 3.426E+01, 1.370E-01/
      DATA (PH2(I,16, 7),I=1,7) /1.040E+01, 5.364E+01,
     1 3.641E+01, 7.090E+00, 2.310E+00, 1.775E+01, 1.663E+00/
      DATA (PH2(I,16, 8),I=1,7) /1.526E-01, 9.646E+03,
     1 1.438E+03, 5.977E+00, 1.492E+00, 1.615E-03, 4.049E-01/
      DATA (PH2(I,16, 9),I=1,7) /1.462E+01, 3.161E+01,
     1 1.611E+01, 8.642E+00, 1.153E-03, 1.869E+01, 3.037E-01/
      DATA (PH2(I,16,10),I=1,7) /3.757E-01, 5.703E-01,
     1 1.460E+02, 1.135E+01, 1.503E+00, 2.222E+02, 4.606E+00/
      DATA (PH2(I,16,11),I=1,7) /1.413E+01, 9.139E+00,
     1 1.656E+03, 3.626E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,16,12),I=1,7) /1.713E-01, 5.072E-04,
     1 1.986E+02, 1.307E+01, 7.880E-01, 9.424E+01, 6.265E-01/
      DATA (PH2(I,16,13),I=1,7) /2.173E+00, 2.606E+00,
     1 6.641E+01, 8.655E+00, 1.863E+00, 1.975E+01, 3.361E+00/
      DATA (PH2(I,16,14),I=1,7) /2.027E+00, 6.666E+00,
     1 5.454E+01, 8.611E+00, 4.109E+00, 1.568E+01, 9.421E+00/
      DATA (PH2(I,16,15),I=1,7) /8.787E+00, 3.136E+02,
     1 3.442E+00, 1.281E+01, 7.354E-01, 2.782E+00, 1.788E-01/
      DATA (PH2(I,16,16),I=1,7) /1.808E+01, 4.564E+04,
     1 1.000E+00, 1.361E+01, 6.385E-01, 9.935E-01, 2.486E-01/
      DATA (PH2(I,18, 1),I=1,7) /1.399E+02, 1.690E+02,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,18, 2),I=1,7) /4.468E+02, 3.108E+01,
     1 3.039E+01, 2.092E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,18, 3),I=1,7) /4.154E+01, 2.135E+01,
     1 4.118E+01, 4.945E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,18, 4),I=1,7) /1.888E+01, 2.571E+04,
     1 1.475E+00, 2.634E+01, 2.909E+01, 2.445E-02, 1.054E-02/
      DATA (PH2(I,18, 5),I=1,7) /1.557E+00, 4.997E-02,
     1 5.031E+02, 8.966E+00, 2.938E-01, 4.552E+02, 6.459E+00/
      DATA (PH2(I,18, 6),I=1,7) /3.209E-01, 2.459E-02,
     1 2.285E+03, 8.810E+00, 6.692E-01, 2.068E+03, 2.113E+01/
      DATA (PH2(I,18, 7),I=1,7) /5.310E+00, 7.018E-01,
     1 1.001E+02, 8.939E+00, 4.987E-01, 1.099E+02, 2.202E-01/
      DATA (PH2(I,18, 8),I=1,7) /1.257E-01, 1.760E+03,
     1 1.579E+03, 6.714E+00, 1.975E+00, 3.286E-03, 3.226E-01/
      DATA (PH2(I,18, 9),I=1,7) /1.040E+01, 8.204E+00,
     1 1.495E+01, 1.115E+01, 9.203E-04, 3.804E+01, 6.390E-01/
      DATA (PH2(I,18,10),I=1,7) /1.926E-01, 8.279E-01,
     1 2.392E+02, 1.121E+01, 1.434E+00, 3.814E+01, 4.649E+00/
      DATA (PH2(I,18,11),I=1,7) /3.884E+00, 3.295E+01,
     1 7.082E+02, 4.645E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,18,12),I=1,7) /2.966E-02, 3.693E+00,
     1 9.951E+03, 7.313E+00, 1.363E-02, 4.383E-04, 2.513E+00/
      DATA (PH2(I,18,13),I=1,7) /5.440E-01, 1.080E+00,
     1 9.419E+02, 7.582E+00, 1.107E+01, 1.700E+02, 1.587E+01/
      DATA (PH2(I,18,14),I=1,7) /1.031E+01, 9.946E+00,
     1 7.444E+01, 6.261E+00, 4.885E-01, 6.406E+00, 3.659E-03/
      DATA (PH2(I,18,15),I=1,7) /6.953E+00, 2.035E+01,
     1 1.400E+01, 9.595E+00, 8.842E-01, 7.501E+00, 1.806E-01/
      DATA (PH2(I,18,16),I=1,7) /1.417E+01, 3.580E+01,
     1 3.776E+01, 5.742E+00, 6.316E-01, 2.384E+00, 1.794E+00/
      DATA (PH2(I,18,17),I=1,7) /2.494E+01, 2.503E+01,
     1 1.272E+02, 4.288E+00, 5.108E-01, 9.299E-01, 7.195E-01/
      DATA (PH2(I,18,18),I=1,7) /1.709E+01, 2.106E+01,
     1 2.645E+02, 4.796E+00, 4.185E-01, 1.688E+00, 8.943E-01/
      DATA (PH2(I,20, 1),I=1,7) /1.729E+02, 1.369E+02,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,20, 2),I=1,7) /6.297E+02, 1.936E+01,
     1 3.921E+01, 1.862E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,20, 3),I=1,7) /9.472E+01, 1.105E+01,
     1 3.818E+01, 4.192E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,20, 4),I=1,7) /2.618E+01, 2.028E+04,
     1 1.456E+00, 2.560E+01, 2.803E+01, 2.402E-02, 9.323E-03/
      DATA (PH2(I,20, 5),I=1,7) /4.293E+00, 1.293E+00,
     1 1.691E+01, 1.438E+01, 3.461E-05, 9.363E-01, 4.589E-02/
      DATA (PH2(I,20, 6),I=1,7) /1.309E+02, 5.513E+01,
     1 3.828E+02, 2.023E+00, 9.084E-02, 1.833E-02, 9.359E-01/
      DATA (PH2(I,20, 7),I=1,7) /9.980E+00, 1.116E+00,
     1 5.918E+01, 9.005E+00, 3.879E+00, 7.104E+01, 5.311E+00/
      DATA (PH2(I,20, 8),I=1,7) /1.008E+01, 1.849E+03,
     1 1.792E+04, 2.868E+00, 2.410E+02, 6.138E-03, 6.931E+01/
      DATA (PH2(I,20, 9),I=1,7) /2.345E+01, 1.227E+01,
     1 1.312E+01, 9.771E+00, 6.842E-04, 2.417E+01, 5.469E-01/
      DATA (PH2(I,20,10),I=1,7) /2.288E-01, 9.384E-01,
     1 2.549E+02, 1.103E+01, 1.390E+00, 2.478E+01, 3.100E+00/
      DATA (PH2(I,20,11),I=1,7) /1.605E+01, 1.437E+01,
     1 6.989E+02, 3.857E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,20,12),I=1,7) /5.520E-02, 2.076E+02,
     1 1.790E+04, 5.893E+00, 1.843E-03, 2.826E-04, 1.657E+00/
      DATA (PH2(I,20,13),I=1,7) /1.366E+00, 6.641E-01,
     1 3.188E+02, 8.138E+00, 2.806E-01, 1.039E+02, 3.329E+00/
      DATA (PH2(I,20,14),I=1,7) /8.080E-01, 4.760E-01,
     1 3.682E+02, 8.634E+00, 5.720E-01, 1.487E+02, 1.283E+00/
      DATA (PH2(I,20,15),I=1,7) /9.515E+00, 7.642E+01,
     1 8.973E+01, 5.141E+00, 2.471E+00, 4.829E+00, 5.824E+00/
      DATA (PH2(I,20,16),I=1,7) /6.882E-01, 1.523E-01,
     1 1.502E+02, 1.061E+01, 8.227E+00, 1.210E+02, 3.876E+00/
      DATA (PH2(I,20,17),I=1,7) /4.255E+00, 7.736E+00,
     1 1.355E+01, 1.236E+01, 1.369E+00, 1.467E+01, 3.298E-02/
      DATA (PH2(I,20,18),I=1,7) /2.436E+01, 3.815E+01,
     1 2.931E+02, 3.944E+00, 3.126E-01, 1.802E+00, 1.233E+00/
      DATA (PH2(I,20,19),I=1,7) /1.553E+01, 1.064E+07,
     1 7.790E-01, 2.130E+01, 6.453E-01, 2.161E-03, 6.706E-02/
      DATA (PH2(I,20,20),I=1,7) /1.278E+01, 5.370E+05,
     1 3.162E-01, 1.242E+01, 4.477E-01, 1.012E-03, 1.851E-02/
      DATA (PH2(I,26, 1),I=1,7) /2.932E+02, 8.099E+01,
     1 3.288E+01, 2.963E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,26, 2),I=1,7) /1.057E+03, 1.195E+01,
     1 5.769E+01, 1.718E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,26, 3),I=1,7) /7.326E+01, 1.276E+01,
     1 4.914E+01, 4.941E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,26, 4),I=1,7) /4.575E+01, 2.580E+04,
     1 1.358E+00, 2.604E+01, 2.723E+01, 3.582E-02, 8.712E-03/
      DATA (PH2(I,26, 5),I=1,7) /9.713E+00, 7.204E-02,
     1 1.853E+02, 8.843E+00, 9.551E-03, 1.702E+02, 4.263E+00/
      DATA (PH2(I,26, 6),I=1,7) /9.243E+00, 1.098E+01,
     1 7.637E+01, 7.962E+00, 1.748E+00, 4.446E+01, 3.512E+00/
      DATA (PH2(I,26, 7),I=1,7) /2.011E+01, 4.455E-01,
     1 4.236E+01, 9.724E+00, 2.757E+00, 6.847E+01, 3.989E+00/
      DATA (PH2(I,26, 8),I=1,7) /7.519E-04, 6.066E-05,
     1 1.606E+06, 8.813E+00, 4.398E+00, 1.915E+06, 3.140E+01/
      DATA (PH2(I,26, 9),I=1,7) /3.190E+01, 2.388E+00,
     1 2.186E+01, 9.589E+00, 2.902E-02, 3.805E+01, 4.805E-01/
      DATA (PH2(I,26,10),I=1,7) /3.444E-01, 1.452E+00,
     1 3.960E+02, 1.013E+01, 1.264E+00, 2.891E+01, 3.404E+00/
      DATA (PH2(I,26,11),I=1,7) /2.873E+01, 1.207E+01,
     1 5.150E+02, 3.846E+00, 0.000E+00, 0.000E+00, 0.000E+00/
      DATA (PH2(I,26,12),I=1,7) /5.555E-02, 2.108E+02,
     1 2.045E+04, 6.033E+00, 1.885E-03, 2.706E-04, 1.628E+00/
      DATA (PH2(I,26,13),I=1,7) /8.509E-01, 1.454E-01,
     1 1.239E+03, 8.066E+00, 4.937E-01, 4.505E+02, 2.504E+00/
      DATA (PH2(I,26,14),I=1,7) /1.317E-01, 2.791E-03,
     1 2.487E+03, 9.791E+00, 6.938E-01, 2.170E+03, 6.852E-03/
      DATA (PH2(I,26,15),I=1,7) /6.295E+00, 1.738E+00,
     1 1.130E+02, 8.037E+00, 3.096E-01, 4.671E+01, 1.425E-01/
      DATA (PH2(I,26,16),I=1,7) /8.284E+00, 3.281E+00,
     1 5.360E+01, 8.571E+00, 3.279E-01, 2.971E+01, 5.220E-01/
      DATA (PH2(I,26,17),I=1,7) /6.886E+01, 6.470E+01,
     1 2.062E+01, 4.111E+00, 2.778E-04, 1.190E-05, 6.570E-03/
      DATA (PH2(I,26,18),I=1,7) /6.741E+00, 2.687E+01,
     1 1.807E+02, 6.290E+00, 2.387E-04, 2.494E+01, 8.251E+00/
      DATA (PH2(I,26,19),I=1,7) /7.098E-02, 1.979E+01,
     1 1.745E+04, 6.750E+00, 2.158E+02, 2.542E+03, 4.672E+02/
      DATA (PH2(I,26,20),I=1,7) /5.059E+00, 2.420E+04,
     1 4.850E+04, 2.374E+00, 2.516E-03, 4.546E-01, 2.683E+01/
      DATA (PH2(I,26,21),I=1,7) /2.656E+00, 5.259E-01,
     1 1.450E+01, 1.632E+01, 1.558E+01, 3.361E+01, 3.743E-03/
      DATA (PH2(I,26,22),I=1,7) /7.256E-01, 1.523E-03,
     1 3.736E+01, 1.767E+01, 5.064E+01, 8.871E+01, 5.280E-02/
      DATA (PH2(I,26,23),I=1,7) /2.544E+01, 3.653E+02,
     1 8.913E+00, 6.538E+00, 5.602E-01, 0.000E+00, 0.000E+00/
      DATA (PH2(I,26,24),I=1,7) /1.698E-01, 6.107E+00,
     1 1.555E+03, 8.055E+00, 8.698E+00, 1.760E+02, 1.847E+01/
      DATA (PH2(I,26,25),I=1,7) /1.761E-01, 4.365E+03,
     1 6.298E+03, 5.204E+00, 1.141E+01, 9.272E+01, 1.075E+02/
      DATA (PH2(I,26,26),I=1,7) /5.461E-02, 3.062E-01,
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
      END" > src/tmp.f
}

function print_meta {
    echo "id	${ID}
type	auto generated
comment	auto generated
model	${model}
tmin	${tmin/d0/}
tmax	${tmax/d0/}
src_abu	${src_abu}
n_pho	${n_pho}" > out/${output_root}_${ID}.txt
}

function run_script {
    h5fc -fopenmp -march=native -O3 src/tmp.f
    time ./a.out
    cat out/tmp*.esc > out/${output_root}_${ID}.esc
}

function run_model {
    model_path="/Users/silver/dat/sne/${model}.h5"
    buffer_size=${#model_path}
    output_root=$(echo "${model}" | tr '[:upper:]' '[:lower:]')
    
    ID="013"
    tmin=160.d0
    tmax=210.d0
    src_abu=2
    print_script
    print_meta
    run_script
    
    ID="014"
    tmin=400.d0
    tmax=450.d0
    print_script
    print_meta
    run_script
    
    ID="017"
    tmin=275.d0
    tmax=325.d0
    print_script
    print_meta
    run_script
    
    ID="015"
    tmin=1.d0
    tmax=1050.d0
    print_script
    print_meta
    run_script
    
    ID="016"
    src_abu=3
    print_script
    print_meta
    run_script
}

function run_stripped_model {
    model_path="/Users/silver/dat/sne/${model}.h5"
    buffer_size=${#model_path}
    output_root=$(echo "${model}" | tr '[:upper:]' '[:lower:]')
    
#    ID="013"
#    tmin=40.d0
#    tmax=60.d0
    src_abu=2
#    print_script
#    print_meta
#    run_script
#    
#    ID="014"
#    tmin=235.d0
#    tmax=265.d0
#    print_script
#    print_meta
#    run_script
#    
#    ID="017"
#    tmin=110.d0
#    tmax=140.d0
#    print_script
#    print_meta
#    run_script
#    
#    ID="018"
#    tmin=85.d0
#    tmax=115.d0
#    print_script
#    print_meta
#    run_script
#    
    ID="011"
    tmin=15.d0
    tmax=25.d0
    print_script
    print_meta
    run_script
#    
#    ID="015"
#    tmin=1.d0
#    tmax=1050.d0
#    print_script
#    print_meta
#    run_script
#    
#    ID="016"
#    src_abu=3
#    print_script
#    print_meta
#    run_script

    ID="010"
    src_abu=1
    tmin=15.d0
    tmax=25.d0
    print_script
    print_meta
    run_script
}

function run_ti_ni {
    model_path="/Users/silver/dat/sne/${model}.h5"
    buffer_size=${#model_path}
    output_root=$(echo "${model}" | tr '[:upper:]' '[:lower:]')
    tmin=1.d0
    tmax=1050.d0
    
    ID="020"
    src_abu=0
    print_script
    print_meta
    run_script
    
    ID="021"
    src_abu=1
    print_script
    print_meta
    run_script
}


################################################################
# Specific for the model
# Type Ib/c # model="Ib"
# Type Ib/c # t_mod=1565943.04402
# Type Ib/c # n_rad=1200
# Type Ib/c # n_tht=90
# Type Ib/c # n_phi=180
# Type Ib/c # n_pho=20000000
# Type Ib/c # run_stripped_model
# Type Ib/c # n_pho=5000000
# Type Ib/c # run_ti_ni
# Type Ib/c # 
# Type Ib/c # model="Ic"
# Type Ib/c # t_mod=1565943.04402
# Type Ib/c # n_rad=1200
# Type Ib/c # n_tht=90
# Type Ib/c # n_phi=180
# Type Ib/c # n_pho=20000000
# Type Ib/c # run_stripped_model
# Type Ib/c # n_pho=5000000
# Type Ib/c # run_ti_ni

# W7 model="W7"
# W7 t_mod=200000.
# W7 n_rad=176
# W7 n_tht=90
# W7 n_phi=180
# W7 n_pho=20000000
# W7 run_stripped_model
# W7 n_pho=5000000
# W7 #run_ti_ni

# model="M157b"
# t_mod=86381.6
# n_rad=1730
# n_tht=90
# n_phi=180
# n_pho=5000000
# run_model
# run_ti_ni
# 
# model="M158b"
# t_mod=86411.2
# n_rad=1789
# run_model
# run_ti_ni
# 
# model="M164a"
# t_mod=86413.6
# n_rad=1727
# run_model
# run_ti_ni
# 
# model="M167b"
# t_mod=86413.6
# n_rad=1723
# run_model
# run_ti_ni
# 
# model="M177a"
# t_mod=86392.6
# n_rad=1745
# run_model
# run_ti_ni
# 
# model="M178a"
# t_mod=86402.6
# n_rad=1775
# run_model
# run_ti_ni

# merger set 2 model="M157b+"
# merger set 2 t_mod=86381.086
# merger set 2 n_rad=1728
# merger set 2 n_tht=90
# merger set 2 n_phi=180
# merger set 2 n_pho=20000000
# merger set 2 run_model
# merger set 2 n_pho=5000000
# merger set 2 run_ti_ni
# merger set 2 
# merger set 2 model="M157b++"
# merger set 2 t_mod=86404.87
# merger set 2 n_rad=1726
# merger set 2 n_pho=20000000
# merger set 2 run_model
# merger set 2 n_pho=5000000
# merger set 2 run_ti_ni
# merger set 2 
# merger set 2 model="M158b-"
# merger set 2 t_mod=86418.42
# merger set 2 n_rad=1820
# merger set 2 n_pho=20000000
# merger set 2 run_model
# merger set 2 n_pho=5000000
# merger set 2 run_ti_ni
# merger set 2 
# merger set 2 model="M164a-"
# merger set 2 t_mod=86390.71
# merger set 2 n_rad=1763
# merger set 2 n_pho=20000000
# merger set 2 run_model
# merger set 2 n_pho=5000000
# merger set 2 run_ti_ni
# merger set 2 
# merger set 2 model="M167b+"
# merger set 2 t_mod=86387.4
# merger set 2 n_rad=1743
# merger set 2 n_pho=20000000
# merger set 2 run_model
# merger set 2 n_pho=5000000
# merger set 2 run_ti_ni
# merger set 2 
# merger set 2 model="M177a+"
# merger set 2 t_mod=86389.32
# merger set 2 n_rad=1743
# merger set 2 n_pho=20000000
# merger set 2 run_model
# merger set 2 n_pho=5000000
# merger set 2 run_ti_ni
# merger set 2 
# merger set 2 model="M178a+"
# merger set 2 t_mod=86417.27
# merger set 2 n_rad=1763
# merger set 2 n_pho=20000000
# merger set 2 run_model
# merger set 2 n_pho=5000000
# merger set 2 run_ti_ni

# merger set 3 model="M157b-1"
# merger set 3 t_mod=86389.9
# merger set 3 n_rad=1728
# merger set 3 n_tht=90
# merger set 3 n_phi=180
# merger set 3 n_pho=20000000
# merger set 3 run_model
# merger set 3 n_pho=5000000
# merger set 3 run_ti_ni
# merger set 3 
# merger set 3 model="M157b-2"
# merger set 3 t_mod=86401.5
# merger set 3 n_rad=1727
# merger set 3 n_pho=20000000
# merger set 3 run_model
# merger set 3 n_pho=5000000
# merger set 3 run_ti_ni
# merger set 3 
# merger set 3 model="M157b-3"
# merger set 3 t_mod=86395
# merger set 3 n_rad=1725
# merger set 3 n_pho=20000000
# merger set 3 run_model
# merger set 3 n_pho=5000000
# merger set 3 run_ti_ni
# merger set 3 
# merger set 3 model="M157b-4"
# merger set 3 t_mod=86397.3
# merger set 3 n_rad=1709
# merger set 3 n_pho=20000000
# merger set 3 run_model
# merger set 3 n_pho=5000000
# merger set 3 run_ti_ni
# merger set 3 
# merger set 3 model="M167b-1"
# merger set 3 t_mod=86406.3
# merger set 3 n_rad=1742
# merger set 3 n_pho=20000000
# merger set 3 run_model
# merger set 3 n_pho=5000000
# merger set 3 run_ti_ni
# merger set 3 
# merger set 3 model="M167b-2"
# merger set 3 t_mod=86400.3
# merger set 3 n_rad=1722
# merger set 3 n_pho=20000000
# merger set 3 run_model
# merger set 3 n_pho=5000000
# merger set 3 run_ti_ni

# # complement old runs model="B15_1Z"
# # complement old runs t_mod=13489751.4508
# n_rad=1200
# n_tht=90
# n_phi=180
# n_pho=5000000
# # complement old runs run_ti_ni
# # complement old runs 
# # complement old runs model="B15_01Z"
# # complement old runs run_ti_ni
# # complement old runs 
# # complement old runs model="B15"
# # complement old runs run_ti_ni
# # complement old runs 
# # complement old runs model="B15_025Z"
# # complement old runs run_ti_ni
# # complement old runs 
# # complement old runs model="HMM"
# # complement old runs t_mod=14832200.0
# # complement old runs run_ti_ni
# # complement old runs 
# # complement old runs model="N20"
# # complement old runs t_mod=12488928.8787
# # complement old runs run_ti_ni
# # complement old runs 
# # complement old runs model="N20_LMC"
# # complement old runs run_ti_ni
# # complement old runs 
# # complement old runs model="L15"
# # complement old runs t_mod=12630166.8468
# # complement old runs run_ti_ni
# # complement old runs 
# model="W15"
# t_mod=12755626.5895
# run_ti_ni
# 
# model="W15_1Z"
# run_ti_ni
# 
# model="IIb"
# t_mod=1565943.04402
# run_ti_ni
# n_pho=20000000
# run_stripped_model
# 
# model="Ib"
# run_stripped_model
# 
# model="Ic"
# run_stripped_model

#model="W7"
#t_mod=200000.
#n_rad=176
#n_tht=90
#n_phi=180
#n_pho=20000000
#run_stripped_model
#n_pho=5000000
##run_ti_ni

# 87A # model="B15_87A"
# 87A # t_mod=13489751.4508
# 87A # n_rad=1200
# 87A # n_tht=90
# 87A # n_phi=180
# 87A # n_pho=20000000
# 87A # run_model
# 87A # n_pho=5000000
# 87A # run_ti_ni

model="M157b-2_min"
t_mod=86401.5
n_rad=1727
n_tht=90
n_phi=180
n_pho=20000000
run_model
n_pho=5000000
run_ti_ni

model="M157b-2_max"
t_mod=86401.5
n_rad=1727
n_tht=90
n_phi=180
n_pho=20000000
run_model
n_pho=5000000
run_ti_ni

model="B15_LMC_min"
t_mod=13489751.4508
n_rad=1200
n_tht=90
n_phi=180
n_pho=20000000
run_model
n_pho=5000000
run_ti_ni

model="B15_LMC_max"
t_mod=13489751.4508
n_rad=1200
n_tht=90
n_phi=180
n_pho=20000000
run_model
n_pho=5000000
run_ti_ni

#model="M157b-2"
#t_mod=86401.5
#n_rad=1727
#run_model

# post processing time python kicks/kicks.py /Users/silver/dat/sne/ M157b
# post processing time python kicks/kicks.py /Users/silver/dat/sne/ M158b
# post processing time python kicks/kicks.py /Users/silver/dat/sne/ M164a
# post processing time python kicks/kicks.py /Users/silver/dat/sne/ M167b
# post processing time python kicks/kicks.py /Users/silver/dat/sne/ M177a
# post processing time python kicks/kicks.py /Users/silver/dat/sne/ M178a

# pp ./run_all.sh m157b
# pp ./run_all.sh m157b+
# pp ./run_all.sh m157b++
# pp ./run_all.sh m158b
# pp ./run_all.sh m158b-
# pp ./run_all.sh m164a
# pp ./run_all.sh m164a-
# pp ./run_all.sh m167b
# pp ./run_all.sh m167b+
# pp ./run_all.sh m177a
# pp ./run_all.sh m177a+
# pp ./run_all.sh m178a
# pp ./run_all.sh m178a+

#python kicks/kicks.py /Users/silver/dat/sne/ M157b+; python kicks/kicks.py /Users/silver/dat/sne/ M157b++; python kicks/kicks.py /Users/silver/dat/sne/ M158b-; python kicks/kicks.py /Users/silver/dat/sne/ M164a-; python kicks/kicks.py /Users/silver/dat/sne/ M167b+; python kicks/kicks.py /Users/silver/dat/sne/ M177a+; python kicks/kicks.py /Users/silver/dat/sne/ M178a+; ./run_all_stripped.sh iib_0.0; ./run_all_stripped.sh iib_1.0; ./run_all.sh M157b+; ./run_all.sh M157b++; ./run_all.sh M158b-; ./run_all.sh M164a-; ./run_all.sh M167b+; ./run_all.sh M177a+; ./run_all.sh M178a+; 

# pp ./run_all.sh b15_1z
# pp ./run_all.sh b15_lmc
# pp ./run_all.sh b15_01z
# pp ./run_all.sh b15_025z
# pp ./run_all.sh n20
# pp ./run_all.sh n20_lmc
# pp ./run_all.sh hmm
# pp ./run_all.sh l15
# pp ./run_all.sh w15
# pp ./run_all.sh w15_1z
# pp 
# pp ./run_all_stripped.sh iib
# pp ./run_all_stripped.sh w7
# pp ./run_all_stripped.sh ib
# pp ./run_all_stripped.sh ic
