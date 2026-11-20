      program belgi_frozen_strengths
C
C     Selection-complete line-strength driver for the frozen Xu/BELGI
C     methanol Hamiltonian.  Every same-species rank-1 pair is evaluated
C     with BELGI's native INTEN routine.  Xu et al. (2008) permanent
C     dipoles are used exactly; all cos(3 gamma) dipoles are fixed to zero.
C
      implicit none

      include 'common_dim.f'
      include 'common_cst.f'
      include 'common_dat.f'
      include 'common_par.f'
      include 'common_expdat.f'
      include 'common_big2.f'
      include 'common_cntl.f'
      include 'common_int.f'

      integer argc, output_jmax, output_kmax, output_vtmax
      integer ii, ios, candidates, retained
      integer count_positive, count_1em16, count_1em14
      integer count_1em12, count_1em10, count_1em8
      character*64 arg, vector_header, dipole_mode
      double precision smu2_min
      logical cross_vibrational_mode

      argc = command_argument_count()
      if (argc .ne. 4 .and. argc .ne. 5) then
         write(0,*)
     &    'usage: belgi_frozen_strengths JMAX KMAX VTMAX SMU2_MIN',
     &    ' [PERMANENT|VCO_A]'
         stop 2
      endif

      call get_command_argument(1, arg)
      read(arg,*,iostat=ios) output_jmax
      if (ios .ne. 0) stop 2
      call get_command_argument(2, arg)
      read(arg,*,iostat=ios) output_kmax
      if (ios .ne. 0) stop 2
      call get_command_argument(3, arg)
      read(arg,*,iostat=ios) output_vtmax
      if (ios .ne. 0) stop 2
      call get_command_argument(4, arg)
      read(arg,*,iostat=ios) smu2_min
      if (ios .ne. 0) stop 2
      dipole_mode = 'PERMANENT'
      if (argc .eq. 5) call get_command_argument(5,dipole_mode)
      cross_vibrational_mode = dipole_mode .eq. 'VCO_A'
      if (.not. cross_vibrational_mode .and.
     &    dipole_mode .ne. 'PERMANENT') stop 2

C     This fixed-size catalog driver is deliberately scoped to the complete
C     Xu/Mekhtiev quantum-number box requested for the project.
      if (output_jmax .lt. 1 .or. output_jmax .gt. 22) stop 2
      if (output_kmax .lt. 0 .or. output_kmax .gt. 14 .or.
     &    output_kmax .gt. output_jmax) stop 2
      if (output_vtmax .lt. 0 .or. output_vtmax .gt. 4) stop 2
      if (smu2_min .lt. 0.d0) stop 2

C     Python validates the named constants and supplies all 205 slots,
C     including explicit zeros, in their canonical data11.f order.
      read(5,'(A)',iostat=ios) vector_header
      if (ios .ne. 0 .or.
     &    vector_header .ne. 'BELGI_FROZEN_CST_205') stop 3
      do ii = 1, cstm
         read(5,*,iostat=ios) cst(1,ii)
         if (ios .ne. 0) stop 3
         cst(2,ii) = cst(1,ii)
      enddo

C     Xu et al. (2008), in Debye.  The three coefficients conventionally
C     denoted with a trailing 3 multiply cos(3 gamma)-type operators.
      if (cross_vibrational_mode) then
C        Owner-reviewed CO-stretch fundamental operator.  The 0.27-D
C        transition moment is placed along methanol's approximate CO/a
C        direction.  INTEN contracts it with the full mixed rotational
C        eigenfunctions; no nominal-Delta-K filter is applied.
         cst(1,102) = 0.27d0
         cst(1,103) = 0.d0
         cst(1,104) = 0.d0
         cst(1,105) = 0.d0
         cst(1,106) = 0.d0
      else
         cst(1,102) = 0.899d0
         cst(1,103) = 0.d0
         cst(1,104) = -1.44d0
         cst(1,105) = 0.d0
         cst(1,106) = 0.d0
      endif
      do ii = 102, 106
         cst(2,ii) = cst(1,ii)
      enddo

C     Prediction-only controls.  INT=1 retains eigenvectors for INTEN;
C     no observations are read and no fit is performed.
      np(1) = 0
      np(2) = 0
      ito = 0
      ktronc = 10
      ktron1 = 11
      conv = 1.d0
      irmw = 1
      jmax = output_jmax
      vtmax = output_vtmax
      int = 1
      unc = 0
      ibug1 = 0
      ibug2 = 0
      ibug3 = 0
      ibug4 = 0
      ibug5 = 0
      ibgtme = 0
      ibgvtd = 0
      evib(1) = 0.d0
      evib(2) = 0.d0
      term = 0.d0
      en0 = 0.d0
      label = 1
C     RINT is not used, but initialize its historical Boltzmann commons so
C     the driver has no undefined state.  SSQ itself is temperature-free.
      tpr = 300.d0
      conv2 = 1.4387768775039338d0
      bltz2 = 1.d0

      call setup(cst)
      call factt

      call emit_levels(output_jmax,output_kmax,output_vtmax)

      write(6,'(A)') 'BEGIN_FROZEN_STRENGTHS_CSV'
      write(6,'(A)')
     & 'species,upper_vt,upper_j,upper_k,upper_parity,'//
     & 'lower_vt,lower_j,lower_k,lower_parity,'//
     & 'upper_energy_raw_cm1,lower_energy_raw_cm1,frequency_mhz,'//
     & 'smu2_debye2,dipole_components,contribution_a_percent,'//
     & 'contribution_b_percent,cancellation_factor'

      candidates = 0
      retained = 0
      count_positive = 0
      count_1em16 = 0
      count_1em14 = 0
      count_1em12 = 0
      count_1em10 = 0
      count_1em8 = 0
      call emit_species('A',cst,output_jmax,output_kmax,
     &     output_vtmax,smu2_min,cross_vibrational_mode,
     &     candidates,retained,
     &     count_positive,count_1em16,count_1em14,count_1em12,
     &     count_1em10,count_1em8)
      call emit_species('E',cst,output_jmax,output_kmax,
     &     output_vtmax,smu2_min,cross_vibrational_mode,
     &     candidates,retained,
     &     count_positive,count_1em16,count_1em14,count_1em12,
     &     count_1em10,count_1em8)
      write(6,'(A)') 'END_FROZEN_STRENGTHS_CSV'

C     The fourth field counts direct native-INTEN evaluations.  No surrogate
C     contraction is used, hence its difference from native INTEN is zero.
      write(6,200) 'FROZEN_STRENGTH_SUMMARY',candidates,retained,
     &             candidates,0.d0
      write(6,201) 'FROZEN_STRENGTH_FLOOR_COUNTS',count_positive,
     &             count_1em16,count_1em14,count_1em12,count_1em10,
     &             count_1em8

 200  format(A,',',I0,',',I0,',',I0,',',ES25.17E3)
 201  format(A,',',I0,',',I0,',',I0,',',I0,',',I0,',',I0)
      end


      subroutine emit_levels(output_jmax,output_kmax,output_vtmax)
      implicit none
      include 'common_dim.f'
      include 'common_big2.f'
      integer output_jmax, output_kmax, output_vtmax
      integer j, k, vt, ipr, ink
      character*1 parity_label
      double precision raw_energy, species_zero_a, species_zero_e

      species_zero_a = ener(1,1,1)
      species_zero_e = ener(6,2,1)
      write(6,'(A)') 'BEGIN_FROZEN_LEVELS_CSV'
      write(6,'(A)')
     & 'species,vt,j,k,parity,energy_raw_cm1,energy_species_cm1'

      do vt = 0, output_vtmax
         do j = 0, output_jmax
            do k = 0, min(j, output_kmax)
               if (k .eq. 0) then
                  ipr = 1
                  parity_label = '+'
                  ink = j**2 + j + k + 1 - ipr*j
                  raw_energy = ener(vt+1,ink,1)
                  write(6,100) 'A',vt,j,k,parity_label,raw_energy,
     &                 raw_energy-species_zero_a
               else
                  do ipr = 0, 1
                     if (ipr .eq. 0) parity_label = '-'
                     if (ipr .eq. 1) parity_label = '+'
                     ink = j**2 + j + k + 1 - ipr*j
                     raw_energy = ener(vt+1,ink,1)
                     write(6,100) 'A',vt,j,k,parity_label,raw_energy,
     &                    raw_energy-species_zero_a
                  enddo
               endif
            enddo
         enddo
      enddo

      do vt = 0, output_vtmax
         do j = 0, output_jmax
            do k = -min(j,output_kmax), min(j,output_kmax)
               ink = j**2 + j + k + 1
               raw_energy = ener(vt+6,ink,1)
               write(6,101) 'E',vt,j,k,raw_energy,
     &              raw_energy-species_zero_e
            enddo
         enddo
      enddo
      write(6,'(A)') 'END_FROZEN_LEVELS_CSV'

 100  format(A1,',',I0,',',I0,',',I0,',',A1,',',
     &       ES25.17E3,',',ES25.17E3)
 101  format(A1,',',I0,',',I0,',',I0,',,',
     &       ES25.17E3,',',ES25.17E3)
      end


      subroutine emit_species(species,cst,output_jmax,output_kmax,
     &     output_vtmax,smu2_min,cross_vibrational_mode,
     &     candidates,retained,
     &     count_positive,count_1em16,count_1em14,count_1em12,
     &     count_1em10,count_1em8)
      implicit none
      include 'common_dim.f'
      include 'common_cst.f'
      include 'common_expdat.f'
      include 'common_big2.f'
      include 'common_quant.f'
      include 'common_int.f'

      integer maxlev
      parameter (maxlev=2500)
      character*1 species, spar(maxlev)
      character*2 components
      character*6 muabc
      integer output_jmax,output_kmax,output_vtmax
      integer candidates,retained,count_positive,count_1em16
      integer count_1em14,count_1em12,count_1em10,count_1em8
      integer nlev,vt,j,k,ipr,ia,ib,ilow,iup
      integer svt(maxlev),sj(maxlev),sk(maxlev),sipr(maxlev)
      integer sivv(maxlev),sink(maxlev),ivvpair(2),inkpair(2)
      integer ictt(6)
      double precision se(maxlev),smu2_min,rint,cf,frequency
      logical cross_vibrational_mode

      nlev = 0
      if (species .eq. 'A') then
         isig = 0
         iisig = 1
         do vt = 0, output_vtmax
            do j = 0, output_jmax
               do k = 0, min(j,output_kmax)
                  if (k .eq. 0) then
                     nlev = nlev + 1
                     if (nlev .gt. maxlev) stop 5
                     ipr = 1
                     svt(nlev)=vt
                     sj(nlev)=j
                     sk(nlev)=k
                     sipr(nlev)=ipr
                     spar(nlev)='+'
                     sivv(nlev)=vt+1
                     sink(nlev)=j**2+j+k+1-ipr*j
                     se(nlev)=ener(sivv(nlev),sink(nlev),1)
                  else
                     do ipr = 0, 1
                        nlev = nlev + 1
                        if (nlev .gt. maxlev) stop 5
                        svt(nlev)=vt
                        sj(nlev)=j
                        sk(nlev)=k
                        sipr(nlev)=ipr
                        if (ipr .eq. 0) spar(nlev)='-'
                        if (ipr .eq. 1) spar(nlev)='+'
                        sivv(nlev)=vt+1
                        sink(nlev)=j**2+j+k+1-ipr*j
                        se(nlev)=ener(sivv(nlev),sink(nlev),1)
                     enddo
                  endif
               enddo
            enddo
         enddo
      else
         isig = 1
         iisig = 2
         do vt = 0, output_vtmax
            do j = 0, output_jmax
               do k = -min(j,output_kmax), min(j,output_kmax)
                  nlev = nlev + 1
                  if (nlev .gt. maxlev) stop 5
                  svt(nlev)=vt
                  sj(nlev)=j
                  sk(nlev)=k
                  sipr(nlev)=0
                  spar(nlev)=' '
                  sivv(nlev)=vt+6
                  sink(nlev)=j**2+j+k+1
                  se(nlev)=ener(sivv(nlev),sink(nlev),1)
               enddo
            enddo
         enddo
      endif

C     Enumerate unordered pairs, enforce only the exact rank-1 rotational
C     rule, and let the full mixed eigenvectors decide whether strength is
C     carried.  Nominal K and vt labels are never used as prefilters.
      do ia = 1, nlev
         do ib = ia, nlev
            if (ia .eq. ib .and. .not. cross_vibrational_mode)
     &         goto 20
            if (abs(sj(ia)-sj(ib)) .gt. 1) goto 20
            if (sj(ia) .eq. 0 .and. sj(ib) .eq. 0) goto 20
            candidates = candidates + 1
            if (se(ia) .lt. se(ib)) then
               ilow = ia
               iup = ib
            else
               ilow = ib
               iup = ia
            endif
            if (se(iup) .lt. se(ilow)) goto 20
            if (se(iup) .eq. se(ilow) .and.
     &          .not. cross_vibrational_mode) goto 20

            isig = 0
            iisig = 1
            if (species .eq. 'E') then
               isig = 1
               iisig = 2
            endif
            nobs(1,1)=sj(ilow)
            nobs(2,1)=sj(iup)
            kobs(1,1)=sk(ilow)
            kobs(2,1)=sk(iup)
            ivvpair(1)=sivv(ilow)
            ivvpair(2)=sivv(iup)
            inkpair(1)=sink(ilow)
            inkpair(2)=sink(iup)
            bltz2=1.d0
            call inten(rint,cst,1,ivvpair,inkpair,muabc,ictt,cf)

            if (ssq .gt. 0.d0) count_positive=count_positive+1
            if (ssq .ge. 1.d-16) count_1em16=count_1em16+1
            if (ssq .ge. 1.d-14) count_1em14=count_1em14+1
            if (ssq .ge. 1.d-12) count_1em12=count_1em12+1
            if (ssq .ge. 1.d-10) count_1em10=count_1em10+1
            if (ssq .ge. 1.d-8) count_1em8=count_1em8+1
            if (ssq .lt. smu2_min .or. ssq .le. 0.d0) goto 20

            retained = retained + 1
            frequency=(se(iup)-se(ilow))*29979.2458d0
            components = '  '
            if (muabc(1:1) .eq. 'A') components(1:1)='A'
            if (muabc(2:2) .eq. 'B') components(2:2)='B'
            write(6,100) species,svt(iup),sj(iup),sk(iup),spar(iup),
     &           svt(ilow),sj(ilow),sk(ilow),spar(ilow),se(iup),
     &           se(ilow),frequency,ssq,components,ictt(1),ictt(2),cf
 20         continue
         enddo
      enddo

 100  format(A1,',',I0,',',I0,',',I0,',',A1,',',
     &       I0,',',I0,',',I0,',',A1,',',ES25.17E3,',',
     &       ES25.17E3,',',ES25.17E3,',',ES25.17E3,',',A2,',',
     &       I0,',',I0,',',ES25.17E3)
      end
