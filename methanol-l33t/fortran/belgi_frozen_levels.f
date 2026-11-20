      program belgi_frozen_levels
C
C     Prediction-only driver for the BELGI Hamiltonian.  The frozen
C     Hamiltonian constants are read directly from a complete canonical
C     205-slot vector on standard input.  No observations or fit are used.
C
      implicit none

      include 'common_dim.f'
      include 'common_cst.f'
      include 'common_dat.f'
      include 'common_par.f'
      include 'common_expdat.f'
      include 'common_big2.f'
      include 'common_cntl.f'

      integer argc, output_jmax, output_kmax, output_vtmax
      integer j, k, vt, ipr, ink, ii, ios
      character*64 arg, vector_header
      character*1 parity_label
      double precision raw_energy, species_zero_a, species_zero_e

      argc = command_argument_count()
      if (argc .ne. 3) then
         write(0,*) 'usage: belgi_frozen_levels JMAX KMAX VTMAX'
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

      if (output_jmax .lt. 0 .or. output_jmax .gt. 46) stop 2
      if (output_kmax .lt. 0 .or.
     &    output_kmax .gt. output_jmax) stop 2
      if (output_vtmax .lt. 0 .or. output_vtmax .gt. 4) stop 2

C     Python has already validated the names and serialized every CST slot,
C     including explicit zeros.  Direct loading avoids DATA11's scratch-file
C     round trip and its historical AMUB3 initialization typo.
      read(5,'(A)',iostat=ios) vector_header
      if (ios .ne. 0 .or.
     &    vector_header .ne. 'BELGI_FROZEN_CST_205') stop 3
      do ii = 1, cstm
         read(5,*,iostat=ios) cst(1,ii)
         if (ios .ne. 0) stop 3
         cst(2,ii) = cst(1,ii)
      enddo

C     Reproduce the controls recorded in the Xu GLOBALFIT output, except
C     that fitting and intensity calculation are deliberately disabled.
      np(1) = 0
      np(2) = 0
      ito = 0
      ktronc = 10
      ktron1 = 11
      conv = 1.d0
      irmw = 1
      jmax = output_jmax
      vtmax = output_vtmax
      int = 0
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
C     Xu's final listing used ENCAL2's high-J continuity pickup scheme.
      label = 1

      call setup(cst)

      species_zero_a = ener(1,1,1)
C     The lowest E-species level is vt=0, J=1, signed K=-1 (INK=2).
      species_zero_e = ener(6,2,1)

      write(6,'(A)') 'BEGIN_FROZEN_LEVELS_CSV'
      write(6,'(A)')
     & 'species,vt,j,k,parity,energy_raw_cm1,energy_species_cm1'

C     A symmetry: non-negative K and explicit parity.  K=0 exists only
C     with positive parity in this representation.
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

C     E symmetry: signed K and no parity label.
      do vt = 0, output_vtmax
         do j = 0, output_jmax
            do k = -min(j, output_kmax), min(j, output_kmax)
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
