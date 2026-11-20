      SUBROUTINE INTEN_CACHED(RINT,CST,K,IVV,INK,MUABC,ICTT,CF)
      implicit none
      include 'common_dim.f'
      include 'common_cst.f'
      include 'common_eigen.f'
      include 'common_tor59.f'
      include 'common_rotor.f'
      include 'common_quant.f'
      include 'common_expdat.f'
      include 'common_int.f'
      integer k,ivv(2),ink(2),ictt(6),jl,ju
      character*6 muabc
      double precision rint,cf,repair_a,repair_b,repair_component_cf
      double precision q,denom
      common/REPAIR_INTEN_COMPONENTS/repair_a,repair_b,
     &     repair_component_cf
      if(cst(1,103).ne.0.d0.or.cst(1,105).ne.0.d0.or.
     &   cst(1,106).ne.0.d0) then
         write(0,*) 'cached INTEN supports declared constant a/b only'
         stop 73
      endif
      jl=nobs(1,k)
      ju=nobs(2,k)
      call repair_cached_components(jl,ju,isig+1,ndimto,
     & ivv(1),ivv(2),ink(1),ink(2),a59,vectpr,cst(1,102),
     & cst(1,104),repair_a,repair_b,cf,repair_component_cf)
      ssq=(repair_a+repair_b)**2
      rint=bltz2*ssq
      q=dble((2*jl+1)*(2*ju+1))
      muabc='      '
      ictt=0
      if(abs(repair_a)/sqrt(q).gt.1.d-10) muabc(1:1)='A'
      if(abs(repair_b)/sqrt(q).gt.1.d-10) muabc(2:2)='B'
      denom=abs(repair_a)+abs(repair_b)
      if(ssq/q.gt.1.d-10.and.denom.gt.0.d0) then
         ictt(1)=nint(100.d0*repair_a/denom)
         ictt(2)=nint(100.d0*repair_b/denom)
      else
         cf=2.d0
      endif
      END
