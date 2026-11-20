! Algebraically exact full signed constant-dipole contractions cached by J pair.
! No state-dependent pruning, parity-based deletion, or nominal label K filter.
module methanol_dipole_cache
  use methanol_parity_repair, only: state_generation
  implicit none
  type dipole_block
    real(8), allocatable :: ampa(:,:),ampb(:,:),absolute_sum(:,:)
  end type
  type(dipole_block), save :: blocks(0:46,0:46,2)
  integer, save :: cached_generation=-1
  real(8), save :: cached_mua=huge(1d0),cached_mub=huge(1d0)
contains
  subroutine reset_cache(mua,mub)
    real(8),intent(in) :: mua,mub
    integer :: j,k,s
    if(cached_generation==state_generation.and.mua==cached_mua.and.mub==cached_mub) return
    do s=1,2
      do j=0,46
        do k=0,46
          if(allocated(blocks(j,k,s)%ampa)) &
            deallocate(blocks(j,k,s)%ampa,blocks(j,k,s)%ampb,blocks(j,k,s)%absolute_sum)
        end do
      end do
    end do
    cached_generation=state_generation;cached_mua=mua;cached_mub=mub
  end subroutine

  subroutine fill_block(jl,ju,sigma,ndim,a,v,mua,mub)
    integer,intent(in) :: jl,ju,sigma,ndim
    real(8),intent(in) :: a(9,41,93,47,2),v(10,2209,93,9),mua,mub
    real(8),allocatable :: da(:,:),db(:,:),vl(:,:),vu(:,:),inter(:,:),absd(:,:),absl(:,:),absu(:,:)
    real(8) :: angular,overlap,phase,qfactor,ad3j
    external :: ad3j
    integer :: nkl,nku,nl,nu,sl,su,kl,ku,delta,t,tp,i,k,col,vt,local
    nkl=2*jl+1;nku=2*ju+1;nl=9*nkl;nu=9*nku;sl=5*nkl;su=5*nku
    allocate(da(nl,nu),db(nl,nu),vl(nl,sl),vu(nu,su),inter(nl,su))
    da=0d0;db=0d0
    qfactor=sqrt(real((2*jl+1)*(2*ju+1),8))
    do kl=-jl,jl
      phase=1d0
      if(mod(jl-kl,2)/=0) phase=-1d0
      do delta=-1,1
        ku=kl+delta
        if(abs(ku)>ju) cycle
        angular=ad3j(2*jl,2,2*ju,-2*kl,-2*delta)*phase*qfactor
        if(delta/=0) angular=angular*delta/sqrt(2d0)
        do t=1,9
          i=(t-1)*nkl+kl+jl+1
          do tp=1,9
            k=(tp-1)*nku+ku+ju+1
            overlap=dot_product(a(t,1:ndim,kl+jl+1,jl+1,sigma),a(tp,1:ndim,ku+ju+1,ju+1,sigma))
            if(delta==0) then
              da(i,k)=angular*overlap*mua
            else
              db(i,k)=angular*overlap*mub
            end if
          end do
        end do
      end do
    end do
    do vt=0,4
      do local=1,nkl
        col=vt*nkl+local
        do t=1,9
          vl((t-1)*nkl+1:t*nkl,col)=v(vt+1+5*(sigma-1),jl*jl+local,1:nkl,t)
        end do
      end do
      do local=1,nku
        col=vt*nku+local
        do t=1,9
          vu((t-1)*nku+1:t*nku,col)=v(vt+1+5*(sigma-1),ju*ju+local,1:nku,t)
        end do
      end do
    end do
    allocate(blocks(jl,ju,sigma)%ampa(sl,su),blocks(jl,ju,sigma)%ampb(sl,su),blocks(jl,ju,sigma)%absolute_sum(sl,su))
    call dgemm('N','N',nl,su,nu,1d0,da,nl,vu,nu,0d0,inter,nl)
    call dgemm('T','N',sl,su,nl,1d0,vl,nl,inter,nl,0d0,blocks(jl,ju,sigma)%ampa,sl)
    call dgemm('N','N',nl,su,nu,1d0,db,nl,vu,nu,0d0,inter,nl)
    call dgemm('T','N',sl,su,nl,1d0,vl,nl,inter,nl,0d0,blocks(jl,ju,sigma)%ampb,sl)
    ! Cowan's basis-contribution cancellation denominator, exactly matching
    ! unpruned native INTEN for the constant a/b operator (disjoint K support).
    absd=abs(da+db);absl=abs(vl);absu=abs(vu)
    call dgemm('N','N',nl,su,nu,1d0,absd,nl,absu,nu,0d0,inter,nl)
    call dgemm('T','N',sl,su,nl,1d0,absl,nl,inter,nl,0d0,blocks(jl,ju,sigma)%absolute_sum,sl)
  end subroutine

  subroutine cached_components(jl,ju,sigma,ndim,ivlo,ivup,inklo,inkup,a,v,mua,mub,ampa,ampb,cowan,component)
    integer,intent(in) :: jl,ju,sigma,ndim,ivlo,ivup,inklo,inkup
    real(8),intent(in) :: a(9,41,93,47,2),v(10,2209,93,9),mua,mub
    real(8),intent(out) :: ampa,ampb,cowan,component
    integer :: il,iu
    real(8) :: denom
    call reset_cache(mua,mub)
    if(.not.allocated(blocks(jl,ju,sigma)%ampa)) call fill_block(jl,ju,sigma,ndim,a,v,mua,mub)
    il=(ivlo-1-5*(sigma-1))*(2*jl+1)+inklo-jl*jl
    iu=(ivup-1-5*(sigma-1))*(2*ju+1)+inkup-ju*ju
    ampa=blocks(jl,ju,sigma)%ampa(il,iu);ampb=blocks(jl,ju,sigma)%ampb(il,iu)
    denom=blocks(jl,ju,sigma)%absolute_sum(il,iu)
    cowan=0d0;component=0d0
    if(denom>0d0) cowan=((ampa+ampb)/denom)**2
    denom=abs(ampa)+abs(ampb)
    if(denom>0d0) component=((ampa+ampb)/denom)**2
  end subroutine
end module

subroutine repair_cached_components(jl,ju,sigma,ndim,ivlo,ivup,inklo,inkup,a,v,mua,mub,ampa,ampb,cowan,component)
  use methanol_dipole_cache
  implicit none
  integer,intent(in) :: jl,ju,sigma,ndim,ivlo,ivup,inklo,inkup
  real(8),intent(in) :: a(9,41,93,47,2),v(10,2209,93,9),mua,mub
  real(8),intent(out) :: ampa,ampb,cowan,component
  call cached_components(jl,ju,sigma,ndim,ivlo,ivup,inklo,inkup,a,v,mua,mub,ampa,ampb,cowan,component)
end subroutine
