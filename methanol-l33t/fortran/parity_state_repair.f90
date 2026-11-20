! Explicit physical-parity resolution and globally bijective A-state export.
! This is an opt-in prediction overlay; the frozen vendor sources are unchanged.
module methanol_parity_repair
  implicit none
  integer, parameter :: maxdim=837, maxk=93, maxtor=41
  integer, save :: previous_j=-1, mapped_j=-1, state_generation=0
  integer, save :: mapping(9,maxk)=0
  real(8), allocatable, save :: previous_vectors(:,:), previous_a(:,:,:)
contains
  subroutine symmetric_eigh(x,e)
    real(8), intent(inout) :: x(:,:)
    real(8), intent(out) :: e(:)
    real(8), allocatable :: w(:)
    real(8) :: query(1)
    integer :: info,n,lw
    n=size(x,1)
    call dsyev('V','U',n,x,n,e,query,-1,info)
    if(info/=0) error stop 'parity repair: LAPACK workspace query failed'
    lw=max(3*n-1,int(query(1)))
    allocate(w(lw))
    call dsyev('V','U',n,x,n,e,w,lw,info)
    if(info/=0) error stop 'parity repair: symmetric eigensolve failed'
  end subroutine

  real(8) function frobenius(x)
    real(8),intent(in) :: x(:,:)
    frobenius=sqrt(sum(x*x))
  end function

  subroutine parity_eigh(h,p,u,e,par,metrics)
    real(8),intent(in) :: h(:,:),p(:,:)
    real(8),intent(out) :: u(:,:),e(:),metrics(8)
    integer,intent(out) :: par(:)
    real(8),allocatable :: q(:,:),pe(:),hp(:,:),block(:,:),be(:),ident(:,:),res(:,:),tmp(:)
    integer :: n,minus,first,last,nb,sgn,i,j,idx,ps
    real(8) :: hscale,ev
    n=size(h,1)
    allocate(q(n,n),pe(n),ident(n,n),res(n,n),tmp(n))
    ident=0d0
    do i=1,n
      ident(i,i)=1d0
    end do
    metrics(1)=maxval(abs(p-transpose(p)))
    metrics(2)=maxval(abs(matmul(p,p)-ident))
    if(metrics(1)>1d-10.or.metrics(2)>1d-10) error stop 'parity repair: invalid physical parity operator'
    q=0.5d0*(p+transpose(p))
    call symmetric_eigh(q,pe)
    if(maxval(abs(abs(pe)-1d0))>1d-10) error stop 'parity repair: nonbinary parity eigenvalues'
    minus=count(pe<0d0)
    hscale=max(1d0,frobenius(h))
    hp=matmul(transpose(q),matmul(h,q))
    metrics(3)=sqrt(2d0)*frobenius(hp(:minus,minus+1:))/hscale
    if(metrics(3)>1d-11) error stop 'parity repair: Hamiltonian violates physical parity'
    do sgn=-1,1,2
      if(sgn==-1) then
        first=1;last=minus
      else
        first=minus+1;last=n
      end if
      nb=last-first+1
      if(nb==0) cycle
      allocate(block(nb,nb),be(nb))
      block=hp(first:last,first:last)
      call symmetric_eigh(block,be)
      u(:,first:last)=matmul(q(:,first:last),block)
      e(first:last)=be
      par(first:last)=sgn
      deallocate(block,be)
    end do
    ! Stable energy ordering is a storage convention, never a parity assignment.
    do i=2,n
      ev=e(i);ps=par(i);tmp=u(:,i);j=i-1
      do while(j>=1)
        if(e(j)<=ev) exit
        e(j+1)=e(j);par(j+1)=par(j);u(:,j+1)=u(:,j);j=j-1
      end do
      e(j+1)=ev;par(j+1)=ps;u(:,j+1)=tmp
    end do
    do i=1,n
      idx=maxloc(abs(u(:,i)),dim=1)
      if(u(idx,i)<0d0) u(:,i)=-u(:,i)
    end do
    metrics(4)=maxval(abs(matmul(transpose(u),u)-ident))
    res=matmul(h,u)
    do i=1,n
      res(:,i)=res(:,i)-e(i)*u(:,i)
    end do
    metrics(5)=frobenius(res)/hscale
    metrics(6)=maxval(sqrt(sum(res*res,dim=1)))
    res=matmul(p,u)
    do i=1,n
      res(:,i)=res(:,i)-par(i)*u(:,i)
    end do
    metrics(7)=maxval(abs(res))
    metrics(8)=maxval(abs(h-transpose(h)))/hscale
    if(metrics(4)>1d-10.or.metrics(5)>1d-11.or.metrics(7)>1d-9) &
      error stop 'parity repair: repaired eigenbasis failed original-problem validation'
  end subroutine

  subroutine maximum_assignment(weight,assignment)
    ! Rectangular Hungarian primal-dual algorithm specialized to a square
    ! full-state assignment. No greedy pickup and no duplicate eigenvectors.
    real(8),intent(in) :: weight(:,:)
    integer,intent(out) :: assignment(:)
    real(8),allocatable :: u(:),v(:),minimum(:)
    integer,allocatable :: matched(:),way(:)
    logical,allocatable :: used(:)
    real(8) :: delta,cur
    integer :: n,i,j,j0,j1,i0
    n=size(weight,1)
    allocate(u(0:n),v(0:n),minimum(0:n),matched(0:n),way(0:n),used(0:n))
    u=0d0;v=0d0;matched=0;way=0
    do i=1,n
      matched(0)=i;j0=0;minimum=huge(1d0);used=.false.
      do
        used(j0)=.true.;i0=matched(j0);delta=huge(1d0);j1=0
        do j=1,n
          if(used(j)) cycle
          cur=-weight(i0,j)-u(i0)-v(j)
          if(cur<minimum(j)) then
            minimum(j)=cur;way(j)=j0
          end if
          if(minimum(j)<delta) then
            delta=minimum(j);j1=j
          end if
        end do
        do j=0,n
          if(used(j)) then
            u(matched(j))=u(matched(j))+delta;v(j)=v(j)-delta
          else
            minimum(j)=minimum(j)-delta
          end if
        end do
        j0=j1
        if(matched(j0)==0) exit
      end do
      do
        j1=way(j0);matched(j0)=matched(j1);j0=j1
        if(j0==0) exit
      end do
    end do
    assignment=0
    do j=1,n
      assignment(matched(j))=j
    end do
  end subroutine

  subroutine native_parity(j,ndim,a,p)
    integer,intent(in) :: j,ndim
    real(8),intent(in) :: a(9,maxtor,maxk)
    real(8),intent(out) :: p(:,:)
    integer :: k,t,tp,q,i,l,nk
    real(8) :: phase
    nk=2*j+1;p=0d0
    do k=-j,j
      phase=1d0
      if(mod(j+k,2)/=0) phase=-1d0
      do t=1,9
        i=(t-1)*nk+k+j+1
        do tp=1,9
          l=(tp-1)*nk-k+j+1
          do q=1,ndim
            p(i,l)=p(i,l)+phase*a(t,q,k+j+1)*a(tp,ndim+1-q,-k+j+1)
          end do
        end do
      end do
    end do
  end subroutine

  subroutine assign_native(j,ndim,a,p,u,e,par)
    integer,intent(in) :: j,ndim,par(:)
    real(8),intent(in) :: a(9,maxtor,maxk),p(:,:),u(:,:),e(:)
    real(8),allocatable :: refs(:,:),weight(:,:),transport(:,:),ident(:,:),ov(:,:)
    integer,allocatable :: target_parity(:),assignment(:),seen(:)
    real(8) :: norm2,gram,prerr,transport_err,second,minimum_overlap,margin,minmargin
    integer :: n,nk,oldnk,vt,sk,k,t,tp,row,oldrow,i,col,idx,ps,unit,envstat
    character(2048) :: auditdir,filename
    n=size(e);nk=2*j+1
    allocate(refs(n,n),target_parity(n),assignment(n),seen(n),transport(9,9),ident(n,n))
    refs=0d0;ident=0d0;transport_err=0d0
    do i=1,n
      ident(i,i)=1d0
    end do
    if(j==0) then
      if(allocated(previous_vectors)) deallocate(previous_vectors,previous_a)
      previous_j=-1
      state_generation=state_generation+1
    else
      if(previous_j/=j-1) error stop 'parity repair: nonsequential A state history'
      oldnk=2*(j-1)+1
    end if
    do vt=0,8
      do sk=-j,j
        row=vt*nk+sk+j+1
        ps=1
        if(mod(j+vt,2)/=0) ps=-1
        if(sk<0) ps=-ps
        target_parity(row)=ps
        if(j>0.and.abs(sk)<j) then
          oldrow=vt*oldnk+sk+j
          do k=-j+1,j-1
            transport=matmul(a(:,1:ndim,k+j+1),transpose(previous_a(:,1:ndim,k+j)))
            do t=1,9
              do tp=1,9
                refs((t-1)*nk+k+j+1,row)=refs((t-1)*nk+k+j+1,row) &
                  +transport(t,tp)*previous_vectors((tp-1)*oldnk+k+j,oldrow)
              end do
            end do
            do t=1,9
              do tp=1,9
                norm2=dot_product(transport(:,t),transport(:,tp))
                if(t==tp) norm2=norm2-1d0
                transport_err=max(transport_err,abs(norm2))
              end do
            end do
          end do
        else
          ! The positive-K torsional basis vector and its exact parity image
          ! span the two new edge-K states. At J=0 this projects K=0 seeds.
          idx=vt*nk+j+j+1
          refs(:,row)=ps*p(:,idx)
          refs(idx,row)=refs(idx,row)+1d0
          norm2=sum(refs(:,row)**2)
          if(norm2<1d-8) error stop 'parity repair: deficient label seed projection'
          refs(:,row)=refs(:,row)/sqrt(norm2)
        end if
      end do
    end do
    gram=maxval(abs(matmul(transpose(refs),refs)-ident))
    ov=matmul(p,refs)
    do i=1,n
      ov(:,i)=ov(:,i)-target_parity(i)*refs(:,i)
    end do
    prerr=maxval(abs(ov))
    if(gram>1d-8.or.prerr>1d-8.or.transport_err>1d-8) &
      error stop 'parity repair: gauge transport or complete references failed'
    ov=matmul(transpose(refs),u)
    weight=ov*ov
    do i=1,n
      do col=1,n
        if(target_parity(i)/=par(col)) weight(i,col)=-1d6
      end do
    end do
    call maximum_assignment(weight,assignment)
    seen=0;minimum_overlap=1d0;minmargin=1d0
    call get_environment_variable('BELGI_REPAIR_AUDIT_DIR',auditdir,status=envstat)
    if(envstat==0.and.len_trim(auditdir)>0) then
      write(filename,'(A,"/J",I2.2,".assignment.csv")') trim(auditdir),j
      open(newunit=unit,file=trim(filename),status='replace',action='write')
      write(unit,'(A)') 'J,vt,signed_K,physical_parity,eigenindex,energy_cm1,overlap_squared,best_alternative_squared,assignment_margin'
    else
      unit=0
    end if
    do vt=0,8
      do sk=-j,j
        row=vt*nk+sk+j+1;col=assignment(row)
        if(col<1.or.col>n) error stop 'parity repair: missing state assignment'
        seen(col)=seen(col)+1
        if(par(col)/=target_parity(row)) error stop 'parity repair: forbidden assignment'
        mapping(vt+1,sk+47)=col
        minimum_overlap=min(minimum_overlap,weight(row,col))
        second=0d0
        do i=1,n
          if(i/=col) second=max(second,weight(row,i))
        end do
        margin=weight(row,col)-second;minmargin=min(minmargin,margin)
        if(unit/=0) write(unit,'(5(I0,","),4(ES25.17E3,:,","))') &
          j,vt,sk,target_parity(row),col,e(col),weight(row,col),second,margin
      end do
    end do
    if(any(seen/=1)) error stop 'parity repair: full-state assignment is not bijective'
    if(unit/=0) close(unit)
    write(0,'(A,I0,5(A,ES12.4E3))') 'PARITY_ASSIGNMENT J=',j, &
      ' gram=',gram,' parity=',prerr,' gauge=',transport_err, &
      ' min_overlap2=',minimum_overlap,' min_margin=',minmargin
    if(allocated(previous_vectors)) deallocate(previous_vectors,previous_a)
    allocate(previous_vectors(n,n),previous_a(9,maxtor,maxk))
    do row=1,n
      previous_vectors(:,row)=u(:,assignment(row))
    end do
    previous_a=a;previous_j=j;mapped_j=j
  end subroutine
end module

subroutine repair_diagonalize(j,isig,ndim,h,e,u,a)
  use methanol_parity_repair
  implicit none
  integer,intent(in) :: j,isig,ndim
  real(8),intent(in) :: h(maxdim,maxdim),a(9,maxtor,maxk)
  real(8),intent(out) :: e(maxdim),u(maxdim,maxdim)
  real(8),allocatable :: hsmall(:,:),p(:,:),usmall(:,:),esmall(:)
  integer,allocatable :: par(:)
  real(8) :: metrics(8)
  integer :: n,i
  n=9*(2*j+1)
  allocate(hsmall(n,n),p(n,n),usmall(n,n),esmall(n),par(n))
  hsmall=h(:n,:n)
  if(isig==0) then
    call native_parity(j,ndim,a,p)
    call parity_eigh(hsmall,p,usmall,esmall,par,metrics)
    call assign_native(j,ndim,a,p,usmall,esmall,par)
    write(0,'(A,I0,A,8(ES12.4E3,1X))') 'PARITY_EIGENSYSTEM J=',j,' metrics=',metrics
  else
    usmall=hsmall
    call symmetric_eigh(usmall,esmall)
  end if
  u(:n,:n)=usmall;e(:n)=esmall
end subroutine

integer function repair_state_index(j,vt,k,pr)
  use methanol_parity_repair
  implicit none
  integer,intent(in) :: j,vt,k
  character,intent(in) :: pr
  integer :: sk
  if(j/=mapped_j.or.vt<0.or.vt>8) error stop 'parity repair: unavailable global state map'
  sk=k
  if(pr=='-') sk=-k
  repair_state_index=mapping(vt+1,sk+47)
  if(repair_state_index<1) error stop 'parity repair: unassigned native state'
end function
