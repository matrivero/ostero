module mod_fortran

!inputs deriv_and_detjac_calc
integer*4 :: num_elements_bulk
real*8, allocatable, dimension(:,:) :: nodes
integer*4, allocatable, dimension(:,:) :: elements_bulk

!inputs vmass_calc
integer*4 :: num_nodes
real*8, allocatable, dimension(:,:,:) :: pgauss
real*8, allocatable, dimension(:,:) :: weight

!inputs assembly
real*8, allocatable, dimension(:) :: displ
real*8, allocatable, dimension(:) :: lame1_mu
real*8, allocatable, dimension(:) :: lame2_lambda
real*8, allocatable, dimension(:) :: plane
logical :: stress_calc_on

!outputs deriv_and_detjac_calc
real*8, allocatable, dimension(:,:) :: det_jac
real*8, allocatable, dimension(:,:,:,:) :: deriv

!output vmass_calc
real*8, allocatable, dimension(:) :: vmass

!outputs assembly
real*8, allocatable, dimension(:,:) :: k_tot
real*8, allocatable, dimension(:) :: r_tot
real*8, allocatable, dimension(:,:) :: strain
real*8, allocatable, dimension(:,:) :: stress

character(len=4) :: model

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init()

allocate(weight(4,2))
weight = 0.0D0
weight(1,1) = 1.0/6.0
weight(2,1) = 1.0/6.0
weight(3,1) = 1.0/6.0

weight(1,2) = 1.0
weight(2,2) = 1.0
weight(3,2) = 1.0
weight(4,2) = 1.0

allocate(pgauss(4,2,2))
pgauss = 0.0D0
pgauss(1,1,1) = 1.0/6.0
pgauss(1,2,1) = 1.0/6.0
pgauss(2,1,1) = 2.0/3.0
pgauss(2,2,1) = 1.0/6.0
pgauss(3,1,1) = 1.0/6.0
pgauss(3,2,1) = 2.0/3.0

pgauss(1,1,2) = -0.577350269189626 
pgauss(1,2,2) = -0.577350269189626 
pgauss(2,1,2) = 0.577350269189626 
pgauss(2,2,2) = -0.577350269189626 
pgauss(3,1,2) = 0.577350269189626 
pgauss(3,2,2) = 0.577350269189626 
pgauss(4,1,2) = -0.577350269189626 
pgauss(4,2,2) = 0.577350269189626 

allocate(plane(num_elements_bulk))
allocate(lame1_mu(num_elements_bulk))
allocate(lame2_lambda(num_elements_bulk))

allocate(nodes(num_nodes,4))
allocate(elements_bulk(num_elements_bulk,9))

allocate(displ(num_nodes*2))

if (model == 'ISOL') print *,"ISOLINEAL MATERIAL MODEL" 
if (model == 'BELY') print *,"BELYTSCHKO / NEO-HOOKEAN MATERIAL MODEL" 
if (model == 'ZIEN') print *,"ZIENKIEWICZ / NEO-HOOKEAN MATERIAL MODEL" 
if (model == 'LAUR') print *,"LAURSEN / NEO-HOOKEAN MATERIAL MODEL" 

end subroutine init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine derivs(pnode,igauss,dh)

implicit none

integer*4, intent(in) :: pnode,igauss
real*8, dimension(4,2), intent(out) :: dh

dh = 0.0D0
if (pnode == 3) then
  dh(1,1) = -1 !dh1/deta
  dh(1,2) = -1 !dh1/dxi
  dh(2,1) = 1 !dh2/deta
  dh(2,2) = 0 !dh2/dxi
  dh(3,1) = 0 !dh3/deta
  dh(3,2) = 1 !dh3/dxi
else if(pnode == 4) then
  dh(1,1) = -(1 - pgauss(igauss,2,2))/4.0 !dh1/eta
  dh(1,2) = -(1 - pgauss(igauss,1,2))/4.0 !dh1/xi
  dh(2,1) = (1 - pgauss(igauss,2,2))/4.0  !dh2/eta
  dh(2,2) = -(1 + pgauss(igauss,1,2))/4.0  !dh2/xi
  dh(3,1) = (1 + pgauss(igauss,2,2))/4.0  !dh3/eta
  dh(3,2) = (1 + pgauss(igauss,1,2))/4.0  !dh3/xi
  dh(4,1) = -(1 + pgauss(igauss,2,2))/4.0  !dh4/deta
  dh(4,2) = (1 - pgauss(igauss,1,2))/4.0 !dh4/dxi
end if

end subroutine derivs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shapefunc(pnode,igauss,gpsha)

implicit none

integer*4, intent(in) :: pnode,igauss
real*8, dimension(4), intent(out) :: gpsha

gpsha = 0.0D0
if (pnode == 3) then
    gpsha(1) = 1 - pgauss(igauss,1,pnode-2) - pgauss(igauss,2,pnode-2) !1-eta-xi 
    gpsha(2) = pgauss(igauss,1,pnode-2) !eta
    gpsha(3) = pgauss(igauss,2,pnode-2) !xi
else if(pnode == 4) then
    gpsha(1) = (1 - pgauss(igauss,1,pnode-2))*(1 - pgauss(igauss,2,pnode-2))/4.0
    gpsha(2) = (1 + pgauss(igauss,1,pnode-2))*(1 - pgauss(igauss,2,pnode-2))/4.0
    gpsha(3) = (1 + pgauss(igauss,1,pnode-2))*(1 + pgauss(igauss,2,pnode-2))/4.0
    gpsha(4) = (1 - pgauss(igauss,1,pnode-2))*(1 + pgauss(igauss,2,pnode-2))/4.0
end if

end subroutine shapefunc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine deriv_and_detjac_calc()

implicit none

real*8, dimension(4,2) :: coord_nodes
real*8, dimension(2,2) :: inv_jac
real*8, dimension(2,2) :: jac
real*8, dimension(4,2) :: dh
real*8 :: det_jac_tmp
integer*4 :: e,i,j
integer*4 :: pnode,ngauss
logical :: inv_flag

allocate(det_jac(num_elements_bulk,4))
allocate(deriv(4,2,num_elements_bulk,4))

do e = 1,num_elements_bulk
  if (elements_bulk(e,2) == 2) then !TRIANGLE ELEMENT
    pnode = 3
    ngauss = 3
    coord_nodes(1,1) = nodes(elements_bulk(e,6),2)
    coord_nodes(1,2) = nodes(elements_bulk(e,6),3)
    coord_nodes(2,1) = nodes(elements_bulk(e,7),2)
    coord_nodes(2,2) = nodes(elements_bulk(e,7),3)
    coord_nodes(3,1) = nodes(elements_bulk(e,8),2)
    coord_nodes(3,2) = nodes(elements_bulk(e,8),3)
  else if (elements_bulk(e,2) == 3) then !QUAD ELEMENT
    pnode = 4
    ngauss = 4
    coord_nodes(1,1) = nodes(elements_bulk(e,6),2)
    coord_nodes(1,2) = nodes(elements_bulk(e,6),3)
    coord_nodes(2,1) = nodes(elements_bulk(e,7),2)
    coord_nodes(2,2) = nodes(elements_bulk(e,7),3)
    coord_nodes(3,1) = nodes(elements_bulk(e,8),2)
    coord_nodes(3,2) = nodes(elements_bulk(e,8),3)
    coord_nodes(4,1) = nodes(elements_bulk(e,9),2)
    coord_nodes(4,2) = nodes(elements_bulk(e,9),3)
  end if

  do i = 1,ngauss !GAUSS POINT LOOP
    call derivs(pnode,i,dh)
    jac = 0.0D0
    do j = 1,pnode
      jac(1,1) = jac(1,1) + coord_nodes(j,1)*dh(j,1)
      jac(1,2) = jac(1,2) + coord_nodes(j,2)*dh(j,1)
      jac(2,1) = jac(2,1) + coord_nodes(j,1)*dh(j,2)
      jac(2,2) = jac(2,2) + coord_nodes(j,2)*dh(j,2)
    end do
    call m22inv(jac, inv_jac, det_jac_tmp, inv_flag)
    if (.not. inv_flag) then
      write(*,*) "Jacobian Matrix is not invertible... bye!"
      stop
    end if
    det_jac(e,i) = abs(det_jac_tmp)
    do j = 1,pnode
      deriv(j,1,e,i) = dh(j,1)*inv_jac(1,1) + dh(j,2)*inv_jac(1,2)
      deriv(j,2,e,i) = dh(j,1)*inv_jac(2,1) + dh(j,2)*inv_jac(2,2)
    end do 
  end do
end do

return

end subroutine deriv_and_detjac_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine vmass_calc()

implicit none

integer*4 :: e,i,inode,ipoin,pnode,ngauss
real*8, dimension(4) :: numer,gpsha
real*8 :: gpvol

allocate(vmass(num_nodes))

do e = 1,num_elements_bulk

  if (elements_bulk(e,2) == 2) then !TRIANGLE ELEMENT
    pnode = 3
    ngauss = 3
  else if (elements_bulk(e,2) == 3) then !QUAD ELEMENT
    pnode = 4
    ngauss = 4
  end if

  numer = 0
  do i = 1,ngauss !GAUSS POINT LOOP
    gpvol = 0.0D0
    gpvol = det_jac(e,i)*weight(i,pnode-2)
    call shapefunc(pnode,i,gpsha)
    do inode = 1,pnode
      numer(inode) = numer(inode) + gpvol * gpsha(inode)
    end do
  end do
  do inode = 1,pnode
    ipoin = elements_bulk(e,5+inode)
    vmass(ipoin) = vmass(ipoin) + numer(inode)
  end do
end do

end subroutine vmass_calc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine m22inv(mat, inv_mat, det, ok_flag)

implicit none

real*8, dimension(2,2), intent(in)  :: mat
real*8, dimension(2,2), intent(out) :: inv_mat
real*8, intent(out) :: det
logical, intent(out) :: ok_flag
real*8, parameter :: eps = 1.0D-10

det = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)

if (abs(det) .le. eps) then
  inv_mat = 0.0D0
  ok_flag = .false.
  return
end if

inv_mat(1,1) = mat(2,2)/det 
inv_mat(1,2) = -mat(1,2)/det
inv_mat(2,1) = -mat(2,1)/det
inv_mat(2,2) = mat(1,1)/det

ok_flag = .true.

return

end subroutine m22inv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assembly()

implicit none

real*8, dimension(8) :: displ_ele
integer*4 :: e,i,j,k,l,m,n,o
real*8, dimension(2,2) :: inv_F
logical :: inv_flag
real*8, dimension(2,2) :: Cau
real*8, dimension(2,2) :: inv_Cau
real*8 :: det_Cau
real*8 :: dummr, traceE
integer*4, dimension(2,2) :: delta_kron
real*8, dimension(2,2) :: S
real*8, dimension(2,2,2,2) :: D
real*8, dimension(2,2,2,2) :: dPdF
integer*4 :: anode,bnode,inode,adofn,bdofn,idofn
integer*4 :: idime,jdime,kdime,ldime
real*8, dimension(8,8) :: k_elem
real*8, dimension(8) :: r_elem
real*8 :: det_F
real*8, dimension(2,2) :: F
real*8, dimension(2,2) :: P
real*8, dimension(2,2) :: Et
real*8, dimension(2,2) :: ST
real*8 :: gpinv
real*8 :: gpvol
real*8 :: xmean
real*8, dimension(4) :: gpsha
integer*4 :: ipoin, ivoigt 
integer*4, dimension(3,2) :: nvgij
integer*4 :: pnode,ngauss 

if (.not. stress_calc_on) then
allocate(k_tot(num_nodes*2,num_nodes*2))
allocate(r_tot(num_nodes*2))
k_tot = 0.0D0
r_tot = 0.0D0
end if

if (stress_calc_on) then
  allocate(stress(3,num_nodes))
  allocate(strain(3,num_nodes))
  stress = 0.0D0
  strain = 0.0D0
end if

delta_kron = 0
do i = 1,2
  delta_kron(i,i) = 1
end do

nvgij = 0
nvgij(1,1) = 1 
nvgij(1,2) = 1
nvgij(2,1) = 2
nvgij(2,2) = 2
nvgij(3,1) = 1
nvgij(3,2) = 2

displ_ele = 0.0D0

do e = 1,num_elements_bulk

  if (elements_bulk(e,2) == 2) then !TRIANGLE ELEMENT
    pnode = 3
    ngauss = 3
    displ_ele(1) = displ(elements_bulk(e,6)*2-1)
    displ_ele(2) = displ(elements_bulk(e,6)*2)
    displ_ele(3) = displ(elements_bulk(e,7)*2-1)
    displ_ele(4) = displ(elements_bulk(e,7)*2)
    displ_ele(5) = displ(elements_bulk(e,8)*2-1)
    displ_ele(6) = displ(elements_bulk(e,8)*2)
  else if (elements_bulk(e,2) == 3) then !QUAD ELEMENT
    pnode = 4
    ngauss = 4
    displ_ele(1) = displ(elements_bulk(e,6)*2-1)
    displ_ele(2) = displ(elements_bulk(e,6)*2)
    displ_ele(3) = displ(elements_bulk(e,7)*2-1)
    displ_ele(4) = displ(elements_bulk(e,7)*2)
    displ_ele(5) = displ(elements_bulk(e,8)*2-1)
    displ_ele(6) = displ(elements_bulk(e,8)*2)
    displ_ele(7) = displ(elements_bulk(e,9)*2-1)
    displ_ele(8) = displ(elements_bulk(e,9)*2)
  end if

  k_elem = 0.0D0
  r_elem = 0.0D0

  do i = 1,ngauss !GAUSS POINT LOOP
    !DEFORMATION GRADIENT TENSOR
    F = 0.0D0
    if (pnode == 3) then
      F(1,1) = 1 + (deriv(1,1,e,i)*displ_ele(1) + deriv(2,1,e,i)*displ_ele(3) + deriv(3,1,e,i)*displ_ele(5))
      F(1,2) = deriv(1,2,e,i)*displ_ele(1) + deriv(2,2,e,i)*displ_ele(3) + deriv(3,2,e,i)*displ_ele(5)
      F(2,1) = deriv(1,1,e,i)*displ_ele(2) + deriv(2,1,e,i)*displ_ele(4) + deriv(3,1,e,i)*displ_ele(6)
      F(2,2) = 1 + (deriv(1,2,e,i)*displ_ele(2) + deriv(2,2,e,i)*displ_ele(4) + deriv(3,2,e,i)*displ_ele(6))
    else if (pnode == 4) then
      F(1,1) = 1 + (deriv(1,1,e,i)*displ_ele(1) + deriv(2,1,e,i)*displ_ele(3) + deriv(3,1,e,i)*displ_ele(5) +&
                    deriv(4,1,e,i)*displ_ele(7))
      F(1,2) = deriv(1,2,e,i)*displ_ele(1) + deriv(2,2,e,i)*displ_ele(3) + deriv(3,2,e,i)*displ_ele(5) +&
                    deriv(4,2,e,i)*displ_ele(7)
      F(2,1) = deriv(1,1,e,i)*displ_ele(2) + deriv(2,1,e,i)*displ_ele(4) + deriv(3,1,e,i)*displ_ele(6) +&
                    deriv(4,1,e,i)*displ_ele(8)
      F(2,2) = 1 + (deriv(1,2,e,i)*displ_ele(2) + deriv(2,2,e,i)*displ_ele(4) + deriv(3,2,e,i)*displ_ele(6) +&
                    deriv(4,2,e,i)*displ_ele(8))
    end if
    call m22inv(F, inv_F, det_F, inv_flag)
    if (det_F < 1.0E-12) then
      write(*,*) "NEGATIVE GPDET IN ELEMENT",e
      write(*,*) "GPDET:",det_F
      write(*,*) "NODES OF THE ELEMENT:",elements_bulk(e,6),elements_bulk(e,7),elements_bulk(e,8)
      stop
    end if
    if (.not. inv_flag) then
      write(*,*) "Deformation gradient tensor is not invertible... bye!"
      stop
    end if

    !CAUCHY DEFORMATION TENSOR
    do k = 1,2
      do j = 1,2
        Cau(j,k) = 0.0D0
        do l = 1,2
          Cau(j,k) = Cau(j,k) + F(l,j)*F(l,k)
        end do
      end do
    end do
    call m22inv(Cau, inv_Cau, det_Cau, inv_flag)
    if (.not. inv_flag) then
      write(*,*) "Cauchy deformation tensor is not invertible... bye!"
      stop
    end if
    
    !GREEN-LAGRANGE STRESS TENSOR
    Et = 0.0D0
    do j = 1,2
      do k = 1,2
        dummr = 0.0D0
        do l = 1,2
          dummr = dummr + F(l,j)*F(l,k)
        Et(j,k) = 0.5*(dummr - delta_kron(j,k))
        end do
      end do
    end do
    traceE = 0.0D0
    do j = 1,2
      traceE = traceE + Et(j,j)
    end do
    
    !!!!!!START ISOLIN MODEL!!!!!!
    if (model == 'ISOL') then
      !SECOND PIOLA-KIRCHHOFF STRESS TENSOR
      S = 0.0D0
      do j = 1,2
        do k = 1,2
          S(j,k) = lame2_lambda(e)*plane(e)*traceE*delta_kron(j,k) + 2.0*lame1_mu(e)*Et(j,k)
        end do
      end do 
    
      !MATERIAL TANGENT MODULI
      D = 0.0D0
      do j = 1,2
        do k = 1,2
          do l = 1,2
            do m = 1,2
              D(j,k,l,m) = lame2_lambda(e)*plane(e)*delta_kron(j,k)*delta_kron(l,m) + & 
                         lame1_mu(e)*(delta_kron(j,l)*delta_kron(k,m) + delta_kron(j,m)*delta_kron(k,l))
            end do 
          end do
        end do
      end do 
    end if
    !!!!!!END ISOLIN MODEL!!!!!!
    
    !!!!!!START BELYTSCHKO MODEL - NEO-HOOKEAN MATERIAL (AS IMPLEMENTED IN ALYA)!!!!!!
    if (model == 'BELY') then 
      !SECOND PIOLA-KIRCHHOFF STRESS TENSOR
      S = 0.0D0
      do j = 1,2
        do k = 1,2
          S(j,k) = S(j,k) + (lame2_lambda(e)*log(det_F)-lame1_mu(e))*inv_Cau(j,k) + lame1_mu(e)*delta_kron(j,k)
        end do
      end do

      !MATERIAL TANGENT MODULI
      D = 0.0D0
      do j = 1,2
        do k = 1,2
          do l = 1,2
            do m = 1,2
              D(j,k,l,m) = lame2_lambda(e)*inv_Cau(j,k)*inv_Cau(l,m) + & 
              (lame1_mu(e)-lame2_lambda(e)*log(det_F))*(inv_Cau(j,l)*inv_Cau(k,m) + &
              inv_Cau(j,m)*inv_Cau(k,l))
            end do
          end do
        end do
      end do
    end if
    !!!!!!END BELYTSCHKO MODEL!!!!!!
    
    !!!!!!START ZIENKIEWICZ MODEL - NEO-HOOKEAN MATERIAL - ZIENKIEWICZ VOL. 2, PAG. 341!!!!!!
    if (model == 'ZIEN') then 
      !SECOND PIOLA-KIRCHHOFF STRESS TENSOR
      S = 0.0D0
      S(1,1) = lame1_mu(e)*(1-inv_Cau(1,1)) + lame2_lambda(e)*det_F*(det_F-1)*inv_Cau(1,1)
      S(1,2) = lame1_mu(e)*(-inv_Cau(1,2)) + lame2_lambda(e)*det_F*(det_F-1)*inv_Cau(1,2)
      S(2,1) = lame1_mu(e)*(-inv_Cau(2,1)) + lame2_lambda(e)*det_F*(det_F-1)*inv_Cau(2,1)
      S(2,2) = lame1_mu(e)*(1-inv_Cau(2,2)) + lame2_lambda(e)*det_F*(det_F-1)*inv_Cau(2,2)

      !MATERIAL TANGENT MODULI
      D = 0.0D0
      do j = 1,2
        do k = 1,2
          do l = 1,2
            do m = 1,2
              D(j,k,l,m) = lame2_lambda(e)*det_F*(2*det_F-1)*inv_Cau(j,k)*inv_Cau(l,m) + & 
              (lame1_mu(e)-lame2_lambda(e)*det_F*(det_F-1))*(inv_Cau(j,l)*inv_Cau(k,m)+ &
              inv_Cau(j,m)*inv_Cau(k,l))
            end do
          end do
        end do
      end do
    end if
    !!!!!!END ZIENKIEWICZ MODEL!!!!!!

    !!!!!!START LAURSEN MODEL - NEO-HOOKEAN MATERIAL - LAURSEN'S CONTACT BOOK, PAG. 34!!!!!!
    if (model == 'LAUR') then 
      !SECOND PIOLA-KIRCHHOFF STRESS TENSOR
      S = 0.0D0
      do j = 1,2
        do k = 1,2
          S(j,k) = lame1_mu(e)*(delta_kron(j,k) - &
                   (inv_F(j,1)*inv_F(k,1) + inv_F(j,2)*inv_F(k,2))) + &
                   lame2_lambda(e)/2*(det_F*det_F - 1)*(inv_F(j,1)*inv_F(k,1)+inv_F(j,2)*inv_F(k,2))
        end do
      end do

      !MATERIAL TANGENT MODULI
      D = 0.0D0
      do j = 1,2
        do k = 1,2
          do l = 1,2
            do m = 1,2
              D(j,k,l,m) = 2*lame1_mu(e)*(1 + (lame2_lambda(e)/(2*lame1_mu(e)))*(1 - &
              det_F*det_F))*( inv_F(j,1)*inv_F(k,1)*inv_F(l,1)*inv_F(m,1) + &
              inv_F(j,2)*inv_F(k,2)*inv_F(l,2)*inv_F(m,2) + inv_F(j,1)*inv_F(k,2)* &
              inv_F(l,1)*inv_F(m,2) + inv_F(j,2)*inv_F(k,1)*inv_F(l,2)*inv_F(m,1) ) + &
              lame2_lambda(e)*det_F*det_F*( inv_F(j,1)*inv_F(k,1)*inv_F(l,1)*inv_F(m,1) + &
              inv_F(j,2)*inv_F(k,2)*inv_F(l,2)*inv_F(m,2) + &
              inv_F(j,1)*inv_F(k,1)*inv_F(l,2)*inv_F(m,2) + inv_F(j,2)*inv_F(k,2)*inv_F(l,1)*inv_F(m,1) )
            end do
          end do
        end do
      end do
    end if
    !!!!!!END LAURSEN MODEL!!!!!!

    !FIRST PIOLA-KIRCHHOFF STRESS TENSOR
    P = 0.0D0
    do j = 1,2
      do k = 1,2
        do l = 1,2
          P(j,k) = P(j,k) + F(j,l)*S(l,k)
        end do
      end do
    end do

    !TANGENT OF FIRST PIOLA-KIRCHHOFF STRESS TENSOR WITH RESPECT TO DEFORMATION GRADIENT
    dPdF = 0.0D0
    do j = 1,2
      do k = 1,2
        do l = 1,2
          do m = 1,2
            dPdF(j,k,l,m) = dPdF(j,k,l,m) + delta_kron(j,l)*S(k,m)
            do n = 1,2
              do o = 1,2
                dPdF(j,k,l,m) = dPdF(j,k,l,m) + D(n,k,o,m)*F(j,n)*F(l,o) 
              end do
            end do
          end do
        end do
      end do
    end do
   
    if (.not. stress_calc_on) then
      !ELEMENTAL TANGENT STIFFNESS MATRIX
      do anode = 1,pnode
        do idime = 1,2
          adofn = (anode-1)*2+idime
          do bnode = 1,pnode
            do kdime = 1,2
              bdofn = (bnode-1)*2+kdime
              do jdime = 1,2
                do ldime = 1,2
                  k_elem(adofn,bdofn) = k_elem(adofn,bdofn) + &
                  dPdF(idime,jdime,kdime,ldime)*deriv(anode,jdime,e,i)*deriv(bnode,ldime,e,i)*weight(i,pnode-2)*det_jac(e,i)
                end do
              end do
            end do
          end do
        end do
      end do

      !ELEMENTAL RESIDUAL VECTOR
      do inode = 1,pnode
        idofn = (inode-1)*2
        do idime = 1,2
          idofn = idofn + 1
          do jdime = 1,2
            r_elem(idofn) = r_elem(idofn) - P(idime,jdime)*deriv(inode,jdime,e,i)*weight(i,pnode-2)*det_jac(e,i)
          end do
        end do
      end do
    end if

    if (stress_calc_on) then
      !ST (SIGMA FOR EACH GAUSS POINT) - ALYA'S VARIABLE IS sigma
      ST = 0.0D0
      gpinv = 1.0/det_F
      do j = 1,2
        do k = 1,2
          do l = 1,2
            ST(k,j) = ST(k,j) + gpinv*F(k,l)*P(j,l)
          end do
        end do
      end do

      !STRAIN AND STRESS CALCULATION (IN VOIGT NOTATION) - ALYA'S VARIABLE IS caust_sld AND green_sld
      gpvol = det_jac(e,i)*weight(i,pnode-2) 
      call shapefunc(pnode,i,gpsha)
      do inode = 1,pnode
        ipoin = elements_bulk(e,5+inode)
        xmean = gpsha(inode)*gpvol
        do ivoigt = 1,3
          jdime = nvgij(ivoigt,1)
          kdime = nvgij(ivoigt,2)
          stress(ivoigt,ipoin) = stress(ivoigt,ipoin) + xmean*ST(jdime,kdime)
          strain(ivoigt,ipoin) = strain(ivoigt,ipoin) + xmean*Et(jdime,kdime)
        end do
      end do
    end if

  end do !end do gauss points

  if (.not. stress_calc_on) then
    !GLOBAL TANGENT STIFFNESS MATRIX
    do i = 1,pnode
      do j = 1,pnode
        k_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)-1) = & !k_tot(u1,u1)
          k_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)-1) + k_elem((2*i)-1,(2*j)-1)
        k_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)) = & !k_tot(u1,v1)
          k_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)) + k_elem((2*i)-1,(2*j))
        k_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)-1) = & !k_tot(v1,u1)
          k_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)-1) + k_elem((2*i),(2*j)-1)
        k_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)) = & !k_tot(v1,v1)
          k_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)) + k_elem((2*i),(2*j))
      end do
    end do

    !GLOBAL RESIDUAL VECTOR
    do i = 1,pnode
      r_tot((elements_bulk(e,(5+i))*2)-1) = r_tot((elements_bulk(e,(5+i))*2)-1) + r_elem((2*i)-1)
      r_tot(elements_bulk(e,(5+i))*2) = r_tot(elements_bulk(e,(5+i))*2) + r_elem(2*i)
    end do
  end if

end do !end do elements

if (stress_calc_on) then
  do ipoin = 1,num_nodes
    do ivoigt = 1,3
      stress(ivoigt,ipoin) = stress(ivoigt,ipoin)/vmass(ipoin)
      strain(ivoigt,ipoin) = strain(ivoigt,ipoin)/vmass(ipoin)
    end do
  end do
end if

end subroutine assembly

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dealloca_global_matrices()

deallocate(k_tot)
deallocate(r_tot)

end subroutine dealloca_global_matrices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dealloca_stress_strain_matrices()

deallocate(stress)
deallocate(strain)

end subroutine dealloca_stress_strain_matrices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dealloca_init()

deallocate(det_jac)
deallocate(deriv)
deallocate(vmass)
deallocate(weight)
deallocate(pgauss)
deallocate(plane)
deallocate(lame1_mu)
deallocate(lame2_lambda)
deallocate(nodes)
deallocate(elements_bulk)
deallocate(displ)

end subroutine dealloca_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_fortran
