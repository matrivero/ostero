!External subroutines for elemental matrix generation and global matrix assembly
!Copyright (C) 2016 Matias Rivero
!
!This program is free software: you can redistribute it and/or modify
!it under the terms of the GNU General Public License as published by
!the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program.  If not, see <http://www.gnu.org/licenses/>.
!matias.rivero@bsc.es

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
real*8, allocatable, dimension(:) :: young
real*8, allocatable, dimension(:) :: poisson
real*8, allocatable, dimension(:) :: density
logical :: stress_calc_on
logical :: gravity_calc_on
real*8 :: grav_magnitude
real*8, allocatable, dimension(:) :: grav_direction
real*8, allocatable, dimension(:) :: damage
real*8, allocatable, dimension(:) :: r_0
real*8, allocatable, dimension(:) :: tau
real*8, allocatable, dimension(:) :: tau_history
real*8, allocatable, dimension(:) :: A_damage

!outputs deriv_and_detjac_calc
real*8, allocatable, dimension(:,:) :: det_jac
real*8, allocatable, dimension(:,:,:,:) :: deriv

!output vmass_calc
real*8, allocatable, dimension(:) :: vmass

!outputs assembly
real*8, allocatable, dimension(:,:) :: k_tot
real*8, allocatable, dimension(:,:) :: m_tot
real*8, allocatable, dimension(:) :: r_tot
real*8, allocatable, dimension(:,:) :: strain
real*8, allocatable, dimension(:,:) :: stress

character(len=4) :: model
character(len=12) :: submodel

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init()

allocate(weight(4,2))
weight = 0.0D0
weight(1,1) = 0.166666666666666D0
weight(2,1) = 0.166666666666666D0
weight(3,1) = 0.166666666666666D0

weight(1,2) = 1.0D0
weight(2,2) = 1.0D0
weight(3,2) = 1.0D0
weight(4,2) = 1.0D0

allocate(pgauss(4,2,2))
pgauss = 0.0D0
pgauss(1,1,1) = 0.166666666666666D0
pgauss(1,2,1) = 0.166666666666666D0
pgauss(2,1,1) = 0.666666666666666D0
pgauss(2,2,1) = 0.166666666666666D0
pgauss(3,1,1) = 0.166666666666666D0
pgauss(3,2,1) = 0.666666666666666D0

pgauss(1,1,2) = -0.577350269189626D0 
pgauss(1,2,2) = -0.577350269189626D0 
pgauss(2,1,2) = 0.577350269189626D0 
pgauss(2,2,2) = -0.577350269189626D0 
pgauss(3,1,2) = 0.577350269189626D0 
pgauss(3,2,2) = 0.577350269189626D0 
pgauss(4,1,2) = -0.577350269189626D0 
pgauss(4,2,2) = 0.577350269189626D0 

allocate(young(num_elements_bulk))
allocate(poisson(num_elements_bulk))
allocate(density(num_elements_bulk))

allocate(damage(num_elements_bulk))
allocate(r_0(num_elements_bulk))
allocate(tau(num_elements_bulk))
allocate(tau_history(num_elements_bulk))
allocate(A_damage(num_elements_bulk))

allocate(nodes(num_nodes,4))
allocate(elements_bulk(num_elements_bulk,9))

allocate(displ(num_nodes*2))

end subroutine init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine derivs(pnode,igauss,dh)

implicit none

integer*4, intent(in) :: pnode,igauss
real*8, dimension(4,2), intent(out) :: dh

dh = 0.0D0
if (pnode == 3) then
  dh(1,1) = -1D0 !dh1/deta
  dh(1,2) = -1D0 !dh1/dxi
  dh(2,1) = 1D0 !dh2/deta
  dh(2,2) = 0D0 !dh2/dxi
  dh(3,1) = 0D0 !dh3/deta
  dh(3,2) = 1D0 !dh3/dxi
else if(pnode == 4) then
  dh(1,1) = -(1D0 - pgauss(igauss,2,2))/4.0D0 !dh1/eta
  dh(1,2) = -(1D0 - pgauss(igauss,1,2))/4.0D0 !dh1/xi
  dh(2,1) = (1D0 - pgauss(igauss,2,2))/4.0D0  !dh2/eta
  dh(2,2) = -(1D0 + pgauss(igauss,1,2))/4.0D0  !dh2/xi
  dh(3,1) = (1D0 + pgauss(igauss,2,2))/4.0D0  !dh3/eta
  dh(3,2) = (1D0 + pgauss(igauss,1,2))/4.0D0  !dh3/xi
  dh(4,1) = -(1D0 + pgauss(igauss,2,2))/4.0D0  !dh4/deta
  dh(4,2) = (1D0 - pgauss(igauss,1,2))/4.0D0 !dh4/dxi
end if

end subroutine derivs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine shapefunc(pnode,igauss,gpsha)

implicit none

integer*4, intent(in) :: pnode,igauss
real*8, dimension(4), intent(out) :: gpsha

gpsha = 0.0D0
if (pnode == 3) then
    gpsha(1) = 1D0 - pgauss(igauss,1,pnode-2) - pgauss(igauss,2,pnode-2) !1-eta-xi 
    gpsha(2) = pgauss(igauss,1,pnode-2) !eta
    gpsha(3) = pgauss(igauss,2,pnode-2) !xi
else if(pnode == 4) then
    gpsha(1) = (1D0 - pgauss(igauss,1,pnode-2))*(1D0 - pgauss(igauss,2,pnode-2))/4.0D0
    gpsha(2) = (1D0 + pgauss(igauss,1,pnode-2))*(1D0 - pgauss(igauss,2,pnode-2))/4.0D0
    gpsha(3) = (1D0 + pgauss(igauss,1,pnode-2))*(1D0 + pgauss(igauss,2,pnode-2))/4.0D0
    gpsha(4) = (1D0 - pgauss(igauss,1,pnode-2))*(1D0 + pgauss(igauss,2,pnode-2))/4.0D0
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
vmass = 0.0D0

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

subroutine m33inv(mat, inv_mat, det, ok_flag)

implicit none

real*8, dimension(3,3), intent(in)  :: mat
real*8, dimension(3,3), intent(out) :: inv_mat
real*8, intent(out)  :: det

real*8, parameter :: eps = 1.0D-10
real*8  :: c11,c12,c13,c21,c22,c23,c31,c32,c33

logical, intent(out) :: ok_flag
! Verificar
 c11 = mat(2,2)*mat(3,3)-mat(3,2)*mat(2,3)
 c12 = mat(1,3)*mat(3,2)-mat(1,2)*mat(3,3)
 c13 = mat(1,2)*mat(2,3)-mat(2,2)*mat(1,3)
 c21 = mat(2,3)*mat(3,1)-mat(3,3)*mat(2,1)
 c22 = mat(1,1)*mat(3,3)-mat(3,1)*mat(1,3)
 c23 = mat(1,3)*mat(2,1)-mat(2,3)*mat(1,1)
 c31 = mat(2,1)*mat(3,2)-mat(3,1)*mat(2,2)
 c32 = mat(1,2)*mat(3,1)-mat(3,2)*mat(1,1)
 c33 = mat(1,1)*mat(2,2)-mat(2,1)*mat(1,2)
 
 det=mat(1,1)*c11+mat(1,2)*c21+mat(1,3)*c31
 
 inv_mat(1,1)=c11/det
 inv_mat(2,1)=c21/det
 inv_mat(3,1)=c31/det! Verificar
 inv_mat(1,2)=c12/det
 inv_mat(2,2)=c22/det
 inv_mat(3,2)=c32/det
 inv_mat(1,3)=c13/det
 inv_mat(2,3)=c23/det
 inv_mat(3,3)=c33/det

if (abs(det) .le. eps) then
  inv_mat = 0.0D0
  ok_flag = .false.
  return
end if

ok_flag = .true.

return

end subroutine m33inv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assembly_nonlinear()

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
real*8, dimension(8,8) :: m_elem
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
real*8 :: lame2_lambda,lame1_mu,plane

if (.not. stress_calc_on) then
  allocate(k_tot(num_nodes*2,num_nodes*2))
  allocate(m_tot(num_nodes*2,num_nodes*2))
  allocate(r_tot(num_nodes*2))
  k_tot = 0.0D0
  m_tot = 0.0D0
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

do e = 1,num_elements_bulk !ELEMENTS LOOP
  
  lame2_lambda = poisson(e)*young(e)/((1+poisson(e))*(1-2*poisson(e))) 
  lame1_mu = young(e)/(2*(1+poisson(e)))
  if ((model == 'ISOL') .and. (submodel == 'PLANE_STRESS')) then
    plane = (1.0-2.0*poisson(e))/(1.0-poisson(e))
  else if ((model == 'ISOL') .and. (submodel == 'PLANE_STRAIN')) then
    plane = 1.0D0
  else if (model == 'ISOL') then
    print *, "You need a submodel for ISOLIN model: PLANE_STRESS or PLANE_STRAIN... bye!"
    stop
  end if
  
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
  m_elem = 0.0D0
  r_elem = 0.0D0

  do i = 1,ngauss !GAUSS POINT LOOP
    !DEFORMATION GRADIENT TENSOR
    F = 0.0D0
    if (pnode == 3) then
      F(1,1) = 1D0 + (deriv(1,1,e,i)*displ_ele(1) + deriv(2,1,e,i)*displ_ele(3) + deriv(3,1,e,i)*displ_ele(5))
      F(1,2) = deriv(1,2,e,i)*displ_ele(1) + deriv(2,2,e,i)*displ_ele(3) + deriv(3,2,e,i)*displ_ele(5)
      F(2,1) = deriv(1,1,e,i)*displ_ele(2) + deriv(2,1,e,i)*displ_ele(4) + deriv(3,1,e,i)*displ_ele(6)
      F(2,2) = 1D0 + (deriv(1,2,e,i)*displ_ele(2) + deriv(2,2,e,i)*displ_ele(4) + deriv(3,2,e,i)*displ_ele(6))
    else if (pnode == 4) then
      F(1,1) = 1D0 + (deriv(1,1,e,i)*displ_ele(1) + deriv(2,1,e,i)*displ_ele(3) + deriv(3,1,e,i)*displ_ele(5) +&
                    deriv(4,1,e,i)*displ_ele(7))
      F(1,2) = deriv(1,2,e,i)*displ_ele(1) + deriv(2,2,e,i)*displ_ele(3) + deriv(3,2,e,i)*displ_ele(5) +&
                    deriv(4,2,e,i)*displ_ele(7)
      F(2,1) = deriv(1,1,e,i)*displ_ele(2) + deriv(2,1,e,i)*displ_ele(4) + deriv(3,1,e,i)*displ_ele(6) +&
                    deriv(4,1,e,i)*displ_ele(8)
      F(2,2) = 1D0 + (deriv(1,2,e,i)*displ_ele(2) + deriv(2,2,e,i)*displ_ele(4) + deriv(3,2,e,i)*displ_ele(6) +&
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
    if (model == 'ISOL')  then
      !SECOND PIOLA-KIRCHHOFF STRESS TENSOR
      S = 0.0D0
      do j = 1,2
        do k = 1,2
          S(j,k) = lame2_lambda*plane*traceE*delta_kron(j,k) + 2.0*lame1_mu*Et(j,k)
        end do
      end do 
    
      !MATERIAL TANGENT MODULI
      D = 0.0D0
      do j = 1,2
        do k = 1,2
          do l = 1,2
            do m = 1,2
              D(j,k,l,m) = lame2_lambda*plane*delta_kron(j,k)*delta_kron(l,m) + & 
                         lame1_mu*(delta_kron(j,l)*delta_kron(k,m) + delta_kron(j,m)*delta_kron(k,l))
            end do 
          end do
        end do
      end do 
    !!!!!!END ISOLIN MODEL!!!!!!
    
    !!!!!!START BELYTSCHKO MODEL - NEO-HOOKEAN MATERIAL (AS IMPLEMENTED IN ALYA)!!!!!!
    else if (model == 'BELY') then 
      !SECOND PIOLA-KIRCHHOFF STRESS TENSOR
      S = 0.0D0
      do j = 1,2
        do k = 1,2
          S(j,k) = S(j,k) + (lame2_lambda*log(det_F)-lame1_mu)*inv_Cau(j,k) + lame1_mu*delta_kron(j,k)
        end do
      end do

      !MATERIAL TANGENT MODULI
      D = 0.0D0
      do j = 1,2
        do k = 1,2
          do l = 1,2
            do m = 1,2
              D(j,k,l,m) = lame2_lambda*inv_Cau(j,k)*inv_Cau(l,m) + & 
              (lame1_mu-lame2_lambda*log(det_F))*(inv_Cau(j,l)*inv_Cau(k,m) + &
              inv_Cau(j,m)*inv_Cau(k,l))
            end do
          end do
        end do
      end do
    !!!!!!END BELYTSCHKO MODEL!!!!!!
    
    !!!!!!START ZIENKIEWICZ MODEL - NEO-HOOKEAN MATERIAL - ZIENKIEWICZ VOL. 2, PAG. 341!!!!!!
    else if (model == 'ZIEN') then 
      !SECOND PIOLA-KIRCHHOFF STRESS TENSOR
      S = 0.0D0
      S(1,1) = lame1_mu*(1-inv_Cau(1,1)) + lame2_lambda*det_F*(det_F-1)*inv_Cau(1,1)
      S(1,2) = lame1_mu*(-inv_Cau(1,2)) + lame2_lambda*det_F*(det_F-1)*inv_Cau(1,2)
      S(2,1) = lame1_mu*(-inv_Cau(2,1)) + lame2_lambda*det_F*(det_F-1)*inv_Cau(2,1)
      S(2,2) = lame1_mu*(1-inv_Cau(2,2)) + lame2_lambda*det_F*(det_F-1)*inv_Cau(2,2)

      !MATERIAL TANGENT MODULI
      D = 0.0D0
      do j = 1,2
        do k = 1,2
          do l = 1,2
            do m = 1,2
              D(j,k,l,m) = lame2_lambda*det_F*(2*det_F-1)*inv_Cau(j,k)*inv_Cau(l,m) + & 
              (lame1_mu-lame2_lambda*det_F*(det_F-1))*(inv_Cau(j,l)*inv_Cau(k,m)+ &
              inv_Cau(j,m)*inv_Cau(k,l))
            end do
          end do
        end do
      end do
    !!!!!!END ZIENKIEWICZ MODEL!!!!!!

    !!!!!!START LAURSEN MODEL - NEO-HOOKEAN MATERIAL - LAURSEN'S CONTACT BOOK, PAG. 34!!!!!!
    else if (model == 'LAUR') then 
      !SECOND PIOLA-KIRCHHOFF STRESS TENSOR
      S = 0.0D0
      do j = 1,2
        do k = 1,2
          S(j,k) = lame1_mu*(delta_kron(j,k) - &
                   (inv_F(j,1)*inv_F(k,1) + inv_F(j,2)*inv_F(k,2))) + &
                   lame2_lambda/2*(det_F*det_F - 1)*(inv_F(j,1)*inv_F(k,1)+inv_F(j,2)*inv_F(k,2))
        end do
      end do

      !MATERIAL TANGENT MODULI
      D = 0.0D0
      do j = 1,2
        do k = 1,2
          do l = 1,2
            do m = 1,2
              D(j,k,l,m) = 2*lame1_mu*(1 + (lame2_lambda/(2*lame1_mu))*(1 - &
              det_F*det_F))*( inv_F(j,1)*inv_F(k,1)*inv_F(l,1)*inv_F(m,1) + &
              inv_F(j,2)*inv_F(k,2)*inv_F(l,2)*inv_F(m,2) + inv_F(j,1)*inv_F(k,2)* &
              inv_F(l,1)*inv_F(m,2) + inv_F(j,2)*inv_F(k,1)*inv_F(l,2)*inv_F(m,1) ) + &
              lame2_lambda*det_F*det_F*( inv_F(j,1)*inv_F(k,1)*inv_F(l,1)*inv_F(m,1) + &
              inv_F(j,2)*inv_F(k,2)*inv_F(l,2)*inv_F(m,2) + &
              inv_F(j,1)*inv_F(k,1)*inv_F(l,2)*inv_F(m,2) + inv_F(j,2)*inv_F(k,2)*inv_F(l,1)*inv_F(m,1) )
            end do
          end do
        end do
      end do
    !!!!!!END LAURSEN MODEL!!!!!!

    else
      print*, "Check model name... bye!"
      stop
    end if

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
    !SEE SECTION 5.4.8 ELASTICITY TENSORS - BELYTSCHKO's BOOK
    dPdF = 0.0D0
    do j = 1,2
      do k = 1,2
        do l = 1,2
          do m = 1,2
            dPdF(j,k,l,m) = dPdF(j,k,l,m) + delta_kron(j,l)*S(k,m) !GEOMETRIC STIFFNESS
            do n = 1,2
              do o = 1,2
                dPdF(j,k,l,m) = dPdF(j,k,l,m) + D(n,k,o,m)*F(j,n)*F(l,o) !MATERIAL TANGENT STIFFNESS
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

      !ELEMENTAL MASS MATRIX
      gpvol = det_jac(e,i)*weight(i,pnode-2) 
      call shapefunc(pnode,i,gpsha)
      do idime = 1,2
        do anode = 1,pnode
          adofn = (anode-1)*2+idime
          do bnode = 1,pnode
              bdofn = (bnode-1)*2+idime
              m_elem(adofn,bdofn) = m_elem(adofn,bdofn) + density(e)*gpsha(anode)*gpsha(bnode)*gpvol
          end do
        end do
      end do

      !ELEMENTAL RESIDUAL VECTOR
      gpvol = det_jac(e,i)*weight(i,pnode-2) 
      call shapefunc(pnode,i,gpsha)
      do inode = 1,pnode
        idofn = (inode-1)*2
        do idime = 1,2
          idofn = idofn + 1
          do jdime = 1,2
            r_elem(idofn) = r_elem(idofn) - P(idime,jdime)*deriv(inode,jdime,e,i)*weight(i,pnode-2)*det_jac(e,i)
          end do
          if (gravity_calc_on) then
            r_elem(idofn) = r_elem(idofn) + (gpsha(inode)*density(e)*gpvol)*(grav_magnitude*grav_direction(idime))
          end if
        end do
      end do
    end if

    if (stress_calc_on) then
      !ST (SIGMA FOR EACH GAUSS POINT) - ALYA'S VARIABLE IS sigma
      ST = 0.0D0
      gpinv = 1.0D0/det_F
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

        m_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)-1) = & !m_tot(u1,u1)
          m_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)-1) + m_elem((2*i)-1,(2*j)-1)
        m_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)) = & !m_tot(u1,v1)
          m_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)) + m_elem((2*i)-1,(2*j))
        m_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)-1) = & !m_tot(v1,u1)
          m_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)-1) + m_elem((2*i),(2*j)-1)
        m_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)) = & !m_tot(v1,v1)
          m_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)) + m_elem((2*i),(2*j))
      end do

      !GLOBAL RESIDUAL VECTOR
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

end subroutine assembly_nonlinear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assembly_linear()

implicit none

real*8, dimension(8) :: displ_ele
real*8, dimension(8,8) :: k_elem
real*8, dimension(8,8) :: m_elem
real*8, dimension(8) :: r_elem
real*8, dimension(3,8) :: B
real*8, dimension(8,3) :: Bt
real*8, dimension(3,3) :: C
real*8, dimension(3,8) :: TMP
real*8, dimension(4) :: gpsha
real*8, dimension(3) :: strain_elem
real*8 :: gpvol, xmean
integer*4 :: e, i, j, k, l, m
integer*4 :: anode, bnode, adofn, bdofn
integer*4 :: pnode, ngauss, inode, ivoigt, ipoin
integer*4 :: idofn, idime

if (.not. stress_calc_on) then
  allocate(k_tot(num_nodes*2,num_nodes*2))
  allocate(m_tot(num_nodes*2,num_nodes*2))
  allocate(r_tot(num_nodes*2))
  k_tot = 0.0D0
  m_tot = 0.0D0
  r_tot = 0.0D0
end if

if (stress_calc_on) then
  allocate(stress(3,num_nodes))
  allocate(strain(3,num_nodes))
  stress = 0.0D0
  strain = 0.0D0
end if

displ_ele = 0.0D0

do e = 1,num_elements_bulk !ELEMENTS LOOP

  if (elements_bulk(e,2) == 2) then !TRIANGLE ELEMENT
    displ_ele(1) = displ(elements_bulk(e,6)*2-1)
    displ_ele(2) = displ(elements_bulk(e,6)*2)
    displ_ele(3) = displ(elements_bulk(e,7)*2-1)
    displ_ele(4) = displ(elements_bulk(e,7)*2)
    displ_ele(5) = displ(elements_bulk(e,8)*2-1)
    displ_ele(6) = displ(elements_bulk(e,8)*2)
    pnode = 3
    ngauss = 3
  else if (elements_bulk(e,2) == 3) then !QUAD ELEMENT
    displ_ele(1) = displ(elements_bulk(e,6)*2-1)
    displ_ele(2) = displ(elements_bulk(e,6)*2)
    displ_ele(3) = displ(elements_bulk(e,7)*2-1)
    displ_ele(4) = displ(elements_bulk(e,7)*2)
    displ_ele(5) = displ(elements_bulk(e,8)*2-1)
    displ_ele(6) = displ(elements_bulk(e,8)*2)
    displ_ele(7) = displ(elements_bulk(e,9)*2-1)
    displ_ele(8) = displ(elements_bulk(e,9)*2)
    pnode = 4
    ngauss = 4
  end if
  strain_elem = 0.0D0

  if (submodel == 'PLANE_STRESS') then 

    !PLANE STRESS, ISOTROPIC MATERIAL (SIGZZ = 0.0)
    !ZIENKIEWICZ, VOL. 1, PAG. 90
    C(1,1) = (young(e)/(1.0D0-poisson(e)**2))*1.0D0 
    C(1,2) = (young(e)/(1.0D0-poisson(e)**2))*poisson(e)
    C(1,3) = (young(e)/(1.0D0-poisson(e)**2))*0.0D0
    C(2,1) = (young(e)/(1.0D0-poisson(e)**2))*poisson(e)
    C(2,2) = (young(e)/(1.0D0-poisson(e)**2))*1.0D0 
    C(2,3) = (young(e)/(1.0D0-poisson(e)**2))*0.0D0
    C(3,1) = (young(e)/(1.0D0-poisson(e)**2))*0.0D0
    C(3,2) = (young(e)/(1.0D0-poisson(e)**2))*0.0D0
    C(3,3) = (young(e)/(1.0D0-poisson(e)**2))*(1-poisson(e))/2.0D0

  else if(submodel == 'PLANE_STRAIN') then  

    !PLANE STRAIN, ISOTROPIC MATERIAL (SIGZZ != 0.0)
    !ZIENKIEWICZ, VOL. 1, PAG. 91
    C(1,1) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*(1-poisson(e))
    C(1,2) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*poisson(e)
    C(1,3) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*0.0D0
    C(2,1) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*poisson(e)
    C(2,2) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*(1-poisson(e))
    C(2,3) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*0.0D0 
    C(3,1) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*0.0D0
    C(3,2) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*0.0D0
    C(3,3) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*(1-2.0D0*poisson(e))/2.0D0

  else

    write(*,*) "In LINEAR treatment you need to specify a submodel: PLANE_STRESS or PLANE_STRAIN... bye!"
    stop

  end if

  C = (1.0-damage(e))*C !DAMAGE CONSTANT - ONLY USED WHEN DAMAGE MODEL IS ON
  
  k_elem = 0.0D0
  m_elem = 0.0D0
  r_elem = 0.0D0

  do i = 1,ngauss !GAUSS POINT LOOP

    do inode = 1,pnode !B assembly

      B(1,2*inode-1) = deriv(inode,1,e,i) !dh_inode/dx
      B(1,2*inode)   = 0.0D0 
      B(2,2*inode-1) = 0.0D0
      B(2,2*inode)   = deriv(inode,2,e,i) !dh_inode/dy
      B(3,2*inode-1) = deriv(inode,2,e,i) !dh_inode/dy
      B(3,2*inode)   = deriv(inode,1,e,i) !dh_inode/dx

    end do !end do pnode

    do k = 1,3
      do m = 1,(2*pnode)
        Bt(m,k) = B(k,m)
      end do
    end do

    TMP = 0.0D0
    do k = 1,3
      do m = 1,(2*pnode)
        do l = 1,3
          TMP(k,m) = TMP(k,m) + C(k,l)*B(l,m)
        end do
      end do
    end do

    if (.not. stress_calc_on) then
      !ELEMENTAL STIFFNESS MATRIX
      do k = 1,(2*pnode)
        do m = 1,(2*pnode)
          do l = 1,3
            k_elem(k,m) = k_elem(k,m) + Bt(k,l)*TMP(l,m)*det_jac(e,i)*weight(i,pnode-2)
          end do
        end do
      end do

      !ELEMENTAL MASS MATRIX
      gpvol = det_jac(e,i)*weight(i,pnode-2) 
      call shapefunc(pnode,i,gpsha)
      do idime = 1,2
        do anode = 1,pnode
          adofn = (anode-1)*2+idime
          do bnode = 1,pnode
              bdofn = (bnode-1)*2+idime
              m_elem(adofn,bdofn) = m_elem(adofn,bdofn) + density(e)*gpsha(anode)*gpsha(bnode)*gpvol
          end do
        end do
      end do

      !ELEMENTAL RESIDUAL VECTOR (GRAVITY CONTRIBUTION)
      if (gravity_calc_on) then
        gpvol = det_jac(e,i)*weight(i,pnode-2) 
        call shapefunc(pnode,i,gpsha)
        do inode = 1,pnode
          idofn = (inode-1)*2
          do idime = 1,2
            idofn = idofn + 1
            r_elem(idofn) = r_elem(idofn) + (gpsha(inode)*density(e)*gpvol)*(grav_magnitude*grav_direction(idime))
          end do
        end do
      end if
    end if

    if (stress_calc_on) then
      gpvol = det_jac(e,i)*weight(i,pnode-2) 
      call shapefunc(pnode,i,gpsha)
      do ivoigt = 1,3
        do inode = 1,pnode
          xmean = gpsha(inode)*gpvol
          strain_elem(ivoigt) = strain_elem(ivoigt) + &
              xmean*B(ivoigt,2*inode-1)*displ_ele(2*inode-1) + xmean*B(ivoigt,2*inode)*displ_ele(2*inode) 
        end do
      end do
    end if

  end do !end do ngauss
 
  !RESIDUAL CALCULATION (R - K*d)
  do i = 1,(2*pnode)
      do m = 1,(2*pnode)
        r_elem(i) = r_elem(i) - k_elem(i,m)*displ_ele(m)
      end do
  end do
  
  if (.not. stress_calc_on) then
    do i = 1,pnode
      !GLOBAL STIFFNESS MATRIX
      do j = 1,pnode
        k_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)-1) = & !k_tot(u1,u1)
          k_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)-1) + k_elem((2*i)-1,(2*j)-1)
        k_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)) = & !k_tot(u1,v1)
          k_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)) + k_elem((2*i)-1,(2*j))
        k_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)-1) = & !k_tot(v1,u1)
          k_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)-1) + k_elem((2*i),(2*j)-1)
        k_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)) = & !k_tot(v1,v1)
          k_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)) + k_elem((2*i),(2*j))

        m_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)-1) = & !m_tot(u1,u1)
          m_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)-1) + m_elem((2*i)-1,(2*j)-1)
        m_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)) = & !m_tot(u1,v1)
          m_tot((elements_bulk(e,(5+i))*2)-1,(elements_bulk(e,(5+j))*2)) + m_elem((2*i)-1,(2*j))
        m_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)-1) = & !m_tot(v1,u1)
          m_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)-1) + m_elem((2*i),(2*j)-1)
        m_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)) = & !m_tot(v1,v1)
          m_tot((elements_bulk(e,(5+i))*2),(elements_bulk(e,(5+j))*2)) + m_elem((2*i),(2*j))
      end do
      
      !GLOBAL RESIDUAL VECTOR
      r_tot((elements_bulk(e,(5+i))*2)-1) = r_tot((elements_bulk(e,(5+i))*2)-1) + r_elem((2*i)-1)
      r_tot(elements_bulk(e,(5+i))*2) = r_tot(elements_bulk(e,(5+i))*2) + r_elem(2*i)
    end do
  end if

  if (stress_calc_on) then
    do inode = 1,pnode
      ipoin = elements_bulk(e,5+inode)
      do ivoigt = 1,3
        strain(ivoigt,ipoin) = strain(ivoigt,ipoin) + strain_elem(ivoigt)
        do k = 1,3
          stress(ivoigt,ipoin) = stress(ivoigt,ipoin) + C(ivoigt,k)*strain_elem(k)
        end do
      end do
    end do
  end if

end do !end do elements

if (stress_calc_on) then
  do ipoin = 1,num_nodes
    do ivoigt = 1,2
      strain(ivoigt,ipoin) = strain(ivoigt,ipoin)/vmass(ipoin)
      stress(ivoigt,ipoin) = stress(ivoigt,ipoin)/vmass(ipoin)
    end do
    strain(3,ipoin) = strain(3,ipoin)/vmass(ipoin)/2
    stress(3,ipoin) = stress(3,ipoin)/vmass(ipoin)
  end do
end if

end subroutine assembly_linear

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calculate_damage()

implicit none

real*8, dimension(3)   :: sig_0
real*8, dimension(3,3) :: C
real*8, dimension(3,3) :: inv_C
real*8, dimension(8)   :: displ_ele
real*8, dimension(3)   :: strain_ele
real*8, dimension(3,8) :: B

real*8    :: r
real*8    :: aux
real*8    :: vol
real*8    :: det_C

integer*4 :: e
integer*4 :: m
integer*4 :: k
integer*4 :: i
integer*4 :: j
integer*4 :: inode
integer*4 :: pnode
integer*4 :: ngauss

logical   :: inv_flag

do e = 1,num_elements_bulk !ELEMENTS LOOP

  if (elements_bulk(e,2) == 2) then !TRIANGLE ELEMENT
    displ_ele(1) = displ(elements_bulk(e,6)*2-1)
    displ_ele(2) = displ(elements_bulk(e,6)*2)
    displ_ele(3) = displ(elements_bulk(e,7)*2-1)
    displ_ele(4) = displ(elements_bulk(e,7)*2)
    displ_ele(5) = displ(elements_bulk(e,8)*2-1)
    displ_ele(6) = displ(elements_bulk(e,8)*2)
    pnode = 3
    ngauss = 3
  else if (elements_bulk(e,2) == 3) then !QUAD ELEMENT
    displ_ele(1) = displ(elements_bulk(e,6)*2-1)
    displ_ele(2) = displ(elements_bulk(e,6)*2)
    displ_ele(3) = displ(elements_bulk(e,7)*2-1)
    displ_ele(4) = displ(elements_bulk(e,7)*2)
    displ_ele(5) = displ(elements_bulk(e,8)*2-1)
    displ_ele(6) = displ(elements_bulk(e,8)*2)
    displ_ele(7) = displ(elements_bulk(e,9)*2-1)
    displ_ele(8) = displ(elements_bulk(e,9)*2)
    pnode = 4
    ngauss = 4
  end if
  
  if (submodel == 'PLANE_STRESS') then 

    !PLANE STRESS, ISOTROPIC MATERIAL (SIGZZ = 0.0)
    !ZIENKIEWICZ, VOL. 1, PAG. 90
    C(1,1) = (young(e)/(1.0D0-poisson(e)**2))*1.0D0 
    C(1,2) = (young(e)/(1.0D0-poisson(e)**2))*poisson(e)
    C(1,3) = (young(e)/(1.0D0-poisson(e)**2))*0.0D0
    C(2,1) = (young(e)/(1.0D0-poisson(e)**2))*poisson(e)
    C(2,2) = (young(e)/(1.0D0-poisson(e)**2))*1.0D0 
    C(2,3) = (young(e)/(1.0D0-poisson(e)**2))*0.0D0
    C(3,1) = (young(e)/(1.0D0-poisson(e)**2))*0.0D0
    C(3,2) = (young(e)/(1.0D0-poisson(e)**2))*0.0D0
    C(3,3) = (young(e)/(1.0D0-poisson(e)**2))*(1-poisson(e))/2.0D0

  else if(submodel == 'PLANE_STRAIN') then  

    !PLANE STRAIN, ISOTROPIC MATERIAL (SIGZZ != 0.0)integer*4 :: k
    !ZIENKIEWICZ, VOL. 1, PAG. 91
    C(1,1) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*(1-poisson(e))
    C(1,2) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*poisson(e)
    C(1,3) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*0.0D0
    C(2,1) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*poisson(e)
    C(2,2) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*(1-poisson(e))
    C(2,3) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*0.0D0 
    C(3,1) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*0.0D0
    C(3,2) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*0.0D0
    C(3,3) = (young(e)/((1.0D0+poisson(e))*(1-2.0D0*poisson(e))))*(1-2.0D0*poisson(e))/2.0D0

  else

    write(*,*) "In LINEAR treatment you need to specify a submodel: PLANE_STRESS or PLANE_STRAIN... bye!"
    stop

  end if
  call m33inv(C, inv_C, det_C, inv_flag)
  if (.not. inv_flag) then
      write(*,*) "Constitutive Matrix is not invertible... bye!"
      stop
  end if
  
  sig_0      = 0
  strain_ele = 0
  vol        = 0

  do i = 1,ngauss !GAUSS POINT LOOP
     do inode = 1,pnode !B assembly
     
        B(1,2*inode-1) = deriv(inode,1,e,i) !dh_inode/dx
        B(1,2*inode)   = 0.0D0 
        B(2,2*inode-1) = 0.0D0
        B(2,2*inode)   = deriv(inode,2,e,i) !dh_inode/dy
        B(3,2*inode-1) = deriv(inode,2,e,i) !dh_inode/dy
        B(3,2*inode)   = deriv(inode,1,e,i) !dh_inode/dx
     
     end do !end do pnode
          
     !esta es una manera gay de calcular la sigma (primero epsilon y despues sigma = C*epsilon) 
     do k = 1,3
        do m = 1,2*pnode
           strain_ele(k) = strain_ele(k) + B(k,m)*displ_ele(m)
        end do
     end do
     do k = 1,3
        do m = 1,2*pnode
           sig_0(k) = sig_0(k) + C(k,m)*strain_ele(m)
        end do
     end do
     sig_0 = sig_0 + sig_0 * det_jac(e,i)*weight(i,pnode-2)
     vol = vol + det_jac(e,i)*weight(i,pnode-2)
  end do
  
  sig_0 = sig_0 / vol
  
  !tau = sqrt(Sig_0 * C_0^-1 * Sig_0)
  aux = 0.0
  do i=1,3
     do j=1,3
        aux = aux + Sig_0(i) * inv_C(i,j) * Sig_0(j)
     end do   
  end do  
  tau(e) = sqrt(aux)

  r = max(r_0(e),tau_history(e),tau(e))
  damage(e) =  1-r_0(e)/r * exp(A_damage(e)*(1-r/r_0(e)));

end do

end subroutine calculate_damage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dealloca_global_matrices()

if (allocated(k_tot)) deallocate(k_tot)
if (allocated(m_tot)) deallocate(m_tot)
if (allocated(r_tot)) deallocate(r_tot)

end subroutine dealloca_global_matrices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dealloca_stress_strain_matrices()

if (allocated(stress)) deallocate(stress)
if (allocated(strain)) deallocate(strain)

end subroutine dealloca_stress_strain_matrices

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dealloca_init()

deallocate(det_jac)
deallocate(deriv)
deallocate(vmass)
deallocate(weight)
deallocate(pgauss)
deallocate(young)
deallocate(poisson)
deallocate(density)
deallocate(nodes)
deallocate(elements_bulk)
deallocate(displ)
deallocate(damage)
deallocate(r_0)
deallocate(tau)
deallocate(tau_history)
deallocate(A_damage)

end subroutine dealloca_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mod_fortran
