Subroutine Writeout_OMPinfo(Nk)
  Use Omp_lib
  Implicit none
  Integer, intent(in) :: Nk
  Integer :: NOMP

  !$ NOMP = omp_get_max_threads()
  !$ Write(*,'("# ++++++++++++++++++++++++++++")') 
  !$ Write(*,'("# Number of OpenMP thread =",i8)')  NOMP
  !$ If (Nk < NOMP) then
  !$   Write(*,'("# CAUTION: OMP parallelization is not efficient because of Nk < NOMP.")')
  !$ End If
  !$ Write(*,'("# ++++++++++++++++++++++++++++")') 
 
  Return
End Subroutine Writeout_OMPinfo
!==========================================================================================
Subroutine compute_sigma(omega, eigval, occ, pmat, ewidth, Nene, Nk, Nb, sigma)
  Implicit None
  Complex(kind(0d0)), parameter :: zI=(0.0d0, 1.0d0)
  Integer, intent(in) :: Nene, Nk, Nb
  Double precision, intent(in) :: omega(Nene), eigval(Nk,Nb), occ(Nk,Nb), ewidth
  Complex(kind(0d0)), intent(in) :: pmat(Nk,Nb,Nb,3)
  Integer :: Nevery, ik, ib, jb, ixyz, jxyz
  Double precision :: docc
  Complex(kind(0d0)), intent(out) :: sigma(3,3,Nene)
  Complex(kind(0d0)) :: moment, ene_denominator(Nene)

  sigma = 0.d0
  Nevery = Nk/10
  Write(*,'("# Following is progress of GenerateSigmaEpsilon.generate function. ")')
  !$   Write(*,'("# CAUTION: This function does not support OMP parallelization. ")')
  Do ik = 1, Nk
    If (mod(ik, Nevery)==0)Then
      Write(*,*) (ik*100.0)/float(Nk), '% is done.'
    End If
    Do ib = 1, Nb
      Do jb = 1, Nb
        If (ib /= jb) Then
          ene_denominator = omega(:) - (eigval(ik,jb) - eigval(ik,ib)) + zI*ewidth
          ene_denominator = ene_denominator*(eigval(ik,jb) - eigval(ik,ib))
          docc = occ(ik,ib) - occ(ik,jb)
          Do ixyz = 1,3
            Do jxyz = 1,3
              moment = pmat(ik,jb,ib,ixyz)*pmat(ik,ib,jb,jxyz)
              sigma(jxyz,ixyz,:) = sigma(jxyz,ixyz,:) + docc*moment/ene_denominator(:)
            End Do
          End Do
        End If 
      End Do
    End Do
  End Do
  Return
End Subroutine compute_Sigma
!==========================================================================================
Subroutine eigh(N,H,E,V)
  Implicit none
  Integer, intent(in) :: N
  Complex(kind(0d0)), intent(in) :: H(1:N, 1:N)
  Double precision, intent(out) :: E(1:N)
  Complex(kind(0d0)), intent(out) :: V(1:N, 1:N)
!For ZHEEV
  Integer :: LWORK_EV 
  Complex(kind(0d0)) :: WORK_EV(2*(2*N-1)) !The argument is just LWORK_EV
  Double precision :: RWORK_EV(3*N - 2)
  Integer :: INFO_EV
  LWORK_EV = 2*(2*N-1)

  V = H
  Call ZHEEV('V','U',N,V,N,E,WORK_EV,LWORK_EV,RWORK_EV,INFO_EV)
    
  Return
End Subroutine eigh
