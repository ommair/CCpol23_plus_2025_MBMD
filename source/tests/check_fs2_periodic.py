#!/usr/bin/env python3
"""Check the actual FS2 routine with gfortran and a serial MPI reduction stub.

No production inputs or MPI installation are needed. Compile/run in a temporary
folder; simulation source and existing build artifacts are left untouched.
"""
from pathlib import Path
import subprocess
import tempfile

SOURCE = Path(__file__).resolve().parents[1] / 'threebody_potential.f90'
PREAMBLE = '''module fs2_test
implicit none
integer, parameter :: nom=3,nsite=5
integer :: nm=3,ns=5,MPIrank=0,MPIcomm_size=1
real(8), parameter :: bohr2a=0.529177249d0,h2kcal=627.510d0
real(8) :: alpha3=6.51292d0,rcutsd3=6.1d0**2,r0_tilda=6.00446d0,delta33=10d0
real(8) :: fxfs2(nom),fyfs2(nom),fzfs2(nom),txfs2(nom),tyfs2(nom),tzfs2(nom)
real(8) :: fxs(nsite,nom),fys(nsite,nom),fzs(nsite,nom),vir
contains
'''
DRIVER = '''end module
program check
use fs2_test
implicit none
real(8) :: r(3,3),s(3,5,3),rbase(3,3),f(3,3),fb(3,3),sb(3,5,3)
real(8) :: e,eb,v,vb,ep,em,h,conv,vcom
integer :: i,a,sgn,geometry
h=1d-5
conv=bohr2a/418.4d0
! Three non-collinear geometries cover short and damped separations.
do geometry=1,3
r(1,:)=[0d0,5.5d0,2d0]*geometry**0.5d0
r(2,:)=[0d0,0.3d0,5d0]*geometry**0.5d0
r(3,:)=[0d0,0.2d0,0.7d0]*geometry**0.5d0
do i=1,3
s(1,:,i)=[0d0,1.43d0,-1.43d0,0.3d0,-0.3d0]
s(2,:,i)=[-0.1d0,1d0,1d0,-0.5d0,-0.5d0]
s(3,:,i)=[0d0,0d0,0d0,0.5d0,-0.5d0]
enddo
rbase=r
call calc(r,40d0,e)
eb=e; vb=vir
fb(1,:)=fxfs2+sum(fxs,dim=1)
fb(2,:)=fyfs2+sum(fys,dim=1)
fb(3,:)=fzfs2+sum(fzs,dim=1)
sb(1,:,:)=fxs; sb(2,:,:)=fys; sb(3,:,:)=fzs
! Shift every molecule independently along each axis, both directions.
do i=1,3
do a=1,3
do sgn=-1,1,2
r=rbase; r(a,i)=r(a,i)+sgn*40d0
call calc(r,40d0,e)
call near(e,eb,1d-10,'periodic energy')
call near(vir,vb,1d-7,'periodic virial')
f(1,:)=fxfs2+sum(fxs,dim=1)
f(2,:)=fyfs2+sum(fys,dim=1)
f(3,:)=fzfs2+sum(fzs,dim=1)
call near(maxval(abs(f-fb)),0d0,1d-7,'periodic molecular forces')
call near(maxval(abs(fxs-sb(1,:,:))),0d0,1d-7,'periodic site x forces')
call near(maxval(abs(fys-sb(2,:,:))),0d0,1d-7,'periodic site y forces')
call near(maxval(abs(fzs-sb(3,:,:))),0d0,1d-7,'periodic site z forces')
enddo
enddo
enddo
! Test derivatives in a representation crossing periodic boundaries.
r=rbase; r(1,1)=r(1,1)+40d0; r(2,2)=r(2,2)-40d0
call calc(r,40d0,e)
vcom=sum(fxs*s(1,:,:)+fys*s(2,:,:)+fzs*s(3,:,:))*conv
v=vir/418.4d0+vcom
f=fb*conv
do i=1,3
do a=1,3
r(a,i)=r(a,i)+h
call calc(r,40d0,ep)
r(a,i)=r(a,i)-2*h
call calc(r,40d0,em)
r(a,i)=r(a,i)+h
call near(f(a,i),-(ep-em)/(2*h),1d-7,'finite difference force')
enddo
enddo
call calc(r*(1+h),40d0*(1+h),ep)
call calc(r*(1-h),40d0*(1-h),em)
call near(v,(ep-em)/(2*h),1d-7,'finite difference rigid virial')
enddo
print *, 'PASS: FS2 periodic energy, site/molecular forces, virial, and finite differences'
contains
subroutine calc(rr,b,energy)
real(8), intent(in) :: rr(3,3),b
real(8), intent(out) :: energy
fxs=0; fys=0; fzs=0; vir=0
call FS2_term(rr(1,:),rr(2,:),rr(3,:),s(1,:,:),s(2,:,:),s(3,:,:),b,energy)
end subroutine
subroutine near(actual,expected,tol,label)
real(8) :: actual,expected,tol
character(*) :: label
if (abs(actual-expected)>tol*max(1d0,abs(expected))) then
print *, 'FAIL: ',label,actual,expected
stop 1
endif
end subroutine
end program
'''


def main():
    text = SOURCE.read_text()
    routine = text[text.index('  subroutine FS2_term'):text.index('subroutine FS3_term')]
    reduction = 'call MPI_ALLREDUCE(MPIbuff1, MPIbuff2, j, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPIierr);'
    assert routine.count(reduction) == 1, 'Update serial reduction stub for changed source'
    routine = routine.replace(reduction, 'MPIbuff2 = MPIbuff1')
    with tempfile.TemporaryDirectory(prefix='fs2_periodic_') as folder:
        path = Path(folder)
        (path / 'check.f90').write_text(PREAMBLE + routine + DRIVER)
        subprocess.run(['gfortran', '-O0', '-g', '-fcheck=all',
                        '-ffree-line-length-none', 'check.f90', '-o', 'check'],
                       cwd=path, check=True)
        subprocess.run([str(path / 'check')], cwd=path, check=True)


if __name__ == '__main__':
    main()
