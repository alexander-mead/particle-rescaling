MODULE cosdef

  TYPE cosmology
     REAL :: om_m, om_b, om_v, h, n, sig8, w, A, gamma, dc, fr0, frn
     INTEGER :: ipow, igro, idc, itk, icos, igrok, irate, iratek, icham
     INTEGER :: np, ng
     REAL, ALLOCATABLE :: kp_tab(:), p_tab(:)
     REAL, ALLOCATABLE :: ag_tab(:), g_tab(:), zg_tab(:)
     REAL, ALLOCATABLE :: afg_tab(:), fg_tab(:), zfg_tab(:)
     REAL, ALLOCATABLE :: m(:), dcol(:)
     !     REAL, ALLOCATABLE :: k_g(:), gk(:), k_r(:), rk(:)
  END TYPE cosmology

  TYPE parameters
     CHARACTER(len=256), ALLOCATABLE :: input_file(:), output_file(:)
     CHARACTER(len=256) :: disp_file
     CHARACTER(len=256) :: pk_initial, growth_initial, rate_initial
     CHARACTER(len=256) :: pk_target, growth_target, rate_target
     REAL :: L_target, sig8_initial, sig8_target, omm_target, omv_target
     REAL, ALLOCATABLE :: z_target(:), knl_target(:)
     INTEGER :: izsd, idisp, ifudge, n_files
  END TYPE parameters

END MODULE cosdef

PROGRAM modcat

  USE, INTRINSIC :: iso_c_binding
  USE cosdef

  IMPLICIT NONE

  !INCLUDE '/usr/include/fftw3.f' !For Linux
  INCLUDE '/usr/local/include/fftw3.f' !For Mac

  REAL :: s, mlim, ms, kmin, kmax, kny, kbox, rsm
  REAL :: ai, at, zi, zt, Li, Lt, mi, mt
  CHARACTER(len=256) :: input_file, output_file
  CHARACTER(len=256) :: pfile
  INTEGER :: n, i, j, k, mesh, minmesh, isim
  INTEGER :: ix, iy, iz
  REAL :: dx, dy, dz, dvx, dvy, dvz
  REAL, ALLOCATABLE :: x(:,:), v(:,:)
  REAL, ALLOCATABLE :: weight(:)
  INTEGER, ALLOCATABLE :: id(:)
  INTEGER*8 :: plan_d, plan_x, plan_y, plan_z, plan_x1, plan_y1, plan_z1
  REAL :: kx, ky, kz, kmod, knl, var_psi, var_psi_theory, var_d_theory, sig8i, sig8f
  REAL :: fudge, vfac
  REAL :: avd, avv
  TYPE(cosmology) :: cosi, cost
  TYPE(parameters) :: params
  REAL, PARAMETER :: pi=3.141592654
  REAL :: powi, powt, ratet, ratei, modfac, modfacv
  DOUBLE COMPLEX, ALLOCATABLE :: df_x(:,:,:), df_y(:,:,:), df_z(:,:,:), d(:,:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: dv_x(:,:,:), dv_y(:,:,:), dv_z(:,:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: psi_x(:,:,:), psi_y(:,:,:), psi_z(:,:,:)
  DOUBLE COMPLEX, ALLOCATABLE :: fvx(:,:,:), fvy(:,:,:), fvz(:,:,:)
  LOGICAL :: lexist

  CALL get_command_argument(1,pfile)
  IF(pfile=='') STOP 'Specify input file'
  INQUIRE(file=pfile,exist=lexist)
  IF(lexist .EQV. .FALSE.) STOP 'This parameters file does not exist'

  !  STOP

  !  CALL get_command_argument(1,infile)
  !  IF(infile=='') STOP 'Specify input file'
  !  INQUIRE(file=infile,exist=lexist)
  !  IF(lexist .EQV. .FALSE.) STOP 'This input file does not exist'

  !  CALL get_command_argument(2,outfile)
  !  IF(outfile=='') STOP 'Specify output file'
  !  INQUIRE(file=outfile,exist=lexist)
  !  IF(lexist .EQV. .TRUE.) STOP 'This output file exists, please delete it'

  WRITE(*,*) 
  WRITE(*,*) 'Hello, you are now in MODPART the particle position modifier'
  WRITE(*,*) '============================================================'
  WRITE(*,*)

  CALL read_params(params,pfile)

!  zt=params%z_target
  Lt=params%L_target

!  WRITE(*,*) 'Input file:', params%input_file
!  WRITE(*,*) 'Output file:', params%output_file
!  WRITE(*,*) 'Target redshift:', zt
  WRITE(*,*) 'Target box size:', Lt
  WRITE(*,*)

  IF(params%izsd==0) WRITE(*,*) 'Size and redshift scaling only (zs)'
  IF(params%izsd==1) WRITE(*,*) 'Full method including displacement fields (zsd)'
  WRITE(*,*)

  cosi%sig8=params%sig8_initial

  cost%sig8=params%sig8_target
  cost%om_m=params%omm_target
  cost%om_v=params%omv_target

  !  WRITE(*,*) 'Initial cosmology'
  !  CALL assign_cosmology(cosi,zi)
  !  ai=1./(1.+zi)

  IF(params%izsd==1) THEN
     CALL input_camb_pk(params%pk_initial,cosi)
     CALL input_growth(params%growth_initial,cosi)
     CALL normalise(cosi)
  END IF

  CALL input_rate(params%rate_initial,cosi)

  DO isim=1,params%n_files

     input_file=params%input_file(isim)
     output_file=params%output_file(isim)
     zt=params%z_target(isim)
     IF(params%izsd==1) knl=params%knl_target(isim)
     
     CALL read_gadget(x,v,id,Li,cosi%om_m,cosi%om_v,mi,ai,zi,n,input_file)
     CALL replace(x,Li)

     IF(params%izsd==1) THEN
        sig8i=sigma(8.,zi,cosi)
        WRITE(*,*) 'Initial cosmology sigma 8 (at your z):', sig8i
        WRITE(*,*)
     END IF

     !  outfile='outdat'
     !  CALL write_gadget(x,v,id,L,om_m,om_v,m,a,z,n,outfile)
     !  STOP

     !  WRITE(*,*) '1 - zs only'
     !  WRITE(*,*) '2 - zsd'
     !  READ(*,*) imode
     !  WRITE(*,*)

     !  WRITE(*,*) 'Simulation particles:', n
     !  WRITE(*,*) 'Cube root of this:', n**(1./3.)
     !  WRITE(*,*)

     !  WRITE(*,*) 'Target cosmology'
     !  CALL assign_cosmology(cost,zt)

     IF(isim==1) THEN

        !Only read in this crud once

        IF(params%izsd==1) THEN
           CALL input_camb_pk(params%pk_target,cost)
           CALL input_growth(params%growth_target,cost)
           CALL normalise(cost)
        END IF

        CALL input_rate(params%rate_target,cost)

     END IF

     at=1./(1.+zt)

     IF(params%izsd==1) THEN
        sig8f=sigma(8.,zt,cost)
        WRITE(*,*) 'Target sigma 8 (at your z):', sig8f
        WRITE(*,*)
     END IF

     !  WRITE(*,*) 'Target box size (Mpc/h):'
     !  READ(*,*) Lt
     !  WRITE(*,*)

     s=Lt/Li

     ms=(s**3.)*cost%om_m/cosi%om_m

     kbox=(2.*pi)/(s*Li)
     vfac=s*(at/ai)*sqrt(hubble2(zt,cost)/hubble2(zi,cosi))*(grow_rate(kbox,zt,cost)/grow_rate(s*kbox,zi,cosi))

     WRITE(*,*) 'Original growth rate at box scale:', grow_rate(s*kbox,zi,cosi)
     WRITE(*,*) 'Target growth rate at box scale:', grow_rate(kbox,zt,cost)

     !Do all the bulk scaling of box quantities
     x=s*x
     v=vfac*v
     mt=ms*mi

     WRITE(*,*) 'Original box length (Mpc/h):', Lt/s
     WRITE(*,*) 'Rescaled box length (Mpc/h):', Lt
     WRITE(*,*) 'Length scale factor:', s
     WRITE(*,*) 'Velocity scale factor:', vfac
     WRITE(*,*) 'Mass scale factor:', ms
     WRITE(*,*)

     IF(params%izsd==0) THEN

        !Just write out Gadget file and exit to do 'zs'
        CALL write_gadget(x,v,id,Lt,cost%om_m,cost%om_v,mt,at,zt,n,output_file)

        DEALLOCATE(x,v,id)

     ELSE IF(params%izsd==1) THEN

        !      WRITE(*,*) '0 - Reconstruct displacement field from evolved positions'x
        !      WRITE(*,*) '1 - Read in displacement field (should be unsmoothed and unscaled)'
        !      READ(*,*) idisp
        !      WRITE(*,*)

        !Now comes all displacement field crap for the zsd part of the method

        !     WRITE(*,*) 'Non-linear wave number (k_nl, defined in the TARGET cosmology):'
        !     READ(*,*) knl
        !     WRITE(*,*)
        
        rsm=1./knl
        WRITE(*,*) 'Smoothing length (r_sm=1/k_nl, TARGET cosmology):', rsm
        WRITE(*,*)

        IF(params%idisp==0) THEN

           !  Maybe it should be greater than this because there will still be some
           !  stuff at the Nyquist here (Gaussian filter). However, results change very little if mesh is increased
           !  Mesh is twice the Nyquist frequency smoothed on scales of k_ny by a Gaussian!
           mesh=CEILING(knl*Lt*2./pi)
           !  mesh=150

           minmesh=16
           IF(mesh<minmesh) mesh=minmesh

           WRITE(*,*) 'Mesh size:', mesh
           WRITE(*,*)

        ELSE IF(params%idisp==1) THEN

           !This reads in a displacement field that should be in the units of the original
           !Simulation (unscaled) at the correct redshift (z) - it should also be unsmoothed

           !        WRITE(*,*) 'Displacement field name:'
           !        READ(*,*) disp_field_name
           !        WRITE(*,*)

           WRITE(*,*) 'Reading in displacement field: ', TRIM(params%disp_file)

           mesh=file_length(params%disp_file)

           mesh=mesh**(1./3.)

           WRITE(*,*) 'Mesh size:', mesh
           WRITE(*,*)

           ALLOCATE(psi_x(mesh,mesh,mesh),psi_y(mesh,mesh,mesh),psi_z(mesh,mesh,mesh))
           psi_x=(0.d0,0.d0)
           psi_y=(0.d0,0.d0)
           psi_z=(0.d0,0.d0)

           OPEN(7,file=params%disp_file)
           DO k=1,mesh
              DO j=1,mesh
                 DO i=1,mesh
                    READ(7,*) psi_x(i,j,k), psi_y(i,j,k), psi_z(i,j,k)
                 END DO
              END DO
           END DO
           CLOSE(7)

           !In this case the displacement field already exists 
           !so needs to be scaled to target cosmology here!

           psi_x=s*psi_x
           psi_y=s*psi_y
           psi_z=s*psi_z

           !Is smoothing this field really necessary?
           CALL smooth(psi_x,mesh,rsm,Lt)
           CALL smooth(psi_y,mesh,rsm,Lt)
           CALL smooth(psi_z,mesh,rsm,Lt)

        END IF

        kny=pi*float(mesh)/Lt

        WRITE(*,*) 'Nyquist frequency of mesh:', kny
        IF(kny<knl) THEN
           WRITE(*,*) 'WARNING: Grid size of insufficient resolution to resolve N-L wavenumber'
           WRITE(*,*) 'Grid needs to be at least:', Lt*knl/pi
           READ(*,*)
        END IF
        WRITE(*,*)


        IF(params%idisp==0) THEN

           ALLOCATE(psi_x(mesh,mesh,mesh),psi_y(mesh,mesh,mesh),psi_z(mesh,mesh,mesh))

           !This reconstructs the displacement fields
           ALLOCATE(weight(n))
           weight=1.
           CALL NGP(x,Lt,weight,d,mesh)
           DEALLOCATE(weight)

           d=d*(float(mesh)**3.)/float(n)

           CALL smooth(d,mesh,rsm,Lt)

           CALL displacement(d,psi_x,psi_y,psi_z,Lt)

           var_d_theory=vard(kbox,rsm,zt,cost)
           WRITE(*,*) 'Over-density field created'
           WRITE(*,*) 'Average:', average(d)
           WRITE(*,*) 'RMS:', sqrt(variance(d))
           WRITE(*,*) 'Theoretical RMS:', sqrt(var_d_theory)
           WRITE(*,*)

           DEALLOCATE(d)

           !This is true, probably something to do with using x rather than q as a coordinate for psi
           !Psi is defined as a function of q so should be considered in Lagrangian space in which the
           !particles have essentially already been moved backwards, I think this accounts for the minus sign
           psi_x=-psi_x
           psi_y=-psi_y
           psi_z=-psi_z

        END IF

        var_psi_theory=varp(kbox,rsm,zt,cost)

        !  name='density.dat'
        !  CALL write_full_field(d,name)

        var_psi=variance(psi_x)+variance(psi_y)+variance(psi_z)

        WRITE(*,*) 'Displacement field created'
        WRITE(*,*) 'Psi RMS/(Mpc/h):', sqrt(var_psi)
        WRITE(*,*) 'Theoretical psi RMS/(Mpc/h):', sqrt(var_psi_theory)
        WRITE(*,*)

        !     WRITE(*,*) 'Fudge displacement field to have correct variance?'
        !     WRITE(*,*) '0 - No'
        !     WRITE(*,*) '1 - Yes'
        !     READ(*,*) ifudge
        !     WRITE(*,*)

        IF(params%ifudge==1) THEN

           WRITE(*,*) 'Fudging displacement fields to have the correct variance'

           fudge=sqrt(var_psi_theory/var_psi)

           psi_x=psi_x*fudge
           psi_y=psi_y*fudge
           psi_z=psi_z*fudge

           WRITE(*,*) 'Done'
           WRITE(*,*)

           var_psi=variance(psi_x)+variance(psi_y)+variance(psi_z)
           WRITE(*,*) 'Modified psi RMS / (Mpc/h):', sqrt(var_psi)
           WRITE(*,*)

        END IF

        !  name='displacement.dat'
        !  CALL write_full_field(sqrt((psi_x**2.)+(psi_y**2.)+(psi_z**2.)),name)

        !  This could be used to create psi from the velocity field rather than from
        !  the over-density, this is potentially cleaner but as of 02/2014 I have not
        !  tested this at all
        !  vfac=100.*af*sqrt(hubble2(zf,cost))*grow_rate(zf,cost)
        !  psi_x=fvx/vfac
        !  psi_y=fvy/vfac
        !  psi_z=fvz/vfac

        WRITE(*,*) 'Fourier transforming displacement fields'
        CALL fft3(psi_x,psi_x,-1)
        CALL fft3(psi_y,psi_y,-1)
        CALL fft3(psi_z,psi_z,-1)
        WRITE(*,*)

        ALLOCATE(df_x(mesh,mesh,mesh),df_y(mesh,mesh,mesh),df_z(mesh,mesh,mesh))
        ALLOCATE(dv_x(mesh,mesh,mesh),dv_y(mesh,mesh,mesh),dv_z(mesh,mesh,mesh))

        df_x=(0.d0,0.d0)
        df_y=(0.d0,0.d0)
        df_z=(0.d0,0.d0)

        dv_x=(0.d0,0.d0)
        dv_y=(0.d0,0.d0)
        dv_z=(0.d0,0.d0)

        WRITE(*,*) 'Creating differential displacement and velocity fields'

        !     varfi=varp(s*kbox,0.,zs,cosi)/3.
        !     varft=varp(kbox,0.,zf,cost)/3.

        !     WRITE(*,*) 'Initial displacement field RMS/(Mpc/h):', sqrt(varfi)
        !     WRITE(*,*) ' Target displacement field RMS/(Mpc/h):', sqrt(varft)

        DO k=1,mesh
           DO j=1,mesh
              DO i=1,mesh

                 CALL k_fft(i,j,k,mesh,kx,ky,kz,kmod,Lt)

                 IF(i==mesh/2+1 .OR. j==mesh/2+1 .OR. k==mesh/2+1) THEN

                    !Set the dpsi and dv Nyquist frequencies to 0

                    df_x(i,j,k)=(0.d0,0.d0)
                    df_y(i,j,k)=(0.d0,0.d0)
                    df_z(i,j,k)=(0.d0,0.d0)

                    dv_x(i,j,k)=(0.d0,0.d0)
                    dv_y(i,j,k)=(0.d0,0.d0)
                    dv_z(i,j,k)=(0.d0,0.d0)

                 ELSE IF(i==1 .AND. j==1 .AND. k==1) THEN

                    !Set the DC mode to 0

                    df_x(i,j,k)=(0.d0,0.d0)
                    df_y(i,j,k)=(0.d0,0.d0)
                    df_z(i,j,k)=(0.d0,0.d0)

                    dv_x(i,j,k)=(0.d0,0.d0)
                    dv_y(i,j,k)=(0.d0,0.d0)
                    dv_z(i,j,k)=(0.d0,0.d0)

                 ELSE

                    !The initial and target power spectra at k
                    powi=plin(s*kmod,zi,cosi)
                    powt=plin(kmod,zt,cost)

                    !The initial and target growth-rate at k
                    !Should be removed for non scale-dependent growth
                    ratei=grow_rate(s*kmod,zi,cosi)
                    ratet=grow_rate(kmod,zt,cost)

                    !The modification factors for the displacement field and growth rate
                    modfac=sqrt(powt/powi)-1.
                    modfacv=(ratet/ratei)*sqrt(powt/powi)-1.

                    !Now assign the dpsi field values
                    df_x(i,j,k)=psi_x(i,j,k)*modfac
                    df_y(i,j,k)=psi_y(i,j,k)*modfac
                    df_z(i,j,k)=psi_z(i,j,k)*modfac

                    !And the dv field values
                    !Extra factors of a, H, f_g that do not depend on k are multiplied below
                    dv_x(i,j,k)=psi_x(i,j,k)*modfacv
                    dv_y(i,j,k)=psi_y(i,j,k)*modfacv
                    dv_z(i,j,k)=psi_z(i,j,k)*modfacv

                 END IF

              END DO
           END DO
        END DO

        DEALLOCATE(psi_x,psi_y,psi_z)

        !Convert from displacement to velocity units for the velocity field!
        dv_x=dv_x*at*sqrt(hubble2(zt,cost))*grow_rate(kbox,zt,cost)
        dv_y=dv_y*at*sqrt(hubble2(zt,cost))*grow_rate(kbox,zt,cost)
        dv_z=dv_z*at*sqrt(hubble2(zt,cost))*grow_rate(kbox,zt,cost)

        WRITE(*,*) 'Fourier differnetial displacement and velocity field calculated'
        WRITE(*,*)

        !Return the dpsi field to real space
        WRITE(*,*) 'Starting FFT for real differential displacement fields'
        CALL fft3(df_x,df_x,1)
        CALL fft3(df_y,df_y,1)
        CALL fft3(df_z,df_z,1)

        !Return the dv field to real space
        WRITE(*,*) 'Starting FFT for real differential velocity fields'
        CALL fft3(dv_x,dv_x,1)
        CALL fft3(dv_y,dv_y,1)
        CALL fft3(dv_z,dv_z,1)

        !Normalise for FFTW conventions
        df_x=df_x/float(mesh)**3.
        df_y=df_y/float(mesh)**3.
        df_z=df_z/float(mesh)**3.
        dv_x=dv_x/float(mesh)**3.
        dv_y=dv_y/float(mesh)**3.
        dv_z=dv_z/float(mesh)**3.

        WRITE(*,*) 'FFT complete'
        WRITE(*,*)

        WRITE(*,*) 'Moving particles according to modulated displacement fields'

        avd=0.d0
        avv=0.d0
        DO i=1,n

           !Cell that particle resides in for NGP binning
           ix=CEILING(x(1,i)*float(mesh)/Lt)
           iy=CEILING(x(2,i)*float(mesh)/Lt)
           iz=CEILING(x(3,i)*float(mesh)/Lt)

           dx=df_x(ix,iy,iz)
           dy=df_y(ix,iy,iz)
           dz=df_z(ix,iy,iz)

           avd=avd+sqrt(dx**2.+dy**2.+dz**2.)

           dvx=dv_x(ix,iy,iz)
           dvy=dv_y(ix,iy,iz)
           dvz=dv_z(ix,iy,iz)

           avv=avv+sqrt(dvx**2.+dvy**2.+dvz**2.)

           x(1,i)=x(1,i)+dx
           x(2,i)=x(2,i)+dy
           x(3,i)=x(3,i)+dz

           v(1,i)=v(1,i)+dvx
           v(2,i)=v(2,i)+dvy
           v(3,i)=v(3,i)+dvz

        END DO
        avd=avd/float(n)
        avv=avv/float(n)

        CALL replace(x,Lt)

        !Deallocate the differential field arrays
        DEALLOCATE(df_x,df_y,df_z,dv_x,dv_y,dv_z)

        WRITE(*,*) 'Average displacement per paticle /(Mpc/h):', avd
        WRITE(*,*) 'Average speed change per particle /(km/s):', avv
        WRITE(*,*) 'Particles moved'
        WRITE(*,*)

        CALL write_gadget(x,v,id,Lt,cost%om_m,cost%om_v,mt,at,zt,n,output_file)
        DEALLOCATE(x,v,id)

     END IF

  END DO

CONTAINS

  SUBROUTINE read_params(params,pfile)

    USE cosdef
    IMPLICIT NONE
    TYPE(parameters) :: params
    CHARACTER(len=256), INTENT(IN) :: pfile
    CHARACTER(len=256), ALLOCATABLE :: names(:)
    CHARACTER(len=256) :: spam
    INTEGER :: i, j, n

    !    params%input_file=''
    !    params%output_file=''
    params%disp_file=''

    params%pk_initial=''
    params%pk_target=''
    params%growth_initial=''
    params%growth_target=''
    params%rate_initial=''
    params%rate_target=''

    !    params%z_target=-1.
    params%L_target=-1.
    params%sig8_initial=-1.
    params%sig8_target=-1.
    params%omm_target=-1.
    params%omv_target=-1.
    !    params%knl_target=-1.

    params%izsd=-1
    params%ifudge=-1
    params%idisp=-1

    n=file_length(pfile)

    WRITE(*,*) 'Parameters file length:', n

    ALLOCATE(names(n))

    OPEN(7,file=pfile)
    DO i=1,n
       READ(7,*) names(i), spam
    END DO
    CLOSE(7)

    DO i=1,n

       IF(names(i)=='n_files') THEN

          OPEN(7,file=pfile)

          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%n_files
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO

          ALLOCATE(params%z_target(params%n_files), params%knl_target(params%n_files))
          ALLOCATE(params%input_file(params%n_files), params%output_file(params%n_files))

          params%z_target=-1.
          params%knl_target=-1.
          params%input_file=''
          params%output_file=''

          WRITE(*,*) 'Number of files to rescale:', params%n_files

          DO j=1,params%n_files
             READ(7,*) spam, params%input_file(j), params%output_file(j), params%z_target(j), params%knl_target(j)
             WRITE(*,*) 'Input:', j, TRIM(params%input_file(j))
          END DO
          WRITE(*,*)

          CLOSE(7)

       END IF

    END DO

    !    IF(ALLOCATED(params%z_initial) .EQV. .FALSE.) STOP 'Error: arrays not correctly allocated'
    IF(ALLOCATED(params%z_target) .EQV. .FALSE.) STOP 'Error: arrays not correctly allocated'
    IF(ALLOCATED(params%knl_target) .EQV. .FALSE.) STOP 'Error: arrays not correctly allocated'
    IF(ALLOCATED(params%input_file) .EQV. .FALSE.) STOP 'Error: arrays not correctly allocated'
    IF(ALLOCATED(params%output_file) .EQV. .FALSE.) STOP 'Error: arrays not correctly allocated'
    IF(params%n_files==-1) STOP 'Error: parameters file does not contain the number of files to rescale'

    DO i=1,n


!!$       IF(names(i)=='input_file') THEN
!!$
!!$          OPEN(7,file=pfile)
!!$          DO j=1,i
!!$             IF(j==i) THEN
!!$                READ(7,*) spam, params%input_file
!!$                EXIT
!!$             ELSE
!!$                READ(7,*) spam
!!$             END IF
!!$          END DO
!!$          CLOSE(7)
!!$
!!$       ELSE IF(names(i)=='output_file') THEN
!!$
!!$          OPEN(7,file=pfile)
!!$          DO j=1,i
!!$             IF(j==i) THEN
!!$                READ(7,*) spam, params%output_file
!!$                EXIT
!!$             ELSE
!!$                READ(7,*) spam
!!$             END IF
!!$          END DO
!!$          CLOSE(7)

       IF(names(i)=='disp_file') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%disp_file
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)

       ELSE IF(names(i)=='pk_initial') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%pk_initial
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)

       ELSE IF(names(i)=='pk_target') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%pk_target
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)    

       ELSE IF(names(i)=='growth_initial') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%growth_initial
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)

       ELSE IF(names(i)=='growth_target') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%growth_target
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)     

       ELSE IF(names(i)=='rate_initial') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%rate_initial
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)

       ELSE IF(names(i)=='rate_target') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%rate_target
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)  

!!$       ELSE IF(names(i)=='z_target') THEN
!!$
!!$          OPEN(7,file=pfile)
!!$          DO j=1,i
!!$             IF(j==i) THEN
!!$                READ(7,*) spam, params%z_target
!!$                EXIT
!!$             ELSE
!!$                READ(7,*) spam
!!$             END IF
!!$          END DO
!!$          CLOSE(7)

       ELSE IF(names(i)=='L_target') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%L_target
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)

       ELSE IF(names(i)=='sig8_initial') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%sig8_initial
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)

       ELSE IF(names(i)=='sig8_target') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%sig8_target
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)

       ELSE IF(names(i)=='omm_target') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%omm_target
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)

       ELSE IF(names(i)=='omv_target') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%omv_target
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)

!!$       ELSE IF(names(i)=='knl_target') THEN
!!$
!!$          OPEN(7,file=pfile)
!!$          DO j=1,i
!!$             IF(j==i) THEN
!!$                READ(7,*) spam, params%knl_target
!!$                EXIT
!!$             ELSE
!!$                READ(7,*) spam
!!$             END IF
!!$          END DO
!!$          CLOSE(7)

       ELSE IF(names(i)=='izsd') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%izsd
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)

       ELSE IF(names(i)=='ifudge') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%ifudge
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)

       ELSE IF(names(i)=='idisp') THEN

          OPEN(7,file=pfile)
          DO j=1,i
             IF(j==i) THEN
                READ(7,*) spam, params%idisp
                EXIT
             ELSE
                READ(7,*) spam
             END IF
          END DO
          CLOSE(7)


       END IF

    END DO

    !    IF(params%input_file=='') STOP 'Error: parameters file does not contain an input file'
    !    IF(params%output_file=='') STOP 'Error: parameters file does not contain an output file'
    IF(params%idisp==1 .AND. params%disp_file=='') STOP 'Error: parameters file does not contain a displacement field file'

    IF(params%izsd==1 .AND. params%pk_initial=='') STOP 'Error: parameters file does not contain an initial P(k) file'
    IF(params%izsd==1 .AND. params%pk_target=='') STOP 'Error: parameters file does not contain a target P(k) file'

    IF(params%izsd==1 .AND. params%growth_initial=='') STOP 'Error: parameters file does not contain an initial g(z) file'
    IF(params%izsd==1 .AND. params%growth_target=='') STOP 'Error: parameters file does not contain a target g(z) file'

    IF(params%rate_initial=='') STOP 'Error: parameters file does not contain an initial fg(z) file'
    IF(params%rate_target=='') STOP 'Error: parameters file does not contain a target fg(z) file'

    !    IF(params%z_target==-1.) STOP 'Error: parameters file does not contain a target redshift'
    IF(params%L_target==-1.) STOP 'Error: parameters file does not contain a target box size'
    !    IF(params%izsd==1 .AND. params%knl_target==-1.) STOP 'Error: parameters file does not contain a target non-linear k'

    IF(params%izsd==1 .AND. params%sig8_initial==-1.) STOP 'Error: parameters file does not contain an initial sigma8'
    IF(params%izsd==1 .AND. params%sig8_target==-1.) STOP 'Error: parameters file does not contain a target sigma8'

    IF(params%omm_target==-1.) STOP 'Error: parameters file does not contain a target Om_m'
    IF(params%omv_target==-1.) STOP 'Error: parameters file does not contain a target Om_v'

    IF(params%izsd==-1) STOP 'Error: parameters file does not specify izsd'
    IF(params%izsd==1 .AND. params%ifudge==-1) STOP 'Error: parameters file does not displacement field fudge'
    IF(params%izsd==1 .AND. params%idisp==-1) STOP 'Error: parameters file does not displacement field option'


  END SUBROUTINE read_params

  SUBROUTINE input_camb_pk(input,cosm)

    USE cosdef
    IMPLICIT NONE
    CHARACTER(len=256) :: input
    INTEGER :: i, n
    REAL, PARAMETER :: pi=3.141592654
    TYPE(cosmology) :: cosm
    LOGICAL :: lexist

    INQUIRE(file=input,exist=lexist)
    IF(lexist .EQV. .FALSE.) STOP 'This input CAMB file does not exist'

    n=file_length(input)

    cosm%np=n

    WRITE(*,*) 'Reading in power spectrum: ', TRIM(input)
    WRITE(*,*) 'Length:', n

    !Allocate the arrays for power spectrum
    ALLOCATE(cosm%kp_tab(n),cosm%p_tab(n))

    OPEN(7,file=input)
    DO i=1,n
       READ(7,*) cosm%kp_tab(i), cosm%p_tab(i)
    END DO
    CLOSE(7)

    WRITE(*,*) 'k_min:', cosm%kp_tab(1)
    WRITE(*,*) 'k_max:', cosm%kp_tab(n)

    WRITE(*,*) 'Converting to D^2(k)'
    !Convert P(k) -> D(k)
    cosm%p_tab=cosm%p_tab*(cosm%kp_tab**3.)/(2.*pi*pi)

    WRITE(*,*) 'Done'
    WRITE(*,*)

  END SUBROUTINE input_camb_pk

  SUBROUTINE input_growth(input,cosm)

    USE cosdef
    IMPLICIT NONE
    CHARACTER(len=256) :: input
    INTEGER :: i, n
    TYPE(cosmology) :: cosm
    LOGICAL :: lexist

    INQUIRE(file=input,exist=lexist)
    IF(lexist .EQV. .FALSE.) STOP 'This input growth function file does not exist'

    n=file_length(input)

    cosm%ng=n

    WRITE(*,*) 'Reading in growth function: ', TRIM(input)
    WRITE(*,*) 'Length:', n

    !Allocate the arrays for growth function
    ALLOCATE(cosm%ag_tab(n),cosm%g_tab(n),cosm%zg_tab(n))

    OPEN(7,file=input)
    DO i=1,n
       READ(7,*) cosm%ag_tab(i), cosm%g_tab(i)
    END DO
    CLOSE(7)

    WRITE(*,*) 'a_min:', cosm%ag_tab(1)
    WRITE(*,*) 'a_max:', cosm%ag_tab(n)

    CALL reverse(cosm%ag_tab)
    CALL reverse(cosm%g_tab)

    cosm%zg_tab=-1.+1./cosm%ag_tab

    WRITE(*,*) 'Done'
    WRITE(*,*)

  END SUBROUTINE input_growth

  SUBROUTINE input_rate(input,cosm)

    USE cosdef
    IMPLICIT NONE
    CHARACTER(len=256) :: input
    INTEGER :: i, n
    TYPE(cosmology) :: cosm
    LOGICAL :: lexist

    INQUIRE(file=input,exist=lexist)
    IF(lexist .EQV. .FALSE.) STOP 'This input growth rate file does not exist'

    n=file_length(input)

    cosm%ng=n

    WRITE(*,*) 'Reading in growth rate: ', TRIM(input)
    WRITE(*,*) 'Length:', n

    !Allocate the arrays for growth function
    ALLOCATE(cosm%afg_tab(n),cosm%fg_tab(n),cosm%zfg_tab(n))

    OPEN(7,file=input)
    DO i=1,n
       READ(7,*) cosm%afg_tab(i), cosm%fg_tab(i)
    END DO
    CLOSE(7)

    WRITE(*,*) 'a_min:', cosm%afg_tab(1)
    WRITE(*,*) 'a_max:', cosm%afg_tab(n)

    CALL reverse(cosm%afg_tab)
    CALL reverse(cosm%fg_tab)

    cosm%zfg_tab=-1.+1./cosm%afg_tab

    WRITE(*,*) 'Done'
    WRITE(*,*)

  END SUBROUTINE input_rate

  SUBROUTINE replace(x,L)

    IMPLICIT NONE
    REAL :: x(:,:), L
    INTEGER :: i, j, n

    WRITE(*,*) 'Replacing particles that may have strayed'

    n=SIZE(x(1,:))

    DO i=1,n
       DO j=1,3
          IF(x(j,i)>L) x(j,i)=x(j,i)-L
          IF(x(j,i)<=0.) x(j,i)=x(j,i)+L
       END DO
    END DO

    WRITE(*,*) 'Done'
    WRITE(*,*)

  END SUBROUTINE replace

  FUNCTION vard(kbox,r,z,cosm)

    USE cosdef
    IMPLICIT NONE
    TYPE(cosmology) :: cosm
    REAL :: kbox, r, z, vard
    INTEGER :: i, n
    REAL, ALLOCATABLE :: integrand(:)
    REAL :: e, g, p

    IF(cosm%ipow==0) THEN

       n=SIZE(cosm%kp_tab)
       ALLOCATE(integrand(n))

       DO i=1,n

          IF(cosm%kp_tab(i)<kbox) THEN

             integrand(i)=0.

          ELSE

             e=exp(-((cosm%kp_tab(i)*r)**2.))
             g=grow(cosm%kp_tab(i),z,cosm)**2.
             p=cosm%p_tab(i)/cosm%kp_tab(i)

             integrand(i)=e*g*p

          END IF

       END DO

       !Should be converted to log k integral
       vard=inttab(cosm%kp_tab,integrand,3)

    ELSE IF(cosm%ipow==1) THEN

       vard=integrate_vard(kbox,100.,r,z,cosm,0.001)

    END IF

  END FUNCTION vard

  FUNCTION varp(kbox,r,z,cosm)

    USE cosdef
    IMPLICIT NONE
    TYPE(cosmology) :: cosm
    REAL :: kbox, r, z, varp
    INTEGER :: i, n
    REAL, ALLOCATABLE :: integrand(:)
    REAL :: e, g, p

    IF(cosm%ipow==0) THEN

       n=SIZE(cosm%kp_tab)
       ALLOCATE(integrand(n))

       DO i=1,n

          IF(cosm%kp_tab(i)<kbox) THEN

             integrand(i)=0.

          ELSE

             e=exp(-((cosm%kp_tab(i)*r)**2.))
             g=grow(cosm%kp_tab(i),z,cosm)**2.
             p=cosm%p_tab(i)/(cosm%kp_tab(i)**3.)

             integrand(i)=e*g*p

          END IF

       END DO

       !Should be converted to log k integral
       varp=inttab(cosm%kp_tab,integrand,3)

    ELSE IF(cosm%ipow==1) THEN

       varp=integrate_varp(kbox,100.,r,z,cosm,0.001)

    END IF

  END FUNCTION varp

  SUBROUTINE write_full_field(f,file_name)

    IMPLICIT NONE
    DOUBLE COMPLEX :: f(:,:,:)
    INTEGER :: i, j, k, m
    CHARACTER(len=256) :: file_name

    WRITE(*,*) 'Writing out full field'

    m=SIZE(f(:,1,1))

    WRITE(*,*) 'Size:', m**3

    OPEN(7,file=file_name)
    DO k=1,m
       DO j=1,m
          DO i=1,m
             WRITE(7,*) real(f(i,j,k))
          END DO
       END DO
    END DO
    CLOSE(7)

    WRITE(*,*) 'Output:', file_name
    WRITE(*,*) 'Done'
    WRITE(*,*)

  END SUBROUTINE write_full_field

  SUBROUTINE read_gadget(pos,vel,id,L,om_m,om_v,m,a,z,n,infile)

    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: infile
    REAL*8 :: massarr(6), z8, a8, L8, om_m8, om_v8
    REAL, ALLOCATABLE, INTENT(OUT) :: pos(:,:), vel(:,:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: id(:)
    REAL, INTENT(OUT) :: om_m, om_v, m, a, z, L
    INTEGER, INTENT(OUT) :: n
    INTEGER :: np(6), np2(6), crap

    WRITE(*,*) 'Reading in Gadget2 particle data: ', TRIM(infile)

    OPEN(7,file=infile,form='unformatted',status='old')
    READ(7) np, massarr, a8, z8, crap, crap, np2, crap, crap, L8, om_m8, om_v8
    CLOSE(7)

    m=massarr(2)
    a=a8
    z=z8
    om_m=om_m8
    om_v=om_v8
    n=np(2)

    L=L8/1000.
    m=m*1.e10

    WRITE(*,*) 'Particle number:', n
    WRITE(*,*) 'Which is:', nint(n**(1./3.)), 'cubed'
    WRITE(*,*) 'Box size (Mpc/h):', L
    WRITE(*,*) 'a:', a
    WRITE(*,*) 'z:', z
    WRITE(*,*) 'Particle mass (Msun):', m
    WRITE(*,*) 'Om_m:', om_m
    WRITE(*,*) 'Om_v:', om_v

    ALLOCATE(pos(3,n),vel(3,n),id(n))

    OPEN(7,file=infile,form='unformatted')
    READ(7)
    READ(7) pos
    READ(7) vel
    READ(7) id
    CLOSE(7)

    !kpc -> Mpc conversion!
    x=x/1000.
    v=v*sqrt(a8)

    WRITE(*,*) 'Finished reading in file'
    WRITE(*,*)

  END SUBROUTINE read_gadget

  SUBROUTINE write_gadget(pos,vel,id,L,om_m,om_v,m,a,z,n,outfile)

    IMPLICIT NONE
    CHARACTER(len=256) :: outfile
    REAL*8 :: massarr(6), z8, a8, L8, crap8, om_m8, om_v8
    REAL, INTENT(INOUT) :: pos(:,:), vel(:,:)
    REAL, INTENT(IN) :: L, a, z, om_m, om_v, m
    INTEGER, INTENT(IN) :: n, id(:)
    INTEGER :: np(6), np2(6), crapi

    WRITE(*,*) 'Outputting particle data in Gadget2 format: ', TRIM(outfile)

    np=0
    np(2)=n
    massarr=0.d0
    a8=a
    z8=z
    crapi=0
    crap8=0.d0
    om_m8=om_m
    om_v8=om_v

    L8=L*1000.d0
    massarr(2)=m/1.e10

    WRITE(*,*) 'Particle number:', n
    WRITE(*,*) 'Which is:', nint(n**(1./3.)), 'cubed'
    WRITE(*,*) 'Box size (Mpc/h):', L
    WRITE(*,*) 'a:', a
    WRITE(*,*) 'z:', z
    WRITE(*,*) 'Particle mass:', m
    WRITE(*,*) 'Om_m:', om_m
    WRITE(*,*) 'Om_v:', om_v

    pos=pos*1000.
    vel=vel/sqrt(a)

    OPEN(7,file=outfile,form='unformatted')
    WRITE(7) np, massarr, a8, z8, crapi, crapi, np2, crapi, crapi, L8, om_m8, om_v8
    WRITE(7) pos
    WRITE(7) vel
    WRITE(7) id
    CLOSE(7)

    !Incase these are to be used again outwith the subroutine
    pos=pos/1000.
    vel=vel*sqrt(a)

    WRITE(*,*) 'Finished writing file'
    WRITE(*,*)

  END SUBROUTINE write_gadget

!!$  SUBROUTINE write_gadget(x,y,z,vx,vy,vz,m,L,om_m,om_v,zr,n,a,outfile)
!!$
!!$    IMPLICIT NONE
!!$    CHARACTER(len=256) :: outfile
!!$    REAL, INTENT(IN) :: x(:), y(:), z(:), vx(:), vy(:), vz(:), m(:), L, a, om_m, om_v, zr
!!$    INTEGER, INTENT(IN) :: n
!!$    REAL*8 :: massarr(6), z8, a8, L8, c8, m8, om_m8, om_v8
!!$    REAL*4, ALLOCATABLE :: pos(:,:), vel(:,:), mass(:)
!!$    INTEGER :: np(6), i, np2(6), ci, crap(15)
!!$
!!$    WRITE(*,*) 'Outputting particle data in Gadget2 format'
!!$    !    WRITE(*,*) m(1)
!!$
!!$    np=0
!!$    np(2)=n
!!$    np2=0
!!$    np2(2)=n
!!$    massarr=0.d0
!!$    massarr(2)=m(1)/1e10
!!$    z8=zr
!!$    a8=a
!!$    L8=L*1000.d0
!!$    ci=0
!!$    c8=0.d0
!!$    crap=0
!!$    om_m8=om_m
!!$    om_v8=om_v
!!$
!!$    ALLOCATE(pos(3,n),vel(3,n))
!!$
!!$    DO i=1,n
!!$       pos(1,i)=x(i)
!!$       pos(2,i)=y(i)
!!$       pos(3,i)=z(i)
!!$       vel(1,i)=vx(i)
!!$       vel(2,i)=vy(i)
!!$       vel(3,i)=vz(i)
!!$    END DO
!!$
!!$    pos=pos*1000.
!!$    vel=vel/sqrt(a)
!!$
!!$    OPEN(7,file=outfile,form='unformatted')
!!$    WRITE(7) np, massarr, a8, z8, ci, ci, np2, ci, ci, L8, om_m8, om_v8, c8, ci, ci, np2, ci, crap
!!$    WRITE(7) pos
!!$    WRITE(7) vel
!!$    WRITE(7) id
!!$    CLOSE(7)
!!$
!!$    DEALLOCATE(pos,vel)
!!$
!!$    WRITE(*,*) 'Finished writing Gadget2 file:', outfile
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE write_gadget

  FUNCTION file_length(file_name)

    IMPLICIT NONE
    CHARACTER(len=256), INTENT(IN) :: file_name
    INTEGER :: n, file_length
    REAL :: data

    OPEN(7,file=file_name)
    n=0
    DO
       n=n+1
       READ(7,*,end=301)
    END DO

    !301 is just the label to jump to when the end of the file is reached

301 CLOSE(7)  

    file_length=n-1

  END FUNCTION file_length

!!$  FUNCTION file_length_string(file_name)
!!$
!!$    IMPLICIT NONE
!!$    CHARACTER(len=256) :: file_name
!!$    INTEGER :: n, file_length
!!$    REAL :: data
!!$
!!$    OPEN(7,file=file_name)
!!$    n=0
!!$    DO
!!$       n=n+1
!!$       READ(7,*,end=301)
!!$    END DO
!!$
!!$    !301 is just the label to jump to when the end of the file is reached
!!$
!!$301 CLOSE(7)  
!!$
!!$    file_length=n-1
!!$
!!$  END FUNCTION file_length_string

  SUBROUTINE range(arr)

    IMPLICIT NONE
    REAL :: arr(:)

    !Writes out the lowest and highest value in the array

    WRITE(*,*) 'Low:', MINVAL(arr)
    WRITE(*,*) 'High:', MAXVAL(arr)
    WRITE(*,*)

  END SUBROUTINE range

  SUBROUTINE normalise(cosm)

    USE cosdef
    IMPLICIT NONE
    TYPE(cosmology) :: cosm
    REAL, ALLOCATABLE :: win(:)
    REAL :: sig8o, sig8f, x
    INTEGER :: n

    WRITE(*,*) 'Normalising input power spectrum'

    !Calculate sigma 8 and compare to desired value
    sig8o=sigma(8.,0.,cosm)
    sig8f=cosm%sig8

    WRITE(*,*) 'Sigma 8 from input file:', sig8o
    WRITE(*,*) 'Target sigma 8:', sig8f

    !Renormalise power spectrum
    cosm%p_tab=cosm%p_tab*(sig8f/sig8o)**2.

    WRITE(*,*) 'Normalisation adjusted'
    WRITE(*,*)

  END SUBROUTINE normalise

!!$  SUBROUTINE normalise(cosm)
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    TYPE(cosmology) :: cosm
!!$    REAL :: sig8o, sig8f
!!$
!!$    cosm%A=1.
!!$    sig8o=sigma(8.,0.,cosm)
!!$    
!!$    WRITE(*,*) 'Normalising'
!!$    WRITE(*,*) 'Original sigma8 (z=0):', sig8o
!!$
!!$    IF(cosm%ipow==0) THEN             
!!$       cosm%p=cosm%p*((cosm%sig8/sig8o)**2.)
!!$    ELSE IF(cosm%ipow==1) THEN
!!$       cosm%A=cosm%sig8/sig8o
!!$    END IF
!!$
!!$    sig8f=sigma(8.,0.,cosm)
!!$
!!$    WRITE(*,*) 'Final sigma8 (z=0):', sig8f
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE normalise

  FUNCTION sigma(r,z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: r, sigma, z
    REAL :: w, y
    REAL, ALLOCATABLE :: integrand(:)
    INTEGER :: i, n
    TYPE(cosmology) :: cosm

    n=SIZE(cosm%kp_tab)

    ALLOCATE(integrand(n))

    DO i=1,n
       y=cosm%kp_tab(i)*r
       w=wk(y)
       integrand(i)=cosm%p_tab(i)*(w**2.)!/(cosm%k(i))
    END DO

    sigma=sqrt(inttab(log(cosm%kp_tab),integrand,3))
    sigma=sigma*grow(0.,z,cosm)

    DEALLOCATE(integrand)

  END FUNCTION sigma

  FUNCTION wk(x)

    IMPLICIT NONE
    REAL :: wk, x

    IF(x<0.) STOP 'WK: Error, argument should be greater than 0'

    !The Fourier transform of a tophat function
    IF(x<0.1) THEN
       wk=1.-(x**2.)/10.
    ELSE
       wk=3.*(sin(x)-x*cos(x))/(x**3.)
    END IF

  END FUNCTION wk

!!$  FUNCTION sigma(r,z,cosm)
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    REAL :: sigma
!!$    REAL, INTENT(IN) :: r, z
!!$    REAL :: y, wk
!!$    REAL, ALLOCATABLE :: integrand(:)
!!$    TYPE(cosmology) :: cosm
!!$    INTEGER :: i, n
!!$
!!$    IF(cosm%ipow==0) THEN
!!$
!!$       n=SIZE(cosm%k)
!!$       ALLOCATE(integrand(n))
!!$   
!!$       !Construct the integrand
!!$       DO i=1,n
!!$          y=cosm%k(i)*r
!!$          wk=tophat_k(y)
!!$          integrand(i)=cosm%p(i)*(wk**2.)!/cosm%k(i)
!!$       END DO
!!$
!!$       sigma=sqrt(inttab(log(cosm%k),integrand,3))
!!$
!!$       DEALLOCATE(integrand)
!!$
!!$       sigma=sigma*grow(0.,z,cosm)
!!$
!!$    ELSE IF(cosm%ipow==1) THEN
!!$
!!$       sigma=sigint(r,z,cosm)
!!$
!!$    END IF
!!$
!!$  END FUNCTION sigma

!!$  FUNCTION tophat_k(x)
!!$
!!$    IMPLICIT NONE
!!$    REAL :: tophat_k
!!$    REAL, INTENT(IN) :: x
!!$
!!$    !Use Maclaurin expansion for x<0.1 to avoid errors in cancelling
!!$    !Otherwise use exact expression
!!$    !x should be >0
!!$
!!$    IF(x<0.1) THEN
!!$       tophat_k=1.-(x**2.)/10.
!!$    ELSE
!!$       tophat_k=3.*(sin(x)-x*cos(x))/(x**3.)
!!$    END IF
!!$
!!$  END FUNCTION tophat_k

!!$  FUNCTION sigint(r,z,cosm)
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    REAL :: sigint, r, z
!!$    REAL*8 :: sum, w_hat, sumold
!!$    REAL :: t, k, y, acc
!!$    INTEGER :: i, nint, j, n_init
!!$    TYPE(cosmology) :: cosm
!!$
!!$    acc=0.0001
!!$    n_init=10
!!$
!!$    DO j=1,100
!!$
!!$       nint=n_init*(2**(j-1))
!!$       sum=0.d0
!!$
!!$       DO i=1,nint-1
!!$
!!$          t=float(i)/float(nint)
!!$          k=(-1.+1./t)*(r**-0.3)
!!$          y=k*r
!!$!          w_hat=3.*(sin(y)-y*cos(y))/(y**3.)
!!$          w_hat=tophat_k(y)
!!$          sum=sum+plin(k,z,cosm)*(w_hat**2.)/(t*(1.-t))
!!$
!!$       END DO
!!$
!!$       sum=sum/float(nint)
!!$
!!$       IF(j>1 .AND. ABS(-1.+sum/sumold)<acc) THEN
!!$          sigint=sqrt(sum)
!!$          EXIT
!!$       ELSE
!!$          sumold=sum
!!$       END IF
!!$
!!$    END DO
!!$
!!$  END FUNCTION sigint

!!$  FUNCTION sigtab(r,z,cosm)
!!$
!!$    IMPLICIT NONE
!!$    REAL :: sigtab, r, z
!!$    REAL*8, ALLOCATABLE :: integrand(:)
!!$    REAL*8 :: sum, w2, y
!!$    INTEGER :: i
!!$    TYPE(cosmology) :: cosm
!!$
!!$    sum=0.d0
!!$
!!$    ALLOCATE(integrand(SIZE(cosm%k)))
!!$
!!$    DO i=1,SIZE(cosm%k)
!!$       y=r*cosm%k(i)
!!$!       w2=3.*(sin(y)-y*cos(y))/(y**3.)
!!$       w2=tophat_k(y)
!!$       w2=w2**2.
!!$       integrand(i)=((grow(cosm%k(i),z,cosm))**2.)*cosm%p(i)*w2/cosm%k(i)
!!$    END DO
!!$
!!$    DO i=1,SIZE(cosm%k)-1
!!$       sum=sum+(integrand(i+1)+integrand(i))*(cosm%k(i+1)-cosm%k(i))/2.
!!$    END DO
!!$
!!$    DEALLOCATE(integrand)
!!$
!!$    sigtab=sqrt(sum)
!!$
!!$  END FUNCTION sigtab

  FUNCTION inttab(x,y,iorder)

    IMPLICIT NONE
    REAL :: inttab
    REAL, INTENT(IN) :: x(:), y(:)
    REAL :: a, b, c, d, h
    REAL :: q1, q2, q3, qi, qf
    REAL :: x1, x2, x3, x4, y1, y2, y3, y4, xi, xf
    REAL*8 :: sum
    INTEGER :: i, n, i1, i2, i3, i4
    INTEGER, INTENT(IN) :: iorder

    n=SIZE(x)

    IF(n .NE. SIZE(y)) STOP 'Tables must be of the same length'

    sum=0.d0

    IF(iorder==1) THEN

       !Sums over all Trapezia (a+b)*h/2
       DO i=1,n-1
          a=y(i+1)
          b=y(i)
          h=x(i+1)-x(i)
          sum=sum+(a+b)*h/2.d0
       END DO

    ELSE IF(iorder==2) THEN

       DO i=1,n-2

          x1=x(i)
          x2=x(i+1)
          x3=x(i+2)

          y1=y(i)
          y2=y(i+1)
          y3=y(i+2)

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          q1=a*(x1**3.)/3.+b*(x1**2.)/2.+c*x1
          q2=a*(x2**3.)/3.+b*(x2**2.)/2.+c*x2
          q3=a*(x3**3.)/3.+b*(x3**2.)/2.+c*x3

          !Takes value for first and last sections but averages over sections where you
          !have two independent estimates of the area
          IF(n==3) THEN
             sum=sum+q3-q1
          ELSE IF(i==1) THEN
             sum=sum+(q2-q1)+(q3-q2)/2.d0
          ELSE IF(i==n-2) THEN
             sum=sum+(q2-q1)/2.d0+(q3-q2)
          ELSE
             sum=sum+(q3-q1)/2.
          END IF

       END DO

    ELSE IF(iorder==3) THEN

       DO i=1,n-1

          !First choose the integers used for defining cubics for each section
          !First and last are different because the section does not lie in the *middle* of a cubic

          IF(i==1) THEN

             i1=1
             i2=2
             i3=3
             i4=4

          ELSE IF(i==n-1) THEN

             i1=n-3
             i2=n-2
             i3=n-1
             i4=n

          ELSE

             i1=i-1
             i2=i
             i3=i+1
             i4=i+2

          END IF

          x1=x(i1)
          x2=x(i2)
          x3=x(i3)
          x4=x(i4)

          y1=y(i1)
          y2=y(i2)
          y3=y(i3)
          y4=y(i4)

          CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

          !These are the limits of the particular section of integral
          xi=x(i)
          xf=x(i+1)

          qi=a*(xi**4.)/4.+b*(xi**3.)/3.+c*(xi**2.)/2.+d*xi
          qf=a*(xf**4.)/4.+b*(xf**3.)/3.+c*(xf**2.)/2.+d*xf

          sum=sum+qf-qi

       END DO

    END IF

    inttab=sum

  END FUNCTION inttab

  SUBROUTINE fit_line(a1,a0,x1,y1,x2,y2)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1
    REAL, INTENT(IN) :: x1, y1, x2, y2

    a1=(y2-y1)/(x2-x1)
    a0=y1-a1*x1

  END SUBROUTINE fit_line

  SUBROUTINE fit_quadratic(a2,a1,a0,x1,y1,x2,y2,x3,y3)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a0, a1, a2
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3

    a2=((y2-y1)/(x2-x1)-(y3-y1)/(x3-x1))/(x2-x3)
    a1=(y2-y1)/(x2-x1)-a2*(x2+x1)
    a0=y1-a2*(x1**2.d0)-a1*x1

  END SUBROUTINE fit_quadratic

  SUBROUTINE fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)

    IMPLICIT NONE
    REAL, INTENT(OUT) :: a, b, c, d
    REAL, INTENT(IN) :: x1, y1, x2, y2, x3, y3, x4, y4
    REAL :: f1, f2, f3

    f1=(y4-y1)/((x4-x2)*(x4-x1)*(x4-x3))
    f2=(y3-y1)/((x3-x2)*(x3-x1)*(x4-x3))
    f3=(y2-y1)/((x2-x1)*(x4-x3))*(1./(x4-x2)-1./(x3-x2))

    a=f1-f2-f3

    f1=(y3-y1)/((x3-x2)*(x3-x1))
    f2=(y2-y1)/((x2-x1)*(x3-x2))
    f3=a*(x3+x2+x1)

    b=f1-f2-f3

    f1=(y4-y1)/(x4-x1)
    f2=a*(x4**2.+x4*x1+x1**2.)
    f3=b*(x4+x1)

    c=f1-f2-f3

    d=y1-a*x1**3.-b*x1**2.-c*x1

  END SUBROUTINE fit_cubic

  FUNCTION plin(k,z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: k, z, plin
    TYPE(cosmology) :: cosm

    plin=exp(find(log(k),log(cosm%kp_tab),log(cosm%p_tab),3,3))
    plin=(grow(0.,z,cosm)**2.)*plin

  END FUNCTION plin

!!$  FUNCTION plin(k,z,cosm)
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    REAL :: plin, k, z
!!$    REAL :: A, n, g
!!$    TYPE(cosmology) :: cosm
!!$
!!$    g=grow(k,z,cosm)
!!$
!!$    IF(cosm%ipow==0) THEN
!!$       plin=(g**2.)*exp(find(log(k),log(cosm%k),log(cosm%p),3,1))
!!$    ELSE IF(cosm%ipow==1) THEN
!!$       A=cosm%A
!!$       n=cosm%n
!!$       plin=(A**2.)*(g**2.)*(Tk(k,cosm)**2.)*(k**(n+3.))
!!$    END IF
!!$
!!$  END FUNCTION plin

  FUNCTION Tk(k,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: Tk, k
    TYPE(cosmology) :: cosm

    IF(cosm%itk==1) THEN
       Tk=Tk_eh(k,cosm)
    ELSE IF(cosm%itk==2) THEN
       Tk=Tk_defw(k,cosm)
    END IF

  END FUNCTION Tk

  FUNCTION Tk_defw(k,cosm)

    IMPLICIT NONE
    REAL :: tk_defw, k
    REAL :: rkeff, q, tk, gamma
    REAL*8 :: q8, tk8
    TYPE(cosmology) :: cosm

    gamma=cosm%gamma

    rkeff=0.172+0.011*log(gamma/0.36)*log(gamma/0.36)
    q=1.e-20 + k/gamma
    q8=1.e-20 + rkeff/gamma
    tk=1./(1.+(6.4*q+(3.0*q)**1.5+(1.7*q)**2)**1.13)**(1./1.13)
    tk8=1./(1.+(6.4*q8+(3.0*q8)**1.5+(1.7*q8)**2)**1.13)**(1./1.13)

    Tk_defw=tk/tk8

  END FUNCTION Tk_defw

  FUNCTION Tk_eh(k,cosm)

    ! the astonishing D.J. Eisenstein & W. Hu fitting formula (ApJ 496 605 [1998])
    ! remember I use k/h, whereas they use pure k, om_m is cdm + baryons

    USE cosdef
    IMPLICIT NONE

    REAL :: Tk_eh, k
    REAL :: om_m, om_b, h
    REAL :: rk, e, thet, b1, b2, zd, ze, rd, re, rke, s, rks
    REAL :: q
    REAL :: y, g, ab
    REAL :: a1, a2, ac
    REAL :: bc
    REAL :: f
    REAL :: c1, c2, tc
    REAL :: bb, bn, ss, tb
    TYPE(cosmology) :: cosm

    om_m=cosm%om_m
    om_b=cosm%om_b
    h=cosm%h

    rk=k*h

    e=exp(1.)

    thet=2.728/2.7
    b1=0.313*(om_m*h*h)**(-0.419)*(1+0.607*(om_m*h*h)**0.674)
    b2=0.238*(om_m*h*h)**0.223
    zd=1291.*(1+b1*(om_b*h*h)**b2)*(om_m*h*h)**0.251/(1.+0.659*(om_m*h*h)**0.828)
    ze=2.50e4*om_m*h*h/thet**4.
    rd=31500.*om_b*h*h/thet**4./zd
    re=31500.*om_b*h*h/thet**4./ze
    rke=7.46e-2*om_m*h*h/thet**2.
    s=(2./3./rke)*sqrt(6./re)*log((sqrt(1.+rd)+sqrt(rd+re))/(1+sqrt(re)))
    rks=1.6*( (om_b*h*h)**0.52 ) * ( (om_m*h*h)**0.73 ) * (1.+(10.4*om_m*h*h)**(-0.95))

    q=rk/13.41/rke

    y=(1.+ze)/(1.+zd)
    g=y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)))
    ab=g*2.07*rke*s/(1.+rd)**(0.75)

    a1=(46.9*om_m*h*h)**0.670*(1+(32.1*om_m*h*h)**(-0.532))
    a2=(12.0*om_m*h*h)**0.424*(1+(45.0*om_m*h*h)**(-0.582))
    ac=(a1**(-om_b/om_m)) * (a2**(-(om_b/om_m)**3.))

    b1=0.944/(1+(458.*om_m*h*h)**(-0.708))
    b2=(0.395*om_m*h*h)**(-0.0266)
    bc=1./(1.+b1*((1.-om_b/om_m)**b2-1.))

    f=1./(1.+(rk*s/5.4)**4.)

    c1=14.2 + 386./(1.+69.9*q**1.08)
    c2=14.2/ac + 386./(1.+69.9*q**1.08)
    tc=f*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c1*q*q) +(1.-f)*log(e+1.8*bc*q)/(log(e+1.8*bc*q)+c2*q*q)

    bb=0.5+(om_b/om_m) + (3.-2.*om_b/om_m)*sqrt((17.2*om_m*h*h)**2.+1.)
    bn=8.41*(om_m*h*h)**0.435
    ss=s/(1.+(bn/rk/s)**3.)**(1./3.)
    tb=log(e+1.8*q)/(log(e+1.8*q)+c1*q*q)/(1+(rk*s/5.2)**2.)
    tb=(tb+ab*exp(-(rk/rks)**1.4)/(1.+(bb/rk/s)**3.))*sin(rk*ss)/rk/ss

    Tk_eh=(om_b/om_m)*tb+(1-om_b/om_m)*tc

  END FUNCTION Tk_eh

  FUNCTION hubble2(z,cosm)

    !This calculates the dimensionless squared hubble parameter squared at redshift z!
    !it ignores contributions from radiation (not accurate at very high z)!

    USE cosdef
    IMPLICIT NONE
    REAL :: hubble2, z
    REAL :: om_m, om_v
    TYPE(cosmology) :: cosm

    om_m=cosm%om_m
    om_v=cosm%om_v

    !hubble2=(om_m*(1.+z)**3.)+om_v*((1.+z)**(3.*(1.+w)))+((1.-om_m-om_v)*(1.+z)**2.)
    hubble2=(om_m*(1.+z)**3.)+om_v+((1.-om_m-om_v)*(1.+z)**2.)

  END FUNCTION hubble2

  FUNCTION omega_m(z,cosm)

    !This calculates omega_m variations with z!

    USE cosdef
    IMPLICIT NONE
    REAL :: omega_m, z
    REAL :: om_m
    TYPE(cosmology) :: cosm

    om_m=cosm%om_m

    omega_m=(om_m*(1.+z)**3.)/hubble2(z,cosm)

  END FUNCTION omega_m

  FUNCTION omega_v(z,cosm)

    !This calculates omega_v variations with z!

    USE cosdef
    IMPLICIT NONE
    REAL :: omega_v, z
    REAL :: om_v
    TYPE(cosmology) :: cosm

    om_v=cosm%om_v

    !    omega_v=om_v*((1.+z)**(3.*(1.+w)))/h2(z)
    omega_v=om_v/hubble2(z,cosm)

  END FUNCTION omega_v

!!$  FUNCTION grow(k,z,cosm)
!!$
!!$    !The growth function, which is a function of k in some MG models
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    REAL :: grow, k, z
!!$    TYPE(cosmology) :: cosm
!!$    
!!$    IF(cosm%igro==0) THEN
!!$       grow=find(z,cosm%z_g,cosm%g,3,1)
!!$    ELSE IF(cosm%igro==1) THEN
!!$       grow=growint(z,cosm)
!!$    ELSE IF(cosm%igro==2) THEN
!!$       IF(z .NE. 0) THEN
!!$          WRITE(*,*) 'z=', z
!!$          STOP 'Error: attempting to calculate scale dependent linear growth at z .ne. 0'
!!$       END IF
!!$       IF(k==0.)THEN
!!$          !Take the first entry in this case
!!$          grow=cosm%gk(1)
!!$       ELSE
!!$          grow=find(log(k),log(cosm%k_g),cosm%gk,3,1)
!!$       END IF
!!$    END IF
!!$
!!$  END FUNCTION grow

  FUNCTION grow(k,z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: grow
    REAL, INTENT(IN) :: z, k
    TYPE(cosmology) :: cosm

    IF(z==0.) THEN
       grow=1.
    ELSE
       grow=find(z,cosm%zg_tab,cosm%g_tab,3,3)
    END IF

  END FUNCTION grow

  FUNCTION grow_rate(k,z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: grow_rate
    REAL, INTENT(IN) :: z, k
    TYPE(cosmology) :: cosm

    grow_rate=find(z,cosm%zfg_tab,cosm%fg_tab,3,3)

  END FUNCTION grow_rate

!!$  FUNCTION growint(z,cosm)
!!$
!!$    !Integrates the Linder approximation for growth rate to get the growth function
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    INTEGER :: i, j, jmax
!!$    REAL :: growint, z
!!$    REAL :: acc, dx, a
!!$    INTEGER :: nint, n_init
!!$    REAL :: x, fac, func, this, gam
!!$    REAL*8 :: sum1, sum2
!!$    TYPE(cosmology) :: cosm
!!$
!!$    IF(z==0.) THEN
!!$
!!$       growint=1.
!!$
!!$    ELSE
!!$
!!$       a=1./(1.+z)
!!$
!!$       sum1=0.d0
!!$       sum2=0.d0
!!$
!!$       n_init=10
!!$       acc=0.001
!!$
!!$       jmax=20
!!$
!!$       DO j=1,jmax
!!$
!!$          nint=n_init*(2**(j-1))
!!$
!!$          DO i=1,nint
!!$
!!$             x=1.+(a-1.)*((float(i)-1)/(float(nint)-1))
!!$
!!$             IF(i==1 .OR. i==nint) THEN
!!$                !multiple of 1 for beginning and end and multiple of 2 for middle points!
!!$                fac=1.
!!$             ELSE
!!$                fac=2.
!!$             END IF
!!$
!!$             !LINDER APPROXIMATION FOR RATE INTEGRATED!
!!$
!!$             !          IF(w<-1.) THEN
!!$             !             gam=0.55+0.02*(1.+w)
!!$             !          ELSE IF(w>-1) THEN
!!$             !             gam=0.55+0.05*(1.+w)
!!$             !          ELSE
!!$             gam=0.55
!!$             !          END IF
!!$
!!$             func=(omega_m(-1.+1./x,cosm)**gam)/x
!!$             sum2=sum2+fac*func
!!$
!!$          END DO
!!$
!!$          dx=(a-1.)/float(nint-1)
!!$          sum2=sum2*dx/2.
!!$
!!$          IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
!!$             growint=exp(sum2)
!!$             EXIT
!!$          ELSE
!!$             sum1=sum2
!!$             sum2=0.d0
!!$          END IF
!!$
!!$       END DO
!!$
!!$    END IF
!!$
!!$  END FUNCTION growint

!!$  SUBROUTINE fft2(in,out,ifb)
!!$
!!$    DOUBLE COMPLEX :: in(:,:), out(:,:)
!!$    INTEGER :: ifb
!!$    INTEGER :: ift
!!$    INTEGER*8 :: plan
!!$
!!$    IF(ifb .NE. 1 .AND. ifb .NE. -1) THEN
!!$       WRITE(*,*) 'Error - need to specify forwards or backwards'
!!$    END IF
!!$
!!$    ift=SIZE(in,1)
!!$
!!$    !    WRITE(*,*) 'Starting FFT - creating plan'
!!$    IF(ifb==-1) THEN
!!$       call dfftw_plan_dft_2d(plan,ift,ift,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
!!$    ELSE IF(ifb==1) THEN
!!$       call dfftw_plan_dft_2d(plan,ift,ift,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
!!$    END IF
!!$
!!$    !This computes the FFT!
!!$    !    WRITE(*,*) 'Executing FFTW'
!!$    call dfftw_execute(plan)
!!$    !    WRITE(*,*) 'FFTW complete'
!!$
!!$    !And this destroys the plan!
!!$    call dfftw_destroy_plan(plan)
!!$    !    WRITE(*,*) 'Plan destroyed'
!!$    !    WRITE(*,*) ''
!!$
!!$  END SUBROUTINE fft2

  SUBROUTINE fft3(in,out,ifb)

    DOUBLE COMPLEX :: in(:,:,:), out(:,:,:)
    INTEGER :: ifb
    INTEGER :: ift
    INTEGER*8 :: plan

    IF(ifb .NE. 1 .AND. ifb .NE. -1) THEN
       WRITE(*,*) 'Error - need to specify forwards or backwards'
    END IF

    ift=SIZE(in,1)      

    !    WRITE(*,*) 'Starting FFT - creating plan'
    IF(ifb==-1) THEN
       call dfftw_plan_dft_3d(plan,ift,ift,ift,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
    ELSE IF(ifb==1) THEN
       call dfftw_plan_dft_3d(plan,ift,ift,ift,in,out,FFTW_BACKWARD,FFTW_ESTIMATE)
    END IF

    !This computes the FFT!
    !    WRITE(*,*) 'Executing FFTW'
    call dfftw_execute(plan)
    !    WRITE(*,*) 'FFTW complete'

    !And this destroys the plan!
    call dfftw_destroy_plan(plan)
    !    WRITE(*,*) 'Plan destroyed'
    !    WRITE(*,*) ''

  END SUBROUTINE fft3

  SUBROUTINE power(d,L,bins,kmin,kmax,fname,np)

    IMPLICIT NONE
    REAL :: kx, ky, kz, modk, kv, a, b, nps, shot
    REAL, INTENT(IN) :: L, kmin, kmax
    INTEGER , INTENT(IN) :: np
    REAL, ALLOCATABLE :: bin(:), pow(:)
    !    INTEGER, ALLOCATABLE :: ninbin(:)
    DOUBLE COMPLEX, INTENT(IN) :: d(:,:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: dk(:,:,:)
    INTEGER :: i, j, k, m, n, bins
    CHARACTER(len=256) :: fname

    ALLOCATE(bin(bins+1),pow(bins))

    pow=0.

    a=log10(kmin)
    b=log10(kmax)

    WRITE(*,*) 'Computing power spectrum'

    DO i=1,bins+1
       bin(i)=a+(b-a)*float(i-1)/float(bins)
    END DO

    WRITE(*,*) 'Bins created'

    bin=10.**bin

    m=SIZE(d,1)

    ALLOCATE(dk(m,m,m))

    WRITE(*,*) 'FFTing density field'

    CALL fft3(d,dk,-1)

    DO i=1,m
       DO j=1,m
          DO k=1,m

             CALL k_fft(i,j,k,m,kx,ky,kz,modk,L)

             DO n=1,bins
                IF(modk>=bin(n) .AND. modk<bin(n+1)) THEN
                   pow(n)=pow(n)+real(dk(i,j,m))**2.+aimag(dk(i,j,m))**2.
                   !                   ninbin(n)=ninbin(n)+1
                   EXIT
                END IF
             END DO

          END DO
       END DO
    END DO

    DEALLOCATE(dk)

    pow=pow/(float(m)**6.)

    nps=(float(np))**(1./3.)

    OPEN(19,file=fname)

    DO i=1,bins

       kv=(bin(i+1)+bin(i))/2.

       pow(i)=pow(i)/(log(bin(i+1)/bin(i)))
       !       shot=(kv**3.)/((2.*pi*pi)*(nps/L)**3.)
       shot=((L**3.)/float(np))*(kv**3.)/(2.*(pi**2.))

       WRITE(19,*) kv, pow(i), shot

    END DO

    DEALLOCATE(bin,pow)

    WRITE(*,*) 'Done computing power, written: ', fname
    WRITE(*,*)

  END SUBROUTINE power

  SUBROUTINE smooth(d,m,r,L)

    IMPLICIT NONE
    DOUBLE COMPLEX :: d(:,:,:)
    REAL :: r, L, kx, ky, kz, modk, kr
    INTEGER :: m, i, j, k

    kr=1./r

    WRITE(*,*) 'Smoothing density field'
    WRITE(*,*) 'Smoothing scale / (Mpc/h):', r
    WRITE(*,*) 'Smoothing k / ((Mpc/h)^-1):', kr

    CALL fft3(d,d,-1)

    DO i=1,m
       DO j=1,m
          DO k=1,m
             CALL k_fft(i,j,k,m,kx,ky,kz,modk,L)
             d(i,j,k)=d(i,j,k)*exp(-((modk*r)**2.)/2.)
          END DO
       END DO
    END DO

    CALL fft3(d,d,1)

    d=d/(float(m)**3.)

    WRITE(*,*) 'Smoothing complete'
    WRITE(*,*)

  END SUBROUTINE smooth

  SUBROUTINE displacement(d,fx,fy,fz,L)

    IMPLICIT NONE

    REAL :: z_thick, xpos, ypos, kx, ky, kz, modk, L, valx, valy, valz, val
    INTEGER :: i, j, k, m_thick, m
    DOUBLE COMPLEX, INTENT(IN) :: d(:,:,:)
    DOUBLE COMPLEX, INTENT(OUT) :: fx(:,:,:), fy(:,:,:), fz(:,:,:)
    DOUBLE COMPLEX, ALLOCATABLE :: ds(:,:,:)

    !This creates a displacement field from an over-density field in the framework of linear theory!

    m=SIZE(d(:,1,1))

    ALLOCATE(ds(m,m,m))

    ds=d

    fx=(0.d0,0.d0)
    fy=(0.d0,0.d0)
    fz=(0.d0,0.d0)

    WRITE(*,*) 'Creating linear displacement field from (should be linear) density field'

    CALL fft3(ds,ds,-1)

    DO i=1,m
       DO j=1,m
          DO k=1,m
             IF(i==1 .AND. j==1 .AND. k==1) THEN
                fx(1,1,1)=(0.d0,0.d0)
                fy(1,1,1)=(0.d0,0.d0)
                fz(1,1,1)=(0.d0,0.d0)
             ELSE
                CALL k_fft(i,j,k,m,kx,ky,kz,modk,L)
                fx(i,j,k)=ds(i,j,k)*(0.d0,-1.d0)*kx/(modk**2.)
                fy(i,j,k)=ds(i,j,k)*(0.d0,-1.d0)*ky/(modk**2.)
                fz(i,j,k)=ds(i,j,k)*(0.d0,-1.d0)*kz/(modk**2.)
             END IF
          END DO
       END DO
    END DO

    DEALLOCATE(ds)

    CALL fft3(fx,fx,1)
    CALL fft3(fy,fy,1)
    CALL fft3(fz,fz,1)

    !    Normalise after FFT
    fx=fx/(float(m)**3.)
    fy=fy/(float(m)**3.)
    fz=fz/(float(m)**3.)

  END SUBROUTINE displacement

!!$  SUBROUTINE write_field_slice(f,z_thick,L,filename)
!!$
!!$    IMPLICIT NONE
!!$
!!$    REAL :: z_thick, xpos, ypos, kx, ky, kz, modk, L, val
!!$    CHARACTER(len=256) :: filename
!!$    INTEGER :: i, j, k, m_thick, m
!!$    DOUBLE COMPLEX, INTENT(IN) :: f(:,:,:)
!!$
!!$    m=SIZE(f(:,1,1))
!!$
!!$    m_thick=CEILING(float(m)*z_thick/L)
!!$
!!$    WRITE(*,*) 'Writing out field'
!!$    WRITE(*,*) 'Thickness in z/(Mpc/h):', z_thick
!!$    WRITE(*,*) 'Tickness in cells:', m_thick
!!$
!!$    WRITE(*,*) 'Writing out field map'
!!$    OPEN(9,file=filename)
!!$    DO j=1,m
!!$       DO i=1,m
!!$
!!$          val=0.
!!$
!!$          DO k=1,m_thick
!!$             val=val+f(i,j,k)
!!$          END DO
!!$
!!$          val=val/float(m_thick)
!!$
!!$          xpos=L*(float(i)-0.5)/float(m)
!!$          ypos=L*(float(j)-0.5)/float(m)
!!$
!!$          WRITE(9,fmt='(3F20.5)') xpos, ypos, val
!!$
!!$       END DO
!!$    END DO
!!$    CLOSE(9)
!!$    WRITE(*,*) 'Field written:', filename
!!$    WRITE(*,*) ''
!!$
!!$  END SUBROUTINE write_field_slice

!!$  SUBROUTINE density_slice(x,x1,x2,y,y1,y2,z,z1,z2,L,s,bias,filename)
!!$
!!$    IMPLICIT NONE
!!$    REAL :: x1, x2, y1, y2, z1, z2, L, zt, s, spam, modk, kx, ky, kz, bias
!!$    CHARACTER(len=12) :: filename
!!$    REAL :: x(:), y(:), z(:), xpos, ypos, d2dr
!!$    INTEGER :: nx, ny, mc, i, j, n, k
!!$    DOUBLE COMPLEX, ALLOCATABLE :: d2d(:,:)
!!$
!!$    WRITE(*,*) 'Mesh size for density plot:'
!!$    READ(*,*) mc
!!$    WRITE(*,*)
!!$
!!$    !    mc=500
!!$    k=0
!!$
!!$    WRITE(*,*) 'Writing out density field'
!!$    WRITE(*,*) 'Thickness in z/(Mpc/h):', z2-z1
!!$
!!$    ALLOCATE(d2d(mc,mc))
!!$
!!$    d2d=0.
!!$
!!$    n=SIZE(x)
!!$
!!$    DO i=1,n
!!$       IF(x1<x(i) .AND. x(i)<=x2 .AND. y1<y(i) .AND. y(i)<=y2 .AND. z1<z(i) .AND. z(i)<=z2) THEN
!!$          nx=CEILING(x(i)*float(mc)/(x2-x1))
!!$          ny=CEILING(y(i)*float(mc)/(y2-y1))
!!$          IF(nx<1) nx=1
!!$          IF(nx>mc) nx=mc
!!$          IF(ny<1) ny=1
!!$          IF(ny>mc) ny=mc
!!$          d2d(nx,ny)=d2d(nx,ny)+1.
!!$       END IF
!!$    END DO
!!$
!!$    !    Normalise density by dividing by average!
!!$    d2d=d2d*(float(mc)**2.)/(float(n)*((z2-z1)/L))
!!$
!!$    !    Lognormal bias
!!$    d2d=d2d**(1./bias)
!!$
!!$    d2d=d2d/average2d(d2d)
!!$
!!$    !    Subtract 1
!!$    d2d=d2d-1.
!!$
!!$    call fft2(d2d,d2d,-1)
!!$
!!$    DO i=1,mc
!!$       DO j=1,mc
!!$
!!$          call k_fft(i,j,k,mc,kx,ky,kz,modk,L)
!!$
!!$          IF(j==1) THEN
!!$             WRITE(*,*) modk
!!$          END IF
!!$
!!$          d2d(i,j)=d2d(i,j)*exp(-((modk*s)**2.)/2.)
!!$
!!$       END DO
!!$    END DO
!!$
!!$    call fft2(d2d,d2d,1)
!!$
!!$    !Normalise after FFT
!!$    d2d=d2d/(float(mc)**2.)
!!$
!!$    WRITE(*,*) 'Writing map'
!!$    OPEN(9,file=filename)
!!$    WRITE(9,fmt='(128F10.5)') d2d
!!$    DO i=1,mc
!!$       DO j=1,mc
!!$          d2dr=real(d2d(i,j))
!!$          xpos=x1+(x2-x1)*(float(i)-0.5)/float(mc)
!!$          ypos=y1+(y2-y1)*(float(j)-0.5)/float(mc)
!!$          WRITE(9,*) xpos, ypos, d2dr
!!$       END DO
!!$    END DO
!!$    CLOSE(9)
!!$    WRITE(*,*) 'Density written:', filename
!!$    WRITE(*,*) ''
!!$
!!$    DEALLOCATE(d2d)
!!$
!!$  END SUBROUTINE density_slice

  SUBROUTINE k_fft(ix,iy,iz,mc,kx,ky,kz,kmod,L)

    IMPLICIT NONE
    INTEGER :: ix, iy, iz, mc
    REAL :: pi, kx, ky, kz, kmod, L

    pi=3.141592654

    kx=float(ix-1)
    ky=float(iy-1)
    kz=float(iz-1)

    IF(ix>mc/2+1) kx=-float(mc-ix+1)
    IF(iy>mc/2+1) ky=-float(mc-iy+1)
    IF(iz>mc/2+1) kz=-float(mc-iz+1)

    kx=kx*2.*pi/L
    ky=ky*2.*pi/L
    kz=kz*2.*pi/L

    kmod=sqrt((kx**2.)+(ky**2.)+(kz**2.))

  END SUBROUTINE k_fft

  !  FUNCTION ran(x1,x2)

  !    IMPLICIT NONE
  !    REAL :: rand, ran
  !    REAL :: x1,x2

  !rand is some inbuilt function!
  !ran=x1+(x2-x1)*(rand(0))

  !  END FUNCTION ran

  !  SUBROUTINE RNG_set

  !    IMPLICIT NONE
  !    INTEGER :: int, timearray(3)
  !    REAL :: rand

  !This fills the time array using the system clock!
  !If called within the same second the numbers will be identical!
  !    CALL itime(timeArray)
  !This then initialises the generator!
  !    int=rand(timeArray(1)+timeArray(2)+timeArray(3))

  !  END SUBROUTINE RNG_set

  FUNCTION integrate(a,b,f,acc)

    !Integrates between a and b until desired accuracy is reached!

    IMPLICIT NONE
    INTEGER :: i, j, jmax
    REAL :: integrate, a, b, acc, dx
    INTEGER :: nint
    REAL :: x, fac, this
    REAL*8 :: sum1, sum2

    INTERFACE
       REAL FUNCTION f(x)
         REAL, INTENT(IN) :: x
       END FUNCTION f
    END INTERFACE

    sum1=0.d0
    sum2=0.d0

    jmax=20

    DO j=1,jmax

       nint=10.*(2.**j)

       DO i=1,nint

          x=a+(b-a)*((float(i)-1)/(float(nint)-1))

          IF(i==1 .OR. i==nint) THEN
             !multiple of 1 for beginning and end and multiple of 2 for middle points!
             fac=1
          ELSE
             fac=2
          END IF

          sum2=sum2+fac*f(x)

       END DO

       dx=((b-a)/(float(nint)-1.))
       sum2=sum2*dx/2.

       IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
          integrate=sum2
          !       WRITE(*,*) nint
          EXIT
       ELSE IF(j==jmax) THEN
          WRITE(*,*) 'Integration timed out'
       ELSE
          sum1=sum2
          sum2=0.d0
       END IF

    END DO

  END FUNCTION integrate

  FUNCTION integrate_vard(a,b,r,z,cosm,acc)

    !Integrates between a and b until desired accuracy is reached!

    USE cosdef
    IMPLICIT NONE
    INTEGER :: i, j, jmax
    REAL :: integrate_vard, a, b, acc, dx
    INTEGER :: nint
    REAL :: fac, this, r, z, fun, x
    REAL*8 :: sum1, sum2
    TYPE(cosmology) :: cosm

    sum1=0.d0
    sum2=0.d0

    jmax=20

    DO j=1,jmax

       nint=10.*(2.**j)

       DO i=1,nint

          x=a+(b-a)*((float(i)-1)/(float(nint)-1))

          IF(i==1 .OR. i==nint) THEN
             !multiple of 1 for beginning and end and multiple of 2 for middle points!
             fac=1.
          ELSE
             fac=2.
          END IF

          fun=exp(-((x*r)**2.))*plin(x,z,cosm)/x
          sum2=sum2+fac*fun

          !          WRITE(*,*) x, r, fac, plin(x,z,cosm)

       END DO

       dx=(b-a)/float(nint-1)
       sum2=sum2*dx/2.

       IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
          integrate_vard=sum2
          !       WRITE(*,*) nint
          EXIT
       ELSE IF(j==jmax) THEN
          WRITE(*,*) 'Integration timed out'
       ELSE
          sum1=sum2
          sum2=0.d0
       END IF

    END DO

  END FUNCTION integrate_vard

  FUNCTION integrate_varp(a,b,r,z,cosm,acc)

    !Integrates between a and b until desired accuracy is reached!

    USE cosdef

    IMPLICIT NONE
    INTEGER :: i, j, jmax
    REAL :: integrate_varp, a, b, acc, dx
    INTEGER :: nint
    REAL :: fac, this, r, z, fun, x
    REAL*8 :: sum1, sum2
    TYPE(cosmology) :: cosm

    sum1=0.d0
    sum2=0.d0

    jmax=20

    DO j=1,jmax

       nint=10.*(2.**j)

       DO i=1,nint

          x=a+(b-a)*((float(i)-1)/(float(nint)-1))

          IF(i==1 .OR. i==nint) THEN
             !multiple of 1 for beginning and end and multiple of 2 for middle points!
             fac=1.
          ELSE
             fac=2.
          END IF

          fun=exp(-((x*r)**2.))*plin(x,z,cosm)/(x**3.)
          sum2=sum2+fac*fun

       END DO

       dx=(b-a)/float(nint-1)
       sum2=sum2*dx/2.

       IF(j .NE. 1 .AND. ABS(-1.+sum2/sum1)<acc) THEN
          integrate_varp=sum2
          !       WRITE(*,*) nint
          EXIT
       ELSE IF(j==jmax) THEN
          WRITE(*,*) 'Integration timed out'
       ELSE
          sum1=sum2
          sum2=0.d0
       END IF

    END DO

  END FUNCTION integrate_varp

  SUBROUTINE assign_cosmology(cosm,z)

    USE cosdef
    IMPLICIT NONE
    INTEGER :: iopt
    REAL :: z
    TYPE(cosmology) :: cosm

    WRITE(*,*) 'Cosmology:'
    WRITE(*,*) '0 - Assign own'
    WRITE(*,*) '1 - Millennium (WMAP1ish)'
    WRITE(*,*) '2 - Multidark (WMAP5ish)'
    WRITE(*,*) '3 - CoDECS (WMAP7ish)'
    WRITE(*,*) '4 - TCDM (defw)'
    WRITE(*,*) '5 - Deus - LCDM'
    WRITE(*,*) '6 - Deus - RPCDM'
    WRITE(*,*) '7 - Deus - SUCDM'
    WRITE(*,*) '8 - TCDM (defw, high sigma 8)'
    WRITE(*,*) '9 - EdS (defw, high sigma 8)'
    WRITE(*,*) '10 - GR (WMAP9)'
    WRITE(*,*) '11 - F4'
    WRITE(*,*) '12 - F5'
    WRITE(*,*) '13 - F6'
    WRITE(*,*) '14 - Vanilla (high sigma 8)'
    READ(*,*) cosm%icos
    WRITE(*,*)

    cosm%ipow=1
    cosm%igro=1
    cosm%idc=1
    cosm%irate=1
    cosm%icham=0

    IF(cosm%icos==0) THEN

       WRITE(*,*) '1 - EH approximation'
       WRITE(*,*) '2 - DEFW'
       READ(*,*) cosm%itk
       WRITE(*,*)

       WRITE(*,*) 'Omega_m:'
       READ(*,*) cosm%om_m
       WRITE(*,*) 'Omega_v:'
       READ(*,*) cosm%om_v
       WRITE(*,*) 'h:'
       READ(*,*) cosm%h
       WRITE(*,*) 'sig8:'
       READ(*,*) cosm%sig8
       WRITE(*,*) 'n:'
       READ(*,*) cosm%n

       IF(cosm%itk==1) THEN
          WRITE(*,*) 'Omega_b:'
          READ(*,*) cosm%om_b
       ELSE IF(cosm%itk==2) THEN
          WRITE(*,*) 'Gamma:'
          READ(*,*) cosm%gamma
       END IF

       WRITE(*,*)

    ELSE IF(cosm%icos==1) THEN

       WRITE(*,*) 'Millennium'

       cosm%om_m=0.25
       cosm%om_b=0.045
       cosm%om_v=0.75
       cosm%h=0.73
       cosm%sig8=0.9
       cosm%n=1.0
       cosm%itk=1

       WRITE(*,*) '1 - EH approximation'
       WRITE(*,*) '2 - Actual Millennium Tk'
       READ(*,*) i
       WRITE(*,*)

       IF(i==2) cosm%ipow=0

    ELSE IF(cosm%icos==2) THEN

       WRITE(*,*) 'Multidark'
       WRITE(*,*)
       STOP 'Not fully supported (edit input tables and files)'

       cosm%om_m=0.27
       cosm%om_b=0.0469
       cosm%om_v=0.73
       cosm%h=0.70
       cosm%sig8=0.82
       cosm%n=0.95
       cosm%itk=1

    ELSE IF(cosm%icos==3) THEN

       WRITE(*,*) 'CoDECS'
       WRITE(*,*)
       STOP 'Not fully supported (edit input tables and files)'

       cosm%om_m=0.271
       cosm%om_b=0.0451
       cosm%om_v=0.729
       cosm%h=0.703
       cosm%sig8=0.809
       cosm%n=0.966
       cosm%itk=1

    ELSE IF(cosm%icos==4) THEN

       WRITE(*,*) 'TCDM'
       WRITE(*,*)

       cosm%itk=2
       cosm%om_m=1.
       cosm%om_v=0.
       cosm%h=0.5
       cosm%sig8=0.51
       cosm%n=1.
       cosm%gamma=0.21

       !       cosm%dc=2.

    ELSE IF(cosm%icos==5) THEN

       WRITE(*,*) 'Deus - LCDM'
       WRITE(*,*)
       STOP 'Not fully supported (edit input tables)'

       cosm%itk=1
       cosm%n=0.96
       cosm%om_b=0.04
       cosm%h=0.72

       cosm%om_v=0.74
       cosm%om_m=0.26
       cosm%sig8=0.79

    ELSE IF(cosm%icos==6) THEN

       WRITE(*,*) 'Deus - RPCDM'
       WRITE(*,*)
       STOP 'Not fully supported (edit input tables)'

       cosm%itk=1
       cosm%n=0.96
       cosm%om_b=0.04
       cosm%h=0.72

       cosm%om_v=0.77
       cosm%om_m=0.23
       cosm%sig8=0.66

    ELSE IF(cosm%icos==7) THEN

       WRITE(*,*) 'Deus - SUCDM'
       WRITE(*,*)
       STOP 'Not fully supported (edit input tables)'

       cosm%itk=1
       cosm%n=0.96
       cosm%om_b=0.04
       cosm%h=0.72

       cosm%om_v=0.75
       cosm%om_m=0.25
       cosm%sig8=0.73

    ELSE IF(cosm%icos==8) THEN

       WRITE(*,*) 'TCDM (high sigma 8)'
       WRITE(*,*)

       cosm%itk=2
       cosm%om_m=1.
       cosm%om_v=0.
       cosm%h=0.5
       cosm%sig8=0.8
       cosm%n=1.
       cosm%gamma=0.21

    ELSE IF(cosm%icos==9) THEN

       WRITE(*,*) 'EdS (high sigma 8)'
       WRITE(*,*)

       cosm%itk=2
       cosm%om_m=1.
       cosm%om_v=0.
       cosm%h=0.5
       cosm%sig8=0.8
       cosm%n=1.
       cosm%gamma=0.5

    ELSE IF(cosm%icos==10) THEN

       WRITE(*,*) 'WMAP9'
       WRITE(*,*)

       cosm%om_m=0.281
       cosm%om_b=0.046
       cosm%om_v=0.719
       cosm%h=0.697
       cosm%sig8=0.82
       cosm%n=0.971
       cosm%itk=1

       cosm%ipow=0

    ELSE IF(cosm%icos==11) THEN

       WRITE(*,*) 'F4'
       WRITE(*,*) 
       WRITE(*,*) '1 - Use fixed dc and scale dependent growth'
       WRITE(*,*) '2 - Use effective dc and LCDM growth'
       READ(*,*) iopt
       WRITE(*,*)

       cosm%om_m=0.281
       cosm%om_b=0.046
       cosm%om_v=0.719
       cosm%h=0.697
       cosm%sig8=0.82
       cosm%n=0.971

       cosm%fr0=-1e-4
       cosm%frn=1.

       cosm%itk=1
       cosm%icham=1

       IF(iopt==2) cosm%idc=0
       cosm%ipow=0
       cosm%irate=2

    ELSE IF(cosm%icos==12) THEN

       WRITE(*,*) 'F5'
       WRITE(*,*) 
       WRITE(*,*) '1 - Use fixed dc and scale dependent growth'
       WRITE(*,*) '2 - Use effective dc and LCDM growth'
       READ(*,*) iopt
       WRITE(*,*)

       cosm%om_m=0.281
       cosm%om_b=0.046
       cosm%om_v=0.719
       cosm%h=0.697
       cosm%sig8=0.82
       cosm%n=0.971

       cosm%fr0=-1e-5
       cosm%frn=1.

       cosm%itk=1
       cosm%icham=1

       IF(iopt==2) cosm%idc=0
       cosm%ipow=0
       cosm%irate=2

    ELSE IF(cosm%icos==13) THEN

       WRITE(*,*) 'F6'
       WRITE(*,*) 
       WRITE(*,*) '1 - Use fixed dc and scale dependent growth'
       WRITE(*,*) '2 - Use effective dc and LCDM growth'
       READ(*,*) iopt
       WRITE(*,*)

       cosm%om_m=0.281
       cosm%om_b=0.046
       cosm%om_v=0.719
       cosm%h=0.697
       cosm%sig8=0.82
       cosm%n=0.971

       cosm%fr0=-1e-6
       cosm%frn=1.

       cosm%itk=1
       cosm%icham=1

       IF(iopt==2) cosm%idc=0
       cosm%ipow=0
       cosm%irate=2

    ELSE IF(cosm%icos==14) THEN

       WRITE(*,*) 'Vanilla (high sigma 8)'
       WRITE(*,*)

       cosm%om_m=0.3
       cosm%om_b=0.05
       cosm%om_v=0.7
       cosm%h=0.7
       cosm%sig8=1.2
       cosm%n=0.97
       cosm%itk=1

       cosm%ipow=0

    END IF

    IF(cosm%icham==1) THEN
       WRITE(*,*) '0 - Ignore enhanced gravity'
       WRITE(*,*) '1 - Include enhancement and screening'
       WRITE(*,*) '2 - Enhance gravity with no screening whatsoever'
       READ(*,*) cosm%icham
       WRITE(*,*)
    END IF

    WRITE(*,*) 'Your chosen cosmology is:'
    WRITE(*,*) 'Om_m', cosm%om_m
    WRITE(*,*) 'Om_v', cosm%om_v
    IF(cosm%itk==1) WRITE(*,*) 'Omega_b:', cosm%om_b
    WRITE(*,*) 'h', cosm%h
    WRITE(*,*) 'sig8', cosm%sig8
    WRITE(*,*) 'n', cosm%n
    IF(cosm%itk==2) WRITE(*,*) 'Gamma:', cosm%gamma
    WRITE(*,*)

    WRITE(*,*) 'Redshift'
    READ(*,*) z
    WRITE(*,*)

    IF(cosm%ipow==0) CALL input_pow_tab(cosm)
    IF(cosm%igro==0 .OR. cosm%igro==2)   CALL input_gro_tab(cosm)
    IF(cosm%irate==0 .OR. cosm%irate==2) CALL input_rate_tab(cosm)
    IF(cosm%idc==0) CALL input_dc_tab(cosm)

    CALL normalise(cosm)

    !Call this after normalisation so that code does GR normalisation!
    IF(iopt==1 .AND. (cosm%icos==11 .OR. cosm%icos==12 .OR. cosm%icos==13)) THEN
       cosm%igro=2
       CALL input_gro_tab(cosm)
    END IF

  END SUBROUTINE assign_cosmology

!!$  SUBROUTINE input_dc_tab(cosm)
!!$
!!$    !Inputs a dc(M) tab for a modified gravity model
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    TYPE(cosmology) :: cosm
!!$    CHARACTER(len=256) :: input
!!$    REAL :: spam
!!$    INTEGER :: i, n
!!$
!!$    IF(cosm%icos==11 .OR. cosm%icos==12 .OR. cosm%icos==13) THEN
!!$
!!$       IF(cosm%icos==11) input='/disk1/am/simulations/mg/dc/deltac_F4.dat'
!!$       IF(cosm%icos==12) input='/disk1/am/simulations/mg/dc/deltac_F5.dat'
!!$       IF(cosm%icos==13) input='/disk1/am/simulations/mg/dc/deltac_F6.dat'
!!$
!!$       !       This does not work for some reason (file header?)
!!$       !       n=file_length(input)
!!$       !       WRITE(*,*) 'n:', n
!!$       !       STOP
!!$
!!$       n=100
!!$
!!$       ALLOCATE(cosm%m(n),cosm%dcol(n))
!!$
!!$       OPEN(7,file=input)
!!$       DO i=0,n
!!$          IF(i==0) THEN
!!$             READ(7,*)
!!$          ELSE
!!$             READ(7,*) spam, cosm%m(i), spam, spam, spam, spam, spam, cosm%dcol(i)
!!$          END IF
!!$       END DO
!!$       CLOSE(7)
!!$
!!$    END IF
!!$
!!$  END SUBROUTINE input_dc_tab

!!$  SUBROUTINE input_pow_tab(cosm)
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    INTEGER :: n, i
!!$    CHARACTER(len=256) :: filebase
!!$    TYPE(cosmology) :: cosm
!!$    REAL :: spam
!!$    LOGICAL :: lexist
!!$
!!$!    WRITE(*,*) '1 - Millennium Tk'
!!$!    WRITE(*,*) '2 - Codecs'
!!$!    WRITE(*,*) '3 - Deus'
!!$!    WRITE(*,*) '4 - QUICC'
!!$!    READ(*,*) ic
!!$!    WRITE(*,*)
!!$
!!$    IF(cosm%icos==1) THEN
!!$
!!$       WRITE(*,*) 'Reading in Millennium Transfer function'
!!$
!!$       filebase='/disk1/am/millennium/millk.dat'
!!$
!!$       n=file_length(filebase)
!!$
!!$       ALLOCATE(cosm%k(n),cosm%p(n))
!!$
!!$       OPEN(7,file=filebase)
!!$       DO i=1,n
!!$          READ(7,*) cosm%k(i), cosm%p(i)
!!$       END DO
!!$       CLOSE(7)
!!$
!!$       !Convert kpc -> Mpc and log scale!
!!$       cosm%k=1000.*(10.**cosm%k)
!!$       cosm%p=10.**cosm%p
!!$
!!$       WRITE(*,*) 'Done'
!!$       WRITE(*,*)
!!$
!!$    ELSE IF(cosm%icos==10 .OR. cosm%icos==11 .OR. cosm%icos==12 .OR. cosm%icos==13 .OR. cosm%icos==14) THEN
!!$   
!!$       IF(cosm%icos==14) THEN
!!$!          filebase='/disk1/am/simulations/mg/linear_power/camb_602.dat'
!!$          filebase='/disk1/am/simulations/mg/linear_power/power_Li_LCDM.dat'
!!$       ELSE IF(cosm%icos==10 .OR. cosm%icos==11 .OR. cosm%icos==12 .OR. cosm%icos==13) THEN
!!$!          filebase='/disk1/am/simulations/mg/linear_power/camb_gr.dat'
!!$          filebase='/disk1/am/simulations/mg/linear_power/power_Li_GR.dat'
!!$       END IF
!!$       
!!$       INQUIRE(file=filebase,exist=lexist)
!!$       IF(lexist==.FALSE.) STOP 'You cock, this does not exist'
!!$       
!!$       n=file_length(filebase)
!!$
!!$       WRITE(*,*) 'Reading in power spectra:', filebase
!!$       WRITE(*,*) 'Length:', n
!!$
!!$       ALLOCATE(cosm%k(n),cosm%p(n))
!!$        
!!$       OPEN(7,file=filebase)
!!$       DO i=1,n
!!$
!!$          !For CAMB spectra
!!$!          READ(7,*) cosm%k(i), cosm%p(i)
!!$!          cosm%p(i)=cosm%p(i)*((cosm%k(i))**3.)/19.7392088
!!$
!!$          !For Baojiu spectra
!!$          READ(7,*) cosm%k(i), spam, cosm%p(i)
!!$
!!$       END DO
!!$       CLOSE(7)
!!$
!!$
!!$       WRITE(*,*) 'Done'
!!$       WRITE(*,*)
!!$
!!$    END IF
!!$
!!$  END SUBROUTINE input_pow_tab

!!$  SUBROUTINE input_gro_tab(cosm)
!!$    
!!$    !Inputs a table of a scale-dependent growth factor
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    TYPE(cosmology) :: cosm
!!$    CHARACTER(len=256) :: input
!!$    REAL :: spam
!!$    INTEGER :: i, n
!!$
!!$    IF(cosm%icos==11 .OR. cosm%icos==12 .OR. cosm%icos==13) THEN
!!$
!!$       input='/disk1/am/simulations/mg/growth/growth_z0.0.txt'
!!$       n=501
!!$
!!$       WRITE(*,*) 'Reading in growth factor table:', input
!!$
!!$       ALLOCATE(cosm%k_g(n),cosm%gk(n))
!!$
!!$       OPEN(7,file=input)
!!$       DO i=1,n
!!$          IF(cosm%icos==11) READ(7,*) cosm%k_g(i), spam, spam, spam, spam, spam, spam, spam, cosm%gk(i)
!!$          IF(cosm%icos==12) READ(7,*) cosm%k_g(i), spam, spam, spam, spam, spam, spam, cosm%gk(i)
!!$          IF(cosm%icos==13) READ(7,*) cosm%k_g(i), spam, spam, spam, spam, spam, cosm%gk(i)
!!$       END DO
!!$       CLOSE(7)
!!$
!!$       WRITE(*,*) 'Done'
!!$       WRITE(*,*)
!!$
!!$    END IF
!!$
!!$  END SUBROUTINE input_gro_tab

!!$  SUBROUTINE input_rate_tab(cosm)
!!$   
!!$    !Inputs a table of a scale-dependent growth rate
!!$ 
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    TYPE(cosmology) :: cosm
!!$    CHARACTER(len=256) :: input
!!$    REAL :: spam
!!$    INTEGER :: i, n
!!$
!!$    IF(cosm%icos==11 .OR. cosm%icos==12 .OR. cosm%icos==13) THEN
!!$
!!$       input='/disk1/am/simulations/mg/growth/growth_z0.0.txt'
!!$       n=501
!!$
!!$       WRITE(*,*) 'Reading in growth rate table:', input
!!$
!!$       ALLOCATE(cosm%k_r(n),cosm%rk(n))
!!$
!!$       OPEN(7,file=input)
!!$       DO i=1,n
!!$          IF(cosm%icos==11) READ(7,*) cosm%k_r(i), spam, spam, spam, cosm%rk(i)
!!$          IF(cosm%icos==12) READ(7,*) cosm%k_r(i), spam, spam, cosm%rk(i)
!!$          IF(cosm%icos==13) READ(7,*) cosm%k_r(i), spam, cosm%rk(i)
!!$!          WRITE(*,*) cosm%k_r(i), cosm%rk(i)
!!$       END DO
!!$       CLOSE(7)
!!$
!!$        WRITE(*,*) 'Done'
!!$        WRITE(*,*)
!!$
!!$    END IF
!!$
!!$  END SUBROUTINE input_rate_tab

  FUNCTION find(x,xtab,ytab,iorder,imeth)

    IMPLICIT NONE
    REAL :: find
    REAL, INTENT(IN) :: x, xtab(:), ytab(:)
    REAL :: a, b, c, d
    REAL :: x1, x2, x3, x4
    REAL :: y1, y2, y3, y4
    INTEGER :: i, n
    INTEGER, INTENT(IN) :: imeth, iorder
    INTEGER :: maxorder, maxmethod

    !This version interpolates if the value is off either end of the array!
    !Care should be chosen to insert x, xtab, ytab as log if this might give better!
    !Results from the interpolation!

    !If the value required is off the table edge the interpolation is always linear

    !imeth = 1 => find x in xtab by crudely searching from x(1) to x(n)
    !imeth = 2 => find x in xtab quickly assuming the table is linearly spaced
    !imeth = 3 => find x in xtab using midpoint splitting (iterations=CEILING(log2(n)))

    !iorder = 1 => linear interpolation
    !iorder = 2 => quadratic interpolation
    !iorder = 3 => cubic interpolation

    n=SIZE(xtab)

    maxorder=3
    maxmethod=3

    IF(xtab(1)>xtab(n)) STOP 'FIND: x table in wrong order'
    IF(n .NE. SIZE(ytab)) STOP 'FIND: Tables not of the same size'
    IF(iorder<1) STOP 'FIND: find order not specified correctly'
    IF(iorder>maxorder) STOP 'FIND: find order not specified correctly'
    IF(imeth<1) STOP 'FIND: Method of finding within a table not specified correctly'
    IF(imeth>maxmethod) STOP 'FIND: Method of finding within a table not specified correctly'

    IF(x<xtab(1)) THEN

       x1=xtab(1)
       x2=xtab(2)

       y1=ytab(1)
       y2=ytab(2)

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b

    ELSE IF(x>xtab(n)) THEN

       x1=xtab(n-1)
       x2=xtab(n)

       y1=ytab(n-1)
       y2=ytab(n)

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b

    ELSE IF(iorder==1) THEN

       IF(n<2) STOP 'FIND: Not enough points in your table for linear interpolation'

       IF(x<=xtab(2)) THEN

          x2=xtab(2)
          x1=xtab(1)

          y2=ytab(2)
          y1=ytab(1)

       ELSE IF (x>=xtab(n-1)) THEN

          x2=xtab(n)
          x1=xtab(n-1)

          y2=ytab(n)
          y1=ytab(n-1)

       ELSE

          IF(imeth==1) i=search_int(x,xtab)
          IF(imeth==2) i=linear_table_integer(x,xtab)
          IF(imeth==3) i=int_split(x,xtab)

          x2=xtab(i+1)
          x1=xtab(i)

          y2=ytab(i+1)
          y1=ytab(i)

       END IF

       CALL fit_line(a,b,x1,y1,x2,y2)
       find=a*x+b

    ELSE IF(iorder==2) THEN

       IF(n<3) STOP 'FIND_QUADRATIC: Not enough points in your table'

       IF(x<=xtab(2) .OR. x>=xtab(n-1)) THEN

          IF(x<=xtab(2)) THEN

             x3=xtab(3)
             x2=xtab(2)
             x1=xtab(1)

             y3=ytab(3)
             y2=ytab(2)
             y1=ytab(1)

          ELSE IF (x>=xtab(n-1)) THEN

             x3=xtab(n)
             x2=xtab(n-1)
             x1=xtab(n-2)

             y3=ytab(n)
             y2=ytab(n-1)
             y1=ytab(n-2)

          END IF

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)

          find=a*(x**2.)+b*x+c

       ELSE

          IF(imeth==1) i=search_int(x,xtab)
          IF(imeth==2) i=linear_table_integer(x,xtab)
          IF(imeth==3) i=int_split(x,xtab)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

          !In this case take the average of two separate quadratic spline values

          find=0.

          CALL fit_quadratic(a,b,c,x1,y1,x2,y2,x3,y3)
          find=find+(a*(x**2.)+b*x+c)/2.

          CALL fit_quadratic(a,b,c,x2,y2,x3,y3,x4,y4)
          find=find+(a*(x**2.)+b*x+c)/2.

       END IF

    ELSE IF(iorder==3) THEN

       IF(n<4) STOP 'FIND_CUBIC: Not enough points in your table'

       IF(x<=xtab(3)) THEN

          x4=xtab(4)
          x3=xtab(3)
          x2=xtab(2)
          x1=xtab(1)

          y4=ytab(4)
          y3=ytab(3)
          y2=ytab(2)
          y1=ytab(1)

       ELSE IF (x>=xtab(n-2)) THEN

          x4=xtab(n)
          x3=xtab(n-1)
          x2=xtab(n-2)
          x1=xtab(n-3)

          y4=ytab(n)
          y3=ytab(n-1)
          y2=ytab(n-2)
          y1=ytab(n-3)

       ELSE

          IF(imeth==1) i=search_int(x,xtab)
          IF(imeth==2) i=linear_table_integer(x,xtab)
          IF(imeth==3) i=int_split(x,xtab)

          x1=xtab(i-1)
          x2=xtab(i)
          x3=xtab(i+1)
          x4=xtab(i+2)

          y1=ytab(i-1)
          y2=ytab(i)
          y3=ytab(i+1)
          y4=ytab(i+2)

       END IF

       CALL fit_cubic(a,b,c,d,x1,y1,x2,y2,x3,y3,x4,y4)
       find=a*x**3.+b*x**2.+c*x+d

    END IF

  END FUNCTION find

  FUNCTION linear_table_integer(x,xtab)

    IMPLICIT NONE
    INTEGER :: linear_table_integer
    REAL, INTENT(IN) :: x, xtab(:)
    INTEGER :: n
    REAL :: x1, x2, xn
    REAL :: acc

    !Returns the integer (table position) below the value of x
    !eg. if x(3)=6. and x(4)=7. and x=6.5 this will return 6
    !Assumes table is organised linearly (care for logs)

    n=SIZE(xtab)
    x1=xtab(1)
    x2=xtab(2)
    xn=xtab(n)

    !Test for linear table
    acc=0.001

    IF(x1>xn) STOP 'LINEAR_TABLE_INTEGER :: table in the wrong order'
    IF(ABS(-1.+float(n-1)*(x2-x1)/(xn-x1))>acc) STOP 'LINEAR_TABLE_INTEGER :: table does not seem to be linear'

    linear_table_integer=1+FLOOR(float(n-1)*(x-x1)/(xn-x1))

  END FUNCTION linear_table_integer

  FUNCTION search_int(x,xtab)

    IMPLICIT NONE
    INTEGER :: search_int
    INTEGER :: i, n
    REAL, INTENT(IN) :: x, xtab(:)

    n=SIZE(xtab)

    IF(xtab(1)>xtab(n)) STOP 'SEARCH_INT: table in wrong order'

    DO i=1,n
       IF(x>=xtab(i) .AND. x<=xtab(i+1)) EXIT
    END DO

    search_int=i

  END FUNCTION search_int

  FUNCTION int_split(x,xtab)

    IMPLICIT NONE
    REAL, INTENT(IN) :: x, xtab(:)
    INTEGER :: i1, i2, imid, n
    INTEGER :: int_split

    !Finds the position of the value in the table by continually splitting it in half

    n=SIZE(xtab)

    IF(xtab(1)>xtab(n)) STOP 'INT_SPLIT: table in wrong order'

    i1=1
    i2=n

    DO

       imid=NINT((i1+i2)/2.)

       IF(x<xtab(imid)) THEN
          i2=imid
       ELSE
          i1=imid
       END IF

       IF(i2==i1+1) EXIT

    END DO

    int_split=i1

  END FUNCTION int_split

  SUBROUTINE reverse(arry)

    IMPLICIT NONE
    INTEGER :: n, i
    REAL, ALLOCATABLE :: hold(:)
    REAL :: arry(:)

    !This reverses the contents of arry!

    n=SIZE(arry)

    ALLOCATE(hold(n))

    hold=arry

    DO i=1,n
       arry(i)=hold(n-i+1)
    END DO

    DEALLOCATE(hold)

  END SUBROUTINE reverse

!!$  SUBROUTINE read_particles(cosm,zz,L,n,iph)
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    REAL :: zz, L
!!$    INTEGER :: n
!!$    INTEGER :: i, icut, n_tot, ii, ims, iph
!!$    REAL :: cut, spam, mas, xxx, yyy, zzz, mmm
!!$    CHARACTER(len=256) :: fname
!!$    LOGICAL :: lexist
!!$    TYPE(cosmology) :: cosm
!!$
!!$    !    WRITE(*,*) 'Input positions (should be that from original cosmology):'
!!$    !    WRITE(*,*) '1 - Multidark haloes'
!!$    !    WRITE(*,*) '2 - CoDECS haloes'
!!$    !    WRITE(*,*) '3 - Deus haloes'
!!$    !    WRITE(*,*) '4 - Millennium haloes'
!!$    !    WRITE(*,*) '5 - Gadget2 data file'
!!$    !    WRITE(*,*) '9 - Use previous (unformatted)'
!!$    !    READ(*,*) icat
!!$
!!$    WRITE(*,*) '0 - Particles'
!!$    WRITE(*,*) '1 - Haloes'
!!$    !    WRITE(*,*) '3 - My halo catalogues'
!!$    READ(*,*) iph
!!$    WRITE(*,*)
!!$
!!$    IF(iph==0) THEN
!!$
!!$       CALL read_gadget(cosm,zz,L,n)
!!$
!!$    ELSE IF(iph==1) THEN
!!$
!!$       CALL read_gadget_catalogue(cosm,zz,L,n)
!!$!       CALL mass_cut(n)
!!$
!!$    END IF
!!$
!!$!    ELSE IF(iph==2 .AND. (icos==5 .OR. icos==6 .OR. icos==7)) THEN
!!$
!!$!       CALL deus_name(icos,zz,L,fname)
!!$
!!$!       INQUIRE(file=fname,exist=lexist)
!!$!       IF(lexist==.FALSE.) STOP 'You cock, this halo catalogue does not exist'
!!$
!!$!       n=file_length(fname)
!!$
!!$!       CALL array_alloc(n)
!!$
!!$!       OPEN(7,file=fname)
!!$!       DO i=1,n
!!$!          READ(7,*) x(i), y(i), z(i), m(i)
!!$!       END DO
!!$!       CLOSE(7)
!!$
!!$!       vx=0.
!!$!       vy=0.
!!$!       vz=0.
!!$
!!$!    ELSE IF(iph==3) THEN
!!$
!!$!       IF(icos==1) THEN
!!$          !          WRITE(*,*) '1 - Read MS halo catalogue'
!!$          !          WRITE(*,*) '2 - Read my halo catalogue'
!!$          !          READ(*,*) ims
!!$          !          WRITE(*,*)
!!$!          CALL read_mill(L)
!!$          !          IF(ims==2) CALL my_cat_name(icos,zz,L,fname)
!!$!       END IF
!!$
!!$!       IF(icos==2) CALL read_mdfile
!!$!       IF(icos==3) CALL read_codecs
!!$!       IF(icos==5 .OR. icos==6 .OR. icos==7) CALL deus_name(icos,zz,L,fname)
!!$!       IF(icos==9) CALL read_unformatted('inpdat')
!!$!       IF(icos==4) CALL my_cat_name(icos,zz,L,fname)
!!$!       IF(icos==8) STOP 'QUICC catalgoues do not exist yet'
!!$
!!$!       inquire(file=fname,exist=lexist)
!!$
!!$!       IF(lexist==.FALSE.) STOP 'You cock, this halo catalogue does not exist'
!!$
!!$!       n=0
!!$!       n_tot=0
!!$
!!$       !cut=0.
!!$!       icut=0
!!$
!!$!       WRITE(*,*) 'Mass cut:'
!!$!       WRITE(*,*) '0 - No cut'
!!$!       WRITE(*,*) '1 - Greater than'
!!$!       WRITE(*,*) '2 - Less than'
!!$!       READ(*,*) icut
!!$
!!$!       IF(icut==1 .OR. icut==2) THEN
!!$!          WRITE(*,*) 'Mass:'
!!$!          READ(*,*) cut
!!$!       END IF
!!$!       WRITE(*,*)
!!$
!!$!       OPEN(7,file=fname)
!!$!       DO
!!$!          READ(7,*, end=401) spam, spam, spam, mas
!!$!          n_tot=n_tot+1
!!$!          IF(icut==1 .AND. mas>cut) n=n+1
!!$!          IF(icut==2 .AND. mas<cut) n=n+1
!!$!       END DO
!!$!401    CLOSE(7)  
!!$
!!$!       IF(icut==0) n=n_tot
!!$
!!$!       IF(icut==0) WRITE(*,*) 'No mass cut'
!!$!       IF(icut==1) WRITE(*,*) 'Only haloes above:', cut
!!$!       IF(icut==2) WRITE(*,*) 'Only haloes below:', cut
!!$!       WRITE(*,*) 'Total particles in file:', n_tot
!!$!       WRITE(*,*) 'Paricles in cut:', n
!!$
!!$!       CALL array_alloc(n)
!!$
!!$!       ii=0
!!$
!!$!       OPEN(7,file=fname)
!!$!       DO i=1,n_tot
!!$!          READ(7,*) xxx, yyy, zzz, mmm
!!$!          IF(icut==0 .OR. (icut==1 .AND. mmm>cut) .OR. (icut==2 .AND. mmm<cut)) THEN
!!$!             ii=ii+1
!!$!             x(ii)=xxx
!!$!             y(ii)=yyy
!!$!             z(ii)=zzz
!!$!             m(ii)=mmm
!!$             !          WRITE(*,*) ii
!!$!          END IF
!!$!       END DO
!!$!       CLOSE(7)
!!$
!!$!       WRITE(*,*) 'Finished reading halo catalogue'
!!$!       WRITE(*,*)
!!$
!!$!    END IF
!!$
!!$  END SUBROUTINE read_particles

!!$  SUBROUTINE vfield(vx,vy,vz,nx,ny,nz,fvx,fvy,fvz,r,L)
!!$
!!$    IMPLICIT NONE
!!$    REAL, INTENT(IN) :: vx(:), vy(:), vz(:), r
!!$    REAL :: kx, ky, kz, kmod, L
!!$    DOUBLE COMPLEX, INTENT(OUT) :: fvx(:,:,:), fvy(:,:,:), fvz(:,:,:)
!!$    INTEGER, INTENT(IN) :: nx(:), ny(:), nz(:)
!!$    INTEGER :: n, m, i, j, k
!!$    INTEGER, ALLOCATABLE :: ncell(:,:,:)
!!$
!!$    !First build velocity field!
!!$
!!$    n=SIZE(vx)
!!$    m=SIZE(fvx(:,1,1))
!!$
!!$    ALLOCATE(ncell(m,m,m))
!!$
!!$    fvx=0.
!!$    fvy=0.
!!$    fvz=0.
!!$
!!$    ncell=0
!!$
!!$    DO i=1,n
!!$       fvx(nx(i),ny(i),nz(i))=fvx(nx(i),ny(i),nz(i))+vx(i)
!!$       fvy(nx(i),ny(i),nz(i))=fvy(nx(i),ny(i),nz(i))+vy(i)
!!$       fvz(nx(i),ny(i),nz(i))=fvz(nx(i),ny(i),nz(i))+vz(i)
!!$       ncell(nx(i),ny(i),nz(i))=ncell(nx(i),ny(i),nz(i))+1
!!$    END DO
!!$
!!$    DO k=1,m
!!$       DO j=1,m
!!$          DO i=1,m
!!$
!!$             IF(ncell(i,j,k) .NE. 0) THEN
!!$
!!$                fvx(i,j,k)=fvx(i,j,k)/ncell(i,j,k)
!!$                fvy(i,j,k)=fvy(i,j,k)/ncell(i,j,k)
!!$                fvz(i,j,k)=fvz(i,j,k)/ncell(i,j,k)
!!$
!!$             END IF
!!$
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    DEALLOCATE(ncell)
!!$
!!$    !Then smooth velocity field!
!!$
!!$    CALL fft3(fvx,fvx,-1)
!!$    CALL fft3(fvy,fvy,-1)
!!$    CALL fft3(fvz,fvz,-1)
!!$
!!$    DO k=1,m
!!$       DO j=1,m
!!$          DO i=1,m
!!$
!!$             CALL k_fft(i,j,k,m,kx,ky,kz,kmod,L)
!!$
!!$             fvx(i,j,k)=fvx(i,j,k)*exp(-((kmod*r)**2.)/2.)
!!$             fvy(i,j,k)=fvy(i,j,k)*exp(-((kmod*r)**2.)/2.)
!!$             fvz(i,j,k)=fvz(i,j,k)*exp(-((kmod*r)**2.)/2.)
!!$
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    CALL fft3(fvx,fvx,1)
!!$    CALL fft3(fvy,fvy,1)
!!$    CALL fft3(fvz,fvz,1)
!!$
!!$    fvx=fvx/(float(m)**3.)
!!$    fvy=fvy/(float(m)**3.)
!!$    fvz=fvz/(float(m)**3.)
!!$
!!$  END SUBROUTINE vfield

  SUBROUTINE NGP(x,L,w,d,m)

    IMPLICIT NONE
    INTEGER :: i, n
    DOUBLE COMPLEX, ALLOCATABLE, INTENT(INOUT) :: d(:,:,:)
    REAL, INTENT(IN) :: x(:,:), L, w(:)
    INTEGER :: ix, iy, iz
    INTEGER, INTENT(IN) :: m

    WRITE(*,*) 'NGP: Binning particles and creating density field'
    WRITE(*,*) 'NGP: Cells:', m
    IF(SIZE(x(:,1)) .NE. 3) STOP 'NGP: Array must be 3 dimensional'

    ALLOCATE(d(m,m,m))

    n=SIZE(x(1,:))
    IF(SIZE(w) .NE. n) STOP 'NGP: Weight function must be defined for all particles'

    DO i=1,n

       ix=CEILING(x(1,i)*float(m)/L)
       iy=CEILING(x(2,i)*float(m)/L)
       iz=CEILING(x(3,i)*float(m)/L)

       IF(ix>m .OR. ix<1) THEN
          WRITE(*,*) 'x:', i, x(1,i)
          STOP 'NGP: Warning, point outside box'
       END IF

       IF(iy>m .OR. iy<1) THEN
          WRITE(*,*) 'y:', i, x(2,i)
          STOP 'NGP: Warning, point outside box'
       END IF

       IF(iz>m .OR. iz<1) THEN
          WRITE(*,*) 'z:', i, x(3,i)
          STOP 'NGP: Warning, point outside box'
       END IF

       d(ix,iy,iz)=d(ix,iy,iz)+w(i)

    END DO

    WRITE(*,*) 'NGP: Binning complete'
    WRITE(*,*)

  END SUBROUTINE NGP

!!$  SUBROUTINE NGP(d,x,y,z,L,nx,ny,nz,m)
!!$
!!$    IMPLICIT NONE
!!$    INTEGER :: i, icon, n
!!$    INTEGER, INTENT(IN) :: m
!!$    REAL, INTENT(IN) :: x(:), y(:), z(:), L
!!$    REAL :: kx, ky, kz, kmod, fcorr, fcx, fcy, fcz, kxh, kyh, kzh
!!$    INTEGER, INTENT(OUT) :: nx(:), ny(:), nz(:)
!!$    DOUBLE COMPLEX, INTENT(OUT) :: d(:,:,:)
!!$
!!$    nx=0
!!$    ny=0
!!$    nz=0
!!$    d=(0.d0,0.d0)
!!$
!!$    n=SIZE(x)
!!$
!!$    WRITE(*,*) 'Binning particles (NGP) and creating density field'
!!$    WRITE(*,*) 'With cells:', m
!!$
!!$    DO i=1,n
!!$
!!$       nx(i)=CEILING(x(i)*float(m)/L)
!!$       ny(i)=CEILING(y(i)*float(m)/L)
!!$       nz(i)=CEILING(z(i)*float(m)/L)
!!$
!!$       IF(nx(i)>m) THEN
!!$          WRITE(*,*) 'Warning, point outside box:', x(i)
!!$          READ(*,*) icon
!!$          nx(i)=1
!!$       ELSE IF(ny(i)>m) THEN
!!$          WRITE(*,*) 'Warning, point outside box:', y(i)
!!$          READ(*,*) icon
!!$          ny(i)=1
!!$       ELSE IF(nz(i)>m) THEN
!!$          WRITE(*,*) 'Warning, point outside box:', z(i)
!!$          READ(*,*) icon
!!$          nz(i)=1
!!$       ELSE IF(nx(i)<1) THEN
!!$          WRITE(*,*) 'Warning, point outside box:', x(i)
!!$          READ(*,*) icon
!!$          nx(i)=m
!!$       ELSE IF(ny(i)<1) THEN
!!$          WRITE(*,*) 'Warning, point outside box:', y(i)
!!$          READ(*,*) icon
!!$          ny(i)=m
!!$       ELSE IF(nz(i)<1) THEN
!!$          WRITE(*,*) 'Warning, point outside box:', z(i)
!!$          READ(*,*) icon
!!$          nz(i)=m
!!$       END IF
!!$
!!$       d(nx(i),ny(i),nz(i))=(1.d0,0.d0)+d(nx(i),ny(i),nz(i))
!!$
!!$    END DO
!!$
!!$    WRITE(*,*) 'Binning complete'
!!$    WRITE(*,*) 'Correcting for binning by sharpening field'
!!$
!!$    !Now correct for binning!
!!$
!!$    CALL fft3(d,d,-1)
!!$
!!$    DO i=1,m
!!$       DO j=1,m
!!$          DO k=1,m
!!$
!!$             call k_fft(i,j,k,m,kx,ky,kz,kmod,L)
!!$
!!$             kxh=L*kx/(2.*float(m))
!!$             kyh=L*ky/(2.*float(m))
!!$             kzh=L*kz/(2.*float(m))
!!$
!!$
!!$             fcx=sinc(kxh)
!!$             fcy=sinc(kyh)
!!$             fcz=sinc(kzh)
!!$
!!$             fcorr=fcx*fcy*fcz
!!$
!!$             d(i,j,k)=d(i,j,k)/fcorr
!!$
!!$          END DO
!!$       END DO
!!$    END DO
!!$
!!$    CALL fft3(d,d,1)
!!$
!!$    d=d/(float(m)**3.)
!!$
!!$    WRITE(*,*) 'Sharpening complete, NGP finished'
!!$    WRITE(*,*)
!!$
!!$  END SUBROUTINE NGP

  FUNCTION sinc(x)

    IMPLICIT NONE
    REAL :: sinc
    REAL, INTENT(IN) :: x

    IF(ABS(x)<0.1) THEN
       sinc=1.-(x**2.)/6.
    ELSE
       sinc=sin(x)/x
    END IF

  END FUNCTION sinc

  FUNCTION variance(arr)

    IMPLICIT NONE
    REAL :: variance
    REAL*8 :: var
    DOUBLE COMPLEX, INTENT(IN) :: arr(:,:,:)
    INTEGER :: n
    INTEGER :: i, j, k

    var=0.d0

    n=size(arr,1)

    !This calculates the variance in an array!

    DO i=1,n
       DO j=1,n
          DO k=1,n
             var=var+abs(arr(i,j,k))**2.
          END DO
       END DO
    END DO

    var=var/(float(n)**3.)
    var=var-average(arr)**2.

    variance=var

  END FUNCTION variance

  FUNCTION average(arr)

    IMPLICIT NONE
    REAL :: average
    REAL*8 :: avg
    DOUBLE COMPLEX, INTENT(IN) :: arr(:,:,:)
    INTEGER :: n
    INTEGER :: i, j, k

    avg=0.d0

    n=size(arr,1)

    DO i=1,n
       DO j=1,n
          DO k=1,n
             avg=avg+arr(i,j,k)
          END DO
       END DO
    END DO

    average=avg/(float(n)**3.)

  END FUNCTION average

  FUNCTION average2d(arr)

    IMPLICIT NONE
    REAL :: average2d
    REAL*8 :: avg
    DOUBLE COMPLEX :: arr(:,:)
    INTEGER :: n
    INTEGER :: i, j

    avg=0.d0

    n=size(arr,1)

    DO i=1,n
       DO j=1,n         
          avg=avg+arr(i,j)
       END DO
    END DO

    average2d=avg/(float(n)**2.)

  END FUNCTION average2d

!!$  SUBROUTINE read_gadget(cosm,z4,L,n)
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    CHARACTER(len=256) :: infile
!!$    REAL*8 :: massarr(6), z8, a8, L8
!!$    REAL*4, ALLOCATABLE :: pos(:,:), vel(:,:)
!!$!    INTEGER, ALLOCATABLE :: id(:)
!!$    REAL :: L, z4, dx
!!$    INTEGER :: np(6), np2(6), crap, n, ibox, i
!!$    TYPE(cosmology) :: cosm
!!$
!!$    IF(cosm%icos==1) THEN
!!$       WRITE(*,*) '1 - 780 Mpc/h box'
!!$       WRITE(*,*) '2 - 156 Mpc/h box'
!!$       READ(*,*) ibox
!!$       WRITE(*,*)
!!$       IF(ibox==1 .AND. z4==0.) infile='../simulations/scaling/mill_780/Data_000'
!!$       IF(ibox==2 .AND. z4==0.) infile='../simulations/scaling/mill_156/Data_000'
!!$!    ELSE IF(cosm%icos==4) THEN
!!$!       WRITE(*,*) '1 - 500 Mpc/h box'
!!$!       WRITE(*,*) '2 - 315 Mpc/h box'
!!$!       READ(*,*) ibox
!!$!       WRITE(*,*)
!!$!       IF(ibox==1) infile='../simulations/scaling/tcdm_500/Data_000.bak'
!!$!       IF(ibox==2) infile='../simulations/scaling/tcdm_old/Data_000.bak'
!!$    ELSE IF(cosm%icos==8) THEN
!!$       WRITE(*,*) '1 - 500 Mpc/h box'
!!$       WRITE(*,*) '2 - 100 Mpc/h box'
!!$       READ(*,*) ibox
!!$       WRITE(*,*)
!!$       IF(ibox==1 .AND. z4==0.00) infile='../simulations/scaling/tcdm_500/Data_001'
!!$       IF(ibox==1 .AND. z4==0.22) infile='../simulations/scaling/tcdm_500/Data_000'
!!$       IF(ibox==2 .AND. z4==0.00) infile='../simulations/scaling/tcdm_100/Data_001'
!!$       IF(ibox==2 .AND. z4==0.22) infile='../simulations/scaling/tcdm_100/Data_000'
!!$    ELSE IF(cosm%icos==10) THEN
!!$       IF(z4==0.) infile='../simulations/mg/GR_512/Data_015'
!!$    ELSE IF(cosm%icos==11) THEN
!!$       IF(z4==0.) infile='../simulations/mg/F4_512/Data_015'
!!$    ELSE IF(cosm%icos==12) THEN
!!$       IF(z4==0.) infile='../simulations/mg/F5_512/Data_015'
!!$    ELSE IF(cosm%icos==13) THEN
!!$       IF(z4==0.) infile='../simulations/mg/F6_512/Data_015'
!!$    ELSE IF(cosm%icos==14) THEN
!!$
!!$!       IF(z4==0.00) infile='../simulations/mg/LCDM_602/Data_004'
!!$!       IF(z4==0.25) infile='../simulations/mg/LCDM_602/Data_003'
!!$!       IF(z4==0.38) infile='../simulations/mg/LCDM_602/Data_002'       
!!$!       IF(z4==0.47) infile='../simulations/mg/LCDM_602/Data_001' 
!!$!       IF(z4==0.51) infile='../simulations/mg/LCDM_602/Data_000'
!!$
!!$       !F4
!!$       IF(z4==0.225) infile='../simulations/mg/LCDM_611/Data_000'
!!$       !F5
!!$       IF(z4==0.381) infile='../simulations/mg/LCDM_602/Data_002'
!!$       !F6
!!$       IF(z4==0.644) infile='../simulations/mg/LCDM_536/Data_000'
!!$       !GR
!!$       IF(z4==0.844) infile='../simulations/mg/LCDM_482/Data_000'
!!$
!!$    END IF
!!$
!!$    WRITE(*,*) 'INPUT FILE:', infile
!!$
!!$    OPEN(7,file=infile,form='unformatted',status='old')
!!$    READ(7) np, massarr, a8, z8, crap, crap, np2, crap, crap, L8
!!$    CLOSE(7)
!!$
!!$    L=L8/1000.
!!$
!!$    WRITE(*,*) 'Particle number:', np(2)
!!$    WRITE(*,*) 'Which is:', nint(np(2)**(1./3.)), 'cubed.'
!!$    WRITE(*,*) 'Box size:', L
!!$    WRITE(*,*) 'a sim:', a8
!!$    WRITE(*,*) 'z sim:', z8
!!$    WRITE(*,*) 'z choice:', z4
!!$    IF(cosm%icos==14) THEN
!!$       z4=z8
!!$       WRITE(*,*) 'Changing choice redshift to actual redshift'
!!$    END IF
!!$    WRITE(*,*)
!!$
!!$    n=np(2)
!!$
!!$    ALLOCATE(pos(3,n),vel(3,n),id(n))
!!$
!!$    CALL array_alloc(n)
!!$
!!$    OPEN(7,file=infile,form='unformatted')
!!$    READ(7)
!!$    READ(7) pos
!!$    READ(7) vel
!!$    READ(7) id
!!$    CLOSE(7)
!!$
!!$    !kpc -> Mpc conversion!
!!$    pos=pos/1000.
!!$    vel=vel*sqrt(a8)
!!$
!!$    dx=0.00001
!!$
!!$    DO i=1,n
!!$
!!$       x(i)=pos(1,i)
!!$       y(i)=pos(2,i)
!!$       z(i)=pos(3,i)
!!$
!!$       vx(i)=vel(1,i)
!!$       vy(i)=vel(2,i)
!!$       vz(i)=vel(3,i)
!!$
!!$       IF(x(i)<dx) x(i)=dx
!!$       IF(y(i)<dx) y(i)=dx
!!$       IF(z(i)<dx) z(i)=dx
!!$       IF(x(i)>L-dx) x(i)=L-dx
!!$       IF(y(i)>L-dx) y(i)=L-dx
!!$       IF(z(i)>L-dx) z(i)=L-dx
!!$
!!$    END DO
!!$
!!$    DEALLOCATE(pos,vel)
!!$
!!$    m=massarr(2)*1e10
!!$
!!$  END SUBROUTINE read_gadget

!!$  FUNCTION grow_rate(k,z,cosm)
!!$
!!$    USE cosdef
!!$    IMPLICIT NONE
!!$    REAL :: grow_rate, z, k
!!$    TYPE(cosmology) :: cosm
!!$
!!$    IF(cosm%irate==0) THEN
!!$       grow_rate=find(z,cosm%z_g,cosm%g,3,1)
!!$    ELSE IF(cosm%irate==1) THEN
!!$       grow_rate=grow_rate_linder(z,cosm)
!!$    ELSE IF(cosm%irate==2) THEN
!!$       IF(z .NE. 0.) THEN
!!$          WRITE(*,*) 'z=', z
!!$          STOP 'Error: attempting to calculate scale dependent growth rate at z .ne. 0'
!!$       END IF
!!$       grow_rate=find(log(k),log(cosm%k_r),cosm%rk,3,1)
!!$    END IF
!!$
!!$  END FUNCTION grow_rate

  FUNCTION grow_rate_linder(z,cosm)

    USE cosdef
    IMPLICIT NONE
    REAL :: grow_rate_linder, z
    TYPE(cosmology) :: cosm
    REAL :: gam


    !          IF(w<-1.) THEN
    !             gam=0.55+0.02*(1.+w)
    !          ELSE IF(w>-1) THEN
    !             gam=0.55+0.05*(1.+w)
    !          ELSE
    gam=0.55
    !          END IF

    grow_rate_linder=omega_m(z,cosm)**gam

  END FUNCTION grow_rate_linder

END PROGRAM modcat
