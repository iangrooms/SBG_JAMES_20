program smooth
  use netcdf
  implicit none

  ! We are reading 3D data, a 2400 x 3600 x 62 grid. 
  integer, parameter :: NX = 3600, NY = 2400, NZ = 62
  integer, parameter :: rec_length = 8*NX*NY
  double precision, allocatable, dimension(:,:,:) :: smoothed, field, fieldT, PD
  double precision, dimension(-179:NX+180,NY,NZ) :: mask
  double precision, dimension(NY,NX,NZ) :: maskT
  double precision, dimension(NX,NY) :: TAREA, TLON, TLAT
  double precision, dimension(NX,NY) :: UAREA, ULON, ULAT
  double precision :: weightsX(361), weightsY(241,121:1182)
  integer :: i, j, k ! Loop indices
  integer :: KMT(NX,NY)
  integer :: stat_ncdf ! netCDF status
  integer, allocatable, dimension(:,:) :: KMU

  ! This will be the netCDF ID for the file and data variable.
  integer :: ncid, varid

  ! Some timing variables
  integer :: tCount0, tCount1, tCountRate


! Read smaller auxiliary fields
  ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to
  ! the file.
  stat_ncdf = nf90_open("current_data", NF90_NOWRITE, ncid) )
  ! Get the varid of the data variable, based on its name.
  ! Read TAREA
  stat_ncdf = nf90_inq_varid(ncid, "TAREA", varid)
  stat_ncdf = nf90_get_var(ncid, varid, TAREA)

  ! Read Longitude (degrees)
  stat_ncdf = nf90_inq_varid(ncid, "TLONG", varid)
  stat_ncdf = nf90_get_var(ncid, varid, TLON)

  ! Read Latitude (degrees)
  stat_ncdf = nf90_inq_varid(ncid, "TLAT", varid)
  stat_ncdf = nf90_get_var(ncid, varid, TLAT)

  ! Read UAREA
  stat_ncdf = nf90_inq_varid(ncid, "UAREA", varid)
  stat_ncdf = nf90_get_var(ncid, varid, UAREA)

  ! Read Longitude (degrees)
  stat_ncdf = nf90_inq_varid(ncid, "ULONG", varid)
  stat_ncdf = nf90_get_var(ncid, varid, ULON)

  ! Read Latitude (degrees)
  stat_ncdf = nf90_inq_varid(ncid, "ULAT", varid)
  stat_ncdf = nf90_get_var(ncid, varid, ULAT)

  ! Read KMT
  stat_ncdf = nf90_inq_varid(ncid, "KMT", varid)
  stat_ncdf = nf90_get_var(ncid, varid, KMT)

  ! Read KMU
  allocate(KMU(NX,NY))
  stat_ncdf = nf90_inq_varid(ncid, "KMU", varid)
  stat_ncdf = nf90_get_var(ncid, varid, KMU)

! Set land mask
  KMT = MIN(KMT,KMU)
  deallocate(KMU)
  mask = 1.
  do k=1,NZ
    do j=1,NY
      do i=1,NX
        if ( KMT(i,j) <= k-1 ) then
          mask(i,j,k) = 0.
        end if
      end do
    end do
    call per_ext(mask(:,:,k))
    maskT(:,:,k) = transpose(mask(1:NX,:,k))
  end do

  call system_clock(tCount0,tCountRate) 
! Read PD
allocate(field(-179:NX+180,NY,NZ))
  stat_ncdf = nf90_inq_varid(ncid, "PD", varid)
  stat_ncdf = nf90_get_var(ncid, varid, field(1:NX,:,:))
  ! Close the file, freeing all resources.
  stat_ncdf = nf90_close(ncid)
  call system_clock(tCount1)
  print *, 'Reading PD: ', real(tCount1-tCount0)/real(tCountRate)

! Set values on land to one immediately above
  do k=2,NZ
    do j=1,NY
      do i=1,NX
        if ( KMT(i,j) <= k-1 ) then
          field(i,j,k) = field(i,j,k-1) 
        end if
      end do
    end do
    call per_ext(field(:,:,k))
  end do
  call per_ext(field(:,:,1))

! Smooth PD
allocate(smoothed(NX,NY,NZ))
  smoothed = field(1:NX,:,:)
! Average over longitude
  call system_clock(tCount0)
!$OMP PARALLEL DO PRIVATE(weightsX,i,j,k)
do j=1,1302
  ! Compute the weights
  ! i=2820--3020 have nonzero lat/lon coordinates at all j
  weightsX(181:360) = TAREA(3000,j)*exp( -(.5*.75**2) * &
        &  (((TLON(2821:3000,j) - TLON(2821,j))/cos(TLAT(3000,j)))**2) )
  weightsX(1:180) = TAREA(3000,j)*exp( -(.5*.75**2) * &
        &  (((TLON(2821:3000,j) - TLON(3000,j))/cos(TLAT(3000,j)))**2) ) 
  ! Loop over depth
  do k=1,NZ 
    ! Loop over longitude
    do i=1,NX
      ! Don't go below the bottom 
      if( KMT(i,j) - k > 1 ) then
        smoothed(i,j,k) = sum(mask((i-180):(i+180),j,k)*weightsX*field((i-180):(i+180),j,k))&
   &/sum(mask((i-180):(i+180),j,k)*weightsX)
        end if
    end do
  end do
end do
!$OMP END PARALLEL DO
!deallocate(field)
  call system_clock(tCount1)
  print *, 'First smoothing loop: ', real(tCount1-tCount0)/real(tCountRate) 

! Average over latitude
allocate(fieldT(NY,NX,NZ))
  do k=1,NZ
    fieldT(:,:,k) = transpose(smoothed(:,:,k))
  end do
deallocate(smoothed)
allocate(smoothed(NY,NX,NZ))

! Precompute weightsY
do j=121,1182
  weightsY(:,j) = TAREA(3000,(j-120):(j+120))*exp( -(.5*.75**2) * &
        &  ((TLAT(3000,(j-120):(j+120)) - TLAT(3000,j))**2) ) 
  weightsY(:,j) = weightsY(:,j)/maxval(weightsY(:,j))
end do

  call system_clock(tCount0)
!$OMP PARALLEL DO PRIVATE(i,j,k)
do i=1,NX
  ! Loop over depth
  do k=1,NZ
    ! Loop over latitude
    do j=121,1182
      ! Skip regions that are too shallow
      if( KMT(i,j) - k > 1 ) then
        smoothed(j,i,k) = sum(maskT((j-120):(j+120),i,k)*weightsY(:,j)*fieldT((j-120):(j+120),i,k))&
   &/sum(maskT((j-120):(j+120),i,k)*weightsY(:,j))
        end if
    end do
  end do
end do
!$OMP END PARALLEL DO
deallocate(fieldT)
  call system_clock(tCount1)
  print *, 'Second smoothing loop: ', real(tCount1-tCount0)/real(tCountRate)

! Write smoothed PD
  open(unit=21,file='PDUSmooth3D.dat',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=rec_length)
  do k=1,NZ
    write(21,REC=k) transpose(smoothed(:,:,k))
  end do
  close(21)
deallocate(smoothed)


  call system_clock(tCount0,tCountRate) 
! Read UVEL
allocate(PD(-179:NX+180,NY,NZ))
  PD = field
  stat_ncdf = nf90_open("current_data", NF90_NOWRITE, ncid)
  stat_ncdf = nf90_inq_varid(ncid, "UVEL", varid)
  stat_ncdf = nf90_get_var(ncid, varid, field(1:NX,:,:))
  ! Close the file, freeing all resources.
  stat_ncdf = nf90_close(ncid)
  call system_clock(tCount1)
  print *, 'Reading UVEL: ', real(tCount1-tCount0)/real(tCountRate)

! Set values on land to one immediately above
  do k=2,NZ
    do j=1,NY
      do i=1,NX
        if ( KMT(i,j) <= k-1 ) then
          field(i,j,k) = field(i,j,k-1)
        end if
      end do
    end do
    call per_ext(field(:,:,k))
  end do
  call per_ext(field(:,:,1))

! Smooth UVEL
allocate(smoothed(NX,NY,NZ))
  smoothed = field(1:NX,:,:)
! Average over longitude
  call system_clock(tCount0)
!$OMP PARALLEL DO PRIVATE(weightsX,i,j,k)
do j=1,1302
  ! Compute the weights
  ! i=2820--3020 have nonzero lat/lon coordinates at all j
  weightsX(181:360) = UAREA(3000,j)*exp( -(.5*.75**2) * &
        &  (((ULON(2821:3000,j) - ULON(2821,j))/cos(ULAT(3000,j)))**2) )
  weightsX(1:180) = UAREA(3000,j)*exp( -(.5*.75**2) * &
        &  (((ULON(2821:3000,j) - ULON(3000,j))/cos(ULAT(3000,j)))**2) ) 
  ! Loop over depth
  do k=1,NZ 
    ! Loop over longitude
    do i=1,NX
      ! Skip regions that are too shallow
      if( KMT(i,j) - k > 1 ) then
        smoothed(i,j,k) = sum(mask((i-180):(i+180),j,k)*weightsX*field((i-180):(i+180),j,k))&
   &/sum(mask((i-180):(i+180),j,k)*weightsX)
        end if
    end do
  end do
end do
!$OMP END PARALLEL DO
!deallocate(field)
  call system_clock(tCount1)
  print *, 'First smoothing loop: ', real(tCount1-tCount0)/real(tCountRate) 

! Average over latitude
allocate(fieldT(NY,NX,NZ))
  do k=1,NZ
    fieldT(:,:,k) = transpose(smoothed(:,:,k))
  end do
deallocate(smoothed)
allocate(smoothed(NY,NX,NZ))

! Precompute weightsY
do j=121,1182
  weightsY(:,j) = UAREA(3000,(j-120):(j+120))*exp( -(.5*.75**2) * &
        &  ((ULAT(3000,(j-120):(j+120)) - ULAT(3000,j))**2) ) 
  weightsY(:,j) = weightsY(:,j)/maxval(weightsY(:,j))
end do

  call system_clock(tCount0)
!$OMP PARALLEL DO PRIVATE(i,j,k)
do i=1,NX
  ! Loop over depth
  do k=1,NZ
    ! Loop over latitude
    do j=121,1182
      ! Skip regions that are too shallow
      if( KMT(i,j) - k > 1 ) then
        smoothed(j,i,k) = sum(maskT((j-120):(j+120),i,k)*weightsY(:,j)*fieldT((j-120):(j+120),i,k))&
   &/sum(maskT((j-120):(j+120),i,k)*weightsY(:,j))
        end if
    end do
  end do
end do
!$OMP END PARALLEL DO
deallocate(fieldT)
  call system_clock(tCount1)
  print *, 'Second smoothing loop: ', real(tCount1-tCount0)/real(tCountRate)

! Write smoothed UVEL
  open(unit=21,file='UVELSmooth3D.dat',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=rec_length)
  do k=1,NZ
    write(21,REC=k) transpose(smoothed(:,:,k))
  end do
  close(21)
deallocate(smoothed)


! Get U*PD
  call system_clock(tCount0,tCountRate)
  field = field*PD
! Smooth U*PD
allocate(smoothed(NX,NY,NZ))
  smoothed = field(1:NX,:,:)
! Average over longitude
  call system_clock(tCount0)
!$OMP PARALLEL DO PRIVATE(weightsX,i,j,k)
do j=1,1302
  ! Compute the weights
  ! i=2820--3020 have nonzero lat/lon coordinates at all j
  weightsX(181:360) = UAREA(3000,j)*exp( -(.5*.75**2) * &
        &  (((ULON(2821:3000,j) - ULON(2821,j))/cos(ULAT(3000,j)))**2) )
  weightsX(1:180) = UAREA(3000,j)*exp( -(.5*.75**2) * &
        &  (((ULON(2821:3000,j) - ULON(3000,j))/cos(ULAT(3000,j)))**2) ) 
  ! Loop over depth
  do k=1,NZ 
    ! Loop over longitude
    do i=1,NX
      ! Skip regions that are too shallow
      if( KMT(i,j) - k > 1 ) then
        smoothed(i,j,k) = sum(mask((i-180):(i+180),j,k)*weightsX*field((i-180):(i+180),j,k))&
   &/sum(mask((i-180):(i+180),j,k)*weightsX)
        end if
    end do
  end do
end do
!$OMP END PARALLEL DO
deallocate(field)
  call system_clock(tCount1)
  print *, 'First smoothing loop: ', real(tCount1-tCount0)/real(tCountRate) 

! Average over latitude
allocate(fieldT(NY,NX,NZ))
  do k=1,NZ
    fieldT(:,:,k) = transpose(smoothed(:,:,k))
  end do
deallocate(smoothed)
allocate(smoothed(NY,NX,NZ))

  call system_clock(tCount0)
!$OMP PARALLEL DO PRIVATE(i,j,k)
do i=1,NX
  ! Loop over depth
  do k=1,NZ
    ! Loop over latitude
    do j=121,1182
      ! Skip regions that are too shallow
      if( KMT(i,j) - k > 1 ) then
        smoothed(j,i,k) = sum(maskT((j-120):(j+120),i,k)*weightsY(:,j)*fieldT((j-120):(j+120),i,k))&
   &/sum(maskT((j-120):(j+120),i,k)*weightsY(:,j))
        end if
    end do
  end do
end do
!$OMP END PARALLEL DO
deallocate(fieldT)
  call system_clock(tCount1)
  print *, 'Second smoothing loop: ', real(tCount1-tCount0)/real(tCountRate)

! Write smoothed U*PD
  open(unit=21,file='UPDSmooth3D.dat',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=rec_length)
  do k=1,NZ
    write(21,REC=k) transpose(smoothed(:,:,k))
  end do
  close(21)
deallocate(smoothed)


  call system_clock(tCount0,tCountRate) 
! Read VVEL
allocate(field(-179:NX+180,NY,NZ))
  stat_ncdf = nf90_open("current_data", NF90_NOWRITE, ncid)
  stat_ncdf = nf90_inq_varid(ncid, "VVEL", varid)
  stat_ncdf = nf90_get_var(ncid, varid, field(1:NX,:,:))
  ! Close the file, freeing all resources.
  stat_ncdf = nf90_close(ncid)
  call system_clock(tCount1)
  print *, 'Reading VVEL: ', real(tCount1-tCount0)/real(tCountRate)

! Set values on land to one immediately above
  do k=2,NZ
    do j=1,NY
      do i=1,NX
        if ( KMT(i,j) <= k-1 ) then
          field(i,j,k) = field(i,j,k-1)
        end if
      end do
    end do
    call per_ext(field(:,:,k))
  end do
  call per_ext(field(:,:,1))

! Smooth VVEL
allocate(smoothed(NX,NY,NZ))
  smoothed = field(1:NX,:,:)
! Average over longitude
  call system_clock(tCount0)
!$OMP PARALLEL DO PRIVATE(weightsX,i,j,k)
do j=1,1302
  ! Compute the weights
  ! i=2820--3020 have nonzero lat/lon coordinates at all j
  weightsX(181:360) = UAREA(3000,j)*exp( -(.5*.75**2) * &
        &  (((ULON(2821:3000,j) - ULON(2821,j))/cos(ULAT(3000,j)))**2) )
  weightsX(1:180) = UAREA(3000,j)*exp( -(.5*.75**2) * &
        &  (((ULON(2821:3000,j) - ULON(3000,j))/cos(ULAT(3000,j)))**2) ) 
  ! Loop over depth
  do k=1,NZ 
    ! Loop over longitude
    do i=1,NX
      ! Skip regions that are too shallow
      if( KMT(i,j) - k > 1 ) then
        smoothed(i,j,k) = sum(mask((i-180):(i+180),j,k)*weightsX*field((i-180):(i+180),j,k))&
   &/sum(mask((i-180):(i+180),j,k)*weightsX)
        end if
    end do
  end do
end do
!$OMP END PARALLEL DO
!deallocate(field)
  call system_clock(tCount1)
  print *, 'First smoothing loop: ', real(tCount1-tCount0)/real(tCountRate) 

! Average over latitude
allocate(fieldT(NY,NX,NZ))
  do k=1,NZ
    fieldT(:,:,k) = transpose(smoothed(:,:,k))
  end do
deallocate(smoothed)
allocate(smoothed(NY,NX,NZ))

  call system_clock(tCount0)
!$OMP PARALLEL DO PRIVATE(i,j,k)
do i=1,NX
  ! Loop over depth
  do k=1,NZ
    ! Loop over latitude
    do j=121,1182
      ! Skip regions that are too shallow
      if( KMT(i,j) - k > 1 ) then
        smoothed(j,i,k) = sum(maskT((j-120):(j+120),i,k)*weightsY(:,j)*fieldT((j-120):(j+120),i,k))&
   &/sum(maskT((j-120):(j+120),i,k)*weightsY(:,j))
        end if
    end do
  end do
end do
!$OMP END PARALLEL DO
deallocate(fieldT)
  call system_clock(tCount1)
  print *, 'Second smoothing loop: ', real(tCount1-tCount0)/real(tCountRate)

! Write smoothed VVEL
  open(unit=21,file='VVELSmooth3D.dat',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=rec_length)
  do k=1,NZ
    write(21,REC=k) transpose(smoothed(:,:,k))
  end do
  close(21)
deallocate(smoothed)


! Get V*PD
  call system_clock(tCount0,tCountRate)
  field = field*PD
deallocate(PD)
! Smooth V*PD
allocate(smoothed(NX,NY,NZ))
  smoothed = field(1:NX,:,:)
! Average over longitude
  call system_clock(tCount0)
!$OMP PARALLEL DO PRIVATE(weightsX,i,j,k)
do j=1,1302
  ! Compute the weights
  ! i=2820--3020 have nonzero lat/lon coordinates at all j
  weightsX(181:360) = UAREA(3000,j)*exp( -(.5*.75**2) * &
        &  (((ULON(2821:3000,j) - ULON(2821,j))/cos(ULAT(3000,j)))**2) )
  weightsX(1:180) = UAREA(3000,j)*exp( -(.5*.75**2) * &
        &  (((ULON(2821:3000,j) - ULON(3000,j))/cos(ULAT(3000,j)))**2) ) 
  ! Loop over depth
  do k=1,NZ 
    ! Loop over longitude
    do i=1,NX
      ! Skip regions that are too shallow
      if( KMT(i,j) - k > 1 ) then
        smoothed(i,j,k) = sum(mask((i-180):(i+180),j,k)*weightsX*field((i-180):(i+180),j,k))&
   &/sum(mask((i-180):(i+180),j,k)*weightsX)
        end if
    end do
  end do
end do
!$OMP END PARALLEL DO
deallocate(field)
  call system_clock(tCount1)
  print *, 'First smoothing loop: ', real(tCount1-tCount0)/real(tCountRate) 

! Average over latitude
allocate(fieldT(NY,NX,NZ))
  do k=1,NZ
    fieldT(:,:,k) = transpose(smoothed(:,:,k))
  end do
deallocate(smoothed)
allocate(smoothed(NY,NX,NZ))

  call system_clock(tCount0)
!$OMP PARALLEL DO PRIVATE(i,j,k)
do i=1,NX
  ! Loop over depth
  do k=1,NZ
    ! Loop over latitude
    do j=121,1182
      ! Skip regions that are too shallow
      if( KMT(i,j) - k > 1 ) then
        smoothed(j,i,k) = sum(maskT((j-120):(j+120),i,k)*weightsY(:,j)*fieldT((j-120):(j+120),i,k))&
   &/sum(maskT((j-120):(j+120),i,k)*weightsY(:,j))
        end if
    end do
  end do
end do
!$OMP END PARALLEL DO
deallocate(fieldT)
  call system_clock(tCount1)
  print *, 'Second smoothing loop: ', real(tCount1-tCount0)/real(tCountRate)

! Write smoothed V*PD
  open(unit=21,file='VPDSmooth3D.dat',access='DIRECT',form='UNFORMATTED',status='UNKNOWN',RECL=rec_length)
  do k=1,NZ
    write(21,REC=k) transpose(smoothed(:,:,k))
  end do
  close(21)
deallocate(smoothed)


contains
  subroutine per_ext(dummy)
    double precision, dimension(-179:NX+180,NY), intent(inout) :: dummy
    dummy(-179:0,:) = dummy(NX-179:NX,:)
    dummy(NX+1:NX+180,:) = dummy(1:180,:)
  end subroutine per_ext
end program smooth
