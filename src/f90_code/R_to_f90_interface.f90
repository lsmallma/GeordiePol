
!
! Next we create the fortran interface subroutine from R.
! This must be seperate from the module as R cannot handle compiling modules. We
! can however link directly from here to the module
! 

subroutine rtof90interface(tmp_file,tmp_char,recal_run_in,group,nos_iter_in   &
                          ,deposition_model_in,domain_extent_x_in             &
                          ,domain_extent_y_in,resolution_in,random_or_read_in &
                          ,ecotone_in,weighted_patch_in,map_fixed_in          &
                          ,nos_species_in,sp_fractions_in,patch_size_in       &
                          ,lake_size_in,lake_centred_in,lake_fixed_in         &
                          ,lake_specified_in,nos_lake_in,nos_water_in         &
                          ,water_locations_in                                 &
                          ,lake_user_x_in,lake_user_y_in,wind_sp_in,fall_speed_in   &
                          ,pollen_radius_in,pollen_productivity_in            &
                          ,i_extent_in,j_extent_in,mean_sp_abundance_at_distance_in &
                          ,sum_pol_out)
                          
  ! be explicit about memory variables as there are other internal variables
  use pollen_deposition_model_mod, only: iterative_loop,idum,recal_run &
                                        ,sum_pol,nos_iter,deposition_model       &
                                        ,domain_extent_x,domain_extent_y         &
                                        ,resolution,random_or_read,ecotone       &
                                        ,weighted_patch,map_fixed,nos_species    &
                                        ,sp_fractions,patch_size,lake_size       &
                                        ,lake_centred,lake_fixed,lake_specified  & 
                                        ,nos_lake,lake_user_x,lake_user_y        &
                                        ,wind_sp,fall_speed,pollen_radius        &
                                        ,pollen_productivity,i_extent,j_extent   &
                                        ,distance_array,sp_hold,resolution2      &
                                        ,mean_sp_abundance_at_distance,nos_water &
                                        ,water_locations

  implicit none

  ! declare input variables
  integer, intent(inout) :: ecotone_in,weighted_patch_in,map_fixed_in            &
                           ,lake_centred_in,lake_fixed_in,lake_specified_in      &
                           ,random_or_read_in,nos_iter_in,deposition_model_in    &
                           ,domain_extent_x_in,domain_extent_y_in,nos_species_in & 
                           ,group,tmp_char,recal_run_in,nos_lake_in,nos_water_in

  integer, dimension(nos_lake_in), intent(inout) :: lake_user_x_in,lake_user_y_in
  integer, dimension(domain_extent_x_in*domain_extent_y_in), intent(inout) :: i_extent_in &
                                                                             ,j_extent_in
  double precision, intent(inout) :: resolution_in,patch_size_in,lake_size_in,wind_sp_in
  double precision, dimension(nos_water_in), intent(inout) :: water_locations_in
  double precision, dimension(nos_species_in), intent(inout) :: sp_fractions_in,fall_speed_in, &
                                                                pollen_radius_in,pollen_productivity_in
  double precision, dimension(domain_extent_x_in, & 
                              domain_extent_y_in,nos_species_in), intent(inout) :: mean_sp_abundance_at_distance_in
  double precision, dimension(nos_species_in,nos_iter_in), intent(inout):: sum_pol_out
  character(tmp_char), intent(inout) :: tmp_file

  ! declare local variables
  integer iter, time1, time2, time3
  character(len=20) :: write_fmt
!if (i == 1) then
!    open(unit=666,file="/home/lsmallma/WORK/R/Scripts/EdHulPOL/src/r_code/out.csv", & 
!         status='replace',action='readwrite' )
!    write(666,*),"1"
!    close(666)
!endif

  ! allocate the output variable
  if (allocated(sum_pol)) deallocate(sum_pol)
  allocate(sum_pol(nos_species_in,nos_iter_in))

  ! allocate variables with the following dimensions: (domain_extent_x,domain_extent_y,nos_species)
  if (allocated(mean_sp_abundance_at_distance)) deallocate (mean_sp_abundance_at_distance)
  allocate(mean_sp_abundance_at_distance(domain_extent_x_in,domain_extent_y_in,nos_species_in))

  ! allocate variables with the following dimensions: (nos_species)
  if (allocated(fall_speed)) deallocate(fall_speed)
  if (allocated(sp_fractions)) deallocate(sp_fractions)
  if (allocated(pollen_productivity)) deallocate(pollen_productivity)
  if (allocated(pollen_radius)) deallocate(pollen_radius)
  allocate(fall_speed(nos_species_in),sp_fractions(nos_species_in) &
          ,pollen_productivity(nos_species_in),pollen_radius(nos_species_in))

  ! allocate variables with the following dimensions:
  ! (domain_extent_x*domain_extent_y)
  if (allocated(i_extent)) deallocate(i_extent)
  if (allocated(j_extent)) deallocate(j_extent)
  if (allocated(distance_array)) deallocate(distance_array)
  allocate(i_extent(domain_extent_x_in*domain_extent_y_in), &
           j_extent(domain_extent_x_in*domain_extent_y_in), &
           distance_array(domain_extent_x_in*domain_extent_y_in))

  nos_lake = nos_lake_in ; nos_water = nos_water_in

  ! allocate water pixel information
  if (allocated(water_locations)) deallocate(water_locations)
  allocate(water_locations(nos_water)) 

  ! allocate variables with the following dimensions: (nos_lake)
  if (allocated(lake_user_x)) deallocate(lake_user_x)
  if (allocated(lake_user_y)) deallocate(lake_user_y)
  allocate(lake_user_x(nos_lake),lake_user_y(nos_lake))

  ! some things to pass along too
  sum_pol = 0.0
  nos_iter = nos_iter_in               ; deposition_model = deposition_model_in
  domain_extent_x = domain_extent_x_in ; domain_extent_y = domain_extent_y_in
  resolution = resolution_in           ; resolution2 = resolution*resolution
  nos_species = nos_species_in         ; sp_fractions = sp_fractions_in       
  patch_size = patch_size_in           ; lake_size = lake_size_in       
  lake_user_x = lake_user_x_in         ; water_locations=nint(water_locations_in)
  lake_user_y = lake_user_y_in         ; wind_sp = wind_sp_in                 
  fall_speed = fall_speed_in           ; pollen_radius = pollen_radius_in
  i_extent = i_extent_in               ; j_extent = j_extent_in 
  pollen_productivity = pollen_productivity_in
  mean_sp_abundance_at_distance = mean_sp_abundance_at_distance_in

  ! variables used in map generation
  if (allocated(sp_hold)) deallocate(sp_hold)
  if (random_or_read_in == 1) then
     allocate(sp_hold(domain_extent_x*domain_extent_y,nos_species))
  endif

  if (recal_run_in == 1) then
      recal_run = .true.
  else
      recal_run = .false.
  endif
  if (random_or_read_in == 1) then
      random_or_read = .true.
  else
      random_or_read = .false.
  endif
  if (ecotone_in == 1) then
      ecotone = .true.
  else
      ecotone = .false.
  endif
  if (weighted_patch_in == 1) then
      weighted_patch = .true.
  else
      weighted_patch = .false.
  endif
  if (map_fixed_in == 1) then
      map_fixed = .true.
  else
      map_fixed = .false.
  endif
  if (lake_centred_in == 1) then
      lake_centred = .true.
  else
      lake_centred = .false.
  endif
  if (lake_fixed_in == 1) then
      lake_fixed = .true.
  else
      lake_fixed = .false.
  end if
  if (lake_specified_in == 1) then
      lake_specified = .true.
  else
      lake_specified = .false.
  endif

  ! seed the random number generator
  ! determine unique (sort of) seed value; based on system time
  call system_clock(time1,time2,time3)
  ! set seed value outside of the function, idum must be a negative number
  idum=time1+time2+time3

  ! if appropriate open up the restart file for writing 
  if (.not.recal_run) then
      open(unit = 777, file=trim(tmp_file), status='unknown' ,action='write' ,position='append')
      write(write_fmt,fmt='(a,i0,a)')'(f0.7,',nos_species-1,'(",",f0.7))'
  endif

  ! loop through all iterations requested
  do iter = 1, nos_iter
     call iterative_loop(iter) 
     ! again if appropriate write to the restart file
     if (.not.recal_run) write(unit=777,fmt=trim(write_fmt))sum_pol(1:nos_species,iter)
  end do

  ! align outputs
  sum_pol_out = sum_pol

  ! close the tmp file
  if (.not.recal_run) close(777)

  ! return back to the user
  return

end subroutine rtof90interface

!
!----------------------------------------------------------------------------
!

