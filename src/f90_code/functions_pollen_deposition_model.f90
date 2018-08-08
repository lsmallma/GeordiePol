
!###
!## Script contains all fortran source code for the deposition model functions
!## If you are not familiar with fortran please leave this alone...
!###

module pollen_deposition_model_mod

implicit none

! assume private
private

! explicit publics
public :: iterative_loop,idum,recal_run &
         ,sum_pol,nos_iter,deposition_model      &
         ,domain_extent_x,domain_extent_y        &
         ,resolution,random_or_read,ecotone      &
         ,weighted_patch,map_fixed,nos_species   &
         ,sp_fractions,patch_size,lake_size      &
         ,lake_centred,lake_fixed,lake_specified &
         ,nos_lake,lake_user_x,lake_user_y       &
         ,wind_sp,fall_speed,pollen_radius       &
         ,pollen_productivity,i_extent,j_extent  &
         ,mean_sp_abundance_at_distance          &
         ,sp_hold,distance_array,resolution2     &
         ,nos_water,water_locations 

! declare needed constants
double precision, parameter :: pi = 3.1415926535           & ! pi 
                            ,pi_2 = 3.1415926535*2         & ! pi*2
                            ,sqrt_pi = sqrt(pi)            &
                            ,dble_one = 1.0, dble_zero = 0 &
                            ,gravity = 9.80665             & ! acceleration due to gravity (m.s-2)
                            ,mu = 1.8e-2                   & ! dynamic viscosity (g.m-1.s-1)
                            ,nine_mu_1 = (9.0*mu)**(-1.0)  &
                            ,rho0 = 2e6                    & ! particle density (g.m-3)
                            ,diffusion_coefficient = 0.12  & ! (0.30) vertical diffusion
                                                             ! coefficient (m^1/8)
                            ,rho = 1.27e6                  & ! air density (g.m-3)
                            ,gamma1 = 0.125                & ! turbulence parameter (2*gamma1)
                            ,two_gamma1 = 2.0*gamma1


! seed value for random number generater
double precision :: idum
! counter optionally needed for the 'which' subroutine
integer :: nos_out
double precision, dimension(:), allocatable :: which_out
! variables needed for the 'subsample' subroutine
double precision, dimension(:), allocatable :: subsample_out &
                                              ,subsample_locations_out

! global variables
double precision, dimension(:,:,:,:), allocatable :: pollen_deposition
double precision, dimension(:,:,:), allocatable :: mean_sp_abundance_at_distance &
                                                  ,distance_from_lake
double precision, dimension(:,:), allocatable :: pollen_load_at_lake_prediction &
                                                ,sp_hold
double precision, dimension(:), allocatable :: distance_array
double precision :: total_area
integer, dimension(:), allocatable :: water_locations
integer :: lake_centre_i,lake_centre_j ! lake pixel coordinates in current use

! model inputs
double precision, dimension(:),allocatable :: fall_speed          & ! pollen fall speed (ms-1)
                                             ,pollen_radius       & ! average pollen radius (m)
                                             ,sp_fractions        & ! cover fraction for each species
                                             ,pollen_productivity   ! relative pollen productivity
double precision :: resolution    & ! spatial resolution (m) 
                   ,resolution2   & ! resolution**2
                   ,patch_size    & ! radius of cleared patches (m)
                   ,each_patch_area_1 & ! (pi * patch_size * patch_size) ** -1
                   ,lake_size     & ! radius of lake (m)
                   ,wind_sp         ! prevailing wind speed (m.s-1)

integer, dimension(:), allocatable :: lake_user_x, lake_user_y  ! lake pixel coordinate
integer, dimension(:), allocatable :: i_extent,j_extent          ! domain description
integer :: deposition_model & ! deposition model to use (1 = 1/d, 2 = 1/d2, 3 = S-P model)
          ,nos_water        & ! number of water pixels
          ,nos_lake         &
          ,nos_iter         & ! number of iterations
          ,nos_species      & ! number of species
          ,domain_extent_x  & ! user defined nos pixels in x 
          ,domain_extent_y    ! user defined nos pixels in y

logical :: ecotone        & ! true = do ecotone experiment
          ,lake_centred   & ! true = lake is in centre of domain
          ,lake_fixed     & ! true = lake stays in same location each iteration
          ,lake_specified & ! true = lake position user provided
          ,weighted_patch & ! true = weight patches to some location
          ,map_fixed        ! true = use same map for all iterations    

logical :: random_or_read ! "read" (false) from file or "random" (true) generation
logical :: recal_run = .false.

! model outputs
double precision, dimension(:,:), allocatable :: sum_pol ! nos_species,nos_iter

contains
!
!--------------------------------------------------------------------------------------------------
!
  !
  !----------------------------------------------------------------------------
  !
  subroutine which(input,length,greaterlessequal,condition_num)

    implicit none

    ! declare inputs
    integer,intent(in) :: length, greaterlessequal
    double precision, intent(in) :: condition_num
    double precision, dimension(length), intent(in) :: input

    ! local variables
    integer :: i,j
    integer, dimension(length) :: nos_out_vector, var_index
    double precision :: condition_num_trunc
    double precision, dimension(length) :: input_trunc
 
    if (allocated(which_out)) deallocate(which_out)
    allocate(which_out(length))

    ! reset output first
    which_out = -9999 ; nos_out = 0 ; nos_out_vector = 0 ; j = 1
    do i = 1, length
       var_index(i) = i
    enddo

    ! restrict the data to three decimal places for the comparison
    input_trunc=dnint(input*1e3) ; input_trunc=input_trunc*1e-3
    condition_num_trunc=dnint(condition_num*1e3) ; condition_num_trunc=condition_num_trunc*1e-3

    ! decide condition type wanted
    if (greaterlessequal == 1) then
       ! now start looping through the datasets and assess the conditions
       where (input_trunc >= condition_num_trunc)
          which_out = var_index ; nos_out_vector = 1
       end where
    else if (greaterlessequal == 2) then
       ! now start looping through the datasets and assess the conditions
       where (input_trunc < condition_num_trunc)
          which_out = var_index ; nos_out_vector = 1
       end where
    else if (greaterlessequal == 3) then
       ! now start looping through the datasets and assess the conditions
       where (input_trunc == condition_num_trunc)
          which_out = var_index ; nos_out_vector = 1
       end where
    else
       print*,"you have not specified a valid greater less or equal"
    endif ! greaterlossequal

    ! sum total number selected
    nos_out = sum(nos_out_vector)

    ! return the solution
    return 

  end subroutine which
  !
  !----------------------------------------------------------------------------
  !
  double precision function randn(option)

    ! from Numerical Receipes p271 Press et al., 1986 2nd Edition Chapter 7,
    ! Random Numbers function returns real random number between 0-1 based on an initial start
    ! point. The start point (default = -1) is reinitialised every time the model runs
    ! providing the same distribution each run
    ! To ensure random numbers each time use the sum of the current system time

    ! modified based on blooms C code to alter range of random numbers

    implicit none
    integer IA,IM,IQ,IR,NTAB,NDIV,option
    double precision AM,EPS,RNMX,const,r1,r2
    parameter(IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32 &
             ,NDIV=1+(IM-1)/NTAB,EPS=1.2e-30,RNMX=1.-EPS)
    INTEGER j,k,iv(NTAB),iy
    SAVE iv,iy
    DATA iv /NTAB*0/
    DATA iy /0/
 
    const = 1.0

    if (option == 0) then
      ! option - is uniform random between 0 and 1
      if (idum < 0. .or. iy == 0) then
          idum=max(-idum,const)
          do j=(NTAB+8),1,-1
             k=idum/IQ
             idum=IA*(idum-k*IQ)-IR*k
             if (idum < 0) idum=idum+IM
             if (j < NTAB) iv(j)=idum
          enddo
          iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum < 0.) idum=idum+IM
      j=1.+iy/NDIV
      iy=iv(j)
      iv(j)=idum

      ! output now
      randn=min(AM*iy, RNMX)
      return

    else
      ! option is uniform distribution between -1 and 1
      ! final line adjusts this to a normal distribution roughly -2 to 2
      if (idum < 0. .or. iy == 0) then
          idum=max(-idum,const)
          do j=(NTAB+8),1,-1
             k=idum/IQ
             idum=IA*(idum-k*IQ)-IR*k
             if (idum < 0.) idum=idum+IM
             if (j < NTAB) iv(j)=idum
          enddo
          iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum < 0.) idum=idum+IM
      j=1.+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      r1=max(min(AM*iy, RNMX),1e-30)

      if (idum < 0. .or. iy == 0) then
          idum=max(-idum,const)
          do j=(NTAB+8),1,-1
             k=idum/IQ
             idum=IA*(idum-k*IQ)-IR*k
             if (idum < 0.) idum=idum+IM
             if (j < NTAB) iv(j)=idum
          enddo
          iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum < 0.) idum=idum+IM
      j=1.+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      r2=max(min(AM*iy, RNMX),1e-30)

      ! output now
      ! this line adjusts the uniform distribution to normal
      randn=sqrt(-2.*log(r1)) * cos(pi_2*r2)
      return

    endif

  end function randn
  !
  !---------------------------------------------------------------------------- 
  !
  subroutine subsample(input,length,required,replace)

    implicit none

    ! declare input variables
    integer, intent(in) :: length, required
    double precision, dimension(length),intent(in) :: input
    logical, intent(in) :: replace

    ! declare local variables
    integer :: i 
    integer :: random_num
    logical, dimension(length) :: already_selected

    ! TLS: DO I REALLY NEED TO ALLOCATE / DEALLOCATE ALL THE TIME, THIS TAKES A
    ! LOT OF TIME...
    if (allocated(subsample_out)) then
       ! we have allocated these variables, are they the correct size or do I
       ! need to deallocate / allocate them again?
       if (size(subsample_out) /= required) then
          deallocate(subsample_out,subsample_locations_out)
          allocate(subsample_out(required),subsample_locations_out(required))
       endif
    else 
       ! have not allocated these so best do that
       allocate(subsample_out(required),subsample_locations_out(required))
    endif

    ! initial values
    i = 1 ; already_selected = .true.
    ! split between replace or not at this point to avoid asking the question
    ! each time
    if (replace) then
       do while (i < required)
          random_num=nint(randn(0)*length)
          ! and allocate its position to the output variable
          subsample_out(i)=input(random_num)
          subsample_locations_out(i)=random_num
          ! rememebering to move the counter along
          i = i + 1
       end do
    else ! replace or not
       do while (i < required)
          random_num=nint(randn(0)*length)
          ! if this is something we have not already sampled then we will take it
          if ( already_selected(random_num) ) then
              ! and allocate its position to the output variable
              subsample_out(i)=input(random_num)
              subsample_locations_out(i)=random_num
              ! rememebering to move the counter along
              i = i + 1
              ! assuming that we do not want to replace the items back in the list
              ! we better make them set to false in the vector we are using to keep
              ! track
              already_selected(random_num) = .false.
          endif
       end do
    endif ! replace or not

    ! return the output now
    return

  end subroutine subsample
  !
  !----------------------------------------------------------------------------
  !
  subroutine generate_species_map(iter)

    implicit none

    ! function calculates species specific fractions for each pixel in the
    ! simulation. The pixels values can be randomly selected, random with some
    ! weighting (yet to be defined) or along an ecotone boundary

    ! declare arguments
    integer, intent(in) :: iter

    ! declare local
    double precision, dimension(domain_extent_x*domain_extent_y) :: tmp,i_extent_dble,j_extent_dble
    double precision, dimension(:), allocatable :: tmpi,tmpj &
                                                  ,patch_area,patch_centre_i,patch_centre_j
    integer, dimension(:), allocatable :: est_nos_patch, patch_locations
    double precision :: patch_done, check_interval, check_interval_max, start,finish
    integer :: sp, i, j, counter, time1, time2, time3

    !# estimate how many forest patches will be needed  (m)
    if (iter == 1) total_area = domain_extent_x*domain_extent_y*resolution2

    ! allocate some stuff
    if (allocated(est_nos_patch) .or. allocated(patch_area)) then
        if (size(patch_area) /= nos_species) deallocate(est_nos_patch,patch_area)
        allocate(est_nos_patch(nos_species),patch_area(nos_species))
    else
        allocate(est_nos_patch(nos_species),patch_area(nos_species))
    endif

    !# how much area in total will become each species
    if (ecotone) then
        patch_area(1)=total_area*(0.5*sp_fractions(1)*1e-2)
        patch_area(2)=total_area*(0.5+(0.5*sp_fractions(2)*1e-2))
        !# multiple by 2 as random allocation could mean that we overlap sometimes
        est_nos_patch=ceiling(patch_area*each_patch_area_1)
    else 
        patch_area=total_area*sp_fractions*1e-2
        !# multiple by 2 as random allocation could mean that we overlap
        !sometimes
        est_nos_patch=ceiling(patch_area*each_patch_area_1)
    end if ! ecotone or not

    !# block out these values in our maps
    sp_hold = 1
    !# assume that sp 1 is blanket coverage first
    sp_hold(:,2:nos_species) = 0
    ! convert extent integers into dble
    i_extent_dble=dble(i_extent) ; j_extent_dble=dble(j_extent)

    ! loop through non-default land covers
    do sp = 2, nos_species

        if (allocated(patch_centre_i)) then
           if (size(patch_centre_i) /= est_nos_patch(sp)) then
               deallocate(patch_centre_i,patch_centre_j)
               allocate(patch_centre_i(est_nos_patch(sp)),patch_centre_j(est_nos_patch(sp)))
           endif
        else
           allocate(patch_centre_i(est_nos_patch(sp)),patch_centre_j(est_nos_patch(sp)))
        endif

        if ( .not.weighted_patch .and. .not.ecotone) then

            ! find required number of random patch locations
            ! i direction first, which returns vector of length, with
            ! locations in original vector that meet the condition specified.
            ! Unused spaces in the vector have value -9999
            call which(input=sp_hold(:,1),length=(domain_extent_x*domain_extent_y),greaterlessequal=3 &
                     ,condition_num=dble_one)
            ! a module level variable generated by the which function contains
            ! information on how many variables are actually valid in the output
            ! and this is used to truncate the vector
            allocate(tmpi(nos_out),tmpj(nos_out))
            tmpi=which_out(1:nos_out) ; tmpj=which_out(1:nos_out)
            tmpj=tmpj/domain_extent_x ; tmpi=tmpi-dble(floor(tmpj)*domain_extent_x)
            tmpj=dble(ceiling(tmpj))
            call subsample(input=tmpi,length=size(tmpi),required=est_nos_patch(sp),replace=.false.)
            patch_centre_i = nint(tmpi(nint(subsample_locations_out)))
            ! j direction next
            patch_centre_j = nint(tmpj(nint(subsample_locations_out)))

            ! after final extraction we need to deallocate the tmpij variables,
            ! while tmp can remain as it will be the same total domain length
            deallocate(tmpi,tmpj)

        else if (.not.weighted_patch .and. ecotone ) then

            ! allocate required patch locations starting with the basic ecotone
            ! assumption of 50 %
            allocate(patch_locations(floor(domain_extent_x*domain_extent_y*0.5)))
            ! then assign them with a positional number
            do i = 1, floor(domain_extent_x*domain_extent_y*0.5)
               patch_locations(i)=i
            enddo
            ! setting initial ecotone assumptions
            sp_hold(patch_locations,1) = 0
            sp_hold(patch_locations,sp) = 1
            ! as patch_locations will be used again deallocate
            deallocate(patch_locations)

            ! find required number of random patch locations
            ! i direction first, which returns vector of length, with
            ! locations in original vector that meet the condition specified.
            ! Unused spaces in the vector have value -9999
            call which(input=sp_hold(:,1),length=(domain_extent_x*domain_extent_y),greaterlessequal=3 &
                      ,condition_num=dble_one)
            ! a module level variable generated by the which function contains
            ! information on how many variables are actually valid in the output
            ! and this is used to truncate the vector
            allocate(tmpi(nos_out),tmpj(nos_out))
            tmpi=which_out(1:nos_out) ; tmpj=which_out(1:nos_out)
            tmpj=tmpj/domain_extent_x ; tmpi=tmpi-dble(floor(tmpj)*domain_extent_x)
            tmpj=dble(ceiling(tmpj))
            call subsample(input=tmpi,length=size(tmpi),required=est_nos_patch(sp),replace=.false.)
            patch_centre_i = nint(tmpi(nint(subsample_locations_out)))
            ! j direction next
            patch_centre_j = nint(tmpj(nint(subsample_locations_out)))

            ! after final extraction we need to deallocate the tmpij variables,
            ! while tmp can remain as it will be the same total domain length
            deallocate(tmpi,tmpj)

        else 

            stop
            ! this should be a weighting approach see R code for additional
            ! ideas

        end if ! ecotone or not in selecting pixels

        i = 1 ; patch_done = 0
        check_interval_max = min(100,nint(est_nos_patch(sp)*0.5))
        check_interval = 2
        do while (patch_done < patch_area(sp)) 
            !# now work out the distance from each patch
            distance_array=sqrt(((i_extent_dble-patch_centre_i(i))*(i_extent_dble-patch_centre_i(i))) &
                               +((j_extent_dble-patch_centre_j(i))*(j_extent_dble-patch_centre_j(i))))*resolution

            !# which of these are within the lake area, assuming circular lake
            where (distance_array < patch_size) 
                   sp_hold(:,1) = 0 ; sp_hold(:,sp) = 1 
            end where

            ! only assess the total cover every so often to avoid a little
            ! computation
            if (mod(i,nint(check_interval)) == 0 .and. check_interval /= dble_one) then
               ! find out the total which has now had current species allocated
               tmp = 0
               where (sp_hold(:,sp) == 1) tmp = 1
               nos_out = sum(tmp)
               ! adjusted number of pixels by resolution
               patch_done = nos_out * resolution2
               ! then update the interval for checking (i.e. we check more often
               ! the closer we get)
               check_interval=(0.99-exp(0.99-(0.99/(patch_done/patch_area(sp))**5)))*check_interval_max
               check_interval=max(dble_one,check_interval)
            endif ! mod(i,check_interval)

            ! keep track of counter
            i = i + 1
            ! if we have ran out of patches we best now generate some more
            if (i == est_nos_patch(sp)) then
               ! first thing re set the counter back to the beginning of the vector
               i = 1

               ! find required number of random patch locations
               ! i direction first, which returns vector of length, with
               ! locations in original vector that meet the condition specified.
               ! Unused spaces in the vector have value -9999
               call which(input=sp_hold(:,1),length=(domain_extent_x*domain_extent_y),greaterlessequal=3 &
                         ,condition_num=dble_one)
               ! a module level variable generated by the which function contains
               ! information on how many variables are actually valid in the output
               ! and this is used to truncate the vector
               allocate(tmpi(nos_out),tmpj(nos_out))
               tmpi=which_out(1:nos_out) ; tmpj=which_out(1:nos_out)
               tmpj=tmpj/domain_extent_x ; tmpi=tmpi-(dble(floor(tmpj))*domain_extent_x)
               tmpj=dble(ceiling(tmpj))
               call subsample(input=tmpi,length=size(tmpi),required=est_nos_patch(sp),replace=.false.)
               patch_centre_i = anint(tmpi(nint(subsample_locations_out)))
               ! j direction next
               patch_centre_j = anint(tmpj(nint(subsample_locations_out)))

               ! after final extraction we need to deallocate the tmpij variables,
               ! while tmp can remain as it will be the same total domain length
               deallocate(tmpi,tmpj)

            endif ! if we have ran out of patches

        end do ! do while not enough patches specified

    enddo !# end species loop

    ! what area is covered by each species type
    tmp = 0
    where (sp_hold(:,1) == 1) tmp = 1
    nos_out = sum(tmp)
    if ( (nos_out*resolution2)/patch_area(1) > 1.15 .or. &
         (nos_out*resolution2)/patch_area(1) < 0.9) then
         print*,"Ratio of actual patch area : desired patch area ", &
               (nos_out*resolution2)/patch_area(1)," of species ",1
         print*,"sp =",1," have - desired ",(nos_out*resolution2)-patch_area(1)," m2"
         print*,"---fatal map error---"
    endif

    !# re-structure for further work and final checks
    if (.not.(allocated(mean_sp_abundance_at_distance))) then
        allocate(mean_sp_abundance_at_distance(domain_extent_x,domain_extent_y,nos_species))
    endif

    do sp = 1, nos_species
        ! then the restructuring bit
        counter = 1
        do j = 1, domain_extent_y
           mean_sp_abundance_at_distance(1:domain_extent_x,j,sp) = sp_hold(counter:(counter+domain_extent_x),sp)
           counter = counter+domain_extent_x
        enddo ! j/y dimension
    enddo ! species

    !# clean up
    deallocate(patch_centre_i,patch_centre_j)

    !# return back to user
    return

  end subroutine generate_species_map !# generate species map
  !
  !----------------------------------------------------------------------------
  !
  subroutine mask_lake_area(iter)

    implicit none    

    !# Function masks out the area covered by the lake in the domain.
    !# This is simply achieved by removing all species coverages in a given area

    ! declare arguments
    integer, intent(in) :: iter

    ! declare local variables
    double precision, dimension(domain_extent_x*domain_extent_y) :: i_extent_dble,j_extent_dble
    integer :: counter, tracker, i, j

    ! make dble to avoid repeated conversion
    i_extent_dble = dble(i_extent) ; j_extent_dble = dble(j_extent)

    if (.not.lake_fixed .or. iter == 1)  then
        !# now work out the distance from the lake
        distance_array=sqrt((((i_extent_dble-lake_centre_i)*(i_extent_dble-lake_centre_i)) &
                           + ((j_extent_dble-lake_centre_j)*(j_extent_dble-lake_centre_j))))*resolution
        !# which of these are within the lake area, assuming circular lake
        call which(input=distance_array,length=size(distance_array),greaterlessequal=2 &
                 ,condition_num=lake_size)
        if (allocated(water_locations)) then
            if (size(water_locations) /= nos_out) then
                deallocate(water_locations) ; allocate(water_locations(nos_out)) 
            endif
        else 
            allocate(water_locations(nos_out))
        endif
        water_locations=which_out(1:nos_out)
        !# keep track of how many
        nos_water=nos_out
    endif

    !# block out these values in our maps
    counter = 1 ; tracker = 1
    do j = 1, domain_extent_y
       do i = 1, domain_extent_x
          ! is our vectorised position in the array the same as the location
          ! specified as the next water position?
          if (counter == water_locations(tracker)) then
              ! is so make water across all species
              mean_sp_abundance_at_distance(i,j,:) = 0
              ! and then move on to the next tracker position in the water
              ! vector
              tracker = tracker + 1
          endif
          ! but always move through the counter position
          counter = counter + 1
       enddo
    end do

    !# return output back to the user
    return

  end subroutine mask_lake_area
  !
  !----------------------------------------------------------------------------
  !
  subroutine lake_distance_map

    implicit none

    !# function determines distance map for each water pixel

    ! declare local variables
    double precision, dimension(domain_extent_x*domain_extent_y) :: i_extent_dble,j_extent_dble
    double precision :: centre_i, centre_j
    integer :: water,j,counter

    ! make dble to avoid repeated conversion
    i_extent_dble = dble(i_extent) ; j_extent_dble = dble(j_extent)
    ! calculate distance arrays for each water pixel in turn
    do water = 1, nos_water

       !# locate new center, i.e. the water pixel
       centre_i=mod(water_locations(water),domain_extent_x)
       centre_j=ceiling(dble(water_locations(water))/dble(domain_extent_y))
       !# now create multiple distance maps for each of the water locations
       distance_array=sqrt(((i_extent_dble-centre_i)*(i_extent_dble-centre_i)) &
                         + ((j_extent_dble-centre_j)*(j_extent_dble-centre_j)))*resolution

       counter = 1
       do j = 1, domain_extent_y
          distance_from_lake(1:domain_extent_x,j,water) = distance_array(counter:(counter+domain_extent_x))
          counter = counter + domain_extent_x
       enddo ! do loop j

    enddo ! water loop
    !# ensure a minimum distance (0.1 m) in distance arrays
    where ( distance_from_lake .lt. 0.1 ) distance_from_lake = 0.1

    !# back to the user
    return

  end subroutine lake_distance_map !# lake_distance_map
  !
  !----------------------------------------------------------------------------
  !
  subroutine pollen_deposition_model
    
    implicit none

    ! function calculate the pollen deposition rate required for the pollen load
    ! modelling

    ! declare local variables
    integer :: dims_a(4),dims(3),i,w
    double precision :: beta_sp

    !# match matrix size with the species distribution information
    dims=shape(mean_sp_abundance_at_distance)
    if (allocated(pollen_deposition)) then
        dims_a = shape(pollen_deposition)
        if (dims_a(1) /= dims(1) .or. dims_a(2) /= dims(2) .or. &
            dims_a(3) /= dims(3) .or. dims_a(4) /= nos_water) then
            deallocate(pollen_deposition)
            allocate(pollen_deposition(dims(1),dims(2),dims(3),nos_water))
        endif
    else
        allocate(pollen_deposition(dims(1),dims(2),dims(3),nos_water))
    endif

    ! provide initial value
    pollen_deposition = 0
    if (deposition_model == 1) then
        !# deposition = 1/distance
        do i = 1, nos_species
           do w = 1, nos_water  
              pollen_deposition(:,:,i,w)=distance_from_lake(:,:,w)**(-dble_one)
           end do
        end do
    else if (deposition_model == 2) then
        !# deposition = 1/distance**2
        do i = 1, nos_species
           do w = 1, nos_water  
              pollen_deposition(:,:,i,w)=distance_from_lake(:,:,w)**(-2)
           end do
        end do
    else if (deposition_model == 3) then
        !# Presentice model (Sutton 1953; Prentice 1988)
        do i = 1, nos_species 
            if (fall_speed(i) == dble_zero) fall_speed(i) = (2*(pollen_radius(i)*pollen_radius(i))*gravity*(rho0-rho))*nine_mu_1
            beta_sp=(4*fall_speed(i))/(two_gamma1*wind_sp*sqrt_pi*diffusion_coefficient)
            pollen_deposition(:,:,i,:)=beta_sp*gamma1*(distance_from_lake**(gamma1-dble_one)) & 
                                      *exp(-beta_sp*(distance_from_lake**gamma1))
        enddo !# species loop
    endif !# all conditions

    !# return to user
    return

  end subroutine pollen_deposition_model !# pollen deposition model
  !
  !----------------------------------------------------------------------------
  !
  subroutine prentice_sugita(sp,out_var)

    ! declare input
    integer :: sp

    ! declare local variables
    integer :: w
    double precision, dimension(nos_water) :: out_var

    !# Prentice-Sugita model (Sugita 1994) modified by Bunting & Middleton (2005); 
    ! where pollen load is linearly related to distance weighted plant abundance
    do w = 1, nos_water 
       out_var(w)=sum((pollen_productivity(sp)*mean_sp_abundance_at_distance(:,:,sp)) &
                      *pollen_deposition(:,:,sp,w))
    end do

    !# return to the user
    return

  end subroutine prentice_sugita !# prentice_sugita
  !
  !----------------------------------------------------------------------------
  !
  subroutine iterative_loop (iter)

      ! this subroutine carries out the iterations commanded by the R f90
      ! interface. All variables should have been loaded into the model memory
      ! prior to coming into this subroutine.

      implicit none

      ! input variables
      integer, intent(in) :: iter ! current iteration 

      ! local variables
      double precision :: forest_productivity,savanna_productivity(9),nos_water_mult
      integer :: i, sp, water

      ! estimate inverse of area for each individual patch for this set up
      if (iter == 1) each_patch_area_1 = (pi * patch_size * patch_size) ** (-dble_one)

      ! determine what the map looks like
      if (random_or_read .and. (iter == 1 .or. .not.map_fixed .or. ecotone)) then

          if (.not.lake_fixed .or. iter == 1) then
              !# determine lake centre location
              if (lake_centred) then
                  lake_centre_i=floor(dble(domain_extent_x)*0.5)
                  lake_centre_j=floor(dble(domain_extent_y)*0.5)
              else if (lake_fixed .and. lake_specified) then
                  lake_centre_i=lake_user_x(1)
                  lake_centre_j=lake_user_y(1)
              else if (.not.lake_fixed .and. lake_specified ) then
                  lake_centre_i=lake_user_x(iter)
                  lake_centre_j=lake_user_y(iter)
              else if (.not.lake_fixed .and. .not.lake_centred .and. .not.lake_specified) then
                  lake_centre_i=nint(randn(0)*domain_extent_x)
                  lake_centre_j=nint(randn(0)*domain_extent_y)
              endif
          endif !# lake fixed or not

          ! ensure we have no carry over from last time
          if (iter == 1 .and. allocated(mean_sp_abundance_at_distance)) then
              deallocate(mean_sp_abundance_at_distance)
          endif
          !# create spatial domain of species info
          call generate_species_map(iter)
          !# update species map for lake area
          call mask_lake_area(iter)

      else if (.not.random_or_read) then

          ! then we should have already done something about this....

      endif ! random_or_read and map_fixed or ecotone

      if (.not.lake_fixed .or. iter == 1) then
          !# calculate distance maps for each water pixels
          !# serial approach
          if (allocated(distance_from_lake)) deallocate(distance_from_lake)
          allocate(distance_from_lake(domain_extent_x,domain_extent_y,nos_water))
          call lake_distance_map
          call pollen_deposition_model

      endif ! fixed lake condition

      if (nos_iter > 1 .and. .not.random_or_read) then
          !# then we assume we are calibrating, we assume a standard calibration
          !range with interval determine by the number of iterations.
          !# This is an assumption currentl and relies on the user subselecting
          !the output for realistic relative combinations of parameters
          forest_productivity=1
          savanna_productivity(1)=0.0675 ; savanna_productivity(2)=0.125
          savanna_productivity(3)=0.25   ; savanna_productivity(4)=0.5
          savanna_productivity(5)=1.0     ; savanna_productivity(6)=2.0
          savanna_productivity(7)=4.0     ; savanna_productivity(8)=8.0
          savanna_productivity(9)=16.0
          ! load specific values
          pollen_productivity(1)=forest_productivity
          pollen_productivity(2)=savanna_productivity(iter)
      endif

      !# Prentice-Sugita model (Sugita 1994) modified by Bunting & Middleton
      !(2005); where pollen load is linearly related to distance weighted plant
      !abundance
      ! declare output field first
      if (iter == 1 .and. allocated(pollen_load_at_lake_prediction)) deallocate(pollen_load_at_lake_prediction)
      if (.not.allocated(pollen_load_at_lake_prediction)) allocate(pollen_load_at_lake_prediction(nos_water,nos_species))
      nos_water_mult = dble(nos_water)**(-dble_one)
      do sp = 1, nos_species
         call prentice_sugita(sp,pollen_load_at_lake_prediction(1:nos_water,sp))
         sum_pol(sp,iter) = sum(pollen_load_at_lake_prediction(1:nos_water,sp)) 
      end do
      ! average across water pixels
      sum_pol = sum_pol * nos_water_mult 

      !# return to user
      return

  end subroutine iterative_loop
  !
  !----------------------------------------------------------------------------
  !
!
!------------------------------------------------------------------------------
!
end module pollen_deposition_model_mod

