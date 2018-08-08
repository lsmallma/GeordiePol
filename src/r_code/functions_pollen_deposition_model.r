
###
## Script contains all source code for the deposition model functions
## If you are not familiar with R please leave this alone...
###

if(use_parallel) {require(parallel)}
require(fields)

generate_species_map<- function (lake_centre_i,lake_centre_j,i_extent,j_extent,domain_extent_x,domain_extent_y,resolution,ecotone,sp_fractions,patch_size,nos_species,weighted_patch) {

    # function generates species maps and rapidly carries out random / random weighted clearance

    # estimate how many forest patches will be needed  (m)
    total_area=domain_extent_x*domain_extent_y*resolution**2
    # how much area in total will become each species
    if (ecotone) {
	patch_area=c(1,2)
  	patch_area[1]=total_area*0.5*sp_fractions[1]*1e-2
  	patch_area[2]=total_area*(0.5+(0.5*sp_fractions[2]*1e-2))
        # multiple by 2 as random allocation could mean that we overlap sometimes
        est_nos_patch=ceiling(patch_area/(pi*patch_size*patch_size))
    } else {
	patch_area=total_area*sp_fractions*1e-2
        # multiple by 2 as random allocation could mean that we overlap sometimes
        est_nos_patch=ceiling(patch_area/(pi*patch_size*patch_size))*2
    }
    # block out these values in our maps
    sp_hold=array(1,dim=c(domain_extent_x*domain_extent_y,nos_species))
    # assume that sp 1 is blanket coverage first
    sp_hold[,2:nos_species]=0

    for (sp in seq(2, nos_species)) {

	if ( weighted_patch == FALSE & ecotone == FALSE) {
	    # find required number of random patch locations
            tmp=sample(which(sp_hold[,1] == 1),est_nos_patch[sp])
            patch_centre_j=tmp/domain_extent_x
            patch_centre_i=tmp-(floor(patch_centre_j)*domain_extent_x)
            patch_centre_j=ceiling(patch_centre_j)

	} else if (weighted_patch == FALSE & ecotone == TRUE ) {
            # apply ecotone first
            patch_locations=1:(floor(dim(sp_hold)[1]*0.5))
            sp_hold[patch_locations,1]=0
            sp_hold[patch_locations,sp]=1

	    # find required number of patch locations
            tmp=sample(which(sp_hold[,1] == 1),est_nos_patch[sp])
            patch_centre_j=tmp/domain_extent_x
            patch_centre_i=tmp-(floor(patch_centre_j)*domain_extent_x)
            patch_centre_j=ceiling(patch_centre_j)
	} else {
	    stop('')
	    # find required number of random-weighted patch locations
	    patch_centre_i=round(pmax(1,pmin(domain_extent_x,rnorm(est_nos_patch,lake_centre_i,domain_extent_x/4))), digits=0)
	    patch_centre_j=round(pmax(1,pmin(domain_extent_y,rnorm(est_nos_patch,lake_centre_j,domain_extent_y/4))), digits=0)
	}

	i = 1 ; patch_done = 0
	while (patch_done < patch_area[sp]) {
	    # now work out the distance from each patch
	    distance_array=sqrt((abs(i_extent-patch_centre_i[i])**2)+(abs(j_extent-patch_centre_j[i])**2))*resolution
	    # which of these are within the lake area, assuming circular lake
	    patch_locations=which(distance_array < patch_size)
	    sp_hold[patch_locations,1]=0
	    sp_hold[patch_locations,sp]=1
	    patch_done = length(which(sp_hold[,sp] == 1))*resolution**2
	    i = i + 1 
            if (i == est_nos_patch[sp]) {
               # we will generate some new variables
               # first reset counter
               i = 1

               tmp=sample(which(sp_hold[,1] == 1),est_nos_patch[sp])
               patch_centre_j=tmp/domain_extent_x
               patch_centre_i=tmp-(floor(patch_centre_j)*domain_extent_x)
               patch_centre_j=ceiling(patch_centre_j)
            
	    } # end of if regenerate
	} # end of while
	
    } # end species loop

    # final checks
    for (sp in seq(1, nos_species)) {
	if ( (length(which(sp_hold[,sp] == 1))*resolution**2)/patch_area[sp] > 1.15 | (length(which(sp_hold[,sp] == 1))*resolution**2)/patch_area[sp] < 0.9) {
	    print(paste("Ratio of actual patch area : desired patch area ",(length(which(sp_hold[,sp] == 1))*resolution**2)/patch_area[sp]," of species ",sp,sep=""))
	    stop('---Fatal map error---')
	}
    }

    # re-structure for further work
    mean_sp_abundance_at_distance=array(sp_hold, dim=c(domain_extent_x,domain_extent_y,nos_species))

    # clean up
    rm(sp_hold,total_area,est_nos_patch,patch_centre_i,patch_centre_j,patch_locations,distance_array)

    # return back to user
    return(mean_sp_abundance_at_distance)

} # generate species map

mask_lake_area<-function(mean_sp_abundance_at_distance,lake_centre_i,lake_centre_j,i_extent,j_extent,lake_size,domain_extent_x,domain_extent_y,resolution,nos_species) {
    
    # Function masks out the area covered by the lake in the domain.
    # This is simply achieved by removing all species coverages in a given area

    # now work out the distance from the lake
    distance_array=sqrt((abs(i_extent-lake_centre_i)**2)+(abs(j_extent-lake_centre_j)**2))*resolution
    # which of these are within the lake area, assuming circular lake
    water_locations=which(distance_array < lake_size)
    # keep track of how many
    nos_water=length(water_locations)
    # block out these values in our maps
    sp_hold=array(1,dim=c(domain_extent_x*domain_extent_y,nos_species))
    for (sp in seq(1,nos_species)) { sp_hold[,sp]=as.vector(mean_sp_abundance_at_distance[,,sp]) }
    sp_hold[water_locations,]=0
    # re-structure for further work
    mean_sp_abundance_at_distance=array(sp_hold, dim=c(domain_extent_x,domain_extent_y,nos_species))

    # clean up
    rm(sp_hold,distance_array)

    # return output back to the user
    return(list(mean_sp_abundance_at_distance=mean_sp_abundance_at_distance,nos_water=nos_water,water_locations=water_locations))
}

lake_distance_map<- function(water,nos_water,water_locations,domain_extent_x,domain_extent_y,resolution,i_extent,j_extent) {

    # function determines distance map for each water pixel

    # locate new center, i.e. the water pixel
    centre_i=water_locations[water]%%domain_extent_x ; centre_j=ceiling(water_locations[water]/domain_extent_y)
    # now create multiple distance maps for each of the water locations
    distance_array=sqrt((abs(i_extent-centre_i)**2)+(abs(j_extent-centre_j)**2))*resolution

    # ensure a minimum distance (1 m) in distance arrays
    distance_array[which(distance_array < 0.1)] = 0.1

    # clean up
    rm(centre_i,centre_j)

    # back to the user
    return(distance_array)

} # lake_distance_map

pollen_deposition_model<-function(nos_water,nos_species,deposition_model,distance_array,mean_sp_abundance_at_distance,prentice_parameters) {

  # unload values
  fall_speed = prentice_parameters$fall_speed       ; gravity = prentice_parameters$gravity
  pollen_radius = prentice_parameters$pollen_radius ; mu = prentice_parameters$mu
  rho0 = prentice_parameters$rho0                   ; diffusion_coefficient = prentice_parameters$diffusion_coefficient
  rho = prentice_parameters$rho                     ; gamma1 = prentice_parameters$gamma1
  wind_sp = prentice_parameters$wind_sp

  # match matrix size with the species distribution information
  pollen_deposition=array(NA, dim=c(dim(mean_sp_abundance_at_distance),nos_water))
  if (deposition_model==1) {
      # deposition = 1/distance
      for (w in seq(1, nos_water)) { pollen_deposition[,,1:nos_species,w]=distance_array[,,w]**-1 }
  } else if (deposition_model==2) {
      # deposition = 1/distance**2
      for (w in seq(1, nos_water)) { pollen_deposition[,,1:nos_species,w]=distance_array[,,w]**-2 }
  } else if (deposition_model==3) {
      # Presentice model (Sutton 1953; Prentice 1988)
      for (i in seq(1, nos_species)) {
	  if (fall_speed[i] == 0) {fall_speed[i]=(2*(pollen_radius**2)*gravity*(rho0-rho))/(9*mu)}
	  beta_sp=(4*fall_speed[i])/((2*gamma1)*wind_sp*(pi**0.5)*diffusion_coefficient)
	  pollen_deposition[,,i,]=beta_sp*gamma1*(distance_array[,,]**(gamma1-1))*exp(-beta_sp*distance_array[,,]**gamma1)
      } # species loop
  } # all conditions

  # return to user
  return(pollen_deposition)

} # pollen deposition model

prentice_sugita<-function (sp,pollen_productivity,mean_sp_abundance_at_distance,pollen_deposition,nos_water) {

    # Prentice-Sugita model (Sugita 1994) modified by Bunting & Middleton (2005); where pollen load is linearly related to distance weighted plant abundance
    pollen_load_at_lake_prediction=array(NA, dim=nos_water)
    for (w in seq(1, nos_water)) {
	pollen_load_at_lake_prediction[w]=sum(pollen_productivity[sp]*mean_sp_abundance_at_distance[,,sp]*pollen_deposition[,,sp,w])
    }

    # return to the user
    return(pollen_load_at_lake_prediction)

} # prentice_sugita

read_map<-function(input_filename,input_type){
  # function is responsible for reading in the designated file and extracting the relevant spatial / map information
  # currently we assume that we have only two species

  if (input_type == "netcdf") {
      require(ncdf)
      map_file=open.ncdf(input_filename)
      if (length(which(names(map_file$var) == "forest_cover")) == 0) {stop("input netcdf file does not contain 'forest_cover' variable")}
      forest_cover=get.var.ncdf(map_file,"forest_cover")
  } else if (input_type == "tiff") {
      require(raster)
      forest_cover_tiff=raster(input_filename)
      # smash up type of object
      dims=dim(forest_cover_tiff)
      forest_cover=as.vector(forest_cover_tiff)
      forest_cover=array(forest_cover, dim=c(dims[2],dims[1]))
      forest_cover=forest_cover[,rev(1:dim(forest_cover)[2])]
      rm(forest_cover_tiff)
  } else {
     stop("have not specified either 'netcdf' or 'tiff' file types for reading input")
  }
  
  # now how big is our spatial map
  domain_extent_x = dim(forest_cover)[1]
  domain_extent_y = dim(forest_cover)[2]
  # where is the water
  water_locations = which(as.vector(forest_cover == -100))
  # how much water is there
  nos_water = length(water_locations)
    
  # generate two species array
  mean_sp_abundance_at_distance=array(forest_cover,dim=c(dim(forest_cover)[1],dim(forest_cover)[2],2))
  # assume data is % cover and we turn it into fractions
  mean_sp_abundance_at_distance=mean_sp_abundance_at_distance*1e-2
  # assume where we have non-zero number that savanna is 1-forest value
  mean_sp_abundance_at_distance[,,2]=1-mean_sp_abundance_at_distance[,,2]
  # mask water areas
  forest_cover=array(NA,dim=c(domain_extent_x*domain_extent_y,2))
  for (sp in seq(1,2)) { forest_cover[,sp]=as.vector(mean_sp_abundance_at_distance[,,sp]) }
  forest_cover[water_locations,]=0
  # re-structure for further work
  mean_sp_abundance_at_distance=array(forest_cover, dim=c(domain_extent_x,domain_extent_y,2))
  
  # that is everything, time to return to the main function
  output=list(domain_extent_x=domain_extent_x,domain_extent_y=domain_extent_y
	     ,water_locations=water_locations,nos_water=nos_water
	     ,mean_sp_abundance_at_distance=mean_sp_abundance_at_distance)
  
  # return values
  return(output)
  
} # read_map

iterative_loop<-function(iter,model_info,i_extent,j_extent) {

      # extract model information
      deposition_model=model_info$deposition_model           ; nos_iter=model_info$nos_iter
      domain_extent_x=model_info$domain_extent_x             ; domain_extent_y=model_info$domain_extent_y
      random_or_read=model_info$random_or_read               ; ecotone=model_info$ecotone
      resolution=model_info$resolution                       ; weighted_patch=model_info$weighted_patch 
      map_fixed=model_info$map_fixed                         ; show_map=model_info$show_map
      nos_species=model_info$nos_species                     ; sp_fractions=model_info$sp_fractions 
      patch_size=model_info$patch_size                       ; lake_size=model_info$lake_size 
      lake_centred=model_info$lake_centred                   ; lake_fixed=model_info$lake_fixed
      lake_specified=model_info$lake_specified               ; lake_user_x=model_info$lake_user_x 
      lake_user_y=model_info$lake_user_y                     ; wind_sp=model_info$wind_sp 
      fall_speed=model_info$fall_speed                       ; pollen_radius=model_info$pollen_radius
      pollen_productivity=model_info$pollen_productivity     ; gamma1=model_info$gamma1
      diffusion_coefficient=model_info$diffusion_coefficient ; gravity=model_info$gravity 
      pollen_deposition=model_info$pollen_deposition         ; lake_centre_i=model_info$lake_centre_i
      lake_centre_j=model_info$lake_centre_j                 ; mean_sp_abundance_at_distance=model_info$mean_sp_abundance_at_distance
      nos_water=model_info$nos_water                         ; water_locations=model_info$water_locations
      rho0=model_info$rho0 ; rho=model_info$rho ; mu=model_info$mu

      if (random_or_read == "random" & (map_fixed == FALSE | ecotone == TRUE)) {

	  if (lake_fixed == FALSE) {
	      # determine lake centre location
	      if (lake_centred) {
		  lake_centre_i=floor(domain_extent_x/2)
		  lake_centre_j=floor(domain_extent_y/2)
	      } else if (lake_fixed == TRUE & lake_specified == TRUE) {
		  lake_centre_i=lake_user_x[1]
		  lake_centre_j=lake_user_y[1]
	      } else if (lake_fixed == FALSE & lake_specified == TRUE) {
		  lake_centre_i=lake_user_x[iter]
		  lake_centre_j=lake_user_y[iter]
	      } else if (lake_fixed == FALSE & lake_centred == FALSE & lake_specified == FALSE) {
		  lake_centre_i=round(runif(1,1,domain_radius*2), digits=0)
		  lake_centre_j=round(runif(1,1,domain_radius*2), digits=0)
	      }
	  } # lake fixed or not


	  # create spatial domain of species info
	  mean_sp_abundance_at_distance=generate_species_map(lake_centre_i,lake_centre_j,i_extent,j_extent,domain_extent_x,domain_extent_y,resolution,ecotone,sp_fractions,patch_size,nos_species,weighted_patch)
	  

	  # update species map for lake area
	  output=mask_lake_area(mean_sp_abundance_at_distance,lake_centre_i,lake_centre_j,i_extent,j_extent,lake_size,domain_extent_x,domain_extent_y,resolution,nos_species)
	  mean_sp_abundance_at_distance=output$mean_sp_abundance_at_distance
	  nos_water=output$nos_water
	  water_locations=output$water_locations

      } else if (random_or_read == "read") {
	  # we don't do anything here as we assume that we have only one map to read in and we already did that in the previous function
      }
      
      if (lake_fixed == FALSE) {
	  # calculate distance maps for each water pixels
	  # serial approach
	  distance_array=lapply(1:nos_water,FUN=lake_distance_map,nos_water=nos_water,water_locations=water_locations,domain_extent_x=domain_extent_x,domain_extent_y=domain_extent_y,resolution=resolution,i_extent=i_extent,j_extent=j_extent)
	  # extract from lists the finished array
	  distance_array=array(unlist(distance_array),dim=c(domain_extent_x,domain_extent_y,nos_water))

	  # important to periodically clean up
	  rm(output) ; gc(reset=TRUE, verbose=FALSE)

	  # calculate pollen deposition
	  prentice_parameters=list(pollen_radius=pollen_radius,fall_speed=fall_speed,gravity=gravity,rho0=rho0,rho=rho,mu=mu,gamma1=gamma1,diffusion_coefficient=diffusion_coefficient,wind_sp=wind_sp)
	  pollen_deposition=pollen_deposition_model(nos_water,nos_species,deposition_model,distance_array,mean_sp_abundance_at_distance,prentice_parameters)
	  rm(prentice_parameters)
      } # lake fixed condition

      if ((iter == 1 | iter == nos_iter) & show_map) {
	  if (iter == 1) {
	      jpeg(file=paste(output_path,"/timeseries_composit_",project_name,".jpg",sep=""), width=7200, height=4000, res=400, quality=100)
	      par(mfrow=c(2,nos_species))
	  }
	  for (sp in seq(1,nos_species)) {
	      image.plot(mean_sp_abundance_at_distance[,,sp], main=paste("Species ",sp," abundance (%)",sep=""))
	  }
	  if (iter == nos_iter) {dev.off()}
      }

      if (nos_iter > 1 & random_or_read == "read") {
	  # then we assume we are calibrating, we assume a standard calibration range with interval determine by the number of iterations.
	  # This is an assumption currentl and relies on the user subselecting the output for realistic relative combinations of parameters
	  forest_productivity=rep(1,times=nos_iter)
	  savanna_productivity=c(0.0675,0.125,0.25,0.5,1,2,4,8,16)
	  pollen_productivity=c(forest_productivity[iter],savanna_productivity[iter])
      }
      # Prentice-Sugita model (Sugita 1994) modified by Bunting & Middleton (2005); where pollen load is linearly related to distance weighted plant abundance
      pollen_load_at_lake_prediction=lapply(1:nos_species,FUN=prentice_sugita,pollen_productivity=pollen_productivity,mean_sp_abundance_at_distance=mean_sp_abundance_at_distance,pollen_deposition=pollen_deposition,nos_water=nos_water)
      pollen_load_at_lake_prediction=t(array(unlist(pollen_load_at_lake_prediction), dim=c(nos_water,nos_species)))
      # number of pollen grains expected; mean across all water pixels
      sum_pol=rowMeans(pollen_load_at_lake_prediction)

      # if this is not a recal run then drop in the intermediate file
      if (model_info$recal_run == FALSE) {write(sum_pol,file=model_info$intermediate_file,ncolumns=nos_species,append=TRUE, sep=",")}

      # garbage collection
      gc(reset=TRUE, verbose=FALSE)

      # return to user
      return(sum_pol)

} # end of iterative loop

call_f90_code<-function(group,model_info,i_extent,j_extent,nositer_group,wdir){

      # call the fortran code
      dyn.load(paste(wdir,"/pollen.so",sep=""))
      # minor adjustements for compatability; f90 cannot recieve logical values, they are transformed into integers, so it makes sense to be explicit here instead
      if (model_info$random_or_read == "random") {random_or_read_in = 1} else {random_or_read_in = 0}
      if (model_info$ecotone) {ecotone_in = 1} else {ecotone_in = 0}
      if (model_info$weighted_patch) {weighted_patch_in = 1 } else {weighted_patch_in = 0}
      if (model_info$map_fixed) {map_fixed_in = 1} else {map_fixed_in = 0}
      if (model_info$lake_centred) {lake_centred_in = 1} else {lake_centred_in = 0}
      if (model_info$lake_fixed) {lake_fixed_in = 1} else {lake_fixed_in = 0}
      if (model_info$lake_specified) {lake_specified_in = 1} else {lake_specified_in = 0}
      if (model_info$recal_run) {recal_run_in = 1} else {recal_run_in = 0}

      # call the actual interface function
      tmp=.Fortran("rtof90interface",tmp_file=as.character(model_info$intermediate_file),tmp_char=as.integer(model_info$outlength)
                                    ,recal_run_in=as.integer(recal_run_in)
                                    ,group=as.integer(group),nos_iter_in=as.integer(nositer_group[group])
                                    ,deposition_model_in=as.integer(model_info$deposition_model)
                                    ,domain_extent_x_in=as.integer(model_info$domain_extent_x)
                                    ,domain_extent_y_in=as.integer(model_info$domain_extent_y)
                                    ,resolution_in=as.double(model_info$resolution)
                                    ,random_or_read_in=as.integer(random_or_read_in)
                                    ,ecotone_in=as.integer(ecotone_in)
                                    ,weighted_patch_in=as.integer(weighted_patch_in)
                                    ,map_fixed_in=as.integer(map_fixed_in)
                                    ,nos_species_in=as.integer(model_info$nos_species)
                                    ,sp_fraction_in=as.double(model_info$sp_fractions)
                                    ,patch_size_in=as.double(model_info$patch_size)
                                    ,lake_size_in=as.double(model_info$lake_size)
                                    ,lake_centred_in=as.integer(lake_centred_in)
                                    ,lake_fixed_in=as.integer(lake_fixed_in)
                                    ,lake_specified_in=as.integer(lake_specified_in)
                                    ,nos_lake_in=as.integer(model_info$nos_lake)
                                    ,nos_water_in=as.integer(model_info$nos_water)
                                    ,water_locations_in=array(model_info$water_locations, dim=model_info$nos_water)
                                    ,lake_user_x_in=as.integer(model_info$lake_user_x)
                                    ,lake_user_y_in=as.integer(model_info$lake_user_y)
                                    ,wind_sp_in=as.double(model_info$wind_sp)
                                    ,fall_speed_in=as.double(model_info$fall_speed)
                                    ,pollen_radius_in=as.double(model_info$pollen_radius)
                                    ,pollen_productivity_in=as.double(model_info$pollen_productivity)
                                    ,i_extent_in=as.integer(i_extent)
                                    ,j_extent_in=as.integer(j_extent)
                                    ,mean_sp_abundance_at_distance_in=array(model_info$mean_sp_abundance_at_distance,dim=c(model_info$domain_extent_x,model_info$domain_extent_y,model_info$nos_species))
                                    ,sum_pol_out=array(0,dim=c(model_info$nos_species,nositer_group[group])))
      # extract and restructure output variable
      sum_pol=tmp$sum_pol_out
      sum_pol=array(sum_pol, dim=c(model_info$nos_species,nositer_group[group]))
      # upload the f90 shared library
      dyn.unload(paste(wdir,"/pollen.so",sep=""))
      # clean up
      rm(tmp) ; gc(reset=TRUE, verbose=FALSE)
      # and return
      return(sum_pol)

} # end of call_f90_code

run_model<-function() {

  stime=proc.time()["elapsed"] ; print(Sys.time())

  # declare some assumed constants, these are the same values as used in f90
  # but remember to update them...
#  pi = 3.1415926535            # ! pi is already built into R
  gravity = 9.80665            #& ! acceleration due to gravity (m.s-2)
  mu = 1.8e-2                  #& ! dynamic viscosity (g.m-1.s-1)
  rho0 = 2e6                   #& ! particle density (g.m-3)
  diffusion_coefficient = 0.12 #& ! (0.30) vertical diffusion coefficient (m^1/8)
  rho = 1.27e6                 #& ! air density (g.m-3)
  gamma1 = 0.125               #  ! turbulence parameter (2*gamma1)

  # check options
  if (ecotone & lake_fixed == FALSE & (nos_iter != length(lake_user_x) | nos_iter != length(lake_user_y))) {
      stop('number of iterations must equal number lakes provided when running moving lake analysis in relation to ecotones')
  }
  if (ecotone) {
      print("NOTE: using ecotone assumes user defined fractions are in ADDITION to assumed 50 % loss")
  }
  if (sum(sp_fractions) != 100) {
      stop('Species fractions provided to not equal 100 %')
  }
  if (random_or_read == "read" & ecotone == TRUE) {
      stop('random_or_read and ecotone options are incompatable')
  }
  if (ecotone == TRUE & nos_species > 2) {
      stop('ecotone function does not currently work with more than two species')
  }
  if (lake_centred & lake_specified) {
      stop('both lake_centred and lake_specified cannot be true...pick one!')
  }
  if (random_or_read == "read" & (lake_fixed == FALSE | map_fixed == FALSE)) {
      stop('lake_fixed and map_fixed must equal TRUE when reading in the map or this doesnt work')
  }
  if (random_or_read == "read" & nos_species > 2) {
      stop('the option to read in a map is currently restricted to two species only, sorry will get to it later')
  }
  recal_run=FALSE
  if (random_or_read == "read" & nos_iter > 1 ) {
      print('when reading a map and iterations > 1 the model assumes this is a calibration attempt, therefore nos_iter will be forced == 9')
      nos_iter = 9 ; recal_run=TRUE
  }

  # if the lake is centred, by default is must also be fixed!
  if (lake_centred & lake_fixed == FALSE) {print("lake_fixed is currently set to FALSE while the lake_centred is TRUE. lake_fixed will be ignored") ; lake_fixed=TRUE}

  # if the lake is not fixed the only time we read in more than one lake is if it is user specified and this must be the same number as iterations
  if (lake_fixed == FALSE & lake_specified == TRUE) {
      nos_lake = nos_iter
  } else {
      nos_lake = 1	
  }

  # check whether this is a restart run, first we need to check if a intermediate file exists
  intermediate_file=paste(output_path,"/",project_name,"_tmp.csv",sep="")
  outlength=nchar(intermediate_file, type = "chars", allowNA = FALSE)
  if (outlength > 254) {
      stop(paste("output_path and project_name combo needs to be ",outlength-254," characters shorter",sep=""))
  }
  outfile=paste(output_path,"/",project_name,"_done.csv",sep="")
  restart = FALSE ; restart_iterations = 0
  if (recal_run == FALSE & file.exists(intermediate_file)){
     # first off this is a restart run
     restart=TRUE
     tmp=read.csv(file=intermediate_file, header=TRUE)
     if (dim(tmp)[1] < 2) {
        # if only just started before then assume this is a fresh run
        restart=FALSE
     } else {
        dims=dim(tmp) ; restart_sum_pol=array(unlist(tmp),dim=dims)
        rm(tmp,dims)     
        restart_iterations=dim(restart_sum_pol)[1]
        # adjust nos_iter to now be the difference we need
        nos_iter = nos_iter-restart_iterations
        # inform the user
        print("...this is a restart run...")
     }
  } else if (recal_run == FALSE & file.exists(intermediate_file) == FALSE) {
     # well if this is not a restart then we better set up the restart file
     write(paste("sp",1:nos_species,sep=""), file=intermediate_file, ncolumns=nos_species, append=FALSE, sep=",")
  } #  is restart run?

     # we do things very differently when we randomly create a map versus reading one in for training
     if (random_or_read == "random") {
         
         # create vectorised versions of the 2D space
         i_extent=rep(1:(domain_extent_x),times=(domain_extent_y)) 
         j_extent=rep(1:(domain_extent_y), each=(domain_extent_x))

         # determine lake centre location
         if (lake_centred) {
             lake_centre_i=floor(domain_extent_x/2)
             lake_centre_j=floor(domain_extent_y/2)
         } else if (lake_fixed == TRUE & lake_specified == TRUE) {
	     lake_centre_i=lake_user_x[1]
	     lake_centre_j=lake_user_y[1]
         } else if (lake_fixed == FALSE & lake_specified == TRUE) {
   	     # this will be iterated through the function itself
	     lake_centre_i=lake_user_x[1]
	     lake_centre_j=lake_user_y[1]
         } else if (lake_fixed == FALSE & lake_centred == FALSE & lake_specified == FALSE) {
   	     lake_centre_i=round(runif(1,1,domain_extent_x), digits=0)
	     lake_centre_j=round(runif(1,1,domain_extent_y), digits=0)
         }

         # create spatial domain of species info
         mean_sp_abundance_at_distance=generate_species_map(lake_centre_i,lake_centre_j,i_extent,j_extent,domain_extent_x,domain_extent_y,resolution
                                                           ,ecotone,sp_fractions,patch_size,nos_species,weighted_patch)

         # update species map for lake area
         output=mask_lake_area(mean_sp_abundance_at_distance,lake_centre_i,lake_centre_j,i_extent,j_extent,lake_size,domain_extent_x,domain_extent_y,resolution,nos_species)
         mean_sp_abundance_at_distance=output$mean_sp_abundance_at_distance
         nos_water=output$nos_water
         water_locations=output$water_locations

     } else if (random_or_read == "read") {

         output=read_map(input_filename,input_type)
         domain_extent_x=output$domain_extent_x
         domain_extent_y=output$domain_extent_y
         mean_sp_abundance_at_distance=output$mean_sp_abundance_at_distance
         nos_water=output$nos_water
         water_locations=output$water_locations      
      
         # set some dummy values for passing in the model
         lake_centre_i = NA ; lake_centre_j = NA

         # create vectorised versions of the 2D space
         i_extent=rep(1:(domain_extent_x),times=(domain_extent_y))
         j_extent=rep(1:(domain_extent_y), each=(domain_extent_x))

     }
     # important to periodically clean up
     rm(output) ; gc(reset=TRUE, verbose=FALSE)

     # calculate distance maps for each water pixels
     # serial approach
     distance_array=lapply(1:nos_water,FUN=lake_distance_map,nos_water=nos_water,water_locations=water_locations
                                                            ,domain_extent_x=domain_extent_x,domain_extent_y=domain_extent_y
                                                            ,resolution=resolution,i_extent=i_extent,j_extent=j_extent)
     # extract from lists the finished array
     distance_array=array(unlist(distance_array),dim=c(domain_extent_x,domain_extent_y,nos_water))

     # calculate pollen deposition
     prentice_parameters=list(pollen_radius=pollen_radius,fall_speed=fall_speed,gravity=gravity,rho0=rho0,rho=rho,mu=mu,gamma1=gamma1
                             ,diffusion_coefficient=diffusion_coefficient,wind_sp=wind_sp)
     pollen_deposition=pollen_deposition_model(nos_water,nos_species,deposition_model,distance_array,mean_sp_abundance_at_distance,prentice_parameters)

     if (show_deposition_curve & output_path != " " & project_name != " ") {
         col_list=c("black","red","blue","green","pink","yellow","purple")
         jpeg(file=paste(output_path,"/",project_name,"_","deposition_distance_relationship.jpg",sep=""), width=7200, height=4000, res=400, quality=100)
         par(mfrow=c(1,1), omi=c(0.2,0.2,0.1,0.1), mai=c(1.0,1.0,0.1,0.1))
         for (sp in seq(1, dim(pollen_deposition)[3])) {
   	     if (sp == 1) {
	         plot(c(1:length(lake_centre_i:dim(pollen_deposition)[1]))*resolution,pollen_deposition[lake_centre_i:dim(pollen_deposition)[1],lake_centre_j,sp,1]/max(as.vector(pollen_deposition[lake_centre_i:dim(pollen_deposition)[1],lake_centre_j,sp,1])), ylim=c(0,1), type="l", lwd=4,ylab="Deposition as fraction of maximum rate", xlab="Distance from lake (m)",cex.axis=2.0,cex.lab=2.0)
   	     } else {
	         lines(c(1:length(lake_centre_i:dim(pollen_deposition)[1]))*resolution,pollen_deposition[lake_centre_i:dim(pollen_deposition)[1],lake_centre_j,sp,1]/max(as.vector(pollen_deposition[lake_centre_i:dim(pollen_deposition)[1],lake_centre_j,sp,1])), col=col_list[sp],lwd=4)
   	     }
	      abline(0.95,0, col="grey", lwd=3)
         }
         dev.off()
         jpeg(file=paste(output_path,"/",project_name,"_","cumulative_deposition_distance_relationship.jpg",sep=""), width=7200, height=4000, res=400, quality=100)
         par(mfrow=c(1,1), omi=c(0.2,0.2,0.1,0.1), mai=c(1.0,1.0,0.1,0.1))
         for (sp in seq(1, dim(pollen_deposition)[3])) {
   	     if (sp == 1) {
	         plot(c(1:length(lake_centre_i:dim(pollen_deposition)[1]))*resolution,(cumsum(pollen_deposition[lake_centre_i:dim(pollen_deposition)[1],lake_centre_j,sp,1])/sum(as.vector(pollen_deposition[lake_centre_i:dim(pollen_deposition)[1],lake_centre_j,sp,1])))*100, ylim=c(0,100), type="l", lwd=4,ylab="Cumulative pollen contribution (%)", xlab="Distance from lake (m)",cex.axis=2.0,cex.lab=2.0)
   	     } else {
	         lines(c(1:length(lake_centre_i:dim(pollen_deposition)[1]))*resolution,(cumsum(pollen_deposition[lake_centre_i:dim(pollen_deposition)[1],lake_centre_j,sp,1])/sum(as.vector(pollen_deposition[lake_centre_i:dim(pollen_deposition)[1],lake_centre_j,sp,1])))*100, col=col_list[sp],lwd=4)
   	     }
	     abline(95,0, col="grey", lwd=3)
         }
         dev.off()
     }

     rm(prentice_parameters)
     # important to periodically clean up
     gc(reset=TRUE, verbose=FALSE)

     # load model information / options into list to be passed to function
     model_info=list(deposition_model=deposition_model,nos_iter=nos_iter
   		    ,domain_extent_x=domain_extent_x,domain_extent_y=domain_extent_y
		    ,resolution=resolution,random_or_read=random_or_read,ecotone=ecotone
                    ,weighted_patch=weighted_patch,map_fixed=map_fixed,show_map=show_map
		    ,nos_species=nos_species,sp_fractions=sp_fractions,patch_size=patch_size
		    ,lake_size=lake_size,lake_centred=lake_centred,lake_fixed=lake_fixed
		    ,lake_specified=lake_specified,lake_user_x=lake_user_x,lake_user_y=lake_user_y
		    ,wind_sp=wind_sp,fall_speed=fall_speed,pollen_radius=pollen_radius
		    ,pollen_productivity=pollen_productivity,gamma1=gamma1
		    ,diffusion_coefficient=diffusion_coefficient,gravity=gravity,rho0=rho0,nos_lake=nos_lake
                    ,rho=rho,mu=mu,pollen_deposition=pollen_deposition,lake_centre_i=lake_centre_i,lake_centre_j=lake_centre_j
		    ,mean_sp_abundance_at_distance=mean_sp_abundance_at_distance,nos_water=nos_water,water_locations=water_locations
                    ,outfile=outfile,intermediate_file=intermediate_file,outlength=outlength,recal_run=recal_run)

     # pass functins needed
     function_info=c("generate_species_map","mask_lake_area","lake_distance_map","prentice_sugita","read_map","pollen_deposition_model")

  # inform the user
  print("Beginning iterative proceedures...")

  if (use_parallel & heavy_lifting != "f90") {
      # use parallel functions
      cl <- makeCluster(nos_cores, type = "PSOCK")
      clusterExport(cl,function_info) 
      sum_pol=parLapply(cl,1:nos_iter,fun=iterative_loop,model_info=model_info,i_extent=i_extent,j_extent=j_extent)
      stopCluster(cl) ; closeAllConnections()
  } else if (heavy_lifting == "f90") {
      old_wd=getwd() ; setwd(paste(old_wd,"/f90_code/",sep=""))
      if (recompile_f90) {
          # for the moment lets assume we are hard-coding the compilation
          system("rm pollen.so") ; system("rm pollen_deposition_model_mod.mod") ; system("rm R_to_f90_interface.o")
          compiler="gfortran"
          system(paste(compiler," -shared functions_pollen_deposition_model.f90 R_to_f90_interface.f90 -o pollen.so -fPIC",sep=""))
      }
      # determine the number of iterations per node
      if (use_parallel & nos_cores > nos_iter) {
          stop('number of cores greater than number of iterations requested please check the control script to ensure that this is not the case')
      }
      # assume serial behaviour is default
      nositer_group=nos_iter
      # call function which submits the f90 code based on number of iterations to cores
      if (use_parallel & nos_cores > 1) {
          # adjust the iteration groupings for each core
          nositer_group=floor(nos_iter/nos_cores)
          nositer_group=rep(nositer_group, times=nos_cores)
          nositer_group[nos_cores]=nositer_group[nos_cores]+nos_iter%%nos_cores
          # create virtual cluster and submit job
          cl <- makeCluster(nos_cores, type = "PSOCK")
          sum_pol=parLapply(cl,1:nos_cores,fun=call_f90_code,model_info=model_info,i_extent=i_extent,j_extent=j_extent,nositer_group=nositer_group,wdir=getwd())
          stopCluster(cl) ; closeAllConnections()
#          cl <- makeForkCluster(nos_cores) # fork after the variables have been set up
#          sum_pol=parLapply(cl,1:nos_cores,fun=call_f90_code,model_info=model_info,i_extent=i_extent,j_extent=j_extent,nositer_group=nositer_group,wdir=getwd())
#          stopCluster(cl) ; closeAllConnections()
      } else {
          # serial f90 useage
          sum_pol=call_f90_code(1,model_info,i_extent,j_extent,nositer_group,getwd())
      } # parallel path or not
      # don't forget to go back to where we began
      if (recompile_f90) { setwd(old_wd) }
  } else {
      # serial approach R
      sum_pol=lapply(1:nos_iter,FUN=iterative_loop,model_info=model_info,i_extent=i_extent,j_extent=j_extent)
  }

  # some restucturing required
  sum_pol=t(array(unlist(sum_pol),dim=c(nos_species,nos_iter)))

  # now add the restart values back in there if they exist
  if (restart) { sum_pol=rbind(restart_sum_pol,sum_pol) }

  if (ecotone == TRUE & lake_fixed == FALSE) {

      # begin writing out report to the user
      report=data.frame( Pollen_load_species_1=(sum_pol[,1]/rowSums(sum_pol))*100, Pollen_load_species_2=(sum_pol[,2]/rowSums(sum_pol))*100 )

  } else if (random_or_read == "read" & nos_iter > 1){	
      
      # recreate the pollen productivies applied
      forest_productivity=rep(1,times=nos_iter)
      savanna_productivity=rep(c(0.0675,0.125,0.25,0.5,1,2,4,8,16),length.out=nos_iter)
      
      # begin writing out report to the user
      report=data.frame( Pollen_productivity_species_1 = forest_productivity, Pollen_productivity_species_2 = savanna_productivity
			,Pollen_load_species_1=(sum_pol[,1]/rowSums(sum_pol))*100, Pollen_load_species_2=(sum_pol[,2]/rowSums(sum_pol))*100 
			,Pollen_raw_species_1=sum_pol[,1], Pollen_raw_species_2=sum_pol[,2])

 } else {
      # begin writing out report to the user
      for (sp in seq(1,nos_species)) {
	    if (sp == 1) {report=data.frame( as.numeric(as.character(round(quantile(sum_pol[,sp]/rowSums(sum_pol), prob=c(0.5))*100,digits=2))),
				      as.numeric(as.character(round(quantile(sum_pol[,sp]/rowSums(sum_pol), prob=c(0.975))*100,digits=2))),
				      as.numeric(as.character(round(quantile(sum_pol[,sp]/rowSums(sum_pol), prob=c(0.025))*100,digits=2))),
				      as.numeric(as.character(round(sd(sum_pol[,sp]/rowSums(sum_pol))*100,digits=2))),
				      as.numeric(as.character(round(min(sum_pol[,sp]/rowSums(sum_pol))*100,digits=2))), 
				      as.numeric(as.character(round(max(sum_pol[,sp]/rowSums(sum_pol))*100,digits=2))))
	    } else {
	    report=rbind(report,c(as.numeric(as.character(round(quantile(sum_pol[,sp]/rowSums(sum_pol), prob=c(0.5))*100,digits=2))),
				  as.numeric(as.character(round(quantile(sum_pol[,sp]/rowSums(sum_pol), prob=c(0.975))*100,digits=2))),
				  as.numeric(as.character(round(quantile(sum_pol[,sp]/rowSums(sum_pol), prob=c(0.025))*100,digits=2))),
				  as.numeric(as.character(round(sd(sum_pol[,sp]/rowSums(sum_pol))*100,digits=2))),
				  as.numeric(as.character(round(min(sum_pol[,sp]/rowSums(sum_pol))*100,digits=2))), 
				  as.numeric(as.character(round(max(sum_pol[,sp]/rowSums(sum_pol))*100,digits=2)))))
	    } # end sp = 1
      } # for species
      colnames(report)<-c("Median_est_%","Upper_95CI_%","Lower_95CI_%","SD_estimate","Minium_est_%","Maximum_est_%")
  }

  # garbage collection
  rm(sum_pol) ; gc(reset=TRUE, verbose=FALSE)  
 
  # write output to target file
  write.table(report, file=outfile, sep=',', col.names=TRUE, row.names=FALSE)
  # and now we can delete the intermediate file
  if (file.exists(intermediate_file)) {tmp = file.remove(intermediate_file)}

  # keep time
  print(paste("Time for model run (minutes) = ",round((proc.time()["elapsed"]-stime)/60, digits=1),sep=""))

  # return to the user
  return(report)


} # end of run_mode 

