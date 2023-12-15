#' Code to manage and run Atlantis simulations using parallel processing and doAzure
#'  
#' @author Alberto Rovellini
#' @date October 2023
#' # https://github.com/Azure/doAzureParallel
#' 

#########################################################
# AR: We need to build Atlantis within the loop. Not ideal, but I cannot get the custom / private docker image to pull and there is no support online
# In the scheme of a 50+ hrs Atlantis run building is not the end of the world 


# set locale to avoid multibyte errors
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
# https://www.r-bloggers.com/web-scraping-and-invalid-multibyte-string/

# devtools::install_github("omegahat/XMLSchema")
# remotes::install_github("sckott/SSOAP") # https://rdrr.io/github/sckott/SSOAP/
# install the Azure parallel packages
# devtools::install_github(c("Azure/rAzureBatch", "Azure/doAzureParallel"))
# install my version of the doAzureParallel package
devtools::install_github("Azure/rAzureBatch")
devtools::install_github("somros/doAzureParallel")
                         
# List of packages for session
# these are the packages that the Controller needs
.packages = c("devtools", "dtplyr","stringi","data.table","tidyverse","stringr","R.utils",
              "magrittr","future","parallel","doSNOW","XMLSchema", "SSOAP",
              "doAzureParallel", "rAzureBatch","here")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

#set paths
# not sure what these should be yet, let's keep going
# workpath <- "~/Dropbox/GOM_Diets/"
# shpath <- "~/Dropbox/GOM_Diets/prms" # this does not appear elsewhere
# savepath <- "~/Atlantis"

#setwd("~/")

# Upgrade files
system("sudo apt-get update; sudo apt-get upgrade -y", wait = TRUE)

# create the prm and run files if needed
# source(here('NOAA_Azure','code','create_prm_sh_files.R')) 
# This needs to be pushed to github

# Fill out your credential config and cluster config files.
# Enter your Azure Batch Account & Azure Storage keys/account-info into your credential config ("credentials.json") and configure your cluster in your cluster config ("cluster.json")
# See https://github.com/Azure/doAzureParallel for details
# Set your credentials - you need to give the R session your credentials to interact with Azure
setCredentials("NOAA_Azure/keys/credentials.json")

# Register the pool. This will create a new pool if your pool hasn't already been provisioned.
cluster <- doAzureParallel::makeCluster("NOAA_Azure/cluster.json")#, resourceFiles = resource_files)
# check the status of your cluster (it seems to hang on the provisioning in the R console sometimes)
getClusterList(filter = NULL)

# Register the pool as your parallel backend
registerDoAzureParallel(cluster)

# Check that your parallel backend has been registered
getDoParWorkers()

# read in lookup table for the indices
f_lookup <- read.csv(here('NOAA_Azure','data','f_lookup_4.csv'), header = TRUE)
runidx <- f_lookup$idx # as many runs as you want to do, would be 96 for us
#runidx <- 1:4 #for debugging

if("package:doSNOW" %in% search()) detach("package:doSNOW", unload=TRUE) 
if("package:future" %in% search()) detach("package:future", unload=TRUE) 
if("package:parallely" %in% search()) detach("package:parallely", unload=TRUE) 
if("package:snow" %in% search()) detach("package:snow", unload=TRUE) 
if("package:parallel" %in% search()) detach("package:parallel", unload=TRUE) 

# azure options
opts <- list(maxTaskRetryCount = 0, # do not retry a task if it fails
             maxDate = Sys.time() + 60 * 60 * 24 * 10, # 10 days
             autoDeleteJob=FALSE)

# Run the parallel loop
atlantis.scenarios <- foreach(idx=runidx, .options.azure=opts, .errorhandling = 'pass', .verbose = TRUE) %dopar% { # move to dopar to run jobs in parallel, do does it on local host
  
  # # Upgrade files
  system("sudo apt-get update; sudo apt-get upgrade -y", wait = TRUE)
  system("sudo apt-get install -y apt-utils", wait = TRUE)
  
  # # Install dependencies needed for Atlantis
  # # Hem's sh file
  # # system("sudo apt-get -y install autoconf automake build-essential libcurl4 libcurl4-openssl-dev curl gdal-bin libgdal-dev libgeos++-dev libgeos-dev flip libcairo2 libcairo2-dev libapparmor1 libhdf5-dev libnetcdf-dev libxml2-dev libproj-dev libssl-dev libv4l-0 libgeotiff5 libglu1-mesa-dev libpoppler-cpp-dev libprotobuf-dev librsvg2-dev libx11-dev lsscsi openjdk-8-jdk python2.7 python3-pip python3-dev proj-bin proj-data protobuf-compiler openssl rpm mesa-common-dev netcdf-bin ntp ntpdate subversion valgrind cdo nco m4 dos2unix gawk gfortran libfribidi-dev libharfbuzz-dev libudunits2-dev libv8-dev gnupg-agent software-properties-common dirmngr --no-install-recommends", wait = TRUE)
  # Hem's original
  system("sudo apt-get install -y subversion build-essential flip autoconf libnetcdf-dev libxml2-dev gawk", wait = TRUE)
   
  #use this to install Proj4 which is required for Atlantis
  #https://gist.github.com/robinkraft/2a8ee4dd7e9ee9126030
  #http://grasswiki.osgeo.org/wiki/Compile_and_Install_Ubuntu#PROJ4
  #uninstall   proj-bin libproj-dev if present
  #look for proj with locate libproj.so
  # sudo apt install mlocate
  system("sudo apt-get install -y subversion make")
  system("cd /tmp; svn co http://svn.osgeo.org/metacrs/proj/branches/4.8/proj/")
  system("cd /tmp/proj/nad; sudo wget http://download.osgeo.org/proj/proj-datumgrid-1.5.zip; unzip -o -q proj-datumgrid-1.5.zip")
  system("cd /tmp/proj/; ./configure  &&  make  &&  sudo make install && sudo ldconfig")

  # Copy Atlantis code
  # I cannot seem to be able to load this from the docker image (cannot get the custom docker repo to pull)
  # Cannot go the Azure blob storage route also for credential issues possibly (see https://github.com/Azure/doAzureParallel/issues/82)
  # The only option left if to svn checkout - this needs credentials so do not commit this to github as it gives access to the code
  # Get code from CSIRO SVN, if need to change version use co -6665 for example
  # do not commit these lines
  # svn_username <- ""
  # svn_pw <- ""
  # 
  # system(paste("svn co -r6665 https://svnserv.csiro.au/svn/ext/atlantis/Atlantis/trunk --username", svn_username, "--password", svn_pw, "--quiet"), wait = TRUE) # this works it seems

  ### v6665 comes with loads of fprints that need to be silenced or will produce tons of output
  # Alternative is to clone a private GitHub repo instead with the code and then build from that
  
  pat <- ""
  system(paste0("git clone https://",pat,"@github.com/somros/for_pw.git"))
  
  # to compile
  system("cd for_pw/trunk/atlantis; aclocal; autoheader; autoconf; automake -a; ./configure; sudo make CFLAGS='-DACCEPT_USE_OF_DEPRECATED_PROJ_API_H -Wno-misleading-indentation -Wno-format -Wno-implicit-fallthrough'; sudo make -d install", wait = TRUE)
  
  # now clone the F test folder from GitHub
  system("git clone https://github.com/somros/AtlantisGOA_F_test_4.git")
  
  # # List of R packages for session
  .packages = c("data.table","dplyr","tidyr","stringr","R.utils","magrittr","here")

  # Install CRAN packages (if not already installed)
  .inst <- .packages %in% installed.packages()
  if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst], repos="https://cloud.r-project.org/")

  # Load packages into session
  lapply(.packages, require, character.only=TRUE)
  
  # now we have all files where we need them. Point to the correct run file
  this_sh <- paste('RunAtlantis_',idx,".sh",sep="") # this is the sh for this run
  #
  # run Atlantis scenario
  system(paste("cd AtlantisGOA_F_test_4/; flip -uv *; sh ", this_sh, sep=""), wait = TRUE)

  # get the files we need
  outFolder <- paste("AtlantisGOA_F_test_4/output",sep="_")
  biomage_filename <- list.files(path = outFolder, pattern = 'AgeBiomIndx.txt')[1]
  
  biomage_path <- file.path(outFolder, biomage_filename, fsep = "/") # biomass by age
  catch_path <- file.path(outFolder, paste0('output_',idx,'Catch.txt'), fsep = '/')
  
  # read files
  biomage <- read.table(biomage_path, sep = ' ', header = T)
  catch <- read.table(catch_path, sep = ' ', header = T)
  
  # package results 
  list_names <- c("sh", paste("biomage",idx,sep="_"), paste("catch",idx,sep="_"))
  results_list <- list(this_sh, biomage,  catch)
  names(results_list) <- list_names
  return(results_list)
  
}

# Stop cluster
doAzureParallel::stopCluster(cluster)

# Stop controller
# system("az vm deallocate --name AtlantisGOAController --no-wait --resource-group Atlantis_GOA")
