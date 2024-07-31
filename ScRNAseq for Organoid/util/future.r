setup_future <- function(object, plan=c('multicore', 'multisession'), workers=4) {
  plan = match.arg(plan)
  library(future)
  options(future.globals.maxSize = object.size(object)*10)
  options(future.resolve.recursive = Inf)
  options(future.wait.timeout = 10)
  plan(strategy=plan, workers=workers)
  rlog::log_info(paste('util/future.r setup_future: Running', plan, 'future with', workers, 'workers.'))

  # I run into the following error with this:
  # [51] "Error in if (any(value < 0L | value >= max_cores)) { : " 
  # [52] "  missing value where TRUE/FALSE needed"

  #if (nbrOfFreeWorkers() < 1) {
  #  plan(strategy='sequential') 
  #  rlog::log_info('setup_future: Running sequential future.')
  #} else {
  #  plan(strategy=plan, workers=workers)
  #  rlog::log_info(paste('setup_future: Running', plan, 'future with', workers, 'workers.'))
  #}
}
