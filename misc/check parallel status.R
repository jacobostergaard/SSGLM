library(misc)
library(SSGLM)
library(parallel)
clean_up()
datalib = "/Users/jacob/Data/" #"Library/Mobile Documents/com~apple~CloudDocs/GitHub/Source/R/packages/SSGLM/data/"
set.seed(1234)


# format(Sys.time(), "%a %b %d %H:%M:%S %Y")
tmp  = Sys.time()
wait = 3600-as.numeric(format(tmp,"%M"))*60+as.numeric(format(tmp,"%S")) # wait untill on the hour
Sys.sleep(wait)

while(as.numeric(format(tmp, "%H")) < 22){
  tmp = Sys.time()
  # if(format(tmp, "%M:%S") == "00:00"){ #on the hour
  #
  # }

  # List current files in datalib
  fns = list.files(datalib)
  fns = fns[!grepl("turtle2",fns)]
  fns = fns[grepl("x2",fns)]

  tmp.info = file.info(paste0(datalib,fns))$ctime       # update times for files
  # fns[order(tmp.info)]
  # tmp.info[grepl("n1_",fns)] <= tmp.info
  idx = as.numeric(tmp)-15*3600 <= as.numeric(tmp.info) # indices of files modified within last 15 hours
  fns = fns[idx]

  pct = round(100*length(fns)/249,2)
  msg = paste0(pct,"% neurons analyzed at ",tmp)

  notify(msg)
  cat("\n",msg)
  Sys.sleep(3600)
}


# fn = paste0(datalib,fns[2])
# load(fn)
# tmpres[[1]]
