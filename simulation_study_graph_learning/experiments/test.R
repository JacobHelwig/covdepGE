# parse command line arguments
args <- (commandArgs(TRUE))
print(args)
if (length(args) > 0){
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
}else{
  stop("No arguments passed")
}
print(ls())
