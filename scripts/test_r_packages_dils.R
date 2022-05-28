
binpath <- paste(Sys.getenv("CONDA_PREFIX"), 'bin', sep="/")
Sys.setenv("R_LIBS_USER"=binpath)

print(Sys.getenv("R_LIBS_USER"))

library("abcrf")

