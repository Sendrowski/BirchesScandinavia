
#binpath <- paste(Sys.getenv("CONDA_PREFIX"), 'bin', sep="/")

install.packages(c("shiny", "shinycssloaders", "shinythemes", "shinydashboard",
    "shinydashboardPlus", "shinyjs", "DT", "shinyWidgets", "shinyhelper", "plotly", "viridis",
    "tidyr", "RColorBrewer", "yaml", "ggpubr", "FactoMineR", "data.table", "ggpubr", "nnet",
    "tidyverse", "abcrf", "matrixStats", "ranger", "RcppArmadillo"),
    dep=T)

library("devtools");
with_libpaths(new=binpath, install_github("nik01010/dashboardthemes"));

#lib=binpath,
#, repos="http://cran.r-project.org"
