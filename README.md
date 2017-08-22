# PaCaOmicsAPP
Application for the analysis and visualisation of PaCaOmics data.

## Launching the application
This application requires to have R installed along with two packages (beeswarm and shiny). Installing instruction are shown at the end of this page if needed.

To launch the application, open R and type in these two lines:
```R
library(shiny)
runGitHub("PaCaOmicsAPP", "RemyNicolle")
```

This will download the application along with the necessary data and launch it.


### Installing R and necessary packages
To install R for windows: [https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)
To install R for Mac, download and run the latest _pkg_ :[https://cran.r-project.org/bin/macosx/](https://cran.r-project.org/bin/macosx/)

Once installed, launch R and type the following commands in R to install the necessary packages to launch the application.
```R
install.packages(c("shiny","beeswarm"))
```

