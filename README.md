# data_viz_shiny_app

This repository contains a Shiny app hosted at https://jfneum.shinyapps.io/data_viz_shiny_app/

The app requires the following inputs:
1. A csv file containing at least the columns "Symbol" (HGNC genes), "adj.p.value" and "logFC"
2. A DisGeNet account (email & password) which can be created here: https://www.disgenet.org/signup/

The app should be possible to run using Docker but one package that needs to be installed from bitbucket and uses an outdated dependency currently crashes the container.
