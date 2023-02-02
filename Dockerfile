FROM rocker/shiny:4.2.2

# Install additional R packages that are available through CRAN.
RUN R -e 'install.packages(c("devtools", "ggplot2", "dplyr", "rbioapi"))'

# Install disgenet2r from repo and, before that, the archived SPARQL library it needs.
RUN R -e 'devtools::install_version("SPARQL", version = "1.16", repos = "https://CRAN.R-project.org")'
RUN R -e 'devtools::install_bitbucket("ibi_group/disgenet2r")'

# Copy scripts to Docker image.
COPY data_viz_shiny_app.Rproj /data_viz_shiny_app/
COPY app.R /data_viz_shiny_app/
COPY example_data.csv /data_viz_shiny_app/

# Set the working directory.
WORKDIR /data_viz_shiny_app

# Run the workflow.
CMD ["R", "-e", "shiny::runApp('app.R')"]