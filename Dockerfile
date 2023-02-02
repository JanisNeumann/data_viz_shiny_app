FROM rocker/shiny:4.2.2

# Install additional R packages that are available through CRAN.
RUN R -e 'install.packages(c("ggplot2", "dplyr", "rbioapi"))'

# Copy scripts to Docker image.
COPY data_viz_shiny_app.Rproj /data_viz_shiny_app/
COPY app.R /data_viz_shiny_app/
COPY example_data.csv /data_viz_shiny_app/

# Set the working directory.
WORKDIR /data_viz_shiny_app

# Run the workflow.
CMD ["R", "-e", "shiny::runApp('app.R')"]