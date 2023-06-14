# RNA-seq analysis

## Run via docker

1. Build the docker images `docker build -t rstudio .`
2. Run the docker container `docker run --rm -p 8787:8787 -v $(pwd)/src/:/home/rstudio/scripts rstudio`
3. Clean dangling images `docker image prune && docker container prune`
4. Go to `http://localhost:8787/` username: `rstudio` password: `gyara2`