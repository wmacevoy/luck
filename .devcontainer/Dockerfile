FROM devcontainer-latex:latest

# Reference implementations of the luck statistics live in scilab/, R/,
# python/, and maxima/.  This image adds the runtimes needed to execute
# them; the LaTeX toolchain comes from the base.

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3 \
        python3-pip \
        python3-numpy \
        python3-scipy \
        r-base \
        r-cran-testthat \
        scilab \
        wxmaxima \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
