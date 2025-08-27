FROM debian:stable

RUN echo base run && ( \
    apt-get update && \
    apt-get install -y --no-install-recommends \
      ca-certificates \
      git \
      make \
      perl \
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* )

RUN echo texlive run && ( apt-get update && \
    apt-get install -y --no-install-recommends \
      texlive \
      asymptote \
      latexmk \
      xindy \
      ghostscript \
    && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* )

RUN echo tufte run && (\
    apt-get update && \
    apt-get install -y --no-install-recommends \
        texlive-fonts-extra \
        texlive-fonts-recommended \
        texlive-latex-base \
        texlive-latex-extra \
        texlive-latex-recommended \
        texlive-plain-generic \
        texlive-science \
        && \
        apt-get clean && \
        rm -rf /var/lib/apt/lists/* )

WORKDIR /data
ENV HOME=/home

VOLUME ["/data","/home"]
