FROM mattjvincent/simplerestrserve:0.2.0
LABEL maintainer="Matthew Vincent <mattjvincent@gmail.com>" \
	  version="0.5.0"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
	    supervisor && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# install the dependencies
RUN R -e 'remotes::install_version("missMDA", version = "1.15")' \
 && R -e 'remotes::install_version("dbplyr", version = "1.4.2")' \
 && R -e 'remotes::install_version("pryr", version = "0.1.4")' \
 && R -e 'remotes::install_version("RSQLite", version = "2.2.0")' \
 && R -e 'remotes::install_version("gtools", version = "3.8.1")'
 
# install qtl2 tools
RUN R -e 'remotes::install_github("rqtl/qtl2@0.20")'

# install intermediate
RUN R -e 'remotes::install_github("churchill-lab/intermediate@v.2.5")'

SHELL ["/bin/bash", "-c"]

ENV INSTALL_PATH /app/qtlapi
RUN mkdir -p $INSTALL_PATH/data $INSTALL_PATH/conf

WORKDIR $INSTALL_PATH

COPY ./src/qtlapi.R qtlapi.R
COPY ./src/restapi.R restapi.R
COPY ./src/run.R run.R
COPY ./conf/supervisor.conf $INSTALL_PATH/conf/supervisor.conf
COPY VERSION .

CMD ["/usr/bin/supervisord", "-n", "-c", "/app/qtlapi/conf/supervisor.conf"]
