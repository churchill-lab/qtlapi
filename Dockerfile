FROM mattjvincent/simplerestrserve:0.4.0
LABEL maintainer="Matthew Vincent <mattjvincent@gmail.com>" \
	  version="0.6.0"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
	    supervisor && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# install the dependencies and qtl2
RUN R -e 'remotes::install_version("missMDA", version = "1.18")' \
 && R -e 'remotes::install_version("dbplyr", version = "2.1.1")' \
 && R -e 'remotes::install_version("pryr", version = "0.1.5")' \
 && R -e 'remotes::install_version("RSQLite", version = "2.2.7")' \
 && R -e 'remotes::install_version("gtools", version = "3.9.2")' \
 && R -e 'remotes::install_version("qtl2", version="0.24")'
 
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
