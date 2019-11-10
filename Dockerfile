FROM mattjvincent/simplerestrserve:0.0.2
LABEL maintainer="Matthew Vincent <mattjvincent@gmail.com>" \
	  version="0.4.4"

# install the dependencies
RUN R -e 'devtools::install_version("missMDA", version = "1.14")' \
 && R -e 'devtools::install_version("dbplyr", version = "1.4.2")' \
 && R -e 'devtools::install_version("pryr", version = "0.1.4")' \
 && R -e 'devtools::install_version("RSQLite", version = "2.1.1")' \
 && R -e 'devtools::install_version("gtools", version = "3.8.1")'
 
# install qtl2 tools
RUN R -e 'devtools::install_github("rqtl/qtl2", ref = "0.16")'

# install intermediate
RUN R -e 'devtools::install_github("churchill-lab/intermediate", ref = "v.2.5")'

SHELL ["/bin/bash", "-c"]

ENV INSTALL_PATH /app/qtlapi
RUN mkdir -p $INSTALL_PATH/data $INSTALL_PATH/conf

WORKDIR $INSTALL_PATH

COPY ./src/qtlapi.R qtlapi.R
COPY ./src/restapi.R restapi.R
COPY ./src/rawapi.R rawapi.R
COPY ./conf/Rserve.conf $INSTALL_PATH/conf/Rserve.conf
COPY ./conf/supervisor.conf $INSTALL_PATH/conf/supervisor.conf
COPY VERSION .

CMD ["/usr/bin/supervisord", "-n", "-c", "/app/qtlapi/conf/supervisor.conf"]
