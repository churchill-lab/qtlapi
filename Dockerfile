FROM mattjvincent/simpleplumber:r.3.4.3_p.0.4.4
LABEL maintainer="Matthew Vincent <mattjvincent@gmail.com>" \
	  version="1.0"

# install the dependencies
RUN R -e 'devtools::install_version("data.table", version = "1.10.4-3")' \
 && R -e 'devtools::install_version("dbplyr", version = "1.1.0")' \
 && R -e 'devtools::install_version("jsonlite", version = "1.5")' \
 && R -e 'devtools::install_version("pryr", version = "0.1.3")' \
 && R -e 'devtools::install_version("RcppEigen", version = "0.3.3.3.1")' \
 && R -e 'devtools::install_version("RSQLite", version = "2.0")' \
 && R -e 'devtools::install_version("yaml", version = "2.1.15")' \
 && R -e 'devtools::install_version("gtools", version = "3.5.0")'

# install qtl2 tools
RUN R -e 'devtools::install_github("rqtl/qtl2", ref = "0.12")'
RUN R -e 'devtools::install_github("churchill-lab/intermediate", ref = "v.1.0")'

EXPOSE 8000

ENV INSTALL_PATH /app/qtlapi
RUN mkdir -p $INSTALL_PATH/data

WORKDIR $INSTALL_PATH

COPY qtlapi.R qtlapi.R
COPY VERSION .

CMD ["/app/qtlapi/qtlapi.R"]
