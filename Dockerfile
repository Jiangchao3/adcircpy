FROM frolvlad/alpine-miniconda3
RUN apk add tini
RUN apk add bash
ADD conda-linux-x86_64.yml /tmp/conda-linux-x86_64.yml
RUN conda update -n base -c defaults conda
RUN conda env create -f /tmp/conda-linux-x86_64.yml
RUN echo "source activate $(head -1 /tmp/conda-linux-x86_64.yml | cut -d' ' -f2)" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /tmp/conda-linux-x86_64.yml | cut -d' ' -f2)/bin:$PATH
ENTRYPOINT ["/sbin/tini", "--"]
CMD ["/bin/bash"]
