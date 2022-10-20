FROM snakemake/snakemake:stable

COPY . /src
RUN pip install snakebids
ENTRYPOINT [ "/src/diffparc/run.py" ]
