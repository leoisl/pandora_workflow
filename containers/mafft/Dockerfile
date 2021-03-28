# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/pandora_workflow_mafft

# This container has mafft and some other dependencies to be used in the pandora container
# to build: sudo docker build . -t leandroishilima/pandora_workflow_mafft

FROM continuumio/miniconda3:4.7.12

RUN conda install -c bioconda bioconda::mafft=7.471 click=7.1.2
