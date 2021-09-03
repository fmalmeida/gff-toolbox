# remove everything
conda build purge-all
conda clean -afy

# build package
conda-build --user falmeida . -c conda-forge -c bioconda -c anaconda -c defaults -c falmeida