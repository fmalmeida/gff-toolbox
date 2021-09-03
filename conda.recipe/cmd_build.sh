# remove everything
conda build purge-all
conda clean -afy

# build package
conda-build --user falmeida . -c conda-forge -c bioconda -c anaconda -c defaults -c falmeida

# convert to other SOs
conda convert -p osx-64 $CONDA_PREFIX/conda-bld/linux-64/gff-toolbox-*.tar.bz2
conda convert -p win-64 $CONDA_PREFIX/conda-bld/linux-64/gff-toolbox-*.tar.bz2

# upload
anaconda upload osx-64/* --force
anaconda upload win-64/* --force

# clean dir
rm -rf osx-64 win-64
