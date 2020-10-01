<p align="left">
<img src="https://github.com/fmalmeida/gff-toolbox/raw/master/images/GFF_ToolBox.png" alt="logo" width="500px"/>
</p>

gff-toolbox is a toolbox of commands that tries to facilitate the work with GFF annotation files.

## Table of contents

* [Requirements](https://github.com/fmalmeida/gff-toolbox#requirements)
* [Installation](https://github.com/fmalmeida/gff-toolbox#installation)
* [Documentation](https://github.com/fmalmeida/gff-toolbox#documentation)
* [How you can colaborate?](https://github.com/fmalmeida/gff-toolbox#collaborating)
* [Citation](https://github.com/fmalmeida/gff-toolbox#citation)

## Requirements

This pipeline has only one requirement:

* Python >= 3.6, and its packages:
    * [pandas](https://pandas.pydata.org/)
    * [biopython](https://biopython.org/)
    * [bcbiogff](https://github.com/chapmanb/bcbb/tree/master/gff)
    * [dna features viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer)
    * [matplotlib](https://matplotlib.org/)
    * [docopt](http://docopt.org/)

## Installation

Installation is super easy and perhaps not required:

```bash
# Download
git clone https://github.com/fmalmeida/gff-toolbox.git
cd gff-toolbox

# Run without installing
python3 gfftoolbox-runner.py -h

# Install and run in any place
python3 setup.py install
gff-toolbox -h
```

## Documentation

The command is very well documented inside its help messages and examples which can be checked as:

* `gff-toolbox -h`

Additionally, users can checkout the execution and results of commands in the [wiki](https://github.com/fmalmeida/gff-toolbox/wiki) provided.

## Collaborating

## Citation

To cite this pipeline users can use the github url. Users are encouraged to cite the python packages used in this pipeline whenever their outputs are valuable.
