<p align="left">
<img src="https://github.com/fmalmeida/gff-toolbox/raw/master/images/GFF_ToolBox.png" alt="logo" width="500px"/>
</p>

gff-toolbox is a toolbox of commands that tries to facilitate the work with GFF annotation files.

## Table of contents

* [Requirements](https://github.com/fmalmeida/gff-toolbox#requirements)
* [Installation](https://github.com/fmalmeida/gff-toolbox#installation)
* [Documentation](https://github.com/fmalmeida/gff-toolbox#documentation)
* [How can you colaborate?](https://github.com/fmalmeida/gff-toolbox#collaborating)
* [Citation](https://github.com/fmalmeida/gff-toolbox#citation)

## Requirements

This pipeline has only two requirement:

* Python >= 3.6, and its packages:
    * [pandas](https://pandas.pydata.org/)
    * [biopython](https://biopython.org/)
    * [bcbiogff](https://github.com/chapmanb/bcbb/tree/master/gff)
    * [dna features viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer)
    * [matplotlib](https://matplotlib.org/)
    * [docopt](http://docopt.org/)
    * [pymongo](https://pypi.org/project/pymongo/)
* [mongo shell](https://docs.mongodb.com/manual/)
    * [In conda](https://anaconda.org/anaconda/mongodb)

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

This is meant to be a collaborative project, which means it is meant to adapt to the community needs. Thus, we encourage users to use it and to collaborate with ideas for different implementations, new commands, additions, etc.

If you have an analysis that you constantly do when working with GFFs and would like to see it implemented in a command-like package to make your life easy, or whenever you feel something can be improved, don't hesitate and **collaborate**.

You can collaborate by:

* flagging an **enhancement issue** discussing your idea in the homepage of the project
* you can fork the repo, create and start the implementation of your own script/command in the project and then submit a **pull request**
    * I'll then check the request, make sure it is in the same format and standards of the already implemented commands and **confirm** the inclusion.
    * Of course, you will be recognized as the developer/creator of that specific implementation.

Checkout more at about forking and contributing to repos at:

* [Chase's tutorial](https://gist.github.com/Chaser324/ce0505fbed06b947d962)
* [github's advises on how to collaborate to projects](https://docs.github.com/en/free-pro-team@latest/github/collaborating-with-issues-and-pull-requests)

## Citation

To cite this pipeline users can use the github url. Users are encouraged to cite the python packages used in this pipeline whenever their outputs are valuable.
