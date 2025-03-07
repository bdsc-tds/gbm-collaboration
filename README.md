<div align="center">
<h1>GBM Collab repository</h1>
Data exploration and analyses of the GBM Mosaic data
</div>

# Table of Contents
- [Installation](#installation)


# Installation
1. Start by initialising Conda.
```bash
conda init bash && bash
```

2. Create the environment.
```bash
conda create -n gbmcollaboration -y python=3.10
conda activate gbmcollaboration
```

3. Clone the repository and install the dependencies.
```bash
git clone https://github.com/owkin/gbm-collaboration.git
cd gbm-collaboration
make install-poetry
make install-all
pre-commit install
```
