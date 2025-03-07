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

4. Optional, fork the repository
If you aim to share your code with your team members, consider doing
a fork or create a copy of the actual repository with a new name.

5. Clone the repository and install the dependencies.
```bash
git clone https://github.com/owkin/gbm-collaboration.git
cd gbm-collaboration
make install-poetry
make install-all
pre-commit install
```

Note that this also makes your environment visible in Jupyter.
