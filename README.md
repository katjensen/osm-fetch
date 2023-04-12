# osm-fetch
___
A repo dedicated to fetching geometries from OpenStreetMap using the Overpass QL API

## Installation
___
### Clone the repo
```console
git clone git@github.com:katjensen/osm-fetch.git
```

### Getting set up
This project uses `poetry` for python package management. If you're new to `poetry`, check out their ["Basic Usage" page](https://python-poetry.org/docs/basic-usage/).

Install Poetry. (Make sure to restart your terminal after installation)
```console
curl -sSL https://install.python-poetry.org | python3 -
```

Create and activate a fresh, new python environment using conda (or pyenv), with Python >=3.8
```console
conda create -n osm-fetch python=3.9
conda activate osm-fetch
```

From the project root directory, you can install this environment, along with pre-commit git hooks, with this command:
```console
make install
```

If you want, you may install `osm-fetch` in your current environment in editable mode
```console
pip install -e .
```

### Some additional tips

To add a package to your project using poetry CLI,

```console
poetry add "[your_package]"
```

You can also manually update `pyproject.toml` to include a new dependency or change dependency versions. This will require you to recreate the `poetry.lock` file in the project, which you can do with the update command below.
```console
poetry update
```
