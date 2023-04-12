import importlib.metadata

from osm_fetch.extract import json_to_geometries

# import package version from root 'pyproject.toml' file
__version__ = importlib.metadata.version("osm-fetch")
