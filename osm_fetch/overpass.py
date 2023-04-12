import logging
from typing import Dict, List, Optional, Union

import rasterio as rio
import requests
from shapely.geometry import Polygon, shape
from tenacity import retry, retry_if_not_exception_type, stop_after_attempt, wait_fixed

from osm_fetch.exceptions import (
    OverpassBadRequest,
    OverpassGatewayTimeout,
    OverpassMoved,
    OverpassTooManyRequests,
)

EPSG_4326 = rio.crs.CRS.from_epsg(4326)  # geo coords

logger = logging.getLogger(__name__)


def write_query(
    bbox: Polygon,
    tag: str,
    crs: rio.crs.CRS = EPSG_4326,
    value_constraints: Optional[Union[list, None]] = None,
    timeout: Optional[int] = 1000,
) -> str:
    """
    Translate a bounding box and OSM tag into an Overpass query
        * see more: https://wiki.openstreetmap.org/wiki/Overpass_API
        * test out queries online: https://overpass-turbo.eu/

    Args:
        bbox: Bounding box polygon
        tag: OSM tag to query (ex: "highway")
        crs: CRS of the provided bbox, default= geo coords
        value_contraints: List of restricted values for the provided OSM tag, or None queries ALL possible
        timeout: Overpass timeout, in seconds

    Returns:
        query: Formatted Overpass query

    TODO: enable more flexible format on input bbox geometry?
    """

    if crs != EPSG_4326:  # make sure its geo coords!
        bbox = shape(rio.warp.transform_geom(src_crs=crs, dst_crs=EPSG_4326, geom=bbox))

    # Shapely bounds = (W, S, E, N) -->  but the overpass ql expects bbox format = (S,W,N,E)
    bbox_tuple = (bbox.bounds[1], bbox.bounds[0], bbox.bounds[3], bbox.bounds[2])

    if value_constraints:
        if len(value_constraints) > 1:
            query = f'["{ tag }"~"{ "|".join(value_constraints) }"]'
        else:
            query = f'["{ tag }"="{ value_constraints[0] }"]'
    else:
        query = f'["{ tag }"]'

    return f"[out:json][timeout:{ timeout }]; nwr{ query }{ bbox_tuple }; out geom qt;"


def dict_to_query(
    bbox: Polygon,
    tag_values_dict: Dict[str, Union[None, List]],
    crs: rio.crs.CRS = EPSG_4326,
    timeout: Optional[int] = 1000,
) -> str:
    """
    Translate a bounding box and a whole dictionary of OSM tags into an Overpass query
        * see more: https://wiki.openstreetmap.org/wiki/Overpass_API
        * test out queries online: https://overpass-turbo.eu/

    Args:
        bbox: Bounding box polygon
        tag: OSM tag to query (ex: "highway")
        crs: CRS of the provided bbox, default= geo coords
        value_contraints: List of restricted values for the provided OSM tag, or None queries ALL possible
        timeout: Overpass timeout, in seconds

    Returns:
        query: Formatted Overpass query

    TODO: enable more flexible format on input bbox geometry?
    """

    if crs != EPSG_4326:  # make sure its geo coords!
        bbox = shape(rio.warp.transform_geom(src_crs=crs, dst_crs=EPSG_4326, geom=bbox))

    # Shapely bounds = (W, S, E, N) -->  but the overpass ql expects bbox format = (S,W,N,E)
    bbox_tuple = (bbox.bounds[1], bbox.bounds[0], bbox.bounds[3], bbox.bounds[2])

    full_query = f"[out:json][timeout:{ timeout }];("

    for tag in tag_values_dict:
        value_constraints = tag_values_dict[tag]
        if value_constraints:
            if len(value_constraints) > 1:
                sub_query = f'["{ tag }"~"{ "|".join(value_constraints) }"]'
            else:
                sub_query = f'["{ tag }"="{ value_constraints[0] }"]'
        else:
            sub_query = f'["{ tag }"]'

        full_query += f"nwr{ sub_query }{ bbox_tuple };"

    full_query += ");out geom qt;"
    return full_query


def request(query: str, endpoint: str) -> dict:
    """Send a request to the Overpass API.

    Args:
        query : Overpass QL query
        endpoint: API endpoint

    Returns
        response : JSON response as a dictionary

    TODO: point to our own local instance ? relying on external server is not ideal!
    """
    response = requests.get(endpoint, params={"data": query})

    if response.status_code == 302:
        raise OverpassMoved
    elif response.status_code == 400:
        raise OverpassBadRequest
    elif response.status_code == 429:
        raise OverpassTooManyRequests
    elif response.status_code == 504:
        raise OverpassGatewayTimeout

    return response.json()


def retry_info(retry_state):
    # consider updating the decorator below (needs more testing):
    #   before_sleep = before_sleep_log(logger, logging.INFO),
    #   after = after_log(logger, logging.INFO),
    logger.info(
        f"OSM request attempt #{retry_state.attempt_number} ended with: {retry_state.outcome}"
    )


@retry(
    retry=retry_if_not_exception_type((OverpassMoved, OverpassBadRequest)),
    stop=stop_after_attempt(10),
    wait=wait_fixed(10),
    before_sleep=retry_info,
)
def request_osm(
    query: str, endpoint: str = "http://overpass-api.de/api/interpreter"
) -> Dict:
    """
    Wrapper for `request()` that attempts retries

    Args:
        query : Overpass QL query
        endpoint: API endpoint, default = "http://overpass-api.de/api/interpreter"

    Returns
        response : JSON response as a dictionary

    TODO: point to our own local instance ? relying on external server is not ideal!
    """
    return request(query, endpoint)
