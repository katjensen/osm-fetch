import logging
import warnings
from typing import Dict, List, Optional, Union

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, MultiPolygon, Point, Polygon
from shapely.geos import TopologicalError
from shapely.ops import linemerge, polygonize

from .lookups import POLYGON_FEATURES

# These are the main functions used for extracting geometries from
#   OpenStreetMap from using the Overpass QL API
#
# Many functions have been adapated from
#    - osmnx: https://github.com/gboeing/osmnx
#    - osmextract: https://github.com/ropensci/osmextract

logger = logging.getLogger(__name__)


def _parse_node_to_coords(element: Dict) -> Dict:
    """
    Parse coordinates from a node in the overpass response.
    The coords are only used to create LineStrings and Polygons.

    Args:
        element: element type "node" from overpass response JSON

    Returns
        coords : dict of latitude/longitude coordinates
    """
    # return the coordinate of a single node element
    coords = {"lat": element["lat"], "lon": element["lon"]}
    return coords


def _parse_node_to_point(element: Dict) -> Dict:
    """
    Parse point from a tagged node in the overpass response.

    Args:
        element : element type "node" from overpass response JSON

    Returns
        point : dict of OSM ID, OSM element type, tags and geometry
    """
    point = dict()
    point["osmid"] = element["id"]
    point["element_type"] = "node"

    if "tags" in element:
        for tag in element["tags"]:
            point[tag] = element["tags"][tag]

    point["geometry"] = Point(element["lon"], element["lat"])
    return point


def _parse_way_to_linestring_or_polygon(
    element: Dict, coords: Dict, polygon_features: Optional[Dict] = POLYGON_FEATURES
) -> Dict:
    """
    Parse open LineString, closed LineString or Polygon from OSM 'way'.
        see: https://wiki.openstreetmap.org/wiki/Overpass_turbo/Polygon_Features
        for more information on which tags should be parsed to polygons

    Args:
        element : element type "way" from overpass response JSON
        coords : dict of node IDs and their latitude/longitude coordinates
        polygon_features : dict for determining whether closed ways are LineStrings or Polygons

    Returns:
        linestring_or_polygon : dict of OSM ID, OSM element type, nodes, tags and geometry
    """
    nodes = element["nodes"]
    geom = element["geometry"]

    linestring_or_polygon = dict()
    linestring_or_polygon["osmid"] = element["id"]
    linestring_or_polygon["element_type"] = "way"
    linestring_or_polygon["nodes"] = nodes

    # un-nest individual tags
    if "tags" in element:
        for tag in element["tags"]:
            linestring_or_polygon[tag] = element["tags"][tag]

    # if the OSM element is an open way (i.e. first and last nodes are not the
    # same) the geometry should be a Shapely LineString
    if element["nodes"][0] != element["nodes"][-1]:
        try:
            geometry = LineString(
                [(coords[node]["lon"], coords[node]["lat"]) for node in nodes]
            )

        except KeyError:  # pragma: no cover
            geometry = LineString([(p["lon"], p["lat"]) for p in geom])

    # if the OSM element is a closed way (i.e. first and last nodes are the
    # same) depending upon the tags the geometry could be a Shapely LineString
    # or Polygon
    elif element["nodes"][0] == element["nodes"][-1]:

        # determine if closed way represents LineString or Polygon
        if _is_closed_way_a_polygon(element):
            # if it is a Polygon
            try:
                geometry = Polygon(
                    [(coords[node]["lon"], coords[node]["lat"]) for node in nodes]
                )
            except KeyError:
                geometry = Polygon([(p["lon"], p["lat"]) for p in geom])
            except ValueError as e:
                # XMLs may include geometries that are incomplete, in which
                # case return an empty geometry
                logger.debug(
                    f"{e} .** The geometry for "
                    f"https://www.openstreetmap.org/{element['type']}/{element['id']} was not created."
                )
                geometry = Polygon()

        else:
            # if it is a LineString
            try:
                geometry = LineString(
                    [(coords[node]["lon"], coords[node]["lat"]) for node in nodes]
                )

            except KeyError:
                geometry = LineString([(p["lon"], p["lat"]) for p in geom])

    linestring_or_polygon["geometry"] = geometry
    return linestring_or_polygon


def _is_closed_way_a_polygon(
    element: Dict, polygon_features: Optional[Dict] = POLYGON_FEATURES
) -> bool:
    """
    Determine whether a closed OSM way represents a Polygon, not a LineString.
        Closed OSM ways may represent LineStrings (e.g. a roundabout or hedge
        round a field) or Polygons (e.g. a building footprint or land use area)
        depending on the tags applied to them.

    The starting assumption is that it is not a polygon, however any polygon
        type tagging will return a polygon unless explicitly tagged with area:no.
        It is possible for a single closed OSM way to have both LineString and
        Polygon type tags (e.g. both barrier=fence and landuse=agricultural).

    We return a single Polygon for elements tagged in this way.
        For more information see: https://wiki.openstreetmap.org/wiki/One_feature,_one_OSM_element)

    Args:
        element: closed element type "way" from overpass response JSON
        polygon_features: dict of tag keys with associated values and blocklist/passlist

    Returns
        is_polygon: True if the tags are for a polygon type geometry
    """
    # polygon_features dict is for determining which ways should become Polygons
    # therefore the starting assumption is that the geometry is a LineString
    is_polygon = False

    # get the element's tags
    element_tags = element.get("tags")

    # if the element doesn't have any tags leave it as a Linestring
    if element_tags is not None:

        # if the element is specifically tagged 'area':'no' -> LineString
        if element_tags.get("area") == "no":
            pass

        # if the element has tags and is not tagged 'area':'no'
        # compare its tags with the polygon_features dict
        else:
            # identify common keys in element's tags and polygon_features dict
            intersecting_keys = element_tags.keys() & polygon_features.keys()

            # for each key in the intersecting keys (if any found)
            for key in intersecting_keys:

                # Get the key's value from the element's tags
                key_value = element_tags.get(key)

                # Determine if the key is for a blocklist or passlist in
                # polygon_features dict
                blocklist_or_passlist = polygon_features.get(key).get("polygon")

                # Get values for the key from the polygon_features dict
                polygon_features_values = polygon_features.get(key).get("values")

                # if all features with that key should be polygons -> Polygon
                if blocklist_or_passlist == "all":
                    is_polygon = True

                # if the key is for a blocklist i.e. tags that should not become Polygons
                elif blocklist_or_passlist == "blocklist":

                    # if the value for that key in the element is not in the blocklist -> Polygon
                    if key_value not in polygon_features_values:
                        is_polygon = True

                # if the key is for a passlist i.e. specific tags should become Polygons
                elif blocklist_or_passlist == "passlist":

                    # if the value for that key in the element is in the passlist -> Polygon
                    if key_value in polygon_features_values:
                        is_polygon = True

    return is_polygon


def _parse_relation_to_multipolygon(element: Dict, geometries: Dict) -> Dict:
    """
    Parse multipolygon from OSM relation (type:MultiPolygon).
        See more information about relations from OSM documentation:
        http://wiki.openstreetmap.org/wiki/Relation

    Args:
        element: element type "relation" from overpass response JSON
        geometries: dict containing all linestrings and polygons generated from OSM ways

    Returns
        multipolygon: dict of tags and geometry for a single multipolygon
    """
    multipolygon = dict()
    multipolygon["osmid"] = element["id"]
    multipolygon["element_type"] = "relation"

    # Parse member 'way' ids
    member_way_refs = [
        member["ref"] for member in element["members"] if member["type"] == "way"
    ]
    multipolygon["ways"] = member_way_refs

    # Add the tags
    if "tags" in element:
        for tag in element["tags"]:
            multipolygon[tag] = element["tags"][tag]

    # Extract the ways from the geometries dict using their unique id.
    # XMLs exported from the openstreetmap.org homepage with a bounding box
    # may include the relation but not the ways outside the bounding box.

    try:
        member_ways = [
            geometries[f"way/{member_way_ref}"] for member_way_ref in member_way_refs
        ]
    except KeyError as e:  # pragma: no cover

        try:
            # TODO: need to revisit how ways are parsed!! this fails in some rare cases where line strings should
            #   be combined to make a polygon
            line_strings_present = False
            for p in element["members"]:
                if len(p["geometry"]) < 3:
                    line_strings_present = True

            if line_strings_present:
                multipolygon["geometry"] = Polygon(
                    [
                        Point(pt["lon"], pt["lat"])
                        for p in element["members"]
                        for pt in p["geometry"]
                        if (p["role"] == "outer")
                    ]
                )  # [::-1]

            else:
                multipolygon["geometry"] = MultiPolygon(
                    [
                        Polygon([Point(pt["lon"], pt["lat"]) for pt in p["geometry"]])
                        for p in element["members"]
                        if (p["role"] == "outer") and (len(p["geometry"]) > 2)
                    ]
                )  # [::-1]

        except (ValueError, KeyError):
            logger.debug(
                f"{e} was not found in `geometries`.\nThe geometry for "
                f"https://www.openstreetmap.org/{element['type']}/{element['id']} was not created."
            )
            multipolygon["geometry"] = MultiPolygon()

        return multipolygon

    # Extract the nodes of those ways
    member_nodes = [[member_way["nodes"] for member_way in member_ways]]
    multipolygon["nodes"] = member_nodes

    # Assemble MultiPolygon component polygons from component LineStrings and Polygons
    outer_polygons, inner_polygons = _assemble_multipolygon_component_polygons(
        element, geometries
    )

    # Subtract inner polygons from outer polygons
    geometry = _subtract_inner_polygons_from_outer_polygons(
        element, outer_polygons, inner_polygons
    )

    multipolygon["geometry"] = geometry
    return multipolygon


def _assemble_multipolygon_component_polygons(
    element: Dict, geometries: Dict
) -> MultiPolygon:
    """
    Assemble a MultiPolygon from its component LineStrings and Polygons.
        The OSM wiki suggests an algorithm for assembling multipolygon geometries
        https://wiki.openstreetmap.org/wiki/Relation:multipolygon/Algorithm.

    This method takes a simpler approach relying on the accurate tagging
        of component ways with 'inner' and 'outer' roles as required on this page
        https://wiki.openstreetmap.org/wiki/Relation:multipolygon.

    Args:
        element: element type "relation" from overpass response JSON
        geometries: dict containing all linestrings and polygons generated from OSM ways

    Returns
        geometry: a single MultiPolygon object
    """
    outer_polygons = []
    inner_polygons = []
    outer_linestrings = []
    inner_linestrings = []

    # get the linestrings and polygons that make up the multipolygon
    for member in element["members"]:
        if member.get("type") == "way":
            # get the member's geometry from linestrings_and_polygons
            linestring_or_polygon = geometries.get(f"way/{member['ref']}")
            # sort it into one of the lists according to its role and geometry
            if (member.get("role") == "outer") and (
                linestring_or_polygon["geometry"].geom_type == "Polygon"
            ):
                outer_polygons.append(linestring_or_polygon["geometry"])
            elif (member.get("role") == "inner") and (
                linestring_or_polygon["geometry"].geom_type == "Polygon"
            ):
                inner_polygons.append(linestring_or_polygon["geometry"])
            elif (member.get("role") == "outer") and (
                linestring_or_polygon["geometry"].geom_type == "LineString"
            ):
                outer_linestrings.append(linestring_or_polygon["geometry"])
            elif (member.get("role") == "inner") and (
                linestring_or_polygon["geometry"].geom_type == "LineString"
            ):
                inner_linestrings.append(linestring_or_polygon["geometry"])

    # Merge outer linestring fragments.
    # Returns a single LineString or MultiLineString collection
    merged_outer_linestrings = linemerge(outer_linestrings)

    # polygonize each linestring separately and append to list of outer polygons
    if merged_outer_linestrings.geom_type == "LineString":
        outer_polygons += polygonize(merged_outer_linestrings)
    elif merged_outer_linestrings.geom_type == "MultiLineString":
        for merged_outer_linestring in list(merged_outer_linestrings):
            outer_polygons += polygonize(merged_outer_linestring)

    # Merge inner linestring fragments.
    # Returns a single LineString or MultiLineString collection
    merged_inner_linestrings = linemerge(inner_linestrings)

    # polygonize each linestring separately and append to list of inner polygons
    if merged_inner_linestrings.geom_type == "LineString":
        inner_polygons += polygonize(merged_inner_linestrings)
    elif merged_inner_linestrings.geom_type == "MultiLineString":
        for merged_inner_linestring in merged_inner_linestrings.geoms:
            inner_polygons += polygonize(merged_inner_linestring)

    if not outer_polygons:
        logger.debug(
            "No outer polygons were created for"
            f" https://www.openstreetmap.org/{element['type']}/{element['id']}"
        )

    return outer_polygons, inner_polygons


def _subtract_inner_polygons_from_outer_polygons(
    element: Dict, outer_polygons: List, inner_polygons: List
) -> Union[Polygon, MultiPolygon]:
    """
    Subtract inner polygons from outer polygons. Creates a Polygon or MultiPolygon with holes.

    Args:
        element: element type "relation" from overpass response JSON
        outer_polygons: list of outer polygons that are part of a multipolygon
        inner_polygons: list of inner polygons that are part of a multipolygon

    Returns
        geometry: a single Polygon or MultiPolygon
    """
    # create a new list to hold the outer polygons with the inner polygons
    # subtracted
    outer_polygons_with_holes = []

    # loop through the outer polygons subtracting the inner polygons and
    # appending to the list
    for outer_polygon in outer_polygons:
        for inner_polygon in inner_polygons:
            if inner_polygon.within(outer_polygon):
                try:
                    outer_polygon = outer_polygon.difference(inner_polygon)
                except TopologicalError:  # pragma: no cover
                    logger.debug(
                        f"relation https://www.openstreetmap.org/relation/{element['id']} caused"
                        " a TopologicalError, trying with zero buffer."
                    )
                    outer_polygon = outer_polygon.buffer(0).difference(
                        inner_polygon.buffer(0)
                    )

        # note: .buffer(0) can return either a Polygon or MultiPolygon
        # if it returns a MultiPolygon we need to extract the component
        # sub Polygons to add to outer_polygons_with_holes
        if outer_polygon.geom_type == "Polygon":
            outer_polygons_with_holes.append(outer_polygon)
        elif outer_polygon.geom_type == "MultiPolygon":
            outer_polygons_with_holes.extend(list(outer_polygon))

    # if only one polygon with holes was created, return that single polygon
    if len(outer_polygons_with_holes) == 1:
        geometry = outer_polygons_with_holes[0]
    # otherwise create a multipolygon from list of outer polygons with holes
    else:
        geometry = MultiPolygon(outer_polygons_with_holes)

    return geometry


def _buffer_invalid_geometries(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Buffer any invalid geometries remaining in the GeoDataFrame.
        Invalid geometries in the GeoDataFrame (which may accurately reproduce
        invalid geometries in OpenStreetMap) can cause the filtering to the query
        polygon and other subsequent geometric operations to fail. This function
        logs the ids of the invalid geometries and applies a buffer of zero to try
        to make them valid.

    Note: the resulting geometries may differ from the originals

    Args:
        gdf: a GeoDataFrame with possibly invalid geometries

    Returns:
        gdf: the GeoDataFrame with .buffer(0) applied to invalid geometries
    """
    # only apply the filters if the GeoDataFrame is not empty
    if not gdf.empty:

        # create a filter for rows with invalid geometries
        invalid_geometry_filter = ~gdf["geometry"].is_valid

        # if there are invalid geometries
        if invalid_geometry_filter.any():
            # get their unique_ids from the index
            invalid_geometry_ids = gdf.loc[invalid_geometry_filter].index.to_list()

            # create a list of their urls and log them
            osm_url = "https://www.openstreetmap.org/"
            invalid_geom_urls = [
                osm_url + unique_id for unique_id in invalid_geometry_ids
            ]
            logger.debug(
                f"{len(invalid_geometry_ids)} invalid geometries"
                f".buffer(0) applied to {invalid_geom_urls}",  # level=lg.WARNING,
            )

            # apply .buffer(0), help fix any invalid geometries?
            gdf.loc[invalid_geometry_filter, "geometry"] = gdf.loc[
                invalid_geometry_filter, "geometry"
            ].buffer(0)

    return gdf


def json_to_geometries(response_json: Dict) -> gpd.GeoDataFrame:
    """
    Parse JSON response from the Overpass API to a GeoDataFrame.

    Args:
        response_json: dict inlcuding a list of JSON responses from from the Overpass API

    Returns:
        gdf: GeoDataFrame of geometries and their associated tags
    """
    try:
        elements = response_json["elements"]
    except KeyError:
        raise Exception("OSM Overpass response json is invalid")

    elements_count = len(elements)

    # if there are no elements in the responses
    if elements_count == 0:
        # create an empty GeoDataFrame
        gdf = gpd.GeoDataFrame()
        gdf["geometry"] = np.nan
        gdf.set_geometry("geometry")
        gdf.crs = "epsg:4326"

        # log a warning
        logger.warning("No data elements: check query tags and location.")
        return gdf

    # else if there were elements in the response
    else:
        logger.debug(
            f"Converting {elements_count} elements in JSON responses to geometries"
        )

        # Dictionaries to hold nodes and complete geometries
        coords = dict()
        geometries = dict()

        # Set to hold the unique IDs of elements that do not have tags
        untagged_element_ids = set()

        # identify which relation types to parse to (multi)polygons
        relation_types = {"boundary", "multipolygon"}

        # extract geometries from the downloaded osm data
        # Parses the JSON of OSM nodes, ways and (multipolygon) relations
        # to dictionaries of coordinates, Shapely Points, LineStrings,
        # Polygons and MultiPolygons
        for element in elements:
            # try:

            # id numbers are only unique within element types
            # create unique id from combination of type and id
            unique_id = f"{element['type']}/{element['id']}"

            # add elements that are not nodes and that are without tags or
            # with empty tags to the untagged_element_ids set (untagged
            # nodes are not added to the geometries dict at all)
            if (element["type"] != "node") and (
                ("tags" not in element) or (not element["tags"])
            ):
                untagged_element_ids.add(unique_id)

            if element["type"] == "node":
                # Parse all nodes to coords
                coords[element["id"]] = _parse_node_to_coords(element=element)

                # If the node has tags and the tags are not empty parse it
                # to a Point. Empty check is necessary for JSONs created
                # from XML where nodes without tags are assigned tags={}
                if "tags" in element and len(element["tags"]) > 0:
                    point = _parse_node_to_point(element=element)
                    geometries[unique_id] = point

            elif element["type"] == "way":
                # Parse all ways to linestrings or polygons
                linestring_or_polygon = _parse_way_to_linestring_or_polygon(
                    element=element, coords=coords
                )
                geometries[unique_id] = linestring_or_polygon

            elif (
                element["type"] == "relation"
                and element.get("tags").get("type") in relation_types
            ):
                # parse relations to (multi)polygons
                multipolygon = _parse_relation_to_multipolygon(
                    element=element, geometries=geometries
                )
                geometries[unique_id] = multipolygon

            else:
                logger.debug(f' .. Error ! Skipping= {element["id"]}')

        logger.debug(f"{len(geometries)} geometries created in the dict")

        # remove untagged elements from the final dict of geometries
        for untagged_element_id in untagged_element_ids:
            geometries.pop(untagged_element_id, None)

        logger.debug(f"{len(untagged_element_ids)} untagged geometries removed")

        # Create GeoDataFrame
        gdf = gpd.GeoDataFrame.from_dict(geometries, orient="index")

        # ensure gdf has a geometry col before assigning crs
        if "geometry" not in gdf.columns:
            # if there is no geometry column, create a null column
            gdf["geometry"] = np.nan
        gdf.set_geometry("geometry")

        gdf.crs = "epsg:4326"  # Set default crs
        gdf = _buffer_invalid_geometries(
            gdf
        )  # Apply .buffer(0) to any invalid geometries

        # bug in geopandas <0.9 raises a TypeError if trying to plot empty
        # geometries but missing geometries (gdf['geometry'] = None) cannot be
        # projected e.g. gdf.to_crs(). Remove rows with empty (e.g. Point())
        # or missing (e.g. None) geometry, and suppress gpd warning caused by
        # calling gdf["geometry"].isna() on GeoDataFrame with empty geometries
        if not gdf.empty:
            warnings.filterwarnings("ignore", "GeoSeries.isna", UserWarning)
            gdf = gdf[~(gdf["geometry"].is_empty | gdf["geometry"].isna())].copy()
            warnings.resetwarnings()

        logger.debug(f"{len(gdf)} geometries in the final GeoDataFrame")

        return gdf.reset_index(drop=True)
