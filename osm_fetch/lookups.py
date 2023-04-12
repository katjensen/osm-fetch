""" Dictionary of tags to determine if closed ways should become polygons

Note: The dictionary is based on the JSON linked to from the following page:
    https://wiki.openstreetmap.org/wiki/Overpass_turbo/Polygon_Features
"""

POLYGON_FEATURES = {
    "building": {"polygon": "all"},
    "highway": {
        "polygon": "passlist",
        "values": ["services", "rest_area", "escape", "elevator"],
    },
    "natural": {
        "polygon": "blocklist",
        "values": ["coastline", "cliff", "ridge", "arete", "tree_row"],
    },
    "landuse": {"polygon": "all"},
    "waterway": {
        "polygon": "passlist",
        "values": ["riverbank", "dock", "boatyard", "dam"],
    },
    "amenity": {"polygon": "all"},
    "leisure": {"polygon": "all"},
    "barrier": {
        "polygon": "passlist",
        "values": ["city_wall", "ditch", "hedge", "retaining_wall", "spikes"],
    },
    "railway": {
        "polygon": "passlist",
        "values": ["station", "turntable", "roundhouse", "platform"],
    },
    "area": {"polygon": "all"},
    "boundary": {"polygon": "all"},
    "man_made": {
        "polygon": "blocklist",
        "values": ["cutline", "embankment", "pipeline"],
    },
    "power": {
        "polygon": "passlist",
        "values": ["plant", "substation", "generator", "transformer"],
    },
    "place": {"polygon": "all"},
    "shop": {"polygon": "all"},
    "aeroway": {"polygon": "blocklist", "values": ["taxiway"]},
    "tourism": {"polygon": "all"},
    "historic": {"polygon": "all"},
    "public_transport": {"polygon": "all"},
    "office": {"polygon": "all"},
    "building:part": {"polygon": "all"},
    "military": {"polygon": "all"},
    "ruins": {"polygon": "all"},
    "area:highway": {"polygon": "all"},
    "craft": {"polygon": "all"},
    "golf": {"polygon": "all"},
    "indoor": {"polygon": "all"},
}


""" Some additional helpful tag restriction definitions """

MANMADE_TAG_RESTRICTIONS = {
    "aeroway": [
        "aerodrome",
        "apron",
        "gate",
        "hangar",
        "helipad",
        "heliport",
        "navigationaid",
        "runway",
        "spaceport",
        "taxiway",
        "terminal",
    ],
    "amenity": ["parking"],
    "building": None,  # grab all!!
    "highway": [
        "motorway",
        "trunk",
        "primary",
        "secondary",
        "tertiary",
        "unclassified",
        "residential",
        "motorway_link",
        "trunk_link",
        "primary_link",
        "secondary_link",
        "tertiary_link",
        "living_street",
        "service",
        "pedestrian",
        "bus_guideway",
        "escape",
        "raceway",
        "busway",
        "rest_area",
        "services",
    ],
    "man_made": None,
    "landuse": [
        "farmyard",
        "construction",
        "industrial",
        "commercial",
        "retail",
        "quarry",
        "railway",
        "garages",
        "residential",
        "brownfield",
    ],
    "railway": [
        "construction",
        "disused",
        "funicular",
        "light_rail",
        "monorail",
        "narrow_gauge",
        "rail",
        "tram",
    ],
}

WATER_TAG_RESTRICTIONS = {
    "landuse": [
        "salt_pond",
        "aquaculture",
        "basin",
        "reservoir",
    ],
    "natural": ["water"],
    "water": None,
    "waterway": None,
}
