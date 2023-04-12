class OverpassBadRequest(Exception):
    """Error 400: Syntax error."""


class OverpassMoved(Exception):
    """Error 302: Moved to new location?"""


class OverpassTooManyRequests(Exception):
    """Error 429: Too many requests."""


class OverpassGatewayTimeout(Exception):
    """Error 504: Too much load."""
