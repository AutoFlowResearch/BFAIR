# Logging utils for the BFAIR modules.

import logging


def get_logger(name):
    """
    Returns a logger with the specified name.

    Paramters
    ---------
    name : str
        Name of the logger.

    Returns
    -------
    logging.Logger
        A logger object.
    """
    logger = logging.getLogger(name)
    if not logger.hasHandlers():
        handler = logging.StreamHandler()
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(logging.Formatter("[%(levelname)s]: %(message)s"))
    return logger
