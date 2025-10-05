# utils/logging_config.py
import logging
import sys
from logging.handlers import RotatingFileHandler

def setup_logging(verbose: bool = False, log_file: str = None):
    """
    Configure logging with console and optional file handler.

    :param verbose: If True, set level to DEBUG
    :param log_file: Optional path to log file
    """
    level = logging.DEBUG if verbose else logging.INFO
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    
    # Root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(level)
    root_logger.addHandler(console_handler)
    
    # File handler if specified
    if log_file:
        file_handler = RotatingFileHandler(log_file, maxBytes=10**6, backupCount=5)
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
        root_logger.addHandler(file_handler)
    
    logging.info("Logging configured")
