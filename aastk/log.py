import logging
import sys
from datetime import datetime
from pathlib import Path
from .util import *


def logger_setup(silent: bool = False, output_dir: str = None):
	timestamp = datetime.now().isoformat(timespec='seconds').replace(":", "-")
	subprocess_log_filename = f"aastk-subprocesses-{timestamp}.log"

	# If output_dir is provided, create the directory if it doesn't exist
	if output_dir:
		subprocess_log_filepath = ensure_path(path=output_dir, target=subprocess_log_filename)
	else:
		subprocess_log_filepath = subprocess_log_filename

	logger = logging.getLogger()
	if silent:
		logger.setLevel(logging.CRITICAL + 1)
	else:
		logger.setLevel(logging.INFO)

	logger.handlers.clear()

	# StreamHandler to print INFO level messages to stdout
	console_handler = logging.StreamHandler(sys.stdout)
	console_handler.setLevel(logging.INFO)
	console_handler.addFilter(lambda record: record.levelno != 99)
	console_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))

	# FileHandler to write custom level messages (for subprocesses) to file
	file_handler = logging.FileHandler(subprocess_log_filepath, delay=True)
	file_handler.setLevel(99)
	file_handler.setFormatter(logging.Formatter('%(asctime)s SUBPROCESS %(message)s'))

	logger.addHandler(console_handler)
	logger.addHandler(file_handler)

	return logger