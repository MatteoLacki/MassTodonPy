import os
from MassTodonPy.Parsers.Paths import parse_path


def create_folder_if_needed(path):
    """Creates folder according to the path, if necessary."""
    file_path, file_name, file_ext = parse_path(path)
    if not os.path.exists(file_path):
        os.makedirs(file_path)
