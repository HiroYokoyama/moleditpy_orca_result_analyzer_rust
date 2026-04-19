import os


def get_default_export_path(base_path, suffix="_analyzed", extension=""):
    """
    Generates a default export path based on the input file path.
    Example: 'job.out' -> 'job_analyzed.csv'
    """
    if not base_path:
        return ""
    dirname = os.path.dirname(base_path)
    filename_base = os.path.splitext(os.path.basename(base_path))[0]
    new_filename = f"{filename_base}{suffix}{extension}"
    return os.path.join(dirname, new_filename)
