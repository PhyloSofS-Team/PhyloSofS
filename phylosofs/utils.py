import shutil
import os


def clear_folder(path):
    "Delete files and folders inside a folder/path."
    if not os.path.isdir(path):
        raise ValueError(str(path) + " should be a folder.")
    for root, dirs, files in os.walk(path):
        for file in files:
            os.unlink(os.path.join(root, file))
        for directory in dirs:
            shutil.rmtree(os.path.join(root, directory))


def makedirifnot(dir):
    """
    It creates a new directory if it doesn't exists.
    """
    if not os.path.exists(dir):
        os.makedirs(dir)
