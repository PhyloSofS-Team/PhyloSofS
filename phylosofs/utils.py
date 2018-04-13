import shutil
import os

def clear_folder(path):
    "Delete files and folders inside a folder/path."
    if not os.path.isdir(path):
        raise ValueError(str(path) + " should be a folder.")
    for root, dirs, files in os.walk(path):
        for f in files:
            os.unlink(os.path.join(root, f))
        for d in dirs:
            shutil.rmtree(os.path.join(root, d))

