import os

initials = "fs"

ROOT_DIR = os.path.dirname(os.path.abspath(__file__)) + "/"
PUB_PATH = ROOT_DIR + "%s_paper/" % initials
FIG_FOLDER = "%s-figures" % initials
PUB_FIG_PATH = PUB_PATH + FIG_FOLDER + "/"
TEMP_FIG_PATH = ROOT_DIR + "temp/"
MODULE_DATA_PATH = ROOT_DIR + "data/"  # no initials, consistent with python science environment from PyCharm

import user_paths as up

if hasattr(up, "PM_STUDY_RESULTS"):
    PM_STUDY_RESULTS = up.PM_STUDY_RESULTS
else:
    PM_STUDY_RESULTS = "<set-path-to-PM_STUDY_RESULTS>"
