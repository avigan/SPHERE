import sys
path = '/Users/avigan/Work/GitHub/pySPHERE/'
if path not in sys.path:
    sys.path.append(path)
import pysphere.IFS as IFS


root_path = '/Users/avigan/data/pySPHERE-test/IFS/'

red = IFS.IFSReduction(root_path)

# one-line reduction
# red.full_reduction()

# manual reduction
# red.init_dataset()
# red.create_static_calibrations()
red.preprocess_science()
red.process_science()
