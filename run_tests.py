import sys
import os

sys.path.append(os.getcwd() + "/..")
import runner

origin_path = os.getcwd()

tags = []
verbose = 0
cache = 0
if len(sys.argv) > 2:
  for i in range(2, len(sys.argv)):
    if sys.argv[i] == "--verbose": # allow verbose by command line argument --verbose
      verbose = 1
    elif sys.argv[i] == "--cache": # allow to run from saved list of tests
      cache = 1
    else:
      tags.append(sys.argv[i])

if len(tags) > 0:
    print(f"Tags to be tested: {tags}, verbose = {verbose}")
else:
    print(f"Tags to be tested: all, verbose = {verbose}")

# Check if the file exists
exclude_file = ".appignore"
if os.path.isfile(exclude_file):
    # Read lines, strip whitespace, and ignore empty lines
    with open(exclude_file, "r") as f:
        exclude_names = [line.strip() for line in f if line.strip()]
else:
    exclude_names = []  # No file, no exclusions

params = {}
params['pytest_py'] = sys.argv[1]
params['tags']      = tags
params['verbose']   = verbose
params['cache']     = cache
params['exclude']   = exclude_names
runner.run(origin_path, params)
