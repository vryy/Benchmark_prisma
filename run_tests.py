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

print("Tags to be tested:", tags, ", verbose =", verbose)

runner.run(origin_path, {'pytest_py': sys.argv[1], 'tags': tags, 'verbose': verbose, 'cache': cache})
