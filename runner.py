import os
import time as time_module
import platform
import subprocess
import shutil

def run(origin_path, params):
   if params['cache'] == 0:
      cpu_info = get_cpu_info()
      print("---------------------------")
      print(f"System information -- OS: {platform.system()}, CPU: {cpu_info}")
      print("---------------------------")
      walk_and_run(origin_path, pytest_py = params['pytest_py'], tags = params['tags'], verbose = params['verbose'], exclude = params['exclude'], dry_run=params['dry_run'])
   else:
      pass

""" Walk the folder and look for test and run them
"""
def walk_and_run(origin_path, pytest_py = "python2", tags = [], verbose = 0, exclude=[], dry_run=False):
   num_tests = 0
   num_passed_tests = 0
   start0 = time_module.time()
   run_tests = {}
   skip_dirs = {"__pycache__", ".git"}
   for root, d_names, f_names in os.walk(origin_path):
      d_names[:] = [d for d in d_names if d not in skip_dirs]
      need_to_run = True
      for p in exclude: # check if the folder is in the exclude list
         exclude_path = os.path.join(origin_path, p)
         if os.path.commonpath([root, exclude_path]) == exclude_path:
            # print(f"{root} is excluded")
            need_to_run = False
      if not need_to_run:
         continue
      for f in f_names:
         if f.startswith("pytest_") and f.endswith(".py"):
            success, elapsed_time = run_file(origin_path, root, f, pytest_py=pytest_py, tags=tags, verbose=verbose, dry_run=dry_run)
            if success:
               num_passed_tests += 1
               fn = os.path.join(root, f)
               run_tests[fn] = elapsed_time
            if elapsed_time > 0.0:
               num_tests += 1
   end0 = time_module.time()
   print("Test completed, %d/%d passed. Total time = %f s." % (num_passed_tests, num_tests, end0 - start0))
   cnt_long_tests = 0
   for fn, elapsed_time in run_tests.items():
      if elapsed_time > 1.0:
         cnt_long_tests += 1
   if cnt_long_tests > 0:
      print("%d tests take more than one second to run" % (cnt_long_tests))
      for fn, elapsed_time in run_tests.items():
         if elapsed_time > 1.0:
            print("  " + fn + ": %f s" % (elapsed_time))
   if len(tags) == 0:
      # save the list of tests
      ifile = open(origin_path+"/updated_list_of_test.txt", 'w')
      for fn, elapsed_time in run_tests.items():
         ifile.write(fn + '\n')
      ifile.close()

""" Run tests saved in cached list
"""
def run_cache(origin_path, pytest_py = "python2", tags = [], verbose = 0):
   ifile = open(origin_path+"/updated_list_of_test.txt", 'r')
   for line in ifile.readlines():
      pass
   ifile.close()

""" Dummy process for dry run
"""
class DummyProcess:
   def __init__(self):
      self.returncode = 0

   def communicate(self):
      return "", ""

""" Run a test file
"""
def run_file(origin_path, root, f, pytest_py = "python2", tags=[], verbose=0, dry_run=False):
   run_with_tag = (len(tags) > 0)

   use_temp_folder_for_test = False
   if os.name == "nt":
      if len(root) > 250:  # deal with the long folder name on Windows
         link_path = os.path.join(origin_path, "ztemp")
         ### try to create symlink # this solution does not work if the local folder contains symlink
         # if os.path.exists(link_path):
         #    try:
         #       shutil.rmtree(link_path)
         #    except Exception:
         #       pass
         # print(f"link_path: {link_path}")
         # cmd = ["cmd", "/c", "mklink", "/J", "\\\\?\\"+root, link_path]
         # subprocess.run(cmd, shell=False, check=True)
         ### try to copy and run in a temporary directory
         use_temp_folder_for_test = True
         shutil.copytree(root, link_path, symlinks=False, dirs_exist_ok=False)
      else:
         link_path = root
   else:
      link_path = root

   # obtain the tags of the test
   proc = subprocess.Popen([pytest_py, f, "print_tag"], cwd=link_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
   tmp_stdout, tmp_stderr = proc.communicate()
   tags_of_test = []
   if proc.returncode == 0:
      tmp = tmp_stdout.decode('ascii').splitlines()
      for t in tmp:
         if ("Tag(s):" in t) or ("Tags" in t):
            words = t.replace(",", " ").split() # tags should be separated by comma
            for i in range(1, len(words)):
               tags_of_test.append(words[i])

   if run_with_tag:
      run_mode = 0
      if ("untested" in tags_of_test) or ("Untested" in tags_of_test):
         run_mode = -1
      for tag in tags_of_test:
         if tag in tags:
            run_mode = 2
   else:
      run_mode = 1
      if ("untested" in tags_of_test) or ("Untested" in tags_of_test):
         run_mode = -1

   success = False
   elapsed_time = 0.0
   if run_mode != 0:
      if run_mode > 0:
         start = time_module.time()
         if not dry_run:
            proc = subprocess.Popen([pytest_py, f, "test"], cwd=link_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
         else:
            proc = DummyProcess()
         tmp_stdout, tmp_stderr = proc.communicate()
         if verbose:
            for line in tmp_stdout.splitlines():
               formatted_line = line.strip()  # Remove trailing newline/whitespace
               print(f"[OUTPUT] {formatted_line}")
         end = time_module.time()
         if not dry_run:
            elapsed_time = end - start
         else:
            elapsed_time = 1e-6
         if proc.returncode == 0:
            # print(fn + " passed, elapsed time = %f s" % (elapsed_time))
            print(f + " -> passed")
            success = True
         else:
            print(f + " -> failed")
            for line in tmp_stderr.splitlines():
               formatted_line = line.strip()  # Remove trailing newline/whitespace
               print(f"[ERROR] {formatted_line}")
         print("  test folder: %s" % (root))
         print("  tags:", end='')
         for tag in tags_of_test:
            print(" " + str(tag), end='')
         print('')
         print("  elapsed time: %f s" % (elapsed_time))
         print('')
      elif run_mode == -1: # ignore the untested ones -> just print the action
         if not run_with_tag:
            print(f + " -> not run (marked as untested)")
            print("  test folder: %s" % (root))
            print("  tags:", end='')
            for tag in tags_of_test:
               print(" " + str(tag), end='')
            print('')

   # clean up
   if use_temp_folder_for_test:
      print(f"{link_path} is deleted")
      shutil.rmtree(link_path)

   return success, elapsed_time

""" Print the tags associated with a test file
"""
def print_tags(origin_path, pytest_py = "python2"):
   for root, d_names, f_names in os.walk(origin_path):
      for f in f_names:
         if f.startswith("pytest_") and f.endswith(".py"):
            # obtain the tags of the test
            proc = subprocess.Popen([pytest_py, f, "print_tag"], cwd=root, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            tmp_stdout, tmp_stderr = proc.communicate()
            tags_of_test = []
            if proc.returncode == 0:
               tmp = tmp_stdout.decode('ascii').splitlines()
               for t in tmp:
                  if ("Tag(s):" in t) or ("Tags" in t):
                     words = t.replace(",", " ").split() # tags should be separated by comma
                     for i in range(1, len(words)):
                        tags_of_test.append(words[i])
            print(root + "/" + f)
            print("proc.returncode: ", proc.returncode)
            print("tags_of_test: ", tags_of_test)
            print("---------------------------------------")

""" Walk the folder and look for test and print the name and path
"""
def print_tests(origin_path, exclude=[]):
   num_tests = 0
   skip_dirs = {"__pycache__", ".git"}
   for root, d_names, f_names in os.walk(origin_path):
      d_names[:] = [d for d in d_names if d not in skip_dirs]
      d_names.sort()
      f_names.sort()
      need_to_run = True
      for p in exclude: # check if the folder is in the exclude list
         exclude_path = os.path.join(origin_path, p)
         if os.path.commonpath([root, exclude_path]) == exclude_path:
            # print(f"{root} is excluded")
            need_to_run = False
      if not need_to_run:
         continue
      for f in f_names:
         if f.startswith("pytest_") and f.endswith(".py"):
            link_path = os.path.join(root, f)
            print(f"{link_path}")
            num_tests += 1

   print("Number of tests: %d" % (num_tests))

""" Get the CPU information
"""
def get_cpu_info():
   # Get CPU model
   cpu_info = "Unknown"

   # Linux
   if platform.system() == "Linux":
       try:
           with open("/proc/cpuinfo") as f:
               for line in f:
                   if "model name" in line:
                       cpu_info = line.strip().split(":")[1].strip()
                       break
       except:
           pass

   # Windows
   elif platform.system() == "Windows":
       try:
           cpu_info = subprocess.check_output("wmic cpu get Name", shell=True)
           cpu_info = cpu_info.decode().split("\n")[1].strip()
       except:
           pass

   # macOS
   elif platform.system() == "Darwin":
       try:
           cpu_info = subprocess.check_output(["sysctl", "-n", "machdep.cpu.brand_string"])
           cpu_info = cpu_info.decode().strip()
       except:
           pass

   return cpu_info
