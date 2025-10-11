import os
import subprocess

pytest_py = "python2"
f = "one_patch.py"
root = os.getcwd()
proc = subprocess.Popen([pytest_py, f, "print_tag"], cwd=root, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
tmp_stdout, tmp_stderr = proc.communicate()
# print("output:")
# print(tmp_stdout)
# print("error:")
# print(tmp_stderr)
print("returncode:")
print(proc.returncode)
