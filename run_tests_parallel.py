import sys
import os
import time as time_module
import subprocess
from datetime import datetime
import socket

import prunner

def send_email(recipient, subject, body):
    # Construct a basic RFC 822 message
    message = f"To: {recipient}\nSubject: {subject}\n\n{body}"

    try:
        # Run msmtp command
        # -t tells msmtp to read the recipient from the message headers
        process = subprocess.Popen(
            ['msmtp', '-t'],
            stdin=subprocess.PIPE,
            text=True
        )
        process.communicate(input=message)

        if process.returncode == 0:
            print(f"Test results sent successfully to {recipient}")
        else:
            print(f"Error: msmtp exited with code {process.returncode}")

    except FileNotFoundError:
        print("Error: msmtp is not installed or not in your PATH.")

def main(params):
    origin_path = params['origin']
    num_cores   = params['num_cores']
    pytest_py   = params['pytest_py']
    tags        = params['tags']
    application = params['application']
    recipient   = params['recipient']
    runner      = params['runner']

    ###

    start_time = time_module.perf_counter()

    test_files, untest_files, untag_files = prunner.collect_tests_parallel(origin_path, \
        pytest_py=pytest_py, \
        max_workers=num_cores, \
        application=application)

    end_time = time_module.perf_counter()

    # print("test files:")
    # for f in test_files:
    #     print(f)
    # print("untest files:")
    # for f in untest_files:
    #     print(f)
    # print("untag files:")
    # for f in untag_files:
    #     print(f)

    print("Number of test files: %d" % (len(test_files)))
    print("Number of untest files: %d" % (len(untest_files)))
    print("Number of untag files: %d" % (len(untag_files)))
    print("Collection time: %.2e s" % (end_time - start_time))

    if len(tags) == 0:
        # if no tag is specified, then all tagged and untagged files are tested
        all_test_files = test_files + untag_files
    else:
        all_test_files = []
        # else only include the test with appropriate tag
        for file_info in test_files:
            has_tag = False
            tags_of_test = file_info[1]
            for tag in tags:
                if tag in tags_of_test:
                    has_tag = True
                    break
            if has_tag:
                all_test_files.append(file_info)
        if len(all_test_files) > 0:
            print("Number of matching tag files: %d" % (len(all_test_files)))

    start_time = time_module.perf_counter()
    log, result = prunner.run(all_test_files, pytest_py=pytest_py, num_proc=num_cores)
    end_time = time_module.perf_counter()

    test_message = log
    test_message += "\n"
    # print(result)
    result.sort(key=lambda x: x[1])

    num_long_tests = 0
    for r in result:
        if r[1] > 1.0:
            num_long_tests += 1

    if num_long_tests > 0:
        test_message += ("%d long tests were performed (taking more than one second to run)\n" % num_long_tests)
        for r in result:
            if r[1] > 1.0:
                test_message += (f"  %s: %.3e s\n" % (r[0], r[1]))

    if len(untest_files) > 0:
        test_message += ("List of untested files:\n")
        for f in untest_files:
            test_message += (f"  %s\n" % (f))
    test_message += ("Test completed, %d/%d passed. Total time = %.3e s.\n" % (len(result), len(all_test_files), end_time - start_time))
    print(test_message)

    if recipient != None:
        now = datetime.now()
        hostname = socket.gethostname()
        if len(result) < len(all_test_files):
            test_subject = "🚨 ALERT! 🚨 "
        else:
            test_subject = ""
        test_subject += "Prisma test result |"
        test_subject += (f" host: {socket.gethostname()} |")
        test_subject += (f" date: {now.year}/{now.month}/{now.day} |")
        test_subject += (f" %d/%d passed |" % (len(result), len(all_test_files)))
        if runner != None:
            test_subject += (f" runner: {runner} |")
        send_email(recipient, test_subject, test_message)

if __name__ == "__main__":
    tags        = []
    verbose     = 0
    cache       = 0
    num_cores   = os.cpu_count() - 1
    application = "all"
    recipient   = None
    runner      = None
    if len(sys.argv) > 2:
      for i in range(2, len(sys.argv)):
        if sys.argv[i] == "--verbose": # allow verbose by command line argument --verbose
          verbose = 1
        elif sys.argv[i] == "--cache": # allow to run from saved list of tests
          cache = 1
        elif "--numcores=" in sys.argv[i]: # specify the number of processes to run tests
          num_cores = int(sys.argv[i].split('=')[-1])
        elif "--application=" in sys.argv[i]: # specify the tests in an application
          application=sys.argv[i].split('=')[-1]
        elif "--recipient=" in sys.argv[i]: # specify the email to sent the test results to
          recipient=sys.argv[i].split('=')[-1]
        elif "--runner=" in sys.argv[i]: # specify the name of the runner to run the test
          runner=sys.argv[i].split('=')[-1]
        else:
          tags.append(sys.argv[i])

    if len(tags) > 0:
        print(f"Tags to be tested: {tags}, application = {application}, verbose = {verbose}")
    else:
        print(f"Tags to be tested: all, application = {application}, verbose = {verbose}")

    # Check if the file exists
    exclude_file = ".appignore"
    if os.path.isfile(exclude_file):
        # Read lines, strip whitespace, and ignore empty lines
        with open(exclude_file, "r") as f:
            exclude_names = [line.strip() for line in f if line.strip()]
    else:
        exclude_names = []  # No file, no exclusions

    params = {}
    params['origin']    = os.getcwd()
    params['num_cores'] = num_cores
    params['pytest_py'] = sys.argv[1]
    params['tags']      = tags
    params['verbose']   = verbose
    params['cache']     = cache
    params['exclude']   = exclude_names
    params['application']   = application
    params['recipient'] = recipient
    params['runner']    = runner
    params['dry_run']   = False # enable this to NOT run the actual test (assume passing)

    ######

    main(params)
