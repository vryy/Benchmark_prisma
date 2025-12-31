import os
import time as time_module
import subprocess
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# def run(test_files, pytest_py="python2", num_proc=4):
#     total_tasks = len(test_files)

#     if total_tasks == 0:
#         print("No tests found to run.")
#         return ""

#     with multiprocessing.Manager() as manager:
#         unified_log_list = manager.list()
#         unified_result_list = manager.list()
#         progress = manager.Value('i', 0)
#         lock = manager.Lock()

#         print(f"\nStarting execution of {total_tasks} files using {num_proc} processes...")

#         with ProcessPoolExecutor(max_workers=num_proc) as executor:
#             # Map the files to the runner function
#             futures = [
#                 executor.submit(run_file, f, unified_log_list, unified_result_list, progress, lock, pytest_py=pytest_py)
#                 for f in test_files
#             ]

#             # Visual Progress Bar
#             while progress.value < total_tasks:
#                 curr = progress.value
#                 percent = int((curr / total_tasks) * 100)
#                 bar = "#" * (percent // 5) + "." * (20 - (percent // 5))
#                 print(f"\rExecuting: [{bar}] {percent}% ({curr}/{total_tasks})", end="", flush=True)
#                 import time
#                 time.sleep(0.2)

#             print(f"\rExecuting: [{"#"*20}] 100% ({total_tasks}/{total_tasks})")

#         # Join all captured logs into one large string
#         return "\n".join(unified_log_list), list(unified_result_list)

def run(test_files, pytest_py="python2", num_proc=4):
    total_tasks = len(test_files)

    if total_tasks == 0:
        print("No tests found to run.")
        return "", []

    with multiprocessing.Manager() as manager:
        unified_log_list = manager.list()
        unified_result_list = manager.list()
        progress = manager.Value('i', 0)
        lock = manager.Lock()

        print(f"\nStarting execution of {total_tasks} files using {num_proc} processes...")

        with ProcessPoolExecutor(max_workers=num_proc) as executor:
            futures = [
                executor.submit(run_file, f, unified_log_list, unified_result_list, progress, lock, pytest_py=pytest_py)
                for f in test_files
            ]

            # --- TQDM Implementation ---
            with tqdm(total=total_tasks, desc="Executing Tests", unit="file") as pbar:
                last_progress = 0
                while last_progress < total_tasks:
                    # Get the current shared value
                    current_progress = progress.value

                    # Calculate the step (how many finished since the last loop)
                    step = current_progress - last_progress

                    if step > 0:
                        pbar.update(step)
                        last_progress = current_progress

                    # Check if all futures are actually done to avoid infinite loop
                    if all(f.done() for f in futures):
                        # Final sync to ensure bar hits 100%
                        pbar.update(total_tasks - last_progress)
                        break
            # ---------------------------

        # Convert managed lists to standard lists before returning
        return "\n".join(list(unified_log_list)), list(unified_result_list)

def run_file(file_info, shared_log, shared_result, progress_counter, lock, pytest_py="python2", dry_run=False):
    # 1. Suppress macOS crash windows for this specific subprocess
    env = os.environ.copy()
    env["NSUnbufferedIO"] = "YES"

    file_path = file_info[0]
    tags_of_test = file_info[1]
    file_name = os.path.basename(file_path)
    dir_path = os.path.dirname(file_path)

    # Initialize output so it's always defined
    output = f"{file_name}"

    try:
        # 2. Run the test
        # Note: we use 'result' here to match the variable name below
        start_time = time_module.perf_counter()

        if not dry_run:
            proc = subprocess.Popen([pytest_py, file_name, "test"], cwd=dir_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
        else:
            proc = DummyProcess()
        tmp_stdout, tmp_stderr = proc.communicate()

        end_time = time_module.perf_counter()

        if proc.returncode == 0:
            output += " -> passed\n"
            output += f"  test folder: {dir_path}\n"
            output += f"  tags:"
            for tag in tags_of_test:
                output += f" {tag}"
            output += f"\n"
            output += f"  elapsed time: %.3e s\n" % (end_time - start_time)
        else:
            output += " -> failed\n"
            output += f"  test folder: {dir_path}\n"
            output += f"  tags:"
            for tag in tags_of_test:
                output += f" {tag}"
            output += f"\n"
            output += f"\n{tmp_stderr}\n"

    except Exception as e:
        output += f"CRITICAL ERROR starting process: {str(e)}\n"

    finally:
        # 3. Always append the result we have gathered
        shared_log.append(output)

        if proc.returncode == 0:
            shared_result.append([file_path, end_time - start_time])

        # 4. Always update the progress bar so it doesn't hang at 99%
        with lock:
            progress_counter.value += 1

def get_file_tag(file_info, shared_lists, progress_counter, lock):
    """
    Function executed by child processes.
    file_info: (root, f, pytest_py)
    """
    root, f, pytest_py = file_info
    test_files, untest_files, untag_files = shared_lists
    link_path = os.path.join(root, f)

    # Run the subprocess
    proc = subprocess.Popen(
        [pytest_py, f, "print_tag"],
        cwd=root,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    tmp_stdout, _ = proc.communicate()

    tags_of_test = []
    if proc.returncode == 0:
        tmp = tmp_stdout.decode('ascii').splitlines()
        for t in tmp:
            if ("Tag(s):" in t) or ("Tags" in t):
                words = t.replace(",", " ").split()
                tags_of_test.extend(words[1:])
    else:
        # if the test is failed to load, then it is not marked for testing
        tags_of_test.append('untested')

    # Thread-safe appending to shared lists
    if len(tags_of_test) == 0:
        untag_files.append([link_path, []])
    elif any(tag.lower() == "untested" for tag in tags_of_test):
        untest_files.append(link_path)
    else:
        test_files.append([link_path, tags_of_test])

    # Update global progress
    with lock:
        progress_counter.value += 1

# def collect_tests_parallel(origin_path, exclude=[], pytest_py="python2", max_workers=4):
#     # 1. First, quickly crawl the directories to find all files (the "Plan")
#     files_to_process = []
#     skip_dirs = {"__pycache__", ".git"}

#     for root, d_names, f_names in os.walk(origin_path):
#         d_names[:] = [d for d in d_names if d not in skip_dirs]

#         # Check exclusion
#         if any(os.path.commonpath([root, os.path.join(origin_path, p)]) == os.path.join(origin_path, p) for p in exclude):
#             continue

#         for f in f_names:
#             if f.startswith("pytest_") and f.endswith(".py"):
#                 files_to_process.append((root, f, pytest_py))

#     total_files = len(files_to_process)
#     if total_files == 0:
#         return [], [], []

#     # 2. Setup shared memory via Manager
#     with multiprocessing.Manager() as manager:
#         test_files = manager.list()
#         untest_files = manager.list()
#         untag_files = manager.list()
#         progress = manager.Value('i', 0)
#         lock = manager.Lock()

#         shared_lists = (test_files, untest_files, untag_files)

#         # 3. Execute in parallel
#         with ProcessPoolExecutor(max_workers=max_workers) as executor:
#             # Start processes
#             futures = [executor.submit(get_file_tag, info, shared_lists, progress, lock)
#                       for info in files_to_process]

#             # 4. Progress Bar Loop
#             while progress.value < total_files:
#                 curr = progress.value
#                 percent = int((curr / total_files) * 100)
#                 bar = "#" * (percent // 5) + "." * (20 - (percent // 5))
#                 print(f"\rCollecting test tags: [{bar}] {percent}% ({curr}/{total_files})", end="", flush=True)
#                 if all(f.done() for f in futures): break
#                 import time
#                 time.sleep(0.1)

#             print(f"\rCollecting test tags: [{"#"*20}] 100% ({total_files}/{total_files})")

#         # Convert managed lists back to standard python lists before returning
#         return list(test_files), list(untest_files), list(untag_files)

def collect_tests_parallel(origin_path, exclude=[], pytest_py="python2", max_workers=4):
    # 1. Quickly crawl the directories to find all files
    files_to_process = []
    skip_dirs = {"__pycache__", ".git"}

    for root, d_names, f_names in os.walk(origin_path):
        d_names[:] = [d for d in d_names if d not in skip_dirs]

        if any(os.path.commonpath([root, os.path.join(origin_path, p)]) == os.path.join(origin_path, p) for p in exclude):
            continue

        for f in f_names:
            if f.startswith("pytest_") and f.endswith(".py"):
                files_to_process.append((root, f, pytest_py))

    total_files = len(files_to_process)
    if total_files == 0:
        return [], [], []

    # 2. Setup shared memory via Manager
    with multiprocessing.Manager() as manager:
        test_files = manager.list()
        untest_files = manager.list()
        untag_files = manager.list()
        progress = manager.Value('i', 0)
        lock = manager.Lock()

        shared_lists = (test_files, untest_files, untag_files)

        # 3. Execute in parallel
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            futures = [executor.submit(get_file_tag, info, shared_lists, progress, lock)
                      for info in files_to_process]

            # 4. TQDM Progress Bar
            with tqdm(total=total_files, desc="Collecting test tags", unit="file") as pbar:
                last_val = 0
                while last_val < total_files:
                    current_val = progress.value
                    step = current_val - last_val
                    if step > 0:
                        pbar.update(step)
                        last_val = current_val

                    if all(f.done() for f in futures):
                        # Final catch-up if needed
                        pbar.update(total_files - last_val)
                        break

        # Convert managed lists back to standard python lists before the manager exits
        final_test = list(test_files)
        final_untest = list(untest_files)
        final_untag = list(untag_files)

    return final_test, final_untest, final_untag
