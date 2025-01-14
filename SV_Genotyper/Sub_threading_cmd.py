import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
def execute_commands(cmd):
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    return stdout.decode('utf-8'), stderr.decode('utf-8'), process.returncode, cmd

def threading_cmd(cmds, log, worker=10):
    with ThreadPoolExecutor(max_workers=worker) as executor:
        futures = [executor.submit(execute_commands, cmd) for cmd in cmds]
        results =[]
        for future in as_completed(futures):
            stdout, stderr, returncode, cmd = future.result()  # Wait on each result future:
            if returncode != 0:
                print(f"Error executing command: {cmd}" , file=log)
                print(stderr, file=log)
            else:
                print(f"Output of command '{cmd}':", file=log)
                print(stdout,file=log)
            results.append((stdout,stderr,returncode,cmd))
    return results

