import subprocess
import sys
import time
import random

def log(message, verbose=True):
    if verbose:
        print(message, flush=True)


def run_cmd(cmd, shell=False, verbose = True):
    log(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}", verbose)
    result = subprocess.run(cmd, shell=shell)
    if result.returncode != 0:
        print(f"Command failed: {cmd}", file=sys.stderr)
        sys.exit(result.returncode)

def run_cmd_raise(cmd, shell=False, verbose = True):
    """Run command and raise RuntimeError on failure instead of exiting.
    
    Use this in library code paths where exceptions are preferred over sys.exit().
    CLI code should catch exceptions and exit with code 1.
    """
    log(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}", verbose)
    result = subprocess.run(cmd, shell=shell)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {result.returncode}: {cmd}")

def run_cmd_retry(cmd, tries=5, base=5, cap=300, verbose=True):
    """Run a command, retrying on failure with exponential backoff + full jitter.

    Full jitter (sleep in [0, base*2**attempt]) de-synchronizes many jobs that fail at once.
    """
    for attempt in range(tries):
        log(f"Running: {' '.join(cmd) if isinstance(cmd, list) else cmd}", verbose)
        if subprocess.run(cmd).returncode == 0:
            return
        if attempt < tries - 1:
            delay = min(cap, random.uniform(0, base * 2 ** attempt))
            log(f"Command failed (attempt {attempt + 1}/{tries}); retrying in {delay:.0f}s", verbose)
            time.sleep(delay)
    raise RuntimeError(f"Command failed after {tries} attempts: {cmd}")


def get_top_n_chromosomes(fai_path, n=2):
    chroms = []
    with open(fai_path) as f:
        for line in f:
            fields = line.strip().split('\t')
            chroms.append((fields[0], int(fields[1])))
    chroms.sort(key=lambda x: -x[1])
    return [c[0] for c in chroms[:n]]

