import re

def get_proc_memory_gb(pid):
    virt, res = None, None
    try:
        with open(f'/proc/{pid}/status', 'r') as fh:
            contents = fh.read()
        virt = int(re.search(r'VmSize:[^\d]+(\d+)', contents).group(1)) / 1e6
        res = int(re.search(r'VmRSS:[^\d]+(\d+)', contents).group(1)) / 1e6
    finally:
        return virt, res
