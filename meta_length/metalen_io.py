import gzip
import os
import shutil

import sys


def process_readline(line, is_python3=sys.version.startswith('3.')):
    if is_python3:
        return str(line, 'utf-8')
    return line

def ensure_dir_existence(dirname):
    if os.path.isfile(dirname):
        os.remove(dirname)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

def recreate_dir(dirname):
    if os.path.exists(dirname):
        shutil.rmtree(dirname)
    os.makedirs(dirname)

def universal_read(file):
    if file.endswith("gz"):
        return gzip.open(file, "r")
    else:
        return open(file, "r")

def online_sys_call(cmd, err_file):
    import shlex
    import subprocess

    if isinstance(cmd, list):
        cmd_list = cmd
    else:
        cmd_list = shlex.split(cmd)

    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=open(err_file, "a"))

    return proc.stdout

def universal_sys_call(cmd, log, out_filename=None, err_filename=None, cwd=None):
    '''
    Runs cmd and redirects stdout to out_filename (if specified), stderr to err_filename (if specified), or to log otherwise
    '''
    import shlex
    import subprocess
    if isinstance(cmd, list):
        cmd_list = cmd
    else:
        cmd_list = shlex.split(cmd)
    if out_filename:
        stdout = open(out_filename, 'w')
    else:
        stdout = subprocess.PIPE
    if err_filename:
        stderr = open(err_filename, 'w')
    else:
        stderr = subprocess.PIPE

    proc = subprocess.Popen(cmd_list, stdout=stdout, stderr=stderr, cwd=cwd)

    if log and (not out_filename or not err_filename):
        while not proc.poll():
            if not out_filename:
                line = process_readline(proc.stdout.readline()).rstrip()
                if line:
                    log.info(line)
            if not err_filename:
                line = process_readline(proc.stderr.readline()).rstrip()
                if line:
                    log.info(line)
            if proc.returncode is not None:
                break

        if not out_filename:
            for line in proc.stdout.readlines():
                if line != '':
                    log.info(process_readline(line).rstrip())
        if not err_filename:
            for line in proc.stderr.readlines():
                if line != '':
                    log.info(process_readline(line).rstrip())
    else:
        proc.wait()

    if out_filename:
        stdout.close()
    if err_filename:
        stderr.close()
    if proc.returncode:
        os.error('system call for: "%s" finished abnormally, err code: %d' % (cmd, proc.returncode), log)
