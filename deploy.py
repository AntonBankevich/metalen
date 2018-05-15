import os
import shutil
import sys


def Deploy(from_path, to_path, filter = lambda x: True):
    if os.path.isfile(from_path):
        if filter(from_path):
            print "Deploying:", from_path
            shutil.copy(from_path, to_path)
    else:
        os.makedirs(to_path)
        for file in os.listdir(from_path):
            Deploy(os.path.join(from_path, file), os.path.join(to_path, file), filter)

if __name__ == "__main__":
    deploy_path = sys.argv[1]
    os.makedirs(deploy_path)
    copyall = ["GPLv2.txt", "LICENSE", "manual.html", "MetaLen.py"]
    copypy = ["py"]
    for path in copyall:
        Deploy(path, os.path.join(deploy_path, path))
    for path in copypy:
        Deploy(path, os.path.join(deploy_path, path), lambda name: name.endswith(".py"))
