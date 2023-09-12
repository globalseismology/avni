#!/usr/bin/env python
"""This script isolates files for a public release"""

import argparse #parsing arguments
import os
import numpy as np
import pdb
import glob

#########################################################
def main():
    parser = argparse.ArgumentParser(description='plot map-view or cross-section plots of 3D Earth models')
    parser.add_argument('-p', '--private', type=str, default=
        'docs/.private',help='List of files to make private')
    parser.add_argument('-u', '--public', type=str, default=
        'docs/.public',help='List of files to make public')
    parser.add_argument('-e', '--exclude', type=str, default=
        '.pyc,.so,.plist',help='Extensions to exclude')
    parser.add_argument('-g', '--github', action='store_true',
        help='Execute git commands for file exlusion')
    arg = parser.parse_args()


    exclude = tuple(arg.exclude.split(','))

    # Get the full paths of public files
    with open(arg.public) as f:
        paths = [line.rstrip('\n') for line in f if not line.startswith('\n') and not line.startswith('#')]
    test_list = []
    for dir_path in paths:
        if os.path.isdir(dir_path):
            for root, dirs, files in os.walk(dir_path):
                for filename in files:
                    if not filename.endswith(exclude): test_list.append(root+'/'+filename)
        else:
            test_list += glob.glob(dir_path, recursive=True)
    print('TEST PUBLIC FILES #'+str(len(test_list)))
    print(test_list)

    # Get the full paths of private files
    with open(arg.private) as f:
        paths = [line.rstrip('\n') for line in f if not line.startswith('\n') and not line.startswith('#')]
    private = []
    for dir_path in paths:
        if os.path.isdir(dir_path):
            for root, dirs, files in os.walk(dir_path):
                for filename in files: private.append(root+'/'+filename)
        else:
            private += glob.glob(dir_path, recursive=True)
    print('PRIVATE FILES #'+str(len(private)))
    print(private)

    # using list comprehension to perform task
    public = [i for i in test_list if i not in private]
    print('PUBLIC FILES #'+str(len(public)))
    print(public)

    # remove files from commit
    if arg.github:
        for root, dirs, files in os.walk('.'):
            for filename in files:
                search = (root+'/'+filename).lstrip('./')
                if search not in public:
                    print('Removing '+search)
                    os.system('git rm --cached '+search)
    return

if __name__== "__main__":
    main()
