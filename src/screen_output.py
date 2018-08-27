#!/usr/bin/python 
#blender build scripts
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def log_out(string):
    print bcolors.OKBLUE+string+bcolors.ENDC


def run_out(string):
    print bcolors.OKGREEN+string+bcolors.ENDC

def err_out(string):
    print bcolors.FAIL+string+bcolors.ENDC

if __name__ == '__main__':
    print bcolors.WARNING+bcolors.BOLD+"Warning: No active remain"+bcolors.ENDC
