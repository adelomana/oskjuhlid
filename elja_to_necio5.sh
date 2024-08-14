#!/bin/bash

# this is a command to transfer files from elja to necio5. Should be run in necio5. OK from IT on 2024.05.27

# rsync source target
# -e specifies a remote shell

time rsync -azPvh --bwlimit=40000  -e "ssh -i ~/.ssh/elja" adrian@elja.hi.is:/hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siCtrl_1 /Users/adrian/tmp/.
time rsync -azPvh --bwlimit=40000  -e "ssh -i ~/.ssh/elja" adrian@elja.hi.is:/hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siCtrl_2 /Users/adrian/tmp/.
time rsync -azPvh --bwlimit=40000  -e "ssh -i ~/.ssh/elja" adrian@elja.hi.is:/hpcdata/Mimir/adrian/research/oskjuhlid/data/raw/siCtrl_3 /Users/adrian/tmp/.
