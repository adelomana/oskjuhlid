#!/bin/bash

rsync -azPvh -e "ssh -i ~/.ssh/elja" adrian@elja.hi.is:/users/home/adrian/projects/oskjuhlid/src .
rm src/*~
