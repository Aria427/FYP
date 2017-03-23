#!/bin/bash
set -e -x

sudo rm /usr/bin/python
sudo ln -s /usr/bin/python2.7 /usr/bin/python

sudo -H pip install ipython[all]

sudo -H pip install smart_open
sudo -H pip install requests[security]