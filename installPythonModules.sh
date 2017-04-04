#!/bin/bash
set -e -x

sudo rm /usr/bin/python
sudo ln -s /usr/bin/python2.7 /usr/bin/python

sudo -H pip install ipython[all]
sudo -H pip install requests[security]

export AWS_ACCESS_KEY_ID=AKIAJLG2HIPH47IWGK4A
export AWS_SECRET_ACCESS_KEY=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