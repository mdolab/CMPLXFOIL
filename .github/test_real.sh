#!/bin/bash
set -e
cd tests
testflo -n 1 -v --coverage --coverpkg pyxlight
