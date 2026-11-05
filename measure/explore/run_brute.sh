#!/bin/bash
PY=/Users/timrozday/Documents/kmer-functional-profiler/.venv/bin/python
$PY brute.py truth/true150.json truth/true150.out.json "$(cat seeds.json)"
$PY brute.py truth/bg60.json truth/bg60.out.json "$(cat seeds.json)"
