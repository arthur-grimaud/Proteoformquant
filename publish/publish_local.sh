#!/usr/bin/env bash
python setup.py sdist bdist_wheel && pip install dist/*.whl --force-reinstall
python clean.py