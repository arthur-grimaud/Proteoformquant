python setup.py sdist bdist_wheel
pip uninstall -y fragannot
pip install --find-links=dist\ fragannot
python clean.py