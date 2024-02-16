from setuptools import find_packages, setup
import os
import glob as gb

# get requirements for installation
lib_folder = os.path.dirname(os.path.realpath(__file__))
requirement_path = lib_folder + '/requirements.txt'
install_requires = [] # Here we'll get: ["gunicorn", "docutils>=0.3", "lxml==0.5a7"]
if os.path.isfile(requirement_path):
    with open(requirement_path) as f:
        install_requires = f.read().splitlines()
        print("install requires:")
        print(install_requires)

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='anisotropicMAIRE',
    packages=find_packages(exclude='pytests'),
    version='1.0.0',
    description='Python library for running the anisotropic version of MAIRE+',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Chris S. W. Davis',
    author_email='ChrisSWDavis@gmail.com',
    license='GNU General Public License v3.0',
    install_requires=install_requires,
    setup_requires=['pytest-runner','wheel'],
    tests_require=['pytest'],
)
