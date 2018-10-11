import os
from setuptools import setup, find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

packages = [
    "molmimic",
    "molmimic.pdb_tools",
    "molmimic.common",
    "molmimic.generate_data",
    "molmimic.keras_model",
    "molmimic.parsers",
    "molmimic.scratch",
    "molmimic.torch_model",
    "molmimic.util",
    "molmimic.visualize",
]

setup(
    name = "molmimic",
    version = "0.0.1",
    author = "Eli Draizen",
    author_email = "edraizen@gmail.com",
    packages=packages,
    long_description=read('README.md'),
    install_requires=read("requirements.txt").splitlines()
)
