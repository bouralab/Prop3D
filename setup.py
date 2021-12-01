import os
from setuptools import setup, find_packages

def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname), "rb") as f:
        reqs = f.read().decode("utf-8")
    return reqs

packages = [
    "molmimic",
    "molmimic.pdb_tools",
    "molmimic.common",
    "molmimic.generate_data",
    "molmimic.parsers",
    "molmimic.parsers.CNS-Templates",
    "molmimic.parsers.superpose",
    "molmimic.torch_model",
    "molmimic.util",
    "molmimic.visualize",
#    "prodigy",
#    "prodigy.lib",
#    "prodigy.data"
]

setup(
    name = "molmimic",
    version = "0.0.1",
    author = "Eli Draizen",
    author_email = "edraizen@gmail.com",
    packages=packages,
    long_description=read('README.md'),
    install_requires=read("requirements.txt").splitlines(),
    include_package_data=True
)
