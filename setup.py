import os
from setuptools import setup, find_packages

def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname), "rb") as f:
        reqs = f.read().decode("utf-8")
    return reqs

packages = [
    "Prop3D",
    "Prop3D.pdb_tools",
    "Prop3D.common",
    "Prop3D.generate_data",
    "Prop3D.parsers",
    "Prop3D.parsers.CNS-Templates",
    "Prop3D.parsers.superpose",
    "Prop3D.util",
    "Prop3D.visualize",
    "Prop3D.ml",
    "Prop3D.ml.datasets",
    "Prop3D.ml.examples"
#    "prodigy",
#    "prodigy.lib",
#    "prodigy.data"
]

setup(
    name = "Prop3D",
    version = "0.0.1",
    author = "Eli Draizen",
    author_email = "edraizen@gmail.com",
    packages=packages,
    long_description=read('README.md'),
    install_requires=read("requirements.txt").splitlines(),
    include_package_data=True,
    zip_safe=False
)
