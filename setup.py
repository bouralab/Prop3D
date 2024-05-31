import os
import warnings
from setuptools import setup, find_packages

def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname), "rb") as f:
        reqs = f.read().decode("utf-8")
    return reqs

from pkg_resources import require, DistributionNotFound, parse_version

def check_provided(distribution, min_version, max_version=None, optional=False):
    # taken from https://github.com/BD2KGenomics/toil-scripts/blob/master/setup.py
    min_version = parse_version(min_version)
    if isinstance(min_version, tuple):
        raise RuntimeError("Setuptools version 8.0 or newer required. Update by running "
                           "'pip install setuptools --upgrade'")
    if max_version is not None:
        max_version = parse_version(max_version)

    messages = []

    toil_missing = 'Cannot find a valid installation of Toil.'
    dist_missing = 'Cannot find an installed copy of the %s distribution, typically provided by Toil.' % distribution
    version_too_low = 'The installed copy of %s is out of date.' % distribution
    version_too_high = 'The installed copy of %s is too new.' % distribution
    required_version = 'Setup requires version %s or higher' % min_version
    required_version += '.' if max_version is None else ', up to but not including %s.' % max_version
    install_toil = 'Installing Toil should fix this problem.'
    upgrade_toil = 'Upgrading Toil should fix this problem.'
    reinstall_dist = 'Uninstalling %s and reinstalling Toil should fix this problem.' % distribution
    reinstall_toil = 'Uninstalling Toil and reinstalling it should fix this problem.'
    footer = ("Setup doesn't install Toil automatically to give you a chance to choose any of the optional extras "
              "that Toil provides. More on installing Toil at http://toil.readthedocs.io/en/latest/installation.html.")
    try:
        # This check will fail if the distribution or any of its dependencies are missing.
        installed_version = parse_version(require(distribution)[0].version)
    except DistributionNotFound:
        installed_version = None
        if not optional:
            messages.extend([toil_missing if distribution == 'toil' else dist_missing, install_toil])
    else:
        if installed_version < min_version:
            messages.extend([version_too_low, required_version,
                             upgrade_toil if distribution == 'toil' else reinstall_dist])
        elif max_version is not None and max_version < installed_version:
            messages.extend([version_too_high, required_version,
                             reinstall_toil if distribution == 'toil' else reinstall_dist])
    if messages:
        messages.append(footer)
        raise RuntimeError(' '.join(messages))
    else:
        return str(installed_version)

try:
    toil_version = check_provided('toil', min_version='3.7.0a1.dev392', max_version='3.7.0a1.dev392')
except RuntimeError as e:
    warnings.warn(f"Toil is not installed, you will not ne able to create datasets with AtomicToil. \n\n {e}")

packages = [
    "Prop3D",
    "Prop3D.pdb_tools",
    "Prop3D.common",
    "Prop3D.custom_featurizers",
    "Prop3D.generate_data",
    "Prop3D.parsers",
    "Prop3D.parsers.CNS-Templates",
    "Prop3D.parsers.superpose",
    "Prop3D.util",
    "Prop3D.visualize",
    "Prop3D.ml",
    "Prop3D.ml.datasets"
]

setup(
    name = "Prop3D",
    version = "0.0.1",
    author = "Eli Draizen",
    author_email = "edraizen@gmail.com",
    packages=packages,
    long_description=read('README.md'),
    install_requires=read("requirements.txt").splitlines(),
    extras_require={"s3": ["boto3", "awscli"]},
    include_package_data=True,
    zip_safe=False
)
