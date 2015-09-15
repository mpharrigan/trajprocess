from setuptools import setup, find_packages

setup(
    name="trajprocess",
    version='0.17',
    packages=find_packages(),
    requires=['numpy', 'mdtraj', 'nose'],
)
