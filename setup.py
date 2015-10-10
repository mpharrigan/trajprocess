from setuptools import setup, find_packages

setup(
    name="trajprocess",
    version='0.32',
    packages=find_packages(),
    requires=['numpy', 'mdtraj', 'nose'],
    zip_safe=False,
    include_package_data=True,
)
