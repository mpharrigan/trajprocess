from setuptools import setup, find_packages

setup(
    name="trajprocess",
    version='2.0.3',
    packages=find_packages(),
    requires=['numpy', 'mdtraj', 'nose'],
    zip_safe=False,
    include_package_data=True,
)
