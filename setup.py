from setuptools import setup, find_packages

setup(
    name="trajprocess",
    version='0.25',
    packages=find_packages(),
    requires=['numpy', 'mdtraj', 'nose'],
    zip_safe=False,
    package_data={
        'trajprocess.tests': ['topol.tpr']
    },
)
