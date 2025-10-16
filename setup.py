# setup.py (create this new file)
from setuptools import setup, find_packages

setup(
    name="metamag",
    version="1.0.0",
    packages=find_packages(),
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'metamag=MetaMAG.cli:main',
        ],
    },
    install_requires=[
        'click',
        'pyyaml',
        'pandas',
        'numpy',
        'biopython',
    ],
)
