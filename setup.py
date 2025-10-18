from setuptools import setup, find_packages

setup(
    name="MetaMAG",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "requests",
        "scikit-learn",
        "plotly",
        "pyyaml",
        "numpy",
        "pandas",
        "matplotlib",
        "seaborn",
        "biopython",
        "beautifulsoup4",
        # Add any other Python dependencies required
    ],
    entry_points={
        'console_scripts': [
            'metamag=MetaMAG.main:main'  # This creates a 'metamag' command
        ]
    },
    author="Maria Satti",
    description="A metagenomics pipeline for sequencing data analysis",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/MetaMAG",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
