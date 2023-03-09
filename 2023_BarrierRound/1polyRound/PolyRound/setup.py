# ©2020-​2021 ETH Zurich, Axel Theorell

from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="PolyRound",
    version="0.1.8",
    description="A python package for rounding polytopes.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitlab.com/csb.ethz/PolyRound",
    author="Axel Theorell",
    author_email="atheorell@ethz.ch",
    license="MIT License",
    packages=[
        "PolyRound",
        "PolyRound.mutable_classes",
        "PolyRound.static_classes",
        "PolyRound.static_classes.rounding",
    ],
    install_requires=[
        "numpy>=1.2",
        "pandas>=1.2",
        "python-dateutil>=2.8",
        "python-libsbml>=5.18",
        "scipy>=1.4",
        "h5py>=2.10",
        "optlang>=1.4",
        "tables>=3.6",
        "cobra>=0.20",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
    python_requires=">=3.6",
)
