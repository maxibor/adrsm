from adrsm import __version__
from setuptools import setup, find_packages

setup(
    name="adrsm",
    version=__version__,
    description="Ancient DNA Read Simulator for Metagenomic ",
    long_description=open("README.md").read(),
    url="https://github.com/maxibor/adrsm",
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.6",
    install_requires=[
        "numpy >=1.18.1",
        "scipy >= 1.3.1",
        "requests >= 2.22",
        "click",
        "xopen",
        "tqdm",
    ],
    packages=find_packages(include=["adrsm", "adrsm.lib"]),
    package_data={"adrsm": ["data/quality/fwd_qual.p", "data/quality/rev_qual.p"]},
    entry_points={"console_scripts": ["adrsm= adrsm.adrsm:cli"]},
)