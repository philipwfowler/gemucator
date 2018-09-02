from setuptools import setup

setup(
    name='genucator',
    version='1.0.0',
    author='Philip W Fowler',
    packages=['genucator'],
    install_requires=[
        "biopython >= 1.70"
    ],
    scripts=["bin/genucator-run.py"],
    license='unknown',
    long_description=open('README.md').read(),
)
