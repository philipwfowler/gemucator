from setuptools import setup

setup(
    name='tbparse',
    version='1.0.0',
    author='Philip W Fowler',
    packages=['tbparse'],
    install_requires=[
        "biopython >= 1.70"
    ],
    scripts=["bin/tbparse-run.py"],
    license='unknown',
    long_description=open('README.md').read(),
)
