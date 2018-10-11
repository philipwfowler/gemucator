from setuptools import setup

setup(
    name='gemucator',
    version='1.0.1',
    author='Philip W Fowler',
    packages=['gemucator'],
    package_data={'gemucator': ['../config/H37Rv.gbk']},
    install_requires=[
        "biopython >= 1.70"
    ],
    scripts=["bin/gemucator-run.py"],
    license='unknown',
    long_description=open('README.md').read(),
    zip_safe=False
)
