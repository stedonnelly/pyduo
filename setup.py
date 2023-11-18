from setuptools import setup, find_packages

setup(
    name='pyduo',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',  # Required dependency
        'pandas',  # Required dependency
        'scipy',  # Required dependency
        'mendeleev',  # Required dependency
        'colorlog',
    ],
    # Additional metadata about your package
    author='Shaun Donnelly',
    author_email='s.t.donnelly@sheffield.ac.uk',
    description='A package designed to interface with the DUO quantum chemistry software developed by Yurchenko et. al.',
    # Include other relevant information
)
