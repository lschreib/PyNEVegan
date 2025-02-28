from setuptools import setup, find_packages

setup(
    name="PyNEVegan",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas"
        # Add other dependencies here
    ],
    author="Lars Schreiber",
    author_email="lars.schreib@googlemail.com",
    description="A Python implementation of R's vegan package for numerical ecology.",
    url="https://github.com/lschreib/PyNEVegan",
)