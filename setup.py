#! usr/bin/env python

"""
install with `pip install -e .`
"""

from setuptools import setup

setup(
    name="hogtie",
    version="0.0.3",
    packages=[],
    entry_points={
        'console_scripts': ['hogtie = hogtie.__main__:main']
    }
)
