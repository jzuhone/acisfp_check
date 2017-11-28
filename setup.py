#!/usr/bin/env python
from setuptools import setup
from acisfp_check import __version__

entry_points = {'console_scripts': 'acisfp_check = acisfp_check.acisfp_check:main'}

url = 'https://github.com/acisops/acisfp_check/tarball/{}'.format(__version__)

setup(name='acisfp_check',
      packages=["acisfp_check"],
      version=__version__,
      description='ACIS Thermal Model for FPTEMP',
      author='John ZuHone',
      author_email='jzuhone@gmail.com',
      url='http://github.com/acisops/acisfp_check',
      download_url=url,
      include_package_data=True,
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.5',
      ],
      entry_points=entry_points,
      )
