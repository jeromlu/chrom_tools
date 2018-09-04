# -*- coding: utf-8 -*-
"""
Created on Sat Aug 25 07:44:11 2018

@author: JEROMLU2
"""

from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='chrom_tools',
      version='0.1',
      description='Some classes and functions for chromatography simulatons',
      long_description =readme(),
      url='no current URL',
      author='Lux',
      author_email='luka.jeromel1@gmail.com',
      license='open license',
      packages=['chrom_tools'],
      install_requires=['numpy','pandas','matplotlib',],
      include_package_data=True,
      zip_safe=False)