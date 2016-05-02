#! /usr/bin/env python

descr = """The 2016 final project for UC Berkeley's Statistics MA Capstone Project
           Spatio-temporal analysis of mouse experiments
        """

DISTNAME = 'mousestyles'
DESCRIPTION = 'Spatio-temporal analysis of mouse experiments'
LONG_DESCRIPTION = descr
MAINTAINER = 'UC Berkeley SPRING 2016 STAT222 class'
URL = 'https://github.com/berkeley-stat222/mousestyles'
LICENSE = 'BSD License'
DOWNLOAD_URL = 'https://github.com/berkeley-stat222/mousestyles'
VERSION = '0.1dev'

INSTALL_REQUIRES = [
                    'numpy>=1.10.1',
                    'pandas>=0.17.0',
                    'pytest>=2.9.1',
                    'scipy>=0.16.1',
                    'scikit_learn>=0.17',
                    'statsmodels>=0.6.1',
                    'Sphinx>=1.4.1',
                    'sphinxcontrib-bibtex>=0.3.3',
                    'numpydoc>=0.6.0',
                    'ghp-import>=0.4.1',
                    ]

TESTS_REQUIRE = [
                'coverage>=4.0.3',
                'coveralls>=1.1',
                'pytest-cov>=2.2.1',
                'pytest-flakes>=1.0.1',
                'pytest-pep8>=1.0.6',
                ]


if __name__ == "__main__":
    from setuptools import setup
    setup(
          name=DISTNAME,
          version=VERSION,
          license=LICENSE,
          description=DESCRIPTION,
          long_description=LONG_DESCRIPTION,
          maintainer=MAINTAINER,
          url=URL,
          download_url=DOWNLOAD_URL,

          classifiers=[
                      'Development Status :: 3 - Alpha',
                      'Environment :: Console',
                      'Intended Audience :: Developers',
                      'Intended Audience :: Science/Research',
                      'License :: OSI Approved :: BSD License',
                      'Programming Language :: Python',
                      'Programming Language :: Python :: 2.7',
                      'Programming Language :: Python :: 3.4',
                      'Programming Language :: Python :: 3.5',
                      'Topic :: Scientific/Engineering',
                      'Operating System :: Microsoft :: Windows',
                      'Operating System :: POSIX',
                      'Operating System :: Unix',
                      'Operating System :: MacOS',
                      ],

          install_requires=INSTALL_REQUIRES,
          tests_require=TESTS_REQUIRE,

          packages=["mousestyles", "mousestyles.tests",
              "mousestyles.data", "mousestyles.data.tests",
              "mousestyles.classification", "mousestyles.classification.tests",
              "mousestyles.path_diversity", "mousestyles.path_diversity.tests",
              "mousestyles.visualization",
			        "mousestyles.ultradian", "mousestyles.ultradian.tests"],
          package_data={'mousestyles.data': ['*.npy', '*/*/*.npy']}
          )
