import os
import codecs
from setuptools import setup

here = os.path.abspath(os.path.dirname(__file__))
with codecs.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='tsconvert',
    description='Tree sequence conversion utilities',
    long_description=long_description,
    url='https://github.com/tskit-dev/tsconvert',
    author='Tskit Developers',
    # TODO setup a tskit developers email address.
    author_email='jerome.kelleher@well.ox.ac.uk',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
    ],
    keywords='tree sequence, newick, nexus',
    packages=['tsconvert'],
    include_package_data=True,
    install_requires=["tskit", "dendropy"],
    entry_points={
        'console_scripts': [
            'tsconvert=tsconvert.__main__:main',
        ],
    },
    project_urls={
        'Bug Reports': 'https://github.com/tskit-dev/tsconvert/issues',
        'Source': 'https://github.com/tsckit-dev/tsconvert',
    },
    license="MIT",
    platforms=["POSIX", "MacOS X"],
    setup_requires=['setuptools_scm'],
    use_scm_version={"write_to": "tsconvert/_version.py"},
)
