from setuptools import setup

setup(
    name="tsconvert",
    description="Tree sequence conversion utilities",
    long_description="Convert various file formats to and from tskit tree sequences",
    url="https://github.com/tskit-dev/tsconvert",
    author="Tskit Developers",
    python_requires=">=3.7",
    author_email="admin@tskit.dev",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Developers",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    keywords="tree sequence, newick, nexus",
    packages=["tsconvert"],
    include_package_data=True,
    install_requires=["tskit", "dendropy"],
    entry_points={
        "console_scripts": [
            "tsconvert=tsconvert.__main__:main",
        ],
    },
    project_urls={
        "Bug Reports": "https://github.com/tskit-dev/tsconvert/issues",
        "Source": "https://github.com/tsckit-dev/tsconvert",
    },
    license="MIT",
    platforms=["POSIX", "MacOS X"],
    setup_requires=["setuptools_scm"],
    use_scm_version={"write_to": "tsconvert/_version.py"},
)
