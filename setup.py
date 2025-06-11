from setuptools import setup, find_packages

setup(
    name="aastk",
    version="0.0.1",
    description="Automated Annotation Support ToolKit",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Marcel Rennig and Daan Speth",
    author_email="daan.speth@univie.ac.at",
    url="https://github.com/dspeth/aastk",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "contourpy>=1.3.2",
        "matplotlib>=3.10.1",
        "numpy>=2.2.5", 
        "openTSNE>=1.0.2",
        "pandas>=2.2.3",
        "PyYAML>=6.0.2",
        "scikit-learn>=1.6.1",
        "scipy>=1.15.3"
    ],
    entry_points={
        "console_scripts": [
            "aastk=aastk.__main__:main",
        ],
    },
    python_requires=">=3.12",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
