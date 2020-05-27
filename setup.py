import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="adsorption_isotherm_fitting",
    version="0.0.4",
    author="Robert F. DeJaco",
    author_email="dejac001@umn.edu",
    description="Fitting Temperature-dependent isotherms",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dejac001/adsorption_isotherm_fitting",
    packages=['isotherm_models'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU Lesser General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "matplotlib==3.2.1", "numpy==1.18.4", "scipy==1.4.1",
        "chem-util==0.0.4", "pandas==1.0.3", "realgas==1.0.2",
        "Pyomo==5.6.9"
    ]
)