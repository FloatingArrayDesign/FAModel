import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="famodel",
    version="0.1.0",
    author="National Renewable Energy Laboratory",
    author_email="matthew.hall@nrel.gov",
    description="A floating array modeling library",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/FloatingArrayDesign/FAModel",
    packages=['famodel',
              'famodel.mooring',
              'famodel.anchors',
              'famodel.anchors.anchors_famodel',
              'famodel.platform',
              'famodel.cables',
              'famodel.turbine',
              'famodel.substation',
              'famodel.failure',
              'famodel.seabed'],
    package_data={'famodel.cables':['cableProps_default.yaml']},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: LGPL",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        "numpy",
        "matplotlib",
        "pyyaml"
    ],
)