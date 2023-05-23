from setuptools import setup

setup(
    name='LOggBe',
    version='1.0.0',
    description='Bacterial core gene evolution analysis',
    install_requires=[
        'iqtree',  # Specify the desired version of iqtree
        'unifrac',
        'biom-format'
    ],
    packages=['LOggBe'],
)
