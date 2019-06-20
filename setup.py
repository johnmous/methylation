from setuptools import setup, find_packages

setup(
    name="methylator",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "methylator = methylator.splitReads:topLevel"
        ]
    }
)
